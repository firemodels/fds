MODULE PRES
 
! Find the perturbation pressure by solving Poisson's Equation
 
USE PRECISION_PARAMETERS
USE MESH_POINTERS

IMPLICIT NONE

PRIVATE
CHARACTER(255), PARAMETER :: presid='$Id$'
CHARACTER(255), PARAMETER :: presrev='$Revision$'
CHARACTER(255), PARAMETER :: presdate='$Date$'
PUBLIC PRESSURE_SOLVER,COMPUTE_A_B,COMPUTE_CORRECTION_PRESSURE,GET_REV_PRES
 
CONTAINS
 
SUBROUTINE PRESSURE_SOLVER(T,NM)
USE POIS, ONLY: H3CZSS,H2CZSS,H2CYSS,H3CSSS
USE COMP_FUNCTIONS, ONLY: SECOND
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE GLOBAL_CONSTANTS
 
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
REAL(EB), POINTER, DIMENSION(:) :: UWP
INTEGER :: I,J,K,IW,IOR,BC_TYPE,NOM,N_INT_CELLS,IIO,JJO,KKO
REAL(EB) :: TRM1,TRM2,TRM3,TRM4,RES,LHSS,RHSS,H_OTHER,DWWDT,DVVDT,DUUDT,HQ2,RFODT,U2,V2,W2,HFAC,H0RR(6),TNOW,DUMMY=0._EB, &
            TSI,TIME_RAMP_FACTOR,H_EXTERNAL
TYPE (VENTS_TYPE), POINTER :: VT
 
IF (SOLID_PHASE_ONLY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   UWP=> UWS
ELSE
   UU => US
   VV => VS
   WW => WS
   UWP=> UW
ENDIF

! Miscellaneous settings for wind and baroclinic cases
 
RFODT = RELAXATION_FACTOR/DT
HFAC  = 1._EB-RHOA/RHO_AVG
H0RR  = 0._EB
IF (U0>=0._EB) H0RR(1) = H0*RHOA/RHO_AVG
IF (U0<=0._EB) H0RR(2) = H0*RHOA/RHO_AVG
IF (V0>=0._EB) H0RR(3) = H0*RHOA/RHO_AVG
IF (V0<=0._EB) H0RR(4) = H0*RHOA/RHO_AVG
IF (W0>=0._EB) H0RR(5) = H0*RHOA/RHO_AVG
IF (W0<=0._EB) H0RR(6) = H0*RHOA/RHO_AVG
IF (EVACUATION_ONLY(NM)) H0RR(1:6) = 0._EB
 
! Apply pressure boundary conditions at external cells.
! If Neumann, BXS, BXF, etc., contain dH/dx(x=XS), dH/dx(x=XF), etc.
! If Dirichlet, BXS, BXF, etc., contain H(x=XS), H(x=XF), etc.
! LBC, MBC and NBC are codes used be Poisson solver to denote type
! of boundary condition at x, y and z boundaries. See Crayfishpak
! manual for details.
 
WALL_CELL_LOOP: DO IW=1,NEWC
 
   I   = IJKW(1,IW)
   J   = IJKW(2,IW)
   K   = IJKW(3,IW)
   IOR = IJKW(4,IW)

   ! Identify the type of pressure BC for each of the six mesh boundaries

   SELECT CASE(IOR)
      CASE( 1)
         IF (LBC==3 .OR. LBC==4 .OR. LBC==6) BC_TYPE = NEUMANN
         IF (LBC==1 .OR. LBC==2 .OR. LBC==5) BC_TYPE = DIRICHLET
      CASE(-1)
         IF (LBC==2 .OR. LBC==3 .OR. LBC==6) BC_TYPE = NEUMANN
         IF (LBC==1 .OR. LBC==4 .OR. LBC==5) BC_TYPE = DIRICHLET
      CASE( 2)
         IF (MBC==3 .OR. MBC==4) BC_TYPE = NEUMANN
         IF (MBC==1 .OR. MBC==2) BC_TYPE = DIRICHLET
      CASE(-2)
         IF (MBC==3 .OR. MBC==2) BC_TYPE = NEUMANN
         IF (MBC==1 .OR. MBC==4) BC_TYPE = DIRICHLET
      CASE( 3)
         IF (NBC==3 .OR. NBC==4) BC_TYPE = NEUMANN
         IF (NBC==1 .OR. NBC==2) BC_TYPE = DIRICHLET
      CASE(-3)
         IF (NBC==3 .OR. NBC==2) BC_TYPE = NEUMANN
         IF (NBC==1 .OR. NBC==4) BC_TYPE = DIRICHLET
   END SELECT

   ! Apply pressure gradients at NEUMANN boundaries: dH/dn = -F_n - d(u_n)/dt

   IF_NEUMANN: IF (BC_TYPE==NEUMANN) THEN
      SELECT CASE(IOR)
         CASE( 1)
            BXS(J,K) = HX(0)   *(-FVX(0,J,K)    + DUWDT(IW))
         CASE(-1)
            BXF(J,K) = HX(IBP1)*(-FVX(IBAR,J,K) - DUWDT(IW))
         CASE( 2)
            BYS(I,K) = HY(0)   *(-FVY(I,0,K)    + DUWDT(IW))
         CASE(-2)
            BYF(I,K) = HY(JBP1)*(-FVY(I,JBAR,K) - DUWDT(IW))
         CASE( 3)
            BZS(I,J) = HZ(0)   *(-FVZ(I,J,0)    + DUWDT(IW))
         CASE(-3)
            BZF(I,J) = HZ(KBP1)*(-FVZ(I,J,KBAR) - DUWDT(IW))
      END SELECT
   ENDIF IF_NEUMANN

   ! Apply pressures at DIRICHLET boundaries, depending on the specific type
 
   IF_DIRICHLET: IF (BC_TYPE==DIRICHLET) THEN

      NOT_OPEN: IF (BOUNDARY_TYPE(IW)/=OPEN_BOUNDARY .AND. BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN

         ! Solid boundary that uses a Dirichlet BC to drive the normal component of velocity towards UWP
 
         SELECT CASE(IOR)
            CASE( 1)
               DUUDT = -RFODT*(UU(0,J,K)   +UWP(IW))
               BXS(J,K) = H(1,J,K)     + 0.5_EB*DX(0)   *(DUUDT+FVX(0,J,K))
            CASE(-1) 
               DUUDT = -RFODT*(UU(IBAR,J,K)-UWP(IW))
               BXF(J,K) = H(IBAR,J,K) - 0.5_EB*DX(IBP1)*(DUUDT+FVX(IBAR,J,K))
            CASE( 2)
               DVVDT = -RFODT*(VV(I,0,K)   +UWP(IW)) 
               BYS(I,K) = H(I,1,K)    + 0.5_EB*DY(0)   *(DVVDT+FVY(I,0,K))
            CASE(-2) 
               DVVDT = -RFODT*(VV(I,JBAR,K)-UWP(IW))
               BYF(I,K) = H(I,JBAR,K) - 0.5_EB*DY(JBP1)*(DVVDT+FVY(I,JBAR,K))
            CASE( 3)
               DWWDT = -RFODT*(WW(I,J,0)   +UWP(IW))
               BZS(I,J) = H(I,J,1)    + 0.5_EB*DZ(0)   *(DWWDT+FVZ(I,J,0))
            CASE(-3) 
               DWWDT = -RFODT*(WW(I,J,KBAR)-UWP(IW))
               BZF(I,J) = H(I,J,KBAR) - 0.5_EB*DZ(KBP1)*(DWWDT+FVZ(I,J,KBAR))
         END SELECT

      ENDIF NOT_OPEN

      ! Interpolated boundary -- set boundary value of H to be average of neighboring cells from previous time step
 
      INTERPOLATED_ONLY: IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) THEN

         NOM     = IJKW(9,IW)
         H_OTHER = 0._EB
         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  H_OTHER = H_OTHER + OMESH(NOM)%H(IIO,JJO,KKO)
               ENDDO
            ENDDO
         ENDDO
         N_INT_CELLS = (IJKW(13,IW)-IJKW(10,IW)+1) * (IJKW(14,IW)-IJKW(11,IW)+1) * (IJKW(15,IW)-IJKW(12,IW)+1)
         H_OTHER = H_OTHER/REAL(N_INT_CELLS,EB)

         SELECT CASE(IOR)
            CASE( 1)
                  BXS(J,K) = 0.5*(H(1,J,K)+H_OTHER)
            CASE(-1) 
                  BXF(J,K) = 0.5*(H(IBAR,J,K)+H_OTHER) 
            CASE( 2) 
                  BYS(I,K) = 0.5*(H(I,1,K)+H_OTHER) 
            CASE(-2) 
                  BYF(I,K) = 0.5*(H(I,JBAR,K)+H_OTHER) 
            CASE( 3) 
                  BZS(I,J) = 0.5*(H(I,J,1)+H_OTHER)
            CASE(-3) 
                  BZF(I,J) = 0.5*(H(I,J,KBAR)+H_OTHER) 
         END SELECT
 
      ENDIF INTERPOLATED_ONLY
 
      ! OPEN (passive opening to exterior of domain) boundary. Apply inflow/outflow BC.

      OPEN: IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
 
         H_EXTERNAL = 0._EB
         IF (VENT_INDEX(IW)>0) THEN
            VT => VENTS(VENT_INDEX(IW))
            TSI = T - T_BEGIN
            TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,DUMMY,VT%PRESSURE_RAMP_INDEX)
            H_EXTERNAL = TIME_RAMP_FACTOR*VT%DYNAMIC_PRESSURE/RHOA
         ENDIF

         SELECT CASE(IOR)
            CASE( 1)
               U2  = UU(0,J,K)**2
               V2  = .25_EB*(VV(1,J,K)+VV(1,J-1,K))**2
               W2  = .25_EB*(WW(1,J,K)+WW(1,J,K-1))**2
               HQ2 = MIN(5000._EB,0.5_EB*(U2+V2+W2))
               IF (UU(0,J,K)<0._EB) THEN
                  BXS(J,K) = H_EXTERNAL + HQ2
               ELSE
                  BXS(J,K) = H_EXTERNAL + H0RR(1) + HQ2*HFAC
               ENDIF
            CASE(-1)
               U2  = UU(IBAR,J,K)**2
               V2  = .25_EB*(VV(IBAR,J,K)+VV(IBAR,J-1,K))**2
               W2  = .25_EB*(WW(IBAR,J,K)+WW(IBAR,J,K-1))**2
               HQ2 = MIN(5000._EB,0.5_EB*(U2+V2+W2))
               IF (UU(IBAR,J,K)>0._EB) THEN
                  BXF(J,K) = H_EXTERNAL + HQ2
               ELSE
                  BXF(J,K) = H_EXTERNAL + H0RR(2) + HQ2*HFAC
               ENDIF
            CASE( 2)
               U2  = .25_EB*(UU(I,1,K)+UU(I-1,1,K))**2
               V2  = VV(I,0,K)**2
               W2  = .25_EB*(WW(I,1,K)+WW(I,1,K-1))**2
               HQ2 = MIN(5000._EB,0.5_EB*(U2+V2+W2))
               IF (VV(I,0,K)<0._EB) THEN
                  BYS(I,K) = H_EXTERNAL + HQ2
               ELSE
                  BYS(I,K) = H_EXTERNAL + H0RR(3) + HQ2*HFAC
               ENDIF
            CASE(-2)
               U2  = .25_EB*(UU(I,JBAR,K)+UU(I-1,JBAR,K))**2
               V2  = VV(I,JBAR,K)**2
               W2  = .25_EB*(WW(I,JBAR,K)+WW(I,JBAR,K-1))**2
               HQ2 = MIN(5000._EB,0.5_EB*(U2+V2+W2))
               IF (VV(I,JBAR,K)>0._EB) THEN
                  BYF(I,K) = H_EXTERNAL + HQ2
               ELSE
                  BYF(I,K) = H_EXTERNAL + H0RR(4) + HQ2*HFAC
               ENDIF
            CASE( 3)
               U2  = .25_EB*(UU(I,J,1)+UU(I-1,J,1))**2
               V2  = .25_EB*(VV(I,J,1)+VV(I,J-1,1))**2
               W2  = WW(I,J,0)**2
               HQ2 = MIN(5000._EB,0.5_EB*(U2+V2+W2))
               IF (WW(I,J,0)<0._EB) THEN
                  BZS(I,J) = H_EXTERNAL + HQ2
               ELSE
                  BZS(I,J) = H_EXTERNAL + H0RR(5) + HQ2*HFAC
               ENDIF
            CASE(-3)
               U2  = .25_EB*(UU(I,J,KBAR)+UU(I-1,J,KBAR))**2
               V2  = .25_EB*(VV(I,J,KBAR)+VV(I,J-1,KBAR))**2
               W2  = WW(I,J,KBAR)**2
               HQ2 = MIN(5000._EB,0.5_EB*(U2+V2+W2))
               IF (WW(I,J,KBAR)>0._EB) THEN
                  BZF(I,J) = H_EXTERNAL + HQ2
               ELSE
                  BZF(I,J) = H_EXTERNAL + H0RR(6) + HQ2*HFAC
               ENDIF
         END SELECT
    
      ENDIF OPEN
 
   ENDIF IF_DIRICHLET
 
ENDDO WALL_CELL_LOOP

! Compute the RHS of the Poisson equation
 
IF (CYLINDRICAL) THEN
   DO K=1,KBAR
      DO I=1,IBAR
         TRM1 = (R(I-1)*FVX(I-1,1,K)-R(I)*FVX(I,1,K))*RDX(I)*RRN(I)
         TRM3 = (FVZ(I,1,K-1)-FVZ(I,1,K))*RDZ(K)
         TRM4 = -DDDT(I,1,K)
         PRHS(I,1,K) = TRM1 + TRM3 + TRM4
      ENDDO
   ENDDO
ENDIF
 
IF ( (IPS<=1 .OR. IPS==4 .OR. IPS==7) .AND. .NOT.CYLINDRICAL) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            TRM1 = (FVX(I-1,J,K)-FVX(I,J,K))*RDX(I)
            TRM2 = (FVY(I,J-1,K)-FVY(I,J,K))*RDY(J) 
            TRM3 = (FVZ(I,J,K-1)-FVZ(I,J,K))*RDZ(K)
            TRM4 = -DDDT(I,J,K)
            PRHS(I,J,K) = TRM1 + TRM2 + TRM3 + TRM4
         ENDDO
      ENDDO
   ENDDO
END IF
 
IF (IPS==2) THEN  ! Switch x and y
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            TRM1 = (FVX(I-1,J,K)-FVX(I,J,K))*RDX(I)
            TRM2 = (FVY(I,J-1,K)-FVY(I,J,K))*RDY(J)
            TRM3 = (FVZ(I,J,K-1)-FVZ(I,J,K))*RDZ(K)
            TRM4 = -DDDT(I,J,K)
            PRHS(J,I,K) = TRM1 + TRM2 + TRM3 + TRM4
         ENDDO
      ENDDO
   ENDDO
   BZST = TRANSPOSE(BZS)
   BZFT = TRANSPOSE(BZF)
ENDIF
 
IF (IPS==3 .OR. IPS==6) THEN  ! Switch x and z
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            TRM1 = (FVX(I-1,J,K)-FVX(I,J,K))*RDX(I)
            TRM2 = (FVY(I,J-1,K)-FVY(I,J,K))*RDY(J)
            TRM3 = (FVZ(I,J,K-1)-FVZ(I,J,K))*RDZ(K)
            TRM4 = -DDDT(I,J,K)
            PRHS(K,J,I) = TRM1 + TRM2 + TRM3 + TRM4
         ENDDO
      ENDDO
   ENDDO
   BXST = TRANSPOSE(BXS)
   BXFT = TRANSPOSE(BXF)
   BYST = TRANSPOSE(BYS)
   BYFT = TRANSPOSE(BYF)
   BZST = TRANSPOSE(BZS)
   BZFT = TRANSPOSE(BZF)
END IF
 
IF (IPS==5) THEN  ! Switch y and z
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            TRM1 = (FVX(I-1,J,K)-FVX(I,J,K))*RDX(I)
            TRM2 = (FVY(I,J-1,K)-FVY(I,J,K))*RDY(J)
            TRM3 = (FVZ(I,J,K-1)-FVZ(I,J,K))*RDZ(K)
            TRM4 = -DDDT(I,J,K)
            PRHS(I,K,J) = TRM1 + TRM2 + TRM3 + TRM4
         ENDDO
      ENDDO
   ENDDO
   BXST = TRANSPOSE(BXS)
   BXFT = TRANSPOSE(BXF)
END IF

! Poisson solve with stretching the 1st coordinate
 
IF (IPS<=1) THEN
   IF (.NOT.TWO_D) CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE,WORK,HX)
   IF (TWO_D .AND. .NOT. CYLINDRICAL) CALL H2CZSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE,WORK,HX)
   IF (TWO_D .AND. CYLINDRICAL) CALL H2CYSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE,WORK)
ENDIF
 
IF (IPS==2) CALL H3CZSS(BYS,BYF,BXS,BXF,BZST,BZFT,ITRN,JTRN,PRHS,POIS_PTB,SAVE,WORK,HY)
 
IF (IPS==3) THEN
   IF (.NOT.TWO_D) CALL H3CZSS(BZST,BZFT,BYST,BYFT,BXST,BXFT,ITRN,JTRN,PRHS,POIS_PTB,SAVE,WORK,HZ)
   IF (TWO_D)CALL H2CZSS(BZS,BZF,BXS,BXF,ITRN,PRHS,POIS_PTB,SAVE,WORK,HZ)
ENDIF
 
! Poisson solve with stretching the 1st and 2nd coordinate
 
IF (IPS==4) CALL H3CSSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE,WORK,HX,HY)
 
IF (IPS==5) THEN
   IF (.NOT.TWO_D) CALL H3CSSS(BXST,BXFT,BZS,BZF,BYS,BYF,ITRN,JTRN,PRHS,POIS_PTB,SAVE,WORK,HX,HZ)
   IF (TWO_D) CALL H2CZSS(BZS,BZF,BXS,BXF,ITRN,PRHS,POIS_PTB,SAVE,WORK,HZ)
ENDIF
 
IF (IPS==6) CALL H3CSSS(BZST,BZFT,BYST,BYFT,BXST,BXFT,ITRN,JTRN,PRHS,POIS_PTB,SAVE,WORK,HZ,HY)
 
IF (IPS==7) CALL H2CZSS(BXS,BXF,BYS,BYF,ITRN,PRHS,POIS_PTB,SAVE,WORK,HX)
 
! Put output of Poisson solver into the H array

SELECT CASE(IPS)
   CASE(:1,4,7)
      FORALL(K=1:KBAR,J=1:JBAR,I=1:IBAR) H(I,J,K) = PRHS(I,J,K)
   CASE(2)  
      FORALL(K=1:KBAR,J=1:JBAR,I=1:IBAR) H(I,J,K) = PRHS(J,I,K)
   CASE(3,6)
      FORALL(K=1:KBAR,J=1:JBAR,I=1:IBAR) H(I,J,K) = PRHS(K,J,I)
   CASE(5) 
      FORALL(K=1:KBAR,J=1:JBAR,I=1:IBAR) H(I,J,K) = PRHS(I,K,J)
END SELECT 

! Apply boundary conditions to H
 
DO K=1,KBAR
   DO J=1,JBAR
      IF (LBC==3 .OR. LBC==4)             H(0,J,K)    = H(1,J,K)    - DXI*BXS(J,K)
      IF (LBC==3 .OR. LBC==2 .OR. LBC==6) H(IBP1,J,K) = H(IBAR,J,K) + DXI*BXF(J,K)
      IF (LBC==1 .OR. LBC==2)             H(0,J,K)    =-H(1,J,K)    + 2._EB*BXS(J,K)
      IF (LBC==1 .OR. LBC==4 .OR. LBC==5) H(IBP1,J,K) =-H(IBAR,J,K) + 2._EB*BXF(J,K)
      IF (LBC==5 .OR. LBC==6)             H(0,J,K) = H(1,J,K)
   ENDDO
ENDDO
 
DO K=1,KBAR
   DO I=1,IBAR
      IF (MBC==3 .OR. MBC==4) H(I,0,K)    = H(I,1,K)    - DETA*BYS(I,K)
      IF (MBC==3 .OR. MBC==2) H(I,JBP1,K) = H(I,JBAR,K) + DETA*BYF(I,K)
      IF (MBC==1 .OR. MBC==2) H(I,0,K)    =-H(I,1,K)    + 2._EB*BYS(I,K)
      IF (MBC==1 .OR. MBC==4) H(I,JBP1,K) =-H(I,JBAR,K) + 2._EB*BYF(I,K)
   ENDDO
ENDDO
 
DO J=1,JBAR
   DO I=1,IBAR
      IF (NBC==3 .OR. NBC==4)  H(I,J,0)    = H(I,J,1)    - DZETA*BZS(I,J)
      IF (NBC==3 .OR. NBC==2)  H(I,J,KBP1) = H(I,J,KBAR) + DZETA*BZF(I,J)
      IF (NBC==1 .OR. NBC==2)  H(I,J,0)    =-H(I,J,1)    + 2._EB*BZS(I,J)
      IF (NBC==1 .OR. NBC==4)  H(I,J,KBP1) =-H(I,J,KBAR) + 2._EB*BZF(I,J)
   ENDDO
ENDDO

! ************************* Check the Solution *************************
 
IF (CHECK_POISSON) THEN     
   POIS_ERR = 0.
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHSS = (R(I-1)*FVX(I-1,J,K)-R(I)*FVX(I,J,K))*RDX(I)*RRN(I) + (FVY(I,J-1,K)-FVY(I,J,K))*RDY(J) + &
                   (FVZ(I,J,K-1)-FVZ(I,J,K))*RDZ(K) - DDDT(I,J,K)
            LHSS = ((H(I+1,J,K)-H(I,J,K))*RDXN(I)*R(I) - (H(I,J,K)-H(I-1,J,K))*RDXN(I-1)*R(I-1) )*RDX(I)*RRN(I) &
                   +((H(I,J+1,K)-H(I,J,K))*RDYN(J) -  (H(I,J,K)-H(I,J-1,K))*RDYN(J-1) )*RDY(J) &
                   +((H(I,J,K+1)-H(I,J,K))*RDZN(K) - (H(I,J,K)-H(I,J,K-1))*RDZN(K-1) )*RDZ(K)
            RES = ABS(RHSS-LHSS)
            POIS_ERR = MAX(RES,POIS_ERR)
         ENDDO
      ENDDO
   ENDDO
ENDIF
 
! **********************************************************************
 
TUSED(5,NM)=TUSED(5,NM)+SECOND()-TNOW
END SUBROUTINE PRESSURE_SOLVER


! Everything below this point is experimental and not currently implemented
 
SUBROUTINE COMPUTE_A_B(A,B,NM)

! Set up linear system of equations for the coarse grid HBAR

USE GLOBAL_CONSTANTS, ONLY: NCGC, SOLID_BOUNDARY, OPEN_BOUNDARY, PREDICTOR
REAL(EB) :: A(NCGC,NCGC),B(NCGC),DUDT_OTHER,DVDT_OTHER,DWDT_OTHER,DA
TYPE (MESH_TYPE), POINTER :: M2
TYPE (OMESH_TYPE), POINTER :: OM
INTEGER :: I,J,K,NM,II,JJ,KK,IIO,JJO,KKO,NOM,IW,IC,JC,KC,N,NO
 
CALL POINT_TO_MESH(NM)
 
K_COARSE: DO KC=1,KBAR2
   J_COARSE: DO JC=1,JBAR2
      I_COARSE: DO IC=1,IBAR2
 
         N = CGI2(IC,JC,KC)
         II = I_LO(IC)-1
         DO KK=K_LO(KC),K_HI(KC)
            DO JJ=J_LO(JC),J_HI(JC)
               DA = DY(JJ)*DZ(KK)
               IF (II>0) THEN
                  NO = CGI(II,JJ,KK)
                  A(N,N)  = A(N,N)  - DA/DX_M(N,NO)
                  A(N,NO) = A(N,NO) + DA/DX_M(N,NO)
                  B(N)    = B(N)    + DA*FVX(II,JJ,KK)
               ELSE 
                  IW = WALL_INDEX(CELL_INDEX(II+1,JJ,KK),-1)
                  NOM = IJKW(9,IW)
                  IF (NOM>0) THEN
                     OM  => OMESH(NOM)
                     M2  => MESHES(NOM)
                     IIO = IJKW(10,IW)
                     JJO = IJKW(11,IW)
                     KKO = IJKW(12,IW)
                     NO  = M2%CGI(IIO,JJO,KKO)
                     A(N,N)  = A(N,N)  - DA/DX_M(N,NO)
                     A(N,NO) = A(N,NO) + DA/DX_M(N,NO)
                     DUDT_OTHER = 0._EB
                     DO KKO=IJKW(12,IW),IJKW(15,IW)
                        DO JJO=IJKW(11,IW),IJKW(14,IW)
                           DO IIO=IJKW(10,IW),IJKW(13,IW)
                              DUDT_OTHER = DUDT_OTHER + OM%DUDT(IIO,JJO,KKO)*MIN(1._EB,M2%DY(JJO)*M2%DZ(KKO)/DA)
                           ENDDO
                        ENDDO
                     ENDDO
                     DUDT_AVG(IW) = 0.5_EB*(DUDT_OTHER-FVX(II,JJ,KK)-RDXN(0)*(H(1,JJ,KK)-H(0,JJ,KK)))
                     B(N)         = B(N) - DA*DUDT_AVG(IW)
                  ENDIF
                  IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) B(N) = B(N) + DUWDT(IW)*DA
                  IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
                     A(N,N) = A(N,N) - 2._EB*DA/DX_M(N,N)
                     B(N)   = B(N)   + (FVX(II,JJ,KK)+(H(II+1,JJ,KK)-H(II,JJ,KK))*RDXN(0))*DA
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
 
         II = I_HI(IC)
         DO KK=K_LO(KC),K_HI(KC)
            DO JJ=J_LO(JC),J_HI(JC)
               DA = DY(JJ)*DZ(KK)
               IF (II<IBAR) THEN
                  NO = CGI(II+1,JJ,KK)
                  A(N,N)  = A(N,N)  - DA/DX_M(N,NO)
                  A(N,NO) = A(N,NO) + DA/DX_M(N,NO)
                  B(N)    = B(N)    - DA*FVX(II,JJ,KK)
               ELSE
                  IW = WALL_INDEX(CELL_INDEX(II,JJ,KK),1)
                  NOM = IJKW(9,IW)
                  IF (NOM>0) THEN
                     OM  => OMESH(NOM)
                     M2  => MESHES(NOM)
                     IIO = IJKW(10,IW)
                     JJO = IJKW(11,IW)
                     KKO = IJKW(12,IW)
                     NO  = M2%CGI(IIO,JJO,KKO)
                     A(N,N)  = A(N,N)  - DA/DX_M(N,NO)
                     A(N,NO) = A(N,NO) + DA/DX_M(N,NO)
                     DUDT_OTHER = 0._EB
                     DO KKO=IJKW(12,IW),IJKW(15,IW)
                        DO JJO=IJKW(11,IW),IJKW(14,IW)
                           DO IIO=IJKW(10,IW),IJKW(13,IW)
                              DUDT_OTHER = DUDT_OTHER + OM%DUDT(IIO-1,JJO,KKO)*MIN(1._EB,M2%DY(JJO)*M2%DZ(KKO)/DA)
                           ENDDO
                        ENDDO
                     ENDDO
                     DUDT_AVG(IW) = 0.5_EB*(DUDT_OTHER-FVX(II,JJ,KK)-RDXN(IBAR)*(H(IBP1,JJ,KK)-H(IBAR,JJ,KK)))
                     B(N)         = B(N) + DA*DUDT_AVG(IW)
                  ENDIF
                  IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) B(N) = B(N) + DUWDT(IW)*DA
                  IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
                     A(N,N) = A(N,N) - 2._EB*DA/DX_M(N,N)
                     B(N)   = B(N)   - (FVX(II,JJ,KK)+(H(II+1,JJ,KK)-H(II,JJ,KK))*RDXN(IBAR))*DA
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
 
         JJ = J_LO(JC)-1
         DO KK=K_LO(KC),K_HI(KC)
            DO II=I_LO(IC),I_HI(IC)
               DA = DX(II)*DZ(KK)
               IF (JJ>0) THEN
                  NO = CGI(II,JJ,KK)
                  A(N,N)  = A(N,N)  - DA/DY_M(N,NO) 
                  A(N,NO) = A(N,NO) + DA/DY_M(N,NO) 
                  B(N)    = B(N)    + DA*FVY(II,JJ,KK)
               ELSE
                  IW = WALL_INDEX(CELL_INDEX(II,JJ+1,KK),-2)
                  NOM = IJKW(9,IW)
                  IF (NOM>0) THEN
                     OM  => OMESH(NOM)
                     M2  => MESHES(NOM)
                     IIO = IJKW(10,IW)
                     JJO = IJKW(11,IW)
                     KKO = IJKW(12,IW)
                     NO  = M2%CGI(IIO,JJO,KKO)
                     A(N,N)  = A(N,N)  - DA/DY_M(N,NO)
                     A(N,NO) = A(N,NO) + DA/DY_M(N,NO)
                     DVDT_OTHER = 0._EB
                     DO KKO=IJKW(12,IW),IJKW(15,IW)
                        DO JJO=IJKW(11,IW),IJKW(14,IW)
                           DO IIO=IJKW(10,IW),IJKW(13,IW)
                              DVDT_OTHER = DVDT_OTHER + OM%DVDT(IIO,JJO,KKO)*MIN(1._EB,M2%DX(IIO)*M2%DZ(KKO)/DA)
                           ENDDO
                        ENDDO
                     ENDDO
                     DVDT_AVG(IW) = 0.5_EB*(DVDT_OTHER-FVY(II,JJ,KK)-RDYN(0)*(H(II,1,KK)-H(II,0,KK)))
                     B(N)         = B(N) - DA*DVDT_AVG(IW)
                  ENDIF
                  IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) B(N) = B(N) + DUWDT(IW)*DA
                  IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
                     A(N,N) = A(N,N) - 2._EB*DA/DY_M(N,N)
                     B(N)   = B(N)   + (FVY(II,JJ,KK)+(H(II,JJ+1,KK)-H(II,JJ,KK))*RDYN(0))*DA
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
 
         JJ = J_HI(JC)
         DO KK=K_LO(KC),K_HI(KC)
            DO II=I_LO(IC),I_HI(IC)
               DA = DX(II)*DZ(KK)
               IF (JJ<JBAR) THEN
                  NO = CGI(II,JJ+1,KK)
                  A(N,N)  = A(N,N)  - DA/DY_M(N,NO)
                  A(N,NO) = A(N,NO) + DA/DY_M(N,NO)
                  B(N)    = B(N)    - DA*FVY(II,JJ,KK)
               ELSE
                  IW = WALL_INDEX(CELL_INDEX(II,JJ,KK),2)
                  NOM = IJKW(9,IW)
                  IF (NOM>0) THEN
                     OM  => OMESH(NOM)
                     M2  => MESHES(NOM)
                     IIO = IJKW(10,IW)
                     JJO = IJKW(11,IW)
                     KKO = IJKW(12,IW)
                     NO  = M2%CGI(IIO,JJO,KKO)
                     A(N,N)  = A(N,N)  - DA/DY_M(N,NO)
                     A(N,NO) = A(N,NO) + DA/DY_M(N,NO)
                     DVDT_OTHER = 0._EB
                     DO KKO=IJKW(12,IW),IJKW(15,IW)
                        DO JJO=IJKW(11,IW),IJKW(14,IW)
                           DO IIO=IJKW(10,IW),IJKW(13,IW)
                              DVDT_OTHER = DVDT_OTHER + OM%DVDT(IIO,JJO-1,KKO)*MIN(1._EB,M2%DX(IIO)*M2%DZ(KKO)/DA)
                           ENDDO
                        ENDDO
                     ENDDO
                     DVDT_AVG(IW)  = 0.5_EB*(DVDT_OTHER-FVY(II,JJ,KK)-RDYN(0)*(H(II,JBP1,KK)-H(II,JBAR,KK)))
                     B(N)          = B(N) + DA*DVDT_AVG(IW)
                  ENDIF
                  IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) B(N) = B(N) + DUWDT(IW)*DA
                  IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
                     A(N,N) = A(N,N) - 2._EB*DA/DY_M(N,N)
                     B(N)   = B(N)   - (FVY(II,JJ,KK)+(H(II,JJ+1,KK)-H(II,JJ,KK))*RDYN(JBAR))*DA
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
 
         ! Lower z

         KK = K_LO(KC)-1
         DO JJ=J_LO(JC),J_HI(JC)
            DO II=I_LO(IC),I_HI(IC)
               DA = DX(II)*DY(JJ)
               IF (KK>0) THEN
                  NO = CGI(II,JJ,KK)
                  A(N,N)  = A(N,N)  - DA/DZ_M(N,NO) 
                  A(N,NO) = A(N,NO) + DA/DZ_M(N,NO) 
                  B(N)    = B(N)    + DA*FVZ(II,JJ,KK)
               ELSE
                  IW = WALL_INDEX(CELL_INDEX(II,JJ,KK+1),-3)
                  NOM = IJKW(9,IW)
                  IF (NOM>0) THEN
                     OM  => OMESH(NOM)
                     M2  => MESHES(NOM)
                     IIO = IJKW(10,IW)
                     JJO = IJKW(11,IW)
                     KKO = IJKW(12,IW)
                     NO  = M2%CGI(IIO,JJO,KKO)
                     A(N,N)  = A(N,N)  - DA/DZ_M(N,NO)
                     A(N,NO) = A(N,NO) + DA/DZ_M(N,NO) 
                     DWDT_OTHER = 0._EB
                     DO KKO=IJKW(12,IW),IJKW(15,IW)
                        DO JJO=IJKW(11,IW),IJKW(14,IW)
                           DO IIO=IJKW(10,IW),IJKW(13,IW)
                              DWDT_OTHER = DWDT_OTHER + OM%DWDT(IIO,JJO,KKO)*MIN(1._EB,M2%DX(IIO)*M2%DY(JJO)/DA)
                           ENDDO
                        ENDDO
                     ENDDO
                     DWDT_AVG(IW) = 0.5_EB*(DWDT_OTHER-FVZ(II,JJ,KK)-RDZN(0)*(H(II,JJ,1)-H(II,JJ,0)))
                     B(N)         = B(N) - DA*DWDT_AVG(IW)
                  ENDIF
                  IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) B(N) = B(N) + DUWDT(IW)*DA
                  IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
                     A(N,N) = A(N,N) - 2._EB*DA/DZ_M(N,N)
                     B(N)   = B(N)   + (FVZ(II,JJ,KK)+(H(II,JJ,KK+1)-H(II,JJ,KK))*RDZN(0))*DA
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
 
         KK = K_HI(KC)
         DO JJ=J_LO(JC),J_HI(JC)
            DO II=I_LO(IC),I_HI(IC)
               DA = DX(II)*DY(JJ)
               IF (KK<KBAR) THEN
                  NO = CGI(II,JJ,KK+1)
                  A(N,N)  = A(N,N)  - DA/DZ_M(N,NO) 
                  A(N,NO) = A(N,NO) + DA/DZ_M(N,NO) 
                  B(N)    = B(N)    - DA*FVZ(II,JJ,KK)
               ELSE
                  IW = WALL_INDEX(CELL_INDEX(II,JJ,KK),3)
                  NOM = IJKW(9,IW)
                  IF (NOM>0) THEN
                     OM  => OMESH(NOM)
                     M2  => MESHES(NOM)
                     IIO = IJKW(10,IW)
                     JJO = IJKW(11,IW)
                     KKO = IJKW(12,IW)
                     NO  = M2%CGI(IIO,JJO,KKO)
                     A(N,N)  = A(N,N)  - DA/DZ_M(N,NO) 
                     A(N,NO) = A(N,NO) + DA/DZ_M(N,NO) 
                     DWDT_OTHER = 0._EB
                     DO KKO=IJKW(12,IW),IJKW(15,IW)
                        DO JJO=IJKW(11,IW),IJKW(14,IW)
                           DO IIO=IJKW(10,IW),IJKW(13,IW)
                              DWDT_OTHER = DWDT_OTHER + OM%DWDT(IIO,JJO,KKO-1)*MIN(1._EB,M2%DX(IIO)*M2%DY(JJO)/DA)
                           ENDDO
                        ENDDO
                     ENDDO
                     DWDT_AVG(IW) = 0.5_EB*(DWDT_OTHER-FVZ(II,JJ,KK)-RDZN(KBAR)*(H(II,JJ,KBP1)-H(II,JJ,KBAR)))
                     B(N)         = B(N) + DA*DWDT_AVG(IW)
                  ENDIF
                  IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) B(N) = B(N) + DUWDT(IW)*DA
                  IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
                     A(N,N) = A(N,N) - 2._EB*DA/DZ_M(N,N)
                     B(N)   = B(N)   - (FVZ(II,JJ,KK)+(H(II,JJ,KK+1)-H(II,JJ,KK))*RDZN(KBAR))*DA
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
 
         DO K=K_LO(KC),K_HI(KC)
            DO J=J_LO(JC),J_HI(JC)
               DO I=I_LO(IC),I_HI(IC)
                  B(N) = B(N) - DDDT(I,J,K)*DX(I)*RC(I)*DY(J)*DZ(K)
               ENDDO
            ENDDO
         ENDDO
 
      ENDDO I_COARSE
   ENDDO J_COARSE
ENDDO K_COARSE
 
END SUBROUTINE COMPUTE_A_B
 

SUBROUTINE COMPUTE_CORRECTION_PRESSURE(B,NM)

! Set up and solve Laplace's Eq for the perturbation pressure HP

USE GLOBAL_CONSTANTS, ONLY: NCGC, SOLID_BOUNDARY, OPEN_BOUNDARY, PREDICTOR, TWO_D, CYLINDRICAL
USE POIS, ONLY: H3CZSS, H2CZSS
REAL(EB) :: B(NCGC),DHDX_F,DHDX_S,DHDY_F,DHDY_S,DHDZ_F,DHDZ_S,DA, &
            AREA_XS,AREA_XF,AREA_YS,AREA_YF,AREA_ZS,AREA_ZF, &
            DUUDT,DVVDT,DWWDT,RDT,B_OTHER, &
            AREA_XS_CLOSED,AREA_XF_CLOSED, &
            AREA_YS_CLOSED,AREA_YF_CLOSED, &
            AREA_ZS_CLOSED,AREA_ZF_CLOSED
INTEGER :: NM,II,JJ,KK,IW,IOR,I,J,K,IIO,JJO,KKO,NOM, &
           IOR_PATCH,II_LOW,II_HIGH,JJ_LOW,JJ_HIGH, KK_LOW,KK_HIGH,IC,JC,KC,N,NO
TYPE (MESH_TYPE), POINTER :: M2
TYPE (OMESH_TYPE), POINTER :: OM
 
CALL POINT_TO_MESH(NM)
 
RDT = 1._EB/DT
 
BXS = 0._EB
BXF = 0._EB
BYS = 0._EB
BYF = 0._EB
BZS = 0._EB
BZF = 0._EB

ORIENT_LOOP:  DO IOR_PATCH=-3,3
   IF (IOR_PATCH==0) CYCLE ORIENT_LOOP
   KC_LOOP: DO KC=1,KBAR2
      IF (IOR_PATCH== 3 .AND. KC/=1)     CYCLE KC_LOOP
      IF (IOR_PATCH==-3 .AND. KC/=KBAR2) CYCLE KC_LOOP
      JC_LOOP: DO JC=1,JBAR2
         IF (IOR_PATCH== 2 .AND. JC/=1)     CYCLE JC_LOOP
         IF (IOR_PATCH==-2 .AND. JC/=JBAR2) CYCLE JC_LOOP
         IC_LOOP: DO IC=1,IBAR2
            IF (IOR_PATCH== 1 .AND. IC/=1)     CYCLE IC_LOOP
            IF (IOR_PATCH==-1 .AND. IC/=IBAR2) CYCLE IC_LOOP
 
            N = CGI2(IC,JC,KC)
 
            SELECT CASE(IOR_PATCH)
 
               CASE(1)
                  II_LOW    = I_LO(IC)-1
                  II_HIGH   = I_LO(IC)-1
                  JJ_LOW    = J_LO(JC)
                  JJ_HIGH   = J_HI(JC)
                  KK_LOW    = K_LO(KC)
                  KK_HIGH   = K_HI(KC)
               CASE(-1)
                  II_LOW    = I_HI(IC)+1
                  II_HIGH   = I_HI(IC)+1
                  JJ_LOW    = J_LO(JC)
                  JJ_HIGH   = J_HI(JC)
                  KK_LOW    = K_LO(KC)
                  KK_HIGH   = K_HI(KC)
               CASE(2)
                  II_LOW    = I_LO(IC)
                  II_HIGH   = I_HI(IC)
                  JJ_LOW    = J_LO(JC)-1
                  JJ_HIGH   = J_LO(JC)-1
                  KK_LOW    = K_LO(KC)
                  KK_HIGH   = K_HI(KC)
               CASE(-2)
                  II_LOW    = I_LO(IC)
                  II_HIGH   = I_HI(IC)
                  JJ_LOW    = J_HI(JC)+1
                  JJ_HIGH   = J_HI(JC)+1
                  KK_LOW    = K_LO(KC)
                  KK_HIGH   = K_HI(KC)
               CASE(3)
                  II_LOW    = I_LO(IC)
                  II_HIGH   = I_HI(IC)
                  JJ_LOW    = J_LO(JC)
                  JJ_HIGH   = J_HI(JC)
                  KK_LOW    = K_LO(KC)-1
                  KK_HIGH   = K_LO(KC)-1
               CASE(-3)
                  II_LOW    = I_LO(IC)
                  II_HIGH   = I_HI(IC)
                  JJ_LOW    = J_LO(JC)
                  JJ_HIGH   = J_HI(JC)
                  KK_LOW    = K_HI(KC)+1
                  KK_HIGH   = K_HI(KC)+1
 
            END SELECT
 
            DHDX_S = 0._EB
            DHDX_F = 0._EB
            DHDY_S = 0._EB 
            DHDY_F = 0._EB 
            DHDZ_S = 0._EB 
            DHDZ_F = 0._EB 
 
            AREA_XS_CLOSED=0._EB
            AREA_XF_CLOSED=0._EB
            AREA_YS_CLOSED=0._EB
            AREA_YF_CLOSED=0._EB
            AREA_ZS_CLOSED=0._EB
            AREA_ZF_CLOSED=0._EB
 
            AREA_XS = 0._EB 
            AREA_XF = 0._EB 
            AREA_YS = 0._EB 
            AREA_YF = 0._EB 
            AREA_ZS = 0._EB 
            AREA_ZF = 0._EB 

            ! Loop over all external wall cells and compute (Sum dH/dn) for each face of mesh

            WALL_CELL_LOOP: DO IW=1,NEWC
 
               II  = IJKW(1,IW)
               JJ  = IJKW(2,IW)
               KK  = IJKW(3,IW)
               IOR = IJKW(4,IW)
               IF (IOR/=IOR_PATCH) CYCLE WALL_CELL_LOOP
               IF (II<II_LOW .OR. II>II_HIGH) CYCLE WALL_CELL_LOOP
               IF (JJ<JJ_LOW .OR. JJ>JJ_HIGH) CYCLE WALL_CELL_LOOP
               IF (KK<KK_LOW .OR. KK>KK_HIGH) CYCLE WALL_CELL_LOOP

               IF_INTERPOLATED_BOUNDARY: IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) THEN
                  NOM = IJKW(9,IW)
                  OM  => OMESH(NOM)
                  M2  => MESHES(NOM)
                  IIO = IJKW(10,IW)
                  JJO = IJKW(11,IW)
                  KKO = IJKW(12,IW)
                  NO = M2%CGI(IIO,JJO,KKO)
                  SELECT CASE(IOR)
                     CASE( 1) 
                        DA      = DY(JJ)*DZ(KK)
                        DHDX_S = DHDX_S + DA* (-DUDT_AVG(IW)-FVX(0,JJ,KK) + (B(N)-B(NO))/DX_M(N,NO) &
                              - (H(1,JJ,KK)-H(0,JJ,KK))*RDXN(0) )
                        AREA_XS = AREA_XS + DA
                     CASE(-1)
                        DA      = DY(JJ)*DZ(KK)
                        DHDX_F = DHDX_F + DA* (-DUDT_AVG(IW)-FVX(IBAR,JJ,KK) + (B(NO)-B(N))/DX_M(N,NO) &
                              - (H(IBP1,JJ,KK)-H(IBAR,JJ,KK))*RDXN(IBAR) )
                        AREA_XF = AREA_XF + DA
                     CASE( 2)
                        DA      = DX(II)*DZ(KK)
                        DHDY_S = DHDY_S + DA* (-DVDT_AVG(IW)-FVY(II,0,KK) + (B(N)-B(NO))/DY_M(N,NO) &
                              - (H(II,1,KK)-H(II,0,KK))*RDYN(0) )
                        AREA_YS = AREA_YS + DA
                     CASE(-2)
                        DA      = DX(II)*DZ(KK)
                        DHDY_F = DHDY_F + DA* (-DVDT_AVG(IW)-FVY(II,JBAR,KK) + (B(NO)-B(N))/DY_M(N,NO) &
                              - (H(II,JBP1,KK)-H(II,JBAR,KK))*RDYN(JBAR) )
                        AREA_YF = AREA_YF + DA
                     CASE( 3)
                        DA      = DX(II)*DY(JJ)
                        DHDZ_S = DHDZ_S + DA* (-DWDT_AVG(IW)-FVZ(II,JJ,0) + (B(N)-B(NO))/DZ_M(N,NO) &
                              - (H(II,JJ,1)-H(II,JJ,0))*RDZN(0) )
                        AREA_ZS = AREA_ZS + DA
                     CASE(-3)
                        DA      = DX(II)*DY(JJ)
                        DHDZ_F = DHDZ_F + DA* (-DWDT_AVG(IW)-FVZ(II,JJ,KBAR) + (B(NO)-B(N))/DZ_M(N,NO)    &
                              - (H(II,JJ,KBP1)-H(II,JJ,KBAR))*RDZN(KBAR) )
                        AREA_ZF = AREA_ZF + DA
                  END SELECT
               ENDIF IF_INTERPOLATED_BOUNDARY
 
               NON_INTERPOLATED: IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY .OR. BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
 
                  SELECT CASE(IOR)
                     CASE(-1)
                        DA      = DY(JJ)*DZ(KK)
                        AREA_XF = AREA_XF + DA
                        SELECT CASE (BOUNDARY_TYPE(IW))
                           CASE (SOLID_BOUNDARY)
                              DUUDT = -FVX(IBAR,JJ,KK) - RDXN(IBAR)*(H(IBP1,JJ,KK)-H(IBAR,JJ,KK))
                              AREA_XF_CLOSED = AREA_XF_CLOSED + DA
                              BXF(JJ,KK) = DUUDT - DUWDT(IW)
                           CASE (OPEN_BOUNDARY) 
                              B_OTHER = -B(N) 
                              DHDX_F = DHDX_F + DA* (  (B_OTHER-B(N))/DX_M(N,N) )
                        END SELECT
                     CASE( 1)
                        DA      = DY(JJ)*DZ(KK)
                        AREA_XS = AREA_XS + DA
                        SELECT CASE (BOUNDARY_TYPE(IW))
                           CASE (SOLID_BOUNDARY)
                              DUUDT = -FVX(0,JJ,KK) - RDXN(0)*(H(1,JJ,KK)-H(0,JJ,KK))
                              AREA_XS_CLOSED = AREA_XS_CLOSED + DA
                              BXS(JJ,KK) = DUUDT + DUWDT(IW)
                           CASE (OPEN_BOUNDARY) 
                              B_OTHER = -B(N)
                              DHDX_S = DHDX_S + DA* (  (B(N)-B_OTHER)/DX_M(N,N)  )
                        END SELECT
                     CASE(-2)
                        DA      = DX(II)*DZ(KK)
                        AREA_YF = AREA_YF + DA
                        SELECT CASE (BOUNDARY_TYPE(IW))
                           CASE (SOLID_BOUNDARY)
                              DVVDT = -FVY(II,JBAR,KK) - RDYN(JBAR)*(H(II,JBP1,KK)-H(II,JBAR,KK))
                              AREA_YF_CLOSED = AREA_YF_CLOSED + DA
                              BYF(II,KK) = DVVDT - DUWDT(IW)
                           CASE (OPEN_BOUNDARY)
                              B_OTHER = -B(N)
                              DHDY_F = DHDY_F + DA* (  (B_OTHER-B(N))/DY_M(N,N) )
                           END SELECT
                     CASE( 2)
                        DA      = DX(II)*DZ(KK)
                        AREA_YS = AREA_YS + DA
                        SELECT CASE (BOUNDARY_TYPE(IW))
                           CASE (SOLID_BOUNDARY)
                              DVVDT = -FVY(II,0,KK) - RDYN(0)*(H(II,1,KK)-H(II,0,KK))
                              AREA_YS_CLOSED = AREA_YS_CLOSED + DA
                              BYS(II,KK) = DVVDT + DUWDT(IW)
                           CASE (OPEN_BOUNDARY)
                              B_OTHER = -B(N)
                              DHDY_S = DHDY_S + DA* (  (B(N)-B_OTHER)/DY_M(N,N)  )
                           END SELECT
                     CASE(-3)
                        DA      = DX(II)*DY(JJ)
                        AREA_ZF = AREA_ZF + DA
                        SELECT CASE (BOUNDARY_TYPE(IW))
                           CASE (SOLID_BOUNDARY)
                              DWWDT = -FVZ(II,JJ,KBAR) - RDZN(KBAR)*(H(II,JJ,KBP1)-H(II,JJ,KBAR))
                              AREA_ZF_CLOSED = AREA_ZF_CLOSED + DA
                              BZF(II,JJ) = DWWDT - DUWDT(IW)
                           CASE (OPEN_BOUNDARY) 
                              B_OTHER = -B(N)
                              DHDZ_F = DHDZ_F + DA* (  (B_OTHER-B(N))/DZ_M(N,N)    )
                           END SELECT
                     CASE( 3)
                        DA      = DX(II)*DY(JJ)
                        AREA_ZS = AREA_ZS + DA
                        SELECT CASE (BOUNDARY_TYPE(IW))
                           CASE (SOLID_BOUNDARY)
                              DWWDT  = -FVZ(II,JJ,0) - RDZN(0)*(H(II,JJ,1)-H(II,JJ,0))
                              AREA_ZS_CLOSED = AREA_ZS_CLOSED + DA
                              BZS(II,JJ) = DWWDT + DUWDT(IW)
                           CASE (OPEN_BOUNDARY) 
                              B_OTHER = -B(N) 
                              DHDZ_S = DHDZ_S + DA* (  (B(N)-B_OTHER)/DZ_M(N,N) )
                           END SELECT
                  END SELECT
               ENDIF NON_INTERPOLATED
 
            ENDDO WALL_CELL_LOOP
 
            ! AREA_XX is only for the OPEN or INTERPOLATED boundary cells

            AREA_XS = AREA_XS - AREA_XS_CLOSED
            AREA_XF = AREA_XF - AREA_XF_CLOSED
            AREA_YS = AREA_YS - AREA_YS_CLOSED
            AREA_YF = AREA_YF - AREA_YF_CLOSED
            AREA_ZS = AREA_ZS - AREA_ZS_CLOSED
            AREA_ZF = AREA_ZF - AREA_ZF_CLOSED
             
            IF (AREA_XS==0._EB) AREA_XS=1._EB
            IF (AREA_XF==0._EB) AREA_XF=1._EB
            IF (AREA_YS==0._EB) AREA_YS=1._EB
            IF (AREA_YF==0._EB) AREA_YF=1._EB
            IF (AREA_ZS==0._EB) AREA_ZS=1._EB
            IF (AREA_ZF==0._EB) AREA_ZF=1._EB

            ! Loop over external wall cells and compute dH/dn and assign to BXS, etc.

            BC_LOOP: DO IW=1,NEWC
 
               II  = IJKW(1,IW)
               JJ  = IJKW(2,IW)
               KK  = IJKW(3,IW)
               IOR = IJKW(4,IW)
               IF (IOR/=IOR_PATCH) CYCLE BC_LOOP
               IF (II<II_LOW .OR. II>II_HIGH) CYCLE BC_LOOP
               IF (JJ<JJ_LOW .OR. JJ>JJ_HIGH) CYCLE BC_LOOP
               IF (KK<KK_LOW .OR. KK>KK_HIGH) CYCLE BC_LOOP
 
               BOUNDARY_SELECT: IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
 
                  SELECT CASE(IOR)
                     CASE( 1) 
                        BXS(JJ,KK) = HX(0)   *BXS(JJ,KK)
                     CASE(-1) 
                        BXF(JJ,KK) = HX(IBP1)*BXF(JJ,KK)
                     CASE( 2) 
                        BYS(II,KK) = HY(0)   *BYS(II,KK)
                     CASE(-2) 
                        BYF(II,KK) = HY(JBP1)*BYF(II,KK)
                     CASE( 3) 
                        BZS(II,JJ) = HZ(0)   *BZS(II,JJ)
                     CASE(-3) 
                        BZF(II,JJ) = HZ(KBP1)*BZF(II,JJ)
                  END SELECT
 
               ELSE BOUNDARY_SELECT
 
                  SELECT CASE(IOR)
                     CASE( 1) 
                        BXS(JJ,KK) = HX(0)   *(BXS(JJ,KK)+DHDX_S/AREA_XS)
                     CASE(-1) 
                        BXF(JJ,KK) = HX(IBP1)*(BXF(JJ,KK)+DHDX_F/AREA_XF)
                     CASE( 2) 
                        BYS(II,KK) = HY(0)   *(BYS(II,KK)+DHDY_S/AREA_YS)
                     CASE(-2) 
                        BYF(II,KK) = HY(JBP1)*(BYF(II,KK)+DHDY_F/AREA_YF)
                     CASE( 3) 
                        BZS(II,JJ) = HZ(0)   *(BZS(II,JJ)+DHDZ_S/AREA_ZS)
                     CASE(-3) 
                        BZF(II,JJ) = HZ(KBP1)*(BZF(II,JJ)+DHDZ_F/AREA_ZF)
                  END SELECT
 
               ENDIF BOUNDARY_SELECT
 
            ENDDO BC_LOOP
 
         ENDDO IC_LOOP
      ENDDO JC_LOOP
   ENDDO KC_LOOP
ENDDO ORIENT_LOOP

! Solve Laplace's Eq to get HP
 
HP   = 0._EB
PRHS = 0._EB
 
IF (.NOT.TWO_D) CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN, PRHS,POIS_PTB,SAVE2,WORK,HX)
 
IF (TWO_D .AND. .NOT. CYLINDRICAL) CALL H2CZSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE2,WORK,HX)
 
IF (ABS(POIS_PTB)>1.E-5_EB) WRITE(LU_ERR,*) ' POIS_PTB=',POIS_PTB, ' MESH=',NM
 
FORALL(I=1:IBAR,J=1:JBAR,K=1:KBAR) HP(I,J,K) = PRHS(I,J,K)
 
! Assign ghost cell values of HP

SET_BC_LOOP: DO IW=1,NEWC
   I   = IJKW(1,IW)
   J   = IJKW(2,IW)
   K   = IJKW(3,IW)
   IOR = IJKW(4,IW)
   SELECT CASE(IOR)
      CASE( 1) 
         HP(0,J,K)    = HP(1,J,K)    - DXI  *BXS(J,K)
      CASE(-1) 
         HP(IBP1,J,K) = HP(IBAR,J,K) + DXI  *BXF(J,K)
      CASE( 2) 
         HP(I,0,K)    = HP(I,1,K)    - DETA *BYS(I,K)
      CASE(-2) 
         HP(I,JBP1,K) = HP(I,JBAR,K) + DETA *BYF(I,K)
      CASE( 3) 
         HP(I,J,0)    = HP(I,J,1)    - DZETA*BZS(I,J)
      CASE(-3) 
         HP(I,J,KBP1) = HP(I,J,KBAR) + DZETA*BZF(I,J)
   END SELECT
ENDDO SET_BC_LOOP

END SUBROUTINE COMPUTE_CORRECTION_PRESSURE



SUBROUTINE GET_REV_pres(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') presrev(INDEX(presrev,':')+1:LEN_TRIM(presrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') presdate

END SUBROUTINE GET_REV_pres

END MODULE PRES
