MODULE PRES
 
! Find the perturbation pressure by solving Poisson's Equation
 
USE PRECISION_PARAMETERS
USE MESH_POINTERS

IMPLICIT NONE

PRIVATE
CHARACTER(255), PARAMETER :: presid='$Id$'
CHARACTER(255), PARAMETER :: presrev='$Revision$'
CHARACTER(255), PARAMETER :: presdate='$Date$'
PUBLIC PRESSURE_SOLVER,COMPUTE_VELOCITY_ERROR,GET_REV_PRES
 
CONTAINS
 
SUBROUTINE PRESSURE_SOLVER(T,NM)

USE POIS, ONLY: H3CZSS,H2CZSS,H2CYSS,H3CSSS
USE COMP_FUNCTIONS, ONLY: SECOND
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE GLOBAL_CONSTANTS
 
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,HP
REAL(EB), POINTER, DIMENSION(:) :: UWP
INTEGER :: I,J,K,IW,IOR,NOM,N_INT_CELLS,IIO,JJO,KKO
REAL(EB) :: TRM1,TRM2,TRM3,TRM4,RES,LHSS,RHSS,H_OTHER,DWWDT,DVVDT,DUUDT,HQ2,RFODT,U2,V2,W2,HFAC,H0RR(6),TNOW,DUMMY=0._EB, &
            TSI,TIME_RAMP_FACTOR,H_EXTERNAL,DX_OTHER,DY_OTHER,DZ_OTHER,P_EXTERNAL
TYPE (VENTS_TYPE), POINTER :: VT
 
IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY) RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   HP => H
   UWP=> UWS
ELSE
   UU => US
   VV => VS
   WW => WS
   HP => HS
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

!$OMP PARALLEL
!$OMP DO PRIVATE(IW,I,J,K,IOR,DUUDT,DVVDT,DWWDT,NOM,H_OTHER,KKO,JJO,IIO) &
!$OMP PRIVATE(N_INT_CELLS,DX_OTHER,DY_OTHER,DZ_OTHER,H_EXTERNAL,VT,TSI,TIME_RAMP_FACTOR) &
!$OMP PRIVATE(U2,V2,W2,HQ2)
WALL_CELL_LOOP: DO IW=1,NEWC
 
   I   = IJKW(1,IW)
   J   = IJKW(2,IW)
   K   = IJKW(3,IW)
   IOR = IJKW(4,IW)

   ! Apply pressure gradients at NEUMANN boundaries: dH/dn = -F_n - d(u_n)/dt

   IF_NEUMANN: IF (PRESSURE_BC_INDEX(IW)==NEUMANN) THEN
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
 
   IF_DIRICHLET: IF (PRESSURE_BC_INDEX(IW)==DIRICHLET) THEN

      NOT_OPEN: IF (BOUNDARY_TYPE(IW)/=OPEN_BOUNDARY .AND. BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) THEN

         ! Solid boundary that uses a Dirichlet BC to drive the normal component of velocity towards UWP
 
         SELECT CASE(IOR)
            CASE( 1)
               IF (PREDICTOR) DUUDT =       RFODT*(-UWP(IW)-        U(0,J,K)           )
               IF (CORRECTOR) DUUDT = 2._EB*RFODT*(-UWP(IW)-0.5_EB*(U(0,J,K)+US(0,J,K)))
               BXS(J,K) = HP(1,J,K)     + 0.5_EB*DX(0)   *(DUUDT+FVX(0,J,K))
            CASE(-1) 
               IF (PREDICTOR) DUUDT =       RFODT*( UWP(IW)-        U(IBAR,J,K)              )
               IF (CORRECTOR) DUUDT = 2._EB*RFODT*( UWP(IW)-0.5_EB*(U(IBAR,J,K)+US(IBAR,J,K)))
               BXF(J,K) = HP(IBAR,J,K) - 0.5_EB*DX(IBP1)*(DUUDT+FVX(IBAR,J,K))
            CASE( 2)
               IF (PREDICTOR) DVVDT =       RFODT*(-UWP(IW)-        V(I,0,K)           ) 
               IF (CORRECTOR) DVVDT = 2._EB*RFODT*(-UWP(IW)-0.5_EB*(V(I,0,K)+VS(I,0,K)))
               BYS(I,K) = HP(I,1,K)    + 0.5_EB*DY(0)   *(DVVDT+FVY(I,0,K))
            CASE(-2) 
               IF (PREDICTOR) DVVDT =       RFODT*( UWP(IW)-        V(I,JBAR,K)              )
               IF (CORRECTOR) DVVDT = 2._EB*RFODT*( UWP(IW)-0.5_EB*(V(I,JBAR,K)+VS(I,JBAR,K)))
               BYF(I,K) = HP(I,JBAR,K) - 0.5_EB*DY(JBP1)*(DVVDT+FVY(I,JBAR,K))
            CASE( 3)
               IF (PREDICTOR) DWWDT =       RFODT*(-UWP(IW)-        W(I,J,0)           )
               IF (CORRECTOR) DWWDT = 2._EB*RFODT*(-UWP(IW)-0.5_EB*(W(I,J,0)+WS(I,J,0)))
               BZS(I,J) = HP(I,J,1)    + 0.5_EB*DZ(0)   *(DWWDT+FVZ(I,J,0))
            CASE(-3) 
               IF (PREDICTOR) DWWDT =       RFODT*( UWP(IW)-        W(I,J,KBAR)              )
               IF (CORRECTOR) DWWDT = 2._EB*RFODT*( UWP(IW)-0.5_EB*(W(I,J,KBAR)+WS(I,J,KBAR)))
               BZF(I,J) = HP(I,J,KBAR) - 0.5_EB*DZ(KBP1)*(DWWDT+FVZ(I,J,KBAR))
         END SELECT

      ENDIF NOT_OPEN

      ! Interpolated boundary -- set boundary value of H to be average of neighboring cells from previous time step
 
      INTERPOLATED_ONLY: IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) THEN

         NOM     = IJKW(9,IW)
         H_OTHER = 0._EB
         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  IF (PREDICTOR) H_OTHER = H_OTHER + OMESH(NOM)%H(IIO,JJO,KKO)
                  IF (CORRECTOR) H_OTHER = H_OTHER + OMESH(NOM)%HS(IIO,JJO,KKO)
               ENDDO
            ENDDO
         ENDDO
         N_INT_CELLS = (IJKW(13,IW)-IJKW(10,IW)+1) * (IJKW(14,IW)-IJKW(11,IW)+1) * (IJKW(15,IW)-IJKW(12,IW)+1)
         H_OTHER = H_OTHER/REAL(N_INT_CELLS,EB)
         
         SELECT CASE(IOR)
            CASE( 1)
               DX_OTHER = MESHES(NOM)%DX(IJKW(10,IW))
               BXS(J,K) = (DX_OTHER*HP(1,J,K) + DX(1)*H_OTHER)/(DX(1)+DX_OTHER) + WALL_WORK1(IW)
            CASE(-1)
               DX_OTHER = MESHES(NOM)%DX(IJKW(10,IW))
               BXF(J,K) = (DX_OTHER*HP(IBAR,J,K) + DX(IBAR)*H_OTHER)/(DX(IBAR)+DX_OTHER)+ WALL_WORK1(IW)
            CASE( 2)
               DY_OTHER = MESHES(NOM)%DY(IJKW(11,IW))
               BYS(I,K) = (DY_OTHER*HP(I,1,K) + DY(1)*H_OTHER)/(DY(1)+DY_OTHER)+ WALL_WORK1(IW)
            CASE(-2)
               DY_OTHER = MESHES(NOM)%DY(IJKW(11,IW))
               BYF(I,K) = (DY_OTHER*HP(I,JBAR,K) + DY(JBAR)*H_OTHER)/(DY(JBAR)+DY_OTHER)+ WALL_WORK1(IW)
            CASE( 3)
               DZ_OTHER = MESHES(NOM)%DZ(IJKW(12,IW))
               BZS(I,J) = (DZ_OTHER*HP(I,J,1) + DZ(1)*H_OTHER)/(DZ(1)+DZ_OTHER)+ WALL_WORK1(IW)
            CASE(-3)
               DZ_OTHER = MESHES(NOM)%DZ(IJKW(12,IW))
               BZF(I,J) = (DZ_OTHER*HP(I,J,KBAR) + DZ(KBAR)*H_OTHER)/(DZ(KBAR)+DZ_OTHER)+ WALL_WORK1(IW)
         END SELECT
 
      ENDIF INTERPOLATED_ONLY
 
      ! OPEN (passive opening to exterior of domain) boundary. Apply inflow/outflow BC.

      OPEN: IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) THEN
 
         H_EXTERNAL = 0._EB
         IF (VENT_INDEX(IW)>0) THEN
            VT => VENTS(VENT_INDEX(IW))
            IF (TW(IW) == T_BEGIN .AND. VT%PRESSURE_RAMP_INDEX >=1) THEN
               TSI = T
            ELSE
               TSI = T - T_BEGIN
            ENDIF
            TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,DUMMY,VT%PRESSURE_RAMP_INDEX)
            !!H_EXTERNAL = TIME_RAMP_FACTOR*VT%DYNAMIC_PRESSURE/RHOA
            IF (.NOT.BAROCLINIC2) H_EXTERNAL = TIME_RAMP_FACTOR*VT%DYNAMIC_PRESSURE/RHO_AVG !! test
            IF (BAROCLINIC2)      P_EXTERNAL = TIME_RAMP_FACTOR*VT%DYNAMIC_PRESSURE         !! also test
         ENDIF

         SELECT CASE(IOR)
            CASE( 1)
               U2  = UU(0,J,K)**2
               V2  = .25_EB*(VV(1,J,K)+VV(1,J-1,K))**2
               W2  = .25_EB*(WW(1,J,K)+WW(1,J,K-1))**2
               HQ2 = MIN(5000._EB,0.5_EB*(U2+V2+W2))
               IF (UU(0,J,K)<0._EB) THEN
                  IF (.NOT.BAROCLINIC2) BXS(J,K) = H_EXTERNAL + HQ2
                  IF (BAROCLINIC2)      BXS(J,K) = P_EXTERNAL/RHO_W(IW) + KRES(0,J,K)
               ELSE
                  !BXS(J,K) = H_EXTERNAL + H0RR(1) + HQ2*HFAC
                  IF (.NOT.BAROCLINIC2) BXS(J,K) = H_EXTERNAL + H0RR(1) + HQ2*(1._EB-RHO_W(IW)/RHO_AVG)
                  IF (BAROCLINIC2)      BXS(J,K) = P_EXTERNAL/RHO_W(IW) + H0
               ENDIF
            CASE(-1)
               U2  = UU(IBAR,J,K)**2
               V2  = .25_EB*(VV(IBAR,J,K)+VV(IBAR,J-1,K))**2
               W2  = .25_EB*(WW(IBAR,J,K)+WW(IBAR,J,K-1))**2
               HQ2 = MIN(5000._EB,0.5_EB*(U2+V2+W2))
               IF (UU(IBAR,J,K)>0._EB) THEN
                  IF (.NOT.BAROCLINIC2) BXF(J,K) = H_EXTERNAL + HQ2
                  IF (BAROCLINIC2)      BXF(J,K) = P_EXTERNAL/RHO_W(IW) + KRES(IBP1,J,K)
               ELSE
                  !BXF(J,K) = H_EXTERNAL + H0RR(2) + HQ2*HFAC
                  IF (.NOT.BAROCLINIC2) BXF(J,K) = H_EXTERNAL + H0RR(2) + HQ2*(1._EB-RHO_W(IW)/RHO_AVG)
                  IF (BAROCLINIC2)      BXF(J,K) = P_EXTERNAL/RHO_W(IW) + H0
               ENDIF
            CASE( 2)
               U2  = .25_EB*(UU(I,1,K)+UU(I-1,1,K))**2
               V2  = VV(I,0,K)**2
               W2  = .25_EB*(WW(I,1,K)+WW(I,1,K-1))**2
               HQ2 = MIN(5000._EB,0.5_EB*(U2+V2+W2))
               IF (VV(I,0,K)<0._EB) THEN
                  IF (.NOT.BAROCLINIC2) BYS(I,K) = H_EXTERNAL + HQ2
                  IF (BAROCLINIC2)      BYS(I,K) = P_EXTERNAL/RHO_W(IW) + KRES(I,0,K)
               ELSE
                  !BYS(I,K) = H_EXTERNAL + H0RR(3) + HQ2*HFAC
                  IF (.NOT.BAROCLINIC2) BYS(I,K) = H_EXTERNAL + H0RR(3) + HQ2*(1._EB-RHO_W(IW)/RHO_AVG)
                  IF (BAROCLINIC2)      BYS(I,K) = P_EXTERNAL/RHO_W(IW) + H0
               ENDIF
            CASE(-2)
               U2  = .25_EB*(UU(I,JBAR,K)+UU(I-1,JBAR,K))**2
               V2  = VV(I,JBAR,K)**2
               W2  = .25_EB*(WW(I,JBAR,K)+WW(I,JBAR,K-1))**2
               HQ2 = MIN(5000._EB,0.5_EB*(U2+V2+W2))
               IF (VV(I,JBAR,K)>0._EB) THEN
                  IF (.NOT.BAROCLINIC2) BYF(I,K) = H_EXTERNAL + HQ2
                  IF (BAROCLINIC2)      BYF(I,K) = P_EXTERNAL/RHO_W(IW) + KRES(I,JBP1,K)
               ELSE
                  !BYF(I,K) = H_EXTERNAL + H0RR(4) + HQ2*HFAC
                  IF (.NOT.BAROCLINIC2) BYF(I,K) = H_EXTERNAL + H0RR(4) + HQ2*(1._EB-RHO_W(IW)/RHO_AVG)
                  IF (BAROCLINIC2)      BYF(I,K) = P_EXTERNAL/RHO_W(IW) + H0
               ENDIF
            CASE( 3)
               U2  = .25_EB*(UU(I,J,1)+UU(I-1,J,1))**2
               V2  = .25_EB*(VV(I,J,1)+VV(I,J-1,1))**2
               W2  = WW(I,J,0)**2
               HQ2 = MIN(5000._EB,0.5_EB*(U2+V2+W2))
               IF (WW(I,J,0)<0._EB) THEN
                  IF (.NOT.BAROCLINIC2) BZS(I,J) = H_EXTERNAL + HQ2
                  IF (BAROCLINIC2)      BZS(I,J) = P_EXTERNAL/RHO_W(IW) + KRES(I,J,0)
               ELSE
                  !BZS(I,J) = H_EXTERNAL + H0RR(5) + HQ2*HFAC
                  IF (.NOT.BAROCLINIC2) BZS(I,J) = H_EXTERNAL + H0RR(5) + HQ2*(1._EB-RHO_W(IW)/RHO_AVG)
                  IF (BAROCLINIC2)      BZS(I,J) = P_EXTERNAL/RHO_W(IW) + H0
               ENDIF
            CASE(-3)
               U2  = .25_EB*(UU(I,J,KBAR)+UU(I-1,J,KBAR))**2
               V2  = .25_EB*(VV(I,J,KBAR)+VV(I,J-1,KBAR))**2
               W2  = WW(I,J,KBAR)**2
               HQ2 = MIN(5000._EB,0.5_EB*(U2+V2+W2))
               IF (WW(I,J,KBAR)>0._EB) THEN
                  IF (.NOT.BAROCLINIC2) BZF(I,J) = H_EXTERNAL + HQ2
                  IF (BAROCLINIC2)      BZF(I,J) = P_EXTERNAL/RHO_W(IW) + KRES(I,J,KBP1)
               ELSE
                  !BZF(I,J) = H_EXTERNAL + H0RR(6) + HQ2*HFAC
                  IF (.NOT.BAROCLINIC2) BZF(I,J) = H_EXTERNAL + H0RR(6) + HQ2*(1._EB-RHO_W(IW)/RHO_AVG)
                  IF (BAROCLINIC2)      BZF(I,J) = P_EXTERNAL/RHO_W(IW) + H0
               ENDIF
         END SELECT
    
      ENDIF OPEN
 
   ENDIF IF_DIRICHLET
 
ENDDO WALL_CELL_LOOP
!$OMP END DO

! Compute the RHS of the Poisson equation
 
SELECT CASE(IPS)

   CASE(:1,4,7)
      IF (CYLINDRICAL) THEN
         !$OMP DO COLLAPSE(2) PRIVATE(K,I,TRM1,TRM3,TRM4)
         DO K=1,KBAR
            DO I=1,IBAR
               TRM1 = (R(I-1)*FVX(I-1,1,K)-R(I)*FVX(I,1,K))*RDX(I)*RRN(I)
               TRM3 = (FVZ(I,1,K-1)-FVZ(I,1,K))*RDZ(K)
               TRM4 = -DDDT(I,1,K)
               PRHS(I,1,K) = TRM1 + TRM3 + TRM4
            ENDDO
         ENDDO
         !$OMP END DO NOWAIT
      ENDIF
      IF (.NOT.CYLINDRICAL) THEN
         !$OMP DO COLLAPSE(3) PRIVATE(TRM1,TRM2,TRM3,TRM4)
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
         !$OMP END DO NOWAIT
      ENDIF
 
   CASE(2)  ! Switch x and y
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,TRM1, TRM2,TRM3,TRM4)
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
      !$OMP END DO NOWAIT
      !$OMP WORKSHARE
      BZST = TRANSPOSE(BZS)
      BZFT = TRANSPOSE(BZF)
      !$OMP END WORKSHARE
 
   CASE(3,6)  ! Switch x and z
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,TRM1,TRM2,TRM3,TRM4)
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
      !$OMP END DO NOWAIT
      !$OMP WORKSHARE
      BXST = TRANSPOSE(BXS)
      BXFT = TRANSPOSE(BXF)
      BYST = TRANSPOSE(BYS)
      BYFT = TRANSPOSE(BYF)
      BZST = TRANSPOSE(BZS)
      BZFT = TRANSPOSE(BZF)
      !$OMP END WORKSHARE
 
   CASE(5)  ! Switch y and z
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,TRM1,TRM2,TRM3,TRM4)
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
      !$OMP END DO NOWAIT
      !$OMP WORKSHARE
      BXST = TRANSPOSE(BXS)
      BXFT = TRANSPOSE(BXF)
      !$OMP END WORKSHARE

END SELECT
!$OMP END PARALLEL

! Call the Poisson solver
 
SELECT CASE(IPS)
   CASE(:1) 
      IF (.NOT.TWO_D) CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
      IF (TWO_D .AND. .NOT. CYLINDRICAL) CALL H2CZSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
      IF (TWO_D .AND.       CYLINDRICAL) CALL H2CYSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE1,WORK)
   CASE(2) 
      CALL H3CZSS(BYS,BYF,BXS,BXF,BZST,BZFT,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HY)
   CASE(3)      
      IF (.NOT.TWO_D) CALL H3CZSS(BZST,BZFT,BYST,BYFT,BXST,BXFT,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HZ)
      IF (TWO_D)      CALL H2CZSS(BZS,BZF,BXS,BXF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HZ)
   CASE(4)  
      CALL H3CSSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HX,HY)
   CASE(5)      
      IF (.NOT.TWO_D) CALL H3CSSS(BXST,BXFT,BZS,BZF,BYS,BYF,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HX,HZ)
      IF (     TWO_D) CALL H2CZSS(BZS,BZF,BXS,BXF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HZ)
   CASE(6) 
      CALL H3CSSS(BZST,BZFT,BYST,BYFT,BXST,BXFT,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HZ,HY)
   CASE(7)  
      CALL H2CZSS(BXS,BXF,BYS,BYF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
END SELECT
 
! Put output of Poisson solver into the H array

!$OMP PARALLEL
SELECT CASE(IPS)
   CASE(:1,4,7)
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               HP(I,J,K) = PRHS(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(2)
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               HP(I,J,K) = PRHS(J,I,K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(3,6)
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               HP(I,J,K) = PRHS(K,J,I)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(5)
      !$OMP DO COLLAPSE(3) PRIVATE(K,J,I)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               HP(I,J,K) = PRHS(I,K,J)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
END SELECT 

! Apply boundary conditions to H

!$OMP DO COLLAPSE(2) PRIVATE(K,J)
DO K=1,KBAR
   DO J=1,JBAR
      IF (LBC==3 .OR. LBC==4)             HP(0,J,K)    = HP(1,J,K)    - DXI*BXS(J,K)
      IF (LBC==3 .OR. LBC==2 .OR. LBC==6) HP(IBP1,J,K) = HP(IBAR,J,K) + DXI*BXF(J,K)
      IF (LBC==1 .OR. LBC==2)             HP(0,J,K)    =-HP(1,J,K)    + 2._EB*BXS(J,K)
      IF (LBC==1 .OR. LBC==4 .OR. LBC==5) HP(IBP1,J,K) =-HP(IBAR,J,K) + 2._EB*BXF(J,K)
      IF (LBC==5 .OR. LBC==6)             HP(0,J,K) = HP(1,J,K)
      IF (LBC==0) THEN
         HP(0,J,K) = HP(IBAR,J,K)
         HP(IBP1,J,K) = HP(1,J,K)
      ENDIF
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2) PRIVATE(K,I)
DO K=1,KBAR
   DO I=1,IBAR
      IF (MBC==3 .OR. MBC==4) HP(I,0,K)    = HP(I,1,K)    - DETA*BYS(I,K)
      IF (MBC==3 .OR. MBC==2) HP(I,JBP1,K) = HP(I,JBAR,K) + DETA*BYF(I,K)
      IF (MBC==1 .OR. MBC==2) HP(I,0,K)    =-HP(I,1,K)    + 2._EB*BYS(I,K)
      IF (MBC==1 .OR. MBC==4) HP(I,JBP1,K) =-HP(I,JBAR,K) + 2._EB*BYF(I,K)
      IF (MBC==0) THEN
         HP(I,0,K) = HP(I,JBAR,K)
         HP(I,JBP1,K) = HP(I,1,K)
      ENDIF
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(2) PRIVATE(J,I)
DO J=1,JBAR
   DO I=1,IBAR
      IF (EVACUATION_ONLY(NM)) CYCLE
      IF (NBC==3 .OR. NBC==4)  HP(I,J,0)    = HP(I,J,1)    - DZETA*BZS(I,J)
      IF (NBC==3 .OR. NBC==2)  HP(I,J,KBP1) = HP(I,J,KBAR) + DZETA*BZF(I,J)
      IF (NBC==1 .OR. NBC==2)  HP(I,J,0)    =-HP(I,J,1)    + 2._EB*BZS(I,J)
      IF (NBC==1 .OR. NBC==4)  HP(I,J,KBP1) =-HP(I,J,KBAR) + 2._EB*BZF(I,J)
      IF (NBC==0) THEN
         HP(I,J,0) = HP(I,J,KBAR)
         HP(I,J,KBP1) = HP(I,J,1)
      ENDIF
   ENDDO
ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

! Optional check of the accuracy of the pressure solver

IF (CHECK_POISSON .AND. .NOT.EVACUATION_ONLY(NM)) THEN
   POIS_ERR = 0.
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHSS = ( R(I-1)*FVX(I-1,J,K) - R(I)*FVX(I,J,K) )*RDX(I)*RRN(I) &
                 + (        FVY(I,J-1,K) -      FVY(I,J,K) )*RDY(J)        &
                 + (        FVZ(I,J,K-1) -      FVZ(I,J,K) )*RDZ(K)        &
                 - DDDT(I,J,K)
            LHSS = ((HP(I+1,J,K)-HP(I,J,K))*RDXN(I)*R(I) - (HP(I,J,K)-HP(I-1,J,K))*RDXN(I-1)*R(I-1) )*RDX(I)*RRN(I) &
                 + ((HP(I,J+1,K)-HP(I,J,K))*RDYN(J)      - (HP(I,J,K)-HP(I,J-1,K))*RDYN(J-1)        )*RDY(J)        &
                 + ((HP(I,J,K+1)-HP(I,J,K))*RDZN(K)      - (HP(I,J,K)-HP(I,J,K-1))*RDZN(K-1)        )*RDZ(K)
            RES = ABS(RHSS-LHSS)
            POIS_ERR = MAX(RES,POIS_ERR)
         ENDDO
      ENDDO
   ENDDO
ENDIF

TUSED(5,NM)=TUSED(5,NM)+SECOND()-TNOW
END SUBROUTINE PRESSURE_SOLVER



SUBROUTINE COMPUTE_VELOCITY_ERROR(NM)

! Check the maximum velocity error at a solid boundary

USE COMP_FUNCTIONS, ONLY: SECOND
USE GLOBAL_CONSTANTS

INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:) :: UWP
INTEGER :: IW,IOR,II,JJ,KK,IIO,JJO,KKO
REAL(EB) :: TNOW,U_NEW,V_NEW,W_NEW,U_NEW_OTHER,V_NEW_OTHER,W_NEW_OTHER,VELOCITY_ERROR
TYPE(OMESH_TYPE), POINTER :: OM

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   UWP=> UWS
ELSE
   UWP=> UW
ENDIF

VELOCITY_ERROR_MAX(NM) = 0._EB
WALL_WORK1 = 0._EB

CHECK_WALL_LOOP: DO IW=1,NWC

   IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY .AND. BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) CYCLE CHECK_WALL_LOOP

   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IOR = IJKW(4,IW)

   ! Update normal component of velocity at the mesh boundary

   IF (PREDICTOR) THEN
      SELECT CASE(IOR)
         CASE( 1)
            U_NEW = U(II,JJ,KK)   - DT*(FVX(II,JJ,KK)   + RDXN(II)  *(H(II+1,JJ,KK)-H(II,JJ,KK)))
         CASE(-1)
            U_NEW = U(II-1,JJ,KK) - DT*(FVX(II-1,JJ,KK) + RDXN(II-1)*(H(II,JJ,KK)-H(II-1,JJ,KK)))
         CASE( 2)
            V_NEW = V(II,JJ,KK)   - DT*(FVY(II,JJ,KK)   + RDYN(JJ)  *(H(II,JJ+1,KK)-H(II,JJ,KK)))
         CASE(-2)
            V_NEW = V(II,JJ-1,KK) - DT*(FVY(II,JJ-1,KK) + RDYN(JJ-1)*(H(II,JJ,KK)-H(II,JJ-1,KK)))
         CASE( 3)
            W_NEW = W(II,JJ,KK)   - DT*(FVZ(II,JJ,KK)   + RDZN(KK)  *(H(II,JJ,KK+1)-H(II,JJ,KK)))
         CASE(-3)
            W_NEW = W(II,JJ,KK-1) - DT*(FVZ(II,JJ,KK-1) + RDZN(KK-1)*(H(II,JJ,KK)-H(II,JJ,KK-1)))
      END SELECT
   ELSE
      SELECT CASE(IOR)
         CASE( 1)
            U_NEW = 0.5_EB*(U(II,JJ,KK)+US(II,JJ,KK)     - DT*(FVX(II,JJ,KK)   + RDXN(II)  *(HS(II+1,JJ,KK)-HS(II,JJ,KK))))
         CASE(-1)
            U_NEW = 0.5_EB*(U(II-1,JJ,KK)+US(II-1,JJ,KK) - DT*(FVX(II-1,JJ,KK) + RDXN(II-1)*(HS(II,JJ,KK)-HS(II-1,JJ,KK))))
         CASE( 2)
            V_NEW = 0.5_EB*(V(II,JJ,KK)+VS(II,JJ,KK)     - DT*(FVY(II,JJ,KK)   + RDYN(JJ)  *(HS(II,JJ+1,KK)-HS(II,JJ,KK))))
         CASE(-2)
            V_NEW = 0.5_EB*(V(II,JJ-1,KK)+VS(II,JJ-1,KK) - DT*(FVY(II,JJ-1,KK) + RDYN(JJ-1)*(HS(II,JJ,KK)-HS(II,JJ-1,KK))))
         CASE( 3)
            W_NEW = 0.5_EB*(W(II,JJ,KK)+WS(II,JJ,KK)     - DT*(FVZ(II,JJ,KK)   + RDZN(KK)  *(HS(II,JJ,KK+1)-HS(II,JJ,KK))))
         CASE(-3)
            W_NEW = 0.5_EB*(W(II,JJ,KK-1)+WS(II,JJ,KK-1) - DT*(FVZ(II,JJ,KK-1) + RDZN(KK-1)*(HS(II,JJ,KK)-HS(II,JJ,KK-1))))
      END SELECT
   ENDIF

   ! At interpolated boundaries, compare updated normal component of velocity with that of the other mesh

   IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) THEN
      OM  => OMESH(IJKW(9,IW))
      IIO = IJKW(10,IW)
      JJO = IJKW(11,IW)
      KKO = IJKW(12,IW)
      IF (PREDICTOR) THEN
         SELECT CASE(IOR)
            CASE( 1)
               U_NEW_OTHER = OM%U(IIO,JJO,KKO)   + DT*OM%DUDT(IIO,JJO,KKO)
            CASE(-1)
               U_NEW_OTHER = OM%U(IIO-1,JJO,KKO) + DT*OM%DUDT(IIO-1,JJO,KKO)
            CASE( 2)
               V_NEW_OTHER = OM%V(IIO,JJO,KKO)   + DT*OM%DVDT(IIO,JJO,KKO)
            CASE(-2)
               V_NEW_OTHER = OM%V(IIO,JJO-1,KKO) + DT*OM%DVDT(IIO,JJO-1,KKO)
            CASE( 3)
               W_NEW_OTHER = OM%W(IIO,JJO,KKO)   + DT*OM%DWDT(IIO,JJO,KKO)
            CASE(-3)
               W_NEW_OTHER = OM%W(IIO,JJO,KKO-1) + DT*OM%DWDT(IIO,JJO,KKO-1)
         END SELECT
      ELSE
         SELECT CASE(IOR)
            CASE( 1)
               U_NEW_OTHER = 0.5_EB*(OM%U(IIO,JJO,KKO)+OM%US(IIO,JJO,KKO)     + DT*OM%DUDT(IIO,JJO,KKO))
            CASE(-1)
               U_NEW_OTHER = 0.5_EB*(OM%U(IIO-1,JJO,KKO)+OM%US(IIO-1,JJO,KKO) + DT*OM%DUDT(IIO-1,JJO,KKO))
            CASE( 2)
               V_NEW_OTHER = 0.5_EB*(OM%V(IIO,JJO,KKO)+OM%VS(IIO,JJO,KKO)     + DT*OM%DVDT(IIO,JJO,KKO))
            CASE(-2)
               V_NEW_OTHER = 0.5_EB*(OM%V(IIO,JJO-1,KKO)+OM%VS(IIO,JJO-1,KKO) + DT*OM%DVDT(IIO,JJO-1,KKO))
            CASE( 3)
               W_NEW_OTHER = 0.5_EB*(OM%W(IIO,JJO,KKO)+OM%WS(IIO,JJO,KKO)     + DT*OM%DWDT(IIO,JJO,KKO))
            CASE(-3)
               W_NEW_OTHER = 0.5_EB*(OM%W(IIO,JJO,KKO-1)+OM%WS(IIO,JJO,KKO-1) + DT*OM%DWDT(IIO,JJO,KKO-1))
         END SELECT
      ENDIF
   ENDIF

   ! At solid boundaries, compare updated normal velocity with specified normal velocity (UWP)

   IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
      SELECT CASE(IOR)
         CASE( 1)
            U_NEW_OTHER = -UWP(IW)
         CASE(-1)
            U_NEW_OTHER =  UWP(IW)
         CASE( 2)
            V_NEW_OTHER = -UWP(IW)
         CASE(-2)
            V_NEW_OTHER =  UWP(IW)
         CASE( 3)
            W_NEW_OTHER = -UWP(IW)
         CASE(-3)
            W_NEW_OTHER =  UWP(IW)
      END SELECT
   ENDIF

   ! Compute velocity difference

   SELECT CASE(ABS(IOR))
      CASE(1)
         VELOCITY_ERROR = U_NEW - U_NEW_OTHER
         IF (IOR>0) WALL_WORK1(IW) = -1.0*VELOCITY_ERROR*DX(II)/(4.*DT)
         IF (IOR<0) WALL_WORK1(IW) =  1.0*VELOCITY_ERROR*DX(II)/(4.*DT)
      CASE(2)
         VELOCITY_ERROR = V_NEW - V_NEW_OTHER
         IF (IOR>0) WALL_WORK1(IW) = -1.0*VELOCITY_ERROR*DY(JJ)/(4.*DT)
         IF (IOR<0) WALL_WORK1(IW) =  1.0*VELOCITY_ERROR*DY(JJ)/(4.*DT)
      CASE(3)
         VELOCITY_ERROR = W_NEW - W_NEW_OTHER
         IF (IOR>0) WALL_WORK1(IW) = -1.0*VELOCITY_ERROR*DZ(KK)/(4.*DT)
         IF (IOR<0) WALL_WORK1(IW) =  1.0*VELOCITY_ERROR*DZ(KK)/(4.*DT)
   END SELECT

   IF (ABS(VELOCITY_ERROR)>VELOCITY_ERROR_MAX(NM)) then
      VELOCITY_ERROR_MAX_I(NM) = II
      VELOCITY_ERROR_MAX_J(NM) = JJ
      VELOCITY_ERROR_MAX_K(NM) = KK
      VELOCITY_ERROR_MAX(NM)   = ABS(VELOCITY_ERROR)
   endif

ENDDO CHECK_WALL_LOOP

TUSED(5,NM)=TUSED(5,NM)+SECOND()-TNOW
END SUBROUTINE COMPUTE_VELOCITY_ERROR



SUBROUTINE GET_REV_pres(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') presrev(INDEX(presrev,':')+1:LEN_TRIM(presrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') presdate

END SUBROUTINE GET_REV_pres


END MODULE PRES
