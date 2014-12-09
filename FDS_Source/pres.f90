MODULE PRES
 
! Find the perturbation pressure by solving Poisson's Equation
 
USE PRECISION_PARAMETERS
USE MESH_POINTERS

IMPLICIT NONE

PRIVATE
CHARACTER(255), PARAMETER :: presid='$Id$'
CHARACTER(255), PARAMETER :: presrev='$Revision$'
CHARACTER(255), PARAMETER :: presdate='$Date$'

PUBLIC PRESSURE_SOLVER,COMPUTE_VELOCITY_ERROR,GET_REV_PRES!,BUILD_SPARSE_MATRIX_LAPLACE
 
CONTAINS
 
SUBROUTINE PRESSURE_SOLVER(T,NM)

USE POIS, ONLY: H3CZSS,H2CZSS,H2CYSS,H3CSSS
USE COMP_FUNCTIONS, ONLY: SECOND
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP, AFILL2
USE GLOBAL_CONSTANTS
USE TRAN, ONLY: GET_IJK
USE TURBULENCE, ONLY: NS_H_EXACT
 
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,HP
INTEGER :: I,J,K,IW,IOR,NOM,N_INT_CELLS,IIO,JJO,KKO,IIX,JJY,KKZ
REAL(EB) :: TRM1,TRM2,TRM3,TRM4,RES,LHSS,RHSS,H_OTHER,TNOW,DUMMY=0._EB, &
            TSI,TIME_RAMP_FACTOR,DX_OTHER,DY_OTHER,DZ_OTHER,P_EXTERNAL, &
            XIO,YJO,ZKO,PP,RR,SS,UBAR,VBAR,WBAR
TYPE (VENTS_TYPE), POINTER :: VT
TYPE (WALL_TYPE), POINTER :: WC
 
IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   HP => H
ELSE
   UU => US
   VV => VS
   WW => WS
   HP => HS
ENDIF

! Apply pressure boundary conditions at external cells.
! If Neumann, BXS, BXF, etc., contain dH/dx(x=XS), dH/dx(x=XF), etc.
! If Dirichlet, BXS, BXF, etc., contain H(x=XS), H(x=XF), etc.
! LBC, MBC and NBC are codes used be Poisson solver to denote type
! of boundary condition at x, y and z boundaries. See Crayfishpak
! manual for details.

WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

   WC => WALL(IW)
   I   = WC%ONE_D%II
   J   = WC%ONE_D%JJ
   K   = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR

   ! Apply pressure gradients at NEUMANN boundaries: dH/dn = -F_n - d(u_n)/dt

   IF_NEUMANN: IF (WC%PRESSURE_BC_INDEX==NEUMANN) THEN
      SELECT CASE(IOR)
         CASE( 1)
            BXS(J,K) = HX(0)   *(-FVX(0,J,K)    + WC%DUWDT)
         CASE(-1)
            BXF(J,K) = HX(IBP1)*(-FVX(IBAR,J,K) - WC%DUWDT)
         CASE( 2)
            BYS(I,K) = HY(0)   *(-FVY(I,0,K)    + WC%DUWDT)
         CASE(-2)
            BYF(I,K) = HY(JBP1)*(-FVY(I,JBAR,K) - WC%DUWDT)
         CASE( 3)
            BZS(I,J) = HZ(0)   *(-FVZ(I,J,0)    + WC%DUWDT)
         CASE(-3)
            BZF(I,J) = HZ(KBP1)*(-FVZ(I,J,KBAR) - WC%DUWDT)
      END SELECT
   ENDIF IF_NEUMANN

   ! Apply pressures at DIRICHLET boundaries, depending on the specific type
 
   IF_DIRICHLET: IF (WC%PRESSURE_BC_INDEX==DIRICHLET) THEN

      NOT_OPEN: IF (WC%BOUNDARY_TYPE/=OPEN_BOUNDARY .AND. WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) THEN

         ! Solid boundary that uses a Dirichlet BC. Assume that the pressure at the boundary (BXS, etc) is the average of the
         ! last computed pressures in the ghost and adjacent gas cells.
 
         SELECT CASE(IOR)
            CASE( 1)
               BXS(J,K) = 0.5_EB*(HP(0,J,K)+HP(1,J,K))
            CASE(-1) 
               BXF(J,K) = 0.5_EB*(HP(IBAR,J,K)+HP(IBP1,J,K))
            CASE( 2)
               BYS(I,K) = 0.5_EB*(HP(I,0,K)+HP(I,1,K))
            CASE(-2) 
               BYF(I,K) = 0.5_EB*(HP(I,JBAR,K)+HP(I,JBP1,K))
            CASE( 3)
               BZS(I,J) = 0.5_EB*(HP(I,J,0)+HP(I,J,1))
            CASE(-3) 
               BZF(I,J) = 0.5_EB*(HP(I,J,KBAR)+HP(I,J,KBP1))
         END SELECT

      ENDIF NOT_OPEN

      ! Interpolated boundary -- set boundary value of H to be average of neighboring cells from previous time step
 
      INTERPOLATED_ONLY: IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
         
         NOM     = WC%NOM
         H_OTHER = 0._EB
         DO KKO=WC%NOM_IB(3),WC%NOM_IB(6)
            DO JJO=WC%NOM_IB(2),WC%NOM_IB(5)
               DO IIO=WC%NOM_IB(1),WC%NOM_IB(4)
                  IF (PREDICTOR) H_OTHER = H_OTHER + OMESH(NOM)%H(IIO,JJO,KKO)
                  IF (CORRECTOR) H_OTHER = H_OTHER + OMESH(NOM)%HS(IIO,JJO,KKO)
               ENDDO
            ENDDO
         ENDDO
         N_INT_CELLS = (WC%NOM_IB(4)-WC%NOM_IB(1)+1) * (WC%NOM_IB(5)-WC%NOM_IB(2)+1) * (WC%NOM_IB(6)-WC%NOM_IB(3)+1)
         H_OTHER = H_OTHER/REAL(N_INT_CELLS,EB)
         
         SELECT CASE(IOR)
            CASE( 1)
               DX_OTHER = MESHES(NOM)%DX(WC%NOM_IB(1))
               BXS(J,K) = (DX_OTHER*HP(1,J,K) + DX(1)*H_OTHER)/(DX(1)+DX_OTHER) + WALL_WORK1(IW)
            CASE(-1)
               DX_OTHER = MESHES(NOM)%DX(WC%NOM_IB(1))
               BXF(J,K) = (DX_OTHER*HP(IBAR,J,K) + DX(IBAR)*H_OTHER)/(DX(IBAR)+DX_OTHER) + WALL_WORK1(IW)
            CASE( 2)
               DY_OTHER = MESHES(NOM)%DY(WC%NOM_IB(2))
               BYS(I,K) = (DY_OTHER*HP(I,1,K) + DY(1)*H_OTHER)/(DY(1)+DY_OTHER) + WALL_WORK1(IW)
            CASE(-2)
               DY_OTHER = MESHES(NOM)%DY(WC%NOM_IB(2))
               BYF(I,K) = (DY_OTHER*HP(I,JBAR,K) + DY(JBAR)*H_OTHER)/(DY(JBAR)+DY_OTHER) + WALL_WORK1(IW)
            CASE( 3)
               DZ_OTHER = MESHES(NOM)%DZ(WC%NOM_IB(3))
               BZS(I,J) = (DZ_OTHER*HP(I,J,1) + DZ(1)*H_OTHER)/(DZ(1)+DZ_OTHER) + WALL_WORK1(IW)
            CASE(-3)
               DZ_OTHER = MESHES(NOM)%DZ(WC%NOM_IB(3))
               BZF(I,J) = (DZ_OTHER*HP(I,J,KBAR) + DZ(KBAR)*H_OTHER)/(DZ(KBAR)+DZ_OTHER) + WALL_WORK1(IW)
         END SELECT
      
      ENDIF INTERPOLATED_ONLY
      
      ! Interpolation for embedded meshes --- set H on the boundary of mesh level 1 to an interpolated value from mesh level 0
      
      EMBEDDED_MESH_IF: IF (PERIODIC_TEST==8 .AND. WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .AND. MESH_LEVEL==1) THEN
         
         NOM = WC%NOM
         
         CALL GET_IJK(XC(I)-TWO_EPSILON_EB,YC(J)-TWO_EPSILON_EB,ZC(K)-TWO_EPSILON_EB,NOM,XIO,YJO,ZKO,IIO,JJO,KKO)
         
         IIX = FLOOR(XIO+.5_EB)
         JJY = FLOOR(YJO+.5_EB)
         KKZ = FLOOR(ZKO+.5_EB)
         
         SELECT CASE(IOR)
            CASE( 1)
               PP = (X(I)  - MESHES(NOM)%XC(IIX)) * MESHES(NOM)%RDX(IIX)
               RR = (YC(J) - MESHES(NOM)%YC(JJY)) * MESHES(NOM)%RDY(JJY)
               SS = (ZC(K) - MESHES(NOM)%ZC(KKZ)) * MESHES(NOM)%RDZ(KKZ)
               
               !IF (PREDICTOR) BXS(J,K) = AFILL2( MESHES(NOM)%H,IIX,JJY,KKZ,PP,RR,SS) + WALL_WORK1(IW)
               !IF (CORRECTOR) BXS(J,K) = AFILL2(MESHES(NOM)%HS,IIX,JJY,KKZ,PP,RR,SS) + WALL_WORK1(IW)
               BXS(J,K) = NS_H_EXACT(X(I),ZC(K),T,MU(I,J,K),WC%RHO_F,2._EB)
            CASE(-1)
               PP = (X(I - 1) - MESHES(NOM)%XC(IIX)) * MESHES(NOM)%RDX(IIX)
               RR = (YC(J)    - MESHES(NOM)%YC(JJY)) * MESHES(NOM)%RDY(JJY)
               SS = (ZC(K)    - MESHES(NOM)%ZC(KKZ)) * MESHES(NOM)%RDZ(KKZ)
               
               !IF (PREDICTOR) BXF(J,K) = AFILL2( MESHES(NOM)%H,IIX,JJY,KKZ,PP,RR,SS) + WALL_WORK1(IW)
               !IF (CORRECTOR) BXF(J,K) = AFILL2(MESHES(NOM)%HS,IIX,JJY,KKZ,PP,RR,SS) + WALL_WORK1(IW)
               BXS(J,K) = NS_H_EXACT(X(I-1),ZC(K),T,MU(I,J,K),WC%RHO_F,2._EB)
            CASE( 2)
               PP = (XC(I) - MESHES(NOM)%XC(IIX)) * MESHES(NOM)%RDX(IIX)
               RR = (Y(J)  - MESHES(NOM)%YC(JJY)) * MESHES(NOM)%RDY(JJY)
               SS = (ZC(K) - MESHES(NOM)%ZC(KKZ)) * MESHES(NOM)%RDZ(KKZ)
               
               !IF (PREDICTOR) BYS(I,K) = AFILL2( MESHES(NOM)%H,IIX,JJY,KKZ,PP,RR,SS) + WALL_WORK1(IW)
               !IF (CORRECTOR) BYS(I,K) = AFILL2(MESHES(NOM)%HS,IIX,JJY,KKZ,PP,RR,SS) + WALL_WORK1(IW)
               BYS(I,K) = 0._EB
            CASE(-2)
               PP = (XC(I)    - MESHES(NOM)%XC(IIX)) * MESHES(NOM)%RDX(IIX)
               RR = (Y(J - 1) - MESHES(NOM)%YC(JJY)) * MESHES(NOM)%RDY(JJY)
               SS = (ZC(K)    - MESHES(NOM)%ZC(KKZ)) * MESHES(NOM)%RDZ(KKZ)
               
               !IF (PREDICTOR) BYF(I,K) = AFILL2( MESHES(NOM)%H,IIX,JJY,KKZ,PP,RR,SS) + WALL_WORK1(IW)
               !IF (CORRECTOR) BYF(I,K) = AFILL2(MESHES(NOM)%HS,IIX,JJY,KKZ,PP,RR,SS) + WALL_WORK1(IW)
               BYF(I,K) = 0._EB
            CASE( 3)
               PP = (XC(I) - MESHES(NOM)%XC(IIX)) * MESHES(NOM)%RDX(IIX)
               RR = (YC(J) - MESHES(NOM)%YC(JJY)) * MESHES(NOM)%RDY(JJY)
               SS = (Z(K)  - MESHES(NOM)%ZC(KKZ)) * MESHES(NOM)%RDZ(KKZ)
               
               !IF (PREDICTOR) BZS(I,J)  = AFILL2( MESHES(NOM)%H,IIX,JJY,KKZ,PP,RR,SS) + WALL_WORK1(IW)
               !IF (CORRECTOR) BZS(I,J)  = AFILL2(MESHES(NOM)%HS,IIX,JJY,KKZ,PP,RR,SS) + WALL_WORK1(IW)
               BZS(I,J) = NS_H_EXACT(XC(I),Z(K),T,MU(I,J,K),WC%RHO_F,2._EB)
            CASE(-3)
               PP = (XC(I)    - MESHES(NOM)%XC(IIX)) * MESHES(NOM)%RDX(IIX)
               RR = (YC(J)    - MESHES(NOM)%YC(JJY)) * MESHES(NOM)%RDY(JJY)
               SS = (Z(K - 1) - MESHES(NOM)%ZC(KKZ)) * MESHES(NOM)%RDZ(KKZ)
               
               !IF (PREDICTOR) BZF(I,J) = AFILL2( MESHES(NOM)%H,IIX,JJY,KKZ,PP,RR,SS) + WALL_WORK1(IW)
               !IF (CORRECTOR) BZF(I,J) = AFILL2(MESHES(NOM)%HS,IIX,JJY,KKZ,PP,RR,SS) + WALL_WORK1(IW)
               BZS(I,J) = NS_H_EXACT(XC(I),Z(K-1),T,MU(I,J,K),WC%RHO_F,2._EB)
         END SELECT
      
      ENDIF EMBEDDED_MESH_IF
      
      ! OPEN (passive opening to exterior of domain) boundary. Apply inflow/outflow BC.

      OPEN_IF: IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
 
         IF (WC%VENT_INDEX>0) THEN
            VT => VENTS(WC%VENT_INDEX)
            IF (ABS(WC%ONE_D%T_IGN-T_BEGIN)<=TWO_EPSILON_EB .AND. VT%PRESSURE_RAMP_INDEX >=1) THEN
               TSI = T
            ELSE
               TSI = T - T_BEGIN
            ENDIF
            TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,DUMMY,VT%PRESSURE_RAMP_INDEX)
            P_EXTERNAL = TIME_RAMP_FACTOR*VT%DYNAMIC_PRESSURE      
         ENDIF

         IF (ANY(MEAN_FORCING)) THEN
            UBAR = U0*EVALUATE_RAMP(T-T_BEGIN,DUMMY,I_RAMP_U0)*EVALUATE_RAMP(ZC(K),DUMMY,I_RAMP_U0_Z)
            VBAR = V0*EVALUATE_RAMP(T-T_BEGIN,DUMMY,I_RAMP_V0)*EVALUATE_RAMP(ZC(K),DUMMY,I_RAMP_V0_Z)
            WBAR = W0*EVALUATE_RAMP(T-T_BEGIN,DUMMY,I_RAMP_W0)*EVALUATE_RAMP(ZC(K),DUMMY,I_RAMP_W0_Z)
            H0 = 0.5_EB*(UBAR**2+VBAR**2+WBAR**2)
         ENDIF

         SELECT CASE(IOR)
            CASE( 1)
               IF (UU(0,J,K)<0._EB) THEN
                  BXS(J,K) = P_EXTERNAL/WC%RHO_F + KRES(1,J,K)
               ELSE
                  BXS(J,K) = P_EXTERNAL/WC%RHO_F + H0
               ENDIF
            CASE(-1)
               IF (UU(IBAR,J,K)>0._EB) THEN
                  BXF(J,K) = P_EXTERNAL/WC%RHO_F + KRES(IBAR,J,K)
               ELSE
                  BXF(J,K) = P_EXTERNAL/WC%RHO_F + H0
               ENDIF
            CASE( 2)
               IF (VV(I,0,K)<0._EB) THEN
                  BYS(I,K) = P_EXTERNAL/WC%RHO_F + KRES(I,1,K)
               ELSE
                  BYS(I,K) = P_EXTERNAL/WC%RHO_F + H0
               ENDIF
            CASE(-2)
               IF (VV(I,JBAR,K)>0._EB) THEN
                  BYF(I,K) = P_EXTERNAL/WC%RHO_F + KRES(I,JBAR,K)
               ELSE
                  BYF(I,K) = P_EXTERNAL/WC%RHO_F + H0
               ENDIF
            CASE( 3)
               IF (WW(I,J,0)<0._EB) THEN
                  BZS(I,J) = P_EXTERNAL/WC%RHO_F + KRES(I,J,1)
               ELSE
                  BZS(I,J) = P_EXTERNAL/WC%RHO_F + H0
               ENDIF
            CASE(-3)
               IF (WW(I,J,KBAR)>0._EB) THEN
                  BZF(I,J) = P_EXTERNAL/WC%RHO_F + KRES(I,J,KBAR)
               ELSE
                  BZF(I,J) = P_EXTERNAL/WC%RHO_F + H0
               ENDIF
         END SELECT
    
      ENDIF OPEN_IF
 
   ENDIF IF_DIRICHLET
 
ENDDO WALL_CELL_LOOP

! Compute the RHS of the Poisson equation

SELECT CASE(IPS)

   CASE(:1,4,7)
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
      IF (.NOT.CYLINDRICAL) THEN
         !$OMP PARALLEL DO PRIVATE(TRM1, TRM2, TRM3, TRM4) SCHEDULE(STATIC)
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
         !$OMP END PARALLEL DO
      ENDIF
 
   CASE(2)  ! Switch x and y
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
 
   CASE(3,6)  ! Switch x and z
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
 
   CASE(5)  ! Switch y and z
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

END SELECT


! In case of ScaRC-method leave routine

IF (PRES_METHOD == 'SCARC') RETURN

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



SELECT CASE(IPS)
   CASE(:1,4,7)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               HP(I,J,K) = PRHS(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   CASE(2)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               HP(I,J,K) = PRHS(J,I,K)
            ENDDO
         ENDDO
      ENDDO
   CASE(3,6)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               HP(I,J,K) = PRHS(K,J,I)
            ENDDO
         ENDDO
      ENDDO
   CASE(5)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               HP(I,J,K) = PRHS(I,K,J)
            ENDDO
         ENDDO
      ENDDO
END SELECT 

! Apply boundary conditions to H

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



! Optional check of the accuracy of the pressure solver

IF (CHECK_POISSON .AND. .NOT.EVACUATION_ONLY(NM)) THEN
   POIS_ERR = 0._EB
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
USE GLOBAL_CONSTANTS, ONLY: PREDICTOR,VELOCITY_ERROR_MAX,SOLID_BOUNDARY,INTERPOLATED_BOUNDARY,HVAC_BOUNDARY,&
                            VELOCITY_ERROR_MAX_I,VELOCITY_ERROR_MAX_J,VELOCITY_ERROR_MAX_K,TUSED

INTEGER, INTENT(IN) :: NM
INTEGER :: IW,IOR,II,JJ,KK,IIO,JJO,KKO,N_INT_CELLS
REAL(EB) :: TNOW,UN_NEW,UN_NEW_OTHER,VELOCITY_ERROR,DUDT,DVDT,DWDT,ITERATIVE_FACTOR
TYPE(OMESH_TYPE), POINTER :: OM
TYPE(MESH_TYPE), POINTER :: M2
TYPE(WALL_TYPE), POINTER :: WC

TNOW=SECOND()
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   ITERATIVE_FACTOR = 0.25_EB
ELSE
   ITERATIVE_FACTOR = 0.50_EB
ENDIF

VELOCITY_ERROR_MAX(NM) = 0._EB
WALL_WORK1 = 0._EB

CHECK_WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

   WC=>WALL(IW)

   IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY        .AND. &
       WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. &
       WC%BOUNDARY_TYPE/=HVAC_BOUNDARY) CYCLE CHECK_WALL_LOOP

   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
      OM => OMESH(WC%NOM)
      M2 => MESHES(WC%NOM)
   ENDIF

   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR

   ! Update normal component of velocity at the mesh boundary

   IF (PREDICTOR) THEN
      SELECT CASE(IOR)
         CASE( 1)
            UN_NEW = U(II,JJ,KK)   - DT*(FVX(II,JJ,KK)   + RDXN(II)  *(H(II+1,JJ,KK)-H(II,JJ,KK)))
         CASE(-1)
            UN_NEW = U(II-1,JJ,KK) - DT*(FVX(II-1,JJ,KK) + RDXN(II-1)*(H(II,JJ,KK)-H(II-1,JJ,KK)))
         CASE( 2)
            UN_NEW = V(II,JJ,KK)   - DT*(FVY(II,JJ,KK)   + RDYN(JJ)  *(H(II,JJ+1,KK)-H(II,JJ,KK)))
         CASE(-2)
            UN_NEW = V(II,JJ-1,KK) - DT*(FVY(II,JJ-1,KK) + RDYN(JJ-1)*(H(II,JJ,KK)-H(II,JJ-1,KK)))
         CASE( 3)
            UN_NEW = W(II,JJ,KK)   - DT*(FVZ(II,JJ,KK)   + RDZN(KK)  *(H(II,JJ,KK+1)-H(II,JJ,KK)))
         CASE(-3)
            UN_NEW = W(II,JJ,KK-1) - DT*(FVZ(II,JJ,KK-1) + RDZN(KK-1)*(H(II,JJ,KK)-H(II,JJ,KK-1)))
      END SELECT
   ELSE
      SELECT CASE(IOR)
         CASE( 1)
            UN_NEW = 0.5_EB*(U(II,JJ,KK)+US(II,JJ,KK)     - DT*(FVX(II,JJ,KK)   + RDXN(II)  *(HS(II+1,JJ,KK)-HS(II,JJ,KK))))
         CASE(-1)
            UN_NEW = 0.5_EB*(U(II-1,JJ,KK)+US(II-1,JJ,KK) - DT*(FVX(II-1,JJ,KK) + RDXN(II-1)*(HS(II,JJ,KK)-HS(II-1,JJ,KK))))
         CASE( 2)
            UN_NEW = 0.5_EB*(V(II,JJ,KK)+VS(II,JJ,KK)     - DT*(FVY(II,JJ,KK)   + RDYN(JJ)  *(HS(II,JJ+1,KK)-HS(II,JJ,KK))))
         CASE(-2)
            UN_NEW = 0.5_EB*(V(II,JJ-1,KK)+VS(II,JJ-1,KK) - DT*(FVY(II,JJ-1,KK) + RDYN(JJ-1)*(HS(II,JJ,KK)-HS(II,JJ-1,KK))))
         CASE( 3)
            UN_NEW = 0.5_EB*(W(II,JJ,KK)+WS(II,JJ,KK)     - DT*(FVZ(II,JJ,KK)   + RDZN(KK)  *(HS(II,JJ,KK+1)-HS(II,JJ,KK))))
         CASE(-3)
            UN_NEW = 0.5_EB*(W(II,JJ,KK-1)+WS(II,JJ,KK-1) - DT*(FVZ(II,JJ,KK-1) + RDZN(KK-1)*(HS(II,JJ,KK)-HS(II,JJ,KK-1))))
      END SELECT
   ENDIF

   ! At interpolated boundaries, compare updated normal component of velocity with that of the other mesh

   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN

      UN_NEW_OTHER = 0._EB

      DO KKO=WC%NOM_IB(3),WC%NOM_IB(6)
         DO JJO=WC%NOM_IB(2),WC%NOM_IB(5)
            DO IIO=WC%NOM_IB(1),WC%NOM_IB(4)
               IF (PREDICTOR) THEN
                  SELECT CASE(IOR)
                     CASE( 1)
                        DUDT = -OM%FVX(IIO,JJO,KKO)   - M2%RDXN(IIO)  *(OM%H(IIO+1,JJO,KKO)-OM%H(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%U(IIO,JJO,KKO)   + DT*DUDT
                     CASE(-1)
                        DUDT = -OM%FVX(IIO-1,JJO,KKO) - M2%RDXN(IIO-1)*(OM%H(IIO,JJO,KKO)-OM%H(IIO-1,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%U(IIO-1,JJO,KKO) + DT*DUDT
                     CASE( 2)
                        DVDT = -OM%FVY(IIO,JJO,KKO)   - M2%RDYN(JJO)  *(OM%H(IIO,JJO+1,KKO)-OM%H(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%V(IIO,JJO,KKO)   + DT*DVDT
                     CASE(-2)
                        DVDT = -OM%FVY(IIO,JJO-1,KKO) - M2%RDYN(JJO-1)*(OM%H(IIO,JJO,KKO)-OM%H(IIO,JJO-1,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%V(IIO,JJO-1,KKO) + DT*DVDT
                     CASE( 3)
                        DWDT = -OM%FVZ(IIO,JJO,KKO)   - M2%RDZN(KKO)  *(OM%H(IIO,JJO,KKO+1)-OM%H(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%W(IIO,JJO,KKO)   + DT*DWDT
                     CASE(-3)
                        DWDT = -OM%FVZ(IIO,JJO,KKO-1) - M2%RDZN(KKO-1)*(OM%H(IIO,JJO,KKO)-OM%H(IIO,JJO,KKO-1))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%W(IIO,JJO,KKO-1) + DT*DWDT
                  END SELECT
               ELSE
                  SELECT CASE(IOR)
                     CASE( 1)
                        DUDT = -OM%FVX(IIO,JJO,KKO)   - M2%RDXN(IIO)  *(OM%HS(IIO+1,JJO,KKO)-OM%HS(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%U(IIO,JJO,KKO)+OM%US(IIO,JJO,KKO)     + DT*DUDT)
                     CASE(-1)
                        DUDT = -OM%FVX(IIO-1,JJO,KKO) - M2%RDXN(IIO-1)*(OM%HS(IIO,JJO,KKO)-OM%HS(IIO-1,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%U(IIO-1,JJO,KKO)+OM%US(IIO-1,JJO,KKO) + DT*DUDT)
                     CASE( 2)
                        DVDT = -OM%FVY(IIO,JJO,KKO)   - M2%RDYN(JJO)  *(OM%HS(IIO,JJO+1,KKO)-OM%HS(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%V(IIO,JJO,KKO)+OM%VS(IIO,JJO,KKO)     + DT*DVDT)
                     CASE(-2)
                        DVDT = -OM%FVY(IIO,JJO-1,KKO) - M2%RDYN(JJO-1)*(OM%HS(IIO,JJO,KKO)-OM%HS(IIO,JJO-1,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%V(IIO,JJO-1,KKO)+OM%VS(IIO,JJO-1,KKO) + DT*DVDT)
                     CASE( 3)
                        DWDT = -OM%FVZ(IIO,JJO,KKO)   - M2%RDZN(KKO)  *(OM%HS(IIO,JJO,KKO+1)-OM%HS(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%W(IIO,JJO,KKO)+OM%WS(IIO,JJO,KKO)     + DT*DWDT)
                     CASE(-3)
                        DWDT = -OM%FVZ(IIO,JJO,KKO-1) - M2%RDZN(KKO-1)*(OM%HS(IIO,JJO,KKO)-OM%HS(IIO,JJO,KKO-1))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%W(IIO,JJO,KKO-1)+OM%WS(IIO,JJO,KKO-1) + DT*DWDT)
                  END SELECT
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      N_INT_CELLS  = (WC%NOM_IB(4)-WC%NOM_IB(1)+1) * (WC%NOM_IB(5)-WC%NOM_IB(2)+1) * (WC%NOM_IB(6)-WC%NOM_IB(3)+1)
      UN_NEW_OTHER = UN_NEW_OTHER/REAL(N_INT_CELLS,EB)

   ENDIF

   ! At solid boundaries, compare updated normal velocity with specified normal velocity 

   IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .OR. WC%BOUNDARY_TYPE==HVAC_BOUNDARY) THEN
      IF (PREDICTOR) THEN
         UN_NEW_OTHER = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UWS
      ELSE
         UN_NEW_OTHER = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%UW
      ENDIF
   ENDIF

   ! Compute velocity difference

   VELOCITY_ERROR = UN_NEW - UN_NEW_OTHER
   WC%VEL_ERR_NEW = VELOCITY_ERROR
   WALL_WORK1(IW) = -SIGN(1._EB,REAL(IOR,EB))*ITERATIVE_FACTOR*VELOCITY_ERROR/(WC%RDN*DT)

   ! If the grid cells in the current mesh are smaller than those of the other mesh, do not include in error tolerance

   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
      IF (OM%NIC_S>OM%NIC_R) CYCLE CHECK_WALL_LOOP
   ENDIF

   ! Save maximum velocity error

   IF (ABS(VELOCITY_ERROR)>VELOCITY_ERROR_MAX(NM)) THEN
      VELOCITY_ERROR_MAX_I(NM) = II
      VELOCITY_ERROR_MAX_J(NM) = JJ
      VELOCITY_ERROR_MAX_K(NM) = KK
      VELOCITY_ERROR_MAX(NM)   = ABS(VELOCITY_ERROR)
   ENDIF

ENDDO CHECK_WALL_LOOP

TUSED(5,NM)=TUSED(5,NM)+SECOND()-TNOW
END SUBROUTINE COMPUTE_VELOCITY_ERROR


! SUBROUTINE BUILD_SPARSE_MATRIX_LAPLACE(NM)
! USE GLOBAL_CONSTANTS, ONLY: FISHPAK_BC_NEUMANN_NEUMANN, FISHPAK_BC_NEUMANN_DIRICHLET, FISHPAK_BC_DIRICHLET_NEUMANN, &
!                             FISHPAK_BC_DIRICHLET_DIRICHLET
! USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR 
! IMPLICIT NONE

! INTEGER, INTENT(IN) :: NM
! INTEGER :: IZERO,I,J,K,I_R,I_S,I_S_D,NX,NY,NZ,NXNY
! REAL(EB) :: DXDX,DYDY,DZDZ,BC(6),BCW,BCE,BCS,BCN,BCB,BCT
! LOGICAL :: RECORD_ROW_INDEX
! INTEGER, ALLOCATABLE, DIMENSION (:) :: A_ROWS_TMP,A_COLUMNS_TMP
! REAL(EB), ALLOCATABLE, DIMENSION (:) :: A_VALUES_TMP
! REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: FULL_A
! TYPE (MESH_TYPE), POINTER :: M=>NULL()

! M=>MESHES(NM)

! M%A_N_ROWS = M%IBAR*M%JBAR*M%KBAR
! M%A_N_COLS = M%A_N_ROWS
! M%A_N_ELEMENTS = 7*M%A_N_ROWS ! high estimate

! ! (debug) allocate memory for full A

! ALLOCATE(FULL_A(M%A_N_ROWS,M%A_N_COLS),STAT=IZERO); CALL ChkMemErr('PRES','FULL_A',IZERO)

! ! Allocate memory for compressed row format

! ALLOCATE(M%A_ROW_INDEX(M%A_N_ROWS+1),STAT=IZERO); CALL ChkMemErr('PRES','A_ROW_INDEX',IZERO)

! ! Allocate temporary memory for sparse matrix A

! ALLOCATE(A_ROWS_TMP(M%A_N_ELEMENTS),STAT=IZERO); CALL ChkMemErr('PRES','A_ROWS_TMP',IZERO)
! ALLOCATE(A_COLUMNS_TMP(M%A_N_ELEMENTS),STAT=IZERO); CALL ChkMemErr('PRES','A_COLUMNS_TMP',IZERO)
! ALLOCATE(A_VALUES_TMP(M%A_N_ELEMENTS),STAT=IZERO); CALL ChkMemErr('PRES','A_VALUES_TMP',IZERO)

! ! Set boundary conditions on exterior of domain

! SELECT CASE(M%LBC)
!    CASE(FISHPAK_BC_NEUMANN_NEUMANN)
!       BC(1)=0._EB
!       BC(2)=0._EB
!    CASE(FISHPAK_BC_NEUMANN_DIRICHLET)
!       BC(1)=0._EB
!       BC(2)=2._EB
!    CASE(FISHPAK_BC_DIRICHLET_NEUMANN)
!       BC(1)=2._EB
!       BC(2)=0._EB
!    CASE(FISHPAK_BC_DIRICHLET_DIRICHLET)
!       BC(1)=2._EB
!       BC(2)=2._EB
! END SELECT

! SELECT CASE(M%MBC)
!    CASE(FISHPAK_BC_NEUMANN_NEUMANN)
!       BC(3)=0._EB
!       BC(4)=0._EB
!    CASE(FISHPAK_BC_NEUMANN_DIRICHLET)
!       BC(3)=0._EB
!       BC(4)=2._EB
!    CASE(FISHPAK_BC_DIRICHLET_NEUMANN)
!       BC(3)=2._EB
!       BC(4)=0._EB
!    CASE(FISHPAK_BC_DIRICHLET_DIRICHLET)
!       BC(3)=2._EB
!       BC(4)=2._EB
! END SELECT

! SELECT CASE(M%NBC)
!    CASE(FISHPAK_BC_NEUMANN_NEUMANN)
!       BC(5)=0._EB
!       BC(6)=0._EB
!    CASE(FISHPAK_BC_NEUMANN_DIRICHLET)
!       BC(5)=0._EB
!       BC(6)=2._EB
!    CASE(FISHPAK_BC_DIRICHLET_NEUMANN)
!       BC(5)=2._EB
!       BC(6)=0._EB
!    CASE(FISHPAK_BC_DIRICHLET_DIRICHLET)
!       BC(5)=2._EB
!       BC(6)=2._EB
! END SELECT

! ! Build sparse matrix A

! ! boundary condition type
! ! -----------------------
! ! bc(1) = west
! ! bc(2) = east
! ! bc(3) = south
! ! bc(4) = north
! ! bc(5) = bottom
! ! bc(6) = top
! ! values: 0 = Neumann, 1 = Dirichlet (ghost), 2 = Dirichlet (face)
! !
! ! Eventually, we will allow BXS(I,J), etc., for a specific type and value of bc for each wall cell.
! !
! ! Consider that we are solving the following discrete equation in 2D:
! !
! ! (p(i+1,j)-2*p(i,j)+p(i-1,j))/dx^2 + (p(i,j+1)-2*p(i,j)+p(i,j-1))^2 = b(i,j)
! !
! ! Example of Neumann bc:
! ! Suppose we have dp/dx = Fx on the right-side bc.  The we have
! ! p(i+1)-p(i,j)= Fx*dx, and the eqn is rewritten as
! !
! ! ( Fx*dx - p(i,j)+p(i-1,j) )/dx^2 + ... = b(i,j).
! !
! ! In this case, the coefficient for p(i+1,j) is 0, and the coefficient
! ! for p(i,j) is -1.  And the source is augmented to be
! ! b(i,j) - Fx/dx.

! DXDX = M%DX(1)**2 ! this setup assumes uniform grid spacing in each direction
! DYDY = M%DY(1)**2
! DZDZ = M%DZ(1)**2

! NX = M%IBAR
! NY = M%JBAR
! NZ = M%KBAR
! NXNY = NX*NY

! ! Build A

! I_R = 0 ! row index, lexicographical ordering, I_R(I,J,K) = (K-1)*NXNY+(J-1)*NX+I
! I_S = 0 ! sparse element index

! ! note the ordering of the coefficients for the discrete Laplacian
! !
! !     phi(i,j,k-1)/dz^2 + phi(i,j-1,k)/dy^2 + phi(i-1,j,k)/dx^2
! !  - (2*phi(i,j,k)/dx^2 + 2*phi(i,j,k)/dy^2 + 2*phi(i,j,k)/dz^2)
! !     phi(i+1,j,k)/dx^2 + phi(i,j+1,k)/dy^2 + phi(i,j,k+1)/dz^2
! !
! ! within a row the order is always
! !
! !    (k-1) (j-1) (i-1) ijk (i+1) (j+1) (k+1)
! !
! ! the matrix is symmetric, so we only store ijk, i+1, j+1, k+1
! ! that is, we store the upper triangular part of A

! DO K=1,NZ
!    DO J=1,NY
!       DO I=1,NX

!          I_R = (K-1)*NXNY+(J-1)*NX+I

!          RECORD_ROW_INDEX = .TRUE.

!          IF (K>1) THEN
!             !I_S = I_S + 1
!             !IF (RECORD_ROW_INDEX) THEN
!             !   M%A_ROW_INDEX(I_R) = I_S
!             !   RECORD_ROW_INDEX = .FALSE.
!             !ENDIF
!             !A_COLUMNS_TMP(I_S) = I_R-NXNY
!             BCB = 1._EB
!             !A_VALUES_TMP(I_S) = BCB/DZDZ
!          ELSE
!             BCB = BC(5)
!          ENDIF

!          IF (J>1) THEN
!             !I_S = I_S + 1
!             !IF (RECORD_ROW_INDEX) THEN
!             !   M%A_ROW_INDEX(I_R) = I_S
!             !   RECORD_ROW_INDEX = .FALSE.
!             !ENDIF
!             !A_COLUMNS_TMP(I_S) = I_R-NX
!             BCS = 1._EB
!             !A_VALUES_TMP(I_S) = BCS/DYDY
!          ELSE
!             BCS = BC(3)
!          ENDIF

!          IF (I>1) THEN
!             !I_S = I_S + 1
!             !IF (RECORD_ROW_INDEX) THEN
!             !   M%A_ROW_INDEX(I_R) = I_S
!             !   RECORD_ROW_INDEX = .FALSE.
!             !ENDIF
!             !A_COLUMNS_TMP(I_S) = I_R-1
!             BCW = 1._EB
!             !A_VALUES_TMP(I_S) = BCW/DXDX
!          ELSE
!             BCW = BC(1)
!          ENDIF

!          ! increment I_S for diagonal coefficient, but wait to compute until bcs are known
!          I_S = I_S + 1
!          IF (RECORD_ROW_INDEX) THEN
!             M%A_ROW_INDEX(I_R) = I_S
!             RECORD_ROW_INDEX = .FALSE.
!          ENDIF
!          I_S_D = I_S

!          IF (I<NX) THEN
!             I_S = I_S + 1
!             IF (RECORD_ROW_INDEX) THEN
!                M%A_ROW_INDEX(I_R) = I_S
!                RECORD_ROW_INDEX = .FALSE.
!             ENDIF
!             A_COLUMNS_TMP(I_S) = I_R+1
!             BCE = 1._EB
!             A_VALUES_TMP(I_S) = BCE/DXDX
!          ELSE
!             BCE = BC(2)
!          ENDIF

!          IF (J<NY) THEN
!             I_S = I_S + 1
!             IF (RECORD_ROW_INDEX) THEN
!                M%A_ROW_INDEX(I_R) = I_S
!                RECORD_ROW_INDEX = .FALSE.
!             ENDIF
!             A_COLUMNS_TMP(I_S) = I_R+NX
!             BCN = 1._EB
!             A_VALUES_TMP(I_S) = BCN/DYDY
!          ELSE
!             BCN = BC(4)
!          ENDIF

!          IF (K<NZ) THEN
!             I_S = I_S + 1
!             IF (RECORD_ROW_INDEX) THEN
!                M%A_ROW_INDEX(I_R) = I_S
!                RECORD_ROW_INDEX = .FALSE.
!             ENDIF
!             A_COLUMNS_TMP(I_S) = I_R+NXNY
!             BCT = 1._EB
!             A_VALUES_TMP(I_S) = BCT/DZDZ
!          ELSE
!             BCT = BC(6)
!          ENDIF

!          ! fill in I_S for diagonal
!          A_COLUMNS_TMP(I_S_D) = I_R
!          A_VALUES_TMP(I_S_D) = -( (BCW+BCE)/DXDX + (BCS+BCN)/DYDY + (BCB+BCT)/DZDZ )

!       ENDDO
!    ENDDO
! ENDDO

! ! Allocate final storage arrays

! M%A_N_ELEMENTS = I_S ! exact element count

! ALLOCATE(M%A_COLUMNS(M%A_N_ELEMENTS),STAT=IZERO); CALL ChkMemErr('PRES','A_COLUMNS',IZERO)
! ALLOCATE(M%A_VALUES(M%A_N_ELEMENTS),STAT=IZERO); CALL ChkMemErr('PRES','A_VALUES',IZERO)

! ! Copy temporary values for final arrays

! M%A_COLUMNS = A_COLUMNS_TMP(1:M%A_N_ELEMENTS)
! M%A_VALUES = A_VALUES_TMP(1:M%A_N_ELEMENTS)

! ! Deallocate temporary arrays

! ! DEALLOCATE(A_COLUMNS_TMP)
! ! DEALLOCATE(A_VALUES_TMP)

! ! Fill in final row index

! M%A_ROW_INDEX(M%A_N_ROWS+1) = M%A_N_ELEMENTS+1

! ! print *,M%A_VALUES
! ! print *,M%A_COLUMNS
! ! print *,M%A_ROW_INDEX
! ! print *,M%A_N_ROWS,M%A_N_ELEMENTS,I_S
! DO I_R=1,M%A_N_ROWS
!    I_S = M%A_ROW_INDEX(I_R)
!    print *,M%A_VALUES(I_S)
! ENDDO

! END SUBROUTINE BUILD_SPARSE_MATRIX_LAPLACE


! SUBROUTINE SPARSE_LU_FACTORIZATION(NM)
! USE GLOBAL_CONSTANTS, ONLY: SCI_KM1, SCI_JM1, SCI_IM1, SCI_IJK, SCI_IP1,SCI_JP1,SCI_KP1, NDIAG
! IMPLICIT NONE

! INTEGER, INTENT(IN) :: NM
! INTEGER :: IZERO,N,N_L,I,J,K,K_S
! TYPE (MESH_TYPE), POINTER :: M=>NULL()

! M=>MESHES(NM)
! N = M%IBAR*M%JBAR*M%KBAR

! ! Go through algorithm once just to count nonzero entries

! ! Heath, Algorithm 2.3, LU Factorization by Gaussian Elimination

! DO K=1,N-1
!    K_S_D = (K-1)*NDIAG
!    K_S = K_S_D + SCI_IJK
!    IF ( ABS(M%A_COEF(K_S))<TWO_EPSILON_EB ) CALL SHUTDOWN('STOP: Error in LU factorization')
!    DO I=K+1,N
!       I_S_D = (I-1)*NDIAG
!       DO J=1:NDIAG
!          I_S = I_S_D + J
!          L_IK = A(I,K)/A(K,K)
!       ENDDO
!    ENDDO

! ENDDO

! ! Allocate 1D arrays

! ! ALLOCATE(M%L_I(N_L),STAT=IZERO);    CALL ChkMemErr('INIT','L_I',IZERO)
! ! ALLOCATE(M%L_J(N_L),STAT=IZERO);    CALL ChkMemErr('INIT','L_J',IZERO)
! ! ALLOCATE(M%L_COEF(N_L),STAT=IZERO); CALL ChkMemErr('INIT','L_COEF',IZERO)

! END SUBROUTINE SPARSE_LU_FACTORIZATION


SUBROUTINE GET_REV_pres(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') presrev(INDEX(presrev,':')+2:LEN_TRIM(presrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') presdate

END SUBROUTINE GET_REV_pres


END MODULE PRES













