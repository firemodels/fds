MODULE PRES

! Find the perturbation pressure by solving Poisson's Equation

USE PRECISION_PARAMETERS
USE MESH_POINTERS

IMPLICIT NONE
PRIVATE

PUBLIC PRESSURE_SOLVER_COMPUTE_RHS,PRESSURE_SOLVER_FFT,PRESSURE_SOLVER_CHECK_RESIDUALS,COMPUTE_VELOCITY_ERROR

CONTAINS

SUBROUTINE PRESSURE_SOLVER_COMPUTE_RHS(T,NM)

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE GLOBAL_CONSTANTS

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,HP,RHOP
INTEGER :: I,J,K,IW,IOR,NOM,N_INT_CELLS,IIO,JJO,KKO
REAL(EB) :: TRM1,TRM2,TRM3,TRM4,H_OTHER,TNOW,DUMMY=0._EB, &
            TSI,TIME_RAMP_FACTOR,DX_OTHER,DY_OTHER,DZ_OTHER,P_EXTERNAL, &
            UBAR,VBAR,WBAR
TYPE (VENTS_TYPE), POINTER :: VT
TYPE (WALL_TYPE), POINTER :: WC
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   HP => H
   RHOP => RHO
ELSE
   UU => US
   VV => VS
   WW => WS
   HP => HS
   RHOP => RHOS
ENDIF

! Apply pressure boundary conditions at external cells.
! If Neumann, BXS, BXF, etc., contain dH/dx(x=XS), dH/dx(x=XF), etc.
! If Dirichlet, BXS, BXF, etc., contain H(x=XS), H(x=XF), etc.
! LBC, MBC and NBC are codes used be Poisson solver to denote type
! of boundary condition at x, y and z boundaries. See Crayfishpak
! manual for details.

WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

   WC => WALL(IW)
   EWC => EXTERNAL_WALL(IW)
   I   = WC%ONE_D%II
   J   = WC%ONE_D%JJ
   K   = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR

   ! Apply pressure gradients at NEUMANN boundaries: dH/dn = -F_n - d(u_n)/dt

   IF_NEUMANN: IF (WC%PRESSURE_BC_INDEX==NEUMANN) THEN
      SELECT CASE(IOR)
         CASE( 1)
            BXS(J,K) = HX(0)   *(-FVX(0,J,K)    + WC%DUNDT)
         CASE(-1)
            BXF(J,K) = HX(IBP1)*(-FVX(IBAR,J,K) - WC%DUNDT)
         CASE( 2)
            BYS(I,K) = HY(0)   *(-FVY(I,0,K)    + WC%DUNDT)
         CASE(-2)
            BYF(I,K) = HY(JBP1)*(-FVY(I,JBAR,K) - WC%DUNDT)
         CASE( 3)
            BZS(I,J) = HZ(0)   *(-FVZ(I,J,0)    + WC%DUNDT)
         CASE(-3)
            BZF(I,J) = HZ(KBP1)*(-FVZ(I,J,KBAR) - WC%DUNDT)
      END SELECT
   ENDIF IF_NEUMANN

   ! Apply pressures at DIRICHLET boundaries, depending on the specific type

   IF_DIRICHLET: IF (WC%PRESSURE_BC_INDEX==DIRICHLET) THEN

      NOT_OPEN: IF (WC%BOUNDARY_TYPE/=OPEN_BOUNDARY .AND. WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) THEN

         ! Solid boundary that uses a Dirichlet BC. Assume that the pressure at the boundary (BXS, etc) is the average of the
         ! last computed pressures in the ghost and adjacent gas cells.

         SELECT CASE(IOR)
            CASE( 1)
               BXS(J,K) = 0.5_EB*(HP(0,J,K)+HP(1,J,K)) + WALL_WORK1(IW)
            CASE(-1)
               BXF(J,K) = 0.5_EB*(HP(IBAR,J,K)+HP(IBP1,J,K)) + WALL_WORK1(IW)
            CASE( 2)
               BYS(I,K) = 0.5_EB*(HP(I,0,K)+HP(I,1,K)) + WALL_WORK1(IW)
            CASE(-2)
               BYF(I,K) = 0.5_EB*(HP(I,JBAR,K)+HP(I,JBP1,K)) + WALL_WORK1(IW)
            CASE( 3)
               BZS(I,J) = 0.5_EB*(HP(I,J,0)+HP(I,J,1)) + WALL_WORK1(IW)
            CASE(-3)
               BZF(I,J) = 0.5_EB*(HP(I,J,KBAR)+HP(I,J,KBP1)) + WALL_WORK1(IW)
         END SELECT

      ENDIF NOT_OPEN

      ! Interpolated boundary -- set boundary value of H to be average of neighboring cells from previous time step

      INTERPOLATED_ONLY: IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN

         NOM     = EWC%NOM
         H_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  IF (PREDICTOR) H_OTHER = H_OTHER + OMESH(NOM)%H(IIO,JJO,KKO)
                  IF (CORRECTOR) H_OTHER = H_OTHER + OMESH(NOM)%HS(IIO,JJO,KKO)
               ENDDO
            ENDDO
         ENDDO
         N_INT_CELLS = (EWC%IIO_MAX-EWC%IIO_MIN+1) * (EWC%JJO_MAX-EWC%JJO_MIN+1) * (EWC%KKO_MAX-EWC%KKO_MIN+1)
         H_OTHER = H_OTHER/REAL(N_INT_CELLS,EB)

         SELECT CASE(IOR)
            CASE( 1)
               DX_OTHER = MESHES(NOM)%DX(EWC%IIO_MIN)
               BXS(J,K) = (DX_OTHER*HP(1,J,K) + DX(1)*H_OTHER)/(DX(1)+DX_OTHER) + WALL_WORK1(IW)
            CASE(-1)
               DX_OTHER = MESHES(NOM)%DX(EWC%IIO_MIN)
               BXF(J,K) = (DX_OTHER*HP(IBAR,J,K) + DX(IBAR)*H_OTHER)/(DX(IBAR)+DX_OTHER) + WALL_WORK1(IW)
            CASE( 2)
               DY_OTHER = MESHES(NOM)%DY(EWC%JJO_MIN)
               BYS(I,K) = (DY_OTHER*HP(I,1,K) + DY(1)*H_OTHER)/(DY(1)+DY_OTHER) + WALL_WORK1(IW)
            CASE(-2)
               DY_OTHER = MESHES(NOM)%DY(EWC%JJO_MIN)
               BYF(I,K) = (DY_OTHER*HP(I,JBAR,K) + DY(JBAR)*H_OTHER)/(DY(JBAR)+DY_OTHER) + WALL_WORK1(IW)
            CASE( 3)
               DZ_OTHER = MESHES(NOM)%DZ(EWC%KKO_MIN)
               BZS(I,J) = (DZ_OTHER*HP(I,J,1) + DZ(1)*H_OTHER)/(DZ(1)+DZ_OTHER) + WALL_WORK1(IW)
            CASE(-3)
               DZ_OTHER = MESHES(NOM)%DZ(EWC%KKO_MIN)
               BZF(I,J) = (DZ_OTHER*HP(I,J,KBAR) + DZ(KBAR)*H_OTHER)/(DZ(KBAR)+DZ_OTHER) + WALL_WORK1(IW)
         END SELECT

      ENDIF INTERPOLATED_ONLY

      ! OPEN (passive opening to exterior of domain) boundary. Apply inflow/outflow BC.

      OPEN_IF: IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN

         VT => VENTS(WC%VENT_INDEX)
         IF (ABS(WC%ONE_D%T_IGN-T_BEGIN)<=TWO_EPSILON_EB .AND. VT%PRESSURE_RAMP_INDEX >=1) THEN
            TSI = T
         ELSE
            TSI = T - T_BEGIN
         ENDIF
         TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,DUMMY,VT%PRESSURE_RAMP_INDEX)
         P_EXTERNAL = TIME_RAMP_FACTOR*VT%DYNAMIC_PRESSURE

         IF (ANY(MEAN_FORCING)) THEN
            UBAR = U0*EVALUATE_RAMP(T,DUMMY,I_RAMP_U0_T)*EVALUATE_RAMP(ZC(K),DUMMY,I_RAMP_U0_Z)
            VBAR = V0*EVALUATE_RAMP(T,DUMMY,I_RAMP_V0_T)*EVALUATE_RAMP(ZC(K),DUMMY,I_RAMP_V0_Z)
            WBAR = W0*EVALUATE_RAMP(T,DUMMY,I_RAMP_W0_T)*EVALUATE_RAMP(ZC(K),DUMMY,I_RAMP_W0_Z)
            H0 = 0.5_EB*(UBAR**2+VBAR**2+WBAR**2)
         ENDIF

         SELECT CASE(IOR)
            CASE( 1)
               IF (UU(0,J,K)<0._EB) THEN
                  BXS(J,K) = P_EXTERNAL/WC%ONE_D%RHO_F + KRES(1,J,K)
               ELSE
                  BXS(J,K) = P_EXTERNAL/WC%ONE_D%RHO_F + H0
               ENDIF
            CASE(-1)
               IF (UU(IBAR,J,K)>0._EB) THEN
                  BXF(J,K) = P_EXTERNAL/WC%ONE_D%RHO_F + KRES(IBAR,J,K)
               ELSE
                  BXF(J,K) = P_EXTERNAL/WC%ONE_D%RHO_F + H0
               ENDIF
            CASE( 2)
               IF (VV(I,0,K)<0._EB) THEN
                  BYS(I,K) = P_EXTERNAL/WC%ONE_D%RHO_F + KRES(I,1,K)
               ELSE
                  BYS(I,K) = P_EXTERNAL/WC%ONE_D%RHO_F + H0
               ENDIF
            CASE(-2)
               IF (VV(I,JBAR,K)>0._EB) THEN
                  BYF(I,K) = P_EXTERNAL/WC%ONE_D%RHO_F + KRES(I,JBAR,K)
               ELSE
                  BYF(I,K) = P_EXTERNAL/WC%ONE_D%RHO_F + H0
               ENDIF
            CASE( 3)
               IF (WW(I,J,0)<0._EB) THEN
                  BZS(I,J) = P_EXTERNAL/WC%ONE_D%RHO_F + KRES(I,J,1)
               ELSE
                  BZS(I,J) = P_EXTERNAL/WC%ONE_D%RHO_F + H0
               ENDIF
            CASE(-3)
               IF (WW(I,J,KBAR)>0._EB) THEN
                  BZF(I,J) = P_EXTERNAL/WC%ONE_D%RHO_F + KRES(I,J,KBAR)
               ELSE
                  BZF(I,J) = P_EXTERNAL/WC%ONE_D%RHO_F + H0
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
         !$OMP PARALLEL DO SIMD PRIVATE(TRM1, TRM2, TRM3, TRM4) SCHEDULE(STATIC)
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
         !$OMP END PARALLEL DO SIMD

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

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
END SUBROUTINE PRESSURE_SOLVER_COMPUTE_RHS


SUBROUTINE PRESSURE_SOLVER_FFT(NM)

USE POIS, ONLY: H3CZSS,H2CZSS,H2CYSS,H3CSSS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE GLOBAL_CONSTANTS

INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
INTEGER :: I,J,K
REAL(EB) :: TNOW

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   HP => H
ELSE
   HP => HS
ENDIF

! Call the Poisson solver

SELECT CASE(IPS)
   CASE(:1)
      IF (.NOT.TWO_D) CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
      IF (TWO_D .AND. .NOT.CYLINDRICAL) CALL H2CZSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
      IF (TWO_D .AND.      CYLINDRICAL) CALL H2CYSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE1,WORK)
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
      IF (LBC==5 .OR. LBC==6)             HP(0,J,K)    = HP(1,J,K)
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

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
END SUBROUTINE PRESSURE_SOLVER_FFT


SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS(NM)

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE GLOBAL_CONSTANTS

INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,RHOP,P,RESIDUAL
INTEGER :: I,J,K
REAL(EB) :: LHSS,RHSS,TNOW

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   HP => H
   RHOP => RHO
ELSE
   HP => HS
   RHOP => RHOS
ENDIF

! Optional check of the accuracy of the separable pressure solution, del^2 H = -del dot F - dD/dt

IF (CHECK_POISSON) THEN
   RESIDUAL => WORK8(1:IBAR,1:JBAR,1:KBAR)
   !$OMP PARALLEL DO PRIVATE(I,J,K,RHSS,LHSS) SCHEDULE(STATIC)
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
            RESIDUAL(I,J,K) = ABS(RHSS-LHSS)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   POIS_ERR = MAXVAL(RESIDUAL)
ENDIF

! Mandatory check of how well the computed pressure satisfies the inseparable Poisson equation:
! LHSS = del dot (1/rho) del p + del K = -del dot F - dD/dt = RHSS

IF (ITERATE_BAROCLINIC_TERM) THEN
   P => WORK7
   P = RHOP*(HP-KRES)
   RESIDUAL => WORK8(1:IBAR,1:JBAR,1:KBAR)
   !$OMP PARALLEL PRIVATE(I,J,K,RHSS,LHSS,NM)
   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHSS = ( R(I-1)*(FVX(I-1,J,K)-FVX_B(I-1,J,K)) - R(I)*(FVX(I,J,K)-FVX_B(I,J,K)) )*RDX(I)*RRN(I) &
                 + (        (FVY(I,J-1,K)-FVY_B(I,J-1,K)) -      (FVY(I,J,K)-FVY_B(I,J,K)) )*RDY(J)        &
                 + (        (FVZ(I,J,K-1)-FVZ_B(I,J,K-1)) -      (FVZ(I,J,K)-FVZ_B(I,J,K)) )*RDZ(K)        &
                 - DDDT(I,J,K)
            LHSS = ((P(I+1,J,K)-P(I,J,K))*RDXN(I)*R(I)    *2._EB/(RHOP(I+1,J,K)+RHOP(I,J,K)) - &
                    (P(I,J,K)-P(I-1,J,K))*RDXN(I-1)*R(I-1)*2._EB/(RHOP(I-1,J,K)+RHOP(I,J,K)))*RDX(I)*RRN(I) &
                 + ((P(I,J+1,K)-P(I,J,K))*RDYN(J)         *2._EB/(RHOP(I,J+1,K)+RHOP(I,J,K)) - &
                    (P(I,J,K)-P(I,J-1,K))*RDYN(J-1)       *2._EB/(RHOP(I,J-1,K)+RHOP(I,J,K)))*RDY(J)        &
                 + ((P(I,J,K+1)-P(I,J,K))*RDZN(K)         *2._EB/(RHOP(I,J,K+1)+RHOP(I,J,K)) - &
                    (P(I,J,K)-P(I,J,K-1))*RDZN(K-1)       *2._EB/(RHOP(I,J,K-1)+RHOP(I,J,K)))*RDZ(K)        &
                 + ((KRES(I+1,J,K)-KRES(I,J,K))*RDXN(I)*R(I) - (KRES(I,J,K)-KRES(I-1,J,K))*RDXN(I-1)*R(I-1) )*RDX(I)*RRN(I) &
                 + ((KRES(I,J+1,K)-KRES(I,J,K))*RDYN(J)      - (KRES(I,J,K)-KRES(I,J-1,K))*RDYN(J-1)        )*RDY(J)        &
                 + ((KRES(I,J,K+1)-KRES(I,J,K))*RDZN(K)      - (KRES(I,J,K)-KRES(I,J,K-1))*RDZN(K-1)        )*RDZ(K)
            RESIDUAL(I,J,K) = ABS(RHSS-LHSS)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
   !$OMP END PARALLEL
   PRESSURE_ERROR_MAX(NM) = MAXVAL(RESIDUAL)
   PRESSURE_ERROR_MAX_LOC(:,NM) = MAXLOC(RESIDUAL)
ENDIF

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
END SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS


SUBROUTINE COMPUTE_VELOCITY_ERROR(DT,NM)

! Check the maximum velocity error at a solid boundary

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE GLOBAL_CONSTANTS, ONLY: PREDICTOR,VELOCITY_ERROR_MAX,SOLID_BOUNDARY,INTERPOLATED_BOUNDARY,VELOCITY_ERROR_MAX_LOC,T_USED,&
                            EXTERNAL_BOUNDARY_CORRECTION,PRES_ON_WHOLE_DOMAIN,PRES_METHOD,FREEZE_VELOCITY,SOLID_PHASE_ONLY

REAL(EB), INTENT(IN) :: DT
INTEGER, INTENT(IN) :: NM
INTEGER :: IW,IOR,II,JJ,KK,IIO,JJO,KKO,N_INT_CELLS,IIO1,IIO2,JJO1,JJO2,KKO1,KKO2
REAL(EB) :: TNOW,UN_NEW,UN_NEW_OTHER,VELOCITY_ERROR,DUDT,DVDT,DWDT,ITERATIVE_FACTOR,DHFCT
TYPE(OMESH_TYPE), POINTER :: OM
TYPE(MESH_TYPE), POINTER :: M2
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
LOGICAL :: GLMAT_ON_WHOLE_DOMAIN

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   ITERATIVE_FACTOR = 0.25_EB
ELSE
   ITERATIVE_FACTOR = 0.50_EB
ENDIF

VELOCITY_ERROR_MAX(NM) = 0._EB
WALL_WORK1 = 0._EB

! Solve Laplace equation for pressure correction, H_PRIME, and add to H or HS.

IF (EXTERNAL_BOUNDARY_CORRECTION) CALL LAPLACE_EXTERNAL_VELOCITY_CORRECTION(DT,NM)

! Logical to define not to apply pressure gradient on external mesh boundaries for GLMAT.
GLMAT_ON_WHOLE_DOMAIN = (PRES_METHOD=='GLMAT') .AND. PRES_ON_WHOLE_DOMAIN

! Loop over wall cells and check velocity error.

CHECK_WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

   WC=>WALL(IW)

   IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY        .AND. &
       WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE CHECK_WALL_LOOP

   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
      EWC=>EXTERNAL_WALL(IW)
      IF (EWC%AREA_RATIO<0.9_EB) CYCLE CHECK_WALL_LOOP
      OM => OMESH(EWC%NOM)
      M2 => MESHES(EWC%NOM)
   ENDIF

   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR

   DHFCT = 1._EB
   IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
      IF ( (.NOT.PRES_ON_WHOLE_DOMAIN) .OR. (GLMAT_ON_WHOLE_DOMAIN .AND.  IW<=N_EXTERNAL_WALL_CELLS) ) DHFCT = 0._EB
   ENDIF

   ! Update normal component of velocity at the mesh boundary

   IF (PREDICTOR) THEN
      SELECT CASE(IOR)
         CASE( 1)
            UN_NEW = U(II,JJ,KK)   - DT*(FVX(II,JJ,KK)   + RDXN(II)  *(H(II+1,JJ,KK)-H(II,JJ,KK))*DHFCT)
         CASE(-1)
            UN_NEW = U(II-1,JJ,KK) - DT*(FVX(II-1,JJ,KK) + RDXN(II-1)*(H(II,JJ,KK)-H(II-1,JJ,KK))*DHFCT)
         CASE( 2)
            UN_NEW = V(II,JJ,KK)   - DT*(FVY(II,JJ,KK)   + RDYN(JJ)  *(H(II,JJ+1,KK)-H(II,JJ,KK))*DHFCT)
         CASE(-2)
            UN_NEW = V(II,JJ-1,KK) - DT*(FVY(II,JJ-1,KK) + RDYN(JJ-1)*(H(II,JJ,KK)-H(II,JJ-1,KK))*DHFCT)
         CASE( 3)
            UN_NEW = W(II,JJ,KK)   - DT*(FVZ(II,JJ,KK)   + RDZN(KK)  *(H(II,JJ,KK+1)-H(II,JJ,KK))*DHFCT)
         CASE(-3)
            UN_NEW = W(II,JJ,KK-1) - DT*(FVZ(II,JJ,KK-1) + RDZN(KK-1)*(H(II,JJ,KK)-H(II,JJ,KK-1))*DHFCT)
      END SELECT
   ELSE
      SELECT CASE(IOR)
         CASE( 1)
            UN_NEW = 0.5_EB*(U(II,JJ,KK)+US(II,JJ,KK)     - DT*(FVX(II,JJ,KK)   + RDXN(II)  *(HS(II+1,JJ,KK)-HS(II,JJ,KK))*DHFCT))
         CASE(-1)
            UN_NEW = 0.5_EB*(U(II-1,JJ,KK)+US(II-1,JJ,KK) - DT*(FVX(II-1,JJ,KK) + RDXN(II-1)*(HS(II,JJ,KK)-HS(II-1,JJ,KK))*DHFCT))
         CASE( 2)
            UN_NEW = 0.5_EB*(V(II,JJ,KK)+VS(II,JJ,KK)     - DT*(FVY(II,JJ,KK)   + RDYN(JJ)  *(HS(II,JJ+1,KK)-HS(II,JJ,KK))*DHFCT))
         CASE(-2)
            UN_NEW = 0.5_EB*(V(II,JJ-1,KK)+VS(II,JJ-1,KK) - DT*(FVY(II,JJ-1,KK) + RDYN(JJ-1)*(HS(II,JJ,KK)-HS(II,JJ-1,KK))*DHFCT))
         CASE( 3)
            UN_NEW = 0.5_EB*(W(II,JJ,KK)+WS(II,JJ,KK)     - DT*(FVZ(II,JJ,KK)   + RDZN(KK)  *(HS(II,JJ,KK+1)-HS(II,JJ,KK))*DHFCT))
         CASE(-3)
            UN_NEW = 0.5_EB*(W(II,JJ,KK-1)+WS(II,JJ,KK-1) - DT*(FVZ(II,JJ,KK-1) + RDZN(KK-1)*(HS(II,JJ,KK)-HS(II,JJ,KK-1))*DHFCT))
      END SELECT
   ENDIF

   ! At interpolated boundaries, compare updated normal component of velocity with that of the other mesh

   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN

      UN_NEW_OTHER = 0._EB

      EWC=>EXTERNAL_WALL(IW)
      IIO1 = EWC%IIO_MIN
      JJO1 = EWC%JJO_MIN
      KKO1 = EWC%KKO_MIN
      IIO2 = EWC%IIO_MAX
      JJO2 = EWC%JJO_MAX
      KKO2 = EWC%KKO_MAX

      PREDICTOR_IF: IF (PREDICTOR) THEN
         IOR_SELECT_1: SELECT CASE(IOR)
            CASE( 1)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DUDT = -OM%FVX(IIO,JJO,KKO)   - M2%RDXN(IIO)  *(OM%H(IIO+1,JJO,KKO)-OM%H(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%U(IIO,JJO,KKO)   + DT*DUDT
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-1)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DUDT = -OM%FVX(IIO-1,JJO,KKO) - M2%RDXN(IIO-1)*(OM%H(IIO,JJO,KKO)-OM%H(IIO-1,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%U(IIO-1,JJO,KKO) + DT*DUDT
                     ENDDO
                  ENDDO
               ENDDO
            CASE( 2)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DVDT = -OM%FVY(IIO,JJO,KKO)   - M2%RDYN(JJO)  *(OM%H(IIO,JJO+1,KKO)-OM%H(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%V(IIO,JJO,KKO)   + DT*DVDT
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-2)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DVDT = -OM%FVY(IIO,JJO-1,KKO) - M2%RDYN(JJO-1)*(OM%H(IIO,JJO,KKO)-OM%H(IIO,JJO-1,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%V(IIO,JJO-1,KKO) + DT*DVDT
                     ENDDO
                  ENDDO
               ENDDO
            CASE( 3)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DWDT = -OM%FVZ(IIO,JJO,KKO)   - M2%RDZN(KKO)  *(OM%H(IIO,JJO,KKO+1)-OM%H(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%W(IIO,JJO,KKO)   + DT*DWDT
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-3)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DWDT = -OM%FVZ(IIO,JJO,KKO-1) - M2%RDZN(KKO-1)*(OM%H(IIO,JJO,KKO)-OM%H(IIO,JJO,KKO-1))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%W(IIO,JJO,KKO-1) + DT*DWDT
                     ENDDO
                  ENDDO
               ENDDO
         END SELECT IOR_SELECT_1
      ELSE PREDICTOR_IF
         IOR_SELECT_2: SELECT CASE(IOR)
            CASE( 1)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DUDT = -OM%FVX(IIO,JJO,KKO)   - M2%RDXN(IIO)  *(OM%HS(IIO+1,JJO,KKO)-OM%HS(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%U(IIO,JJO,KKO)+OM%US(IIO,JJO,KKO)     + DT*DUDT)
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-1)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DUDT = -OM%FVX(IIO-1,JJO,KKO) - M2%RDXN(IIO-1)*(OM%HS(IIO,JJO,KKO)-OM%HS(IIO-1,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%U(IIO-1,JJO,KKO)+OM%US(IIO-1,JJO,KKO) + DT*DUDT)
                     ENDDO
                  ENDDO
               ENDDO
            CASE( 2)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DVDT = -OM%FVY(IIO,JJO,KKO)   - M2%RDYN(JJO)  *(OM%HS(IIO,JJO+1,KKO)-OM%HS(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%V(IIO,JJO,KKO)+OM%VS(IIO,JJO,KKO)     + DT*DVDT)
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-2)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DVDT = -OM%FVY(IIO,JJO-1,KKO) - M2%RDYN(JJO-1)*(OM%HS(IIO,JJO,KKO)-OM%HS(IIO,JJO-1,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%V(IIO,JJO-1,KKO)+OM%VS(IIO,JJO-1,KKO) + DT*DVDT)
                     ENDDO
                  ENDDO
               ENDDO
            CASE( 3)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DWDT = -OM%FVZ(IIO,JJO,KKO)   - M2%RDZN(KKO)  *(OM%HS(IIO,JJO,KKO+1)-OM%HS(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%W(IIO,JJO,KKO)+OM%WS(IIO,JJO,KKO)     + DT*DWDT)
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-3)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DWDT = -OM%FVZ(IIO,JJO,KKO-1) - M2%RDZN(KKO-1)*(OM%HS(IIO,JJO,KKO)-OM%HS(IIO,JJO,KKO-1))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%W(IIO,JJO,KKO-1)+OM%WS(IIO,JJO,KKO-1) + DT*DWDT)
                     ENDDO
                  ENDDO
               ENDDO
         END SELECT IOR_SELECT_2
      ENDIF PREDICTOR_IF

      N_INT_CELLS  = (EWC%IIO_MAX-EWC%IIO_MIN+1) * (EWC%JJO_MAX-EWC%JJO_MIN+1) * (EWC%KKO_MAX-EWC%KKO_MIN+1)
      UN_NEW_OTHER = UN_NEW_OTHER/REAL(N_INT_CELLS,EB)

   ENDIF

   ! At solid boundaries, compare updated normal velocity with specified normal velocity

   IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
      IF (PREDICTOR) THEN
         UN_NEW_OTHER = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%U_NORMAL_S
      ELSE
         UN_NEW_OTHER = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%U_NORMAL
      ENDIF
   ENDIF

   ! Compute velocity difference

   VELOCITY_ERROR = UN_NEW - UN_NEW_OTHER
   WC%VEL_ERR_NEW = VELOCITY_ERROR
   WALL_WORK1(IW) = -SIGN(1._EB,REAL(IOR,EB))*ITERATIVE_FACTOR*VELOCITY_ERROR/(WC%ONE_D%RDN*DT)

   ! If the grid cells in the current mesh are smaller than those of the other mesh, do not include in error tolerance

   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
      IF (OM%NIC_R>OM%NIC_S) CYCLE CHECK_WALL_LOOP
   ENDIF

   ! Save maximum velocity error

   IF (ABS(VELOCITY_ERROR)>VELOCITY_ERROR_MAX(NM)) THEN
      VELOCITY_ERROR_MAX_LOC(1,NM) = II
      VELOCITY_ERROR_MAX_LOC(2,NM) = JJ
      VELOCITY_ERROR_MAX_LOC(3,NM) = KK
      VELOCITY_ERROR_MAX(NM)       = ABS(VELOCITY_ERROR)
   ENDIF

ENDDO CHECK_WALL_LOOP

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
END SUBROUTINE COMPUTE_VELOCITY_ERROR


SUBROUTINE LAPLACE_EXTERNAL_VELOCITY_CORRECTION(DT,NM)
USE POIS, ONLY: H3CZSS,H2CZSS
USE GLOBAL_CONSTANTS, ONLY: PREDICTOR,TWO_D,SOLID_BOUNDARY,OPEN_BOUNDARY,INTERPOLATED_BOUNDARY,CHECK_POISSON,DIRICHLET,NEUMANN

REAL(EB), INTENT(IN) :: DT
INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K,IOR,II,JJ,KK,IW
REAL(EB) :: AXS_SUM,AXF_SUM,AYS_SUM,AYF_SUM,AZS_SUM,AZF_SUM,BXS_SUM,BXF_SUM,BYS_SUM,BYF_SUM,BZS_SUM,BZF_SUM,&
            UN_NEW,UN_NEW_OTHER,LHSS,DT_LOC
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
TYPE(WALL_TYPE), POINTER :: WC

CALL POINT_TO_MESH(NM)

! Build Neumann boundary conditions from external wall velocity errors

AXS_SUM=0._EB
AXF_SUM=0._EB
AYS_SUM=0._EB
AYF_SUM=0._EB
AZS_SUM=0._EB
AZF_SUM=0._EB

BXS_SUM=0._EB
BXF_SUM=0._EB
BYS_SUM=0._EB
BYF_SUM=0._EB
BZS_SUM=0._EB
BZF_SUM=0._EB

BXS=0._EB
BXF=0._EB
BYS=0._EB
BYF=0._EB
BZS=0._EB
BZF=0._EB

IF (PREDICTOR) THEN
   DT_LOC=DT
ELSE
   DT_LOC=0.5_EB*DT
ENDIF

EXTERNAL_SOLID_BOUNDARY_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

   WC=>WALL(IW)

   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR

   ! Sum OPEN flow area for compatibility correction used later

   IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
      SELECT CASE(IOR)
         CASE( 1)
            AXS_SUM  = AXS_SUM + DY(JJ)*DZ(KK)
         CASE(-1)
            AXF_SUM  = AXF_SUM + DY(JJ)*DZ(KK)
         CASE( 2)
            AYS_SUM  = AYS_SUM + DX(II)*DZ(KK)
         CASE(-2)
            AYF_SUM  = AYF_SUM + DX(II)*DZ(KK)
         CASE( 3)
            AZS_SUM  = AZS_SUM + DX(II)*DY(JJ)
         CASE(-3)
            AZF_SUM  = AZF_SUM + DX(II)*DY(JJ)
      END SELECT
      CYCLE EXTERNAL_SOLID_BOUNDARY_LOOP
   ENDIF

   IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE EXTERNAL_SOLID_BOUNDARY_LOOP

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

   ! At solid boundaries, compare updated normal velocity with specified normal velocity

   IF (PREDICTOR) THEN
      UN_NEW_OTHER = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%U_NORMAL_S
   ELSE
      UN_NEW_OTHER = -SIGN(1._EB,REAL(IOR,EB))*WC%ONE_D%U_NORMAL
   ENDIF

   LAPLACE_BC_INDEX_SELECT: SELECT CASE(WC%LAPLACE_BC_INDEX)

      CASE(NEUMANN) LAPLACE_BC_INDEX_SELECT

         ! Compute velocity error and apply as Neumann BC to Laplace

         SELECT CASE(IOR)
            CASE( 1)
               BXS(JJ,KK) = (UN_NEW - UN_NEW_OTHER)/DT_LOC
               BXS_SUM = BXS_SUM + BXS(JJ,KK)*DY(JJ)*DZ(KK)
            CASE(-1)
               BXF(JJ,KK) = (UN_NEW - UN_NEW_OTHER)/DT_LOC
               BXF_SUM = BXF_SUM + BXF(JJ,KK)*DY(JJ)*DZ(KK)
            CASE( 2)
               BYS(II,KK) = (UN_NEW - UN_NEW_OTHER)/DT_LOC
               BYS_SUM = BYS_SUM + BYS(II,KK)*DX(II)*DZ(KK)
            CASE(-2)
               BYF(II,KK) = (UN_NEW - UN_NEW_OTHER)/DT_LOC
               BYF_SUM = BYF_SUM + BYF(II,KK)*DX(II)*DZ(KK)
            CASE( 3)
               BZS(II,JJ) = (UN_NEW - UN_NEW_OTHER)/DT_LOC
               BZS_SUM = BZS_SUM + BZS(II,JJ)*DX(II)*DY(JJ)
            CASE(-3)
               BZF(II,JJ) = (UN_NEW - UN_NEW_OTHER)/DT_LOC
               BZF_SUM = BZF_SUM + BZF(II,JJ)*DX(II)*DY(JJ)
         END SELECT

      CASE(DIRICHLET) LAPLACE_BC_INDEX_SELECT

         ! Trick Dirichlet bc into approximating Neumann bc

         SELECT CASE(IOR)
            CASE( 1)
               BXS(JJ,KK) = 0._EB ! (UN_NEW_OTHER - UN_NEW)/(2._EB*DT_LOC) * DXN(II)
            CASE(-1)
               BXF(JJ,KK) = 0._EB ! (UN_NEW - UN_NEW_OTHER)/(2._EB*DT_LOC) * DXN(II-1)
            CASE( 2)
               BYS(II,KK) = 0._EB ! (UN_NEW_OTHER - UN_NEW)/(2._EB*DT_LOC) * DYN(JJ)
            CASE(-2)
               BYF(II,KK) = 0._EB ! (UN_NEW - UN_NEW_OTHER)/(2._EB*DT_LOC) * DYN(JJ-1)
            CASE( 3)
               BZS(II,JJ) = 0._EB ! (UN_NEW_OTHER - UN_NEW)/(2._EB*DT_LOC) * DZN(KK)
            CASE(-3)
               BZF(II,JJ) = 0._EB ! (UN_NEW - UN_NEW_OTHER)/(2._EB*DT_LOC) * DZN(KK-1)
         END SELECT

   END SELECT LAPLACE_BC_INDEX_SELECT

ENDDO EXTERNAL_SOLID_BOUNDARY_LOOP

! Now loop over OPEN boundary and set Neumann BC

EXTERNAL_OPEN_BOUNDARY_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

   WC=>WALL(IW)

   IF (WC%BOUNDARY_TYPE/=OPEN_BOUNDARY .OR. WC%LAPLACE_BC_INDEX/=NEUMANN) CYCLE EXTERNAL_OPEN_BOUNDARY_LOOP

   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR

   SELECT CASE(IOR)
      CASE( 1)
         BXS(JJ,KK) = -BXS_SUM/AXS_SUM ! there should be no div by zero problem, if ASX_SUM=0 there is a problem elsewhere
      CASE(-1)
         BXF(JJ,KK) = -BXF_SUM/AXF_SUM
      CASE( 2)
         BYS(II,KK) = -BYS_SUM/AYS_SUM
      CASE(-2)
         BYF(II,KK) = -BYF_SUM/AYF_SUM
      CASE( 3)
         BZS(II,JJ) = -BZS_SUM/AZS_SUM
      CASE(-3)
         BZF(II,JJ) = -BZF_SUM/AZF_SUM
   END SELECT

ENDDO EXTERNAL_OPEN_BOUNDARY_LOOP

! Initialize solution and Laplace RHS (zero)

HP=>H_PRIME; HP=0._EB
PRHS=0._EB

! Solve Laplace equation

IF (.NOT.TWO_D) CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,LAPLACE_PTB,SAVE2,WORK,HX)
HP(1:IBAR,1:JBAR,1:KBAR) = PRHS

! Apply boundary conditions to H'

DO K=1,KBAR
   DO J=1,JBAR
      IF (LBC2==3 .OR. LBC2==4) HP(0,J,K)    = HP(1,J,K)    - DXI*BXS(J,K)
      IF (LBC2==3 .OR. LBC2==2) HP(IBP1,J,K) = HP(IBAR,J,K) + DXI*BXF(J,K)
      IF (LBC2==1 .OR. LBC2==2) HP(0,J,K)    =-HP(1,J,K)    + 2._EB*BXS(J,K)
      IF (LBC2==1 .OR. LBC2==4) HP(IBP1,J,K) =-HP(IBAR,J,K) + 2._EB*BXF(J,K)
   ENDDO
ENDDO

DO K=1,KBAR
   DO I=1,IBAR
      IF (MBC2==3 .OR. MBC2==4) HP(I,0,K)    = HP(I,1,K)    - DETA*BYS(I,K)
      IF (MBC2==3 .OR. MBC2==2) HP(I,JBP1,K) = HP(I,JBAR,K) + DETA*BYF(I,K)
      IF (MBC2==1 .OR. MBC2==2) HP(I,0,K)    =-HP(I,1,K)    + 2._EB*BYS(I,K)
      IF (MBC2==1 .OR. MBC2==4) HP(I,JBP1,K) =-HP(I,JBAR,K) + 2._EB*BYF(I,K)
   ENDDO
ENDDO

DO J=1,JBAR
   DO I=1,IBAR
      IF (NBC2==3 .OR. NBC2==4) HP(I,J,0)    = HP(I,J,1)    - DZETA*BZS(I,J)
      IF (NBC2==3 .OR. NBC2==2) HP(I,J,KBP1) = HP(I,J,KBAR) + DZETA*BZF(I,J)
      IF (NBC2==1 .OR. NBC2==2) HP(I,J,0)    =-HP(I,J,1)    + 2._EB*BZS(I,J)
      IF (NBC2==1 .OR. NBC2==4) HP(I,J,KBP1) =-HP(I,J,KBAR) + 2._EB*BZF(I,J)
   ENDDO
ENDDO

! ************************* Check the Solution *************************

IF (CHECK_POISSON) THEN

   LAPLACE_ERR = 0._EB
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            LHSS = ( (HP(I+1,J,K)-HP(I,J,K))*RDXN(I) - (HP(I,J,K)-HP(I-1,J,K))*RDXN(I-1) )*RDX(I) &
                 + ( (HP(I,J+1,K)-HP(I,J,K))*RDYN(J) - (HP(I,J,K)-HP(I,J-1,K))*RDYN(J-1) )*RDY(J) &
                 + ( (HP(I,J,K+1)-HP(I,J,K))*RDZN(K) - (HP(I,J,K)-HP(I,J,K-1))*RDZN(K-1) )*RDZ(K)
            LAPLACE_ERR = MAX(ABS(LHSS),LAPLACE_ERR)
         ENDDO
      ENDDO
   ENDDO

   ! IF (PREDICTOR) THEN
   !    WRITE(0,*) 'PREDICTOR, MESH:',NM
   ! ELSE
   !    WRITE(0,*) 'CORRECTOR, MESH:',NM
   ! ENDIF
   ! WRITE(0,*) 'LAPLACE PTB/ERROR:',LAPLACE_PTB,LAPLACE_ERR

ENDIF

! Enforce continuity of HP(=0) at interpolated boundaries

INTERPOLATED_BOUNDARY_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

   WC=>WALL(IW)

   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE INTERPOLATED_BOUNDARY_LOOP

   II  = WC%ONE_D%II
   JJ  = WC%ONE_D%JJ
   KK  = WC%ONE_D%KK
   IOR = WC%ONE_D%IOR

   SELECT CASE(IOR)
      CASE( 1)
         HP(0,JJ,KK)    = -HP(1,JJ,KK)
      CASE(-1)
         HP(IBP1,JJ,KK) = -HP(IBAR,JJ,KK)
      CASE( 2)
         HP(II,0,KK)    = -HP(II,1,KK)
      CASE(-2)
         HP(II,JBP1,KK) = -HP(II,JBAR,KK)
      CASE( 3)
         HP(II,JJ,0)    = -HP(II,JJ,1)
      CASE(-3)
         HP(II,JJ,KBP1) = -HP(II,JJ,KBAR)
   END SELECT

ENDDO INTERPOLATED_BOUNDARY_LOOP

! Add correction to H or HS field

IF (PREDICTOR) THEN
   H  = H  + HP
ELSE
   HS = HS + HP
ENDIF

END SUBROUTINE LAPLACE_EXTERNAL_VELOCITY_CORRECTION


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


END MODULE PRES


! ---------------------------------- GLOBALMATRIX_SOLVER --------------------------------------------

MODULE GLOBMAT_SOLVER

! Module that contains global matrix vector builds for Poisson equation on gas-cells only when PRES_ON_WHOLE_DOMAIN=.FALSE.
! Builds Matrices and RHS entries per MPI process in parallel.
! Calls MKL sparse cluster solver or Pardiso for the time being.

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MESH_POINTERS

USE COMPLEX_GEOMETRY, ONLY : IBM_CGSC,IBM_FGSC, IBM_UNKH, IBM_NCVARS, GET_H_CUTFACES,        &
                             GET_BOUNDFACE_GEOM_INFO_H, ADD_INPLACE_NNZ_H_WHLDOM,   &
                             NUNKH_LOC, NUNKH_TOT, UNKH_IND, NUNKH_LOCAL, NUNKH_TOTAL, NM_START, &
                             NNZ_ROW_H, TOT_NNZ_H, NNZ_D_MAT_H, D_MAT_H, JD_MAT_H, IA_H,       &
                             JA_H, A_H, H_MATRIX_INDEFINITE, F_H, X_H, PT_H, IPARM, COPY_CC_UNKH_TO_HS, &
                             COPY_CC_HS_TO_UNKH

#ifdef WITH_MKL
USE MKL_CLUSTER_SPARSE_SOLVER
#endif /* WITH_MKL */

IMPLICIT NONE

! These definitions are the same as geom.f90:
INTEGER,  PARAMETER :: NGUARD= 2 ! Two layers of guard-cells.
INTEGER,  PARAMETER :: FCELL = 1 ! Right face index.

! Media definition parameters, same numerical values as in geom.f90:
INTEGER,  PARAMETER :: IS_GASPHASE  = -1
INTEGER,  PARAMETER :: IS_CUTCFE    =  0
INTEGER,  PARAMETER :: IS_SOLID     =  1
INTEGER,  PARAMETER :: IS_UNDEFINED =-11

! Cartesian Cell centered variables, case CC_IBM=.FALSE.:
INTEGER,  PARAMETER :: IS_CGSC   = 1 ! Face media type: IS_GASPHASE, IS_SOLID or IS_CUTCFE.
INTEGER,  PARAMETER :: IS_UNKH   = 2 ! H unknown number.
INTEGER,  PARAMETER :: IS_NCVARS = 2 ! Number of face variables in MESHES(NM)%CCVAR.

INTEGER, SAVE :: ILO_CELL,IHI_CELL,JLO_CELL,JHI_CELL,KLO_CELL,KHI_CELL
INTEGER, SAVE :: ILO_FACE,IHI_FACE,JLO_FACE,JHI_FACE,KLO_FACE,KHI_FACE
INTEGER, SAVE :: NXB, NYB, NZB

! Cartesian Cell centered variables, actual case initialized as CC_IBM=.FALSE.:
INTEGER :: CGSC=IS_CGSC, UNKH=IS_UNKH, NCVARS=IS_NCVARS

! Pardiso or Sparse cluster solver message level:
INTEGER, SAVE :: MSGLVL = 0  ! 0 no messages, 1 print statistical information

! Factor to drop DY in cylindrical axisymmetric coordinates.
REAL(EB), SAVE :: CYL_FCT

!#define SINGLE_PRECISION_PSN_SOLVE
#ifdef SINGLE_PRECISION_PSN_SOLVE
REAL(FB), ALLOCATABLE, DIMENSION(:) :: F_H_FB, X_H_FB, A_H_FB
#endif

! Timing variable:
REAL(EB):: TNOW

PRIVATE

PUBLIC GLMAT_SOLVER_SETUP_H,GLMAT_SOLVER_H,COPY_H_OMESH_TO_MESH,FINISH_GLMAT_SOLVER_H,PRESSURE_SOLVER_CHECK_RESIDUALS_U

CONTAINS

! --------------------------- GLMAT_SOLVER_H -------------------------------------

SUBROUTINE GLMAT_SOLVER_H

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE COMPLEX_GEOMETRY, ONLY : GET_CUTCELL_FH,GET_CUTCELL_HP,GET_CC_IROW
USE MPI

! Local Variables:
INTEGER :: MAXFCT, MNUM, MTYPE, PHASE, NRHS, ERROR
#ifdef WITH_MKL
INTEGER :: PERM(1)
#endif
INTEGER :: NM, IW, IIG, JJG, KKG, IOR, IROW, I, J, K, ICC, I_ZONE
TYPE (WALL_TYPE), POINTER :: WC=>NULL()
REAL(EB) :: IDX, AF, VAL
REAL(EB), POINTER, DIMENSION(:,:,:)   :: HP
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: SUM_FH, SUM_XH ! For use in cases with all compartments without OPEN boundaries.
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: MEAN_FH, MEAN_XH
INTEGER :: IERR

! CHARACTER(30) :: FILE_NAME
! INTEGER :: ICC, IERR

IF (FREEZE_VELOCITY) RETURN ! Fixed velocity soln. i.e. PERIODIC_TEST=102 => FREEZE_VELOCITY=.TRUE.
IF (SOLID_PHASE_ONLY) RETURN
TNOW=CURRENT_TIME()

! Solve:
NRHS   =  1
MAXFCT =  1
MNUM   =  1
ERROR  =  0 ! initialize error flag

! Define rhs F_H, here we use Source and BCs populated on PRESSURE_SOLVER:
F_H(1:NUNKH_LOCAL) = 0._EB
X_H(1:NUNKH_LOCAL) = 0._EB

! Main Mesh Loop:
MESH_LOOP_1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION_ONLY(NM) .OR. EVACUATION_SKIP(NM)) CYCLE MESH_LOOP_1
   CALL POINT_TO_MESH(NM)

   ! First Source on Cartesian cells with IBM_UNKH > 0:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,UNKH) <= 0) CYCLE ! Gasphase Cartesian cells.
            ! Row number:
            IROW = CCVAR(I,J,K,UNKH) - UNKH_IND(NM_START) ! Local numeration.
            ! Add to F_H: If CYL_FCT=0. -> Cartesian coordinates volume (RC(I)=1.). If CYL_FCT=1. -> Cylindrical coords volume.
            F_H(IROW) = F_H(IROW) + PRHS(I,J,K) * ((1._EB-CYL_FCT)*DY(J) + CYL_FCT*RC(I))*DX(I)*DZ(K)
         ENDDO
      ENDDO
   ENDDO

   IF (CC_IBM) CALL GET_CUTCELL_FH(NM) ! Note: CYL_FCT not used for cut-cells.

   ! Then External BCs:
   WALL_CELL_LOOP_1: DO IW=1,N_EXTERNAL_WALL_CELLS

      WC => WALL(IW)

      IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE

      ! NEUMANN boundaries:
      IF_NEUMANN: IF (WC%PRESSURE_BC_INDEX==NEUMANN) THEN

         ! Gasphase cell indexes:
         IIG   = WC%ONE_D%IIG
         JJG   = WC%ONE_D%JJG
         KKG   = WC%ONE_D%KKG
         IOR   = WC%ONE_D%IOR

         ! Define cell size, normal to WC:
         SELECT CASE (IOR)
         CASE(-1) ! -IAXIS oriented, high face of IIG cell.
            AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG  )) * DZ(KKG)
            VAL = -BXF(JJG,KKG)*AF
         CASE( 1) ! +IAXIS oriented, low face of IIG cell.
            AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG-1)) * DZ(KKG)
            VAL =  BXS(JJG,KKG)*AF
         CASE(-2) ! -JAXIS oriented, high face of JJG cell.
            AF  =  DX(IIG)*DZ(KKG)
            VAL = -BYF(IIG,KKG)*AF
         CASE( 2) ! +JAXIS oriented, low face of JJG cell.
            AF  =  DX(IIG)*DZ(KKG)
            VAL =  BYS(IIG,KKG)*AF
         CASE(-3) ! -KAXIS oriented, high face of KKG cell.
            AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG)
            VAL = -BZF(IIG,JJG)*AF
         CASE( 3) ! +KAXIS oriented, low face of KKG cell.
            AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG)
            VAL =  BZS(IIG,JJG)*AF
         END SELECT

         ! Row number:
         IROW = CCVAR(IIG,JJG,KKG,UNKH) - UNKH_IND(NM_START) ! Local numeration.
         IF (IROW <= 0 .AND. CC_IBM) THEN
            CALL GET_CC_IROW(IIG,JJG,KKG,IROW)
            IF (IROW <= 0) CYCLE
         ENDIF

         IF(IROW==IS_UNDEFINED) &
            WRITE(LU_ERR,*) 'CELL W IBM_UNDEFINED IN UNKH=',IIG,JJG,KKG,IROW,CCVAR(IIG,JJG,KKG,CGSC)

         ! Add to F_H:
         F_H(IROW) = F_H(IROW) + VAL

      ENDIF IF_NEUMANN

      ! DIRICHLET boundaries:
      IF_DIRICHLET: IF (WC%PRESSURE_BC_INDEX==DIRICHLET) THEN

         IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
             WC%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
             WC%BOUNDARY_TYPE==SOLID_BOUNDARY) CYCLE    ! No need for these, that's the whole point of a
                                                        ! global solve.

         ! Gasphase cell indexes:
         IIG   = WC%ONE_D%IIG
         JJG   = WC%ONE_D%JJG
         KKG   = WC%ONE_D%KKG
         IOR   = WC%ONE_D%IOR

         ! Define cell size, normal to WC:
         SELECT CASE (IOR)
         CASE(-1) ! -IAXIS oriented, high face of IIG cell.
            IDX = 1._EB / DXN(IIG)
            AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG  )) * DZ(KKG)
            VAL = -2._EB*IDX*AF*BXF(JJG,KKG)
         CASE( 1) ! +IAXIS oriented, low face of IIG cell.
            IDX = 1._EB / DXN(IIG-1)
            AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG-1)) * DZ(KKG)
            VAL = -2._EB*IDX*AF*BXS(JJG,KKG)
         CASE(-2) ! -JAXIS oriented, high face of JJG cell.
            IDX = 1._EB / DYN(JJG)
            AF  =  DX(IIG)*DZ(KKG)
            VAL = -2._EB*IDX*AF*BYF(IIG,KKG)
         CASE( 2) ! +JAXIS oriented, low face of JJG cell.
            IDX = 1._EB / DYN(JJG-1)
            AF  =  DX(IIG)*DZ(KKG)
            VAL = -2._EB*IDX*AF*BYS(IIG,KKG)
         CASE(-3) ! -KAXIS oriented, high face of KKG cell.
            IDX = 1._EB / DZN(KKG)
            AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG)
            VAL = -2._EB*IDX*AF*BZF(IIG,JJG)
         CASE( 3) ! +KAXIS oriented, low face of KKG cell.
            IDX = 1._EB / DZN(KKG-1)
            AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG)
            VAL = -2._EB*IDX*AF*BZS(IIG,JJG)
         END SELECT

         ! Row number:
         IROW = CCVAR(IIG,JJG,KKG,UNKH) - UNKH_IND(NM_START) ! Local numeration.
         IF (IROW <= 0 .AND. CC_IBM) CALL GET_CC_IROW(IIG,JJG,KKG,IROW)
         IF (IROW <= 0) CYCLE
         ! Add to F_H:
         F_H(IROW) = F_H(IROW) + VAL

      ENDIF IF_DIRICHLET

   ENDDO WALL_CELL_LOOP_1

   ! Here we include pressure boundary conditions due to OBST and GEOM surfaces:
   ! If non-reacting SOLIDS, no need to do anything.


ENDDO MESH_LOOP_1

IF ( H_MATRIX_INDEFINITE ) THEN
   MTYPE  = -2 ! symmetric indefinite
ELSE ! positive definite
   MTYPE  =  2
ENDIF

IF (H_MATRIX_INDEFINITE) THEN
   ALLOCATE(SUM_FH(1:3,0:N_ZONE),SUM_XH(1:3,0:N_ZONE)); SUM_FH = 0._EB;  SUM_XH = 0._EB
   ALLOCATE(MEAN_FH(0:N_ZONE), MEAN_XH(0:N_ZONE));     MEAN_FH = 0._EB; MEAN_XH = 0._EB
   WHOLE_DOM_IF1 : IF(.NOT.PRES_ON_WHOLE_DOMAIN) THEN
      ! Sum source F_H by Pressure Zone:
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION_ONLY(NM) .OR. EVACUATION_SKIP(NM)) CYCLE
         CALL POINT_TO_MESH(NM)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (CCVAR(I,J,K,UNKH)<=0 .OR. PRESSURE_ZONE(I,J,K)<=0) CYCLE ! Gasphase Cartesian cells.
                  ! Row number:
                  IROW = CCVAR(I,J,K,UNKH) - UNKH_IND(NM_START) ! Local numeration.
                  ! Sum FH:
                  SUM_FH(1,PRESSURE_ZONE(I,J,K)) = SUM_FH(1,PRESSURE_ZONE(I,J,K)) + F_H(IROW)
                  SUM_FH(2,PRESSURE_ZONE(I,J,K)) = SUM_FH(2,PRESSURE_ZONE(I,J,K)) + 1._EB
               ENDDO
            ENDDO
         ENDDO
         ! Add cut-cell region contribution:
         DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
            I = CUT_CELL(ICC)%IJK(IAXIS)
            J = CUT_CELL(ICC)%IJK(JAXIS)
            K = CUT_CELL(ICC)%IJK(KAXIS)
            IF (PRESSURE_ZONE(I,J,K)<=0) CYCLE
            IROW     = MESHES(NM)%CUT_CELL(ICC)%UNKH(1) - UNKH_IND(NM_START) ! Local numeration.
            SUM_FH(1,PRESSURE_ZONE(I,J,K)) = SUM_FH(1,PRESSURE_ZONE(I,J,K)) + F_H(IROW)
            SUM_FH(2,PRESSURE_ZONE(I,J,K)) = SUM_FH(2,PRESSURE_ZONE(I,J,K)) + 1._EB
         ENDDO
      ENDDO

      IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,SUM_FH(1:2,0:N_ZONE),2*(N_ZONE+1),&
                                                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

      ! Compute arithmetic mean by pressure zone:
      DO I_ZONE=0,N_ZONE
         MEAN_FH(I_ZONE) = SUM_FH(1,I_ZONE)/(SUM_FH(2,I_ZONE)+TWO_EPSILON_EB)
      ENDDO
      ! Write out:
      ! IF (MYID==0) THEN
      !    DO I_ZONE=0,N_ZONE
      !    WRITE(LU_ERR,*) PREDICTOR,'INDEFINITE POISSON MATRIX, I_ZONE, MEAN(RHS), SUM(RHS)=',&
      !                    I_ZONE,MEAN_FH(I_ZONE),SUM_FH(1:2,I_ZONE)
      !    ENDDO
      ! ENDIF

      ! Substract Mean:
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION_ONLY(NM) .OR. EVACUATION_SKIP(NM)) CYCLE
         CALL POINT_TO_MESH(NM)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (CCVAR(I,J,K,UNKH)<=0 .OR. PRESSURE_ZONE(I,J,K)<=0) CYCLE ! Gasphase Cartesian cells.
                  ! Row number:
                  IROW = CCVAR(I,J,K,UNKH) - UNKH_IND(NM_START) ! Local numeration.
                  F_H(IROW) = F_H(IROW) - MEAN_FH(PRESSURE_ZONE(I,J,K))
               ENDDO
            ENDDO
         ENDDO
         ! Add cut-cell region contribution:
         DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
            I = CUT_CELL(ICC)%IJK(IAXIS)
            J = CUT_CELL(ICC)%IJK(JAXIS)
            K = CUT_CELL(ICC)%IJK(KAXIS)
            IF (PRESSURE_ZONE(I,J,K)<=0) CYCLE
            IROW     = MESHES(NM)%CUT_CELL(ICC)%UNKH(1) - UNKH_IND(NM_START) ! Local numeration.
            F_H(IROW) = F_H(IROW) - MEAN_FH(PRESSURE_ZONE(I,J,K))
         ENDDO
      ENDDO
   ELSE WHOLE_DOM_IF1
      SUM_FH(1:2,0) = SUM(F_H(1:NUNKH_LOCAL))
      IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(SUM_FH(1,0),SUM_FH(2,0),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
      MEAN_FH(0) = SUM_FH(2,0)/REAL(NUNKH_TOTAL,EB)
      ! IF (MYID==0) WRITE(LU_ERR,*) 'INDEFINITE POISSON MATRIX, MEAN(RHS), SUM(RHS)=',MEAN_FH(0),SUM_FH(2,0)
      ! Substract Mean:
      F_H(:) = F_H(:) - MEAN_FH(0)
   ENDIF WHOLE_DOM_IF1
ENDIF

! WRITE(LU_ERR,*) 'SUM_FH=',SUM(F_H),H_MATRIX_INDEFINITE

!.. Back substitution and iterative refinement
IPARM(8) =  0 ! max numbers of iterative refinement steps
PHASE    = 33 ! only solving
!   CALL PARDISO(PT_H, MAXFCT, MNUM, MTYPE, PHASE, NUNKH_TOTAL, &
!              A_H, IA_H, JA_H, PERM, NRHS, IPARM, MSGLVL, F_H, X_H, ERROR)
#ifdef WITH_MKL
#ifdef SINGLE_PRECISION_PSN_SOLVE
F_H_FB(1:NUNKH_LOCAL) = REAL(F_H(1:NUNKH_LOCAL),FB)
X_H_FB(1:NUNKH_LOCAL) = 0._FB
CALL CLUSTER_SPARSE_SOLVER(PT_H, MAXFCT, MNUM, MTYPE, PHASE, NUNKH_TOTAL, &
             A_H_FB, IA_H, JA_H, PERM, NRHS, IPARM, MSGLVL, F_H_FB, X_H_FB, MPI_COMM_WORLD, ERROR)
X_H(1:NUNKH_LOCAL) = REAL(X_H_FB(1:NUNKH_LOCAL),EB)
#else
CALL CLUSTER_SPARSE_SOLVER(PT_H, MAXFCT, MNUM, MTYPE, PHASE, NUNKH_TOTAL, &
             A_H, IA_H, JA_H, PERM, NRHS, IPARM, MSGLVL, F_H, X_H, MPI_COMM_WORLD, ERROR)
#endif
#endif
IF (ERROR /= 0) &
WRITE(0,*) 'GLMAT_SOLVER_H: The following ERROR was detected: ', ERROR

IF (H_MATRIX_INDEFINITE) THEN
   WHOLE_DOM_IF2 : IF(.NOT.PRES_ON_WHOLE_DOMAIN) THEN
      ! Sum H by Pressure Zone:
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION_ONLY(NM) .OR. EVACUATION_SKIP(NM)) CYCLE
         CALL POINT_TO_MESH(NM)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (CCVAR(I,J,K,UNKH)<=0 .OR. PRESSURE_ZONE(I,J,K)<=0) CYCLE ! Gasphase Cartesian cells.
                  ! Row number:
                  IROW = CCVAR(I,J,K,UNKH) - UNKH_IND(NM_START) ! Local numeration.
                  ! Sum FH:
                  SUM_XH(1,PRESSURE_ZONE(I,J,K)) = SUM_XH(1,PRESSURE_ZONE(I,J,K)) + X_H(IROW)
                  SUM_XH(2,PRESSURE_ZONE(I,J,K)) = SUM_XH(2,PRESSURE_ZONE(I,J,K)) + 1._EB
               ENDDO
            ENDDO
         ENDDO
         ! Add cut-cell region contribution:
         DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
            I = CUT_CELL(ICC)%IJK(IAXIS)
            J = CUT_CELL(ICC)%IJK(JAXIS)
            K = CUT_CELL(ICC)%IJK(KAXIS)
            IF (PRESSURE_ZONE(I,J,K)<=0) CYCLE
            IROW     = MESHES(NM)%CUT_CELL(ICC)%UNKH(1) - UNKH_IND(NM_START) ! Local numeration.
            SUM_XH(1,PRESSURE_ZONE(I,J,K)) = SUM_XH(1,PRESSURE_ZONE(I,J,K)) + X_H(IROW)
            SUM_XH(2,PRESSURE_ZONE(I,J,K)) = SUM_XH(2,PRESSURE_ZONE(I,J,K)) + 1._EB
         ENDDO
      ENDDO

      IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,SUM_XH(1:2,0:N_ZONE),2*(N_ZONE+1),&
                                                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

      ! Compute arithmetic mean by pressure zone:
      DO I_ZONE=0,N_ZONE
         MEAN_XH(I_ZONE) = SUM_XH(1,I_ZONE)/(SUM_XH(2,I_ZONE)+TWO_EPSILON_EB)
      ENDDO
      ! Write out:
      ! IF (MYID==0) THEN
      !    DO I_ZONE=0,N_ZONE
      !    WRITE(LU_ERR,*) PREDICTOR,'INDEFINITE POISSON MATRIX, I_ZONE, MEAN(H), SUM(H)=',I_ZONE,MEAN_XH(I_ZONE),SUM_XH(1:2,I_ZONE)
      !    ENDDO
      ! ENDIF

      ! Substract Mean:
      DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         IF (EVACUATION_ONLY(NM) .OR. EVACUATION_SKIP(NM)) CYCLE
         CALL POINT_TO_MESH(NM)
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (CCVAR(I,J,K,UNKH)<=0 .OR. PRESSURE_ZONE(I,J,K)<=0) CYCLE ! Gasphase Cartesian cells.
                  ! Row number:
                  IROW = CCVAR(I,J,K,UNKH) - UNKH_IND(NM_START) ! Local numeration.
                  X_H(IROW) = X_H(IROW) - MEAN_XH(PRESSURE_ZONE(I,J,K))
               ENDDO
            ENDDO
         ENDDO
         ! Add cut-cell region contribution:
         DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
            I = CUT_CELL(ICC)%IJK(IAXIS)
            J = CUT_CELL(ICC)%IJK(JAXIS)
            K = CUT_CELL(ICC)%IJK(KAXIS)
            IF (PRESSURE_ZONE(I,J,K)<=0) CYCLE
            IROW     = MESHES(NM)%CUT_CELL(ICC)%UNKH(1) - UNKH_IND(NM_START) ! Local numeration.
            X_H(IROW) = X_H(IROW) - MEAN_XH(PRESSURE_ZONE(I,J,K))
         ENDDO
      ENDDO
   ELSE WHOLE_DOM_IF2
      SUM_XH(1:2,0) = SUM(X_H(1:NUNKH_LOCAL))
      IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(SUM_XH(1,0),SUM_XH(2,0),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
      MEAN_XH(0) = SUM_XH(2,0)/REAL(NUNKH_TOTAL,EB)
      ! IF (MYID==0) WRITE(LU_ERR,*) 'INDEFINITE POISSON MATRIX, MEAN(H), SUM(H)=',MEAN_XH(0),SUM_XH(2,0)
      ! Substract Mean:
      X_H(:) = X_H(:) - MEAN_XH(0)
   ENDIF WHOLE_DOM_IF2
   DEALLOCATE(SUM_FH,SUM_XH,MEAN_FH,MEAN_XH)
ENDIF

! WRITE(LU_ERR,*) 'SUM_XH=',SUM(X_H),SUM(A_H(1:IA_H(NUNKH_LOCAL+1)))
!
! IF (CORRECTOR) THEN
!    DO NM=1,NMESHES
!       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!       IF(MYID/=PROCESS(NM))CYCLE
!       CALL POINT_TO_MESH(NM)
!       WRITE(FILE_NAME,'(A,I2.2,A,I2.2,A)') "FHXH_",N_MPI_PROCESSES,'_',NMESHES,".dat"
!       IF(NM==1)THEN
!          OPEN(unit=33, file=TRIM(FILE_NAME), status='unknown')
!       ELSE
!          OPEN(unit=33, file=TRIM(FILE_NAME), status='old',position='append')
!       ENDIF
!       DO K=1,KBAR
!          DO J=1,JBAR
!             DO I=1,IBAR
!                IF(CCVAR(I,J,K,1)==1) CYCLE ! IBM_SOLID
!                IF(CCVAR(I,J,K,UNKH) > 0) THEN ! Gasphase Cartesian cells.
!                   IROW = CCVAR(I,J,K,UNKH) - UNKH_IND(NM_START)
!                ELSEIF (CCVAR(I,J,K,1)==0) THEN
!                   ICC=CCVAR(I,J,K,4)
!                   IROW= CUT_CELL(ICC)%UNKH(1) - UNKH_IND(NM_START)
!                ENDIF
!                WRITE(33,'(4I8,5F24.18)') NM,I,J,K,XC(I),YC(J),ZC(K),F_H(IROW),X_H(IROW)
!             ENDDO
!          ENDDO
!       ENDDO
!       CLOSE(33)
!    ENDDO
! ENDIF


! Dump result back to mesh containers:
MESH_LOOP_2 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION_ONLY(NM) .OR. EVACUATION_SKIP(NM)) CYCLE MESH_LOOP_2
   CALL POINT_TO_MESH(NM)

   IF (PREDICTOR) THEN
      HP => H
   ELSE
      HP => HS
   ENDIF

   ! First Source on Cartesian cells with IBM_UNKH > 0:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,UNKH) <= 0) CYCLE
            ! Row number:
            IROW = CCVAR(I,J,K,UNKH) - UNKH_IND(NM_START) ! Local numeration.
            ! Assign to HP:
            HP(I,J,K) = -X_H(IROW)
         ENDDO
      ENDDO
   ENDDO

   IF (CC_IBM) CALL GET_CUTCELL_HP(NM,HP)

   ! Fill external boundary conditions for Mesh, if necesary:
   WALL_CELL_LOOP_2: DO IW=1,N_EXTERNAL_WALL_CELLS

      WC => WALL(IW)

      ! NEUMANN boundaries:
      IF_NEUMANN2: IF (WC%PRESSURE_BC_INDEX==NEUMANN) THEN

         ! Gasphase cell indexes:
         I   = WC%ONE_D%II
         J   = WC%ONE_D%JJ
         K   = WC%ONE_D%KK
         IIG   = WC%ONE_D%IIG
         JJG   = WC%ONE_D%JJG
         KKG   = WC%ONE_D%KKG
         IOR   = WC%ONE_D%IOR

         ! Define cell size, normal to WC:
         SELECT CASE (IOR)
            CASE(-1) ! -IAXIS oriented, high face of IIG cell.
               HP(I,J,K) = HP(IIG,JJG,KKG) + DXN(IIG)*BXF(J,K)
            CASE( 1) ! +IAXIS oriented, low face of IIG cell.
               HP(I,J,K) = HP(IIG,JJG,KKG) - DXN(IIG-1)*BXS(J,K)
            CASE(-2) ! -JAXIS oriented, high face of JJG cell.
               HP(I,J,K) = HP(IIG,JJG,KKG) + DYN(JJG)*BYF(I,K)
            CASE( 2) ! +JAXIS oriented, low face of JJG cell.
               HP(I,J,K) = HP(IIG,JJG,KKG) - DYN(JJG-1)*BYS(I,K)
            CASE(-3) ! -KAXIS oriented, high face of KKG cell.
               HP(I,J,K) = HP(IIG,JJG,KKG) + DZN(KKG)*BZF(I,J)
            CASE( 3) ! +KAXIS oriented, low face of KKG cell.
               HP(I,J,K) = HP(IIG,JJG,KKG) - DZN(KKG-1)*BZS(I,J)
         END SELECT

      ENDIF IF_NEUMANN2

      ! DIRICHLET boundaries:
      IF_DIRICHLET2: IF (WC%PRESSURE_BC_INDEX==DIRICHLET) THEN

         IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
             WC%BOUNDARY_TYPE==        NULL_BOUNDARY ) CYCLE ! No need for these, that's the whole point of a
                                                             ! global solve.

         ! Gasphase cell indexes:
         I   = WC%ONE_D%II
         J   = WC%ONE_D%JJ
         K   = WC%ONE_D%KK
         IIG   = WC%ONE_D%IIG
         JJG   = WC%ONE_D%JJG
         KKG   = WC%ONE_D%KKG
         IOR   = WC%ONE_D%IOR

         ! Define cell size, normal to WC:
         IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
            SELECT CASE (IOR) ! Set Homogeneous Neumann in external SOLID_BOUNDARY.
               CASE(-1) ! -IAXIS oriented, high face of IIG cell.
                  HP(I,J,K) =HP(IIG,JJG,KKG)
               CASE( 1) ! +IAXIS oriented, low face of IIG cell.
                  HP(I,J,K) =HP(IIG,JJG,KKG)
               CASE(-2) ! -JAXIS oriented, high face of JJG cell.
                  HP(I,J,K) =HP(IIG,JJG,KKG)
               CASE( 2) ! +JAXIS oriented, low face of JJG cell.
                  HP(I,J,K) =HP(IIG,JJG,KKG)
               CASE(-3) ! -KAXIS oriented, high face of KKG cell.
                  HP(I,J,K) =HP(IIG,JJG,KKG)
               CASE( 3) ! +KAXIS oriented, low face of KKG cell.
                  HP(I,J,K) =HP(IIG,JJG,KKG)
            END SELECT
         ELSE
            SELECT CASE (IOR)
               CASE(-1) ! -IAXIS oriented, high face of IIG cell.
                  HP(I,J,K) =-HP(IIG,JJG,KKG) + 2._EB*BXF(J,K)
               CASE( 1) ! +IAXIS oriented, low face of IIG cell.
                  HP(I,J,K) =-HP(IIG,JJG,KKG) + 2._EB*BXS(J,K)
               CASE(-2) ! -JAXIS oriented, high face of JJG cell.
                  HP(I,J,K) =-HP(IIG,JJG,KKG) + 2._EB*BYF(I,K)
               CASE( 2) ! +JAXIS oriented, low face of JJG cell.
                  HP(I,J,K) =-HP(IIG,JJG,KKG) + 2._EB*BYS(I,K)
               CASE(-3) ! -KAXIS oriented, high face of KKG cell.
                  HP(I,J,K) =-HP(IIG,JJG,KKG) + 2._EB*BZF(I,J)
               CASE( 3) ! +KAXIS oriented, low face of KKG cell.
                  HP(I,J,K) =-HP(IIG,JJG,KKG) + 2._EB*BZS(I,J)
            END SELECT
         ENDIF

      ENDIF IF_DIRICHLET2

   ENDDO WALL_CELL_LOOP_2

ENDDO MESH_LOOP_2

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW

RETURN
END SUBROUTINE GLMAT_SOLVER_H

! ------------------------- GLMAT_SOLVER_SETUP ----------------------------------

SUBROUTINE GLMAT_SOLVER_SETUP_H(STAGE_FLAG)

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME

INTEGER, INTENT(IN) :: STAGE_FLAG

! Local Variables:
LOGICAL :: SUPPORTED_MESH=.TRUE.

IF (FREEZE_VELOCITY)  RETURN ! Fixed velocity soln. i.e. PERIODIC_TEST=102 => FREEZE_VELOCITY=.TRUE.
IF (SOLID_PHASE_ONLY) RETURN
TNOW=CURRENT_TIME()
SELECT CASE(STAGE_FLAG)
CASE(1)

    ! Factor to drop DY(J) in cylindrical coordinates. Soln assumes DTheta=1.
    CYL_FCT = 0._EB; IF (CYLINDRICAL) CYL_FCT = 1._EB

   ! Check for unsupported mesh configurations:
   CALL CHECK_UNSUPPORTED_MESH(SUPPORTED_MESH)
   IF (.NOT.SUPPORTED_MESH) RETURN

   ITERATE_PRESSURE = .TRUE.  ! Although there is no need to do pressure iterations to drive down velocity error
                              ! on wall cells (i.e. the solution should give the right unique dH/dxn), leave it
                              ! .TRUE. to write out velocity error diagnostics.

   ! Test for CC_IBM, define CGSC and UNKH locations in CCVAR:
   IF (CC_IBM) THEN
      CGSC   = IBM_CGSC
      UNKH   = IBM_UNKH
      NCVARS = IBM_NCVARS
   ENDIF

   ! 1. Define unknown numbers for Pressure:
   CALL SET_CCVAR_CGSC_H ! This checks if CC_IBM has been defined, and defines SOLID and GASPHASE cells accordingly,
                         ! in entry CGSC of CCVAR.

   ! Copy CCVAR CGSC to HS:
   CALL COPY_CCVAR_IN_HS(CGSC)
   ! Now go back to main.f90 and Fill Guard-cells for HS.

CASE(2)

   ! Dump back to CCVAR CGSC values from HS.
   CALL COPY_HS_IN_CCVAR(CGSC)

   CALL GET_MATRIX_INDEXES_H ! Fills Guard-cells for UNKH

   ! Copy CCVAR UNKH to HS:
   CALL COPY_CCVAR_IN_HS(UNKH)

   ! Now go back to main.f90 and Fill Guard-cells for HS.

CASE(3)

   ! Dump back to CCVAR UNKH values from HS.
   CALL COPY_HS_IN_CCVAR(UNKH)

   ! 2. For each GASPHASE (cut or regular) face, find global numeration of the volumes
   ! that share it, store a list of areas and centroids for diffussion operator in FV form.
   CALL GET_H_REGFACES

   ! 3. WALL faces have already been populated.

   ! 4. IBM_GASPHASE cut-faces:
   IF(CC_IBM) CALL GET_H_CUTFACES

   ! 5. Exchange information at block boundaries for IBM_RCFACE_H, CUT_FACE
   ! fields on each mesh:
   CALL GET_BOUNDFACE_GEOM_INFO_H

   ! 6. Get nonzeros graph of the Poisson matrix, defined as:
   !    - NNZ_D_MAT_H(1:NUNKH_LOCAL) Number of nonzeros on per matrix row.
   !    - JD_MAT_H(1:NNZ_ROW_H,1:NUNKH_LOCAL) Column location of nonzeros, global numeration.
   CALL GET_MATRIXGRAPH_H_WHLDOM ! Define the Graph of the Matrix for Gasphase cells on whole domain.

   ! 7. Build discrete Laplace operator matrix:
   CALL GET_H_MATRIX

   ! 8. Make changes to H_MATRIX due to boundary conditions (i.e. WALL faces or IBM_INBOUNDARY faces
   ! with DIRICHLET boundary condition):
   CALL GET_BCS_H_MATRIX

   ! 9. Pass D_MAT_H, NNZ_D_MAT_H, JD_MAT_H to CSR format and invoque LU solver:
   CALL GET_H_MATRIX_LUDCMP

END SELECT

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW

RETURN
END SUBROUTINE GLMAT_SOLVER_SETUP_H


! ---------------------------- CHECK_UNSUPPORTED_MESH -------------------------------

SUBROUTINE CHECK_UNSUPPORTED_MESH(SUPPORTED_MESH)

USE MPI
USE GLOBAL_CONSTANTS, ONLY : N_MPI_PROCESSES
USE TRAN, ONLY : TRANS

LOGICAL, INTENT(OUT) :: SUPPORTED_MESH

INTEGER :: NM,TRN_ME(2),IERR
REAL(EB):: DX_P(IAXIS:KAXIS),MIN_XS(3),MAX_XF(3),LX,LY,LZ
INTEGER :: COUNT
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: MESH_GRAPH,DSETS
LOGICAL, ALLOCATABLE, DIMENSION(:)   :: COUNTED
INTEGER, ALLOCATABLE, DIMENSION(:)   :: DIRI_SET,MESH_LIST
TYPE (WALL_TYPE), POINTER :: WC=>NULL()

INTEGER :: NOM,IW,NMLOC,NSETS,ISET,PIVOT,PIVOT_LOC,MESHES_LEFT,CTMSH_LO,CTMSH_HI

SUPPORTED_MESH = .TRUE.

! 1. Stretched grids which is untested:
TRN_ME(1:2) = 0
MESH_LOOP_TRN : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (EVACUATION_ONLY(NM) .OR. EVACUATION_SKIP(NM)) CYCLE MESH_LOOP_TRN
   TRN_ME(1) = TRN_ME(1) + TRANS(NM)%NOCMAX
ENDDO MESH_LOOP_TRN
TRN_ME(2)=TRN_ME(1)
IF (N_MPI_PROCESSES > 1) CALL MPI_ALLREDUCE(TRN_ME(1),TRN_ME(2),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
IF (TRN_ME(2) > 0) THEN ! There is a TRNX, TRNY or TRNZ line defined for stretched grids. Not Unsupported.
   IF (MYID == 0) WRITE(LU_ERR,*) 'GLMAT Setup Error : Stretched grids currently unsupported.'
   SUPPORTED_MESH = .FALSE.
   STOP_STATUS = SETUP_STOP
   RETURN
ENDIF

IF (NMESHES == 1) RETURN

! 2. Two different cell sizes in mesh (i.e. different refinement levels):
NM = 1
IF (MYID==PROCESS(NM)) THEN
   CALL POINT_TO_MESH(NM)
   DX_P(IAXIS) = DX(1)
   DX_P(JAXIS) = DY(1)
   DX_P(KAXIS) = DZ(1)
ENDIF
IF (N_MPI_PROCESSES > 1) CALL MPI_BCAST(DX_P,3,MPI_DOUBLE_PRECISION,PROCESS(NM),MPI_COMM_WORLD,IERR)
! Find domain sizes to define relative epsilon:
MIN_XS(1:3) = (/ MESHES(NM)%XS, MESHES(NM)%YS, MESHES(NM)%ZS /)
MAX_XF(1:3) = (/ MESHES(NM)%XF, MESHES(NM)%YF, MESHES(NM)%ZF /)
DO NM=2,NMESHES
   MIN_XS(1) = MIN(MIN_XS(1),MESHES(NM)%XS)
   MIN_XS(2) = MIN(MIN_XS(2),MESHES(NM)%YS)
   MIN_XS(3) = MIN(MIN_XS(3),MESHES(NM)%ZS)
   MAX_XF(1) = MAX(MAX_XF(1),MESHES(NM)%XF)
   MAX_XF(2) = MAX(MAX_XF(2),MESHES(NM)%YF)
   MAX_XF(3) = MAX(MAX_XF(3),MESHES(NM)%ZF)
ENDDO
LX = MAX(MAX_XF(1)-MIN_XS(1),1._EB)
LY = MAX(MAX_XF(2)-MIN_XS(2),1._EB)
LZ = MAX(MAX_XF(3)-MIN_XS(3),1._EB)
TRN_ME(1:2) = 0
MESH_LOOP_CELL : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (EVACUATION_ONLY(NM) .OR. EVACUATION_SKIP(NM)) CYCLE MESH_LOOP_CELL
   CALL POINT_TO_MESH(NM)
   IF(ABS(DX_P(IAXIS)-DX(1)) > 10._EB*TWO_EPSILON_EB*LX) TRN_ME(1) = TRN_ME(1) + 1
   IF(ABS(DX_P(JAXIS)-DY(1)) > 10._EB*TWO_EPSILON_EB*LY) TRN_ME(1) = TRN_ME(1) + 1
   IF(ABS(DX_P(KAXIS)-DZ(1)) > 10._EB*TWO_EPSILON_EB*LZ) TRN_ME(1) = TRN_ME(1) + 1
ENDDO MESH_LOOP_CELL
TRN_ME(2)=TRN_ME(1)
IF (N_MPI_PROCESSES > 1) CALL MPI_ALLREDUCE(TRN_ME(1),TRN_ME(2),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
IF (TRN_ME(2) > 0) THEN ! Meshes at different refinement levels. Not Unsupported.
   IF (MYID == 0) WRITE(LU_ERR,*) 'GLMAT Setup Error : Meshes at different refinement levels unsupported.'
   SUPPORTED_MESH = .FALSE.
   STOP_STATUS = SETUP_STOP
   RETURN
ENDIF

! 3. Two (or more) disjoint domains, where at least one has all Neumann BCs and one has some Dirichlet bcs.
! This is a topological problem that would require different Matrix types (i.e. one positive definite and one
! indefinite), which would require separate solutions.
! A possible approach to look at is to solve the whole system as indefinite, and then substract a constant in
! zones with Dirichlet condition, s.t. the value of H is zero in open boundaries.

! 1. Build global lists of other connected meshes:
ALLOCATE(MESH_GRAPH(1:6,NMESHES)); MESH_GRAPH(:,:) = 0
MESH_LOOP_GRAPH : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   COUNT=0
   IF (EVACUATION_ONLY(NM) .OR. EVACUATION_SKIP(NM)) CYCLE MESH_LOOP_GRAPH
   DO NOM=1,NMESHES
      IF(MESHES(NM)%CONNECTED_MESH(NOM))THEN
         COUNT=COUNT+1
         MESH_GRAPH(COUNT,NM) = NOM
      ENDIF
   ENDDO
ENDDO MESH_LOOP_GRAPH
IF (N_MPI_PROCESSES > 1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,MESH_GRAPH,6*NMESHES,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)

! 2. Build sets of disjoint meshes:
! Number of sets:
ALLOCATE(COUNTED(NMESHES));  COUNTED(1:NMESHES)  =.FALSE.
ALLOCATE(DSETS(2,NMESHES));  DSETS(1:2,1:NMESHES)= 0
ALLOCATE(MESH_LIST(NMESHES));MESH_LIST(1:NMESHES)= 0
NSETS    = 1
CTMSH_LO = 1
CTMSH_HI = CTMSH_LO
PIVOT    = 1
COUNTED(PIVOT)       = .TRUE.
DSETS(LOW_IND,NSETS) =  CTMSH_LO
MESH_LIST(CTMSH_LO)  =  PIVOT
MESHES_LEFT          = NMESHES-1
DISJ_LOOP : DO

   DO NMLOC=1,6
      PIVOT_LOC = MESH_GRAPH(NMLOC,PIVOT)
      IF(PIVOT_LOC==0) CYCLE ! Cycle is other mesh not present.
      IF(COUNTED(PIVOT_LOC)) CYCLE ! Cycle if mesh has been already added to mesh list.
      ! Add mesh to list:
      CTMSH_HI = CTMSH_HI + 1
      COUNTED(PIVOT_LOC)   = .TRUE.
      MESH_LIST(CTMSH_HI)  =  PIVOT_LOC
      MESHES_LEFT          = MESHES_LEFT-1
   ENDDO

   ! No more new meshes on the set, increase NSETS by one and start again:
   IF (CTMSH_LO == CTMSH_HI) THEN
      DSETS(HIGH_IND,NSETS) = CTMSH_HI
      IF (MESHES_LEFT==0) EXIT DISJ_LOOP ! Done with all meshes, finish and exit.
      DO NM=1,NMESHES
         IF (EVACUATION_ONLY(NM) .OR. EVACUATION_SKIP(NM)) THEN
            MESHES_LEFT          = MESHES_LEFT-1
            CYCLE DISJ_LOOP
         ENDIF
         IF(.NOT.COUNTED(NM))THEN
            NSETS=NSETS+1
            CTMSH_LO = CTMSH_LO + 1
            CTMSH_HI = CTMSH_LO
            PIVOT    = NM
            COUNTED(PIVOT)       = .TRUE.
            DSETS(LOW_IND,NSETS) = CTMSH_LO
            MESH_LIST(CTMSH_LO)  = PIVOT
            MESHES_LEFT          = MESHES_LEFT-1
            CYCLE DISJ_LOOP ! This is such that we don't increase CTMSH_LO when starting a new set.
         ENDIF
      ENDDO
   ENDIF

   ! Increase CTMSH_LO by one:
   CTMSH_LO = CTMSH_LO + 1
   PIVOT    = MESH_LIST(CTMSH_LO)

ENDDO DISJ_LOOP

! If only one set, topology supported, return:
IF (NSETS==1) THEN
   DEALLOCATE(MESH_GRAPH,DSETS,MESH_LIST,COUNTED)
   RETURN
ENDIF

! 3. Check for each set of disjoint meshes that all of them have at least one Pressure Dirichlet boundary condition
! (open boundary). If not, return SETUP stop:
ALLOCATE(DIRI_SET(NSETS)); DIRI_SET(:)=0
SETS_LOOP : DO ISET=1,NSETS
   DO NMLOC=DSETS(LOW_IND,ISET),DSETS(HIGH_IND,ISET)
      NM=MESH_LIST(NMLOC)
      IF (MYID/=PROCESS(NM)) CYCLE
      IF (EVACUATION_ONLY(NM) .OR. EVACUATION_SKIP(NM)) CYCLE SETS_LOOP
      CALL POINT_TO_MESH(NM)

      ! Now for Mesh NM test for Dirichlet Pressure external BCs. Assume External Dirichlet is only related to
      ! OPEN_BOUNDARY condition.
      DO IW=1,N_EXTERNAL_WALL_CELLS
         WC => WALL(IW)
         IF (WC%PRESSURE_BC_INDEX==DIRICHLET .AND. WC%BOUNDARY_TYPE==OPEN_BOUNDARY) DIRI_SET(ISET) = 1
      ENDDO
   ENDDO
ENDDO SETS_LOOP

IF(N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,DIRI_SET,NSETS,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)

! IF (MYID==0) THEN
!    WRITE(LU_ERR,*) ' '
!    WRITE(LU_ERR,*) 'NSETS=',NSETS
!    DO ISET=1,NSETS
!       WRITE(LU_ERR,*) 'ISET=',ISET,DSETS(HIGH_IND,ISET)-DSETS(LOW_IND,ISET)+1,DIRI_SET(ISET)
!    ENDDO
! ENDIF

! Finally do test:
IF (ANY(DIRI_SET(1:NSETS) == 0)) THEN
   IF (MYID==0) WRITE(LU_ERR,*) 'GLMAT Setup Error : Unsupported disjoint domains present on the model.'
   DEALLOCATE(MESH_GRAPH,DSETS,MESH_LIST,COUNTED,DIRI_SET)
   SUPPORTED_MESH = .FALSE.
   STOP_STATUS = SETUP_STOP
   RETURN
ENDIF

DEALLOCATE(MESH_GRAPH,DSETS,MESH_LIST,COUNTED,DIRI_SET)
RETURN
END SUBROUTINE CHECK_UNSUPPORTED_MESH

! ----------------------------- COPY_H_OMESH_TO_MESH --------------------------------

SUBROUTINE COPY_H_OMESH_TO_MESH

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME

! Local Variables:
INTEGER  :: NM,NOM,II,JJ,KK,IOR,IW,IIO,JJO,KKO
TYPE (OMESH_TYPE), POINTER :: OM=>NULL()
TYPE (WALL_TYPE), POINTER :: WC=>NULL()
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()
LOGICAL :: FLG

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN
TNOW=CURRENT_TIME()
! Loop:
PREDCORR_LOOP : IF (PREDICTOR) THEN

   MESH_LOOP_1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

      IF (EVACUATION_ONLY(NM)) CYCLE MESH_LOOP_1
      CALL POINT_TO_MESH(NM)

      ! Loop over all cell edges

      EXTERNAL_WALL_LOOP_1 : DO IW=1,N_EXTERNAL_WALL_CELLS

         WC=>WALL(IW)
         EWC=>EXTERNAL_WALL(IW)
         IF (PRES_ON_WHOLE_DOMAIN) THEN
            ! Matrix connectivities kept in cases of BOUNDARY_TYPE=INTERPOLATED_BOUNDARY, NULL_BOUNDARY, or
            ! SOLID_BOUNDARY where there is an OMESH, case of OBSTS right in the boundary of meshes.
            FLG = WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY
            FLG = FLG .OR. (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. EWC%NOM > 0)
            IF (.NOT.FLG) CYCLE EXTERNAL_WALL_LOOP_1
         ELSE
            ! Case of solving for H only on the gas phase, connectivity kept only for INTERPOLATED_BOUNDARY.
            IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP_1
         ENDIF

         II  = WC%ONE_D%II
         JJ  = WC%ONE_D%JJ
         KK  = WC%ONE_D%KK
         IOR = WC%ONE_D%IOR
         NOM = EWC%NOM
         ! Here if NOM==0 means it is an OBST laying on an external boundary -> CYCLE
         IF(NOM < 1) CYCLE
         OM => OMESH(NOM)

         ! This assumes all meshes at the same level of refinement:
         KKO=EWC%KKO_MIN
         JJO=EWC%JJO_MIN
         IIO=EWC%IIO_MIN

         SELECT CASE(IOR)
         CASE( IAXIS)
            H(II,JJ,KK) = OM%H(IIO,JJO,KKO)
         CASE(-IAXIS)
            H(II,JJ,KK) = OM%H(IIO,JJO,KKO)
         CASE( JAXIS)
            H(II,JJ,KK) = OM%H(IIO,JJO,KKO)
         CASE(-JAXIS)
            H(II,JJ,KK) = OM%H(IIO,JJO,KKO)
         CASE( KAXIS)
            H(II,JJ,KK) = OM%H(IIO,JJO,KKO)
         CASE(-KAXIS)
            H(II,JJ,KK) = OM%H(IIO,JJO,KKO)
         END SELECT

      ENDDO EXTERNAL_WALL_LOOP_1

   ENDDO MESH_LOOP_1

ELSE ! PREDCORR_LOOP

   MESH_LOOP_2 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

      IF (EVACUATION_ONLY(NM)) CYCLE MESH_LOOP_2
      CALL POINT_TO_MESH(NM)

      ! Loop over all cell edges

      EXTERNAL_WALL_LOOP_2 : DO IW=1,N_EXTERNAL_WALL_CELLS

         WC=>WALL(IW)
         EWC=>EXTERNAL_WALL(IW)
         IF (PRES_ON_WHOLE_DOMAIN) THEN
            ! Matrix connectivities kept in cases of BOUNDARY_TYPE=INTERPOLATED_BOUNDARY, NULL_BOUNDARY, or
            ! SOLID_BOUNDARY where there is an OMESH, case of OBSTS right in the boundary of meshes.
            FLG = WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY
            FLG = FLG .OR. (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. EWC%NOM > 0)
            IF (.NOT.FLG) CYCLE EXTERNAL_WALL_LOOP_2
         ELSE
            ! Case of solving for H only on the gas phase, connectivity kept only for INTERPOLATED_BOUNDARY.
            IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP_2
         ENDIF

         II  = WC%ONE_D%II
         JJ  = WC%ONE_D%JJ
         KK  = WC%ONE_D%KK
         IOR = WC%ONE_D%IOR
         NOM = EWC%NOM
         ! Here if NOM==0 means it is an OBST laying on an external boundary -> CYCLE
         IF(NOM < 1) CYCLE
         OM => OMESH(NOM)

         ! This assumes all meshes at the same level of refinement:
         KKO=EWC%KKO_MIN
         JJO=EWC%JJO_MIN
         IIO=EWC%IIO_MIN

         SELECT CASE(IOR)
         CASE( IAXIS)
            HS(II,JJ,KK) = OM%HS(IIO,JJO,KKO)
         CASE(-IAXIS)
            HS(II,JJ,KK) = OM%HS(IIO,JJO,KKO)
         CASE( JAXIS)
            HS(II,JJ,KK) = OM%HS(IIO,JJO,KKO)
         CASE(-JAXIS)
            HS(II,JJ,KK) = OM%HS(IIO,JJO,KKO)
         CASE( KAXIS)
            HS(II,JJ,KK) = OM%HS(IIO,JJO,KKO)
         CASE(-KAXIS)
            HS(II,JJ,KK) = OM%HS(IIO,JJO,KKO)
         END SELECT

      ENDDO EXTERNAL_WALL_LOOP_2

   ENDDO MESH_LOOP_2

ENDIF PREDCORR_LOOP

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW

RETURN
END SUBROUTINE COPY_H_OMESH_TO_MESH

! ------------------------------- COPY_HS_IN_CCVAR ----------------------------------

SUBROUTINE COPY_HS_IN_CCVAR(VAR_CC)

INTEGER, INTENT(IN) :: VAR_CC

! Local Variables:
INTEGER  :: NM,NOM,II,JJ,KK,IOR,IW,IIO,JJO,KKO
TYPE (OMESH_TYPE), POINTER :: OM=>NULL()
TYPE (WALL_TYPE), POINTER :: WC=>NULL()
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()
LOGICAL :: FLG

MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION_ONLY(NM)) CYCLE MESH_LOOP
   CALL POINT_TO_MESH(NM)

   !
   IF (CC_IBM .AND. VAR_CC==UNKH) THEN

      CALL COPY_CC_HS_TO_UNKH(NM)

   ELSE

      ! Loop over external wall cells:
      EXTERNAL_WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

         WC=>WALL(IW)
         EWC=>EXTERNAL_WALL(IW)
         IF (PRES_ON_WHOLE_DOMAIN) THEN
            ! Matrix connectivities kept in cases of BOUNDARY_TYPE=INTERPOLATED_BOUNDARY, NULL_BOUNDARY, or
            ! SOLID_BOUNDARY where there is an OMESH, case of OBSTS right in the boundary of meshes.
            FLG = WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY
            FLG = FLG .OR. (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. EWC%NOM > 0)
            IF (.NOT.FLG) CYCLE EXTERNAL_WALL_LOOP
         ELSE
            ! Case of solving for H only on the gas phase, connectivity kept only for INTERPOLATED_BOUNDARY.
            IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP
         ENDIF

         II  = WC%ONE_D%II
         JJ  = WC%ONE_D%JJ
         KK  = WC%ONE_D%KK
         IOR = WC%ONE_D%IOR
         NOM = EWC%NOM
         ! Here if NOM==0 means it is an OBST laying on an external boundary -> CYCLE
         IF(NOM < 1) CYCLE
         OM => OMESH(NOM)

         ! This assumes all meshes at the same level of refinement:
         KKO=EWC%KKO_MIN
         JJO=EWC%JJO_MIN
         IIO=EWC%IIO_MIN

         CCVAR(II,JJ,KK,VAR_CC) = INT(OM%HS(IIO,JJO,KKO))

         IF (VAR_CC == UNKH) OM%HS(IIO,JJO,KKO) = 0._EB

      ENDDO EXTERNAL_WALL_LOOP

   ENDIF

   IF (VAR_CC == UNKH) THEN
      HS(:,:,:)    = 0._EB
   ENDIF

ENDDO MESH_LOOP

RETURN
END SUBROUTINE COPY_HS_IN_CCVAR


! ------------------------------- COPY_CCVAR_IN_HS ----------------------------------

SUBROUTINE COPY_CCVAR_IN_HS(VAR_CC)

INTEGER, INTENT(IN) :: VAR_CC

! Local Variables:
INTEGER :: NM

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (EVACUATION_ONLY(NM)) CYCLE
   CALL POINT_TO_MESH(NM)
   HS(0:IBP1,0:JBP1,0:KBP1) = REAL(CCVAR(0:IBP1,0:JBP1,0:KBP1,VAR_CC),EB)

   ! Now cut-cells add their single Cartesian UNKH value in HS:
   IF(CC_IBM .AND. VAR_CC==UNKH) CALL COPY_CC_UNKH_TO_HS(NM)

ENDDO

RETURN
END SUBROUTINE COPY_CCVAR_IN_HS

! ------------------------------- GET_H_MATRIX_LUDCMP -------------------------------

SUBROUTINE GET_H_MATRIX_LUDCMP

USE MPI

! Local Variables:
INTEGER :: INNZ, IROW, JCOL
#ifdef WITH_MKL
INTEGER :: PHASE, PERM(1)
INTEGER :: I, IPROC
#endif
!.. All other variables
INTEGER MAXFCT, MNUM, MTYPE, NRHS, ERROR
#ifdef WITH_MKL
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: MB_FACTOR
INTEGER :: IERR
#endif

! Define parameters:
NRHS   = 1
MAXFCT = 1
MNUM   = 1

! Set level MSG to 1 for factorization:
IF(GLMAT_VERBOSE) MSGLVL = 1

! Define control parameter vector iparm:
ALLOCATE(IPARM(64)); IPARM(:) = 0

IPARM(1) = 1   ! no solver default
! Pardiso: IPARM(2) = 2   ! fill-in reordering from METIS
#ifdef WITH_MKL
IF (N_MPI_PROCESSES > 4) THEN ! Typical number of computing cores inside one chip.
   IPARM(2) =10   ! 10 = MPI Parallel fill-in reordering from METIS. If 3 = OpenMP parallel reordering in Master Node.
ELSE              ! Note IPARM(2)=10 has a bug which has been fixed from Intel MKL 2018 update 2 onwards.
   IPARM(2) = 3
ENDIF
#endif
IPARM(4) = 0   ! no iterative-direct algorithm
IPARM(5) = 0   ! no user fill-in reducing permutation
IPARM(6) = 0   ! =0 solution on the first n components of x
IPARM(8) = 2   ! numbers of iterative refinement steps
IPARM(10) = 13 ! perturb the pivot elements with 1E-13
IPARM(11) = 1  ! use nonsymmetric permutation and scaling MPS  !!!!! was 1
IPARM(13) = 1  ! maximum weighted matching algorithm is switched-off
               !(default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
IPARM(14) = 0  ! Output: number of perturbed pivots
IPARM(18) = 0 !-1 ! Output: number of nonzeros in the factor LU
IPARM(19) = 0 !-1 ! Output: Mflops for LU factorization
IPARM(20) = 0  ! Output: Numbers of CG Iterations

IPARM(21) = 1  ! 1x1 diagonal pivoting for symmetric indefinite matrices.

IPARM(24) = 0

IPARM(27) = 1 ! Check matrix

IPARM(40) = 2 ! Matrix, solution and rhs provided in distributed assembled matrix input format.

ERROR     = 0 ! initialize error flag

! Each MPI process builds its local set of rows.
! Matrix blocks defined on CRS distributed format.
! Total number of nonzeros for JD_MAT_H, D_MAT_H:
TOT_NNZ_H = sum( NNZ_D_MAT_H(1:NUNKH_LOCAL) )

! Allocate A_H IA_H and JA_H matrices, considering all matrix coefficients:
ALLOCATE ( A_H(TOT_NNZ_H) , IA_H(NUNKH_LOCAL+1) , JA_H(TOT_NNZ_H) )

! Store upper triangular part of symmetric D_MAT_H in CSR format:
INNZ = 0
DO IROW=1,NUNKH_LOCAL
   IA_H(IROW) = INNZ + 1
   DO JCOL=1,NNZ_D_MAT_H(IROW)
      IF ( JD_MAT_H(JCOL,IROW) < UNKH_IND(NM_START)+IROW ) CYCLE ! Only upper Triangular part.
      INNZ = INNZ + 1
      A_H(INNZ)  =  D_MAT_H(JCOL,IROW)
      JA_H(INNZ) = JD_MAT_H(JCOL,IROW)
   ENDDO
ENDDO
IA_H(NUNKH_LOCAL+1) = INNZ + 1

! OPEN(unit=20,file="Matrix_H.txt",action="write",status="replace")
! DO IROW=1,NUNKH_LOCAL
!    DO JCOL=1,NNZ_D_MAT_H(IROW)
!       WRITE(20,'(2I6,F18.12)') IROW,JD_MAT_H(JCOL,IROW),D_MAT_H(JCOL,IROW)
!    ENDDO
! ENDDO
! ! WRITE(20,'(A)') 'EOF'
! CLOSE(20)
! WRITE(0,*) 'H Matrix file written...'
! PAUSE

! Here each process defines de beginning and end rows in global numeration, for the equations
! it has assembled:
IPARM(41) = UNKH_IND(NM_START) + 1
IPARM(42) = UNKH_IND(NM_START) + NUNKH_LOCAL

IF ( H_MATRIX_INDEFINITE ) THEN
   MTYPE  = -2 ! symmetric indefinite
ELSE ! positive definite
   MTYPE  =  2
ENDIF

ALLOCATE( X_H(NUNKH_LOCAL) , F_H(NUNKH_LOCAL) ) ! JUST ZERO FOR NOW.
F_H(:) = 0._EB
X_H(:) = 0._EB

NUNKH_TOTAL = sum(NUNKH_TOT(1:NMESHES))

ALLOCATE(PT_H(64))

! PARDISO:
! Initialize solver pointer for H matrix solves:
! DO I=1,64
!   PT_H(I)%DUMMY = 0
! ENDDO

! Reorder and Symbolic factorization:
! PHASE = 11
! CALL PARDISO (PT_H, MAXFCT, MNUM, MTYPE, PHASE, NUNKH_TOTAL, &
!     A_H, IA_H, JA_H, PERM, NRHS, IPARM, MSGLVL, F_H, X_H, ERROR)
!
! IF (ERROR /= 0) THEN
!    IF (MYID==0) &
!    WRITE(LU_ERR,'(A,I5)') 'GET_H_MATRIX_LUDCMP PARDISO Sym Factor: The following ERROR was detected: ', ERROR
!    ! Some error - stop flag for CALL STOP_CHECK(1).
!    STOP_STATUS = SETUP_STOP
!    RETURN
! END IF

! Numerical Factorization.
! PHASE = 22 ! only factorization
! CALL PARDISO (PT_H, MAXFCT, MNUM, MTYPE, PHASE, NUNKH_TOTAL, &
!   A_H, IA_H, JA_H, PERM, NRHS, IPARM, MSGLVL, F_H, X_H, ERROR)
!
! IF (ERROR /= 0) THEN
!    IF (MYID==0) &
!    WRITE(LU_ERR,'(A,I5)') 'GET_H_MATRIX_LUDCMP PARDISO Num Factor: The following ERROR was detected: ', ERROR
!    ! Some error - stop flag for CALL STOP_CHECK(1).
!    STOP_STATUS = SETUP_STOP
!    RETURN
! ENDIF

! Define 4 byte A_H, F_H and X_H:
#ifdef SINGLE_PRECISION_PSN_SOLVE
IPARM(28) = 1 ! Single Precision solve.
ALLOCATE(F_H_FB(1:NUNKH_LOCAL)); F_H_FB(1:NUNKH_LOCAL) = 0._FB
ALLOCATE(X_H_FB(1:NUNKH_LOCAL)); X_H_FB(1:NUNKH_LOCAL) = 0._FB
ALLOCATE( A_H_FB(TOT_NNZ_H) );   A_H_FB(1:TOT_NNZ_H)   = REAL(A_H(1:TOT_NNZ_H),FB)
#endif

#ifdef WITH_MKL
! Initialize solver pointer for H matrix solves:
DO I=1,64
  PT_H(I)%DUMMY = 0
ENDDO

! Reorder and Symbolic factorization:
PHASE = 11
#ifdef __INTEL_COMPILER
   CALL KMP_SET_WARNINGS_OFF()
#endif
#ifdef SINGLE_PRECISION_PSN_SOLVE
CALL CLUSTER_SPARSE_SOLVER (PT_H, MAXFCT, MNUM, MTYPE, PHASE, NUNKH_TOTAL, &
    A_H_FB, IA_H, JA_H, PERM, NRHS, IPARM, MSGLVL, F_H_FB, X_H_FB, MPI_COMM_WORLD, ERROR)
#else
CALL CLUSTER_SPARSE_SOLVER (PT_H, MAXFCT, MNUM, MTYPE, PHASE, NUNKH_TOTAL, &
    A_H, IA_H, JA_H, PERM, NRHS, IPARM, MSGLVL, F_H, X_H, MPI_COMM_WORLD, ERROR)
#endif
#ifdef __INTEL_COMPILER
   CALL KMP_SET_WARNINGS_ON()
#endif

IF (ERROR /= 0) THEN
   IF (MYID==0) THEN
   WRITE(LU_ERR,'(A,I5)') 'GET_H_MATRIX_LUDCMP CLUSTER_SOLVER Sym Factor: The following ERROR was detected: ', ERROR
   IF(ERROR == -4) THEN
      WRITE(LU_ERR,'(A,A)') 'This error is probably due to having one or more sealed compartments ',&
      'besides a pressure zone with open boundary. Currently this situation is not supported.'
   ELSEIF(ERROR == -2) THEN
      WRITE(LU_ERR,'(A)') 'Not enough physical memory in your system for factoring the Poisson Matrix.'
   ENDIF
   ENDIF
   ! Some error - stop flag for CALL STOP_CHECK(1).
   STOP_STATUS = SETUP_STOP
   RETURN
END IF

! Define array of MB required by process, gather data to Master:
ALLOCATE(MB_FACTOR(2,0:N_MPI_PROCESSES-1)); MB_FACTOR(:,:)=0
MB_FACTOR(1,MYID) = MAX(IPARM(15),IPARM(16)+IPARM(17))/1000
MB_FACTOR(2,MYID) = NUNKH_LOCAL
IF(N_MPI_PROCESSES > 1) &
CALL MPI_ALLREDUCE(MPI_IN_PLACE, MB_FACTOR, 2*N_MPI_PROCESSES, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
! Write to output file:
IF(MYID==0) THEN
   IPROC=MAXLOC(MB_FACTOR(1,0:N_MPI_PROCESSES-1),DIM=1) - 1 ! MaxLoc defines which element in the array, not index.
   WRITE(LU_OUTPUT,*) '   MPI Process, H unknowns =',IPROC,MB_FACTOR(2,IPROC), &
      ', Peak Factorization Memory Required (MB)=',MB_FACTOR(1,IPROC)
ENDIF
! Here do Memory test trying to allocate an array?
! Difficult to do, nodes have varying numbers of MPI_PROCESSES and RAM.
DEALLOCATE(MB_FACTOR)

! Numerical Factorization.
PHASE = 22 ! only factorization
#ifdef SINGLE_PRECISION_PSN_SOLVE
CALL CLUSTER_SPARSE_SOLVER (PT_H, MAXFCT, MNUM, MTYPE, PHASE, NUNKH_TOTAL, &
    A_H_FB, IA_H, JA_H, PERM, NRHS, IPARM, MSGLVL, F_H_FB, X_H_FB, MPI_COMM_WORLD, ERROR)
#else
CALL CLUSTER_SPARSE_SOLVER (PT_H, MAXFCT, MNUM, MTYPE, PHASE, NUNKH_TOTAL, &
  A_H, IA_H, JA_H, PERM, NRHS, IPARM, MSGLVL, F_H, X_H, MPI_COMM_WORLD, ERROR)
#endif

IF (ERROR /= 0) THEN
   IF (MYID==0) THEN
   WRITE(LU_ERR,'(A,I5)') 'GET_H_MATRIX_LUDCMP CLUSTER_SOLVER Num Factor: The following ERROR was detected: ', ERROR
   IF(ERROR == -4) THEN
      WRITE(LU_ERR,'(A,A)') 'This error is probably due to having one or more sealed compartments ',&
      ' besides a compartment with/without open boundary. Currently only one pressure zone is supported.'
   ELSEIF(ERROR == -2) THEN
      WRITE(LU_ERR,'(A)') 'Not enough physical memory in your system for factoring the Poisson Matrix.'
   ENDIF
   ENDIF
   ! Some error - stop flag for CALL STOP_CHECK(1).
   STOP_STATUS = SETUP_STOP
   RETURN
ENDIF

#else

IF (MYID==0) WRITE(LU_ERR,'(A)') &
'Error: MKL Library compile flag was not defined for GLMAT as pressure solver.'
! Some error - stop flag for CALL STOP_CHECK(1).
STOP_STATUS = SETUP_STOP
RETURN

#endif

! Set level MSG to 0 for solution:
IF(GLMAT_VERBOSE) MSGLVL = 0

END SUBROUTINE GET_H_MATRIX_LUDCMP

! -------------------------------- GET_BCS_H_MATRIX ---------------------------------

SUBROUTINE GET_BCS_H_MATRIX

USE MPI
USE COMPLEX_GEOMETRY, ONLY : GET_CC_UNKH

! Local Variables:
INTEGER :: NM,NM1
INTEGER :: JLOC,JCOL,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND)
INTEGER :: DIRI_SUM,IERR
REAL(EB):: AF,IDX,BIJ
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
INTEGER :: IIG,JJG,KKG,II,JJ,KK,IW
REAL(EB), POINTER, DIMENSION(:,:) :: D_MAT_HP
INTEGER, ALLOCATABLE, DIMENSION(:) :: H_MAT_IVEC, H_MAT_IVEC_AUX

ALLOCATE( H_MAT_IVEC(1:NMESHES) ); H_MAT_IVEC = 0

! Allocate D_MAT_H:
D_MAT_HP => D_MAT_H
NM1 = NM_START

! Main Mesh Loop:
MESH_LOOP_1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION_ONLY(NM)) CYCLE MESH_LOOP_1
   CALL POINT_TO_MESH(NM)

   WALL_LOOP_1 : DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC => WALL(IW)

      ! Only OPEN_BOUNDARY leads to a Dirichlet BC for H when we solve the problem on the whole
      ! unstructured domain. Everything else leads to Neuman BCs on H, no need to modify D_MAT_HP.
      IF ( WC%BOUNDARY_TYPE/=OPEN_BOUNDARY ) CYCLE

      IIG = WC%ONE_D%IIG; JJG = WC%ONE_D%JJG; KKG = WC%ONE_D%KKG
      II  = WC%ONE_D%II;  JJ  = WC%ONE_D%JJ;  KK  = WC%ONE_D%KK
      ! Unknowns on related cells:
      IF(CCVAR(IIG,JJG,KKG,CGSC)==IS_GASPHASE .OR. PRES_ON_WHOLE_DOMAIN) THEN
         IND(LOW_IND)  = CCVAR(IIG,JJG,KKG,UNKH)  ! internal cell.
      ELSEIF(CCVAR(IIG,JJG,KKG,CGSC)==IS_CUTCFE) THEN
         CALL GET_CC_UNKH(IIG,JJG,KKG,IND(LOW_IND))
      ELSE
         CYCLE ! Solid cell and .NOT.PRES_ON_WHOLE_DOMAIN
      ENDIF

      IND_LOC(LOW_IND) = IND(LOW_IND) - UNKH_IND(NM1) ! All row indexes must refer to ind_loc.
      SELECT CASE(WC%ONE_D%IOR)
      CASE( IAXIS)
         AF = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG-1)) * DZ(KKG);            IDX= 1._EB/DXN(IIG-1)
      CASE(-IAXIS)
         AF = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG  )) * DZ(KKG);            IDX= 1._EB/DXN(IIG  )
      CASE( JAXIS)
         AF = DX(IIG)*DZ(KKG);            IDX= 1._EB/DYN(JJG-1)
      CASE(-JAXIS)
         AF = DX(IIG)*DZ(KKG);            IDX= 1._EB/DYN(JJG  )
      CASE( KAXIS)
         AF = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG);            IDX= 1._EB/DZN(KKG-1)
      CASE(-KAXIS)
         AF = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG);            IDX= 1._EB/DZN(KKG  )
      END SELECT

      ! Now add to Adiff corresponding coeff:
      BIJ = IDX*AF
      ! Find diagonal column number:
      JCOL = -1
      DO JLOC = 1,NNZ_D_MAT_H(IND_LOC(LOW_IND))
         IF (IND(LOW_IND) == JD_MAT_H(JLOC,IND_LOC(LOW_IND))) THEN
            JCOL = JLOC
            EXIT
         ENDIF
      ENDDO
      ! Add diagonal coefficient due to DIRICHLET BC:
      D_MAT_HP(JCOL,IND_LOC(LOW_IND)) = D_MAT_HP(JCOL,IND_LOC(LOW_IND)) + 2._EB*BIJ

      ! Add to mesh dirichlet bc counter
      H_MAT_IVEC(NM) = H_MAT_IVEC(NM) + 1

   ENDDO WALL_LOOP_1

   ! Dirchlet bcs in IBM_INBOUNDARY faces:
   ! Might not be needed.

ENDDO MESH_LOOP_1


! Is the resulting Matrix Indefinite?
! Here all reduce with sum among MPI processes:
ALLOCATE( H_MAT_IVEC_AUX(1:NMESHES) )
! CALL MPI_ALLREDUCE:
IF (N_MPI_PROCESSES > 1) THEN
 H_MAT_IVEC_AUX = H_MAT_IVEC
 CALL MPI_ALLREDUCE(H_MAT_IVEC_AUX, H_MAT_IVEC, NMESHES, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
ENDIF
DIRI_SUM = SUM(H_MAT_IVEC(1:NMESHES))

H_MATRIX_INDEFINITE = .TRUE. ! Set this value regarding H BCs.
IF ( DIRI_SUM > 0 ) H_MATRIX_INDEFINITE = .FALSE. ! At least one Dirichlet BC, matrix positive definite.
DEALLOCATE(H_MAT_IVEC_AUX)

DEALLOCATE(H_MAT_IVEC)

RETURN
END SUBROUTINE GET_BCS_H_MATRIX

! --------------------------------- GET_H_MATRIX ------------------------------------

SUBROUTINE GET_H_MATRIX

USE COMPLEX_GEOMETRY, ONLY : GET_H_MATRIX_CC

! Local Variables:
INTEGER :: NM,NM1,NREG
INTEGER :: LOW_FACE,HIGH_FACE,X1AXIS,X2AXIS,X3AXIS,IFACE
REAL(EB), POINTER, DIMENSION(:,:) :: D_MAT_HP
REAL(EB), POINTER, DIMENSION(:)   :: DX1,DX2,DX3
INTEGER :: I,J,K,I1,I2,I3,IIP,JJP,KKP,IIM,JJM,KKM
INTEGER :: ILOC,JLOC,IROW,JCOL,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND)
REAL(EB):: AF,IDX,BIJ,KFACE(1:2,1:2)
TYPE(IBM_REGFACE_TYPE), POINTER, DIMENSION(:) :: REGFACE_H=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()
INTEGER :: IIG,JJG,KKG,II,JJ,KK,IOR,LOCROW,IW
INTEGER :: WC_JD(1:2,1:2)
LOGICAL :: FLG

! Allocate D_MAT_H:
ALLOCATE( D_MAT_H(1:NNZ_ROW_H,1:NUNKH_LOCAL) )
D_MAT_H(:,:)  = 0._EB
D_MAT_HP => D_MAT_H
NM1 = NM_START

! Main Mesh Loop:
MESH_LOOP_1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION_ONLY(NM)) CYCLE MESH_LOOP_1
   CALL POINT_TO_MESH(NM)

   ! X direction bounds:
   ILO_FACE = 0                    ! Low mesh boundary face index.
   IHI_FACE = IBAR                 ! High mesh boundary face index.
   ILO_CELL = ILO_FACE + FCELL     ! First internal cell index. See notes.
   IHI_CELL = IHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Y direction bounds:
   JLO_FACE = 0                    ! Low mesh boundary face index.
   JHI_FACE = JBAR                 ! High mesh boundary face index.
   JLO_CELL = JLO_FACE + FCELL     ! First internal cell index. See notes.
   JHI_CELL = JHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Z direction bounds:
   KLO_FACE = 0                    ! Low mesh boundary face index.
   KHI_FACE = KBAR                 ! High mesh boundary face index.
   KLO_CELL = KLO_FACE + FCELL     ! First internal cell index. See notes.
   KHI_CELL = KHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Regular Faces
   AXIS_LOOP_1 : DO X1AXIS=IAXIS,KAXIS

      NREG = MESHES(NM)%NREGFACE_H(X1AXIS)

      SELECT CASE(X1AXIS)
      CASE(IAXIS)
         REGFACE_H => REGFACE_IAXIS_H
         IIM = FCELL-1; JJM = 0; KKM = 0
         IIP =   FCELL; JJP = 0; KKP = 0
         LOW_FACE=ILO_FACE; HIGH_FACE=IHI_FACE
         X2AXIS=JAXIS; X3AXIS=KAXIS
         DX1 => DXN
         DX2 => DY
         DX3 => DZ
      CASE(JAXIS)
         REGFACE_H => REGFACE_JAXIS_H
         IIM = 0; JJM = FCELL-1; KKM = 0
         IIP = 0; JJP =   FCELL; KKP = 0
         LOW_FACE=JLO_FACE; HIGH_FACE=JHI_FACE
         X2AXIS=KAXIS; X3AXIS=IAXIS
         DX1 => DYN
         DX2 => DZ
         DX3 => DX
      CASE(KAXIS)
         REGFACE_H => REGFACE_KAXIS_H
         IIM = 0; JJM = 0; KKM = FCELL-1
         IIP = 0; JJP = 0; KKP =   FCELL
         LOW_FACE=KLO_FACE; HIGH_FACE=KHI_FACE
         X2AXIS=IAXIS; X3AXIS=JAXIS
         DX1 => DZN
         DX2 => DX
         DX3 => DY
      END SELECT

      IFACE_LOOP_1 : DO IFACE=1,NREG

         I  = REGFACE_H(IFACE)%IJK(IAXIS)
         J  = REGFACE_H(IFACE)%IJK(JAXIS)
         K  = REGFACE_H(IFACE)%IJK(KAXIS)
         I1 = REGFACE_H(IFACE)%IJK(X1AXIS)
         I2 = REGFACE_H(IFACE)%IJK(X2AXIS)
         I3 = REGFACE_H(IFACE)%IJK(X3AXIS)

         ! Unknowns on related cells:
         IND(LOW_IND)  = CCVAR(I+IIM,J+JJM,K+KKM,UNKH)
         IND(HIGH_IND) = CCVAR(I+IIP,J+JJP,K+KKP,UNKH)

         IND_LOC(LOW_IND) = IND(LOW_IND) - UNKH_IND(NM1) ! All row indexes must refer to ind_loc.
         IND_LOC(HIGH_IND)= IND(HIGH_IND)- UNKH_IND(NM1)

         ! Face Area and inv DX1:
         IF (CYLINDRICAL) THEN
            SELECT CASE(X1AXIS)
            CASE(IAXIS); AF  = R(I) *DX3(I3)
            CASE(KAXIS); AF  = RC(I)*DX2(I2)
            END SELECT
         ELSE
            AF  = DX2(I2)*DX3(I3)
         ENDIF
         IDX =    1._EB/DX1(I1)

         ! Now add to Adiff corresponding coeff:
         BIJ = IDX*AF

         ! Cols    ind(1)          ind(2)
         KFACE(1,1)= BIJ; KFACE(1,2)=-BIJ ! Row ind(1)
         KFACE(2,1)=-BIJ; KFACE(2,2)= BIJ ! Row ind(2)

         DO ILOC = LOW_IND,HIGH_IND                   ! Local row number in Kface
            DO JLOC = LOW_IND,HIGH_IND                ! Local col number in Kface, JD
               IROW = IND_LOC(ILOC)                   ! Unknown number.
               JCOL = REGFACE_H(IFACE)%JD(ILOC,JLOC); ! Local position of coef in D_MAT_H
               ! Add coefficient:
               D_MAT_HP(JCOL,IROW) = D_MAT_HP(JCOL,IROW) + KFACE(ILOC,JLOC)
            ENDDO
         ENDDO

      ENDDO IFACE_LOOP_1

      NULLIFY(REGFACE_H)

   ENDDO AXIS_LOOP_1

   ! Next, Wall faces of type INTERPOLATED_BOUNDARY or PERIODIC_BOUNDARY:
   ! Here We have to do something about WALL cells that are also cut-faces, who wins? Make cut-faces take precedence.
   LOCROW = LOW_IND
   WALL_LOOP_1 : DO IW=1,N_EXTERNAL_WALL_CELLS

      WC => WALL(IW)
      EWC=>EXTERNAL_WALL(IW)
      FLG = WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY .OR. WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY
      IF(PRES_ON_WHOLE_DOMAIN) &
      FLG = FLG .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. EWC%NOM > 0)
      IF ( .NOT.FLG ) CYCLE
      ! Here if NOM==0 means it is an OBST laying on an external boundary -> CYCLE
      IF(EWC%NOM < 1) CYCLE

      WC_JD(1,1) = WC%JD11_INDEX
      WC_JD(1,2) = WC%JD12_INDEX
      WC_JD(2,1) = WC%JD21_INDEX
      WC_JD(2,2) = WC%JD22_INDEX

      IIG = WC%ONE_D%IIG; JJG = WC%ONE_D%JJG; KKG = WC%ONE_D%KKG
      II  = WC%ONE_D%II;  JJ  = WC%ONE_D%JJ;  KK  = WC%ONE_D%KK

      IOR = WC%ONE_D%IOR
      ! Check if CC_IBM -> If either IIG,JJG,KKG or II,JJ,KK cell is type IS_CUTCFE or IS_SOLID cycle:
      IF ( .NOT.PRES_ON_WHOLE_DOMAIN .AND. CC_IBM ) THEN
         IF(CCVAR(II ,JJ ,KK ,IBM_CGSC) /= IS_GASPHASE) CYCLE
         IF(CCVAR(IIG,JJG,KKG,IBM_CGSC) /= IS_GASPHASE) CYCLE
      ENDIF

      ! Unknowns on related cells:
      IND(LOW_IND)  = CCVAR(IIG,JJG,KKG,UNKH)  ! internal.
      IND(HIGH_IND) = CCVAR(II,JJ,KK,UNKH)     ! guard-cell.

      IND_LOC(LOW_IND) = IND(LOW_IND) - UNKH_IND(NM1) ! All row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= IND(HIGH_IND)- UNKH_IND(NM1)

      SELECT CASE(IOR)
      CASE( IAXIS)
         AF = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG-1)) * DZ(KKG);         IDX= 1._EB/DXN(IIG-1)
      CASE(-IAXIS)
         AF = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG  )) * DZ(KKG);         IDX= 1._EB/DXN(IIG)
      CASE( JAXIS)
         AF = DX(IIG)*DZ(KKG);         IDX= 1._EB/DYN(JJG-1)
      CASE(-JAXIS)
         AF = DX(IIG)*DZ(KKG);         IDX= 1._EB/DYN(JJG)
      CASE( KAXIS)
         AF = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG);         IDX= 1._EB/DZN(KKG-1)
      CASE(-KAXIS)
         AF = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG);         IDX= 1._EB/DZN(KKG)
      END SELECT

      ! Now add to Adiff corresponding coeff:
      BIJ = IDX*AF

      ! Cols    ind(1)          ind(2)
      KFACE(1,1)= BIJ; KFACE(1,2)=-BIJ ! Row ind(1)
      KFACE(2,1)=-BIJ; KFACE(2,2)= BIJ ! Row ind(2)

      ILOC = LOCROW                             ! Local row number in Kface, only for cell IIG,JJG,KKG.
      DO JLOC = LOW_IND,HIGH_IND                ! Local col number in Kface, JD
         IROW = IND_LOC(ILOC)                   ! Unknown number.
         JCOL = WC_JD(ILOC,JLOC)                ! Local position of coef in D_MAT_H
         ! Add coefficient:
         D_MAT_HP(JCOL,IROW) = D_MAT_HP(JCOL,IROW) + KFACE(ILOC,JLOC)
      ENDDO

   ENDDO WALL_LOOP_1

   ! Contribution to Laplacian matrix from Cut-cells:
   IF ( CC_IBM ) CALL GET_H_MATRIX_CC(NM,NM1,D_MAT_HP)

ENDDO MESH_LOOP_1

RETURN
END SUBROUTINE GET_H_MATRIX


! --------------------------- GET_MATRIXGRAPH_H_WHLDOM ------------------------------

SUBROUTINE GET_MATRIXGRAPH_H_WHLDOM

USE COMPLEX_GEOMETRY, ONLY : GET_CC_MATRIXGRAPH_H, ADD_INPLACE_NNZ_H_WHLDOM
USE MPI

! Local Variables:
INTEGER :: NM
INTEGER :: X1AXIS,IFACE,I,I1,J,K,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND)
INTEGER :: LOCROW,IIND,NII,ILOC
INTEGER :: NREG,IIM,JJM,KKM,IIP,JJP,KKP,LOW_FACE,HIGH_FACE,JLOC,IW,II,JJ,KK,IIG,JJG,KKG
LOGICAL :: INLIST
TYPE(IBM_REGFACE_TYPE), POINTER, DIMENSION(:) :: REGFACE_H=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()
INTEGER :: WC_JD(1:2,1:2)
LOGICAL :: FLG

NUNKH_LOCAL = sum(NUNKH_LOC(1:NMESHES)) ! Filled in GET_MATRIX_INDEXES_H, only nonzeros are for meshes
                                        ! that belong to this process.

! Write number of pressure unknowns to output:
IF (MYID==0) THEN
   WRITE(LU_OUTPUT,'(A)') '   Using GLMAT as pressure solver. List of H unknown numbers per proc:'
ENDIF

! Allocate NNZ_D_MAT_H, JD_MAT_H:
ALLOCATE( NNZ_D_MAT_H(1:NUNKH_LOCAL) )
ALLOCATE( JD_MAT_H(1:NNZ_ROW_H,1:NUNKH_LOCAL) ) ! Contains on first index nonzeros per local row.
NNZ_D_MAT_H(:) = 0
JD_MAT_H(:,:)  = HUGE(I)

! Find NM_START: first mesh that belongs to the processor.
NM_START = LOWER_MESH_INDEX

! First run over all regular and gasphase cut faces and insert-add lists of
! related unknows per unknown:
MESH_LOOP_1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION_ONLY(NM)) CYCLE MESH_LOOP_1
   CALL POINT_TO_MESH(NM)

   ! X direction bounds:
   ILO_FACE = 0                    ! Low mesh boundary face index.
   IHI_FACE = IBAR                 ! High mesh boundary face index.
   ILO_CELL = ILO_FACE + FCELL     ! First internal cell index. See notes.
   IHI_CELL = IHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Y direction bounds:
   JLO_FACE = 0                    ! Low mesh boundary face index.
   JHI_FACE = JBAR                 ! High mesh boundary face index.
   JLO_CELL = JLO_FACE + FCELL     ! First internal cell index. See notes.
   JHI_CELL = JHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Z direction bounds:
   KLO_FACE = 0                    ! Low mesh boundary face index.
   KHI_FACE = KBAR                 ! High mesh boundary face index.
   KLO_CELL = KLO_FACE + FCELL     ! First internal cell index. See notes.
   KHI_CELL = KHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Regular Faces
   AXIS_LOOP_1 : DO X1AXIS=IAXIS,KAXIS

      NREG = MESHES(NM)%NREGFACE_H(X1AXIS)

      SELECT CASE(X1AXIS)
      CASE(IAXIS)
         REGFACE_H => REGFACE_IAXIS_H
         IIM = FCELL-1; JJM = 0; KKM = 0
         IIP =   FCELL; JJP = 0; KKP = 0
         LOW_FACE=ILO_FACE; HIGH_FACE=IHI_FACE
      CASE(JAXIS)
         REGFACE_H => REGFACE_JAXIS_H
         IIM = 0; JJM = FCELL-1; KKM = 0
         IIP = 0; JJP =   FCELL; KKP = 0
         LOW_FACE=JLO_FACE; HIGH_FACE=JHI_FACE
      CASE(KAXIS)
         REGFACE_H => REGFACE_KAXIS_H
         IIM = 0; JJM = 0; KKM = FCELL-1
         IIP = 0; JJP = 0; KKP =   FCELL
         LOW_FACE=KLO_FACE; HIGH_FACE=KHI_FACE
      END SELECT


      IFACE_LOOP_1 : DO IFACE=1,NREG

         I  = REGFACE_H(iface)%IJK(IAXIS)
         J  = REGFACE_H(iface)%IJK(JAXIS)
         K  = REGFACE_H(iface)%IJK(KAXIS)
         I1 = REGFACE_H(iface)%IJK(X1AXIS)

         ! Unknowns on related cells:
         IND(LOW_IND)  = CCVAR(I+IIM,J+JJM,K+KKM,UNKH)
         IND(HIGH_IND) = CCVAR(I+IIP,J+JJP,K+KKP,UNKH)

         IND_LOC(LOW_IND) = IND(LOW_IND) - UNKH_IND(NM_START) ! All row indexes must refer to ind_loc.
         IND_LOC(HIGH_IND)= IND(HIGH_IND)- UNKH_IND(NM_START)

         CALL ADD_INPLACE_NNZ_H_WHLDOM(LOW_IND,HIGH_IND,IND,IND_LOC)

      ENDDO IFACE_LOOP_1

      NULLIFY(REGFACE_H)

   ENDDO AXIS_LOOP_1

   ! Next, Wall faces of type INTERPOLATED_BOUNDARY or PERIODIC_BOUNDARY:
   ! Here We have to do something about WALL cells that are also cut-faces, who wins? Make cut-faces take precedence.
   LOCROW = LOW_IND
   WALL_LOOP_1 : DO IW=1,N_EXTERNAL_WALL_CELLS

     WC => WALL(IW)
     EWC=>EXTERNAL_WALL(IW)
     FLG = WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY .OR. WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY
     IF(PRES_ON_WHOLE_DOMAIN) &
     FLG = FLG .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. EWC%NOM > 0)
     IF ( .NOT.FLG ) CYCLE
     ! Here if NOM==0 means it is an OBST laying on an external boundary -> CYCLE
     IF(EWC%NOM < 1) CYCLE

     IIG = WC%ONE_D%IIG; JJG = WC%ONE_D%JJG; KKG = WC%ONE_D%KKG
     II  = WC%ONE_D%II;  JJ  = WC%ONE_D%JJ;  KK  = WC%ONE_D%KK

     ! Check if CC_IBM -> If either IIG,JJG,KKG or II,JJ,KK cell is type IS_CUTCFE or IS_SOLID cycle:
     IF (  .NOT.PRES_ON_WHOLE_DOMAIN .AND. CC_IBM ) THEN
        IF(CCVAR(II ,JJ ,KK ,IBM_CGSC) /= IS_GASPHASE) CYCLE
        IF(CCVAR(IIG,JJG,KKG,IBM_CGSC) /= IS_GASPHASE) CYCLE
     ENDIF

     ! Unknowns on related cells:
     IND(LOW_IND)  = CCVAR(IIG,JJG,KKG,UNKH)  ! internal.
     IND(HIGH_IND) = CCVAR(II,JJ,KK,UNKH)     ! guard-cell.

     IND_LOC(LOW_IND) = IND(LOW_IND) - UNKH_IND(NM_START) ! All row indexes must refer to ind_loc.
     IND_LOC(HIGH_IND)= IND(HIGH_IND)- UNKH_IND(NM_START)

     ! Same code as inner loop in ADD_INPLACE_NNZ_H_WHLDOM in geom.f90, might be changed by call to it using
     ! LOCROW_1=LOCROW_2=LOCROW
     ! CALL ADD_INPLACE_NNZ_H_WHLDOM(LOCROW,LOCROW,IND,IND_LOC)
     DO IIND=LOW_IND,HIGH_IND

         NII = NNZ_D_MAT_H(IND_LOC(LOCROW))

         ! Check that column index hasn't been already counted:
         INLIST = .FALSE.
         DO ILOC=1,NII
             IF ( IND(IIND) == JD_MAT_H(ILOC,IND_LOC(LOCROW)) ) THEN
                 INLIST = .TRUE.
                 EXIT
             ENDIF
         ENDDO
         IF (INLIST) CYCLE

         ! Now add in place:
         NII = NII + 1
         DO ILOC=1,NII
             IF ( JD_MAT_H(ILOC,IND_LOC(LOCROW)) > IND(IIND) ) EXIT
         ENDDO
         DO JLOC=NII,ILOC+1,-1
             JD_MAT_H(JLOC,IND_LOC(LOCROW)) = JD_MAT_H(JLOC-1,IND_LOC(LOCROW))
         ENDDO
         NNZ_D_MAT_H(IND_LOC(LOCROW))   = NII
         JD_MAT_H(ILOC,IND_LOC(LOCROW)) = IND(IIND)
      ENDDO

   ENDDO WALL_LOOP_1

   ! Finally Add nonzeros corresponding to IBM_RCFACE_H, CUT_FACE
   IF (CC_IBM) CALL GET_CC_MATRIXGRAPH_H(NM,NM_START,.TRUE.)

ENDDO MESH_LOOP_1

! Now add local column location to Faces data structures:
MESH_LOOP_2 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION_ONLY(NM)) CYCLE MESH_LOOP_2
   CALL POINT_TO_MESH(NM)

   ! X direction bounds:
   ILO_FACE = 0                    ! Low mesh boundary face index.
   IHI_FACE = IBAR                 ! High mesh boundary face index.
   ILO_CELL = ILO_FACE + FCELL     ! First internal cell index. See notes.
   IHI_CELL = IHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Y direction bounds:
   JLO_FACE = 0                    ! Low mesh boundary face index.
   JHI_FACE = JBAR                 ! High mesh boundary face index.
   JLO_CELL = JLO_FACE + FCELL     ! First internal cell index. See notes.
   JHI_CELL = JHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Z direction bounds:
   KLO_FACE = 0                    ! Low mesh boundary face index.
   KHI_FACE = KBAR                 ! High mesh boundary face index.
   KLO_CELL = KLO_FACE + FCELL     ! First internal cell index. See notes.
   KHI_CELL = KHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Regular Faces, loop is similar to before:
   AXIS_LOOP_2 : DO X1AXIS=IAXIS,KAXIS

      NREG = MESHES(NM)%NREGFACE_H(X1AXIS)

      SELECT CASE(X1AXIS)
      CASE(IAXIS)
         REGFACE_H => REGFACE_IAXIS_H
         IIM = FCELL-1; JJM = 0; KKM = 0
         IIP =   FCELL; JJP = 0; KKP = 0
         LOW_FACE=ILO_FACE; HIGH_FACE=IHI_FACE
      CASE(JAXIS)
         REGFACE_H => REGFACE_JAXIS_H
         IIM = 0; JJM = FCELL-1; KKM = 0
         IIP = 0; JJP =   FCELL; KKP = 0
         LOW_FACE=JLO_FACE; HIGH_FACE=JHI_FACE
      CASE(KAXIS)
         REGFACE_H => REGFACE_KAXIS_H
         IIM = 0; JJM = 0; KKM = FCELL-1
         IIP = 0; JJP = 0; KKP =   FCELL
         LOW_FACE=KLO_FACE; HIGH_FACE=KHI_FACE
      END SELECT

      IFACE_LOOP_2 : DO IFACE=1,NREG

         I  = REGFACE_H(iface)%IJK(IAXIS)
         J  = REGFACE_H(iface)%IJK(JAXIS)
         K  = REGFACE_H(iface)%IJK(KAXIS)
         I1 = REGFACE_H(iface)%IJK(X1AXIS)

         ! Unknowns on related cells:
         IND(LOW_IND)  = CCVAR(I+IIM,J+JJM,K+KKM,UNKH)
         IND(HIGH_IND) = CCVAR(I+IIP,J+JJP,K+KKP,UNKH)

         IND_LOC(LOW_IND) = IND(LOW_IND) - UNKH_IND(NM_START) ! All row indexes must refer to ind_loc.
         IND_LOC(HIGH_IND)= IND(HIGH_IND)- UNKH_IND(NM_START)

         REGFACE_H(IFACE)%JD(1:2,1:2) = IS_UNDEFINED

         DO LOCROW = LOW_IND,HIGH_IND
            DO IIND=LOW_IND,HIGH_IND
               NII = NNZ_D_MAT_H(IND_LOC(LOCROW))
               DO ILOC=1,NII
                  IF ( IND(IIND) == JD_MAT_H(ILOC,IND_LOC(LOCROW)) ) THEN
                     REGFACE_H(IFACE)%JD(LOCROW,IIND) = ILOC;
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

      ENDDO IFACE_LOOP_2

      NULLIFY(REGFACE_H)

   ENDDO AXIS_LOOP_2

   ! Now Wall faces column locations:
   LOCROW = LOW_IND
   WALL_LOOP_2 : DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

      WC => WALL(IW)

      WC_JD(1:2,1:2) = IS_UNDEFINED

      IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE
      IF(IW <= N_EXTERNAL_WALL_CELLS) THEN
         ! Here if NOM==0 means it is an OBST laying on an external boundary -> CYCLE
         EWC=>EXTERNAL_WALL(IW); IF(EWC%NOM < 1) CYCLE
      ENDIF

      IIG = WC%ONE_D%IIG; JJG = WC%ONE_D%JJG; KKG = WC%ONE_D%KKG
      II  = WC%ONE_D%II;  JJ  = WC%ONE_D%JJ;  KK  = WC%ONE_D%KK

      ! Check if CC_IBM -> If either IIG,JJG,KKG or II,JJ,KK cell is type IS_CUTCFE or IS_SOLID cycle:
      IF (  .NOT.PRES_ON_WHOLE_DOMAIN .AND. CC_IBM ) THEN
         IF(CCVAR(II ,JJ ,KK ,IBM_CGSC) /= IS_GASPHASE) CYCLE
         IF(CCVAR(IIG,JJG,KKG,IBM_CGSC) /= IS_GASPHASE) CYCLE
      ENDIF

      ! Unknowns on related cells:
      IND(LOW_IND)  = CCVAR(IIG,JJG,KKG,UNKH)  ! internal.
      IND(HIGH_IND) = CCVAR(II,JJ,KK,UNKH)     ! guard-cell.

      IND_LOC(LOW_IND) = IND(LOW_IND) - UNKH_IND(NM_START) ! All row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= IND(HIGH_IND)- UNKH_IND(NM_START)

      DO IIND=LOW_IND,HIGH_IND
         NII = NNZ_D_MAT_H(IND_LOC(LOCROW))
         DO ILOC=1,NII
            IF ( IND(IIND) == JD_MAT_H(ILOC,IND_LOC(LOCROW)) ) THEN
               WC_JD(LOCROW,IIND) = ILOC
               EXIT
            ENDIF
         ENDDO
      ENDDO
      WC%JD11_INDEX = WC_JD(1,1)
      WC%JD12_INDEX = WC_JD(1,2)
      WC%JD21_INDEX = WC_JD(2,1)
      WC%JD22_INDEX = WC_JD(2,2)

   ENDDO WALL_LOOP_2

   ! Finally cut-face:
   IF (CC_IBM) CALL GET_CC_MATRIXGRAPH_H(NM,NM_START,.FALSE.)

ENDDO MESH_LOOP_2


RETURN
END SUBROUTINE GET_MATRIXGRAPH_H_WHLDOM

! ---------------------------- GET_H_REGFACES -----------------------------------

SUBROUTINE GET_H_REGFACES

USE COMPLEX_GEOMETRY, ONLY : GET_RCFACES_H

! Local Variables:
INTEGER :: NM
INTEGER :: ILO,IHI,JLO,JHI,KLO,KHI
INTEGER :: I,J,K,II,IREG,X1AXIS
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IJKBUFFER
INTEGER :: IW, IIG, JJG, KKG, IOR, N_INTERNAL_WALL_CELLS_AUX
LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: LOG_INTWC
TYPE(WALL_TYPE), POINTER :: WC=>NULL()

! Mesh loop:
MAIN_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION_ONLY(NM)) CYCLE MAIN_MESH_LOOP
   CALL POINT_TO_MESH(NM)

   ! Mesh sizes:
   NXB=IBAR; NYB=JBAR; NZB=KBAR

   ! X direction bounds:
   ILO_FACE = 0                    ! Low mesh boundary face index.
   IHI_FACE = IBAR                 ! High mesh boundary face index.
   ILO_CELL = ILO_FACE + FCELL     ! First internal cell index. See notes.
   IHI_CELL = IHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Y direction bounds:
   JLO_FACE = 0                    ! Low mesh boundary face index.
   JHI_FACE = JBAR                 ! High mesh boundary face index.
   JLO_CELL = JLO_FACE + FCELL     ! First internal cell index. See notes.
   JHI_CELL = JHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Z direction bounds:
   KLO_FACE = 0                    ! Low mesh boundary face index.
   KHI_FACE = KBAR                 ! High mesh boundary face index.
   KLO_CELL = KLO_FACE + FCELL     ! First internal cell index. See notes.
   KHI_CELL = KHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Set starting number of regular faces for NM to zero:
   MESHES(NM)%NREGFACE_H(IAXIS:KAXIS) = 0

   ! 1. Regular GASPHASE faces connected to Gasphase cells:
   ALLOCATE(IJKBUFFER(IAXIS:KAXIS,1:(NXB+1)*(NYB+1)*(NZB+1)))
   ! Check internal SOLID_BOUNDARY faces:
   ALLOCATE(LOG_INTWC(ILO_FACE:IHI_FACE,JLO_FACE:JHI_FACE,KLO_FACE:KHI_FACE,IAXIS:KAXIS)); LOG_INTWC(:,:,:,:) = .FALSE.

   N_INTERNAL_WALL_CELLS_AUX = 0
   IF(.NOT. PRES_ON_WHOLE_DOMAIN) N_INTERNAL_WALL_CELLS_AUX = N_INTERNAL_WALL_CELLS
   WALL_LOOP_1 : DO IW=N_EXTERNAL_WALL_CELLS+1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS_AUX
      WC => WALL(IW)
      IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE
      IIG = WC%ONE_D%IIG; JJG = WC%ONE_D%JJG; KKG = WC%ONE_D%KKG

      IOR = WC%ONE_D%IOR
      SELECT CASE(IOR)
      CASE( IAXIS)
         LOG_INTWC(IIG-1,JJG  ,KKG  ,IAXIS) = .TRUE.
      CASE(-IAXIS)
         LOG_INTWC(IIG  ,JJG  ,KKG  ,IAXIS) = .TRUE.
      CASE( JAXIS)
         LOG_INTWC(IIG  ,JJG-1,KKG  ,JAXIS) = .TRUE.
      CASE(-JAXIS)
         LOG_INTWC(IIG  ,JJG  ,KKG  ,JAXIS) = .TRUE.
      CASE( KAXIS)
         LOG_INTWC(IIG  ,JJG  ,KKG-1,KAXIS) = .TRUE.
      CASE(-KAXIS)
         LOG_INTWC(IIG  ,JJG  ,KKG  ,KAXIS) = .TRUE.
      END SELECT
   ENDDO WALL_LOOP_1

   ! axis = IAXIS:
   X1AXIS = IAXIS
   ILO = ILO_FACE+1; IHI = IHI_FACE-1 ! ILO_FACE and IHI_FACE faces are defined in WALL
   JLO = JLO_CELL;   JHI = JHI_CELL
   KLO = KLO_CELL;   KHI = KHI_CELL
   ! Loop on Cartesian cells, define cut cells and solid cells CGSC:
   ! First count for allocation:
   IREG = 0
   DO K=KLO,KHI
      DO J=JLO,JHI
         DO I=ILO,IHI
            IF (LOG_INTWC(I,J,K,X1AXIS))        CYCLE
            IF (CCVAR(I+FCELL-1,J,K,UNKH) <= 0) CYCLE
            IF (CCVAR(I+FCELL  ,J,K,UNKH) <= 0) CYCLE
            IREG = IREG + 1
            IJKBUFFER(IAXIS:KAXIS,IREG) = (/ I, J, K /)
         ENDDO
      ENDDO
   ENDDO
   MESHES(NM)%NREGFACE_H(X1AXIS) = IREG
   NULLIFY(REGFACE_IAXIS_H) ! Nullify pointer to mesh variable M%REGFACE_IAXIS_H, we are about to allocate it.
   IF(ALLOCATED(MESHES(NM)%REGFACE_IAXIS_H)) DEALLOCATE(MESHES(NM)%REGFACE_IAXIS_H)
   ALLOCATE(MESHES(NM)%REGFACE_IAXIS_H(IREG))
   DO II=1,IREG
      MESHES(NM)%REGFACE_IAXIS_H(II)%IJK(IAXIS:KAXIS) = IJKBUFFER(IAXIS:KAXIS,II)
   ENDDO

   ! axis = JAXIS:
   X1AXIS = JAXIS
   ILO = ILO_CELL;   IHI = IHI_CELL
   JLO = JLO_FACE+1; JHI = JHI_FACE-1 ! JLO_FACE and JHI_FACE faces are defined in WALL
   KLO = KLO_CELL;   KHI = KHI_CELL
   ! Loop on Cartesian cells, define cut cells and solid cells CGSC:
   ! First count for allocation:
   IREG = 0
   DO K=KLO,KHI
      DO J=JLO,JHI
         DO I=ILO,IHI
            IF (LOG_INTWC(I,J,K,X1AXIS))        CYCLE
            IF (CCVAR(I,J+FCELL-1,K,UNKH) <= 0) CYCLE
            IF (CCVAR(I,J+FCELL  ,K,UNKH) <= 0) CYCLE
            IREG = IREG + 1
            IJKBUFFER(IAXIS:KAXIS,IREG) = (/ I, J, K /)
         ENDDO
      ENDDO
   ENDDO
   MESHES(NM)%NREGFACE_H(X1AXIS) = IREG
   NULLIFY(REGFACE_JAXIS_H) ! Nullify pointer to mesh variable M%REGFACE_JAXIS_H, we are about to allocate it.
   IF(ALLOCATED(MESHES(NM)%REGFACE_JAXIS_H)) DEALLOCATE(MESHES(NM)%REGFACE_JAXIS_H)
   ALLOCATE(MESHES(NM)%REGFACE_JAXIS_H(IREG))
   DO II=1,IREG
      MESHES(NM)%REGFACE_JAXIS_H(II)%IJK(IAXIS:KAXIS) = IJKBUFFER(IAXIS:KAXIS,II)
   ENDDO

   ! axis = KAXIS:
   X1AXIS = KAXIS
   ILO = ILO_CELL;   IHI = IHI_CELL
   JLO = JLO_CELL;   JHI = JHI_CELL
   KLO = KLO_FACE+1; KHI = KHI_FACE-1 ! KLO_FACE and KHI_FACE faces are defined in WALL
   ! Loop on Cartesian cells, define cut cells and solid cells CGSC:
   ! First count for allocation:
   IREG = 0
   DO K=KLO,KHI
      DO J=JLO,JHI
         DO I=ILO,IHI
            IF (LOG_INTWC(I,J,K,X1AXIS))        CYCLE
            IF (CCVAR(I,J,K+FCELL-1,UNKH) <= 0) CYCLE
            IF (CCVAR(I,J,K+FCELL  ,UNKH) <= 0) CYCLE
            IREG = IREG + 1
            IJKBUFFER(IAXIS:KAXIS,IREG) = (/ I, J, K /)
         ENDDO
      ENDDO
   ENDDO
   MESHES(NM)%NREGFACE_H(X1AXIS) = IREG
   NULLIFY(REGFACE_KAXIS_H) ! Nullify pointer to mesh variable M%REGFACE_KAXIS_H, we are about to allocate it.
   IF(ALLOCATED(MESHES(NM)%REGFACE_KAXIS_H)) DEALLOCATE(MESHES(NM)%REGFACE_KAXIS_H)
   ALLOCATE(MESHES(NM)%REGFACE_KAXIS_H(IREG))
   DO II=1,IREG
      MESHES(NM)%REGFACE_KAXIS_H(II)%IJK(IAXIS:KAXIS) = IJKBUFFER(IAXIS:KAXIS,II)
   ENDDO

   ! 2. Lists of Regular Gasphase faces, connected to one regular gasphase and one cut-cell:
   IF (CC_IBM) CALL GET_RCFACES_H(NM)

   DEALLOCATE(IJKBUFFER,LOG_INTWC)

ENDDO MAIN_MESH_LOOP


RETURN
END SUBROUTINE GET_H_REGFACES

! ------------------------- GET_MATRIX_INDEXES_H --------------------------------

SUBROUTINE GET_MATRIX_INDEXES_H

USE COMPLEX_GEOMETRY, ONLY : NUMBER_UNKH_CUTCELLS
USE MPI

! Local Variables:
INTEGER :: NM
INTEGER :: ILO,IHI,JLO,JHI,KLO,KHI
INTEGER :: I,J,K,IERR
INTEGER, PARAMETER :: IMPADD = 1
INTEGER, PARAMETER :: SHFTM(1:3,1:6) = RESHAPE((/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))

! Define local number of cut-cell:
IF (ALLOCATED(NUNKH_LOC)) DEALLOCATE(NUNKH_LOC)
ALLOCATE(NUNKH_LOC(1:NMESHES)); NUNKH_LOC = 0

! Cell numbers for Poisson equation:
MAIN_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION_ONLY(NM)) CYCLE MAIN_MESH_LOOP
   CALL POINT_TO_MESH(NM)

   ! Mesh sizes:
   NXB=IBAR; NYB=JBAR; NZB=KBAR

   ! X direction bounds:
   ILO_FACE = 0                    ! Low mesh boundary face index.
   IHI_FACE = IBAR                 ! High mesh boundary face index.
   ILO_CELL = ILO_FACE + FCELL     ! First internal cell index. See notes.
   IHI_CELL = IHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Y direction bounds:
   JLO_FACE = 0                    ! Low mesh boundary face index.
   JHI_FACE = JBAR                 ! High mesh boundary face index.
   JLO_CELL = JLO_FACE + FCELL     ! First internal cell index. See notes.
   JHI_CELL = JHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Z direction bounds:
   KLO_FACE = 0                    ! Low mesh boundary face index.
   KHI_FACE = KBAR                 ! High mesh boundary face index.
   KLO_CELL = KLO_FACE + FCELL     ! First internal cell index. See notes.
   KHI_CELL = KHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Initialize unknown numbers for H and Z:
   ! We assume SET_CUTCELLS_3D has been called and CCVAR has been allocated:

   CCVAR(:,:,:,UNKH) = IS_UNDEFINED

   ! For Pressure Number regular GASPHASE cells:
   ILO = ILO_CELL; IHI = IHI_CELL
   JLO = JLO_CELL; JHI = JHI_CELL
   KLO = KLO_CELL; KHI = KHI_CELL

   IF_PRES_ON_WHOLE_DOMAIN : IF(PRES_ON_WHOLE_DOMAIN) THEN ! Classic IBM.
      ! Loop on Cartesian cells, define cut cells and solid cells CGSC:
      DO K=KLO,KHI
         DO J=JLO,JHI
            DO I=ILO,IHI
               NUNKH_LOC(NM) = NUNKH_LOC(NM) + 1
               CCVAR(I,J,K,UNKH) = NUNKH_LOC(NM)
            ENDDO
         ENDDO
      ENDDO

   ELSE

      ! Loop on Cartesian cells, define cut cells and solid cells CGSC:
      DO K=KLO,KHI
         DO J=JLO,JHI
            DO I=ILO,IHI
               IF (  CCVAR(I,J,K,CGSC) == IS_GASPHASE ) THEN
                  NUNKH_LOC(NM) = NUNKH_LOC(NM) + 1
                  CCVAR(I,J,K,UNKH) = NUNKH_LOC(NM)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ! Number MESH local unknowns on cut-cells (fully unstructured) or their under underlaying Cartesian cells
      ! (Cartesian unstructured):
      IF(CC_IBM) CALL NUMBER_UNKH_CUTCELLS(.TRUE.,NM,NUNKH_LOC)

   ENDIF IF_PRES_ON_WHOLE_DOMAIN

ENDDO MAIN_MESH_LOOP

! Define total number of unknowns and global unknow index start per MESH:
IF (ALLOCATED(NUNKH_TOT)) DEALLOCATE(NUNKH_TOT)
ALLOCATE(NUNKH_TOT(1:NMESHES)); NUNKH_TOT = 0
IF (N_MPI_PROCESSES > 1) THEN
   CALL MPI_ALLREDUCE(NUNKH_LOC, NUNKH_TOT, NMESHES, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
ELSE
   NUNKH_TOT = NUNKH_LOC
ENDIF
! Define global start indexes for each mesh:
IF (ALLOCATED(UNKH_IND)) DEALLOCATE(UNKH_IND)
ALLOCATE(UNKH_IND(1:NMESHES)); UNKH_IND = 0
DO NM=2,NMESHES
   UNKH_IND(NM) = UNKH_IND(NM-1) + NUNKH_TOT(NM-1)
ENDDO

! Add initial index UNKX_ind to mesh blocks (regular + cut-cells):
MAIN_MESH_LOOP2 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION_ONLY(NM)) CYCLE MAIN_MESH_LOOP2
   CALL POINT_TO_MESH(NM)

   ! Mesh sizes:
   NXB=IBAR; NYB=JBAR; NZB=KBAR

   ! X direction bounds:
   ILO_FACE = 0                    ! Low mesh boundary face index.
   IHI_FACE = IBAR                 ! High mesh boundary face index.
   ILO_CELL = ILO_FACE + FCELL     ! First internal cell index. See notes.
   IHI_CELL = IHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Y direction bounds:
   JLO_FACE = 0                    ! Low mesh boundary face index.
   JHI_FACE = JBAR                 ! High mesh boundary face index.
   JLO_CELL = JLO_FACE + FCELL     ! First internal cell index. See notes.
   JHI_CELL = JHI_FACE + FCELL - 1 ! Last internal cell index.

   ! Z direction bounds:
   KLO_FACE = 0                    ! Low mesh boundary face index.
   KHI_FACE = KBAR                 ! High mesh boundary face index.
   KLO_CELL = KLO_FACE + FCELL     ! First internal cell index. See notes.
   KHI_CELL = KHI_FACE + FCELL - 1 ! Last internal cell index.

   ! For Pressure Global Number regular GASPHASE cells:
   ILO = ILO_CELL; IHI = IHI_CELL
   JLO = JLO_CELL; JHI = JHI_CELL
   KLO = KLO_CELL; KHI = KHI_CELL

   IF_PRES_ON_WHOLE_DOMAIN2 : IF(PRES_ON_WHOLE_DOMAIN) THEN ! Classic IBM.

      ! Loop on all Cartesian cells:
      DO K=KLO,KHI
         DO J=JLO,JHI
            DO I=ILO,IHI
               CCVAR(I,J,K,UNKH) = CCVAR(I,J,K,UNKH) + UNKH_IND(NM)
            ENDDO
         ENDDO
      ENDDO

   ELSE

      ! Loop on Cartesian cells, GASPHASE:
      DO K=KLO,KHI
         DO J=JLO,JHI
            DO I=ILO,IHI
               IF (  CCVAR(I,J,K,CGSC) == IS_GASPHASE ) THEN
                  CCVAR(I,J,K,UNKH) = CCVAR(I,J,K,UNKH) + UNKH_IND(NM)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ! Number MESH global unknowns on cut-cells (fully unstructured) or their under underlaying Cartesian cells
      ! (Cartesian unstructured):
      IF(CC_IBM) CALL NUMBER_UNKH_CUTCELLS(.FALSE.,NM,NUNKH_LOC)

   ENDIF IF_PRES_ON_WHOLE_DOMAIN2


ENDDO MAIN_MESH_LOOP2


RETURN
END SUBROUTINE GET_MATRIX_INDEXES_H


! ------------------------- SET_CCVAR_CGSC_H ---------------------------------

SUBROUTINE SET_CCVAR_CGSC_H


! Local Variables
INTEGER :: NM
INTEGER :: ISTR, IEND, JSTR, JEND, KSTR, KEND
INTEGER :: I,J,K

! First Loop, in case CC_IBM hasn't been defined +> the CCVAR ARRAY hasn't been allocated:
IF (.NOT.CC_IBM) THEN
   FIRST_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (EVACUATION_ONLY(NM)) CYCLE FIRST_MESH_LOOP

      ! Mesh sizes:
      NXB=MESHES(NM)%IBAR
      NYB=MESHES(NM)%JBAR
      NZB=MESHES(NM)%KBAR

      ! X direction bounds:
      ILO_FACE = 0                    ! Low mesh boundary face index.
      IHI_FACE = MESHES(NM)%IBAR      ! High mesh boundary face index.
      ILO_CELL = ILO_FACE + FCELL     ! First internal cell index. See notes.
      IHI_CELL = IHI_FACE + FCELL - 1 ! Last internal cell index.
      ISTR     = ILO_FACE - NGUARD    ! Allocation start x arrays.
      IEND     = IHI_FACE + NGUARD    ! Allocation end x arrays.

      ! Y direction bounds:
      JLO_FACE = 0                    ! Low mesh boundary face index.
      JHI_FACE = MESHES(NM)%JBAR      ! High mesh boundary face index.
      JLO_CELL = JLO_FACE + FCELL     ! First internal cell index. See notes.
      JHI_CELL = JHI_FACE + FCELL - 1 ! Last internal cell index.
      JSTR     = JLO_FACE - NGUARD    ! Allocation start y arrays.
      JEND     = JHI_FACE + NGUARD    ! Allocation end y arrays.

      ! Z direction bounds:
      KLO_FACE = 0                    ! Low mesh boundary face index.
      KHI_FACE = MESHES(NM)%KBAR      ! High mesh boundary face index.
      KLO_CELL = KLO_FACE + FCELL     ! First internal cell index. See notes.
      KHI_CELL = KHI_FACE + FCELL - 1 ! Last internal cell index.
      KSTR     = KLO_FACE - NGUARD    ! Allocation start z arrays.
      KEND     = KHI_FACE + NGUARD    ! Allocation end z arrays.

      ! Cartesian cells:
      IF (.NOT. ALLOCATED(MESHES(NM)%CCVAR)) &
      ALLOCATE(MESHES(NM)%CCVAR(ISTR:IEND,JSTR:JEND,KSTR:KEND,NCVARS))
      MESHES(NM)%CCVAR = 0
      MESHES(NM)%CCVAR(:,:,:,CGSC) = IS_UNDEFINED
      MESHES(NM)%CCVAR(ILO_CELL:IHI_CELL,JLO_CELL:JHI_CELL,KLO_CELL:KHI_CELL,CGSC) = IS_GASPHASE ! Set all internal
                                                                                     ! Block cells to GASPHASE.

   ENDDO FIRST_MESH_LOOP
ENDIF

! At this point, if CC_IBM is being used we have SOLID and CUTCELLS from GEOM objects.
! Now add OBST SOLID cells to CCVAR CGSC array:
SECND_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   IF (EVACUATION_ONLY(NM)) CYCLE SECND_MESH_LOOP
   CALL POINT_TO_MESH(NM)

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CCVAR(I,J,K,CGSC) = IS_SOLID
         ENDDO
      ENDDO
   ENDDO

ENDDO SECND_MESH_LOOP

RETURN
END SUBROUTINE SET_CCVAR_CGSC_H


! --------------------- PRESSURE_SOLVER_CHECK_RESIDUALS_U ----------------------

SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS_U(NM)

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE GLOBAL_CONSTANTS

INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,RHOP,P,RESIDUAL
INTEGER :: I,J,K
REAL(EB) :: LHSS,RHSS,IMFCT,JMFCT,KMFCT,IPFCT,JPFCT,KPFCT

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   HP => H
   RHOP => RHO
ELSE
   HP => HS
   RHOP => RHOS
ENDIF

! Optional check of the accuracy of the separable pressure solution, del^2 H = -del dot F - dD/dt

IF (CHECK_POISSON) THEN
   RESIDUAL => WORK8(1:IBAR,1:JBAR,1:KBAR); RESIDUAL = 0._EB
   !$OMP PARALLEL DO PRIVATE(I,J,K,RHSS,LHSS,IMFCT,JMFCT,KMFCT,IPFCT,JPFCT,KPFCT) SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF(SOLID(CELL_INDEX(I,J,K))) CYCLE
            IMFCT = 1._EB; JMFCT = 1._EB; KMFCT = 1._EB; IPFCT = 1._EB; JPFCT = 1._EB; KPFCT = 1._EB
            ! If surrounding wall_cell is type SOLID_BOUNDARY set FCT gradient factor to zero:
            IF (WALL(WALL_INDEX(CELL_INDEX(I,J,K),-1))%BOUNDARY_TYPE==SOLID_BOUNDARY) IMFCT = 0._EB
            IF (WALL(WALL_INDEX(CELL_INDEX(I,J,K), 1))%BOUNDARY_TYPE==SOLID_BOUNDARY) IPFCT = 0._EB
            IF (WALL(WALL_INDEX(CELL_INDEX(I,J,K),-2))%BOUNDARY_TYPE==SOLID_BOUNDARY) JMFCT = 0._EB
            IF (WALL(WALL_INDEX(CELL_INDEX(I,J,K), 2))%BOUNDARY_TYPE==SOLID_BOUNDARY) JPFCT = 0._EB
            IF (WALL(WALL_INDEX(CELL_INDEX(I,J,K),-3))%BOUNDARY_TYPE==SOLID_BOUNDARY) KMFCT = 0._EB
            IF (WALL(WALL_INDEX(CELL_INDEX(I,J,K), 3))%BOUNDARY_TYPE==SOLID_BOUNDARY) KPFCT = 0._EB
            IF (CC_IBM) THEN
               IF(CCVAR(I,J,K,IBM_CGSC)==IS_SOLID) CYCLE
               IF(FCVAR(I-1,J  ,K  ,IBM_FGSC,IAXIS)==IS_SOLID) IMFCT = 0._EB
               IF(FCVAR(I  ,J  ,K  ,IBM_FGSC,IAXIS)==IS_SOLID) IPFCT = 0._EB
               IF(FCVAR(I  ,J-1,K  ,IBM_FGSC,JAXIS)==IS_SOLID) JMFCT = 0._EB
               IF(FCVAR(I  ,J  ,K  ,IBM_FGSC,JAXIS)==IS_SOLID) JPFCT = 0._EB
               IF(FCVAR(I  ,J  ,K-1,IBM_FGSC,KAXIS)==IS_SOLID) KMFCT = 0._EB
               IF(FCVAR(I  ,J  ,K  ,IBM_FGSC,KAXIS)==IS_SOLID) KPFCT = 0._EB
            ENDIF
            RHSS = ( R(I-1)*FVX(I-1,J,K) - R(I)*FVX(I,J,K) )*RDX(I)*RRN(I) &
                 + (        FVY(I,J-1,K) -      FVY(I,J,K) )*RDY(J)        &
                 + (        FVZ(I,J,K-1) -      FVZ(I,J,K) )*RDZ(K)        &
                 - DDDT(I,J,K)
            LHSS = ((HP(I+1,J,K)-HP(I,J,K))*RDXN(I)*R(I)*IPFCT - (HP(I,J,K)-HP(I-1,J,K))*RDXN(I-1)*R(I-1)*IMFCT )*RDX(I)*RRN(I) &
                 + ((HP(I,J+1,K)-HP(I,J,K))*RDYN(J)*JPFCT      - (HP(I,J,K)-HP(I,J-1,K))*RDYN(J-1)*JMFCT        )*RDY(J)        &
                 + ((HP(I,J,K+1)-HP(I,J,K))*RDZN(K)*KPFCT      - (HP(I,J,K)-HP(I,J,K-1))*RDZN(K-1)*KMFCT        )*RDZ(K)
            RESIDUAL(I,J,K) = ABS(RHSS-LHSS)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   POIS_ERR = MAXVAL(RESIDUAL)
ENDIF

! Mandatory check of how well the computed pressure satisfies the inseparable Poisson equation:
! LHSS = del dot (1/rho) del p + del K = -del dot F - dD/dt = RHSS

IF (ITERATE_BAROCLINIC_TERM) THEN
   P => WORK7
   P = RHOP*(HP-KRES)
   RESIDUAL => WORK8(1:IBAR,1:JBAR,1:KBAR); RESIDUAL = 0._EB
   !$OMP PARALLEL PRIVATE(I,J,K,RHSS,LHSS,IMFCT,JMFCT,KMFCT,IPFCT,JPFCT,KPFCT,NM)
   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF(SOLID(CELL_INDEX(I,J,K))) CYCLE
            IMFCT = 1._EB; JMFCT = 1._EB; KMFCT = 1._EB; IPFCT = 1._EB; JPFCT = 1._EB; KPFCT = 1._EB
            ! If surrounding wall_cell is type SOLID_BOUNDARY set FCT gradient factor to zero:
            IF (WALL(WALL_INDEX(CELL_INDEX(I,J,K),-1))%BOUNDARY_TYPE==SOLID_BOUNDARY) IMFCT = 0._EB
            IF (WALL(WALL_INDEX(CELL_INDEX(I,J,K), 1))%BOUNDARY_TYPE==SOLID_BOUNDARY) IPFCT = 0._EB
            IF (WALL(WALL_INDEX(CELL_INDEX(I,J,K),-2))%BOUNDARY_TYPE==SOLID_BOUNDARY) JMFCT = 0._EB
            IF (WALL(WALL_INDEX(CELL_INDEX(I,J,K), 2))%BOUNDARY_TYPE==SOLID_BOUNDARY) JPFCT = 0._EB
            IF (WALL(WALL_INDEX(CELL_INDEX(I,J,K),-3))%BOUNDARY_TYPE==SOLID_BOUNDARY) KMFCT = 0._EB
            IF (WALL(WALL_INDEX(CELL_INDEX(I,J,K), 3))%BOUNDARY_TYPE==SOLID_BOUNDARY) KPFCT = 0._EB
            IF (CC_IBM) THEN
               IF(CCVAR(I,J,K,IBM_CGSC)==IS_SOLID) CYCLE
               IF(FCVAR(I-1,J  ,K  ,IBM_FGSC,IAXIS)==IS_SOLID) IMFCT = 0._EB
               IF(FCVAR(I  ,J  ,K  ,IBM_FGSC,IAXIS)==IS_SOLID) IPFCT = 0._EB
               IF(FCVAR(I  ,J-1,K  ,IBM_FGSC,JAXIS)==IS_SOLID) JMFCT = 0._EB
               IF(FCVAR(I  ,J  ,K  ,IBM_FGSC,JAXIS)==IS_SOLID) JPFCT = 0._EB
               IF(FCVAR(I  ,J  ,K-1,IBM_FGSC,KAXIS)==IS_SOLID) KMFCT = 0._EB
               IF(FCVAR(I  ,J  ,K  ,IBM_FGSC,KAXIS)==IS_SOLID) KPFCT = 0._EB
            ENDIF
            RHSS = ( R(I-1)*(FVX(I-1,J,K)-FVX_B(I-1,J,K)) - R(I)*(FVX(I,J,K)-FVX_B(I,J,K)) )*RDX(I)*RRN(I) &
                 + (        (FVY(I,J-1,K)-FVY_B(I,J-1,K)) -      (FVY(I,J,K)-FVY_B(I,J,K)) )*RDY(J)        &
                 + (        (FVZ(I,J,K-1)-FVZ_B(I,J,K-1)) -      (FVZ(I,J,K)-FVZ_B(I,J,K)) )*RDZ(K)        &
                 - DDDT(I,J,K)
            LHSS = ((P(I+1,J,K)-P(I,J,K))*RDXN(I)*R(I)    *2._EB/(RHOP(I+1,J,K)+RHOP(I,J,K))*IPFCT - &
                    (P(I,J,K)-P(I-1,J,K))*RDXN(I-1)*R(I-1)*2._EB/(RHOP(I-1,J,K)+RHOP(I,J,K))*IMFCT)*RDX(I)*RRN(I) &
                 + ((P(I,J+1,K)-P(I,J,K))*RDYN(J)         *2._EB/(RHOP(I,J+1,K)+RHOP(I,J,K))*JPFCT - &
                    (P(I,J,K)-P(I,J-1,K))*RDYN(J-1)       *2._EB/(RHOP(I,J-1,K)+RHOP(I,J,K))*JMFCT)*RDY(J)        &
                 + ((P(I,J,K+1)-P(I,J,K))*RDZN(K)         *2._EB/(RHOP(I,J,K+1)+RHOP(I,J,K))*KPFCT - &
                    (P(I,J,K)-P(I,J,K-1))*RDZN(K-1)       *2._EB/(RHOP(I,J,K-1)+RHOP(I,J,K))*KMFCT)*RDZ(K)        &
            + ((KRES(I+1,J,K)-KRES(I,J,K))*RDXN(I)*R(I)*IPFCT - (KRES(I,J,K)-KRES(I-1,J,K))*RDXN(I-1)*R(I-1)*IMFCT )*RDX(I)*RRN(I) &
            + ((KRES(I,J+1,K)-KRES(I,J,K))*RDYN(J)*JPFCT      - (KRES(I,J,K)-KRES(I,J-1,K))*RDYN(J-1)*JMFCT        )*RDY(J)        &
            + ((KRES(I,J,K+1)-KRES(I,J,K))*RDZN(K)*KPFCT      - (KRES(I,J,K)-KRES(I,J,K-1))*RDZN(K-1)*KMFCT        )*RDZ(K)
            RESIDUAL(I,J,K) = ABS(RHSS-LHSS)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
   !$OMP END PARALLEL
   PRESSURE_ERROR_MAX(NM) = MAXVAL(RESIDUAL)
   PRESSURE_ERROR_MAX_LOC(:,NM) = MAXLOC(RESIDUAL)
ENDIF

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
END SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS_U


! --------------------------- FINISH_GLMAT_SOLVER_H --------------------------------

SUBROUTINE FINISH_GLMAT_SOLVER_H

USE MPI

! Local variables:
INTEGER :: MAXFCT, MNUM, MTYPE, PHASE, NRHS, ERROR, MSGLVL
#ifdef WITH_MKL
INTEGER :: PERM(1)
#endif

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

! Solve:
NRHS   =  1
MAXFCT =  1
MNUM   =  1
ERROR  =  0 ! initialize error flag
MSGLVL =  0 ! print statistical information
IF ( H_MATRIX_INDEFINITE ) THEN
   MTYPE  = -2 ! symmetric indefinite
ELSE ! positive definite
   MTYPE  =  2
ENDIF

! Finalize Pardiso:
PHASE = -1
! PARDISO:
! CALL PARDISO(PT_H, MAXFCT, MNUM, MTYPE, PHASE, NUNKH_TOTAL, &
!      A_H, IA_H, JA_H, PERM, NRHS, IPARM, MSGLVL, F_H, X_H, ERROR)
#ifdef WITH_MKL
CALL CLUSTER_SPARSE_SOLVER(PT_H, MAXFCT, MNUM, MTYPE, PHASE, NUNKH_TOTAL, &
     A_H, IA_H, JA_H, PERM, NRHS, IPARM, MSGLVL, F_H, X_H, MPI_COMM_WORLD, ERROR)
#endif /* WITH_MKL */

RETURN
END SUBROUTINE FINISH_GLMAT_SOLVER_H


END MODULE GLOBMAT_SOLVER
