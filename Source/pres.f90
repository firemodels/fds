MODULE PRES

! Find the perturbation pressure by solving Poisson's Equation

USE PRECISION_PARAMETERS
USE MESH_VARIABLES

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC PRESSURE_SOLVER_COMPUTE_RHS,PRESSURE_SOLVER_FFT,TUNNEL_POISSON_SOLVER,PRESSURE_SOLVER_CHECK_RESIDUALS, &
       COMPUTE_VELOCITY_ERROR

CONTAINS

SUBROUTINE PRESSURE_SOLVER_COMPUTE_RHS(T,DT,NM)

USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE GLOBAL_CONSTANTS

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,HP,RHOP
INTEGER :: I,J,K,IW,IOR,NOM
REAL(EB) :: TRM1,TRM2,TRM3,TRM4,TNOW, &
            TSI,TIME_RAMP_FACTOR,DX_OTHER,DY_OTHER,DZ_OTHER,P_EXTERNAL,VEL_EDDY,H0
TYPE (VENTS_TYPE), POINTER :: VT
TYPE (WALL_TYPE), POINTER :: WC
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE (BOUNDARY_PROP1_TYPE), POINTER :: B1
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

!$OMP PARALLEL

! Apply pressure boundary conditions at external cells.
! If Neumann, BXS, BXF, etc., contain dH/dx(x=XS), dH/dx(x=XF), etc.
! If Dirichlet, BXS, BXF, etc., contain H(x=XS), H(x=XF), etc.
! LBC, MBC and NBC are codes used be Poisson solver to denote type
! of boundary condition at x, y and z boundaries. See Crayfishpak
! manual for details.

!$OMP DO PRIVATE(IW,WC,EWC,BC,B1,I,J,K,IOR,NOM,DX_OTHER,DY_OTHER,DZ_OTHER,VT,TSI) &
!$OMP&   PRIVATE(TIME_RAMP_FACTOR,P_EXTERNAL,VEL_EDDY,H0)
WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

   WC => WALL(IW)
   EWC => EXTERNAL_WALL(IW)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   I   = BC%II
   J   = BC%JJ
   K   = BC%KK
   IOR = BC%IOR

   ! Apply pressure gradients at NEUMANN boundaries: dH/dn = -F_n - d(u_n)/dt

   IF_NEUMANN: IF (EWC%PRESSURE_BC_TYPE==NEUMANN) THEN

      SELECT CASE(IOR)
         CASE( 1)
            BXS(J,K) = HX(0)   *(-FVX(0,J,K)    + EWC%DUNDT)
         CASE(-1)
            BXF(J,K) = HX(IBP1)*(-FVX(IBAR,J,K) - EWC%DUNDT)
         CASE( 2)
            BYS(I,K) = HY(0)   *(-FVY(I,0,K)    + EWC%DUNDT)
         CASE(-2)
            BYF(I,K) = HY(JBP1)*(-FVY(I,JBAR,K) - EWC%DUNDT)
         CASE( 3)
            BZS(I,J) = HZ(0)   *(-FVZ(I,J,0)    + EWC%DUNDT)
         CASE(-3)
            BZF(I,J) = HZ(KBP1)*(-FVZ(I,J,KBAR) - EWC%DUNDT)
      END SELECT
   ENDIF IF_NEUMANN

   ! Apply pressures at DIRICHLET boundaries, depending on the specific type

   IF_DIRICHLET: IF (EWC%PRESSURE_BC_TYPE==DIRICHLET) THEN

      NOT_OPEN: IF (WC%BOUNDARY_TYPE/=OPEN_BOUNDARY .AND. WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) THEN

         ! Solid boundary that uses a Dirichlet BC. Assume that the pressure at the boundary (BXS, etc) is the average of the
         ! last computed pressures in the ghost and adjacent gas cells.

         SELECT CASE(IOR)
            CASE( 1) ; BXS(J,K) = 0.5_EB*(HP(0,J,K)   +HP(1,J,K))    + WALL_WORK1(IW)
            CASE(-1) ; BXF(J,K) = 0.5_EB*(HP(IBAR,J,K)+HP(IBP1,J,K)) + WALL_WORK1(IW)
            CASE( 2) ; BYS(I,K) = 0.5_EB*(HP(I,0,K)   +HP(I,1,K))    + WALL_WORK1(IW)
            CASE(-2) ; BYF(I,K) = 0.5_EB*(HP(I,JBAR,K)+HP(I,JBP1,K)) + WALL_WORK1(IW)
            CASE( 3) ; BZS(I,J) = 0.5_EB*(HP(I,J,0)   +HP(I,J,1))    + WALL_WORK1(IW)
            CASE(-3) ; BZF(I,J) = 0.5_EB*(HP(I,J,KBAR)+HP(I,J,KBP1)) + WALL_WORK1(IW)
         END SELECT

      ENDIF NOT_OPEN

      ! Interpolated boundary -- set boundary value of H to be average of neighboring cells from previous time step
      ! HP from the neighboring mesh NOM has already been copied to the external cells of mesh NM in NO_FLUX.

      INTERPOLATED_ONLY: IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN

         NOM = EWC%NOM

         SELECT CASE(IOR)
            CASE( 1)
               DX_OTHER = MESHES(NOM)%DX(EWC%IIO_MIN)
               BXS(J,K) = (DX_OTHER*HP(1,J,K) + DX(1)*HP(0,J,K))/(DX(1)+DX_OTHER) + WALL_WORK1(IW)
            CASE(-1)
               DX_OTHER = MESHES(NOM)%DX(EWC%IIO_MIN)
               BXF(J,K) = (DX_OTHER*HP(IBAR,J,K) + DX(IBAR)*HP(IBP1,J,K))/(DX(IBAR)+DX_OTHER) + WALL_WORK1(IW)
            CASE( 2)
               DY_OTHER = MESHES(NOM)%DY(EWC%JJO_MIN)
               BYS(I,K) = (DY_OTHER*HP(I,1,K) + DY(1)*HP(I,0,K))/(DY(1)+DY_OTHER) + WALL_WORK1(IW)
            CASE(-2)
               DY_OTHER = MESHES(NOM)%DY(EWC%JJO_MIN)
               BYF(I,K) = (DY_OTHER*HP(I,JBAR,K) + DY(JBAR)*HP(I,JBP1,K))/(DY(JBAR)+DY_OTHER) + WALL_WORK1(IW)
            CASE( 3)
               DZ_OTHER = MESHES(NOM)%DZ(EWC%KKO_MIN)
               BZS(I,J) = (DZ_OTHER*HP(I,J,1) + DZ(1)*HP(I,J,0))/(DZ(1)+DZ_OTHER) + WALL_WORK1(IW)
            CASE(-3)
               DZ_OTHER = MESHES(NOM)%DZ(EWC%KKO_MIN)
               BZF(I,J) = (DZ_OTHER*HP(I,J,KBAR) + DZ(KBAR)*HP(I,J,KBP1))/(DZ(KBAR)+DZ_OTHER) + WALL_WORK1(IW)
         END SELECT

      ENDIF INTERPOLATED_ONLY

      ! OPEN (passive opening to exterior of domain) boundary. Apply inflow/outflow BC.

      OPEN_IF: IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN

         B1 => BOUNDARY_PROP1(WC%B1_INDEX)
         VT => VENTS(WC%VENT_INDEX)
         IF (ABS(B1%T_IGN-T_BEGIN)<=TWO_EPSILON_EB .AND. VT%PRESSURE_RAMP_INDEX >=1) THEN
            TSI = T
         ELSE
            TSI = T - T_BEGIN
         ENDIF
         TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,VT%PRESSURE_RAMP_INDEX)
         P_EXTERNAL = TIME_RAMP_FACTOR*VT%DYNAMIC_PRESSURE

         ! Synthetic eddy method for OPEN inflow boundaries

         VEL_EDDY = 0._EB
         IF (VT%N_EDDY>0) THEN
            SELECT CASE(ABS(VT%IOR))
               CASE(1); VEL_EDDY = VT%U_EDDY(J,K)
               CASE(2); VEL_EDDY = VT%V_EDDY(I,K)
               CASE(3); VEL_EDDY = VT%W_EDDY(I,J)
            END SELECT
         ENDIF

         ! Wind inflow boundary conditions

         IF (INITIAL_SPEED>0._EB) THEN
            H0 = 0._EB
         ELSE
            H0 = 0.5_EB*(U0**2+V0**2+W0**2)
         ENDIF

         IF (OPEN_WIND_BOUNDARY) THEN
            SELECT CASE(IOR)
               CASE( 1); H0 = HP(1,J,K)    + 0.5_EB/(DT*RDXN(0)   )*(U_WIND(K) + VEL_EDDY - UU(0,   J,K))
               CASE(-1); H0 = HP(IBAR,J,K) - 0.5_EB/(DT*RDXN(IBAR))*(U_WIND(K) + VEL_EDDY - UU(IBAR,J,K))
               CASE( 2); H0 = HP(I,1,K)    + 0.5_EB/(DT*RDYN(0)   )*(V_WIND(K) + VEL_EDDY - VV(I,0,   K))
               CASE(-2); H0 = HP(I,JBAR,K) - 0.5_EB/(DT*RDYN(JBAR))*(V_WIND(K) + VEL_EDDY - VV(I,JBAR,K))
               CASE( 3); H0 = HP(I,J,1)    + 0.5_EB/(DT*RDZN(0)   )*(W_WIND(K) + VEL_EDDY - WW(I,J,0   ))
               CASE(-3); H0 = HP(I,J,KBAR) - 0.5_EB/(DT*RDZN(KBAR))*(W_WIND(K) + VEL_EDDY - WW(I,J,KBAR))
            END SELECT
         ENDIF

         SELECT CASE(IOR)
            CASE( 1)
               IF (UU(0,J,K)<0._EB) THEN
                  BXS(J,K) = P_EXTERNAL/B1%RHO_F + KRES(1,J,K)
               ELSE
                  BXS(J,K) = P_EXTERNAL/B1%RHO_F + H0
               ENDIF
            CASE(-1)
               IF (UU(IBAR,J,K)>0._EB) THEN
                  BXF(J,K) = P_EXTERNAL/B1%RHO_F + KRES(IBAR,J,K)
               ELSE
                  BXF(J,K) = P_EXTERNAL/B1%RHO_F + H0
               ENDIF
            CASE( 2)
               IF (VV(I,0,K)<0._EB) THEN
                  BYS(I,K) = P_EXTERNAL/B1%RHO_F + KRES(I,1,K)
               ELSE
                  BYS(I,K) = P_EXTERNAL/B1%RHO_F + H0
               ENDIF
            CASE(-2)
               IF (VV(I,JBAR,K)>0._EB) THEN
                  BYF(I,K) = P_EXTERNAL/B1%RHO_F + KRES(I,JBAR,K)
               ELSE
                  BYF(I,K) = P_EXTERNAL/B1%RHO_F + H0
               ENDIF
            CASE( 3)
               IF (WW(I,J,0)<0._EB) THEN
                  BZS(I,J) = P_EXTERNAL/B1%RHO_F + KRES(I,J,1)
               ELSE
                  BZS(I,J) = P_EXTERNAL/B1%RHO_F + H0
               ENDIF
            CASE(-3)
               IF (WW(I,J,KBAR)>0._EB) THEN
                  BZF(I,J) = P_EXTERNAL/B1%RHO_F + KRES(I,J,KBAR)
               ELSE
                  BZF(I,J) = P_EXTERNAL/B1%RHO_F + H0
               ENDIF
         END SELECT

      ENDIF OPEN_IF

   ENDIF IF_DIRICHLET

ENDDO WALL_CELL_LOOP
!$OMP END DO

! Compute the RHS of the Poisson equation

SELECT CASE(IPS)

   CASE(:1,4,7)
      IF (CYLINDRICAL) THEN
         !$OMP DO PRIVATE(TRM1,TRM3,TRM4)
         DO K=1,KBAR
            DO I=1,IBAR
               TRM1 = (R(I-1)*FVX(I-1,1,K)-R(I)*FVX(I,1,K))*RDX(I)*RRN(I)
               TRM3 = (FVZ(I,1,K-1)-FVZ(I,1,K))*RDZ(K)
               TRM4 = -DDDT(I,1,K)
               PRHS(I,1,K) = TRM1 + TRM3 + TRM4
            ENDDO
         ENDDO
         !$OMP END DO
      ENDIF
      IF (.NOT.CYLINDRICAL) THEN
         !$OMP DO PRIVATE(TRM1,TRM2,TRM3,TRM4)
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
         !$OMP END DO

      ENDIF

   CASE(2)  ! Switch x and y
      !$OMP DO PRIVATE(TRM1,TRM2,TRM3,TRM4)
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
      !$OMP END DO

   CASE(3,6)  ! Switch x and z
      !$OMP DO PRIVATE(TRM1,TRM2,TRM3,TRM4)
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
      !$OMP END DO

   CASE(5)  ! Switch y and z
      !$OMP DO PRIVATE(TRM1,TRM2,TRM3,TRM4)
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
      !$OMP END DO

END SELECT

!$OMP END PARALLEL

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
END SUBROUTINE PRESSURE_SOLVER_COMPUTE_RHS


SUBROUTINE PRESSURE_SOLVER_FFT(NM)

USE MESH_POINTERS
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
      IF (.NOT.TWO_D) THEN
         CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
      ELSE
         IF (.NOT.CYLINDRICAL) CALL H2CZSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
         IF (     CYLINDRICAL) CALL H2CYSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE1,WORK)
      ENDIF
   CASE(2)
      BZST = TRANSPOSE(BZS)
      BZFT = TRANSPOSE(BZF)
      CALL H3CZSS(BYS,BYF,BXS,BXF,BZST,BZFT,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HY)
   CASE(3)
      IF (.NOT.TWO_D) THEN
         BXST = TRANSPOSE(BXS)
         BXFT = TRANSPOSE(BXF)
         BYST = TRANSPOSE(BYS)
         BYFT = TRANSPOSE(BYF)
         BZST = TRANSPOSE(BZS)
         BZFT = TRANSPOSE(BZF)
         CALL H3CZSS(BZST,BZFT,BYST,BYFT,BXST,BXFT,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HZ)
      ELSE
         CALL H2CZSS(BZS,BZF,BXS,BXF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HZ)
      ENDIF
   CASE(4)
      CALL H3CSSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HX,HY)
   CASE(5)
      IF (.NOT.TWO_D) THEN
         BXST = TRANSPOSE(BXS)
         BXFT = TRANSPOSE(BXF)
         CALL H3CSSS(BXST,BXFT,BZS,BZF,BYS,BYF,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HX,HZ)
      ELSE
         CALL H2CZSS(BZS,BZF,BXS,BXF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HZ)
      ENDIF
   CASE(6)
      BXST = TRANSPOSE(BXS)
      BXFT = TRANSPOSE(BXF)
      BYST = TRANSPOSE(BYS)
      BYFT = TRANSPOSE(BYF)
      BZST = TRANSPOSE(BZS)
      BZFT = TRANSPOSE(BZF)
      CALL H3CSSS(BZST,BZFT,BYST,BYFT,BXST,BXFT,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HZ,HY)
   CASE(7)
      CALL H2CZSS(BXS,BXF,BYS,BYF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
END SELECT

!$OMP PARALLEL

SELECT CASE(IPS)
   CASE(:1,4,7)
      !$OMP DO
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               HP(I,J,K) = PRHS(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(2)
      !$OMP DO
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               HP(I,J,K) = PRHS(J,I,K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(3,6)
      !$OMP DO
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               HP(I,J,K) = PRHS(K,J,I)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(5)
      !$OMP DO
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

!$OMP DO
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
!$OMP END DO

!$OMP DO
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
!$OMP END DO

!$OMP DO
DO J=1,JBAR
   DO I=1,IBAR
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
!$OMP END DO

!$OMP END PARALLEL

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
END SUBROUTINE PRESSURE_SOLVER_FFT


!> \brief Solve a special 1-D Poisson equation for a tunnel to be used as a preconditioner for the 3-D Poisson solver
!> \details For details, refer to the Appendix in the FDS Technical Reference Guide entitled "A Special Preconditioning
!> Scheme for Solving the Poisson Equation in Tunnels."

SUBROUTINE TUNNEL_POISSON_SOLVER

USE MPI_F08
USE GLOBAL_CONSTANTS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
REAL(EB) :: RR,BXS_BAR,BXF_BAR,BXS_BAR_LEFT,BXF_BAR_RIGHT
REAL(EB), POINTER, DIMENSION(:) :: H_BAR_P,DUDT_BAR_P
INTEGER :: IERR,II,NM,I,J,K
REAL(EB) :: TNOW
TYPE (MESH_TYPE), POINTER :: M
LOGICAL :: SINGULAR_CASE

TNOW=CURRENT_TIME()

IF (PREDICTOR) THEN
   H_BAR_P => H_BAR
   DUDT_BAR_P => DUDT_BAR
ELSE
   H_BAR_P => H_BAR_S
   DUDT_BAR_P => DUDT_BAR_S
ENDIF

! For each mesh, compute the diagonal, off-diagonal, and right hand side terms of the tri-diagonal linear system of equations

MESH_LOOP_1: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   M => MESHES(NM)

   TP_RDXN(I_OFFSET(NM):I_OFFSET(NM)+M%IBAR) = M%RDXN(0:M%IBAR)
   IF (NM>1)       TP_RDXN(I_OFFSET(NM))        = 2._EB/(MESHES(NM-1)%DX(MESHES(NM-1)%IBAR)+M%DX(1))
   IF (NM<NMESHES) TP_RDXN(I_OFFSET(NM)+M%IBAR) = 2._EB/(MESHES(NM+1)%DX(1)                +M%DX(M%IBAR))

   DO I=1,M%IBAR
      II = I_OFFSET(NM) + I  ! Spatial index of the entire tunnel, not just this mesh
      TP_CC(II) = 0._EB
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            SELECT CASE(M%IPS)
               CASE DEFAULT ; TP_CC(II) = TP_CC(II) + M%PRHS(I,J,K)*M%DY(J)*M%DZ(K)
               CASE(2)      ; TP_CC(II) = TP_CC(II) + M%PRHS(J,I,K)*M%DY(J)*M%DZ(K)
               CASE(3,6)    ; TP_CC(II) = TP_CC(II) + M%PRHS(K,J,I)*M%DY(J)*M%DZ(K)
               CASE(5)      ; TP_CC(II) = TP_CC(II) + M%PRHS(I,K,J)*M%DY(J)*M%DZ(K)
            END SELECT
         ENDDO
      ENDDO
      TP_CC(II) = TP_CC(II)/((M%YF-M%YS)*(M%ZF-M%ZS))  ! RHS linear system of equations
      TP_DD(II) = -M%RDX(I)*(TP_RDXN(II)+TP_RDXN(II-1))  ! Diagonal of tri-diagonal matrix
      TP_AA(II) =  M%RDX(I)*TP_RDXN(II)    ! Upper band of matrix
      TP_BB(II) =  M%RDX(I)*TP_RDXN(II-1)  ! Lower band of matrix
      SELECT CASE(M%IPS)
         CASE DEFAULT ; M%PRHS(I,1:M%JBAR,1:M%KBAR) = M%PRHS(I,1:M%JBAR,1:M%KBAR) - TP_CC(II)  ! New RHS of the 3-D Poisson eq
         CASE(2)      ; M%PRHS(1:M%JBAR,I,1:M%KBAR) = M%PRHS(1:M%JBAR,I,1:M%KBAR) - TP_CC(II)  ! New RHS of the 3-D Poisson eq
         CASE(3,6)    ; M%PRHS(1:M%KBAR,1:M%JBAR,I) = M%PRHS(1:M%KBAR,1:M%JBAR,I) - TP_CC(II)  ! New RHS of the 3-D Poisson eq
         CASE(5)      ; M%PRHS(I,1:M%KBAR,1:M%JBAR) = M%PRHS(I,1:M%KBAR,1:M%JBAR) - TP_CC(II)  ! New RHS of the 3-D Poisson eq
      END SELECT
   ENDDO

   ! Subtract average BCs (BXS_BAR, BXF_BAR) from the 3-D BCs (BXS and BXF) for all meshes, including tunnel ends.

   BXS_BAR = 0._EB
   BXF_BAR = 0._EB
   DO K=1,M%KBAR
     DO J=1,M%JBAR
         BXS_BAR = BXS_BAR + M%BXS(J,K)*M%DY(J)*M%DZ(K)
         BXF_BAR = BXF_BAR + M%BXF(J,K)*M%DY(J)*M%DZ(K)
      ENDDO
   ENDDO
   BXS_BAR = BXS_BAR/((M%YF-M%YS)*(M%ZF-M%ZS))  ! Left boundary condition, bar(b)_x,1
   BXF_BAR = BXF_BAR/((M%YF-M%YS)*(M%ZF-M%ZS))  ! Right boundary condition, bar(b)_x,2

   M%BXS = M%BXS - BXS_BAR  ! This new BXS (b_x,1(j,k)) will be used for the 3-D pressure solve
   M%BXF = M%BXF - BXF_BAR  ! This new BXF (b_x,2(j,k)) will be used for the 3-D pressure solve

   ! Apply boundary conditions at end of tunnel to the matrix components

   IF (NM==1) THEN
      BXS_BAR_LEFT = BXS_BAR
      IF (M%LBC==FISHPAK_BC_NEUMANN_NEUMANN .OR. M%LBC==FISHPAK_BC_NEUMANN_DIRICHLET) THEN  ! Neumann BC
         TP_CC(1) = TP_CC(1) + M%DXI*BXS_BAR_LEFT*TP_BB(1)
         TP_DD(1) = TP_DD(1) + TP_BB(1)
      ELSE  ! Dirichlet BC
         TP_CC(1) = TP_CC(1) - 2._EB*BXS_BAR_LEFT*TP_BB(1)
         TP_DD(1) = TP_DD(1) - TP_BB(1)
      ENDIF
   ENDIF

   IF (NM==NMESHES) THEN
      BXF_BAR_RIGHT = BXF_BAR
      IF (M%LBC==FISHPAK_BC_NEUMANN_NEUMANN .OR. M%LBC==FISHPAK_BC_DIRICHLET_NEUMANN) THEN  ! Neumann BC
         TP_CC(TUNNEL_NXP) = TP_CC(TUNNEL_NXP) - M%DXI*BXF_BAR_RIGHT*TP_AA(TUNNEL_NXP)
         TP_DD(TUNNEL_NXP) = TP_DD(TUNNEL_NXP) + TP_AA(TUNNEL_NXP)
      ELSE  ! Dirichet BC
         TP_CC(TUNNEL_NXP) = TP_CC(TUNNEL_NXP) - 2._EB*BXF_BAR_RIGHT*TP_AA(TUNNEL_NXP)
         TP_DD(TUNNEL_NXP) = TP_DD(TUNNEL_NXP) - TP_AA(TUNNEL_NXP)
      ENDIF
   ENDIF

ENDDO MESH_LOOP_1

IF (MY_RANK>0) THEN  ! MPI processes greater than 0 send their matrix components to MPI process 0

   CALL MPI_GATHERV(TP_AA(DISPLS_TP(MY_RANK)+1),COUNTS_TP(MY_RANK),MPI_DOUBLE_PRECISION,TP_AA,COUNTS_TP,DISPLS_TP,&
                    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHERV(TP_BB(DISPLS_TP(MY_RANK)+1),COUNTS_TP(MY_RANK),MPI_DOUBLE_PRECISION,TP_BB,COUNTS_TP,DISPLS_TP,&
                    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHERV(TP_CC(DISPLS_TP(MY_RANK)+1),COUNTS_TP(MY_RANK),MPI_DOUBLE_PRECISION,TP_CC,COUNTS_TP,DISPLS_TP,&
                    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHERV(TP_DD(DISPLS_TP(MY_RANK)+1),COUNTS_TP(MY_RANK),MPI_DOUBLE_PRECISION,TP_DD,COUNTS_TP,DISPLS_TP,&
                    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

ELSE  ! MPI process 0 receives matrix components and solves tri-diagonal linear system of equations.

   CALL MPI_GATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,TP_AA,COUNTS_TP,DISPLS_TP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,TP_BB,COUNTS_TP,DISPLS_TP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,TP_CC,COUNTS_TP,DISPLS_TP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,TP_DD,COUNTS_TP,DISPLS_TP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

   TRIDIAGONAL_SOLVER_1: DO I=2,TUNNEL_NXP
      RR    = TP_BB(I)/TP_DD(I-1)
      TP_DD(I) = TP_DD(I) - RR*TP_AA(I-1)
      TP_CC(I) = TP_CC(I) - RR*TP_CC(I-1)
   ENDDO TRIDIAGONAL_SOLVER_1
   IF (ABS(TP_DD(TUNNEL_NXP))>TWO_EPSILON_EB) THEN
      TP_CC(TUNNEL_NXP) = TP_CC(TUNNEL_NXP)/TP_DD(TUNNEL_NXP)
      SINGULAR_CASE = .FALSE.
   ELSE  ! Singular matrix when both sides of tunnel have Neumann BC
      TP_CC(TUNNEL_NXP) = 1._EB
      SINGULAR_CASE = .TRUE.
   ENDIF
   TRIDIAGONAL_SOLVER_2: DO I=TUNNEL_NXP-1,1,-1
      TP_CC(I) = (TP_CC(I) - TP_AA(I)*TP_CC(I+1))/TP_DD(I)
   ENDDO TRIDIAGONAL_SOLVER_2

   ! In the SINGULAR_CASE, the solution is unique up to a constant. Set the constant so that the average value is zero.

   IF (SINGULAR_CASE) TP_CC(1:TUNNEL_NXP) = TP_CC(1:TUNNEL_NXP) - SUM(TP_CC(1:TUNNEL_NXP))/REAL(TUNNEL_NXP,EB)

ENDIF

! The solution to the tri-diagonal linear system is TP_CC. Broadcast this to all the MPI processes.

CALL MPI_BCAST(TP_CC(1),TUNNEL_NXP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

! Contruct the 1-D solution H_BAR and add boundary conditions at the ends of the tunnel.

H_BAR_P(1:TUNNEL_NXP) = TP_CC(1:TUNNEL_NXP)

! Apply boundary conditions at the ends of the tunnel.

MESH_LOOP_2: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   M => MESHES(NM)

   IF (NM==1) THEN
      IF (M%LBC==FISHPAK_BC_NEUMANN_NEUMANN .OR. M%LBC==FISHPAK_BC_NEUMANN_DIRICHLET) THEN  ! Neumann BC
         H_BAR_P(0) = H_BAR_P(1) - M%DXI*BXS_BAR_LEFT
      ELSE  ! Dirichlet BC
         H_BAR_P(0) = -H_BAR_P(1) + 2._EB*BXS_BAR_LEFT
      ENDIF
   ENDIF

   IF (NM==NMESHES) THEN
      IF (M%LBC==FISHPAK_BC_NEUMANN_NEUMANN .OR. M%LBC==FISHPAK_BC_DIRICHLET_NEUMANN) THEN  ! Neumann BC
         H_BAR_P(TUNNEL_NXP+1) = H_BAR_P(TUNNEL_NXP) + M%DXI*BXF_BAR_RIGHT
      ELSE  ! Dirichlet BC
         H_BAR_P(TUNNEL_NXP+1) = -H_BAR_P(TUNNEL_NXP) + 2._EB*BXF_BAR_RIGHT
      ENDIF
   ENDIF

ENDDO MESH_LOOP_2

! Create the array DUDT_BAR

DO I=0,TUNNEL_NXP
   DUDT_BAR_P(I) = -TP_RDXN(I)*(H_BAR_P(I+1)-H_BAR_P(I))
ENDDO

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
END SUBROUTINE TUNNEL_POISSON_SOLVER


SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS(NM)

USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE GLOBAL_CONSTANTS

INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,RHOP,P,RESIDUAL
REAL(EB), POINTER, DIMENSION(:) :: DUDT_BAR_P
INTEGER :: I,J,K
REAL(EB) :: LHSS,RHSS,TNOW

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

TNOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   HP => H
   RHOP => RHO
   IF (TUNNEL_PRECONDITIONER) DUDT_BAR_P => DUDT_BAR
ELSE
   HP => HS
   RHOP => RHOS
   IF (TUNNEL_PRECONDITIONER) DUDT_BAR_P => DUDT_BAR_S
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
            IF (TUNNEL_PRECONDITIONER) LHSS = LHSS - (DUDT_BAR_P(I_OFFSET(NM)+I) - DUDT_BAR_P(I_OFFSET(NM)+I-1))*RDX(I)
            RESIDUAL(I,J,K) = ABS(RHSS-LHSS)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   POIS_ERR = MAXVAL(RESIDUAL)
ENDIF

! Mandatory check of how well the computed pressure satisfies the inseparable Poisson equation:
! LHSS = del dot ((1/rho) del p + del K) = -del dot F - dD/dt = RHSS

IF (ITERATE_BAROCLINIC_TERM) THEN

   P => WORK7
   RESIDUAL => WORK8(1:IBAR,1:JBAR,1:KBAR)

   !$OMP PARALLEL

   !$OMP DO SCHEDULE(STATIC)
   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
            P(I,J,K) = RHOP(I,J,K)*(HP(I,J,K)-KRES(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO

   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(I,J,K,RHSS,LHSS)
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
            IF (TUNNEL_PRECONDITIONER) LHSS = LHSS - (DUDT_BAR_P(I_OFFSET(NM)+I) - DUDT_BAR_P(I_OFFSET(NM)+I-1))*RDX(I)
            RESIDUAL(I,J,K) = ABS(RHSS-LHSS)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO

   !$OMP END PARALLEL

   PRESSURE_ERROR_MAX(NM) = MAXVAL(RESIDUAL)
   PRESSURE_ERROR_MAX_LOC(:,NM) = MAXLOC(RESIDUAL)
   IF (STORE_PRESSURE_POISSON_RESIDUAL) PP_RESIDUAL(1:IBAR,1:JBAR,1:KBAR)=RESIDUAL(1:IBAR,1:JBAR,1:KBAR)

ENDIF

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
END SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS


SUBROUTINE COMPUTE_VELOCITY_ERROR(DT,NM)

! Check the maximum velocity error at a solid boundary

USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE GLOBAL_CONSTANTS, ONLY: PREDICTOR,VELOCITY_ERROR_MAX,SOLID_BOUNDARY,INTERPOLATED_BOUNDARY,VELOCITY_ERROR_MAX_LOC,T_USED,&
                            PRES_FLAG,FREEZE_VELOCITY,SOLID_PHASE_ONLY,GLMAT_FLAG,UGLMAT_FLAG,ULMAT_FLAG,TUNNEL_PRECONDITIONER,&
                            DUDT_BAR,DUDT_BAR_S,I_OFFSET

REAL(EB), INTENT(IN) :: DT
INTEGER, INTENT(IN) :: NM
INTEGER :: IW,IOR,II,JJ,KK,IIO,JJO,KKO,N_INT_CELLS,IIO1,IIO2,JJO1,JJO2,KKO1,KKO2
REAL(EB) :: TNOW,UN_NEW,UN_NEW_OTHER,VELOCITY_ERROR,DUDT,DVDT,DWDT,ITERATIVE_FACTOR,DHFCT
TYPE(OMESH_TYPE), POINTER :: OM
TYPE(MESH_TYPE), POINTER :: M2
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC

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

! Loop over wall cells and check velocity error.

CHECK_WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

   WC=>WALL(IW)

   IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY        .AND. &
       WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE CHECK_WALL_LOOP

   IF (WC%CUT_FACE_INDEX>0) CYCLE CHECK_WALL_LOOP


   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
      EWC=>EXTERNAL_WALL(IW)
      IF (EWC%AREA_RATIO<0.9_EB) CYCLE CHECK_WALL_LOOP
      OM => OMESH(EWC%NOM)
      M2 => MESHES(EWC%NOM)
   ENDIF

   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)

   II  = BC%II
   JJ  = BC%JJ
   KK  = BC%KK
   IOR = BC%IOR

   DHFCT = 1._EB
   IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
      SELECT CASE(PRES_FLAG)
      CASE(UGLMAT_FLAG,ULMAT_FLAG); DHFCT=0._EB
      CASE(GLMAT_FLAG); IF (IW<=N_EXTERNAL_WALL_CELLS) DHFCT=0._EB
      END SELECT
   ENDIF

   ! Update normal component of velocity at the mesh boundary

   IF (PREDICTOR) THEN
      SELECT CASE(IOR)
         CASE( 1)
            UN_NEW = U(II,JJ,KK)   - DT*(FVX(II,JJ,KK)   + RDXN(II)  *(H(II+1,JJ,KK)-H(II,JJ,KK))*DHFCT)
            IF (TUNNEL_PRECONDITIONER) UN_NEW = UN_NEW + DT*DUDT_BAR(I_OFFSET(NM)+II)*DHFCT
         CASE(-1)
            UN_NEW = U(II-1,JJ,KK) - DT*(FVX(II-1,JJ,KK) + RDXN(II-1)*(H(II,JJ,KK)-H(II-1,JJ,KK))*DHFCT)
            IF (TUNNEL_PRECONDITIONER) UN_NEW = UN_NEW + DT*DUDT_BAR(I_OFFSET(NM)+II-1)*DHFCT
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
            UN_NEW =0.5_EB*(U(II,JJ,KK)+US(II,JJ,KK)    -DT*(FVX(II,JJ,KK)  +RDXN(II)  *(HS(II+1,JJ,KK)-HS(II,JJ,KK))*DHFCT))
            IF (TUNNEL_PRECONDITIONER) UN_NEW = UN_NEW + 0.5*DT*DUDT_BAR_S(I_OFFSET(NM)+II)*DHFCT
         CASE(-1)
            UN_NEW =0.5_EB*(U(II-1,JJ,KK)+US(II-1,JJ,KK)-DT*(FVX(II-1,JJ,KK)+RDXN(II-1)*(HS(II,JJ,KK)-HS(II-1,JJ,KK))*DHFCT))
            IF (TUNNEL_PRECONDITIONER) UN_NEW = UN_NEW + 0.5*DT*DUDT_BAR_S(I_OFFSET(NM)+II-1)*DHFCT
         CASE( 2)
            UN_NEW =0.5_EB*(V(II,JJ,KK)+VS(II,JJ,KK)    -DT*(FVY(II,JJ,KK)  +RDYN(JJ)  *(HS(II,JJ+1,KK)-HS(II,JJ,KK))*DHFCT))
         CASE(-2)
            UN_NEW =0.5_EB*(V(II,JJ-1,KK)+VS(II,JJ-1,KK)-DT*(FVY(II,JJ-1,KK)+RDYN(JJ-1)*(HS(II,JJ,KK)-HS(II,JJ-1,KK))*DHFCT))
         CASE( 3)
            UN_NEW =0.5_EB*(W(II,JJ,KK)+WS(II,JJ,KK)    -DT*(FVZ(II,JJ,KK)  +RDZN(KK)  *(HS(II,JJ,KK+1)-HS(II,JJ,KK))*DHFCT))
         CASE(-3)
            UN_NEW =0.5_EB*(W(II,JJ,KK-1)+WS(II,JJ,KK-1)-DT*(FVZ(II,JJ,KK-1)+RDZN(KK-1)*(HS(II,JJ,KK)-HS(II,JJ,KK-1))*DHFCT))
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
                        IF (TUNNEL_PRECONDITIONER) DUDT = DUDT + DUDT_BAR(I_OFFSET(EWC%NOM)+IIO)
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%U(IIO,JJO,KKO)   + DT*DUDT
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-1)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DUDT = -OM%FVX(IIO-1,JJO,KKO) - M2%RDXN(IIO-1)*(OM%H(IIO,JJO,KKO)-OM%H(IIO-1,JJO,KKO))
                        IF (TUNNEL_PRECONDITIONER) DUDT = DUDT + DUDT_BAR(I_OFFSET(EWC%NOM)+IIO-1)
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
                        IF (TUNNEL_PRECONDITIONER) DUDT = DUDT + DUDT_BAR_S(I_OFFSET(EWC%NOM)+IIO)
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%U(IIO,JJO,KKO)+OM%US(IIO,JJO,KKO)     + DT*DUDT)
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-1)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DUDT = -OM%FVX(IIO-1,JJO,KKO) - M2%RDXN(IIO-1)*(OM%HS(IIO,JJO,KKO)-OM%HS(IIO-1,JJO,KKO))
                        IF (TUNNEL_PRECONDITIONER) DUDT = DUDT + DUDT_BAR_S(I_OFFSET(EWC%NOM)+IIO-1)
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
         UN_NEW_OTHER = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL_S
      ELSE
         UN_NEW_OTHER = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL
      ENDIF
   ENDIF

   ! Compute velocity difference

   VELOCITY_ERROR = UN_NEW - UN_NEW_OTHER
   B1%VEL_ERR_NEW = VELOCITY_ERROR
   WALL_WORK1(IW) = -SIGN(1._EB,REAL(IOR,EB))*ITERATIVE_FACTOR*VELOCITY_ERROR/(B1%RDN*DT)

   ! Save maximum velocity error

   IF (ABS(VELOCITY_ERROR)>VELOCITY_ERROR_MAX(NM)) THEN
      VELOCITY_ERROR_MAX_LOC(1,NM) = II
      VELOCITY_ERROR_MAX_LOC(2,NM) = JJ
      VELOCITY_ERROR_MAX_LOC(3,NM) = KK
      SELECT CASE(IOR)
         CASE(-1) ; VELOCITY_ERROR_MAX_LOC(1,NM) = II-1
         CASE(-2) ; VELOCITY_ERROR_MAX_LOC(2,NM) = JJ-1
         CASE(-3) ; VELOCITY_ERROR_MAX_LOC(3,NM) = KK-1
      END SELECT
      VELOCITY_ERROR_MAX(NM)       = ABS(VELOCITY_ERROR)
   ENDIF

ENDDO CHECK_WALL_LOOP

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW

END SUBROUTINE COMPUTE_VELOCITY_ERROR


END MODULE PRES


! ---------------------------------- LOCAL MATRIX SOLVER --------------------------------------------

MODULE LOCMAT_SOLVER

! Unstructured Poisson solver with Pardiso by MESH
! Using this solver eliminates penetration errors on solid boundaries,
! but still requires iteration to reduce mesh-to-mesh velocity errors.
! This solver allows for coarse-fine mesh interfaces (GLMAT does not).

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MESH_POINTERS
#ifdef WITH_MKL
USE MKL_PARDISO
#endif
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME

IMPLICIT NONE (TYPE,EXTERNAL)

INTEGER, PARAMETER :: IS_UNDEFINED  =-11
INTEGER, PARAMETER :: NNZ_STENCIL_H = 15 ! 7 Point stencil + linked cells.
INTEGER, ALLOCATABLE, DIMENSION(:)   :: NNZ_H_MAT
INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  JD_H_MAT
REAL(EB),ALLOCATABLE, DIMENSION(:,:) ::   D_H_MAT

! PARDISO solver control parameters:
INTEGER, PARAMETER :: SYMM_INDEFINITE       =-2
INTEGER, PARAMETER :: SYMM_POSITIVE_DEFINITE= 2
INTEGER, ALLOCATABLE :: IPARM( : )

! Message level:
INTEGER, SAVE ::  MSGLVL = 0

! Factor to drop DY in cylindrical axisymmetric coordinates.
REAL(EB), SAVE :: CYL_FCT

! Timing variable:
REAL(EB):: TNOW

PRIVATE

PUBLIC ULMAT_SOLVER, ULMAT_SOLVER_SETUP, FINISH_ULMAT_SOLVER

CONTAINS

SUBROUTINE ULMAT_SOLVER_SETUP(NM)

USE COMPLEX_GEOMETRY, ONLY : CC_GASPHASE,CC_CGSC
USE CC_SCALARS, ONLY : GET_H_CUTFACES
USE MEMORY_FUNCTIONS, ONLY : CHKMEMERR
INTEGER, INTENT(IN) :: NM

! Local Variables:
INTEGER :: I,J,K,IPZ,IOPZ,ICC,JCC,IW,IOR,ZBTYPE_LAST(-3:3),WALL_BTYPE,NZIM,IPZIM,IZERO,JDIM
INTEGER, POINTER :: IBAR=>NULL(),JBAR=>NULL(),KBAR=>NULL(),IBP1=>NULL(),JBP1=>NULL(),KBP1=>NULL(),&
                    ITRN=>NULL(),JTRN=>NULL(),KTRN=>NULL()
TYPE(ZONE_MESH_TYPE), POINTER :: ZM
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (WALL_TYPE), POINTER :: WC=>NULL()
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
INTEGER, PARAMETER :: NULL_BTYPE=0,DIRICHLET_BTYPE=1,NEUMANN_BTYPE=2,PERIODIC_BTYPE=3

IF (FREEZE_VELOCITY)  RETURN ! Fixed velocity soln. i.e. PERIODIC_TEST=102 => FREEZE_VELOCITY=.TRUE.
IF (SOLID_PHASE_ONLY) RETURN
TNOW=CURRENT_TIME()

! If MKL library not present stop.
#ifndef WITH_MKL
IF (MY_RANK==0) WRITE(LU_ERR,'(A)') &
'Error: MKL Library compile flag was not defined for ULMAT as pressure solver.'
! Some error - stop flag for CALL STOP_CHECK(1).
STOP_STATUS = SETUP_STOP
RETURN
#endif

! Factor to drop DY(J) in cylindrical coordinates. Soln assumes DTheta=1.
CYL_FCT = 0._EB; IF (CYLINDRICAL) CYL_FCT = 1._EB

M    =>MESHES(NM)
IBAR =>M%IBAR
JBAR =>M%JBAR
KBAR =>M%KBAR
IBP1 =>M%IBP1
JBP1 =>M%JBP1
KBP1 =>M%KBP1
ITRN =>M%ITRN
JTRN =>M%JTRN
KTRN =>M%KTRN

IF (.NOT.ALLOCATED(M%ZONE_MESH)) THEN
   ALLOCATE(M%ZONE_MESH(0:N_ZONE),STAT=IZERO)
   CALL ChkMemErr('INIT','ZONE_MESH',IZERO)
ENDIF

! Identify connected ZONEs

! translate logical array to local and fill diagonal
CONNECTED_ZONES_LOC = 0
DO IOPZ=0,N_ZONE
   DO IPZ=0,N_ZONE
      IF (IPZ==IOPZ .OR. CONNECTED_ZONES(IPZ,IOPZ,NM)) CONNECTED_ZONES_LOC(IPZ,IOPZ) = 1
   ENDDO
ENDDO

! establish sets of connected zones
DO IPZ=1,N_ZONE
   CONNECTED_ZONES_LOC = MATMUL(CONNECTED_ZONES_LOC,CONNECTED_ZONES_LOC)
ENDDO
CONNECTED_ZONES_LOC = MIN(1,CONNECTED_ZONES_LOC)

! select the parent zone as the first in the row
DO IPZ=0,N_ZONE
   ZM=>MESHES(NM)%ZONE_MESH(IPZ)
   ZM%CONNECTED_ZONE_PARENT = MINLOC(CONNECTED_ZONES_LOC(IPZ,:), DIM=1, MASK = CONNECTED_ZONES_LOC(IPZ,:)/=0) - 1
ENDDO

! Test if FFT solver can be used for this MESH

NZIM=0
ZONE_MESH_LOOP: DO IPZ=0,N_ZONE
   ZM=>M%ZONE_MESH(IPZ)
   ZM%USE_FFT=.TRUE.
   IF (ZM%CONNECTED_ZONE_PARENT/=IPZ) CYCLE ZONE_MESH_LOOP

   ! Test for multiple zones in MESH
   NZIM=NZIM+1
   IPZIM=IPZ
   IF (NZIM>1) ZM%USE_FFT=.FALSE.

   ! Test for internal wall cells
   IF (M%N_INTERNAL_WALL_CELLS>0) ZM%USE_FFT=.FALSE.

   ! Test for internal CFACES
   IF (M%N_INTERNAL_CFACE_CELLS>0) ZM%USE_FFT=.FALSE.

   ! Test external wall loop for inhomogeneous boundary types
   IF (ZM%USE_FFT) THEN
      ZBTYPE_LAST=NULL_BTYPE
      WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS
         WC=>M%WALL(IW)
         EWC=>M%EXTERNAL_WALL(IW)
         BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
         IOR = BC%IOR
         WALL_BTYPE=NULL_BTYPE
         SELECT CASE(WC%BOUNDARY_TYPE)
            CASE(SOLID_BOUNDARY,MIRROR_BOUNDARY)      ; WALL_BTYPE=NEUMANN_BTYPE
            CASE(OPEN_BOUNDARY,INTERPOLATED_BOUNDARY) ; WALL_BTYPE=DIRICHLET_BTYPE
            CASE(PERIODIC_BOUNDARY)                   ; WALL_BTYPE=PERIODIC_BTYPE
         END SELECT
         IF (ZBTYPE_LAST(IOR)==NULL_BTYPE) THEN
            ZBTYPE_LAST(IOR)=WALL_BTYPE
         ELSEIF (WALL_BTYPE/=ZBTYPE_LAST(IOR)) THEN
            ZM%USE_FFT=.FALSE.
            EXIT WALL_LOOP
         ENDIF

         ! Test for obstructions on an open boundary or mirror set to dirichlet
         IF ( (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .OR. WC%BOUNDARY_TYPE==MIRROR_BOUNDARY) &
            .AND. EWC%PRESSURE_BC_TYPE==DIRICHLET) THEN
            ZM%USE_FFT=.FALSE.
            EXIT WALL_LOOP
         ENDIF

      ENDDO WALL_LOOP
   ENDIF

   USE_ULMAT_IF: IF (.NOT.ZM%USE_FFT) THEN
      M%IPS=0 ! which implies... (see INITIALIZE_POISSON_SOLVER, here we must reallocate PRHS correctly for ULMAT)
      ITRN = IBP1
      IF (JBAR>1) JTRN = JBP1
      IF (JBAR==1) JTRN = 1
      KTRN = KBP1

      ! pressure periodic boundary conditions
      IF (FISHPAK_BC(1)==FISHPAK_BC_PERIODIC) ITRN=IBAR
      IF (FISHPAK_BC(2)==FISHPAK_BC_PERIODIC) JTRN=JBAR
      IF (FISHPAK_BC(3)==FISHPAK_BC_PERIODIC) KTRN=KBAR

      IF (ALLOCATED(M%PRHS)) DEALLOCATE(M%PRHS)
      IF (ALLOCATED(M%BXS)) DEALLOCATE(M%BXS)
      IF (ALLOCATED(M%BXF)) DEALLOCATE(M%BXF)
      IF (ALLOCATED(M%BYS)) DEALLOCATE(M%BYS)
      IF (ALLOCATED(M%BYF)) DEALLOCATE(M%BYF)
      IF (ALLOCATED(M%BZS)) DEALLOCATE(M%BZS)
      IF (ALLOCATED(M%BZF)) DEALLOCATE(M%BZF)

      ALLOCATE(M%PRHS(ITRN,JTRN,KTRN),STAT=IZERO) ; CALL ChkMemErr('INIT ULMAT','PRHS',IZERO)
      IF (JBAR>1 ) JDIM = JBP1
      IF (JBAR==1) JDIM = 1
      ALLOCATE(M%BXS(JDIM,KBP1),STAT=IZERO) ; CALL ChkMemErr('INIT ULMAT','BXS',IZERO)
      ALLOCATE(M%BXF(JDIM,KBP1),STAT=IZERO) ; CALL ChkMemErr('INIT ULMAT','BXF',IZERO)
      ALLOCATE(M%BYS(IBP1,KBP1),STAT=IZERO) ; CALL ChkMemErr('INIT ULMAT','BYS',IZERO)
      ALLOCATE(M%BYF(IBP1,KBP1),STAT=IZERO) ; CALL ChkMemErr('INIT ULMAT','BYF',IZERO)
      ALLOCATE(M%BZS(IBP1,JDIM),STAT=IZERO) ; CALL ChkMemErr('INIT ULMAT','BZS',IZERO)
      ALLOCATE(M%BZF(IBP1,JDIM),STAT=IZERO) ; CALL ChkMemErr('INIT ULMAT','BZF',IZERO)

      M%PRHS = 0._EB
      M%BXS  = 0._EB
      M%BXF  = 0._EB
      M%BYS  = 0._EB
      M%BYF  = 0._EB
      M%BZS  = 0._EB
      M%BZF  = 0._EB
   ENDIF USE_ULMAT_IF

ENDDO ZONE_MESH_LOOP

! FFT solver will be used for this mesh, no further setup required
IF (NZIM==1 .AND. M%ZONE_MESH(IPZIM)%USE_FFT) RETURN

! If mesh solver is ULMAT, initialize:
! 3. Initialize:
! 3.a Add index per zone in MUNKH array, the test goes by PRESSURE_ZONE(I,J,K), and MUNKH(I,J,K).
!     Similar to GET_MATRIX_INDEXES_H in GLOBMAT_SOLVER. Count number of unknowns ZM%NUNKH.
IF (.NOT.ALLOCATED(M%MUNKH)) ALLOCATE(M%MUNKH(0:IBAR+1,0:JBAR+1,0:KBAR+1))
M%MUNKH=IS_UNDEFINED
CALL POINT_TO_MESH(NM)
ZONE_MESH(:)%NCVLH=0
ZONE_MESH(:)%NUNKH=0
ZONE_MESH(:)%ICVL=0
ZONE_MESH(:)%IROW=0
ZONE_MESH_LOOP_2: DO IPZ=0,N_ZONE
   ZM=>ZONE_MESH(IPZ)
   IF (ZM%CONNECTED_ZONE_PARENT/=IPZ) CYCLE ZONE_MESH_LOOP_2
   ! Count NCVLH, NUNKH:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF(ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
            IF(CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
            IF (CC_IBM) THEN; IF(CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE) CYCLE; ENDIF
            ZM%NCVLH = ZM%NCVLH + 1
            ZM%NUNKH = ZM%NUNKH + 1
         ENDDO
      ENDDO
   ENDDO
   ZM%NCVLH_CART=ZM%NCVLH
   ZM%NUNKH_CART=ZM%NUNKH
   DO ICC=1,M%N_CUTCELL_MESH
      I = CUT_CELL(ICC)%IJK(IAXIS); J = CUT_CELL(ICC)%IJK(JAXIS); K = CUT_CELL(ICC)%IJK(KAXIS)
      IF(ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      IF(CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      ZM%NCVLH = ZM%NCVLH + 1
      IF(ONE_UNKH_PER_CUTCELL) THEN
         DO JCC=1,CUT_CELL(ICC)%NCELL; ZM%NUNKH = ZM%NUNKH + 1; ENDDO
      ELSE
         ZM%NUNKH = ZM%NUNKH + 1
      ENDIF
   ENDDO
   IF (ALLOCATED(ZM%MESH_IJK)) DEALLOCATE(ZM%MESH_IJK)
   ALLOCATE(ZM%MESH_IJK(1:3,1:ZM%NCVLH))
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF(ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
            IF(CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
            IF (CC_IBM) THEN; IF(CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE) CYCLE; ENDIF
            ZM%ICVL                   = ZM%ICVL + 1
            ZM%MESH_IJK(1:3,ZM%ICVL)  = (/ I, J, K/)
            ZM%IROW                   = ZM%IROW + 1
            MUNKH(I,J,K)              = ZM%IROW
         ENDDO
      ENDDO
   ENDDO
   DO ICC=1,M%N_CUTCELL_MESH
      I = CUT_CELL(ICC)%IJK(IAXIS); J = CUT_CELL(ICC)%IJK(JAXIS); K = CUT_CELL(ICC)%IJK(KAXIS)
      IF(ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      IF(CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      ZM%ICVL                  = ZM%ICVL + 1
      ZM%MESH_IJK(1:3,ZM%ICVL) = (/ I, J, K/)
      IF(ONE_UNKH_PER_CUTCELL) THEN
         DO JCC=1,CUT_CELL(ICC)%NCELL
            ZM%IROW                 = ZM%IROW + 1
            CUT_CELL(ICC)%UNKH(JCC) = ZM%IROW
         ENDDO
      ELSE
         ZM%IROW                                   = ZM%IROW + 1
         CUT_CELL(ICC)%UNKH(1:CUT_CELL(ICC)%NCELL) = ZM%IROW
      ENDIF
   ENDDO
ENDDO ZONE_MESH_LOOP_2

! 3.b Build REGFACE_H, RCF_H arrays. These face arrays are defined per mesh and axis and have
!     an integer field PRES_ZONE that provides the pressure zone the face is immersed in.
CALL ULMAT_GET_H_REGFACES(NM)
IF(CC_IBM) CALL GET_H_CUTFACES(ONE_NM=NM)

! Define Pardiso solver control parameters:
CALL ULMAT_DEFINE_IPARM

! Build zone matrix, apply BCs to it and factorize with PARDISO:
CALL POINT_TO_MESH(NM)
ZONE_MESH_LOOP_4: DO IPZ=0,N_ZONE
   ZM=>ZONE_MESH(IPZ)
   IF (ZM%CONNECTED_ZONE_PARENT/=IPZ .OR. ZM%NUNKH<1)  CYCLE ZONE_MESH_LOOP_4

   ! 3.c Per pressure zone (parent) add PRES_ZONE value to REG, RC and cut-faces, etc.
   !     Get nonzeros graph of the Poisson matrix, defined as:
   !    - NNZ_H_MAT(1:NUNKH) Number of nonzeros on per matrix row.
   !    - JD_H_MAT(1:NNZ_STENCIL_H,1:NUNKH) Column location of nonzeros.
   CALL ULMAT_MATRIXGRAPH_H(NM,IPZ) ! Define the Graph of the Matrix for Gasphase cells

   ! 3.d Build discrete Laplace operator matrix H_MAT:
   CALL ULMAT_H_MATRIX(NM,IPZ)

   ! 3.e Make changes to H_MAT due to boundary conditions (i.e. WALL faces with DIRICHLET boundary condition):
   CALL ULMAT_BCS_H_MATRIX(NM,IPZ)

   ! 3.f Pass H_MAT, NNZ_H_MAT, JD_H_MAT to CSR format and invoque LU solver:
   CALL ULMAT_H_MATRIX_LUDCMP(NM,IPZ)

ENDDO ZONE_MESH_LOOP_4

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW

RETURN
END SUBROUTINE ULMAT_SOLVER_SETUP


! ------------------------------- ULMAT_SOLVER -----------------------------------

SUBROUTINE ULMAT_SOLVER(NM,T,DT)

USE PRES, ONLY : PRESSURE_SOLVER_FFT
USE CC_SCALARS, ONLY : GET_PRES_CFACE_BCS

INTEGER, INTENT(IN) :: NM
REAL(EB),INTENT(IN) :: T,DT

! Local variables
INTEGER :: IPZ
TYPE(ZONE_MESH_TYPE), POINTER :: ZM


IF (FREEZE_VELOCITY .OR. SOLID_PHASE_ONLY) RETURN
CALL POINT_TO_MESH(NM)

! Pressure Boundary conditions due to CFACES change BXS, BXF, BYS, BYF.. in external CFACES, and
IF(CC_IBM) CALL GET_PRES_CFACE_BCS(NM,T,DT)

! Loop over zones within MESH NM and solve the unstructured Poisson problem directly.

ZONE_MESH_LOOP: DO IPZ=0,N_ZONE

   ZM=>ZONE_MESH(IPZ)
   IF (ZM%CONNECTED_ZONE_PARENT/=IPZ) CYCLE ZONE_MESH_LOOP
   IF (ZM%USE_FFT) THEN
      CALL PRESSURE_SOLVER_FFT(NM)
      RETURN
   ELSE
      IF(ZM%NUNKH<1) CYCLE ZONE_MESH_LOOP
      CALL ULMAT_SOLVE_ZONE(NM,IPZ)
   ENDIF

ENDDO ZONE_MESH_LOOP

! STOP_STATUS=USER_STOP ! on testing
END SUBROUTINE ULMAT_SOLVER

! -------------------------- ULMAT_SOLVE_ZONE -----------------------------------

SUBROUTINE ULMAT_SOLVE_ZONE(NM,IPZ)

USE COMPLEX_GEOMETRY, ONLY : CC_IDCC,CC_IDRC
USE CC_SCALARS, ONLY : GET_FN_DIVERGENCE_CUTCELL,GET_H_GUARD_CUTCELL,GRADH_ON_CARTESIAN

INTEGER, INTENT(IN) :: NM, IPZ

! Local Variables:
INTEGER :: NRHS,MAXFCT,MNUM,ERROR,I,J,K,ICC,JCC,IIG,JJG,KKG,IOR,IW,IROW,NCELL,ICFACE,IFACE,JFACE,ICVL,ILH,JLH,KLH,IRC
REAL(EB):: SUM_FH(1:2),MEAN_FH,SUM_XH(1:2),MEAN_XH,DIV_FN_VOL,DIV_FN,IDX,AF,VAL,BCV
TYPE(ZONE_MESH_TYPE), POINTER :: ZM
TYPE (WALL_TYPE),  POINTER :: WC
TYPE (EXTERNAL_WALL_TYPE),  POINTER :: EWC
TYPE (CFACE_TYPE), POINTER :: CFA
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP=>NULL()
#ifdef WITH_MKL
INTEGER :: PHASE, PERM(1)
#endif

TNOW=CURRENT_TIME()

! Solve:
NRHS   =  1
MAXFCT =  1
MNUM   =  1
ERROR  =  0 ! initialize error flag

ZM=>ZONE_MESH(IPZ)

! Build FH(:):
ZM%F_H(:) = 0._EB
! First Source on Cartesian cells:
DO ICVL=1,ZM%NCVLH_CART ! Regular Cartesian cells.
   I=ZM%MESH_IJK(IAXIS,ICVL); J=ZM%MESH_IJK(JAXIS,ICVL); K=ZM%MESH_IJK(KAXIS,ICVL)
   IROW = MUNKH(I,J,K)
   ZM%F_H(IROW) = ZM%F_H(IROW) + PRHS(I,J,K) * ((1._EB-CYL_FCT)*DY(J) + CYL_FCT*RC(I))*DX(I)*DZ(K)
ENDDO
! Then source from cut-cells:
DO ICVL=ZM%NCVLH_CART+1,ZM%NCVLH ! Cut-cells.
   I=ZM%MESH_IJK(IAXIS,ICVL); J=ZM%MESH_IJK(JAXIS,ICVL); K=ZM%MESH_IJK(KAXIS,ICVL)
   ICC=CCVAR(I,J,K,CC_IDCC);  NCELL=CUT_CELL(ICC)%NCELL
   ! Here we add div(F) in the cut-cell and DDDT:
   IF(ONE_UNKH_PER_CUTCELL) THEN
      DO JCC=1,NCELL
         CALL GET_FN_DIVERGENCE_CUTCELL(ICC,JCC,DIV_FN,SUBSTRACT_BAROCLINIC=.FALSE.)
         DIV_FN_VOL = DIV_FN*CUT_CELL(ICC)%VOLUME(JCC)
         ! Add to F_H:
         IROW = CUT_CELL(ICC)%UNKH(JCC)
         ZM%F_H(IROW) = ZM%F_H(IROW) - (CUT_CELL(ICC)%DDDTVOL(JCC) + DIV_FN_VOL)
      ENDDO
   ELSE
      DIV_FN_VOL = 0._EB
      DO JCC=1,NCELL
         CALL GET_FN_DIVERGENCE_CUTCELL(ICC,JCC,DIV_FN,SUBSTRACT_BAROCLINIC=.FALSE.)
         DIV_FN_VOL = DIV_FN_VOL + DIV_FN*CUT_CELL(ICC)%VOLUME(JCC)
      ENDDO
      ! Add to F_H:
      IROW = CUT_CELL(ICC)%UNKH(1)
      ZM%F_H(IROW) = ZM%F_H(IROW) - (CUT_CELL(ICC)%DDDTVOL(1) + DIV_FN_VOL)
   ENDIF
ENDDO

! Finally Boundary condition corrections to RHS:
! Compute FV in boundary and external CFACEs with DIRICHLET external BCs:
CFACE_LOOP : DO ICFACE=1,N_EXTERNAL_CFACE_CELLS
   CFA => CFACE(ICFACE)
   ! Here case where SOLID and OPEN or interpolated are mixed on a boundary:
   IF( CFA%BOUNDARY_TYPE==NULL_BOUNDARY .OR. CFA%BOUNDARY_TYPE==SOLID_BOUNDARY) CYCLE CFACE_LOOP
   IFACE= CFA%CUT_FACE_IND1
   JFACE= CFA%CUT_FACE_IND2
   WC  => WALL(CUT_FACE(IFACE)%IWC)
   EWC  => EXTERNAL_WALL(CUT_FACE(IFACE)%IWC)
   BC  => BOUNDARY_COORD(WC%BC_INDEX)
   ! DIRICHLET boundaries:
   IF_CFACE_DIRICHLET: IF (EWC%PRESSURE_BC_TYPE==DIRICHLET) THEN
      ! Gasphase cell indexes:
      IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG
      IF(ZONE_MESH(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE CFACE_LOOP
      IROW = MUNKH(IIG,JJG,KKG)
      IF(CC_IBM) THEN
         IF (IROW <= 0) THEN
            ICC = CCVAR(IIG,JJG,KKG,CC_IDCC); IF(ICC<1) CYCLE CFACE_LOOP
            ! Note: this only works with single pressure unknown per cartesian cell.
            IROW = CUT_CELL(ICC)%UNKH(1)
         ENDIF
      ENDIF
      IOR   = BC%IOR
      ! Define centroid to centroid distance, normal to WC:
      IF(.NOT.GRADH_ON_CARTESIAN) THEN
         IDX=1._EB/(CUT_FACE(IFACE)%XCENHIGH(ABS(IOR),JFACE)-CUT_FACE(IFACE)%XCENLOW(ABS(IOR),JFACE))
      ELSE
         SELECT CASE (IOR)
         CASE(-1); IDX = RDXN(IIG)
         CASE( 1); IDX = RDXN(IIG-1)
         CASE(-2); IDX = RDYN(JJG)
         CASE( 2); IDX = RDYN(JJG-1)
         CASE(-3); IDX = RDZN(KKG)
         CASE( 3); IDX = RDZN(KKG-1)
         END SELECT
      ENDIF
      ! Add to F_H:
      ZM%F_H(IROW) = ZM%F_H(IROW) + (-2._EB*IDX * CFA%AREA * CFA%PRES_BXN)
   ENDIF IF_CFACE_DIRICHLET
ENDDO CFACE_LOOP

! Finally add External Wall cell BCs:
WALL_CELL_LOOP_1 : DO IW=1,N_EXTERNAL_WALL_CELLS
   WC => WALL(IW)
   EWC => EXTERNAL_WALL(IW)
   ! Drop if NULL or this is a cut-face. Dealt with external CFACE.
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. WC%CUT_FACE_INDEX>0) CYCLE WALL_CELL_LOOP_1
   BC  => BOUNDARY_COORD(WC%BC_INDEX)
   ! Gasphase cell indexes:
   IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG; IOR = BC%IOR
   IF(ZONE_MESH(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE WALL_CELL_LOOP_1
   ! NEUMANN boundaries:
   IF_NEUMANN_1: IF (EWC%PRESSURE_BC_TYPE==NEUMANN) THEN

      IROW = MUNKH(IIG,JJG,KKG)
      IF(CC_IBM) THEN
         IF (IROW <= 0) THEN
            ICC = CCVAR(IIG,JJG,KKG,CC_IDCC); IF(ICC<1) CYCLE WALL_CELL_LOOP_1
            ! Note: this only works with single pressure unknown per cartesian cell.
            IROW = CUT_CELL(ICC)%UNKH(1)
         ENDIF
      ENDIF
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
      ! Add to F_H:
      ZM%F_H(IROW) = ZM%F_H(IROW) + VAL
   ENDIF IF_NEUMANN_1

   ! DIRICHLET boundaries:
   IF_DIRICHLET_1: IF (EWC%PRESSURE_BC_TYPE==DIRICHLET) THEN
      ! Here case where SOLID and OPEN or interpolated are mixed on a boundary:
      IF( WC%BOUNDARY_TYPE==SOLID_BOUNDARY .OR. WC%BOUNDARY_TYPE==MIRROR_BOUNDARY) CYCLE WALL_CELL_LOOP_1
      IROW = MUNKH(IIG,JJG,KKG)
      IF(CC_IBM) THEN
         IF (IROW <= 0) THEN
            ICC = CCVAR(IIG,JJG,KKG,CC_IDCC); IF(ICC<1) CYCLE WALL_CELL_LOOP_1
            ! Note: this only works with single pressure unknown per cartesian cell.
            IROW = CUT_CELL(ICC)%UNKH(1)
         ENDIF
      ENDIF
      ! Define cell size, normal to WC:
      ILH    = 0; JLH = 0; KLH = 0
      SELECT CASE (IOR)
      CASE(-1) ! -IAXIS oriented, high face of IIG cell.
         IDX = RDXN(IIG+ILH); BCV = BXF(JJG,KKG)
         AF  = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG  )) * DZ(KKG)
      CASE( 1) ! +IAXIS oriented, low face of IIG cell.
         ILH = -1; IDX = RDXN(IIG+ILH); BCV = BXS(JJG,KKG)
         AF  = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG-1)) * DZ(KKG)
      CASE(-2) ! -JAXIS oriented, high face of JJG cell.
         IDX = RDYN(JJG+JLH); BCV = BYF(IIG,KKG)
         AF  = DX(IIG)*DZ(KKG)
      CASE( 2) ! +JAXIS oriented, low face of JJG cell.
         JLH = -1; IDX = RDYN(JJG+JLH); BCV = BYS(IIG,KKG)
         AF  = DX(IIG)*DZ(KKG)
      CASE(-3) ! -KAXIS oriented, high face of KKG cell.
         IDX = RDZN(KKG+KLH); BCV = BZF(IIG,JJG)
         AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG)
      CASE( 3) ! +KAXIS oriented, low face of KKG cell.
         KLH = -1; IDX = RDZN(KKG+KLH); BCV = BZS(IIG,JJG)
         AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG)
      END SELECT
      ! Address case of RC face in the boundary:
      IF (CC_IBM .AND. .NOT.GRADH_ON_CARTESIAN) THEN
         IRC = FCVAR(IIG+ILH,JJG+JLH,KKG+KLH,CC_IDRC,ABS(BC%IOR))
         IF(IRC > 0) IDX = 1._EB / ( RC_FACE(IRC)%XCEN(ABS(BC%IOR),HIGH_IND) - RC_FACE(IRC)%XCEN(ABS(BC%IOR),LOW_IND) )
      ENDIF
      ! Add to F_H:
      ZM%F_H(IROW) = ZM%F_H(IROW) + (-2._EB*IDX*AF*BCV)
   ENDIF IF_DIRICHLET_1
ENDDO WALL_CELL_LOOP_1

! For indefinite matrices substract mean of source F_H:
H_INDEFINITE_IF_1 : IF (ZM%MTYPE==SYMM_INDEFINITE ) THEN
   SUM_FH(1:2) = 0._EB; MEAN_FH = 0._EB
   ! Get Arithmetic Mean of F_H:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (MUNKH(I,J,K)<=0 .OR. ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE ! Gasphase Cartesian cells.
            SUM_FH(1) = SUM_FH(1) + ZM%F_H(MUNKH(I,J,K))
            SUM_FH(2) = SUM_FH(2) + 1._EB
         ENDDO
      ENDDO
   ENDDO
   ! Add cut-cell region contribution:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      I = CUT_CELL(ICC)%IJK(IAXIS); J = CUT_CELL(ICC)%IJK(JAXIS); K = CUT_CELL(ICC)%IJK(KAXIS)
      IF (CUT_CELL(ICC)%UNKH(1)<=0 .OR. ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      IF(ONE_UNKH_PER_CUTCELL) THEN
         DO JCC=1,CUT_CELL(ICC)%NCELL
            SUM_FH(1) = SUM_FH(1) + ZM%F_H(CUT_CELL(ICC)%UNKH(JCC))
            SUM_FH(2) = SUM_FH(2) + 1._EB
         ENDDO
      ELSE
         SUM_FH(1) = SUM_FH(1) + ZM%F_H(CUT_CELL(ICC)%UNKH(1))
         SUM_FH(2) = SUM_FH(2) + 1._EB
      ENDIF
   ENDDO
   MEAN_FH = SUM_FH(1)/SUM_FH(2)
   ! IF (MY_RANK==0) WRITE(LU_ERR,*) PREDICTOR,'INDEFINITE POISSON MATRIX, IPZ, MEAN(RHS), SUM(RHS)=',&
   !                                 IPZ,MEAN_FH,SUM_FH(1:2)
   ! Substract Mean:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (MUNKH(I,J,K)<=0 .OR. ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE ! Gasphase Cartesian cells.
            ZM%F_H(MUNKH(I,J,K)) = ZM%F_H(MUNKH(I,J,K)) - MEAN_FH
         ENDDO
      ENDDO
   ENDDO
   ! Add cut-cell region contribution:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      I = CUT_CELL(ICC)%IJK(IAXIS); J = CUT_CELL(ICC)%IJK(JAXIS); K = CUT_CELL(ICC)%IJK(KAXIS)
      IF (CUT_CELL(ICC)%UNKH(1)<=0 .OR. ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      IF(ONE_UNKH_PER_CUTCELL) THEN
         DO JCC=1,CUT_CELL(ICC)%NCELL
            ZM%F_H(CUT_CELL(ICC)%UNKH(JCC)) = ZM%F_H(CUT_CELL(ICC)%UNKH(JCC)) - MEAN_FH
         ENDDO
      ELSE
         ZM%F_H(CUT_CELL(ICC)%UNKH(1)) = ZM%F_H(CUT_CELL(ICC)%UNKH(1)) - MEAN_FH
      ENDIF
   ENDDO
ENDIF H_INDEFINITE_IF_1

!.. Back substitution and iterative refinement
IPARM(8) =  0 ! max numbers of iterative refinement steps
#ifdef WITH_MKL
PHASE    = 33 ! only solving
CALL PARDISO(ZM%PT_H, MAXFCT, MNUM, ZM%MTYPE, PHASE, ZM%NUNKH, &
             ZM%A_H, ZM%IA_H, ZM%JA_H, PERM, NRHS, IPARM, MSGLVL, ZM%F_H, ZM%X_H, ERROR)
#endif
IF (ERROR /= 0) WRITE(0,*) 'ULMAT_SOLVER: The following ERROR was detected: ', ERROR

! For indefinite matrices, substract mean of solution X_H:
H_INDEFINITE_IF_2 : IF (ZM%MTYPE==SYMM_INDEFINITE ) THEN
   SUM_XH(1:2) = 0._EB; MEAN_XH = 0._EB
   ! Get Arithmetic Mean of H:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (MUNKH(I,J,K)<=0 .OR. ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE ! Gasphase Cartesian cells.
            SUM_XH(1) = SUM_XH(1) + ZM%X_H(MUNKH(I,J,K))
            SUM_XH(2) = SUM_XH(2) + 1._EB
         ENDDO
      ENDDO
   ENDDO
   ! Add cut-cell region contribution:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      I = CUT_CELL(ICC)%IJK(IAXIS); J = CUT_CELL(ICC)%IJK(JAXIS); K = CUT_CELL(ICC)%IJK(KAXIS)
      IF (CUT_CELL(ICC)%UNKH(1)<=0 .OR. ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      IF(ONE_UNKH_PER_CUTCELL) THEN
         DO JCC=1,CUT_CELL(ICC)%NCELL
            SUM_XH(1) = SUM_XH(1) + ZM%X_H(CUT_CELL(ICC)%UNKH(JCC))
            SUM_XH(2) = SUM_XH(2) + 1._EB
         ENDDO
      ELSE
         SUM_XH(1) = SUM_XH(1) + ZM%X_H(CUT_CELL(ICC)%UNKH(1))
         SUM_XH(2) = SUM_XH(2) + 1._EB
      ENDIF
   ENDDO
   MEAN_XH = SUM_XH(1)/SUM_XH(2)
   ! IF (MY_RANK==0) WRITE(LU_ERR,*) PREDICTOR,'INDEFINITE POISSON MATRIX, IPZ, MEAN(H),SUM(H)=',&
   !                                 IPZ,MEAN_XH,SUM_XH(1:2)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (MUNKH(I,J,K)<=0 .OR. ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE ! Gasphase Cartesian cells.
            ZM%X_H(MUNKH(I,J,K)) = ZM%X_H(MUNKH(I,J,K)) - MEAN_XH
         ENDDO
      ENDDO
   ENDDO
   ! Add cut-cell region contribution:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      I = CUT_CELL(ICC)%IJK(IAXIS); J = CUT_CELL(ICC)%IJK(JAXIS); K = CUT_CELL(ICC)%IJK(KAXIS)
      IF (CUT_CELL(ICC)%UNKH(1)<=0 .OR. ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      IF(ONE_UNKH_PER_CUTCELL) THEN
         DO JCC=1,CUT_CELL(ICC)%NCELL
            ZM%X_H(CUT_CELL(ICC)%UNKH(JCC)) = ZM%X_H(CUT_CELL(ICC)%UNKH(JCC)) - MEAN_XH
         ENDDO
      ELSE
         ZM%X_H(CUT_CELL(ICC)%UNKH(1)) = ZM%X_H(CUT_CELL(ICC)%UNKH(1)) - MEAN_XH
      ENDIF
   ENDDO
ENDIF H_INDEFINITE_IF_2

! WRITE(LU_ERR,*) 'SUM_XH=',SUM(ZM%X_H),SUM(ZM%A_H(1:ZM%IA_H(ZM%NUNKH+1)))

! Dump result back to mesh containers:
IF (PREDICTOR) THEN
   HP => H
ELSE
   HP => HS
ENDIF

! First Source on Cartesian cells with CC_UNKH > 0:
DO IROW=1,ZM%NUNKH_CART ! Regular Cartesian cells.
   I=ZM%MESH_IJK(IAXIS,IROW); J=ZM%MESH_IJK(JAXIS,IROW); K=ZM%MESH_IJK(KAXIS,IROW)
   HP(I,J,K) = -ZM%X_H(MUNKH(I,J,K))
ENDDO
IF (PREDICTOR) THEN
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      I = CUT_CELL(ICC)%IJK(IAXIS); J = CUT_CELL(ICC)%IJK(JAXIS); K = CUT_CELL(ICC)%IJK(KAXIS)
      IF (CUT_CELL(ICC)%UNKH(1)<=0 .OR. ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      IF(ONE_UNKH_PER_CUTCELL) THEN
         DO JCC=1,CUT_CELL(ICC)%NCELL
            CUT_CELL(ICC)%H(JCC) = -ZM%X_H(CUT_CELL(ICC)%UNKH(JCC))
         ENDDO
      ELSE
         CUT_CELL(ICC)%H(1:MESHES(NM)%CUT_CELL(ICC)%NCELL) = -ZM%X_H(CUT_CELL(ICC)%UNKH(1))
      ENDIF
      HP(I,J,K) = -ZM%X_H(CUT_CELL(ICC)%UNKH(1))
   ENDDO
ELSE
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      I = CUT_CELL(ICC)%IJK(IAXIS); J = CUT_CELL(ICC)%IJK(JAXIS); K = CUT_CELL(ICC)%IJK(KAXIS)
      IF (CUT_CELL(ICC)%UNKH(1)<=0 .OR. ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      IF(ONE_UNKH_PER_CUTCELL) THEN
         DO JCC=1,CUT_CELL(ICC)%NCELL
            CUT_CELL(ICC)%HS(JCC) = -ZM%X_H(CUT_CELL(ICC)%UNKH(JCC))
         ENDDO
      ELSE
         CUT_CELL(ICC)%HS(1:MESHES(NM)%CUT_CELL(ICC)%NCELL) = -ZM%X_H(CUT_CELL(ICC)%UNKH(1))
      ENDIF
      HP(I,J,K) = -ZM%X_H(CUT_CELL(ICC)%UNKH(1))
   ENDDO
ENDIF

! Fill external boundary conditions for Mesh, if necessary:
WALL_CELL_LOOP_2 : DO IW=1,N_EXTERNAL_WALL_CELLS
   WC => WALL(IW)
   EWC => EXTERNAL_WALL(IW)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   ! Gasphase cell indexes:
   IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG; IOR = BC%IOR
   IF (ZONE_MESH(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE WALL_CELL_LOOP_2
   ! NEUMANN boundaries:
   IF_NEUMANN_2 : IF (EWC%PRESSURE_BC_TYPE==NEUMANN) THEN
      I   = BC%II;  J   = BC%JJ;  K   = BC%KK
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
   ENDIF IF_NEUMANN_2

   ! DIRICHLET boundaries:
   IF_DIRICHLET_2 : IF (EWC%PRESSURE_BC_TYPE==DIRICHLET) THEN
      I   = BC%II;  J   = BC%JJ;  K   = BC%KK
      ! Define cell size, normal to WC:
      IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .OR. WC%BOUNDARY_TYPE==MIRROR_BOUNDARY) THEN
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
   ENDIF IF_DIRICHLET_2
ENDDO WALL_CELL_LOOP_2

IF(CC_IBM) CALL GET_H_GUARD_CUTCELL(IPZ,HP)

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW

RETURN
END SUBROUTINE ULMAT_SOLVE_ZONE

! -------------------------ULMAT_GET_H_REGFACES ---------------------------------

SUBROUTINE ULMAT_GET_H_REGFACES(NM)

USE CC_SCALARS, ONLY : GET_RCFACES_H

INTEGER, INTENT(IN) :: NM

! Local Variables:
INTEGER :: ILO,IHI,JLO,JHI,KLO,KHI
INTEGER :: I,J,K,II,IREG,X1AXIS
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IJKBUFFER
INTEGER :: IW, IIG, JJG, KKG, IOR
LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: LOG_INTWC
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(MESH_TYPE), POINTER :: M
TYPE(ZONE_MESH_TYPE), POINTER :: ZM1,ZM2

M=>MESHES(NM)

! 1. Regular GASPHASE faces connected to Gasphase cells:
IF (ALLOCATED(IJKBUFFER)) DEALLOCATE(IJKBUFFER)
ALLOCATE(IJKBUFFER(IAXIS:KAXIS,1:(IBAR+1)*(JBAR+1)*(KBAR+1)))

! Check internal SOLID_BOUNDARY faces:
IF (ALLOCATED(LOG_INTWC)) DEALLOCATE(LOG_INTWC)
ALLOCATE(LOG_INTWC(0:IBAR,0:JBAR,0:KBAR,IAXIS:KAXIS)); LOG_INTWC(:,:,:,:) = .FALSE.
WALL_LOOP_1 : DO IW=N_EXTERNAL_WALL_CELLS+1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC => WALL(IW)
   IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IIG = BC%IIG
   JJG = BC%JJG
   KKG = BC%KKG
   IOR = BC%IOR
   SELECT CASE(IOR)
   CASE( 1); LOG_INTWC(IIG-1,JJG  ,KKG  ,IAXIS) = .TRUE.
   CASE(-1); LOG_INTWC(IIG  ,JJG  ,KKG  ,IAXIS) = .TRUE.
   CASE( 2); LOG_INTWC(IIG  ,JJG-1,KKG  ,JAXIS) = .TRUE.
   CASE(-2); LOG_INTWC(IIG  ,JJG  ,KKG  ,JAXIS) = .TRUE.
   CASE( 3); LOG_INTWC(IIG  ,JJG  ,KKG-1,KAXIS) = .TRUE.
   CASE(-3); LOG_INTWC(IIG  ,JJG  ,KKG  ,KAXIS) = .TRUE.
   END SELECT
ENDDO WALL_LOOP_1

! Regular faces in axis = IAXIS : External boundary faces not counted.
X1AXIS = IAXIS
ILO = 1; IHI = IBAR-1
JLO = 1; JHI = JBAR
KLO = 1; KHI = KBAR
! First count for allocation:
IREG = 0
DO K=KLO,KHI
   DO J=JLO,JHI
      DO I=ILO,IHI
         IF (LOG_INTWC(I,J,K,X1AXIS))   CYCLE ! Wall cell in X face.
         ZM1=>M%ZONE_MESH(PRESSURE_ZONE(I,J,K))
         ZM2=>M%ZONE_MESH(PRESSURE_ZONE(I+1,J,K))
         IF ( ZM1%CONNECTED_ZONE_PARENT /= ZM2%CONNECTED_ZONE_PARENT ) CYCLE ! Solid face in X face.
         IF (MUNKH(I  ,J,K) <= 0) CYCLE ! No UNKH defined at low cell.
         IF (MUNKH(I+1,J,K) <= 0) CYCLE ! No UNKH defined at high cell.
         IREG = IREG + 1
         IJKBUFFER(IAXIS:KAXIS,IREG) = (/ I, J, K /)
      ENDDO
   ENDDO
ENDDO
M%NREGFACE_H(X1AXIS) = IREG
NULLIFY(REGFACE_IAXIS_H) ! Nullify pointer to mesh variable M%REGFACE_IAXIS_H, we are about to allocate it.
IF(ALLOCATED(M%REGFACE_IAXIS_H)) DEALLOCATE(M%REGFACE_IAXIS_H)
ALLOCATE(M%REGFACE_IAXIS_H(IREG))
DO II=1,IREG
   M%REGFACE_IAXIS_H(II)%IJK(IAXIS:KAXIS) = IJKBUFFER(IAXIS:KAXIS,II)
ENDDO

! Regular faces in axis = JAXIS : External boundary faces not counted.
X1AXIS = JAXIS
ILO = 1; IHI = IBAR
JLO = 1; JHI = JBAR-1
KLO = 1; KHI = KBAR
! First count for allocation:
IREG = 0
DO K=KLO,KHI
   DO J=JLO,JHI
      DO I=ILO,IHI
         IF (LOG_INTWC(I,J,K,X1AXIS)) CYCLE
         ZM1=>M%ZONE_MESH(PRESSURE_ZONE(I,J,K))
         ZM2=>M%ZONE_MESH(PRESSURE_ZONE(I,J+1,K))
         IF ( ZM1%CONNECTED_ZONE_PARENT /= ZM2%CONNECTED_ZONE_PARENT ) CYCLE ! Solid face in Y face.
         IF (MUNKH(I,J  ,K) <= 0) CYCLE
         IF (MUNKH(I,J+1,K) <= 0) CYCLE
         IREG = IREG + 1
         IJKBUFFER(IAXIS:KAXIS,IREG) = (/ I, J, K /)
      ENDDO
   ENDDO
ENDDO
M%NREGFACE_H(X1AXIS) = IREG
NULLIFY(REGFACE_JAXIS_H) ! Nullify pointer to mesh variable M%REGFACE_JAXIS_H, we are about to allocate it.
IF(ALLOCATED(M%REGFACE_JAXIS_H)) DEALLOCATE(M%REGFACE_JAXIS_H)
ALLOCATE(M%REGFACE_JAXIS_H(IREG))
DO II=1,IREG
   M%REGFACE_JAXIS_H(II)%IJK(IAXIS:KAXIS) = IJKBUFFER(IAXIS:KAXIS,II)
ENDDO

! Regular faces in axis = KAXIS : External boundary faces not counted.
X1AXIS = KAXIS
ILO = 1; IHI = IBAR
JLO = 1; JHI = JBAR
KLO = 1; KHI = KBAR-1
! Loop on Cartesian cells, define cut cells and solid cells CGSC:
! First count for allocation:
IREG = 0
DO K=KLO,KHI
   DO J=JLO,JHI
      DO I=ILO,IHI
         IF (LOG_INTWC(I,J,K,X1AXIS))        CYCLE
         ZM1=>M%ZONE_MESH(PRESSURE_ZONE(I,J,K))
         ZM2=>M%ZONE_MESH(PRESSURE_ZONE(I,J,K+1))
         IF ( ZM1%CONNECTED_ZONE_PARENT /= ZM2%CONNECTED_ZONE_PARENT ) CYCLE ! Solid face in Z face.
         IF (MUNKH(I,J,K  ) <= 0) CYCLE
         IF (MUNKH(I,J,K+1) <= 0) CYCLE
         IREG = IREG + 1
         IJKBUFFER(IAXIS:KAXIS,IREG) = (/ I, J, K /)
      ENDDO
   ENDDO
ENDDO
M%NREGFACE_H(X1AXIS) = IREG
NULLIFY(REGFACE_KAXIS_H) ! Nullify pointer to mesh variable M%REGFACE_KAXIS_H, we are about to allocate it.
IF(ALLOCATED(M%REGFACE_KAXIS_H)) DEALLOCATE(M%REGFACE_KAXIS_H)
ALLOCATE(M%REGFACE_KAXIS_H(IREG))
DO II=1,IREG
   M%REGFACE_KAXIS_H(II)%IJK(IAXIS:KAXIS) = IJKBUFFER(IAXIS:KAXIS,II)
ENDDO

! 2. Lists of Regular Gasphase faces, connected to one regular gasphase and one cut-cell:
IF (CC_IBM) CALL GET_RCFACES_H(NM)

DEALLOCATE(IJKBUFFER,LOG_INTWC)

END SUBROUTINE ULMAT_GET_H_REGFACES

! ------------------------- ULMAT_MATRIXGRAPH_H ---------------------------------

SUBROUTINE ULMAT_MATRIXGRAPH_H(NM,IPZ)

USE COMPLEX_GEOMETRY, ONLY : CC_GASPHASE,CC_IDCC
USE CC_SCALARS, ONLY : GET_CC_MATRIXGRAPH_H

INTEGER, INTENT(IN) :: NM,IPZ

! Local Variables:
INTEGER :: X1AXIS,IFACE,I,J,K,ICF,IND(LOW_IND:HIGH_IND)
INTEGER :: LOCROW,LOCROW_1,LOCROW_2,IIND,NII,ILOC,NUNKH,IROW,ICC
INTEGER :: NREG,IIM,JJM,KKM,IIP,JJP,KKP,LOW_FACE,HIGH_FACE,IW,II,JJ,KK,IIG,JJG,KKG
TYPE(CC_REGFACE_TYPE), POINTER, DIMENSION(:) :: REGFACE_H=>NULL()
TYPE(CC_RCFACE_TYPE), POINTER :: RCF=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
INTEGER :: WC_JD(1:2,1:2)
TYPE(ZONE_MESH_TYPE), POINTER :: ZM

ZM    =>ZONE_MESH(IPZ)
NUNKH = ZM%NUNKH

! Allocate NNZ_H_MAT, JD_H_MAT:
IF (ALLOCATED(NNZ_H_MAT)) DEALLOCATE(NNZ_H_MAT)
ALLOCATE( NNZ_H_MAT(1:NUNKH) ); NNZ_H_MAT(:) = 0

IF (ALLOCATED(JD_H_MAT)) DEALLOCATE(JD_H_MAT)
ALLOCATE( JD_H_MAT(1:NNZ_STENCIL_H,1:NUNKH) ); JD_H_MAT(:,:) = HUGE(I) ! Contains on first index nonzeros per local row.

! 1. First define PRES_ZONE=IPZ for all faces:
! Regular faces:
DO X1AXIS=IAXIS,KAXIS
   SELECT CASE(X1AXIS)
   CASE(IAXIS); REGFACE_H=>MESHES(NM)%REGFACE_IAXIS_H
   CASE(JAXIS); REGFACE_H=>MESHES(NM)%REGFACE_JAXIS_H
   CASE(KAXIS); REGFACE_H=>MESHES(NM)%REGFACE_KAXIS_H
   END SELECT
   DO IFACE=1,MESHES(NM)%NREGFACE_H(X1AXIS)
      I =REGFACE_H(IFACE)%IJK(IAXIS) ! Low Side cell
      J =REGFACE_H(IFACE)%IJK(JAXIS)
      K =REGFACE_H(IFACE)%IJK(KAXIS)
      IF(ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT==IPZ) REGFACE_H(IFACE)%PRES_ZONE=IPZ
   ENDDO
ENDDO

! RC faces:
DO IFACE=1,MESHES(NM)%CC_NRCFACE_H
   RCF => RC_FACE(MESHES(NM)%RCF_H(IFACE));
   I   = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS); K = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1)
   IF(ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT==IPZ) RCF%PRES_ZONE=IPZ
ENDDO

! Cut faces:
DO ICF=1,MESHES(NM)%N_CUTFACE_MESH
   I = CUT_FACE(ICF)%IJK(IAXIS)
   J = CUT_FACE(ICF)%IJK(JAXIS)
   K = CUT_FACE(ICF)%IJK(KAXIS)
   IF(ZONE_MESH(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT==IPZ) CUT_FACE(ICF)%PRES_ZONE=IPZ
ENDDO


! 2. Then proceed to build NNZ_H_MAT, JD_H_MAT arrays:
! Regular Faces
AXIS_LOOP_1 : DO X1AXIS=IAXIS,KAXIS
   NREG = MESHES(NM)%NREGFACE_H(X1AXIS)
   SELECT CASE(X1AXIS)
   CASE(IAXIS)
      REGFACE_H => MESHES(NM)%REGFACE_IAXIS_H
      IIM =   0; JJM = 0; KKM = 0
      IIP =   1; JJP = 0; KKP = 0
      LOW_FACE=0; HIGH_FACE=IBAR
   CASE(JAXIS)
      REGFACE_H => MESHES(NM)%REGFACE_JAXIS_H
      IIM = 0; JJM =   0; KKM = 0
      IIP = 0; JJP =   1; KKP = 0
      LOW_FACE=0; HIGH_FACE=JBAR
   CASE(KAXIS)
      REGFACE_H => MESHES(NM)%REGFACE_KAXIS_H
      IIM = 0; JJM = 0; KKM =   0
      IIP = 0; JJP = 0; KKP =   1
      LOW_FACE=0; HIGH_FACE=KBAR
   END SELECT

   IFACE_LOOP_1 : DO IFACE=1,NREG
      IF(REGFACE_H(IFACE)%PRES_ZONE/=IPZ) CYCLE IFACE_LOOP_1
      I  = REGFACE_H(IFACE)%IJK(IAXIS)
      J  = REGFACE_H(IFACE)%IJK(JAXIS)
      K  = REGFACE_H(IFACE)%IJK(KAXIS)
      ! Unknowns on related cells:
      IND(LOW_IND)  = MUNKH(I+IIM,J+JJM,K+KKM)
      IND(HIGH_IND) = MUNKH(I+IIP,J+JJP,K+KKP)
      CALL ADD_INPLACE_NNZ_H(LOW_IND,HIGH_IND,IND)
   ENDDO IFACE_LOOP_1
   NULLIFY(REGFACE_H)
ENDDO AXIS_LOOP_1

! Here check wall cells/external CFACEs of Type INTERPOLATED or OPEN_BOUNDARY and add in place-one sided where
! no nonzeros have been counted.
WALL_LOOP_1 : DO IW=1,N_EXTERNAL_WALL_CELLS
   WC => WALL(IW)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) CYCLE WALL_LOOP_1
   IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG
   IF(ZONE_MESH(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE WALL_LOOP_1 ! Cycle if pressure zone not IPZ.
   II  = BC%II;  JJ  = BC%JJ;  KK  = BC%KK
   IROW = MUNKH(IIG,JJG,KKG)
   IF(CC_IBM) THEN
      IF (IROW <= 0) THEN
         ICC = CCVAR(IIG,JJG,KKG,CC_IDCC); IF(ICC<1) CYCLE WALL_LOOP_1
         ! Note: this only works with single pressure unknown per cartesian cell.
         IROW = CUT_CELL(ICC)%UNKH(1)
      ENDIF
   ENDIF
   ! Unknowns on related cells:
   IND(LOW_IND)   = IROW  ! internal.
   IF(NNZ_H_MAT(IROW)>0) CYCLE WALL_LOOP_1
   NNZ_H_MAT(IROW) = 1
   JD_H_MAT(1,IROW)= IROW
ENDDO WALL_LOOP_1

! Finally Add nonzeros corresponding to RC_FACE, CUT_FACE
CC_IF_1 : IF (CC_IBM) THEN
   ! Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
   RCFACE_LOOP_1 : DO IFACE=1,MESHES(NM)%CC_NRCFACE_H
      RCF => RC_FACE(MESHES(NM)%RCF_H(IFACE)); IF(RCF%PRES_ZONE/=IPZ) CYCLE RCFACE_LOOP_1
      I   = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS); K = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1)
      ! Unknowns on related cells:
      IND(LOW_IND)  = RCF%UNKH(LOW_IND)
      IND(HIGH_IND) = RCF%UNKH(HIGH_IND)
      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND
      LOCROW_2 = HIGH_IND
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            IF ( I == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( I == IBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(JAXIS)
            IF ( J == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( J == JBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(KAXIS)
            IF ( K == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( K == KBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
      ENDSELECT
      CALL ADD_INPLACE_NNZ_H(LOCROW_1,LOCROW_2,IND) ! Add to matrix arrays.
   ENDDO RCFACE_LOOP_1

   CF_LOOP_1 : DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
      IF ( CUT_FACE(ICF)%STATUS/=CC_GASPHASE .OR. CUT_FACE(ICF)%IWC>0) CYCLE CF_LOOP_1
      IF ( CUT_FACE(ICF)%PRES_ZONE/=IPZ) CYCLE CF_LOOP_1
      I = CUT_FACE(ICF)%IJK(IAXIS)
      J = CUT_FACE(ICF)%IJK(JAXIS)
      K = CUT_FACE(ICF)%IJK(KAXIS)
      X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)
      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND; LOCROW_2 = HIGH_IND
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            IF ( I == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( I == IBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(JAXIS)
            IF ( J == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( J == JBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(KAXIS)
            IF ( K == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( K == KBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
      ENDSELECT
      DO IFACE=1,CUT_FACE(ICF)%NFACE
         IND(LOW_IND)  = CUT_FACE(ICF)%UNKH(LOW_IND,IFACE)
         IND(HIGH_IND) = CUT_FACE(ICF)%UNKH(HIGH_IND,IFACE)
         CALL ADD_INPLACE_NNZ_H(LOCROW_1,LOCROW_2,IND)
      ENDDO
   ENDDO CF_LOOP_1
ENDIF CC_IF_1

! 3. Finally bring back indexes associated with faces to their JD(:,:) matrix.
! Regular Faces, loop is similar to before:
AXIS_LOOP_2 : DO X1AXIS=IAXIS,KAXIS
   NREG = MESHES(NM)%NREGFACE_H(X1AXIS)
   SELECT CASE(X1AXIS)
   CASE(IAXIS)
      REGFACE_H => MESHES(NM)%REGFACE_IAXIS_H
      IIM =   0; JJM = 0; KKM = 0
      IIP =   1; JJP = 0; KKP = 0
      LOW_FACE=0; HIGH_FACE=IBAR
   CASE(JAXIS)
      REGFACE_H => MESHES(NM)%REGFACE_JAXIS_H
      IIM = 0; JJM =   0; KKM = 0
      IIP = 0; JJP =   1; KKP = 0
      LOW_FACE=0; HIGH_FACE=JBAR
   CASE(KAXIS)
      REGFACE_H => MESHES(NM)%REGFACE_KAXIS_H
      IIM = 0; JJM = 0; KKM =   0
      IIP = 0; JJP = 0; KKP =   1
      LOW_FACE=0; HIGH_FACE=KBAR
   END SELECT
   IFACE_LOOP_2 : DO IFACE=1,NREG
      IF(REGFACE_H(IFACE)%PRES_ZONE/=IPZ) CYCLE IFACE_LOOP_2
      I  = REGFACE_H(IFACE)%IJK(IAXIS)
      J  = REGFACE_H(IFACE)%IJK(JAXIS)
      K  = REGFACE_H(IFACE)%IJK(KAXIS)
      ! Unknowns on related cells:
      IND(LOW_IND)  = MUNKH(I+IIM,J+JJM,K+KKM)
      IND(HIGH_IND) = MUNKH(I+IIP,J+JJP,K+KKP)
      REGFACE_H(IFACE)%JD(1:2,1:2) = IS_UNDEFINED
      DO LOCROW = LOW_IND,HIGH_IND
         DO IIND=LOW_IND,HIGH_IND
            NII = NNZ_H_MAT(IND(LOCROW))
            DO ILOC=1,NII
               IF ( IND(IIND) == JD_H_MAT(ILOC,IND(LOCROW)) ) THEN
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
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP_2
   IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG
   IF(ZONE_MESH(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE WALL_LOOP_2 ! Cycle if pressure zone not IPZ.
   II  = BC%II;  JJ  = BC%JJ;  KK  = BC%KK
   IROW = MUNKH(IIG,JJG,KKG)
   IF(CC_IBM) THEN
      IF (IROW <= 0) THEN
         ICC = CCVAR(IIG,JJG,KKG,CC_IDCC); IF(ICC<1) CYCLE WALL_LOOP_2
         ! Note: this only works with single pressure unknown per cartesian cell.
         IROW = CUT_CELL(ICC)%UNKH(1)
      ENDIF
   ENDIF
   ! Unknowns on related cells:
   IND(LOW_IND)   = IROW  ! internal.

   WC_JD(1:2,1:2) = IS_UNDEFINED
   DO IIND=LOW_IND,HIGH_IND
      NII = NNZ_H_MAT(IND(LOCROW))
      DO ILOC=1,NII
         IF ( IND(IIND) == JD_H_MAT(ILOC,IND(LOCROW)) ) THEN
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


CC_IF_2 : IF (CC_IBM) THEN
   ! Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
   RCFACE_LOOP_2 : DO IFACE=1,MESHES(NM)%CC_NRCFACE_H
      RCF => RC_FACE(MESHES(NM)%RCF_H(IFACE)); IF(RCF%PRES_ZONE/=IPZ) CYCLE RCFACE_LOOP_2
      I   = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS); K = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1)
      ! Unknowns on related cells:
      IND(LOW_IND)  = RCF%UNKH(LOW_IND)
      IND(HIGH_IND) = RCF%UNKH(HIGH_IND)
      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND
      LOCROW_2 = HIGH_IND
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            IF ( I == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( I == IBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(JAXIS)
            IF ( J == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( J == JBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(KAXIS)
            IF ( K == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( K == KBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
      ENDSELECT
      RCF%JDH(1:2,1:2) = 0
      ! Add to global matrix arrays:
      DO LOCROW=LOCROW_1,LOCROW_2
         DO IIND=LOW_IND,HIGH_IND
            NII = NNZ_H_MAT(IND(LOCROW))
            DO ILOC=1,NII
               IF ( IND(IIND) == JD_H_MAT(ILOC,IND(LOCROW)) ) THEN
                   RCF%JDH(LOCROW,IIND) = ILOC
                   EXIT
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO RCFACE_LOOP_2

   CF_LOOP_3 : DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
      IF ( CUT_FACE(ICF)%STATUS/=CC_GASPHASE .OR. CUT_FACE(ICF)%IWC>0) CYCLE CF_LOOP_3
      IF ( CUT_FACE(ICF)%PRES_ZONE/=IPZ) CYCLE CF_LOOP_3
      I = CUT_FACE(ICF)%IJK(IAXIS)
      J = CUT_FACE(ICF)%IJK(JAXIS)
      K = CUT_FACE(ICF)%IJK(KAXIS)
      X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)
      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND; LOCROW_2 = HIGH_IND
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            IF ( I == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( I == IBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(JAXIS)
            IF ( J == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( J == JBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(KAXIS)
            IF ( K == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( K == KBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
      ENDSELECT
      CUT_FACE(ICF)%JDH(:,:,:) = 0
      DO IFACE=1,CUT_FACE(ICF)%NFACE
         IND(LOW_IND)  = CUT_FACE(ICF)%UNKH(LOW_IND,IFACE)
         IND(HIGH_IND) = CUT_FACE(ICF)%UNKH(HIGH_IND,IFACE)
         ! Add to global matrix arrays:
         DO LOCROW=LOCROW_1,LOCROW_2
            DO IIND=LOW_IND,HIGH_IND
               NII = NNZ_H_MAT(IND(LOCROW))
               DO ILOC=1,NII
                  IF ( IND(IIND) == JD_H_MAT(ILOC,IND(LOCROW)) ) THEN
                        CUT_FACE(ICF)%JDH(LOCROW,IIND,IFACE) = ILOC
                        EXIT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO CF_LOOP_3
ENDIF CC_IF_2

RETURN
END SUBROUTINE ULMAT_MATRIXGRAPH_H


! ------------------------ ADD_INPLACE_NNZ_H ----------------------------

SUBROUTINE ADD_INPLACE_NNZ_H(LOCROW_1,LOCROW_2,IND)

INTEGER, INTENT(IN) :: LOCROW_1,LOCROW_2,IND(LOW_IND:HIGH_IND)

! Local Variables:
INTEGER LOCROW, IIND, NII, ILOC, JLOC
LOGICAL INLIST

LOCROW_LOOP : DO LOCROW=LOCROW_1,LOCROW_2
   DO IIND=LOW_IND,HIGH_IND
      NII = NNZ_H_MAT(IND(LOCROW))
      ! Check that column index hasn't been already counted:
      INLIST = .FALSE.
      DO ILOC=1,NII
         IF ( IND(IIND) == JD_H_MAT(ILOC,IND(LOCROW)) ) THEN
            INLIST = .TRUE.
            EXIT
         ENDIF
      ENDDO
      IF ( INLIST ) CYCLE

      ! Now add in place:
      NII = NII + 1
      DO ILOC=1,NII
          IF ( JD_H_MAT(ILOC,IND(LOCROW)) > IND(IIND) ) EXIT
      ENDDO
      DO JLOC=NII,ILOC+1,-1
          JD_H_MAT(JLOC,IND(LOCROW)) = JD_H_MAT(JLOC-1,IND(LOCROW))
      ENDDO
      NNZ_H_MAT(IND(LOCROW))   = NII
      JD_H_MAT(ILOC,IND(LOCROW)) = IND(IIND)
   ENDDO
ENDDO LOCROW_LOOP

RETURN
END SUBROUTINE ADD_INPLACE_NNZ_H


! ---------------------------------- ULMAT_H_MATRIX ----------------------------------
SUBROUTINE ULMAT_H_MATRIX(NM,IPZ)

USE COMPLEX_GEOMETRY, ONLY : CC_GASPHASE
USE CC_SCALARS, ONLY : GRADH_ON_CARTESIAN

INTEGER, INTENT(IN) :: NM,IPZ

! Local Variables:
INTEGER :: X1AXIS,X2AXIS,X3AXIS,IFACE,I,J,K,I1,I2,I3,ICF,IND(LOW_IND:HIGH_IND),IROW,JLOC,JCOL
INTEGER :: LOCROW_1,LOCROW_2,ILOC,NUNKH
INTEGER :: NREG,IIM,JJM,KKM,IIP,JJP,KKP,LOW_FACE,HIGH_FACE
REAL(EB):: AF,IDX,BIJ,KFACE(2,2)
TYPE(CC_REGFACE_TYPE), POINTER, DIMENSION(:) :: REGFACE_H=>NULL()
TYPE(CC_RCFACE_TYPE), POINTER :: RCF=>NULL()
TYPE(ZONE_MESH_TYPE), POINTER :: ZM
REAL(EB), POINTER, DIMENSION(:)   :: DX1,DX2,DX3

ZM=>ZONE_MESH(IPZ)
NUNKH=ZM%NUNKH

! Allocate D_H_MAT:
IF (ALLOCATED(D_H_MAT)) DEALLOCATE(D_H_MAT)
ALLOCATE( D_H_MAT(1:NNZ_STENCIL_H,1:NUNKH) ); D_H_MAT(:,:)  = 0._EB

! Regular Faces
AXIS_LOOP_1 : DO X1AXIS=IAXIS,KAXIS
   NREG = MESHES(NM)%NREGFACE_H(X1AXIS)
   SELECT CASE(X1AXIS)
   CASE(IAXIS)
      REGFACE_H => MESHES(NM)%REGFACE_IAXIS_H
      IIM = 0; JJM = 0; KKM = 0
      IIP = 1; JJP = 0; KKP = 0
      LOW_FACE=0; HIGH_FACE=IBAR
      X2AXIS=JAXIS; X3AXIS=KAXIS
      DX1 => DXN
      DX2 => DY
      DX3 => DZ
   CASE(JAXIS)
      REGFACE_H => MESHES(NM)%REGFACE_JAXIS_H
      IIM = 0; JJM =   0; KKM = 0
      IIP = 0; JJP =   1; KKP = 0
      LOW_FACE=0; HIGH_FACE=JBAR
      X2AXIS=KAXIS; X3AXIS=IAXIS
      DX1 => DYN
      DX2 => DZ
      DX3 => DX
   CASE(KAXIS)
      REGFACE_H => MESHES(NM)%REGFACE_KAXIS_H
      IIM = 0; JJM = 0; KKM =   0
      IIP = 0; JJP = 0; KKP =   1
      LOW_FACE=0; HIGH_FACE=KBAR
      X2AXIS=IAXIS; X3AXIS=JAXIS
      DX1 => DZN
      DX2 => DX
      DX3 => DY
   END SELECT
   IFACE_LOOP_1 : DO IFACE=1,NREG
      IF(REGFACE_H(IFACE)%PRES_ZONE/=IPZ) CYCLE IFACE_LOOP_1
      I  = REGFACE_H(IFACE)%IJK(IAXIS)
      J  = REGFACE_H(IFACE)%IJK(JAXIS)
      K  = REGFACE_H(IFACE)%IJK(KAXIS)
      I1 = REGFACE_H(IFACE)%IJK(X1AXIS)
      I2 = REGFACE_H(IFACE)%IJK(X2AXIS)
      I3 = REGFACE_H(IFACE)%IJK(X3AXIS)
      ! Unknowns on related cells:
      IND(LOW_IND)  = MUNKH(I+IIM,J+JJM,K+KKM)
      IND(HIGH_IND) = MUNKH(I+IIP,J+JJP,K+KKP)
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
            IROW = IND(ILOC)                   ! Unknown number.
            JCOL = REGFACE_H(IFACE)%JD(ILOC,JLOC); ! Local position of coef in D_MAT_H
            ! Add coefficient:
            D_H_MAT(JCOL,IROW) = D_H_MAT(JCOL,IROW) + KFACE(ILOC,JLOC)
         ENDDO
      ENDDO
   ENDDO IFACE_LOOP_1
   NULLIFY(REGFACE_H)
ENDDO AXIS_LOOP_1


! Contribution to Laplacian matrix from Cut-cells:
CC_IF : IF ( CC_IBM ) THEN
   ! Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
   RCFACE_LOOP : DO IFACE=1,MESHES(NM)%CC_NRCFACE_H
      RCF => RC_FACE(MESHES(NM)%RCF_H(IFACE)); IF(RCF%PRES_ZONE/=IPZ) CYCLE RCFACE_LOOP
      I   = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS); K = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1)
      ! Unknowns on related cells:
      IND(LOW_IND)  = RCF%UNKH(LOW_IND)
      IND(HIGH_IND) = RCF%UNKH(HIGH_IND)
      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND; LOCROW_2 = HIGH_IND
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            AF  = DY(J)*DZ(K); IDX = RDXN(I)
            IF ( I == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( I == IBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(JAXIS)
            AF  = DX(I)*DZ(K); IDX = RDYN(J)
            IF ( J == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( J == JBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(KAXIS)
            AF  = DX(I)*DY(J); IDX = RDZN(K)
            IF ( K == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( K == KBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
      ENDSELECT
      IF(.NOT.GRADH_ON_CARTESIAN) IDX = 1._EB / ( RCF%XCEN(X1AXIS,HIGH_IND) - RCF%XCEN(X1AXIS,LOW_IND) )
      ! Now add to Adiff corresponding coeff:
      BIJ   = IDX*AF
      !    Cols 1,2: ind(LOW_IND) ind(HIGH_IND), Rows 1,2: ind_loc(LOW_IND) ind_loc(HIGH_IND)
      KFACE(1,1) = BIJ; KFACE(2,1) =-BIJ; KFACE(1,2) =-BIJ; KFACE(2,2) = BIJ
      DO ILOC=LOCROW_1,LOCROW_2      ! Local row number in Kface
         DO JLOC=LOW_IND,HIGH_IND    ! Local col number in Kface, JD
             IROW=IND(ILOC)          ! Process Local Unknown number.
             JCOL=RCF%JDH(ILOC,JLOC) ! Local position of coef in D_MAT_H
             ! Add coefficient:
             D_H_MAT(JCOL,IROW) = D_H_MAT(JCOL,IROW) + KFACE(ILOC,JLOC)
         ENDDO
      ENDDO
   ENDDO RCFACE_LOOP

   ! Now Gasphase CUT_FACES:
   CF_LOOP_1 : DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
         IF ( CUT_FACE(ICF)%STATUS/=CC_GASPHASE .OR. CUT_FACE(ICF)%IWC>0) CYCLE CF_LOOP_1
         IF ( CUT_FACE(ICF)%PRES_ZONE/=IPZ) CYCLE CF_LOOP_1
         I = CUT_FACE(ICF)%IJK(IAXIS)
         J = CUT_FACE(ICF)%IJK(JAXIS)
         K = CUT_FACE(ICF)%IJK(KAXIS)
         X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)
      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND; LOCROW_2 = HIGH_IND
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            IDX = RDXN(I)
            IF ( I == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( I == IBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(JAXIS)
            IDX = RDYN(J)
            IF ( J == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( J == JBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(KAXIS)
            IDX = RDZN(K)
            IF ( K == 0    ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( K == KBAR ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
      ENDSELECT
      DO IFACE=1,MESHES(NM)%CUT_FACE(ICF)%NFACE
         IND(LOW_IND)  = CUT_FACE(ICF)%UNKH(LOW_IND,IFACE)
         IND(HIGH_IND) = CUT_FACE(ICF)%UNKH(HIGH_IND,IFACE)
         AF = CUT_FACE(ICF)%AREA(IFACE)
         IF(.NOT.GRADH_ON_CARTESIAN) IDX= 1._EB/ ( CUT_FACE(ICF)%XCENHIGH(X1AXIS,IFACE) - &
                                                   CUT_FACE(ICF)%XCENLOW(X1AXIS, IFACE) )
         ! Now add to Adiff corresponding coeff:
         BIJ   = IDX*AF
         !    Cols 1,2: ind(LOW_IND) ind(HIGH_IND), Rows 1,2: ind_loc(LOW_IND) ind_loc(HIGH_IND)
         KFACE(1,1) = BIJ; KFACE(2,1) =-BIJ; KFACE(1,2) =-BIJ; KFACE(2,2) = BIJ
         DO ILOC=LOCROW_1,LOCROW_2 ! Local row number in Kface
            DO JLOC=LOW_IND,HIGH_IND ! Local col number in Kface, JD
                  IROW=IND(ILOC)
                  JCOL=CUT_FACE(ICF)%JDH(ILOC,JLOC,IFACE)
                  ! Add coefficient:
                  D_H_MAT(JCOL,IROW) = D_H_MAT(JCOL,IROW) + KFACE(ILOC,JLOC)
            ENDDO
         ENDDO
      ENDDO
   ENDDO CF_LOOP_1

ENDIF CC_IF

RETURN
END SUBROUTINE ULMAT_H_MATRIX


! -------------------------------- ULMAT_BCS_H_MATRIX ----------------------------------
SUBROUTINE ULMAT_BCS_H_MATRIX(NM,IPZ)

USE COMPLEX_GEOMETRY, ONLY : CC_IDCC, CC_IDRC
USE CC_SCALARS, ONLY : GET_CFACE_OPEN_BC_COEF,GRADH_ON_CARTESIAN
INTEGER, INTENT(IN) :: NM,IPZ

! Local Variables:
INTEGER :: DUM,H_MAT_IVEC
INTEGER :: JLOC,JCOL,IND(LOW_IND:HIGH_IND),ILH,JLH,KLH,IRC
REAL(EB):: AF,IDX,BIJ
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(ZONE_MESH_TYPE), POINTER :: ZM
INTEGER :: IIG,JJG,KKG,II,JJ,KK,IW,IROW,ICC

DUM=NM
ZM=>ZONE_MESH(IPZ)

! Dirichlet condition counter, to define if matrix is positive definite or indefinite.
H_MAT_IVEC = 0

WALL_LOOP_1 : DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
   WC => WALL(IW)
   ! Only OPEN_BOUNDARY or INTERPOLATED_BOUNDARY leads to a Dirichlet BC for H.
   ! Everything else leads to Neuman BCs on H, no need to modify D_MAT_HP.
   IF ( .NOT.(WC%BOUNDARY_TYPE==OPEN_BOUNDARY .OR. WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) ) CYCLE WALL_LOOP_1
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG
   IF (ZONE_MESH(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE WALL_LOOP_1
   II  = BC%II;  JJ  = BC%JJ;  KK  = BC%KK
   ! Unknowns on related cells:
   IROW = MUNKH(IIG,JJG,KKG)
   IF(CC_IBM) THEN
      IF (IROW <= 0) THEN
         ICC = CCVAR(IIG,JJG,KKG,CC_IDCC); IF(ICC<1) CYCLE WALL_LOOP_1
         ! Note: this only works with single pressure unknown per cartesian cell.
         IROW = CUT_CELL(ICC)%UNKH(1)
      ENDIF
   ENDIF
   IND(LOW_IND)  = IROW
   ILH           = 0; JLH = 0; KLH = 0
   SELECT CASE(BC%IOR)
   CASE( IAXIS)
      AF  = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG-1)) * DZ(KKG); ILH = -1; IDX= RDXN(IIG+ILH)
   CASE(-IAXIS)
      AF  = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG  )) * DZ(KKG);           IDX= RDXN(IIG+ILH)
   CASE( JAXIS)
      AF  = DX(IIG)*DZ(KKG); JLH= -1; IDX= RDYN(JJG+JLH)
   CASE(-JAXIS)
      AF  = DX(IIG)*DZ(KKG);          IDX= RDYN(JJG+JLH)
   CASE( KAXIS)
      AF  = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG); KLH = -1; IDX= RDZN(KKG+KLH)
   CASE(-KAXIS)
      AF  = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG);           IDX= RDZN(KKG+KLH)
   END SELECT
   IF (CC_IBM .AND. .NOT.GRADH_ON_CARTESIAN) THEN
      IRC = FCVAR(IIG+ILH,JJG+JLH,KKG+KLH,CC_IDRC,ABS(BC%IOR))
      IF(IRC > 0) IDX = 1._EB / ( RC_FACE(IRC)%XCEN(ABS(BC%IOR),HIGH_IND) - RC_FACE(IRC)%XCEN(ABS(BC%IOR),LOW_IND) )
   ENDIF

   ! Now add to Adiff corresponding coeff:
   BIJ = IDX*AF
   ! Case of unstructured projection:
   IF(WC%CUT_FACE_INDEX>0) CALL GET_CFACE_OPEN_BC_COEF(WC%CUT_FACE_INDEX,BC%IOR,IDX,BIJ)

   ! Find diagonal column number:
   JCOL = -1
   DO JLOC = 1,NNZ_H_MAT(IND(LOW_IND))
      IF (IND(LOW_IND) == JD_H_MAT(JLOC,IND(LOW_IND))) THEN
         JCOL = JLOC
         EXIT
      ENDIF
   ENDDO
   ! Add diagonal coefficient due to DIRICHLET BC:
   D_H_MAT(JCOL,IND(LOW_IND)) = D_H_MAT(JCOL,IND(LOW_IND)) + 2._EB*BIJ

   ! Add to mesh dirichlet bc counter
   H_MAT_IVEC = H_MAT_IVEC + 1

ENDDO WALL_LOOP_1

! Now change PARDISO Matrix type of matrix either positive definite of indefinite:
ZM%MTYPE   = SYMM_INDEFINITE ! Initialize to real symmetric indefinite.
IF (H_MAT_IVEC>0) ZM%MTYPE   = SYMM_POSITIVE_DEFINITE ! Real symmetric positive definite.

RETURN
END SUBROUTINE ULMAT_BCS_H_MATRIX

! -------------------------------- ULMAT_DEFINE_IPARM ------------------------------------

SUBROUTINE ULMAT_DEFINE_IPARM

!..
!.. SET UP PARDISO CONTROL PARAMETER
!..
IF(ALLOCATED(IPARM)) RETURN
ALLOCATE(IPARM(64)); IPARM(:) = 0

IPARM(1)  = 1  ! no solver default
IPARM(2)  = 2  ! fill-in reordering from METIS
IPARM(4)  = 0  ! no iterative-direct algorithm
IPARM(5)  = 0  ! no user fill-in reducing permutation
IPARM(6)  = 0  ! =0 solution on the first n components of x
IPARM(8)  = 0  ! numbers of iterative refinement steps
IPARM(10) =13  ! perturb the pivot elements with 1E-13
IPARM(11) = 1  ! use nonsymmetric permutation and scaling MPS
IPARM(13) = 1  ! maximum weighted matching algorithm is switched-off (default for symmetric).
              ! Try IPARM(13) = 1 in case of inappropriate accuracy
IPARM(14) = 0  ! Output: number of perturbed pivots
IPARM(18) = 0  ! -1 Output: number of nonzeros in the factor LU
IPARM(19) = 0  ! -1 Output: Mflops for LU factorization
IPARM(20) = 0  ! Output: Numbers of CG Iterations
IPARM(21) = 1  ! 1x1 diagonal pivoting for symmetric indefinite matrices.
IPARM(24) = 0
IPARM(27) = 1  ! Check matrix

RETURN
END SUBROUTINE ULMAT_DEFINE_IPARM

! ------------------------------- ULMAT_H_MATRIX_LUDCMP ----------------------------------
SUBROUTINE ULMAT_H_MATRIX_LUDCMP(NM,IPZ)

INTEGER, INTENT(IN) :: NM,IPZ

! Local Variables:
INTEGER :: INNZ, IROW, JCOL
#ifdef WITH_MKL
INTEGER :: PHASE, PERM(1)
INTEGER :: I
#endif
!.. All other variables
INTEGER MAXFCT, MNUM, NRHS, ERROR, TOT_NNZ_H
TYPE(ZONE_MESH_TYPE), POINTER :: ZM

INNZ=NM
ZM=>ZONE_MESH(IPZ)
! Define parameters:
NRHS   = 1
MAXFCT = 1
MNUM   = 1

! Set level MSG to 1 for factorization:
IF(CHECK_POISSON) THEN
   MSGLVL = 1
   IF(MY_RANK==0) WRITE(LU_ERR,*) 'ULMAT : PARDISO factorization for MESH,ZONE=',NM,IPZ,ZM%NUNKH
ENDIF
ERROR     = 0 ! initialize error flag

! Each MPI process builds its local set of rows.
! Matrix blocks defined on CRS distributed format.
! Total number of nonzeros for JD_MAT_H, D_MAT_H:
TOT_NNZ_H = SUM( NNZ_H_MAT(1:ZM%NUNKH) )

! Allocate A_H IA_H and JA_H matrices, considering all matrix coefficients:
IF (ALLOCATED(ZM%A_H))  DEALLOCATE(ZM%A_H)
IF (ALLOCATED(ZM%IA_H)) DEALLOCATE(ZM%IA_H)
IF (ALLOCATED(ZM%JA_H)) DEALLOCATE(ZM%JA_H)
ALLOCATE ( ZM%A_H(TOT_NNZ_H) , ZM%IA_H(ZM%NUNKH+1) , ZM%JA_H(TOT_NNZ_H) )

! Store upper triangular part of symmetric D_MAT_H in CSR format for the ZONE_MESH (ZM%IA_H,ZM%JA_H,ZM%A_H):
INNZ = 0
DO IROW=1,ZM%NUNKH
   ZM%IA_H(IROW) = INNZ + 1
   DO JCOL=1,NNZ_H_MAT(IROW)
      IF ( JD_H_MAT(JCOL,IROW) < IROW ) CYCLE ! Only upper Triangular part.
      INNZ = INNZ + 1
      ZM%A_H(INNZ)  =  D_H_MAT(JCOL,IROW)
      ZM%JA_H(INNZ) = JD_H_MAT(JCOL,IROW)
   ENDDO
ENDDO
ZM%IA_H(ZM%NUNKH+1) = INNZ + 1

! OPEN(unit=20,file="IJKUNKH_H_ULMAT.txt",action="write",status="replace")
! DO IROW=1,ZM%NUNKH
!    WRITE(20,'(4I6)') ZM%MESH_IJK(IAXIS:KAXIS,IROW),IROW
! ENDDO
! CLOSE(20)
! WRITE(0,*) 'IJKUNKH file written...'
! STOP

! OPEN(unit=20,file="Matrix_H_ULMAT.txt",action="write",status="replace")
! DO IROW=1,ZM%NUNKH
!    DO JCOL=1,NNZ_H_MAT(IROW)
!       WRITE(20,'(2I6,F18.12)') IROW,JD_H_MAT(JCOL,IROW),D_H_MAT(JCOL,IROW)
!    ENDDO
! ENDDO
! ! WRITE(20,'(A)') 'EOF'
! CLOSE(20)
! WRITE(0,*) 'H Matrix file written...'
! STOP

! Deallocate NNZ_H_MAT, D_H_MAT, JD_H_MAT:
DEALLOCATE(NNZ_H_MAT, D_H_MAT, JD_H_MAT)

! Allocate Solution and RHS vectors:
IF (ALLOCATED(ZM%X_H)) DEALLOCATE(ZM%X_H)
IF (ALLOCATED(ZM%F_H)) DEALLOCATE(ZM%F_H)
ALLOCATE( ZM%X_H(ZM%NUNKH) , ZM%F_H(ZM%NUNKH) ); ZM%F_H(:) = 0._EB; ZM%X_H(:) = 0._EB

! PARDISO:
! Initialize solver pointer for H matrix solves:

IF (.NOT.ALLOCATED(ZM%PT_H)) ALLOCATE(ZM%PT_H(64))

#ifdef WITH_MKL

DO I=1,64
  ZM%PT_H(I)%DUMMY = 0
ENDDO

! Reorder and Symbolic factorization:
PHASE = 11
CALL PARDISO (ZM%PT_H, MAXFCT, MNUM, ZM%MTYPE, PHASE, ZM%NUNKH, &
              ZM%A_H, ZM%IA_H, ZM%JA_H, PERM, NRHS, IPARM, MSGLVL, ZM%F_H, ZM%X_H, ERROR)

IF (ERROR /= 0) THEN
   IF (MY_RANK==0) THEN
   WRITE(LU_ERR,'(A,I5)') 'ULMAT_H_MATRIX_LUDCMP PARDISO Sym Factor: The following ERROR was detected: ', ERROR
   ! Some error - stop flag for CALL STOP_CHECK(1).
   IF(ERROR==-2) WRITE(LU_ERR,'(A)') 'Insufficient Memory for Poisson Matrix Factorization.'
   ENDIF
   STOP_STATUS = SETUP_STOP
   RETURN
END IF

! Numerical Factorization.
PHASE = 22 ! only factorization
CALL PARDISO (ZM%PT_H, MAXFCT, MNUM, ZM%MTYPE, PHASE, ZM%NUNKH, &
              ZM%A_H, ZM%IA_H, ZM%JA_H, PERM, NRHS, IPARM, MSGLVL, ZM%F_H, ZM%X_H, ERROR)

IF (ERROR /= 0) THEN
   IF (MY_RANK==0) THEN
   WRITE(LU_ERR,'(A,I5)') 'ULMAT_H_MATRIX_LUDCMP PARDISO Num Factor: The following ERROR was detected: ', ERROR
   ! Some error - stop flag for CALL STOP_CHECK(1).
   IF(ERROR==-2) WRITE(LU_ERR,'(A)') 'Insufficient Memory for Poisson Matrix Factorization.'
   ENDIF
   STOP_STATUS = SETUP_STOP
   RETURN
ENDIF

#endif

! Set level MSG to 0 for solution:
IF(CHECK_POISSON) MSGLVL = 0

RETURN
END SUBROUTINE ULMAT_H_MATRIX_LUDCMP


SUBROUTINE FINISH_ULMAT_SOLVER(NM)

INTEGER, INTENT(IN) :: NM

! Local variables:
INTEGER :: IPZ,MAXFCT,MNUM,NRHS,ERROR,MSGLVL
TYPE(ZONE_MESH_TYPE), POINTER :: ZM
#ifdef WITH_MKL
INTEGER :: PHASE,PERM(1)
#endif

IF (FREEZE_VELOCITY .OR. SOLID_PHASE_ONLY) RETURN

! Solve:
NRHS   =  1
MAXFCT =  1
MNUM   =  1
ERROR  =  0 ! initialize error flag
MSGLVL =  0 ! print statistical information

CALL POINT_TO_MESH(NM)

! Loop over zones within MESH NM and deallocate ZONE_MESH Pardiso internal arrays.
ZONE_MESH_LOOP: DO IPZ=0,N_ZONE
   ZM=>ZONE_MESH(IPZ)
   IF(ZM%USE_FFT .OR. .NOT.ALLOCATED(ZM%PT_H)) CYCLE ZONE_MESH_LOOP
   ! Finalize Pardiso:
#ifdef WITH_MKL
   PHASE = -1 ! Free memory.
   CALL PARDISO(ZM%PT_H, MAXFCT, MNUM, ZM%MTYPE, PHASE, ZM%NUNKH, &
                ZM%A_H, ZM%IA_H, ZM%JA_H, PERM, NRHS, IPARM, MSGLVL, ZM%F_H, ZM%X_H, ERROR)
#endif /* WITH_MKL */
   IF (ALLOCATED(ZM%A_H))  DEALLOCATE(ZM%A_H)
   IF (ALLOCATED(ZM%IA_H)) DEALLOCATE(ZM%IA_H)
   IF (ALLOCATED(ZM%JA_H)) DEALLOCATE(ZM%JA_H)
   IF (ALLOCATED(ZM%X_H))  DEALLOCATE(ZM%X_H)
   IF (ALLOCATED(ZM%F_H))  DEALLOCATE(ZM%F_H)
ENDDO ZONE_MESH_LOOP

RETURN
END SUBROUTINE FINISH_ULMAT_SOLVER

END MODULE LOCMAT_SOLVER


! ---------------------------------- GLOBAL MATRIX SOLVER --------------------------------------------

MODULE GLOBMAT_SOLVER

! Module that contains global matrix vector builds for Poisson equation on gas-cells only when PRES_ON_WHOLE_DOMAIN=.FALSE.
! Builds Matrices and RHS entries per MPI process in parallel.
! Calls MKL sparse cluster solver or Pardiso for the time being.

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MESH_POINTERS

USE COMPLEX_GEOMETRY, ONLY : CALL_FOR_GLMAT, CC_CGSC,CC_FGSC, CC_UNKH, CC_NCVARS,         &
                             NM_START,IPARM,NNZ_ROW_H,CALL_FROM_GLMAT_SETUP
USE CC_SCALARS, ONLY :   GET_H_CUTFACES, GET_BOUNDFACE_GEOM_INFO_H, ADD_INPLACE_NNZ_H_WHLDOM, &
                         COPY_CC_HS_TO_UNKH, COPY_CC_UNKH_TO_HS

#ifdef WITH_MKL
USE MKL_CLUSTER_SPARSE_SOLVER
#endif /* WITH_MKL */

IMPLICIT NONE (TYPE,EXTERNAL)

! These definitions are the same as geom.f90:
INTEGER,  PARAMETER :: NGUARD= 2 ! Two layers of guard-cells.

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

! Cartesian Cell centered variables, actual case initialized as CC_IBM=.FALSE.:
INTEGER :: CGSC=IS_CGSC, UNKH=IS_UNKH, NCVARS=IS_NCVARS

! Define CC pointers:
TYPE(CC_CUTCELL_TYPE), POINTER :: CC=>NULL()

! Pardiso or Sparse cluster solver message level:
INTEGER, SAVE :: MSGLVL = 0  ! 0 no messages, 1 print statistical information

! Factor to drop DY in cylindrical axisymmetric coordinates.
REAL(EB), SAVE :: CYL_FCT

! Pressure zone loops index:
INTEGER :: IPZ

! Handle for ZONE_SOLVER array entries:
TYPE(ZONE_SOLVE_TYPE), POINTER :: ZSL

!#define SINGLE_PRECISION_PSN_SOLVE

! Timing variable:
REAL(EB):: TNOW

PRIVATE

PUBLIC GLMAT_SOLVER_SETUP,GLMAT_SOLVER,COPY_H_OMESH_TO_MESH,FINISH_GLMAT_SOLVER,PRESSURE_SOLVER_CHECK_RESIDUALS_U

CONTAINS

! --------------------------- GLMAT_SOLVER -------------------------------------

SUBROUTINE GLMAT_SOLVER(T,DT)

USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE CC_SCALARS, ONLY : GET_CUTCELL_HP,GET_PRES_CFACE_BCS,GET_FH_FROM_PRHS_AND_BCS
USE MPI_F08

REAL(EB), INTENT(IN) :: T,DT

! Local Variables:
INTEGER :: MAXFCT, MNUM, NRHS, ERROR
#ifdef WITH_MKL
INTEGER :: PERM(1), PHASE
#endif
INTEGER :: NM, IW, IIG, JJG, KKG, IOR, IROW, I, J, K, ICC
TYPE (WALL_TYPE), POINTER :: WC=>NULL()
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
REAL(EB), POINTER, DIMENSION(:,:,:)   :: HP
REAL(EB) :: SUM_FH(2), SUM_XH(2), MEAN_FH, MEAN_XH
INTEGER :: IERR

! INTEGER  :: JCOL
! REAL(EB) :: LHS
! CHARACTER(30) :: FILE_NAME
! INTEGER :: ICC, IERR

IF (CC_IBM) CALL_FOR_GLMAT = .TRUE.
! Fixed velocity soln. i.e. PERIODIC_TEST=102 => FREEZE_VELOCITY=.TRUE.
IF (FREEZE_VELOCITY .OR. SOLID_PHASE_ONLY) RETURN
TNOW=CURRENT_TIME()

! Solve:
NRHS   =  1
MAXFCT =  1
MNUM   =  1
ERROR  =  0 ! initialize error flag

! Pressure BCs:
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   ! Pressure Boundary conditions due to CFACES change BXS, BXF, BYS, BYF.. in external CFACES, and
   CALL GET_PRES_CFACE_BCS(NM,T,DT)
ENDDO

IPZ_LOOP : DO IPZ=0,N_ZONE

   ZSL => ZONE_SOLVE(IPZ)

   IF (ZSL%NUNKH_TOTAL==0) CYCLE

   ! Define rhs F_H, here we use Source and BCs populated on PRESSURE_SOLVER:
   ZSL%F_H = 0._EB
   ZSL%X_H = 0._EB

   ! Main Mesh Loop:
   MESH_LOOP_1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (ZSL%NUNKH_LOCAL==0) CYCLE
      CALL POINT_TO_MESH(NM)
      ! Build FH(:):
      CALL GET_FH_FROM_PRHS_AND_BCS(NM,DT,CYL_FCT,UNKH,ZSL%NUNKH_LOCAL,IPZ,ZSL%F_H)
   ENDDO MESH_LOOP_1

   IF (ZSL%MTYPE==-2) THEN
      SUM_FH = 0._EB; MEAN_FH = 0._EB
      WHOLE_DOM_IF1 : IF(.NOT.PRES_ON_WHOLE_DOMAIN) THEN
         ! Sum source F_H by Pressure Zone:
         DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
            CALL POINT_TO_MESH(NM)
            DO K=1,KBAR
               DO J=1,JBAR
                  DO I=1,IBAR
                     IF (CCVAR(I,J,K,UNKH)<=0 .OR. ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
                     ! Row number:
                     IROW = CCVAR(I,J,K,UNKH) - ZSL%UNKH_IND(NM_START) ! Local numeration.
                     ! Sum FH:
                     SUM_FH(1) = SUM_FH(1) + ZSL%F_H(IROW)
                     SUM_FH(2) = SUM_FH(2) + 1._EB
                  ENDDO
               ENDDO
            ENDDO
            ! Add cut-cell region contribution:
            DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
               CC => CUT_CELL(ICC); I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
               IF (ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
               IROW      = CC%UNKH(1)- ZSL%UNKH_IND(NM_START) ! Local numeration.
               SUM_FH(1) = SUM_FH(1) + ZSL%F_H(IROW)
               SUM_FH(2) = SUM_FH(2) + 1._EB
            ENDDO
         ENDDO
         IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,SUM_FH(1),2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
         ! Compute arithmetic mean by pressure zone:
         MEAN_FH = SUM_FH(1)/(SUM_FH(2)+TWO_EPSILON_EB)
         ! Substract Mean:
         ZSL%F_H = ZSL%F_H - MEAN_FH
      ELSE WHOLE_DOM_IF1
         SUM_FH(1) = SUM(ZSL%F_H(1:ZSL%NUNKH_LOCAL))
         IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(SUM_FH(1),SUM_FH(2),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
         MEAN_FH = SUM_FH(2)/REAL(ZSL%NUNKH_TOTAL,EB)
         ! Substract Mean:
         ZSL%F_H = ZSL%F_H - MEAN_FH
      ENDIF WHOLE_DOM_IF1
   ENDIF
   ! WRITE(LU_ERR,*) 'SUM_FH=',SUM(F_H),H_MATRIX_INDEFINITE

#ifdef WITH_MKL
   ! Dump local low an high rows assembled by this process in IPARM:
   IPARM(41) = ZSL%LOWER_ROW
   IPARM(42) = ZSL%UPPER_ROW
   !.. Back substitution and iterative refinement
   IPARM(8) =  0 ! max numbers of iterative refinement steps
   PHASE    = 33 ! only solving
#ifdef SINGLE_PRECISION_PSN_SOLVE
   ZSL%F_H_FB(1:ZSL%NUNKH_LOCAL) = REAL(ZSL%F_H(1:ZSL%NUNKH_LOCAL),FB)
   ZSL%X_H_FB(1:ZSL%NUNKH_LOCAL) = 0._FB
   CALL CLUSTER_SPARSE_SOLVER(ZSL%PT_H, MAXFCT, MNUM, ZSL%MTYPE, PHASE, ZSL%NUNKH_TOTAL, &
   ZSL%A_H_FB, ZSL%IA_H, ZSL%JA_H, PERM, NRHS, IPARM, MSGLVL, ZSL%F_H_FB, ZSL%X_H_FB, MPI_COMM_WORLD, ERROR)
   ZSL%X_H(1:ZSL%NUNKH_LOCAL) = REAL(ZSL%X_H_FB(1:ZSL%NUNKH_LOCAL),EB)
#else
   CALL CLUSTER_SPARSE_SOLVER(ZSL%PT_H, MAXFCT, MNUM, ZSL%MTYPE, PHASE, ZSL%NUNKH_TOTAL, &
   ZSL%A_H, ZSL%IA_H, ZSL%JA_H, PERM, NRHS, IPARM, MSGLVL, ZSL%F_H, ZSL%X_H, MPI_COMM_WORLD, ERROR)
#endif
IF (ERROR /= 0 .AND. MY_RANK==0) WRITE(LU_ERR,*) 'GLMAT_SOLVER: The following ERROR was detected: ', ERROR
#endif

   IF (ZSL%MTYPE==-2) THEN
      SUM_XH = 0._EB; MEAN_XH = 0._EB
      WHOLE_DOM_IF2 : IF(.NOT.PRES_ON_WHOLE_DOMAIN) THEN
         ! Sum H by Pressure Zone:
         DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
            CALL POINT_TO_MESH(NM)
            DO K=1,KBAR
               DO J=1,JBAR
                  DO I=1,IBAR
                     IF (CCVAR(I,J,K,UNKH)<=0 .OR. ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
                     IROW = CCVAR(I,J,K,UNKH) - ZSL%UNKH_IND(NM_START) ! Local numeration.
                     SUM_XH(1) = SUM_XH(1) + ZSL%X_H(IROW)
                     SUM_XH(2) = SUM_XH(2) + 1._EB
                  ENDDO
               ENDDO
            ENDDO
            ! Add cut-cell region contribution:
            DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
               CC => CUT_CELL(ICC); I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
               IF (ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
               IROW      = CC%UNKH(1)- ZSL%UNKH_IND(NM_START) ! Local numeration.
               SUM_XH(1) = SUM_XH(1) + ZSL%X_H(IROW)
               SUM_XH(2) = SUM_XH(2) + 1._EB
            ENDDO
         ENDDO
         IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,SUM_XH(1),2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
         ! Compute arithmetic mean by pressure zone:
         MEAN_XH = SUM_XH(1)/(SUM_XH(2)+TWO_EPSILON_EB)
         ! Substract Mean:
         ZSL%X_H = ZSL%X_H - MEAN_XH
      ELSE WHOLE_DOM_IF2
         SUM_XH(1) = SUM(ZSL%X_H(1:ZSL%NUNKH_LOCAL))
         IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(SUM_XH(1),SUM_XH(2),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
         MEAN_XH = SUM_XH(2)/REAL(ZSL%NUNKH_TOTAL,EB)
         ! Substract Mean:
         ZSL%X_H = ZSL%X_H - MEAN_XH
      ENDIF WHOLE_DOM_IF2
   ENDIF
   ! WRITE(LU_ERR,*) 'SUM_XH=',SUM(X_H),SUM(A_H(1:IA_H(NUNKH_LOCAL+1)))

   ! Dump result back to mesh containers:
   MESH_LOOP_2 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      IF (PREDICTOR) THEN
         HP => H
      ELSE
         HP => HS
      ENDIF
      ! First Cartesian cells with CC_UNKH > 0:
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CCVAR(I,J,K,UNKH) <= 0 .OR. ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
               IROW = CCVAR(I,J,K,UNKH) - ZSL%UNKH_IND(NM_START) ! Local numeration.
               ! Assign to HP:
               HP(I,J,K) = -ZSL%X_H(IROW)
            ENDDO
         ENDDO
      ENDDO
      IF (CC_IBM) CALL GET_CUTCELL_HP(NM,IPZ,HP)

      ! Fill external boundary conditions for Mesh, if necesary:
      WALL_CELL_LOOP_2: DO IW=1,N_EXTERNAL_WALL_CELLS
         WC => WALL(IW)
         EWC => EXTERNAL_WALL(IW)
         BC => BOUNDARY_COORD(WC%BC_INDEX)

         IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG
         IF (ZONE_SOLVE(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
         I = BC%II; J = BC%JJ; K = BC%KK; IOR = BC%IOR

         ! NEUMANN boundaries:
         IF_NEUMANN2: IF (EWC%PRESSURE_BC_TYPE==NEUMANN) THEN
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
         IF_DIRICHLET2: IF (EWC%PRESSURE_BC_TYPE==DIRICHLET) THEN
            IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                WC%BOUNDARY_TYPE==        NULL_BOUNDARY ) CYCLE ! No need for these, that's the whole point of a
                                                                ! global solve.
            ! Define cell size, normal to WC:
            IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .OR. WC%BOUNDARY_TYPE==MIRROR_BOUNDARY) THEN
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
ENDDO IPZ_LOOP

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW

RETURN
END SUBROUTINE GLMAT_SOLVER

! ------------------------- GLMAT_SOLVER_SETUP ----------------------------------

SUBROUTINE GLMAT_SOLVER_SETUP(STAGE_FLAG)

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME

INTEGER, INTENT(IN) :: STAGE_FLAG

! Local Variables:
LOGICAL :: SUPPORTED_MESH=.TRUE.

IF (FREEZE_VELOCITY)  RETURN ! Fixed velocity soln. i.e. PERIODIC_TEST=102 => FREEZE_VELOCITY=.TRUE.
IF (SOLID_PHASE_ONLY) RETURN
TNOW=CURRENT_TIME()
SELECT CASE(STAGE_FLAG)
CASE(1)

    CALL_FROM_GLMAT_SETUP = .TRUE.

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
      CGSC   = CC_CGSC
      UNKH   = CC_UNKH
      NCVARS = CC_NCVARS
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

   ! 4. CC_GASPHASE cut-faces:
   IF(CC_IBM) CALL GET_H_CUTFACES

   ! 5. Exchange information at block boundaries for RC_FACE, CUT_FACE
   ! fields on each mesh:
   CALL GET_BOUNDFACE_GEOM_INFO_H

   ! 6. Get nonzeros graph of the Poisson matrix, defined as:
   !    - NNZ_D_MAT_H(1:NUNKH_LOCAL) Number of nonzeros on per matrix row.
   !    - JD_MAT_H(1:NNZ_ROW_H,1:NUNKH_LOCAL) Column location of nonzeros, global numeration.
   CALL GET_MATRIXGRAPH_H_WHLDOM ! Define the Graph of the Matrix for Gasphase cells on whole domain.

   ! 7. Build discrete Laplace operator matrix:
   CALL GET_H_MATRIX

   ! 8. Make changes to H_MATRIX due to boundary conditions (i.e. WALL faces or CC_INBOUNDARY faces
   ! with DIRICHLET boundary condition):
   CALL GET_BCS_H_MATRIX

   ! 9. Pass D_MAT_H, NNZ_D_MAT_H, JD_MAT_H to CSR format and invoque LU solver:
   CALL GET_H_MATRIX_LUDCMP

   CALL_FROM_GLMAT_SETUP = .FALSE.

END SELECT

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW

RETURN
END SUBROUTINE GLMAT_SOLVER_SETUP


! ---------------------------- CHECK_UNSUPPORTED_MESH -------------------------------

SUBROUTINE CHECK_UNSUPPORTED_MESH(SUPPORTED_MESH)

USE MPI_F08
USE MESH_POINTERS
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
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()

INTEGER :: NOM,IW,NMLOC,NSETS,ISET,PIVOT,PIVOT_LOC,MESHES_LEFT,CTMSH_LO,CTMSH_HI

SUPPORTED_MESH = .TRUE.

! 1. Stretched grids which is untested:
GLMAT_IF_1 : IF(PRES_FLAG==GLMAT_FLAG) THEN
TRN_ME(1:2) = 0
MESH_LOOP_TRN : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   TRN_ME(1) = TRN_ME(1) + TRANS(NM)%NOCMAX
ENDDO MESH_LOOP_TRN
TRN_ME(2)=TRN_ME(1)
IF (N_MPI_PROCESSES > 1) CALL MPI_ALLREDUCE(TRN_ME(1),TRN_ME(2),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
IF (TRN_ME(2) > 0) THEN ! There is a TRNX, TRNY or TRNZ line defined for stretched grids. Not Unsupported.
   IF (MY_RANK == 0) WRITE(LU_ERR,*) 'GLMAT Setup Error : Stretched grids currently unsupported.'
   SUPPORTED_MESH = .FALSE.
   STOP_STATUS = SETUP_STOP
   RETURN
ENDIF
ENDIF GLMAT_IF_1

IF (NMESHES == 1) RETURN

! 2. Two different cell sizes in mesh (i.e. different refinement levels):
NM = 1
IF (MY_RANK==PROCESS(NM)) THEN
   CALL POINT_TO_MESH(NM)
   DX_P(IAXIS) = DX(1)
   DX_P(JAXIS) = DY(1)
   DX_P(KAXIS) = DZ(1)
ENDIF
IF (N_MPI_PROCESSES > 1) CALL MPI_BCAST(DX_P(1),3,MPI_DOUBLE_PRECISION,PROCESS(NM),MPI_COMM_WORLD,IERR)
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
   CALL POINT_TO_MESH(NM)
   IF(ABS(DX_P(IAXIS)-DX(1)) > 10._EB*TWO_EPSILON_EB*LX) TRN_ME(1) = TRN_ME(1) + 1
   IF(ABS(DX_P(JAXIS)-DY(1)) > 10._EB*TWO_EPSILON_EB*LY) TRN_ME(1) = TRN_ME(1) + 1
   IF(ABS(DX_P(KAXIS)-DZ(1)) > 10._EB*TWO_EPSILON_EB*LZ) TRN_ME(1) = TRN_ME(1) + 1
ENDDO MESH_LOOP_CELL
TRN_ME(2)=TRN_ME(1)
IF (N_MPI_PROCESSES > 1) CALL MPI_ALLREDUCE(TRN_ME(1),TRN_ME(2),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
IF (TRN_ME(2) > 0) THEN ! Meshes at different refinement levels. Not Unsupported.
   IF (MY_RANK == 0) WRITE(LU_ERR,*) 'GLMAT Setup Error : Meshes at different refinement levels unsupported.'
   SUPPORTED_MESH = .FALSE.
   STOP_STATUS = SETUP_STOP
   RETURN
ENDIF

! 3. Two (or more) disjoint domains, where at least one has all Neumann BCs and one has some Dirichlet bcs.
! This is a topological problem that would require different Matrix types (i.e. one positive definite and one
! indefinite), which would require separate solutions.
! A possible approach to look at is to solve the whole system as indefinite, and then substract a constant in
! zones with Dirichlet condition, s.t. the value of H is zero in open boundaries.
GLMAT_IF_2 : IF(PRES_FLAG==GLMAT_FLAG) THEN
! 1. Build global lists of other connected meshes:
ALLOCATE(MESH_GRAPH(1:6,NMESHES)); MESH_GRAPH(:,:) = 0
MESH_LOOP_GRAPH : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   COUNT=0
   DO NOM=1,NMESHES
      IF(MESHES(NM)%CONNECTED_MESH(NOM))THEN
         COUNT=COUNT+1
         MESH_GRAPH(COUNT,NM) = NOM
      ENDIF
   ENDDO
ENDDO MESH_LOOP_GRAPH
IF (N_MPI_PROCESSES > 1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,MESH_GRAPH(1,1),6*NMESHES,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)

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
      IF (MY_RANK/=PROCESS(NM)) CYCLE
      CALL POINT_TO_MESH(NM)

      ! Now for Mesh NM test for Dirichlet Pressure external BCs. Assume External Dirichlet is only related to
      ! OPEN_BOUNDARY condition.
      DO IW=1,N_EXTERNAL_WALL_CELLS
         WC => WALL(IW)
         EWC => EXTERNAL_WALL(IW)
         IF (EWC%PRESSURE_BC_TYPE==DIRICHLET .AND. WC%BOUNDARY_TYPE==OPEN_BOUNDARY) DIRI_SET(ISET) = 1
      ENDDO
   ENDDO
ENDDO SETS_LOOP

IF(N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,DIRI_SET(1),NSETS,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)

! IF (MY_RANK==0) THEN
!    WRITE(LU_ERR,*) ' '
!    WRITE(LU_ERR,*) 'NSETS=',NSETS
!    DO ISET=1,NSETS
!       WRITE(LU_ERR,*) 'ISET=',ISET,DSETS(HIGH_IND,ISET)-DSETS(LOW_IND,ISET)+1,DIRI_SET(ISET)
!    ENDDO
! ENDIF

! Finally do test:
IF (ANY(DIRI_SET(1:NSETS) == 0)) THEN
   IF (MY_RANK==0) WRITE(LU_ERR,*) 'GLMAT Setup Error : Unsupported disjoint domains present on the model. Consider using ULMAT.'
   DEALLOCATE(MESH_GRAPH,DSETS,MESH_LIST,COUNTED,DIRI_SET)
   SUPPORTED_MESH = .FALSE.
   STOP_STATUS = SETUP_STOP
   RETURN
ENDIF
DEALLOCATE(MESH_GRAPH,DSETS,MESH_LIST,COUNTED,DIRI_SET)
ENDIF GLMAT_IF_2

RETURN
END SUBROUTINE CHECK_UNSUPPORTED_MESH

! ----------------------------- COPY_H_OMESH_TO_MESH --------------------------------

SUBROUTINE COPY_H_OMESH_TO_MESH

USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE CC_SCALARS, ONLY: GET_H_GUARD_CUTCELL
! Local Variables:
INTEGER  :: NM,NOM,II,JJ,KK,IOR,IW,IIO,JJO,KKO
TYPE (OMESH_TYPE), POINTER :: OM=>NULL()
TYPE (WALL_TYPE), POINTER :: WC=>NULL()
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()
LOGICAL :: FLG

IF (CC_IBM) CALL_FOR_GLMAT = .FALSE.

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN
TNOW=CURRENT_TIME()
! Loop:
PREDCORR_LOOP : IF (PREDICTOR) THEN

   MESH_LOOP_1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

      CALL POINT_TO_MESH(NM)

      ! Loop over all cell edges

      EXTERNAL_WALL_LOOP_1 : DO IW=1,N_EXTERNAL_WALL_CELLS

         WC=>WALL(IW)
         BC=>BOUNDARY_COORD(WC%BC_INDEX)
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

         II  = BC%II
         JJ  = BC%JJ
         KK  = BC%KK
         IOR = BC%IOR
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

      IF(CC_IBM) CALL GET_H_GUARD_CUTCELL(0,H)

   ENDDO MESH_LOOP_1

ELSE ! PREDCORR_LOOP

   MESH_LOOP_2 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

      CALL POINT_TO_MESH(NM)

      ! Loop over all cell edges

      EXTERNAL_WALL_LOOP_2 : DO IW=1,N_EXTERNAL_WALL_CELLS

         WC=>WALL(IW)
         BC=>BOUNDARY_COORD(WC%BC_INDEX)
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

         II  = BC%II
         JJ  = BC%JJ
         KK  = BC%KK
         IOR = BC%IOR
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

      IF(CC_IBM) CALL GET_H_GUARD_CUTCELL(0,HS)
   ENDDO MESH_LOOP_2

ENDIF PREDCORR_LOOP

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW

RETURN
END SUBROUTINE COPY_H_OMESH_TO_MESH

! ------------------------------- COPY_HS_IN_CCVAR ----------------------------------

SUBROUTINE COPY_HS_IN_CCVAR(VAR_CC)

USE MESH_POINTERS
INTEGER, INTENT(IN) :: VAR_CC

! Local Variables:
INTEGER  :: NM,NOM,II,JJ,KK,IOR,IW,IIO,JJO,KKO
TYPE (OMESH_TYPE), POINTER :: OM=>NULL()
TYPE (WALL_TYPE), POINTER :: WC=>NULL()
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()
LOGICAL :: FLG

MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)

   !
   IF (CC_IBM .AND. VAR_CC==UNKH) THEN

      CALL COPY_CC_HS_TO_UNKH(NM)

   ELSE

      ! Loop over external wall cells:
      EXTERNAL_WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

         WC=>WALL(IW)
         BC=>BOUNDARY_COORD(WC%BC_INDEX)
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

         II  = BC%II
         JJ  = BC%JJ
         KK  = BC%KK
         IOR = BC%IOR
         NOM = EWC%NOM
         ! Here if NOM==0 means it is an OBST laying on an external boundary -> CYCLE
         IF(NOM < 1) CYCLE EXTERNAL_WALL_LOOP
         OM => OMESH(NOM)

         ! This assumes all meshes at the same level of refinement:
         KKO=EWC%KKO_MIN
         JJO=EWC%JJO_MIN
         IIO=EWC%IIO_MIN

         CCVAR(II,JJ,KK,VAR_CC) = INT(OM%HS(IIO,JJO,KKO))

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

USE MESH_POINTERS
INTEGER, INTENT(IN) :: VAR_CC

! Local Variables:
INTEGER :: NM

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   HS(0:IBP1,0:JBP1,0:KBP1) = REAL(CCVAR(0:IBP1,0:JBP1,0:KBP1,VAR_CC),EB)

   ! Now cut-cells add their single Cartesian UNKH value in HS:
   IF(CC_IBM .AND. VAR_CC==UNKH) CALL COPY_CC_UNKH_TO_HS(NM)

ENDDO

RETURN
END SUBROUTINE COPY_CCVAR_IN_HS

! ------------------------------- GET_H_MATRIX_LUDCMP -------------------------------

SUBROUTINE GET_H_MATRIX_LUDCMP
#ifdef WITH_MKL
USE MPI_F08
#endif
! Local Variables:
INTEGER :: INNZ, IROW, JCOL
#ifdef WITH_MKL
INTEGER :: PHASE, PERM(1)
INTEGER :: I, IPROC
#endif
!.. All other variables
INTEGER MAXFCT, MNUM, NRHS, ERROR
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

ERROR     = 0 ! initialize error flag

#ifdef WITH_MKL
CALL SET_CLUSTER_SOLVER_IPARM
#endif

IPZ_LOOP : DO IPZ=0,N_ZONE

   ZSL => ZONE_SOLVE(IPZ); IF(ZSL%NUNKH_TOTAL==0) CYCLE

   ! Each MPI process builds its local set of rows.
   ! Matrix blocks defined on CRS distributed format.
   ! Total number of nonzeros for ZSL%JD_MAT_H, ZSL%D_MAT_H:
   ZSL%TOT_NNZ_H = 1; IF(ZSL%NUNKH_LOCAL>0) ZSL%TOT_NNZ_H = SUM( ZSL%NNZ_D_MAT_H(1:ZSL%NUNKH_LOCAL) )

   ! Allocate F_H ans H_H for this Process and IPZ:
   ALLOCATE( ZSL%X_H(MAX(ZSL%NUNKH_LOCAL,1)) , ZSL%F_H(MAX(ZSL%NUNKH_LOCAL,1)) ); ZSL%F_H=0._EB; ZSL%X_H=0._EB

   !--- This matrix definitoin used with MKL cluster solver -----
   ! Allocate A_H IA_H and JA_H matrices, considering all matrix coefficients:
   ALLOCATE ( ZSL%A_H(ZSL%TOT_NNZ_H) , ZSL%IA_H(MAX(ZSL%NUNKH_LOCAL,1)+1) , ZSL%JA_H(ZSL%TOT_NNZ_H) )
   ! Store upper triangular part of symmetric D_MAT_H in CSR format:
   IF(ZSL%NUNKH_LOCAL>0) THEN
      ! Here each process defines de beginning and end rows in global numeration, for the equations
      ! it has assembled:
      ZSL%LOWER_ROW = ZSL%UNKH_IND(NM_START) + 1
      ZSL%UPPER_ROW = ZSL%UNKH_IND(NM_START) + ZSL%NUNKH_LOCAL

      INNZ = 0
      DO IROW=1,ZSL%NUNKH_LOCAL
         ZSL%IA_H(IROW) = INNZ + 1
         DO JCOL=1,ZSL%NNZ_D_MAT_H(IROW)
            IF ( ZSL%JD_MAT_H(JCOL,IROW) < ZSL%UNKH_IND(NM_START)+IROW ) CYCLE ! Only upper Triangular part.
            INNZ = INNZ + 1
            ZSL%A_H(INNZ)  =  ZSL%D_MAT_H(JCOL,IROW)
            ZSL%JA_H(INNZ) = ZSL%JD_MAT_H(JCOL,IROW)
         ENDDO
      ENDDO
      ZSL%IA_H(ZSL%NUNKH_LOCAL+1) = INNZ + 1
   ELSE
      ZSL%LOWER_ROW = MAX(1,ZSL%UNKH_IND(NM_START))
      ZSL%UPPER_ROW = MAX(1,ZSL%UNKH_IND(NM_START))
      ZSL%A_H       = 0._EB     ! Add a zero coefficient in A(ZSL%UNKH_IND(NM_START),ZSL%UNKH_IND(NM_START)).
      ZSL%IA_H(1:2) = (/1,2/)
      ZSL%JA_H(1)   = MAX(1,ZSL%UNKH_IND(NM_START))
   ENDIF

   ! Define 4 byte A_H, F_H and X_H:
#ifdef SINGLE_PRECISION_PSN_SOLVE
   IPARM(28) = 1 ! Single Precision solve.
   ALLOCATE(ZSL%F_H_FB(1:MAX(ZSL%NUNKH_LOCAL,1))); ZSL%F_H_FB = 0._FB
   ALLOCATE(ZSL%X_H_FB(1:MAX(ZSL%NUNKH_LOCAL,1))); ZSL%X_H_FB = 0._FB
   ALLOCATE( ZSL%A_H_FB(ZSL%TOT_NNZ_H) );   ZSL%A_H_FB(1:ZSL%TOT_NNZ_H)   = REAL(ZSL%A_H(1:ZSL%TOT_NNZ_H),FB)
#endif

#ifdef WITH_MKL
   ! Lower and uppper rows handled by this process:
   IPARM(41) = ZSL%LOWER_ROW
   IPARM(42) = ZSL%UPPER_ROW

   ! Initialize solver pointer for H matrix solves:
   ALLOCATE(ZSL%PT_H(64))
   DO I=1,64
      ZSL%PT_H(I)%DUMMY = 0
   ENDDO

   ! Reorder and Symbolic factorization:
   PHASE = 11
#ifdef SINGLE_PRECISION_PSN_SOLVE
   CALL CLUSTER_SPARSE_SOLVER (ZSL%PT_H, MAXFCT, MNUM, ZSL%MTYPE, PHASE, ZSL%NUNKH_TOTAL, &
   ZSL%A_H_FB, ZSL%IA_H, ZSL%JA_H, PERM, NRHS, IPARM, MSGLVL, ZSL%F_H_FB, ZSL%X_H_FB, MPI_COMM_WORLD, ERROR)
#else
   CALL CLUSTER_SPARSE_SOLVER (ZSL%PT_H, MAXFCT, MNUM, ZSL%MTYPE, PHASE, ZSL%NUNKH_TOTAL, &
   ZSL%A_H, ZSL%IA_H, ZSL%JA_H, PERM, NRHS, IPARM, MSGLVL, ZSL%F_H, ZSL%X_H, MPI_COMM_WORLD, ERROR)
#endif

   IF (ERROR /= 0) THEN
      IF (MY_RANK==0) THEN
      WRITE(LU_ERR,'(A,2I5)') 'GET_H_MATRIX_LUDCMP CLUSTER_SOLVER Sym Factor: The following ERROR was detected: ', ERROR,IPZ
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
   MB_FACTOR(1,MY_RANK) = MAX(IPARM(15),IPARM(16)+IPARM(17))/1000
   MB_FACTOR(2,MY_RANK) = ZSL%NUNKH_LOCAL
   IF(N_MPI_PROCESSES > 1) &
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, MB_FACTOR(1,0), 2*N_MPI_PROCESSES, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
   ! Write to output file:
   IF(MY_RANK==0 .AND. GLMAT_VERBOSE) THEN
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
   CALL CLUSTER_SPARSE_SOLVER (ZSL%PT_H, MAXFCT, MNUM, ZSL%MTYPE, PHASE, ZSL%NUNKH_TOTAL, &
   ZSL%A_H_FB, ZSL%IA_H, ZSL%JA_H, PERM, NRHS, IPARM, MSGLVL, ZSL%F_H_FB, ZSL%X_H_FB, MPI_COMM_WORLD, ERROR)
#else
   CALL CLUSTER_SPARSE_SOLVER (ZSL%PT_H, MAXFCT, MNUM, ZSL%MTYPE, PHASE, ZSL%NUNKH_TOTAL, &
   ZSL%A_H, ZSL%IA_H, ZSL%JA_H, PERM, NRHS, IPARM, MSGLVL, ZSL%F_H, ZSL%X_H, MPI_COMM_WORLD, ERROR)
#endif

   IF (ERROR /= 0) THEN
      IF (MY_RANK==0) THEN
      WRITE(LU_ERR,'(A,2I5)') 'GET_H_MATRIX_LUDCMP CLUSTER_SOLVER Num Factor: The following ERROR was detected: ', ERROR,IPZ
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

   IF (MY_RANK==0) WRITE(LU_ERR,'(A)') &
   'Error: MKL Library compile flag was not defined for GLMAT as pressure solver.'
   ! Some error - stop flag for CALL STOP_CHECK(1).
   STOP_STATUS = SETUP_STOP
   RETURN

#endif

ENDDO IPZ_LOOP

! Set level MSG to 0 for solution:
IF(GLMAT_VERBOSE) MSGLVL = 0

END SUBROUTINE GET_H_MATRIX_LUDCMP

! -------------------------------- SET_CLUSTER_SOLVER_IPARM ---------------------------------
#ifdef WITH_MKL
SUBROUTINE SET_CLUSTER_SOLVER_IPARM

! Define control parameter vector iparm:
IF(ALLOCATED(IPARM)) DEALLOCATE(IPARM)
ALLOCATE(IPARM(64)); IPARM(:) = 0

IPARM(1) = 1   ! no solver default
! Pardiso: IPARM(2) = 2   ! fill-in reordering from METIS
IF (N_MPI_PROCESSES > 4) THEN ! Typical number of computing cores inside one chip.
   IPARM(2) =10   ! 10 = MPI Parallel fill-in reordering from METIS. If 3 = OpenMP parallel reordering in Master Node.
ELSE              ! Note IPARM(2)=10 has a bug which has been fixed from Intel MKL 2018 update 2 onwards.
   IPARM(2) = 3
ENDIF
IPARM(4) = 0   ! no iterative-direct algorithm
IPARM(5) = 0   ! no user fill-in reducing permutation
IPARM(6) = 0   ! =0 solution on the first n components of x
IPARM(8) = 2   ! numbers of iterative refinement steps
IPARM(10) = 13 ! perturb the pivot elements with 1E-13
IPARM(11) = 1  ! use nonsymmetric permutation and scaling MPS  !!!!! was 1
IPARM(13) = 1  ! maximum weighted matching algorithm is switched-off
               !(default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
IPARM(14) = 0  ! Output: number of perturbed pivots
IPARM(18) = 0  !-1 ! Output: number of nonzeros in the factor LU
IPARM(19) = 0  !-1 ! Output: Mflops for LU factorization
IPARM(20) = 0  ! Output: Numbers of CG Iterations

IPARM(21) = 1  ! 1x1 diagonal pivoting for symmetric indefinite matrices.

IPARM(24) = 0

IPARM(27) = 1  ! Check matrix

IPARM(40) = 2  ! Matrix, solution and rhs provided in distributed assembled matrix input format.

END SUBROUTINE SET_CLUSTER_SOLVER_IPARM
#endif

! -------------------------------- GET_BCS_H_MATRIX ---------------------------------

SUBROUTINE GET_BCS_H_MATRIX

USE MPI_F08
USE MESH_POINTERS
USE COMPLEX_GEOMETRY, ONLY : CC_IDRC
USE CC_SCALARS, ONLY : GET_CC_UNKH, GET_CFACE_OPEN_BC_COEF, GRADH_ON_CARTESIAN

! Local Variables:
INTEGER :: NM,NM1,JLOC,JCOL,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND),IERR,IIG,JJG,KKG,II,JJ,KK,IW,ILH,JLH,KLH,IRC
REAL(EB):: AF,IDX,BIJ
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
REAL(EB), POINTER, DIMENSION(:,:) :: D_MAT_HP
INTEGER :: H_MAT_IVEC

IPZ_LOOP : DO IPZ=0,N_ZONE
   ZSL => ZONE_SOLVE(IPZ)

   ! Allocate D_MAT_H:
   D_MAT_HP => ZSL%D_MAT_H
   NM1 = NM_START

   H_MAT_IVEC = 0
   ! Main Mesh Loop:
   MESH_LOOP_1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      WALL_LOOP_1 : DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC => WALL(IW)
         BC => BOUNDARY_COORD(WC%BC_INDEX)
         ! Only OPEN_BOUNDARY leads to a Dirichlet BC for H when we solve the problem on the whole
         ! unstructured domain. Everything else leads to Neuman BCs on H, no need to modify D_MAT_HP.
         IF ( WC%BOUNDARY_TYPE/=OPEN_BOUNDARY ) CYCLE
         IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG
         IF((.NOT.CC_IBM .AND. CCVAR(IIG,JJG,KKG,UNKH)<1) &
             .OR. ZONE_SOLVE(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
         II  = BC%II;  JJ  = BC%JJ;  KK  = BC%KK
         ! Unknowns on related cells:
         IF(CCVAR(IIG,JJG,KKG,CGSC)==IS_GASPHASE .OR. PRES_ON_WHOLE_DOMAIN) THEN
            IND(LOW_IND)  = CCVAR(IIG,JJG,KKG,UNKH)  ! internal cell.
         ELSEIF(CCVAR(IIG,JJG,KKG,CGSC)==IS_CUTCFE) THEN
            CALL GET_CC_UNKH(IIG,JJG,KKG,IND(LOW_IND))
         ELSE
            CYCLE ! Solid cell and .NOT.PRES_ON_WHOLE_DOMAIN
         ENDIF

         IND_LOC(LOW_IND) = IND(LOW_IND) - ZSL%UNKH_IND(NM1) ! All row indexes must refer to ind_loc.
         ILH           = 0; JLH = 0; KLH = 0
         SELECT CASE(BC%IOR)
         CASE( IAXIS)
            AF  = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG-1)) * DZ(KKG); ILH = -1; IDX= RDXN(IIG+ILH)
         CASE(-IAXIS)
            AF  = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG  )) * DZ(KKG);           IDX= RDXN(IIG+ILH)
         CASE( JAXIS)
            AF  = DX(IIG)*DZ(KKG); JLH= -1; IDX= RDYN(JJG+JLH)
         CASE(-JAXIS)
            AF  = DX(IIG)*DZ(KKG);          IDX= RDYN(JJG+JLH)
         CASE( KAXIS)
            AF  = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG); KLH = -1; IDX= RDZN(KKG+KLH)
         CASE(-KAXIS)
            AF  = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG);           IDX= RDZN(KKG+KLH)
         END SELECT
         IF (CC_IBM .AND. .NOT.GRADH_ON_CARTESIAN) THEN
            IRC = FCVAR(IIG+ILH,JJG+JLH,KKG+KLH,CC_IDRC,ABS(BC%IOR))
            IF(IRC > 0) IDX = 1._EB / ( RC_FACE(IRC)%XCEN(ABS(BC%IOR),HIGH_IND) - RC_FACE(IRC)%XCEN(ABS(BC%IOR),LOW_IND) )
         ENDIF
         ! Now add to Adiff corresponding coeff:
         BIJ = IDX*AF
         ! Case of unstructured projection:
         IF(WC%CUT_FACE_INDEX>0) CALL GET_CFACE_OPEN_BC_COEF(WC%CUT_FACE_INDEX,BOUNDARY_COORD(WC%BC_INDEX)%IOR,IDX,BIJ)

         ! Find diagonal column number:
         JCOL = -1
         DO JLOC = 1,ZSL%NNZ_D_MAT_H(IND_LOC(LOW_IND))
            IF (IND(LOW_IND) == ZSL%JD_MAT_H(JLOC,IND_LOC(LOW_IND))) THEN
               JCOL = JLOC
               EXIT
            ENDIF
         ENDDO
         ! Add diagonal coefficient due to DIRICHLET BC:
         D_MAT_HP(JCOL,IND_LOC(LOW_IND)) = D_MAT_HP(JCOL,IND_LOC(LOW_IND)) + 2._EB*BIJ
         ! Add to mesh dirichlet bc counter
         H_MAT_IVEC = H_MAT_IVEC + 1
      ENDDO WALL_LOOP_1
   ENDDO MESH_LOOP_1

   ! Is the resulting Matrix Indefinite?
   ! Here all reduce with sum among MPI processes:
   IF (N_MPI_PROCESSES > 1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,H_MAT_IVEC,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
   IF ( H_MAT_IVEC > 0 ) THEN
      ZSL%MTYPE = 2 ! At least one Dirichlet BC, matrix positive definite.
   ELSE
      ZSL%MTYPE =-2 ! Indefinite.
   ENDIF
ENDDO IPZ_LOOP

RETURN
END SUBROUTINE GET_BCS_H_MATRIX

! --------------------------------- GET_H_MATRIX ------------------------------------

SUBROUTINE GET_H_MATRIX

USE MESH_POINTERS
USE CC_SCALARS, ONLY : GET_H_MATRIX_CC

! Local Variables:
INTEGER :: NM,NM1,NREG
INTEGER :: LOW_FACE,HIGH_FACE,X1AXIS,X2AXIS,X3AXIS,IFACE
REAL(EB), POINTER, DIMENSION(:,:) :: D_MAT_HP
REAL(EB), POINTER, DIMENSION(:)   :: DX1,DX2,DX3
INTEGER :: I,J,K,I1,I2,I3,IIP,JJP,KKP,IIM,JJM,KKM
INTEGER :: ILOC,JLOC,IROW,JCOL,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND)
REAL(EB):: AF,IDX,BIJ,KFACE(1:2,1:2)
TYPE(CC_REGFACE_TYPE), POINTER, DIMENSION(:) :: RGF=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()
INTEGER :: IIG,JJG,KKG,II,JJ,KK,IOR,LOCROW,IW
INTEGER :: WC_JD(1:2,1:2)
LOGICAL :: FLG

IPZ_LOOP : DO IPZ=0,N_ZONE

   ZSL => ZONE_SOLVE(IPZ)

   ! Allocate D_MAT_H:
   ALLOCATE( ZSL%D_MAT_H(1:NNZ_ROW_H,1:ZSL%NUNKH_LOCAL) ); ZSL%D_MAT_H = 0._EB
   D_MAT_HP => ZSL%D_MAT_H
   NM1 = NM_START

   ! Main Mesh Loop:
   MESH_LOOP_1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      ! X direction bounds:
      ILO_FACE = 0             ! Low mesh boundary face index.
      IHI_FACE = IBAR          ! High mesh boundary face index.
      ILO_CELL = ILO_FACE + 1  ! First internal cell index. See notes.
      IHI_CELL = IHI_FACE      ! Last internal cell index.
      ! Y direction bounds:
      JLO_FACE = 0             ! Low mesh boundary face index.
      JHI_FACE = JBAR          ! High mesh boundary face index.
      JLO_CELL = JLO_FACE + 1  ! First internal cell index. See notes.
      JHI_CELL = JHI_FACE      ! Last internal cell index.
      ! Z direction bounds:
      KLO_FACE = 0             ! Low mesh boundary face index.
      KHI_FACE = KBAR          ! High mesh boundary face index.
      KLO_CELL = KLO_FACE + 1  ! First internal cell index. See notes.
      KHI_CELL = KHI_FACE      ! Last internal cell index.

      ! Regular Faces
      AXIS_LOOP_1 : DO X1AXIS=IAXIS,KAXIS
         NREG = MESHES(NM)%NREGFACE_H(X1AXIS)
         IIM = 0; JJM = 0; KKM = 0
         IIP = 0; JJP = 0; KKP = 0
         SELECT CASE(X1AXIS)
         CASE(IAXIS)
            RGF => REGFACE_IAXIS_H
            IIP = 1; LOW_FACE=ILO_FACE; HIGH_FACE=IHI_FACE
            X2AXIS=JAXIS; X3AXIS=KAXIS
            DX1 => DXN
            DX2 => DY
            DX3 => DZ
         CASE(JAXIS)
            RGF => REGFACE_JAXIS_H
            JJP = 1; LOW_FACE=JLO_FACE; HIGH_FACE=JHI_FACE
            X2AXIS=KAXIS; X3AXIS=IAXIS
            DX1 => DYN
            DX2 => DZ
            DX3 => DX
         CASE(KAXIS)
            RGF => REGFACE_KAXIS_H
            KKP =  1; LOW_FACE=KLO_FACE; HIGH_FACE=KHI_FACE
            X2AXIS=IAXIS; X3AXIS=JAXIS
            DX1 => DZN
            DX2 => DX
            DX3 => DY
         END SELECT

         IFACE_LOOP_1 : DO IFACE=1,NREG
            I  = RGF(IFACE)%IJK(IAXIS);  J  = RGF(IFACE)%IJK(JAXIS);  K  = RGF(IFACE)%IJK(KAXIS)
            IF(ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
            I1 = RGF(IFACE)%IJK(X1AXIS); I2 = RGF(IFACE)%IJK(X2AXIS); I3 = RGF(IFACE)%IJK(X3AXIS)
            ! Unknowns on related cells:
            IND(LOW_IND)     = CCVAR(I+IIM,J+JJM,K+KKM,UNKH)
            IND(HIGH_IND)    = CCVAR(I+IIP,J+JJP,K+KKP,UNKH)
            IND_LOC(LOW_IND) = IND(LOW_IND) - ZSL%UNKH_IND(NM1) ! All row indexes must refer to ind_loc.
            IND_LOC(HIGH_IND)= IND(HIGH_IND)- ZSL%UNKH_IND(NM1)
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
                  JCOL = RGF(IFACE)%JD(ILOC,JLOC); ! Local position of coef in D_MAT_H
                  ! Add coefficient:
                  D_MAT_HP(JCOL,IROW) = D_MAT_HP(JCOL,IROW) + KFACE(ILOC,JLOC)
               ENDDO
            ENDDO
         ENDDO IFACE_LOOP_1
         NULLIFY(RGF)
      ENDDO AXIS_LOOP_1

      ! Next, Wall faces of type INTERPOLATED_BOUNDARY or PERIODIC_BOUNDARY:
      ! Here We have to do something about WALL cells that are also cut-faces, who wins? Make cut-faces take precedence.
      LOCROW = LOW_IND
      WALL_LOOP_1 : DO IW=1,N_EXTERNAL_WALL_CELLS

         WC => WALL(IW)
         BC => BOUNDARY_COORD(WC%BC_INDEX)
         EWC=>EXTERNAL_WALL(IW)
         FLG = WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY .OR. WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY
         IF(PRES_ON_WHOLE_DOMAIN) &
         FLG = FLG .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. EWC%NOM > 0)
         IF ( .NOT.FLG .OR. EWC%NOM<1) CYCLE ! Here if NOM==0 means it is an OBST laying on an external boundary -> CYCLE
         IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG
         IF(ZONE_SOLVE(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE

         WC_JD(1,1) = WC%JD11_INDEX
         WC_JD(1,2) = WC%JD12_INDEX
         WC_JD(2,1) = WC%JD21_INDEX
         WC_JD(2,2) = WC%JD22_INDEX

         II = BC%II; JJ = BC%JJ; KK = BC%KK; IOR = BC%IOR
         ! Check if CC_IBM -> If either IIG,JJG,KKG or II,JJ,KK cell is type IS_CUTCFE or IS_SOLID cycle:
         IF ( .NOT.PRES_ON_WHOLE_DOMAIN .AND. CC_IBM ) THEN
            IF(CCVAR(II ,JJ ,KK ,CC_CGSC) /= IS_GASPHASE) CYCLE
            IF(CCVAR(IIG,JJG,KKG,CC_CGSC) /= IS_GASPHASE) CYCLE
         ENDIF

         ! Unknowns on related cells:
         IND(LOW_IND)     = CCVAR(IIG,JJG,KKG,UNKH)  ! internal.
         IND(HIGH_IND)    = CCVAR(II,JJ,KK,UNKH)     ! guard-cell.
         IND_LOC(LOW_IND) = IND(LOW_IND) - ZSL%UNKH_IND(NM1) ! All row indexes must refer to ind_loc.
         IND_LOC(HIGH_IND)= IND(HIGH_IND)- ZSL%UNKH_IND(NM1)
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

      ! Contribution to Laplacian matrix from RC and cut-faces:
      IF ( CC_IBM ) CALL GET_H_MATRIX_CC(NM,NM1,IPZ,D_MAT_HP)

   ENDDO MESH_LOOP_1
ENDDO IPZ_LOOP
RETURN
END SUBROUTINE GET_H_MATRIX


! --------------------------- GET_MATRIXGRAPH_H_WHLDOM ------------------------------

SUBROUTINE GET_MATRIXGRAPH_H_WHLDOM

USE MESH_POINTERS
USE CC_SCALARS, ONLY : GET_CC_MATRIXGRAPH_H, ADD_INPLACE_NNZ_H_WHLDOM
USE MPI_F08

! Local Variables:
INTEGER :: NM
INTEGER :: X1AXIS,IFACE,I,I1,J,K,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND)
INTEGER :: LOCROW,IIND,NII,ILOC
INTEGER :: NREG,IIM,JJM,KKM,IIP,JJP,KKP,LOW_FACE,HIGH_FACE,JLOC,IW,II,JJ,KK,IIG,JJG,KKG,ICF
LOGICAL :: INLIST
TYPE(CC_REGFACE_TYPE), POINTER, DIMENSION(:) :: RGF=>NULL()
TYPE(CC_RCFACE_TYPE), POINTER :: RCF=>NULL()
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()
INTEGER :: WC_JD(1:2,1:2)
LOGICAL :: FLG

! Write number of pressure unknowns to output:
IF (MY_RANK==0 .AND. GLMAT_VERBOSE) THEN
   WRITE(LU_OUTPUT,'(A)') '   Using GLMAT as pressure solver. List of H unknown numbers per proc:'
ENDIF

! Dump PRES_ZONE on each REG, RC and Cut-face of my MESHES:
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   ! Regular faces:
   DO X1AXIS=IAXIS,KAXIS
      SELECT CASE(X1AXIS)
      CASE(IAXIS); RGF=>MESHES(NM)%REGFACE_IAXIS_H
      CASE(JAXIS); RGF=>MESHES(NM)%REGFACE_JAXIS_H
      CASE(KAXIS); RGF=>MESHES(NM)%REGFACE_KAXIS_H
      END SELECT
      DO IFACE=1,MESHES(NM)%NREGFACE_H(X1AXIS)
         I = RGF(IFACE)%IJK(IAXIS); J = RGF(IFACE)%IJK(JAXIS); K = RGF(IFACE)%IJK(KAXIS)
         RGF(IFACE)%PRES_ZONE=ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT
      ENDDO
   ENDDO
   ! RC faces:
   DO IFACE=1,MESHES(NM)%CC_NRCFACE_H
      RCF => RC_FACE(MESHES(NM)%RCF_H(IFACE));
      I   = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS); K = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1)
      RCF%PRES_ZONE=ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT
   ENDDO
   ! Cut faces:
   DO ICF=1,MESHES(NM)%N_CUTFACE_MESH
      I = CUT_FACE(ICF)%IJK(IAXIS)
      J = CUT_FACE(ICF)%IJK(JAXIS)
      K = CUT_FACE(ICF)%IJK(KAXIS)
      CUT_FACE(ICF)%PRES_ZONE=ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT
   ENDDO
ENDDO

! Now Loop by Pressure Zone:
IPZ_LOOP : DO IPZ=0,N_ZONE

   ZSL => ZONE_SOLVE(IPZ)
   ZSL%NUNKH_LOCAL = SUM(ZSL%NUNKH_LOC(LOWER_MESH_INDEX:UPPER_MESH_INDEX))
   ! Allocate NNZ_D_MAT_H, JD_MAT_H:
   ALLOCATE( ZSL%NNZ_D_MAT_H(1:ZSL%NUNKH_LOCAL) );          ZSL%NNZ_D_MAT_H(:) = 0
   ALLOCATE( ZSL%JD_MAT_H(1:NNZ_ROW_H,1:ZSL%NUNKH_LOCAL) ); ZSL%JD_MAT_H(:,:)  = HUGE(I)

   ! Define NM_START: first mesh that belongs to the processor.
   NM_START = LOWER_MESH_INDEX

   ! First run over all regular and gasphase cut faces and insert-add lists of related unknows per unknown:
   MESH_LOOP_1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      ! X direction bounds:
      ILO_FACE = 0            ! Low mesh boundary face index.
      IHI_FACE = IBAR         ! High mesh boundary face index.
      ILO_CELL = ILO_FACE + 1 ! First internal cell index. See notes.
      IHI_CELL = IHI_FACE     ! Last internal cell index.
      ! Y direction bounds:
      JLO_FACE = 0            ! Low mesh boundary face index.
      JHI_FACE = JBAR         ! High mesh boundary face index.
      JLO_CELL = JLO_FACE + 1 ! First internal cell index. See notes.
      JHI_CELL = JHI_FACE     ! Last internal cell index.
      ! Z direction bounds:
      KLO_FACE = 0            ! Low mesh boundary face index.
      KHI_FACE = KBAR         ! High mesh boundary face index.
      KLO_CELL = KLO_FACE + 1 ! First internal cell index. See notes.
      KHI_CELL = KHI_FACE     ! Last internal cell index.

      ! Regular Faces
      AXIS_LOOP_1 : DO X1AXIS=IAXIS,KAXIS
         NREG = MESHES(NM)%NREGFACE_H(X1AXIS)
         IIM = 0; JJM = 0; KKM = 0
         IIP = 0; JJP = 0; KKP = 0
         SELECT CASE(X1AXIS)
         CASE(IAXIS)
            RGF => REGFACE_IAXIS_H
            IIP = 1; LOW_FACE=ILO_FACE; HIGH_FACE=IHI_FACE
         CASE(JAXIS)
            RGF => REGFACE_JAXIS_H
            JJP = 1; LOW_FACE=JLO_FACE; HIGH_FACE=JHI_FACE
         CASE(KAXIS)
            RGF => REGFACE_KAXIS_H
            KKP = 1; LOW_FACE=KLO_FACE; HIGH_FACE=KHI_FACE
         END SELECT

         IFACE_LOOP_1 : DO IFACE=1,NREG
            I = RGF(IFACE)%IJK(IAXIS); J = RGF(IFACE)%IJK(JAXIS); K = RGF(IFACE)%IJK(KAXIS); I1 = RGF(IFACE)%IJK(X1AXIS)
            IF(ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
            ! Unknowns on related cells:
            IND(LOW_IND)     = CCVAR(I+IIM,J+JJM,K+KKM,UNKH)
            IND(HIGH_IND)    = CCVAR(I+IIP,J+JJP,K+KKP,UNKH)
            IND_LOC(LOW_IND) = IND(LOW_IND) - ZSL%UNKH_IND(NM_START) ! All row indexes must refer to ind_loc.
            IND_LOC(HIGH_IND)= IND(HIGH_IND)- ZSL%UNKH_IND(NM_START)

            CALL ADD_INPLACE_NNZ_H_WHLDOM(LOW_IND,HIGH_IND,IND,IND_LOC,IPZ)
         ENDDO IFACE_LOOP_1
         NULLIFY(RGF)
      ENDDO AXIS_LOOP_1

      ! Next, Wall faces of type INTERPOLATED_BOUNDARY or PERIODIC_BOUNDARY:
      ! Here We have to do something about WALL cells that are also cut-faces, who wins? Make cut-faces take precedence.
      LOCROW = LOW_IND
      WALL_LOOP_1 : DO IW=1,N_EXTERNAL_WALL_CELLS
         WC => WALL(IW)
         BC => BOUNDARY_COORD(WC%BC_INDEX)
         EWC=>EXTERNAL_WALL(IW)
         FLG = WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY .OR. WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY
         IF(PRES_ON_WHOLE_DOMAIN) &
         FLG = FLG .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. EWC%NOM > 0)
         IF ( .NOT.FLG .OR. EWC%NOM<1) CYCLE ! Here if NOM==0 means it is an OBST laying on an external boundary -> CYCLE
         IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG; II = BC%II; JJ = BC%JJ; KK = BC%KK
         IF(ZONE_SOLVE(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
         ! Check if CC_IBM -> If either IIG,JJG,KKG or II,JJ,KK cell is type IS_CUTCFE or IS_SOLID cycle:
         IF (CC_IBM) THEN
            IF(CCVAR(II ,JJ ,KK ,CC_CGSC) /= IS_GASPHASE) CYCLE
            IF(CCVAR(IIG,JJG,KKG,CC_CGSC) /= IS_GASPHASE) CYCLE
         ENDIF
         ! Unknowns on related cells:
         IND(LOW_IND)     = CCVAR(IIG,JJG,KKG,UNKH)  ! internal.
         IND(HIGH_IND)    = CCVAR(II,JJ,KK,UNKH)     ! guard-cell.
         IND_LOC(LOW_IND) = IND(LOW_IND) - ZSL%UNKH_IND(NM_START) ! All row indexes must refer to ind_loc.
         IND_LOC(HIGH_IND)= IND(HIGH_IND)- ZSL%UNKH_IND(NM_START)
         ! Same code ADD_INPLACE_NNZ_H_WHLDOM(LOCROW,LOCROW,IND,IND_LOC):
         DO IIND=LOW_IND,HIGH_IND
            NII = ZSL%NNZ_D_MAT_H(IND_LOC(LOCROW))
            ! Check that column index hasn't been already counted:
            INLIST = .FALSE.
            DO ILOC=1,NII
               IF ( IND(IIND) == ZSL%JD_MAT_H(ILOC,IND_LOC(LOCROW)) ) THEN
                  INLIST = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (INLIST) CYCLE
            ! Now add in place:
            NII = NII + 1
            DO ILOC=1,NII
               IF ( ZSL%JD_MAT_H(ILOC,IND_LOC(LOCROW)) > IND(IIND) ) EXIT
            ENDDO
            DO JLOC=NII,ILOC+1,-1
               ZSL%JD_MAT_H(JLOC,IND_LOC(LOCROW)) = ZSL%JD_MAT_H(JLOC-1,IND_LOC(LOCROW))
            ENDDO
            ZSL%NNZ_D_MAT_H(IND_LOC(LOCROW))   = NII
            ZSL%JD_MAT_H(ILOC,IND_LOC(LOCROW)) = IND(IIND)
         ENDDO
      ENDDO WALL_LOOP_1

      ! Finally Add nonzeros corresponding to RC_FACE, CUT_FACE
      IF (CC_IBM) CALL GET_CC_MATRIXGRAPH_H(NM,NM_START,IPZ,.TRUE.)

   ENDDO MESH_LOOP_1

   ! Now add local column location to Faces data structures:
   MESH_LOOP_2 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      ! X direction bounds:
      ILO_FACE = 0            ! Low mesh boundary face index.
      IHI_FACE = IBAR         ! High mesh boundary face index.
      ILO_CELL = ILO_FACE + 1 ! First internal cell index. See notes.
      IHI_CELL = IHI_FACE     ! Last internal cell index.
      ! Y direction bounds:
      JLO_FACE = 0            ! Low mesh boundary face index.
      JHI_FACE = JBAR         ! High mesh boundary face index.
      JLO_CELL = JLO_FACE + 1 ! First internal cell index. See notes.
      JHI_CELL = JHI_FACE     ! Last internal cell index.
      ! Z direction bounds:
      KLO_FACE = 0            ! Low mesh boundary face index.
      KHI_FACE = KBAR         ! High mesh boundary face index.
      KLO_CELL = KLO_FACE + 1 ! First internal cell index. See notes.
      KHI_CELL = KHI_FACE     ! Last internal cell index.
      ! Regular Faces, loop is similar to before:
      AXIS_LOOP_2 : DO X1AXIS=IAXIS,KAXIS
         NREG = MESHES(NM)%NREGFACE_H(X1AXIS)
         IIM = 0; JJM = 0; KKM = 0
         IIP = 0; JJP = 0; KKP = 0
         SELECT CASE(X1AXIS)
         CASE(IAXIS)
            RGF => REGFACE_IAXIS_H
            IIP = 1; LOW_FACE=ILO_FACE; HIGH_FACE=IHI_FACE
         CASE(JAXIS)
            RGF => REGFACE_JAXIS_H
            JJP = 1; LOW_FACE=JLO_FACE; HIGH_FACE=JHI_FACE
         CASE(KAXIS)
            RGF => REGFACE_KAXIS_H
            KKP = 1; LOW_FACE=KLO_FACE; HIGH_FACE=KHI_FACE
         END SELECT

         IFACE_LOOP_2 : DO IFACE=1,NREG
            I = RGF(IFACE)%IJK(IAXIS); J = RGF(IFACE)%IJK(JAXIS); K = RGF(IFACE)%IJK(KAXIS); I1 = RGF(IFACE)%IJK(X1AXIS)
            IF(ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
            ! Unknowns on related cells:
            IND(LOW_IND)     = CCVAR(I+IIM,J+JJM,K+KKM,UNKH)
            IND(HIGH_IND)    = CCVAR(I+IIP,J+JJP,K+KKP,UNKH)
            IND_LOC(LOW_IND) = IND(LOW_IND) - ZSL%UNKH_IND(NM_START) ! All row indexes must refer to ind_loc.
            IND_LOC(HIGH_IND)= IND(HIGH_IND)- ZSL%UNKH_IND(NM_START)

            RGF(IFACE)%JD(1:2,1:2) = IS_UNDEFINED
            DO LOCROW = LOW_IND,HIGH_IND
               DO IIND=LOW_IND,HIGH_IND
                  NII = ZSL%NNZ_D_MAT_H(IND_LOC(LOCROW))
                  DO ILOC=1,NII
                     IF ( IND(IIND) == ZSL%JD_MAT_H(ILOC,IND_LOC(LOCROW)) ) THEN
                        RGF(IFACE)%JD(LOCROW,IIND) = ILOC;
                        EXIT
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO IFACE_LOOP_2
         NULLIFY(RGF)
      ENDDO AXIS_LOOP_2

      ! Now Wall faces column locations:
      LOCROW = LOW_IND
      WALL_LOOP_2 : DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC => WALL(IW)
         BC => BOUNDARY_COORD(WC%BC_INDEX)
         WC_JD(1:2,1:2) = IS_UNDEFINED
         IF (.NOT.PRES_ON_WHOLE_DOMAIN .AND. WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE
         IF(IW <= N_EXTERNAL_WALL_CELLS) THEN
            ! Here if NOM==0 means it is an OBST laying on an external boundary -> CYCLE
            EWC=>EXTERNAL_WALL(IW); IF(EWC%NOM < 1) CYCLE
         ENDIF
         IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG; II = BC%II; JJ = BC%JJ; KK = BC%KK
         IF((.NOT.CC_IBM .AND. CCVAR(IIG,JJG,KKG,UNKH)<1) &
            .OR. ZONE_SOLVE(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
         ! Check if CC_IBM -> If either IIG,JJG,KKG or II,JJ,KK cell is type IS_CUTCFE or IS_SOLID cycle:
         IF (CC_IBM) THEN
            IF(CCVAR(II ,JJ ,KK ,CC_CGSC) /= IS_GASPHASE) CYCLE
            IF(CCVAR(IIG,JJG,KKG,CC_CGSC) /= IS_GASPHASE) CYCLE
         ENDIF
         ! Unknowns on related cells:
         IND(LOW_IND)     = CCVAR(IIG,JJG,KKG,UNKH)  ! internal.
         IND(HIGH_IND)    = CCVAR(II,JJ,KK,UNKH)     ! guard-cell.
         IND_LOC(LOW_IND) = IND(LOW_IND) - ZSL%UNKH_IND(NM_START) ! All row indexes must refer to ind_loc.
         IND_LOC(HIGH_IND)= IND(HIGH_IND)- ZSL%UNKH_IND(NM_START)
         DO IIND=LOW_IND,HIGH_IND
            NII = ZSL%NNZ_D_MAT_H(IND_LOC(LOCROW))
            DO ILOC=1,NII
               IF ( IND(IIND) == ZSL%JD_MAT_H(ILOC,IND_LOC(LOCROW)) ) THEN
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
      IF (CC_IBM) CALL GET_CC_MATRIXGRAPH_H(NM,NM_START,IPZ,.FALSE.)
   ENDDO MESH_LOOP_2

ENDDO IPZ_LOOP

RETURN
END SUBROUTINE GET_MATRIXGRAPH_H_WHLDOM

! ---------------------------- GET_H_REGFACES -----------------------------------

SUBROUTINE GET_H_REGFACES

USE MESH_POINTERS
USE CC_SCALARS, ONLY : GET_RCFACES_H

! Local Variables:
INTEGER :: NM
INTEGER :: ILO,IHI,JLO,JHI,KLO,KHI
INTEGER :: I,J,K,II,IREG,X1AXIS
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IJKBUFFER
INTEGER :: IW, IIG, JJG, KKG, IOR, N_INTERNAL_WALL_CELLS_AUX
LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: LOG_INTWC
TYPE(WALL_TYPE), POINTER :: WC=>NULL()
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC

! Mesh loop:
MAIN_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)

   ! X direction bounds:
   ILO_FACE = 0            ! Low mesh boundary face index.
   IHI_FACE = IBAR         ! High mesh boundary face index.
   ILO_CELL = ILO_FACE + 1 ! First internal cell index. See notes.
   IHI_CELL = IHI_FACE     ! Last internal cell index.

   ! Y direction bounds:
   JLO_FACE = 0            ! Low mesh boundary face index.
   JHI_FACE = JBAR         ! High mesh boundary face index.
   JLO_CELL = JLO_FACE + 1 ! First internal cell index. See notes.
   JHI_CELL = JHI_FACE     ! Last internal cell index.

   ! Z direction bounds:
   KLO_FACE = 0            ! Low mesh boundary face index.
   KHI_FACE = KBAR         ! High mesh boundary face index.
   KLO_CELL = KLO_FACE + 1 ! First internal cell index. See notes.
   KHI_CELL = KHI_FACE     ! Last internal cell index.

   ! Set starting number of regular faces for NM to zero:
   MESHES(NM)%NREGFACE_H(IAXIS:KAXIS) = 0

   ! 1. Regular GASPHASE faces connected to Gasphase cells:
   ALLOCATE(IJKBUFFER(IAXIS:KAXIS,1:(IBAR+1)*(JBAR+1)*(KBAR+1)))
   ! Check internal SOLID_BOUNDARY faces:
   ALLOCATE(LOG_INTWC(ILO_FACE:IHI_FACE,JLO_FACE:JHI_FACE,KLO_FACE:KHI_FACE,IAXIS:KAXIS)); LOG_INTWC(:,:,:,:) = .FALSE.

   N_INTERNAL_WALL_CELLS_AUX = 0
   IF(.NOT. PRES_ON_WHOLE_DOMAIN) N_INTERNAL_WALL_CELLS_AUX = N_INTERNAL_WALL_CELLS
   WALL_LOOP_1 : DO IW=N_EXTERNAL_WALL_CELLS+1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS_AUX
      WC => WALL(IW)
      BC => BOUNDARY_COORD(WC%BC_INDEX)
      IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE
      IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG

      IOR = BC%IOR
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
            IF (LOG_INTWC(I,J,K,X1AXIS))  CYCLE
            IF (CCVAR(I  ,J,K,UNKH) <= 0) CYCLE
            IF (CCVAR(I+1,J,K,UNKH) <= 0) CYCLE
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
            IF (LOG_INTWC(I,J,K,X1AXIS))  CYCLE
            IF (CCVAR(I,J  ,K,UNKH) <= 0) CYCLE
            IF (CCVAR(I,J+1,K,UNKH) <= 0) CYCLE
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
            IF (LOG_INTWC(I,J,K,X1AXIS))  CYCLE
            IF (CCVAR(I,J,K  ,UNKH) <= 0) CYCLE
            IF (CCVAR(I,J,K+1,UNKH) <= 0) CYCLE
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

USE MESH_POINTERS
USE CC_SCALARS, ONLY : NUMBER_UNKH_CUTCELLS
USE MPI_F08

! Local Variables:
INTEGER :: NM, IOPZ
INTEGER :: I,J,K,IERR
INTEGER, ALLOCATABLE, DIMENSION(:) :: NUNKH_TOT

IF(ALLOCATED(ZONE_SOLVE)) DEALLOCATE(ZONE_SOLVE)
ALLOCATE(ZONE_SOLVE(0:N_ZONE))

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   ! Initialize unknown numbers for H:
   MESHES(NM)%CCVAR(:,:,:,UNKH) = IS_UNDEFINED
   ! Identify connected ZONEs
   ! translate logical array to local and fill diagonal
   CONNECTED_ZONES_LOC = 0
   DO IOPZ=0,N_ZONE
      DO IPZ=0,N_ZONE
         IF (IPZ==IOPZ .OR. CONNECTED_ZONES(IPZ,IOPZ,NM)) CONNECTED_ZONES_LOC(IPZ,IOPZ) = 1
      ENDDO
   ENDDO
   ! establish sets of connected zones
   DO IPZ=1,N_ZONE
      CONNECTED_ZONES_LOC = MATMUL(CONNECTED_ZONES_LOC,CONNECTED_ZONES_LOC)
   ENDDO
   CONNECTED_ZONES_LOC = MIN(1,CONNECTED_ZONES_LOC)
   ! select the parent zone as the first in the row
   DO IPZ=0,N_ZONE
      ZSL => ZONE_SOLVE(IPZ)
      IF(.NOT.PRES_ON_WHOLE_DOMAIN) &
      ZSL%CONNECTED_ZONE_PARENT = MINLOC(CONNECTED_ZONES_LOC(IPZ,:), DIM=1, MASK = CONNECTED_ZONES_LOC(IPZ,:)/=0) - 1
   ENDDO
ENDDO

PRES_ZONE_LOOP : DO IPZ=0,N_ZONE

   ZSL => ZONE_SOLVE(IPZ)

   ! Define local number of H unknowns for this pressure zone:
   IF (ALLOCATED(ZSL%NUNKH_LOC)) DEALLOCATE(ZSL%NUNKH_LOC)
   ALLOCATE(ZSL%NUNKH_LOC(LOWER_MESH_INDEX:UPPER_MESH_INDEX)); ZSL%NUNKH_LOC = 0

   IF (ALLOCATED(ZSL%UNKH_IND)) DEALLOCATE(ZSL%UNKH_IND)
   ALLOCATE(ZSL%UNKH_IND(LOWER_MESH_INDEX:UPPER_MESH_INDEX)); ZSL%UNKH_IND = 0

   ! Unknown numbers for Poisson equation:
   MAIN_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      IF_PRES_ON_WHOLE_DOMAIN : IF(PRES_ON_WHOLE_DOMAIN) THEN ! Classic IBM.
         ! GLMAT : Loop on Cartesian cells, define cut cells and solid cells CGSC:
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  ZSL%NUNKH_LOC(NM) = ZSL%NUNKH_LOC(NM) + 1
                  CCVAR(I,J,K,UNKH) = ZSL%NUNKH_LOC(NM)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         ! UGLMAT : Loop on Cartesian cells, define cut cells and solid cells CGSC:
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF ( CCVAR(I,J,K,CGSC) /= IS_GASPHASE ) CYCLE
                  IF ( ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT /= IPZ ) CYCLE
                  ZSL%NUNKH_LOC(NM) = ZSL%NUNKH_LOC(NM) + 1
                  CCVAR(I,J,K,UNKH) = ZSL%NUNKH_LOC(NM)
               ENDDO
            ENDDO
         ENDDO
         ! Number MESH local unknowns on cut-cells:
         IF(CC_IBM) CALL NUMBER_UNKH_CUTCELLS(.TRUE.,NM,IPZ,ZSL%NUNKH_LOC)
      ENDIF IF_PRES_ON_WHOLE_DOMAIN
   ENDDO MAIN_MESH_LOOP

   ! Define total number of unknowns and global unknow index start per MESH:
   ALLOCATE(NUNKH_TOT(1:NMESHES)); NUNKH_TOT = 0
   NUNKH_TOT(LOWER_MESH_INDEX:UPPER_MESH_INDEX) = ZSL%NUNKH_LOC(LOWER_MESH_INDEX:UPPER_MESH_INDEX)
   IF(N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(MPI_IN_PLACE, NUNKH_TOT(1), NMESHES, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
   ! Define global start indexes for each mesh:
   I = 0
   DO NM=2,NMESHES
      I = I + NUNKH_TOT(NM-1); IF(NM<LOWER_MESH_INDEX .OR. NM>UPPER_MESH_INDEX) CYCLE
      ZSL%UNKH_IND(NM) = I
   ENDDO
   ZSL%NUNKH_TOTAL=SUM(NUNKH_TOT(1:NMESHES)); DEALLOCATE(NUNKH_TOT)

   ! Add initial index UNKX_ind to mesh blocks (regular + cut-cells):
   MAIN_MESH_LOOP2 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      IF_PRES_ON_WHOLE_DOMAIN2 : IF(PRES_ON_WHOLE_DOMAIN) THEN ! Classic IBM.
         ! Loop on all Cartesian cells:
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  CCVAR(I,J,K,UNKH) = CCVAR(I,J,K,UNKH) + ZSL%UNKH_IND(NM)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         ! Loop on Cartesian cells, GASPHASE:
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF ( CCVAR(I,J,K,CGSC) /= IS_GASPHASE ) CYCLE
                  IF ( ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT /= IPZ ) CYCLE
                  CCVAR(I,J,K,UNKH) = CCVAR(I,J,K,UNKH) + ZSL%UNKH_IND(NM)
               ENDDO
            ENDDO
         ENDDO
         ! Number MESH global unknowns on cut-cells (fully unstructured) or their under underlaying Cartesian cells
         ! (Cartesian unstructured):
         IF(CC_IBM) CALL NUMBER_UNKH_CUTCELLS(.FALSE.,NM,IPZ,ZSL%NUNKH_LOC)
      ENDIF IF_PRES_ON_WHOLE_DOMAIN2
   ENDDO MAIN_MESH_LOOP2

ENDDO PRES_ZONE_LOOP

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

      ! X direction bounds:
      ILO_FACE = 0                 ! Low mesh boundary face index.
      IHI_FACE = MESHES(NM)%IBAR   ! High mesh boundary face index.
      ILO_CELL = ILO_FACE + 1      ! First internal cell index. See notes.
      IHI_CELL = IHI_FACE          ! Last internal cell index.
      ISTR     = ILO_FACE - NGUARD ! Allocation start x arrays.
      IEND     = IHI_FACE + NGUARD ! Allocation end x arrays.

      ! Y direction bounds:
      JLO_FACE = 0                 ! Low mesh boundary face index.
      JHI_FACE = MESHES(NM)%JBAR   ! High mesh boundary face index.
      JLO_CELL = JLO_FACE + 1      ! First internal cell index. See notes.
      JHI_CELL = JHI_FACE          ! Last internal cell index.
      JSTR     = JLO_FACE - NGUARD ! Allocation start y arrays.
      JEND     = JHI_FACE + NGUARD ! Allocation end y arrays.

      ! Z direction bounds:
      KLO_FACE = 0                 ! Low mesh boundary face index.
      KHI_FACE = MESHES(NM)%KBAR   ! High mesh boundary face index.
      KLO_CELL = KLO_FACE + 1      ! First internal cell index. See notes.
      KHI_CELL = KHI_FACE          ! Last internal cell index.
      KSTR     = KLO_FACE - NGUARD ! Allocation start z arrays.
      KEND     = KHI_FACE + NGUARD ! Allocation end z arrays.

      ! Cartesian cells:
      IF (.NOT. ALLOCATED(MESHES(NM)%CCVAR)) ALLOCATE(MESHES(NM)%CCVAR(ISTR:IEND,JSTR:JEND,KSTR:KEND,NCVARS))
      MESHES(NM)%CCVAR = 0; MESHES(NM)%CCVAR(:,:,:,CGSC) = IS_UNDEFINED
      MESHES(NM)%CCVAR(ILO_CELL:IHI_CELL,JLO_CELL:JHI_CELL,KLO_CELL:KHI_CELL,CGSC) = IS_GASPHASE ! Set all internal
                                                                                     ! Block cells to GASPHASE.
   ENDDO FIRST_MESH_LOOP
ENDIF

! At this point, if CC_IBM is being used we have SOLID and CUTCELLS from GEOM objects.
! Now add OBST SOLID cells to CCVAR CGSC array:
SECND_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CELL(CELL_INDEX(I,J,K))%SOLID) CCVAR(I,J,K,CGSC) = IS_SOLID
         ENDDO
      ENDDO
   ENDDO
ENDDO SECND_MESH_LOOP

RETURN
END SUBROUTINE SET_CCVAR_CGSC_H


! --------------------- PRESSURE_SOLVER_CHECK_RESIDUALS_U ----------------------

SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS_U(NM)

USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE GLOBAL_CONSTANTS
USE PRES, ONLY : PRESSURE_SOLVER_CHECK_RESIDUALS
USE CC_SCALARS, ONLY : UNSTRUCTURED_POISSON_RESIDUAL, UNSTRUCTURED_POISSON_RESIDUAL_RC, &
                       COMPUTE_LINKED_CUTFACE_BAROCLINIC


INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,RHOP,P,RESIDUAL
INTEGER :: I,J,K,IPZ,IC
REAL(EB) :: LHSS,RHSS,IMFCT,JMFCT,KMFCT,IPFCT,JPFCT,KPFCT

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

TNOW=CURRENT_TIME()
IF (PRES_FLAG==ULMAT_FLAG) THEN
   DO IPZ=0,N_ZONE
      IF (MESHES(NM)%ZONE_MESH(IPZ)%USE_FFT) THEN
         CALL PRESSURE_SOLVER_CHECK_RESIDUALS(NM)
         T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
         RETURN
      ENDIF
   ENDDO
ENDIF

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
            IF(CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
            IMFCT = 1._EB; JMFCT = 1._EB; KMFCT = 1._EB; IPFCT = 1._EB; JPFCT = 1._EB; KPFCT = 1._EB
            IF (CC_IBM) THEN
               IF(CCVAR(I,J,K,CC_CGSC)==IS_SOLID) CYCLE

               IF(CCVAR(I,J,K,CC_CGSC)==IS_CUTCFE) THEN
                  CALL UNSTRUCTURED_POISSON_RESIDUAL(I,J,K,HP,RHOP,P,RESIDUAL(I,J,K),DO_SEPARABLE=.TRUE.);    CYCLE
               ELSEIF(ANY(CCVAR((/I-1,I+1/),J,K,CC_CGSC)==IS_CUTCFE) .OR. &
                      ANY(CCVAR(I,(/J-1,J+1/),K,CC_CGSC)==IS_CUTCFE) .OR. &
                      ANY(CCVAR(I,J,(/K-1,K+1/),CC_CGSC)==IS_CUTCFE)) THEN
                  CALL UNSTRUCTURED_POISSON_RESIDUAL_RC(I,J,K,HP,RHOP,P,RESIDUAL(I,J,K),DO_SEPARABLE=.TRUE.); CYCLE
               ENDIF

               IF(FCVAR(I-1,J  ,K  ,CC_FGSC,IAXIS)==IS_SOLID) IMFCT = 0._EB
               IF(FCVAR(I  ,J  ,K  ,CC_FGSC,IAXIS)==IS_SOLID) IPFCT = 0._EB
               IF(FCVAR(I  ,J-1,K  ,CC_FGSC,JAXIS)==IS_SOLID) JMFCT = 0._EB
               IF(FCVAR(I  ,J  ,K  ,CC_FGSC,JAXIS)==IS_SOLID) JPFCT = 0._EB
               IF(FCVAR(I  ,J  ,K-1,CC_FGSC,KAXIS)==IS_SOLID) KMFCT = 0._EB
               IF(FCVAR(I  ,J  ,K  ,CC_FGSC,KAXIS)==IS_SOLID) KPFCT = 0._EB
            ENDIF
            ! If surrounding wall_cell is type SOLID_BOUNDARY set FCT gradient factor to zero:
            IC = CELL_INDEX(I,J,K)
            IF (WALL(CELL(IC)%WALL_INDEX(-1))%BOUNDARY_TYPE==SOLID_BOUNDARY) IMFCT = 0._EB
            IF (WALL(CELL(IC)%WALL_INDEX( 1))%BOUNDARY_TYPE==SOLID_BOUNDARY) IPFCT = 0._EB
            IF (WALL(CELL(IC)%WALL_INDEX(-2))%BOUNDARY_TYPE==SOLID_BOUNDARY) JMFCT = 0._EB
            IF (WALL(CELL(IC)%WALL_INDEX( 2))%BOUNDARY_TYPE==SOLID_BOUNDARY) JPFCT = 0._EB
            IF (WALL(CELL(IC)%WALL_INDEX(-3))%BOUNDARY_TYPE==SOLID_BOUNDARY) KMFCT = 0._EB
            IF (WALL(CELL(IC)%WALL_INDEX( 3))%BOUNDARY_TYPE==SOLID_BOUNDARY) KPFCT = 0._EB
            RHSS = ( R(I-1)*FVX(I-1,J,K) - R(I)*FVX(I,J,K) )*RDX(I)*RRN(I) &
                 + (        FVY(I,J-1,K) -      FVY(I,J,K) )*RDY(J)        &
                 + (        FVZ(I,J,K-1) -      FVZ(I,J,K) )*RDZ(K)        &
                 - DDDT(I,J,K)
            LHSS = ((HP(I+1,J,K)-HP(I,J,K))*RDXN(I)  *R(I)  *IPFCT -                                                &
                    (HP(I,J,K)-HP(I-1,J,K))*RDXN(I-1)*R(I-1)*IMFCT                                  )*RDX(I)*RRN(I) &
                 + ((HP(I,J+1,K)-HP(I,J,K))*RDYN(J)*JPFCT-(HP(I,J,K)-HP(I,J-1,K))*RDYN(J-1)*JMFCT   )*RDY(J)        &
                 + ((HP(I,J,K+1)-HP(I,J,K))*RDZN(K)*KPFCT-(HP(I,J,K)-HP(I,J,K-1))*RDZN(K-1)*KMFCT   )*RDZ(K)
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
   IF(CC_IBM) CALL COMPUTE_LINKED_CUTFACE_BAROCLINIC(NM,HP,RHOP,P)
   !$OMP PARALLEL PRIVATE(I,J,K,RHSS,LHSS,IMFCT,JMFCT,KMFCT,IPFCT,JPFCT,KPFCT)
   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF(CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
            IMFCT = 1._EB; JMFCT = 1._EB; KMFCT = 1._EB; IPFCT = 1._EB; JPFCT = 1._EB; KPFCT = 1._EB
            IF (CC_IBM) THEN
               IF(CCVAR(I,J,K,CC_CGSC)==IS_SOLID) CYCLE

               IF(CCVAR(I,J,K,CC_CGSC)==IS_CUTCFE) THEN
               CALL UNSTRUCTURED_POISSON_RESIDUAL(I,J,K,HP,RHOP,P,RESIDUAL(I,J,K),DO_SEPARABLE=.FALSE.);    CYCLE
               ELSEIF (ANY(CCVAR((/I-1,I+1/),J,K,CC_CGSC)==IS_CUTCFE) .OR. &
                       ANY(CCVAR(I,(/J-1,J+1/),K,CC_CGSC)==IS_CUTCFE) .OR. &
                       ANY(CCVAR(I,J,(/K-1,K+1/),CC_CGSC)==IS_CUTCFE)) THEN
               CALL UNSTRUCTURED_POISSON_RESIDUAL_RC(I,J,K,HP,RHOP,P,RESIDUAL(I,J,K),DO_SEPARABLE=.FALSE.); CYCLE
               ENDIF

               IF(FCVAR(I-1,J  ,K  ,CC_FGSC,IAXIS)==IS_SOLID) IMFCT = 0._EB
               IF(FCVAR(I  ,J  ,K  ,CC_FGSC,IAXIS)==IS_SOLID) IPFCT = 0._EB
               IF(FCVAR(I  ,J-1,K  ,CC_FGSC,JAXIS)==IS_SOLID) JMFCT = 0._EB
               IF(FCVAR(I  ,J  ,K  ,CC_FGSC,JAXIS)==IS_SOLID) JPFCT = 0._EB
               IF(FCVAR(I  ,J  ,K-1,CC_FGSC,KAXIS)==IS_SOLID) KMFCT = 0._EB
               IF(FCVAR(I  ,J  ,K  ,CC_FGSC,KAXIS)==IS_SOLID) KPFCT = 0._EB
            ENDIF
            ! If surrounding wall_cell is type SOLID_BOUNDARY set FCT gradient factor to zero:
            IMFCT = GRADIENT_WEIGHT(NM, I, J, K, -1)
            IPFCT = GRADIENT_WEIGHT(NM, I, J, K,  1)
            JMFCT = GRADIENT_WEIGHT(NM, I, J, K, -2)
            JPFCT = GRADIENT_WEIGHT(NM, I, J, K,  2)
            KMFCT = GRADIENT_WEIGHT(NM, I, J, K, -3)
            KPFCT = GRADIENT_WEIGHT(NM, I, J, K,  3)
            RHSS =(R(I-1)*(FVX(I-1,J,K)-FVX_B(I-1,J,K)*IMFCT) - R(I)*(FVX(I,J,K)-FVX_B(I,J,K)*IPFCT) )*RDX(I)*RRN(I)&
                 +(       (FVY(I,J-1,K)-FVY_B(I,J-1,K)*JMFCT) -      (FVY(I,J,K)-FVY_B(I,J,K)*JPFCT) )*RDY(J)       &
                 +(       (FVZ(I,J,K-1)-FVZ_B(I,J,K-1)*KMFCT) -      (FVZ(I,J,K)-FVZ_B(I,J,K)*KPFCT) )*RDZ(K)       &
                 -DDDT(I,J,K)
            LHSS = ((P(I+1,J,K)-P(I,J,K))*RDXN(I)*R(I)    *2._EB/(RHOP(I+1,J,K)+RHOP(I,J,K))*IPFCT -                &
                    (P(I,J,K)-P(I-1,J,K))*RDXN(I-1)*R(I-1)*2._EB/(RHOP(I-1,J,K)+RHOP(I,J,K))*IMFCT)*RDX(I)*RRN(I)   &
                 + ((P(I,J+1,K)-P(I,J,K))*RDYN(J)         *2._EB/(RHOP(I,J+1,K)+RHOP(I,J,K))*JPFCT -                &
                    (P(I,J,K)-P(I,J-1,K))*RDYN(J-1)       *2._EB/(RHOP(I,J-1,K)+RHOP(I,J,K))*JMFCT)*RDY(J)          &
                 + ((P(I,J,K+1)-P(I,J,K))*RDZN(K)         *2._EB/(RHOP(I,J,K+1)+RHOP(I,J,K))*KPFCT -                &
                    (P(I,J,K)-P(I,J,K-1))*RDZN(K-1)       *2._EB/(RHOP(I,J,K-1)+RHOP(I,J,K))*KMFCT)*RDZ(K)          &
            + ((KRES(I+1,J,K)-KRES(I,J,K))*RDXN(I)  *R(I)  *IPFCT -                                                 &
               (KRES(I,J,K)-KRES(I-1,J,K))*RDXN(I-1)*R(I-1)*IMFCT )*RDX(I)*RRN(I)                                   &
            + ((KRES(I,J+1,K)-KRES(I,J,K))*RDYN(J)*JPFCT      -(KRES(I,J,K)-KRES(I,J-1,K))*RDYN(J-1)*JMFCT)*RDY(J)  &
            + ((KRES(I,J,K+1)-KRES(I,J,K))*RDZN(K)*KPFCT      -(KRES(I,J,K)-KRES(I,J,K-1))*RDZN(K-1)*KMFCT)*RDZ(K)
            RESIDUAL(I,J,K) = ABS(RHSS-LHSS)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
   !$OMP END PARALLEL
   PRESSURE_ERROR_MAX(NM) = MAXVAL(RESIDUAL)
   PRESSURE_ERROR_MAX_LOC(:,NM) = MAXLOC(RESIDUAL)
   IF (STORE_PRESSURE_POISSON_RESIDUAL) PP_RESIDUAL(1:IBAR,1:JBAR,1:KBAR)=RESIDUAL(1:IBAR,1:JBAR,1:KBAR)
ENDIF

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
END SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS_U

! Determine if the gradient towards a solid cell must be considered or not
! This is only the case if the solid cell is an internal cell or adjacent to a mesh interface

REAL(EB) FUNCTION GRADIENT_WEIGHT(NM, I,J,K, IOR0)
INTEGER, INTENT(IN) :: NM, I, J, K, IOR0
INTEGER::  IW, NOM
LOGICAL:: INTERNAL_WALL_CELL

GRADIENT_WEIGHT = 1.0_EB
IW = CELL(CELL_INDEX(I,J,K))%WALL_INDEX(IOR0) ;  IF (IW == 0) RETURN

INTERNAL_WALL_CELL = (IW > N_EXTERNAL_WALL_CELLS)
NOM = 0 ;  IF (.NOT.INTERNAL_WALL_CELL) NOM = MESHES(NM)%EXTERNAL_WALL(IW)%NOM
IF (WALL(IW)%BOUNDARY_TYPE == SOLID_BOUNDARY .AND. (INTERNAL_WALL_CELL .OR. NOM/=0)) GRADIENT_WEIGHT = 0.0_EB

END FUNCTION GRADIENT_WEIGHT

! --------------------------- FINISH_GLMAT_SOLVER --------------------------------

SUBROUTINE FINISH_GLMAT_SOLVER

USE MPI_F08

! Local variables:
INTEGER :: MAXFCT, MNUM, PHASE, NRHS, ERROR, MSGLVL
#ifdef WITH_MKL
INTEGER :: PERM(1)
#endif

IF (SOLID_PHASE_ONLY .OR. FREEZE_VELOCITY) RETURN

! Solve:
NRHS   =  1
MAXFCT =  1
MNUM   =  1
ERROR  =  0 ! initialize error flag
MSGLVL =  0 ! print statistical information

! Finalize Pardiso:
PHASE = -1

DO IPZ=0,N_ZONE
   ZSL => ZONE_SOLVE(IPZ); IF(ZSL%NUNKH_TOTAL==0) CYCLE

#ifdef WITH_MKL
   CALL CLUSTER_SPARSE_SOLVER(ZSL%PT_H, MAXFCT, MNUM, ZSL%MTYPE, PHASE, ZSL%NUNKH_TOTAL, &
   ZSL%A_H, ZSL%IA_H, ZSL%JA_H, PERM, NRHS, IPARM, MSGLVL, ZSL%F_H, ZSL%X_H, MPI_COMM_WORLD, ERROR)
#endif /* WITH_MKL */
ENDDO

DEALLOCATE(ZONE_SOLVE)

RETURN
END SUBROUTINE FINISH_GLMAT_SOLVER


END MODULE GLOBMAT_SOLVER
