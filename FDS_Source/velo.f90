MODULE VELO

! Module computes the velocity flux terms, baroclinic torque correction terms, and performs the CFL Check

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: SECOND

IMPLICIT NONE

PRIVATE
CHARACTER(255), PARAMETER :: veloid='$Id$'
CHARACTER(255), PARAMETER :: velorev='$Revision$'
CHARACTER(255), PARAMETER :: velodate='$Date$'

PUBLIC COMPUTE_VELOCITY_FLUX,VELOCITY_PREDICTOR,VELOCITY_CORRECTOR,NO_FLUX,GET_REV_velo,MATCH_VELOCITY,VELOCITY_BC
PRIVATE VELOCITY_FLUX,VELOCITY_FLUX_ISOTHERMAL,VELOCITY_FLUX_CYLINDRICAL
 
CONTAINS
 
SUBROUTINE COMPUTE_VELOCITY_FLUX(T,NM,FUNCTION_CODE)

REAL(EB), INTENT(IN) :: T
REAL(EB) :: TNOW
INTEGER, INTENT(IN) :: NM,FUNCTION_CODE

IF (SOLID_PHASE_ONLY) RETURN

TNOW = SECOND()

SELECT CASE(FUNCTION_CODE)
   CASE(1)
      IF (DNS .AND. ISOTHERMAL .AND. N_SPECIES==0) THEN
         CALL VELOCITY_FLUX_ISOTHERMAL(NM)
      ELSE
         IF (PREDICTOR .OR. COMPUTE_VISCOSITY_TWICE) CALL COMPUTE_VISCOSITY(NM)
      ENDIF
   CASE(2)
      IF (PREDICTOR .OR. COMPUTE_VISCOSITY_TWICE) CALL VISCOSITY_BC(NM)
      IF (.NOT.CYLINDRICAL) CALL VELOCITY_FLUX(T,NM)
      IF (     CYLINDRICAL) CALL VELOCITY_FLUX_CYLINDRICAL(T,NM)
!!    CALL VELOCITY_FLUX_CONSERVATIVE ! experimental, currently only works for NMESHES=1
END SELECT

TUSED(4,NM) = TUSED(4,NM) + SECOND() - TNOW
END SUBROUTINE COMPUTE_VELOCITY_FLUX



SUBROUTINE COMPUTE_VISCOSITY(NM)

! Compute turblent eddy viscosity from constant coefficient Smagorinsky model

USE PHYSICAL_FUNCTIONS, ONLY: GET_MU
INTEGER, INTENT(IN) :: NM
REAL(EB) :: DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ,SS,S12,S13,S23,DELTA,CS,MU_SUM
INTEGER :: I,J,K,ITMP,N
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP
 
CALL POINT_TO_MESH(NM)
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   RHOP => RHO
   IF (N_SPECIES > 0) YYP => YY
ELSE
   UU => US
   VV => VS
   WW => WS
   RHOP => RHOS
   IF (N_SPECIES > 0) YYP => YYS
ENDIF

! Compute eddy viscosity using Smagorinsky model

IF (LES) THEN
   CS = CSMAG
   IF (EVACUATION_ONLY(NM)) CS = 0.9_EB
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (TWO_D) THEN
               DELTA = SQRT(DX(I)*DZ(K))
            ELSE
               DELTA = (DX(I)*DY(J)*DZ(K))**ONTH
            ENDIF
            DUDX = RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
            DVDY = RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
            DWDZ = RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
            DUDY = 0.25_EB*RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K)+UU(I-1,J+1,K)-UU(I-1,J-1,K))
            DUDZ = 0.25_EB*RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1)) 
            DVDX = 0.25_EB*RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K)+VV(I+1,J-1,K)-VV(I-1,J-1,K))
            DVDZ = 0.25_EB*RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1)+VV(I,J-1,K+1)-VV(I,J-1,K-1))
            DWDX = 0.25_EB*RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
            DWDY = 0.25_EB*RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K)+WW(I,J+1,K-1)-WW(I,J-1,K-1))
            S12 = 0.5_EB*(DUDY+DVDX)
            S13 = 0.5_EB*(DUDZ+DWDX)
            S23 = 0.5_EB*(DVDZ+DWDY)
            SS = SQRT(2._EB*(DUDX**2 + DVDY**2 + DWDZ**2 + 2._EB*(S12**2 + S13**2 + S23**2)))
            ITMP = 0.1_EB*TMP(I,J,K)
            MU(I,J,K) = SPECIES(0)%MU(ITMP) + RHOP(I,J,K)*(CS*DELTA)**2*SS
         ENDDO
      ENDDO
   ENDDO
ENDIF
   
! Compute viscosity for DNS using primitive species rather than mixture fraction model

IF (DNS .AND. .NOT.MIXTURE_FRACTION) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            ITMP = 0.1_EB*TMP(I,J,K)
            MU_SUM = SPECIES(0)%MU(ITMP)
            DO N=1,N_SPECIES
               MU_SUM = MU_SUM + YYP(I,J,K,N)*(SPECIES(N)%MU(ITMP)-SPECIES(0)%MU(ITMP))
            ENDDO
            MU(I,J,K) = MU_SUM
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Compute viscosity for DNS using mixture fraction model

IF (DNS .AND. MIXTURE_FRACTION) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            ITMP = 0.1_EB*TMP(I,J,K)
            CALL GET_MU(YY(I,J,K,I_Z_MIN:I_Z_MAX),Y_SUM(I,J,K),MU(I,J,K),ITMP)                  
            MU_SUM = MU(I,J,K)*(1-Y_SUM(I,J,K))
            DO N=1,N_SPECIES
               IF(SPECIES(N)%MODE/=MIXTURE_FRACTION_SPECIES) MU_SUM = MU_SUM + YYP(I,J,K,N)*(SPECIES(N)%MU(ITMP)-MU(I,J,K))
            ENDDO
            MU(I,J,K) = MU_SUM
         ENDDO
      ENDDO
   ENDDO
ENDIF 

END SUBROUTINE COMPUTE_VISCOSITY



SUBROUTINE VISCOSITY_BC(NM)

! Specify ghost cell values of the viscosity array MU

INTEGER, INTENT(IN) :: NM
REAL(EB) :: C_RATIO,MU_OTHER
INTEGER :: IIG,JJG,KKG,II,JJ,KK,IW,IIO,JJO,KKO,NOM,N_INT_CELLS

CALL POINT_TO_MESH(NM)

! Mirror viscosity into solids and exterior boundary cells
 
IF (CSMAG==0._EB .OR. CSMAG_WALL>CSMAG) THEN
   C_RATIO = 1._EB
ELSE
   C_RATIO = (CSMAG_WALL/CSMAG)**2
ENDIF

WALL_LOOP: DO IW=1,NWC
   IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_LOOP
   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   SELECT CASE(BOUNDARY_TYPE(IW))
      CASE(SOLID_BOUNDARY)
         MU(IIG,JJG,KKG) = C_RATIO*MU(IIG,JJG,KKG)
         MU(II,JJ,KK) = MU(IIG,JJG,KKG)
      CASE(OPEN_BOUNDARY,MIRROR_BOUNDARY)
         MU(II,JJ,KK) = MU(IIG,JJG,KKG)
      CASE(INTERPOLATED_BOUNDARY)
         NOM = IJKW(9,IW)
         MU_OTHER = 0._EB
         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  MU_OTHER = MU_OTHER + OMESH(NOM)%MU(IIO,JJO,KKO)
               ENDDO
            ENDDO
         ENDDO
         N_INT_CELLS = (IJKW(13,IW)-IJKW(10,IW)+1) * (IJKW(14,IW)-IJKW(11,IW)+1) * (IJKW(15,IW)-IJKW(12,IW)+1)
         MU_OTHER = MU_OTHER/REAL(N_INT_CELLS,EB)
         MU(II,JJ,KK) = MU_OTHER
   END SELECT
ENDDO WALL_LOOP
    
MU(   0,0:JBP1,   0) = MU(   1,0:JBP1,1)
MU(IBP1,0:JBP1,   0) = MU(IBAR,0:JBP1,1)
MU(IBP1,0:JBP1,KBP1) = MU(IBAR,0:JBP1,KBAR)
MU(   0,0:JBP1,KBP1) = MU(   1,0:JBP1,KBAR)
MU(0:IBP1,   0,   0) = MU(0:IBP1,   1,1)
MU(0:IBP1,JBP1,0)    = MU(0:IBP1,JBAR,1)
MU(0:IBP1,JBP1,KBP1) = MU(0:IBP1,JBAR,KBAR)
MU(0:IBP1,0,KBP1)    = MU(0:IBP1,   1,KBAR)
MU(0,   0,0:KBP1)    = MU(   1,   1,0:KBP1)
MU(IBP1,0,0:KBP1)    = MU(IBAR,   1,0:KBP1)
MU(IBP1,JBP1,0:KBP1) = MU(IBAR,JBAR,0:KBP1)
MU(0,JBP1,0:KBP1)    = MU(   1,JBAR,0:KBP1)

! Special periodic boundary conditions

IF (PERIODIC_BC .AND. NMESHES==1) THEN
   MU(0,:,:) = MU(IBAR,:,:)
   MU(:,0,:) = MU(:,JBAR,:)
   MU(:,:,0) = MU(:,:,KBAR)

   MU(IBP1,:,:) = MU(1,:,:)
   MU(:,JBP1,:) = MU(:,1,:)
   MU(:,:,KBP1) = MU(:,:,1)
ENDIF

END SUBROUTINE VISCOSITY_BC



SUBROUTINE VELOCITY_FLUX(T,NM)

! Compute convective and diffusive terms of the momentum equations

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
INTEGER, INTENT(IN) :: NM
REAL(EB) :: T,MUX,MUY,MUZ,UP,UM,VP,VM,WP,WM,VTRM, &
            DTXYDY,DTXZDZ,DTYZDZ,DTXYDX,DTXZDX,DTYZDY, &
            DUDX,DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY, &
            VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX,PMDT,MPDT, &
            AH,RRHO,GX,GY,GZ,TXXP,TXXM,TYYP,TYYM,TZZP,TZZM,DTXXDX,DTYYDY,DTZZDZ, &
            EPSUP,EPSUM,EPSVP,EPSVM,EPSWP,EPSWM,DUMMY=0._EB
INTEGER :: II,JJ,KK,I,J,K,IE,IEC
REAL(EB), POINTER, DIMENSION(:,:,:) :: TXY,TXZ,TYZ,OMX,OMY,OMZ,UU,VV,WW,RHOP,DP
 
CALL POINT_TO_MESH(NM)
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   DP => D  
   RHOP => RHO
ELSE
   UU => US
   VV => VS
   WW => WS
   DP => DS
   RHOP => RHOS
ENDIF

TXY => WORK1
TXZ => WORK2
TYZ => WORK3
OMX => WORK4
OMY => WORK5
OMZ => WORK6
 
! Compute vorticity and stress tensor components

DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         DUDY = RDYN(J)*(UU(I,J+1,K)-UU(I,J,K))
         DVDX = RDXN(I)*(VV(I+1,J,K)-VV(I,J,K))
         DUDZ = RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
         DWDX = RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
         DVDZ = RDZN(K)*(VV(I,J,K+1)-VV(I,J,K))
         DWDY = RDYN(J)*(WW(I,J+1,K)-WW(I,J,K))
         OMX(I,J,K) = DWDY - DVDZ
         OMY(I,J,K) = DUDZ - DWDX
         OMZ(I,J,K) = DVDX - DUDY
         MUX = 0.25_EB*(MU(I,J+1,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I,J+1,K+1))
         MUY = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I+1,J,K+1))
         MUZ = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J+1,K)+MU(I+1,J+1,K))
         TXY(I,J,K) = MUZ*(DVDX + DUDY)
         TXZ(I,J,K) = MUY*(DUDZ + DWDX)
         TYZ(I,J,K) = MUX*(DVDZ + DWDY)
      ENDDO
   ENDDO
ENDDO

! Correct vorticity and stress tensor components at solid edges

EDGE_LOOP: DO IE=1,N_EDGES
   IF (EDGE_TYPE(IE,1)==NULL_EDGE .AND. EDGE_TYPE(IE,2)==NULL_EDGE) CYCLE EDGE_LOOP
   II  = IJKE(1,IE)
   JJ  = IJKE(2,IE)
   KK  = IJKE(3,IE)
   IEC = IJKE(4,IE)
   SELECT CASE(IEC)
      CASE(1)
         OMX(II,JJ,KK) = OME_E(IE)
         TYZ(II,JJ,KK) = TAU_E(IE)
      CASE(2)
         OMY(II,JJ,KK) = OME_E(IE)
         TXZ(II,JJ,KK) = TAU_E(IE)
      CASE(3)
         OMZ(II,JJ,KK) = OME_E(IE)
         TXY(II,JJ,KK) = TAU_E(IE)
   END SELECT
ENDDO EDGE_LOOP

! Compute gravity components
 
GX = EVALUATE_RAMP(T,DUMMY,I_RAMP_GX)*GVEC(1)
GY = EVALUATE_RAMP(T,DUMMY,I_RAMP_GY)*GVEC(2)
GZ = EVALUATE_RAMP(T,DUMMY,I_RAMP_GZ)*GVEC(3)
 
! Upwind/Downwind bias factors
 
IF (PREDICTOR) THEN
   PMDT =  0.5_EB*DT
   MPDT = -0.5_EB*DT
ELSE
   PMDT = -0.5_EB*DT
   MPDT =  0.5_EB*DT
ENDIF
 
! Compute x-direction flux term FVX
 
DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         WP    = WW(I,J,K)   + WW(I+1,J,K)
         WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
         VP    = VV(I,J,K)   + VV(I+1,J,K)
         VM    = VV(I,J-1,K) + VV(I+1,J-1,K)
         EPSWP = 1._EB + WP*MPDT*RDZN(K)
         EPSWM = 1._EB + WM*PMDT*RDZN(K-1)
         EPSVP = 1._EB + VP*MPDT*RDYN(J)
         EPSVM = 1._EB + VM*PMDT*RDYN(J-1)
         WOMY  = EPSWP*WP*OMY(I,J,K) + EPSWM*WM*OMY(I,J,K-1)
         VOMZ  = EPSVP*VP*OMZ(I,J,K) + EPSVM*VM*OMZ(I,J-1,K)
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
         AH    = RHO_0(K)*RRHO - 1._EB   
         DVDY  = (VV(I+1,J,K)-VV(I+1,J-1,K))*RDY(J)
         DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*RDZ(K)
         TXXP  = MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2._EB*(DVDY+DWDZ) )
         DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
         DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
         TXXM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DVDY+DWDZ) )
         DTXXDX= RDXN(I)*(TXXP      -TXXM)
         DTXYDY= RDY(J) *(TXY(I,J,K)-TXY(I,J-1,K))
         DTXZDZ= RDZ(K) *(TXZ(I,J,K)-TXZ(I,J,K-1))
         VTRM  = RRHO*(DTXXDX + DTXYDY + DTXZDZ)
         FVX(I,J,K) = 0.25_EB*(WOMY - VOMZ) + GX*AH - VTRM 
      ENDDO 
   ENDDO   
ENDDO   
 
! Compute y-direction flux term FVY

DO K=1,KBAR
   DO J=0,JBAR
      DO I=1,IBAR
         UP    = UU(I,J,K)   + UU(I,J+1,K)
         UM    = UU(I-1,J,K) + UU(I-1,J+1,K)
         WP    = WW(I,J,K)   + WW(I,J+1,K)
         WM    = WW(I,J,K-1) + WW(I,J+1,K-1)
         EPSUP = 1._EB + UP*MPDT*RDXN(I)
         EPSUM = 1._EB + UM*PMDT*RDXN(I-1)
         EPSWP = 1._EB + WP*MPDT*RDZN(K)
         EPSWM = 1._EB + WM*PMDT*RDZN(K-1)
         WOMX  = EPSWP*WP*OMX(I,J,K) + EPSWM*WM*OMX(I,J,K-1)
         UOMZ  = EPSUP*UP*OMZ(I,J,K) + EPSUM*UM*OMZ(I-1,J,K)
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
         AH    = RHO_0(K)*RRHO - 1._EB
         DUDX  = (UU(I,J+1,K)-UU(I-1,J+1,K))*RDX(I)
         DWDZ  = (WW(I,J+1,K)-WW(I,J+1,K-1))*RDZ(K)
         TYYP  = MU(I,J+1,K)*( FOTH*DP(I,J+1,K) - 2._EB*(DUDX+DWDZ) )
         DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
         DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
         TYYM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DWDZ) )
         DTXYDX= RDX(I) *(TXY(I,J,K)-TXY(I-1,J,K))
         DTYYDY= RDYN(J)*(TYYP      -TYYM)
         DTYZDZ= RDZ(K) *(TYZ(I,J,K)-TYZ(I,J,K-1))
         VTRM  = RRHO*(DTXYDX + DTYYDY + DTYZDZ)
         FVY(I,J,K) = 0.25_EB*(UOMZ - WOMX) + GY*AH - VTRM 
      ENDDO
   ENDDO   
ENDDO   
 
! Compute z-direction flux term FVZ
 
DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         UP    = UU(I,J,K)   + UU(I,J,K+1)
         UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
         VP    = VV(I,J,K)   + VV(I,J,K+1)
         VM    = VV(I,J-1,K) + VV(I,J-1,K+1)
         EPSUP = 1._EB + UP*MPDT*RDXN(I)
         EPSUM = 1._EB + UM*PMDT*RDXN(I-1)
         EPSVP = 1._EB + VP*MPDT*RDYN(J)
         EPSVM = 1._EB + VM*PMDT*RDYN(J-1)
         UOMY  = EPSUP*UP*OMY(I,J,K) + EPSUM*UM*OMY(I-1,J,K)
         VOMX  = EPSVP*VP*OMX(I,J,K) + EPSVM*VM*OMX(I,J-1,K)
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
         AH    = 0.5_EB*(RHO_0(K)+RHO_0(K+1))*RRHO - 1._EB
         DUDX  = (UU(I,J,K+1)-UU(I-1,J,K+1))*RDX(I)
         DVDY  = (VV(I,J,K+1)-VV(I,J-1,K+1))*RDY(J)
         TZZP  = MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2._EB*(DUDX+DVDY) )
         DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
         DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
         TZZM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DVDY) )
         DTXZDX= RDX(I) *(TXZ(I,J,K)-TXZ(I-1,J,K))
         DTYZDY= RDY(J) *(TYZ(I,J,K)-TYZ(I,J-1,K))
         DTZZDZ= RDZN(K)*(TZZP      -TZZM)
         VTRM  = RRHO*(DTXZDX + DTYZDY + DTZZDZ)
         FVZ(I,J,K) = 0.25_EB*(VOMX - UOMY) + GZ*AH - VTRM       
      ENDDO
   ENDDO   
ENDDO   
 
! Baroclinic torque correction
 
IF (BAROCLINIC) CALL BAROCLINIC_CORRECTION
 
! Adjust FVX, FVY and FVZ at solid, internal obstructions for no flux
 
CALL NO_FLUX
 
END SUBROUTINE VELOCITY_FLUX



SUBROUTINE VELOCITY_FLUX_CONSERVATIVE !! RJM
 
! Compute RHS of momentum equation in conservative form for periodic bcs

REAL(EB) :: MUX,MUY,MUZ,RHOX,RHOY,RHOZ,DIVV, &
            UUP,VVP,WWP,UUN,UUT,VVE,VVT,WWE,WWN, &
            DUDX,DUDY,DUDZ, &
            DVDX,DVDY,DVDZ, &
            DWDX,DWDY,DWDZ, &
            SS,S12,S13,S23,DELTA 
INTEGER :: I,J,K,IP1,JP1,KP1,IM1,JM1,KM1,NM,ITMP
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP,TXX,TYY,TZZ,TXY,TXZ,TYZ

NM = 1 ! periodic bcs only work for a single mesh
CALL POINT_TO_MESH(NM)
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   RHOP => RHO
ELSE
   UU => US
   VV => VS
   WW => WS
   RHOP => RHOS
ENDIF

TXX => WORK1
TYY => WORK2
TZZ => WORK3
TXY => WORK4
TXZ => WORK5
TYZ => WORK6

TXX = 0._EB
TYY = 0._EB
TZZ = 0._EB
TXY = 0._EB
TXZ = 0._EB
TYZ = 0._EB

! Smagorinsky model
IF (LES) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
      
            IP1 = I+1
            JP1 = J+1
            KP1 = K+1
            IF (IP1>IBAR) IP1=IP1-IBAR
            IF (JP1>JBAR) JP1=JP1-JBAR
            IF (KP1>KBAR) KP1=KP1-KBAR
         
            IM1 = I-1
            JM1 = J-1
            KM1 = K-1
            IF (IM1<1) IM1=IM1+IBAR
            IF (JM1<1) JM1=JM1+JBAR
            IF (KM1<1) KM1=KM1+KBAR
      
            DELTA = (DX(I)*DY(J)*DZ(K))**ONTH
         
            DUDX = RDX(I)*(UU(I,J,K)-UU(IM1,J,K))
            DVDY = RDY(J)*(VV(I,J,K)-VV(I,JM1,K))
            DWDZ = RDZ(K)*(WW(I,J,K)-WW(I,J,KM1))
            DUDY = 0.25_EB*RDY(J)*(UU(I,JP1,K)-UU(I,JM1,K)+UU(IM1,JP1,K)-UU(IM1,JM1,K))
            DVDX = 0.25_EB*RDX(I)*(VV(IP1,J,K)-VV(IM1,J,K)+VV(IP1,JM1,K)-VV(IM1,JM1,K))
            DWDX = 0.25_EB*RDX(I)*(WW(IP1,J,K)-WW(IM1,J,K)+WW(IP1,J,KM1)-WW(IM1,J,KM1))
            DWDY = 0.25_EB*RDY(J)*(WW(I,JP1,K)-WW(I,JM1,K)+WW(I,JP1,KM1)-WW(I,JM1,KM1))
            DUDZ = 0.25_EB*RDZ(K)*(UU(I,J,KP1)-UU(I,J,KM1)+UU(IM1,J,KP1)-UU(IM1,J,KM1)) 
            DVDZ = 0.25_EB*RDZ(K)*(VV(I,J,KP1)-VV(I,J,KM1)+VV(I,JM1,KP1)-VV(I,JM1,KM1))
         
            S12 = 0.5_EB*(DUDY+DVDX)
            S13 = 0.5_EB*(DUDZ+DWDX)
            S23 = 0.5_EB*(DVDZ+DWDY)
            SS = SQRT(2._EB*(DUDX**2 + DVDY**2 + DWDZ**2 + 2._EB*(S12**2 + S13**2 + S23**2)))
         
            ITMP = 0.1_EB*TMP(I,J,K)
            MU(I,J,K) = SPECIES(0)%MU(ITMP) + RHOP(I,J,K)*(CSMAG*DELTA)**2*SS
         ENDDO
      ENDDO
   ENDDO
ENDIF


! Compute stress tensor components
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
      
         IP1 = I+1
         JP1 = J+1
         KP1 = K+1
         IF (IP1>IBAR) IP1=IP1-IBAR
         IF (JP1>JBAR) JP1=JP1-JBAR
         IF (KP1>KBAR) KP1=KP1-KBAR
         
         IM1 = I-1
         JM1 = J-1
         KM1 = K-1
         IF (IM1<1) IM1=IM1+IBAR
         IF (JM1<1) JM1=JM1+JBAR
         IF (KM1<1) KM1=KM1+KBAR
         
         ! average face velocities for advective terms
         UUP = 0.5_EB*(UU(IM1,J,K)+UU(I,J,K))
         VVP = 0.5_EB*(VV(I,JM1,K)+VV(I,J,K))
         WWP = 0.5_EB*(WW(I,J,KM1)+WW(I,J,K))        
         UUN = 0.5_EB*(UU(I,J,K)+UU(I,JP1,K))
         UUT = 0.5_EB*(UU(I,J,K)+UU(I,J,KP1))
         VVE = 0.5_EB*(VV(I,J,K)+VV(IP1,J,K))
         VVT = 0.5_EB*(VV(I,J,K)+VV(I,J,KP1))
         WWE = 0.5_EB*(WW(I,J,K)+WW(IP1,J,K))
         WWN = 0.5_EB*(WW(I,J,K)+WW(I,JP1,K))
         
         ! velocity gradient tensor
         DUDX = RDX(I)*(UU(I,J,K)-UU(IM1,J,K))
         DVDY = RDY(J)*(VV(I,J,K)-VV(I,JM1,K))
         DWDZ = RDZ(K)*(WW(I,J,K)-WW(I,J,KM1))
         DUDY = RDYN(J)*(UU(I,JP1,K)-UU(I,J,K))
         DUDZ = RDZN(K)*(UU(I,J,KP1)-UU(I,J,K))
         DVDX = RDXN(I)*(VV(IP1,J,K)-VV(I,J,K))
         DVDZ = RDZN(K)*(VV(I,J,KP1)-VV(I,J,K))
         DWDX = RDXN(I)*(WW(IP1,J,K)-WW(I,J,K))
         DWDY = RDYN(J)*(WW(I,JP1,K)-WW(I,J,K))
         
         ! average viscosity to cell vertex
         MUX = 0.25_EB*(MU(I,JP1,K)+MU(I,J,K)+MU(I,J,KP1)+MU(I,JP1,KP1))
         MUY = 0.25_EB*(MU(IP1,J,K)+MU(I,J,K)+MU(I,J,KP1)+MU(IP1,J,KP1))
         MUZ = 0.25_EB*(MU(IP1,J,K)+MU(I,J,K)+MU(I,JP1,K)+MU(IP1,JP1,K))
         
         ! average density to cell vertex
         RHOX = 0.25_EB*(RHOP(I,JP1,K)+RHOP(I,J,K)+RHOP(I,J,KP1)+RHOP(I,JP1,KP1))
         RHOY = 0.25_EB*(RHOP(IP1,J,K)+RHOP(I,J,K)+RHOP(I,J,KP1)+RHOP(IP1,J,KP1))
         RHOZ = 0.25_EB*(RHOP(IP1,J,K)+RHOP(I,J,K)+RHOP(I,JP1,K)+RHOP(IP1,JP1,K))
         
         ! divergence
         DIVV = DUDX+DVDY+DWDZ !! if DIVV\=0, there is a problem
         
         ! component-normal advective and viscous flux
         ! stored at cell center
         TXX(I,J,K) = RHOP(I,J,K)*UUP*UUP - 2._EB*MU(I,J,K)*(DUDX-ONTH*DIVV)
         TYY(I,J,K) = RHOP(I,J,K)*VVP*VVP - 2._EB*MU(I,J,K)*(DVDY-ONTH*DIVV)
         TZZ(I,J,K) = RHOP(I,J,K)*WWP*WWP - 2._EB*MU(I,J,K)*(DWDZ-ONTH*DIVV)
         
         ! off-diagonal tensor components
         ! stored at cell vertices
         TXY(I,J,K) = RHOZ*UUN*VVE - MUZ*(DUDY + DVDX)
         TXZ(I,J,K) = RHOY*UUT*WWE - MUY*(DUDZ + DWDX)
         TYZ(I,J,K) = RHOX*VVT*WWN - MUX*(DVDZ + DWDY)  
      ENDDO
   ENDDO
ENDDO

! Compute force term
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         
         IP1 = I+1
         JP1 = J+1
         KP1 = K+1
         IF (IP1>IBAR) IP1=IP1-IBAR
         IF (JP1>JBAR) JP1=JP1-JBAR
         IF (KP1>KBAR) KP1=KP1-KBAR
         
         IM1 = I-1
         JM1 = J-1
         KM1 = K-1
         IF (IM1<1) IM1=IM1+IBAR
         IF (JM1<1) JM1=JM1+JBAR
         IF (KM1<1) KM1=KM1+KBAR
      
         FVX(I,J,K) = RDXN(I)*(TXX(IP1,J,K)-TXX(I,J,K)) &
                    + RDY(J) *(TXY(I,J,K)-TXY(I,JM1,K)) &
                    + RDZ(K) *(TXZ(I,J,K)-TXZ(I,J,KM1))
                    
         FVY(I,J,K) = RDX(I) *(TXY(I,J,K)-TXY(IM1,J,K)) &
                    + RDYN(J)*(TYY(I,JP1,K)-TYY(I,J,K)) &
                    + RDZ(K) *(TYZ(I,J,K)-TYZ(I,J,KM1))
                    
         FVZ(I,J,K) = RDX(I) *(TXZ(I,J,K)-TXZ(IM1,J,K)) &
                    + RDY(J) *(TYZ(I,J,K)-TYZ(I,JM1,K)) &
                    + RDZN(K)*(TZZ(I,J,KP1)-TZZ(I,J,K))
      ENDDO
   ENDDO
ENDDO

! Baroclinic torque correction
IF (BAROCLINIC) CALL BAROCLINIC_CORRECTION

! apply periodicity
FVX(0,:,:) = FVX(IBAR,:,:)
FVY(:,0,:) = FVY(:,JBAR,:)
FVZ(:,:,0) = FVZ(:,:,KBAR)
 
END SUBROUTINE VELOCITY_FLUX_CONSERVATIVE
 
 
SUBROUTINE VELOCITY_FLUX_ISOTHERMAL(NM)
 
! Compute the velocity flux at cell edges (ISOTHERMAL DNS ONLY)
 
REAL(EB) :: UP,UM,VP,VM,WP,WM,VTRM,VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX, &
            DVDZ,DVDX,DWDY,DWDX,DUDZ,DUDY,PMDT,MPDT, &
            EPSUP,EPSUM,EPSVP,EPSVM,EPSWP,EPSWM
INTEGER :: I,J,K,II,JJ,KK,IE,IEC
REAL(EB), POINTER, DIMENSION(:,:,:) :: OMX,OMY,OMZ,UU,VV,WW
INTEGER, INTENT(IN) :: NM
 
CALL POINT_TO_MESH(NM)
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
ELSE
   UU => US
   VV => VS
   WW => WS
ENDIF
 
OMX => WORK4
OMY => WORK5
OMZ => WORK6
 
! Compute vorticity and stress tensor components
 
DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         DUDY = RDYN(J)*(UU(I,J+1,K)-UU(I,J,K))
         DVDX = RDXN(I)*(VV(I+1,J,K)-VV(I,J,K))
         DUDZ = RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
         DWDX = RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
         DVDZ = RDZN(K)*(VV(I,J,K+1)-VV(I,J,K))
         DWDY = RDYN(J)*(WW(I,J+1,K)-WW(I,J,K))
         OMX(I,J,K) = DWDY - DVDZ
         OMY(I,J,K) = DUDZ - DWDX
         OMZ(I,J,K) = DVDX - DUDY
      ENDDO
   ENDDO
ENDDO
 
! Correct vorticity and stress tensor components at solid edges
 
EDGE_LOOP: DO IE=1,N_EDGES
   IF (EDGE_TYPE(IE,1)==NULL_EDGE .AND. EDGE_TYPE(IE,2)==NULL_EDGE) CYCLE EDGE_LOOP
   II  = IJKE(1,IE)
   JJ  = IJKE(2,IE)
   KK  = IJKE(3,IE)
   IEC = IJKE(4,IE)
   SELECT CASE(IEC)
      CASE(1)
         OMX(II,JJ,KK) = OME_E(IE)
      CASE(2)
         OMY(II,JJ,KK) = OME_E(IE)
      CASE(3)
         OMZ(II,JJ,KK) = OME_E(IE)
   END SELECT
ENDDO EDGE_LOOP
 
! Upwind/Downwind bias factors
 
IF (PREDICTOR) THEN
   PMDT =  0.5_EB*DT
   MPDT = -0.5_EB*DT
ELSE
   PMDT = -0.5_EB*DT
   MPDT =  0.5_EB*DT
ENDIF
 
! Compute x-direction flux term FVX
 
DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         WP    = WW(I,J,K)   + WW(I+1,J,K)
         WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
         VP    = VV(I,J,K)   + VV(I+1,J,K)
         VM    = VV(I,J-1,K) + VV(I+1,J-1,K)
         EPSWP = 1._EB + WP*MPDT*RDZN(K)
         EPSWM = 1._EB + WM*PMDT*RDZN(K-1)
         EPSVP = 1._EB + VP*MPDT*RDYN(J)
         EPSVM = 1._EB + VM*PMDT*RDYN(J-1)
         WOMY  = EPSWP*WP*OMY(I,J,K) + EPSWM*WM*OMY(I,J,K-1)
         VOMZ  = EPSVP*VP*OMZ(I,J,K) + EPSVM*VM*OMZ(I,J-1,K)
         VTRM  = RREDZ(K)*(OMY(I,J,K)-OMY(I,J,K-1)) - RREDY(J)*(OMZ(I,J,K)-OMZ(I,J-1,K))
         FVX(I,J,K) = 0.25_EB*(WOMY - VOMZ) - VTRM
      ENDDO
   ENDDO
ENDDO
! Compute y-direction flux term FVY
 
DO K=1,KBAR
   DO J=0,JBAR
      DO I=1,IBAR
         UP    = UU(I,J,K)   + UU(I,J+1,K)
         UM    = UU(I-1,J,K) + UU(I-1,J+1,K)
         WP    = WW(I,J,K)   + WW(I,J+1,K)
         WM    = WW(I,J,K-1) + WW(I,J+1,K-1)
         EPSUP = 1._EB + UP*MPDT*RDXN(I)
         EPSUM = 1._EB + UM*PMDT*RDXN(I-1)
         EPSWP = 1._EB + WP*MPDT*RDZN(K)
         EPSWM = 1._EB + WM*PMDT*RDZN(K-1)
         WOMX  = EPSWP*WP*OMX(I,J,K) + EPSWM*WM*OMX(I,J,K-1)
         UOMZ  = EPSUP*UP*OMZ(I,J,K) + EPSUM*UM*OMZ(I-1,J,K)
         VTRM  = RREDX(I)*(OMZ(I,J,K)-OMZ(I-1,J,K)) - RREDZ(K)*(OMX(I,J,K)-OMX(I,J,K-1))
         FVY(I,J,K) = 0.25_EB*(UOMZ - WOMX) - VTRM
      ENDDO
   ENDDO
ENDDO
 
! Compute z-direction flux term FVZ
 
DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         UP    = UU(I,J,K)   + UU(I,J,K+1)
         UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
         VP    = VV(I,J,K)   + VV(I,J,K+1)
         VM    = VV(I,J-1,K) + VV(I,J-1,K+1)
         EPSUP = 1._EB + UP*MPDT*RDXN(I)
         EPSUM = 1._EB + UM*PMDT*RDXN(I-1)
         EPSVP = 1._EB + VP*MPDT*RDYN(J)
         EPSVM = 1._EB + VM*PMDT*RDYN(J-1)
         UOMY  = EPSUP*UP*OMY(I,J,K) + EPSUM*UM*OMY(I-1,J,K)
         VOMX  = EPSVP*VP*OMX(I,J,K) + EPSVM*VM*OMX(I,J-1,K)
         VTRM  = RREDY(J)*(OMX(I,J,K)-OMX(I,J-1,K)) - RREDX(I)*(OMY(I,J,K)-OMY(I-1,J,K))
         FVZ(I,J,K) = 0.25_EB*(VOMX - UOMY) - VTRM
      ENDDO
   ENDDO
ENDDO
 
! Adjust FVX, FVY and FVZ at solid, internal obstructions for no flux
 
CALL NO_FLUX
 
END SUBROUTINE VELOCITY_FLUX_ISOTHERMAL
 
 
SUBROUTINE VELOCITY_FLUX_CYLINDRICAL(T,NM)

! Compute convective and diffusive terms for 2D axisymmetric

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP 
REAL(EB) :: T,DMUDX
INTEGER :: I0
INTEGER, INTENT(IN) :: NM
REAL(EB) :: MUY,UP,UM,WP,WM,VTRM,DTXZDZ,DTXZDX,DUDX,DWDZ,DUDZ,DWDX,WOMY,UOMY,PMDT,MPDT, &
            AH,RRHO,GX,GZ,TXXP,TXXM,TZZP,TZZM,DTXXDX,DTZZDZ,EPSUP,EPSUM,EPSWP,EPSWM,DUMMY=0._EB
INTEGER :: II,JJ,KK,I,J,K,IE,IEC
REAL(EB), POINTER, DIMENSION(:,:,:) :: TXZ,OMY,UU,WW,RHOP,DP
 
CALL POINT_TO_MESH(NM)
 
IF (PREDICTOR) THEN
   UU => U
   WW => W
   DP => D  
   RHOP => RHO
ELSE
   UU => US
   WW => WS
   DP => DS
   RHOP => RHOS
ENDIF
 
TXZ => WORK2
OMY => WORK5
 
! Compute vorticity and stress tensor components
 
DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         DUDZ = RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
         DWDX = RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
         OMY(I,J,K) = DUDZ - DWDX
         MUY = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I+1,J,K+1))
         TXZ(I,J,K) = MUY*(DUDZ + DWDX)
      ENDDO
   ENDDO
ENDDO
 
! Correct vorticity and stress tensor components at solid edges
 
EDGE_LOOP: DO IE=1,N_EDGES
   IF (EDGE_TYPE(IE,1)==NULL_EDGE .AND. EDGE_TYPE(IE,2)==NULL_EDGE) CYCLE EDGE_LOOP
   II  = IJKE(1,IE)
   JJ  = IJKE(2,IE)
   KK  = IJKE(3,IE)
   IEC = IJKE(4,IE)
   SELECT CASE(IEC)
      CASE(2)
         OMY(II,JJ,KK) = OME_E(IE)
         TXZ(II,JJ,KK) = TAU_E(IE)
   END SELECT
ENDDO EDGE_LOOP
 
! Compute gravity components
 
GX  = 0._EB
GZ  = EVALUATE_RAMP(T,DUMMY,I_RAMP_GZ)*GVEC(3)
 
! Upwind/Downwind bias factors
 
IF (PREDICTOR) THEN
   PMDT =  0.5_EB*DT
   MPDT = -0.5_EB*DT
ELSE
   PMDT = -0.5_EB*DT
   MPDT =  0.5_EB*DT
ENDIF
 
! Compute r-direction flux term FVX
 
IF (XS==0) THEN 
   I0 = 1
ELSE
   I0 = 0
ENDIF
 
J = 1
DO K= 1,KBAR
   DO I=I0,IBAR
      WP    = WW(I,J,K)   + WW(I+1,J,K)
      WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
      EPSWP = 1._EB + WP*MPDT*RDZN(K)
      EPSWM = 1._EB + WM*PMDT*RDZN(K-1)
      WOMY  = EPSWP*WP*OMY(I,J,K) + EPSWM*WM*OMY(I,J,K-1)
      RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
      AH    = RHO_0(K)*RRHO - 1._EB   
      DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*RDZ(K)
      TXXP  = MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2._EB*DWDZ )
      DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
      TXXM  = MU(I,J,K)  *( FOTH*DP(I,J,K) -2._EB*DWDZ )
      DTXXDX= RDXN(I)*(TXXP      -TXXM)
      DTXZDZ= RDZ(K) *(TXZ(I,J,K)-TXZ(I,J,K-1))
      DMUDX = (MU(I+1,J,K)-MU(I,J,K))*RDXN(I)
      VTRM  = RRHO*( DTXXDX + DTXZDZ - 2._EB*UU(I,J,K)*DMUDX/R(I) ) 
      FVX(I,J,K) = 0.25_EB*WOMY + GX*AH - VTRM 
   ENDDO
ENDDO   
 
! Compute z-direction flux term FVZ
 
J = 1
DO K=0,KBAR
   DO I=1,IBAR
      UP    = UU(I,J,K)   + UU(I,J,K+1)
      UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
      EPSUP = 1._EB + UP*MPDT*RDXN(I)
      EPSUM = 1._EB + UM*PMDT*RDXN(I-1)
      UOMY  = EPSUP*UP*OMY(I,J,K) + EPSUM*UM*OMY(I-1,J,K)
      RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
      AH    = 0.5_EB*(RHO_0(K)+RHO_0(K+1))*RRHO - 1._EB
      DUDX  = (R(I)*UU(I,J,K+1)-R(I-1)*UU(I-1,J,K+1))*RDX(I)*RRN(I)
      TZZP  = MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2._EB*DUDX )
      DUDX  = (R(I)*UU(I,J,K)-R(I-1)*UU(I-1,J,K))*RDX(I)*RRN(I)
      TZZM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*DUDX )
      DTXZDX= RDX(I) *(R(I)*TXZ(I,J,K)-R(I-1)*TXZ(I-1,J,K))*RRN(I)
      DTZZDZ= RDZN(K)*(TZZP      -TZZM)
      VTRM  = RRHO*(DTXZDX + DTZZDZ)
      FVZ(I,J,K) = -0.25_EB*UOMY + GZ*AH - VTRM 
   ENDDO
ENDDO   
 
! Baroclinic torque correction terms
 
IF (BAROCLINIC) CALL BAROCLINIC_CORRECTION
 
! Adjust FVX and FVZ at solid, internal obstructions for no flux
 
CALL NO_FLUX
 
END SUBROUTINE VELOCITY_FLUX_CYLINDRICAL
 
 
SUBROUTINE NO_FLUX

! Set FVX,FVY,FVZ inside internal blockages to maintain no flux

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP 
REAL(EB) :: RFODT
INTEGER  :: IC2,IC1,N,I,J,K,IW,II,JJ,KK,IOR
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
TYPE (OBSTRUCTION_TYPE), POINTER :: OB
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
ELSE
   UU => US
   VV => VS
   WW => WS
ENDIF
 
RFODT = RELAXATION_FACTOR/DT
 
! Drive velocity components inside solid obstructions towards zero
 
OBST_LOOP: DO N=1,N_OBST
   OB=>OBSTRUCTION(N)
   DO K=OB%K1+1,OB%K2
      DO J=OB%J1+1,OB%J2
         LOOP1: DO I=OB%I1  ,OB%I2
            IC1 = CELL_INDEX(I,J,K)
            IC2 = CELL_INDEX(I+1,J,K)
            IF (SOLID(IC1) .AND. SOLID(IC2)) FVX(I,J,K) = -RDXN(I)*(H(I+1,J,K)-H(I,J,K)) + RFODT*UU(I,J,K)
         ENDDO LOOP1
      ENDDO 
   ENDDO 
   DO K=OB%K1+1,OB%K2
      DO J=OB%J1  ,OB%J2
         LOOP2: DO I=OB%I1+1,OB%I2
            IC1 = CELL_INDEX(I,J,K)
            IC2 = CELL_INDEX(I,J+1,K)
            IF (SOLID(IC1) .AND. SOLID(IC2)) FVY(I,J,K) = -RDYN(J)*(H(I,J+1,K)-H(I,J,K)) + RFODT*VV(I,J,K)
         ENDDO LOOP2
      ENDDO 
   ENDDO 
   DO K=OB%K1  ,OB%K2
      DO J=OB%J1+1,OB%J2
         LOOP3: DO I=OB%I1+1,OB%I2
            IC1 = CELL_INDEX(I,J,K)
            IC2 = CELL_INDEX(I,J,K+1)
            IF (SOLID(IC1) .AND. SOLID(IC2)) FVZ(I,J,K) = -RDZN(K)*(H(I,J,K+1)-H(I,J,K)) + RFODT*WW(I,J,K)
         ENDDO LOOP3
      ENDDO 
   ENDDO 
ENDDO OBST_LOOP
 
! Add normal velocity to FVX, etc. for surface cells
 
WALL_LOOP: DO IW=1,NWC
 
IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) CYCLE WALL_LOOP
IF (BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP
II  = IJKW(1,IW)
JJ  = IJKW(2,IW)
KK  = IJKW(3,IW)
IOR = IJKW(4,IW)
SELECT CASE (BOUNDARY_TYPE(IW))
   CASE (SOLID_BOUNDARY,POROUS_BOUNDARY)
      SELECT CASE(IOR)
         CASE( 1) 
            FVX(II,JJ,KK)   = -RDXN(II)  *(H(II+1,JJ,KK)-H(II,JJ,KK)) + RFODT*(UU(II,JJ,KK)   + UWS(IW))
         CASE(-1) 
            FVX(II-1,JJ,KK) = -RDXN(II-1)*(H(II,JJ,KK)-H(II-1,JJ,KK)) + RFODT*(UU(II-1,JJ,KK) - UWS(IW))
         CASE( 2) 
            FVY(II,JJ,KK)   = -RDYN(JJ)  *(H(II,JJ+1,KK)-H(II,JJ,KK)) + RFODT*(VV(II,JJ,KK)   + UWS(IW))
         CASE(-2)
            FVY(II,JJ-1,KK) = -RDYN(JJ-1)*(H(II,JJ,KK)-H(II,JJ-1,KK)) + RFODT*(VV(II,JJ-1,KK) - UWS(IW))
         CASE( 3) 
            FVZ(II,JJ,KK)   = -RDZN(KK)  *(H(II,JJ,KK+1)-H(II,JJ,KK)) + RFODT*(WW(II,JJ,KK)   + UWS(IW))
         CASE(-3) 
            FVZ(II,JJ,KK-1) = -RDZN(KK-1)*(H(II,JJ,KK)-H(II,JJ,KK-1)) + RFODT*(WW(II,JJ,KK-1) - UWS(IW))
      END SELECT
   CASE (MIRROR_BOUNDARY)
      SELECT CASE(IOR)
         CASE( 1)
            FVX(II  ,JJ,KK) = 0._EB
         CASE(-1)
            FVX(II-1,JJ,KK) = 0._EB
         CASE( 2)
            FVY(II  ,JJ,KK) = 0._EB
         CASE(-2)
            FVY(II,JJ-1,KK) = 0._EB
         CASE( 3)
            FVZ(II  ,JJ,KK) = 0._EB
         CASE(-3)
            FVZ(II,JJ,KK-1) = 0._EB
      END SELECT
END SELECT 
 
ENDDO WALL_LOOP
 
END SUBROUTINE NO_FLUX
 
 

SUBROUTINE VELOCITY_PREDICTOR(NM,STOP_STATUS)

! Estimates the velocity components at the next time step

REAL(EB) :: TNOW
INTEGER  :: STOP_STATUS,I,J,K
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: HT

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY) THEN
   CALL CHECK_STABILITY(NM)
   RETURN
ENDIF

TNOW=SECOND() 
CALL POINT_TO_MESH(NM)

IF (PRESSURE_CORRECTION) THEN
   HT => WORK1
   HT = H + HP
ELSE
   HT => H
ENDIF

DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         US(I,J,K) = U(I,J,K) - DT*( FVX(I,J,K) + RDXN(I)*(HT(I+1,J,K)-HT(I,J,K)) )
      ENDDO
   ENDDO
ENDDO

DO K=1,KBAR
   DO J=0,JBAR
      DO I=1,IBAR
         VS(I,J,K) = V(I,J,K) - DT*( FVY(I,J,K) + RDYN(J)*(HT(I,J+1,K)-HT(I,J,K)) )
      ENDDO
   ENDDO
ENDDO

DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         WS(I,J,K) = W(I,J,K) - DT*( FVZ(I,J,K) + RDZN(K)*(HT(I,J,K+1)-HT(I,J,K)) )
      ENDDO
   ENDDO
ENDDO

! Check the stability criteria, and if the time step is too small, send back a signal to kill the job
 
DTOLD = DT
CALL CHECK_STABILITY(NM)
 
IF (DT<DTINT*1.E-4) STOP_STATUS = INSTABILITY_STOP
 
TUSED(4,NM)=TUSED(4,NM)+SECOND()-TNOW
END SUBROUTINE VELOCITY_PREDICTOR
 
 

SUBROUTINE VELOCITY_CORRECTOR(NM)

! Correct the velocity components

REAL(EB) :: TNOW
INTEGER  :: I,J,K
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: HT
 
IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY) RETURN

TNOW=SECOND() 
CALL POINT_TO_MESH(NM)

IF (PRESSURE_CORRECTION) THEN
   HT => WORK1
   HT = H + HP
ELSE
   HT => H
ENDIF

DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         U(I,J,K) = .5_EB*( U(I,J,K) + US(I,J,K) - DT*(FVX(I,J,K) + RDXN(I)*(HT(I+1,J,K)-HT(I,J,K))) )
      ENDDO
   ENDDO
ENDDO

DO K=1,KBAR
   DO J=0,JBAR
      DO I=1,IBAR
         V(I,J,K) = .5_EB*( V(I,J,K) + VS(I,J,K) - DT*(FVY(I,J,K) + RDYN(J)*(HT(I,J+1,K)-HT(I,J,K))) )
      ENDDO
   ENDDO
ENDDO

DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         W(I,J,K) = .5_EB*( W(I,J,K) + WS(I,J,K) - DT*(FVZ(I,J,K) + RDZN(K)*(HT(I,J,K+1)-HT(I,J,K))) )
      ENDDO
   ENDDO
ENDDO

TUSED(4,NM)=TUSED(4,NM)+SECOND()-TNOW
END SUBROUTINE VELOCITY_CORRECTOR
 

 
SUBROUTINE VELOCITY_BC(T,NM)

! Assert tangential velocity boundary conditions

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP 
REAL(EB), INTENT(IN) :: T
REAL(EB) :: BC(2),MUA,FVT(2,2),UP,UM,VP,VM,WP,WM,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY,TSI,WGT,TNOW,PROF(2),RAMP_T
INTEGER  :: I,J,K,IBC(2),NOM(2),IIO(2),JJO(2),KKO(2),IE,II,JJ,KK,IEC,IOR(2),IWM,IWP,ICMM,ICMP,ICPM,ICPP,IC,N_IOR
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,U_Y,U_Z,V_X,V_Z,W_X,W_Y,OM_UU,OM_VV,OM_WW
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE (OMESH_TYPE), POINTER :: OM

IF (SOLID_PHASE_ONLY) RETURN

TNOW = SECOND()

! Assign local names to variables

CALL POINT_TO_MESH(NM)

! Point to the appropriate velocity field

IF (PREDICTOR) THEN
   UU => US
   VV => VS
   WW => WS
ELSE
   UU => U
   VV => V
   WW => W
ENDIF

! Set the boundary velocity place holder to some large negative number

IF (CORRECTOR) THEN
   UVW_GHOST = -1.E6_EB
   U_Y => WORK1
   U_Z => WORK2
   V_X => WORK3
   V_Z => WORK4
   W_X => WORK5
   W_Y => WORK6
   U_Y = -1.E6_EB
   U_Z = -1.E6_EB
   V_X = -1.E6_EB
   V_Z = -1.E6_EB
   W_X = -1.E6_EB
   W_Y = -1.E6_EB
ENDIF

! Loop over all cell edges and determine the appropriate velocity BCs

EDGE_LOOP: DO IE=1,N_EDGES

   IF (EDGE_TYPE(IE,1)==NULL_EDGE .AND. EDGE_TYPE(IE,2)==NULL_EDGE) CYCLE EDGE_LOOP

   II     = IJKE( 1,IE)
   JJ     = IJKE( 2,IE)
   KK     = IJKE( 3,IE)
   IEC    = IJKE( 4,IE)
   IBC(1) = IJKE( 5,IE)
   IOR(1) = IJKE( 6,IE)
   IBC(2) = IJKE( 7,IE)
   IOR(2) = IJKE( 8,IE)
   NOM(1) = IJKE( 9,IE)
   IIO(1) = IJKE(10,IE)
   JJO(1) = IJKE(11,IE)
   KKO(1) = IJKE(12,IE)
   NOM(2) = IJKE(13,IE)
   IIO(2) = IJKE(14,IE)
   JJO(2) = IJKE(15,IE)
   KKO(2) = IJKE(16,IE)

   DO N_IOR=1,2
      SF  => SURFACE(IBC(N_IOR))
      BC(N_IOR)  = SF%SLIP_FACTOR
      IF (BC(N_IOR)>1.5_EB) THEN
         IF (SF%T_IGN==T_BEGIN .AND. SF%RAMP_INDEX(TIME_VELO)>=1) THEN
            TSI = T
         ELSE
            TSI=T-SF%T_IGN
         ENDIF
         PROF(N_IOR) = 1._EB
         IF (SF%PROFILE==ATMOSPHERIC) PROF(N_IOR)=((ZC(KK)-GROUND_LEVEL)/SF%Z0)**SF%PLE
         RAMP_T = EVALUATE_RAMP(TSI,SF%TAU(TIME_VELO),SF%RAMP_INDEX(TIME_VELO))
         FVT(N_IOR,1) = RAMP_T*SF%VEL_T(1)
         FVT(N_IOR,2) = RAMP_T*SF%VEL_T(2)
      ENDIF
   ENDDO
 
   COMPONENT: SELECT CASE(IEC)
 
      CASE(1) COMPONENT  ! Treat edges that point in the x direction (omega_x, tau_x) 

         ICMM = CELL_INDEX(II,MAX(1,JJ)     ,MAX(1,KK)     )
         ICMP = CELL_INDEX(II,MAX(1,JJ)     ,MIN(KBAR,KK+1))
         ICPM = CELL_INDEX(II,MIN(JBAR,JJ+1),MAX(1,KK)     )
         ICPP = CELL_INDEX(II,MIN(JBAR,JJ+1),MIN(KBAR,KK+1))
         VP   = VV(II,JJ,KK+1)
         VM   = VV(II,JJ,KK)
         WP   = WW(II,JJ+1,KK)
         WM   = WW(II,JJ,KK)

         IF (NOM(2)==0 .OR. EDGE_TYPE(IE,2)/=INTERPOLATED_EDGE) THEN
            SELECT CASE(IOR(2))
               CASE(-3)
                  IWM = WALL_INDEX(ICMM, 3) 
                  IWP = WALL_INDEX(ICPM, 3) 
                  IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .OR. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
                     IF (BC(2)<1.5_EB) VP = BC(2)*VM
                     IF (BC(2)>1.5_EB) VP = FVT(2,2)*PROF(2)
                  ENDIF
               CASE( 3)
                  IWM = WALL_INDEX(ICMP,-3) 
                  IWP = WALL_INDEX(ICPP,-3) 
                  IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .OR. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
                     IF (BC(2)<1.5_EB) VM = BC(2)*VP
                     IF (BC(2)>1.5_EB) VM = FVT(2,2)*PROF(2)
                  ENDIF
            END SELECT
         ELSE
            OM => OMESH(ABS(NOM(2)))
            IF (PREDICTOR) THEN
               OM_VV => OM%VS
            ELSE
               OM_VV => OM%V
            ENDIF
            WGT = EDGE_INTERPOLATION_FACTOR(IE,1)
            IF (NOM(2)<0) VM = WGT*OM_VV(IIO(2),JJO(2),KKO(2)) + (1._EB-WGT)*OM_VV(IIO(2),JJO(2)-1,KKO(2))
            IF (NOM(2)>0) VP = WGT*OM_VV(IIO(2),JJO(2),KKO(2)) + (1._EB-WGT)*OM_VV(IIO(2),JJO(2)-1,KKO(2))
         ENDIF

         IF (NOM(1)==0 .OR. EDGE_TYPE(IE,1)/=INTERPOLATED_EDGE) THEN
            SELECT CASE(IOR(1))
               CASE(-2)
                  IWM = WALL_INDEX(ICMM, 2) 
                  IWP = WALL_INDEX(ICMP, 2) 
                  IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .OR. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
                     IF (BC(1)<1.5_EB) WP = BC(1)*WM
                     IF (BC(1)>1.5_EB) WP = FVT(1,2)
                  ENDIF
               CASE( 2)
                  IWM = WALL_INDEX(ICPM,-2) 
                  IWP = WALL_INDEX(ICPP,-2) 
                  IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .OR. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
                     IF (BC(1)<1.5_EB) WM = BC(1)*WP
                     IF (BC(1)>1.5_EB) WM = FVT(1,2)
                  ENDIF
            END SELECT
         ELSE
            OM => OMESH(ABS(NOM(1))) 
            IF (PREDICTOR) THEN
               OM_WW => OM%WS
            ELSE
               OM_WW => OM%W
            ENDIF
            WGT = EDGE_INTERPOLATION_FACTOR(IE,2) 
            IF (NOM(1)<0) WM = WGT*OM_WW(IIO(1),JJO(1),KKO(1)) + (1._EB-WGT)*OM_WW(IIO(1),JJO(1),KKO(1)-1)
            IF (NOM(1)>0) WP = WGT*OM_WW(IIO(1),JJO(1),KKO(1)) + (1._EB-WGT)*OM_WW(IIO(1),JJO(1),KKO(1)-1) 
         ENDIF

         MUA = .25_EB*( MU(II,JJ,KK) + MU(II,JJ+1,KK) + MU(II,JJ+1,KK+1) + MU(II,JJ,KK+1) )
         DVDZ = RDZN(KK)*(VP-VM)
         DWDY = RDYN(JJ)*(WP-WM)
         OME_E(IE) = DWDY - DVDZ
         TAU_E(IE) = MUA*(DVDZ + DWDY)
         IF (JJ==0)    WW(II,JJ,KK)   = WM
         IF (JJ==JBAR) WW(II,JJ+1,KK) = WP
         IF (KK==0)    VV(II,JJ,KK)   = VM
         IF (KK==KBAR) VV(II,JJ,KK+1) = VP
         IF (CORRECTOR .AND. JJ>0 .AND. JJ<JBAR .AND. KK>0 .AND. KK<KBAR) THEN
            W_Y(II,JJ,KK) = 0.5_EB*(WM+WP)
            V_Z(II,JJ,KK) = 0.5_EB*(VM+VP)
         ENDIF
 
      CASE(2) COMPONENT  ! Treat edges that point in the y direction (omega_y, tau_y)

         ICMM = CELL_INDEX(MAX(1,II)     ,JJ,MAX(1,KK)     )
         ICMP = CELL_INDEX(MAX(1,II)     ,JJ,MIN(KBAR,KK+1))
         ICPM = CELL_INDEX(MIN(IBAR,II+1),JJ,MAX(1,KK)     )
         ICPP = CELL_INDEX(MIN(IBAR,II+1),JJ,MIN(KBAR,KK+1))
         UP   = UU(II,JJ,KK+1)
         UM   = UU(II,JJ,KK)
         WP   = WW(II+1,JJ,KK)
         WM   = WW(II,JJ,KK)

         IF (NOM(1)==0 .OR. EDGE_TYPE(IE,1)/=INTERPOLATED_EDGE) THEN
            SELECT CASE(IOR(1))
               CASE(-3)
                  IWM = WALL_INDEX(ICMM, 3) 
                  IWP = WALL_INDEX(ICPM, 3) 
                  IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .OR. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
                     IF (BC(1)<1.5_EB) UP = BC(1)*UM
                     IF (BC(1)>1.5_EB) UP = FVT(1,1)*PROF(1)
                  ENDIF
               CASE( 3)
                  IWM = WALL_INDEX(ICMP,-3) 
                  IWP = WALL_INDEX(ICPP,-3) 
                  IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .OR. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
                     IF (BC(1)<1.5_EB) UM = BC(1)*UP
                     IF (BC(1)>1.5_EB) UM = FVT(1,1)*PROF(1)
                  ENDIF
            END SELECT
         ELSE
            OM => OMESH(ABS(NOM(1)))
            IF (PREDICTOR) THEN
               OM_UU => OM%US
            ELSE
               OM_UU => OM%U
            ENDIF
            WGT = EDGE_INTERPOLATION_FACTOR(IE,1)
            IF (NOM(1)<0) UM = WGT*OM_UU(IIO(1),JJO(1),KKO(1)) + (1._EB-WGT)*OM_UU(IIO(1)-1,JJO(1),KKO(1))
            IF (NOM(1)>0) UP = WGT*OM_UU(IIO(1),JJO(1),KKO(1)) + (1._EB-WGT)*OM_UU(IIO(1)-1,JJO(1),KKO(1))
         ENDIF

         IF (NOM(2)==0 .OR. EDGE_TYPE(IE,2)/=INTERPOLATED_EDGE) THEN
            SELECT CASE(IOR(2))
               CASE(-1)
                  IWM = WALL_INDEX(ICMM, 1) 
                  IWP = WALL_INDEX(ICMP, 1) 
                  IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .OR. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
                     IF (BC(2)<1.5_EB) WP = BC(2)*WM
                     IF (BC(2)>1.5_EB) WP = FVT(2,2)
                  ENDIF
               CASE( 1)
                  IWM = WALL_INDEX(ICPM,-1) 
                  IWP = WALL_INDEX(ICPP,-1) 
                  IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .OR. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
                     IF (BC(2)<1.5_EB) WM = BC(2)*WP
                     IF (BC(2)>1.5_EB) WM = FVT(2,2)
                  ENDIF
            END SELECT
         ELSE
            OM => OMESH(ABS(NOM(2)))
            IF (PREDICTOR) THEN
               OM_WW => OM%WS
            ELSE
               OM_WW => OM%W
            ENDIF
            WGT = EDGE_INTERPOLATION_FACTOR(IE,2)
            IF (NOM(2)<0) WM = WGT*OM_WW(IIO(2),JJO(2),KKO(2)) + (1._EB-WGT)*OM_WW(IIO(2),JJO(2),KKO(2)-1)
            IF (NOM(2)>0) WP = WGT*OM_WW(IIO(2),JJO(2),KKO(2)) + (1._EB-WGT)*OM_WW(IIO(2),JJO(2),KKO(2)-1)
         ENDIF
         
         MUA = .25_EB*( MU(II,JJ,KK) + MU(II+1,JJ,KK) + MU(II+1,JJ,KK+1) + MU(II,JJ,KK+1) )
         DUDZ = RDZN(KK)*(UP-UM)
         DWDX = RDXN(II)*(WP-WM)
         OME_E(IE) = DUDZ - DWDX
         TAU_E(IE) = MUA*(DUDZ + DWDX)
         IF (II==0)    WW(II,JJ,KK)   = WM
         IF (II==IBAR) WW(II+1,JJ,KK) = WP
         IF (KK==0)    UU(II,JJ,KK)   = UM
         IF (KK==KBAR) UU(II,JJ,KK+1) = UP
         IF (CORRECTOR .AND. II>0 .AND. II<IBAR .AND. KK>0 .AND. KK<KBAR) THEN
            U_Z(II,JJ,KK) = 0.5_EB*(UM+UP) 
            W_X(II,JJ,KK) = 0.5_EB*(WM+WP) 
         ENDIF

      CASE(3) COMPONENT   ! Treat edges that point in the z direction (omega_z, tau_z)

         ICMM = CELL_INDEX(MAX(1,II)     ,MAX(1,JJ)     ,KK)
         ICMP = CELL_INDEX(MAX(1,II)     ,MIN(JBAR,JJ+1),KK)
         ICPM = CELL_INDEX(MIN(IBAR,II+1),MAX(1,JJ)     ,KK)
         ICPP = CELL_INDEX(MIN(IBAR,II+1),MIN(JBAR,JJ+1),KK)
         UP   = UU(II,JJ+1,KK)
         UM   = UU(II,JJ,KK)
         VP   = VV(II+1,JJ,KK)
         VM   = VV(II,JJ,KK)

         IF (NOM(2)==0 .OR. EDGE_TYPE(IE,2)/=INTERPOLATED_EDGE) THEN
            SELECT CASE(IOR(2))
               CASE(-2)
                  IWM = WALL_INDEX(ICMM, 2) 
                  IWP = WALL_INDEX(ICPM, 2) 
                  IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .OR. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
                     IF (BC(2)<1.5_EB) UP = BC(2)*UM
                     IF (BC(2)>1.5_EB) UP = FVT(2,1)*PROF(2)
                  ENDIF
               CASE( 2)
                  IWM = WALL_INDEX(ICMP,-2) 
                  IWP = WALL_INDEX(ICPP,-2) 
                  IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .OR. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
                     IF (BC(2)<1.5_EB) UM = BC(2)*UP
                     IF (BC(2)>1.5_EB) UM = FVT(2,1)*PROF(2)
                  ENDIF
            END SELECT
         ELSE
            OM => OMESH(ABS(NOM(2)))
              IF (PREDICTOR) THEN
                 OM_UU => OM%US
              ELSE
                 OM_UU => OM%U
              ENDIF
              WGT = EDGE_INTERPOLATION_FACTOR(IE,1)
              IF (NOM(2)<0) UM = WGT*OM_UU(IIO(2),JJO(2),KKO(2)) + (1._EB-WGT)*OM_UU(IIO(2)-1,JJO(2),KKO(2))
              IF (NOM(2)>0) UP = WGT*OM_UU(IIO(2),JJO(2),KKO(2)) + (1._EB-WGT)*OM_UU(IIO(2)-1,JJO(2),KKO(2))
         ENDIF

         IF (NOM(1)==0 .OR. EDGE_TYPE(IE,1)/=INTERPOLATED_EDGE) THEN
            SELECT CASE(IOR(1))
               CASE(-1)
                  IWM = WALL_INDEX(ICMM, 1) 
                  IWP = WALL_INDEX(ICMP, 1) 
                  IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .OR. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
                     IF (BC(1)<1.5_EB) VP = BC(1)*VM
                     IF (BC(1)>1.5_EB) VP = FVT(1,1)*PROF(1)
                  ENDIF
               CASE( 1)
                  IWM = WALL_INDEX(ICPM,-1) 
                  IWP = WALL_INDEX(ICPP,-1) 
                  IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .OR. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
                     IF (BC(1)<1.5_EB) VM = BC(1)*VP
                     IF (BC(1)>1.5_EB) VM = FVT(1,1)*PROF(1)
                  ENDIF
            END SELECT
         ELSE
            OM => OMESH(ABS(NOM(1))) 
            IF (PREDICTOR) THEN
               OM_VV => OM%VS
            ELSE
               OM_VV => OM%V
            ENDIF
            WGT = EDGE_INTERPOLATION_FACTOR(IE,2) 
            IF (NOM(1)<0) VM = WGT*OM_VV(IIO(1),JJO(1),KKO(1)) + (1._EB-WGT)*OM_VV(IIO(1),JJO(1)-1,KKO(1)) 
            IF (NOM(1)>0) VP = WGT*OM_VV(IIO(1),JJO(1),KKO(1)) + (1._EB-WGT)*OM_VV(IIO(1),JJO(1)-1,KKO(1)) 
         ENDIF
   
         MUA = .25_EB*( MU(II,JJ,KK) + MU(II+1,JJ,KK) + MU(II+1,JJ+1,KK) + MU(II,JJ+1,KK) )
         DVDX = RDXN(II)*(VP-VM)
         DUDY = RDYN(JJ)*(UP-UM)
         OME_E(IE) = DVDX - DUDY
         TAU_E(IE) = MUA*(DVDX + DUDY)
         IF (II==0)    VV(II,JJ,KK)   = VM
         IF (II==IBAR) VV(II+1,JJ,KK) = VP
         IF (JJ==0)    UU(II,JJ,KK)   = UM
         IF (JJ==JBAR) UU(II,JJ+1,KK) = UP
         IF (CORRECTOR .AND. II>0 .AND. II<IBAR .AND. JJ>0 .AND. JJ<JBAR) THEN
            U_Y(II,JJ,KK) = 0.5_EB*(UM+UP) 
            V_X(II,JJ,KK) = 0.5_EB*(VM+VP) 
         ENDIF

   END SELECT COMPONENT

   ! For SAWTOOTH=.FALSE., zero out vorticity and strain at boundary

   IF (IOR(1)==0 .AND. IOR(2)==0) THEN
      OME_E(IE) = 0._EB
      TAU_E(IE) = 0._EB
   ENDIF

ENDDO EDGE_LOOP

! Store cell node averages of the velocity components in UVW_GHOST for use in Smokeview

IF (CORRECTOR) THEN
   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            IC = CELL_INDEX(I,J,K) 
            IF (IC==0) CYCLE
            IF (U_Y(I,J,K)  >-1.E5_EB) UVW_GHOST(IC,1) = U_Y(I,J,K) 
            IF (U_Z(I,J,K)  >-1.E5_EB) UVW_GHOST(IC,1) = U_Z(I,J,K) 
            IF (V_X(I,J,K)  >-1.E5_EB) UVW_GHOST(IC,2) = V_X(I,J,K) 
            IF (V_Z(I,J,K)  >-1.E5_EB) UVW_GHOST(IC,2) = V_Z(I,J,K) 
            IF (W_X(I,J,K)  >-1.E5_EB) UVW_GHOST(IC,3) = W_X(I,J,K) 
            IF (W_Y(I,J,K)  >-1.E5_EB) UVW_GHOST(IC,3) = W_Y(I,J,K) 
         ENDDO
      ENDDO
   ENDDO
ENDIF

TUSED(4,NM)=TUSED(4,NM)+SECOND()-TNOW
END SUBROUTINE VELOCITY_BC 
 
 
 
SUBROUTINE MATCH_VELOCITY(NM)

! Force normal component of velocity to match at interpolated boundaries

INTEGER  :: NOM,II,JJ,KK,IOR,IW,IIO,JJO,KKO
INTEGER, INTENT(IN) :: NM
REAL(EB) :: UU_AVG,VV_AVG,WW_AVG,TNOW,DA_OTHER,UU_OTHER,VV_OTHER,WW_OTHER
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,OM_UU,OM_VV,OM_WW
TYPE (OMESH_TYPE), POINTER :: OM
TYPE (MESH_TYPE), POINTER :: M2

IF (SOLID_PHASE_ONLY) RETURN
IF (NMESHES==1 .AND. PERIODIC_TEST==0) RETURN

TNOW = SECOND()

! Assign local variable names

CALL POINT_TO_MESH(NM)

! Point to the appropriate velocity field

IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
ELSE
   UU => US
   VV => VS
   WW => WS
ENDIF

! Zero out D_CORR

IF (PREDICTOR) DS_CORR = 0._EB
IF (CORRECTOR) D_CORR  = 0._EB

! Loop over all cell edges and determine the appropriate velocity BCs

EXTERNAL_WALL_LOOP: DO IW=1,NEWC

   IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP

   II  = IJKW( 1,IW)
   JJ  = IJKW( 2,IW)
   KK  = IJKW( 3,IW)
   IOR = IJKW( 4,IW)
   NOM = IJKW( 9,IW)
   OM => OMESH(NOM)
   M2 => MESHES(NOM)
   DA_OTHER = 0._EB

   SELECT CASE(ABS(IOR))
      CASE(1)
         IF (PREDICTOR) OM_UU => OM%U
         IF (CORRECTOR) OM_UU => OM%US
         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  DA_OTHER = DA_OTHER + M2%DY(JJO)*M2%DZ(KKO)
               ENDDO
            ENDDO
         ENDDO
      CASE(2)
         IF (PREDICTOR) OM_VV => OM%V
         IF (CORRECTOR) OM_VV => OM%VS
         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  DA_OTHER = DA_OTHER + M2%DX(IIO)*M2%DZ(KKO)
               ENDDO
            ENDDO
         ENDDO
      CASE(3)
         IF (PREDICTOR) OM_WW => OM%W
         IF (CORRECTOR) OM_WW => OM%WS
         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  DA_OTHER = DA_OTHER + M2%DX(IIO)*M2%DY(JJO)
               ENDDO
            ENDDO
         ENDDO
   END SELECT

   SELECT CASE(IOR)
      CASE( 1)
         UU_OTHER = 0._EB
         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  UU_OTHER = UU_OTHER + OM_UU(IIO,JJO,KKO)*M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         UU_AVG = 0.5_EB*(UU(0,JJ,KK) + UU_OTHER)
         IF (PREDICTOR) DS_CORR(IW) = (UU_AVG-UU(0,JJ,KK))*RDX(1)
         IF (CORRECTOR) D_CORR(IW) = DS_CORR(IW) + 0.5*(UU_AVG-UU(0,JJ,KK))*RDX(1)
         UVW_SAVE(IW) = UU(0,JJ,KK)
         UU(0,JJ,KK)    = UU_AVG

      CASE(-1)
         UU_OTHER = 0._EB
         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  UU_OTHER = UU_OTHER + OM_UU(IIO-1,JJO,KKO)*M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         UU_AVG = 0.5_EB*(UU(IBAR,JJ,KK) + UU_OTHER)
         IF (PREDICTOR) DS_CORR(IW) = -(UU_AVG-UU(IBAR,JJ,KK))*RDX(IBAR)
         IF (CORRECTOR) D_CORR(IW) = DS_CORR(IW) - 0.5*(UU_AVG-UU(IBAR,JJ,KK))*RDX(IBAR)
         UVW_SAVE(IW) = UU(IBAR,JJ,KK)
         UU(IBAR,JJ,KK) = UU_AVG

      CASE( 2)
         VV_OTHER = 0._EB
         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  VV_OTHER = VV_OTHER + OM_VV(IIO,JJO,KKO)*M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         VV_AVG = 0.5_EB*(VV(II,0,KK) + VV_OTHER)
         IF (PREDICTOR) DS_CORR(IW) = (VV_AVG-VV(II,0,KK))*RDY(1)
         IF (CORRECTOR) D_CORR(IW) = DS_CORR(IW) + 0.5*(VV_AVG-VV(II,0,KK))*RDY(1)
         UVW_SAVE(IW) = VV(II,0,KK)
         VV(II,0,KK)    = VV_AVG

      CASE(-2)
         VV_OTHER = 0._EB
         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  VV_OTHER = VV_OTHER + OM_VV(IIO,JJO-1,KKO)*M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         VV_AVG = 0.5_EB*(VV(II,JBAR,KK) + VV_OTHER)
         IF (PREDICTOR) DS_CORR(IW) = -(VV_AVG-VV(II,JBAR,KK))*RDY(JBAR)
         IF (CORRECTOR) D_CORR(IW) = DS_CORR(IW) - 0.5*(VV_AVG-VV(II,JBAR,KK))*RDY(JBAR)
         UVW_SAVE(IW) = VV(II,JBAR,KK)
         VV(II,JBAR,KK) = VV_AVG

      CASE( 3)
         WW_OTHER = 0._EB
         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  WW_OTHER = WW_OTHER + OM_WW(IIO,JJO,KKO)*M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         WW_AVG = 0.5_EB*(WW(II,JJ,0) + WW_OTHER)
         IF (PREDICTOR) DS_CORR(IW) = (WW_AVG-WW(II,JJ,0))*RDZ(1)
         IF (CORRECTOR) D_CORR(IW) = DS_CORR(IW) + 0.5*(WW_AVG-WW(II,JJ,0))*RDZ(1)
         UVW_SAVE(IW) = WW(II,JJ,0)
         WW(II,JJ,0)    = WW_AVG

      CASE(-3)
         WW_OTHER = 0._EB
         DO KKO=IJKW(12,IW),IJKW(15,IW)
            DO JJO=IJKW(11,IW),IJKW(14,IW)
               DO IIO=IJKW(10,IW),IJKW(13,IW)
                  WW_OTHER = WW_OTHER + OM_WW(IIO,JJO,KKO-1)*M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         WW_AVG = 0.5_EB*(WW(II,JJ,KBAR) + WW_OTHER)
         IF (PREDICTOR) DS_CORR(IW) = -(WW_AVG-WW(II,JJ,KBAR))*RDZ(KBAR)
         IF (CORRECTOR) D_CORR(IW) = DS_CORR(IW) - 0.5*(WW_AVG-WW(II,JJ,KBAR))*RDZ(KBAR)
         UVW_SAVE(IW) = WW(II,JJ,KBAR)
         WW(II,JJ,KBAR) = WW_AVG
   END SELECT

ENDDO EXTERNAL_WALL_LOOP

TUSED(4,NM)=TUSED(4,NM)+SECOND()-TNOW
END SUBROUTINE MATCH_VELOCITY



SUBROUTINE CHECK_STABILITY(NM)
 
! Checks the Courant and Von Neumann stability criteria, and if necessary, reduces the time step accordingly
 
REAL(EB) :: UODX,VODY,WODZ,UVW,UVWMAX,R_DX2,MU_MAX,MUTRM,DMAX,RDMAX
INTEGER  :: NM,I,J,K
 
CHANGE_TIME_STEP(NM) = .FALSE.
UVWMAX = 0._EB
VN     = 0._EB
MUTRM  = 1.E-9_EB
R_DX2  = 1.E-9_EB
 
! Determine max CFL number from all grid cells

SELECT_VELOCITY_NORM: SELECT CASE (CFL_VELOCITY_NORM)

   CASE(0)
      DO K=0,KBAR
         DO J=0,JBAR
            DO I=0,IBAR
               UODX = ABS(US(I,J,K))*RDXN(I)
               VODY = ABS(VS(I,J,K))*RDYN(J)
               WODZ = ABS(WS(I,J,K))*RDZN(K)
               UVW  = MAX(UODX,VODY,WODZ)
               IF (UVW>=UVWMAX) THEN
                  UVWMAX = UVW 
                  ICFL=I
                  JCFL=J
                  KCFL=K
               ENDIF
            ENDDO 
         ENDDO   
      ENDDO
      
   CASE(1)
      DO K=0,KBAR
         DO J=0,JBAR
            DO I=0,IBAR
               UVW = (ABS(US(I,J,K)) + ABS(VS(I,J,K)) + ABS(WS(I,J,K)))**2
               IF (UVW>=UVWMAX) THEN
                  UVWMAX = UVW 
                  ICFL=I
                  JCFL=J
                  KCFL=K
               ENDIF
            ENDDO 
         ENDDO   
      ENDDO
      
   CASE(2)
      DO K=0,KBAR
         DO J=0,JBAR
            DO I=0,IBAR
               UVW = US(I,J,K)**2 + VS(I,J,K)**2 + WS(I,J,K)**2
               IF (UVW>=UVWMAX) THEN
                  UVWMAX = UVW 
                  ICFL=I
                  JCFL=J
                  KCFL=K
               ENDIF
            ENDDO 
         ENDDO   
      ENDDO
      
END SELECT SELECT_VELOCITY_NORM

! Find minimum time step allowed by divergence constraint
RDMAX = HUGE(1._EB)
IF (CFL_VELOCITY_NORM>0) THEN
   UVWMAX = SQRT(UVWMAX)*MAX(RDXN(ICFL),RDYN(JCFL),RDZN(KCFL))
   DMAX = UVWMAX + MAXVAL(D)
   IF (DMAX>0._EB) RDMAX = CFL_MAX/DMAX
ENDIF

CFL = DT*UVWMAX
 
! Determine max Von Neumann Number for fine grid calcs
 
PARABOLIC_IF: IF (DNS .OR. CELL_SIZE<0.005_EB) THEN
 
   INCOMPRESSIBLE_IF: IF (ISOTHERMAL .AND. N_SPECIES==0) THEN
      IF (TWO_D) THEN
         R_DX2 = 1._EB/DXMIN**2 + 1._EB/DZMIN**2
      ELSE
         R_DX2 = 1._EB/DXMIN**2 + 1._EB/DYMIN**2 + 1._EB/DZMIN**2
      ENDIF
      MUTRM = MAX(RPR,RSC)*SPECIES(0)%MU(NINT(0.1_EB*TMPA))/RHOA
   ELSE INCOMPRESSIBLE_IF
      MU_MAX = 0._EB   
      DO K=1,KBAR
         DO J=1,JBAR
            IILOOP: DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE IILOOP
               IF (MU(I,J,K)>=MU_MAX) THEN
                  MU_MAX = MU(I,J,K)
                  I_VN=I
                  J_VN=J
                  K_VN=K
               ENDIF
            ENDDO IILOOP
         ENDDO  
      ENDDO  
      IF (TWO_D) THEN
         R_DX2 = RDX(I_VN)**2 + RDZ(K_VN)**2
      ELSE
         R_DX2 = RDX(I_VN)**2 + RDY(J_VN)**2 + RDZ(K_VN)**2
      ENDIF
      MUTRM = MAX(RPR,RSC)*MU_MAX/RHOS(I_VN,J_VN,K_VN)
   ENDIF INCOMPRESSIBLE_IF
 
   VN = DT*2._EB*R_DX2*MUTRM
 
ENDIF PARABOLIC_IF
 
! Adjust time step size if necessary
 
IF ((CFL<CFL_MAX.AND.VN<VN_MAX.AND.DT<RDMAX) .OR. LOCK_TIME_STEP) THEN
   DTNEXT = DT
   IF (CFL<=CFL_MIN .AND. VN<VN_MIN .AND. .NOT.LOCK_TIME_STEP) DTNEXT = MIN(1.1_EB*DT,DTINT)
ELSE
   IF (UVWMAX==0._EB) UVWMAX = 1._EB
   DT = 0.9_EB*MIN( CFL_MAX/UVWMAX , VN_MAX/(2._EB*R_DX2*MUTRM), RDMAX )
   CHANGE_TIME_STEP(NM) = .TRUE.
ENDIF
 
END SUBROUTINE CHECK_STABILITY
 
 
SUBROUTINE BAROCLINIC_CORRECTION
 
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,HQS,RTRM,RHOP
REAL(EB) :: RHO_AVG_OLD,RRAT,U2,V2,W2
INTEGER  :: I,J,K
 
HQS  => WORK1
RTRM => WORK2
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   RHOP=>RHO
ELSE
   UU => US
   VV => VS
   WW => WS
   RHOP=>RHOS
ENDIF
 
RHO_AVG_OLD = RHO_AVG
RHO_AVG = 2._EB*RHO_LOWER*RHO_UPPER/(RHO_LOWER+RHO_UPPER)
RRAT = RHO_AVG_OLD/RHO_AVG
RTRM = (1._EB-RHO_AVG/RHOP)*RRAT
 
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         U2 = 0.25_EB*(UU(I,J,K)+UU(I-1,J,K))**2
         V2 = 0.25_EB*(VV(I,J,K)+VV(I,J-1,K))**2
         W2 = 0.25_EB*(WW(I,J,K)+WW(I,J,K-1))**2
         HQS(I,J,K) = 0.5_EB*(U2+V2+W2)
      ENDDO
   ENDDO
ENDDO
 
DO K=1,KBAR
   DO J=1,JBAR
      U2 = (1.5_EB*UU(0,J,K)-0.5_EB*UU(1,J,K))**2
      V2 = 0.25_EB*(VV(1,J,K)+VV(1,J-1,K))**2
      W2 = 0.25_EB*(WW(1,J,K)+WW(1,J,K-1))**2
      HQS(0,J,K) = MIN(0.5_EB*(U2+V2+W2),10._EB*DX(1)+HQS(1,J,K))
      U2 = (1.5_EB*UU(IBAR,J,K)-0.5_EB*UU(IBM1,J,K))**2
      V2 = 0.25_EB*(VV(IBAR,J,K)+VV(IBAR,J-1,K))**2
      W2 = 0.25_EB*(WW(IBAR,J,K)+WW(IBAR,J,K-1))**2
      HQS(IBP1,J,K) = MIN(0.5_EB*(U2+V2+W2),10._EB*DX(IBAR)+HQS(IBAR,J,K))
   ENDDO
ENDDO
 
IF (.NOT.TWO_D) THEN
   DO K=1,KBAR
      DO I=1,IBAR
         U2 = 0.25_EB*(UU(I,1,K)+UU(I-1,1,K))**2
         V2 = (1.5_EB*VV(I,0,K)-0.5_EB*VV(I,1,K))**2
         W2 = 0.25_EB*(WW(I,1,K)+WW(I,1,K-1))**2
         HQS(I,0,K) = MIN(0.5_EB*(U2+V2+W2),10._EB*DY(1)+HQS(I,1,K))
         U2 = 0.25_EB*(UU(I,JBAR,K)+UU(I-1,JBAR,K))**2
         V2 = (1.5_EB*VV(I,JBAR,K)-0.5_EB*VV(I,JBM1,K))**2
         W2 = 0.25_EB*(WW(I,JBAR,K)+WW(I,JBAR,K-1))**2
         HQS(I,JBP1,K) = MIN(0.5_EB*(U2+V2+W2),10._EB*DY(JBAR)+HQS(I,JBAR,K))
      ENDDO
   ENDDO
ENDIF
 
DO J=1,JBAR
   DO I=1,IBAR
      U2 = 0.25_EB*(UU(I,J,1)+UU(I-1,J,1))**2
      V2 = 0.25_EB*(VV(I,J,1)+VV(I,J-1,1))**2
      W2 = (1.5_EB*WW(I,J,0)-0.5_EB*WW(I,J,1))**2
      HQS(I,J,0) = MIN(0.5_EB*(U2+V2+W2),10._EB*DZ(1)+HQS(I,J,1))
      U2 = 0.25_EB*(UU(I,J,KBAR)+UU(I-1,J,KBAR))**2
      V2 = 0.25_EB*(VV(I,J,KBAR)+VV(I,J-1,KBAR))**2
      W2 = (1.5_EB*WW(I,J,KBAR)-0.5_EB*WW(I,J,KBM1))**2
      HQS(I,J,KBP1) = MIN(0.5_EB*(U2+V2+W2),10._EB*DZ(KBAR)+HQS(I,J,KBAR))
   ENDDO
ENDDO
 
DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         FVX(I,J,K) = FVX(I,J,K) - 0.5_EB*(RTRM(I+1,J,K)+RTRM(I,J,K))*(H(I+1,J,K)-H(I,J,K)-HQS(I+1,J,K)+HQS(I,J,K))*RDXN(I)
      ENDDO
   ENDDO
ENDDO
 
IF (.NOT.TWO_D) THEN
   DO K=1,KBAR
      DO J=0,JBAR
         DO I=1,IBAR
            FVY(I,J,K) = FVY(I,J,K) - 0.5_EB*(RTRM(I,J+1,K)+RTRM(I,J,K))*(H(I,J+1,K)-H(I,J,K)-HQS(I,J+1,K)+HQS(I,J,K))*RDYN(J)
         ENDDO
      ENDDO
   ENDDO
ENDIF
 
DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         FVZ(I,J,K) = FVZ(I,J,K) - 0.5_EB*(RTRM(I,J,K+1)+RTRM(I,J,K))*(H(I,J,K+1)-H(I,J,K)-HQS(I,J,K+1)+HQS(I,J,K))*RDZN(K)
      ENDDO
   ENDDO
ENDDO
 
END SUBROUTINE BAROCLINIC_CORRECTION

SUBROUTINE GET_REV_velo(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') velorev(INDEX(velorev,':')+1:LEN_TRIM(velorev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') velodate

END SUBROUTINE GET_REV_velo
 
END MODULE VELO
