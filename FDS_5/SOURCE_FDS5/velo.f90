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
 
SUBROUTINE COMPUTE_VELOCITY_FLUX(T,NM)

REAL(EB), INTENT(IN) :: T
REAL(EB) :: TNOW
INTEGER, INTENT(IN) :: NM

IF (SOLID_PHASE_ONLY) RETURN

TNOW = SECOND()

IF (DNS .AND. ISOTHERMAL .AND. N_SPECIES==0) THEN
   CALL VELOCITY_FLUX_ISOTHERMAL(NM)
ELSE
   IF (.NOT.CYLINDRICAL) CALL VELOCITY_FLUX(T,NM)
   IF (     CYLINDRICAL) CALL VELOCITY_FLUX_CYLINDRICAL(T,NM)
ENDIF

TUSED(4,NM) = TUSED(4,NM) + SECOND() - TNOW
END SUBROUTINE COMPUTE_VELOCITY_FLUX

SUBROUTINE VELOCITY_FLUX(T,NM)
! Compute convective and diffusive terms
USE PHYSICAL_FUNCTIONS, ONLY: GET_MU2
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
INTEGER, INTENT(IN) :: NM
REAL(EB) :: T,MUX,MUY,MUZ,UP,UM,VP,VM,WP,WM,VTRM, &
            DTXYDY,DTXZDZ,DTYZDZ,DTXYDX,DTXZDX,DTYZDY, &
            DUDX,DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY, &
            VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX,PMDT,MPDT, &
            S2,C,CDXDYDZTT,AH,RRHO,GX,GY,GZ,Z_2, &
            TXXP,TXXM,TYYP,TYYM,TZZP,TZZM,DTXXDX,DTYYDY,DTZZDZ, &
            EPSUP,EPSUM,EPSVP,EPSVM,EPSWP,EPSWM,MU_SUM,DUMMY=0._EB
INTEGER :: II,JJ,KK,I,J,K,IW,IIG,JJG,KKG,ITMP,N,IE,IEC
REAL(EB), POINTER, DIMENSION(:,:,:) :: TXY,TXZ,TYZ,OMX,OMY,OMZ, UU,VV,WW,RHOP,DP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP
 
CALL POINT_TO_MESH(NM)
 
IF (PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   DP => D  
   RHOP => RHO
   IF (N_SPECIES > 0) YYP => YY
ELSE
   UU => US
   VV => VS
   WW => WS
   DP => DS
   RHOP => RHOS
   IF (N_SPECIES > 0) YYP => YYS
ENDIF

TXY => WORK1
TXZ => WORK2
TYZ => WORK3
OMX => WORK4
OMY => WORK5
OMZ => WORK6
 
CALC_MU: IF (PREDICTOR) THEN    
 
   LES_VS_DNS: IF (LES) THEN   ! Smagorinsky model (LES)

      C = CSMAG**2
      IF (EVACUATION_ONLY(NM)) C = (0.9_EB)**2
      KLOOP: DO K=1,KBAR
         JLOOP: DO J=1,JBAR
            ILOOP: DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ILOOP
               IF (.NOT.TWO_D) THEN
                  CDXDYDZTT = C*(DX(I)*DY(J)*DZ(K))**TWTH
               ELSE
                  CDXDYDZTT = C*DX(I)*DZ(K)
               ENDIF
               DUDX = RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
               DVDY = RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
               DWDZ = RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
               DUDY = 0.25_EB*RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K)+UU(I-1,J+1,K)-UU(I-1,J-1,K))
               DVDX = 0.25_EB*RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K)+VV(I+1,J-1,K)-VV(I-1,J-1,K))
               DWDX = 0.25_EB*RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
               DWDY = 0.25_EB*RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K)+WW(I,J+1,K-1)-WW(I,J-1,K-1))
               DUDZ = 0.25_EB*RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1)) 
               DVDZ = 0.25_EB*RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1)+VV(I,J-1,K+1)-VV(I,J-1,K-1)) 
               S2   = 2._EB*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ ) + (DUDY+DVDX)**2 + (DUDZ+DWDX)**2 + &
                      (DVDZ+DWDY)**2-TWTH*DP(I,J,K)**2
               S2   = MAX(0._EB,S2)
               ITMP = 0.1_EB*TMP(I,J,K)
               MU(I,J,K) = MAX(SPECIES(0)%MU(ITMP), RHOP(I,J,K)*CDXDYDZTT*SQRT(S2))
            ENDDO ILOOP
         ENDDO JLOOP
      ENDDO KLOOP
 
   ELSE LES_VS_DNS ! DNS viscosity
 
      MIXTURE_FRACTION_DNS: IF (.NOT.MIXTURE_FRACTION) THEN
         DO K=1,KBAR
            DO J=1,JBAR
               IVLOOP: DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE IVLOOP
                  ITMP = 0.1_EB*TMP(I,J,K)
                  MU_SUM = SPECIES(0)%MU(ITMP)
                  DO N=1,N_SPECIES
                     MU_SUM = MU_SUM + YYP(I,J,K,N)*(SPECIES(N)%MU(ITMP)-SPECIES(0)%MU(ITMP))
                  ENDDO
                  MU(I,J,K) = MU_SUM
               ENDDO IVLOOP
            ENDDO
         ENDDO
      ELSE MIXTURE_FRACTION_DNS
         Z_2 = 0._EB
         DO K=1,KBAR
            DO J=1,JBAR
               IVLOOP2: DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE IVLOOP2
                  ITMP = 0.1_EB*TMP(I,J,K)
                  IF(CO_PRODUCTION) THEN
                     CALL GET_MU2(YY(I,J,K,I_FUEL),YY(I,J,K,I_PROG_CO),YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),MU(I,J,K),ITMP)
                  ELSE
                     CALL GET_MU2(YY(I,J,K,I_FUEL),Z_2,YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),MU(I,J,K),ITMP)                  
                  ENDIF
                  MU_SUM = MU(I,J,K)*(1-Y_SUM(I,J,K))
                  DO N=1,N_SPECIES
                     IF(SPECIES(N)%MODE/=MIXTURE_FRACTION_SPECIES) MU_SUM = MU_SUM + YYP(I,J,K,N)*(SPECIES(N)%MU(ITMP)-MU(I,J,K))
                  ENDDO
                  MU(I,J,K) = MU_SUM
               ENDDO IVLOOP2
            ENDDO
         ENDDO
      ENDIF MIXTURE_FRACTION_DNS
 
   ENDIF LES_VS_DNS
 
! Mirror viscosity into solids
 
   WALL_LOOP: DO IW=1,NWC
      IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_LOOP
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      MU(II,JJ,KK) = MU(IIG,JJG,KKG)
   ENDDO WALL_LOOP
    
   MU(0,0:JBP1,0)       = MU(   1,0:JBP1,1)
   MU(IBP1,0:JBP1,0)    = MU(IBAR,0:JBP1,1)
   MU(IBP1,0:JBP1,KBP1) = MU(IBAR,0:JBP1,KBAR)
   MU(0,0:JBP1,KBP1)    = MU(   1,0:JBP1,KBAR)
   MU(0:IBP1,0,0)       = MU(0:IBP1,   1,1)
   MU(0:IBP1,JBP1,0)    = MU(0:IBP1,JBAR,1)
   MU(0:IBP1,JBP1,KBP1) = MU(0:IBP1,JBAR,KBAR)
   MU(0:IBP1,0,KBP1)    = MU(0:IBP1,   1,KBAR)
   MU(0,0,0:KBP1)       = MU(   1,   1,0:KBP1)
   MU(IBP1,0,0:KBP1)    = MU(IBAR,   1,0:KBP1)
   MU(IBP1,JBP1,0:KBP1) = MU(IBAR,JBAR,0:KBP1)
   MU(0,JBP1,0:KBP1)    = MU(   1,JBAR,0:KBP1)
 
ENDIF CALC_MU
 
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
 
GX  = EVALUATE_RAMP(T,DUMMY,I_RAMP_GX)*GVEC(1)
GY  = EVALUATE_RAMP(T,DUMMY,I_RAMP_GY)*GVEC(2)
GZ  = EVALUATE_RAMP(T,DUMMY,I_RAMP_GZ)*GVEC(3)
 
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
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP 
USE PHYSICAL_FUNCTIONS, ONLY: GET_MU2
! Compute convective and diffusive terms for 2D axisymmetric
 
REAL(EB) :: T,DMUDX
INTEGER :: I0
INTEGER, INTENT(IN) :: NM
REAL(EB) :: MUY,UP,UM,WP,WM,VTRM, &
            DTXZDZ,DTXZDX, &
            DUDX,DWDZ,DUDZ,DWDX,WOMY,UOMY,PMDT,MPDT, &
            S2,C,CDXDYDZTT,AH,RRHO,GX,GZ,&
            TXXP,TXXM,TZZP,TZZM,DTXXDX,DTZZDZ, &
            EPSUP,EPSUM,EPSWP,EPSWM,MU_SUM,DUMMY=0._EB
INTEGER :: II,JJ,KK,I,J,K,IW,IIG,JJG,KKG,ITMP,N,IE,IEC
REAL(EB), POINTER, DIMENSION(:,:,:) :: TXZ,OMY,UU,WW,RHOP,DP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP
 
CALL POINT_TO_MESH(NM)
 
IF (PREDICTOR) THEN
   UU => U
   WW => W
   DP => D  
   RHOP => RHO
   IF (N_SPECIES>0) YYP => YY
ELSE
   UU => US
   WW => WS
   DP => DS
   RHOP => RHOS
   IF (N_SPECIES>0) YYP => YYS
ENDIF
 
TXZ => WORK2
OMY => WORK5
 
CALC_MU: IF (PREDICTOR) THEN    
 
   LES_VS_DNS: IF (LES) THEN   ! Smagorinsky model (LES)

      C = CSMAG**2
      J = 1
      KLOOP: DO K=1,KBAR
         ILOOP: DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ILOOP
            CDXDYDZTT = C*DX(I)*DZ(K)
            DUDX = RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
            DWDZ = RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
            DWDX = 0.25_EB*RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
            DUDZ = 0.25_EB*RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1))
            S2   = 2._EB*(DUDX*DUDX     +  DWDZ*DWDZ ) +  (DUDZ+DWDX)**2 -  TWTH*DP(I,J,K)**2
            S2   = MAX(0._EB,S2)
            ITMP = 0.1_EB*TMP(I,J,K)
            MU(I,J,K) = MAX(SPECIES(0)%MU(ITMP), RHOP(I,J,K)*CDXDYDZTT*SQRT(S2))
         ENDDO ILOOP
      ENDDO KLOOP

   ELSE LES_VS_DNS           ! DNS viscosity

      MIXTURE_FRACTION_DNS: IF (.NOT.MIXTURE_FRACTION) THEN
         J = 1
         DO K=1,KBAR
            IVLOOP: DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE IVLOOP
               ITMP = 0.1_EB*TMP(I,J,K)
               MU_SUM = SPECIES(0)%MU(ITMP)
               DO N=1,N_SPECIES
                  MU_SUM = MU_SUM + YYP(I,J,K,N)*(SPECIES(N)%MU(ITMP)-SPECIES(0)%MU(ITMP))
               ENDDO
               MU(I,J,K) = MU_SUM
            ENDDO IVLOOP
         ENDDO
      ELSE MIXTURE_FRACTION_DNS
         J = 1
         DO K=1,KBAR
            IVLOOP2: DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE IVLOOP2
               ITMP = 0.1_EB*TMP(I,J,K)
               IF(CO_PRODUCTION) THEN
                  CALL GET_MU2(YY(I,J,K,I_FUEL),YY(I,J,K,I_PROG_CO),YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),MU(I,J,K),ITMP)
               ELSE
                  CALL GET_MU2(YY(I,J,K,I_FUEL),0._EB,YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),MU(I,J,K),ITMP)                  
               ENDIF
            ENDDO IVLOOP2
         ENDDO
      ENDIF MIXTURE_FRACTION_DNS
 
   ENDIF LES_VS_DNS
 
! Mirror viscosity into solids
 
   DO IW=1,NWC
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      MU(II,JJ,KK) = MU(IIG,JJG,KKG)
   ENDDO
 
ENDIF CALC_MU
 
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
 
!!! IF (.NOT.OBSTRUCTION(OBST_INDEX_W(IW))%SAWTOOTH) CYCLE WALL_LOOP  ! Allow velocity flux through surface of no SAWTOOTH obsts.

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
REAL(EB) :: BC,MUA,FVT,UP,UM,VP,VM,WP,WM,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY,TSI,WGT,TNOW,PROF
INTEGER  :: I,J,K,IBC,NOM1,NOM2,IIO1,IIO2,JJO1,JJO2,KKO1,KKO2,IE,II,JJ,KK,IEC,IOR,IWM,IWP,ICMM,ICMP,ICPM,ICPP,IC
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

   II   = IJKE( 1,IE)
   JJ   = IJKE( 2,IE)
   KK   = IJKE( 3,IE)
   IEC  = IJKE( 4,IE)
   IBC  = IJKE( 5,IE)
   IOR  = IJKE( 6,IE)
   NOM1 = IJKE( 7,IE)
   IIO1 = IJKE( 8,IE)
   JJO1 = IJKE( 9,IE)
   KKO1 = IJKE(10,IE)
   NOM2 = IJKE(11,IE)
   IIO2 = IJKE(12,IE)
   JJO2 = IJKE(13,IE)
   KKO2 = IJKE(14,IE)
   SF  => SURFACE(IBC)
   BC  = SF%SLIP_FACTOR
   IF (BC>1.5_EB) THEN
      IF (SF%T_IGN==T_BEGIN .AND. SF%RAMP_INDEX(TIME_VELO)>=1) THEN
         TSI = T
      ELSE
         TSI=T-SF%T_IGN
      ENDIF
      PROF = 1._EB
      FVT = EVALUATE_RAMP(TSI,SF%TAU(TIME_VELO),SF%RAMP_INDEX(TIME_VELO))
      IF (SF%PROFILE==ATMOSPHERIC) PROF=((ZC(KK)-GROUND_LEVEL)/SF%Z0)**SF%PLE
   ENDIF
 
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

         IF (NOM1==NM) THEN
            IWM = WALL_INDEX(ICMM, 3) 
            IWP = WALL_INDEX(ICPM, 3) 
            IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
               IF (BC<1.5_EB) VP = BC*VM
               IF (BC>1.5_EB) VP = FVT*SF%VEL_T(2)*PROF
            ENDIF
         ENDIF
         IF (NOM1==-NM) THEN
            IWM = WALL_INDEX(ICMP,-3) 
            IWP = WALL_INDEX(ICPP,-3) 
            IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
               IF (BC<1.5_EB) VM = BC*VP
               IF (BC>1.5_EB) VM = FVT*SF%VEL_T(2)*PROF
            ENDIF
         ENDIF

         IF (NOM2==NM) THEN
            IWM = WALL_INDEX(ICMM, 2) 
            IWP = WALL_INDEX(ICMP, 2) 
            IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
               IF (BC<1.5_EB) WP = BC*WM
               IF (BC>1.5_EB) WP = FVT*SF%VEL_T(2)
            ENDIF
         ENDIF
         IF (NOM2==-NM) THEN
            IWM = WALL_INDEX(ICPM,-2) 
            IWP = WALL_INDEX(ICPP,-2) 
            IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
               IF (BC<1.5_EB) WM = BC*WP
               IF (BC>1.5_EB) WM = FVT*SF%VEL_T(2)
            ENDIF
         ENDIF

         IF (ABS(NOM1)/=NM .AND. NOM1/=0) THEN
            OM => OMESH(ABS(NOM1)) 
            IF (PREDICTOR) THEN
               OM_VV => OM%VS
            ELSE
               OM_VV => OM%V
            ENDIF
            WGT = EDGE_INTERPOLATION_FACTOR(IE,1) 
            IF (NOM1<0) VM = WGT*OM_VV(IIO1,JJO1,KKO1) + (1._EB-WGT)*OM_VV(IIO1,JJO1-1,KKO1) 
            IF (NOM1>0) VP = WGT*OM_VV(IIO1,JJO1,KKO1) + (1._EB-WGT)*OM_VV(IIO1,JJO1-1,KKO1) 
         ENDIF
         IF (ABS(NOM2)/=NM .AND. NOM2/=0) THEN
            OM => OMESH(ABS(NOM2)) 
            IF (PREDICTOR) THEN
               OM_WW => OM%WS
            ELSE
               OM_WW => OM%W
            ENDIF
            WGT = EDGE_INTERPOLATION_FACTOR(IE,2) 
            IF (NOM2<0) WM = WGT*OM_WW(IIO2,JJO2,KKO2) + (1._EB-WGT)*OM_WW(IIO2,JJO2,KKO2-1)
            IF (NOM2>0) WP = WGT*OM_WW(IIO2,JJO2,KKO2) + (1._EB-WGT)*OM_WW(IIO2,JJO2,KKO2-1) 
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

         IF (NOM1==NM) THEN
            IWM = WALL_INDEX(ICMM, 3) 
            IWP = WALL_INDEX(ICPM, 3) 
            IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
               IF (BC<1.5_EB) UP = BC*UM
               IF (BC>1.5_EB) UP = FVT*SF%VEL_T(1)*PROF
            ENDIF
         ENDIF
         IF (NOM1==-NM) THEN
            IWM = WALL_INDEX(ICMP,-3) 
            IWP = WALL_INDEX(ICPP,-3) 
            IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
               IF (BC<1.5_EB) UM = BC*UP
               IF (BC>1.5_EB) UM = FVT*SF%VEL_T(1)*PROF
            ENDIF
         ENDIF

         IF (NOM2==NM) THEN
            IWM = WALL_INDEX(ICMM, 1) 
            IWP = WALL_INDEX(ICMP, 1) 
            IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
               IF (BC<1.5_EB) WP = BC*WM
               IF (BC>1.5_EB) WP = FVT*SF%VEL_T(2)
            ENDIF
         ENDIF
         IF (NOM2==-NM) THEN
            IWM = WALL_INDEX(ICPM,-1) 
            IWP = WALL_INDEX(ICPP,-1) 
            IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
               IF (BC<1.5_EB) WM = BC*WP
               IF (BC>1.5_EB) WM = FVT*SF%VEL_T(2)
            ENDIF
         ENDIF

         IF (ABS(NOM1)/=NM .AND. NOM1/=0) THEN
            OM => OMESH(ABS(NOM1)) 
            IF (PREDICTOR) THEN
               OM_UU => OM%US
            ELSE
               OM_UU => OM%U
            ENDIF
            WGT = EDGE_INTERPOLATION_FACTOR(IE,1) 
            IF (NOM1<0) UM = WGT*OM_UU(IIO1,JJO1,KKO1) + (1._EB-WGT)*OM_UU(IIO1-1,JJO1,KKO1)
            IF (NOM1>0) UP = WGT*OM_UU(IIO1,JJO1,KKO1) + (1._EB-WGT)*OM_UU(IIO1-1,JJO1,KKO1)
         ENDIF
         IF (ABS(NOM2)/=NM .AND. NOM2/=0) THEN
            OM => OMESH(ABS(NOM2)) 
            IF (PREDICTOR) THEN
               OM_WW => OM%WS
            ELSE
               OM_WW => OM%W
            ENDIF
            WGT = EDGE_INTERPOLATION_FACTOR(IE,2) 
            IF (NOM2<0) WM = WGT*OM_WW(IIO2,JJO2,KKO2) + (1._EB-WGT)*OM_WW(IIO2,JJO2,KKO2-1) 
            IF (NOM2>0) WP = WGT*OM_WW(IIO2,JJO2,KKO2) + (1._EB-WGT)*OM_WW(IIO2,JJO2,KKO2-1) 
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

         IF (NOM1==NM) THEN
            IWM = WALL_INDEX(ICMM, 2) 
            IWP = WALL_INDEX(ICPM, 2) 
            IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
               IF (BC<1.5_EB) UP = BC*UM
               IF (BC>1.5_EB) UP = FVT*SF%VEL_T(1)*PROF
            ENDIF
         ENDIF
         IF (NOM1==-NM) THEN
            IWM = WALL_INDEX(ICMP,-2) 
            IWP = WALL_INDEX(ICPP,-2) 
            IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
               IF (BC<1.5_EB) UM = BC*UP
               IF (BC>1.5_EB) UM = FVT*SF%VEL_T(1)*PROF
            ENDIF
         ENDIF

         IF (NOM2==NM) THEN
            IWM = WALL_INDEX(ICMM, 1) 
            IWP = WALL_INDEX(ICMP, 1) 
            IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
               IF (BC<1.5_EB) VP = BC*VM
               IF (BC>1.5_EB) VP = FVT*SF%VEL_T(1)*PROF
            ENDIF
         ENDIF
         IF (NOM2==-NM) THEN
            IWM = WALL_INDEX(ICPM,-1) 
            IWP = WALL_INDEX(ICPP,-1) 
            IF (BOUNDARY_TYPE(IWM)/=NULL_BOUNDARY .AND. BOUNDARY_TYPE(IWP)/=NULL_BOUNDARY) THEN
               IF (BC<1.5_EB) VM = BC*VP
               IF (BC>1.5_EB) VM = FVT*SF%VEL_T(1)*PROF
            ENDIF
         ENDIF

         IF (ABS(NOM1)/=NM .AND. NOM1/=0) THEN
            OM => OMESH(ABS(NOM1)) 
            IF (PREDICTOR) THEN
               OM_UU => OM%US
            ELSE
               OM_UU => OM%U
            ENDIF
            WGT = EDGE_INTERPOLATION_FACTOR(IE,1) 
            IF (NOM1<0) UM = WGT*OM_UU(IIO1,JJO1,KKO1) + (1._EB-WGT)*OM_UU(IIO1-1,JJO1,KKO1)
            IF (NOM1>0) UP = WGT*OM_UU(IIO1,JJO1,KKO1) + (1._EB-WGT)*OM_UU(IIO1-1,JJO1,KKO1)
         ENDIF
         IF (ABS(NOM2)/=NM .AND. NOM2/=0) THEN
            OM => OMESH(ABS(NOM2)) 
            IF (PREDICTOR) THEN
               OM_VV => OM%VS
            ELSE
               OM_VV => OM%V
            ENDIF
            WGT = EDGE_INTERPOLATION_FACTOR(IE,2) 
            IF (NOM2<0) VM = WGT*OM_VV(IIO2,JJO2,KKO2) + (1._EB-WGT)*OM_VV(IIO2,JJO2-1,KKO2) 
            IF (NOM2>0) VP = WGT*OM_VV(IIO2,JJO2,KKO2) + (1._EB-WGT)*OM_VV(IIO2,JJO2-1,KKO2) 
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

   IF (IOR==0) THEN
      OME_E(IE) = 0._EB
      TAU_E(IE) = 0._EB
   ENDIF

ENDDO EDGE_LOOP

! Store cell node averages of the velocity components in UVW_GHOST for use in Smokeview

IF (CORRECTOR) THEN
   DO K=0,KBAR
      DO J=0,JBAR
         I_LOOP: DO I=0,IBAR
            IC = CELL_INDEX(I,J,K) 
            IF (IC==0) CYCLE I_LOOP
            IF (U_Y(I,J,K)  >-1.E5_EB) UVW_GHOST(IC,1) = U_Y(I,J,K) 
            IF (U_Y(I,J,K+1)>-1.E5_EB) UVW_GHOST(IC,1) = U_Y(I,J,K+1) 
            IF (U_Z(I,J,K)  >-1.E5_EB) UVW_GHOST(IC,1) = U_Z(I,J,K) 
            IF (U_Z(I,J+1,K)>-1.E5_EB) UVW_GHOST(IC,1) = U_Z(I,J+1,K) 
            IF (V_X(I,J,K)  >-1.E5_EB) UVW_GHOST(IC,2) = V_X(I,J,K) 
            IF (V_X(I,J,K+1)>-1.E5_EB) UVW_GHOST(IC,2) = V_X(I,J,K+1) 
            IF (V_Z(I,J,K)  >-1.E5_EB) UVW_GHOST(IC,2) = V_Z(I,J,K) 
            IF (V_Z(I+1,J,K)>-1.E5_EB) UVW_GHOST(IC,2) = V_Z(I+1,J,K) 
            IF (W_X(I,J,K)  >-1.E5_EB) UVW_GHOST(IC,3) = W_X(I,J,K) 
            IF (W_X(I,J+1,K)>-1.E5_EB) UVW_GHOST(IC,3) = W_X(I,J+1,K) 
            IF (W_Y(I,J,K)  >-1.E5_EB) UVW_GHOST(IC,3) = W_Y(I,J,K) 
            IF (W_Y(I+1,J,K)>-1.E5_EB) UVW_GHOST(IC,3) = W_Y(I+1,J,K) 
       !!!  UVW_GHOST(IC,1) = 0.25_EB*(U_Y(I,J,K)+U_Y(I,J,K+1)+U_Z(I,J,K)+U_Z(I,J+1,K)) 
       !!!  UVW_GHOST(IC,2) = 0.25_EB*(V_X(I,J,K)+V_X(I,J,K+1)+V_Z(I,J,K)+V_Z(I+1,J,K)) 
       !!!  UVW_GHOST(IC,3) = 0.25_EB*(W_X(I,J,K)+W_X(I,J+1,K)+W_Y(I,J,K)+W_Y(I+1,J,K)) 
         ENDDO I_LOOP
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

IF (SOLID_PHASE_ONLY .OR. NMESHES==1) RETURN

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
         WW(II,JJ,KBAR) = WW_AVG
   END SELECT

ENDDO EXTERNAL_WALL_LOOP

TUSED(4,NM)=TUSED(4,NM)+SECOND()-TNOW
END SUBROUTINE MATCH_VELOCITY



SUBROUTINE CHECK_STABILITY(NM)
 
! Checks the Courant and Von Neumann stability criteria, and if necessary, reduces the time step accordingly
 
REAL(EB) :: UODX,VODY,WODZ,UVW,UVWMAX,R_DX2,MU_MAX,MUTRM
INTEGER  :: NM,I,J,K
 
CHANGE_TIME_STEP(NM) = .FALSE.
UVWMAX = 0._EB
VN     = 0._EB
MUTRM  = 1.E-9_EB
R_DX2  = 1.E-9_EB
 
! Determine max CFL number from all grid cells
 
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
 
CFL = DT*UVWMAX
 
! Determine max Von Neumann Number for fine grid calcs
 
PARABOLIC_IF: IF (DNS .OR. DXMIN<0.005_EB) THEN
 
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
 
IF (CFL<CFL_MAX .AND. VN<VN_MAX) THEN
   DTNEXT = DT
   IF (CFL<=CFL_MIN .AND. VN<VN_MIN) DTNEXT = MIN(1.1_EB*DT,DTINT)
ELSE
   IF (UVWMAX==0._EB) UVWMAX = 1._EB
   DT = 0.9_EB*MIN( CFL_MAX/UVWMAX , VN_MAX/(2._EB*R_DX2*MUTRM) )
   CHANGE_TIME_STEP(NM) = .TRUE.
ENDIF
 
END SUBROUTINE CHECK_STABILITY
 
 
SUBROUTINE BAROCLINIC_CORRECTION
 
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,HQS,RTRM,RHOP
REAL(EB) :: RHO_AVG_OLD,RMIN,RMAX,RRAT,U2,V2,W2
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
RMIN =  1000._EB
RMAX = -1000._EB
DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (.NOT.SOLID(CELL_INDEX(I,J,K))) THEN
               RMIN = MIN(RHOP(I,J,K),RMIN)
               RMAX = MAX(RHOP(I,J,K),RMAX)
            ENDIF
         ENDDO
      ENDDO
ENDDO
 
RHO_AVG = 2._EB*RMIN*RMAX/(RMIN+RMAX)
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
