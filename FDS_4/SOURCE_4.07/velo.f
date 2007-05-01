      MODULE VELO
C
      USE PREC
      USE CONS
      USE PACKER
C
      IMPLICIT NONE
      PRIVATE
      REAL(EB) MUX,MUY,MUZ,UP,UM,VP,VM,WP,WM,VTRM,RHS,
     .        DTXYDY,DTXZDZ,DTYZDZ,DTXYDX,DTXZDX,DTYZDY,
     .        DUDX,DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY,
     .        VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX,PMDT,MPDT,
     .        S2,C,CDXDYDZTT,AH,RRHO,GX,GY,GZ,
     .        TXXP,TXXM,TYYP,TYYM,TZZP,TZZM,DTXXDX,DTYYDY,DTZZDZ,
     .        EPSUP,EPSUM,EPSVP,EPSVM,EPSWP,EPSWM,MU_SUM,
     .        U2,V2,W2,RMIN,RMAX,RHO_AVG_OLD,RRAT
      INTEGER II,JJ,KK,I,J,K,IW,IIG,JJG,KKG,ITMP,N,IE,IOR,IYY
      REAL(EB), POINTER, DIMENSION(:,:,:) :: TXY,TXZ,TYZ,OMX,OMY,OMZ,
     .                                       UU,VV,WW,RHOP,DP,HQS,RTRM
      REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP
      PUBLIC VELOCITY_FLUX,VELOCITY_FLUX_ISOTHERMAL,VELOCITY_PREDICTOR,
     .       VELOCITY_CORRECTOR,VELOCITY_FLUX_CYLINDRICAL,NO_FLUX
      TYPE (OBSTRUCTION_TYPE), POINTER :: OB
      TYPE (VENTS_TYPE), POINTER :: VT
C
C
      CONTAINS
C
C
      SUBROUTINE VELOCITY_FLUX(T,NM)
C
C Compute convective and diffusive terms
C
      REAL(EB) T
      INTEGER, INTENT(IN) :: NM
C
      CALL UNPACK_VAR(NM)
C
      TNOW=SECOND()
C
      IF (PREDICTOR) THEN
         UU => U
         VV => V
         WW => W
         DP => D  
       RHOP => RHO
         IF (NSPEC.GT.0) YYP => YY
         ELSE
         UU => US
         VV => VS
         WW => WS
         DP => DS
       RHOP => RHOS
         IF (NSPEC.GT.0) YYP => YYS
         ENDIF
C
      TXY => WORK1
      TXZ => WORK2
      TYZ => WORK3
      OMX => WORK4
      OMY => WORK5
      OMZ => WORK6
C
      CALC_MU: IF (PREDICTOR) THEN    
C
      LES_VS_DNS: IF (LES) THEN   ! Smagorinsky model (LES)
C
      C = CSMAG**2
      IF (EVACUATION_ONLY(NM)) C = (0.9)**2
C
      KLOOP: DO K=1,KBAR
      JLOOP: DO J=1,JBAR
      ILOOP: DO I=1,IBAR
      IF (SOLID(ICA(I,J,K))) CYCLE ILOOP
      IF (.NOT.TWO_D) THEN
         CDXDYDZTT = C*(DX(I)*DY(J)*DZ(K))**TWTH
         ELSE
         CDXDYDZTT = C*DX(I)*DZ(K)
         ENDIF
      DUDX = RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
      DVDY = RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
      DWDZ = RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
      DUDY = 0.25*RDY(J)
     .       *(UU(I,J+1,K)-UU(I,J-1,K)+UU(I-1,J+1,K)-UU(I-1,J-1,K))
      DVDX = 0.25*RDX(I) 
     .       *(VV(I+1,J,K)-VV(I-1,J,K)+VV(I+1,J-1,K)-VV(I-1,J-1,K))
      DWDX = 0.25*RDX(I) 
     .       *(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
      DWDY = 0.25*RDY(J) 
     .       *(WW(I,J+1,K)-WW(I,J-1,K)+WW(I,J+1,K-1)-WW(I,J-1,K-1))
      DUDZ = 0.25*RDZ(K) 
     .       *(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1)) 
      DVDZ = 0.25*RDZ(K) 
     .       *(VV(I,J,K+1)-VV(I,J,K-1)+VV(I,J-1,K+1)-VV(I,J-1,K-1)) 
      S2   = 2.*(DUDX*DUDX     +  DVDY*DVDY     +  DWDZ*DWDZ ) 
     .       +  (DUDY+DVDX)**2 + (DUDZ+DWDX)**2 + (DVDZ+DWDY)**2
     .       -  TWTH*DP(I,J,K)**2
      S2   = MAX(0._EB,S2)
      ITMP = 0.1*TMP(I,J,K)
      MU(I,J,K) = MAX(MU_SPEC(0,ITMP),RHOP(I,J,K)*CDXDYDZTT*SQRT(S2))
      ENDDO ILOOP
      ENDDO JLOOP
      ENDDO KLOOP
C
      ELSE   ! DNS viscosity
C
      MIXTURE_FRACTION_DNS: IF (.NOT.MIXTURE_FRACTION) THEN
C
      DO K=1,KBAR
      DO J=1,JBAR
      IVLOOP: DO I=1,IBAR
      IF (SOLID(ICA(I,J,K))) CYCLE IVLOOP
      ITMP = 0.1*TMP(I,J,K)
      MU_SUM = MU_SPEC(0,ITMP)
      DO N=1,NSPEC
      MU_SUM = MU_SUM + YYP(I,J,K,N)*(MU_SPEC(N,ITMP)-MU_SPEC(0,ITMP))
      ENDDO
      MU(I,J,K) = MU_SUM
      ENDDO IVLOOP
      ENDDO
      ENDDO
C
      ELSE MIXTURE_FRACTION_DNS
C
      DO K=1,KBAR
      DO J=1,JBAR
      IVLOOP2: DO I=1,IBAR
      IF (SOLID(ICA(I,J,K))) CYCLE IVLOOP2
      ITMP = 0.1*TMP(I,J,K)
      IYY  = MAX(0,NINT(YYP(I,J,K,IFUEL)*100.))
      IYY  = MIN(100,IYY)
      MU(I,J,K) = MU_SPEC(IYY,ITMP)
      ENDDO IVLOOP2
      ENDDO
      ENDDO
C
      ENDIF MIXTURE_FRACTION_DNS
C
      ENDIF LES_VS_DNS
C
C Mirror viscosity into solids
C
      BLOOP: DO IW=1,NWC
      IF (IV(IW).EQ.0) CYCLE BLOOP
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      MU(II,JJ,KK) = MU(IIG,JJG,KKG)
      ENDDO BLOOP
C
      ENDIF CALC_MU
C
C Compute vorticity and stress tensor components
C
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
      MUX = 0.25*(MU(I,J+1,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I,J+1,K+1))
      MUY = 0.25*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I+1,J,K+1))
      MUZ = 0.25*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J+1,K)+MU(I+1,J+1,K))
      TXY(I,J,K) = MUZ*(DVDX + DUDY)
      TXZ(I,J,K) = MUY*(DUDZ + DWDX)
      TYZ(I,J,K) = MUX*(DVDZ + DWDY)
      ENDDO
      ENDDO
      ENDDO
C
C Correct vorticity and stress tensor components at solid edges
C
      EDGE_LOOP: DO IE=1,NEDGES
      II  = IJKE(1,IE)
      JJ  = IJKE(2,IE)
      KK  = IJKE(3,IE)
      IOR = IJKE(4,IE)
      SELECT CASE(IOR)
      CASE(1) ; OMX(II,JJ,KK) = OME_E(IE) ; TYZ(II,JJ,KK) = TAU_E(IE)
      CASE(2) ; OMY(II,JJ,KK) = OME_E(IE) ; TXZ(II,JJ,KK) = TAU_E(IE)
      CASE(3) ; OMZ(II,JJ,KK) = OME_E(IE) ; TXY(II,JJ,KK) = TAU_E(IE)
      END SELECT
      ENDDO EDGE_LOOP
C
C Compute gravity components
C
      GX  = FGX(T)*GVEC(1)
      GY  = FGY(T)*GVEC(2)
      GZ  = FGZ(T)*GVEC(3)
C
C Upwind/Downwind bias factors
C
      IF (PREDICTOR) THEN
      PMDT =  0.5*DT
      MPDT = -0.5*DT
      ELSE
      PMDT = -0.5*DT
      MPDT =  0.5*DT
      ENDIF
C
C Compute x-direction flux terms
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=0,IBAR
      WP    = WW(I,J,K)   + WW(I+1,J,K)
      WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
      VP    = VV(I,J,K)   + VV(I+1,J,K)
      VM    = VV(I,J-1,K) + VV(I+1,J-1,K)
      EPSWP = 1. + WP*MPDT*RDZN(K)
      EPSWM = 1. + WM*PMDT*RDZN(K-1)
      EPSVP = 1. + VP*MPDT*RDYN(J)
      EPSVM = 1. + VM*PMDT*RDYN(J-1)
      WOMY  = EPSWP*WP*OMY(I,J,K) + EPSWM*WM*OMY(I,J,K-1)
      VOMZ  = EPSVP*VP*OMZ(I,J,K) + EPSVM*VM*OMZ(I,J-1,K)
      RRHO  = 2./(RHOP(I,J,K)+RHOP(I+1,J,K))
      AH    = RHOA*RRHO - 1.   
      DVDY  = (VV(I+1,J,K)-VV(I+1,J-1,K))*RDY(J)
      DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*RDZ(K)
      TXXP  = MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2.*(DVDY+DWDZ) )
      DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
      DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
      TXXM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2.*(DVDY+DWDZ) )
      DTXXDX= RDXN(I)*(TXXP      -TXXM)
      DTXYDY= RDY(J) *(TXY(I,J,K)-TXY(I,J-1,K))
      DTXZDZ= RDZ(K) *(TXZ(I,J,K)-TXZ(I,J,K-1))
      VTRM  = RRHO*(DTXXDX + DTXYDY + DTXZDZ)
      FVX(I,J,K) = 0.25*(WOMY - VOMZ) + GX*AH - VTRM 
      ENDDO 
      ENDDO   
      ENDDO   
C
C Compute y-direction flux terms
C
      DO K=1,KBAR
      DO J=0,JBAR
      DO I=1,IBAR
      UP    = UU(I,J,K)   + UU(I,J+1,K)
      UM    = UU(I-1,J,K) + UU(I-1,J+1,K)
      WP    = WW(I,J,K)   + WW(I,J+1,K)
      WM    = WW(I,J,K-1) + WW(I,J+1,K-1)
      EPSUP = 1. + UP*MPDT*RDXN(I)
      EPSUM = 1. + UM*PMDT*RDXN(I-1)
      EPSWP = 1. + WP*MPDT*RDZN(K)
      EPSWM = 1. + WM*PMDT*RDZN(K-1)
      WOMX  = EPSWP*WP*OMX(I,J,K) + EPSWM*WM*OMX(I,J,K-1)
      UOMZ  = EPSUP*UP*OMZ(I,J,K) + EPSUM*UM*OMZ(I-1,J,K)
      RRHO  = 2./(RHOP(I,J,K)+RHOP(I,J+1,K))
      AH    = RHOA*RRHO - 1.
      DUDX  = (UU(I,J+1,K)-UU(I-1,J+1,K))*RDX(I)
      DWDZ  = (WW(I,J+1,K)-WW(I,J+1,K-1))*RDZ(K)
      TYYP  = MU(I,J+1,K)*( FOTH*DP(I,J+1,K) - 2.*(DUDX+DWDZ) )
      DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
      DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
      TYYM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2.*(DUDX+DWDZ) )
      DTXYDX= RDX(I) *(TXY(I,J,K)-TXY(I-1,J,K))
      DTYYDY= RDYN(J)*(TYYP      -TYYM)
      DTYZDZ= RDZ(K) *(TYZ(I,J,K)-TYZ(I,J,K-1))
      VTRM  = RRHO*(DTXYDX + DTYYDY + DTYZDZ)
      FVY(I,J,K) = 0.25*(UOMZ - WOMX) + GY*AH - VTRM 
      ENDDO
      ENDDO   
      ENDDO   
C
C Compute z-direction flux terms
C
      DO K=0,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      UP    = UU(I,J,K)   + UU(I,J,K+1)
      UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
      VP    = VV(I,J,K)   + VV(I,J,K+1)
      VM    = VV(I,J-1,K) + VV(I,J-1,K+1)
      EPSUP = 1. + UP*MPDT*RDXN(I)
      EPSUM = 1. + UM*PMDT*RDXN(I-1)
      EPSVP = 1. + VP*MPDT*RDYN(J)
      EPSVM = 1. + VM*PMDT*RDYN(J-1)
      UOMY  = EPSUP*UP*OMY(I,J,K) + EPSUM*UM*OMY(I-1,J,K)
      VOMX  = EPSVP*VP*OMX(I,J,K) + EPSVM*VM*OMX(I,J-1,K)
      RRHO  = 2./(RHOP(I,J,K)+RHOP(I,J,K+1))
      AH    = RHOA*RRHO - 1.
      DUDX  = (UU(I,J,K+1)-UU(I-1,J,K+1))*RDX(I)
      DVDY  = (VV(I,J,K+1)-VV(I,J-1,K+1))*RDY(J)
      TZZP  = MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2.*(DUDX+DVDY) )
      DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
      DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
      TZZM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2.*(DUDX+DVDY) )
      DTXZDX= RDX(I) *(TXZ(I,J,K)-TXZ(I-1,J,K))
      DTYZDY= RDY(J) *(TYZ(I,J,K)-TYZ(I,J-1,K))
      DTZZDZ= RDZN(K)*(TZZP      -TZZM)
      VTRM  = RRHO*(DTXZDX + DTYZDY + DTZZDZ)
      FVZ(I,J,K) = 0.25*(VOMX - UOMY) + GZ*AH - VTRM 
      ENDDO
      ENDDO   
      ENDDO   
C
C Baroclinic torque correction
C
      IF (BAROCLINIC) CALL BAROCLINIC_CORRECTION
C
C Adjust FVX, FVY and FVZ at solid, internal obstructions for no flux
C
      CALL NO_FLUX(T)
C
      TUSED(4,NM) = TUSED(4,NM) + SECOND() - TNOW
      END SUBROUTINE VELOCITY_FLUX
C
C
      SUBROUTINE VELOCITY_FLUX_ISOTHERMAL(T,NM)
C
C Compute the velocity flux at cell edges (ISOTHERMAL DNS ONLY)
C
      REAL(EB) UP,UM,VP,VM,WP,WM,VTRM,T,
     .        VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX,
     .        DVDZ,DVDX,DWDY,DWDX,DUDZ,DUDY,FROUDEGZ
      INTEGER I,J,K
      REAL(EB), POINTER, DIMENSION(:,:,:) :: OMX,OMY,OMZ,UU,VV,WW
      REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP
      INTEGER, INTENT(IN) :: NM
C
      CALL UNPACK_VAR(NM)
C
      TNOW=SECOND()
C
      IF (PREDICTOR) THEN
         UU => U
         VV => V
         WW => W
         ELSE
         UU => US
         VV => VS
         WW => WS
         ENDIF
C
      OMX => WORK4
      OMY => WORK5
      OMZ => WORK6
C
C Compute vorticity and stress tensor components
C
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
C
C Correct vorticity and stress tensor components at solid edges
C
      EDGE_LOOP: DO IE=1,NEDGES
      II  = IJKE(1,IE)
      JJ  = IJKE(2,IE)
      KK  = IJKE(3,IE)
      IOR = IJKE(4,IE)
      SELECT CASE(IOR)
      CASE(1) ; OMX(II,JJ,KK) = OME_E(IE)
      CASE(2) ; OMY(II,JJ,KK) = OME_E(IE)
      CASE(3) ; OMZ(II,JJ,KK) = OME_E(IE)
      END SELECT
      ENDDO EDGE_LOOP
C
C Upwind/Downwind bias factors
C
      IF (PREDICTOR) THEN
      PMDT =  0.5*DT
      MPDT = -0.5*DT
      ELSE
      PMDT = -0.5*DT
      MPDT =  0.5*DT
      ENDIF
C
C Compute x-direction flux terms
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=0,IBAR
      WP    = WW(I,J,K)   + WW(I+1,J,K)
      WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
      VP    = VV(I,J,K)   + VV(I+1,J,K)
      VM    = VV(I,J-1,K) + VV(I+1,J-1,K)
      EPSWP = 1. + WP*MPDT*RDZN(K)
      EPSWM = 1. + WM*PMDT*RDZN(K-1)
      EPSVP = 1. + VP*MPDT*RDYN(J)
      EPSVM = 1. + VM*PMDT*RDYN(J-1)
      WOMY  = EPSWP*WP*OMY(I,J,K) + EPSWM*WM*OMY(I,J,K-1)
      VOMZ  = EPSVP*VP*OMZ(I,J,K) + EPSVM*VM*OMZ(I,J-1,K)
      VTRM  = RREDZ(K)*(OMY(I,J,K)-OMY(I,J,K-1)) -
     .        RREDY(J)*(OMZ(I,J,K)-OMZ(I,J-1,K))
      FVX(I,J,K) = 0.25*(WOMY - VOMZ) - VTRM
      ENDDO
      ENDDO
      ENDDO
C
C Compute y-direction flux terms
C
      DO K=1,KBAR
      DO J=0,JBAR
      DO I=1,IBAR
      UP    = UU(I,J,K)   + UU(I,J+1,K)
      UM    = UU(I-1,J,K) + UU(I-1,J+1,K)
      WP    = WW(I,J,K)   + WW(I,J+1,K)
      WM    = WW(I,J,K-1) + WW(I,J+1,K-1)
      EPSUP = 1. + UP*MPDT*RDXN(I)
      EPSUM = 1. + UM*PMDT*RDXN(I-1)
      EPSWP = 1. + WP*MPDT*RDZN(K)
      EPSWM = 1. + WM*PMDT*RDZN(K-1)
      WOMX  = EPSWP*WP*OMX(I,J,K) + EPSWM*WM*OMX(I,J,K-1)
      UOMZ  = EPSUP*UP*OMZ(I,J,K) + EPSUM*UM*OMZ(I-1,J,K)
      VTRM  = RREDX(I)*(OMZ(I,J,K)-OMZ(I-1,J,K)) -
     .        RREDZ(K)*(OMX(I,J,K)-OMX(I,J,K-1))
      FVY(I,J,K) = 0.25*(UOMZ - WOMX) - VTRM
      ENDDO
      ENDDO
      ENDDO
C
C Compute z-direction flux terms
C
      DO K=0,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      UP    = UU(I,J,K)   + UU(I,J,K+1)
      UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
      VP    = VV(I,J,K)   + VV(I,J,K+1)
      VM    = VV(I,J-1,K) + VV(I,J-1,K+1)
      EPSUP = 1. + UP*MPDT*RDXN(I)
      EPSUM = 1. + UM*PMDT*RDXN(I-1)
      EPSVP = 1. + VP*MPDT*RDYN(J)
      EPSVM = 1. + VM*PMDT*RDYN(J-1)
      UOMY  = EPSUP*UP*OMY(I,J,K) + EPSUM*UM*OMY(I-1,J,K)
      VOMX  = EPSVP*VP*OMX(I,J,K) + EPSVM*VM*OMX(I,J-1,K)
      VTRM  = RREDY(J)*(OMX(I,J,K)-OMX(I,J-1,K)) -
     .        RREDX(I)*(OMY(I,J,K)-OMY(I-1,J,K))
      FVZ(I,J,K) = 0.25*(VOMX - UOMY) - VTRM
      ENDDO
      ENDDO
      ENDDO
C
C Special gravity term for INCOMPRESSIBLE case (FROUDE scaling)
C
      IF (NSPEC.GT.0) THEN
      IF (PREDICTOR) YYP => YY
      IF (CORRECTOR) YYP => YYS
      FROUDEGZ = 0.5*GVEC(3)*(MWN(1)-MWN(0))/MWN(0)
      DO K=0,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      FVZ(I,J,K) = FVZ(I,J,K) - FROUDEGZ*(YYP(I,J,K,1)+YYP(I,J,K+1,1))
      ENDDO
      ENDDO
      ENDDO
      ENDIF
C
C Adjust FVX, FVY and FVZ at solid, internal obstructions for no flux
C
      CALL NO_FLUX(T)
C
      TUSED(4,NM) = TUSED(4,NM) + SECOND() - TNOW
      END SUBROUTINE VELOCITY_FLUX_ISOTHERMAL
C
C
      SUBROUTINE VELOCITY_FLUX_CYLINDRICAL(T,NM)
C
C Compute convective and diffusive terms for 2D axisymmetric
C
      REAL(EB) :: T,DMUDX
      INTEGER :: I0
      INTEGER, INTENT(IN) :: NM
C
      CALL UNPACK_VAR(NM)
C
      TNOW=SECOND()
C
      IF (PREDICTOR) THEN
         UU => U
         WW => W
         DP => D  
       RHOP => RHO
         IF (NSPEC.GT.0) YYP => YY
         ELSE
         UU => US
         WW => WS
         DP => DS
       RHOP => RHOS
         IF (NSPEC.GT.0) YYP => YYS
         ENDIF
C
      TXZ => WORK2
      OMY => WORK5
C
      CALC_MU: IF (PREDICTOR) THEN    
C
      LES_VS_DNS: IF (LES) THEN   ! Smagorinsky model (LES)
C
      C = CSMAG**2
C
      J = 1
      KLOOP: DO K=1,KBAR
      ILOOP: DO I=1,IBAR
      IF (SOLID(ICA(I,J,K))) CYCLE ILOOP
      CDXDYDZTT = C*DX(I)*DZ(K)
      DUDX = RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
      DWDZ = RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
      DWDX = 0.25*RDX(I)
     .       *(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
      DUDZ = 0.25*RDZ(K)
     .       *(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1))
      S2   = 2.*(DUDX*DUDX     +  DWDZ*DWDZ )
     .       +  (DUDZ+DWDX)**2 -  TWTH*DP(I,J,K)**2
      S2   = MAX(0._EB,S2)
      ITMP = 0.1*TMP(I,J,K)
      MU(I,J,K) = MAX(MU_SPEC(0,ITMP),RHOP(I,J,K)*CDXDYDZTT*SQRT(S2))
      ENDDO ILOOP
      ENDDO KLOOP
C
      ELSE LES_VS_DNS           ! DNS viscosity
C
      MIXTURE_FRACTION_DNS: IF (.NOT.MIXTURE_FRACTION) THEN
C
      J = 1
      DO K=1,KBAR
      IVLOOP: DO I=1,IBAR
      IF (SOLID(ICA(I,J,K))) CYCLE IVLOOP
      ITMP = 0.1*TMP(I,J,K)
      MU_SUM = MU_SPEC(0,ITMP)
      DO N=1,NSPEC
      MU_SUM = MU_SUM + YYP(I,J,K,N)*(MU_SPEC(N,ITMP)-MU_SPEC(0,ITMP))
      ENDDO
      MU(I,J,K) = MU_SUM
      ENDDO IVLOOP
      ENDDO
C
      ELSE MIXTURE_FRACTION_DNS
C
      J = 1
      DO K=1,KBAR
      IVLOOP2: DO I=1,IBAR
      IF (SOLID(ICA(I,J,K))) CYCLE IVLOOP2
      ITMP = 0.1*TMP(I,J,K)
      IYY  = MAX(0,NINT(YYP(I,J,K,IFUEL)*100.))
      IYY  = MIN(100,IYY)
      MU(I,J,K) = MU_SPEC(IYY,ITMP)
      ENDDO IVLOOP2
      ENDDO
C
      ENDIF MIXTURE_FRACTION_DNS
C
      ENDIF LES_VS_DNS
C
C Mirror viscosity into solids
C
      DO IW=1,NWC
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      MU(II,JJ,KK) = MU(IIG,JJG,KKG)
      ENDDO
C
      ENDIF CALC_MU
C
C Compute vorticity and stress tensor components
C
      DO K=0,KBAR
      DO J=0,JBAR
      DO I=0,IBAR
      DUDZ = RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
      DWDX = RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
      OMY(I,J,K) = DUDZ - DWDX
      MUY = 0.25*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I+1,J,K+1))
      TXZ(I,J,K) = MUY*(DUDZ + DWDX)
      ENDDO
      ENDDO
      ENDDO
C
C Correct vorticity and stress tensor components at solid edges
C
      EDGE_LOOP: DO IE=1,NEDGES
      II  = IJKE(1,IE)
      JJ  = IJKE(2,IE)
      KK  = IJKE(3,IE)
      IOR = IJKE(4,IE)
      SELECT CASE(IOR)
      CASE(2) ; OMY(II,JJ,KK) = OME_E(IE) ; TXZ(II,JJ,KK) = TAU_E(IE)
      END SELECT
      ENDDO EDGE_LOOP
C
C Compute gravity components
C
      GX  = 0.
      GY  = 0.
      GZ  = FGZ(T)*GVEC(3)
C
C Upwind/Downwind bias factors
C
      IF (PREDICTOR) THEN
      PMDT =  0.5*DT
      MPDT = -0.5*DT
      ELSE
      PMDT = -0.5*DT
      MPDT =  0.5*DT
      ENDIF
C
C Compute r-direction flux terms
C
      IF (XS.EQ.0) THEN ; I0 = 1 ; ELSE ; I0 = 0 ; ENDIF
C
      J = 1
      DO K= 1,KBAR
      DO I=I0,IBAR
      WP    = WW(I,J,K)   + WW(I+1,J,K)
      WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
      EPSWP = 1. + WP*MPDT*RDZN(K)
      EPSWM = 1. + WM*PMDT*RDZN(K-1)
      WOMY  = EPSWP*WP*OMY(I,J,K) + EPSWM*WM*OMY(I,J,K-1)
      RRHO  = 2./(RHOP(I,J,K)+RHOP(I+1,J,K))
      AH    = RHOA*RRHO - 1.   
      DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*RDZ(K)
      TXXP  = MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2.*DWDZ )
      DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
      TXXM  = MU(I,J,K)  *( FOTH*DP(I,J,K) -2.*DWDZ )
      DTXXDX= RDXN(I)*(TXXP      -TXXM)
      DTXZDZ= RDZ(K) *(TXZ(I,J,K)-TXZ(I,J,K-1))
      DMUDX = (MU(I+1,J,K)-MU(I,J,K))*RDXN(I)
      VTRM  = RRHO*( DTXXDX + DTXZDZ - 2.*UU(I,J,K)*DMUDX/R(I) ) 
      FVX(I,J,K) = 0.25*WOMY + GX*AH - VTRM 
      ENDDO
      ENDDO   
C
C Compute z-direction flux terms
C
      J = 1
      DO K=0,KBAR
      DO I=1,IBAR
      UP    = UU(I,J,K)   + UU(I,J,K+1)
      UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
      EPSUP = 1. + UP*MPDT*RDXN(I)
      EPSUM = 1. + UM*PMDT*RDXN(I-1)
      UOMY  = EPSUP*UP*OMY(I,J,K) + EPSUM*UM*OMY(I-1,J,K)
      RRHO  = 2./(RHOP(I,J,K)+RHOP(I,J,K+1))
      AH    = RHOA*RRHO - 1.
      DUDX  = (R(I)*UU(I,J,K+1)-R(I-1)*UU(I-1,J,K+1))*RDX(I)*RRN(I)
      TZZP  = MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2.*DUDX )
      DUDX  = (R(I)*UU(I,J,K)-R(I-1)*UU(I-1,J,K))*RDX(I)*RRN(I)
      TZZM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2.*DUDX )
      DTXZDX= RDX(I) *(R(I)*TXZ(I,J,K)-R(I-1)*TXZ(I-1,J,K))*RRN(I)
      DTZZDZ= RDZN(K)*(TZZP      -TZZM)
      VTRM  = RRHO*(DTXZDX + DTZZDZ)
      FVZ(I,J,K) = -0.25*UOMY + GZ*AH - VTRM 
      ENDDO
      ENDDO   
C
C Baroclinic torque correction terms
C
      IF (BAROCLINIC) CALL BAROCLINIC_CORRECTION
C
C Adjust FVX and FVZ at solid, internal obstructions for no flux
C
      CALL NO_FLUX(T)
C
      TUSED(4,NM) = TUSED(4,NM) + SECOND() - TNOW
      END SUBROUTINE VELOCITY_FLUX_CYLINDRICAL
C
C
      SUBROUTINE NO_FLUX(T)
C
C Set FVX,FVY,FVZ inside internal blockages to maintain no flux
C
      REAL(EB) RFODT,VEL,T,FVT,VEL0
      INTEGER IC2,IC1,IBC,ITAU
C
      IF (PREDICTOR) THEN
         UU => U
         VV => V
         WW => W
         ELSE
         UU => US
         VV => VS
         WW => WS
         ENDIF
C
      RFODT = RF/DT
C
C Blow air through VENTs with no solid backing
C
      VENT_LOOP: DO N=1,NV
C
      VT => VENTS(N)
C
      IF (VT%INDEX.NE.4)         CYCLE VENT_LOOP
      IF (VT%T_OPEN .LT.100000.) CYCLE VENT_LOOP
      IF (VT%T_CLOSE.EQ.999999.) CYCLE VENT_LOOP
      IF (VT%HEAT_INDEX_ACTIVATE.GT.0) CYCLE VENT_LOOP
C
      SELECT CASE(ABS(VT%IOR))
         CASE(1) ; VEL0 = U0
         CASE(2) ; VEL0 = V0
         CASE(3) ; VEL0 = W0
         END SELECT
      IBC = VT%IBC
      IF (TAUV(IBC).LT.1000000.) ITAU = IBC
      IF (TAUV(IBC).GT.1000000.) ITAU = -(TAUV(IBC)-1000000.)
      IF (TAUV(IBC).EQ.      0.) ITAU = 999
      FVT = FV(T-TIGNS(IBC),ITAU)
      VEL = VEL0
      IF (VELS(IBC) .NE.-999.) VEL = VEL0+FVT*(VELS(IBC)-VEL0)
      IF (VFLUX(IBC).NE.-999.) VEL = VEL0+FVT*(VFLUX(IBC)/VT%AREA-VEL0)
C
      SELECT CASE(ABS(VT%IOR))
C
      CASE(1)
      DO K=VT%K1+1,VT%K2
      DO J=VT%J1+1,VT%J2
      I = VT%I1
      FVX(I,J,K) = -RDXN(I)*(H(I+1,J,K)-H(I,J,K))+RFODT*(UU(I,J,K)-VEL)
      ENDDO
      ENDDO
C
      CASE(2)
      DO K=VT%K1+1,VT%K2
      DO I=VT%I1+1,VT%I2
      J = VT%J1
      FVY(I,J,K) = -RDYN(J)*(H(I,J+1,K)-H(I,J,K))+RFODT*(VV(I,J,K)-VEL)
      ENDDO
      ENDDO
C
      CASE(3)
      DO J=VT%J1+1,VT%J2
      DO I=VT%I1+1,VT%I2
      K = VT%K1
      FVZ(I,J,K) = -RDZN(K)*(H(I,J,K+1)-H(I,J,K))+RFODT*(WW(I,J,K)-VEL)
      ENDDO
      ENDDO
C
      END SELECT
C
      ENDDO VENT_LOOP
C
C Stop (or control) air flow at solid obstructions
C
      OBST_LOOP: DO N=1,NB
      OB=>OBSTRUCTION(N)
C
             DO K=OB%K1+1,OB%K2
             DO J=OB%J1+1,OB%J2
      LOOP1: DO I=OB%I1  ,OB%I2
      IC1 = ICA(I,J,K)
      IC2 = ICA(I+1,J,K)
      IF (SOLID(IC1) .AND. SOLID(IC2))
     .FVX(I,J,K) = -RDXN(I)*(H(I+1,J,K)-H(I,J,K)) + RFODT*UU(I,J,K)
      ENDDO LOOP1
      ENDDO 
      ENDDO 
C
             DO K=OB%K1+1,OB%K2
             DO J=OB%J1  ,OB%J2
      LOOP2: DO I=OB%I1+1,OB%I2
      IC1 = ICA(I,J,K)
      IC2 = ICA(I,J+1,K)
      IF (SOLID(IC1) .AND. SOLID(IC2))
     .FVY(I,J,K) = -RDYN(J)*(H(I,J+1,K)-H(I,J,K)) + RFODT*VV(I,J,K)
      ENDDO LOOP2
      ENDDO 
      ENDDO 
C
             DO K=OB%K1  ,OB%K2
             DO J=OB%J1+1,OB%J2
      LOOP3: DO I=OB%I1+1,OB%I2
      IC1 = ICA(I,J,K)
      IC2 = ICA(I,J,K+1)
      IF (SOLID(IC1) .AND. SOLID(IC2))
     .FVZ(I,J,K) = -RDZN(K)*(H(I,J,K+1)-H(I,J,K)) + RFODT*WW(I,J,K)
      ENDDO LOOP3
      ENDDO 
      ENDDO 
C
      ENDDO OBST_LOOP
C
C Add normal velocity to FVX, etc. for surface cells
C
      WALL_LOOP: DO IW=1,NWC
C
      IF (IV(IW).NE.1 .AND. IV(IW).NE.3) CYCLE WALL_LOOP
C
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IOR = IJKW(4,IW)
C
      BOUNDARY_TYPE: SELECT CASE(IV(IW))
C
      CASE(1)
C
      SELECT CASE(IOR)
      CASE( 1) 
      FVX(II,JJ,KK)   = -RDXN(II)  *(H(II+1,JJ,KK)-H(II,JJ,KK)) + 
     .                   RFODT*(UU(II,JJ,KK) + UWS(IW))
      CASE(-1) 
      FVX(II-1,JJ,KK) = -RDXN(II-1)*(H(II,JJ,KK)-H(II-1,JJ,KK)) + 
     .                   RFODT*(UU(II-1,JJ,KK) - UWS(IW))
      CASE( 2) 
      FVY(II,JJ,KK)   = -RDYN(JJ)  *(H(II,JJ+1,KK)-H(II,JJ,KK)) + 
     .                   RFODT*(VV(II,JJ,KK) + UWS(IW))
      CASE(-2)
      FVY(II,JJ-1,KK) = -RDYN(JJ-1)*(H(II,JJ,KK)-H(II,JJ-1,KK)) + 
     .                   RFODT*(VV(II,JJ-1,KK) - UWS(IW))
      CASE( 3) 
      FVZ(II,JJ,KK)   = -RDZN(KK)  *(H(II,JJ,KK+1)-H(II,JJ,KK)) + 
     .                   RFODT*(WW(II,JJ,KK) + UWS(IW))
      CASE(-3) 
      FVZ(II,JJ,KK-1) = -RDZN(KK-1)*(H(II,JJ,KK)-H(II,JJ,KK-1)) + 
     .                   RFODT*(WW(II,JJ,KK-1) - UWS(IW))
      END SELECT
C
      CASE(3)
C
      SELECT CASE(IOR)
      CASE( 1) ; FVX(II  ,JJ,KK) = 0.
      CASE(-1) ; FVX(II-1,JJ,KK) = 0.
      CASE( 2) ; FVY(II  ,JJ,KK) = 0.
      CASE(-2) ; FVY(II,JJ-1,KK) = 0.
      CASE( 3) ; FVZ(II  ,JJ,KK) = 0.
      CASE(-3) ; FVZ(II,JJ,KK-1) = 0.
      END SELECT
C
      END SELECT BOUNDARY_TYPE
C
      ENDDO WALL_LOOP
C
C
      END SUBROUTINE NO_FLUX
C
C
      SUBROUTINE VELOCITY_PREDICTOR(T,NM,ISTOP)
C
      REAL(EB) T
      INTEGER ISTOP
      INTEGER, INTENT(IN) :: NM
C
      CALL UNPACK_VAR(NM)
C
C Estimates the velocity components at the next time step (Predictor)
C
      TNOW=SECOND() 
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=0,IBAR
      RHS = -FVX(I,J,K) - RDXN(I)*(H(I+1,J,K)-H(I,J,K))
      US(I,J,K) = U(I,J,K) + DT*RHS
      ENDDO 
      ENDDO 
      ENDDO 
C
      DO K=1,KBAR
      DO J=0,JBAR
      DO I=1,IBAR
      RHS = -FVY(I,J,K) - RDYN(J)*(H(I,J+1,K)-H(I,J,K))
      VS(I,J,K) = V(I,J,K) + DT*RHS
      ENDDO 
      ENDDO 
      ENDDO 
C
      DO K=0,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      RHS = -FVZ(I,J,K) - RDZN(K)*(H(I,J,K+1)-H(I,J,K))
      WS(I,J,K) = W(I,J,K) + DT*RHS
      ENDDO 
      ENDDO 
      ENDDO 
C
C Check the stability criteria
C
      DTOLD = DT
      CALL CHECK_STABILITY
C
      IF (DT.LT.DTINT*1.E-4) THEN
         ISTOP = 1
         RETURN
         ENDIF
C
      IF (NEW_TIME_STEP) RETURN
C
      CALL VELOCITY_BC(T)
C
      TUSED(4,NM)=TUSED(4,NM)+SECOND()-TNOW
      END SUBROUTINE VELOCITY_PREDICTOR
C
C
      SUBROUTINE VELOCITY_CORRECTOR(T,NM)
C
      REAL(EB) T
      INTEGER, INTENT(IN) :: NM
C
      CALL UNPACK_VAR(NM)
C
C Correct the velocity components
C
      TNOW=SECOND() 
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=0,IBAR
      RHS = -FVX(I,J,K) - RDXN(I)*(H(I+1,J,K)-H(I,J,K))
      U(I,J,K) = .5*(U(I,J,K) + US(I,J,K) + DT*RHS)
      ENDDO 
      ENDDO 
      ENDDO 
C
      DO K=1,KBAR
      DO J=0,JBAR
      DO I=1,IBAR
      RHS = -FVY(I,J,K) - RDYN(J)*(H(I,J+1,K)-H(I,J,K))
      V(I,J,K) = .5*(V(I,J,K) + VS(I,J,K) + DT*RHS)
      ENDDO 
      ENDDO 
      ENDDO 
C
      DO K=0,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      RHS = -FVZ(I,J,K) - RDZN(K)*(H(I,J,K+1)-H(I,J,K))
      W(I,J,K) = .5*(W(I,J,K) + WS(I,J,K) + DT*RHS)
      ENDDO 
      ENDDO 
      ENDDO 
C
      CALL VELOCITY_BC(T)
C
      TUSED(4,NM)=TUSED(4,NM)+SECOND()-TNOW
      END SUBROUTINE VELOCITY_CORRECTOR
C
C
      SUBROUTINE VELOCITY_BC(T)
C
C Assert tangential velocity boundary conditions
C
      REAL(EB) BC,MUA,T,FVT
      INTEGER IBC,ITAU,NOM
C
      IF (PREDICTOR) THEN
         UU => US
         VV => VS
         WW => WS
         ELSE
         UU => U
         VV => V
         WW => W
         ENDIF
C
      EDGE_LOOP: DO IE=1,NEDGES
C
      II  = IJKE(1,IE)
      JJ  = IJKE(2,IE)
      KK  = IJKE(3,IE)
      IOR = IJKE(4,IE)
      IBC = IJKE(5,IE)
      NOM = IJKE(7,IE)
      BC  = BCV(IBC)
C
      IF (BC.EQ.2.) THEN
      IF (TAUV(IBC).LT.1000000.) ITAU = IBC
      IF (TAUV(IBC).GT.1000000.) ITAU = -(TAUV(IBC)-1000000.)
      IF (TAUV(IBC).EQ.      0.) ITAU = 999
      FVT = FV(T-TIGNS(IBC),ITAU)
      ENDIF
C
      SELECT CASE(IOR)
C
      CASE(1)
      VP = VV(II,JJ,KK+1)
      VM = VV(II,JJ,KK)
      WP = WW(II,JJ+1,KK)
      WM = WW(II,JJ,KK)
      SELECT CASE(IJKE(6,IE))
      CASE( 0) 
         OME_E(IE) = 0. 
         TAU_E(IE) = 0. 
         CYCLE EDGE_LOOP
      CASE(-2) 
         IF (BC.NE.2.) THEN ; WP = BC*WM ; ELSE 
                              WP = FVT*VEL_TS(2,IBC) ; ENDIF
         MUA = .5*(MU(II,JJ,KK)+MU(II,JJ,KK+1))
         IF (IBC.EQ.OPEN_INDEX .OR. IBC.EQ.MIRROR_INDEX) 
     .      WW(II,JJ+1,KK) = WP
         IF (IBC.EQ.INTERPOLATED_INDEX) THEN
            WP = OMESH(NOM)%W(IJKE(8,IE),IJKE(9,IE),IJKE(10,IE))
            WW(II,JJ+1,KK) = WP
            ENDIF
      CASE( 2) 
         IF (BC.NE.2.) THEN ; WM = BC*WP ; ELSE
                              WM = FVT*VEL_TS(2,IBC) ; ENDIF
         MUA = .5*(MU(II,JJ+1,KK)+MU(II,JJ+1,KK+1))
         IF (IBC.EQ.OPEN_INDEX .OR. IBC.EQ.MIRROR_INDEX) 
     .      WW(II,JJ,KK) = WM
         IF (IBC.EQ.INTERPOLATED_INDEX) THEN
            WM = OMESH(NOM)%W(IJKE(8,IE),IJKE(9,IE),IJKE(10,IE))
            WW(II,JJ,KK) = WM
            ENDIF
      CASE(-3) 
         IF (BC.NE.2.) THEN ; VP = BC*VM ; ELSE
                              VP = FVT*VEL_TS(2,IBC) ; ENDIF
         MUA = .5*(MU(II,JJ,KK)+MU(II,JJ+1,KK))
         IF (IBC.EQ.OPEN_INDEX .OR. IBC.EQ.MIRROR_INDEX) 
     .      VV(II,JJ,KK+1) = VP
         IF (IBC.EQ.INTERPOLATED_INDEX) THEN
            VP = OMESH(NOM)%V(IJKE(8,IE),IJKE(9,IE),IJKE(10,IE))
            VV(II,JJ,KK+1) = VP
            ENDIF
      CASE( 3) 
         IF (BC.NE.2.) THEN ; VM = BC*VP ; ELSE
                              VM = FVT*VEL_TS(2,IBC) ; ENDIF 
         MUA = .5*(MU(II,JJ,KK+1)+MU(II,JJ+1,KK+1))
         IF (IBC.EQ.OPEN_INDEX .OR. IBC.EQ.MIRROR_INDEX) 
     .      VV(II,JJ,KK) = VM
         IF (IBC.EQ.INTERPOLATED_INDEX) THEN
            VM = OMESH(NOM)%V(IJKE(8,IE),IJKE(9,IE),IJKE(10,IE))
            VV(II,JJ,KK) = VM
            ENDIF
      END SELECT
      DVDZ = RDZN(KK)*(VP-VM)
      DWDY = RDYN(JJ)*(WP-WM)
      OME_E(IE) = DWDY - DVDZ
      TAU_E(IE) = MUA*(DVDZ + DWDY)
C
      CASE(2)
      UP = UU(II,JJ,KK+1)
      UM = UU(II,JJ,KK)
      WP = WW(II+1,JJ,KK)
      WM = WW(II,JJ,KK)
      SELECT CASE(IJKE(6,IE))
      CASE( 0) 
         OME_E(IE) = 0. 
         TAU_E(IE) = 0. 
         CYCLE EDGE_LOOP
      CASE(-1) 
         IF (BC.NE.2.) THEN ; WP = BC*WM ; ELSE
                              WP = FVT*VEL_TS(2,IBC) ; ENDIF
         MUA = .5*(MU(II,JJ,KK)+MU(II,JJ,KK+1))
         IF (IBC.EQ.OPEN_INDEX .OR. IBC.EQ.MIRROR_INDEX) 
     .      WW(II+1,JJ,KK) = WP
         IF (IBC.EQ.INTERPOLATED_INDEX) THEN
            WP = OMESH(NOM)%W(IJKE(8,IE),IJKE(9,IE),IJKE(10,IE))
            WW(II+1,JJ,KK) = WP
            ENDIF
      CASE( 1) 
         IF (BC.NE.2.) THEN ; WM = BC*WP ; ELSE
                              WM = FVT*VEL_TS(2,IBC) ; ENDIF
         MUA = .5*(MU(II+1,JJ,KK)+MU(II+1,JJ,KK+1))
         IF (IBC.EQ.OPEN_INDEX .OR. IBC.EQ.MIRROR_INDEX) 
     .      WW(II,JJ,KK) = WM
         IF (IBC.EQ.INTERPOLATED_INDEX) THEN
            WM = OMESH(NOM)%W(IJKE(8,IE),IJKE(9,IE),IJKE(10,IE))
            WW(II,JJ,KK) = WM
            ENDIF
      CASE(-3) 
         IF (BC.NE.2.) THEN ; UP = BC*UM ; ELSE
                              UP = FVT*VEL_TS(1,IBC) ; ENDIF
         MUA = .5*(MU(II,JJ,KK)+MU(II+1,JJ,KK))
         IF (IBC.EQ.OPEN_INDEX .OR. IBC.EQ.MIRROR_INDEX) 
     .      UU(II,JJ,KK+1) = UP
         IF (IBC.EQ.INTERPOLATED_INDEX) THEN
            UP = OMESH(NOM)%U(IJKE(8,IE),IJKE(9,IE),IJKE(10,IE))
            UU(II,JJ,KK+1) = UP
            ENDIF
      CASE( 3) 
         IF (BC.NE.2.) THEN ; UM = BC*UP ; ELSE
                              UM = FVT*VEL_TS(1,IBC) ; ENDIF
         MUA = .5*(MU(II,JJ,KK+1)+MU(II+1,JJ,KK+1))
         IF (IBC.EQ.OPEN_INDEX .OR. IBC.EQ.MIRROR_INDEX) 
     .      UU(II,JJ,KK) = UM
         IF (IBC.EQ.INTERPOLATED_INDEX) THEN
            UM = OMESH(NOM)%U(IJKE(8,IE),IJKE(9,IE),IJKE(10,IE))
            UU(II,JJ,KK) = UM
            ENDIF
      END SELECT
      DUDZ = RDZN(KK)*(UP-UM)
      DWDX = RDXN(II)*(WP-WM)
      OME_E(IE) = DUDZ - DWDX
      TAU_E(IE) = MUA*(DUDZ + DWDX)
C
      CASE(3)
      UP = UU(II,JJ+1,KK)
      UM = UU(II,JJ,KK)
      VP = VV(II+1,JJ,KK)
      VM = VV(II,JJ,KK)
      SELECT CASE(IJKE(6,IE))
      CASE( 0) 
         OME_E(IE) = 0. 
         TAU_E(IE) = 0. 
         CYCLE EDGE_LOOP
      CASE(-1) 
         IF (BC.NE.2.) THEN ; VP = BC*VM ; ELSE
                              VP = FVT*VEL_TS(1,IBC) ; ENDIF
         MUA = .5*(MU(II,JJ,KK)+MU(II,JJ+1,KK))
         IF (IBC.EQ.OPEN_INDEX .OR. IBC.EQ.MIRROR_INDEX) 
     .      VV(II+1,JJ,KK) = VP
         IF (IBC.EQ.INTERPOLATED_INDEX) THEN
            VP = OMESH(NOM)%V(IJKE(8,IE),IJKE(9,IE),IJKE(10,IE))
            VV(II+1,JJ,KK) = VP
            ENDIF
      CASE( 1) 
         IF (BC.NE.2.) THEN ; VM = BC*VP ; ELSE
                              VM = FVT*VEL_TS(1,IBC) ; ENDIF
         MUA = .5*(MU(II+1,JJ,KK)+MU(II+1,JJ+1,KK))
         IF (IBC.EQ.OPEN_INDEX .OR. IBC.EQ.MIRROR_INDEX) 
     .      VV(II,JJ,KK) = VM
         IF (IBC.EQ.INTERPOLATED_INDEX) THEN
            VM = OMESH(NOM)%V(IJKE(8,IE),IJKE(9,IE),IJKE(10,IE))
            VV(II,JJ,KK) = VM
            ENDIF
      CASE(-2) 
         IF (BC.NE.2.) THEN ; UP = BC*UM ; ELSE
                              UP = FVT*VEL_TS(1,IBC) ; ENDIF
         MUA = .5*(MU(II,JJ,KK)+MU(II+1,JJ,KK))
         IF (IBC.EQ.OPEN_INDEX .OR. IBC.EQ.MIRROR_INDEX) 
     .      UU(II,JJ+1,KK) = UP
         IF (IBC.EQ.INTERPOLATED_INDEX) THEN
            UP = OMESH(NOM)%U(IJKE(8,IE),IJKE(9,IE),IJKE(10,IE))
            UU(II,JJ+1,KK) = UP
            ENDIF
      CASE( 2) 
         IF (BC.NE.2.) THEN ; UM = BC*UP ; ELSE
                              UM = FVT*VEL_TS(1,IBC) ; ENDIF
         MUA = .5*(MU(II,JJ+1,KK)+MU(II+1,JJ+1,KK))
         IF (IBC.EQ.OPEN_INDEX .OR. IBC.EQ.MIRROR_INDEX) 
     .      UU(II,JJ,KK) = UM
         IF (IBC.EQ.INTERPOLATED_INDEX) THEN
            UM = OMESH(NOM)%U(IJKE(8,IE),IJKE(9,IE),IJKE(10,IE))
            UU(II,JJ,KK) = UM
            ENDIF
      END SELECT
      DVDX = RDXN(II)*(VP-VM)
      DUDY = RDYN(JJ)*(UP-UM)
      OME_E(IE) = DVDX - DUDY
      TAU_E(IE) = MUA*(DVDX + DUDY)
C
      END SELECT
      ENDDO EDGE_LOOP
C
      END SUBROUTINE VELOCITY_BC 
C
C
      SUBROUTINE CHECK_STABILITY
C
C Checks the Courant/Von Neumann stability criterion, and if necessary, 
C reduces the time step accordingly
C
      REAL(EB) UODX,VODY,WODZ,UVW,UVWMAX,R_DX2,MU_MAX,MUTRM
C
      NEW_TIME_STEP = .FALSE.
      UVWMAX = 0.
      VN     = 0.
      MUTRM  = 1.E-9
      R_DX2  = 1.E-9
C
C     Determine max CFL number
C
      DO K=0,KBAR
      DO J=0,JBAR
      DO I=0,IBAR
      UODX = ABS(US(I,J,K))*RDXN(I)
      VODY = ABS(VS(I,J,K))*RDYN(J)
      WODZ = ABS(WS(I,J,K))*RDZN(K)
      UVW  = MAX(UODX,VODY,WODZ)
      IF (UVW.GE.UVWMAX) THEN
         UVWMAX = UVW 
         ICFL=I
         JCFL=J
         KCFL=K
         ENDIF
      ENDDO 
      ENDDO   
      ENDDO   
C
      CFL = DT*UVWMAX
C
C     Determine max Von Neumann Number
C
      PARABOLIC_IF: IF (DNS .OR. DXMIN.LT.0.005) THEN
C
      INCOMPRESSIBLE_IF: IF (INCOMPRESSIBLE) THEN
C
      IF (TWO_D) THEN
      R_DX2 = 1./DXMIN**2 + 1./DZMIN**2
      ELSE
      R_DX2 = 1./DXMIN**2 + 1./DYMIN**2 + 1./DZMIN**2
      ENDIF
C
      MUTRM = MAX(RPR,RSC)*MU_SPEC(0,NINT(0.1*TMPA))/RHOA
C
      ELSE INCOMPRESSIBLE_IF
C
      MU_MAX = 0.   
      DO K=1,KBAR
      DO J=1,JBAR
      IILOOP: DO I=1,IBAR
      IF (SOLID(ICA(I,J,K))) CYCLE IILOOP
      IF (MU(I,J,K).GE.MU_MAX) THEN
         MU_MAX = MU(I,J,K)
         I_VN=I
         J_VN=J
         K_VN=K
         ENDIF
      ENDDO IILOOP
      ENDDO  
      ENDDO  
C
      IF (TWO_D) THEN
         R_DX2 = RDX(I_VN)**2 + RDZ(K_VN)**2
         ELSE
         R_DX2 = RDX(I_VN)**2 + RDY(J_VN)**2 + RDZ(K_VN)**2
         ENDIF
      MUTRM = MAX(RPR,RSC)*MU_MAX/RHOS(I_VN,J_VN,K_VN)
C
      ENDIF INCOMPRESSIBLE_IF
C
      VN = 2.*DT*R_DX2*MUTRM
C
      ENDIF PARABOLIC_IF
C
C     Adjust time step size if necessary
C
      IF (CFL.LT.1.0 .AND. VN.LT.1.0) THEN
C
      DTNEXT = DT
      IF (CFL.LE.0.8 .AND. VN.LT.0.8) DTNEXT = MIN(1.1*DT,DTINT)
C
      ELSE
C
      IF (UVWMAX.EQ.0.) UVWMAX = 1.
      DT    = 0.9*MIN( 1.0/UVWMAX , 1.0/(2.*R_DX2*MUTRM) )
      NEW_TIME_STEP = .TRUE.
C
      ENDIF
C
      END SUBROUTINE CHECK_STABILITY
C
C
      SUBROUTINE BAROCLINIC_CORRECTION
C
      HQS  => WORK1
      RTRM => WORK2
C
      IF (PREDICTOR) THEN
         UU => U
         VV => V
         WW => W
         ELSE
         UU => US
         VV => VS
         WW => WS
         ENDIF
C
      RHO_AVG_OLD = RHO_AVG
      RMIN =  1000.
      RMAX = -1000.
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      IF (.NOT.SOLID(ICA(I,J,K))) THEN
      RMIN = MIN(RHOP(I,J,K),RMIN)
      RMAX = MAX(RHOP(I,J,K),RMAX)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      RHO_AVG = 2.*RMIN*RMAX/(RMIN+RMAX)
      RRAT = RHO_AVG_OLD/RHO_AVG
      RTRM = (1.-RHO_AVG/RHOP)*RRAT
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      U2 = 0.25*(UU(I,J,K)+UU(I-1,J,K))**2
      V2 = 0.25*(VV(I,J,K)+VV(I,J-1,K))**2
      W2 = 0.25*(WW(I,J,K)+WW(I,J,K-1))**2
      HQS(I,J,K) = 0.5*(U2+V2+W2)
      ENDDO
      ENDDO
      ENDDO
C
      DO K=1,KBAR
      DO J=1,JBAR
      U2 = (1.5*UU(0,J,K)-0.5*UU(1,J,K))**2
      V2 = 0.25*(VV(1,J,K)+VV(1,J-1,K))**2
      W2 = 0.25*(WW(1,J,K)+WW(1,J,K-1))**2
      HQS(0,J,K) = MIN(0.5*(U2+V2+W2),10*DX(1)+HQS(1,J,K))
      U2 = (1.5*UU(IBAR,J,K)-0.5*UU(IBM1,J,K))**2
      V2 = 0.25*(VV(IBAR,J,K)+VV(IBAR,J-1,K))**2
      W2 = 0.25*(WW(IBAR,J,K)+WW(IBAR,J,K-1))**2
      HQS(IBP1,J,K) = MIN(0.5*(U2+V2+W2),10*DX(IBAR)+HQS(IBAR,J,K))
      ENDDO
      ENDDO
C
      IF (.NOT.TWO_D) THEN
      DO K=1,KBAR
      DO I=1,IBAR
      U2 = 0.25*(UU(I,1,K)+UU(I-1,1,K))**2
      V2 = (1.5*VV(I,0,K)-0.5*VV(I,1,K))**2
      W2 = 0.25*(WW(I,1,K)+WW(I,1,K-1))**2
      HQS(I,0,K) = MIN(0.5*(U2+V2+W2),10*DY(1)+HQS(I,1,K))
      U2 = 0.25*(UU(I,JBAR,K)+UU(I-1,JBAR,K))**2
      V2 = (1.5*VV(I,JBAR,K)-0.5*VV(I,JBM1,K))**2
      W2 = 0.25*(WW(I,JBAR,K)+WW(I,JBAR,K-1))**2
      HQS(I,JBP1,K) = MIN(0.5*(U2+V2+W2),10*DY(JBAR)+HQS(I,JBAR,K))
      ENDDO
      ENDDO
      ENDIF
C
      DO J=1,JBAR
      DO I=1,IBAR
      U2 = 0.25*(UU(I,J,1)+UU(I-1,J,1))**2
      V2 = 0.25*(VV(I,J,1)+VV(I,J-1,1))**2
      W2 = (1.5*WW(I,J,0)-0.5*WW(I,J,1))**2
      HQS(I,J,0) = MIN(0.5*(U2+V2+W2),10*DZ(1)+HQS(I,J,1))
      U2 = 0.25*(UU(I,J,KBAR)+UU(I-1,J,KBAR))**2
      V2 = 0.25*(VV(I,J,KBAR)+VV(I,J-1,KBAR))**2
      W2 = (1.5*WW(I,J,KBAR)-0.5*WW(I,J,KBM1))**2
      HQS(I,J,KBP1) = MIN(0.5*(U2+V2+W2),10*DZ(KBAR)+HQS(I,J,KBAR))
      ENDDO
      ENDDO
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=0,IBAR
      FVX(I,J,K) = FVX(I,J,K) - 0.5*(RTRM(I+1,J,K)+RTRM(I,J,K))*
     .        (H(I+1,J,K)-H(I,J,K)-HQS(I+1,J,K)+HQS(I,J,K))*RDXN(I)
      ENDDO
      ENDDO
      ENDDO
C
      IF (.NOT.TWO_D) THEN
      DO K=1,KBAR
      DO J=0,JBAR
      DO I=1,IBAR
      FVY(I,J,K) = FVY(I,J,K) - 0.5*(RTRM(I,J+1,K)+RTRM(I,J,K))*
     .        (H(I,J+1,K)-H(I,J,K)-HQS(I,J+1,K)+HQS(I,J,K))*RDYN(J)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
C
      DO K=0,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      FVZ(I,J,K) = FVZ(I,J,K) - 0.5*(RTRM(I,J,K+1)+RTRM(I,J,K))*
     .        (H(I,J,K+1)-H(I,J,K)-HQS(I,J,K+1)+HQS(I,J,K))*RDZN(K)
      ENDDO
      ENDDO
      ENDDO
C
      END SUBROUTINE BAROCLINIC_CORRECTION
C
C
      END MODULE VELO
