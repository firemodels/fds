      MODULE MASS
C
C Compute the mass equation differences 
C
      USE PREC
      USE CONS
      USE PACKER
      USE PYRO
      USE RAD 
C
      IMPLICIT NONE
      PRIVATE
      REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP
      REAL(EB), POINTER, DIMENSION(:,:,:) :: UDRHODX,VDRHODY,WDRHODZ
      REAL(EB), POINTER, DIMENSION(:,:,:) :: EPSX,EPSY,EPSZ,UU,VV,WW,
     .                                       RHOP,DP
      REAL(EB) H_DNS,H_NATURAL,H_FORCED,VELCON,DT_BC,U2,V2,W2
      INTEGER IYY,ITMP
      PUBLIC MASS_FLUX,DENSITY
C
C
      CONTAINS
C
      SUBROUTINE MASS_FLUX(NM)
C
      INTEGER, INTENT(IN) :: NM
      REAL(EB) FXYZ,PMDT,UDRHODN
      INTEGER I,J,K,N,II,JJ,KK,IIG,JJG,KKG,IW,IOR
      REAL(EB), POINTER, DIMENSION(:) :: UWP
C
      IF (EVACUATION_ONLY(NM)) RETURN
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
         UWP  => UW
         PMDT = DT
         ELSE
         UU => US
         VV => VS
         WW => WS
         DP => DS
         RHOP => RHOS
         UWP  => UWS
         PMDT = -DT
         ENDIF
C
C Define local CFL numbers
C
      EPSX => WORK1
      EPSY => WORK2
      EPSZ => WORK3
C
      DO K=0,KBAR
      DO J=0,JBAR
      DO I=0,IBAR
      EPSX(I,J,K) = PMDT*UU(I,J,K)*RDXN(I)
      EPSY(I,J,K) = PMDT*VV(I,J,K)*RDYN(J)
      EPSZ(I,J,K) = PMDT*WW(I,J,K)*RDZN(K)
      ENDDO
      ENDDO
      ENDDO
C
C Compute spatial differences for density equation
C
      NOT_ISOTHERMAL_IF: IF (.NOT.ISOTHERMAL) THEN
C
      UDRHODX => WORK4
      VDRHODY => WORK5
      WDRHODZ => WORK6
C
      DO K=0,KBAR
      DO J=0,JBAR
      DO I=0,IBAR
      UDRHODX(I,J,K) = UU(I,J,K)*(RHOP(I+1,J,K)-RHOP(I,J,K))*RDXN(I)
      VDRHODY(I,J,K) = VV(I,J,K)*(RHOP(I,J+1,K)-RHOP(I,J,K))*RDYN(J)
      WDRHODZ(I,J,K) = WW(I,J,K)*(RHOP(I,J,K+1)-RHOP(I,J,K))*RDZN(K)
      ENDDO
      ENDDO
      ENDDO
C
      WLOOP: DO IW=1,NWC
      IF (IV(IW).EQ.0) CYCLE WLOOP
      II  = IJKW(1,IW) ; IIG = IJKW(6,IW)
      JJ  = IJKW(2,IW) ; JJG = IJKW(7,IW)
      KK  = IJKW(3,IW) ; KKG = IJKW(8,IW)
      IOR = IJKW(4,IW)
      UDRHODN = UWP(IW)*(RHO_W(IW)-RHOP(IIG,JJG,KKG))*RDN(IW)
      SELECT CASE(IOR)
      CASE( 1) ; UDRHODX(II,JJ,KK)   = UDRHODN
      CASE(-1) ; UDRHODX(II-1,JJ,KK) = UDRHODN
      CASE( 2) ; VDRHODY(II,JJ,KK)   = UDRHODN
      CASE(-2) ; VDRHODY(II,JJ-1,KK) = UDRHODN
      CASE( 3) ; WDRHODZ(II,JJ,KK)   = UDRHODN
      CASE(-3) ; WDRHODZ(II,JJ,KK-1) = UDRHODN
      END SELECT
      ENDDO WLOOP
C
      DO K=1,KBAR
      DO J=1,JBAR
      FLUX_LOOP: DO I=1,IBAR
      FXYZ   = .5*(UDRHODX(I,J,K)  *(1.-EPSX(I,J,K))   + 
     .             UDRHODX(I-1,J,K)*(1.+EPSX(I-1,J,K)) +
     .             VDRHODY(I,J,K)  *(1.-EPSY(I,J,K))   + 
     .             VDRHODY(I,J-1,K)*(1.+EPSY(I,J-1,K)) +
     .             WDRHODZ(I,J,K)  *(1.-EPSZ(I,J,K))   + 
     .             WDRHODZ(I,J,K-1)*(1.+EPSZ(I,J,K-1)) )
      FRHO(I,J,K) = FXYZ + RHOP(I,J,K)*DP(I,J,K)
C
      ENDDO FLUX_LOOP
      ENDDO
      ENDDO
C
      ENDIF NOT_ISOTHERMAL_IF
C
C Compute the species equation differences
C
      IF (NSPEC.GT.0) THEN
      IF (PREDICTOR) YYP => YY
      IF (CORRECTOR) YYP => YYS
      UDRHODX => WORK4
      VDRHODY => WORK5
      WDRHODZ => WORK6
      ENDIF
C
      SPECIESLOOP: DO N=1,NSPEC
C
      DO K=0,KBAR
      DO J=0,JBAR
      DO I=0,IBAR
      UDRHODX(I,J,K) = UU(I,J,K)*(RHOP(I+1,J,K)*YYP(I+1,J,K,N)-
     .                            RHOP(I,J,K)  *YYP(I,J,K,N)  )*RDXN(I)
      VDRHODY(I,J,K) = VV(I,J,K)*(RHOP(I,J+1,K)*YYP(I,J+1,K,N)-
     .                            RHOP(I,J,K)  *YYP(I,J,K,N)  )*RDYN(J)
      WDRHODZ(I,J,K) = WW(I,J,K)*(RHOP(I,J,K+1)*YYP(I,J,K+1,N)-
     .                            RHOP(I,J,K)  *YYP(I,J,K,N)  )*RDZN(K)
      ENDDO
      ENDDO
      ENDDO
C
C Correct U d(RHO*Y)/dx etc. on boundaries
C
      WLOOP2: DO IW=1,NWC
      IF (IV(IW).EQ.0) CYCLE WLOOP2
      II  = IJKW(1,IW) ; IIG = IJKW(6,IW)
      JJ  = IJKW(2,IW) ; JJG = IJKW(7,IW)
      KK  = IJKW(3,IW) ; KKG = IJKW(8,IW)
      IOR = IJKW(4,IW)
      UDRHODN = UWP(IW)*(RHO_W(IW)*YY_W(IW,N)-
     .                   RHOP(IIG,JJG,KKG)*YYP(IIG,JJG,KKG,N))*RDN(IW)
      SELECT CASE(IOR)
      CASE( 1) ; UDRHODX(II,JJ,KK)   = UDRHODN
      CASE(-1) ; UDRHODX(II-1,JJ,KK) = UDRHODN
      CASE( 2) ; VDRHODY(II,JJ,KK)   = UDRHODN
      CASE(-2) ; VDRHODY(II,JJ-1,KK) = UDRHODN
      CASE( 3) ; WDRHODZ(II,JJ,KK)   = UDRHODN
      CASE(-3) ; WDRHODZ(II,JJ,KK-1) = UDRHODN
      END SELECT
      ENDDO WLOOP2
C
C Construct weighted averages of U d(RHO*Y)/dx etc. in grid cells
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      FXYZ   = .5*(UDRHODX(I,J,K)  *(1.-EPSX(I,J,K))   + 
     .             UDRHODX(I-1,J,K)*(1.+EPSX(I-1,J,K)) +
     .             VDRHODY(I,J,K)  *(1.-EPSY(I,J,K))   + 
     .             VDRHODY(I,J-1,K)*(1.+EPSY(I,J-1,K)) +
     .             WDRHODZ(I,J,K)  *(1.-EPSZ(I,J,K))   + 
     .             WDRHODZ(I,J,K-1)*(1.+EPSZ(I,J,K-1)) )
      FYY(I,J,K,N) = FYY(I,J,K,N) + FXYZ + 
     .               RHOP(I,J,K)*YYP(I,J,K,N)*DP(I,J,K) 
      ENDDO
      ENDDO
      ENDDO
C
      ENDDO SPECIESLOOP
C
      TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
C
      END SUBROUTINE MASS_FLUX
C
C
      SUBROUTINE DENSITY(T,NM)
C
C Update the density
C
      REAL(EB) T,WFAC,DTRATIO,OMDTRATIO
      INTEGER I,J,K,N,II,JJ,KK,IW,IOR
      REAL(EB), POINTER, DIMENSION(:,:,:) :: RCON_MF
      INTEGER, INTENT(IN) :: NM
C
      IF (EVACUATION_ONLY(NM)) RETURN
C
      CALL UNPACK_VAR(NM)
C
      TNOW=SECOND()
C
      PREDICTOR_CORRECTOR: SELECT CASE(PREDICTOR)
C
      CASE(.TRUE.)     ! Predictor Step
C
      IF (.NOT.NEW_TIME_STEP) THEN
C
      DO N=1,NSPEC
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      YYS(I,J,K,N) = RHO(I,J,K)*YY(I,J,K,N) - DT*FYY(I,J,K,N)
      ENDDO 
      ENDDO
      ENDDO
      ENDDO
C
      ELSE
C
      DTRATIO   = DT/DTOLD
      OMDTRATIO = 1. - DTRATIO
C
      DO N=1,NSPEC
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      YYS(I,J,K,N) = OMDTRATIO*RHO(I,J,K) *YY(I,J,K,N) + 
     .                 DTRATIO*RHOS(I,J,K)*YYS(I,J,K,N)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      ENDIF
C
C Compute density
C
      IF (.NOT.ISOTHERMAL) THEN
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      RHOS(I,J,K) = RHO(I,J,K)-DT*FRHO(I,J,K)
      ENDDO 
      ENDDO
      ENDDO
C
      ELSE
C
      RHOS = P0S/(TMPA*RCON(0))
      DO N=1,NSPEC
      WFAC = 1. - RCON(N)/RCON(0)
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      RHOS(I,J,K) = RHOS(I,J,K) + WFAC*YYS(I,J,K,N)
      ENDDO 
      ENDDO
      ENDDO
      ENDDO
C
      ENDIF
C
C Correct densities above or below clip limits
C
      CALL CHECK_DENSITY
C
C Extract mass fraction from RHO * YY
C
      DO N=1,NSPEC
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      YYS(I,J,K,N) = YYS(I,J,K,N)/RHOS(I,J,K)
      ENDDO 
      ENDDO
      ENDDO
      ENDDO
C
C Correct mass fractions above or below clip limits
C
      CALL CHECK_MASS_FRACTION
C
C Predict background pressure at next time step
C
      P0S = P0 + DP0DT*DT
C
C Compute mixture molecular weight
C
      IF (NSPEC.GT.0 .AND. .NOT.MIXTURE_FRACTION) THEN
      RSUM = RCON(0)
      DO N=1,NSPEC
      WFAC = RCON(N) - RCON(0)
      RSUM(:,:,:) = RSUM(:,:,:) + WFAC*YYS(:,:,:,N)
      ENDDO
      ENDIF
C
      IF (MIXTURE_FRACTION) THEN
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      IYY = MAX(0,NINT(YYS(I,J,K,IFUEL)*10000.))
      IYY = MIN(10000,IYY)
      RSUM(I,J,K) = RSUM_MF(IYY)
      ENDDO
      ENDDO
      ENDDO
C
      IF (NSPEC.GT.1) THEN
         RCON_MF => WORK2
         RCON_MF =  RSUM
         ENDIF
C
      SPEC_LOOP: DO N=1,NSPEC
      IF (N.EQ.IFUEL) CYCLE SPEC_LOOP
      RSUM(:,:,:)=RSUM(:,:,:) + YYS(:,:,:,N)*(RCON(N)-RCON_MF(:,:,:))
      ENDDO SPEC_LOOP
C
      ENDIF
C
C Compute predicted temperature
C
      IF (.NOT.ISOTHERMAL) THEN
C
      IF (NSPEC.EQ.0) THEN
         TMP = P0S/(RCON(0)*RHOS)
      ELSE
         TMP = P0S/(RSUM*RHOS)
      ENDIF
C
      TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))
C
      ENDIF
C
C
C *** Corrector Step ***
C
      CASE(.FALSE.)    ! Corrector Step
C
      DO N=1,NSPEC
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      YY(I,J,K,N) = .5*(RHO(I,J,K) *YY(I,J,K,N)  +
     .                  RHOS(I,J,K)*YYS(I,J,K,N) -
     .                  DT*FYY(I,J,K,N) ) 
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      IF (.NOT.ISOTHERMAL) THEN
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      RHO(I,J,K)= .5*(RHO(I,J,K)+RHOS(I,J,K)-DT*FRHO(I,J,K))
      ENDDO
      ENDDO
      ENDDO
C
      ELSE
C
      RHO = P0/(RCON(0)*TMPA)
      DO N=1,NSPEC
      WFAC = 1. - RCON(N)/RCON(0)
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      RHO(I,J,K) = RHO(I,J,K) + WFAC*YY(I,J,K,N)
      ENDDO 
      ENDDO
      ENDDO
      ENDDO
C
      ENDIF
C
C Correct densities above or below clip limits
C
      CALL CHECK_DENSITY
C
C Extract Y_n from rho*Y_n
C
      DO N=1,NSPEC
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      YY(I,J,K,N) = YY(I,J,K,N)/RHO(I,J,K)
      ENDDO 
      ENDDO
      ENDDO
      ENDDO
C
C Correct mass fractions above or below clip limits
C
      CALL CHECK_MASS_FRACTION
C
C Update background pressure
C
      P0 = .5*(P0 + P0S + DP0DTS*DT)
C
C Compute mixture molecular weight
C
      IF (NSPEC.GT.0 .AND. .NOT. MIXTURE_FRACTION) THEN
      RSUM = RCON(0)
      DO N=1,NSPEC
      WFAC = RCON(N) - RCON(0)
      RSUM(:,:,:) = RSUM(:,:,:) + WFAC*YY(:,:,:,N)
      ENDDO
      ENDIF
C
      IF (MIXTURE_FRACTION) THEN
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      IYY = MAX(0,NINT(YY(I,J,K,IFUEL)*10000.))
      IYY = MIN(10000,IYY)
      RSUM(I,J,K) = RSUM_MF(IYY)
      ENDDO
      ENDDO
      ENDDO
C
      IF (NSPEC.GT.1) THEN
         RCON_MF => WORK2
         RCON_MF =  RSUM
         ENDIF
C
      SPEC_LOOP_2: DO N=1,NSPEC
      IF (N.EQ.IFUEL) CYCLE SPEC_LOOP_2
      RSUM(:,:,:)=RSUM(:,:,:) + YY(:,:,:,N)*(RCON(N)-RCON_MF(:,:,:))
      ENDDO SPEC_LOOP_2
C
      ENDIF
C
C Compute temperature
C
      IF (.NOT.ISOTHERMAL) THEN
C
      IF (NSPEC.EQ.0) THEN
         TMP = P0/(RCON(0)*RHO)
      ELSE
         TMP = P0/(RSUM*RHO)
      ENDIF
C
      TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))
C
      ENDIF
C
      END SELECT PREDICTOR_CORRECTOR
C
C Apply thermal, material and density boundary conditions
C
      PREDICTOR_OR_CORRECTOR: IF (PREDICTOR) THEN
C
      IF (.NOT.ISOTHERMAL) CALL TEMPERATURE_BC(T,NM)
      IF (NSPEC.GT.0) CALL SPECIES_BC(T)
      CALL DENSITY_BC
C
      ELSE PREDICTOR_OR_CORRECTOR
C
      IF (COMBUSTION_MODEL.EQ.3) CALL COMBUSTION(NM)
C
      IF (.NOT.ISOTHERMAL) CALL TEMPERATURE_BC(T,NM)
      IF (NSPEC.GT.0) CALL SPECIES_BC(T)
      CALL DENSITY_BC
C
      IF (COMBUSTION_MODEL.EQ.2) CALL COMBUSTION(NM)
C
      IF (RADIATION) THEN
         CALL RADIATION_FVM(NM)
      ELSE
         QR = -CHI_R*Q
      ENDIF
C
      ENDIF PREDICTOR_OR_CORRECTOR
C
      TUSED(3,NM)=TUSED(3,NM)+SECOND()-TNOW
C
C
      CONTAINS
C
      SUBROUTINE CHECK_DENSITY
C
C Redistribute mass from cells below the lower density cut off
C
      REAL(EB) SUM,ISUM,CONST,RHOMI,RHOPI,RHOMJ,RHOPJ,RHOMK,RHOPK,
     .         RHO00
      INTEGER IC
      LOGICAL LC(-3:3)
C
      IF (PREDICTOR) THEN
         RHOP=>RHOS
      ELSE
         RHOP=>RHO
      ENDIF
C
      DO K=1,KBAR
      DO J=1,JBAR
      CHECK_LOOP: DO I=1,IBAR
      IF (RHOP(I,J,K).GE.RHOMAX) RHOP(I,J,K) = RHOMAX
      IF (RHOP(I,J,K).GE.RHOMIN) CYCLE CHECK_LOOP
      IC = ICA(I,J,K)
      IF (SOLID(IC)) THEN
         RHOP(I,J,K) = MAX(RHOP(I,J,K),RHOMIN)
         CYCLE CHECK_LOOP
         ENDIF
      SUM = 0.
      ISUM = 0.
      LC = .FALSE.
      RHO00 = RHOP(I,J,K)
      RHOMI = RHOP(I-1,J,K)
      RHOPI = RHOP(I+1,J,K)
      RHOMJ = RHOP(I,J-1,K)
      RHOPJ = RHOP(I,J+1,K)
      RHOMK = RHOP(I,J,K-1)
      RHOPK = RHOP(I,J,K+1)
      IF (IWA(IC,-1).EQ.0 .AND. RHOMI.GT.RHOMIN) THEN
         SUM = SUM + RHOMI
         ISUM = ISUM + 1.
         LC(-1) = .TRUE.
         ENDIF
      IF (IWA(IC, 1).EQ.0 .AND. RHOPI.GT.RHOMIN) THEN
         SUM = SUM + RHOPI
         ISUM = ISUM + 1.
         LC( 1) = .TRUE.
         ENDIF
      IF (IWA(IC,-2).EQ.0 .AND. RHOMJ.GT.RHOMIN) THEN
         SUM = SUM + RHOMJ
         ISUM = ISUM + 1.
         LC(-2) = .TRUE.
         ENDIF
      IF (IWA(IC, 2).EQ.0 .AND. RHOPJ.GT.RHOMIN) THEN
         SUM = SUM + RHOPJ
         ISUM = ISUM + 1.
         LC( 2) = .TRUE.
         ENDIF
      IF (IWA(IC,-3).EQ.0 .AND. RHOMK.GT.RHOMIN) THEN
         SUM = SUM + RHOMK
         ISUM = ISUM + 1.
         LC(-3) = .TRUE.
         ENDIF
      IF (IWA(IC, 3).EQ.0 .AND. RHOPK.GT.RHOMIN) THEN
         SUM = SUM + RHOPK
         ISUM = ISUM + 1.
         LC( 3) = .TRUE.
         ENDIF
      IF (ISUM.EQ.0) THEN
         RHOP(I,J,K) = RHOMIN
         CYCLE CHECK_LOOP
         ENDIF
      CONST = (RHOMIN-RHO00)/(SUM-ISUM*RHO00)
      IF (LC(-1)) RHOP(I-1,J,K) = MAX(RHOMIN,RHOMI+CONST*(RHO00-RHOMI))
      IF (LC( 1)) RHOP(I+1,J,K) = MAX(RHOMIN,RHOPI+CONST*(RHO00-RHOPI))
      IF (LC(-2)) RHOP(I,J-1,K) = MAX(RHOMIN,RHOMJ+CONST*(RHO00-RHOMJ))
      IF (LC( 2)) RHOP(I,J+1,K) = MAX(RHOMIN,RHOPJ+CONST*(RHO00-RHOPJ))
      IF (LC(-3)) RHOP(I,J,K-1) = MAX(RHOMIN,RHOMK+CONST*(RHO00-RHOMK))
      IF (LC( 3)) RHOP(I,J,K+1) = MAX(RHOMIN,RHOPK+CONST*(RHO00-RHOPK))
      RHOP(I,J,K) = RHOMIN
      ENDDO CHECK_LOOP
      ENDDO
      ENDDO
C
      END SUBROUTINE CHECK_DENSITY
C
C
      SUBROUTINE CHECK_MASS_FRACTION
C
      REAL(EB) SUM,ISUM,CONST,RHYMI,RHYPI,RHYMJ,RHYPJ,RHYMK,RHYPK,
     .         RHY0
      INTEGER IC,N
      LOGICAL LC(-3:3)
C
      IF (PREDICTOR) THEN
         RHOP=>RHOS
         YYP =>YYS
      ELSE
         RHOP=>RHO
         YYP =>YY
      ENDIF
C
      DO N=1,NSPEC
      DO K=1,KBAR
      DO J=1,JBAR
      CHECK_LOOP: DO I=1,IBAR
      IC = ICA(I,J,K)
      IF (SOLID(IC)) CYCLE CHECK_LOOP
      IF (MIXTURE_FRACTION .AND. N.EQ.IFUEL) THEN
         YYP(I,J,K,N) = MAX(YYMIN(N),YYP(I,J,K,N))
         YYP(I,J,K,N) = MIN(YYMAX(N),YYP(I,J,K,N))
         CYCLE CHECK_LOOP
         ENDIF
      IF (YYP(I,J,K,N).GE.YYMAX(N)) YYP(I,J,K,N) = YYMAX(N)
      IF (YYP(I,J,K,N).GE.YYMIN(N)) CYCLE CHECK_LOOP
      SUM = 0.
      ISUM = 0.
      LC = .FALSE.
      RHY0  = RHOP(I,J,K)  *YYP(I,J,K,N)
      RHYMI = RHOP(I-1,J,K)*YYP(I-1,J,K,N)
      RHYPI = RHOP(I+1,J,K)*YYP(I+1,J,K,N)
      RHYMJ = RHOP(I,J-1,K)*YYP(I,J-1,K,N)
      RHYPJ = RHOP(I,J+1,K)*YYP(I,J+1,K,N)
      RHYMK = RHOP(I,J,K-1)*YYP(I,J,K-1,N)
      RHYPK = RHOP(I,J,K+1)*YYP(I,J,K+1,N)
      IF (IWA(IC,-1).EQ.0 .AND. RHYMI.GT.0.) THEN
         SUM = SUM + RHYMI
         ISUM = ISUM + 1.
         LC(-1) = .TRUE.
         ENDIF
      IF (IWA(IC, 1).EQ.0 .AND. RHYPI.GT.0.) THEN
         SUM = SUM + RHYPI
         ISUM = ISUM + 1.
         LC( 1) = .TRUE.
         ENDIF
      IF (IWA(IC,-2).EQ.0 .AND. RHYMJ.GT.0.) THEN
         SUM = SUM + RHYMJ
         ISUM = ISUM + 1.
         LC(-2) = .TRUE.
         ENDIF
      IF (IWA(IC, 2).EQ.0 .AND. RHYPJ.GT.0.) THEN
         SUM = SUM + RHYPJ
         ISUM = ISUM + 1.
         LC( 2) = .TRUE.
         ENDIF
      IF (IWA(IC,-3).EQ.0 .AND. RHYMK.GT.0.) THEN
         SUM = SUM + RHYMK
         ISUM = ISUM + 1.
         LC(-3) = .TRUE.
         ENDIF
      IF (IWA(IC, 3).EQ.0 .AND. RHYPK.GT.0.) THEN
         SUM = SUM + RHYPK
         ISUM = ISUM + 1.
         LC( 3) = .TRUE.
         ENDIF
      IF (ISUM.EQ.0) THEN
         YYP(I,J,K,N) = YYMIN(N)
         CYCLE CHECK_LOOP
         ENDIF
      CONST = (YYMIN(N)-RHY0)/(SUM-ISUM*RHY0)
      IF(LC(-1)) YYP(I-1,J,K,N)=(RHYMI+CONST*(RHY0-RHYMI))/RHOP(I-1,J,K)
      IF(LC( 1)) YYP(I+1,J,K,N)=(RHYPI+CONST*(RHY0-RHYPI))/RHOP(I+1,J,K)
      IF(LC(-2)) YYP(I,J-1,K,N)=(RHYMJ+CONST*(RHY0-RHYMJ))/RHOP(I,J-1,K)
      IF(LC( 2)) YYP(I,J+1,K,N)=(RHYPJ+CONST*(RHY0-RHYPJ))/RHOP(I,J+1,K)
      IF(LC(-3)) YYP(I,J,K-1,N)=(RHYMK+CONST*(RHY0-RHYMK))/RHOP(I,J,K-1)
      IF(LC( 3)) YYP(I,J,K+1,N)=(RHYPK+CONST*(RHY0-RHYPK))/RHOP(I,J,K+1)
      YYP(I,J,K,N) = YYMIN(N)
      ENDDO CHECK_LOOP
      ENDDO
      ENDDO
      ENDDO
C
      END SUBROUTINE CHECK_MASS_FRACTION
C
C
      SUBROUTINE TEMPERATURE_BC(T,NM)
C
C Thermal boundary conditions
C
      REAL(EB) T,TSI,TMP_G,DTMP,TMP_OTHER,FLUX,CP_TERM,RHOWAL
      INTEGER  IOR,II,JJ,KK,IBC,IIG,JJG,KKG
      INTEGER, INTENT(IN) :: NM
C
      IF (PREDICTOR) THEN
         UU => U
         VV => V
         WW => W
         RHOP => RHOS
         ELSE
         UU => US
         VV => VS
         WW => WS
         RHOP => RHO
         ENDIF
C
C Loop through all boundary cells and apply heat transfer conditions
C
      HEAT_FLUX_LOOP: DO IW=1,NWC
C
      IF (IV(IW).EQ.0) CYCLE HEAT_FLUX_LOOP
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      IOR = IJKW(4,IW)
      IBC = IJKW(5,IW)
C
      HEAT_CONDUCTION: SELECT CASE(ABS(MHC(IBC)))
C
      CASE(0) HEAT_CONDUCTION  ! Adiabatic walls 
C
      TMP_F(IW) = TMP(IIG,JJG,KKG)
C
      IF (IV(IW).EQ.2) THEN
      SELECT CASE(IOR)
      CASE( 1) ; IF (UU(II,JJ,KK).GE.0.)   TMP_F(IW) = TMP_O(IW)
      CASE(-1) ; IF (UU(II-1,JJ,KK).LE.0.) TMP_F(IW) = TMP_O(IW)
      CASE( 2) ; IF (VV(II,JJ,KK).GE.0.)   TMP_F(IW) = TMP_O(IW)
      CASE(-2) ; IF (VV(II,JJ-1,KK).LE.0.) TMP_F(IW) = TMP_O(IW)
      CASE( 3) ; IF (WW(II,JJ,KK).GE.0.)   TMP_F(IW) = TMP_O(IW)
      CASE(-3) ; IF (WW(II,JJ,KK-1).LE.0.) TMP_F(IW) = TMP_O(IW)
      END SELECT
      TMP(II,JJ,KK) = TMP_F(IW)
      ENDIF
C
      TMP_W(IW) = TMP_F(IW)
C
      CASE(1) HEAT_CONDUCTION ! Prescribed temperature/heat flux
C
         TMP_G = TMP(IIG,JJG,KKG)
         TSI   = T-TW(IW)
C
         SELECT CASE(MHC(IBC))
C
         CASE(1)    ! Prescribed Temperature
            TMP_F(IW) = TMPA + FQ(TSI,ITAUQ(IW))*(TMP_P(IBC)-TMPA)
            DTMP = TMP_G - TMP_F(IW)
            IF (DNS) THEN
               HEAT_TRANS_COEF(IW) = 2.*KW(IW)*RDN(IW)
            ELSE
               SELECT CASE(ABS(IOR))
               CASE(1)
               V2 = 0.25*(VV(IIG,JJG,KKG)+VV(IIG,JJG-1,KKG))**2
               W2 = 0.25*(WW(IIG,JJG,KKG)+WW(IIG,JJG,KKG-1))**2
               VELCON = (V2+W2)**0.4
               H_NATURAL = HCV*ABS(DTMP)**ONTH
               CASE(2)
               U2 = 0.25*(UU(IIG,JJG,KKG)+UU(IIG-1,JJG,KKG))**2
               W2 = 0.25*(WW(IIG,JJG,KKG)+WW(IIG,JJG,KKG-1))**2
               VELCON = (U2+W2)**0.4
               H_NATURAL = HCV*ABS(DTMP)**ONTH
               CASE(3)
               U2 = 0.25*(UU(IIG,JJG,KKG)+UU(IIG-1,JJG,KKG))**2
               V2 = 0.25*(VV(IIG,JJG,KKG)+VV(IIG,JJG-1,KKG))**2
               VELCON = (U2+V2)**0.4
               H_NATURAL = HCH*ABS(DTMP)**ONTH
               END SELECT
               H_FORCED  = C_FORCED*VELCON*RHOP(IIG,JJG,KKG)**0.8
               ITMP = 0.1*TMP_G
               H_DNS = MU_SPEC(0,ITMP)*CP_GAMMA*RPR*2.*RDN(IW)
               HEAT_TRANS_COEF(IW) = MAX(H_DNS,H_FORCED,H_NATURAL)
            ENDIF
            QCONF(IW) = HEAT_TRANS_COEF(IW)*DTMP
C
         CASE(-1)   ! Prescribed Heat Flux
            TMP_F(IW) = TMPA + FQ(TSI,ITAUQ(IW))*(TMP_P(IBC)-TMPA)
            FLUX = FQ(TSI,ITAUQ(IW))*HEATFLUX(IBC)*AREA_ADJUST(IW)
            QCONF(IW) = -FLUX
         END SELECT
C
         RHOWAL    = 0.5*(RHOP(IIG,JJG,KKG)+RHO_W(IW))
         CP_TERM   = MAX(0._EB,-CP_GAMMA*UW(IW)*RHOWAL)
         TMP_W(IW) = ( (RDN(IW)*KW(IW)-0.5*CP_TERM)*TMP_G +
     .        CP_TERM*TMP_F(IW)-QCONF(IW) )/(0.5*CP_TERM+RDN(IW)*KW(IW))
         TMP_W(IW) = MAX(TMPMIN,TMP_W(IW))
C
      CASE(6) HEAT_CONDUCTION   ! Interpolated Boundary
C
      TMP_G = TMP(IIG,JJG,KKG)
      TMP_OTHER =
     .   OMESH(IJKW(9,IW))%TMP(IJKW(10,IW),IJKW(11,IW),IJKW(12,IW))
C
      IF (CELL_VOLUME_RATIO(IW).LT.0.5 .OR. 
     .    CELL_VOLUME_RATIO(IW).GT.2.0) THEN
C
      TMP_W(IW) = TMP_G
C
      SELECT CASE(IOR)
      CASE( 1) ; IF (UU(II,JJ,KK).GE.0.)   TMP_W(IW) = TMP_OTHER
      CASE(-1) ; IF (UU(II-1,JJ,KK).LE.0.) TMP_W(IW) = TMP_OTHER
      CASE( 2) ; IF (VV(II,JJ,KK).GE.0.)   TMP_W(IW) = TMP_OTHER
      CASE(-2) ; IF (VV(II,JJ-1,KK).LE.0.) TMP_W(IW) = TMP_OTHER
      CASE( 3) ; IF (WW(II,JJ,KK).GE.0.)   TMP_W(IW) = TMP_OTHER
      CASE(-3) ; IF (WW(II,JJ,KK-1).LE.0.) TMP_W(IW) = TMP_OTHER
      END SELECT
C
      ELSE
C
      TMP_W(IW) = TMP_OTHER
C
      ENDIF
C
      TMP_F(IW) = TMP_W(IW)
      QCONF(IW) = KW(IW)*(TMP_G-TMP_W(IW))*RDN(IW)
      TMP(II,JJ,KK) = TMP_W(IW)
C
      END SELECT HEAT_CONDUCTION
C
C Record wall temperature in the ghost cell
C
      IF (SOLID(ICA(II,JJ,KK))) 
     .   TMP(II,JJ,KK) = MAX(100._EB,MIN(4900._EB,TMP_W(IW)))
C
      ENDDO HEAT_FLUX_LOOP
C
C Thermally-thin and thermally-thick boundary conditions (MHC=2,3)
C
      IF (CORRECTOR) THEN
      WALL_COUNTER = WALL_COUNTER + 1
      IF (WALL_COUNTER.EQ.WALL_INCREMENT) THEN
      DT_BC   = T - WALLCLK
      WALLCLK = T
      CALL PYROLYSIS(T)
      WALL_COUNTER = 0
      ENDIF
      ENDIF
C
C Pyrolysis Models (MHC=4,5)
C
      IF (PREDICTOR .AND. PAPERMODEL) CALL SOLID_PHASE_REACTIONS(T,NM)
C
      END SUBROUTINE TEMPERATURE_BC
C
C
      SUBROUTINE SPECIES_BC(T)
C
      REAL(EB) T,YYWAL,YY_G,DENOM,YY_OTHER(1:20),
     .         RHO_G,UN,DD,EPSB,MFT,TSI
      INTEGER IBC,IIG,JJG,KKG,IOR,IC,IWB
C
      IF (PREDICTOR) THEN
         UU => U
         VV => V
         WW => W
         RHOP => RHOS
         YYP => YYS
         ELSE
         UU => US
         VV => VS
         WW => WS
         RHOP => RHO
         YYP => YY
         ENDIF 
C
C Loop through the wall cells, apply mass boundary conditions
C
      MASS_FLUX_LOOP: DO IW=1,NWC
C
      IF (IV(IW).EQ.0) CYCLE MASS_FLUX_LOOP
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IOR = IJKW(4,IW)
      IBC = IJKW(5,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
C
      MASS_TRANSFER: SELECT CASE(MMT(IBC))
C
      CASE(1) MASS_TRANSFER  ! Solid boundary, no mass flux
C
      YY_W(IW,1:NSPEC) = MAX(0.0_EB,YYP(IIG,JJG,KKG,1:NSPEC))
C
      IF (IV(IW).EQ.2) THEN
      SELECT CASE(IOR)
      CASE( 1) 
      IF (UU(II,JJ,KK).GT.0.)   YY_W(IW,1:NSPEC) = YY0(1:NSPEC)
      CASE(-1)
      IF (UU(II-1,JJ,KK).LT.0.) YY_W(IW,1:NSPEC) = YY0(1:NSPEC)
      CASE( 2) 
      IF (VV(II,JJ,KK).GT.0.)   YY_W(IW,1:NSPEC) = YY0(1:NSPEC)
      CASE(-2)
      IF (VV(II,JJ-1,KK).LT.0.) YY_W(IW,1:NSPEC) = YY0(1:NSPEC)
      CASE( 3)
      IF (WW(II,JJ,KK).GT.0.)   YY_W(IW,1:NSPEC) = YY0(1:NSPEC)
      CASE(-3)
      IF (WW(II,JJ,KK-1).LT.0.) YY_W(IW,1:NSPEC) = YY0(1:NSPEC)
      END SELECT
      YYP(II,JJ,KK,1:NSPEC)=YY_W(IW,1:NSPEC)
      ENDIF
C
      IF (LEAKY(IBC) .AND. UWS(IW).LT.0.) THEN
      RHO_G = RHOP(IIG,JJG,KKG)
      UN = -UWS(IW)
      IF (PREDICTOR) EPSB = -.5*UN**2*DT*RDN(IW)
      IF (CORRECTOR) EPSB =  .5*UN**2*DT*RDN(IW)
      SPECIES_LOOP_1: DO N=1,NSPEC
      DD    = RHODW(IW,N)*RDN(IW)
      YY_G  = YYP(IIG,JJG,KKG,N)
      DENOM = DD + (.5*UN+EPSB)*RHO_W(IW)
      YY_W(IW,N) = ( MASSFLUX(IW,N)*ADJUST_BURN_RATE(IBC,N) +
     .    YY_G*(DD + (EPSB-.5*UN)*RHO_G) ) / DENOM
      ENDDO SPECIES_LOOP_1
      ENDIF
C
      CASE(2) MASS_TRANSFER   ! Mass fraction prescribed
C
      IF (PREDICTOR) TSI = T + DT - TW(IW)
      IF (CORRECTOR) TSI = T      - TW(IW)
      DO N=1,NSPEC
      YYWAL = YY0(N) + FZ(TSI,ITAUMF(IW,N),N)*(MASSFRACS(IBC,N)-YY0(N))
      IF (DNS) YY_W(IW,N) = 2.*YYWAL - YYP(IIG,JJG,KKG,N)
      IF (LES) YY_W(IW,N) =    YYWAL 
      ENDDO
C
      CASE(3:5) MASS_TRANSFER    ! Mass flux prescribed
C
      IF (T.LT.TW(IW)) THEN
         YY_W(IW,1:NSPEC) = YYP(IIG,JJG,KKG,1:NSPEC)
         IF (SOLID(ICA(II,JJ,KK))) 
     .      YYP(II,JJ,KK,1:NSPEC)=YY_W(IW,1:NSPEC)
         IF (PREDICTOR) UWS(IW) = 0. 
         MASSFLUX(IW,1:NSPEC) = 0.
         CYCLE MASS_FLUX_LOOP
         ENDIF
C
C Compute fuel mass flux
C
      MFT = 0.
C
      MASSFLUX_LOOP: DO N=0,NSPEC
C
      IF (MASSFLUXS(IBC,N).GT.0.) THEN ! use user-prescribed flux
         TSI = T - TW(IW)
         MASSFLUX(IW,N) = FZ(TSI,ITAUMF(IW,N),N)*MASSFLUXS(IBC,N)
         IF (TMP_I(IBC).LT.2000. .AND. TMP_F(IW).LT.TMP_I(IBC) .AND.
     .       N.EQ.IFUEL .AND. T.GT.TW(IW)) THEN  ! Too cold, don't burn.
            TW(IW) = TW(IW) + DT
            MASSFLUX(IW,N) = 0.
            ENDIF
         IF (EW(IW).GT.0. .AND. N.EQ.IFUEL)
     .      MASSFLUX(IW,N) = MASSFLUX(IW,N)*EXP(-EW(IW))
      ENDIF
C
      MASSFLUX(IW,N) = MASSFLUX(IW,N)*AREA_ADJUST(IW)
      MFT = MFT + MASSFLUX(IW,N)*ADJUST_BURN_RATE(IBC,N)
C
      ENDDO MASSFLUX_LOOP
C
      IF (MASS_LOSS(IW).GT.SURFACE_DENSITY_S(IBC)) THEN
         MFT = 0.
         MASSFLUX(IW,:) = 0.
         ENDIF
C
C Subtract mass loss rate from fuel counters
C
      CONSUME_FUEL: IF (CORRECTOR .AND. IFUEL.GT.0) THEN
C
      IF (DENSITY_F(IBC).GT.0.) THEN
         IC = ICA(II,JJ,KK) 
         IF (SOLID(IC)) CELL_MASS(IC) = 
     .                  CELL_MASS(IC) - MASSFLUX(IW,IFUEL)*AW(IW)*DT
         ENDIF
C
      IF (SURFACE_DENSITY_S(IBC).GT.0. .AND. IBACK(IBC).EQ.3) THEN
         IC = ICA(II,JJ,KK)
         IF (SOLID(IC)) THEN
            IWB=IWA(IC,-IOR)
         ELSE
            IWB=IWA(IC,IOR)
         ENDIF
         IF (IWB.GT.0) MASS_LOSS(IWB) = 
     .                 MASS_LOSS(IWB) + MASSFLUX(IW,IFUEL)*DT
         ENDIF
C
      MASS_LOSS(IW) = MASS_LOSS(IW) + MASSFLUX(IW,IFUEL)*DT
C
      ENDIF CONSUME_FUEL
C
C Get normal velocity
C
      RHO_G = RHOP(IIG,JJG,KKG)
      UN    = 2.*MFT/(RHO_W(IW)+RHO_G)
      IF (PREDICTOR) UWS(IW) = -UN
C
      IF (PREDICTOR) EPSB = -.5*UN**2*DT*RDN(IW)
      IF (CORRECTOR) EPSB =  .5*UN**2*DT*RDN(IW)
C
      SPECIES_LOOP: DO N=1,NSPEC
      DD    = RHODW(IW,N)*RDN(IW)
      YY_G  = YYP(IIG,JJG,KKG,N)
      DENOM = DD + (.5*UN+EPSB)*RHO_W(IW)
      YY_W(IW,N) = ( MASSFLUX(IW,N)*ADJUST_BURN_RATE(IBC,N) + 
     .    YY_G*(DD + (EPSB-.5*UN)*RHO_G) ) / DENOM
      ENDDO SPECIES_LOOP
C
C
      CASE(6) MASS_TRANSFER   ! Interpolated Boundary
C
      YY_OTHER(1:NSPEC) = 
     . OMESH(IJKW(9,IW))%YY(IJKW(10,IW),IJKW(11,IW),IJKW(12,IW),1:NSPEC)
C
      IF (CELL_VOLUME_RATIO(IW).LT.0.5 .OR. 
     .    CELL_VOLUME_RATIO(IW).GT.2.0) THEN
C
      YY_W(IW,1:NSPEC) = YYP(IIG,JJG,KKG,1:NSPEC) 
C
      SELECT CASE(IOR)
      CASE( 1)
      IF (UU(II,JJ,KK).GT.0.)   YY_W(IW,1:NSPEC) = YY_OTHER(1:NSPEC)
      CASE(-1)
      IF (UU(II-1,JJ,KK).LT.0.) YY_W(IW,1:NSPEC) = YY_OTHER(1:NSPEC)
      CASE( 2)
      IF (VV(II,JJ,KK).GT.0.)   YY_W(IW,1:NSPEC) = YY_OTHER(1:NSPEC)
      CASE(-2)
      IF (VV(II,JJ-1,KK).LT.0.) YY_W(IW,1:NSPEC) = YY_OTHER(1:NSPEC)
      CASE( 3)
      IF (WW(II,JJ,KK).GT.0.)   YY_W(IW,1:NSPEC) = YY_OTHER(1:NSPEC)
      CASE(-3)
      IF (WW(II,JJ,KK-1).LT.0.) YY_W(IW,1:NSPEC) = YY_OTHER(1:NSPEC)
      END SELECT
C
      ELSE
C
      YY_W(IW,1:NSPEC) = YY_OTHER(1:NSPEC)
C
      ENDIF
C
      YYP(II,JJ,KK,1:NSPEC) = YY_W(IW,1:NSPEC)
C
      END SELECT MASS_TRANSFER
C
C Set ghost cell value of wall mass fraction
C
      IF (SOLID(ICA(II,JJ,KK))) 
     .    YYP(II,JJ,KK,1:NSPEC)=YY_W(IW,1:NSPEC)
C
      ENDDO MASS_FLUX_LOOP
C
      END SUBROUTINE SPECIES_BC
C
C
      SUBROUTINE DENSITY_BC
C
C Compute density at wall from wall temperatures and mass fractions 
C
      REAL(EB) P0X,RCON_W
C
      IF (PREDICTOR) THEN ; P0X = P0S  ; RHOP => RHOS ; ENDIF
      IF (CORRECTOR) THEN ; P0X = P0   ; RHOP => RHO  ; ENDIF
C
      WALL_CELL_LOOP: DO IW=1,NWC
C
      IF (IV(IW).EQ.0) CYCLE WALL_CELL_LOOP
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
C
      IF (.NOT.MIXTURE_FRACTION .AND. NSPEC.GT.0) THEN
      RSUM_W(IW) = RCON(0)
      DO N=1,NSPEC
      WFAC = RCON(N) - RCON(0)
      RSUM_W(IW) = RSUM_W(IW) + WFAC*YY_W(IW,N)
      ENDDO
      ENDIF
C
      IF (MIXTURE_FRACTION) THEN
      IYY = MAX(0,NINT(YY_W(IW,IFUEL)*10000.))
      IYY = MIN(10000,IYY)
      RSUM_W(IW) = RSUM_MF(IYY)
      RCON_W = RSUM_W(IW)
      SPEC_LOOP: DO N=1,NSPEC
      IF (N.EQ.IFUEL) CYCLE SPEC_LOOP
      RSUM_W(IW) = RSUM_W(IW) + YY_W(IW,N)*(RCON(N)-RCON_W)
      ENDDO SPEC_LOOP
      ENDIF

      IF (NSPEC.EQ.0) THEN
         RHO_W(IW) = P0X/(RCON(0)*TMP_W(IW))
      ELSE
         RHO_W(IW) = P0X/(RSUM_W(IW)*TMP_W(IW))
      ENDIF
C
      IF (SOLID(ICA(II,JJ,KK)) .OR. 
     .    IV(IW).EQ.2 .OR. IV(IW).EQ.4) THEN
         RHOP(II,JJ,KK) = RHO_W(IW)
         IF (NSPEC.GT.0) RSUM(II,JJ,KK) = RSUM_W(IW)
         ENDIF
C
      ENDDO WALL_CELL_LOOP
C
      END SUBROUTINE DENSITY_BC
C
C
      SUBROUTINE PYROLYSIS(T)
C
      REAL(EB) DTMP,QNETF,QNETB,QDXKF,QDXKB,RR,TMP_G,T,
     .         DDT(NWP_MAX),DTS,RFACF,RFACB,RFACF2,RFACB2,
     .         HVRG,PPCLAUS,PPSURF,YY_S,RSUM_S,TMP_G_B,QRADB,RHOWAL,
     .         CP_TERM,TSI, PFRONT_R,QPYR_RELAX,A_RATE
      REAL(EB) DT_F, LTMP, UTMP, FTMP, LFP, MDOT, FVELO, FDIS, DT_N,
     .         C_P_W, K_V, K_C, C_P_V_BAR, C_P_C_BAR, C_P_V_A,
     .         C_P_G_BAR, C_P_W_BAR, RHO_A, RHO_C, H_V_W, DXKF, DXKB,
     .         C_P_C_A,C_D_R,C_SS
      INTEGER IBC,IIG,JJG,KKG,IIB,JJB,KKB,IWB,IC,NN, FI,I
      TYPE(WALL_TYPE), POINTER :: WC
C
      C_P_W    = 4190.        ! Specific Heat of Water
      H_V_W    = 2259.*1000.  ! Heat of Vap of Water
C
C Loop through all the boundary cells, doing thin or thick heat tran.
C
      WALL_CELL_LOOP: DO IW=1,NWC
C
      IF (IV(IW).NE.1) CYCLE WALL_CELL_LOOP
C
      IBC = IJKW(5,IW)
C
      IF (WMPUA(IW).GT.0. .AND. T.GT.TW(IW))
     .    EW(IW) = EW(IW) + ECOEF(IBC)*WMPUA(IW)*DT_BC
C
      IF (MHC(IBC).LT.2 .OR. MHC(IBC).GT.3) CYCLE WALL_CELL_LOOP
C 
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      NN  = NWP(IBC)
      WC  => WALL(IW)
C
      IOR = IJKW(4,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
C
      TMP_G = TMP(IIG,JJG,KKG)
      DTMP = TMP_G - TMP_F(IW)
C
C Heat fluxes from radiation and convection at surface front
C
      IF (DNS) THEN
C
      HEAT_TRANS_COEF(IW) = 2.*KW(IW)*RDN(IW)
C
      ELSE
C
      SELECT CASE(ABS(IOR))
         CASE(1)
         V2 = 0.25*(V(IIG,JJG,KKG)+V(IIG,JJG-1,KKG))**2
         W2 = 0.25*(W(IIG,JJG,KKG)+W(IIG,JJG,KKG-1))**2
         VELCON = (V2+W2)**0.4
         H_NATURAL = HCV*ABS(DTMP)**ONTH
         CASE(2)
         U2 = 0.25*(U(IIG,JJG,KKG)+U(IIG-1,JJG,KKG))**2
         W2 = 0.25*(W(IIG,JJG,KKG)+W(IIG,JJG,KKG-1))**2
         VELCON = (U2+W2)**0.4
         H_NATURAL = HCV*ABS(DTMP)**ONTH
         CASE(3)
         U2 = 0.25*(U(IIG,JJG,KKG)+U(IIG-1,JJG,KKG))**2
         V2 = 0.25*(V(IIG,JJG,KKG)+V(IIG,JJG-1,KKG))**2
         VELCON = (U2+V2)**0.4
         H_NATURAL = HCH*ABS(DTMP)**ONTH
         END SELECT
      H_FORCED  = C_FORCED*VELCON*RHO(IIG,JJG,KKG)**0.8
      ITMP = 0.1*TMP_G
      H_DNS = MU_SPEC(0,ITMP)*CP_GAMMA*RPR*2.*RDN(IW)
      HEAT_TRANS_COEF(IW) = MAX(H_FORCED,H_NATURAL,H_DNS)
C
      ENDIF
C
      QCONF(IW) = HEAT_TRANS_COEF(IW)*DTMP
C
C Heat losses from convection and radiation out of back of surface
C
      SELECT CASE(IBACK(IBC))
C
      CASE(:1)  ! Non-insulated backing to an ambient void
C
         DTMP = TMP_VOID(IBC) - TMP_B(IW)
         IF (DNS) THEN
            QCONB(IW) = 2.*KW(IW)*DTMP*RDN(IW)
         ELSE
            IF (ABS(IOR).EQ.3) THEN
               QCONB(IW) = HCH*DTMP*ABS(DTMP)**ONTH
            ELSE
               QCONB(IW) = HCV*DTMP*ABS(DTMP)**ONTH
            ENDIF
         ENDIF
C
         QRADB = -SIGMA*E_WALL(IW)*(TMP_B(IW)**4-TMP_VOID(IBC)**4)
C
      CASE(2)  ! Insulated backing
C
         QCONB(IW) = 0.
         QRADB     = 0.
C
      CASE(3)  ! Exposed backing
C
         IC = ICA(II,JJ,KK)
C
         IF (SOLID(IC)) THEN
            IWB=IWA(IC,-IOR)
         ELSE
            IWB=IWA(IC,IOR)
         ENDIF
C
         IF (IWB.GT.0) THEN
            IIB = IJKW(6,IWB)
            JJB = IJKW(7,IWB)
            KKB = IJKW(8,IWB)
            TMP_G_B  = TMP(IIB,JJB,KKB)
            QRADB    = QRAD(IWB) 
         ELSE
            TMP_G_B  = TMPA
            QRADB    = -SIGMA*E_WALL(IW)*(TMP_B(IW)**4-TMPA4)
         ENDIF
C
         DTMP = TMP_G_B - TMP_B(IW)
         IF (DNS) THEN
            QCONB(IW) = 2.*KW(IW)*DTMP*RDN(IW)
         ELSE
            IF (ABS(IOR).EQ.3) THEN
               QCONB(IW) = HCH*DTMP*ABS(DTMP)**ONTH 
            ELSE
               QCONB(IW) = HCV*DTMP*ABS(DTMP)**ONTH
            ENDIF
         ENDIF
C
      END SELECT
C
C Sum up net heat fluxes, front and back
C
      QNETF = QRAD(IW) + QCONF(IW)
      QNETB = QRADB    + QCONB(IW)
C
C Take away energy flux due to water evaporation
C
      IF (NLP.GT.0) QNETF = QNETF - WCPUA(IW)
C
C Calculate temperature-dependent thermal properties for 1-D solver
C
      IF (MHC(IBC).EQ.2) THEN
C
      IF (CDR(IBC).GT.0.) THEN
         C_D_R = CDR(IBC)
         ELSE
         ITMP  = NINT(TMP_F(IW))
         C_SS  = FF2(ITMP,-NINT(C_PV(IBC)))*1000.
         C_D_R = C_SS*WALL_THICKNESS(IBC)*DENSITY_S(IBC)
         ENDIF
C
      ENDIF
C
      IF (MHC(IBC).EQ.3) THEN
C
      IF (XKS(IBC).GT.0.) THEN  
         K_S = XKS(IBC)
         ELSE
         DO I=0,NN+1
         ITMP   = NINT(WC%TMP_S(I))
         K_S(I) = FF2(ITMP,-NINT(XKS(IBC)))
         ENDDO
         ENDIF
      IF (C_PV(IBC).GT.0.) THEN
         C_P_V = C_PV(IBC)
         C_P_V_A = C_PV(IBC)
         ELSE
         DO I=1,NN
         ITMP   = NINT(WC%TMP_S(I))
         C_P_V(I) = FF2(ITMP,-NINT(C_PV(IBC)))*1000.
         C_P_V_A  = FF2(NINT(TMPA),-NINT(C_PV(IBC)))*1000.
         ENDDO
         ENDIF

      RHOCBAR(1:NN) = DENSITY_S(IBC)*C_P_V(1:NN)
C
      ENDIF
C
C     Compute the burning rate MASSFLUX(IW,IFUEL) and the energy
C     consumed by the pyrolysis process QPYR(IW)
C
      CHOOSE_PYRO: SELECT CASE(PYROLYSIS_MODEL(IBC))
C
      CASE('GAS_BURNER')
C
      IF (H_V(IBC).GT.0. .AND.
     .    MASS_LOSS(IW).LT.SURFACE_DENSITY_S(IBC)) THEN
         QPYR(IW) = MASSFLUX(IW,IFUEL)*H_V(IBC)
      ELSE
         QPYR(IW) = 0.
      ENDIF
C
      CASE('THERMOPLASTIC')
C
      IF (MASS_LOSS(IW).GT.SURFACE_DENSITY_S(IBC)) THEN
         IF (IBACK(IBC).LE.0) THEN
            IJKW(5,IW) = -IBACK(IBC)
            MASS_LOSS(IW) = 0.
            ENDIF
         MASSFLUX(IW,IFUEL) = 0.
         QPYR(IW) = 0.
      ELSE
         MASSFLUX(IW,IFUEL) = A_SOLID(IBC)*DENSITY_S(IBC)*
     .                        EXP(-E_SOLID(IBC)/(R0*TMP_F(IW)))
         QPYR(IW)           = MASSFLUX(IW,IFUEL)*H_V(IBC)
      ENDIF
C
      CASE('ABLATION')
C
      MASSFLUX(IW,IFUEL) = 0.
      IF (MASS_LOSS(IW).GT.SURFACE_DENSITY_S(IBC)) THEN
         IF (IBACK(IBC).LE.0) THEN
            IJKW(5,IW) = -IBACK(IBC)
            MASS_LOSS(IW) = 0.
            ENDIF
         QPYR(IW) = 0.
      ELSE
         A_RATE = A_SOLID(IBC)*DENSITY_S(IBC)*
     .              EXP(-E_SOLID(IBC)/(R0*TMP_F(IW)))
         MASS_LOSS(IW) = MASS_LOSS(IW) + A_RATE*DT_BC
         QPYR(IW) = H_V(IBC)*A_RATE
      ENDIF
C
      CASE('LIQUID')
C
      HVRG    = MW_FUEL*H_V(IBC)/R0
      PPCLAUS = MIN(1._EB,EXP(HVRG*(1./TMP_I(IBC)-1./TMP_F(IW))))
      YY_S    = 0.5*(YY_W(IW,IFUEL)+YY(IIG,JJG,KKG,IFUEL))
      RSUM_S  = 0.5*(RSUM_W(IW)+RSUM(IIG,JJG,KKG))
      PPSURF  = MIN(1._EB,YY_S*R0/(MW_FUEL*RSUM_S))
      IF (PPSURF.NE.PPCLAUS) THEN
         IF (QPYR(IW).EQ.0.AND.PPSURF.LT.PPCLAUS) THEN
            QPYR(IW) = 2.E-5*MW_FUEL*H_V(IBC)
            TW(IW)=T
           ELSE
            QPYR(IW) = QPYR(IW)*MIN(1.02_EB,MAX(0.98_EB,PPCLAUS/PPSURF))
            IF (QPYR(IW).LT.0.1) QPYR(IW)=0.
           ENDIF
         ENDIF
      IF (TMP_F(IW).GT.TMP_I(IBC)) QPYR(IW) = MAX(QPYR(IW),1.02*
     .         (QNETF - 2.*(TMP_I(IBC)-WC%TMP_S(1))/
     .         DXF(IBC)/K_S(1)))
C
      MASSFLUX(IW,IFUEL) = QPYR(IW)/H_V(IBC)
C
      IF (MASS_LOSS(IW).GT.SURFACE_DENSITY_S(IBC)) THEN
         MASSFLUX(IW,IFUEL) = 0.
         QPYR(IW) = 0.
         ENDIF
C
      CASE('CHAR')
C
      QPYR(IW)           = 0.
      MASSFLUX(IW,IFUEL) = 0.
      DRHOWDT(:)         = 0.
      DRHOMDT(:)         = 0.
C
C Temperature-dependent char specific heat and conductivity
C
      IF (C_PC(IBC).GT.0.) THEN
         C_P_C = C_PC(IBC)
         C_P_C_A = C_PC(IBC)
         ELSE
         DO I=1,NN
         ITMP   = NINT(WC%TMP_S(I))
         C_P_C(I) = FF2(ITMP,-NINT(C_PC(IBC)))*1000.
         C_P_C_A  = FF2(NINT(TMPA),-NINT(C_PC(IBC)))*1000.
         ENDDO
         ENDIF
C
      IF (XKS_C(IBC).GT.0.) THEN
         K_S_C = XKS_C(IBC)
         ELSE
         DO I=0,NN+1
         ITMP   = NINT(WC%TMP_S(I))
         K_S_C(I) = FF2(ITMP,-NINT(XKS_C(IBC)))
         ENDDO
         ENDIF
C
C Calculate material properties for charring solid
C
      DO I=1,NN
         LFP      = MIN(1._EB,MAX(0._EB,PFRONT_I(IW) - REAL(I,EB)))
         RHO_S(I) = DENSITY_S(IBC) - LFP*(DENSITY_S(IBC)-RHO_CHAR(IBC))
         RHO_C    = RHO_CHAR(IBC)*LFP
         RHO_A    = RHO_S(I) - RHO_C
         LFP      = MIN(1._EB,MAX(0._EB,MFRONT_I(IW) - REAL(I,EB)))
         RHO_M(I) = RHO_MOIS(IBC)*(1. - LFP)
         K_V      = K_S(I) + 0.396*RHO_M(I)/DENSITY_S(IBC)
         K_C      = K_S_C(I)
         K_S(I)   = K_V*RHO_A/DENSITY_S(IBC)+K_C*RHO_C/RHO_CHAR(IBC)
         RHOCBAR(I) = RHO_A*C_P_V(I) + RHO_C*C_P_C(I) + RHO_M(I)*C_P_W
         FTMP     = 0.5*(WC%TMP_S(I)+TMPA)
         C_P_V_BAR = 0.5*(C_P_V(I) + C_P_V_A)
         C_P_C_BAR = 0.5*(C_P_C(I) + C_P_C_A)
         C_P_G_BAR  = 66.8*FTMP**0.5 - 136.
         C_P_W_BAR  = C_P_W
         C_S(I)     = (DENSITY_S(IBC)*C_P_V_BAR 
     .                - RHO_CHAR(IBC)*C_P_C_BAR)
     .                /(DENSITY_S(IBC)-RHO_CHAR(IBC)) - C_P_G_BAR
         D_S(I)     = C_P_W_BAR - C_P_G_BAR
      ENDDO
C
      K_S(0)    = K_S(1)
      K_S(NN+1) = K_S(NN)
C     
C     Calculate moisture evaporation rate
C
      IF (IWATER.GT.0) MASSFLUX(IW,IWATER) = 0.
      FI     = MAX(1,FLOOR(MFRONT_I(IW)))
      DT_F   = DT_BC
      DO WHILE ((DT_F.GT.0.0).AND.(FI.LE.NN))
         LTMP          = 0.5*(WC%TMP_S(FI-1) + WC%TMP_S(FI))
         UTMP          = 0.5*(WC%TMP_S(FI)+WC%TMP_S(FI+1))
         LFP           = MFRONT_I(IW) - REAL(FI,EB)
         FTMP          = (LFP*UTMP + (1.-LFP)*LTMP)
         IF (FTMP.GT.TMP_E(IBC)) THEN
         MDOT          = MAX(0._EB,0.0001*(FTMP-TMP_E(IBC)))
         FVELO         = MDOT/(RHO_MOIS(IBC))
         FDIS          = FVELO*DT_F
         FDIS          = MIN(FDIS,(REAL(FI,EB)+1.-MFRONT_I(IW))
     .                           *(X_S(FI,IBC)-X_S(FI-1,IBC)))
         DT_N          = FDIS/(FVELO+1.0E-8)
         DRHOMDT(FI)   = -(DT_N/DT_BC)*MDOT*RDX_W(FI,IBC)
         DT_F          = DT_F-DT_N
         IF (SURF_GEOM(IBC).EQ.'CYLINDRICAL') THEN
            PFRONT_R     = WALL_THICKNESS(IBC) - X_S(FI-1,IBC)
            DRHOMDT(FI)  = DRHOMDT(FI)*PFRONT_R
            MFRONT_I(IW) = MFRONT_I(IW) + FDIS*RDX_W(FI,IBC)*PFRONT_R
            IF (IWATER.GT.0) MASSFLUX(IW,IWATER) = MASSFLUX(IW,IWATER) + 
     .        (DT_N/DT_BC)*MDOT*(PFRONT_R/WALL_THICKNESS(IBC))
         ELSE
            MFRONT_I(IW) = MFRONT_I(IW) + FDIS*RDX_W(FI,IBC)
            IF (IWATER.GT.0) MASSFLUX(IW,IWATER) = 
     .                       MASSFLUX(IW,IWATER) + (DT_N/DT_BC)*MDOT
         ENDIF
         FI            = FI + 1
         ELSE
         DT_F = 0.
         ENDIF
      ENDDO
C
C Calculate pyrolysis front movement
C
      FI      = MAX(1,FLOOR(PFRONT_I(IW)))
      DT_F    = DT_BC
      DO WHILE ((DT_F.GT.EPSILON(DT_F)).AND.(FI.LE.NN))
         LFP           = PFRONT_I(IW) - REAL(FI,EB)
         LTMP          = 0.5*(WC%TMP_S(FI-1) + WC%TMP_S(FI))
         UTMP          = 0.5*(WC%TMP_S(FI)+WC%TMP_S(FI+1))
         FTMP          = (LFP*UTMP + (1.-LFP)*LTMP)
         FVELO         = A_SOLID(IBC)*EXP(-E_SOLID(IBC)/(R0*FTMP))
         MDOT          = FVELO*(DENSITY_S(IBC)-RHO_CHAR(IBC))
         FDIS          = FVELO*DT_F
         FDIS          = MIN(FDIS,(REAL(FI,EB)+1.-PFRONT_I(IW))
     .                           *(X_S(FI,IBC)-X_S(FI-1,IBC)))
         DT_N          = FDIS/FVELO
         DRHOWDT(FI)   = -(DT_N/DT_BC)*MDOT*RDX_W(FI,IBC)
         DT_F          = DT_F-DT_N
         IF (SURF_GEOM(IBC).EQ.'CYLINDRICAL') THEN
            PFRONT_R     = WALL_THICKNESS(IBC) - X_S(FI-1,IBC)
            DRHOWDT(FI)  = DRHOWDT(FI)*PFRONT_R
            PFRONT_I(IW) = PFRONT_I(IW) + FDIS*RDX_W(FI,IBC)*PFRONT_R
            MASSFLUX(IW,IFUEL) = MASSFLUX(IW,IFUEL) + 
     .        (DT_N/DT_BC)*MDOT*(PFRONT_R/WALL_THICKNESS(IBC))
         ELSE
            PFRONT_I(IW) = PFRONT_I(IW) + FDIS*RDX_W(FI,IBC)
            MASSFLUX(IW,IFUEL) = MASSFLUX(IW,IFUEL) + (DT_N/DT_BC)*MDOT
         ENDIF
         FI            = FI + 1
      ENDDO
C
      IF (MASSFLUXS(IBC,IFUEL).GT.0.) THEN
         TSI = T-TW(IW)
         MASSFLUX(IW,IFUEL) = 
     .      FZ(TSI,ITAUMF(IW,IFUEL),IFUEL)*MASSFLUXS(IBC,IFUEL)
         ENDIF
C
      END SELECT CHOOSE_PYRO
C
C Limit the burning rate if specified by user
C
      IF (MASSFLUXS_MAX(IBC).LT.MASSFLUX(IW,IFUEL)) THEN
         MASSFLUX(IW,IFUEL) = MASSFLUXS_MAX(IBC)
         QPYR(IW) = MASSFLUXS_MAX(IBC)*H_V(IBC)
         ENDIF
C
C Subtract off energy used for pyrolysis from total
C
      QNETF = QNETF - QPYR(IW)
C
C Thermally-thick heat transfer calculation           
C
      IF (MHC(IBC).EQ.3) THEN   
C
      WC => WALL(IW)
C
C Set up arrays for thermally-thick solver
C
      DTS    = -.25*DT_BC
      DXKF   = DXF(IBC)/K_S(1)
      DXKB   = DXB(IBC)/K_S(NN)
C
      DO I=1,NN
      BBS(I)=DTS*(K_S(I)+K_S(I-1))*
     .            RDXN_W(I-1,IBC)*RDX_W(I,IBC)/RHOCBAR(I)
      AAS(I)=DTS*(K_S(I)+K_S(I+1))*
     .            RDXN_W(I,IBC)  *RDX_W(I,IBC)/RHOCBAR(I)
      ENDDO
C
      DDS(1:NN) = 1.      - AAS(1:NN) - BBS(1:NN)
      DDS(1)    = DDS(1)  + BBS(1)
      DDS(NN)   = DDS(NN) + AAS(NN)
C
      DO I=1,NN
      CCS(I) = WC%TMP_S(I) - AAS(I)*(WC%TMP_S(I+1)-WC%TMP_S(I)) +
     .                       BBS(I)*(WC%TMP_S(I)-WC%TMP_S(I-1)) 
      ENDDO
C
      SELECT CASE(PYROLYSIS_MODEL(IBC))
C
      CASE DEFAULT
      QPYR_RELAX = 0.
C
      CASE('THERMOPLASTIC')
      QPYR_RELAX = QPYR(IW)*E_SOLID(IBC)/(R0*TMP_F(IW)**2)
C
      CASE('ABLATION')
      QPYR_RELAX = QPYR(IW)*E_SOLID(IBC)/(R0*TMP_F(IW)**2)
C
      CASE('LIQUID')
      QPYR_RELAX = QPYR(IW)*MW_FUEL*H_V(IBC)/(R0*TMP_F(IW)**2)
C
      CASE('CHAR')
      QPYR_RELAX = 0.
      DO I=1,NN
      DTMP = WC%TMP_S(I) - TMPA
      CCS(I) = CCS(I)
     .       + DRHOWDT(I)*(H_V(IBC)-C_S(I)*DTMP)*DT_BC/RHOCBAR(I)
     .       + DRHOMDT(I)*(H_V_W   -D_S(I)*DTMP)*DT_BC/RHOCBAR(I)
      ENDDO
C
      END SELECT
C
      RFACF     = 0.5*DXKF*(4.*E_WALL(IW)*SIGMA*TMP_F(IW)**3 +
     .            QPYR_RELAX + HEAT_TRANS_COEF(IW) )
      RFACB     = DXKB*2.*E_WALL(IW)*SIGMA*TMP_B(IW)**3
      RFACF2    = (1.-RFACF)/(1.+RFACF)
      RFACB2    = (1.-RFACB)/(1.+RFACB)
C
      QDXKF     = (QNETF*DXKF + 2.*RFACF*TMP_F(IW))/(1.+RFACF)
      QDXKB     = (QNETB*DXKB + 2.*RFACB*TMP_B(IW))/(1.+RFACB)
C
      CCS(1)    = CCS(1)  - BBS(1) *QDXKF
      CCS(NN)   = CCS(NN) - AAS(NN)*QDXKB
      DDT       = DDS(:)
      DDT(1)    = DDT(1)  - BBS(1) *(1.-RFACF2)
      DDT(NN)   = DDT(NN) - AAS(NN)*(1.-RFACB2)
C
      DO I=2,NN
      RR = BBS(I)/DDT(I-1)
      DDT(I)    = DDT(I)  - RR*AAS(I-1)
      CCS(I)    = CCS(I)  - RR*CCS(I-1)
      ENDDO
      CCS(NN)   = CCS(NN)/DDT(NN)
      DO I=NN-1,1,-1
      CCS(I)    = (CCS(I) - AAS(I)*CCS(I+1))/DDT(I)
      ENDDO
C
      WC%TMP_S(1:NN) = MAX(TMPMIN,CCS(1:NN))
      WC%TMP_S(0)    = MAX(TMPMIN,WC%TMP_S(1) *RFACF2+QDXKF)
      WC%TMP_S(NN+1) = MAX(TMPMIN,WC%TMP_S(NN)*RFACB2+QDXKB)
      TMP_F(IW)      = WC%TMP_S(1)
ccc   TMP_F(IW)      = 0.5*(WC%TMP_S(0) +WC%TMP_S(1))
      TMP_F(IW)      = MIN(TMPMAX,MAX(TMPMIN,TMP_F(IW)))
      TMP_B(IW)      = 0.5*(WC%TMP_S(NN)+WC%TMP_S(NN+1))
      TMP_B(IW)      = MIN(TMPMAX,MAX(TMPMIN,TMP_B(IW)))
C
      ENDIF
C
C Thermally-thin heat transfer calculation
C
      IF (MHC(IBC).EQ.2) THEN
C
      TMP_F(IW)    = TMP_F(IW) + DT_BC*(QNETF+QNETB)/C_D_R
      TMP_F(IW)    = MIN(TMPMAX,MAX(TMPMIN,TMP_F(IW)))
cc    IF (H_V(IBC).EQ.0.)
cc   .TMP_F(IW)    = MIN(TMP_I(IBC),TMP_F(IW))
      TMP_B(IW)    = TMP_F(IW)
C
      ENDIF
C
C If the wall temperature exceeds the ignition temperature, burn it
C
      IF ( TMP_F(IW).GE.TMP_I(IBC) .AND. T.LT.TW(IW) ) TW(IW) = T
C
C Determine internal wall temperature for temperature BC
C
      RHOWAL  = 0.5*(RHO(IIG,JJG,KKG)+RHO_W(IW))
      CP_TERM = MAX(0._EB,-CP_GAMMA*UW(IW)*RHOWAL)
      TMP_W(IW) = ( (RDN(IW)*KW(IW)-0.5*CP_TERM)*TMP_G +
     .      CP_TERM*TMP_F(IW)-QCONF(IW) )/(0.5*CP_TERM+RDN(IW)*KW(IW))
      TMP_W(IW) = MAX(TMPMIN,TMP_W(IW))
C
      IF (SOLID(ICA(II,JJ,KK))) 
     .   TMP(II,JJ,KK) = MAX(100._EB,MIN(4900._EB,TMP_W(IW)))
C
      ENDDO WALL_CELL_LOOP
C
      END SUBROUTINE PYROLYSIS
C
C
      END SUBROUTINE DENSITY
C
C
      SUBROUTINE COMBUSTION(NM)
C
C Distribute particle heat to nearest four grid cells.
C
      REAL(EB) TMPD,YFU0,BFAC,ETRM,YO2MIN,YFUMIN,YO20,YO2N,YFUN,
     .         DYF,YO2,YFU,DYFDT,YO2S,YFUS,
     .         WFAC,RWFAC,HFAC,BWOW,DTT,VC,RHO_P,RHO_M,
     .         RHOD_N,YYNEW,Q_TOTAL,Q_EXCESS,QFAC,
     .         Z_P,Z_M,ZFLUX,DYODZ1,Z_MAX,FLUX_MAX,Q_STAR,Z_F_OLD,
     .         DNMIN,YO2_MAX,TMP_MIN,ZZZ_MIN,Y_O2_CORR,Y_EXTRA,
     .         MFLOW,DFLOW,TFLOW,RHO_A,Z_A,DA,NNN,
     .         TNOW_COMB,Q_MAX
      INTEGER NODETS,N,I,J,K,II,JJ,KK,IIG,JJG,KKG,IW,NN,IC1,IC2,
     .        IMIN,JMIN,KMIN,IBC
      REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOD,QT
      INTEGER, INTENT(IN) :: NM
C
      TNOW_COMB = SECOND()
C
C Select model for combustion
C
      SELECT CASE(COMBUSTION_MODEL)
C
      CASE(2)   ! Infinite-rate Reaction (LES, usually)
C
      IF (PREDICTOR) THEN
         UU=>U
         VV=>V
         WW=>W
      ELSE
         UU=>US
         VV=>VS
         WW=>WS
      ENDIF
C
      RHOD => WORK4
      QT   => WORK5
      QT   =  0.      
C
C Determine effective Z_F
C
      IF (AUTOMATIC_Z) THEN
C
      FLUX_MAX = 0.
      Z_F_OLD  = Z_F_EFF
C
      WALL_LOOP: DO IW=1,NWC
      IF (IV(IW).NE.1) CYCLE WALL_LOOP
      FLUX_MAX = MAX( FLUX_MAX , MASSFLUX(IW,IFUEL) )
      ENDDO WALL_LOOP
C
      Z_MAX = 0.
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      IF (.NOT.SOLID(ICA(I,J,K))) Z_MAX = MAX(YY(I,J,K,IFUEL),Z_MAX)
      ENDDO
      ENDDO
      ENDDO
C
      IF (Z_MAX.EQ.0. .AND. FLUX_MAX.EQ.0.) RETURN
C
      DNMIN = (DXMIN*DYMIN*DZMIN)**ONTH
      Q_STAR = DELTAH_FUEL*FLUX_MAX/
     .         (RHOA*CP_GAMMA*TMPA*(GRAV*DNMIN)**.5)
      IF (FLUX_MAX.EQ.0.) THEN
      Z_F_EFF = Z_F
      ELSE
      Z_MAX   = MAX(Z_MAX,Z_F)
      NNN     = Z_MAX/TANH(13.*FLUX_MAX)
      Z_F_EFF = Z_F*MIN(1._EB,(0.6/Z_CONSTANT)*NNN**0.8*Q_STAR**0.4)
      ENDIF
C
      IF (Z_F_OLD.GT.0.) Z_F_EFF = 0.1*Z_F_EFF + 0.9*Z_F_OLD
C
      ENDIF
C
C Set a few parameters for the HRR calculation
C
      DYODZ1 = -(MW_O2*NU_O2/(MW_FUEL*NU_FUEL))*Y_F_INLET/(1.-Z_F_EFF)
C
C Compute rho*D for either DNS or LES
C
      IF (DNS) THEN
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      ITMP = 0.1*TMP(I,J,K)
      IYY  = MAX(0,NINT(YY(I,J,K,IFUEL)*100.))
      IYY  = MIN(100,IYY)
      RHOD(I,J,K) = RHO(I,J,K)*D_SPEC(IYY,ITMP)
      ENDDO
      ENDDO
      ENDDO
C
      ELSE
C
      RHOD = MU*RSC
C
      ENDIF
C
C Compute Q from fuel flux across x-face boundaries
C
      MFLOW = 0.
      DFLOW = 0.
      TFLOW = 0.
C
      DO K=1,KBAR
      DO J=1,JBAR
      LOOP1: DO I=1,IBM1
      IC1 = ICA(I,J,K) ; IC2 = ICA(I+1,J,K)
      IF (SOLID(IC1) .OR. SOLID(IC2)) CYCLE LOOP1
      IW = IWA(IC1,1) ; IF (IV(IW).EQ.1) CYCLE LOOP1
      Z_P = YY(I+1,J,K,IFUEL)
      Z_M = YY(I,J,K,IFUEL)
      IF (Z_P.GE.Z_F_EFF .EQV. Z_M.GE.Z_F_EFF) CYCLE LOOP1
      RHO_P  = RHO(I+1,J,K)
      RHO_M  = RHO(I,J,K)
      RHO_A  = 0.5*(RHO_P+RHO_M)
      Z_A    = 0.5*(RHO_P*Z_P+RHO_M*Z_M)/RHO_A
      RHOD_N = 0.5*(RHOD(I+1,J,K)+RHOD(I,J,K))
      ZFLUX  = RDXN(I)*( RHOD_N*(Z_P-Z_M) )
     .                  -UU(I,J,K)*RHO_A*(Z_A-Z_F_EFF) 
      IF (Z_P.LT.Z_M) ZFLUX = -ZFLUX
      IF (RHO_P.GT.RHO_M) THEN
         QT(I+1,J,K) = QT(I+1,J,K) + ZFLUX*RDX(I+1)
      ELSE
         QT(I,J,K)   = QT(I,J,K)   + ZFLUX*RDX(I)
      ENDIF
c     DA = DY(J)*DZ(K)
c     IF (Z_P.LT.Z_M) THEN
c        MFLOW = MFLOW + UU(I,J,K)*RHO_A*DA
c        DFLOW = DFLOW + ZFLUX*DA
c        TFLOW = TFLOW + (ZFLUX + UU(I,J,K)*RHO_A*Z_F_EFF)*DA
c     ELSE
c        MFLOW = MFLOW - UU(I,J,K)*RHO_A*DA
c        DFLOW = DFLOW + ZFLUX*DA
c        TFLOW = TFLOW + (ZFLUX - UU(I,J,K)*RHO_A*Z_F_EFF)*DA
c     ENDIF
      ENDDO LOOP1
      ENDDO
      ENDDO
C
C Compute Q from fuel flux across y-face boundaries
C
      DO K=1,KBAR
      DO J=1,JBM1
      LOOP2: DO I=1,IBAR
      IC1 = ICA(I,J,K) ; IC2 = ICA(I,J+1,K)
      IF (SOLID(IC1) .OR. SOLID(IC2)) CYCLE LOOP2
      IW = IWA(IC1,2) ; IF (IV(IW).EQ.1) CYCLE LOOP2
      Z_P = YY(I,J+1,K,IFUEL)
      Z_M = YY(I,J,K,IFUEL)
      IF (Z_P.GE.Z_F_EFF .EQV. Z_M.GE.Z_F_EFF) CYCLE LOOP2
      RHO_P  = RHO(I,J+1,K)
      RHO_M  = RHO(I,J,K)
      RHO_A  = 0.5*(RHO_P+RHO_M)
      Z_A    = 0.5*(RHO_P*Z_P+RHO_M*Z_M)/RHO_A
      RHOD_N = 0.5*(RHOD(I,J+1,K)+RHOD(I,J,K))
      ZFLUX  = RDYN(J)*( RHOD_N*(Z_P-Z_M) )
     .                  -VV(I,J,K)*RHO_A*(Z_A-Z_F_EFF) 
      IF (Z_P.LT.Z_M) ZFLUX = -ZFLUX
      IF (RHO_P.GT.RHO_M) THEN
         QT(I,J+1,K) = QT(I,J+1,K) + ZFLUX*RDY(J+1)
      ELSE
         QT(I,J,K)   = QT(I,J,K)   + ZFLUX*RDY(J)
      ENDIF
c     DA = DX(I)*DZ(K)
c     IF (Z_P.LT.Z_M) THEN
c        MFLOW = MFLOW + VV(I,J,K)*RHO_A*DA
c        DFLOW = DFLOW + ZFLUX*DA
c        TFLOW = TFLOW + (ZFLUX + VV(I,J,K)*RHO_A*Z_F_EFF)*DA
c     ELSE
c        MFLOW = MFLOW - VV(I,J,K)*RHO_A*DA
c        DFLOW = DFLOW + ZFLUX*DA
c        TFLOW = TFLOW + (ZFLUX - VV(I,J,K)*RHO_A*Z_F_EFF)*DA
c     ENDIF
      ENDDO LOOP2
      ENDDO
      ENDDO
C
C Compute Q from fuel flux across z-face boundaries
C
      DO K=1,KBM1
      DO J=1,JBAR
      LOOP3: DO I=1,IBAR
      IC1 = ICA(I,J,K) ; IC2 = ICA(I,J,K+1)
      IF (SOLID(IC1) .OR. SOLID(IC2)) CYCLE LOOP3
      IW = IWA(IC1,3) ; IF (IV(IW).EQ.1) CYCLE LOOP3
      Z_P = YY(I,J,K+1,IFUEL)
      Z_M = YY(I,J,K,IFUEL)
      IF (Z_P.GE.Z_F_EFF .EQV. Z_M.GE.Z_F_EFF) CYCLE LOOP3
      RHO_P  = RHO(I,J,K+1)
      RHO_M  = RHO(I,J,K)
      RHO_A  = 0.5*(RHO_P+RHO_M)
      Z_A    = 0.5*(RHO_P*Z_P+RHO_M*Z_M)/RHO_A
      RHOD_N = 0.5*(RHOD(I,J,K+1)+RHOD(I,J,K))
      ZFLUX  = RDZN(K)*( RHOD_N*(Z_P-Z_M) )
     .                  -WW(I,J,K)*RHO_A*(Z_A-Z_F_EFF) 
      IF (Z_P.LT.Z_M) ZFLUX = -ZFLUX
      IF (RHO_P.GT.RHO_M) THEN
         QT(I,J,K+1) = QT(I,J,K+1) + ZFLUX*RDZ(K+1)
      ELSE
         QT(I,J,K)   = QT(I,J,K)   + ZFLUX*RDZ(K)
      ENDIF
c     DA = DX(I)*DY(J)
c     IF (Z_P.LT.Z_M) THEN
c        MFLOW = MFLOW + WW(I,J,K)*RHO_A*DA
c        DFLOW = DFLOW + ZFLUX*DA
c        TFLOW = TFLOW + (ZFLUX + WW(I,J,K)*RHO_A*Z_F_EFF)*DA
c     ELSE
c        MFLOW = MFLOW - WW(I,J,K)*RHO_A*DA
c        DFLOW = DFLOW + ZFLUX*DA
c        TFLOW = TFLOW + (ZFLUX - WW(I,J,K)*RHO_A*Z_F_EFF)*DA
c     ENDIF
      ENDDO LOOP3
      ENDDO
      ENDDO
C
C Compute Q near solid wall boundaries
C
      WLOOP: DO IW=1,NWC
      IF (IV(IW).NE.1) CYCLE WLOOP
      IF (MASSFLUX(IW,IFUEL).LE.0.) CYCLE WLOOP
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      Z_P = YY(IIG,JJG,KKG,IFUEL)
      IF (Z_P.GE.Z_F_EFF) CYCLE WLOOP
      IBC = IJKW(5,IW)
      QT(IIG,JJG,KKG) = QT(IIG,JJG,KKG) + 
     .  (1.-Z_F_EFF)*MASSFLUX(IW,IFUEL)*
     .    ADJUST_BURN_RATE(IBC,IFUEL)*RDN(IW)
c     MFLOW = MFLOW + MASSFLUX(IW,IFUEL)*AW(IW)
c     TFLOW = TFLOW + MASSFLUX(IW,IFUEL)*AW(IW)
      ENDDO WLOOP
C
c     write(6000,*) T,MFLOW,DFLOW,TFLOW
C
C Compute heat release rate from fuel flux
C
      QFAC = -DYODZ1*EPUMO2
      QT = QT*QFAC
C
C Adjust heat release rate to account for total fuel flux
C
      Q_EXCESS = 0.
      Q_TOTAL  = 0.
C
      DO K=1,KBAR
      DO J=1,JBAR
      ILOOP: DO I=1,IBAR
      IF (QT(I,J,K).EQ.0.) CYCLE ILOOP
      VC = RC(I)*DX(I)*DY(J)*DZ(K)
      IF (QT(I,J,K).GT.Q_UPPER) THEN
         Q_EXCESS = Q_EXCESS + (QT(I,J,K)-Q_UPPER)*VC
         QT(I,J,K) = Q_UPPER
         ENDIF
      IF (QT(I,J,K).LT.0.) THEN
         Q_EXCESS = Q_EXCESS + (QT(I,J,K)-0.5*Q_UPPER)*VC
         QT(I,J,K) = 0.5*Q_UPPER
         ENDIF
      Q_TOTAL = Q_TOTAL + QT(I,J,K)*VC
      ENDDO ILOOP
      ENDDO
      ENDDO
C
      IF (Q_TOTAL.GT.0.) THEN
         QFAC = MAX(0._EB,(Q_TOTAL+Q_EXCESS)/Q_TOTAL)
      ELSE
         QFAC = 1.
      ENDIF
C
      QT = QT*QFAC
C
      Q = 0.8*Q + 0.2*QT
ccc   Q = QT
C
C     Gas phase suppression routine
C
      IF_SUPPRESSION: IF (SUPPRESSION) THEN
C
      DO K=1,KBAR
      DO J=1,JBAR
      LOOP4: DO I=1,IBAR
      IF (Q(I,J,K).LE.0.) CYCLE LOOP4
C
      TMP_MIN = TMPMAX
      ZZZ_MIN = 1.
      Q_MAX   = 0.
C
      DO KK=MAX(1,K-2),MIN(KBAR,K+2)
      DO JJ=MAX(1,J-2),MIN(JBAR,J+2)
      DO II=MAX(1,I-2),MIN(IBAR,I+2)
      Q_MAX = MAX(Q_MAX,Q(II,JJ,KK))
      IF (YY(II,JJ,KK,IFUEL).LT.ZZZ_MIN .AND. 
     .    .NOT.SOLID(ICA(II,JJ,KK)) ) THEN
         ZZZ_MIN = YY(II,JJ,KK,IFUEL)
         TMP_MIN = TMP(II,JJ,KK)
         IMIN = II
         JMIN = JJ
         KMIN = KK
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      IYY = MIN(10000,MAX(0,NINT(ZZZ_MIN*10000.)))
      Y_EXTRA = 0.
      DO NN=1,NSPEC
      IF (NN.NE.IFUEL) Y_EXTRA = Y_EXTRA + YY(IMIN,JMIN,KMIN,NN)
      ENDDO
C
      YO2_MAX   = (1.-Y_EXTRA)*Y_STATE(IYY,2)
      Y_O2_CORR = Y_O2_LL*(TMP_CRIT-TMP_MIN)/(TMP_CRIT-TMPA)
C
      IF (YO2_MAX.LT.Y_O2_CORR) Q(I,J,K) = 0.
C
      ENDDO LOOP4
      ENDDO
      ENDDO
C
      ENDIF IF_SUPPRESSION
C
C
C
      CASE(3)   ! Finite-rate Global Arrhenius Reaction (DNS, usually)
C
      Q      = 0.
C
      HFAC   = DELTAH_FUEL/DT
      WFAC   = MWN(IOXYGEN)*NUN(IOXYGEN)/(MWN(IFUEL)*NUN(IFUEL))
      RWFAC  = 1./WFAC
      BWOW   = (0.001)**(XNF+XNO-1.)*BOF/
     .         (MWN(IOXYGEN)**XNO*MWN(IFUEL)**(XNF-1.))
      NODETS = 20
C     NODETS = 100
      DTT    = DT/REAL(NODETS,EB)
C
C Cycle through the cells to do the reaction rate
C
      DO K=1,KBAR
      DO J=1,JBAR
      ILOOP3: DO I=1,IBAR
C
      IF (SOLID(ICA(I,J,K))) CYCLE ILOOP3
      YO20  = YY(I,J,K,IOXYGEN)
      YFU0  = YY(I,J,K,IFUEL)
      YO2MIN = MAX(1.E-10_EB,YO20+Q_UPPER*DT*MPUE(IOXYGEN)/RHO(I,J,K))
      YFUMIN = MAX(1.E-10_EB,YFU0+Q_UPPER*DT*MPUE(IFUEL)  /RHO(I,J,K))
      IF (YO20.LE.YO2MIN .OR. YFU0.LE.YFUMIN) CYCLE ILOOP3
      TMPD = TMP(I,J,K)
      IF (TMPD.LT.TMP_LOWER) CYCLE ILOOP3
      ETRM = EXP(-E_GAS/(R0*TMPD))
      BFAC = -BWOW*ETRM*RHO(I,J,K)**(XNF+XNO-1.)
      YFU  = YFU0
      YO2  = YO20
C
      ODELOOP: DO II=1,NODETS
      IF (YFU.LE.YFUMIN .OR. YO2.LE.YO2MIN) EXIT ODELOOP
      DYFDT= BFAC*YFU**XNF*YO2**XNO
      YFUS = YFU + DTT*DYFDT
      YO2S = YO2 + DTT*WFAC*DYFDT
      IF (YO2S.LE.YO2MIN) THEN
         YFU = MAX(0.0_EB,YFU-YO2*RWFAC)
         EXIT ODELOOP
         ENDIF
      IF (YFUS.LE.YFUMIN) THEN
         YFU = YFUMIN
         EXIT ODELOOP
         ENDIF
      DYFDT= BFAC*YFUS**XNF*YO2S**XNO
      YFUN = .5*(YFU+YFUS+DTT*DYFDT)
      YO2N = .5*(YO2+YO2S+DTT*WFAC*DYFDT)
      IF (YO2N.LE.YO2MIN) THEN
         YFU = MAX(0.0_EB,0.5*((YFU+YFUS)-(YO2+YO2S)*RWFAC))
         EXIT ODELOOP
         ENDIF
      IF (YFUN.LE.YFUMIN) THEN
         YFU = YFUMIN
         EXIT ODELOOP
         ENDIF
      YFU  = YFUN
      YO2  = YO2N
      ENDDO ODELOOP
C
      DYF   = YFU - YFU0
      DYF   = MAX(-YFU0,DYF)
      Q(I,J,K)  = MIN(Q_UPPER,-DYF*RHO(I,J,K)*HFAC)
C
      ENDDO ILOOP3
      ENDDO
      ENDDO
C
C Change species mass fractions
C
      DO N=1,NSPEC
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      YYNEW = YY(I,J,K,N) + DT*MPUE(N)*Q(I,J,K)/RHO(I,J,K)
      YY(I,J,K,N) = MAX(YYMIN(N),YYNEW)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      END SELECT
C
      TUSED(10,NM)=TUSED(10,NM)+SECOND()-TNOW_COMB
C
      END SUBROUTINE COMBUSTION
C
C
      END MODULE MASS
