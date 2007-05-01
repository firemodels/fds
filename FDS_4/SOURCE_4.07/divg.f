      MODULE DIVG              
C
      USE PREC
      USE CONS
      USE PACKER
C
      IMPLICIT NONE
      PRIVATE
      INTEGER I,J,K
      PUBLIC DIVERGENCE,CHECK_DIVERGENCE
C
      CONTAINS
C
C
      SUBROUTINE DIVERGENCE(T,NM)
C
C Compute the divergence of the flow D, and the time derivative of D
C
      INTEGER, INTENT(IN) :: NM
      REAL(EB), POINTER, DIMENSION(:,:,:) :: KDTDX,KDTDY,KDTDZ,DP,KP,
     .         RDDYDX,RDDYDY,RDDYDZ,RHOD,RHOP,TRDDYDX,TRDDYDY,TRDDYDZ,
     .         D_OLD,RTRM
      REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YYP
      REAL(EB) DELKDELT,VC,RP0,DTDX,DTDY,DTDZ,K_SUM,CP_SUM,
     .         USUM,DSUM,PSUM,GFAC,XIF,YIF,ZIF,U_LEAK,
     .         HDIFF,DYDX,DYDY,DYDZ,T,RDT,RDD,AREA_V,AREA_L,TSI,LEAK
      SAVE AREA_V,AREA_L
      INTEGER IW,N,IOR,II,JJ,KK,IIG,JJG,KKG,ITMP,IBC,
     .        IYY,NM_IN,II_IN,JJ_IN,KK_IN,IC
C
      TNOW=SECOND()
C
      CALL UNPACK_VAR(NM)
C
      RDT = 1./DT
C
      SELECT CASE(PREDICTOR)
      CASE(.TRUE.)  ; DP => DS   ; RP0=1./P0
      CASE(.FALSE.) ; DP => DDDT ; RP0=1./P0S
      END SELECT
C
      IF (NSPEC.GT.0) THEN
      RDDYDX  => WORK1
      RDDYDY  => WORK2
      RDDYDZ  => WORK3
      SELECT CASE(PREDICTOR)
      CASE(.TRUE.)  ; YYP => YYS ; RHOP => RHOS
      CASE(.FALSE.) ; YYP => YY  ; RHOP => RHO
      END SELECT
      ENDIF
C
C Zero out divergence to start
C
      DP  = 0.
      IF (NSPEC.GT.0) FYY = 0.
C
C Add species diffusion terms to divergence expression and
C compute diffusion term for species equations
C
      SPECIES_LOOP: DO N=1,NSPEC
C
C Compute rho*D
C
      RHOD => WORK4
C
      LES_VS_DNS: IF (DNS) THEN
C
      DNS_RHOD: IF (.NOT.MIXTURE_FRACTION) THEN
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      ITMP = 0.1*TMP(I,J,K)
      RHOD(I,J,K) = RHOP(I,J,K)*D_SPEC(N,ITMP)
      ENDDO 
      ENDDO
      ENDDO
C
      ELSE DNS_RHOD
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      ITMP = 0.1*TMP(I,J,K)
      IYY  = MIN(100,MAX(0,NINT(YYP(I,J,K,IFUEL)*100.)))
      RHOD(I,J,K) = RHOP(I,J,K)*D_SPEC(IYY,ITMP)
      ENDDO
      ENDDO
      ENDDO
C
      ENDIF DNS_RHOD
C
      ELSE LES_VS_DNS
C
      RHOD = MU*RSC
C
      ENDIF LES_VS_DNS
C
C Compute rho*D*dY/dx
C
      DO K=0,KBAR
      DO J=0,JBAR
      DO I=0,IBAR
      DYDX = (YYP(I+1,J,K,N)-YYP(I,J,K,N))*RDXN(I)
      RDDYDX(I,J,K) = .5*(RHOD(I+1,J,K)+RHOD(I,J,K))*DYDX
      DYDY = (YYP(I,J+1,K,N)-YYP(I,J,K,N))*RDYN(J)
      RDDYDY(I,J,K) = .5*(RHOD(I,J+1,K)+RHOD(I,J,K))*DYDY
      DYDZ = (YYP(I,J,K+1,N)-YYP(I,J,K,N))*RDZN(K)
      RDDYDZ(I,J,K) = .5*(RHOD(I,J,K+1)+RHOD(I,J,K))*DYDZ
      ENDDO
      ENDDO
      ENDDO
C
C Correct rho*D*dY/dx at boundaries and store rho*D at boundaries
C
      WLOOP: DO IW=1,NWC
      IF (IV(IW).EQ.0) CYCLE WLOOP
      IIG = IJKW(6,IW) 
      JJG = IJKW(7,IW) 
      KKG = IJKW(8,IW) 
      RHODW(IW,N) = RHOD(IIG,JJG,KKG)
      RDD = RHODW(IW,N)*(YYP(IIG,JJG,KKG,N)-YY_W(IW,N))*RDN(IW)
      FYY(IIG,JJG,KKG,N) = FYY(IIG,JJG,KKG,N) + RDD*RDN(IW)
      IOR = IJKW(4,IW)
      SELECT CASE(IOR) 
      CASE( 1) ; RDDYDX(IIG-1,JJG,KKG) = 0.
      CASE(-1) ; RDDYDX(IIG,JJG,KKG)   = 0.
      CASE( 2) ; RDDYDY(IIG,JJG-1,KKG) = 0.
      CASE(-2) ; RDDYDY(IIG,JJG,KKG)   = 0.
      CASE( 3) ; RDDYDZ(IIG,JJG,KKG-1) = 0.
      CASE(-3) ; RDDYDZ(IIG,JJG,KKG)   = 0.
      END SELECT
      ENDDO WLOOP
C
C Compute species transport terms for divergence expression
C
      SPECIES_DIFFUSION: IF (.NOT.MIXTURE_FRACTION) THEN
C
      TRDDYDX => WORK4
      TRDDYDY => WORK5
      TRDDYDZ => WORK6
C
      DO K=0,KBAR
      DO J=0,JBAR
      DO I=0,IBAR
      ITMP = .05*(TMP(I+1,J,K)+TMP(I,J,K))
      HDIFF = HH(N,ITMP)-HH(0,ITMP)
      TRDDYDX(I,J,K) = HDIFF*RDDYDX(I,J,K)
      ITMP = .05*(TMP(I,J+1,K)+TMP(I,J,K))
      HDIFF = HH(N,ITMP)-HH(0,ITMP)
      TRDDYDY(I,J,K) = HDIFF*RDDYDY(I,J,K)
      ITMP = .05*(TMP(I,J,K+1)+TMP(I,J,K))
      HDIFF = HH(N,ITMP)-HH(0,ITMP)
      TRDDYDZ(I,J,K) = HDIFF*RDDYDZ(I,J,K)
      ENDDO
      ENDDO
      ENDDO
C
      WLOOP2: DO IW=1,NWC
      IF (IV(IW).EQ.0) CYCLE WLOOP2
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      IOR  = IJKW(4,IW)
      ITMP = 0.1*TMP_F(IW)
      HDIFF = HH(N,ITMP)-HH(0,ITMP)
      RDD = RHODW(IW,N)*(YYP(IIG,JJG,KKG,N)-YY_W(IW,N))*RDN(IW)
      SELECT CASE(IOR)
      CASE( 1) ; TRDDYDX(IIG-1,JJG,KKG) = 0.
      CASE(-1) ; TRDDYDX(IIG,JJG,KKG)   = 0.
      CASE( 2) ; TRDDYDY(IIG,JJG-1,KKG) = 0.
      CASE(-2) ; TRDDYDY(IIG,JJG,KKG)   = 0.
      CASE( 3) ; TRDDYDZ(IIG,JJG,KKG-1) = 0.
      CASE(-3) ; TRDDYDZ(IIG,JJG,KKG)   = 0.
      END SELECT
      DP(IIG,JJG,KKG) = DP(IIG,JJG,KKG) - RDN(IW)*HDIFF*RDD
      ENDDO WLOOP2
C
      CYLINDER: SELECT CASE(CYLINDRICAL)
C
      CASE(.FALSE.) CYLINDER  ! 3D or 2D Cartesian Coords
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      DP(I,J,K) = DP(I,J,K) +
     .            (TRDDYDX(I,J,K)-TRDDYDX(I-1,J,K))*RDX(I) +
     .            (TRDDYDY(I,J,K)-TRDDYDY(I,J-1,K))*RDY(J) +
     .            (TRDDYDZ(I,J,K)-TRDDYDZ(I,J,K-1))*RDZ(K)
      ENDDO
      ENDDO
      ENDDO
C
      CASE(.TRUE.) CYLINDER  ! 2D Cylindrical Coords
C
      J = 1
      DO K=1,KBAR
      DO I=1,IBAR
      DP(I,J,K) = DP(I,J,K) +
     . (R(I)*TRDDYDX(I,J,K)-R(I-1)*TRDDYDX(I-1,J,K))*RDX(I)*RRN(I) +
     .      (TRDDYDZ(I,J,K)-       TRDDYDZ(I,J,K-1))*RDZ(K)
      ENDDO
      ENDDO
C
      END SELECT CYLINDER
C
      ENDIF SPECIES_DIFFUSION
C
C Compute del dot rho D del Y_i or del dot rho D del Z
C
      CYLINDER2: SELECT CASE(CYLINDRICAL)
C
      CASE(.FALSE.) CYLINDER2  ! 3D or 2D Cartesian Coords
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      FYY(I,J,K,N) = FYY(I,J,K,N) +
     .               (RDDYDX(I-1,J,K)-RDDYDX(I,J,K))*RDX(I) +
     .               (RDDYDY(I,J-1,K)-RDDYDY(I,J,K))*RDY(J) +
     .               (RDDYDZ(I,J,K-1)-RDDYDZ(I,J,K))*RDZ(K)
      ENDDO
      ENDDO 
      ENDDO
C
      CASE(.TRUE.) CYLINDER2  ! 2D Cylindrical Coords
C
      J=1
      DO K=1,KBAR
      DO I=1,IBAR
      FYY(I,J,K,N) = FYY(I,J,K,N) +
     . (R(I-1)*RDDYDX(I-1,J,K)-R(I)*RDDYDX(I,J,K))*RDX(I)*RRN(I) +
     .        (RDDYDZ(I,J,K-1)-     RDDYDZ(I,J,K))*RDZ(K)
      ENDDO
      ENDDO
C
      END SELECT CYLINDER2
C
      ENDDO SPECIES_LOOP
C
C Compute RHS of Temperature Equation
C
      ENERGY: IF (.NOT.ISOTHERMAL) THEN
C
      KDTDX => WORK1
      KDTDY => WORK2
      KDTDZ => WORK3
      KP    => WORK4
C
C Compute thermal conductivity k (KP)
C
      K_DNS_OR_LES: IF (DNS) THEN
C
      IF (.NOT.MIXTURE_FRACTION) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      ITMP = 0.1*TMP(I,J,K)
      K_SUM = K_SPEC(0,ITMP)
      DO N=1,NSPEC
      K_SUM = K_SUM + YYP(I,J,K,N)*(K_SPEC(N,ITMP)-K_SPEC(0,ITMP))
      ENDDO
      KP(I,J,K) = K_SUM
      ENDDO
      ENDDO
      ENDDO
      ENDIF
C
      IF (MIXTURE_FRACTION) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      ITMP = 0.1*TMP(I,J,K)
      IYY  = MAX(0,NINT(YYP(I,J,K,IFUEL)*100.))
      IYY  = MIN(100,IYY)
      KP(I,J,K) = MU(I,J,K)*CP(IYY,ITMP)*RPR
      ENDDO
      ENDDO
      ENDDO
      ENDIF
C
      BOUNDARY_LOOP: DO IW=1,NEWC
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      KP(II,JJ,KK) = KP(IIG,JJG,KKG)
      ENDDO BOUNDARY_LOOP
C
      ELSE K_DNS_OR_LES
C
      KP = MU*CPOPR
C
      ENDIF K_DNS_OR_LES
C
C Compute k*dT/dx, etc
C
      DO K=0,KBAR
      DO J=0,JBAR
      DO I=0,IBAR
      DTDX = (TMP(I+1,J,K)-TMP(I,J,K))*RDXN(I)
      KDTDX(I,J,K) = .5*(KP(I+1,J,K)+KP(I,J,K))*DTDX
      DTDY = (TMP(I,J+1,K)-TMP(I,J,K))*RDYN(J)
      KDTDY(I,J,K) = .5*(KP(I,J+1,K)+KP(I,J,K))*DTDY
      DTDZ = (TMP(I,J,K+1)-TMP(I,J,K))*RDZN(K)
      KDTDZ(I,J,K) = .5*(KP(I,J,K+1)+KP(I,J,K))*DTDZ
      ENDDO
      ENDDO
      ENDDO
C
C Correct thermal gradient (k dT/dn) at boundaries
C
      CORRECTION_LOOP: DO IW=1,NWC
      IF (IV(IW).EQ.0) CYCLE CORRECTION_LOOP
      II  = IJKW(1,IW) 
      JJ  = IJKW(2,IW) 
      KK  = IJKW(3,IW) 
      IIG = IJKW(6,IW) 
      JJG = IJKW(7,IW) 
      KKG = IJKW(8,IW) 
      KW(IW) = KP(IIG,JJG,KKG)
      IOR = IJKW(4,IW)
      SELECT CASE(IOR)
      CASE( 1) ; KDTDX(II,JJ,KK)   = 0.
      CASE(-1) ; KDTDX(II-1,JJ,KK) = 0.
      CASE( 2) ; KDTDY(II,JJ,KK)   = 0.
      CASE(-2) ; KDTDY(II,JJ-1,KK) = 0.
      CASE( 3) ; KDTDZ(II,JJ,KK)   = 0.
      CASE(-3) ; KDTDZ(II,JJ,KK-1) = 0.
      END SELECT
      DP(IIG,JJG,KKG) = DP(IIG,JJG,KKG) - QCONF(IW)*RDN(IW)
      ENDDO CORRECTION_LOOP
C
C Compute (q + del dot k del T) and add to the divergence
C
      CYLINDER3: SELECT CASE(CYLINDRICAL)
C
      CASE(.FALSE.) CYLINDER3   ! 3D or 2D Cartesian
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      ITMP = 0.1*TMP(I,J,K)
      DELKDELT = (KDTDX(I,J,K)-KDTDX(I-1,J,K))*RDX(I) +
     .           (KDTDY(I,J,K)-KDTDY(I,J-1,K))*RDY(J) +
     .           (KDTDZ(I,J,K)-KDTDZ(I,J,K-1))*RDZ(K)
      DP(I,J,K) = DP(I,J,K) + DELKDELT + Q(I,J,K) + QR(I,J,K)
      ENDDO 
      ENDDO
      ENDDO
C
      CASE(.TRUE.) CYLINDER3   ! 2D Cylindrical
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      ITMP = 0.1*TMP(I,J,K)
      DELKDELT =
     . (R(I)*KDTDX(I,J,K)-R(I-1)*KDTDX(I-1,J,K))*RDX(I)*RRN(I) +
     .      (KDTDZ(I,J,K)-       KDTDZ(I,J,K-1))*RDZ(K)
      DP(I,J,K) = DP(I,J,K) + DELKDELT + Q(I,J,K) + QR(I,J,K)
      ENDDO 
      ENDDO
      ENDDO
C
      END SELECT CYLINDER3
C
      ENDIF ENERGY
C
C Correct divergence term to account for temperature dependent CP
C
      RTRM => WORK1
C
      IF (MIXTURE_FRACTION) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      ITMP = 0.1*TMP(I,J,K)
      IYY  = MAX(0,NINT(YYP(I,J,K,IFUEL)*100.))
      IYY  = MIN(100,IYY)
      RTRM(I,J,K) = RP0*RCP(IYY,ITMP)*RSUM(I,J,K)
      DP(I,J,K) = RTRM(I,J,K)*DP(I,J,K)
      ENDDO
      ENDDO 
      ENDDO
      ENDIF
C
      IF (.NOT.MIXTURE_FRACTION .AND. NSPEC.GT.0) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      ITMP = 0.1*TMP(I,J,K)
      CP_SUM = CP(0,ITMP)
      DO N=1,NSPEC
      CP_SUM = CP_SUM + YYP(I,J,K,N)*(CP(N,ITMP)-CP(0,ITMP))
      ENDDO
      RTRM(I,J,K) = RP0*RSUM(I,J,K)/CP_SUM
      DP(I,J,K) = RTRM(I,J,K)*DP(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
C
      IF (.NOT.MIXTURE_FRACTION .AND. NSPEC.EQ.0) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      ITMP = 0.1*TMP(I,J,K)
      CP_SUM = CP(0,ITMP)
      RTRM(I,J,K) = RP0*RCON(0)/CP_SUM
      DP(I,J,K) = RTRM(I,J,K)*DP(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
C
C Add water vapor if sprinklers are on
C
      IF (NLP.GT.0) THEN
C
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      DP(I,J,K) = DP(I,J,K) + (RP0-RTRM(I,J,K))*D_VAP(I,J,K)
      ENDDO
      ENDDO
      ENDDO
C
      ENDIF
C
C Compute normal component of velocity at boundaries
C
      IF (PREDICTOR) THEN
C
      SEALED = .TRUE.
      AREA_V = 0.
      AREA_L = 0.
C
      WALL_LOOP: DO IW=1,NWC
      IOR = IJKW(4,IW)
      WALL_CELL_TYPE: SELECT CASE(IV(IW))
      CASE(0) ; UWS(IW) = 0.
      CASE(1)
         IF (MMT(IJKW(5,IW)).GE.3) CYCLE WALL_LOOP
         IBC = IJKW(5,IW)
         IF (LEAKY(IBC)) AREA_L = AREA_L + AW(IW)
         TSI = T+DT-TW(IW)
         IF (TSI.LT.0.) THEN
            UWS(IW) = 0.
            CYCLE WALL_LOOP
            ENDIF
         SELECT CASE(IOR) 
         CASE( 1); UWS(IW) =-U0 + FV(TSI,ITAUV(IW))*(UW0(IW)+U0)
         CASE(-1); UWS(IW) = U0 + FV(TSI,ITAUV(IW))*(UW0(IW)-U0)
         CASE( 2); UWS(IW) =-V0 + FV(TSI,ITAUV(IW))*(UW0(IW)+V0)
         CASE(-2); UWS(IW) = V0 + FV(TSI,ITAUV(IW))*(UW0(IW)-V0)
         CASE( 3); UWS(IW) =-W0 + FV(TSI,ITAUV(IW))*(UW0(IW)+W0)
         CASE(-3); UWS(IW) = W0 + FV(TSI,ITAUV(IW))*(UW0(IW)-W0)
         END SELECT
      CASE(2)
         SEALED = .FALSE.
         AREA_V = AREA_V + AW(IW)
         II = IJKW(1,IW)
         JJ = IJKW(2,IW)
         KK = IJKW(3,IW)
         SELECT CASE(IOR)
         CASE( 1) ; UWS(IW) = -U(II,JJ,KK)
         CASE(-1) ; UWS(IW) =  U(II-1,JJ,KK)
         CASE( 2) ; UWS(IW) = -V(II,JJ,KK)
         CASE(-2) ; UWS(IW) =  V(II,JJ-1,KK)
         CASE( 3) ; UWS(IW) = -W(II,JJ,KK)
         CASE(-3) ; UWS(IW) =  W(II,JJ,KK-1)
         END SELECT
      CASE(4)
         NM_IN = IJKW(9,IW)
         II_IN = IJKW(10,IW)
         JJ_IN = IJKW(11,IW)
         KK_IN = IJKW(12,IW)
         SEALED = .FALSE.
         SELECT CASE(ABS(IOR))
         CASE(1) 
         XIF = INTERPOLATION_FACTOR(IW,1)
         UWS(IW) = -SIGN(1,IOR)*
     .             (OMESH(NM_IN)%U(II_IN,JJ_IN,KK_IN)  * XIF
     .             +OMESH(NM_IN)%U(II_IN-1,JJ_IN,KK_IN)*(1.-XIF))
         CASE(2) 
         YIF = INTERPOLATION_FACTOR(IW,2)
         UWS(IW) = -SIGN(1,IOR)*
     .             (OMESH(NM_IN)%V(II_IN,JJ_IN,KK_IN)  * YIF
     .             +OMESH(NM_IN)%V(II_IN,JJ_IN-1,KK_IN)*(1.-YIF))
         CASE(3) 
         ZIF = INTERPOLATION_FACTOR(IW,3)
         UWS(IW) = -SIGN(1,IOR)*
     .             (OMESH(NM_IN)%W(II_IN,JJ_IN,KK_IN)  * ZIF
     .             +OMESH(NM_IN)%W(II_IN,JJ_IN,KK_IN-1)*(1.-ZIF))
         END SELECT
      END SELECT WALL_CELL_TYPE
      ENDDO WALL_LOOP
C
      DUWDT(1:NEWC) = RDT*(UWS(1:NEWC)-UW(1:NEWC))
C
      ELSE    ! Corrector Step
C  
      UW = UWS
C
      ENDIF
C
C Sum up divergence integrals for use in dP0/dt calculation
C
      USUM = 0.
      DSUM  = 0.
      PSUM  = 0.
C
      IF (SEALED) THEN
C
      LEAK=LEAK_AREA*SIGN(1._EB,P0-PINF)*SQRT(2.*ABS(P0-PINF)/RHOA)
      IF (P0.GE.PINF) LEAK = MIN(LEAK,(P0-PINF)*RDT*VOL/(GAMMA*P0))
      IF (P0.LT.PINF) LEAK = MAX(LEAK,(P0-PINF)*RDT*VOL/(GAMMA*P0))
      U_LEAK = 0.
      IF (AREA_L.GT.0.) U_LEAK = LEAK/AREA_L
C
      WLOOP3: DO IW=1,NWC
      IF (IV(IW).NE.1) CYCLE WLOOP3
      IBC = IJKW(5,IW)
      IF (LEAKY(IBC) .AND. PREDICTOR) THEN
         UWS(IW) = UWS(IW)+U_LEAK
         IF (IW.LE.NEWC) DUWDT(IW) = RDT*(UWS(IW)-UW(IW))
         ENDIF
      USUM = USUM + UWS(IW)*AW(IW)   ! int u dot dS
      ENDDO WLOOP3
C
      DO K=1,KBAR
      DO J=1,JBAR
      INNER: DO I=1,IBAR
      IF (SOLID(ICA(I,J,K))) CYCLE INNER
      VC   = DX(I)*RC(I)*DY(J)*DZ(K)
      DSUM = DSUM + VC*DP(I,J,K)
      PSUM = PSUM + VC*(RP0-RTRM(I,J,K))
      ENDDO INNER
      ENDDO 
      ENDDO
C
      ENDIF
C
C Compute change in background pressure
C
      IF (PREDICTOR) THEN
C
         IF (SEALED) THEN
            DP0DTS = (DSUM - USUM)/PSUM
            ELSE
            DP0DTS = -GAMMA*(AREA_V/VOL)*SIGN(P0S,P0S-PINF)*
     .                SQRT(2.*ABS(P0S-PINF)/RHOA)
            IF (P0S.GT.PINF) DP0DTS = MAX(DP0DTS,(PINF-P0S)*RDT)
            IF (P0S.LT.PINF) DP0DTS = MIN(DP0DTS,(PINF-P0S)*RDT)
            ENDIF
C
      ELSE   ! Corrector Step
C
         IF (SEALED) THEN
            DP0DT = (DSUM - USUM)/PSUM
            ELSE
            DP0DT = -GAMMA*(AREA_V/VOL)*SIGN(P0,P0-PINF)*
     .               SQRT(2.*ABS(P0-PINF)/RHOA)
            IF (P0.GT.PINF) DP0DT = MAX(DP0DT,(PINF-P0)*RDT)
            IF (P0.LT.PINF) DP0DT = MIN(DP0DT,(PINF-P0)*RDT)
            ENDIF
C
      ENDIF
C
C Add pressure term to divergence
C
      IF (SEALED) THEN
      IF (PREDICTOR) DP = DP + (RTRM-RP0)*DP0DTS
      IF (CORRECTOR) DP = DP + (RTRM-RP0)*DP0DT
      ENDIF
C
C Atmospheric Stratification Case
C
      IF (DT0DZ.NE.-GRAV/CP_GAMMA) THEN
      GFAC = (DT0DZ+GRAV/CP_GAMMA)/TMPA
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      DP(I,J,K) = DP(I,J,K) - 0.5*(W(I,J,K)+W(I,J,K-1))*GFAC
      ENDDO 
      ENDDO 
      ENDDO
      ENDIF
C
C Zero out divergence in solid cells
C
      SOLID_LOOP: DO IC=1,NDBC
      IF (.NOT.SOLID(IC)) CYCLE SOLID_LOOP
      I = I_CELL(IC)
      J = J_CELL(IC)
      K = K_CELL(IC)
      DP(I,J,K) = 0.
      ENDDO SOLID_LOOP
C
C Set divergence in boundary cells
C
      BC_LOOP: DO IW=1,NWC
      IF (IV(IW).EQ.0) CYCLE BC_LOOP
      II = IJKW(1,IW)
      JJ = IJKW(2,IW)
      KK = IJKW(3,IW)
      SELECT CASE(IV(IW))
      CASE(1)
      IF (.NOT.SOLID(ICA(II,JJ,KK))) CYCLE BC_LOOP
      IOR = IJKW(4,IW)
      SELECT CASE(IOR)
      CASE( 1);DP(II,JJ,KK)=DP(II,JJ,KK)-UWS(IW)*RDX(II)*RRN(II)*R(II)
      CASE(-1);DP(II,JJ,KK)=DP(II,JJ,KK)-UWS(IW)*RDX(II)*RRN(II)*R(II-1)
      CASE( 2);DP(II,JJ,KK)=DP(II,JJ,KK)-UWS(IW)*RDY(JJ)
      CASE(-2);DP(II,JJ,KK)=DP(II,JJ,KK)-UWS(IW)*RDY(JJ)
      CASE( 3);DP(II,JJ,KK)=DP(II,JJ,KK)-UWS(IW)*RDZ(KK)
      CASE(-3);DP(II,JJ,KK)=DP(II,JJ,KK)-UWS(IW)*RDZ(KK)
      END SELECT
      CASE(2:4)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      DP(II,JJ,KK) = DP(IIG,JJG,KKG)
      END SELECT
      ENDDO BC_LOOP
C
C Compute time derivative of the divergence
C
      IF (PREDICTOR) THEN
C
      DDDT = (DS-D)*RDT
C
      ELSE
C
      D_OLD => WORK1
C
      D_OLD = DDDT
      DDDT  = (2.*DDDT-DS-D)*RDT
      D     = D_OLD
C
      ENDIF
C
      TUSED(2,NM)=TUSED(2,NM)+SECOND()-TNOW
      END SUBROUTINE DIVERGENCE
C
C
      SUBROUTINE CHECK_DIVERGENCE(NM)
C
C Computes maximum velocity divergence
C
      INTEGER, INTENT(IN) :: NM
      REAL(EB) DIV,RES
C
      TNOW=SECOND()
C
      CALL UNPACK_VAR(NM)
C
      RESMAX = 0.
      DIVMX  = -10000.
      DIVMN  =  10000.
      IMX    = 0
      JMX    = 0
      KMX    = 0
C
      DO K=1,KBAR
      DO J=1,JBAR
      LOOP1: DO I=1,IBAR
      IF (SOLID(ICA(I,J,K))) CYCLE LOOP1
      SELECT CASE(CYLINDRICAL)
      CASE(.FALSE.)
      DIV = (U(I,J,K)-U(I-1,J,K))*RDX(I) + 
     .      (V(I,J,K)-V(I,J-1,K))*RDY(J) +
     .      (W(I,J,K)-W(I,J,K-1))*RDZ(K)
      CASE(.TRUE.)
      DIV = (R(I)*U(I,J,K)-R(I-1)*U(I-1,J,K))*RDX(I)*RRN(I) + 
     .      (W(I,J,K)-W(I,J,K-1))*RDZ(K)
      END SELECT
      RES = ABS(DIV-D(I,J,K))
      IF (ABS(RES).GE.RESMAX) THEN
      RESMAX = ABS(RES)
      IRM=I
      JRM=J
      KRM=K
      ENDIF
      RESMAX = MAX(RES,RESMAX)
      IF (DIV.GE.DIVMX) THEN
         DIVMX = DIV
         IMX=I
         JMX=J
         KMX=K
         ENDIF
      IF (DIV.LT.DIVMN) THEN
         DIVMN = DIV
         IMN=I
         JMN=J
         KMN=K
         ENDIF
      ENDDO LOOP1
      ENDDO
      ENDDO
C
      TUSED(2,NM)=TUSED(2,NM)+SECOND()-TNOW
      END SUBROUTINE CHECK_DIVERGENCE
C
C
      END MODULE DIVG
