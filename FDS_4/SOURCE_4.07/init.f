      MODULE INIT      
C
      USE PREC
      USE VARS
      USE PACKER
      USE CONS
      USE TRAN
      USE RAD
      USE POIS
C
      IMPLICIT NONE
      PRIVATE
      INTEGER IZERO
      CHARACTER(80) MESSAGE
      PUBLIC INITIALIZE_MESH_VARIABLES,INITIALIZE_GLOBAL_VARIABLES,
     .       OPEN_AND_CLOSE
C
C
      CONTAINS
C
C
      SUBROUTINE INITIALIZE_MESH_VARIABLES(NM)
C
      INTEGER N,I,J,K,II,JJ,KK,IPTS,JPTS,KPTS,
     .        NEDGES_DIM,IW,IC,IBC,IOR
      INTEGER, INTENT(IN) :: NM
      REAL(EB) :: MU_N
      INTEGER, POINTER :: IBP1, JBP1, KBP1,
     .                    IBAR, JBAR, KBAR,
     .                    NDWC, NEDGES, NWC
      REAL(EB),POINTER :: XS,XF,YS,YF,ZS,ZF
      TYPE (MESH_TYPE), POINTER :: M
      TYPE (OBSTRUCTION_TYPE), POINTER :: OB
C
      M => MESH(NM)
C
      IBP1 =>M%IBP1
      JBP1 =>M%JBP1
      KBP1 =>M%KBP1
      IBAR =>M%IBAR
      JBAR =>M%JBAR
      KBAR =>M%KBAR
      NDWC =>M%NDWC
       NWC =>M%NWC
      NEDGES=>M%NEDGES
      XS=>M%XS ; YS=>M%YS ; ZS=>M%ZS
      XF=>M%XF ; YF=>M%YF ; ZF=>M%ZF
C
      ALLOCATE(  M%RHO(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','RHO',IZERO)
      ALLOCATE( M%RHOS(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','RHOS',IZERO)
      ALLOCATE(  M%TMP(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','TMP',IZERO)
      ALLOCATE( M%FRHO(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','FRHO',IZERO)
C
      ALLOCATE( M%U(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','U',IZERO)
      ALLOCATE( M%V(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','V',IZERO)
      ALLOCATE( M%W(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','W',IZERO)
      ALLOCATE(M%US(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','US',IZERO)
      ALLOCATE(M%VS(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','VS',IZERO)
      ALLOCATE(M%WS(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','WS',IZERO)
      ALLOCATE(M%FVX(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','FVX',IZERO)
      ALLOCATE(M%FVY(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','FVY',IZERO)
      ALLOCATE(M%FVZ(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','FVZ',IZERO)
C
      ALLOCATE(   M%H(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','H',IZERO)
      ALLOCATE(M%DDDT(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','DDDT',IZERO)
      ALLOCATE(   M%D(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','D',IZERO)
      ALLOCATE(  M%DS(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','DS',IZERO)
      ALLOCATE( M%MU(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','MU',IZERO)
C
      ALLOCATE(   M%Q(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','Q',IZERO)
C
C Allocate species arrays
C
      IF (NSPEC.GT.0) THEN
      ALLOCATE( M%YY(0:IBP1,0:JBP1,0:KBP1,NSPEC),STAT=IZERO)
      CALL ChkMemErr('INIT','YY',IZERO)
      ALLOCATE(M%YYS(0:IBP1,0:JBP1,0:KBP1,NSPEC),STAT=IZERO)
      CALL ChkMemErr('INIT','YYS',IZERO)
      ALLOCATE(M%FYY(0:IBP1,0:JBP1,0:KBP1,NSPEC),STAT=IZERO)
      CALL ChkMemErr('INIT','FYY',IZERO)
      ALLOCATE(M%RSUM(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','RSUM',IZERO)
      ENDIF
C
C Allocate water mass arrays if sprinklers are present
C
      IF (DROPLET_FILE) THEN
      ALLOCATE(M%AVG_DROP_DEN(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','AVG_DROP_DEN',IZERO) ; M%AVG_DROP_DEN=0.
      ALLOCATE(M%AVG_DROP_TMP(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','AVG_DROP_TMP',IZERO) ; M%AVG_DROP_TMP=TMPM
      ALLOCATE(M%AVG_DROP_RAD(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','AVG_DROP_RAD',IZERO) ; M%AVG_DROP_RAD=0.
      ALLOCATE(M%QR_W(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','QR_W',IZERO) ; M%QR_W = 0.
      ALLOCATE(M%D_VAP(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','D_VAP',IZERO) ; M%D_VAP = 0.
      ENDIF
C
C If radiation absorption desired allocate arrays
C
      ALLOCATE(M%QR(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','QR',IZERO)
      ALLOCATE(M%KAPPA(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','KAPPA',IZERO) 
      ALLOCATE(M%UII(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','UII',IZERO)
C
C Work arrays
C
      ALLOCATE(M%WORK1(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','WORK1',IZERO)
      ALLOCATE(M%WORK2(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','WORK2',IZERO)
      ALLOCATE(M%WORK3(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','WORK3',IZERO)
      ALLOCATE(M%WORK4(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','WORK4',IZERO)
      ALLOCATE(M%WORK5(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','WORK5',IZERO)
      ALLOCATE(M%WORK6(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','WORK6',IZERO)
      ALLOCATE(M%WORK7(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO)
      CALL ChkMemErr('INIT','WORK7',IZERO)
C
C Boundary file patch counter
C
      ALLOCATE(M%INC(-3:3,0:M%NB))
      CALL ChkMemErr('INIT','INC',IZERO)
C
C Initialize all variables
C
      M%DTOLD  = M%DT
      M%DTNEXT = M%DT
      M%DTINT  = M%DT
      M%NEW_TIME_STEP = .FALSE.
C
      M%RHO   = RHOA
      M%RHOS  = M%RHO
      M%RHO_AVG = RHOA
      M%TMP   = TMPA
      M%FRHO  = 0.
      M%U     = U0
      M%V     = V0
      M%W     = W0
      M%US    = U0
      M%VS    = V0
      M%WS    = W0
      IF (NOISE) CALL INITIAL_NOISE
      M%FVX   = 0.
      M%FVY   = 0.
      M%FVZ   = 0.
      M%H     = H0
      M%DDDT  = 0.
      M%D     = 0.
      M%DS    = 0.
      M%Q     = 0.
C
C Background pressure
C
      M%P0     = PINF
      M%P0S    = M%P0
      M%DP0DT  = 0.
      M%DP0DTS = M%DP0DT
C
C Upper bounds on HRR
C
      IF (TWO_D) THEN
      M%Q_UPPER = HRRPUA_SHEET*1000./(M%DXMIN*M%DZMIN)**0.5
      ELSE
      M%Q_UPPER = HRRPUA_SHEET*1000./(M%DXMIN*M%DYMIN*M%DZMIN)**ONTH
      ENDIF
C
C Viscosity
C
      MU_N = MU_SPEC(0,NINT(0.1*TMPA))
      M%MU = MU_N
C
      IF (DNS .AND. (ISOTHERMAL .OR. INCOMPRESSIBLE)) THEN
C
      ALLOCATE(M%RREDX(M%IBAR),STAT=IZERO)
      CALL ChkMemErr('READ','RREDX',IZERO)
      ALLOCATE(M%RREDY(M%JBAR),STAT=IZERO)
      CALL ChkMemErr('READ','RREDY',IZERO)
      ALLOCATE(M%RREDZ(M%KBAR),STAT=IZERO)
      CALL ChkMemErr('READ','RREDZ',IZERO)
C
      DO I=1,M%IBAR
      M%RREDX(I) = (MU_N/RHOA)*M%RDX(I)
      ENDDO
      DO J=1,M%JBAR
      M%RREDY(J) = (MU_N/RHOA)*M%RDY(J)
      ENDDO
      DO K=1,M%KBAR
      M%RREDZ(K) = (MU_N/RHOA)*M%RDZ(K)
      ENDDO
C
      ENDIF
C
      IF (NSPEC.GT.0) THEN
         M%RSUM = RSUM0
         DO N=1,NSPEC
         M%YY(:,:,:,N)  = YY0(N)
         M%YYS(:,:,:,N) = YY0(N)
         M%FYY(:,:,:,N) = 0.
         ENDDO
         ENDIF
C
C Set initial value blocks
C
      DO N=1,NIB
      DO K=0,KBP1
      DO J=0,JBP1
      DO I=0,IBP1
      IF (M%XC(I).GT.XB1(N) .AND. M%XC(I).LT.XB2(N) .AND.
     .    M%YC(J).GT.YB1(N) .AND. M%YC(J).LT.YB2(N) .AND. 
     .    M%ZC(K).GT.ZB1(N) .AND. M%ZC(K).LT.ZB2(N)) THEN
          SELECT CASE(INIT_INDEX(N))
          CASE(-1) ; M%TMP(I,J,K)  = INIT_VALUE(N)
c                    M%RHO(I,J,K)  = M%P0/(M%TMP(I,J,K)*RSUM0)
          CASE(-2) ; M%RHO(I,J,K)  = INIT_VALUE(N)
c                    M%TMP(I,J,K)  = M%P0/(M%RHO(I,J,K)*RSUM0)
          CASE(1:) ; M%YY(I,J,K,INIT_INDEX(N)) = INIT_VALUE(N)
          END SELECT
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
C Radiation
C
      M%QR     = 0.
      M%KAPPA  = 0.
      M%UII    = 4.*SIGMA*TMPA4

      M%WORK1 = 0.
      M%WORK2 = 0.
      M%WORK3 = 0.
      M%WORK4 = 0.
      M%WORK5 = 0.
      M%WORK6 = 0.
      M%WORK7 = 0.
C
C Stoichiometric Z
C
      IF (MIXTURE_FRACTION) M%Z_F_EFF = Z_F
C
C Allocate grid cell arrays
C
      ALLOCATE(M%CELL_MASS(M%NDBC),STAT=IZERO)
      CALL ChkMemErr('INIT','CELL_MASS',IZERO); M%CELL_MASS = 1.E10
C
C Designate each boundary cell with a reference number for wall BC's
C
      NWC  = 0
      NDWC = 0
C
C Determine the number of wall cells to allocate
C
      OBST_LOOP: DO N=1,M%NB
      OB=>M%OBSTRUCTION(N)
C
      IF ( BURNAWAY(OB%IBC( 1)) .OR.
     .     BURNAWAY(OB%IBC(-1)) .OR.
     .     BURNAWAY(OB%IBC( 2)) .OR.
     .     BURNAWAY(OB%IBC(-2)) .OR.
     .     BURNAWAY(OB%IBC( 3)) .OR.
     .     BURNAWAY(OB%IBC(-3)) ) THEN
C
      NDWC = NDWC +
     .       3*(OB%I2-OB%I1+1)*(OB%J2-OB%J1+1)*(OB%K2-OB%K1+1)
C
      ELSE IF (OB%T_REMOVE.LT.100000. .OR. 
     .         OB%HEAT_INDEX_REMOVE.GT.0) THEN
C
      NDWC = NDWC +
     .       3*(OB%I2-OB%I1+1)*(OB%J2-OB%J1+1)*(OB%K2-OB%K1+1)
C
      ELSE
C
      DO K=OB%K1+1,OB%K2
      DO J=OB%J1+1,OB%J2
      IC = M%ICA(OB%I1  ,J,K)
      IF (.NOT.M%SOLID(IC)) NDWC = NDWC + 1
      IC = M%ICA(OB%I2+1,J,K)
      IF (.NOT.M%SOLID(IC)) NDWC = NDWC + 1
      ENDDO 
      ENDDO
C
      DO K=OB%K1+1,OB%K2
      DO I=OB%I1+1,OB%I2
      IC = M%ICA(I,OB%J1  ,K)
      IF (.NOT.M%SOLID(IC)) NDWC = NDWC + 1
      IC = M%ICA(I,OB%J2+1,K)
      IF (.NOT.M%SOLID(IC)) NDWC = NDWC + 1
      ENDDO 
      ENDDO
C
      DO J=OB%J1+1,OB%J2
      DO I=OB%I1+1,OB%I2
      IC = M%ICA(I,J,OB%K1  )
      IF (.NOT.M%SOLID(IC)) NDWC = NDWC + 1
      IC = M%ICA(I,J,OB%K2+1)
      IF (.NOT.M%SOLID(IC)) NDWC = NDWC + 1
      ENDDO 
      ENDDO
C
      ENDIF
C
      ENDDO OBST_LOOP
C
      NDWC = NDWC + M%NEWC
C
      ALLOCATE(M%WALL(1:NDWC))
      CALL ChkMemErr('INIT','WALL',IZERO)
C
      ALLOCATE(M%TMP_F(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','TMP_F',IZERO) ; M%TMP_F = TMPA
      ALLOCATE(M%TMP_B(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','TMP_B',IZERO) ; M%TMP_B = TMPA
      ALLOCATE(M%TMP_W(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','TMP_W',IZERO) ; M%TMP_W = TMPA
      ALLOCATE(M%RHO_W(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','RHO_W',IZERO) ; M%RHO_W = RHOA
      ALLOCATE(M%RSUM_W(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','RSUM_W',IZERO) ; M%RSUM_W = RSUM0
      ALLOCATE(M%YY_W(NDWC,NSPEC),STAT=IZERO)
      CALL ChkMemErr('INIT','YY_W',IZERO)  
      DO IW=1,NDWC
      M%YY_W(IW,1:NSPEC) = YY0(1:NSPEC)
      ENDDO
C
      ALLOCATE(M%PFRONT_I(NDWC),STAT=IZERO)
      CALL ChkMemErr('READ','PFRONT_I',IZERO) ; M%PFRONT_I = 1.
      ALLOCATE(M%MFRONT_I(NDWC),STAT=IZERO)
      CALL ChkMemErr('READ','MFRONT_I',IZERO) ; M%MFRONT_I = 1.
      ALLOCATE(M%E_WALL(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','E_WALL',IZERO) ; M%E_WALL = 1.
      ALLOCATE(M%QPYR(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','QPYR',IZERO) ; M%QPYR = 0.
      ALLOCATE(M%QRAD(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','QRAD',IZERO) ; M%QRAD = 0.
      ALLOCATE(M%QCONF(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','QCONF',IZERO) ; M%QCONF = 0.
      ALLOCATE(M%QCONB(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','QCONB',IZERO) ; M%QCONB = 0.
      ALLOCATE(M%HEAT_TRANS_COEF(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','HEAT_TRANS_COEF',IZERO) 
      M%HEAT_TRANS_COEF = 0.
      ALLOCATE(M%XW(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','XW',IZERO)
      ALLOCATE(M%YW(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','YW',IZERO)
      ALLOCATE(M%ZW(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','ZW',IZERO)
      ALLOCATE(M%TW(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','TW',IZERO)  ; M%TW = 0.
      ALLOCATE(M%EW(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','EW',IZERO)  ; M%EW = 0.
      ALLOCATE(M%KW(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','KW',IZERO)  ; M%KW = 1.
      ALLOCATE(M%RHODW(NDWC,NSPEC),STAT=IZERO)
      CALL ChkMemErr('INIT','RHODW',IZERO) ; M%RHODW = 1.
      ALLOCATE(M%AREA_ADJUST(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','AREA_ADJUST',IZERO)  ; M%AREA_ADJUST = 1.
      ALLOCATE(M%MASSFLUX(NDWC,0:NSPEC),STAT=IZERO)
      CALL ChkMemErr('INIT','MASSFLUX',IZERO)  ; M%MASSFLUX = 0.
      ALLOCATE(M%RDN(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','RDN',IZERO) ; M%RDN = 1.
      ALLOCATE(M%AW(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','AW',IZERO)  ; M%AW = 0.
      ALLOCATE(M%RAW(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','RAW',IZERO) ; M%RAW = 0.
      ALLOCATE(M%NPPCW(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','NPPCW',IZERO) ; M%NPPCW = 1
      ALLOCATE(M%ITAUV(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','ITAUV',IZERO) ; M%ITAUV = 0
      ALLOCATE(M%ITAUQ(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','ITAUQ',IZERO) ; M%ITAUQ = 0
      ALLOCATE(M%ITAUMF(NDWC,0:NSPEC),STAT=IZERO)
      CALL ChkMemErr('INIT','ITAUMF',IZERO) ; M%ITAUMF = 0
      ALLOCATE(M%UW0(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','UW0',IZERO)    ; M%UW0 = 0.
      ALLOCATE(M%UW(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','UW',IZERO)    ; M%UW = 0.
      ALLOCATE(M%UWS(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','UWS',IZERO)   ; M%UWS = 0.
      ALLOCATE(M%CELL_PARTICLE_CLASS(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','CELL_PARTICLE_CLASS',IZERO) 
      M%CELL_PARTICLE_CLASS = -1
      ALLOCATE(M%OBST_INDEX(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','OBST_INDEX',IZERO) 
      M%OBST_INDEX = 0
      ALLOCATE(M%DUWDT(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','DUWDT',IZERO) ; M%DUWDT = 0.
      ALLOCATE(M%IJKW(12,NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','IJKW',IZERO)  ; M%IJKW = 0
      ALLOCATE(M%INTERPOLATION_FACTOR(M%NEWC,3),STAT=IZERO)
      CALL ChkMemErr('INIT','INTERPOLATION_FACTOR',IZERO)  
      ALLOCATE(M%CELL_VOLUME_RATIO(M%NEWC),STAT=IZERO)
      CALL ChkMemErr('INIT','CELL_VOLUME_RATIO',IZERO)  
      ALLOCATE(M%IV(0:NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','IV',IZERO) ; M%IV = 0 
      ALLOCATE(M%MASS_LOSS(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','MASS_LOSS',IZERO); M%MASS_LOSS = 0.
      ALLOCATE(M%IWA(0:M%NDBC,-3:3),STAT=IZERO)
      CALL ChkMemErr('INIT','IWA',IZERO) ; M%IWA = 0
      ALLOCATE(M%IEA(0:M%NDBC,1:12),STAT=IZERO)
      CALL ChkMemErr('INIT','IEA',IZERO) ; M%IEA = 0
C
C Surface water arrays
C
      IF (ACCUMULATE_WATER) THEN
      ALLOCATE(M%AWMPUA(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','AWMPUA',IZERO) ; M%AWMPUA = 0.
      ENDIF
      ALLOCATE(M%WMPUA(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','WMPUA',IZERO) ; M%WMPUA = 0.
      ALLOCATE(M%WCPUA(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','WCPUA',IZERO) ; M%WCPUA = 0.
C
C Surface work arrays
C
      ALLOCATE(M%WALL_WORK1(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','WALL_WORK1',IZERO) 
      ALLOCATE(M%WALL_WORK2(NDWC),STAT=IZERO)
      CALL ChkMemErr('INIT','WALL_WORK2',IZERO) 
C
C Set up boundary arrays for external walls
C
      DO K=1,KBAR
      DO J=1,JBAR
      I   = 0
      IBC = IBCDEF
      IF (M%SOLID(M%ICA(1,J,K))) IBC = 0
      IOR = 1
      NWC = NWC + 1
      IW  = NWC
      CALL INIT_WALL_CELL(NM,I,J,K,0,IW,IOR,IBC)
      ENDDO
      ENDDO
      DO K=1,KBAR
      DO J=1,JBAR
      I   = IBP1
      IBC = IBCDEF
      IF (M%SOLID(M%ICA(IBAR,J,K))) IBC = 0
      IOR = -1
      NWC = NWC + 1
      IW  = NWC
      CALL INIT_WALL_CELL(NM,I,J,K,0,IW,IOR,IBC)
      ENDDO
      ENDDO
C
      DO K=1,KBAR
      DO I=1,IBAR
      J   = 0
      IBC = IBCDEF
      IF (M%SOLID(M%ICA(I,1,K))) IBC = 0
      IOR = 2
      NWC = NWC + 1
      IW  = NWC
      CALL INIT_WALL_CELL(NM,I,J,K,0,IW,IOR,IBC)
      ENDDO
      ENDDO
      DO K=1,KBAR
      DO I=1,IBAR
      J   = JBP1
      IBC = IBCDEF
      IF (M%SOLID(M%ICA(I,JBAR,K))) IBC = 0
      IOR = -2
      NWC = NWC + 1
      IW  = NWC
      CALL INIT_WALL_CELL(NM,I,J,K,0,IW,IOR,IBC)
      ENDDO
      ENDDO
C
      DO J=1,JBAR
      DO I=1,IBAR
      K   = 0
      IBC = IBCDEF
      IF (M%SOLID(M%ICA(I,J,1))) IBC = 0
      IOR = 3
      NWC = NWC + 1
      IW  = NWC
      CALL INIT_WALL_CELL(NM,I,J,K,0,IW,IOR,IBC)
      ENDDO
      ENDDO
      DO J=1,JBAR
      DO I=1,IBAR
      K   = KBP1
      IBC = IBCDEF
      IF (M%SOLID(M%ICA(I,J,KBAR))) IBC = 0
      IOR = -3
      NWC = NWC + 1
      IW  = NWC
      CALL INIT_WALL_CELL(NM,I,J,K,0,IW,IOR,IBC)
      ENDDO
      ENDDO
C
C Set up boundary arrays for internal walls
C
      NBLOOP: DO N=1,M%NB
      OB=>M%OBSTRUCTION(N)
C
      IF (OB%HIDDEN) CYCLE NBLOOP
C
      DO          K=OB%K1+1,OB%K2
      LOOP251: DO J=OB%J1+1,OB%J2
      I   = OB%I1+1
      IC  = M%ICA(I-1,J,K)
      IF (M%SOLID(IC) .OR. I.EQ.1) CYCLE LOOP251
      IOR = -1 
      IBC = OB%IBC(IOR)
      IW  = M%IWA(IC,-IOR)
      IF (IW.EQ.0) THEN
         NWC = NWC + 1
         IW  = NWC
         ENDIF
      CALL INIT_WALL_CELL(NM,I,J,K,N,IW,IOR,IBC)
      ENDDO LOOP251
      ENDDO
C
      DO          K=OB%K1+1,OB%K2
      LOOP252: DO J=OB%J1+1,OB%J2
      I   = OB%I2
      IC  = M%ICA(I+1,J,K)
      IF (M%SOLID(IC) .OR. I.EQ.M%IBAR) CYCLE LOOP252
      IOR = 1
      IBC = OB%IBC(IOR)
      IW  = M%IWA(IC,-IOR)
      IF (IW.EQ.0) THEN
         NWC = NWC + 1
         IW  = NWC
         ENDIF
      CALL INIT_WALL_CELL(NM,I,J,K,N,IW,IOR,IBC)
      ENDDO LOOP252
      ENDDO
C
      DO          K=OB%K1+1,OB%K2
      LOOP253: DO I=OB%I1+1,OB%I2
      J   = OB%J1+1
      IC  = M%ICA(I,J-1,K)
      IF (M%SOLID(IC) .OR. J.EQ.1) CYCLE LOOP253
      IOR = -2
      IBC = OB%IBC(IOR)
      IW  = M%IWA(IC,-IOR)
      IF (IW.EQ.0) THEN
         NWC = NWC + 1
         IW  = NWC
         ENDIF
      CALL INIT_WALL_CELL(NM,I,J,K,N,IW,IOR,IBC)
      ENDDO LOOP253
      ENDDO   
C
      DO          K=OB%K1+1,OB%K2
      LOOP254: DO I=OB%I1+1,OB%I2
      J   = OB%J2
      IC  = M%ICA(I,J+1,K)
      IF (M%SOLID(IC) .OR. J.EQ.M%JBAR) CYCLE LOOP254
      IOR = 2
      IBC = OB%IBC(IOR)
      IW  = M%IWA(IC,-IOR)
      IF (IW.EQ.0) THEN
         NWC = NWC + 1
         IW  = NWC
         ENDIF
      CALL INIT_WALL_CELL(NM,I,J,K,N,IW,IOR,IBC)
      ENDDO LOOP254
      ENDDO   
C
      DO          J=OB%J1+1,OB%J2
      LOOP255: DO I=OB%I1+1,OB%I2
      K   = OB%K1+1
      IC  = M%ICA(I,J,K-1)
      IF (M%SOLID(IC) .OR. K.EQ.1) CYCLE LOOP255
      IOR = -3
      IBC = OB%IBC(IOR)
      IW  = M%IWA(IC,-IOR)
      IF (IW.EQ.0) THEN
         NWC = NWC + 1
         IW  = NWC
         ENDIF
      CALL INIT_WALL_CELL(NM,I,J,K,N,IW,IOR,IBC)
      ENDDO LOOP255
      ENDDO   
C
      DO          J=OB%J1+1,OB%J2
      LOOP256: DO I=OB%I1+1,OB%I2
      K   = OB%K2
      IC  = M%ICA(I,J,K+1)
      IF (M%SOLID(IC) .OR. K.EQ.M%KBAR) CYCLE LOOP256
      IOR = 3 
      IBC = OB%IBC(IOR)
      IW  = M%IWA(IC,-IOR)
      IF (IW.EQ.0) THEN
         NWC = NWC + 1
         IW  = NWC
         ENDIF
      CALL INIT_WALL_CELL(NM,I,J,K,N,IW,IOR,IBC)
      ENDDO LOOP256
      ENDDO   
C
      ENDDO NBLOOP
C
      M%WALLCLK  = 0.
      M%WALL_COUNTER = 0
C
C Allocate arrays for storing velocity boundary condition info
C
      NEDGES_DIM = 4*(IBP1*JBP1+IBP1*KBP1+JBP1*KBP1)
      DO N=1,M%NB
      OB=>M%OBSTRUCTION(N)
      IPTS = OB%I2-OB%I1
      JPTS = OB%J2-OB%J1
      KPTS = OB%K2-OB%K1
      NEDGES_DIM = NEDGES_DIM + 4*(IPTS*JPTS+IPTS*KPTS+JPTS*KPTS)
      ENDDO
C
      ALLOCATE(M%IJKE(10,NEDGES_DIM),STAT=IZERO)
      CALL ChkMemErr('INIT','IJKE',IZERO)   ; M%IJKE  = 0
      ALLOCATE(M%OME_E(NEDGES_DIM),STAT=IZERO)
      CALL ChkMemErr('INIT','OME_E',IZERO)  ; M%OME_E = 0.
      ALLOCATE(M%TAU_E(NEDGES_DIM),STAT=IZERO)
      CALL ChkMemErr('INIT','TAU_E',IZERO)  ; M%TAU_E = 0.
C
      M%NLP = 0
      M%NLPDIM = 1000
      IF (DROPLET_FILE) THEN
      ALLOCATE(M%DROPLET(M%NLPDIM),STAT=IZERO)
      CALL ChkMemErr('INIT','DROPLET',IZERO)
      ENDIF
C
C Allocate array to hold character strings for Smokeview file
C
      M%N_STRINGS     =   0
      M%N_STRINGS_MAX = 100
      ALLOCATE(M%STRING(M%N_STRINGS_MAX),STAT=IZERO)
      CALL ChkMemErr('INIT','STRING',IZERO)
C
C Set up arrays to hold velocity boundary condition info
C
      CALL INITIALIZE_EDGES
C
C Initialize Pressure solver
C
      CALL INITIALIZE_POISSON_SOLVER
C
C Set up thermocouple arrays
C
      CALL INITIALIZE_THCP
C
C Allocate and Initialize Mesh-Dependent Radiation Arrays
C
      M%ANGLE_INC_COUNTER = 0
      M%RAD_CALL_COUNTER  = 0
C
      ALLOCATE(M%UIID(0:M%IBP1,0:M%JBP1,0:M%KBP1,1:UIIDIM),STAT=IZERO)
      CALL ChkMemErr('INIT','UIID',IZERO)
C
      DO IW=1,M%NDWC
      IF (M%IV(IW).NE.2) THEN
         ALLOCATE(M%WALL(IW)%ILW(NRA,NSB),STAT=IZERO)
         CALL ChkMemErr('INIT','ILW',IZERO)
         M%WALL(IW)%ILW  = SIGMA*TMPA4*RPI
         ENDIF
      ENDDO
C
      M%UIID = 0.
C
C Initialize Mesh Exchange
C
      CALL INITIALIZE_INTERPOLATION
C
      CONTAINS
C
C
      SUBROUTINE INITIALIZE_EDGES
C
C Set up edge arrays for velocity boundary conditions
C
      INTEGER I,J,K,IW,IIG,JJG,KKG,NOM,
     .        II,JJ,KK,IOR,IC,ICG,IC1,IC2,IBC,NN,IW1,IW2
      LOGICAL :: SOLID_ONLY
C
      CALL UNPACK_VAR(NM)
C
      NEDGES = 0
C
      WALL_CELL_LOOP: DO IW=1,NWC
C
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IOR = IJKW(4,IW)
      IBC = IJKW(5,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      IC  = ICA(II,JJ,KK)
      ICG = ICA(IIG,JJG,KKG)
C
      SELECT CASE(IOR)
C
      CASE( 1)
      SOLID_ONLY = .FALSE.
      CALL GET_OBST(NM,II,JJ,KK,IOR,SOLID_ONLY,NN)
      OB=>OBSTRUCTION(NN)
      IC1 = ICA(II,JJ,KK+1)
      IC2 = ICA(IIG,JJG,KKG+1)
      IW1 = IWA(ICG,-IOR)
      IW2 = IWA(IC2,-IOR)
      NOM = IJKW(9,IW)
      IF (IW1.GT.0 .AND. IW2.GT.0) THEN
         NEDGES = NEDGES + 1
         IJKE(1,NEDGES) = II
         IJKE(2,NEDGES) = JJ
         IJKE(3,NEDGES) = KK
         IJKE(4,NEDGES) = 2
         IJKE(5,NEDGES) = IBC
         IF (IJKW(5,IW1).EQ.INTERPOLATED_INDEX .OR.
     .       IJKW(5,IW2).EQ.INTERPOLATED_INDEX) THEN
             IF (NOM.GT.0) IJKE(5,NEDGES) = INTERPOLATED_INDEX
             ENDIF
         IJKE(6,NEDGES) = IOR
         IF (.NOT.OB%SAWTOOTH) IJKE(6,NEDGES) = 0
         IEA(IC ,8) = NEDGES
         IEA(IC1,6) = NEDGES
         IEA(IC2,5) = NEDGES
         IEA(ICG,7) = NEDGES
         IF (NOM.GT.0) THEN
            IJKE( 7,NEDGES) = NOM
            IJKE( 8,NEDGES) = IJKW(10,IW)
            IJKE( 9,NEDGES) = IJKW(11,IW)
            IJKE(10,NEDGES) = IJKW(12,IW)
            ENDIF
         ENDIF
      IC1 = ICA(II,JJ+1,KK)
      IC2 = ICA(IIG,JJG+1,KKG)
      IW1 = IWA(ICG,-IOR)
      IW2 = IWA(IC2,-IOR)
      NOM = IJKW(9,IW)
      IF (IW1.GT.0 .AND. IW2.GT.0) THEN
         NEDGES = NEDGES + 1
         IJKE(1,NEDGES) = II
         IJKE(2,NEDGES) = JJ
         IJKE(3,NEDGES) = KK
         IJKE(4,NEDGES) = 3
         IJKE(5,NEDGES) = IBC
         IF (IJKW(5,IW1).EQ.INTERPOLATED_INDEX .OR.
     .       IJKW(5,IW2).EQ.INTERPOLATED_INDEX) THEN
             IF (NOM.GT.0) IJKE(5,NEDGES) = INTERPOLATED_INDEX
             ENDIF
         IJKE(6,NEDGES) = IOR
         IF (.NOT.OB%SAWTOOTH) IJKE(6,NEDGES) = 0
         IEA(IC ,12) = NEDGES
         IEA(IC1,10) = NEDGES
         IEA(IC2, 9) = NEDGES
         IEA(ICG,11) = NEDGES
         IF (NOM.GT.0) THEN
            IJKE( 7,NEDGES) = NOM
            IJKE( 8,NEDGES) = IJKW(10,IW)
            IJKE( 9,NEDGES) = IJKW(11,IW)
            IJKE(10,NEDGES) = IJKW(12,IW)
            ENDIF
         ENDIF
C
      CASE(-1)
      SOLID_ONLY = .FALSE.
      CALL GET_OBST(NM,II-1,JJ,KK,IOR,SOLID_ONLY,NN)
      OB=>OBSTRUCTION(NN)
      IC1 = ICA(II,JJ,KK+1)
      IC2 = ICA(IIG,JJG,KKG+1)
      IW1 = IWA(ICG,-IOR)
      IW2 = IWA(IC2,-IOR)
      NOM = IJKW(9,IW)
      IF (IW1.GT.0 .AND. IW2.GT.0) THEN
         NEDGES = NEDGES + 1
         IJKE(1,NEDGES) = II-1
         IJKE(2,NEDGES) = JJ
         IJKE(3,NEDGES) = KK
         IJKE(4,NEDGES) = 2
         IJKE(5,NEDGES) = IBC
         IF (IJKW(5,IW1).EQ.INTERPOLATED_INDEX .OR.
     .       IJKW(5,IW2).EQ.INTERPOLATED_INDEX) THEN
             IF (NOM.GT.0) IJKE(5,NEDGES) = INTERPOLATED_INDEX
             ENDIF
         IJKE(6,NEDGES) = IOR
         IF (.NOT.OB%SAWTOOTH) IJKE(6,NEDGES) = 0
         IEA(IC ,7) = NEDGES
         IEA(IC1,5) = NEDGES
         IEA(IC2,6) = NEDGES
         IEA(ICG,8) = NEDGES
         IF (NOM.GT.0) THEN
            IJKE( 7,NEDGES) = NOM
            IJKE( 8,NEDGES) = IJKW(10,IW)
            IJKE( 9,NEDGES) = IJKW(11,IW)
            IJKE(10,NEDGES) = IJKW(12,IW)
            ENDIF
         ENDIF
      IC1 = ICA(II,JJ+1,KK)
      IC2 = ICA(IIG,JJG+1,KKG)
      IW1 = IWA(ICG,-IOR)
      IW2 = IWA(IC2,-IOR)
      NOM = IJKW(9,IW)
      IF (IW1.GT.0 .AND. IW2.GT.0) THEN
         NEDGES = NEDGES + 1
         IJKE(1,NEDGES) = II-1
         IJKE(2,NEDGES) = JJ
         IJKE(3,NEDGES) = KK
         IJKE(4,NEDGES) = 3
         IJKE(5,NEDGES) = IBC
         IF (IJKW(5,IW1).EQ.INTERPOLATED_INDEX .OR.
     .       IJKW(5,IW2).EQ.INTERPOLATED_INDEX) THEN
             IF (NOM.GT.0) IJKE(5,NEDGES) = INTERPOLATED_INDEX
             ENDIF
         IJKE(6,NEDGES) = IOR
         IF (.NOT.OB%SAWTOOTH) IJKE(6,NEDGES) = 0
         IEA(IC ,11) = NEDGES
         IEA(IC1, 9) = NEDGES
         IEA(IC2,10) = NEDGES
         IEA(ICG,12) = NEDGES
         IF (NOM.GT.0) THEN
            IJKE( 7,NEDGES) = NOM
            IJKE( 8,NEDGES) = IJKW(10,IW)
            IJKE( 9,NEDGES) = IJKW(11,IW)
            IJKE(10,NEDGES) = IJKW(12,IW)
            ENDIF
         ENDIF
C
      CASE( 2)
      SOLID_ONLY = .FALSE.
      CALL GET_OBST(NM,II,JJ,KK,IOR,SOLID_ONLY,NN)
      OB=>OBSTRUCTION(NN)
      IC1 = ICA(II,JJ,KK+1)
      IC2 = ICA(IIG,JJG,KKG+1)
      IW1 = IWA(ICG,-IOR)
      IW2 = IWA(IC2,-IOR)
      NOM = IJKW(9,IW)
      IF (IW1.GT.0 .AND. IW2.GT.0) THEN
         NEDGES = NEDGES + 1
         IJKE(1,NEDGES) = II
         IJKE(2,NEDGES) = JJ
         IJKE(3,NEDGES) = KK
         IJKE(4,NEDGES) = 1
         IJKE(5,NEDGES) = IBC
         IF (IJKW(5,IW1).EQ.INTERPOLATED_INDEX .OR.
     .       IJKW(5,IW2).EQ.INTERPOLATED_INDEX) THEN
             IF (NOM.GT.0) IJKE(5,NEDGES) = INTERPOLATED_INDEX
             ENDIF
         IJKE(6,NEDGES) = IOR
         IF (.NOT.OB%SAWTOOTH) IJKE(6,NEDGES) = 0
         IEA(IC ,4) = NEDGES
         IEA(IC1,2) = NEDGES
         IEA(IC2,1) = NEDGES
         IEA(ICG,3) = NEDGES
         IF (NOM.GT.0) THEN
            IJKE( 7,NEDGES) = NOM
            IJKE( 8,NEDGES) = IJKW(10,IW)
            IJKE( 9,NEDGES) = IJKW(11,IW)
            IJKE(10,NEDGES) = IJKW(12,IW)
            ENDIF
         ENDIF
      IC1 = ICA(II+1,JJ,KK)
      IC2 = ICA(IIG+1,JJG,KKG)
      IW1 = IWA(ICG,-IOR)
      IW2 = IWA(IC2,-IOR)
      NOM = IJKW(9,IW)
      IF (IW1.GT.0 .AND. IW2.GT.0) THEN
         NEDGES = NEDGES + 1
         IJKE(1,NEDGES) = II
         IJKE(2,NEDGES) = JJ
         IJKE(3,NEDGES) = KK
         IJKE(4,NEDGES) = 3
         IJKE(5,NEDGES) = IBC
         IF (IJKW(5,IW1).EQ.INTERPOLATED_INDEX .OR.
     .       IJKW(5,IW2).EQ.INTERPOLATED_INDEX) THEN
             IF (NOM.GT.0) IJKE(5,NEDGES) = INTERPOLATED_INDEX
             ENDIF
         IJKE(6,NEDGES) = IOR
         IF (.NOT.OB%SAWTOOTH) IJKE(6,NEDGES) = 0
         IEA(IC ,12) = NEDGES
         IEA(IC1,11) = NEDGES
         IEA(IC2, 9) = NEDGES
         IEA(ICG,10) = NEDGES
         IF (NOM.GT.0) THEN
            IJKE( 7,NEDGES) = NOM
            IJKE( 8,NEDGES) = IJKW(10,IW)
            IJKE( 9,NEDGES) = IJKW(11,IW)
            IJKE(10,NEDGES) = IJKW(12,IW)
            ENDIF
         ENDIF
C
      CASE(-2)
      SOLID_ONLY = .FALSE.
      CALL GET_OBST(NM,II,JJ-1,KK,IOR,SOLID_ONLY,NN)
      OB=>OBSTRUCTION(NN)
      IC1 = ICA(II,JJ,KK+1)
      IC2 = ICA(IIG,JJG,KKG+1)
      IW1 = IWA(ICG,-IOR)
      IW2 = IWA(IC2,-IOR)
      NOM = IJKW(9,IW)
      IF (IW1.GT.0 .AND. IW2.GT.0) THEN
         NEDGES = NEDGES + 1
         IJKE(1,NEDGES) = II
         IJKE(2,NEDGES) = JJ-1
         IJKE(3,NEDGES) = KK
         IJKE(4,NEDGES) = 1
         IJKE(5,NEDGES) = IBC
         IF (IJKW(5,IW1).EQ.INTERPOLATED_INDEX .OR.
     .       IJKW(5,IW2).EQ.INTERPOLATED_INDEX) THEN
             IF (NOM.GT.0) IJKE(5,NEDGES) = INTERPOLATED_INDEX
             ENDIF
         IJKE(6,NEDGES) = IOR
         IF (.NOT.OB%SAWTOOTH) IJKE(6,NEDGES) = 0
         IEA(IC ,3) = NEDGES
         IEA(IC1,1) = NEDGES
         IEA(IC2,2) = NEDGES
         IEA(ICG,4) = NEDGES
         IF (NOM.GT.0) THEN
            IJKE( 7,NEDGES) = NOM
            IJKE( 8,NEDGES) = IJKW(10,IW)
            IJKE( 9,NEDGES) = IJKW(11,IW)
            IJKE(10,NEDGES) = IJKW(12,IW)
            ENDIF
         ENDIF
      IC1 = ICA(II+1,JJ,KK)
      IC2 = ICA(IIG+1,JJG,KKG)
      IW1 = IWA(ICG,-IOR)
      IW2 = IWA(IC2,-IOR)
      NOM = IJKW(9,IW)
      IF (IW1.GT.0 .AND. IW2.GT.0) THEN
         NEDGES = NEDGES + 1
         IJKE(1,NEDGES) = II
         IJKE(2,NEDGES) = JJ-1
         IJKE(3,NEDGES) = KK
         IJKE(4,NEDGES) = 3
         IJKE(5,NEDGES) = IBC
         IF (IJKW(5,IW1).EQ.INTERPOLATED_INDEX .OR.
     .       IJKW(5,IW2).EQ.INTERPOLATED_INDEX) THEN
             IF (NOM.GT.0) IJKE(5,NEDGES) = INTERPOLATED_INDEX
             ENDIF
         IJKE(6,NEDGES) = IOR
         IF (.NOT.OB%SAWTOOTH) IJKE(6,NEDGES) = 0
         IEA(IC ,10) = NEDGES
         IEA(IC1, 9) = NEDGES
         IEA(IC2,11) = NEDGES
         IEA(ICG,12) = NEDGES
         IF (NOM.GT.0) THEN
            IJKE( 7,NEDGES) = NOM
            IJKE( 8,NEDGES) = IJKW(10,IW)
            IJKE( 9,NEDGES) = IJKW(11,IW)
            IJKE(10,NEDGES) = IJKW(12,IW)
            ENDIF
         ENDIF
C
      CASE( 3)
      SOLID_ONLY = .FALSE.
      CALL GET_OBST(NM,II,JJ,KK,IOR,SOLID_ONLY,NN)
      OB=>OBSTRUCTION(NN)
      IC1 = ICA(II+1,JJ,KK)
      IC2 = ICA(IIG+1,JJG,KKG)
      IW1 = IWA(ICG,-IOR)
      IW2 = IWA(IC2,-IOR)
      NOM = IJKW(9,IW)
      IF (IW1.GT.0 .AND. IW2.GT.0) THEN
         NEDGES = NEDGES + 1
         IJKE(1,NEDGES) = II
         IJKE(2,NEDGES) = JJ
         IJKE(3,NEDGES) = KK
         IJKE(4,NEDGES) = 2
         IJKE(5,NEDGES) = IBC
         IF (IJKW(5,IW1).EQ.INTERPOLATED_INDEX .OR.
     .       IJKW(5,IW2).EQ.INTERPOLATED_INDEX) THEN
             IF (NOM.GT.0) IJKE(5,NEDGES) = INTERPOLATED_INDEX
             ENDIF
         IJKE(6,NEDGES) = IOR
         IF (.NOT.OB%SAWTOOTH) IJKE(6,NEDGES) = 0
         IEA(IC ,8) = NEDGES
         IEA(IC1,7) = NEDGES
         IEA(IC2,5) = NEDGES
         IEA(ICG,6) = NEDGES
         IF (NOM.GT.0) THEN
            IJKE( 7,NEDGES) = NOM
            IJKE( 8,NEDGES) = IJKW(10,IW)
            IJKE( 9,NEDGES) = IJKW(11,IW)
            IJKE(10,NEDGES) = IJKW(12,IW)
            ENDIF
         ENDIF
      IC1 = ICA(II,JJ+1,KK)
      IC2 = ICA(IIG,JJG+1,KKG)
      IW1 = IWA(ICG,-IOR)
      IW2 = IWA(IC2,-IOR)
      NOM = IJKW(9,IW)
      IF (IW1.GT.0 .AND. IW2.GT.0) THEN
         NEDGES = NEDGES + 1
         IJKE(1,NEDGES) = II
         IJKE(2,NEDGES) = JJ
         IJKE(3,NEDGES) = KK
         IJKE(4,NEDGES) = 1
         IJKE(5,NEDGES) = IBC
         IF (IJKW(5,IW1).EQ.INTERPOLATED_INDEX .OR.
     .       IJKW(5,IW2).EQ.INTERPOLATED_INDEX) THEN
             IF (NOM.GT.0) IJKE(5,NEDGES) = INTERPOLATED_INDEX
             ENDIF
         IJKE(6,NEDGES) = IOR
         IF (.NOT.OB%SAWTOOTH) IJKE(6,NEDGES) = 0
         IEA(IC ,4) = NEDGES
         IEA(IC1,3) = NEDGES
         IEA(IC2,1) = NEDGES
         IEA(ICG,2) = NEDGES
         IF (NOM.GT.0) THEN
            IJKE( 7,NEDGES) = NOM
            IJKE( 8,NEDGES) = IJKW(10,IW)
            IJKE( 9,NEDGES) = IJKW(11,IW)
            IJKE(10,NEDGES) = IJKW(12,IW)
            ENDIF
         ENDIF
C
      CASE(-3)
      SOLID_ONLY = .FALSE.
      CALL GET_OBST(NM,II,JJ,KK-1,IOR,SOLID_ONLY,NN)
      OB=>OBSTRUCTION(NN)
      IC1 = ICA(II+1,JJ,KK)
      IC2 = ICA(IIG+1,JJG,KKG)
      IW1 = IWA(ICG,-IOR)
      IW2 = IWA(IC2,-IOR)
      NOM = IJKW(9,IW)
      IF (IW1.GT.0 .AND. IW2.GT.0) THEN
         NEDGES = NEDGES + 1
         IJKE(1,NEDGES) = II
         IJKE(2,NEDGES) = JJ
         IJKE(3,NEDGES) = KK-1
         IJKE(4,NEDGES) = 2
         IJKE(5,NEDGES) = IBC
         IF (IJKW(5,IW1).EQ.INTERPOLATED_INDEX .OR.
     .       IJKW(5,IW2).EQ.INTERPOLATED_INDEX) THEN
             IF (NOM.GT.0) IJKE(5,NEDGES) = INTERPOLATED_INDEX
             ENDIF
         IJKE(6,NEDGES) = IOR
         IF (.NOT.OB%SAWTOOTH) IJKE(6,NEDGES) = 0
         IEA(IC ,6) = NEDGES
         IEA(IC1,5) = NEDGES
         IEA(IC2,7) = NEDGES
         IEA(ICG,8) = NEDGES
         IF (NOM.GT.0) THEN
            IJKE( 7,NEDGES) = NOM
            IJKE( 8,NEDGES) = IJKW(10,IW)
            IJKE( 9,NEDGES) = IJKW(11,IW)
            IJKE(10,NEDGES) = IJKW(12,IW)
            ENDIF
         ENDIF
      IC1 = ICA(II,JJ+1,KK)
      IC2 = ICA(IIG,JJG+1,KKG)
      IW1 = IWA(ICG,-IOR)
      IW2 = IWA(IC2,-IOR)
      NOM = IJKW(9,IW)
      IF (IW1.GT.0 .AND. IW2.GT.0) THEN
         NEDGES = NEDGES + 1
         IJKE(1,NEDGES) = II
         IJKE(2,NEDGES) = JJ
         IJKE(3,NEDGES) = KK-1
         IJKE(4,NEDGES) = 1
         IJKE(5,NEDGES) = IBC
         IF (IJKW(5,IW1).EQ.INTERPOLATED_INDEX .OR.
     .       IJKW(5,IW2).EQ.INTERPOLATED_INDEX) THEN
             IF (NOM.GT.0) IJKE(5,NEDGES) = INTERPOLATED_INDEX
             ENDIF
         IJKE(6,NEDGES) = IOR
         IF (.NOT.OB%SAWTOOTH) IJKE(6,NEDGES) = 0
         IEA(IC ,2) = NEDGES
         IEA(IC1,1) = NEDGES
         IEA(IC2,3) = NEDGES
         IEA(ICG,4) = NEDGES
         IF (NOM.GT.0) THEN
            IJKE( 7,NEDGES) = NOM
            IJKE( 8,NEDGES) = IJKW(10,IW)
            IJKE( 9,NEDGES) = IJKW(11,IW)
            IJKE(10,NEDGES) = IJKW(12,IW)
            ENDIF
         ENDIF
C
      END SELECT
C
      ENDDO WALL_CELL_LOOP
C
C Outer edges of domain
C
      DO I=1,IBAR
      NEDGES = NEDGES + 1
      IJKE(1,NEDGES) = I
      IJKE(2,NEDGES) = 0
      IJKE(3,NEDGES) = 0
      IJKE(4,NEDGES) = 1
      IJKE(5,NEDGES) = OPEN_INDEX
      IJKE(6,NEDGES) = 0
      IEA(ICA(I,1,1),1) = NEDGES
      ENDDO
C
      DO I=1,IBAR
      NEDGES = NEDGES + 1
      IJKE(1,NEDGES) = I
      IJKE(2,NEDGES) = JBAR
      IJKE(3,NEDGES) = 0
      IJKE(4,NEDGES) = 1
      IJKE(5,NEDGES) = OPEN_INDEX
      IJKE(6,NEDGES) = 0
      IEA(ICA(I,JBAR,1),2) = NEDGES
      ENDDO
C
      DO I=1,IBAR
      NEDGES = NEDGES + 1
      IJKE(1,NEDGES) = I
      IJKE(2,NEDGES) = 0
      IJKE(3,NEDGES) = KBAR
      IJKE(4,NEDGES) = 1
      IJKE(5,NEDGES) = OPEN_INDEX
      IJKE(6,NEDGES) = 0
      IEA(ICA(I,1,KBAR),3) = NEDGES
      ENDDO
C
      DO I=1,IBAR
      NEDGES = NEDGES + 1
      IJKE(1,NEDGES) = I
      IJKE(2,NEDGES) = JBAR
      IJKE(3,NEDGES) = KBAR
      IJKE(4,NEDGES) = 1
      IJKE(5,NEDGES) = OPEN_INDEX
      IJKE(6,NEDGES) = 0
      IEA(ICA(I,JBAR,KBAR),4) = NEDGES
      ENDDO
C
C
      DO J=1,JBAR
      NEDGES = NEDGES + 1
      IJKE(1,NEDGES) = 0
      IJKE(2,NEDGES) = J
      IJKE(3,NEDGES) = 0
      IJKE(4,NEDGES) = 2
      IJKE(5,NEDGES) = OPEN_INDEX
      IJKE(6,NEDGES) = 0
      IEA(ICA(1,J,1),5) = NEDGES
      ENDDO
C
      DO J=1,JBAR
      NEDGES = NEDGES + 1
      IJKE(1,NEDGES) = 0
      IJKE(2,NEDGES) = J
      IJKE(3,NEDGES) = KBAR
      IJKE(4,NEDGES) = 2
      IJKE(5,NEDGES) = OPEN_INDEX
      IJKE(6,NEDGES) = 0
      IEA(ICA(1,J,KBAR),7) = NEDGES
      ENDDO
C
      DO J=1,JBAR
      NEDGES = NEDGES + 1
      IJKE(1,NEDGES) = IBAR
      IJKE(2,NEDGES) = J
      IJKE(3,NEDGES) = 0
      IJKE(4,NEDGES) = 2
      IJKE(5,NEDGES) = OPEN_INDEX
      IJKE(6,NEDGES) = 0
      IEA(ICA(IBAR,J,1),6)=NEDGES
      ENDDO
C
      DO J=1,JBAR
      NEDGES = NEDGES + 1
      IJKE(1,NEDGES) = IBAR
      IJKE(2,NEDGES) = J
      IJKE(3,NEDGES) = KBAR
      IJKE(4,NEDGES) = 2
      IJKE(5,NEDGES) = OPEN_INDEX
      IJKE(6,NEDGES) = 0
      IEA(ICA(IBAR,J,KBAR),8)=NEDGES
      ENDDO
C
      DO K=1,KBAR
      NEDGES = NEDGES + 1
      IJKE(1,NEDGES) = 0
      IJKE(2,NEDGES) = 0
      IJKE(3,NEDGES) = K
      IJKE(4,NEDGES) = 3
      IJKE(5,NEDGES) = OPEN_INDEX
      IJKE(6,NEDGES) = 0
      IEA(ICA(1,1,K),9)=NEDGES
      ENDDO
C
      DO K=1,KBAR
      NEDGES = NEDGES + 1
      IJKE(1,NEDGES) = IBAR
      IJKE(2,NEDGES) = 0
      IJKE(3,NEDGES) = K
      IJKE(4,NEDGES) = 3
      IJKE(5,NEDGES) = OPEN_INDEX
      IJKE(6,NEDGES) = 0
      IEA(ICA(IBAR,1,K),10)=NEDGES
      ENDDO
C
      DO K=1,KBAR
      NEDGES = NEDGES + 1
      IJKE(1,NEDGES) = 0
      IJKE(2,NEDGES) = JBAR
      IJKE(3,NEDGES) = K
      IJKE(4,NEDGES) = 3
      IJKE(5,NEDGES) = OPEN_INDEX
      IJKE(6,NEDGES) = 0
      IEA(ICA(1,JBAR,K),11)=NEDGES
      ENDDO
C
      DO K=1,KBAR
      NEDGES = NEDGES + 1
      IJKE(1,NEDGES) = IBAR
      IJKE(2,NEDGES) = JBAR
      IJKE(3,NEDGES) = K
      IJKE(4,NEDGES) = 3
      IJKE(5,NEDGES) = OPEN_INDEX
      IJKE(6,NEDGES) = 0
      IEA(ICA(IBAR,JBAR,K),12)=NEDGES
      ENDDO
C
      IF (EVACUATION_ONLY(NM)) IJKE(6,:) = 0
C
      END SUBROUTINE INITIALIZE_EDGES
C
C
      SUBROUTINE INITIALIZE_POISSON_SOLVER
C
      REAL(EB) XLM,XMU
      INTEGER N,IERR
      INTEGER, POINTER :: ITRN,JTRN,KTRN,LBC,MBC,NBC
      INTEGER, POINTER, DIMENSION(:) :: NOC
      TYPE (VENTS_TYPE), POINTER :: VT
C
C Allocate major arrays
C
      ITRN =>M%ITRN ; JTRN =>M%JTRN ; KTRN =>M%KTRN
       LBC =>M%LBC  ;  MBC =>M%MBC  ;  NBC =>M%NBC
      NOC=>TRANS(NM)%NOC
C
      IF (NOC(1).EQ.0 .AND. NOC(2).EQ.0 .AND. NOC(3).EQ.0) M%IPS=0
      IF (NOC(1).NE.0 .AND. NOC(2).EQ.0 .AND. NOC(3).EQ.0) M%IPS=1
      IF (NOC(1).EQ.0 .AND. NOC(2).NE.0 .AND. NOC(3).EQ.0) M%IPS=2
      IF (NOC(1).EQ.0 .AND. NOC(2).EQ.0 .AND. NOC(3).NE.0) M%IPS=3
      IF (NOC(1).NE.0 .AND. NOC(2).NE.0 .AND. NOC(3).EQ.0) M%IPS=4
      IF (NOC(1).NE.0 .AND. NOC(2).EQ.0 .AND. NOC(3).NE.0) M%IPS=5
      IF (NOC(1).EQ.0 .AND. NOC(2).NE.0 .AND. NOC(3).NE.0) M%IPS=6
      IF (EVACUATION_ONLY(NM)                            ) M%IPS=7
      IF (NOC(1).NE.0 .AND. NOC(2).NE.0 .AND. NOC(3).NE.0) 
     .  CALL SHUTDOWN('ERROR: Stretch at most 2 coordinate directions')
C
      IF (M%IPS.LE.1 .OR. M%IPS.EQ.4) THEN
         ITRN = IBP1
         IF (JBAR.GT.1) JTRN = JBP1
         IF (JBAR.EQ.1) JTRN = 1
         KTRN = KBP1
         ENDIF
C
      IF (M%IPS.EQ.2) THEN
         ITRN = JBP1
         JTRN = IBP1
         KTRN = KBP1
         ALLOCATE(M%BZST(JBP1,IBP1),STAT=IZERO)
         CALL ChkMemErr('INIT','BZST',IZERO)
         ALLOCATE(M%BZFT(JBP1,IBP1),STAT=IZERO)
         CALL ChkMemErr('INIT','BZFT',IZERO)
         ENDIF
C
      IF (M%IPS.EQ.3 .OR. M%IPS.EQ.6) THEN
         ITRN = KBP1
         IF (JBAR.GT.1) JTRN = JBP1
         IF (JBAR.EQ.1) JTRN = 1
         KTRN = IBP1
         ALLOCATE(M%BXST(KBP1,JTRN),STAT=IZERO)
         CALL ChkMemErr('INIT','BXST',IZERO)
         ALLOCATE(M%BXFT(KBP1,JTRN),STAT=IZERO)
         CALL ChkMemErr('INIT','BXFT',IZERO)
         ALLOCATE(M%BYST(KBP1,IBP1),STAT=IZERO)
         CALL ChkMemErr('INIT','BYST',IZERO)
         ALLOCATE(M%BYFT(KBP1,IBP1),STAT=IZERO)
         CALL ChkMemErr('INIT','BYFT',IZERO)
         ALLOCATE(M%BZST(JTRN,IBP1),STAT=IZERO)
         CALL ChkMemErr('INIT','BZST',IZERO)
         ALLOCATE(M%BZFT(JTRN,IBP1),STAT=IZERO)
         CALL ChkMemErr('INIT','BZFT',IZERO)
         ENDIF
C
      IF (M%IPS.EQ.5) THEN
         ITRN = IBP1
         JTRN = KBP1
         KTRN = JBP1
         ALLOCATE(M%BXST(KBP1,JBP1),STAT=IZERO)
         CALL ChkMemErr('INIT','BXST',IZERO)
         ALLOCATE(M%BXFT(KBP1,JBP1),STAT=IZERO)
         CALL ChkMemErr('INIT','BXFT',IZERO)
         ENDIF
C
         IF (M%IPS.EQ.7) THEN
         ITRN = IBP1
         JTRN = JBP1
         KTRN = 1
         ENDIF
C
      IF (M%IPS.LE.3 .OR. M%IPS.EQ.7) THEN
         M%LSAVE = (ITRN+1)*JTRN*KTRN+7*ITRN+5*JTRN+6*KTRN+56
         M%LWORK = (ITRN+1)*JTRN*KTRN
         ELSE
         N_LOOP: DO N=1,50
         IF ((JTRN+1).LE.2**N) EXIT N_LOOP
         ENDDO N_LOOP
         M%LSAVE = KTRN*(6*N*(2**N)+2*N+19)+8*ITRN+7*JTRN+38
         M%LWORK = JTRN*(ITRN*(KTRN+1)+1)
         ENDIF
C
      ALLOCATE(M%SAVE(-3:M%LSAVE),STAT=IZERO)
      CALL ChkMemErr('INIT','SAVE',IZERO)
      ALLOCATE(M%WORK(M%LWORK),STAT=IZERO)
      CALL ChkMemErr('INIT','WORK',IZERO)
      ALLOCATE(M%PRHS(ITRN,JTRN,KTRN),STAT=IZERO)
      CALL ChkMemErr('INIT','PRHS',IZERO)
      IF (KBAR.GT.1) THEN
      IF (JBAR.GT.1) ALLOCATE(M%BXS(JBP1,KBP1),STAT=IZERO)
      IF (JBAR.EQ.1) ALLOCATE(M%BXS(1,KBP1)   ,STAT=IZERO)
      ELSE
                     ALLOCATE(M%BXS(JBP1,1)   ,STAT=IZERO)
      ENDIF
      CALL ChkMemErr('INIT','BXS',IZERO)
      IF (KBAR.GT.1) THEN
      IF (JBAR.GT.1) ALLOCATE(M%BXF(JBP1,KBP1),STAT=IZERO)
      IF (JBAR.EQ.1) ALLOCATE(M%BXF(1,KBP1)   ,STAT=IZERO)
      ELSE
                     ALLOCATE(M%BXF(JBP1,1)   ,STAT=IZERO)
      ENDIF
      CALL ChkMemErr('INIT','BXF',IZERO)
      IF (KBAR.GT.1) THEN
      ALLOCATE(M%BYS(IBP1,KBP1),STAT=IZERO)
      ELSE
      ALLOCATE(M%BYS(IBP1,1),STAT=IZERO)
      ENDIF
      CALL ChkMemErr('INIT','BYS',IZERO)
      IF (KBAR.GT.1) THEN
      ALLOCATE(M%BYF(IBP1,KBP1),STAT=IZERO)
      ELSE
      ALLOCATE(M%BYF(IBP1,1),STAT=IZERO)
      ENDIF
      CALL ChkMemErr('INIT','BYF',IZERO)
      IF (JBAR.GT.1) ALLOCATE(M%BZS(IBP1,JBP1),STAT=IZERO)
      IF (JBAR.EQ.1) ALLOCATE(M%BZS(IBP1,1)   ,STAT=IZERO)
      CALL ChkMemErr('INIT','BZS',IZERO)
      IF (JBAR.GT.1) ALLOCATE(M%BZF(IBP1,JBP1),STAT=IZERO)
      IF (JBAR.EQ.1) ALLOCATE(M%BZF(IBP1,1)   ,STAT=IZERO)
      CALL ChkMemErr('INIT','BZF',IZERO)
C
      M%SAVE  = 0.
      M%WORK  = 0.
      M%PRHS  = 0.
      M%BXS   = 0.
      M%BXF   = 0.
      M%BYS   = 0.
      M%BYF   = 0.
      M%BZS   = 0.
      M%BZF   = 0.
C
C Initialize pressure solver   
C
      XLM = 0.         ! No Helmholtz equation
      XMU = 0.         ! No Helmholtz equation
      LBC = 3
      MBC = 3
      NBC = 3
C
      NVLOOP: DO N=1,M%NV
C
      VT => M%VENTS(N)
C
      IF (VT%INDEX.NE.2) CYCLE NVLOOP
C
         IF (VT%I1.EQ.0 .AND. VT%I2.EQ.0) THEN
            IF (LBC.EQ.3) LBC = 2
            IF (LBC.EQ.4) LBC = 1
            ENDIF
         IF (VT%I1.EQ.M%IBAR .AND. VT%I2.EQ.M%IBAR) THEN
            IF (LBC.EQ.3) LBC = 4
            IF (LBC.EQ.2) LBC = 1
            ENDIF
C
         IF (VT%J1.EQ.0 .AND. VT%J2.EQ.0) THEN
            IF (MBC.EQ.3) MBC = 2
            IF (MBC.EQ.4) MBC = 1
            ENDIF
         IF (VT%J1.EQ.M%JBAR .AND. VT%J2.EQ.M%JBAR) THEN
            IF (MBC.EQ.3) MBC = 4
            IF (MBC.EQ.2) MBC = 1
            ENDIF
C
         IF (VT%K1.EQ.0 .AND. VT%K2.EQ.0) THEN
            IF (NBC.EQ.3) NBC = 2
            IF (NBC.EQ.4) NBC = 1
            ENDIF
         IF (VT%K1.EQ.M%KBAR .AND. VT%K2.EQ.M%KBAR) THEN
            IF (NBC.EQ.3) NBC = 4
            IF (NBC.EQ.2) NBC = 1
            ENDIF
C
      ENDDO NVLOOP
C
      DO IW=1,M%NEWC
      IF (M%IJKW(9,IW).GT.0) THEN
         SELECT CASE(M%IJKW(4,IW))
         CASE( 1)
         IF (LBC.EQ.3) LBC = 2
         IF (LBC.EQ.4) LBC = 1
         CASE(-1)
         IF (LBC.EQ.3) LBC = 4
         IF (LBC.EQ.2) LBC = 1
         CASE( 2)
         IF (MBC.EQ.3) MBC = 2
         IF (MBC.EQ.4) MBC = 1
         CASE(-2)
         IF (MBC.EQ.3) MBC = 4
         IF (MBC.EQ.2) MBC = 1
         CASE( 3)
         IF (NBC.EQ.3) NBC = 2
         IF (NBC.EQ.4) NBC = 1
         CASE(-3)
         IF (NBC.EQ.3) NBC = 4
         IF (NBC.EQ.2) NBC = 1
         END SELECT
      ENDIF
      ENDDO
C
C User over-rides of Poisson boundary conditions
C
      IF (PBC(1,NM).EQ.0 .AND. PBC(2,NM).EQ.0) LBC = 1
      IF (PBC(1,NM).EQ.0 .AND. PBC(2,NM).EQ.1) LBC = 2
      IF (PBC(1,NM).EQ.1 .AND. PBC(2,NM).EQ.1) LBC = 3
      IF (PBC(1,NM).EQ.1 .AND. PBC(2,NM).EQ.0) LBC = 4
      IF (PBC(3,NM).EQ.0 .AND. PBC(4,NM).EQ.0) MBC = 1
      IF (PBC(3,NM).EQ.0 .AND. PBC(4,NM).EQ.1) MBC = 2
      IF (PBC(3,NM).EQ.1 .AND. PBC(4,NM).EQ.1) MBC = 3
      IF (PBC(3,NM).EQ.1 .AND. PBC(4,NM).EQ.0) MBC = 4
      IF (PBC(5,NM).EQ.0 .AND. PBC(6,NM).EQ.0) NBC = 1
      IF (PBC(5,NM).EQ.0 .AND. PBC(6,NM).EQ.1) NBC = 2
      IF (PBC(5,NM).EQ.1 .AND. PBC(6,NM).EQ.1) NBC = 3
      IF (PBC(5,NM).EQ.1 .AND. PBC(6,NM).EQ.0) NBC = 4
C
C Poisson solver with stretching in the 1st coordinate
C
      IF (M%IPS.LE.1)  THEN
      IF (JBAR.GT.1) THEN
      CALL H3CZIS(XS,XF,IBAR,LBC,YS,YF,JBAR,MBC,ZS,ZF,KBAR,NBC,M%HX,
     .            XLM,ITRN,JTRN,IERR,M%SAVE)
      ELSE
      IF (.NOT.CYLINDRICAL)
     .CALL H2CZIS(XS,XF,IBAR,LBC,ZS,ZF,KBAR,NBC,M%HX,XLM,
     .            ITRN,IERR,M%SAVE)
      IF (CYLINDRICAL) THEN
      IF (XS.EQ.0. .AND. LBC.EQ.1) LBC = 5
      IF (XS.EQ.0. .AND. LBC.EQ.2) LBC = 6
      IF (XS.EQ.0. .AND. LBC.EQ.3) LBC = 6
      IF (XS.EQ.0. .AND. LBC.EQ.4) LBC = 5
      CALL H2CYIS(XS,XF,IBAR,LBC,ZS,ZF,KBAR,NBC,XLM,XMU,ITRN,
     .            IERR,M%SAVE)
      ENDIF
      ENDIF
      ENDIF
      IF (M%IPS.EQ.2) 
     .CALL H3CZIS(YS,YF,JBAR,MBC,XS,XF,IBAR,LBC,ZS,ZF,KBAR,NBC,M%HY,
     .            XLM,ITRN,JTRN,IERR,M%SAVE)
      IF (M%IPS.EQ.3) THEN
      IF (.NOT.TWO_D)
     .CALL H3CZIS(ZS,ZF,KBAR,NBC,YS,YF,JBAR,MBC,XS,XF,IBAR,LBC,M%HZ,
     .            XLM,ITRN,JTRN,IERR,M%SAVE)
      IF (TWO_D)
     .CALL H2CZIS(ZS,ZF,KBAR,NBC,XS,XF,IBAR,LBC,M%HZ,XLM,
     .            ITRN,IERR,M%SAVE)
      ENDIF
C
C Poisson solver with stretching in the 1st and 2nd coordinate
C
      IF (M%IPS.EQ.4) 
     .CALL H3CSIS(XS,XF,IBAR,LBC,YS,YF,JBAR,MBC,ZS,ZF,KBAR,NBC,
     .            XLM,ITRN,JTRN,IERR,M%SAVE,M%WORK,M%HX,M%HY)
      IF (M%IPS.EQ.5) THEN
      IF (.NOT.TWO_D)
     .CALL H3CSIS(XS,XF,IBAR,LBC,ZS,ZF,KBAR,NBC,YS,YF,JBAR,MBC,
     .            XLM,ITRN,JTRN,IERR,M%SAVE,M%WORK,M%HX,M%HZ)
      IF (TWO_D)
     .CALL H2CZIS(ZS,ZF,KBAR,NBC,XS,XF,IBAR,LBC,M%HZ,XLM,
     .            ITRN,IERR,M%SAVE)
      ENDIF
      IF (M%IPS.EQ.6) 
     .CALL H3CSIS(ZS,ZF,KBAR,NBC,YS,YF,JBAR,MBC,XS,XF,IBAR,LBC,
     .            XLM,ITRN,JTRN,IERR,M%SAVE,M%WORK,M%HZ,M%HY)
      IF (M%IPS.EQ.7)
     .CALL H2CZIS(XS,XF,IBAR,LBC,YS,YF,JBAR,MBC,M%HX,XLM,
     .            ITRN,IERR,M%SAVE)
C
C Check for errors with Poisson solver initialization
C
      IF (IERR.NE.0) THEN
         WRITE(MESSAGE,'(A,I2,A,I3)') 
     .   'ERROR: Poisson initialization error, Number=',IERR,
     .   ', Mesh=',NM
         CALL SHUTDOWN(MESSAGE)
         ENDIF
C
      END SUBROUTINE INITIALIZE_POISSON_SOLVER
C
C
      SUBROUTINE INITIAL_NOISE
C
C Generate random noise at the start of the simulation
C
      REAL(EB) VFAC,RN
      INTEGER I,J,K
C
      IF (EVACUATION_ONLY(NM)) RETURN
C
      VFAC = 0.005
C
      DO K=0,M%KBAR
      DO J=0,M%JBAR
      LOOP_1: DO I=1,M%IBAR
      IF (M%SOLID(M%ICA(I,J,K))   .OR. 
     .    M%SOLID(M%ICA(I,J,K+1)) .OR.
     .    M%SOLID(M%ICA(I,J+1,K)) .OR. 
     .    M%SOLID(M%ICA(I,J+1,K+1))) 
     .    CYCLE LOOP_1
      CALL RANDOM_NUMBER(RN)
      RN = VFAC*(-1. + 2.*RN)*M%DXMIN
      M%W(I,J,K)   = M%W(I,J,K)   - RN*M%RDY(J)
      M%W(I,J+1,K) = M%W(I,J+1,K) + RN*M%RDY(J+1)
      M%V(I,J,K)   = M%V(I,J,K)   + RN*M%RDZ(K)
      M%V(I,J,K+1) = M%V(I,J,K+1) - RN*M%RDZ(K+1)
      ENDDO LOOP_1
      ENDDO
      ENDDO
C
      DO K=0,M%KBAR
      DO J=1,M%JBAR
      LOOP_2: DO I=0,M%IBAR
      IF (M%SOLID(M%ICA(I,J,K))   .OR. 
     .    M%SOLID(M%ICA(I,J,K+1)) .OR.
     .    M%SOLID(M%ICA(I+1,J,K)) .OR. 
     .    M%SOLID(M%ICA(I+1,J,K+1))) 
     .    CYCLE LOOP_2
      CALL RANDOM_NUMBER(RN)
      RN = VFAC*(-1. + 2.*RN)*M%DXMIN
      M%W(I,J,K)   = M%W(I,J,K)   - RN*M%RDX(I)*M%R(I)*M%RRN(I)
      M%W(I+1,J,K) = M%W(I+1,J,K) + RN*M%RDX(I+1)*M%R(I)*M%RRN(I+1)
      M%U(I,J,K)   = M%U(I,J,K)   + RN*M%RDZ(K)
      M%U(I,J,K+1) = M%U(I,J,K+1) - RN*M%RDZ(K+1)
      ENDDO LOOP_2
      ENDDO
      ENDDO
C
      DO K=1,M%KBAR
      DO J=0,M%JBAR
      LOOP_3: DO I=0,M%IBAR
      IF (M%SOLID(M%ICA(I,J,K))   .OR. 
     .    M%SOLID(M%ICA(I,J+1,K)) .OR.
     .    M%SOLID(M%ICA(I+1,J,K)) .OR. 
     .    M%SOLID(M%ICA(I+1,J+1,K))) 
     .    CYCLE LOOP_3
      CALL RANDOM_NUMBER(RN)
      RN = VFAC*(-1. + 2.*RN)*M%DXMIN
      M%V(I,J,K)   = M%V(I,J,K)   - RN*M%RDX(I)
      M%V(I+1,J,K) = M%V(I+1,J,K) + RN*M%RDX(I+1)
      M%U(I,J,K)   = M%U(I,J,K)   + RN*M%RDY(J)
      M%U(I,J+1,K) = M%U(I,J+1,K) - RN*M%RDY(J+1)
      ENDDO LOOP_3
      ENDDO
      ENDDO
C
      END SUBROUTINE INITIAL_NOISE
C
C
      SUBROUTINE INITIALIZE_THCP
C
      INTEGER III
      TYPE (THERMOCOUPLE_TYPE), POINTER :: TC
C
      TC_LOOP: DO N=1,NTC
      TC => THERMOCOUPLE(N)
C
      IF (NM.NE.TC%MESH) CYCLE TC_LOOP
      IF (TC%INDEX.GT.0) CYCLE TC_LOOP
C
      II  = GINV(TC%X-M%XS,1,NM)*M%RDXI   + 1.
      JJ  = GINV(TC%Y-M%YS,2,NM)*M%RDETA  + 1.
      KK  = GINV(TC%Z-M%ZS,3,NM)*M%RDZETA + 1.
      IOR = TC%IOR
      IC  = M%ICA(II,JJ,KK)
C
      IF (M%SOLID(IC)) THEN
         SELECT CASE(IOR)
         CASE(-1) ; IF (II.GT.0)      II = II-1
         CASE( 1) ; IF (II.LT.M%IBP1) II = II+1
         CASE(-2) ; IF (JJ.GT.0)      JJ = JJ-1
         CASE( 2) ; IF (JJ.LT.M%JBP1) JJ = JJ+1
         CASE(-3) ; IF (KK.GT.0)      KK = KK-1
         CASE( 3) ; IF (KK.LT.M%KBP1) KK = KK+1
         END SELECT
      ENDIF
C
      IC  = M%ICA(II,JJ,KK)
      IW  = M%IWA(IC,-IOR)
C
      IF (IW.LE.0) THEN
         SELECT CASE(IOR)
         CASE(-1) ; IF (II.GT.0)      IC = M%ICA(II-1,JJ,KK)
         CASE( 1) ; IF (II.LT.M%IBP1) IC = M%ICA(II+1,JJ,KK)
         CASE(-2) ; IF (JJ.GT.0)      IC = M%ICA(II,JJ-1,KK)
         CASE( 2) ; IF (JJ.LT.M%JBP1) IC = M%ICA(II,JJ+1,KK)
         CASE(-3) ; IF (KK.GT.0)      IC = M%ICA(II,JJ,KK-1)
         CASE( 3) ; IF (KK.LT.M%KBP1) IC = M%ICA(II,JJ,KK+1)
         END SELECT
         IW = M%IWA(IC,-IOR)
      ENDIF
C
      IF (IW.GT.0) THEN
         TC%IW = IW
         IF (TC%INDEX.EQ.-6) THEN
         IBC = M%IJKW(5,IW)
         IF (MHC(IBC).NE.3) THEN
            TC%I_DEPTH = 1
            ELSE
            TC%I_DEPTH = NWP(IBC)
            DO III=NWP(IBC),1,-1
            IF (TC%DEPTH.LE.X_S(III,IBC)) TC%I_DEPTH = III
            ENDDO
            ENDIF
         ENDIF
      ELSE
         WRITE(MESSAGE,'(A,I4,A)') 'Reposition THCP No.',TC%ORDINAL,
     .   '. FDS cannot determine which boundary cell to assign'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
C
      ENDDO TC_LOOP
C
      END SUBROUTINE INITIALIZE_THCP
C
C
      SUBROUTINE INITIALIZE_INTERPOLATION
C
C Create arrays by which info is to exchanged across meshes
C
      INTEGER NOM,I,J,K
      TYPE (MESH_TYPE), POINTER :: M2
C
      IF (NM.EQ.1) RETURN
C
      ALLOCATE(M%INTERPOLATED_MESH(1:M%IBAR,1:M%JBAR,1:M%KBAR),
     .         STAT=IZERO)
      CALL ChkMemErr('INIT','INTERPOLATED_MESH',IZERO)  
      M%INTERPOLATED_MESH = 0
C
      DO K=1,M%KBAR
      DO J=1,M%JBAR
      DO I=1,M%IBAR
C
      OTHER_MESH_LOOP: DO NOM=1,NM-1
      M2=>MESH(NOM)
      IF (M%X(I-1).GE.M2%XS .AND. M%X(I).LE.M2%XF .AND.
     .    M%Y(J-1).GE.M2%YS .AND. M%Y(J).LE.M2%YF .AND.
     .    M%Z(K-1).GE.M2%ZS .AND. M%Z(K).LE.M2%ZF) THEN 
         M%INTERPOLATED_MESH(I,J,K) = NOM
         EXIT OTHER_MESH_LOOP
         ENDIF
      ENDDO OTHER_MESH_LOOP
C
      ENDDO
      ENDDO
      ENDDO
C
      END SUBROUTINE INITIALIZE_INTERPOLATION
C
C
      END SUBROUTINE INITIALIZE_MESH_VARIABLES
C
C
C
      SUBROUTINE INITIALIZE_GLOBAL_VARIABLES
C
C Initialize time, printout and plot clocks
C
      ALLOCATE(PINCLK(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','PINCLK',IZERO) 
      ALLOCATE(SPINCLK(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','SPINCLK',IZERO) 
      ALLOCATE(PARCLK(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','PARCLK',IZERO) 
      ALLOCATE(ISOCLK(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','ISOCLK',IZERO) 
      ALLOCATE(BFCLK(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','BFCLK',IZERO) 
      ALLOCATE(SFCLK(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','SFCLK',IZERO) 
      ALLOCATE(CORCLK(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','CORCLK',IZERO) 
      ALLOCATE(PLTCLK(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','PLTCLK',IZERO) 
C
      ICYC     = 0
      PARCLK   = 0.
      PINCLK   = 0.
      SPINCLK  = 0.
      SPRKCLK  = 0.      ; IF (NSPR.EQ.0) SPRKCLK = 1.E10
      HEATCLK  = 0.      ; IF (NHD .EQ.0) HEATCLK = 1.E10
      SMOKECLK = 0.      ; IF (NSD .EQ.0) SMOKECLK= 1.E10
      TCCLK    = 0.      ; IF (NTC .EQ.0) TCCLK   = 1.E10
      PLTCLK   = WPLT
      ISOCLK   = 0.      ; IF (NIF.EQ.0)  ISOCLK  = 1.E10
      SFCLK    = 0.   
      BFCLK    = 0.      ; IF (NBF.EQ.0)  BFCLK   = 1.E10
      CORCLK   = DTCORE
      HRRCLK   = 0.
      MINTCLK  = 0.      ; IF (NSPEC.EQ.0) MINTCLK = 1.E10
C
      ALLOCATE(HRR(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','HRR',IZERO) 
      ALLOCATE(RHRR(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','RHRR',IZERO) 
      ALLOCATE(CHRR(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','CHRR',IZERO) 
      ALLOCATE(FHRR(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','FHRR',IZERO) 
      ALLOCATE(MLR(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','MLR',IZERO) 
C
      ALLOCATE(HRR_SUM(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','HRR_SUM',IZERO)  ; HRR_SUM=0.
      ALLOCATE(RHRR_SUM(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','RHRR_SUM',IZERO) ; RHRR_SUM=0.
      ALLOCATE(CHRR_SUM(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','CHRR_SUM',IZERO) ; CHRR_SUM=0.
      ALLOCATE(FHRR_SUM(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','FHRR_SUM',IZERO) ; FHRR_SUM=0.
      ALLOCATE(MLR_SUM(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','MLR_SUM',IZERO)  ; MLR_SUM=0.
      ALLOCATE(HRR_COUNT(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','MLR_COUNT',IZERO)  ; HRR_COUNT=0.
C
      ALLOCATE(MINT(0:20,NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','MINT',IZERO) 
      ALLOCATE(MINT_SUM(0:20,NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','MINT_SUM',IZERO) ; MINT_SUM=0.
      ALLOCATE(MINT_COUNT(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','MINT_COUNT',IZERO) ; MINT_COUNT=0.
C
      ALLOCATE(I_MIN(NMESHES,NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','I_MIN',IZERO) ; I_MIN = -10
      ALLOCATE(I_MAX(NMESHES,NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','I_MAX',IZERO) ; I_MAX = -10
      ALLOCATE(J_MIN(NMESHES,NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','J_MIN',IZERO) ; J_MIN = -10
      ALLOCATE(J_MAX(NMESHES,NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','J_MAX',IZERO) ; J_MAX = -10
      ALLOCATE(K_MIN(NMESHES,NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','K_MIN',IZERO) ; K_MIN = -10
      ALLOCATE(K_MAX(NMESHES,NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','K_MAX',IZERO) ; K_MAX = -10
      ALLOCATE(NIC(NMESHES,NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','NIC',IZERO) ; NIC = 0
C
      ALLOCATE(T_PER_STEP(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','T_PER_STEP',IZERO) ; T_PER_STEP = 0.
      ALLOCATE(T_ACCUM(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','T_ACCUM',IZERO) ; T_ACCUM = 0.
      ALLOCATE(NTCYC(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','NTCYC',IZERO) ; NTCYC = 0
      ALLOCATE(NCYC(NMESHES),STAT=IZERO)
      CALL ChkMemErr('INIT','NCYC',IZERO) ; NCYC = 0
C
      END SUBROUTINE INITIALIZE_GLOBAL_VARIABLES
C
C
      SUBROUTINE INIT_WALL_CELL(NM,I,J,K,I_OBST,IW,IOR,IBC)
C
      INTEGER NM,NOM,ICO,IBC,IOR,NWP_MAX
      REAL(EB) PX,PY,PZ,X1,X2,Y1,Y2,Z1,Z2,T_ACTIVATE,XI,YJ,ZK,
     .         XIN,YIN,ZIN,DIST
      INTEGER N,NN,I_OBST,I,J,K,IBCX,IIG,JJG,KKG,IW,II,JJ,KK,IC,ICG
      TYPE (MESH_TYPE), POINTER :: M,MM
      TYPE (WALL_TYPE), POINTER :: WC
      TYPE (OBSTRUCTION_TYPE), POINTER :: OBX
      TYPE (VENTS_TYPE), POINTER :: VT
C
      M=>MESH(NM)
C
C Compute boundary cell physical coords (XW,YW,ZW) and area (AW)
C
      IBCX         = IBC
c     M%IJKW(9,IW) = NM
      M%IJKW(1,IW) = I
      M%IJKW(2,IW) = J
      M%IJKW(3,IW) = K
      M%IJKW(4,IW) = IOR
      M%IJKW(5,IW) = IBCX
      M%OBST_INDEX(IW) = I_OBST
C
      IF (ABS(IOR).EQ.1) THEN
         IF (IOR.EQ. 1) THEN
            M%XW(IW)     = M%X(I)
            M%IJKW(6,IW) = I+1
            M%RDN(IW)    = M%RDXN(I) 
            M%AW(IW)     = M%R(I)*M%DY(J)*M%DZ(K)
            M%UW(IW)     = -U0
            ENDIF
         IF (IOR.EQ.-1) THEN
            M%XW(IW)     = M%X(I-1)
            M%IJKW(6,IW) = I-1
            M%RDN(IW)    = M%RDXN(I-1) 
            M%AW(IW)     = M%R(I-1)*M%DY(J)*M%DZ(K)
            M%UW(IW)     = U0
            ENDIF
         M%IJKW(7,IW) = J
         M%IJKW(8,IW) = K
         M%YW(IW) = 0.5*(M%Y(J)+M%Y(J-1))  
         M%ZW(IW) = 0.5*(M%Z(K)+M%Z(K-1))
         ENDIF
      IF (ABS(IOR).EQ.2) THEN
         IF (IOR.EQ. 2) THEN
            M%YW(IW)     = M%Y(J)
            M%IJKW(7,IW) = J+1
            M%RDN(IW)    = M%RDYN(J) 
            M%UW(IW)     = -V0
            ENDIF
         IF (IOR.EQ.-2) THEN
            M%YW(IW)     = M%Y(J-1)
            M%IJKW(7,IW) = J-1
            M%RDN(IW)    = M%RDYN(J-1) 
            M%UW(IW)     = V0
            ENDIF
         M%IJKW(6,IW) = I
         M%IJKW(8,IW) = K
         M%XW(IW) = 0.5*(M%X(I)+M%X(I-1))
         M%ZW(IW) = 0.5*(M%Z(K)+M%Z(K-1))
         M%AW(IW) = M%DX(I)*M%DZ(K)
         ENDIF
      IF (ABS(IOR).EQ.3) THEN
         IF (IOR.EQ. 3) THEN
            M%ZW(IW)     = M%Z(K)
            M%IJKW(8,IW) = K+1
            M%RDN(IW)    = M%RDZN(K) 
            M%UW(IW)     = -W0
            ENDIF
         IF (IOR.EQ.-3) THEN
            M%ZW(IW)     = M%Z(K-1)
            M%IJKW(8,IW) = K-1
            M%RDN(IW)    = M%RDZN(K-1) 
            M%UW(IW)     = W0
            ENDIF
         M%IJKW(6,IW) = I
         M%IJKW(7,IW) = J
         M%XW(IW) = 0.5*(M%X(I)+M%X(I-1))
         M%YW(IW) = 0.5*(M%Y(J)+M%Y(J-1))
         M%AW(IW) = M%DX(I)*M%RC(I)*M%DY(J)
         ENDIF
C
C Assign codes for boundary functions
C
      IF (M%AW(IW).GT.0.) M%RAW(IW) = 1./M%AW(IW)
C
      IIG = M%IJKW(6,IW)
      JJG = M%IJKW(7,IW)
      KKG = M%IJKW(8,IW)
C
      M%IV(IW) = 1
C
      IC  = M%ICA(I,J,K)
      ICG = M%ICA(IIG,JJG,KKG)
      IF (M%SOLID(IC)) M%IWA(IC,IOR) = IW
      IF (M%SOLID(IC) .AND. M%SOLID(ICG)) M%IV(IW) = 0
      M%IWA(ICG,-IOR) = IW
C
cc    IF (M%OBSTRUCTION(I_OBST)%HIDDEN) M%IV(IW) = 0
C
C Check for neighboring meshes
C
      CHECK_MESHES: IF (IW.LE.M%NEWC .AND.
     .                  .NOT.EVACUATION_ONLY(NM)) THEN
C
      XIN = M%XW(IW) 
      YIN = M%YW(IW) 
      ZIN = M%ZW(IW) 
C
      FIND_MESH: DO NOM=1,NMESHES
      IF (NOM.EQ.NM) CYCLE FIND_MESH
      IF (EVACUATION_ONLY(NOM)) CYCLE FIND_MESH
      MM=>MESH(NOM)
      IF (XIN.LT.MM%XS) CYCLE FIND_MESH
      IF (XIN.GT.MM%XF) CYCLE FIND_MESH
      IF (YIN.LT.MM%YS) CYCLE FIND_MESH
      IF (YIN.GT.MM%YF) CYCLE FIND_MESH
      IF (ZIN.LT.MM%ZS) CYCLE FIND_MESH
      IF (ZIN.GT.MM%ZF) CYCLE FIND_MESH
      SELECT CASE(IOR)
      CASE( 1) ; XIN = XIN - 0.05*M%DX(0)
      CASE(-1) ; XIN = XIN + 0.05*M%DX(M%IBP1)
      CASE( 2) ; YIN = YIN - 0.05*M%DY(0)
      CASE(-2) ; YIN = YIN + 0.05*M%DY(M%JBP1)
      CASE( 3) ; ZIN = ZIN - 0.05*M%DZ(0)
      CASE(-3) ; ZIN = ZIN + 0.05*M%DZ(M%KBP1)
      END SELECT
      IF (XIN.LT.MM%XS) CYCLE FIND_MESH
      IF (XIN.GT.MM%XF) CYCLE FIND_MESH
      IF (YIN.LT.MM%YS) CYCLE FIND_MESH
      IF (YIN.GT.MM%YF) CYCLE FIND_MESH
      IF (ZIN.LT.MM%ZS) CYCLE FIND_MESH
      IF (ZIN.GT.MM%ZF) CYCLE FIND_MESH
      XI = MM%CELLSI(NINT((XIN-MM%XS)*MM%RDXINT)) + 1.
      YJ = MM%CELLSJ(NINT((YIN-MM%YS)*MM%RDYINT)) + 1.
      ZK = MM%CELLSK(NINT((ZIN-MM%ZS)*MM%RDZINT)) + 1.
      II = FLOOR(XI)
      JJ = FLOOR(YJ)
      KK = FLOOR(ZK)
      M%INTERPOLATION_FACTOR(IW,1) = XI-II
      M%INTERPOLATION_FACTOR(IW,2) = YJ-JJ
      M%INTERPOLATION_FACTOR(IW,3) = ZK-KK
      M%CELL_VOLUME_RATIO(IW) = M%DX(IIG)*M%RC(IIG)*M%DY(JJG)*M%DZ(KKG)/
     .                         MM%DX(II)*MM%RC(II)*MM%DY(JJ)*MM%DZ(KK)
      M%IJKW(9,IW)  = NOM
      M%IJKW(10,IW) = II 
      M%IJKW(11,IW) = JJ 
      M%IJKW(12,IW) = KK 
      IF (I_OBST.EQ.0 .AND. .NOT.M%SOLID(ICG)) THEN
cc    IF (I_OBST.EQ.0) THEN
         M%IV(IW) = 4
         IBCX     = INTERPOLATED_INDEX
         M%IJKW(5,IW) = IBCX
         ENDIF
      ICO = MM%ICA(II,JJ,KK)
      IF (.NOT.MM%SOLID(ICO) .AND. .NOT.M%SOLID(ICG)) 
cc    IF (.NOT.MM%SOLID(ICO)) 
     .   M%SOLID(M%ICA(I,J,K)) = .FALSE.
      EXIT FIND_MESH
      ENDDO FIND_MESH
C
      ENDIF CHECK_MESHES
C
C Assign internal values of temp, density, and mass fraction
C
      IF (NSPEC.GT.0) THEN
      M%RSUM_W(IW)= M%RSUM(IIG,JJG,KKG)
      M%RSUM(I,J,K) = M%RSUM(IIG,JJG,KKG)
      M%YY_W(IW,1:NSPEC)  = M%YY(IIG,JJG,KKG,1:NSPEC)
      M%YY(I,J,K,1:NSPEC) = M%YY(IIG,JJG,KKG,1:NSPEC)
      ENDIF
C
      M%RHO_W(IW) = M%RHO(IIG,JJG,KKG)
      M%TMP_W(IW) = M%TMP(IIG,JJG,KKG)
      M%TMP(I,J,K) = M%TMP(IIG,JJG,KKG)
      M%RHO(I,J,K) = M%RHO(IIG,JJG,KKG)
C
C Assign various other quantities to the cell
C
      IF (I_OBST.GT.0) THEN
         OBX=>M%OBSTRUCTION(I_OBST)
         M%AREA_ADJUST(IW) = OBX%AREA_0(ABS(IOR))/OBX%AREA(ABS(IOR))
         IF (M%AREA_ADJUST(IW).EQ.0.) M%AREA_ADJUST(IW) = 1.
         ENDIF
C
C Prescribe exit velocity for surface cell
C
      IF (VELS(IBCX).NE.-999.) THEN
         M%UW0(IW) = VELS(IBCX)
         ELSE
         M%UW0(IW) = 0.
         ENDIF
C
      IF (I_OBST.GT.0 .AND. VFLUX(IBCX).NE.-999.) THEN
         OBX=>M%OBSTRUCTION(I_OBST)
         M%UW0(IW) = VFLUX(IBCX)/OBX%AREA(ABS(IOR))
         ENDIF
C
      T_ACTIVATE = -1.
C
C Check if there is a vent embedded in the surface
C
      VLOOP: DO N=1,M%NV
C
      VT => M%VENTS(N)
C
      IF (I_OBST.GT.0 .AND. VT%INDEX.EQ.2) CYCLE VLOOP
      IF (VT%IOR.NE.IOR) CYCLE VLOOP
C
      IF (IBCX.EQ.INTERPOLATED_INDEX) CYCLE VLOOP
      IF ((VT%INDEX.EQ.2 .OR. VT%INDEX.EQ.3) .AND.
     .     M%IJKW(9,IW).GT.0) CYCLE VLOOP
C
      IF (ABS(IOR).EQ.1) THEN
      IF (IOR.EQ. 1 .AND. I.NE.VT%I1  ) CYCLE VLOOP
      IF (IOR.EQ.-1 .AND. I.NE.VT%I1+1) CYCLE VLOOP
      IF (J.LT.VT%J1+1 .OR. J.GT.VT%J2) CYCLE VLOOP
      IF (K.LT.VT%K1+1 .OR. K.GT.VT%K2) CYCLE VLOOP
      ENDIF
C
      IF (ABS(IOR).EQ.2) THEN
      IF (IOR.EQ. 2 .AND. J.NE.VT%J1  ) CYCLE VLOOP
      IF (IOR.EQ.-2 .AND. J.NE.VT%J1+1) CYCLE VLOOP
      IF (I.LT.VT%I1+1 .OR. I.GT.VT%I2) CYCLE VLOOP
      IF (K.LT.VT%K1+1 .OR. K.GT.VT%K2) CYCLE VLOOP
      ENDIF
C
      IF (ABS(IOR).EQ.3) THEN
      IF (IOR.EQ. 3 .AND. K.NE.VT%K1  ) CYCLE VLOOP
      IF (IOR.EQ.-3 .AND. K.NE.VT%K1+1) CYCLE VLOOP
      IF (I.LT.VT%I1+1 .OR. I.GT.VT%I2) CYCLE VLOOP
      IF (J.LT.VT%J1+1 .OR. J.GT.VT%J2) CYCLE VLOOP
      ENDIF
C
      IBCX = VT%IBC
      M%AREA_ADJUST(IW) = VT%AREA_0/VT%AREA
      IF (M%AREA_ADJUST(IW).EQ.0.) M%AREA_ADJUST(IW) = 1.
C
C Set the velocity at each surface cell
C
      IF (VELS(IBCX).NE.-999.) THEN
         M%UW0(IW) = VELS(IBCX)
         ELSE
         M%UW0(IW) = 0.
         ENDIF
C
      IF (VFLUX(IBCX).NE.-999.) M%UW0(IW) = VFLUX(IBCX)/VT%AREA
C
C Special velocity profiles
C
      IF (IPROF(IBCX).EQ.1) THEN              ! Parabolic profile
         SELECT CASE(ABS(IOR))
         CASE(1)
            Y1 = M%Y(VT%J1)
            Y2 = M%Y(VT%J2)
            Z1 = M%Z(VT%K1)
            Z2 = M%Z(VT%K2)
            PY = 4.*(M%YC(J)-Y1)*(Y2-M%YC(J))/(Y2-Y1)**2
            PZ = 4.*(M%ZC(K)-Z1)*(Z2-M%ZC(K))/(Z2-Z1)**2
            M%UW0(IW) = M%UW0(IW)*PY*PZ
         CASE(2)
            X1 = M%X(VT%I1)
            X2 = M%X(VT%I2)
            Z1 = M%Z(VT%K1)
            Z2 = M%Z(VT%K2)
            PX = 4.*(M%XC(I)-X1)*(X2-M%XC(I))/(X2-X1)**2
            PZ = 4.*(M%ZC(K)-Z1)*(Z2-M%ZC(K))/(Z2-Z1)**2
            M%UW0(IW) = M%UW0(IW)*PX*PZ
         CASE(3)
            X1 = M%X(VT%I1)
            X2 = M%X(VT%I2)
            IF (CYLINDRICAL .AND. X1.EQ.0.) X1 = -X2
            Y1 = M%Y(VT%J1)
            Y2 = M%Y(VT%J2)
            PX = 4.*(M%XC(I)-X1)*(X2-M%XC(I))/(X2-X1)**2
            PY = 4.*(M%YC(J)-Y1)*(Y2-M%YC(J))/(Y2-Y1)**2
            IF (CYLINDRICAL) THEN
               M%UW0(IW) = M%UW0(IW)*PX
               ELSE
               M%UW0(IW) = M%UW0(IW)*PX*PY
               ENDIF
         END SELECT
      ENDIF
C
      IF (IPROF(IBCX).EQ.2) THEN     ! Atmospheric profile
         Z1 = M%Z(VT%K1)
         M%UW0(IW) = M%UW0(IW)*((M%ZC(K)-Z1)/Z0S(IBCX))**PLES(IBCX)
         ENDIF
C
      IF (IPROF(IBCX).EQ.3) THEN     ! 1D-Parabolic profile
         SELECT CASE(ABS(IOR))
         CASE(1)
            Y1 = M%Y(VT%J1)
            Y2 = M%Y(VT%J2)
            PY = 4.*(M%YC(J)-Y1)*(Y2-M%YC(J))/(Y2-Y1)**2
            M%UW0(IW) = M%UW0(IW)*PY
         CASE(2)
            X1 = M%X(VT%I1)
            X2 = M%X(VT%I2)
            PX = 4.*(M%XC(I)-X1)*(X2-M%XC(I))/(X2-X1)**2
            M%UW0(IW) = M%UW0(IW)*PX
         CASE(3)
            X1 = M%X(VT%I1)
            X2 = M%X(VT%I2)
            Y1 = M%Y(VT%J1)
            Y2 = M%Y(VT%J2)
            PX = 4.*(M%XC(I)-X1)*(X2-M%XC(I))/(X2-X1)**2
            PY = 4.*(M%YC(J)-Y1)*(Y2-M%YC(J))/(Y2-Y1)**2
            M%UW0(IW) = M%UW0(IW)*PX*PY
         END SELECT
      ENDIF
C
C Check if fire spreads radially
C
      IF (VT%X0.GT.-900.) THEN
         DIST = SQRT((M%XW(IW)-VT%X0)**2 +
     .               (M%YW(IW)-VT%Y0)**2 +
     .               (M%ZW(IW)-VT%Z0)**2)
         T_ACTIVATE = DIST/VT%FIRE_SPREAD_RATE
         ENDIF
C
C Miscellaneous settings
C
      IF (VT%INDEX.EQ.2) THEN  ! Open boundary stuff
         M%IV(IW) = 2
         M%QPYR(IW) = VT%TMP_OUTSIDE  ! QPYR is just a surrogate
         ENDIF
C
      IF (VT%INDEX.EQ.3) M%IV(IW) = 3  ! Mirror boundary
C
      IF (VT%CLOSED) THEN 
         M%IV(IW)   = 1
         T_ACTIVATE = 1000000.
         ENDIF
C
      M%IJKW(5,IW) = IBCX
C
      ENDDO VLOOP
C
C Compute mass of solid cells in cases where object burns away
C
      OBSTRUCTION_IF: IF (I_OBST.GT.0 .AND. BURNAWAY(IBCX)) THEN
      OBX=>M%OBSTRUCTION(I_OBST)
      IF ( DENSITY_F(IBCX).GT.0. .AND. 
     .     M%SOLID(IC) .AND.
     .     MINVAL(OBX%DIMENSIONS(:)).GT.0 .AND.
     .     M%CELL_MASS(IC).GT.1000000.) 
     .M%CELL_MASS(IC) = DENSITY_F(IBCX)*M%DX(I)*M%RC(I)*M%DY(J)*M%DZ(K)
      ENDIF OBSTRUCTION_IF
C
C Set ignition time of each boundary cell
C
      IF (T_ACTIVATE.LT.0.) THEN
         M%TW(IW) = TIGNS(IBCX)
      ELSE
         M%TW(IW) = T_ACTIVATE
      ENDIF
C
C Initialize solid temperature
C
      WC => M%WALL(IW)
      NWP_MAX = NWP(IBCX)
      IF (IBACK(IBCX).LE.0) NWP_MAX = MAX(NWP_MAX,NWP(-IBACK(IBCX)))
      ALLOCATE(WC%TMP_S(0:NWP_MAX+1),STAT=IZERO)
      CALL ChkMemErr('INIT','TMP_S',IZERO)
C
      WC%TMP_S(:)   = TMP_P0(IBCX)
      M%TMP_F(IW)   = TMP_P0(IBCX)
      M%TMP_B(IW)   = TMP_P0(IBCX)
      M%TMP_W(IW)   = TMP_P0(IBCX)
C
C Set ramp-up indicators
C
      IF (TAUQ(IBCX).LT.1000000.) M%ITAUQ(IW) = IBCX
      IF (TAUQ(IBCX).GT.1000000.) M%ITAUQ(IW) = -(TAUQ(IBCX)-1000000.)
      IF (TAUQ(IBCX).EQ.      0.) M%ITAUQ(IW) = 999
      IF (TAUV(IBCX).LT.1000000.) M%ITAUV(IW) = IBCX
      IF (TAUV(IBCX).GT.1000000.) M%ITAUV(IW) = -(TAUV(IBCX)-1000000.)
      IF (TAUV(IBCX).EQ.      0.) M%ITAUV(IW) = 999
      DO NN=0,NSPEC
      IF (TAUMF(IBCX,NN).LT.1000000.) M%ITAUMF(IW,NN) = IBCX
      IF (TAUMF(IBCX,NN).GT.1000000.) M%ITAUMF(IW,NN) =
     .                                -(TAUMF(IBCX,NN)-1000000.)
      IF (TAUMF(IBCX,NN).EQ.      0.) M%ITAUMF(IW,NN) = 999
      ENDDO
C
C Misc
C
      M%E_WALL(IW) = E_WALL_S(IBCX)
      IF (SURF_PARTICLE_CLASS(IBCX).GE.0) 
     .   M%CELL_PARTICLE_CLASS(IW) = SURF_PARTICLE_CLASS(IBCX)
      M%NPPCW(IW) = NPPCS(IBCX)      ! Number of particles per cell
C
      END SUBROUTINE INIT_WALL_CELL
C
C
      SUBROUTINE OPEN_AND_CLOSE(T,NM)
C
C Check to see if a cell is to be removed
C
      REAL(EB) T
      INTEGER N,II,JJ,KK,IIG,JJG,KKG,NWC0,IW,IOR,IC,IBC,I,J,K
      INTEGER, INTENT(IN) :: NM
      LOGICAL :: REMOVE_IT,SOLID_ONLY
      TYPE (OBSTRUCTION_TYPE), POINTER :: OB
      TYPE (VENTS_TYPE), POINTER :: VT
C
      CALL UNPACK_VAR(NM)
C
      NWC0 = NWC
C
      WALL_CELL_LOOP: DO IW=1,NWC0
      IF (IV(IW).NE.1) CYCLE WALL_CELL_LOOP
      IBC = IJKW(5,IW)
      IF (.NOT.BURNAWAY(IBC)) CYCLE WALL_CELL_LOOP
C
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IC  = ICA(II,JJ,KK)
C
      IF (CELL_MASS(IC).LE.0.) THEN
         IIG = IJKW(6,IW)
         JJG = IJKW(7,IW)
         KKG = IJKW(8,IW)
         TMP(II,JJ,KK) = TMP(IIG,JJG,KKG)
         RHO(II,JJ,KK) = RHO(IIG,JJG,KKG)
         IF (NSPEC.GT.0) THEN
         RSUM(II,JJ,KK) = RSUM(IIG,JJG,KKG)
         YY(II,JJ,KK,1:NSPEC) = YY(IIG,JJG,KKG,1:NSPEC)
         ENDIF
         IF (IW.GT.NEWC) CALL REMOVE_OBST(NM,II,II,JJ,JJ,KK,KK)
         CELL_PARTICLE_CLASS(IW) = -1
         ENDIF
C
      IF (MASS_LOSS(IW).GT.SURFACE_DENSITY_S(IBC) .AND.
     .    IBACK(IBC).EQ.3 .AND. IW.GT.NEWC) THEN
         IOR = IJKW(4,IW)
         IF (SOLID(IC)) THEN
            CALL REMOVE_OBST(NM,II,II,JJ,JJ,KK,KK)
         ELSE
            SELECT CASE(IOR)
            CASE( 1) ; CALL REMOVE_OBST(NM,II+1,II,JJ,JJ,KK,KK)
            CASE(-1) ; CALL REMOVE_OBST(NM,II,II-1,JJ,JJ,KK,KK)
            CASE( 2) ; CALL REMOVE_OBST(NM,II,II,JJ+1,JJ,KK,KK)
            CASE(-2) ; CALL REMOVE_OBST(NM,II,II,JJ,JJ-1,KK,KK)
            CASE( 3) ; CALL REMOVE_OBST(NM,II,II,JJ,JJ,KK+1,KK)
            CASE(-3) ; CALL REMOVE_OBST(NM,II,II,JJ,JJ,KK,KK-1)
            END SELECT
         ENDIF
         CELL_PARTICLE_CLASS(IW) = -1
         ENDIF
C
      IF (CELL_MASS(IC).LE.0. .AND. IW.GT.NEWC) THEN
         IOR = IJKW(4,IW)
         SOLID_ONLY = .FALSE.
         CALL GET_OBST(NM,II,JJ,KK,IOR,SOLID_ONLY,N)
         OB=>OBSTRUCTION(N)
         IF (N.GT.0 .AND. MINVAL(OB%DIMENSIONS(:)).GT.0) THEN
            REMOVE_IT = .TRUE.
            DO K=OB%K1+1,OB%K2
            DO J=OB%J1+1,OB%J2
            DO I=OB%I1+1,OB%I2
            IF (SOLID(ICA(I,J,K))) REMOVE_IT = .FALSE.
            ENDDO
            ENDDO
            ENDDO
            IF (REMOVE_IT) OB%T_REMOVE = T-1.
            ENDIF
         ENDIF
C
      ENDDO WALL_CELL_LOOP
C
C Check to see if an entire obstacle is to be removed
C
      BLOOP: DO N=1,NB
      OB=>OBSTRUCTION(N)
      IF (T.LT.OB%T_REMOVE) CYCLE BLOOP
      OB%HIDDEN = .TRUE.
      CALL REMOVE_OBST(NM,OB%I1+1,OB%I2,OB%J1+1,OB%J2,OB%K1+1,OB%K2)
      OB%T_REMOVE = 1000000.
      IF (N_STRINGS+2.GT.N_STRINGS_MAX) THEN
         CALL RE_ALLOCATE_STRINGS(NM)
         STRING => MESH(NM)%STRING
         ENDIF
      N_STRINGS = N_STRINGS + 1
      WRITE(STRING(N_STRINGS),'(A,I3)') 'HIDE_OBST',NM
      N_STRINGS = N_STRINGS + 1
      WRITE(STRING(N_STRINGS),'(I6,F10.2)') N,T
      ENDDO BLOOP
C
C Check to see if an entire obstacle is to be created
C
      BLOOP2: DO N=1,NB
      OB=>OBSTRUCTION(N)
      IF (T.LT.OB%T_CREATE) CYCLE BLOOP2
      OB%HIDDEN = .FALSE.
      CALL CREATE_OBST(NM,OB%I1+1,OB%I2,OB%J1+1,OB%J2,OB%K1+1,OB%K2,N)
      OB%T_CREATE = 1000000.
      IF (N_STRINGS+2.GT.N_STRINGS_MAX) THEN
         CALL RE_ALLOCATE_STRINGS(NM)
         STRING => MESH(NM)%STRING
         ENDIF
      N_STRINGS = N_STRINGS + 1
      WRITE(STRING(N_STRINGS),'(A,I3)') 'SHOW_OBST',NM
      N_STRINGS = N_STRINGS + 1
      WRITE(STRING(N_STRINGS),'(I6,F10.2)') N,T
      ENDDO BLOOP2
C
C Check to see if a vent should be opened
C
      VLOOP: DO N=1,NV
C
      VT => VENTS(N)
C
      IF ( T.LT.VT%T_OPEN .AND. 
     .     (P0.LT.P0_MAX .OR. VT%INDEX.NE.2) ) CYCLE VLOOP
      IF (VT%T_CLOSE.EQ.999999.) CYCLE VLOOP
      SEARCH: DO IW=1,NWC
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IOR = IJKW(4,IW)
      IF (IOR.NE.VT%IOR) CYCLE SEARCH
      SELECT CASE(ABS(IOR))
      CASE(1)
      IF (IOR.EQ. 1 .AND.    II.NE.VT%I1)   CYCLE SEARCH
      IF (IOR.EQ.-1 .AND.    II.NE.VT%I1+1) CYCLE SEARCH
      IF (JJ.LT.VT%J1+1 .OR. JJ.GT.VT%J2)   CYCLE SEARCH
      IF (KK.LT.VT%K1+1 .OR. KK.GT.VT%K2)   CYCLE SEARCH
      CASE(2)
      IF (IOR.EQ. 2 .AND.    JJ.NE.VT%J1)   CYCLE SEARCH
      IF (IOR.EQ.-2 .AND.    JJ.NE.VT%J1+1) CYCLE SEARCH
      IF (II.LT.VT%I1+1 .OR. II.GT.VT%I2)   CYCLE SEARCH
      IF (KK.LT.VT%K1+1 .OR. KK.GT.VT%K2)   CYCLE SEARCH
      CASE(3)
      IF (IOR.EQ. 3 .AND.    KK.NE.VT%K1)   CYCLE SEARCH
      IF (IOR.EQ.-3 .AND.    KK.NE.VT%K1+1) CYCLE SEARCH
      IF (II.LT.VT%I1+1 .OR. II.GT.VT%I2)   CYCLE SEARCH
      IF (JJ.LT.VT%J1+1 .OR. JJ.GT.VT%J2)   CYCLE SEARCH
      END SELECT
      IF (VT%INDEX.NE.2) THEN
         TW(IW) = T
         CYCLE SEARCH
         ENDIF
      CALL BLKCLL(NM,II,II,JJ,JJ,KK,KK,0)
      IV(IW) = 2
      IJKW(5,IW) = OPEN_INDEX
      ENDDO SEARCH
      VT%T_OPEN = 1000000.
      IF (N_STRINGS+2.GT.N_STRINGS_MAX) THEN
         CALL RE_ALLOCATE_STRINGS(NM)
         STRING => MESH(NM)%STRING
         ENDIF
      N_STRINGS = N_STRINGS + 1
      WRITE(STRING(N_STRINGS),'(A,I3)') 'OPEN_VENT',NM
      N_STRINGS = N_STRINGS + 1
      WRITE(STRING(N_STRINGS),'(I6,F10.2)') N,T
      ENDDO VLOOP
C
C Check to see if a vent should be closed
C
      VLOOP2: DO N=1,NV
C
      VT => VENTS(N)
C
      IF (T.LT.VT%T_CLOSE) CYCLE VLOOP2
      SEARCH2: DO IW=1,NWC
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IOR = IJKW(4,IW)
      IF (IOR.NE.VT%IOR) CYCLE SEARCH2
      SELECT CASE(ABS(IOR))
      CASE(1)
      IF (IOR.EQ. 1 .AND.    II.NE.VT%I1)   CYCLE SEARCH2
      IF (IOR.EQ.-1 .AND.    II.NE.VT%I1+1) CYCLE SEARCH2
      IF (JJ.LT.VT%J1+1 .OR. JJ.GT.VT%J2)   CYCLE SEARCH2
      IF (KK.LT.VT%K1+1 .OR. KK.GT.VT%K2)   CYCLE SEARCH2
      CASE(2)
      IF (IOR.EQ. 2 .AND.    JJ.NE.VT%J1)   CYCLE SEARCH2
      IF (IOR.EQ.-2 .AND.    JJ.NE.VT%J1+1) CYCLE SEARCH2
      IF (II.LT.VT%I1+1 .OR. II.GT.VT%I2)   CYCLE SEARCH2
      IF (KK.LT.VT%K1+1 .OR. KK.GT.VT%K2)   CYCLE SEARCH2
      CASE(3)
      IF (IOR.EQ. 3 .AND.    KK.NE.VT%K1)   CYCLE SEARCH2
      IF (IOR.EQ.-3 .AND.    KK.NE.VT%K1+1) CYCLE SEARCH2
      IF (II.LT.VT%I1+1 .OR. II.GT.VT%I2)   CYCLE SEARCH2
      IF (JJ.LT.VT%J1+1 .OR. JJ.GT.VT%J2)   CYCLE SEARCH2
      END SELECT
      IF (VT%INDEX.NE.2) THEN
         TW(IW) = 1000000.
         CYCLE SEARCH2
         ENDIF
      CALL BLKCLL(NM,II,II,JJ,JJ,KK,KK,1)
      IV(IW) = 1
      IF (RADIATION) THEN
         NULLIFY(MESH(NM)%WALL(IW)%ILW)
         ALLOCATE(MESH(NM)%WALL(IW)%ILW(NRA,NSB))
         MESH(NM)%WALL(IW)%ILW = SIGMA*TMPA4/PI
         ENDIF
      ENDDO SEARCH2
      VT%T_CLOSE =  999999.
      VT%INDEX   = 1
      IF (N_STRINGS+2.GT.N_STRINGS_MAX) THEN
         CALL RE_ALLOCATE_STRINGS(NM)
         STRING => MESH(NM)%STRING
         ENDIF
      N_STRINGS = N_STRINGS + 1
      WRITE(STRING(N_STRINGS),'(A,I3)') 'CLOSE_VENT',NM
      N_STRINGS = N_STRINGS + 1
      WRITE(STRING(N_STRINGS),'(I6,F10.2)') N,T
      ENDDO VLOOP2
C
C
      END SUBROUTINE OPEN_AND_CLOSE
C
C
      SUBROUTINE REMOVE_OBST(NM,I1,I2,J1,J2,K1,K2)
C
      INTEGER I1,I2,J1,J2,K1,K2,I,J,K,NN,IE,IW,IBC,IOR,IC
      INTEGER, INTENT(IN) :: NM
      LOGICAL :: SOLID_ONLY
      TYPE (OBSTRUCTION_TYPE), POINTER :: OBO
C
      CALL UNPACK_VAR(NM)
C
C Change the volume of gas in the entire domain
C
      DO K=K1,K2
      DO J=J1,J2
      DO I=I1,I2
      IC = ICA(I,J,K)
      IF (SOLID(IC)) VOL = VOL + RC(I)*DX(I)*DY(J)*DZ(K)
      ENDDO
      ENDDO
      ENDDO
C
C Remove the blockage
C
      CALL BLKCLL(NM,I1,I2,J1,J2,K1,K2,0)
C
C Nullify wall cells on blockage that is to be removed
C
      DO K=K1,K2
      DO J=J1,J2
      DO I=I1,I2
      DO NN=-3,3
      IW = IWA(ICA(I,J,K),NN)
C      
      IF ( IW.LE.0 ) THEN 
          CONTINUE
      ELSE IF ( IJKW(9,IW) .EQ. 0 ) THEN 
          IV(IW) = 0
      ENDIF
      
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      DO K=K1,K2
      DO J=J1,J2
      IW = IWA(ICA(I1-1,J,K),1)

      IF ( IW.LE.0 ) THEN 
          CONTINUE
      ELSE IF ( IJKW(9,IW) .EQ. 0 ) THEN 
          IV(IW) = 0
      ENDIF
     
      IW = IWA(ICA(I2+1,J,K),-1)
C      
      IF ( IW.LE.0 ) THEN 
          CONTINUE
      ELSE IF ( IJKW(9,IW) .EQ. 0 ) THEN 
          IV(IW) = 0
      ENDIF
     
      ENDDO
      ENDDO
C
      DO K=K1,K2
      DO I=I1,I2
      IW = IWA(ICA(I,J1-1,K),2)
C      
      IF ( IW.LE.0 ) THEN 
          CONTINUE
      ELSE IF ( IJKW(9,IW) .EQ. 0 ) THEN 
          IV(IW) = 0
      ENDIF
     
      IW = IWA(ICA(I,J2+1,K),-2)
C      
      IF ( IW.LE.0 ) THEN 
          CONTINUE
      ELSE IF ( IJKW(9,IW) .EQ. 0 ) THEN 
          IV(IW) = 0
      ENDIF
     
      ENDDO
      ENDDO
C
      DO J=J1,J2
      DO I=I1,I2
      IW = IWA(ICA(I,J,K1-1),3)
C      
      IF ( IW.LE.0 ) THEN 
          CONTINUE
      ELSE IF ( IJKW(9,IW) .EQ. 0 ) THEN 
          IV(IW) = 0
      ENDIF
     
      IW = IWA(ICA(I,J,K2+1),-3)
C      
      IF ( IW.LE.0 ) THEN 
          CONTINUE
      ELSE IF ( IJKW(9,IW) .EQ. 0 ) THEN 
          IV(IW) = 0
      ENDIF
     
      ENDDO
      ENDDO
C
C Nullify block edges on blockage that is to be removed
C
      DO K=K1-1,K2
      DO J=J1-1,J2
      DO I=I1  ,I2
      IE = IEA(ICA(I,J,K),4)
      IF (IE.GT.0 .AND. J.NE.0 .AND. J.NE.JBAR .AND. 
     .                  K.NE.0 .AND. K.NE.KBAR) IJKE(4,IE) = 0
C      
      IF ( IE .LE. 0 ) THEN 
          CONTINUE
      ELSE IF ( IJKE(7,IE) .GT. 0 ) THEN 
          IJKE(5,IE)=INTERPOLATED_INDEX
      ENDIF
     
      ENDDO
      ENDDO
      ENDDO
C
      DO K=K1-1,K2
      DO J=J1  ,J2
      DO I=I1-1,I2
      IE = IEA(ICA(I,J,K),8)
      IF (IE.GT.0 .AND. I.NE.0 .AND. I.NE.IBAR .AND.
     .                  K.NE.0 .AND. K.NE.KBAR) IJKE(4,IE) = 0

C     APG changed IF stmt that was illegaly using side-effects
C      IF (IE.GT.0 .AND. IJKE(7,IE).GT.0) IJKE(5,IE)=INTERPOLATED_INDEX
C      
      IF ( IE .LE. 0 ) THEN 
          CONTINUE
      ELSE IF ( IJKE(7,IE) .GT. 0 ) THEN 
          IJKE(5,IE)=INTERPOLATED_INDEX
      ENDIF
     
      ENDDO
      ENDDO
      ENDDO
C
      DO K=K1  ,K2
      DO J=J1-1,J2
      DO I=I1-1,I2
      IE = IEA(ICA(I,J,K),12)
      IF (IE.GT.0 .AND. I.NE.0 .AND. I.NE.IBAR .AND.
     .                  J.NE.0 .AND. J.NE.JBAR) IJKE(4,IE) = 0
  
C     APG changed IF stmt that was illegaly using side-effects
C      IF (IE.GT.0 .AND. IJKE(7,IE).GT.0) IJKE(5,IE)=INTERPOLATED_INDEX
C      
      IF ( IE .LE. 0 ) THEN 
          CONTINUE
      ELSE IF ( IJKE(7,IE) .GT. 0 ) THEN 
          IJKE(5,IE)=INTERPOLATED_INDEX
      ENDIF
     
      ENDDO
      ENDDO
      ENDDO
C
C Initialize exposed surfaces
C
      LOOP1: DO K=K1,K2
      LOOP2: DO J=J1,J2
      I   = I1
      IC  = ICA(I-1,J,K)
      IF (I.EQ.1) THEN
         IW = IWA(IC,1) ; IF (IW.EQ.0) CYCLE LOOP2
         IF (IV(IW).EQ.0) IV(IW) = 1
         IF (IJKW(9,IW).GT.0) THEN
            IV(IW) = 4
            IJKW(5,IW) = INTERPOLATED_INDEX
            ENDIF
         CYCLE LOOP2
         ENDIF
      IF (.NOT.SOLID(IC)) CYCLE LOOP2
      IOR = 1
      NWC = NWC + 1
      IW  = NWC
      SOLID_ONLY = .TRUE.
      CALL GET_OBST(NM,I-1,J,K,IOR,SOLID_ONLY,NN)
      OBO=>OBSTRUCTION(NN)
      IBC = OBO%IBC(IOR)
      CALL INIT_WALL_CELL(NM,I-1,J,K,NN,IW,IOR,IBC)
      ENDDO LOOP2
      ENDDO LOOP1
C
      LOOP3: DO K=K1,K2
      LOOP4: DO J=J1,J2
      I   = I2
      IC  = ICA(I+1,J,K)
      IF (I.EQ.IBAR) THEN
         IW = IWA(IC,-1) ; IF (IW.EQ.0) CYCLE LOOP4
         IF (IV(IW).EQ.0) IV(IW) = 1
         IF (IJKW(9,IW).GT.0) THEN
            IV(IW) = 4
            IJKW(5,IW) = INTERPOLATED_INDEX
            ENDIF
         CYCLE LOOP4
         ENDIF
      IF (.NOT.SOLID(IC)) CYCLE LOOP4
      IOR = -1
      NWC = NWC + 1
      IW  = NWC
      SOLID_ONLY = .TRUE.
      CALL GET_OBST(NM,I,J,K,IOR,SOLID_ONLY,NN)
      OBO=>OBSTRUCTION(NN)
      IBC = OBO%IBC(IOR)
      CALL INIT_WALL_CELL(NM,I+1,J,K,NN,IW,IOR,IBC)
      ENDDO LOOP4
      ENDDO LOOP3
C
      LOOP5: DO K=K1,K2
      LOOP6: DO I=I1,I2
      J   = J1
      IC  = ICA(I,J-1,K)
      IF (J.EQ.1) THEN
         IW = IWA(IC,2) ; IF (IW.EQ.0) CYCLE LOOP6
         IF (IV(IW).EQ.0) IV(IW) = 1
         IF (IJKW(9,IW).GT.0) THEN
            IV(IW) = 4
            IJKW(5,IW) = INTERPOLATED_INDEX
            ENDIF
         CYCLE LOOP6
         ENDIF
      IF (.NOT.SOLID(IC)) CYCLE LOOP6
      IOR = 2
      NWC = NWC + 1
      IW  = NWC
      SOLID_ONLY = .TRUE.
      CALL GET_OBST(NM,I,J-1,K,IOR,SOLID_ONLY,NN)
      OBO=>OBSTRUCTION(NN)
      IBC = OBO%IBC(IOR)
      CALL INIT_WALL_CELL(NM,I,J-1,K,NN,IW,IOR,IBC)
      ENDDO LOOP6
      ENDDO LOOP5
C
      LOOP7: DO K=K1,K2
      LOOP8: DO I=I1,I2
      J   = J2
      IC  = ICA(I,J+1,K)
      IF (J.EQ.JBAR) THEN
         IW = IWA(IC,-2) ; IF (IW.EQ.0) CYCLE LOOP8
         IF (IV(IW).EQ.0) IV(IW) = 1
         IF (IJKW(9,IW).GT.0) THEN
            IV(IW) = 4
            IJKW(5,IW) = INTERPOLATED_INDEX
            ENDIF
         CYCLE LOOP8
         ENDIF
      IF (.NOT.SOLID(IC)) CYCLE LOOP8
      IOR = -2
      NWC = NWC + 1
      IW  = NWC
      SOLID_ONLY = .TRUE.
      CALL GET_OBST(NM,I,J,K,IOR,SOLID_ONLY,NN)
      OBO=>OBSTRUCTION(NN)
      IBC = OBO%IBC(IOR)
      CALL INIT_WALL_CELL(NM,I,J+1,K,NN,IW,IOR,IBC)
      ENDDO LOOP8
      ENDDO LOOP7
C
      LOOP9:  DO J=J1,J2
      LOOP10: DO I=I1,I2
      K   = K1
      IC  = ICA(I,J,K-1)
      IF (K.EQ.1) THEN
         IW = IWA(IC,3) ; IF (IW.EQ.0) CYCLE LOOP10
         IF (IV(IW).EQ.0) IV(IW) = 1
         IF (IJKW(9,IW).GT.0) THEN
            IV(IW) = 4
            IJKW(5,IW) = INTERPOLATED_INDEX
            ENDIF
         CYCLE LOOP10
         ENDIF
      IF (.NOT.SOLID(IC)) CYCLE LOOP10
      IOR = 3
      NWC = NWC + 1
      IW  = NWC
      SOLID_ONLY = .TRUE.
      CALL GET_OBST(NM,I,J,K-1,IOR,SOLID_ONLY,NN)
      OBO=>OBSTRUCTION(NN)
      IBC = OBO%IBC(IOR)
      CALL INIT_WALL_CELL(NM,I,J,K-1,NN,IW,IOR,IBC)
      ENDDO LOOP10
      ENDDO LOOP9
C
      LOOP11: DO J=J1,J2
      LOOP12: DO I=I1,I2
      K   = K2
      IC  = ICA(I,J,K+1)
      IF (K.EQ.KBAR) THEN
         IW = IWA(IC,-3) ; IF (IW.EQ.0) CYCLE LOOP12
         IF (IV(IW).EQ.0) IV(IW) = 1
         IF (IJKW(9,IW).GT.0) THEN
            IV(IW) = 4
            IJKW(5,IW) = INTERPOLATED_INDEX
            ENDIF
         CYCLE LOOP12
         ENDIF
      IF (.NOT.SOLID(IC)) CYCLE LOOP12
      IOR = -3
      NWC = NWC + 1
      IW  = NWC
      SOLID_ONLY = .TRUE.
      CALL GET_OBST(NM,I,J,K,IOR,SOLID_ONLY,NN)
      OBO=>OBSTRUCTION(NN)
      IBC = OBO%IBC(IOR)
      CALL INIT_WALL_CELL(NM,I,J,K+1,NN,IW,IOR,IBC)
      ENDDO LOOP12
      ENDDO LOOP11
C
      END SUBROUTINE REMOVE_OBST
C
C
      SUBROUTINE CREATE_OBST(NM,I1,I2,J1,J2,K1,K2,NN)
C
      INTEGER I1,I2,J1,J2,K1,K2,I,J,K,NN,IW,IBC,IOR,IC
      INTEGER, INTENT(IN) :: NM
      TYPE (OBSTRUCTION_TYPE), POINTER :: OB
C
      CALL UNPACK_VAR(NM)
      OB=>OBSTRUCTION(NN)
C
C Change the volume of gas in the entire domain
C
      DO K=K1,K2
      DO J=J1,J2
      DO I=I1,I2
      IC = ICA(I,J,K)
      IF (.NOT.SOLID(IC)) VOL = VOL - RC(I)*DX(I)*DY(J)*DZ(K)
      ENDDO
      ENDDO
      ENDDO
C
C Create the blockage 
C
      CALL BLKCLL(NM,I1,I2,J1,J2,K1,K2,1)
C
C Initialize exposed surfaces
C
      LOOP1: DO K=K1,K2
      LOOP2: DO J=J1,J2
      I   = I1
      IC  = ICA(I-1,J,K)
      IF (SOLID(IC)) THEN
         IW = IWA(IC,1)
         IF (IV(IW).EQ.1) IV(IW) = 0
         CYCLE LOOP2
         ENDIF
      IF (I.EQ.1) CYCLE LOOP2
      IF (IWA(IC,1).GT.0) THEN
         IW = IWA(IC,1)
      ELSE
         NWC = NWC + 1
         IW  = NWC
      ENDIF
      IOR = -1
      IBC = OB%IBC(IOR)
      CALL INIT_WALL_CELL(NM,I,J,K,NN,IW,IOR,IBC)
      ENDDO LOOP2
      ENDDO LOOP1
C
      LOOP3: DO K=K1,K2
      LOOP4: DO J=J1,J2
      I   = I2
      IC  = ICA(I+1,J,K)
      IF (SOLID(IC)) THEN
         IW = IWA(IC,-1)
         IF (IV(IW).EQ.1) IV(IW) = 0
         CYCLE LOOP4
         ENDIF
      IF (I.EQ.IBAR) CYCLE LOOP4
      IF (IWA(IC,-1).GT.0) THEN
         IW = IWA(IC,-1)
      ELSE
         NWC = NWC + 1
         IW  = NWC
      ENDIF
      IOR = 1
      IBC = OB%IBC(IOR)
      CALL INIT_WALL_CELL(NM,I,J,K,NN,IW,IOR,IBC)
      ENDDO LOOP4
      ENDDO LOOP3
C
      LOOP5: DO K=K1,K2
      LOOP6: DO I=I1,I2
      J   = J1
      IC  = ICA(I,J-1,K)
      IF (SOLID(IC)) THEN
         IW = IWA(IC,2)
         IF (IV(IW).EQ.1) IV(IW) = 0
         CYCLE LOOP6
         ENDIF
      IF (J.EQ.1) CYCLE LOOP6
      IF (IWA(IC,2).GT.0) THEN
         IW = IWA(IC,2)
      ELSE
         NWC = NWC + 1
         IW  = NWC
      ENDIF
      IOR = -2
      IBC = OB%IBC(IOR)
      CALL INIT_WALL_CELL(NM,I,J,K,NN,IW,IOR,IBC)
      ENDDO LOOP6
      ENDDO LOOP5
C
      LOOP7: DO K=K1,K2
      LOOP8: DO I=I1,I2
      J   = J2
      IC  = ICA(I,J+1,K)
      IF (SOLID(IC)) THEN
         IW = IWA(IC,-2)
         IF (IV(IW).EQ.1) IV(IW) = 0
         CYCLE LOOP8
         ENDIF
      IF (J.EQ.JBAR) CYCLE LOOP8
      IF (IWA(IC,-2).GT.0) THEN
         IW = IWA(IC,-2)
      ELSE
         NWC = NWC + 1
         IW  = NWC
      ENDIF
      IOR = 2
      IBC = OB%IBC(IOR)
      CALL INIT_WALL_CELL(NM,I,J,K,NN,IW,IOR,IBC)
      ENDDO LOOP8
      ENDDO LOOP7
C
      LOOP9:  DO J=J1,J2
      LOOP10: DO I=I1,I2
      K   = K1
      IC  = ICA(I,J,K-1)
      IF (SOLID(IC)) THEN
         IW = IWA(IC,3)
         IF (IV(IW).EQ.1) IV(IW) = 0
         CYCLE LOOP10
         ENDIF
      IF (K.EQ.1) CYCLE LOOP10
      IF (IWA(IC,3).GT.0) THEN
         IW = IWA(IC,3)
      ELSE
         NWC = NWC + 1
         IW  = NWC
      ENDIF
      IOR = -3
      IBC = OB%IBC(IOR)
      CALL INIT_WALL_CELL(NM,I,J,K,NN,IW,IOR,IBC)
      ENDDO LOOP10
      ENDDO LOOP9
C
      LOOP11: DO J=J1,J2
      LOOP12: DO I=I1,I2
      K   = K2
      IC  = ICA(I,J,K+1)
      IF (SOLID(IC)) THEN
         IW = IWA(IC,-3)
         IF (IV(IW).EQ.1) IV(IW) = 0
         CYCLE LOOP12
         ENDIF
      IF (K.EQ.KBAR) CYCLE LOOP12
      IF (IWA(IC,-3).GT.0) THEN
         IW = IWA(IC,-3)
      ELSE
         NWC = NWC + 1
         IW  = NWC
      ENDIF
      IOR = 3
      IBC = OB%IBC(IOR)
      CALL INIT_WALL_CELL(NM,I,J,K,NN,IW,IOR,IBC)
      ENDDO LOOP12
      ENDDO LOOP11
C
      END SUBROUTINE CREATE_OBST
C
C
      END MODULE INIT
