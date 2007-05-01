      MODULE SPRK
C
      USE PREC
      USE CONS
      USE PACKER
      USE TRAN
C
      IMPLICIT NONE
C
      PRIVATE
      PUBLIC INSERT_DROPLETS,INSERT_PARTICLES,
     .       TRACK_DROPLETS,CHECK_SPRINKLERS,
     .       CHECK_HEAT_DETECTORS,SMOKE_DETECTORS,
     .       INITIALIZE_DROPLETS,INITIALIZE_TREES
      TYPE (DROPLET_TYPE), POINTER :: DR
      TYPE (LAGRANGIAN_TYPE), POINTER :: LP
C
C
      CONTAINS
C
C
      SUBROUTINE INITIALIZE_DROPLETS(NM)
C
C Insert droplets into the domain at the start of calculation
C
      REAL(EB) RN,MASS_SUM,LL,UL,BIN_SIZE
      REAL(EB) VOL1, VOL2, X1, X2, Y1, Y2, Z1, Z2
      INTEGER I,J,II,JJ,KK,IL,IU,STRATUM,IPC,NTOT
      INTEGER, INTENT(IN) :: NM
C
      TNOW=SECOND()
C
      CALL UNPACK_VAR(NM)
C
      PART_CLASS_LOOP: DO IPC=1,NPC
C
      LP=>LAGRANGIAN(IPC)
C
C     Check if droplets evaporate into a SPECIES
C
      SPECIES_LOOP: DO I=1,NSPEC
         IF (LP%SPECIES=='null') CYCLE SPECIES_LOOP
         IF (SPECIES_ID(I)==LP%SPECIES) THEN
            LP%SPECIES_INDEX = I
            EXIT SPECIES_LOOP
         ENDIF
      ENDDO SPECIES_LOOP
C
C     If particles have size distribution, initialize here
C
      IF_SIZE_DISTRIBUTION: IF (LP%DIAMETER.GT.0.) THEN
C
      CALL DROPLET_SIZE_DISTRIBUTION(LP%DIAMETER,
     .     LP%R_CDF(:),LP%CDF(:),NDC,LP%GAMMA,LP%SIGMA)
C
      BIN_SIZE = LP%R_CDF(NDC)/REAL(NSTRATA,EB)
C
      STRATIFY: DO I=1,NSTRATA
         LL = (I-1)*BIN_SIZE
         UL =  I   *BIN_SIZE
         LL_LOOP: DO J=1,NDC
            IF (LP%R_CDF(J).GT.LL) THEN
               IL = J-1 ; LP%IL_CDF(I) = J-1
               EXIT LL_LOOP
               ENDIF
            ENDDO LL_LOOP
         UL_LOOP: DO J=NDC,1,-1
            IF (LP%R_CDF(J).LE.UL) THEN
               IU = J ; LP%IU_CDF(I) = J
               EXIT UL_LOOP
               ENDIF
            ENDDO UL_LOOP
C
         LP%W_CDF(I) = LP%CDF(IU) - LP%CDF(IL)
      ENDDO STRATIFY
C
      ENDIF IF_SIZE_DISTRIBUTION
C
C     If there is an initial number of droplets/particles, initialize
C
      IF (LP%TREE)           CYCLE PART_CLASS_LOOP
      IF (LP%N_INITIAL.EQ.0) CYCLE PART_CLASS_LOOP
C
      IF (LP%X1 .EQ. 0.0 .AND. LP%X2 .EQ. 0.0 .AND.
     .    LP%Y1 .EQ. 0.0 .AND. LP%Y2 .EQ. 0.0 .AND.
     .    LP%Z1 .EQ. 0.0 .AND. LP%Z2 .EQ. 0.0 ) THEN
         X1 = XS ; X2 = XF
         Y1 = YS ; Y2 = YF
         Z1 = ZS ; Z2 = ZF
         VOL2 = (XF - XS) * (YF - YS) * (ZF - ZS)
         VOL1 = VOL2
         ELSE
         IF (LP%X1.GT.XF .OR. LP%X2.LT.XS .OR.
     .       LP%Y1.GT.YF .OR. LP%Y2.LT.YS .OR.
     .       LP%Z1.GT.ZF .OR. LP%Z2.LT.ZS) THEN
            CYCLE PART_CLASS_LOOP
            ENDIF
         X1 = MAX(LP%X1,XS) ; X2 = MIN(LP%X2, XF)
         Y1 = MAX(LP%Y1,YS) ; Y2 = MIN(LP%Y2, YF)
         Z1 = MAX(LP%Z1,ZS) ; Z2 = MIN(LP%Z2, ZF)
         VOL2 = (LP%X2 - LP%X1) * (LP%Y2 - LP%Y1) * (LP%Z2 - LP%Z1)
         VOL1 = (X2 - X1) * (Y2 - Y1) * (Z2 - Z1)
         ENDIF
C     
      NTOT = INT(REAL(LP%N_INITIAL,EB)*VOL1/VOL2)
C
      MASS_SUM = 0.
      INITIALIZATION_LOOP: DO I=1,LP%N_INITIAL
      NLP = NLP + 1
      IF (NLP.GT.NLPDIM) THEN
         CALL RE_ALLOCATE_DROPLETS(1,NM,0)
         DROPLET=>MESH(NM)%DROPLET
         ENDIF
      DR=>DROPLET(NLP)
      BLK_LOOP:  DO
      CALL RANDOM_NUMBER(RN)
      DR%X = X1 + RN*(X2-X1)
      CALL RANDOM_NUMBER(RN)
      DR%Y = Y1 + RN*(Y2-Y1)
      CALL RANDOM_NUMBER(RN)
      DR%Z = Z1 + RN*(Z2-Z1)
      II = CELLSI(FLOOR((DR%X-XS)*RDXINT)) + 1.
      JJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT)) + 1.
      KK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT)) + 1.
      IF (.NOT.SOLID(ICA(II,JJ,KK))) EXIT BLK_LOOP
      ENDDO BLK_LOOP
      DR%U = 0.
      DR%V = 0.
      DR%W = 0.
      DR%TMP = LP%TMP_INITIAL
      DR%T   = 0.
      DR%R   = 0.
      DR%PWT = 1.
      DR%IOR = 0  
      DR%CLASS = IPC
      IF (MOD(NLP,LP%SAMPLING).EQ.0) THEN
         DR%SHOW = .TRUE.    
      ELSE
         DR%SHOW = .FALSE.    
      ENDIF
C
      IF (LP%DIAMETER.GT.0.) THEN
      STRATUM = MOD(NLP,NSTRATA) + 1
      IL = LP%IL_CDF(STRATUM)
      IU = LP%IU_CDF(STRATUM)
      CALL RANDOM_CHOICE(LP%CDF(IL:IU),LP%R_CDF(IL:IU),IU-IL,DR%R)
      DR%PWT = LP%W_CDF(STRATUM)
      MASS_SUM = MASS_SUM + DR%PWT*LP%FTPR*DR%R**3
      ENDIF
C
      ENDDO INITIALIZATION_LOOP
C
C     Adjust particle weighting factor PWT so that desired 
C     MASS_PER_VOLUME is achieved
C
      IF (LP%DIAMETER.GT.0.)
     .DROPLET(NLP-LP%N_INITIAL+1:NLP)%PWT = 
     .DROPLET(NLP-LP%N_INITIAL+1:NLP)%PWT*
     .        LP%MASS_PER_VOLUME*VOL1/MASS_SUM
C
      ENDDO PART_CLASS_LOOP
C
      TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
      END SUBROUTINE INITIALIZE_DROPLETS
C
C
      SUBROUTINE INITIALIZE_TREES(NM)
C
      REAL(EB) DELTAZ_BIN,CANOPY_LENGTH,CANOPY_VOLUME,COSINE,TANGENT,
     .         XC_MAX_ZC_MAX,XC_MAX_ZC_MIN,XC_MIN,XC_MAX,
     .         YC_MIN,YC_MAX,ZBIN_VOLUME,ZC_MIN,ZC_MAX,CANOPY_WIDTH
      INTEGER NCT,NNLP,NLP_TREE,NZB,NZBINS,P_PER_ZBIN,IPC
      REAL(EB) RN,MASS_SUM
      INTEGER I,II,JJ,KK,IL,IU,STRATUM
      INTEGER, INTENT(IN) :: NM
C
      CALL UNPACK_VAR(NM)
C
      TREE_LOOP: DO NCT=1,N_TREES
C
      IF (TREE_MESH(NCT).NE.NM) CYCLE TREE_LOOP
      IPC = TREE_PARTICLE_CLASS(NCT)
      LP=>LAGRANGIAN(IPC)
C
C     Build a conical tree
C
      CANOPY_WIDTH  = CANOPY_W(NCT)
      CANOPY_LENGTH = TREE_H(NCT) - CANOPY_B_H(NCT)
      TANGENT = 0.5*CANOPY_W(NCT)/CANOPY_LENGTH
      COSINE = CANOPY_LENGTH/
     .         SQRT(CANOPY_LENGTH**2 + (0.5*CANOPY_WIDTH)**2)
      NZBINS = 1./(1.-COSINE)
      DELTAZ_BIN = CANOPY_LENGTH/REAL(NZBINS,EB)
      CANOPY_VOLUME = PI*CANOPY_WIDTH**2*CANOPY_LENGTH/12.
C
      NLP_TREE = 0
      ZBIN_LOOP: DO NZB=1,NZBINS
       ZC_MIN = CANOPY_B_H(NCT) + (NZB-1)*DELTAZ_BIN
       ZC_MAX = ZC_MIN + DELTAZ_BIN
       XC_MAX_ZC_MIN = (CANOPY_LENGTH-(ZC_MIN-CANOPY_B_H(NCT)))*TANGENT
       XC_MAX_ZC_MAX = (CANOPY_LENGTH-(ZC_MAX-CANOPY_B_H(NCT)))*TANGENT
       ZBIN_VOLUME = PI*DELTAZ_BIN/3.*(
     .               XC_MAX_ZC_MIN**2 + XC_MAX_ZC_MAX**2 +
     .               XC_MAX_ZC_MIN*XC_MAX_ZC_MAX )
       P_PER_ZBIN  = LP%N_INITIAL*ZBIN_VOLUME/CANOPY_VOLUME
C
       NISP_LOOP1: DO NNLP=1,P_PER_ZBIN
        NLP  = NLP + 1
        NLP_TREE = NLP_TREE + 1
        IF (NLP.GT.NLPDIM) THEN
          CALL RE_ALLOCATE_DROPLETS(1,NM,0)
          DROPLET=>MESH(NM)%DROPLET
          ENDIF
        DR=>DROPLET(NLP)
        BLK_LOOP:  DO
         CALL RANDOM_NUMBER(RN)
         DR%Z = ZC_MIN + RN*(ZC_MAX-ZC_MIN)
         CALL RANDOM_NUMBER(RN)
         XC_MAX =
     .    (CANOPY_LENGTH-(DR%Z-CANOPY_B_H(NCT)))*TANGENT
         XC_MIN = -XC_MAX
         DR%X = XC_MIN + RN*(XC_MAX-XC_MIN)
         CALL RANDOM_NUMBER(RN)
         YC_MAX = SQRT(XC_MAX**2 - DR%X**2)
         YC_MIN = -YC_MAX
         DR%Y = YC_MIN + RN*(YC_MAX-YC_MIN)
         DR%X = X_TREE(NCT) + DR%X
         DR%Y = Y_TREE(NCT) + DR%Y
         DR%Z = Z_TREE(NCT) + DR%Z
         II = CELLSI(FLOOR((DR%X - XS)*RDXINT)) + 1.
         JJ = CELLSJ(FLOOR((DR%Y - YS)*RDYINT)) + 1.
         KK = CELLSK(FLOOR((DR%Z - ZS)*RDZINT)) + 1.
         IF (.NOT.SOLID(ICA(II,JJ,KK))) EXIT BLK_LOOP
        ENDDO BLK_LOOP
       ENDDO NISP_LOOP1
      ENDDO ZBIN_LOOP
C
C Define physical characteristics of fuel elements
C
      MASS_SUM = 0.
      NISP_LOOP2: DO I=NLP-NLP_TREE+1,NLP
       DR=>DROPLET(I)
       DR%U = 0.
       DR%V = 0.
       DR%W = 0.
       DR%TMP = LP%TMP_INITIAL
       DR%T   = 0.
       DR%IOR = 0
       DR%A_X = 0.
       DR%A_Y = 0.
       DR%A_Z = 0.
       DR%CLASS = IPC
       IF (MOD(I,LP%SAMPLING).EQ.0) THEN
          DR%SHOW = .TRUE.
       ELSE
          DR%SHOW = .FALSE.
       ENDIF
       STRATUM = MOD(I,NSTRATA) + 1
       IL = LP%IL_CDF(STRATUM)
       IU = LP%IU_CDF(STRATUM)
       CALL RANDOM_CHOICE(LP%CDF(IL:IU),LP%R_CDF(IL:IU),IU-IL,DR%R)
       DR%PWT = LP%W_CDF(STRATUM)
       MASS_SUM = MASS_SUM + DR%PWT*LP%FTPR*DR%R**3
      ENDDO NISP_LOOP2
C
C     Adjust particle weighting factor PWT so that desired 
C     MASS_PER_VOLUME is achieved
C
      DROPLET(NLP-NLP_TREE+1:NLP)%PWT =
     . DROPLET(NLP-NLP_TREE+1:NLP)%PWT*
     . LP%MASS_PER_VOLUME*CANOPY_VOLUME/MASS_SUM
C
      ENDDO TREE_LOOP
C
      END SUBROUTINE INITIALIZE_TREES
C
C
      SUBROUTINE INSERT_DROPLETS(T,NM)
C
C Insert sprinkler particles into the domain every DTSPAR seconds.
C
      REAL(EB) T,PHI_RN,RN,FLOW_RATE,DM_ADJUSTED,
     .         THETA_RN,SPHI,CPHI,MASS_SUM,
     .         STHETA,CTHETA,PWT0,TT,PP,
     .         DROPLET_SPEED,XI,YJ,ZK,SHIFT1,SHIFT2,XTMP,YTMP,ZTMP,VLEN,
     .         TRIGT1,TRIGT2,LL,UL,BIN_SIZE
      INTEGER I,J,KS,II,JJ,KK,IND,ITHETA,IPHI,IC,IL,IU,STRATUM,IPC,
     .        DROP_SUM
      INTEGER, INTENT(IN) :: NM
      TYPE (SPRINKLER_TYPE), POINTER :: S
      TYPE (SPRINKLER_HEAD_TYPE), POINTER :: SH
C
      CALL UNPACK_VAR(NM)
C
      TNOW=SECOND()
C
      SPRINKLER_INSERT_LOOP: DO KS=1,NSPR
C
      SH => SPRINKLER_HEAD(KS)
      IF (SH%MESH.NE.NM)    CYCLE SPRINKLER_INSERT_LOOP
      IF (SH%ACT_CODE.LE.0) CYCLE SPRINKLER_INSERT_LOOP
      IF (SH%T.GT.T)        CYCLE SPRINKLER_INSERT_LOOP
C
      IPC = SH%DROPLET_CLASS
      LP=>LAGRANGIAN(IPC)
      IND = SH%INDEX
      S => SPRINKLER(IND)
C
C Compute sprinkler droplet Cumulative Distribution Function (CDF_DROP)
C
      IF (S%OPERATING_PRESSURE.GT.0.) THEN
C
      IF (PRESSURE.LT.0.) THEN
         PIPE_PRESSURE = S%OPERATING_PRESSURE
         ELSE
         PIPE_PRESSURE = PRESSURE
         ENDIF
C
      S%P_FACTOR = SQRT(PIPE_PRESSURE/S%OPERATING_PRESSURE)
C
      PHI_LOOP:   DO IPHI=1,S%NPHI(3)
      THETA_LOOP: DO ITHETA=1,S%NTHETA(3)
      DM_ADJUSTED = S%DROPLET_DIAMETER(ITHETA,IPHI)*
     .   (S%OPERATING_PRESSURE/PIPE_PRESSURE)**ONTH
      CALL DROPLET_SIZE_DISTRIBUTION(DM_ADJUSTED,S%RDROP(:,ITHETA,IPHI),
     .     S%CDF_DROP(:,ITHETA,IPHI),NDC,S%GAMMA,S%SIGMA)
      BIN_SIZE = S%RDROP(NDC,ITHETA,IPHI)/REAL(NSTRATA,EB)
      STRATIFY2: DO I=1,NSTRATA
         LL = (I-1)*BIN_SIZE
         UL =  I   *BIN_SIZE
         LL_LOOP2: DO J=1,NDC
            IF (S%RDROP(J,ITHETA,IPHI).GT.LL) THEN
               IL = J-1 ; S%IL_DROP(I,ITHETA,IPHI) = J-1
               EXIT LL_LOOP2
               ENDIF
            ENDDO LL_LOOP2
         UL_LOOP2: DO J=NDC,1,-1
            IF (S%RDROP(J,ITHETA,IPHI).LE.UL) THEN
               IU = J ; S%IU_DROP(I,ITHETA,IPHI) = J
               EXIT UL_LOOP2
               ENDIF
            ENDDO UL_LOOP2
C
         S%W_DROP(I,ITHETA,IPHI) = S%CDF_DROP(IU,ITHETA,IPHI) -
     .                             S%CDF_DROP(IL,ITHETA,IPHI)
      ENDDO STRATIFY2
      ENDDO THETA_LOOP
      ENDDO PHI_LOOP
C
      S%OPERATING_PRESSURE = -1.
C
      ENDIF
C
C Direction initialization stuff
C
      TRIGT1 = ACOS(-SH%DIR(3))
      IF (SH%DIR(2).EQ.0.) THEN
      TRIGT2 = ACOS(1.)
      ELSE
      TRIGT2 = ACOS(ABS(SH%DIR(1))/
     .             SQRT(SH%DIR(1)**2+SH%DIR(2)**2))
      ENDIF
C
C Droplet insertion loop
C
      MASS_SUM = 0.
      DROP_SUM = 0
C
      DROPLET_INSERT_LOOP: DO I=1,LP%N_INSERT
C
      IF (NLP+1.GT.MAXIMUM_DROPLETS) EXIT DROPLET_INSERT_LOOP
C
      NLP = NLP+1
      IF (NLP.GT.NLPDIM) THEN
         CALL RE_ALLOCATE_DROPLETS(1,NM,0)
         DROPLET=>MESH(NM)%DROPLET
         ENDIF
      DR=>DROPLET(NLP)
C
C Set droplet temperature and phase. Phase indicates
C whether the droplet is in the air (0) or stuck on an object (not 0).
C
      DR%TMP = S%WATER_TEMPERATURE
      DR%T   = T 
      DR%IOR = 0      
      DR%A_X    = 0.
      DR%A_Y    = 0.
      DR%A_Z    = 0.
      DR%CLASS = IPC
      IF (MOD(NLP,LP%SAMPLING).EQ.0) THEN
         DR%SHOW = .TRUE.    
      ELSE
         DR%SHOW = .FALSE.
      ENDIF
C
C Randomly choose theta and phi
C
      CHOOSE_COORDS: DO
C
      CALL RANDOM_NUMBER(RN)
      II = S%I_TABLE(NINT(RN*N_TE))
      JJ = S%J_TABLE(NINT(RN*N_TE))
      CALL RANDOM_NUMBER(RN)
      THETA_RN = S%THETA(II-1) + RN*S%DTHETA
      CALL RANDOM_NUMBER(RN)
      PHI_RN   = S%PHI(JJ-1)   + RN*S%DPHI
C
C  Adjust for rotation of head by rotating about z-axis
C
      PHI_RN = PHI_RN + SH%DIR(4)
C
C  Adjust for tilt of sprinkler pipe
C
      SPHI   = SIN(PHI_RN)
      CPHI   = COS(PHI_RN)
      STHETA = SIN(THETA_RN)
      CTHETA = COS(THETA_RN)
      XTMP   = STHETA*CPHI
      YTMP   = STHETA*SPHI
      ZTMP   = -CTHETA
C
C First rotate about y-axis away from x-axis
C
      VLEN   = SQRT(XTMP**2+ZTMP**2)
      SHIFT1 = ACOS(ABS(XTMP)/VLEN)
      SELECT CASE(INT(SIGN(1._EB,ZTMP)))
      CASE (-1); IF (XTMP.LT.0) SHIFT1 = PI-SHIFT1
      CASE ( 1)
      SELECT CASE(INT(SIGN(1._EB,XTMP)))
      CASE (-1) ; SHIFT1 = SHIFT1+PI
      CASE ( 1) ; SHIFT1 = TWOPI - SHIFT1
      END SELECT
      END SELECT
C
      SHIFT1 = SHIFT1 + TRIGT1
      XTMP = VLEN * COS(SHIFT1)
      ZTMP = -VLEN * SIN(SHIFT1)
C
C Second rotate about z-axis away from x-axis
C
      VLEN   = SQRT(XTMP**2+YTMP**2)
      SHIFT1 = ACOS(ABS(XTMP)/VLEN)
      SELECT CASE(INT(SIGN(1._EB,YTMP)))
      CASE ( 1); IF (XTMP.LT.0) SHIFT1 = PI-SHIFT1
      CASE (-1)
      SELECT CASE(INT(SIGN(1._EB,XTMP)))
      CASE (-1) ; SHIFT1 = SHIFT1+PI
      CASE ( 1) ; SHIFT1 = TWOPI - SHIFT1
      END SELECT
      END SELECT
C
      SHIFT2 = TRIGT2
      SELECT CASE(INT(SIGN(1._EB,SH%DIR(1))))
      CASE (-1); IF (SH%DIR(2).GT.0) SHIFT2 = TWOPI - SHIFT2
      CASE ( 1)
      SELECT CASE(INT(SIGN(1._EB,SH%DIR(2))))
      CASE (-1) ; SHIFT2 = PI-SHIFT2
      CASE ( 1) ; SHIFT2 = SHIFT2+ PI
      END SELECT
      END SELECT
      SHIFT1=SHIFT1+SHIFT2
C
      XTMP = VLEN * COS(SHIFT1)
      YTMP = VLEN * SIN(SHIFT1)
C
      TT = THETA_RN*S%NTHETA(1)/PI
      PP = (PHI_RN-SH%DIR(4))*S%NPHI(1)/TWOPI
      JJ = INT(PP)+1
      II = INT(TT)+1
      DROPLET_SPEED = S%P_FACTOR*S%DROPLET_VELOCITY(II,JJ)
C
C Compute initial position and velocity of droplets
C
      DR%U = DROPLET_SPEED*XTMP
      DR%V = DROPLET_SPEED*YTMP
      DR%W = DROPLET_SPEED*ZTMP
      DR%X = SH%X + S%OFFSET*XTMP
      DR%Y = SH%Y + S%OFFSET*YTMP
      DR%Z = SH%Z + S%OFFSET*ZTMP
C
      IF (TWO_D) THEN
         DR%V = 0.
         DR%Y = SH%Y
         ENDIF
C
      IF (DR%X.LE.XS .OR. DR%X.GE.XF) CYCLE CHOOSE_COORDS
      IF (DR%Y.LE.YS .OR. DR%Y.GE.YF) CYCLE CHOOSE_COORDS
      IF (DR%Z.LE.ZS .OR. DR%Z.GE.ZF) CYCLE CHOOSE_COORDS
C
      XI = CELLSI(FLOOR((DR%X-XS)*RDXINT))
      YJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
      ZK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
      II = FLOOR(XI+1.)
      JJ = FLOOR(YJ+1.)
      KK = FLOOR(ZK+1.)
      IC = ICA(II,JJ,KK)
C
      IF (.NOT.SOLID(IC)) EXIT CHOOSE_COORDS
C
      ENDDO CHOOSE_COORDS
C
C Randomly choose droplet size according to Cumulative Distribution
C Function (CDF)
C
      TT = THETA_RN*S%NTHETA(3)/PI
      PP = (PHI_RN-SH%DIR(4))*S%NPHI(3)/TWOPI
      JJ = INT(PP)+1
      II = INT(TT)+1
      STRATUM = MOD(NLP,NSTRATA) + 1
      IL = S%IL_DROP(STRATUM,II,JJ)
      IU = S%IU_DROP(STRATUM,II,JJ)
      CALL RANDOM_CHOICE(S%CDF_DROP(IL:IU,II,JJ),
     .                   S%RDROP(IL:IU,II,JJ),
     .                   IU-IL,DR%R)
      DR%PWT = S%W_DROP(STRATUM,II,JJ)
C
C Sum up mass of water being introduced
C
      MASS_SUM = MASS_SUM + DR%PWT*LP%FTPR*DR%R**3
      DROP_SUM = DROP_SUM + 1
C
      ENDDO DROPLET_INSERT_LOOP
C
C Compute weighting factor for the droplets just inserted
C
      IF (DROP_SUM.GT.0) THEN
      FLOW_RATE = S%K_FACTOR*SQRT(PIPE_PRESSURE)*(LP%DENSITY/1000.)/60.
      PWT0 = FLOW_RATE*MAX(T-SH%T,0.01_EB)/MASS_SUM
      DROPLET(NLP-DROP_SUM+1:NLP)%PWT = 
     .DROPLET(NLP-DROP_SUM+1:NLP)%PWT*PWT0
      ENDIF
C
      SH%T = T    
C
      ENDDO SPRINKLER_INSERT_LOOP
C
C Update the sprinkler droplet insert clock
C
      SPINCLK(NM) = SPINCLK(NM) + DTSPAR
C
      TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
      END SUBROUTINE INSERT_DROPLETS
C
C
      SUBROUTINE INSERT_PARTICLES(T,NM)
C
C Particle initialization
C
      REAL(EB), INTENT(IN) :: T
      INTEGER, INTENT(IN) :: NM
      INTEGER I,II,JJ,KK,IIG,JJG,KKG,IW,IC,IOR,IPC,IL,IU,STRATUM
      REAL(EB) :: RN
      REAL(EB), PARAMETER :: VENT_OFFSET=0.1
C
      TNOW=SECOND()
C
      CALL UNPACK_VAR(NM)
C
C Loop through all boundary cells and assign particles if needed
C
      WALL_CELL_LOOP: DO IW=1,NWC
C
      IF (PINCLK(NM).LT.TW(IW)) CYCLE WALL_CELL_LOOP
      IF (CELL_PARTICLE_CLASS(IW).LT.0) CYCLE WALL_CELL_LOOP
      IF (UW(IW).GE.-0.0001) CYCLE WALL_CELL_LOOP
C
      II = IJKW(1,IW)
      JJ = IJKW(2,IW)
      KK = IJKW(3,IW)
      IC = ICA(II,JJ,KK)
      IF (.NOT.SOLID(IC)) CYCLE WALL_CELL_LOOP
C
      IF (NM.GT.1) THEN
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      IF (INTERPOLATED_MESH(IIG,JJG,KKG).GT.0) CYCLE WALL_CELL_LOOP
      ENDIF
C
      IOR = IJKW(4,IW)
C
      PARTICLE_INSERT_LOOP: DO I=1,NPPCW(IW)
C
      IF (NLP+1.GT.MAXIMUM_DROPLETS) EXIT PARTICLE_INSERT_LOOP
C
      NLP = NLP+1
      IF (NLP.GT.NLPDIM) THEN
         CALL RE_ALLOCATE_DROPLETS(1,NM,0)
         DROPLET=>MESH(NM)%DROPLET
         ENDIF
      DR=>DROPLET(NLP)
C
      IF (ABS(IOR).EQ.1) THEN
         IF (IOR.EQ. 1) DR%X = X(II)   + VENT_OFFSET*DX(II+1)
         IF (IOR.EQ.-1) DR%X = X(II-1) - VENT_OFFSET*DX(II-1)
         CALL RANDOM_NUMBER(RN)
         DR%Y = Y(JJ-1) + DY(JJ)*RN
         CALL RANDOM_NUMBER(RN)
         DR%Z = Z(KK-1) + DZ(KK)*RN
         ENDIF
      IF (ABS(IOR).EQ.2) THEN
         IF (IOR.EQ. 2) DR%Y = Y(JJ)   + VENT_OFFSET*DY(JJ+1)
         IF (IOR.EQ.-2) DR%Y = Y(JJ-1) - VENT_OFFSET*DY(JJ-1)
         CALL RANDOM_NUMBER(RN)
         DR%X = X(II-1) + DX(II)*RN
         CALL RANDOM_NUMBER(RN)
         DR%Z = Z(KK-1) + DZ(KK)*RN
         ENDIF
      IF (ABS(IOR).EQ.3) THEN
         IF (IOR.EQ. 3) DR%Z = Z(KK)   + VENT_OFFSET*DZ(KK+1)
         IF (IOR.EQ.-3) DR%Z = Z(KK-1) - VENT_OFFSET*DZ(KK-1)
         CALL RANDOM_NUMBER(RN)
         DR%X = X(II-1) + DX(II)*RN
         CALL RANDOM_NUMBER(RN)
         DR%Y = Y(JJ-1) + DY(JJ)*RN
         ENDIF
C
C Give particles an initial velocity (if desired)
C
      SELECT CASE(IOR)
      CASE( 1) ; DR%U = -UW(IW)
      CASE(-1) ; DR%U =  UW(IW)
      CASE( 2) ; DR%V = -UW(IW)
      CASE(-2) ; DR%V =  UW(IW)
      CASE( 3) ; DR%W = -UW(IW)
      CASE(-3) ; DR%W =  UW(IW)
      END SELECT
C
C Save the insertion time (TP) and scalar property (SP) for the particle
C
      DR%A_X    = 0.
      DR%A_Y    = 0.
      DR%A_Z    = 0.
      IPC       = CELL_PARTICLE_CLASS(IW)
      LP=>LAGRANGIAN(IPC)
      DR%CLASS  = IPC
      IF (MOD(NLP,LP%SAMPLING).EQ.0) THEN
         DR%SHOW = .TRUE.
      ELSE
         DR%SHOW = .FALSE.
      ENDIF
C
      DR%R   = 0.0
      DR%TMP = LP%TMP_INITIAL
      DR%T   = T
      DR%PWT = 1.
      DR%IOR = 0.
C
      IF (LP%DIAMETER.GT.0) THEN
      STRATUM = MOD(NLP,NSTRATA) + 1
      IL = LP%IL_CDF(STRATUM)
      IU = LP%IU_CDF(STRATUM)
      CALL RANDOM_CHOICE(LP%CDF(IL:IU),LP%R_CDF(IL:IU),IU-IL,DR%R)
      DR%PWT = LP%W_CDF(STRATUM)
      ENDIF
C
      ENDDO PARTICLE_INSERT_LOOP
      ENDDO WALL_CELL_LOOP
C
      PINCLK(NM) = PINCLK(NM) + DTPAR
C
      TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
      END SUBROUTINE INSERT_PARTICLES
C
C
      SUBROUTINE RANDOM_CHOICE(CDF,VAR,NPTS,CHOICE)
C
      INTEGER,  INTENT(IN)  :: NPTS
      REAL(EB), INTENT(IN)  :: CDF(0:NPTS),VAR(0:NPTS)
      REAL(EB), INTENT(OUT) :: CHOICE
      INTEGER   IT
      REAL(EB)  CFRAC,RN,A,B
C
      CALL RANDOM_NUMBER(RN)
C
C     Scale random number
C
      A = MINVAL(CDF)
      B = MAXVAL(CDF)
      RN = A + (B-A)*RN
C
      CDFLOOP: DO IT=1,NPTS
      IF (CDF(IT).GT.RN) THEN
      CFRAC  = (RN-CDF(IT-1))/(CDF(IT)-CDF(IT-1))
      CHOICE = VAR(IT-1) + (VAR(IT)-VAR(IT-1))*CFRAC
      EXIT CDFLOOP
      ENDIF
      ENDDO CDFLOOP
C
      END SUBROUTINE RANDOM_CHOICE
C
C
      SUBROUTINE TRACK_DROPLETS(T,NM)
C
C Momentum and heat transfer from sprinkler droplets
C
      REAL(EB), POINTER, DIMENSION(:,:,:) :: FVXS,FVYS,FVZS,
     .                            DROP_DEN,DROP_RAD,DROP_TMP
      REAL(EB), POINTER, DIMENSION(:) :: WMPUAOLD,WCPUAOLD
C
      REAL(EB) NUSSELT,KA,QREL,SFAC,TMPD,RRD,HV,
     .        RVC,RDS,RDC,UREL,VREL,WREL,TMP_G,WGT,RN,THETA_RN,
     .        FLUXMIN,OMWGT,FLUXMAX,QUSE,HD,TMPDN,RD_NEW,RD,
     .        SH_FAC,NU_FAC,RSP0,T,PR,CDRAG,MVAP,
     .        RSP1,ZK,RDT,YJ,XI,RE_D,MU_AIR,
     .        TMPS,HS,HLF,DTLF,DTSP,DTMIN,UBAR,VBAR,WBAR,
     .        BFAC,GRVT1,GRVT2,GRVT3,AUREL,AVREL,AWREL,CONST,
     .        Y_W,UVW,QREF,QRADD,QWALLC,QWALLR,TDAVE,DEN_ADD
      INTEGER ICN,I,J,K,IIN,JJN,KKN,II,JJ,KK,IIX,JJY,KKZ,IW,N,NITER,NH,
     .        IWN,IBC,IC,IYY,IWP1,IWM1,IWP2,IWM2,IWP3,IWM3,IOR_OLD,
     .        IPC,IGAS
      REAL(EB) SC_AIR,D_AIR,DHOR,SHER,XVAPC,XVAP,QCONV,QVAP,
     .         MXVAP,MDROP,RHOG,MW_RATIO,MW_DROP,
     .         CP_DROP,Z1,Z2,Z3,Z4,Z5,MGAS,ATOP,ABOT,CR2,WALL,
     .         FILMD
      INTEGER, INTENT(IN) :: NM
C
      TNOW=SECOND()
C
      CALL UNPACK_VAR(NM)
C
      IF (NLP.EQ.0) RETURN
C
      TRANSPORT: IF (CORRECTOR) THEN
C
      D_VAP = 0.
C
      DROP_DEN => WORK4
      DROP_RAD => WORK5
      DROP_TMP => WORK6
      WMPUAOLD => WALL_WORK1
      WCPUAOLD => WALL_WORK2
C
      DROP_DEN = 0. 
      DROP_TMP = 0. 
      DROP_RAD = 0. 
C
      WMPUAOLD = WMPUA
      WCPUAOLD = WCPUA
      WMPUA    = 0.
      WCPUA    = 0.
C
C Miscellaneous constants
C
      RSP0   = .001     ! Droplet radius when droplet hits solid (m)
      RSP1   = .001     ! Droplet radius when droplet drips off  (m)
      HS     = 3000.    ! Heat transfer coef from solid to drop (W/m2/K)
      HLF    = 1000.    ! Leidenfrost heat transfer coef (W/m2/K)
      PR     = 0.7      ! Prandtl number of air
      FILMD  = 0.002    ! Thickness of continous water film
C
C Definitions
C
      NU_FAC   = 0.6*PR**ONTH
      GRVT1    = -FGX(T)*GVEC(1) 
      GRVT2    = -FGY(T)*GVEC(2) 
      GRVT3    = -FGZ(T)*GVEC(3) 
      RDT      = 1./DT
      CR2      = 2.**ONTH
      RVC      = 0.
      TDAVE    = 0.
C
C Constants for new evaporation routine
C
      D_AIR  = 3.33E-5
      SC_AIR = 0.6
      SH_FAC = 0.6*SC_AIR**ONTH
C
C Loop through sprinkler droplets
C
      DROPLET_LOOP: DO I=1,NLP  
      DR=>DROPLET(I)
C
      IPC = DR%CLASS
      LP=>LAGRANGIAN(IPC)
C
      XI = CELLSI(FLOOR((DR%X-XS)*RDXINT))
      YJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
      ZK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
      II  = FLOOR(XI+1.)
      JJ  = FLOOR(YJ+1.)
      KK  = FLOOR(ZK+1.)
      IIX = FLOOR(XI+.5)
      JJY = FLOOR(YJ+.5)
      KKZ = FLOOR(ZK+.5)
C
      UBAR = AFILL(U(II-1,JJY,KKZ),U(II,JJY,KKZ),
     .             U(II-1,JJY+1,KKZ),U(II,JJY+1,KKZ),
     .             U(II-1,JJY,KKZ+1),U(II,JJY,KKZ+1),
     .             U(II-1,JJY+1,KKZ+1),U(II,JJY+1,KKZ+1),
     .             XI-II+1,YJ-JJY+.5,ZK-KKZ+.5)
      VBAR = AFILL(V(IIX,JJ-1,KKZ),V(IIX+1,JJ-1,KKZ),
     .             V(IIX,JJ,KKZ),V(IIX+1,JJ,KKZ),
     .             V(IIX,JJ-1,KKZ+1),V(IIX+1,JJ-1,KKZ+1),
     .             V(IIX,JJ,KKZ+1),V(IIX+1,JJ,KKZ+1),
     .             XI-IIX+.5,YJ-JJ+1,ZK-KKZ+.5)
      WBAR = AFILL(W(IIX,JJY,KK-1),W(IIX+1,JJY,KK-1),
     .             W(IIX,JJY+1,KK-1),W(IIX+1,JJY+1,KK-1),
     .             W(IIX,JJY,KK),W(IIX+1,JJY,KK),
     .             W(IIX,JJY+1,KK),W(IIX+1,JJY+1,KK),
     .             XI-IIX+.5,YJ-JJY+.5,ZK-KK+1)
C
      IF (LP%MASSLESS) THEN
      DR%U = UBAR
      DR%V = VBAR
      DR%W = WBAR
      DR%X = DR%X + DR%U*DT
      DR%Y = DR%Y + DR%V*DT
      DR%Z = DR%Z + DR%W*DT
      CYCLE DROPLET_LOOP
      ENDIF
C
      RVC = RDX(II)*RDY(JJ)*RDZ(KK)
      IC  = ICA(II,JJ,KK)
      RD  = DR%R
      RDS = RD*RD
      RDC = RD*RDS
      RRD = 1./RD
C
C Transfer droplet momentum to gas
C
      UREL  = DR%U - UBAR
      VREL  = DR%V - VBAR
      WREL  = DR%W - WBAR
      QREL  = MAX(1.E-6_EB,SQRT(UREL*UREL + VREL*VREL + WREL*WREL))
      TMP_G = MAX(TMPMIN,TMP(II,JJ,KK))
      RHOG  = RHO(II,JJ,KK)
      MU_AIR = MU_SPEC(0,MIN(500,NINT(0.1*TMP_G)))
      RE_D  = RHOG*QREL*2.*RD/MU_AIR
C
      DRAG_CALC: IF (DR%IOR.EQ.0 .AND. RE_D.GT.0) THEN
         CDRAG   = DRAG(RE_D)
         SFAC    = DR%PWT*RDS*PIO2*QREL*CDRAG
         DR%A_X  = SFAC*UREL*RVC
         DR%A_Y  = SFAC*VREL*RVC
         DR%A_Z  = SFAC*WREL*RVC
         IF (.NOT.LP%STATIC) THEN
         CONST   = 8.*LP%DENSITY*RD/(3.*RHOG*CDRAG*QREL)
         BFAC    = EXP(-DT/CONST)
         AUREL   = CONST*GRVT1
         AVREL   = CONST*GRVT2
         AWREL   = CONST*GRVT3
         DR%U  = UBAR + (UREL+AUREL)*BFAC - AUREL
         DR%V  = VBAR + (VREL+AVREL)*BFAC - AVREL
         DR%W  = WBAR + (WREL+AWREL)*BFAC - AWREL
         ENDIF
      ELSE
         DR%A_X  = 0.
         DR%A_Y  = 0.
         DR%A_Z  = 0.
      ENDIF DRAG_CALC
C
C Determine radiation for droplet
C
      QRADD = 0.
      IF (AVG_DROP_DEN(II,JJ,KK).GT.0.) 
     .   QRADD = QR_W(II,JJ,KK)*LP%FTPR*AVG_DROP_RAD(II,JJ,KK)/
     .           AVG_DROP_DEN(II,JJ,KK)
C
C Decide how many time steps to use in tracking droplet
C
      DTMIN = DT
      UVW = MAX( ABS(DR%U)*RDX(II),
     .           ABS(DR%V)*RDY(JJ),
     .           ABS(DR%W)*RDZ(KK) )
      IF (UVW.NE.0.) DTMIN = MIN(DTMIN,1./UVW)
C
      NITER = 1
      DTSP  = DT
      NLOOP: DO N=0,3
      IF (DTMIN.LT.DT*0.5**N) THEN
         NITER = 2**(N+1)
         DTSP =  DT*0.5**(N+1)
         ENDIF
      ENDDO NLOOP
C
      IF (LP%STATIC) THEN
         NITER = 0
         IWN   = 0
         IIN   = II
         JJN   = JJ
         KKN   = KK
         ICN   = IC
         ENDIF
C
      SUB_TIME_STEP_ITERATIONS: DO N=1,NITER
C
      IF (N.GT.1) THEN
      XI = CELLSI(FLOOR((DR%X-XS)*RDXINT))
      YJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
      ZK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
      II = FLOOR(XI+1.)
      JJ = FLOOR(YJ+1.)
      KK = FLOOR(ZK+1.)
      IC = ICA(II,JJ,KK)
      ENDIF
C
C Update droplet position
C
      DR%X = DR%X + DR%U*DTSP
      DR%Y = DR%Y + DR%V*DTSP
      DR%Z = DR%Z + DR%W*DTSP
C
C Droplet hits the floor
C
      IF (POROUS_FLOOR .AND. DR%Z.LT.ZS) THEN
         IC = ICA(II,JJ,0)
         IW = IWA(IC,3)
         IF (IV(IW).EQ.1 .AND. ACCUMULATE_WATER) THEN
            AWMPUA(IW) = AWMPUA(IW) + DR%PWT*LP%FTPR*RDC*RAW(IW)
            ENDIF
         CYCLE DROPLET_LOOP
         ENDIF
C
      IF (.NOT.POROUS_FLOOR .AND. DR%Z.LT.ZS) THEN
         IC = ICA(II,JJ,0)
         IW = IWA(IC,3)
         IF (IV(IW).EQ.1) THEN
            DR%Z = DR%Z - DR%W*DTSP
            DR%IOR = 3
            CALL RANDOM_NUMBER(RN)
            THETA_RN = TWOPI*RN
            DR%U = DROP_HORIZONTAL_VELOCITY*COS(THETA_RN)
            DR%V = DROP_HORIZONTAL_VELOCITY*SIN(THETA_RN)
            DR%W = 0.
            DR%PWT = DR%PWT*(DR%R/RSP0)**3
            DR%R = RSP0
            ENDIF
         ENDIF
C
C Where is the droplet?
C
      XI  = CELLSI(FLOOR((DR%X-XS)*RDXINT))
      YJ  = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
      ZK  = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
      IIN = FLOOR(XI+1.)
      JJN = FLOOR(YJ+1.)
      KKN = FLOOR(ZK+1.)
      IF (IIN.LT.0 .OR. IIN.GT.IBP1) CYCLE DROPLET_LOOP
      IF (JJN.LT.0 .OR. JJN.GT.JBP1) CYCLE DROPLET_LOOP
      IF (KKN.LT.0 .OR. KKN.GT.KBP1) CYCLE DROPLET_LOOP
      ICN = ICA(IIN,JJN,KKN)
C
      IF (IC.EQ.0 .OR. ICN.EQ.0) CYCLE SUB_TIME_STEP_ITERATIONS
C
      IF (.NOT.SOLID(ICN)) THEN
      IF (DR%X.LT.XS .OR. DR%X.GT.XF) CYCLE DROPLET_LOOP
      IF (DR%Y.LT.YS .OR. DR%Y.GT.YF) CYCLE DROPLET_LOOP
      IF (DR%Z.LT.ZS .OR. DR%Z.GT.ZF) CYCLE DROPLET_LOOP
      ENDIF
C
C Don't process porous droplet, if it re-emerges, make it gas phase
C
      IF (DR%IOR.EQ.4) THEN
         IF (.NOT.SOLID(ICN)) DR%IOR = 0
         CYCLE SUB_TIME_STEP_ITERATIONS
         ENDIF
C
C If droplet hits an obstacle, change its properties
C
      AIR_TO_SOLID: IF (II.NE.IIN .OR. JJ.NE.JJN .OR. KK.NE.KKN) THEN
C
      IOR_OLD = DR%IOR
C
C     Check droplet exiting old cell
C
      IWP1 = IWA(IC, 1) 
      IWM1 = IWA(IC,-1)
      IWP2 = IWA(IC, 2)
      IWM2 = IWA(IC,-2)
      IWP3 = IWA(IC, 3)
      IWM3 = IWA(IC,-3)
C
      NH = 0
      IF (KKN.GT.KK.AND.IV(IWP3).EQ.1) THEN; DR%IOR=-3; NH=NH+1; ENDIF
      IF (KKN.LT.KK.AND.IV(IWM3).EQ.1) THEN; DR%IOR= 3; NH=NH+1; ENDIF
      IF (IIN.GT.II.AND.IV(IWP1).EQ.1) THEN; DR%IOR=-1; NH=NH+1; ENDIF
      IF (IIN.LT.II.AND.IV(IWM1).EQ.1) THEN; DR%IOR= 1; NH=NH+1; ENDIF
      IF (JJN.GT.JJ.AND.IV(IWP2).EQ.1) THEN; DR%IOR=-2; NH=NH+1; ENDIF
      IF (JJN.LT.JJ.AND.IV(IWM2).EQ.1) THEN; DR%IOR= 2; NH=NH+1; ENDIF
C
      IWN = IWA(IC,-DR%IOR)
C
      IF (IWN.EQ.0) THEN
C
      IWP1 = IWA(ICN, 1)
      IWM1 = IWA(ICN,-1)
      IWP2 = IWA(ICN, 2)
      IWM2 = IWA(ICN,-2)
      IWP3 = IWA(ICN, 3)
      IWM3 = IWA(ICN,-3)
C
      NH = 0
      IF (KKN.GT.KK.AND.IV(IWM3).EQ.1) THEN; DR%IOR=-3; NH=NH+1; ENDIF
      IF (KKN.LT.KK.AND.IV(IWP3).EQ.1) THEN; DR%IOR= 3; NH=NH+1; ENDIF
      IF (IIN.GT.II.AND.IV(IWM1).EQ.1) THEN; DR%IOR=-1; NH=NH+1; ENDIF
      IF (IIN.LT.II.AND.IV(IWP1).EQ.1) THEN; DR%IOR= 1; NH=NH+1; ENDIF
      IF (JJN.GT.JJ.AND.IV(IWM2).EQ.1) THEN; DR%IOR=-2; NH=NH+1; ENDIF
      IF (JJN.LT.JJ.AND.IV(IWP2).EQ.1) THEN; DR%IOR= 2; NH=NH+1; ENDIF
C
      IWN = IWA(ICN,DR%IOR)
C
      ENDIF
C
C Check if droplet has crossed no solid planes or too many
C
      HIT_SOLID: IF (NH.NE.0) THEN
C
      IF (IWN.EQ.0) CYCLE SUB_TIME_STEP_ITERATIONS
C
C Move droplet out of solid and place it near surface
C
      SELECT CASE(DR%IOR)
      CASE( 1) ; DR%X = XW(IWN) + 0.1*DX(II)
      CASE(-1) ; DR%X = XW(IWN) - 0.1*DX(II)
      CASE( 2) ; DR%Y = YW(IWN) + 0.1*DY(JJ)
      CASE(-2) ; DR%Y = YW(IWN) - 0.1*DY(JJ)
      CASE( 3) ; DR%Z = ZW(IWN) + 0.1*DZ(KK)
      CASE(-3) ; DR%Z = ZW(IWN) - 0.1*DZ(KK)
      END SELECT
C
      SELECT CASE(ABS(DR%IOR))
      CASE(1)
      XI = CELLSI(FLOOR((DR%X-XS)*RDXINT))
      IIN = FLOOR(XI+1.)
      CASE(2)
      YJ = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
      JJN = FLOOR(YJ+1.)
      CASE(3)
      ZK = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
      KKN = FLOOR(ZK+1.)
      END SELECT
C
      ICN = ICA(IIN,JJN,KKN)
C
      IF (IOR_OLD.EQ.DR%IOR) CYCLE SUB_TIME_STEP_ITERATIONS
C
C Reweight the droplets
C
      IF (DR%IOR.NE.0) THEN
         DR%PWT = DR%PWT*(DR%R/RSP0)**3
         DR%R = RSP0
         ENDIF
C
C Choose a direction for the droplets to move
C
      DIRECTION: SELECT CASE(ABS(DR%IOR))
C
      CASE (1:2) DIRECTION  
      DR%U = 0.
      DR%V = 0.
      DR%W = -DROP_VERTICAL_VELOCITY 
C
      CASE(3) DIRECTION
      IBC = IJKW(5,IWN)
      CALL RANDOM_NUMBER(RN)
      IF (RN.LT.POROS(IBC)) THEN
         DR%IOR = 4
         DR%Z = ZW(IWN) - 0.1*DZ(KK)
         DR%U = 0.
         DR%V = 0.
         DR%W = -DROP_VERTICAL_VELOCITY
         CYCLE DROPLET_LOOP
         ELSE
         CALL RANDOM_NUMBER(RN)
         THETA_RN = TWOPI*RN
         DR%U = DROP_HORIZONTAL_VELOCITY*COS(THETA_RN)
         DR%V = DROP_HORIZONTAL_VELOCITY*SIN(THETA_RN)
         DR%W = 0.
         ENDIF
C
      END SELECT DIRECTION
C
      ENDIF HIT_SOLID
      ENDIF AIR_TO_SOLID 
C
C Check if attached droplets are still attached
C
      SELECT CASE(DR%IOR)
      CASE( 1)
         IW = IWA(ICN,-1)
         IF (IV(IW).NE.1) DR%U = -DROP_HORIZONTAL_VELOCITY
      CASE(-1)
         IW = IWA(ICN, 1)
         IF (IV(IW).NE.1) DR%U =  DROP_HORIZONTAL_VELOCITY
      CASE( 2)
         IW = IWA(ICN,-2)
         IF (IV(IW).NE.1) DR%V = -DROP_HORIZONTAL_VELOCITY
      CASE(-2)
         IW = IWA(ICN, 2)
         IF (IV(IW).NE.1) DR%V =  DROP_HORIZONTAL_VELOCITY
      CASE( 3)
         IW = IWA(ICN,-3)
         IF (IV(IW).NE.1) THEN
            DR%U = -DR%U
            DR%V = -DR%V
            DR%Z = DR%Z - 0.2*DZ(KK)
            ENDIF
      CASE(-3)
         IF (.NOT.SOLID(ICN)) DR%IOR = 0
      END SELECT
C
      IF (DR%IOR.NE.0 .AND. IV(IW).NE.1) THEN
          DR%IOR = 0
          IWN = 0
      ELSE
          IWN = IWA(ICN,-DR%IOR)
      ENDIF
C
      ENDDO SUB_TIME_STEP_ITERATIONS
C
C
C Compute mass and energy transfer between droplet and gas 
C
      IF (LP%SPECIES_INDEX.EQ.0) CYCLE DROPLET_LOOP
C
      HEAT_AND_MASS_TRANSFER: IF (DR%IOR.NE.4 .AND. DR%R.GT.0.) THEN
C
      IGAS = LP%SPECIES_INDEX
      MW_DROP = MWN(IGAS)
C
      MW_RATIO = MWN(0)/MW_DROP
      CP_DROP  = GAMMA*R0/((GAMMA-1.)*MW_DROP)
C
      RD   = DR%R
      RDS  = RD*RD
      RDC  = RD*RDS
      RRD  = 1./RD
      TMPD = DR%TMP
C
      IF (IGAS==IWATER) THEN 
         HV = LP%H_V_0 + 2400.*(LP%TMP_V-TMPD) ! Incropera, 4th, p 846
      ELSE
         HV = LP%H_V_0
      ENDIF
C
      DHOR = HV*MW_DROP/R0
C
      XVAP  = MIN(1._EB,EXP(DHOR*(1./LP%TMP_V-1./TMPD)))
      XVAP  = XVAP/(MW_RATIO + (1.-MW_RATIO)*XVAP)
      TMPDN = TMPD
      WALL  = 1
      ABOT  = 0.
      ATOP  = 4.*PI*RDS
      QUSE  = 0.
      QWALLR= 0.
      QWALLC= 0.
      MDROP = LP%FTPR*RDC       ! Mass of drop
      WGT   = DR%PWT
      RD_NEW = RD
      MVAP  = 0.
C
C Set variables for heat transfer on solid
C
      SET_SOLID: IF (DR%IOR.NE.0 .AND. IWN.GT.0) THEN
C
      IW = IWN
C
C Droplet dimension adjust assuming hemisphere
C
      WALL = 0.5
      RD   = CR2*RD
      RD_NEW = RD
      RDS  = RD*RD
      RDC  = RD*RDS
      RRD  = 1./RD
      ABOT = PI*RDS
      ATOP = 2.*PI*RDS
      IF (WMPUAOLD(IW)/LP%DENSITY.GT.FILMD) THEN
         ABOT = MDROP/WMPUAOLD(IW)
         ATOP = ABOT
         ENDIF
C
      IIN  = IJKW(6,IW)
      JJN  = IJKW(7,IW)
      KKN  = IJKW(8,IW)
      IBC  = IJKW(5,IW)
      RVC  = RDX(IIN)*RDY(JJN)*RDZ(KKN)
C
      TMPS   = TMP_F(IW)
      TMP_G  = TMP(IIN,JJN,KKN)
      RHOG   = RHO(IIN,JJN,KKN)
      RE_D   = RHOG*DROP_HORIZONTAL_VELOCITY*2.*RD/MU_AIR
C
      ELSE SET_SOLID
C
      TMP_G = TMP(IIN,JJN,KKN)
      RHOG = RHO(IIN,JJN,KKN)
C
      ENDIF SET_SOLID  ! Solid Heat Transfer
C
      KA      = CPOPR*MU_AIR
      NUSSELT = 2. + NU_FAC*SQRT(RE_D)
      HD      = 0.5*NUSSELT*KA*RRD
      MDROP   = WALL*LP%FTPR*RDC      ! Mass of drop
      MGAS    = RHOG/RVC           ! Mass of gas cell
      QCONV   = CP_GAMMA*LP%C_P*MDROP*MGAS/
     .          (LP%C_P*MDROP*WGT+CP_GAMMA*MGAS) ! Max conv. per drop
      QCONV   = MIN(QCONV,HD*ATOP*DT)*(TMP_G-TMPD)
      QCONV   = MAX(QCONV,LP%C_P*MDROP*(TMPM-TMPD))
      QRADD   = QRADD*RDS*DT
C
      SHER    = 2. + SH_FAC*SQRT(RE_D)
      Y_W     = YY(IIN,JJN,KKN,IGAS)
C
      IF (MIXTURE_FRACTION) THEN
      IF (IGAS.EQ.IWATER) THEN  ! Add water vapor from combustion
         IYY   = MAX(0,MIN(10000,NINT(YY(IIN,JJN,KKN,IFUEL)*10000.)))
         XVAPC = Y_W + Y_STATE(IYY,4)*(1.-Y_W)
         ELSE
         XVAPC = Y_W
         ENDIF
      ELSE
         XVAPC = Y_W
      ENDIF
C
      IF (WALL.LT.1.) THEN
         DTLF = 2.*(LP%TMP_V-TMPM)
         IF (TMPS-LP%TMP_V.GT.DTLF) THEN  ! Leidenfrost droplet
            QWALLC = ABOT*DT*HLF*(TMPS-TMPD)
         ELSE
            QWALLC = ABOT*DT*HS*(TMPS-TMPD)
         ENDIF
         QWALLR = ABOT*DT*QRAD(IW)
      ENDIF
C
      QUSE = QWALLR+QWALLC
      QREF = QWALLR+QWALLC+QCONV+QRADD
C
C Mass loss due to vapor pressure difference
C
      VAPOR_PRESSURE:  IF (XVAP.GT.XVAPC.AND.TMPD.GT.LP%TMP_MELT) THEN
C
      Z1 = LP%C_P*(TMPD-LP%TMP_MELT)
      Z2 = XVAP - 1.
      Z3 = XVAP*(QREF*WGT-MGAS*(HV+Z1))+Z1*(MGAS*XVAPC+MDROP*WGT*Z2)
      Z4 = 2*WGT*(HV*XVAP+Z1*Z2)
      Z5 = 2*MGAS*(QREF*XVAP+MDROP*(XVAP-XVAPC)*Z1)*Z4+Z3**2
C
      IF (Z5.LT.0.) THEN
         MXVAP=MDROP
      ELSE
         MXVAP=MAX(0._EB,(Z3+SQRT(Z5))/Z4)
      ENDIF
C
      Z1     = LP%C_P*(TMPD-LP%TMP_MELT)
      IF (XVAP.EQ.1.) THEN
        MVAP = MDROP
        ELSE
        MVAP = SHER*RHOG*D_AIR*WALL*TWOPI*RD*LOG((XVAPC-1.)/(XVAP-1.))*
     .         DT
        ENDIF
      MVAP   = MAX(MIN(MVAP,MDROP,MXVAP),0._EB)
      QVAP   = HV*MVAP
      RD_NEW = (MAX(0._EB,RDC-MVAP/LP%FTPR/WALL))**ONTH
C
      IF (QVAP.GE.0. .AND. QVAP.GT.QREF+(MDROP-MVAP)*Z1) THEN
         QVAP=QREF+(MDROP-MVAP)*Z1
         MVAP=MAX(QVAP/HV,0._EB)
         ENDIF
      RD_NEW=(MAX(0._EB,RDC-MVAP/LP%FTPR/WALL))**ONTH
C
      IF (MVAP.LT.MDROP) TMPDN = TMPD + (QREF-QVAP)/
     .                                  (LP%C_P*(MDROP-MVAP))
C
      IF (TMPDN.GT.LP%TMP_V) THEN !Adjust to avoid superheating drop
         QVAP = MIN(QVAP+LP%C_P*(MDROP-MVAP)*(TMPDN-LP%TMP_V),
     .              HV*MDROP)
         TDAVE = LP%TMP_V
         TMPDN = LP%TMP_V
         MVAP  = QVAP/HV
         RD_NEW = (MAX(0._EB,RDC-MVAP/LP%FTPR/WALL))**ONTH
      ELSE
         TDAVE=0.5*(TMPDN+TMPD)
      ENDIF
C
      QVAP = QVAP+LP%C_P*(MDROP-MVAP)*(TMPDN-TMPD)
C
      ELSE VAPOR_PRESSURE ! Heating due to convection
C
      MVAP  = 0.
      QVAP  = QREF
      TMPDN = TMPD+QREF/(MDROP*LP%C_P)
      IF (TMPDN.GT.LP%TMP_V) THEN
         TMPDN  = LP%TMP_V  !Set drop temp to boiling
         QVAP   = MDROP*LP%C_P*(LP%TMP_V-TMPD)
         MVAP   = MIN((QREF-QVAP)/HV,MDROP)
         QVAP   = MVAP*HV+QVAP
         RD_NEW = (MAX(0._EB,RDC-MVAP/LP%FTPR/WALL))**ONTH
         TDAVE  = LP%TMP_V
      ELSEIF (TMPDN.LE.LP%TMP_MELT) THEN
         TMPDN  = LP%TMP_MELT  ! Set drop to freezing
         QVAP   = MDROP*LP%C_P*(TMPDN-TMPD)
      ENDIF
C
      ENDIF VAPOR_PRESSURE
C
      IF (QREF.GT.QVAP.AND.(ABS(QCONV)+ABS(QWALLC)).NE.0) THEN
         QCONV =QCONV -ABS(QCONV) /(ABS(QCONV)+ABS(QWALLC))*(QREF-QVAP)
         QWALLC=QWALLC-ABS(QWALLC)/(ABS(QCONV)+ABS(QWALLC))*(QREF-QVAP)
         ENDIF
C
      QUSE = QWALLR+QWALLC
C
      IF (WALL.LT.1.) THEN
         RDC       = RD_NEW**3
         WCPUA(IW) = WCPUA(IW) + WGT*RDT*QUSE*RAW(IW)
         WMPUA(IW) = WMPUA(IW) + WGT*WALL*LP%FTPR*RDC*RAW(IW)
         RD_NEW    = RD_NEW/CR2
         ENDIF
C
      DR%R  = RD_NEW
      RDS   = RD_NEW*RD_NEW
      RDC   = RD_NEW*RDS
      MVAP  = MVAP*WGT
C
C Decrease temperature due to droplet heating/vaporization
C This is energy to raise droplet vapor to new gas temperature
C
      TMP(IIN,JJN,KKN)=(CP_GAMMA*MGAS*TMP_G+CP_DROP*MVAP*TDAVE-
     .                  QCONV*WGT)/(CP_GAMMA*MGAS+CP_DROP*MVAP)
      TMP(IIN,JJN,KKN) = MIN(TMPMAX,MAX(TMPMIN,TMP(IIN,JJN,KKN)))
C
C Save water vapor production rate for diverence expression
C
      D_VAP(IIN,JJN,KKN) = D_VAP(IIN,JJN,KKN) +
     .    R0*RDT*(
     .    RHOG*((1-Y_W)/MWN(0)+Y_W/MW_DROP)*(TMP(IIN,JJN,KKN)-TMP_G) +
     .    RVC*MVAP*TMP(IIN,JJN,KKN)/MW_DROP) 
C
C Adjust mass of evaporated liquid to account for different Heat of
C Combustion between fuel droplet and gas
C 
      MVAP = LP%ADJUST_EVAPORATION*MVAP
C
C Add water vapor or fuel gas to the cell
C
      YY(IIN,JJN,KKN,IGAS)= MIN(1._EB,(MVAP+MGAS*Y_W)/(MVAP+MGAS))
C
C Add new mass from vaporized water droplet 
C
      RHO(IIN,JJN,KKN) = RHO(IIN,JJN,KKN) + MVAP*RVC
C
      DR%TMP = TMPDN
C
      IF (DR%R.LE.0.) CYCLE DROPLET_LOOP
C
      ENDIF HEAT_AND_MASS_TRANSFER
C
C Assign water mass to the cell for airborne drops
C
      IF (DR%IOR.EQ.0) THEN
      DEN_ADD = WGT*LP%FTPR*RDC*RVC
      DROP_DEN(IIN,JJN,KKN) = DROP_DEN(IIN,JJN,KKN) + DEN_ADD
      DROP_TMP(IIN,JJN,KKN) = DROP_TMP(IIN,JJN,KKN) + DEN_ADD*TMPDN
      DROP_RAD(IIN,JJN,KKN) = DROP_RAD(IIN,JJN,KKN) + DEN_ADD*RD_NEW
      ENDIF
C
      ENDDO DROPLET_LOOP
C
C Weight the new water mass array
C
      WGT   = .5
      OMWGT = 1.-WGT
C
      DROP_RAD = DROP_RAD/(DROP_DEN+1.E-10)
      DROP_TMP = DROP_TMP/(DROP_DEN+1.E-10)
C
      AVG_DROP_RAD = (AVG_DROP_DEN*AVG_DROP_RAD + DROP_DEN*DROP_RAD)
     .              /(AVG_DROP_DEN + DROP_DEN + 1.0E-10)
      AVG_DROP_TMP = (AVG_DROP_DEN*AVG_DROP_TMP + DROP_DEN*DROP_TMP)
     .              /(AVG_DROP_DEN + DROP_DEN + 1.0E-10)
      AVG_DROP_TMP = MAX(TMPM,AVG_DROP_TMP)
      AVG_DROP_DEN = WGT*AVG_DROP_DEN + OMWGT*DROP_DEN
      WHERE (AVG_DROP_DEN.LT.0.0001 .AND. DROP_DEN.EQ.0.)
     .   AVG_DROP_DEN = 0.0
C
      WMPUA = WGT*WMPUAOLD + OMWGT*WMPUA
      WCPUA = WGT*WCPUAOLD + OMWGT*WCPUA
C
C     Remove out-of-bounds droplets
C
      CALL REMOVE_DROPLETS(T,NM)
C
      ENDIF TRANSPORT
C
C Add droplet momentum as a force term in momentum equation
C
      FVXS  => WORK1
      FVYS  => WORK2
      FVZS  => WORK3
      FVXS  = 0.
      FVYS  = 0.
      FVZS  = 0.
C
      SUM_MOMENTUM_LOOP: DO I=1,NLP
C
      DR=>DROPLET(I)
      IPC=DR%CLASS
      LP=>LAGRANGIAN(IPC)
      IF (DR%IOR.NE.0) CYCLE SUM_MOMENTUM_LOOP
      IF (LP%MASSLESS) CYCLE SUM_MOMENTUM_LOOP
C
      XI  = CELLSI(FLOOR((DR%X-XS)*RDXINT))
      YJ  = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
      ZK  = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
      II  = FLOOR(XI+1.)
      JJ  = FLOOR(YJ+1.)
      KK  = FLOOR(ZK+1.)
      IF (SOLID(ICA(II,JJ,KK))) CYCLE SUM_MOMENTUM_LOOP
      IIX = FLOOR(XI+.5)
      JJY = FLOOR(YJ+.5)
      KKZ = FLOOR(ZK+.5)
C
      IC = ICA(IIX,JJ,KK)
      IW = IWA(IC,1)
      IF (IV(IW).EQ.0) FVXS(IIX,JJ,KK) = FVXS(IIX,JJ,KK) - DR%A_X
C
      IC = ICA(II,JJY,KK)
      IW = IWA(IC,2) 
      IF (IV(IW).EQ.0) FVYS(II,JJY,KK) = FVYS(II,JJY,KK) - DR%A_Y
C
      IC = ICA(II,JJ,KKZ)
      IW = IWA(IC,3) 
      IF (IV(IW).EQ.0) FVZS(II,JJ,KKZ) = FVZS(II,JJ,KKZ) - DR%A_Z
C
      ENDDO SUM_MOMENTUM_LOOP
C
      FLUXMAX =  100.
      FLUXMIN = -FLUXMAX
C
      DO K=0,KBAR
      DO J=0,JBAR
      DO I=0,IBAR
      IF (FVXS(I,J,K).GT.FLUXMAX) FVXS(I,J,K) = FLUXMAX
      IF (FVXS(I,J,K).LT.FLUXMIN) FVXS(I,J,K) = FLUXMIN
      FVX(I,J,K) = FVX(I,J,K) + FVXS(I,J,K)
      IF (FVYS(I,J,K).GT.FLUXMAX) FVYS(I,J,K) = FLUXMAX
      IF (FVYS(I,J,K).LT.FLUXMIN) FVYS(I,J,K) = FLUXMIN
      FVY(I,J,K) = FVY(I,J,K) + FVYS(I,J,K)
      IF (FVZS(I,J,K).GT.FLUXMAX) FVZS(I,J,K) = FLUXMAX
      IF (FVZS(I,J,K).LT.FLUXMIN) FVZS(I,J,K) = FLUXMIN
      FVZ(I,J,K) = FVZ(I,J,K) + FVZS(I,J,K)
      ENDDO
      ENDDO
      ENDDO
C
      TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
      END SUBROUTINE TRACK_DROPLETS
C
C
C
      SUBROUTINE REMOVE_DROPLETS(T,NM)
C
C Remove droplets that do not lie in any mesh
C
      INTEGER :: IKILL,I,NM,IPC
      REAL(EB) :: T
      REAL(EB), PARAMETER :: RDMIN=1.0E-6   ! Minimum droplet size (m)
C
      IKILL = 0
      DROP_LOOP: DO I=1,NLP
      WEED_LOOP: DO
      DR=>DROPLET(I) ; IPC=DR%CLASS ; LP=>LAGRANGIAN(IPC)
      IF (I.GT.NLP-IKILL) EXIT DROP_LOOP
      IF (DR%R.LE.RDMIN .AND. .NOT.LP%MASSLESS) THEN 
         CALL REPLACE ; CYCLE WEED_LOOP ; ENDIF
      IF (T-DR%T.GT.LP%LIFETIME) THEN
         CALL REPLACE ; CYCLE WEED_LOOP ; ENDIF
      IF (DR%X.GT.XS .AND. DR%X.LT.XF .AND.
     .    DR%Y.GT.YS .AND. DR%Y.LT.YF .AND.
     .    DR%Z.GT.ZS .AND. DR%Z.LT.ZF)
     .    CYCLE DROP_LOOP
      CALL REPLACE
      ENDDO WEED_LOOP
      ENDDO DROP_LOOP
C
      NLP = NLP - IKILL
C
      CONTAINS
C
C
      SUBROUTINE REPLACE
C
      INTEGER OM,NOM
      TYPE (MESH_TYPE), POINTER :: M
      TYPE (OMESH_TYPE), POINTER :: M2
C
      NOM = 0
      SEARCH_LOOP: DO OM=1,NMESHES
      IF (NIC(NM,OM).EQ.0) CYCLE SEARCH_LOOP
      M=>MESH(OM)
      IF (DR%X.GT.M%XS .AND. DR%X.LT.M%XF .AND.
     .    DR%Y.GT.M%YS .AND. DR%Y.LT.M%YF .AND.
     .    DR%Z.GT.M%ZS .AND. DR%Z.LT.M%ZF) THEN
          NOM = OM
          EXIT SEARCH_LOOP
          ENDIF
      ENDDO SEARCH_LOOP
C
      IF (NOM.NE.0) THEN
      M2=>MESH(NM)%OMESH(NOM)
      M2%N_DROP_ORPHANS = M2%N_DROP_ORPHANS + 1
      IF (M2%N_DROP_ORPHANS.GT.M2%N_DROP_ORPHANS_DIM) 
     .    CALL RE_ALLOCATE_DROPLETS(2,NM,NOM)
      M2%DROPLET(M2%N_DROP_ORPHANS) = DROPLET(I)
      ENDIF
C
      DROPLET(I) = DROPLET(NLP-IKILL)
      IKILL = IKILL + 1
C
      END SUBROUTINE REPLACE
C
C
      END SUBROUTINE REMOVE_DROPLETS
C
C
      SUBROUTINE CHECK_SPRINKLERS(T,NM)
C
C Checks if sprinklers are to go off
C
      REAL(EB) TMP_G,VEL2,VEL,RHS,VELSR,T,WATER_VOL_FRAC
      INTEGER II,JJ,KK,N,IND
      REAL(EB), PARAMETER :: C_DIMARZO = 6.0E6
      INTEGER, INTENT(IN) :: NM
      TYPE (SPRINKLER_HEAD_TYPE), POINTER :: SH
      TYPE (SPRINKLER_TYPE), POINTER :: S
C
      TNOW=SECOND()
C
      CALL UNPACK_VAR(NM)
C
      SPRINKLER_LOOP: DO N=1,NSPR
C
      SH => SPRINKLER_HEAD(N)
      IF (NM.NE.SH%MESH) CYCLE SPRINKLER_LOOP
C
      IF (T.GT.SH%T_DEACT) THEN
         SH%ACT_CODE = -1 
         CYCLE SPRINKLER_LOOP
         ENDIF
C
      IF (SH%ACT_CODE.NE.0) CYCLE SPRINKLER_LOOP
C
C Determine the temperature, velocity, water mass fraction at sprinkler
C
      II = SH%I
      JJ = SH%J
      KK = SH%K
      TMP_G = TMP(II,JJ,KK)
      VEL2  = 0.25*( (U(II,JJ,KK)+U(II-1,JJ,KK))**2 +
     .               (V(II,JJ,KK)+V(II,JJ-1,KK))**2 +
     .               (W(II,JJ,KK)+W(II,JJ,KK-1))**2 )
      VEL   = SQRT(VEL2)
      VELSR = SQRT(VEL)
      WATER_VOL_FRAC = AVG_DROP_DEN(II,JJ,KK)/LAGRANGIAN(NPC)%DENSITY
C
C Update the link temperature TMP_L using a Runge-Kutta scheme
C
      IND = SH%INDEX
C
      S => SPRINKLER(IND)
C
      IF (PREDICTOR) THEN
      RHS      = ( VELSR     *(TMP_G-SH%TMP_L) 
     .            -S%C_FACTOR*(SH%TMP_L -TMPA) 
     .            -C_DIMARZO*VEL*WATER_VOL_FRAC )/S%RTI
      SH%TMP_L_S = MAX(TMPA,SH%TMP_L + DT*RHS)
      ELSE
      RHS     = ( VELSR     *(TMP_G-SH%TMP_L_S) 
     .           -S%C_FACTOR*(SH%TMP_L_S -TMPA)
     .           -C_DIMARZO*VEL*WATER_VOL_FRAC )/S%RTI
      SH%TMP_L = MAX(TMPA,.5*(SH%TMP_L + SH%TMP_L_S + DT*RHS))
      ENDIF
C
      ENDDO SPRINKLER_LOOP
C
      TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
      END SUBROUTINE CHECK_SPRINKLERS
C
C
      SUBROUTINE CHECK_HEAT_DETECTORS(T,NM)
C
C Checks if heat detectors are to go off
C
      REAL(EB) TMP_G,VEL2,VEL,RHS,VELSR,T
      INTEGER II,JJ,KK,N
      INTEGER, INTENT(IN) :: NM
      TYPE (HEAT_DETECTOR_TYPE), POINTER :: HD
C
      TNOW=SECOND()
C
      CALL UNPACK_VAR(NM)
C
      DETECTOR_LOOP: DO N=1,NHD
C
      HD => HEAT_DETECTOR(N)
C
      IF (HD%MESH.NE.NM)         CYCLE DETECTOR_LOOP
      IF (HD%TMP_ACT.GT.100000.) CYCLE DETECTOR_LOOP
C
C Determine the temperature, velocity at heat detector
C
      II = HD%I
      JJ = HD%J
      KK = HD%K
      TMP_G = TMP(II,JJ,KK)
      VEL2  = 0.25*( (U(II,JJ,KK)+U(II-1,JJ,KK))**2 +
     .               (V(II,JJ,KK)+V(II,JJ-1,KK))**2 +
     .               (W(II,JJ,KK)+W(II,JJ,KK-1))**2 )
      VEL   = SQRT(VEL2)
      VELSR = SQRT(VEL)
C
C Update the link temperature TMPL using a Runge-Kutta scheme
C
      IF (PREDICTOR) THEN
      RHS         = VELSR*(TMP_G-HD%TMP_L)/HD%RTI
      HD%TMP_L_S  = HD%TMP_L + DT*RHS
      ELSE
      RHS         = VELSR*(TMP_G-HD%TMP_L_S)/HD%RTI
      HD%TMP_L    = 0.5*(HD%TMP_L + HD%TMP_L_S + DT*RHS)
      ENDIF
C
      ENDDO DETECTOR_LOOP
C
      TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
      END SUBROUTINE CHECK_HEAT_DETECTORS
C
C
      SUBROUTINE SMOKE_DETECTORS(T,NM)
C
C     01/20/2004 Smoke Detector modeling routine
C
      REAL(EB) YSMOKE_E1,VEL2,VEL,ARHS,T,XI,YJ,ZK,DT_CURVE
     .        ,SUM,TAS,TBS,Y_ODE,T_TAR1,T_TAR2
      REAL(EB) YYHAT,IFAC,Y_STATE_INT,Y_EXTRA
      REAL(EB) OLD_SS,FUNC_A,FUNC_B,FUNC_XDX,UPDATE_SS,XDX,DEL,EPS
     .        ,C_ELECTRICAL
      INTEGER  IT_NUB,N_N,JMAX,NN,N,ITSX,IXI,IYY1,IYY2
      REAL(EB),POINTER, DIMENSION(:,:,:) :: FF_SMOKE
      INTEGER  N_IERT,I_XDX,II,I,J,K,I1,J1,K1
      INTEGER, INTENT(IN) :: NM
      TYPE (SMOKE_DETECTOR_TYPE),POINTER :: SD
C
      TNOW=SECOND()
C
      CALL UNPACK_VAR(NM)
C
      FF_SMOKE => WORK1
C
C Smoke mass fraction:
C
      DO K=0,KBP1
      DO J=0,JBP1
      DO I=0,IBP1
      YYHAT = YY(I,J,K,IFUEL)
      YYHAT = MIN(1._EB,MAX(0._EB,YYHAT))
      IYY1  = FLOOR(YYHAT*10000.)
      IFAC  = YYHAT*10000. - IYY1
      IYY2  = IYY1+1.0
      IF(IYY2.GT.10000) IYY2 = 10000
      Y_EXTRA = 0.
      DO NN=1,NSPEC
      IF (NN.NE.IFUEL) Y_EXTRA = Y_EXTRA + YY(I,J,K,NN)
      ENDDO
      Y_STATE_INT  = (1.-IFAC)*Y_STATE(IYY1,8) +
     .                   IFAC *Y_STATE(IYY2,8)
      FF_SMOKE(I,J,K) = Y_STATE_INT*RHO(I,J,K)*(1.-Y_EXTRA)*1.E6
      ENDDO
      ENDDO
      ENDDO
C
C Modeling of smoke detectors   
C
      SMOKE_DETECTOR_LOOP: DO N=1,NSD

      SD => SMOKE_DETECTOR(N)
C    
      IF (SD%MESH.NE.NM)            CYCLE SMOKE_DETECTOR_LOOP
cc    IF (SD%YSMOKE_ACT.GT.100000.) CYCLE SMOKE_DETECTOR_LOOP
C
C Determine the Mass fraction of smoke, velocity at
C the smoke detector locations.
C
      XI = CELLSI(FLOOR((SD%X-XS)*RDXINT)) + .5
      YJ = CELLSJ(FLOOR((SD%Y-YS)*RDYINT)) + .5
      ZK = CELLSK(FLOOR((SD%Z-ZS)*RDZINT)) + .5
      I1 = FLOOR(XI)
      J1 = FLOOR(YJ)
      K1 = FLOOR(ZK)
C
      YSMOKE_E1 = AFILL(FF_SMOKE(I1,J1,K1),FF_SMOKE(I1+1,J1,K1),
     .           FF_SMOKE(I1,J1+1,K1),FF_SMOKE(I1+1,J1+1,K1),
     .           FF_SMOKE(I1,J1,K1+1),FF_SMOKE(I1+1,J1,K1+1),
     .           FF_SMOKE(I1,J1+1,K1+1),FF_SMOKE(I1+1,J1+1,K1+1),
     .           XI-I1,YJ-J1,ZK-K1)
C
      VEL2  = 0.25*(  (U(I1,J1,K1)+U(I1-1,J1,K1))**2 
     .              + (V(I1,J1,K1)+V(I1,J1-1,K1))**2 
     .              + (W(I1,J1,K1)+W(I1,J1,K1-1))**2 )
      VEL   = SQRT(VEL2)
      VEL   = MAX(VEL,1.0E-10_EB)
C
C Converting the mass fractoin of smoke into obscuration(%/m) 
C at the detector location. 
C
      PREDICTOR_IF: IF (PREDICTOR) THEN
C
      N_IERT= FLOOR(T*10)
C
c Converting the mass fraction of smoke into obscuration(%/m)  
c
      SD%YSMOK_E(N_IERT)  = YSMOKE_E1
      SD%TAR_CURVE(N_IERT)= SD%ALPHA_C*VEL**SD%BETA_C
      Y_ODE  = (MASS_EXTINCTION_COEFFICIENT/1000.)*
     .         (SD%YSMOK_E(N_IERT)/1000.)
      SD%YSMOKE_OUT  = (1-EXP(-Y_ODE))*100
C
C Update the mass fraction of smoke at the inside of detector.
C
      DT_CURVE   = SD%ALPHA_E*VEL**SD%BETA_E
      OLD_SS     = -1.0E30
      JMAX       = 3
      EPS        = 1.0E-6
      UPDATE_SS  = 0.0
      TAS        = 0.0
      TBS        = T
      CALL TAU_SUB(T_TAR1,DT_CURVE,T,N,NM,N_IERT)
      ARHS       = EXP(-T_TAR1)
C
      TBS_GT_DT_CURVE: IF (TBS.GT.DT_CURVE) THEN
C
      INTLOOP: DO J=1,JMAX
C
      J_EQ_1: IF (J.EQ.1) THEN
C
      IF (TAS.EQ.0.0) FUNC_A = 0.0
      IF (TBS.LT.DT_CURVE) CYCLE INTLOOP
      N_N = FLOOR((TBS - DT_CURVE)*10)
C
      DO II=1,N_IERT
      IF (II.EQ.N_N)
     .FUNC_B  = EXP(T_TAR1)*SD%YSMOK_E(II)/SD%TAR_CURVE(N_IERT)
      ENDDO
      UPDATE_SS  = 0.5*(TBS - TAS)*(FUNC_A+FUNC_B)
C
      ELSE

      IT_NUB     = 2**(J-2)
      DEL        = (TBS-TAS)/IT_NUB
      XDX        = TAS + 0.5*DEL
      CALL TAU_SUB(T_TAR2,DT_CURVE,XDX,N,NM,N_IERT)
      SUM        = 0.0
      IF(XDX.LT.DT_CURVE) CYCLE INTLOOP
      I_XDX      = FLOOR(XDX*10)
      ITSX       = FLOOR((XDX-DT_CURVE)*10)
      DO II      = 1,IT_NUB
      DO IXI     = 1,I_XDX
      IF(IXI.EQ.ITSX)
     .  FUNC_XDX = EXP(T_TAR2)*SD%YSMOK_E(IXI)/SD%TAR_CURVE(I_XDX)
      ENDDO
      SUM        = SUM + FUNC_XDX
      XDX        = XDX + DEL
      ENDDO
      UPDATE_SS  = 0.5*(UPDATE_SS +(TBS-TAS)*SUM/IT_NUB)
C
      ENDIF J_EQ_1
C
      OLD_SS      = UPDATE_SS
c
      ENDDO INTLOOP
      ENDIF TBS_GT_DT_CURVE
C
      SD%YSMOKE  = ARHS*UPDATE_SS
C
C Converting the smoke at inside detector into obscuration(%/m):
C
      C_ELECTRICAL=1.0
      Y_ODE     = C_ELECTRICAL*(MASS_EXTINCTION_COEFFICIENT/1000.)*
     .            (SD%YSMOKE/1000.)
      Y_ODE     = (1-EXP(-Y_ODE))*100
      SD%YSMOKE_IN = MIN(Y_ODE,SD%YSMOKE_OUT) 
C
      ENDIF PREDICTOR_IF
c
      ENDDO SMOKE_DETECTOR_LOOP
C
      TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
      END SUBROUTINE SMOKE_DETECTORS
C
C
      SUBROUTINE TAU_SUB(UPDATE_S,X1,T,N,NM,N_IER_TAU)
C
C Computation of the intergration of the mixing time TAU from 0 to T. 
C
      REAL(EB) OLD_S,FUN_A,FUN_B,FUN_XDX,UPDATE_S,XDX1,DEL1,EPS
     .        ,TA,TB,SUM1,T,X1
      INTEGER  IT_NU,N_IER_TAU,JMAX,II,JJ,ITSX,N
      PARAMETER(EPS=1.0E-6)
      INTEGER, INTENT(IN) :: NM
      TYPE (SMOKE_DETECTOR_TYPE),POINTER :: SD
C
      TNOW=SECOND()
C
      SD => SMOKE_DETECTOR(N)
C
      OLD_S     = -1.0E30
      JMAX      = 20
      UPDATE_S  = 0.0
      TA        = 0.0
      TB        = T
      N_IER_TAU= FLOOR(T*10)
C
      INTLOOP1: DO JJ  = 1,JMAX
      IF(JJ.EQ.1) THEN
      FUN_A = 0.0
      IF(SD%TAR_CURVE(N_IER_TAU).EQ.0.0) CYCLE INTLOOP1 
      FUN_B = 1/SD%TAR_CURVE(N_IER_TAU)
      UPDATE_S  = 0.5*(TB - TA)*(FUN_A+FUN_B)
      ELSE
      IT_NU     = 2**(JJ-2)
      DEL1      = (TB-TA)/IT_NU
      XDX1      = TA + 0.5*DEL1
      SUM1      = 0.0
      IF(XDX1.LT.X1) CYCLE INTLOOP1
      ITSX      = FLOOR((XDX1)*10)
      DO II     = 1,IT_NU
      FUN_XDX = 1.0/SD%TAR_CURVE(ITSX)
      SUM1      = SUM1 + FUN_XDX
      XDX1      = XDX1 + DEL1
      ENDDO
      UPDATE_S  = 0.5*(UPDATE_S +(TB-TA)*SUM1/IT_NU)
      ENDIF
C
      IF(JJ.GT.5) THEN
      IF(ABS(UPDATE_S - OLD_S).LT.EPS*ABS(OLD_S).OR.
     .   (UPDATE_S.EQ.0..AND.OLD_S.EQ.0.)) GOTO 101
      ENDIF
      OLD_S      = UPDATE_S
      ENDDO INTLOOP1
 101  CONTINUE
C
      TUSED(8,NM)=TUSED(8,NM)+SECOND()-TNOW
      END SUBROUTINE TAU_SUB
C
C
      END MODULE SPRK
