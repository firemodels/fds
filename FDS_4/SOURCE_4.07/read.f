      MODULE READ 
C
      USE PREC
      USE VARS
      USE CONS
      USE TRAN
      USE PACKER
C
      IMPLICIT NONE
C
      PRIVATE
      PUBLIC READ_DATA
      CHARACTER(60) FYI
      CHARACTER(30) QUANTITY,MAKE,LABEL,CB
      CHARACTER(30) COLOR_QUANTITY
      CHARACTER(80) MESSAGE
      CHARACTER(6)  PHASE
      CHARACTER(26) PART_ID,SURF_ID,SURF_IDS(3),SURF_ID6(6),
     .              ID,AKA,AKA_NAME(0:200),
     .              RAMP_MF(0:20),RAMP_Q,RAMP_V,RAMP_G,SURF_DEFAULT,
     .              BACKGROUND_SPECIES,REACTION,RAMP_KS,RAMP_C_P,
     .              RAMP_KS_CHAR,RAMP_C_P_CHAR,EVAC_SURF_DEFAULT
      CHARACTER(20) PROFILE,BACKING,PARTICLE_COLOR
      LOGICAL SUCCESS,BNDF_FACE(-3:3),BNDF_BLOCK,EX,
     .        THICKEN_OBSTRUCTIONS,PARTICLES
      LOGICAL PAPER_MODEL,BNDF_DEFAULT,ADIABATIC,BAD,
     .        BURN_AWAY,LEAKING,OUTLINE
      REAL(EB) KS,TAU_MF(0:20)
      REAL(EB) XYZ(3),XB(6),VEL_T(2),TEXTURE_ORIGIN(3)
      CHARACTER(60) TEXTURE_MAP
      REAL(EB) TM,TAU_Q,TAU_V,TIGN,F,HRRPUA,EMISSIVITY,TEXTURE_WIDTH,
     .         TEXTURE_HEIGHT,
     .         E_COEFFICIENT,RADIATIVE_FRACTION,VOLUME_FLUX,
     .         C_HORIZONTAL,C_VERTICAL,
     .         ALPHA,C_P,TMPWAL,TMPWAL0,DELTA,VEL,C_DELTA_RHO,VBC,
     .         TMPIGN,TMPEVAP,VBC0,DTSAM,PBX,PBY,PBZ,
     .         T,MASS_FLUX(0:20),MASS_FRACTION(0:20),XI,ETA,ZETA,
     .         Z0,PLE,DELTAH,SURFACE_DENSITY,
     .         MW_BACKGROUND,MW,DENSITY,DUMMY,
     .         THERMAL_CONDUCTIVITY,VISCOSITY,HEAT_FLUX
      REAL(EB) DIFFUSION_COEFFICIENT,MU_USER(0:20),K_USER(0:20),
     .         D_USER(0:20),POROSITY,
     .         HEAT_OF_VAPORIZATION,HEAT_OF_COMBUSTION,
     .         HEAT_OF_ABLATION,ABLATION_TEMPERATURE,ABLATION_RATE,
     .         BURNING_RATE_MAX,X_O2_LL,IN_DEPTH_COEFFICIENT,
     .         GEN_TO_YLD,HUMIDITY,MW_MIN,MW_MAX,
     .         ORIENTATION(3),ROTATION,RGB(3),RGB4(4),RADIUS,
     .         MASS_FLUX_CRITICAL,A,
     .         E,MOISTURE_FRACTION,FUEL_FRACTION,
     .         CHAR_DENSITY, C_P_CHAR, KS_CHAR,
     .         CRITICAL_FLAME_TEMPERATURE,DX_SOLID,
     .         EXTERNAL_FLUX,TMP_BACK
      INTEGER NSPRO,NTCO,II,NSFO,NNN,NR,WALL_POINTS,
     .        NPPC,IOR,NSPC,ND,NN,
     .        I1,I2,J1,J2,K1,K2,N,I,J,K,IZERO,
     .        IOS,NVO,LUDUM,ITER,NFILES,NITER
      TYPE (MESH_TYPE), POINTER :: M
      TYPE(OBSTRUCTION_TYPE), POINTER :: OB,OB2
      TYPE (VENTS_TYPE), POINTER :: VT
      TYPE(LAGRANGIAN_TYPE), POINTER :: LP
C
C
      CONTAINS
C
C
      SUBROUTINE READ_DATA
C
      CALL READ_HEAD
      CALL READ_GRID
      CALL READ_PDIM
      CALL READ_TRAN
      CALL READ_TIME
      CALL READ_MISC
      CALL READ_OBST
      CALL READ_VENT
      CALL READ_PART
      CALL READ_TREE
      CALL READ_SPEC
      CALL READ_SURF
      CALL READ_INIT
      CALL READ_RAMP
      CALL READ_SPRK
      CALL READ_HEAT
      CALL READ_SMOD
      CALL READ_THCP
      CALL READ_PL3D
      CALL READ_SLCF
      CALL READ_ISOF
      CALL READ_BNDF
C
      END SUBROUTINE READ_DATA
C
C
C
      SUBROUTINE READ_HEAD
C
      NAMELIST /HEAD/ TITLE,CHID,FYI
C
      CHID    = 'output'
      TITLE   = '      '
C
      HEAD_LOOP: DO
      CALL CHECKREAD('HEAD',LU5,IOS) ; IF (IOS.EQ.1) EXIT HEAD_LOOP
      READ(LU5,HEAD,END=13,ERR=14,IOSTAT=IOS)
   14 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with HEAD line')
      ENDDO HEAD_LOOP
   13 REWIND(LU5)
C
      CLOOP: DO I=1,39
      IF (CHID(I:I).EQ.'.')
     .   CALL SHUTDOWN('ERROR: No periods allowed in CHID')
      IF (CHID(I:I).EQ.' ') EXIT CLOOP
      ENDDO CLOOP
C
      INQUIRE(FILE=TRIM(CHID)//'.stop',EXIST=EX)
      IF (EX) THEN
         WRITE(MESSAGE,'(A,A,A)') "ERROR: Remove the file, ",
     .   TRIM(CHID)//'.stop',", from the current directory"
         CALL SHUTDOWN(MESSAGE)
         ENDIF
C
      END SUBROUTINE READ_HEAD
C
C
      SUBROUTINE READ_GRID
C
      INTEGER :: IBAR,JBAR,KBAR,NM,POISSON_BC(6)
      LOGICAL :: EVACUATION, EVAC_HUMANS
      NAMELIST /GRID/ IBAR,JBAR,KBAR,FYI,ID,SYNCHRONIZE,
     .     EVACUATION,EVAC_HUMANS,POISSON_BC
      TYPE (MESH_TYPE), POINTER :: M
C
      NMESHES = 0
C
      GRID_LOOP: DO
      CALL CHECKREAD('GRID',LU5,IOS) ; IF (IOS.EQ.1) EXIT GRID_LOOP
      READ(LU5,GRID,END=15,ERR=16,IOSTAT=IOS)
      NMESHES = NMESHES + 1
   16 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with GRID line')
      ENDDO GRID_LOOP
   15 REWIND(LU5)
C
      ALLOCATE(MESH(NMESHES),STAT=IZERO)
      CALL ChkMemErr('READ','MESH',IZERO)
      ALLOCATE(MESH_NAME(NMESHES),STAT=IZERO)
      CALL ChkMemErr('READ','MESH_NAME',IZERO)
      ALLOCATE(TUSED(N_TIMERS,NMESHES),STAT=IZERO)
      CALL ChkMemErr('READ','TUSED',IZERO)
      ALLOCATE(SYNC_TIME_STEP(NMESHES),STAT=IZERO)
      CALL ChkMemErr('READ','SYNC_TIME_STEP',IZERO)
      SYNC_TIME_STEP = .FALSE.
      ALLOCATE(EVACUATION_ONLY(NMESHES),STAT=IZERO)
      CALL ChkMemErr('READ','EVACUATION_ONLY',IZERO)
      EVACUATION_ONLY = .FALSE.
      ALLOCATE(EVACUATION_GRID(NMESHES),STAT=IZERO)
      CALL ChkMemErr('READ','EVACUATION_GRID',IZERO)
      EVACUATION_GRID = .FALSE.
      ALLOCATE(PBC(6,NMESHES),STAT=IZERO)
      CALL ChkMemErr('READ','PBC',IZERO)
C
C Read and store GRID dimensions
C
      MESH_LOOP: DO NM=1,NMESHES
      IBAR=10 ; JBAR=10 ; KBAR=10
      TWO_D = .FALSE.
      ID = 'null'
      SYNCHRONIZE = .FALSE.
      EVACUATION  = .FALSE.
      EVAC_HUMANS = .FALSE.
      POISSON_BC  = -1
      WRITE(MESH_NAME(NM),'(A,I3)') 'MESH',NM
      CALL CHECKREAD('GRID',LU5,IOS) ; IF (IOS.EQ.1) EXIT MESH_LOOP
      READ(LU5,GRID,END=115)
      M => MESH(NM)
      M%IBAR = IBAR
      M%JBAR = JBAR
      M%KBAR = KBAR
      M%NEWC = 2*IBAR*JBAR+2*IBAR*KBAR+2*JBAR*KBAR
      IF (SYNCHRONIZE) SYNC_TIME_STEP(NM)  = .TRUE.
      IF (EVACUATION)  EVACUATION_ONLY(NM) = .TRUE.
      IF (EVAC_HUMANS) EVACUATION_GRID(NM) = .TRUE.
      IF (JBAR.EQ.1) TWO_D = .TRUE.
      IF (TWO_D .AND. JBAR.NE.1) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: JBAR must be 1 for all grids'
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      IF (ID.NE.'null') MESH_NAME(NM) = ID
      PBC(:,NM) = POISSON_BC(:)
      ENDDO MESH_LOOP
  115 REWIND(LU5)
C
C Start the timing arrays
C
      TUSED      = 0.
      TUSED(1,:) = SECOND()
C
      END SUBROUTINE READ_GRID
C
C
      SUBROUTINE READ_PDIM
C
      REAL(EB) :: XBAR,YBAR,ZBAR,XBAR0,YBAR0,ZBAR0
      INTEGER :: NM
C
      REAL(EB) :: RBAR0,RBAR
      NAMELIST /PDIM/ XBAR,YBAR,ZBAR,XBAR0,YBAR0,ZBAR0,FYI,RBAR0,RBAR
C
      MESH_LOOP: DO NM=1,NMESHES
C
      RBAR0 = 0. ; RBAR = -1.
      XBAR0 = 0. ; XBAR =  1.
      YBAR0 = 0. ; YBAR =  1.
      ZBAR0 = 0. ; ZBAR =  1.
      CYLINDRICAL = .FALSE.
C
      CALL CHECKREAD('PDIM',LU5,IOS) ; IF (IOS.EQ.1) EXIT MESH_LOOP
      READ(LU5,PDIM,END=19,ERR=20,IOSTAT=IOS) 
   20 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with PDIM line')
C
      M => MESH(NM)
C
      IF (RBAR.GT.0.) THEN
         CYLINDRICAL = .TRUE.
         XBAR        = RBAR
         XBAR0       = RBAR0
         ENDIF
C
      M%XS    = XBAR0
      M%XF    = XBAR
      M%YS    = YBAR0
      M%YF    = YBAR
      M%ZS    = ZBAR0
      M%ZF    = ZBAR
      M%DXI   = (XBAR-XBAR0)/REAL(M%IBAR,EB)
      M%DETA  = (YBAR-YBAR0)/REAL(M%JBAR,EB)
      M%DZETA = (ZBAR-ZBAR0)/REAL(M%KBAR,EB)
      M%RDXI  = 1./M%DXI
      M%RDETA = 1./M%DETA
      M%RDZETA= 1./M%DZETA
      M%IBM1  = M%IBAR-1
      M%JBM1  = M%JBAR-1
      M%KBM1  = M%KBAR-1
      M%IBP1  = M%IBAR+1
      M%JBP1  = M%JBAR+1
      M%KBP1  = M%KBAR+1
C
      ENDDO MESH_LOOP
   19 REWIND(LU5)
C
      END SUBROUTINE READ_PDIM
C
C
      SUBROUTINE READ_TRAN
C
C Compute the polynomial transform function for the vertical coordinate
C
      REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: A,XX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ND
      REAL(EB)  PC,CC,COEF
      INTEGER  IEXP,IC,IDERIV,N,K,IERROR,IOS,I,MESH_NUMBER,
     .         NIPX,NIPY,NIPZ,NIPXS,NIPYS,NIPZS,NIPXF,NIPYF,NIPZF,NM
      TYPE (MESH_TYPE), POINTER :: M
      TYPE (TRAN_TYPE), POINTER :: T
      NAMELIST /TRNX/ IDERIV,CC,PC,FYI,MESH_NUMBER
      NAMELIST /TRNY/ IDERIV,CC,PC,FYI,MESH_NUMBER
      NAMELIST /TRNZ/ IDERIV,CC,PC,FYI,MESH_NUMBER
C
C Scan the input file, counting the number of NAMELIST entries
C
      ALLOCATE(TRANS(NMESHES))
C
      MESH_LOOP: DO NM=1,NMESHES
      M => MESH(NM)
      T => TRANS(NM)
C
      DO N=1,3
      T%NOC(N) = 0
      TRNLOOP: DO
      IF (N.EQ.1) THEN
         CALL CHECKREAD('TRNX',LU5,IOS) ; IF (IOS.EQ.1) EXIT TRNLOOP
         MESH_NUMBER = 1
         READ(LU5,NML=TRNX,END=17,ERR=18,IOSTAT=IOS)
         IF (MESH_NUMBER.NE.NM) CYCLE TRNLOOP
         ENDIF
      IF (N.EQ.2) THEN
         CALL CHECKREAD('TRNY',LU5,IOS) ; IF (IOS.EQ.1) EXIT TRNLOOP
         MESH_NUMBER = 1
         READ(LU5,NML=TRNY,END=17,ERR=18,IOSTAT=IOS)
         IF (MESH_NUMBER.NE.NM) CYCLE TRNLOOP
         ENDIF
      IF (N.EQ.3) THEN
         CALL CHECKREAD('TRNZ',LU5,IOS) ; IF (IOS.EQ.1) EXIT TRNLOOP
         MESH_NUMBER = 1
         READ(LU5,NML=TRNZ,END=17,ERR=18,IOSTAT=IOS)
         IF (MESH_NUMBER.NE.NM) CYCLE TRNLOOP
         ENDIF
      T%NOC(N) = T%NOC(N) + 1
   18 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with TRN* line')
      ENDDO TRNLOOP
   17 REWIND(LU5)
      ENDDO
C
      T%NOCMAX = MAX(T%NOC(1),T%NOC(2),T%NOC(3))
      ALLOCATE(A(T%NOCMAX+1,T%NOCMAX+1))
      ALLOCATE(XX(T%NOCMAX+1,3))
      ALLOCATE(ND(T%NOCMAX+1,3))
      ALLOCATE(T%C1(0:T%NOCMAX+1,3))
      T%C1               = 0.
      T%C1(1,1:3)        = 1.
      ALLOCATE(T%C2(0:T%NOCMAX+1,3))
      ALLOCATE(T%C3(0:T%NOCMAX+1,3))
      ALLOCATE(T%CCSTORE(T%NOCMAX,3))
      ALLOCATE(T%PCSTORE(T%NOCMAX,3))
      ALLOCATE(T%IDERIVSTORE(T%NOCMAX,3))
C
      T%ITRAN  = 0
C
      DO IC=1,3
      NLOOP:  DO N=1,T%NOC(IC)
      IDERIV = -1
      IF (IC.EQ.1) THEN
         LOOP1: DO
         CALL CHECKREAD('TRNX',LU5,IOS) ; IF (IOS.EQ.1) EXIT NLOOP
         MESH_NUMBER = 1
         READ(LU5,TRNX,END=1,ERR=2)
         IF (MESH_NUMBER.EQ.NM) EXIT LOOP1
         ENDDO LOOP1
         ENDIF
      IF (IC.EQ.2) THEN
         LOOP2: DO
         CALL CHECKREAD('TRNY',LU5,IOS) ; IF (IOS.EQ.1) EXIT NLOOP
         MESH_NUMBER = 1
         READ(LU5,TRNY,END=1,ERR=2)
         IF (MESH_NUMBER.EQ.NM) EXIT LOOP2
         ENDDO LOOP2
         ENDIF
      IF (IC.EQ.3) THEN
         LOOP3: DO
         CALL CHECKREAD('TRNZ',LU5,IOS) ; IF (IOS.EQ.1) EXIT NLOOP
         MESH_NUMBER = 1
         READ(LU5,TRNZ,END=1,ERR=2)
         IF (MESH_NUMBER.EQ.NM) EXIT LOOP3
         ENDDO LOOP3
         ENDIF
      T%CCSTORE(N,IC) = CC
      T%PCSTORE(N,IC) = PC
      T%IDERIVSTORE(N,IC) = IDERIV
      IF (IDERIV.GE.0) T%ITRAN(IC) = 1
      IF (IDERIV.LT.0) T%ITRAN(IC) = 2
    2 ENDDO NLOOP
    1 REWIND(LU5)
      ENDDO 
C
      ICLOOP: DO IC=1,3
C
      IF (T%ITRAN(IC).EQ.1) THEN
C
C If ITRAN=1, do a polynomial transformation 
C
      ND(1,IC)  = 0
      SELECT CASE(IC)
      CASE(1)
         XX(1,IC)    = M%XF-M%XS
         T%C1(1,IC)  = M%XF-M%XS
      CASE(2)
         XX(1,IC)    = M%YF-M%YS
         T%C1(1,IC)  = M%YF-M%YS
      CASE(3)
         XX(1,IC)    = M%ZF-M%ZS
         T%C1(1,IC)  = M%ZF-M%ZS
      END SELECT
 
      NNLOOP:  DO N=2,T%NOC(IC)+1
      IDERIV = T%IDERIVSTORE(N-1,IC)
      IF (IC.EQ.1) CC = T%CCSTORE(N-1,IC)-M%XS
      IF (IC.EQ.2) CC = T%CCSTORE(N-1,IC)-M%YS
      IF (IC.EQ.3) CC = T%CCSTORE(N-1,IC)-M%ZS
      IF (IC.EQ.1 .AND. IDERIV.EQ.0) PC = T%PCSTORE(N-1,IC)-M%XS
      IF (IC.EQ.2 .AND. IDERIV.EQ.0) PC = T%PCSTORE(N-1,IC)-M%YS
      IF (IC.EQ.3 .AND. IDERIV.EQ.0) PC = T%PCSTORE(N-1,IC)-M%ZS
      IF (IC.EQ.1 .AND. IDERIV.GT.0) PC = T%PCSTORE(N-1,IC)
      IF (IC.EQ.2 .AND. IDERIV.GT.0) PC = T%PCSTORE(N-1,IC)
      IF (IC.EQ.3 .AND. IDERIV.GT.0) PC = T%PCSTORE(N-1,IC)
        ND(N,IC) = IDERIV
        XX(N,IC) = CC
      T%C1(N,IC) = PC
      ENDDO NNLOOP
C
      DO K=1,T%NOC(IC)+1
      DO N=1,T%NOC(IC)+1
      COEF = IFAC(K,ND(N,IC))
      IEXP = K-ND(N,IC)
      IF (IEXP.LT.0) A(N,K) = 0.
      IF (IEXP.EQ.0) A(N,K) = COEF
      IF (IEXP.GT.0) A(N,K) = COEF*XX(N,IC)**IEXP
      ENDDO
      ENDDO
C
      IERROR = 0
      CALL GAUSSJ(A,T%NOC(IC)+1,T%NOCMAX+1,T%C1(1:T%NOCMAX+1,IC),
     .            1,1,IERROR)
      IF (IERROR.NE.0)
     .   CALL SHUTDOWN('ERROR: Problem with grid transformation')
C
      ENDIF
C
      IF (T%ITRAN(IC).EQ.2) THEN
C
C If ITRAN=2, do a linear transformation
C
      T%C1(0,IC) = 0.
      T%C2(0,IC) = 0.
      DO N=1,T%NOC(IC)
      IF (IC.EQ.1) CC = T%CCSTORE(N,IC)-M%XS
      IF (IC.EQ.2) CC = T%CCSTORE(N,IC)-M%YS
      IF (IC.EQ.3) CC = T%CCSTORE(N,IC)-M%ZS
      IF (IC.EQ.1) PC = T%PCSTORE(N,IC)-M%XS
      IF (IC.EQ.2) PC = T%PCSTORE(N,IC)-M%YS
      IF (IC.EQ.3) PC = T%PCSTORE(N,IC)-M%ZS
      T%C1(N,IC) = CC
      T%C2(N,IC) = PC
      ENDDO
C
      SELECT CASE(IC)
      CASE(1)
         T%C1(T%NOC(1)+1,1) = M%XF-M%XS
         T%C2(T%NOC(1)+1,1) = M%XF-M%XS
      CASE(2)
         T%C1(T%NOC(2)+1,2) = M%YF-M%YS
         T%C2(T%NOC(2)+1,2) = M%YF-M%YS
      CASE(3)
         T%C1(T%NOC(3)+1,3) = M%ZF-M%ZS
         T%C2(T%NOC(3)+1,3) = M%ZF-M%ZS
      END SELECT
C
      DO N=1,T%NOC(IC)+1
      T%C3(N,IC) = (T%C2(N,IC)-T%C2(N-1,IC))/(T%C1(N,IC)-T%C1(N-1,IC))
      ENDDO
C
      ENDIF 
C
      ENDDO ICLOOP
C
      DEALLOCATE(A)
      DEALLOCATE(XX)
      DEALLOCATE(ND)
C
C Set up grid stretching arrays
C
      ALLOCATE(M%R(0:M%IBAR),STAT=IZERO)
      CALL ChkMemErr('READ','R',IZERO)
      ALLOCATE(M%RC(0:M%IBAR+1),STAT=IZERO)
      CALL ChkMemErr('READ','RC',IZERO) ; M%RC = 1.
      ALLOCATE(M%RRN(0:M%IBP1),STAT=IZERO)
      CALL ChkMemErr('READ','RRN',IZERO) ; M%RRN = 1.
      ALLOCATE(M%X(0:M%IBAR),STAT=IZERO)
      CALL ChkMemErr('READ','X',IZERO)
      ALLOCATE(M%XC(0:M%IBP1),STAT=IZERO)
      CALL ChkMemErr('READ','XC',IZERO)
      ALLOCATE(M%HX(0:M%IBP1),STAT=IZERO)
      CALL ChkMemErr('READ','HX',IZERO)
      ALLOCATE(M%DX(0:M%IBP1),STAT=IZERO)
      CALL ChkMemErr('READ','DX',IZERO)
      ALLOCATE(M%RDX(0:M%IBP1),STAT=IZERO)
      CALL ChkMemErr('READ','RDX',IZERO)
      ALLOCATE(M%DXN(0:M%IBAR),STAT=IZERO)
      CALL ChkMemErr('READ','DXN',IZERO)
      ALLOCATE(M%RDXN(0:M%IBAR),STAT=IZERO)
      CALL ChkMemErr('READ','RDXN',IZERO)
      ALLOCATE(M%Y(0:M%JBAR),STAT=IZERO)
      CALL ChkMemErr('READ','Y',IZERO)
      ALLOCATE(M%YC(0:M%JBP1),STAT=IZERO)
      CALL ChkMemErr('READ','YC',IZERO)
      ALLOCATE(M%HY(0:M%JBP1),STAT=IZERO)
      CALL ChkMemErr('READ','HY',IZERO)
      ALLOCATE(M%DY(0:M%JBP1),STAT=IZERO)
      CALL ChkMemErr('READ','DY',IZERO)
      ALLOCATE(M%RDY(0:M%JBP1),STAT=IZERO)
      CALL ChkMemErr('READ','RDY',IZERO)
      ALLOCATE(M%DYN(0:M%JBAR),STAT=IZERO)
      CALL ChkMemErr('READ','DYN',IZERO)
      ALLOCATE(M%RDYN(0:M%JBAR),STAT=IZERO)
      CALL ChkMemErr('READ','RDYN',IZERO)
      ALLOCATE(M%Z(0:M%KBAR),STAT=IZERO)
      CALL ChkMemErr('READ','Z',IZERO)
      ALLOCATE(M%ZC(0:M%KBP1),STAT=IZERO)
      CALL ChkMemErr('READ','ZC',IZERO)
      ALLOCATE(M%HZ(0:M%KBP1),STAT=IZERO)
      CALL ChkMemErr('READ','HZ',IZERO)
      ALLOCATE(M%DZ(0:M%KBP1),STAT=IZERO)
      CALL ChkMemErr('READ','DZ',IZERO)
      ALLOCATE(M%RDZ(0:M%KBP1),STAT=IZERO)
      CALL ChkMemErr('READ','RDZ',IZERO)
      ALLOCATE(M%DZN(0:M%KBAR),STAT=IZERO)
      CALL ChkMemErr('READ','DZN',IZERO)
      ALLOCATE(M%RDZN(0:M%KBAR),STAT=IZERO)
      CALL ChkMemErr('READ','RDZN',IZERO)
C
C Define X grid stretching terms
C
      M%DXMIN = 1000.
      DO I=1,M%IBAR
      XI    = (REAL(I,EB)-.5)*M%DXI
      M%HX(I) = GP(XI,1,NM)
      M%DX(I) = M%HX(I)*M%DXI
      M%DXMIN = MIN(M%DXMIN,M%DX(I))
      IF (M%HX(I).LE.0.) THEN
         WRITE(MESSAGE,'(A,I2)') 
     .       'ERROR: x transformation not monotonic, mesh ',NM
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      M%RDX(I) = 1./M%DX(I)
      ENDDO
C
      M%HX(0)    = M%HX(1)
      M%HX(M%IBP1) = M%HX(M%IBAR)
      M%DX(0)    = M%DX(1)
      M%DX(M%IBP1) = M%DX(M%IBAR)
      M%RDX(0)    = 1./M%DX(1)
      M%RDX(M%IBP1) = 1./M%DX(M%IBAR)
C
      DO I=0,M%IBAR
      XI     = I*M%DXI
      M%X(I) = M%XS + G(XI,1,NM)
      IF (CYLINDRICAL) THEN ; M%R(I) = M%X(I)
                       ELSE ; M%R(I) = 1. ; ENDIF
      M%DXN(I)  = 0.5*(M%DX(I)+M%DX(I+1))
      M%RDXN(I) = 1./M%DXN(I)
      ENDDO
      M%X(0)      = M%XS
      M%X(M%IBAR) = M%XF
C
      DO I=1,M%IBAR
      M%XC(I) = 0.5*(M%X(I)+M%X(I-1))
      ENDDO
      M%XC(0)      = M%XS - 0.5*M%DX(0)
      M%XC(M%IBP1) = M%XF + 0.5*M%DX(M%IBP1)
C
      IF (CYLINDRICAL) THEN  
      DO I=1,M%IBAR
      M%RRN(I) = 2./(M%R(I)+M%R(I-1))
      M%RC(I)  = 0.5*(M%R(I)+M%R(I-1))
      ENDDO
      M%RRN(0)    = M%RRN(1)
      M%RRN(M%IBP1) = M%RRN(M%IBAR)
      ENDIF
C
C Define Y grid stretching terms
C
      M%DYMIN = 1000.
      DO J=1,M%JBAR
      ETA   = (REAL(J,EB)-.5)*M%DETA
      M%HY(J) = GP(ETA,2,NM)
      M%DY(J) = M%HY(J)*M%DETA
      M%DYMIN = MIN(M%DYMIN,M%DY(J))
      IF (M%HY(J).LE.0.) THEN
         WRITE(MESSAGE,'(A,I2)') 
     .       'ERROR: y transformation not monotonic, mesh ',NM
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      M%RDY(J) = 1./M%DY(J)
      ENDDO
C
      M%HY(0)    = M%HY(1)
      M%HY(M%JBP1) = M%HY(M%JBAR)
      M%DY(0)    = M%DY(1)
      M%DY(M%JBP1) = M%DY(M%JBAR)
      M%RDY(0)    = 1./M%DY(1)
      M%RDY(M%JBP1) = 1./M%DY(M%JBAR)
C
      DO J=0,M%JBAR
      ETA     = J*M%DETA
      M%Y(J)    = M%YS + G(ETA,2,NM)
      M%DYN(J)  = 0.5*(M%DY(J)+M%DY(J+1))
      M%RDYN(J) = 1./M%DYN(J)
      ENDDO
C
      M%Y(0)      = M%YS
      M%Y(M%JBAR) = M%YF
C
      DO J=1,M%JBAR
      M%YC(J) = 0.5*(M%Y(J)+M%Y(J-1))
      ENDDO
      M%YC(0)      = M%YS - 0.5*M%DY(0)
      M%YC(M%JBP1) = M%YF + 0.5*M%DY(M%JBP1)
C
C Define Z grid stretching terms
C
      M%DZMIN = 1000.
      DO K=1,M%KBAR
      ZETA  = (REAL(K,EB)-.5)*M%DZETA
      M%HZ(K) = GP(ZETA,3,NM)
      M%DZ(K) = M%HZ(K)*M%DZETA
      M%DZMIN = MIN(M%DZMIN,M%DZ(K))
      IF (M%HZ(K).LE.0.) THEN
         WRITE(MESSAGE,'(A,I2)') 
     .       'ERROR: z transformation not monotonic, mesh ',NM
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      M%RDZ(K) = 1./M%DZ(K)
      ENDDO
C
      M%HZ(0)    = M%HZ(1)
      M%HZ(M%KBP1) = M%HZ(M%KBAR)
      M%DZ(0)    = M%DZ(1)
      M%DZ(M%KBP1) = M%DZ(M%KBAR)
      M%RDZ(0)    = 1./M%DZ(1)
      M%RDZ(M%KBP1) = 1./M%DZ(M%KBAR)
C
      DO K=0,M%KBAR
      ZETA      = K*M%DZETA
      M%Z(K)    = M%ZS + G(ZETA,3,NM)
      M%DZN(K)  = 0.5*(M%DZ(K)+M%DZ(K+1))
      M%RDZN(K) = 1./M%DZN(K)
      ENDDO
C
      M%Z(0)      = M%ZS
      M%Z(M%KBAR) = M%ZF
C
      DO K=1,M%KBAR
      M%ZC(K) = 0.5*(M%Z(K)+M%Z(K-1))
      ENDDO
      M%ZC(0)      = M%ZS - 0.5*M%DZ(0)
      M%ZC(M%KBP1) = M%ZF + 0.5*M%DZ(M%KBP1)
C
C Set up arrays that will return coordinate positions
C
      NIPX   = 100*M%IBAR
      NIPY   = 100*M%JBAR
      NIPZ   = 100*M%KBAR
      NIPXS  = NINT(NIPX*M%DX(0)/(M%XF-M%XS))
      NIPXF  = NINT(NIPX*M%DX(M%IBP1)/(M%XF-M%XS))
      NIPYS  = NINT(NIPY*M%DY(0)/(M%YF-M%YS))
      NIPYF  = NINT(NIPY*M%DY(M%JBP1)/(M%YF-M%YS))
      NIPZS  = NINT(NIPZ*M%DZ(0)/(M%ZF-M%ZS))
      NIPZF  = NINT(NIPZ*M%DZ(M%KBP1)/(M%ZF-M%ZS))
      M%RDXINT = REAL(NIPX,EB)/(M%XF-M%XS)
      M%RDYINT = REAL(NIPY,EB)/(M%YF-M%YS)
      M%RDZINT = REAL(NIPZ,EB)/(M%ZF-M%ZS)
C
      ALLOCATE(M%CELLSI(-NIPXS:NIPX+NIPXF),STAT=IZERO)
      CALL ChkMemErr('READ','CELLSI',IZERO)
      ALLOCATE(M%CELLSJ(-NIPYS:NIPY+NIPYF),STAT=IZERO)
      CALL ChkMemErr('READ','CELLSJ',IZERO)
      ALLOCATE(M%CELLSK(-NIPZS:NIPZ+NIPZF),STAT=IZERO)
      CALL ChkMemErr('READ','CELLSK',IZERO)
C
      DO I=-NIPXS,NIPX+NIPXF
      M%CELLSI(I) = GINV(REAL(I,EB)/M%RDXINT,1,NM)*M%RDXI
      M%CELLSI(I) = MAX(M%CELLSI(I),-0.9_EB)
      M%CELLSI(I) = MIN(M%CELLSI(I),REAL(M%IBAR)+0.9_EB)
      ENDDO
      DO J=-NIPYS,NIPY+NIPYF
      M%CELLSJ(J) = GINV(REAL(J,EB)/M%RDYINT,2,NM)*M%RDETA
      M%CELLSJ(J) = MAX(M%CELLSJ(J),-0.9_EB)
      M%CELLSJ(J) = MIN(M%CELLSJ(J),REAL(M%JBAR)+0.9_EB)
      ENDDO
      DO K=-NIPZS,NIPZ+NIPZF
      M%CELLSK(K) = GINV(REAL(K,EB)/M%RDZINT,3,NM)*M%RDZETA
      M%CELLSK(K) = MAX(M%CELLSK(K),-0.9_EB)
      M%CELLSK(K) = MIN(M%CELLSK(K),REAL(M%KBAR)+0.9_EB)
      ENDDO
C
      ENDDO MESH_LOOP
C
C
      CONTAINS
C
      INTEGER FUNCTION IFAC(II,N)
      INTEGER II,N
      IFAC = 1
      DO I=II-N+1,II
      IFAC = IFAC*I
      ENDDO
      END FUNCTION IFAC

      END SUBROUTINE READ_TRAN
C
C
      SUBROUTINE READ_TIME
C
      REAL(EB) :: DT,VEL_CHAR
      INTEGER :: NM
      NAMELIST /TIME/ DT,TWFIN,FYI,WALL_INCREMENT,SYNCHRONIZE,
     .                EVAC_DT_FLOWFIELD,EVAC_DT_STEADY_STATE
      TYPE (MESH_TYPE), POINTER :: M
C
      DT             =-1.
      TWFIN          = 1.0
      WALL_INCREMENT = 2
      SET_UP         = .FALSE.
      SYNCHRONIZE    = .FALSE.
      EVAC_DT_FLOWFIELD = 0.01
      EVAC_DT_STEADY_STATE = 0.05
C
      TIME_LOOP: DO
      CALL CHECKREAD('TIME',LU5,IOS) ; IF (IOS.EQ.1) EXIT TIME_LOOP
      READ(LU5,TIME,END=21,ERR=22,IOSTAT=IOS)
   22 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with TIME line')
      ENDDO TIME_LOOP
   21 REWIND(LU5)
C
      IF (TWFIN.LE.0.) SET_UP = .TRUE.
C
      IF (SYNCHRONIZE) SYNC_TIME_STEP = .TRUE.
      IF (ANY(SYNC_TIME_STEP)) SYNCHRONIZE = .TRUE.
C
      MESH_LOOP: DO NM=1,NMESHES
      M=>MESH(NM)
      IF (DT.GT.0.) THEN
         M%DT = DT
      ELSE
         VEL_CHAR = 0.2*SQRT(10.*(M%ZF-M%ZS))
         M%DT     = (M%DXMIN*M%DYMIN*M%DZMIN)**(1./3.)/VEL_CHAR
      ENDIF
      IF (EVACUATION_ONLY(NM)) THEN
         SYNC_TIME_STEP(NM) = .FALSE.
         M%DT = EVAC_DT_FLOWFIELD
         ENDIF
      ENDDO MESH_LOOP
C
      END SUBROUTINE READ_TIME
C
C
      SUBROUTINE READ_MISC
C
      REAL(EB) MAX_OVER_PRESSURE,DTSAM_PART,GAUGE_TEMPERATURE
      NAMELIST /MISC/ PR,SC,TMPA,TMPO,GVEC,RF,NFRAMES,FYI,
     .                CSMAG,RAMP_G,BAROCLINIC,
     .                DT0DZ,ISOTHERMAL,INCOMPRESSIBLE,
     .                PINF,DATABASE,SURF_DEFAULT,EVAC_SURF_DEFAULT,
     .                DENSITY,DATABASE_DIRECTORY,RENDER_FILE,
     .                C_FORCED,C_VERTICAL,C_HORIZONTAL,RESTART,
     .                DTCORE,BACKGROUND_SPECIES,MW,LES,DNS,
     .                VISCOSITY,THERMAL_CONDUCTIVITY,NOISE,
     .                RADIATION,GAMMA,BNDF_DEFAULT,REACTION,
     .                MAX_OVER_PRESSURE,AUTOMATIC_Z,U0,V0,W0,HUMIDITY,
     .                POROUS_FLOOR,SUPPRESSION,CHECK_POISSON,
     .                TEXTURE_ORIGIN,NSTRATA,SMOKE3D,
     .                THICKEN_OBSTRUCTIONS,LEAK_AREA,
     .                DROP_VERTICAL_VELOCITY,
     .                DROP_HORIZONTAL_VELOCITY,RHO_SOOT,
     .                DTSAM_PART,NPPS,DTPAR,DTSPAR,DEBUG,TIMING,
     .                GAUGE_TEMPERATURE,CHARACTERISTIC_VELOCITY,
     .                FLUSH_FILE_BUFFERS,MAXIMUM_DROPLETS,
     .                EVAC_PRESSURE_ITERATIONS,EVAC_TIME_ITERATIONS,
     .                TEXTURE_DIRECTORY
C
C Physical constants
C
      R0      = 8314.3    ! Universal Gas Constant (J/kmol/K)
      R1      = 1.9862E-3 ! Universal Gas Constant (kcal/mol/K)
      TMPA    = 20.       ! Ambient temperature (C)
      TMPO    = -10000. 
      GRAV    = 9.81      ! Acceleration of gravity (m/s**2)
      GAMMA   = 1.4       ! Heat capacity ratio for air
      PINF    = 101325.   ! Ambient pressure (Pa)
      TMPM    = 273.15    ! Melting temperature of water (K)
      SIGMA   = 5.67E-8   ! Stefan-Boltzmann constant (W/m**2/K**4)
      HUMIDITY= -1.       ! Relative Humidity
      RHO_SOOT= 1850.     ! Density of soot particle (kg/m3)
C
C Empirical constants
C
      C_VERTICAL   = 1.31 ! Vertical free convection (Holman, Table 7-2)
      C_HORIZONTAL = 1.52 ! Horizontal free convection 
      C_FORCED     = 0.037 ! Forced convection coefficient
C
C Miscellanious parameters
C
      PI      = 4.*ATAN(1.0_EB)
      RPI     = 1./PI
      TWOPI   = 2.*PI
      PIO2    = PI/2.
      ONTH    = 1./3.
      THFO    = 3./4.
      ONSI    = 1./6.
      TWTH    = 2./3.
      FOTH    = 4./3.
      RFPI    = 1./(4.*PI)
C
C Numerical Parameters
C
      U0 = 0.
      V0 = 0.
      W0 = 0.
      BACKGROUND_SPECIES = 'AIR'
      VISCOSITY = -1.
      THERMAL_CONDUCTIVITY = -1.
      MU_USER = -1.
      K_USER  = -1.
      MW      = 0.   
      DENSITY = -1.
      RESTART = .FALSE.
      RADIATION      = .TRUE.
      SUPPRESSION    = .TRUE.
      CHECK_POISSON  = .FALSE.
      BAROCLINIC     = .FALSE.
      NOISE          = .TRUE.
      ISOTHERMAL     = .FALSE.  
      INCOMPRESSIBLE = .FALSE.
      BNDF_DEFAULT   = .TRUE.
      LES            = .TRUE.
      DNS            = .FALSE.
      COMBUSTION_MODEL = 1
      AUTOMATIC_Z    = .TRUE.
      MAX_OVER_PRESSURE = 1000000.
      POROUS_FLOOR = .TRUE.
      TEXTURE_ORIGIN(1) = 0.
      TEXTURE_ORIGIN(2) = 0.
      TEXTURE_ORIGIN(3) = 0.
      LEAK_AREA = 0.
      GAUGE_TEMPERATURE = -1.
      CHARACTERISTIC_VELOCITY = 1.
      MAXIMUM_DROPLETS = 500000
C
C EVACuation parameters
C
      EVAC_PRESSURE_ITERATIONS = 50
      EVAC_TIME_ITERATIONS     = 50
C
C LES parameters
C
      CSMAG                = 0.20     ! Smagorinsky constant
      PR                   = -1.0     ! Turbulent Prandtl number
      SC                   = -1.0     ! Turbulent Schmidt number
C
C Misc
C
      DATABASE             = 'null'
      DATABASE_DIRECTORY   = 'null'
      TEXTURE_DIRECTORY    = 'null'
      RAMP_G               = 'null'
      SURF_DEFAULT         = 'INERT'
      EVAC_SURF_DEFAULT    = 'INERT'
      RENDER_FILE          = 'null'
      REACTION             = 'null'
      DTCORE               = 1000000.
      GVEC(1)              = 0.        ! x-component of gravity 
      GVEC(2)              = 0.        ! y-component of gravity 
      GVEC(3)              = -GRAV     ! z-component of gravity 
      DT0DZ                = -1000000.
      NFRAMES              = 1000      ! Number of output data sets
      DTSAM_PART           = -1.
      NPPS                 = 100000    ! Number Particles Per Set
      DTPAR                = 0.05      ! Particle Insertion Interval
      DTSPAR               = 0.05      ! Droplet Insertion Interval
      RF                   = 1.00      ! Relaxation factor for no-flux
      DROPLET_FILE         = .FALSE.
      NSTRATA              = 7         ! Number bins for drop dist.
      IF (.NOT.TWO_D) SMOKE3D = .TRUE.
      IF (     TWO_D) SMOKE3D = .FALSE.
      THICKEN_OBSTRUCTIONS = .FALSE.
      DROP_VERTICAL_VELOCITY   = 0.20  
      DROP_HORIZONTAL_VELOCITY = 0.50  
      DEBUG                = .FALSE.
      TIMING               = .FALSE.
      FLUSH_FILE_BUFFERS   = .TRUE.
C
      MISC_LOOP: DO 
      CALL CHECKREAD('MISC',LU5,IOS) ; IF (IOS.EQ.1) EXIT MISC_LOOP
      READ(LU5,MISC,END=23,ERR=24,IOSTAT=IOS)
   24 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with MISC line')
      ENDDO MISC_LOOP
   23 REWIND(LU5)
C
C Look to see if DATABASE exists
C
      IF (DATABASE.EQ.'null' .AND. DATABASE_DIRECTORY.NE.'null') THEN
      DATABASE = TRIM(DATABASE_DIRECTORY)//'database4.data'
      ENDIF
C
      NFILES = 1
      IF (DATABASE.NE.'null') THEN
      INQUIRE(FILE=DATABASE,EXIST=EX)
      IF (.NOT.EX) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: DATABASE '//TRIM(DATABASE)//
     .                        ' not found'
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      NFILES = 2
      OPEN(LU80,FILE=DATABASE,FORM='FORMATTED',STATUS='OLD',ERR=190,
     .     ACTION='READ',IOSTAT=IOS)
  190 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem opening DATABASE')
      ENDIF
C
      H0    = 0.5*(U0**2+V0**2+W0**2)
      TMPA  = TMPA + TMPM
      TMPA4 = TMPA**4
C
      IF (TMPO.GT.-1000.) THEN
         TMPO = TMPO + TMPM
      ELSE 
         TMPO = TMPA
      ENDIF
C
C Humidity (40% by default, but limited for high temps)
C
      IF (HUMIDITY.LT.0.)
     .HUMIDITY=0.4*MIN(1._EB,EXP(-4890.5*(1/TMPA-1/373.15))/0.0812)
C
C Miscellaneous
C
      IF (DENSITY.GT.0.) MW = TMPA*DENSITY*R0/PINF
      MW_BACKGROUND = MW
      MU_USER(0) = VISCOSITY
      K_USER(0)  = THERMAL_CONDUCTIVITY
      HCH    = C_HORIZONTAL
      HCV    = C_VERTICAL
      SURFNAME(0) = SURF_DEFAULT
      IF (ISOTHERMAL) RADIATION = .FALSE.
      C_FORCED = C_FORCED*(1012.)*(1.8E-5)**0.2 / (0.7)**(2./3.)
      TEX_ORI = TEXTURE_ORIGIN
      IF (GAUGE_TEMPERATURE.LT.0.) THEN
         TMP_GAUGE = TMPA
      ELSE
         TMP_GAUGE = GAUGE_TEMPERATURE + TMPM
      ENDIF
C
      IF (DTSAM_PART.LT.0.) THEN
      WPAR = TWFIN/REAL(NFRAMES,EB)
      ELSE
      WPAR = DTSAM_PART
      ENDIF
C
C Background Pressure Stuff
C
      P0_MAX = PINF + MAX_OVER_PRESSURE*PINF
C
C Gravity ramp
C
      IRAMPG = 0
      NRAMP  = 0
      IF (RAMP_G.NE.'null') THEN
      NRAMP = NRAMP + 1
      RAMPID(NRAMP) = RAMP_G
      RAMPTYPE(NRAMP) = 'TIME'
      IRAMPG = NRAMP
      ENDIF
C
C Prandtl and Schmidt numbers
C
      IF (DNS) THEN
         BAROCLINIC = .TRUE.
         LES = .FALSE.
         IF (PR.LT.0.) PR = 0.7
         IF (SC.LT.0.) SC = 1.0
         ELSE
         IF (PR.LT.0.) PR = 0.5
         IF (SC.LT.0.) SC = 0.5
         ENDIF
C
      RSC = 1./SC
      RPR = 1./PR
C
C Check for a restart file
C
      APPEND = .FALSE.
      IF (RESTART) APPEND = .TRUE.
C
C Min and Max values of species and temperature
C
      TMPMIN = MIN(TMPA,TMPO)
      TMPMAX = 3000.
      YYMIN  = 0.
      YYMAX  = 1.
C
      END SUBROUTINE READ_MISC
C
C
      SUBROUTINE READ_SPEC
C
      REAL(EB) EPSIJ,SIGMA2,AMW,MASS_FRACTION_0,MW,OMEGA,TSTAR,
     .         EPSK(0:5),SIG(0:5),SIGMALJ,EPSILONKLJ,NU,MU_N,K_N,D_N0,
     .         ZZ,DZZ,TE,WGT,EPSK_N,SIG_N,MW_N
      REAL(EB) CP_N2,CP_F,CP_CO2,CP_O2,CP_H2O,XVAP,MAXIMUM_VISIBILITY
      REAL(EB) RCON_MF
      INTEGER IPC,NN
      LOGICAL READ_MIX
      INTEGER  ITMP,IYY,NSPEC_READ
      CHARACTER(26) FUEL,STATE_SPECIES(10)
      NAMELIST /SURF/ VBC,TMPIGN,TMPWAL,TMPWAL0,ALPHA,C_P,KS,DELTA,
     .                C_DELTA_RHO,MASS_FRACTION,VEL,VEL_T,NPPC,
     .                E_COEFFICIENT,HEAT_FLUX,WALL_POINTS,
     .                TIGN,TAU_Q,TAU_V,RAMP_Q,SURFACE_DENSITY,TAU_MF,
     .                RAMP_MF,PART_ID,RAMP_V,VOLUME_FLUX,RAMP_KS,
     .                RAMP_KS_CHAR,RAMP_C_P_CHAR,
     .                PROFILE,PLE,Z0,ID,AKA,MASS_FLUX,RAMP_C_P,
     .                FYI,POROSITY,PAPER_MODEL,
     .                BACKING,TMP_BACK,HRRPUA,EMISSIVITY,PHASE,RADIUS,
     .                HEAT_OF_VAPORIZATION,DENSITY,IN_DEPTH_COEFFICIENT,
     .                HEAT_OF_ABLATION,ABLATION_TEMPERATURE,
     .                ABLATION_RATE,ADIABATIC,HEAT_OF_COMBUSTION,
     .                BURNING_RATE_MAX,PARTICLES,PARTICLE_COLOR,
     .                TEXTURE_MAP,TEXTURE_WIDTH,TEXTURE_HEIGHT,RGB,RGB4,
     .                MASS_FLUX_CRITICAL,A,
     .                E,MOISTURE_FRACTION,FUEL_FRACTION,
     .                CHAR_DENSITY, C_P_CHAR, KS_CHAR,BURN_AWAY,LEAKING,
     .                DX_SOLID,EXTERNAL_FLUX,TMPEVAP
      NAMELIST /SPEC/ MASS_FRACTION_0,MW,FYI,ID,SIGMALJ,
     .                EPSILONKLJ,NU,DENSITY,THERMAL_CONDUCTIVITY,
     .                VISCOSITY,DIFFUSION_COEFFICIENT
      NAMELIST /REAC/ E,BOF,XNO,XNF,DELTAH,FYI,FUEL,EPUMO2,ID,
     .                RADIATIVE_FRACTION,HRRPUA_SHEET,TMP_LOWER,
     .                Y_O2_INFTY,Y_F_INLET,NU_FUEL,NU_O2,
     .                NU_H2O,NU_CO2,MW_FUEL,FUEL_N2,
     .                H2_YIELD,SOOT_YIELD,CO_YIELD,
     .                MASS_EXTINCTION_COEFFICIENT,DTSAM,
     .                X_O2_LL,Z_CONSTANT,VISIBILITY_FACTOR,
     .                MAXIMUM_VISIBILITY,
     .                CRITICAL_FLAME_TEMPERATURE
C
C Scan all active SURF lines to see if mixture fraction combustion
C model is desired
C
      IWATER  = 0
      IFUEL   = 0
      IOXYGEN = 0     
      ICO2    = 0
      MIXTURE_FRACTION = .FALSE.
C
      IF (FUEL_EVAPORATION) MIXTURE_FRACTION = .TRUE.
C
      N = -1
      SURF_LOOP: DO 
      N = N+1
      IF (N.GT.NBT) EXIT SURF_LOOP
C
      AKA_NAME(N) = 'null'
      IF (SURFNAME(N).EQ.'INERT') CYCLE SURF_LOOP
C
      FILE_LOOP: DO ITER=1,NFILES
C
      IF (ITER.EQ.1) LUDUM = LU5
      IF (ITER.EQ.2) LUDUM = LU80
      REWIND(LUDUM)
C
      LOOP1: DO
      CALL CHECKREAD('SURF',LUDUM,IOS) ; IF (IOS.EQ.1) CYCLE FILE_LOOP
      CALL SET_SURF_DEFAULTS
      READ(LUDUM,SURF,ERR=34,IOSTAT=IOS)
C
      NN = -1
      FIND_BACKING: DO
      NN = NN+1
      IF (NN.GT.NBT) EXIT FIND_BACKING
      IF (BACKING.EQ.SURFNAME(NN)) EXIT FIND_BACKING
      IF (NN.EQ.NBT .AND. BACKING.NE.'VOID' .AND. BACKING.NE.'EXPOSED'
     .    .AND. BACKING.NE.'INSULATED') THEN
         NBT = NBT+1
         SURFNAME(NBT) = BACKING
         EXIT FIND_BACKING
         ENDIF
      ENDDO FIND_BACKING

      IF (ID.EQ.SURFNAME(N) .OR. ID.EQ.AKA_NAME(N)) THEN
         IF (AKA.NE.'null') THEN
             AKA_NAME(N) = AKA
             CYCLE LOOP1
             ENDIF
         IF (HRRPUA.GT.0. .OR. HEAT_OF_VAPORIZATION.GT.0.) THEN
             MIXTURE_FRACTION = .TRUE.
             ENDIF 
         EXIT FILE_LOOP
         ENDIF
   34 IF (IOS.GT.0) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Problem with SURF '//
     .       TRIM(SURFNAME(N))//' or some SURF line preceding it'
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO LOOP1
C
      ENDDO FILE_LOOP
C
      IF (ID.NE.SURFNAME(N) .AND. ID.NE.AKA_NAME(N)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SURFNAME(N))//
     .                     ' not found'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
C
      REWIND(LU5)
      IF (DATABASE.NE.'null') REWIND(LU80)
C
      ENDDO SURF_LOOP
C
C Check SPEC lines for errors and determine number of extra species
C
      ALLOCATE(SPECIES_ID(0:5),STAT=IZERO)
      CALL ChkMemErr('READ','SPECIES_ID',IZERO)
      SPECIES_ID(0) = BACKGROUND_SPECIES
      SPECIES_ID(1:5) = 'nullspecies'
C
      NSPEC = 0
      READ_MIX = .FALSE.
C
      SPECLOOP: DO
      CALL CHECKREAD('SPEC',LU5,IOS) ; IF (IOS.EQ.1) EXIT SPECLOOP
      ID = 'nullspecies'
      NU = 0.
      READ(LU5,NML=SPEC,END=29,ERR=30,IOSTAT=IOS)
C
      IF (ID.EQ.'nullspecies') THEN
         WRITE(MESSAGE,'(A,I2,A)') 'ERROR: Species',N,
     .                             'needs a name (ID=...)'
         CALL SHUTDOWN(MESSAGE)
         ENDIF
C
      NSPEC = NSPEC + 1
      SPECIES_ID(NSPEC) = ID
C
      IF (NU.NE.0.) MIXTURE_FRACTION = .FALSE.
      IF (ID.EQ.'MIXTURE_FRACTION') THEN
         READ_MIX         = .TRUE.
         MIXTURE_FRACTION = .TRUE.
         IFUEL            = NSPEC
         ENDIF
      IF (ID.EQ.'WATER VAPOR') IWATER = NSPEC
      IF (ID.EQ.'OXYGEN') IOXYGEN = NSPEC
      IF (ID.EQ.'CARBON DIOXIDE') ICO2 = NSPEC
C
   30 IF (IOS.GT.0) THEN
         WRITE(MESSAGE,'(A,I2)') 
     .       'ERROR: Problem with SPECies number',NSPEC+1
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO SPECLOOP
   29 REWIND(LU5)
C
      NSPEC_READ = NSPEC
C
C Check if MIXTURE_FRACTION needs to be added to SPEC list
C
      IF (MIXTURE_FRACTION .AND. IFUEL.EQ.0) THEN
         NSPEC = NSPEC + 1
         SPECIES_ID(NSPEC) = 'MIXTURE_FRACTION'
         IFUEL = NSPEC
         ENDIF
C
C Check if water vapor should be included in the calculation
C
      IF (WATER_EVAPORATION .AND. IWATER.EQ.0) THEN
         NSPEC = NSPEC + 1
         SPECIES_ID(NSPEC) = 'WATER VAPOR'
         IWATER = NSPEC
         ENDIF
C
C Allocate species-related arrays
C
      ALLOCATE(NUN(0:NSPEC),STAT=IZERO)
      CALL ChkMemErr('READ','NUN',IZERO)
C
      ALLOCATE(YY0(0:NSPEC),STAT=IZERO)
      CALL ChkMemErr('READ','YY0',IZERO)
      ALLOCATE(RCON(0:NSPEC),STAT=IZERO)
      CALL ChkMemErr('READ','RCON',IZERO)
C
      IF (MIXTURE_FRACTION) THEN
      ALLOCATE(CP(0:100,500),STAT=IZERO)
      CALL ChkMemErr('READ','CP',IZERO)
      ALLOCATE(RCP(0:100,500),STAT=IZERO)
      CALL ChkMemErr('READ','RCP',IZERO)
      ALLOCATE(HH(0:100,0:500),STAT=IZERO)
      CALL ChkMemErr('READ','HH',IZERO)
      ALLOCATE(RCON_STATE(0:10),STAT=IZERO)
      CALL ChkMemErr('READ','RCON_STATE',IZERO)
      ALLOCATE(MWN(0:10),STAT=IZERO)
      CALL ChkMemErr('READ','MWN',IZERO)
      ELSE
      ALLOCATE(CP(0:NSPEC,500),STAT=IZERO)
      CALL ChkMemErr('READ','CP',IZERO)
      ALLOCATE(RCP(0:NSPEC,500),STAT=IZERO)
      CALL ChkMemErr('READ','RCP',IZERO)
      ALLOCATE(HH(0:NSPEC,0:500),STAT=IZERO)
      CALL ChkMemErr('READ','HH',IZERO)
      ALLOCATE(MWN(0:NSPEC),STAT=IZERO)
      CALL ChkMemErr('READ','MWN',IZERO)
      ENDIF
C
      YY0     = 0.
      MWN     = 0.
      MWN(0)  = MW_BACKGROUND
      EPSK    = 0.
      SIG     = 0.
      NUN     = 0.
C
C Initialize variables for MIX_FRAC and WATER_VAPOR, if necessary
C
      IF (MIXTURE_FRACTION .AND. .NOT.READ_MIX) THEN
         YY0(IFUEL) = 0.
         MWN(IFUEL) = 0.
         EPSK(IFUEL) = 0.
         SIG(IFUEL)  = 0.
         NUN(IFUEL)  = 0.
         ENDIF
C
      IF (IWATER.GT.0) THEN
         LP=>LAGRANGIAN(NPC)
         XVAP  = MIN(1._EB,EXP(LP%H_V_0*18./R0*(1./LP%TMP_V-1./
     .                     MIN(TMPA,LP%TMP_V))))
         YY0(IWATER) = HUMIDITY*XVAP/(29./18. + (1.-29./18.)*XVAP)
         MWN(IWATER) = 0.
         EPSK(IWATER) = 0.
         SIG(IWATER)  = 0.
         NUN(IWATER)  = 0.
         ENDIF
C
C Read SPEC info for species other than MIXTURE_FRACTION and WATER VAPOR
C
      SPEC_LOOP: DO N=1,NSPEC_READ
C
      MASS_FRACTION_0 = -1.
      MW  = 0.
      DENSITY = -1.
      EPSILONKLJ = 0.
      SIGMALJ    = 0.
      NU         = 0.
      THERMAL_CONDUCTIVITY  = -1.
      VISCOSITY     = -1.
      DIFFUSION_COEFFICIENT = -1.
C
      CALL CHECKREAD('SPEC',LU5,IOS) ; IF (IOS.EQ.1) CYCLE SPEC_LOOP
      READ(LU5,SPEC,END=31,ERR=32)
C
      IF (MASS_FRACTION_0.LT.0.) THEN
         YY0(N) = 0. 
         IF (ID.EQ.'OXYGEN') YY0(N) = 0.23
         ELSE
         YY0(N) = MASS_FRACTION_0
         ENDIF
C
      IF (DENSITY.GT.0.) MW = TMPA*DENSITY*R0/PINF
      MU_USER(N) = VISCOSITY
      K_USER(N) = THERMAL_CONDUCTIVITY
      D_USER(N) = DIFFUSION_COEFFICIENT
      MWN(N) = MW
      EPSK(N) = EPSILONKLJ
      SIG(N)  = SIGMALJ
      NUN(N)  = NU
   32 ENDDO SPEC_LOOP
   31 REWIND(LU5)
C
C Set default REACtion parameters
C
      CALL SET_REAC_DEFAULTS
C
C Search for user specified REACtion parameters
C
      REAC_LOOP: DO ITER=1,NFILES
C
      IF (ITER.EQ.2 .AND. REACTION.EQ.'null') EXIT REAC_LOOP
C
      IF (ITER.EQ.1) LUDUM = LU5
      IF (ITER.EQ.2) LUDUM = LU80
      REWIND(LUDUM)
C
      SCANLOOP: DO
      CALL CHECKREAD('REAC',LUDUM,IOS) ; IF (IOS.EQ.1) CYCLE REAC_LOOP
      CALL SET_REAC_DEFAULTS
      READ(LUDUM,REAC,ERR=434,IOSTAT=IOS)
      IF (ID.EQ.REACTION) EXIT REAC_LOOP
  434 IF (IOS.GT.0) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Problem with REAC '//
     .              TRIM(REACTION)//' or some REAC line preceding it'
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO SCANLOOP
C
      ENDDO REAC_LOOP
C
      IF (ID.NE.REACTION .AND. REACTION.NE.'null') THEN
      WRITE(MESSAGE,'(A,A,A)') 
     .    'ERROR: REAC ',TRIM(REACTION),' not found'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
C
      REWIND(LU5)
      IF (DATABASE.NE.'null') REWIND(LU80)
C
C Adjust units of input parameters
C
      IF (DELTAH.GT.0.) DELTAH_FUEL = DELTAH*1000.
      EPUMO2  = EPUMO2*1000.
      TMP_LOWER = TMP_LOWER + TMPM
      TMP_CRIT = CRITICAL_FLAME_TEMPERATURE + TMPM
      DTHRR   = DTSAM
      DTMINT  = DTSAM
      Y_O2_LL = X_O2_LL*32./(X_O2_LL*32. + (1.-X_O2_LL)*28.)
      EC_LL   = VISIBILITY_FACTOR/MAXIMUM_VISIBILITY
      E_GAS   = E*1000.
C
C Set up Mixture Fraction variables
C
      IF (MIXTURE_FRACTION) CALL STATE_RELATIONSHIPS
C
C Set chi_r
C
      CHI_R = RADIATIVE_FRACTION
C
C Get molecular weights
C
      DO N=0,NSPEC
      CALL GAS_PROPS(SPECIES_ID(N),SIG(N),EPSK(N),MWN(N))
      ENDDO
C
C Compute average molecular weight
C
      YY0(0)  = 1.
      DO N=0,NSPEC
      RCON(N) = R0/MWN(N)
      ENDDO
C
      IF (MIXTURE_FRACTION) THEN
C
      RCON_STATE(1) = R0/MW_FUEL
      RCON_STATE(2) = R0/MW_O2
      RCON_STATE(3) = R0/MW_N2
      RCON_STATE(4) = R0/MW_H2O
      RCON_STATE(5) = R0/MW_CO2
      RCON_STATE(6) = R0/MW_CO
      RCON_STATE(7) = R0/MW_H2
      RCON_STATE(8) = R0/MW_C
      RSUM0         = R0/MW_AVG(0)
C
      MW_MIN = MW_AVG(0)
      MW_MAX = MW_AVG(0)
      DO I=1,10000
      MW_MIN = MIN(MW_MIN,MW_AVG(I))
      MW_MAX = MAX(MW_MAX,MW_AVG(I))
      ENDDO
C
      IF (NSPEC.GT.1) THEN
         RCON_MF =  RSUM0
         EXTRA_SPECIES_LOOP: DO N=1,NSPEC
         IF (N.EQ.IFUEL) CYCLE EXTRA_SPECIES_LOOP
         RSUM0 = RSUM0 + YY0(N)*(RCON(N)-RCON_MF)
         MW_MIN = MIN(MW_MIN,MWN(N))
         MW_MAX = MAX(MW_MAX,MWN(N))
         ENDDO EXTRA_SPECIES_LOOP
      ENDIF
C
      RHOA = PINF/(TMPA*RSUM0)
C
      ENDIF
C
      IF (.NOT.MIXTURE_FRACTION) THEN
C
      DO N=1,NSPEC
      YY0(0) = YY0(0) - YY0(N)
      ENDDO
C
      MW_MIN = MWN(0)
      MW_MAX = MWN(0)
      RSUM0  = YY0(0)*R0/MWN(0)
      DO N=1,NSPEC
      MW_MIN = MIN(MW_MIN,MWN(N))
      MW_MAX = MAX(MW_MAX,MWN(N))
      RSUM0 = RSUM0 + YY0(N)*R0/MWN(N)
      ENDDO
C
      RHOA = PINF/(TMPA*RSUM0)
C
      ENDIF
C
C Compute constant-temperature specific heats
C
      CP_GAMMA = RCON(0)*GAMMA/(GAMMA-1.)
      IF (DT0DZ.EQ.-1000000.) DT0DZ= -GRAV/CP_GAMMA
      CPOPR = CP_GAMMA/PR
C
C Compute variable-temperature specific heats for specific species
C
      IF (.NOT.MIXTURE_FRACTION) THEN
C
      SMOKE3D = .FALSE.
C
      SPECIES_LOOP: DO N=0,NSPEC
C
      HH(N,0) = 0.
C
      TEMPERATURE_LOOP: DO J=1,500
      TE = MAX(J*0.01_EB,0.301_EB)
C
      SELECT CASE(SPECIES_ID(N))
      CASE DEFAULT   !  Use Molecular Weight
      CP(N,J) = 0.001*RCON(N)*GAMMA/(GAMMA-1.)
c     CASE('AIR')  ! Nitrogen
c     CP(N,J) = 0.931857+0.293529*TE-0.0705765*TE**2+0.00568836*TE**3+
c    .          0.001589693/TE**2
      CASE('NITROGEN')
      CP(N,J) = 0.931857+0.293529*TE-0.0705765*TE**2+0.00568836*TE**3+
     .          0.001589693/TE**2
      CASE('OXYGEN')
      CP(N,J) = 0.926844+0.191789*TE-0.0370788*TE**2+0.00299313*TE**3-
     .          0.00686447/TE**2
      CASE('METHANE')
      IF (J.LE.130) THEN
      CP(N,J) = -0.0439393+6.77983*TE-2.6576*TE**2+0.366424*TE**3+
     .          0.0424103/TE**2
      ELSE
      CP(N,J) = 5.36326+0.704042*TE-0.132134*TE**2+0.00863688*TE**3-
     .          1.65139/TE**2
      ENDIF
      CASE('CARBON DIOXIDE')
      IF (J.LE.120) THEN
      CP(N,J) = 0.568122+1.25425*TE-0.765713*TE**2+0.180645*TE**3-
     .          0.00310541/TE**2
      ELSE
      CP(N,J) = 1.32196+0.0618199*TE-0.0111884*TE**2+0.00082818*TE**3-
     .          0.146529/TE**2
      ENDIF
      CASE('WATER VAPOR')
      IF (J.LE.170) THEN
      CP(N,J) = 1.67178+0.379583*TE+0.377413*TE**2-0.140804*TE**3+
     .          0.00456328/TE**2
      ELSE
      CP(N,J) = 2.33135+0.479003*TE-0.00833212*TE**2+0.00545106*TE**3-
     .          0.619869/TE**2
      ENDIF
      CASE('HELIUM')
      CP(N,J) = 0.001*RCON(N)*GAMMA/(GAMMA-1.)
      END SELECT
      CP(N,J) = CP(N,J)*1000.            ! J/kg/K
      HH(N,J) = HH(N,J-1) + CP(N,J)*10.  ! J/kg
      RCP(N,J) = 1./CP(N,J)
      ENDDO TEMPERATURE_LOOP
C
      ENDDO SPECIES_LOOP
C
      ENDIF
C
C Compute temperature-dependent specific heat based on mixture fraction
C
      VARIABLE_CP: IF (MIXTURE_FRACTION) THEN
C
      HH(:,0) = 0.
      T_LOOP: DO J=1,500
C
      TE  = MAX(0.01_EB*J,0.301_EB)
C
      CP_O2 = 0.926844+0.191789*TE-0.0370788*TE**2+0.00299313*TE**3-
     .        0.00686447/TE**2
C
      CP_N2 = 0.931857+0.293529*TE-0.0705765*TE**2+0.00568836*TE**3+
     .        0.001589693/TE**2
C
      IF (J.LE.130) THEN   ! Ethylene
      CP_F  = -0.0439393 + 6.77983*TE - 2.6576*TE**2 + 0.366424*TE**3 +
     .         0.0424103/TE**2
      ELSE
      CP_F  = 5.36326+0.704042*TE-0.132134*TE**2+0.00863688*TE**3-
     .        1.65139/TE**2
      ENDIF
C
      IF (J.LE.120) THEN
      CP_CO2 = 0.568122+1.25425*TE-0.765713*TE**2+0.180645*TE**3-
     .         0.00310541/TE**2
      ELSE
      CP_CO2 = 1.32196+0.0618199*TE-0.0111884*TE**2+0.00082818*TE**3-
     .         0.146529/TE**2
      ENDIF
      IF (J.LE.170) THEN
      CP_H2O = 1.67178+0.379583*TE+0.377413*TE**2-0.140804*TE**3+
     .         0.00456328/TE**2
      ELSE
      CP_H2O = 2.33135+0.479003*TE-0.00833212*TE**2+0.00545106*TE**3-
     .         0.619869/TE**2
      ENDIF
C
      Z_LOOP: DO I=0,100
      IYY = 100*I
      CP(I,J) = ( Y_STATE(IYY,1)*CP_F   + 
     .            Y_STATE(IYY,2)*CP_O2  +
     .            Y_STATE(IYY,3)*CP_N2  + 
     .            Y_STATE(IYY,4)*CP_H2O +
     .            Y_STATE(IYY,5)*CP_CO2 )*1000.
      RCP(I,J) = 1./CP(I,J)
C
      HH(I,J) = HH(I,J-1) + CP(I,J)*10.
C
      ENDDO Z_LOOP
C
      ENDDO T_LOOP
C
      ENDIF VARIABLE_CP
C
C For mixture fraction model, get the fuel heat of combustion
C
      IF (MIXTURE_FRACTION) THEN
         COMBUSTION_MODEL = 2
         DELTAH_FUEL = EPUMO2*MW_O2*NU_O2/(MW_FUEL*NU_FUEL)
         DO IPC=1,NPC
         LP=>LAGRANGIAN(IPC)
         IF (LP%DELTAH.GT.0.) 
     .   LP%ADJUST_EVAPORATION = LP%DELTAH/DELTAH_FUEL
         ENDDO
         ENDIF
C
C For finite rate reaction, set parameters
C
      DO N=1,NSPEC
      IF (FUEL.EQ.SPECIES_ID(N)) IFUEL = N
      ENDDO
C
      IF (IFUEL.GT.0 .AND. IOXYGEN.GT.0) THEN
         COMBUSTION_MODEL = 3
         MW_FUEL          = MWN(IFUEL)
         IF (NUN(IFUEL)  .EQ.0.) NUN(IFUEL)   = -1.
         IF (NUN(IFUEL)  .GT.0.) NUN(IFUEL)   = -NUN(IFUEL)
         IF (NUN(IOXYGEN).GT.0.) NUN(IOXYGEN) = -NUN(IOXYGEN)
         IF (NUN(IOXYGEN).EQ.0.) THEN
            WRITE(MESSAGE,'(A)') 
     .     'ERROR: Specify a stoichiometric coefficient for oxygen'
            CALL SHUTDOWN(MESSAGE)
            ENDIF
         EPUMO2 = DELTAH_FUEL*MWN(IFUEL)*NUN(IFUEL)/
     .                       (MWN(IOXYGEN)*NUN(IOXYGEN))
      ENDIF
C
C Compute viscosities and thermal conductivities for species 0-N
C
      IF (.NOT.MIXTURE_FRACTION) THEN
C
      ALLOCATE(MU_SPEC(0:NSPEC,1:500),STAT=IZERO)
      CALL ChkMemErr('READ','MU_SPEC',IZERO)
      ALLOCATE(K_SPEC(0:NSPEC,1:500),STAT=IZERO)
      CALL ChkMemErr('READ','K_SPEC',IZERO)
C
      DO N=0,NSPEC
      SIGMA2 = SIG(N)**2
      DO ITMP=1,500
      TSTAR = ITMP*10/EPSK(N)
      OMEGA = 1.16145*TSTAR**(-0.14874) + 
     .        0.52487*EXP(-0.77320*TSTAR) +
     .        2.16178*EXP(-2.43787*TSTAR)
      MU_N = 26.69E-7*(MWN(N)*ITMP*10)**0.5/(SIGMA2*OMEGA)
      IF (MU_USER(N).GE.0.) MU_N = MU_USER(N)*(ITMP*10/TMPA)**0.75
      K_N = MU_N*CP(N,ITMP)/PR
      IF (K_USER(N).GE.0.)  K_N = K_USER(N)*(ITMP*10/TMPA)**0.75
      MU_SPEC(N,ITMP) = MU_N
      K_SPEC(N,ITMP)  = K_N
      ENDDO
      ENDDO
C
      ALLOCATE(D_SPEC(0:NSPEC,1:500),STAT=IZERO)
      CALL ChkMemErr('READ','D_SPEC',IZERO)
C
      DO N=1,NSPEC
      SIGMA2 = (0.5*(SIG(N)+SIG(0)))**2
      EPSIJ  = SQRT(EPSK(N)*EPSK(0))
      AMW    = SQRT( (MWN(N)+MWN(0))/(2.*MWN(N)*MWN(0)) )
      DO ITMP=1,500
      TSTAR = ITMP*10/EPSIJ
      OMEGA = 1.06036/TSTAR**(0.15610) +
     .        0.19300/EXP(0.47635*TSTAR) +
     .        1.03587/EXP(1.52996*TSTAR) +
     .        1.76474/EXP(3.89411*TSTAR)
      D_N0 = 0.00266E-4*AMW*(10*ITMP)**1.5/(SIGMA2*OMEGA)
      IF (D_USER(N).GE.0.) D_N0 = D_USER(N)*(ITMP*10/TMPA)**1.75
      D_SPEC(N,ITMP) = D_N0
      ENDDO
      ENDDO
C
      ENDIF
C
      IF (MIXTURE_FRACTION) THEN
C
      ALLOCATE(MU_SPEC(0:100,1:500),STAT=IZERO)
      CALL ChkMemErr('READ','MU_SPEC',IZERO) ; MU_SPEC = 0.
      ALLOCATE(K_SPEC(0:100,1:500),STAT=IZERO)
      CALL ChkMemErr('READ','K_SPEC',IZERO)  ; K_SPEC = 0.
      ALLOCATE(D_SPEC(0:100,1:500),STAT=IZERO)
      CALL ChkMemErr('READ','D_SPEC',IZERO)    ; D_SPEC = 0.
C
      STATE_SPECIES(1) = 'ETHYLENE'
      STATE_SPECIES(2) = 'OXYGEN'
      STATE_SPECIES(3) = 'NITROGEN'
      STATE_SPECIES(4) = 'WATER VAPOR'
      STATE_SPECIES(5) = 'CARBON DIOXIDE'
C
      DO N=1,5
      SIG_N = -1. ; EPSK_N = -1. ; MW_N = -.1
      CALL GAS_PROPS(STATE_SPECIES(N),SIG_N,EPSK_N,MW_N)
      SIGMA2 = SIG_N**2
      DO ITMP=1,500
      TSTAR = ITMP*10/EPSK_N
      OMEGA = 1.16145*TSTAR**(-0.14874) +
     .        0.52487*EXP(-0.77320*TSTAR) +
     .        2.16178*EXP(-2.43787*TSTAR)
      MU_N = 26.69E-7*(MW_N*ITMP*10)**0.5/(SIGMA2*OMEGA)
      DO IYY=0,100
      K_N  = MU_N*CP(IYY,ITMP)/PR
      WGT  = Y_STATE(100*IYY,N)*MW_AVG(100*IYY)/MW_N
      MU_SPEC(IYY,ITMP) = MU_SPEC(IYY,ITMP) + WGT*MU_N
      K_SPEC(IYY,ITMP)  = K_SPEC(IYY,ITMP)  + WGT*K_N
      ENDDO
      ENDDO
C
      SIGMA2 = (0.5*(SIG_N+SIG(0)))**2
      EPSIJ  = SQRT(EPSK_N*EPSK(0))
      AMW    = SQRT( (MW_N+MWN(0))/(2.*MW_N*MWN(0)) )
      DO ITMP=1,500
      TSTAR = ITMP*10/EPSIJ
      OMEGA = 1.06036/TSTAR**(0.15610) +
     .        0.19300/EXP(0.47635*TSTAR) +
     .        1.03587/EXP(1.52996*TSTAR) +
     .        1.76474/EXP(3.89411*TSTAR)
      D_N0 = 0.00266E-4*AMW*(10*ITMP)**1.5/(SIGMA2*OMEGA)
      DO IYY=0,100
      D_SPEC(IYY,ITMP) = D_SPEC(IYY,ITMP) + Y_STATE(100*IYY,N)*D_N0
      ENDDO
      ENDDO
      ENDDO
C
      ENDIF
C
C Assign Mass Per Unit Energy
C
      IF (IOXYGEN.GT.0) THEN
      ALLOCATE(MPUE(NSPEC))
      ALLOCATE(EPUM(NSPEC))
      CALL ChkMemErr('READ','MPUE',IZERO)
      MPUE(1:NSPEC) = NUN(1:NSPEC)*MWN(1:NSPEC)/
     .          ABS(EPUMO2*MWN(IOXYGEN)*NUN(IOXYGEN))
      EPUM = 1./MPUE
      ENDIF
C
C Define all possible output quantities
C
      CALL DEFINE_OUTPUT_QUANTITIES
C
      CONTAINS
C 
C
      SUBROUTINE GAS_PROPS(SPECIES,SIGMA,EPSOK,MW)
C
C Molecular weight (g/mol) and Lennard-Jones properties
C
      REAL(EB) SIGMA,EPSOK,MW,SIGMAIN,EPSOKIN,MWIN
      CHARACTER(26) SPECIES
C
      SIGMAIN = SIGMA
      EPSOKIN = EPSOK
      MWIN    = MW
C
      SELECT CASE(SPECIES)
      CASE('AIR')             ; SIGMA=3.711 ; EPSOK= 78.6  ; MW=29.
      CASE('CARBON MONOXIDE') ; SIGMA=3.690 ; EPSOK= 91.7  ; MW=28.
      CASE('CARBON DIOXIDE')  ; SIGMA=3.941 ; EPSOK=195.2  ; MW=44.
      CASE('ETHYLENE')        ; SIGMA=4.163 ; EPSOK=224.7  ; MW=28.
      CASE('HELIUM')          ; SIGMA=2.551 ; EPSOK= 10.22 ; MW= 4.
      CASE('HYDROGEN')        ; SIGMA=2.827 ; EPSOK= 59.7  ; MW= 2.
      CASE('METHANE')         ; SIGMA=3.758 ; EPSOK=148.6  ; MW=16.
      CASE('NITROGEN')        ; SIGMA=3.798 ; EPSOK= 71.4  ; MW=28.
      CASE('OXYGEN')          ; SIGMA=3.467 ; EPSOK=106.7  ; MW=32.
      CASE('PROPANE')         ; SIGMA=5.118 ; EPSOK=237.1  ; MW=44.
      CASE('WATER VAPOR')     ; SIGMA=2.641 ; EPSOK=809.1  ; MW=18.
      CASE DEFAULT            ; SIGMA=3.711 ; EPSOK= 78.6  ; MW=29.
      END SELECT  
C
      IF (SIGMAIN.GT.0.) SIGMA = SIGMAIN
      IF (EPSOKIN.GT.0.) EPSOK = EPSOKIN
      IF (MWIN   .GT.0.) MW    = MWIN 
C
      IF (SPECIES.EQ.'MIXTURE_FRACTION') MW = MW_FUEL
C
      END SUBROUTINE GAS_PROPS
C
C
      SUBROUTINE STATE_RELATIONSHIPS
C
C Calculate state relations, Y_STATE(Z), and average molecular weight.
C Y_STATE(0:10000,I): I=1-Fuel,2-O2,3-N2,4-H2O,5-CO2,6-CO,7-H2,8-SOOT
C
      REAL(EB) Z_TERM,S_TERM,TOTAL_MASS,NU_N2_AIR,NU_N2_DILUENT,ETA,CHI,
     .         Y_N2_INLET
C
C Normalize ideal stoichiometric coefficients
C
      IF (NU_FUEL.NE.1.) THEN
         NU_O2   = NU_O2/NU_FUEL
         NU_H2O  = NU_H2O/NU_FUEL
         NU_CO2  = NU_CO2/NU_FUEL
         NU_FUEL = 1.
         ENDIF
C
      Y_N2_INFTY = 1.-Y_O2_INFTY
C
C If nitrogen is in the fuel, correct MW_FUEL and Y_F_INLET
C
      Y_F_INLET  = Y_F_INLET*(1.-MW_N2*FUEL_N2/MW_FUEL)
      MW_FUEL    = MW_FUEL - MW_N2*FUEL_N2
      Y_N2_INLET = 1.-Y_F_INLET 
C
C Apply Faeth's correlation to get CO yield from soot yield
C
      IF (CO_YIELD.LT.0.) THEN
      GEN_TO_YLD = NU_CO2*MW_C/MW_FUEL
      CO_YIELD   = GEN_TO_YLD*(0.0014 + 0.37*SOOT_YIELD/GEN_TO_YLD)
      ENDIF
C
C Adjust ideal stoichiometric coefficients to account for soot, CO, etc.
C
      NU_SOOT = SOOT_YIELD*MW_FUEL/MW_C
      NU_H2   = H2_YIELD*MW_FUEL/MW_H2
      NU_CO   = CO_YIELD*MW_FUEL/MW_CO
      NU_CO2  = NU_CO2 - NU_SOOT - NU_CO
      NU_H2O  = NU_H2O - NU_H2
      NU_O2   = NU_O2 - NU_SOOT - 0.5*(NU_H2+NU_CO)
      NU_N2_AIR   = NU_O2*MW_O2*Y_N2_INFTY/(Y_O2_INFTY*MW_N2)
      NU_N2_DILUENT = NU_FUEL*MW_FUEL*Y_N2_INLET/(Y_F_INLET*MW_N2)
C
C Stoichiometric value of mixture fraction
C
      Z_F  = Y_O2_INFTY/
     .(Y_O2_INFTY+Y_F_INLET*NU_O2*MW_O2/MW_FUEL)
C
      IF (Y_F_INLET.LT.(1.-MW_N2*FUEL_N2/MW_FUEL)) 
     .         AUTOMATIC_Z = .FALSE.
      IF (DNS) AUTOMATIC_Z = .FALSE.
C
C Allocate state relation arrays
C
      ALLOCATE (MW_AVG(0:10000))
      CALL ChkMemErr('READ','MW_AVG',IZERO) ; MW_AVG=0.
      ALLOCATE (Y_STATE(0:10000,1:8))
      CALL ChkMemErr('READ','Y_STATE',IZERO) ; Y_STATE=0.
      ALLOCATE (RSUM_MF(0:10000))
      CALL ChkMemErr('READ','RSUM_MF',IZERO) ; RSUM_MF=0.
C
C Compute Y_STATE and MW_AVG
C
      Y_STATE(0,1) = 0.
      Y_STATE(0,2) = Y_O2_INFTY
      Y_STATE(0,3) = Y_N2_INFTY
      MW_AVG(0) = 1/(Y_N2_INFTY/MW_N2+Y_O2_INFTY/MW_O2)
C
      Y_STATE(10000,1) = Y_F_INLET
      Y_STATE(10000,2) = 0.
      Y_STATE(10000,3) = Y_N2_INLET
      MW_AVG(10000) = 1/(Y_F_INLET/MW_FUEL+Y_N2_INLET/MW_N2)
C
      DZZ=1./10000.
      CHI=1.0
C
      DO I=1,9999
      ZZ=I*DZZ
C
      S_TERM = NU_O2*MW_O2/MW_FUEL
      Z_TERM = ZZ*(S_TERM*Y_F_INLET+Y_O2_INFTY) - Y_O2_INFTY
      ETA = ((S_TERM-Z_TERM)*MW_FUEL-Z_TERM*NU_N2_DILUENT*MW_N2)/
     .      (NU_O2*MW_O2*(Z_TERM+1.)+NU_N2_AIR*MW_N2*Z_TERM)
      TOTAL_MASS = MW_FUEL + ETA*(NU_O2*MW_O2 + NU_N2_AIR*MW_N2)+
     .             NU_N2_DILUENT*MW_N2
C
      Y_STATE(I,1) = MW_FUEL*MAX(0._EB,1.-ETA)/TOTAL_MASS 
     .             + MW_FUEL*MIN(1._EB,ETA)*(1.-CHI)/TOTAL_MASS
      Y_STATE(I,2) = MW_O2*NU_O2*MAX(0._EB,ETA-1.)/TOTAL_MASS
     .             + MW_O2*NU_O2*MIN(1._EB,ETA)*(1.-CHI)/TOTAL_MASS
      Y_STATE(I,3) = MW_N2  *(NU_N2_AIR*ETA+NU_N2_DILUENT)/TOTAL_MASS
      Y_STATE(I,4) = MW_H2O *NU_H2O *MIN(1._EB,ETA)*CHI/TOTAL_MASS
      Y_STATE(I,5) = MW_CO2 *NU_CO2 *MIN(1._EB,ETA)*CHI/TOTAL_MASS
      Y_STATE(I,6) = MW_CO  *NU_CO  *MIN(1._EB,ETA)*CHI/TOTAL_MASS
      Y_STATE(I,7) = MW_H2  *NU_H2  *MIN(1._EB,ETA)*CHI/TOTAL_MASS
      Y_STATE(I,8) = MW_C   *NU_SOOT*MIN(1._EB,ETA)*CHI/TOTAL_MASS

      MW_AVG(I) = 1./(Y_STATE(I,1)/MW_FUEL +
     .                Y_STATE(I,2)/MW_O2   +
     .                Y_STATE(I,3)/MW_N2   +
     .                Y_STATE(I,4)/MW_H2O  +
     .                Y_STATE(I,5)/MW_CO2  +
     .                Y_STATE(I,6)/MW_CO   +
     .                Y_STATE(I,7)/MW_H2   +
     .                Y_STATE(I,8)/MW_C)
C
      ENDDO
C
C Compute R0/MW_AVG
C
      RSUM_MF = R0/MW_AVG
C
      END SUBROUTINE STATE_RELATIONSHIPS
C
C
      SUBROUTINE SET_REAC_DEFAULTS
C
      ID   = 'null'
      BOF  = 0.          ! Pre-exponential Factor (cm**3/mol-s)
      E    = 0.          ! Activation Energy (kJ/kmol)
      XNO  = 1.
      XNF  = 1.
      FUEL = 'nullspecies'
      DELTAH = -1.
      DELTAH_FUEL = 40000.  ! Energy per unit mass fuel consumed (kJ/kg)
      HRRPUA_SHEET = 200.   ! Max HRR per unit area of flame sheet
      TMP_LOWER  = 325.  ! Lower Temp for One Step Reaction
      EPUMO2     = 13100.   ! Energy per unit mass oxygen (kJ/kg)
      IF (LES) THEN
         RADIATIVE_FRACTION = 0.35
         SOOT_YIELD = 0.01
         ELSE
         RADIATIVE_FRACTION = 0.0
         SOOT_YIELD = 0.
         ENDIF
      DTSAM      = TWFIN/REAL(NFRAMES)
      Y_O2_INFTY = 0.23
      X_O2_LL    = 0.15
      CRITICAL_FLAME_TEMPERATURE = 1427. 
      Y_F_INLET  = 1.
      NU_O2      = 5.       ! Default fuel is propane
      NU_CO2     = 3.
      NU_H2O     = 4.
      NU_FUEL    = 1.
      MW_FUEL    = 44.
      MW_O2      = 32.
      MW_CO2     = 44.
      MW_C       = 12.
      MW_N2      = 28.
      MW_H2O     = 18.
      MW_CO      = 28.
      MW_H2      = 2.
      FUEL_N2    = 0.
      CO_YIELD   = -1.
      H2_YIELD   =  0.
      MASS_EXTINCTION_COEFFICIENT = 7600.     ! m2/kg
      VISIBILITY_FACTOR = 3.
      MAXIMUM_VISIBILITY = 30.
      Z_CONSTANT = 1.0
C
      END SUBROUTINE SET_REAC_DEFAULTS
C
      END SUBROUTINE READ_SPEC
C
C
      SUBROUTINE READ_PART
C
      INTEGER :: NUMBER_INITIAL_DROPLETS,NIP,NISP,NSPINS,NPDIM,
     .           NSPDIM,SAMPLING_FACTOR,DROPLETS_PER_SECOND,
     .           NPSAM,NSPSAM
      REAL(EB) :: SPECIFIC_HEAT,VAPORIZATION_TEMPERATURE,
     .            MELTING_TEMPERATURE,MASS_PER_VOLUME,DIAMETER,
     .            GAMMA_D,WMPUV,DTSAM,AGE,DELAY,INITIAL_TEMPERATURE
      CHARACTER(30) :: SPEC_ID,HEAT_ACTIVATE
      LOGICAL :: MASSLESS,STATIC,FUEL,WATER,TREE,FUEL_DROPLETS
      NAMELIST /PART/ NUMBER_INITIAL_DROPLETS,QUANTITY,FYI,
     .                DROPLETS_PER_SECOND,MASS_PER_VOLUME,
     .                SAMPLING_FACTOR,NPSAM,NSPSAM,
     .                ID,STATIC,MASSLESS,FUEL,WATER,TREE,
     .                DENSITY,VAPORIZATION_TEMPERATURE,
     .                SPECIFIC_HEAT,HEAT_OF_VAPORIZATION,
     .                MELTING_TEMPERATURE,DIAMETER,GAMMA_D,WMPUV,
     .                HEAT_OF_COMBUSTION,DTPAR,DTSPAR,
     .                NIP,NISP,NSPINS,DTSAM,AGE,NPDIM,NSPDIM,
     .                FUEL_DROPLETS,SPEC_ID,INITIAL_TEMPERATURE,XB
      EQUIVALENCE (SAMPLING_FACTOR,NPSAM,NSPSAM)
      EQUIVALENCE (NUMBER_INITIAL_DROPLETS,NISP,NIP)
      EQUIVALENCE (WMPUV,MASS_PER_VOLUME)
      EQUIVALENCE (FUEL,FUEL_DROPLETS)
C
      REAL(EB) T_ACTIVATE,T_DEACTIVATE
      NAMELIST /SPRK/ XYZ,MAKE,T_ACTIVATE,T_DEACTIVATE,FYI,
     .                ORIENTATION,ROTATION,PART_ID,LABEL,DELAY,
     .                HEAT_ACTIVATE
C
      WATER_EVAPORATION = .FALSE.
      FUEL_EVAPORATION  = .FALSE.
      EVAPORATION       = .FALSE.
C
C Look for sprinklers, if found WATER VAPOR will be a species
C
      NSPR = 0
      SPRLOOP: DO
      PART_ID = 'null'
      CALL CHECKREAD('SPRK',LU5,IOS) ; IF (IOS.EQ.1) EXIT SPRLOOP
      READ(LU5,SPRK,END=7,ERR=8,IOSTAT=IOS)
      NSPR = NSPR + 1
      IF (PART_ID.EQ.'null') WATER_EVAPORATION = .TRUE.
    8 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with SPRK line')
      ENDDO SPRLOOP
    7 REWIND(LU5)
C
C Determine total number of PART lines in the input file
C
      NPC = 0
      COUNT_PART_LOOP: DO
      CALL CHECKREAD('PART',LU5,IOS) 
      IF (IOS.EQ.1) EXIT COUNT_PART_LOOP
      READ(LU5,NML=PART,END=219,ERR=220,IOSTAT=IOS)
      NPC = NPC + 1
  220 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with PART line')
      ENDDO COUNT_PART_LOOP
  219 REWIND(LU5)
C
C Add two special PART classes, sprinkler droplets and tracers
C
      NPC = NPC + 2
C
C Allocate quantities for PART classes
C
      IPART = 0
C
      ALLOCATE(LAGRANGIAN(NPC),STAT=IZERO)
      CALL ChkMemErr('READ','LAGRANGIAN',IZERO) 

      LAGRANGIAN(1:NPC)%ADJUST_EVAPORATION = 1. 
      LAGRANGIAN(1:NPC)%SPECIES_INDEX = 0
      LAGRANGIAN(1:NPC)%SPECIES = ' '
      LAGRANGIAN(1:NPC)%COLOR_INDEX = 0
C
      READ_PART_LOOP: DO N=1,NPC
C
      LP=>LAGRANGIAN(N)
C
      DENSITY                  = 1000.   ! kg/m3
      MASS_PER_VOLUME          = 1.      ! kg/m3
      VAPORIZATION_TEMPERATURE = 100.0   ! C
      INITIAL_TEMPERATURE      = TMPA - TMPM ! C
      MELTING_TEMPERATURE      = TMPM - TMPM ! C
      SPECIFIC_HEAT            = 4.184   ! kJ/kg-K
      HEAT_OF_VAPORIZATION     = 2259.   ! kJ/kg
      HEAT_OF_COMBUSTION       = -1.     ! kJ/kg
      DROPLETS_PER_SECOND      = 1000
      DIAMETER                 = 100.    ! mu-m
      GAMMA_D                  = 2.4
      AGE                      = 1.E6    ! s
C
      ID                       = 'null'
      QUANTITY                 = 'null'
      SPEC_ID                  = 'null'
      SAMPLING_FACTOR          = -1      
      NUMBER_INITIAL_DROPLETS  = 0
      XB                       = 0.
      FUEL                     = .FALSE.
      WATER                    = .FALSE.
      STATIC                   = .FALSE.
      MASSLESS                 = .FALSE.
      TREE                     = .FALSE.
C
C Read the PART line
C
      IF (N.LT.NPC-1) THEN
      CALL CHECKREAD('PART',LU5,IOS) ; IF (IOS.EQ.1) EXIT READ_PART_LOOP
      READ(LU5,PART,END=25,IOSTAT=IOS)
C
      IF (NUMBER_INITIAL_DROPLETS.GT.0 .AND. RESTART)
     .    NUMBER_INITIAL_DROPLETS =  0
      ENDIF
C
      IF (N.EQ.NPC-1) THEN
      MASSLESS = .TRUE.
      ID       = 'tracer particles'
      ENDIF
C
      IF (SAMPLING_FACTOR.LE.0) THEN
      IF (MASSLESS) THEN
         SAMPLING_FACTOR = 1
         ELSE
         SAMPLING_FACTOR = 10
         ENDIF
         ENDIF
C
      IF (N.EQ.NPC) THEN
      WATER = .TRUE.
      ID    = 'water droplets'
      SPEC_ID = 'WATER VAPOR'
      AGE   = 60.
      ENDIF
C
      IF (TREE) FUEL   = .TRUE.
      IF (TREE) STATIC = .TRUE.
C
      IF (MASSLESS) DIAMETER = 0.
C
      LP%CLASS_NAME = ID
C
      IF (QUANTITY.EQ.'TEMPERATURE')       IPART = -1
      IF (QUANTITY.EQ.'HRRPUV')            IPART = -2
      IF (QUANTITY.EQ.'DIAMETER')          IPART = -3
      IF (QUANTITY.EQ.'VELOCITY')          IPART = -4
      IF (QUANTITY.EQ.'DROPLET_PHASE')     IPART = -5
      IF (QUANTITY.EQ.'AGE')               IPART = -6
      IF (QUANTITY.EQ.'BLACK')    LP%COLOR_INDEX = 0
      IF (QUANTITY.EQ.'YELLOW')   LP%COLOR_INDEX = 1
      IF (QUANTITY.EQ.'BLUE')     LP%COLOR_INDEX = 2
      IF (QUANTITY.EQ.'RED')      LP%COLOR_INDEX = 3
      IF (QUANTITY.EQ.'GREEN')    LP%COLOR_INDEX = 4
      IF (QUANTITY.EQ.'MAGENTA')  LP%COLOR_INDEX = 5
      IF (QUANTITY.EQ.'CYAN')     LP%COLOR_INDEX = 6
      IF (QUANTITY.EQ.'WHITE')    LP%COLOR_INDEX = 7
C
      IF (QUANTITY.EQ.'null' .AND. WATER) LP%COLOR_INDEX = 2
      IF (QUANTITY.EQ.'null' .AND. TREE)  LP%COLOR_INDEX = 4
C
      LP%SAMPLING = SAMPLING_FACTOR
C
      IF (NUMBER_INITIAL_DROPLETS.GT.0) DROPLET_FILE  = .TRUE.
C
      IF (NUMBER_INITIAL_DROPLETS.GT.0) THEN
         DO I=1,5,2
         IF (XB(I).GT.XB(I+1)) THEN
            DUMMY   = XB(I)
            XB(I)   = XB(I+1)
            XB(I+1) = DUMMY
            ENDIF
         ENDDO
         ENDIF
C
      LP%X1 = XB(1)
      LP%X2 = XB(2)
      LP%Y1 = XB(3)
      LP%Y2 = XB(4)
      LP%Z1 = XB(5)
      LP%Z2 = XB(6)
C
      LP%N_INITIAL = NUMBER_INITIAL_DROPLETS
      LP%N_INSERT  = DROPLETS_PER_SECOND*DTSPAR
C
      IF (DIAMETER.GT.0.) THEN
      ALLOCATE(LP%CDF(0:NDC),STAT=IZERO)
      CALL ChkMemErr('READ','CDF',IZERO)
      ALLOCATE(LP%R_CDF(0:NDC),STAT=IZERO)
      CALL ChkMemErr('READ','R_CDF',IZERO)
      ALLOCATE(LP%IL_CDF(NSTRATA),STAT=IZERO)
      CALL ChkMemErr('READ','IL_CDF',IZERO)
      ALLOCATE(LP%IU_CDF(NSTRATA),STAT=IZERO)
      CALL ChkMemErr('READ','IU_CDF',IZERO)
      ALLOCATE(LP%W_CDF(NSTRATA),STAT=IZERO)
      CALL ChkMemErr('READ','W_CDF',IZERO)
      ENDIF
C
C     Assign property data to LAGRANGIAN class
C
      LP%DIAMETER = DIAMETER*1.E-6
      LP%GAMMA    = GAMMA_D
      LP%SIGMA    = 1.15/GAMMA_D
      LP%TMP_INITIAL = INITIAL_TEMPERATURE + TMPM
      LP%C_P      = SPECIFIC_HEAT*1000.
      LP%H_V_0    = HEAT_OF_VAPORIZATION*1000.
      LP%DELTAH   = HEAT_OF_COMBUSTION*1000.
      LP%TMP_V    = VAPORIZATION_TEMPERATURE + TMPM
      LP%DENSITY  = DENSITY
      LP%FTPR     = FOTH*PI*DENSITY
      LP%MASS_PER_VOLUME = MASS_PER_VOLUME
      LP%TMP_MELT = MELTING_TEMPERATURE + TMPM
      LP%MASSLESS = MASSLESS
      LP%LIFETIME = AGE
      LP%TREE     = TREE
      LP%FUEL     = FUEL
      LP%WATER    = WATER
      LP%STATIC   = STATIC
      IF (FUEL)  SPEC_ID = 'MIXTURE_FRACTION'
      IF (WATER) SPEC_ID = 'WATER VAPOR'
      LP%SPECIES  = SPEC_ID
C
      IF (N.LT.NPC-1 .AND. FUEL)  FUEL_EVAPORATION = .TRUE.
      IF (N.LT.NPC-1 .AND. WATER) WATER_EVAPORATION = .TRUE.
C
      ENDDO READ_PART_LOOP
   25 REWIND(LU5)
C
      IF (FUEL_EVAPORATION .OR. WATER_EVAPORATION) EVAPORATION=.TRUE.
C
      END SUBROUTINE READ_PART
C
C
      SUBROUTINE READ_TREE
C
      INTEGER :: IPC,N_TREES_0,NM,NN,N
      REAL(EB) :: CANOPY_WIDTH,CANOPY_BASE_HEIGHT,TREE_HEIGHT
      NAMELIST /TREE/ XYZ,CANOPY_WIDTH,CANOPY_BASE_HEIGHT,
     .                     TREE_HEIGHT,PART_ID
C
C Read the TREE lines to determine how many cone shaped trees
C there will be
C
      N_TREES = 0
      COUNT_TREE_LOOP: DO
      CALL CHECKREAD('TREE',LU5,IOS) 
      IF (IOS.EQ.1) EXIT COUNT_TREE_LOOP
      READ(LU5,NML=TREE,END=11,ERR=12,IOSTAT=IOS)
      N_TREES = N_TREES + 1
   12 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with TREE line')
      ENDDO COUNT_TREE_LOOP
   11 REWIND(LU5)
C
C Sequentially read the CONE_TREE namelist to get shape and size
C parameters for each tree.
C
      IF (N_TREES.EQ.0) RETURN
C
      ALLOCATE(CANOPY_W(N_TREES),STAT=IZERO)
      CALL ChkMemErr('READ','CANOPY_W',IZERO)
      ALLOCATE(CANOPY_B_H(N_TREES),STAT=IZERO)
      CALL ChkMemErr('READ','CANOPY_B_H',IZERO)
      ALLOCATE(TREE_H(N_TREES),STAT=IZERO)
      CALL ChkMemErr('READ','TREE_H',IZERO)
      ALLOCATE(X_TREE(N_TREES),STAT=IZERO)
      CALL ChkMemErr('READ','X_TREE',IZERO)
      ALLOCATE(Y_TREE(N_TREES),STAT=IZERO)
      CALL ChkMemErr('READ','Y_TREE',IZERO)
      ALLOCATE(Z_TREE(N_TREES),STAT=IZERO)
      CALL ChkMemErr('READ','Z_TREE',IZERO)
      ALLOCATE(TREE_PARTICLE_CLASS(N_TREES),STAT=IZERO)
      CALL ChkMemErr('READ','TREE_PARTICLE_CLASS',IZERO)
      ALLOCATE(TREE_MESH(N_TREES),STAT=IZERO)
      CALL ChkMemErr('READ','TREE_MESH',IZERO)
C
      N_TREES_0 = N_TREES
      N = 0
C
      CONE_LOOP: DO NN=1,N_TREES_0
      N = N + 1
C
      PART_ID = 'null'
C
      CALL CHECKREAD('TREE',LU5,IOS) ; IF (IOS.EQ.1) EXIT CONE_LOOP
      READ(LU5,TREE,END=25,IOSTAT=IOS)
C
      MESH_LOOP: DO NM=1,NMESHES
      IF (XYZ(1).GE.MESH(NM)%XS .AND. XYZ(1).LE.MESH(NM)%XF .AND.
     .    XYZ(2).GE.MESH(NM)%YS .AND. XYZ(2).LE.MESH(NM)%YF .AND.
     .    XYZ(3).GE.MESH(NM)%ZS .AND. XYZ(3).LE.MESH(NM)%ZF) THEN
         TREE_MESH(N) = NM
         EXIT MESH_LOOP
         ENDIF
      IF (NM.EQ.NMESHES) THEN
         N    = N-1
         N_TREES = N_TREES - 1
         CYCLE CONE_LOOP
         ENDIF
      ENDDO MESH_LOOP
C
      IF (PART_ID.EQ.'null')
     .   CALL SHUTDOWN('ERROR: Specify PART_ID of tree')
C
      CANOPY_W(N)   = CANOPY_WIDTH
      CANOPY_B_H(N) = CANOPY_BASE_HEIGHT
      TREE_H(N)     = TREE_HEIGHT
      X_TREE(N) = XYZ(1)
      Y_TREE(N) = XYZ(2)
      Z_TREE(N) = XYZ(3)
C
      DO IPC=1,NPC
      LP=>LAGRANGIAN(IPC)
      IF (LP%CLASS_NAME.EQ.PART_ID)
     .   TREE_PARTICLE_CLASS(N) = IPC
      ENDDO
      DROPLET_FILE=.TRUE.
C
      ENDDO CONE_LOOP
   25 REWIND(LU5)
C
      END SUBROUTINE READ_TREE
C
C
      SUBROUTINE READ_SURF
C
      REAL(EB) :: MINIMUM_DENSITY,MAXIMUM_DENSITY,
     .            MINIMUM_MASS_FRACTION,MAXIMUM_MASS_FRACTION,
     .            MINIMUM_TEMPERATURE,MAXIMUM_TEMPERATURE
      INTEGER :: IBC,IPC,NN
      LOGICAL :: BURNING
      REAL(EB) :: WIDTH,SFAC,LS,DLDS,CONST,X1,X2,X0
      NAMELIST /SURF/ VBC,TMPIGN,TMPWAL,TMPWAL0,ALPHA,C_P,KS,DELTA,
     .                C_DELTA_RHO,MASS_FRACTION,VEL,VEL_T,NPPC,
     .                E_COEFFICIENT,HEAT_FLUX,WALL_POINTS,
     .                TIGN,TAU_Q,TAU_V,RAMP_Q,SURFACE_DENSITY,TAU_MF,
     .                RAMP_MF,PART_ID,RAMP_V,VOLUME_FLUX,RAMP_KS,
     .                RAMP_KS_CHAR,RAMP_C_P_CHAR,
     .                PROFILE,PLE,Z0,ID,AKA,MASS_FLUX,RAMP_C_P,
     .                FYI,POROSITY,PAPER_MODEL,
     .                BACKING,TMP_BACK,HRRPUA,EMISSIVITY,PHASE,RADIUS,
     .                HEAT_OF_VAPORIZATION,DENSITY,IN_DEPTH_COEFFICIENT,
     .                HEAT_OF_ABLATION,ABLATION_TEMPERATURE,
     .                ABLATION_RATE,ADIABATIC,HEAT_OF_COMBUSTION,
     .                BURNING_RATE_MAX,PARTICLES,PARTICLE_COLOR,
     .                TEXTURE_MAP,TEXTURE_WIDTH,TEXTURE_HEIGHT,RGB,RGB4,
     .                MASS_FLUX_CRITICAL,A, 
     .                E,MOISTURE_FRACTION,FUEL_FRACTION,
     .                CHAR_DENSITY, C_P_CHAR, KS_CHAR,BURN_AWAY,LEAKING,
     .                DX_SOLID,EXTERNAL_FLUX,TMPEVAP
      NAMELIST /CLIP/ MINIMUM_DENSITY,MAXIMUM_DENSITY,FYI,
     .                MINIMUM_MASS_FRACTION,MAXIMUM_MASS_FRACTION,
     .                MINIMUM_TEMPERATURE,MAXIMUM_TEMPERATURE
C
C No clipping for mixture fraction because of HRR computation
C
      IF (MIXTURE_FRACTION) YYMIN(IFUEL) = -0.5
C
C Miscellaneous control parameters
C
      PAPERMODEL   = .FALSE.
      NWP_MAX      = 0
      IBCDEF       = 0
C
C Get default surface properties
C
      CALL SET_SURF_DEFAULTS
C
      ALLOCATE(BCV(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','BCV',IZERO) ; BCV = VBC
      ALLOCATE(TMP_I(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','TMP_I',IZERO) ; TMP_I = TMPIGN+TMPM
      ALLOCATE(TMP_E(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','TMP_E',IZERO) ; TMP_E = TMPEVAP+TMPM
      ALLOCATE(E_SOLID(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','E_SOLID',IZERO) 
      ALLOCATE(A_SOLID(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','A_SOLID',IZERO) 
      ALLOCATE(TMP_P(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','TMP_P',IZERO) ; TMP_P = TMPWAL+TMPM
      ALLOCATE(TMP_P0(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','TMP_P0',IZERO) ; TMP_P0 = TMPWAL0+TMPM
      ALLOCATE(TIGNS(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','TIGNS',IZERO) ; TIGNS = TIGN
      ALLOCATE(C_PV(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','C_PV',IZERO) ; C_PV = C_P
      ALLOCATE(C_PC(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','C_PC',IZERO) ; C_PC = C_P_CHAR
      ALLOCATE(TMP_VOID(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','TMP_VOID',IZERO) ; TMP_VOID = TMPA
      ALLOCATE(XKS(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','XKS',IZERO) ; XKS=KS
      ALLOCATE(XKS_C(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','XKS_C',IZERO) ; XKS_C=KS_CHAR
      ALLOCATE(VELS(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','VELS',IZERO) ; VELS = VEL
      ALLOCATE(VEL_TS(2,0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','VEL_TS',IZERO) ; VEL_TS(:,0:NBT+3) = -999.
      ALLOCATE(CDR(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','CDR',IZERO) ; CDR = C_DELTA_RHO*1000.
      ALLOCATE(PCI(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','PCI',IZERO) ; PCI = 0
      ALLOCATE(POROS(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','POROS',IZERO) ; POROS = 0.
      ALLOCATE(MASSFRACS(0:NBT+3,0:NSPEC),STAT=IZERO)
      CALL ChkMemErr('READ','MASSFRACS',IZERO) ; MASSFRACS = -1.
      ALLOCATE(MASSFLUXS(0:NBT+3,0:NSPEC),STAT=IZERO)
      CALL ChkMemErr('READ','MASSFLUXS',IZERO) ; MASSFLUXS = 0.
      ALLOCATE(MASSFLUXS_MAX(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','MASSFLUXS_MAX',IZERO) ; MASSFLUXS_MAX = 0.
      ALLOCATE(H_V(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','H_V',IZERO) ; H_V=1000*HEAT_OF_VAPORIZATION
      ALLOCATE(ADJUST_BURN_RATE(0:NBT+3,0:NSPEC),STAT=IZERO)
      CALL ChkMemErr('READ','ADJUST_BURN_RATE',IZERO) 
      ADJUST_BURN_RATE=1.
      ALLOCATE(HEATFLUX(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','HEATFLUX',IZERO) ; HEATFLUX = HEAT_FLUX
      ALLOCATE(TAUMF(0:NBT+3,0:NSPEC),STAT=IZERO)
      CALL ChkMemErr('READ','TAUMF',IZERO) ; TAUMF = 1.
      ALLOCATE(DX_W(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','DX_W',IZERO) 
      ALLOCATE(DXI_W(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','DXI_W',IZERO) ; DXI_W = -1.
      ALLOCATE(MHC(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','MHC',IZERO) ; MHC = 1
      ALLOCATE(NWP(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','NWP',IZERO) ; NWP = WALL_POINTS
      ALLOCATE(IBACK(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','IBACK',IZERO) ; IBACK = 1
      ALLOCATE(MMT(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','MMT',IZERO) ; MMT = 1
      ALLOCATE(ECOEF(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','ECOEF',IZERO) ; ECOEF = 0.
      ALLOCATE(TAUQ(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','TAUQ',IZERO) ; TAUQ = TAU_Q
      ALLOCATE(TAUV(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','TAUV',IZERO) ; TAUV = TAU_V
      ALLOCATE(LEAKY(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','LEAKY',IZERO) ; LEAKY = .FALSE.
      ALLOCATE(DENSITY_S(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','DENSITY_S',IZERO)
      ALLOCATE(DENSITY_F(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','DENSITY_F',IZERO) ; DENSITY_F = 0.
      ALLOCATE(BURNAWAY(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','BURNAWAY',IZERO) ; BURNAWAY=.FALSE.
      ALLOCATE(WALL_THICKNESS(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','WALL_THICKNESS',IZERO) 
      ALLOCATE(SURFACE_DENSITY_S(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','SURFACE_DENSITY_S',IZERO) 
      ALLOCATE(E_WALL_S(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','E_WALL_S',IZERO) ; E_WALL_S = 1.
      ALLOCATE(SURF_TYPE(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','SURF_TYPE',IZERO) ; SURF_TYPE = 0
      ALLOCATE(SURF_RGB(0:NBT+3,4),STAT=IZERO)
      CALL ChkMemErr('READ','SURF_RGB',IZERO) 
      SURF_RGB(:,1)=1.0 ; SURF_RGB(:,2)=0.8 ; SURF_RGB(:,3)=0.4 
      SURF_RGB(:,4)=1.0
      ALLOCATE(TEX_MAP(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','TEX_MAP',IZERO) ; TEX_MAP = 'null'
      ALLOCATE(TEX_WIDTH(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','TEX_WIDTH',IZERO) ; TEX_WIDTH = 1.
      ALLOCATE(TEX_HEIGHT(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','TEX_HEIGHT',IZERO) ; TEX_HEIGHT = 1.
      ALLOCATE(PYROLYSIS_MODEL(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','PYROLYSIS_MODEL',IZERO)  
      PYROLYSIS_MODEL = 'null'
      ALLOCATE(SURF_GEOM(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','SURF_GEOM',IZERO) ; SURF_GEOM = 'CARTESIAN'
      ALLOCATE(SURF_PARTICLE_CLASS(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','SURF_PARTICLE_CLASS',IZERO)
      SURF_PARTICLE_CLASS = -1
      ALLOCATE(NPPCS(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','NPPCS',IZERO) ; NPPCS = NPPC
      ALLOCATE(IPROF(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','IPROF',IZERO) ; IPROF = 0
      ALLOCATE(PLES(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','PLES',IZERO) ; PLES = PLE
      ALLOCATE(Z0S(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','Z0S',IZERO) ; Z0S = Z0 
      ALLOCATE(VFLUX(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','VFLUX',IZERO) ; VFLUX = VOLUME_FLUX
      ALLOCATE(RHO_CHAR(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','RHO_CHAR',IZERO) ; RHO_CHAR = CHAR_DENSITY
      ALLOCATE(RHO_MOIS(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','RHO_MOIS',IZERO)
      RHO_MOIS = MOISTURE_FRACTION
      ALLOCATE(QEXT(0:NBT+3),STAT=IZERO)
      CALL ChkMemErr('READ','QEXT',IZERO) ; QEXT = 0.
C
C     Read the SURF lines
C
      BCLOOP: DO N=0,NBT
C
      AKA_NAME(N) = 'null'
      IF (SURFNAME(N).EQ.SURF_DEFAULT) IBCDEF = N
      IF (SURFNAME(N).EQ.'INERT') CYCLE BCLOOP
C
      SURF_LOOP: DO ITER=1,NFILES
C
      IF (ITER.EQ.1) LUDUM = LU5
      IF (ITER.EQ.2) LUDUM = LU80
      REWIND(LUDUM)
C
      LOOP1: DO
      CALL CHECKREAD('SURF',LUDUM,IOS) ; IF (IOS.EQ.1) CYCLE SURF_LOOP
      CALL SET_SURF_DEFAULTS
      READ(LUDUM,SURF,ERR=34,IOSTAT=IOS) 
      IF (ID.EQ.SURFNAME(N) .OR. ID.EQ.AKA_NAME(N)) THEN
         IF (AKA.NE.'null') THEN
             AKA_NAME(N) = AKA
             CYCLE LOOP1
             ENDIF
         EXIT SURF_LOOP
         ENDIF
   34 IF (IOS.GT.0) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Problem with SURF '//
     .     TRIM(SURFNAME(N))//' or some SURF line preceding it'
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO LOOP1
C
      ENDDO SURF_LOOP
C
      IF (ID.NE.SURFNAME(N) .AND. ID.NE.AKA_NAME(N)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SURFNAME(N))//
     .                     ' not found'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
C
      REWIND(LU5)
      IF (DATABASE.NE.'null') REWIND(LU80)
C
C     Count up the RAMPs
C
      QUANLOOP: DO I=1,6+NSPEC
      ID = 'null'
      IF (I.EQ.1 .AND. RAMP_Q.NE.'null') ID = RAMP_Q
      IF (I.EQ.2 .AND. RAMP_V.NE.'null') ID = RAMP_V
      IF (I.EQ.3 .AND. RAMP_KS.NE.'null') ID = RAMP_KS
      IF (I.EQ.4 .AND. RAMP_C_P.NE.'null') ID = RAMP_C_P
      IF (I.EQ.5 .AND. RAMP_KS_CHAR.NE.'null') ID = RAMP_KS_CHAR
      IF (I.EQ.6 .AND. RAMP_C_P_CHAR.NE.'null') ID = RAMP_C_P_CHAR
      IF (I.GT.6) THEN
         IF (RAMP_MF(I-6).NE.'null') ID = RAMP_MF(I-6)
         ENDIF
      IF (ID.EQ.'null') CYCLE QUANLOOP
      DO NR=1,NRAMP
      IF (ID.EQ.RAMPID(NR)) CYCLE QUANLOOP
      ENDDO
      NRAMP = NRAMP + 1
      RAMPID(NRAMP) = ID
      RAMPTYPE(NRAMP) = 'TIME'
      IF (I.GE.3 .AND. I.LE.6) RAMPTYPE(NRAMP)='TEMPERATURE'
      ENDDO QUANLOOP
C
C     Particle Information
C
      IF (PARTICLES) PART_ID = 'tracer particles'
C
      IF (PART_ID.NE.'null') THEN
         DO IPC=1,NPC
         LP=>LAGRANGIAN(IPC)
         IF (LP%CLASS_NAME.EQ.PART_ID) 
     .       SURF_PARTICLE_CLASS(N) = IPC
         ENDDO
         DROPLET_FILE=.TRUE.
         ENDIF
C
C     Ramps for velocity and heat/pyrolysis
C
      TAUQ(N)  = TAU_Q
      TAUV(N)  = TAU_V
      IF (RAMP_Q.NE.'null') THEN
         RAMPLOOP4: DO NR=1,NRAMP
         IF (RAMP_Q.EQ.RAMPID(NR)) EXIT RAMPLOOP4
         ENDDO RAMPLOOP4
         TAUQ(N) = 1000000 + NR 
         ENDIF
      IF (RAMP_V.NE.'null') THEN
         RAMPLOOP5: DO NR=1,NRAMP
         IF (RAMP_V.EQ.RAMPID(NR)) EXIT RAMPLOOP5
         ENDDO RAMPLOOP5
         TAUV(N) = 1000000 + NR
         ENDIF
C
C     Determine if the surface is combustible/burning
C
      BURNING  = .FALSE.
      IF (HRRPUA.GT.0. .OR. HEAT_OF_VAPORIZATION.GT.0.) BURNING=.TRUE.
C
      BURNING_IF: IF (BURNING) THEN
C
      IF (HRRPUA.GT.0.) THEN
         IF (HEAT_OF_COMBUSTION.GT.0.) THEN
            MASS_FLUX(IFUEL) = HRRPUA/HEAT_OF_COMBUSTION
            ELSE
            MASS_FLUX(IFUEL) = HRRPUA*1000./(DELTAH_FUEL*Y_F_INLET)
            ENDIF
         ENDIF
C
      IF (HEAT_OF_COMBUSTION.GT.0. .AND. DELTAH_FUEL.GT.0.) 
     .   ADJUST_BURN_RATE(N,IFUEL) = 
     .       (1000.*HEAT_OF_COMBUSTION)/(Y_F_INLET*DELTAH_FUEL)
C
      IF (HEAT_OF_VAPORIZATION.GT.0.) THEN
         H_V(N) = 1000.*HEAT_OF_VAPORIZATION
         MASSFLUXS_MAX(N) = BURNING_RATE_MAX
         TAU_Q = 0.
         ENDIF
C
      TAU_MF(IFUEL)  = TAU_Q
      RAMP_MF(IFUEL) = RAMP_Q
C
      ENDIF BURNING_IF
C
C     Species Arrays
C
      SPECIES_LOOP: DO NSPC=0,NSPEC
      IF (VEL.LT.0. .AND. VEL.NE.-999. .AND. 
     .    MASS_FRACTION(NSPC).LT.0. .AND. .NOT.BURNING) THEN
         MMT(N) = 2
         MASS_FRACTION(NSPC) = YY0(NSPC)
         ENDIF
      IF (VOLUME_FLUX.LT.0. .AND. VOLUME_FLUX.NE.-999. .AND.
     .    MASS_FRACTION(NSPC).LT.0. .AND. .NOT.BURNING) THEN
         MMT(N) = 2
         MASS_FRACTION(NSPC) = YY0(NSPC)
         ENDIF
      IF (MASS_FRACTION(NSPC).GE.0.) THEN
         MMT(N) = 2
         MASSFRACS(N,NSPC)  = MASS_FRACTION(NSPC)
         ENDIF
      IF (BURNING .OR. MASS_FLUX(NSPC).NE.0.) THEN
         MMT(N) = 3
         MASSFLUXS(N,NSPC) = MASS_FLUX(NSPC)
         ENDIF
      TAUMF(N,NSPC) = TAU_MF(NSPC)
      IF (RAMP_MF(NSPC).NE.'null') THEN
         RAMPLOOP6: DO NR=1,NRAMP
         IF (RAMP_MF(NSPC).EQ.RAMPID(NR)) EXIT RAMPLOOP6
         ENDDO RAMPLOOP6
         TAUMF(N,NSPC) = 1000000 + NR
         ENDIF
      ENDDO SPECIES_LOOP
C
C     Velocity Information
C
      BCV(N)        = VBC
      IF (VEL_T(1).NE.-999.) BCV(N) = 2.
      VEL_TS(:,N) = VEL_T
      VELS(N)  = VEL
      VFLUX(N) = VOLUME_FLUX
      IF (RADIUS.GT.0.)           IBACK(N) = 2
      IF (BACKING.EQ.'INSULATED') IBACK(N) = 2
      IF (BACKING.EQ.'EXPOSED')   IBACK(N) = 3
      DO NN=0,NBT
      IF (BACKING.EQ.SURFNAME(NN)) IBACK(N) = -NN
      ENDDO
      TMP_VOID(N) = TMP_BACK + TMPM
C
      IF (PROFILE.EQ.'PARABOLIC')   IPROF(N) = 1
      IF (PROFILE.EQ.'ATMOSPHERIC') THEN
         IPROF(N) = 2
         PLES(N)  = PLE
         Z0S(N)   = Z0
         ENDIF 
      IF (PROFILE.EQ.'1D-PARABOLIC') IPROF(N) = 3
C
      LEAKY(N) = LEAKING
C
C     Thermal properties. MHC is Method of Heat Conduction
C
      IF (ADIABATIC) THEN
         MHC(N) = 0
         EMISSIVITY = 0.0
         ENDIF
C
      IF (HEAT_FLUX.NE.0.) THEN
         HEATFLUX(N) = HEAT_FLUX*1000.
         MHC(N) = -1
         ENDIF
C
      C_P = C_P * 1000.
C
      IF_THICK: IF (KS.GT.0.) THEN
         MHC(N) = 3
         IF (ALPHA .GT. 0.) THEN
            IF (C_P    .GT.0.) DENSITY = KS/(ALPHA*C_P)
            IF (DENSITY.GT.0.) C_P     = KS/(ALPHA*DENSITY)
            ENDIF
         IF (C_P.LT.0.) C_P = 1000.  ! If user did not set ALPHA or C_P
      ENDIF IF_THICK
C
      DENSITY_S(N)  = DENSITY*(1.-MOISTURE_FRACTION)
      POROS(N)      = POROSITY
      QEXT(N)       = 1000.*EXTERNAL_FLUX
      TMP_P(N)      = TMPWAL + TMPM
      TMP_P0(N)     = TMPWAL0+ TMPM
      TMPMIN        = MIN(TMPMIN,TMP_P(N),TMP_P0(N))
      ECOEF(N)      = E_COEFFICIENT
      TMP_I(N)      = TMPIGN + TMPM
      TMP_E(N)      = TMPEVAP + TMPM
C
      XKS(N)      = KS
      IF (RAMP_KS.NE.'null') THEN
         MHC(N) = 3
         RAMPLOOP7: DO NR=1,NRAMP
         IF (RAMP_KS.EQ.RAMPID(NR)) EXIT RAMPLOOP7
         ENDDO RAMPLOOP7
         XKS(N) = -NR
         ENDIF
C
      XKS_C(N)      = KS_CHAR
      IF (RAMP_KS_CHAR.NE.'null') THEN
         MHC(N) = 3
         RAMPLOOP9: DO NR=1,NRAMP
         IF (RAMP_KS_CHAR.EQ.RAMPID(NR)) EXIT RAMPLOOP9
         ENDDO RAMPLOOP9
         XKS_C(N) = -NR
         ENDIF
C
      C_PV(N)      = C_P
      IF (RAMP_C_P.NE.'null') THEN
         RAMPLOOP8: DO NR=1,NRAMP
         IF (RAMP_C_P.EQ.RAMPID(NR)) EXIT RAMPLOOP8
         ENDDO RAMPLOOP8
         C_PV(N) = -NR
         ENDIF
C
      C_PC(N)      = C_P_CHAR*1000.
      IF (RAMP_C_P_CHAR.NE.'null') THEN
         RAMPLOOP10: DO NR=1,NRAMP
         IF (RAMP_C_P_CHAR.EQ.RAMPID(NR)) EXIT RAMPLOOP10
         ENDDO RAMPLOOP10
         C_PC(N) = -NR
         ENDIF
C
      IF (MHC(N).EQ.3) THEN
         IF (WALL_POINTS.EQ.0) THEN
            NWP(N) = 20
         ELSE
            NWP(N) = WALL_POINTS
         ENDIF
         IF (RADIUS.GT.0.) DELTA = RADIUS
         IF (RADIUS.GT.0.) SURF_GEOM(N) = 'CYLINDRICAL'
         DXI_W(N) = DELTA/REAL(NWP(N),EB)
         DX_W(N) = DX_SOLID
         ENDIF
C
      WALL_THICKNESS(N) = DELTA
      CDR(N) = C_DELTA_RHO*1000.
C
      IF (CDR(N).GT.0.) MHC(N) = 2
      IF (MHC(N).NE.3 .AND. DELTA.GT.0. .AND. DENSITY.GT.0. .AND.
     .   (C_P.GT.0. .OR. RAMP_C_P.NE.'null') ) THEN
        MHC(N) = 2
        IF (C_P.GT.0.) CDR(N) = C_P*DELTA*DENSITY
        ENDIF
C
C     Pyrolysis Models
C
      IF (BURN_AWAY) BURNAWAY(N)=.TRUE.
      DENSITY_F(N) = DENSITY_S(N)*FUEL_FRACTION
      A_SOLID(N)   = A
C
      IF (PAPER_MODEL) PHASE = 'PAPER'
C
      IF (PHASE.EQ.'null') THEN
         PYROLYSIS_MODEL(N) = 'GAS_BURNER'
         IF (MHC(N).GE.2 .AND. H_V(N).GT.0.)  
     .           PYROLYSIS_MODEL(N) = 'THERMOPLASTIC'
         ELSE
         PYROLYSIS_MODEL(N) = PHASE
         ENDIF
C
      IF (HEAT_OF_ABLATION.GT.0.) THEN
         H_V(N) = HEAT_OF_ABLATION*1000.
         PYROLYSIS_MODEL(N) = 'ABLATION'
         ENDIF
C
      PYRO_MODEL: SELECT CASE(PYROLYSIS_MODEL(N))
C
      CASE('GAS_BURNER')
         RHO_CHAR(N) = 0.
         IF (SURFACE_DENSITY.GT.0.) THEN
            SURFACE_DENSITY_S(N) = SURFACE_DENSITY
            ELSE
            SURFACE_DENSITY_S(N) = 1000000.      
            ENDIF
C
      CASE('THERMOPLASTIC')
         IF (SURFACE_DENSITY.GT.0.) THEN
            SURFACE_DENSITY_S(N) = SURFACE_DENSITY
            ELSE
            SURFACE_DENSITY_S(N) = DENSITY_F(N)*WALL_THICKNESS(N)
            ENDIF
         RHO_CHAR(N) = 0.
         IF (E.LT.0.) THEN
            E_SOLID(N) = -R0*TMP_I(N)*
     .       LOG( MASS_FLUX_CRITICAL/(A_SOLID(N)*DENSITY_S(N)) )
            ELSE
            E_SOLID(N) = E*1000.
            ENDIF
C
      CASE('ABLATION')
         RHO_CHAR(N) = 0.
         IF (SURFACE_DENSITY.GT.0.) THEN
            SURFACE_DENSITY_S(N) = SURFACE_DENSITY
            ELSE
            SURFACE_DENSITY_S(N) = 1000000.
            ENDIF
         IF (E.LT.0.) THEN
            E_SOLID(N) = -R0*(ABLATION_TEMPERATURE+TMPM)*
     .       LOG( ABLATION_RATE/(A_SOLID(N)*DENSITY_S(N)) )
            ELSE
            E_SOLID(N) = E*1000.
            ENDIF
C
      CASE('LIQUID')
         IF (SURFACE_DENSITY.GT.0.) THEN
            SURFACE_DENSITY_S(N) = SURFACE_DENSITY
            ELSE
            SURFACE_DENSITY_S(N) = DENSITY_F(N)*WALL_THICKNESS(N)
            ENDIF
         RHO_CHAR(N) = 0.
C
      CASE('CHAR')
         RHO_MOIS(N)  = MOISTURE_FRACTION*DENSITY
         RHO_CHAR(N)  = CHAR_DENSITY
         DENSITY_F(N) = DENSITY_F(N)*(DENSITY_S(N)-CHAR_DENSITY)/
     .                                DENSITY_S(N)
         SURFACE_DENSITY_S(N) = DENSITY_F(N)*DELTA
         IF (E.LT.0.) THEN
            E_SOLID(N) = -R0*TMP_I(N)*LOG( MASS_FLUX_CRITICAL/
     .                  (A_SOLID(N)*(DENSITY_S(N)-RHO_CHAR(N))) )
            ELSE
            E_SOLID(N) = E*1000.
            ENDIF
C
      CASE('PAPER')
         MHC(N) = 4
         MMT(N) = 4
         PAPERMODEL= .TRUE.
         TAUV(N)   = 0.
         TAUQ(N)   = 0.
         TAUMF(N,IFUEL)   = 0.
         TAUMF(N,IOXYGEN) = 0.
C
      END SELECT PYRO_MODEL
C
C     Ignition Time
C
      IF (TIGN.GE.0.) THEN
         TIGNS(N) = TIGN
         ELSE
         TIGNS(N) = 0.
         IF (MHC(N).GT.1) TIGNS(N) = 1000000.
         IF (H_V(N).GT.0. .AND. MASSFLUXS(N,IFUEL).EQ.0.) TIGNS(N) = 0.
         ENDIF
C
C     Miscellaneous Surface Info
C
      NWP_MAX   = MAX(NWP_MAX,NWP(N))
      NPPCS(N)  = NPPC
      E_WALL_S(N) = EMISSIVITY
      TEX_MAP(N) = TEXTURE_MAP
      IF (TEXTURE_MAP.NE.'null') SURF_TYPE(N) = 1
      TEX_WIDTH(N) = TEXTURE_WIDTH
      TEX_HEIGHT(N) = TEXTURE_HEIGHT
      IF (RGB4(1).GT.-0.5) THEN
      SURF_RGB(N,:) = RGB4
      ELSE
      SURF_RGB(N,1) = RGB(1)
      SURF_RGB(N,2) = RGB(2)
      SURF_RGB(N,3) = RGB(3)
      SURF_RGB(N,4) = 1.0    
      ENDIF
C
      ENDDO BCLOOP
      REWIND(LU5)
C
C Check for user-defined mins and maxes.
C
      MINIMUM_DENSITY = -999.
      MAXIMUM_DENSITY = -999.
      MINIMUM_TEMPERATURE = -999.
      MAXIMUM_TEMPERATURE = -999.
      MINIMUM_MASS_FRACTION = -999.
      MAXIMUM_MASS_FRACTION = -999.
C
      CLIP_LOOP: DO
      CALL CHECKREAD('CLIP',LU5,IOS) ; IF (IOS.EQ.1) EXIT CLIP_LOOP
      READ(LU5,CLIP,END=431,ERR=432,IOSTAT=IOS)
  432 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with CLIP line')
      ENDDO CLIP_LOOP
  431 REWIND(LU5)
C
      IF (MINIMUM_TEMPERATURE.GT.-TMPM) 
     .   TMPMIN = MINIMUM_TEMPERATURE + TMPM
C
      IF (MAXIMUM_TEMPERATURE.GT.-TMPM) 
     .   TMPMAX = MAXIMUM_TEMPERATURE + TMPM
C
      IF (MINIMUM_DENSITY.GT.0.) THEN
      RHOMIN = MINIMUM_DENSITY
      ELSE
      RHOMIN = 0.1*RHOA
      ENDIF
C
      IF (MAXIMUM_DENSITY.GT.0.) THEN
      RHOMAX = MAXIMUM_DENSITY
      ELSE
      RHOMAX = 1.2*PINF*MW_MAX/(R0*TMPMIN)
      ENDIF
C
      IF (MINIMUM_MASS_FRACTION.GT.-1.) YYMIN = MINIMUM_MASS_FRACTION
      IF (MAXIMUM_MASS_FRACTION.GT.-1.) YYMAX = MAXIMUM_MASS_FRACTION
C
C Define open boundary and mirror boundary surface characteristics
C
      OPEN_INDEX = NBT+1
      MIRROR_INDEX = NBT+2
      INTERPOLATED_INDEX = NBT+3
C
      BCV(OPEN_INDEX) = 1. ; SURFNAME(OPEN_INDEX) = 'OPEN'
      MHC(OPEN_INDEX) = 0
      MMT(OPEN_INDEX) = 1
      SURF_TYPE(OPEN_INDEX) = 2
C
      BCV(MIRROR_INDEX) = 1. ; SURFNAME(MIRROR_INDEX) = 'MIRROR'
      MHC(MIRROR_INDEX) = 0
      MMT(MIRROR_INDEX) = 1
      SURF_TYPE(MIRROR_INDEX) = -2
C
      BCV(INTERPOLATED_INDEX) = 1.  
      SURFNAME(INTERPOLATED_INDEX) = 'INTERPOLATED'
      MHC(INTERPOLATED_INDEX) = 6
      MMT(INTERPOLATED_INDEX) = 6
C
      NBT = NBT + 3
C
C Set up 1-D grids and arrays for thermally-thick calcs
C
      ALLOCATE(AAS(NWP_MAX),STAT=IZERO)
      CALL ChkMemErr('INIT','AAS',IZERO) ; AAS = 0.
      ALLOCATE(CCS(NWP_MAX),STAT=IZERO)
      CALL ChkMemErr('INIT','CCS',IZERO) ; CCS = 0.
      ALLOCATE(BBS(NWP_MAX),STAT=IZERO)
      CALL ChkMemErr('INIT','BBS',IZERO) ; BBS = 0.
      ALLOCATE(DDS(NWP_MAX),STAT=IZERO)
      CALL ChkMemErr('INIT','DDS',IZERO) ; DDS = 0.
      ALLOCATE(RDX_W(NWP_MAX,0:NBT),STAT=IZERO)
      CALL ChkMemErr('INIT','RDX_W',IZERO)
      ALLOCATE(RDXN_W(0:NWP_MAX,0:NBT),STAT=IZERO)
      CALL ChkMemErr('INIT','RDXN_W',IZERO)
      ALLOCATE(X_S(0:NWP_MAX+1,0:NBT),STAT=IZERO)
      CALL ChkMemErr('INIT','X_S',IZERO) ; X_S = 0.
      ALLOCATE(DXF(0:NBT),STAT=IZERO)
      CALL ChkMemErr('INIT','DXF',IZERO)
      ALLOCATE(DXB(0:NBT),STAT=IZERO)
      CALL ChkMemErr('INIT','DXB',IZERO)
      ALLOCATE(RHO_M(NWP_MAX),STAT=IZERO)
      CALL ChkMemErr('INIT','RHO_M',IZERO)
      ALLOCATE(RHO_S(NWP_MAX),STAT=IZERO)
      CALL ChkMemErr('INIT','RHO_S',IZERO)
      ALLOCATE(DRHOWDT(NWP_MAX),STAT=IZERO)
      CALL ChkMemErr('INIT','DRHOWDT',IZERO) ; DRHOWDT = 0.
      ALLOCATE(DRHOMDT(NWP_MAX),STAT=IZERO)
      CALL ChkMemErr('INIT','DRHOMDT',IZERO) ; DRHOMDT = 0.
      ALLOCATE(K_S(0:NWP_MAX+1),STAT=IZERO)
      CALL ChkMemErr('INIT','K_S',IZERO)
      ALLOCATE(K_S_C(0:NWP_MAX+1),STAT=IZERO)
      CALL ChkMemErr('INIT','K_S_C',IZERO)
      ALLOCATE(C_P_V(0:NWP_MAX+1),STAT=IZERO)
      CALL ChkMemErr('INIT','C_P_V',IZERO)
      ALLOCATE(C_P_C(0:NWP_MAX+1),STAT=IZERO)
      CALL ChkMemErr('INIT','C_P_C',IZERO)
      ALLOCATE(RHOCBAR(NWP_MAX),STAT=IZERO)
      CALL ChkMemErr('INIT','RHOCBAR',IZERO)
      ALLOCATE(C_S(NWP_MAX),STAT=IZERO)
      CALL ChkMemErr('INIT','C_S',IZERO)
      ALLOCATE(D_S(NWP_MAX),STAT=IZERO)
      CALL ChkMemErr('INIT','D_S',IZERO)
C
C Set up wall conduction arrays
C
      BLOOP: DO IBC=0,NBT
      IF (MHC(IBC).NE.3) CYCLE BLOOP
C
C     Find stretch factor SFAC to ensure desired first grid cell size
C     Solve DX_W=WIDTH*(EXP(SFAC/NWP)-1)/(EXP(SFAC)-1) for SFAC
C
      WIDTH = DXI_W(IBC)*NWP(IBC)
      SFAC  = 4.
      IF (DX_W(IBC).LT.0.001) SFAC = 6.
      IF (DX_W(IBC).LT.0.0001) SFAC = 8.
C
      IF (DXI_W(IBC).LE.DX_W(IBC)) THEN
      SFAC=0.001   ! Just use nearly-linear x=xi
      ELSE
      DO ITER=1,50 ! Iterate to find SFAC using Newton's Method
      LS   = NWP(IBC)*WIDTH*(EXP(SFAC/NWP(IBC))-1.)-
     .       NWP(IBC)*DX_W(IBC)*(EXP(SFAC)-1.)
      DLDS = WIDTH*(EXP(SFAC/NWP(IBC))-1.)-NWP(IBC)*DX_W(IBC)*EXP(SFAC)
      SFAC = SFAC - LS/DLDS
      ENDDO
      ENDIF
C
      CONST = WIDTH/(EXP(SFAC)-1.)
C
      DO I=1,NWP(IBC)
      X1 = CONST*(EXP(SFAC*(I-1)/NWP(IBC))-1.)
      X2 = CONST*(EXP(SFAC*(I  )/NWP(IBC))-1.)
      RDX_W(I,IBC) = 1./(X2-X1)
      IF (SURF_GEOM(IBC).EQ.'CYLINDRICAL')
     .   RDX_W(I,IBC) = RDX_W(I,IBC)/(WIDTH-0.5*(X1+X2))
      ENDDO
      DO I=0,NWP(IBC)
      X0 = CONST*(EXP(SFAC*(I-1)/NWP(IBC))-1.)
      X1 = CONST*(EXP(SFAC*(I  )/NWP(IBC))-1.)
      X2 = CONST*(EXP(SFAC*(I+1)/NWP(IBC))-1.)
      X0 = 0.5*(X0+X1)
      X2 = 0.5*(X1+X2)
      X_S(I,IBC) = X1
      RDXN_W(I,IBC) = 1./(X2-X0)
      IF (I.EQ.0)        DXF(IBC) = (X2-X0)
      IF (I.EQ.NWP(IBC)) DXB(IBC) = (X2-X0)
      IF (SURF_GEOM(IBC).EQ.'CYLINDRICAL')
     .   RDXN_W(I,IBC) = RDXN_W(I,IBC)*(WIDTH-X1)
      ENDDO
C
      ENDDO BLOOP
C
      END SUBROUTINE READ_SURF
C
C
      SUBROUTINE SET_SURF_DEFAULTS
C
      ID                 = 'null'
      AKA                = 'null'
      HRRPUA             = 0.
      HEAT_OF_VAPORIZATION = 0.
      HEAT_OF_COMBUSTION = -1.
      HEAT_OF_ABLATION   = -1.
      ABLATION_TEMPERATURE = 1000.
      ABLATION_RATE      = 0.05
      BURNING_RATE_MAX   = 0.1
      EMISSIVITY         = 0.9
      BURN_AWAY          = .FALSE.
      RGB(1)=1. ; RGB(2)=0.8 ; RGB(3)=0.4
      RGB4(1)=-1.;RGB4(2)=-1.; RGB4(3)=-1.;RGB4(4)=1.
      TEXTURE_MAP        = 'null'
      TEXTURE_WIDTH      = 1.
      TEXTURE_HEIGHT     = 1.
      PHASE              = 'null'  
      RADIUS             = -1.
      MASS_FLUX          = 0.
      MASS_FLUX_CRITICAL = 0.02
c      A_PYR              = 2.6E8
      A                  = 6.5E5    ! m/s
      E                  = -1.
      HEAT_FLUX          = 0.
      WALL_POINTS        = 0
      MASS_FRACTION      = -1.
      TAU_MF             =  1.
      IF (LES) VBC0      =  0.5 
      IF (DNS) VBC0      = -1.0
      VBC                = VBC0
      TMPIGN             = 5000.
      TMPEVAP            = 100.
      VEL_T              = -999.
      VEL                = -999.
      VOLUME_FLUX        = -999.
      TIGN               = -1.
      TMPWAL             = TMPA-TMPM 
      TMPWAL0            = TMPA-TMPM 
      TMP_BACK           = TMPA-TMPM
      ALPHA              = -1.
      C_P                = -1.
      C_P_CHAR           = 0.7
      KS                 = -1.
      KS_CHAR            = 0.1
      DELTA              = 0.1
      C_DELTA_RHO        = -1.
      POROSITY           = 0.
      E_COEFFICIENT      = 0.
      TAU_Q              = 1.
      TAU_V              = 1.
      RAMP_Q             = 'null'
      RAMP_V             = 'null'
      RAMP_KS            = 'null'
      RAMP_KS_CHAR       = 'null'
      RAMP_C_P           = 'null'
      RAMP_C_P_CHAR      = 'null'
      RAMP_MF            = 'null'
      SURFACE_DENSITY    = -1.     
      DENSITY            = 450.
      CHAR_DENSITY       = 120.
      MOISTURE_FRACTION  = 0.
      FUEL_FRACTION      = 1.
      PART_ID            = 'null'
      NPPC               = 1
      PROFILE            = 'null'
      PLE                = 0.3
      Z0                 = 10.
      PAPER_MODEL        = .FALSE.
      BACKING            = 'VOID'
      ADIABATIC          = .FALSE.
      LEAKING            = .FALSE.
      PARTICLES          = .FALSE.
      DX_SOLID           = 0.0001
      EXTERNAL_FLUX      = 0.
C
      END SUBROUTINE SET_SURF_DEFAULTS
C
C
      SUBROUTINE READ_RAMP
C
      REAL(EB) F3(3),TT_MAX
      INTEGER IC
      NAMELIST /RAMP/ T,F,ID,FYI,F3
C
      IF (NRAMP.EQ.0) RETURN
C
      ALLOCATE(NTT(NRAMP)) ; NTT = 0
      ALLOCATE(TT(50000,NRAMP))
      ALLOCATE(FF(50000,NRAMP))
      ALLOCATE(FF3(50000,3)) ; FF3 = 0.
      ALLOCATE(FF2(10001,NRAMP)) ; FF2 = 0.
      ALLOCATE(FF2G(10001,3)) ; FF2G = 0.
C
      REWIND(LU5)
      IF (DATABASE.NE.'null') REWIND(LU80)
C
      RAMPLOOP: DO N=1,NRAMP
C
      SEARCH_LOOP: DO ITER=1,NFILES
      IF (ITER.EQ.1) LUDUM = LU5
      IF (ITER.EQ.2) LUDUM = LU80
      REWIND(LUDUM)
C
      NTTLOOP: DO 
      CALL CHECKREAD('RAMP',LUDUM,IOS) ; IF (IOS.EQ.1) EXIT NTTLOOP
      READ(LUDUM,NML=RAMP,ERR=56,IOSTAT=IOS)
      IF (ID.NE.RAMPID(N)) CYCLE NTTLOOP
      NTT(N) = NTT(N) + 1
      TT(NTT(N),N) = T
      IF (N.NE.IRAMPG) FF(NTT(N),N) = F
      IF (N.EQ.IRAMPG) THEN
         FF3(NTT(N),1) = F3(1)
         FF3(NTT(N),2) = F3(2)
         FF3(NTT(N),3) = F3(3)
         ENDIF
   56 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with RAMP line')
      ENDDO NTTLOOP
C
      REWIND(LUDUM)
      ENDDO SEARCH_LOOP
C
      IF (NTT(N).EQ.0) THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: RAMP ',TRIM(RAMPID(N)),
     .                            ' not found'
         CALL SHUTDOWN(MESSAGE)
         ENDIF
C
      ENDDO RAMPLOOP
C
C Extend ramp function
C
      DO N=1,NRAMP
      IF (RAMPTYPE(N).EQ.'TEMPERATURE') 
     .   TT(1:NTT(N),N) = TT(1:NTT(N),N) + TMPM
      IF (TT(1,N).GT.0. .AND. N.NE.IRAMPG) THEN
         DO I=NTT(N)+1,2,-1
         TT(I,N) = TT(I-1,N)
         FF(I,N) = FF(I-1,N)
         ENDDO
         TT(1,N) = 0.
         FF(1,N) = FF(2,N)
         NTT(N)  = NTT(N)+1
         ENDIF
      IF (RAMPTYPE(N).EQ.'TEMPERATURE') TT_MAX = 10000.
      IF (RAMPTYPE(N).EQ.'TIME')        TT_MAX = TWFIN
      IF (TT(NTT(N),N).LT.TT_MAX) THEN
         NTT(N) = NTT(N) + 1
         TT(NTT(N),N) = TT_MAX
         IF (N.NE.IRAMPG) FF(NTT(N),N) = FF(NTT(N)-1,N)
         IF (N.EQ.IRAMPG) FF3(NTT(N),:) = FF3(NTT(N)-1,:)
         ENDIF
      ENDDO
C
C Set up interpolated ramp values in FF2
C
      FF2  = 0.
      FF2G = 0.
      DTRAMP  = TT_MAX/10000.
      RDTRAMP = 1./DTRAMP
      DO N=1,NRAMP
      FF2(1,N)    = FF(1,N)
      FF2G(1,1:3) = FF3(1,1:3)
      DO I=2,10001
      TM = DTRAMP*(I-1)
      TLOOP: DO II=1,NTT(N)-1
      IF (TM.GE.TT(II,N) .AND. TM.LT.TT(II+1,N)) THEN
         IF (IRAMPG.NE.N) THEN
         FF2(I,N) = FF(II,N) +
     .      (TM-TT(II,N))*(FF(II+1,N)-FF(II,N))/(TT(II+1,N)-TT(II,N))
         ELSE 
         DO IC=1,3
         FF2G(I,IC) = FF3(II,IC) + (TM-TT(II,N))*
     .                (FF3(II+1,IC)-FF3(II,IC))/(TT(II+1,N)-TT(II,N))
         ENDDO
         ENDIF
         EXIT TLOOP
         ENDIF
      ENDDO TLOOP
      ENDDO
      ENDDO
C
C Close the Database file
C
      IF (DATABASE.NE.'null') CLOSE(LU80)
C
      END SUBROUTINE READ_RAMP
C
C
      SUBROUTINE READ_OBST
C
      INTEGER :: NM,NOM,NBO,NNN,IC
      CHARACTER(26) :: HEAT_REMOVE,HEAT_CREATE
      CHARACTER(60) :: MESH_ID
      CHARACTER(11) COLOR,BLOCK_COLOR
      REAL(EB) :: T_REMOVE,T_CREATE,T_OPEN
      EQUIVALENCE(T_OPEN,T_REMOVE)
      EQUIVALENCE(COLOR,BLOCK_COLOR)
      LOGICAL :: SAWTOOTH,EMBEDDED,THICKEN,PERMIT_HOLE,EVACUATION
      NAMELIST /OBST/ XB,SURF_ID,SURF_IDS,SURF_ID6,FYI,BLOCK_COLOR,
     .                T_REMOVE,T_CREATE,T_OPEN,BNDF_FACE,BNDF_BLOCK,
     .                SAWTOOTH,RGB,RGB4,TEXTURE_ORIGIN,THICKEN,
     .                OUTLINE,HEAT_REMOVE,HEAT_CREATE,COLOR,
     .                PERMIT_HOLE,EVACUATION,MESH_ID
C
      MESH_LOOP: DO NM=1,NMESHES
      M=>MESH(NM)
      CALL UNPACK_VAR(NM)
C
      NBT = 0
C
C Evacuation boundary type added to list if defined
C
      IF (EVAC_SURF_DEFAULT.NE.SURF_DEFAULT) THEN
        NBT = NBT + 1
        SURFNAME(NBT) = EVAC_SURF_DEFAULT
      END IF
C
C Read OBST lines and look for SURF's and add to list
C
      NB = 0
      BLOOP: DO
      SURF_ID  = 'null'
      SURF_IDS = 'null'
      SURF_ID6 = 'null'
      CALL CHECKREAD('OBST',LU5,IOS) ; IF (IOS.EQ.1) EXIT BLOOP
      READ(LU5,NML=OBST,END=1,ERR=2,IOSTAT=IOS)
      NB = NB + 1
      IF (SURF_ID.NE.'null') THEN
         DO N=0,NBT
         IF (SURF_ID.EQ.SURFNAME(N)) CYCLE BLOOP
         ENDDO
         NBT = NBT + 1
         SURFNAME(NBT) = SURF_ID
         ENDIF
      NNLOOP: DO NN=1,3
      IF (SURF_IDS(NN).NE.'null') THEN
         DO N=0,NBT
         IF (SURF_IDS(NN).EQ.SURFNAME(N)) CYCLE NNLOOP
         ENDDO
         NBT = NBT + 1
         SURFNAME(NBT) = SURF_IDS(NN)
         ENDIF
      ENDDO NNLOOP
      NNLOOP2: DO NN=1,6
      IF (SURF_ID6(NN).NE.'null') THEN
         DO N=0,NBT
         IF (SURF_ID6(NN).EQ.SURFNAME(N)) CYCLE NNLOOP2
         ENDDO
         NBT = NBT + 1
         SURFNAME(NBT) = SURF_ID6(NN)
         ENDIF
      ENDDO NNLOOP2
    2 IF (IOS.GT.0) THEN
         WRITE(MESSAGE,'(A,I5)') 
     .   'ERROR: Problem with OBSTruction number',NB+1
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO BLOOP
    1 REWIND(LU5)
C
      ALLOCATE(M%OBSTRUCTION(0:NB),STAT=IZERO)
      CALL ChkMemErr('READ','OBSTRUCTION',IZERO)
      OBSTRUCTION=>M%OBSTRUCTION
C
      OBSTRUCTION(0)   %IBC = 0
      OBSTRUCTION(0:NB)%HEAT_INDEX_CREATE = 0
      OBSTRUCTION(0:NB)%HEAT_INDEX_REMOVE = 0
      OBSTRUCTION(0:NB)%T_REMOVE = 1000000.
      OBSTRUCTION(0:NB)%T_CREATE = 1000000.
      OBSTRUCTION(0:NB)%HIDDEN   = .FALSE.
      OBSTRUCTION(0:NB)%BCI      = -1     
      OBSTRUCTION(0:NB)%BTI      = -1     
      OBSTRUCTION(0)   %RGB      = 0.     
      OBSTRUCTION(0)   %RGB(4)   = 1.     
      OBSTRUCTION(0)   %TEXTURE  = 0.     
      OBSTRUCTION(0)   %SHOW_BNDF = BNDF_DEFAULT
      OBSTRUCTION(0:NB)%SAWTOOTH  = .TRUE.      
      OBSTRUCTION(0:NB)%PERMIT_HOLE = .TRUE.      
      OBSTRUCTION(0)   %DIMENSIONS = 0           
C
      N   = 0
      NBO = NB
C
      OBSTLOOP: DO NN=1,NBO
      N        = N + 1
      SURF_ID  = 'null'
      SURF_IDS = 'null'
      SURF_ID6 = 'null'
      BLOCK_COLOR = 'null'
      MESH_ID     = 'null'
      RGB         = -1.
      RGB4        = -1.
      T_REMOVE    = 1000000.
      T_CREATE    = 1000000.
      BNDF_FACE   = BNDF_DEFAULT
      BNDF_BLOCK  = BNDF_DEFAULT
      SAWTOOTH    = .TRUE.
      THICKEN     = THICKEN_OBSTRUCTIONS
      OUTLINE     = .FALSE.
      TEXTURE_ORIGIN = -999.
      HEAT_REMOVE = 'null'
      HEAT_CREATE = 'null'
      PERMIT_HOLE = .TRUE.
      IF (.NOT.EVACUATION_ONLY(NM)) EVACUATION = .FALSE.
      IF (     EVACUATION_ONLY(NM)) EVACUATION = .TRUE.
      IF (     EVACUATION_ONLY(NM)) THICKEN    = .TRUE.
      IF (     EVACUATION_ONLY(NM)) SAWTOOTH   = .FALSE.
C
      CALL CHECKREAD('OBST',LU5,IOS) ; IF (IOS.EQ.1) EXIT OBSTLOOP
      READ(LU5,OBST,END=35,ERR=36)
C
C Evacuation criteria
C
      IF (MESH_ID.NE.MESH_NAME(NM) .AND. MESH_ID.NE.'null') THEN
          N = N-1
          NB= NB-1
          CYCLE OBSTLOOP
          ENDIF
C
      IF ((.NOT.EVACUATION .AND. EVACUATION_ONLY(NM)) .OR.
     .    (EVACUATION .AND. .NOT.EVACUATION_ONLY(NM))) THEN
          N = N-1
          NB= NB-1
          CYCLE OBSTLOOP
          ENDIF
C
C Reorder coords if necessary
C
      DO I=1,5,2
      IF (XB(I).GT.XB(I+1)) THEN
         DUMMY   = XB(I)
         XB(I)   = XB(I+1)
         XB(I+1) = DUMMY
         ENDIF
      ENDDO
C
      XB(1) = MAX(XB(1),XS)
      XB(2) = MIN(XB(2),XF)
      XB(3) = MAX(XB(3),YS)
      XB(4) = MIN(XB(4),YF)
      XB(5) = MAX(XB(5),ZS)
      XB(6) = MIN(XB(6),ZF)
      IF (XB(1).GT.XF .OR. XB(2).LT.XS .OR.
     .    XB(3).GT.YF .OR. XB(4).LT.YS .OR.
     .    XB(5).GT.ZF .OR. XB(6).LT.ZS) THEN
          N = N-1
          NB= NB-1
          CYCLE OBSTLOOP
          ENDIF
C
C Begin processing of OBSTruction
C
      OB=>OBSTRUCTION(N)
C
      OB%X1 = XB(1)
      OB%X2 = XB(2)
      OB%Y1 = XB(3)
      OB%Y2 = XB(4)
      OB%Z1 = XB(5)
      OB%Z2 = XB(6)
C
      OB%I1 = NINT( GINV(XB(1)-XS,1,NM)*RDXI   ) 
      OB%I2 = NINT( GINV(XB(2)-XS,1,NM)*RDXI   )
      OB%J1 = NINT( GINV(XB(3)-YS,2,NM)*RDETA  ) 
      OB%J2 = NINT( GINV(XB(4)-YS,2,NM)*RDETA  )
      OB%K1 = NINT( GINV(XB(5)-ZS,3,NM)*RDZETA ) 
      OB%K2 = NINT( GINV(XB(6)-ZS,3,NM)*RDZETA )
C
C     If desired, thicken small obstructions
C
      IF (THICKEN .AND. OB%I1.EQ.OB%I2) THEN
         OB%I1 = GINV(.5*(XB(1)+XB(2))-XS,1,NM)*RDXI
         OB%I2 = MIN(OB%I1+1,IBAR)
         ENDIF
      IF (THICKEN .AND. OB%J1.EQ.OB%J2) THEN
         OB%J1 = GINV(.5*(XB(3)+XB(4))-YS,2,NM)*RDETA
         OB%J2 = MIN(OB%J1+1,JBAR)
         ENDIF
      IF (THICKEN .AND. OB%K1.EQ.OB%K2) THEN
         OB%K1 = GINV(.5*(XB(5)+XB(6))-ZS,3,NM)*RDZETA
         OB%K2 = MIN(OB%K1+1,KBAR)
         ENDIF
C
C     Throw out obstructions that are too small
C
      IF ((OB%I1.EQ.OB%I2 .AND. OB%J1.EQ.OB%J2) .OR.
     .    (OB%I1.EQ.OB%I2 .AND. OB%K1.EQ.OB%K2) .OR.
     .    (OB%J1.EQ.OB%J2 .AND. OB%K1.EQ.OB%K2)) THEN
         N = N-1
         NB= NB-1
         CYCLE OBSTLOOP
         ENDIF
C
C     Check to see if obstacle is completely embedded in another
C
      EMBEDDED = .FALSE.
      EMBED_LOOP: DO NNN=1,N-1
      OB2=>OBSTRUCTION(NNN)
      IF (OB%I1.GE.OB2%I1 .AND. OB%I2.LE.OB2%I2 .AND.
     .    OB%J1.GE.OB2%J1 .AND. OB%J2.LE.OB2%J2 .AND.
     .    OB%K1.GE.OB2%K1 .AND. OB%K2.LE.OB2%K2) THEN
         EMBEDDED = .TRUE.
         EXIT EMBED_LOOP
         ENDIF
      ENDDO EMBED_LOOP
C
      IF (EMBEDDED .AND. T_CREATE.GT.100000. .AND. 
     .                   T_REMOVE.GT.100000. .AND.
     .                   PERMIT_HOLE) THEN
          N = N-1
          NB= NB-1
          CYCLE OBSTLOOP
          ENDIF
C
C     Save boundary condition info for obstacles
C
      OB%IBC(:) = 0
C
      IF (EVACUATION_ONLY(NM)) SURF_ID=EVAC_SURF_DEFAULT
C
      DO NNN=0,NBT
      IF (SURF_ID    .EQ.SURFNAME(NNN)) OB%IBC(:)    = NNN
      IF (SURF_IDS(1).EQ.SURFNAME(NNN)) OB%IBC(3)    = NNN
      IF (SURF_IDS(2).EQ.SURFNAME(NNN)) OB%IBC(-2:2) = NNN
      IF (SURF_IDS(3).EQ.SURFNAME(NNN)) OB%IBC(-3)   = NNN
      IF (SURF_ID6(1).EQ.SURFNAME(NNN)) OB%IBC(-1)   = NNN
      IF (SURF_ID6(2).EQ.SURFNAME(NNN)) OB%IBC( 1)   = NNN
      IF (SURF_ID6(3).EQ.SURFNAME(NNN)) OB%IBC(-2)   = NNN
      IF (SURF_ID6(4).EQ.SURFNAME(NNN)) OB%IBC( 2)   = NNN
      IF (SURF_ID6(5).EQ.SURFNAME(NNN)) OB%IBC(-3)   = NNN
      IF (SURF_ID6(6).EQ.SURFNAME(NNN)) OB%IBC( 3)   = NNN
      ENDDO
C
C     Creation and removal logic
C
      OB%HEAT_REMOVE = HEAT_REMOVE 
      OB%HEAT_CREATE = HEAT_CREATE 
C
      IF (T_CREATE.LT.0.) THEN
         OB%HEAT_INDEX_CREATE = NINT(ABS(T_CREATE))
      ELSE
         OB%T_CREATE = T_CREATE
      ENDIF
C
      IF (T_REMOVE.LT.0.) THEN
         OB%HEAT_INDEX_REMOVE = NINT(ABS(T_REMOVE))
      ELSE
         OB%T_REMOVE = T_REMOVE
      ENDIF
C
C     Choose obstruction color index
C
      SELECT CASE(BLOCK_COLOR)
      CASE('WHITE')     ; OB%BCI =  0
      CASE('YELLOW')    ; OB%BCI =  1
      CASE('BLUE')      ; OB%BCI =  2
      CASE('RED')       ; OB%BCI =  3
      CASE('GREEN')     ; OB%BCI =  4
      CASE('MAGENTA')   ; OB%BCI =  5
      CASE('CYAN')      ; OB%BCI =  6
      CASE('BLACK')     ; OB%BCI =  7
      CASE('INVISIBLE') ; OB%BCI = -2
      CASE DEFAULT      ; OB%BCI = -1
      END SELECT
C
      IF (RGB(1).GE.0. .AND. RGB(2).GE.0. .AND. RGB(3).GE.0.) THEN
         OB%RGB(1) = RGB(1)
         OB%RGB(2) = RGB(2)
         OB%RGB(3) = RGB(3)
         OB%RGB(4) = 1.0   
         OB%BCI    = -3
         ENDIF
C
      IF (RGB4(1).GE.0.) THEN
         OB%RGB(1:4) = RGB4(1:4)
         OB%BCI    = -3
         ENDIF
C
C     Miscellaneous assignments
C
      OB%TEXTURE(:) = TEXTURE_ORIGIN(:)  ! Origin of texture map
      OB%ORDINAL = NN  ! Order of OBST in original input file
      OB%PERMIT_HOLE = PERMIT_HOLE
C
C     Make obstruction invisible if it's within a finer mesh
C
      DO NOM=1,NM-1
      IF (XB(1).GE.MESH(NOM)%XS .AND. XB(2).LE.MESH(NOM)%XF .AND.
     .    XB(3).GE.MESH(NOM)%YS .AND. XB(4).LE.MESH(NOM)%YF .AND.
     .    XB(5).GE.MESH(NOM)%ZS .AND. XB(6).LE.MESH(NOM)%ZF) OB%BCI=-2
      ENDDO
C
C     Prevent drawing of boundary info if desired
C
      IF (BNDF_DEFAULT) THEN
      OB%SHOW_BNDF(:) = BNDF_FACE(:)
      IF (.NOT.BNDF_BLOCK) OB%SHOW_BNDF(:) = .FALSE.
      ELSE
      OB%SHOW_BNDF(:) = BNDF_FACE(:)
      IF (BNDF_BLOCK) OB%SHOW_BNDF(:) = .TRUE.
      ENDIF
C
C     Smooth obstacles if desired
C
      IF (.NOT.SAWTOOTH) THEN
         OB%BTI = 3
         OB%SAWTOOTH = .FALSE.
         ENDIF
C
C     In Smokeview, draw the outline of the obstruction
C
      IF (OUTLINE) OB%BTI = 2
C
   36 ENDDO OBSTLOOP
   35 REWIND(LU5)
C
      ENDDO MESH_LOOP
C
C Read HOLEs and cut out blocks
C
      CALL READ_HOLE
C
C Go through all meshes, recording which cells are solid
C
      MESH_LOOP2: DO NM=1,NMESHES
      M=>MESH(NM)
      CALL UNPACK_VAR(NM)
C
C     Compute areas of obstruction faces, both actual (AB0) and
C     FDS approximated (AB)
C
      DO N=1,NB
      OB=>OBSTRUCTION(N)
      OB%AREA_0(1) = (OB%Y2-OB%Y1)*(OB%Z2-OB%Z1)
      OB%AREA_0(2) = (OB%X2-OB%X1)*(OB%Z2-OB%Z1)
      OB%AREA_0(3) = (OB%X2-OB%X1)*(OB%Y2-OB%Y1)
      OB%AREA(1)   = (Y(OB%J2)-Y(OB%J1))*(Z(OB%K2)-Z(OB%K1))
      OB%AREA(2)   = (X(OB%I2)-X(OB%I1))*(Z(OB%K2)-Z(OB%K1))
      OB%AREA(3)   = (X(OB%I2)-X(OB%I1))*(Y(OB%J2)-Y(OB%J1))
      OB%DIMENSIONS(1) = OB%I2 - OB%I1
      OB%DIMENSIONS(2) = OB%J2 - OB%J1
      OB%DIMENSIONS(3) = OB%K2 - OB%K1
      IF (OB%T_CREATE.LT.100000.   .OR. 
     .    OB%HEAT_CREATE.NE.'null' .OR. 
     .    OB%HEAT_INDEX_CREATE.GT.0)  OB%HIDDEN=.TRUE.
      IF (OB%T_REMOVE.LT.OB%T_CREATE) OB%HIDDEN=.FALSE.
      ENDDO
C
C     Create main blockage index array (ICA)
C
      ALLOCATE(M%ICA(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO) 
      CALL ChkMemErr('READ','ICA',IZERO) 
      ICA=>M%ICA
C
      ICA  = 0 
      NDBC = 0
C
      DO K=0,KBP1
      DO J=0,JBP1
      DO I=0,1
      IF (ICA(I,J,K).EQ.0) THEN
         NDBC = NDBC + 1
         ICA(I,J,K) = NDBC
         ENDIF
      ENDDO
      DO I=IBAR,IBP1
      IF (ICA(I,J,K).EQ.0) THEN
         NDBC = NDBC + 1
         ICA(I,J,K) = NDBC
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      DO K=0,KBP1
      DO I=0,IBP1
      DO J=0,1
      IF (ICA(I,J,K).EQ.0) THEN
         NDBC = NDBC + 1
         ICA(I,J,K) = NDBC
         ENDIF
      ENDDO
      DO J=JBAR,JBP1
      IF (ICA(I,J,K).EQ.0) THEN
         NDBC = NDBC + 1
         ICA(I,J,K) = NDBC
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      DO J=0,JBP1
      DO I=0,IBP1
      DO K=0,1
      IF (ICA(I,J,K).EQ.0) THEN
         NDBC = NDBC + 1
         ICA(I,J,K) = NDBC
         ENDIF
      ENDDO
      DO K=KBAR,KBP1
      IF (ICA(I,J,K).EQ.0) THEN
         NDBC = NDBC + 1
         ICA(I,J,K) = NDBC
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      DO N=1,NB
      OB=>OBSTRUCTION(N)
      DO K=OB%K1,OB%K2+1
      DO J=OB%J1,OB%J2+1
      DO I=OB%I1,OB%I2+1
      IF (ICA(I,J,K).EQ.0) THEN
         NDBC = NDBC + 1
         ICA(I,J,K) = NDBC
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
C     Store in SOLID which cells are solid and which are not
C
      ALLOCATE(M%SOLID(0:NDBC),STAT=IZERO) 
      CALL ChkMemErr('READ','SOLID',IZERO) ; M%SOLID = .FALSE.
      SOLID=>M%SOLID
C
      CALL BLKCLL(NM,   0,   0,   1,JBAR,   1,KBAR,1)
      CALL BLKCLL(NM,IBP1,IBP1,   1,JBAR,   1,KBAR,1)
      IF (TWO_D) THEN
      CALL BLKCLL(NM,   0,IBP1,   0,   0,   0,KBP1,1)
      CALL BLKCLL(NM,   0,IBP1,JBP1,JBP1,   0,KBP1,1)
      ELSE
      CALL BLKCLL(NM,   1,IBAR,   0,   0,   1,KBAR,1)
      CALL BLKCLL(NM,   1,IBAR,JBP1,JBP1,   1,KBAR,1)
      ENDIF
      CALL BLKCLL(NM,   1,IBAR,   1,JBAR,   0,   0,1)
      CALL BLKCLL(NM,   1,IBAR,   1,JBAR,KBP1,KBP1,1)
C
      DO N=1,NB
      OB=>OBSTRUCTION(N)
      IF (.NOT.OB%HIDDEN)
     .   CALL BLKCLL(NM,OB%I1+1,OB%I2,OB%J1+1,OB%J2,OB%K1+1,OB%K2,1)
      ENDDO
C
      ALLOCATE(M%I_CELL(NDBC),STAT=IZERO) 
      CALL ChkMemErr('READ','I_CELL',IZERO) ; M%I_CELL = -1
      ALLOCATE(M%J_CELL(NDBC),STAT=IZERO) 
      CALL ChkMemErr('READ','J_CELL',IZERO) ; M%J_CELL = -1
      ALLOCATE(M%K_CELL(NDBC),STAT=IZERO) 
      CALL ChkMemErr('READ','K_CELL',IZERO) ; M%K_CELL = -1
      I_CELL=>M%I_CELL ; J_CELL=>M%J_CELL ; K_CELL=>M%K_CELL
C
      DO K=0,KBP1
      DO J=0,JBP1
      DO I=0,IBP1
      IC = ICA(I,J,K)
      IF (IC.GT.0) THEN
         I_CELL(IC) = I
         J_CELL(IC) = J
         K_CELL(IC) = K
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      ENDDO MESH_LOOP2
C
      END SUBROUTINE READ_OBST
C
C
      SUBROUTINE READ_HOLE
C
      CHARACTER(26) :: HEAT_REMOVE,HEAT_CREATE
      CHARACTER(60) :: MESH_ID
      LOGICAL :: EVACUATION
      INTEGER :: NM,NH,NN,NDO
      REAL(EB) :: X1,X2,Y1,Y2,Z1,Z2,T_REMOVE,T_CREATE
      NAMELIST /HOLE/ XB,FYI,T_REMOVE,T_CREATE,RGB,RGB4,
     .                HEAT_REMOVE,HEAT_CREATE,EVACUATION,MESH_ID
      TYPE(OBSTRUCTION_TYPE), DIMENSION(0:6) :: TEMP_OBST
C
      NH  = 0
      COUNT_LOOP: DO
      CALL CHECKREAD('HOLE',LU5,IOS) ; IF (IOS.EQ.1) EXIT COUNT_LOOP
      READ(LU5,NML=HOLE,END=1,ERR=2,IOSTAT=IOS)
      NH = NH + 1
    2 IF (IOS.GT.0) THEN
         WRITE(MESSAGE,'(A,I5)') 
     .   'ERROR: Problem with HOLE number',NH+1
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO COUNT_LOOP
    1 REWIND(LU5)
C
      READ_HOLE_LOOP: DO N=1,NH
C
      T_REMOVE = 1000000.
      T_CREATE = 0.
      HEAT_REMOVE = 'null'
      HEAT_CREATE = 'null'
      MESH_ID  = 'null'
      RGB      = -1.
      RGB4     = -1.
      EVACUATION = .FALSE.
C
      CALL CHECKREAD('HOLE',LU5,IOS) ; IF (IOS.EQ.1) EXIT READ_HOLE_LOOP
      READ(LU5,HOLE)
C
      DO I=1,5,2
      IF (XB(I).GT.XB(I+1)) THEN
         DUMMY   = XB(I)
         XB(I)   = XB(I+1)
         XB(I+1) = DUMMY
         ENDIF
      ENDDO
C
      MESH_LOOP: DO NM=1,NMESHES
      M=>MESH(NM)
      CALL UNPACK_VAR(NM)
C
C Evacuation criteria
C
      IF (MESH_ID.NE.'null' .AND. MESH_ID.NE.MESH_NAME(NM))
     .   CYCLE MESH_LOOP
C
      IF ((.NOT.EVACUATION .AND. EVACUATION_ONLY(NM)) .OR.
     .    (EVACUATION .AND. .NOT.EVACUATION_ONLY(NM))) 
     .   CYCLE MESH_LOOP
C
C Check if hole is contained within the current mesh
C
      X1 = XB(1)
      X2 = XB(2)
      Y1 = XB(3)
      Y2 = XB(4)
      Z1 = XB(5)
      Z2 = XB(6)
C
      IF (X1.GE.XF .OR. X2.LE.XS .OR.
     .    Y1.GT.YF .OR. Y2.LE.YS .OR.
     .    Z1.GT.ZF .OR. Z2.LE.ZS) CYCLE MESH_LOOP
C
      X1 = MAX(X1,XS-0.001*DX(0))
      X2 = MIN(X2,XF+0.001*DX(IBP1))
      Y1 = MAX(Y1,YS-0.001*DY(0))
      Y2 = MIN(Y2,YF+0.001*DY(JBP1))
      Z1 = MAX(Z1,ZS-0.001*DZ(0))
      Z2 = MIN(Z2,ZF+0.001*DZ(KBP1))
C
      I1 = NINT( GINV(XB(1)-XS,1,NM)*RDXI   ) 
      I2 = NINT( GINV(XB(2)-XS,1,NM)*RDXI   )
      J1 = NINT( GINV(XB(3)-YS,2,NM)*RDETA  ) 
      J2 = NINT( GINV(XB(4)-YS,2,NM)*RDETA  )
      K1 = NINT( GINV(XB(5)-ZS,3,NM)*RDZETA ) 
      K2 = NINT( GINV(XB(6)-ZS,3,NM)*RDZETA )
C
      NN=0
      OBST_LOOP: DO
      NN=NN+1
      IF (NN.GT.NB) EXIT OBST_LOOP
      OB=>OBSTRUCTION(NN)
      IF (.NOT.OB%PERMIT_HOLE) CYCLE OBST_LOOP
C
C     TEMP_OBST(0) is the intersection of HOLE and OBST
C
      TEMP_OBST(0)    = OBSTRUCTION(NN)
C
      TEMP_OBST(0)%I1 = MAX(I1,OB%I1)
      TEMP_OBST(0)%I2 = MIN(I2,OB%I2)
      TEMP_OBST(0)%J1 = MAX(J1,OB%J1)
      TEMP_OBST(0)%J2 = MIN(J2,OB%J2)
      TEMP_OBST(0)%K1 = MAX(K1,OB%K1)
      TEMP_OBST(0)%K2 = MIN(K2,OB%K2)
C
      TEMP_OBST(0)%X1 = MAX(X1,OB%X1)
      TEMP_OBST(0)%X2 = MIN(X2,OB%X2)
      TEMP_OBST(0)%Y1 = MAX(Y1,OB%Y1)
      TEMP_OBST(0)%Y2 = MIN(Y2,OB%Y2)
      TEMP_OBST(0)%Z1 = MAX(Z1,OB%Z1)
      TEMP_OBST(0)%Z2 = MIN(Z2,OB%Z2)
C
C     Ignore OBSTs that do not intersect with HOLE or are merely
C     sliced by the hole.
C
      IF (TEMP_OBST(0)%I2-TEMP_OBST(0)%I1.LT.0 .OR.
     .    TEMP_OBST(0)%J2-TEMP_OBST(0)%J1.LT.0 .OR.
     .    TEMP_OBST(0)%K2-TEMP_OBST(0)%K1.LT.0) CYCLE OBST_LOOP
C
      IF (TEMP_OBST(0)%I2-TEMP_OBST(0)%I1.EQ.0) THEN
         IF (OB%I1.LT.TEMP_OBST(0)%I1 .OR. 
     .       OB%I2.GT.TEMP_OBST(0)%I2) CYCLE OBST_LOOP
      ENDIF
      IF (TEMP_OBST(0)%J2-TEMP_OBST(0)%J1.EQ.0) THEN
         IF (OB%J1.LT.TEMP_OBST(0)%J1 .OR. 
     .       OB%J2.GT.TEMP_OBST(0)%J2) CYCLE OBST_LOOP
      ENDIF
      IF (TEMP_OBST(0)%K2-TEMP_OBST(0)%K1.EQ.0) THEN
         IF (OB%K1.LT.TEMP_OBST(0)%K1 .OR. 
     .       OB%K2.GT.TEMP_OBST(0)%K2) CYCLE OBST_LOOP
      ENDIF
C
      IF (TEMP_OBST(0)%X2.LE.X1 .OR. TEMP_OBST(0)%X1.GE.X2 .OR.
     .    TEMP_OBST(0)%Y2.LE.Y1 .OR. TEMP_OBST(0)%Y1.GE.Y2 .OR.
     .    TEMP_OBST(0)%Z2.LE.Z1 .OR. TEMP_OBST(0)%Z1.GE.Z2)     
     . CYCLE OBST_LOOP
C
C     Start counting new OBSTs that need to be created
C
      NDO=0
C
      IF (OB%I1.LT.I1 .AND. I1.LT.OB%I2) THEN
      NDO=NDO+1
      TEMP_OBST(NDO)=OBSTRUCTION(NN)
      TEMP_OBST(NDO)%I1 = OB%I1
      TEMP_OBST(NDO)%I2 = I1
      TEMP_OBST(NDO)%X1 = OB%X1
      TEMP_OBST(NDO)%X2 = X1
      ENDIF
C
      IF (OB%I1.LT.I2 .AND. I2.LT.OB%I2) THEN
      NDO=NDO+1
      TEMP_OBST(NDO)=OBSTRUCTION(NN)
      TEMP_OBST(NDO)%I1 = I2
      TEMP_OBST(NDO)%I2 = OB%I2
      TEMP_OBST(NDO)%X1 = X2 
      TEMP_OBST(NDO)%X2 = OB%X2
      ENDIF
C
      IF (OB%J1.LT.J1 .AND. J1.LT.OB%J2) THEN
      NDO=NDO+1
      TEMP_OBST(NDO)=OBSTRUCTION(NN)
      TEMP_OBST(NDO)%I1 = MAX(I1,OB%I1)
      TEMP_OBST(NDO)%I2 = MIN(I2,OB%I2)
      TEMP_OBST(NDO)%X1 = MAX(X1,OB%X1)
      TEMP_OBST(NDO)%X2 = MIN(X2,OB%X2)
      TEMP_OBST(NDO)%J1 = OB%J1
      TEMP_OBST(NDO)%J2 = J1
      TEMP_OBST(NDO)%Y1 = OB%Y1
      TEMP_OBST(NDO)%Y2 = Y1
      ENDIF
C
      IF (OB%J1.LT.J2 .AND. J2.LT.OB%J2) THEN
      NDO=NDO+1
      TEMP_OBST(NDO)=OBSTRUCTION(NN)
      TEMP_OBST(NDO)%I1 = MAX(I1,OB%I1)
      TEMP_OBST(NDO)%I2 = MIN(I2,OB%I2)
      TEMP_OBST(NDO)%X1 = MAX(X1,OB%X1)
      TEMP_OBST(NDO)%X2 = MIN(X2,OB%X2)
      TEMP_OBST(NDO)%J1 = J2    
      TEMP_OBST(NDO)%J2 = OB%J2
      TEMP_OBST(NDO)%Y1 = Y2
      TEMP_OBST(NDO)%Y2 = OB%Y2
      ENDIF
C
      IF (OB%K1.LT.K1 .AND. K1.LT.OB%K2) THEN
      NDO=NDO+1
      TEMP_OBST(NDO)=OBSTRUCTION(NN)
      TEMP_OBST(NDO)%I1 = MAX(I1,OB%I1)
      TEMP_OBST(NDO)%I2 = MIN(I2,OB%I2)
      TEMP_OBST(NDO)%X1 = MAX(X1,OB%X1)
      TEMP_OBST(NDO)%X2 = MIN(X2,OB%X2)
      TEMP_OBST(NDO)%J1 = MAX(J1,OB%J1)
      TEMP_OBST(NDO)%J2 = MIN(J2,OB%J2)
      TEMP_OBST(NDO)%Y1 = MAX(Y1,OB%Y1)
      TEMP_OBST(NDO)%Y2 = MIN(Y2,OB%Y2)
      TEMP_OBST(NDO)%K1 = OB%K1
      TEMP_OBST(NDO)%K2 = K1
      TEMP_OBST(NDO)%Z1 = OB%Z1
      TEMP_OBST(NDO)%Z2 = Z1
      ENDIF
C
      IF (OB%K1.LT.K2 .AND. K2.LT.OB%K2) THEN
      NDO=NDO+1
      TEMP_OBST(NDO)=OBSTRUCTION(NN)
      TEMP_OBST(NDO)%I1 = MAX(I1,OB%I1)
      TEMP_OBST(NDO)%I2 = MIN(I2,OB%I2)
      TEMP_OBST(NDO)%X1 = MAX(X1,OB%X1)
      TEMP_OBST(NDO)%X2 = MIN(X2,OB%X2)
      TEMP_OBST(NDO)%J1 = MAX(J1,OB%J1)
      TEMP_OBST(NDO)%J2 = MIN(J2,OB%J2)
      TEMP_OBST(NDO)%Y1 = MAX(Y1,OB%Y1)
      TEMP_OBST(NDO)%Y2 = MIN(Y2,OB%Y2)
      TEMP_OBST(NDO)%K1 = K2
      TEMP_OBST(NDO)%K2 = OB%K2
      TEMP_OBST(NDO)%Z1 = Z2
      TEMP_OBST(NDO)%Z2 = OB%Z2
      ENDIF
C
C     Maintain ordinal rank of original obstruction, but negate it.
C     This will be a code for Smokeview.
C
      TEMP_OBST(:)%ORDINAL = -OB%ORDINAL
C
C     Re-allocate space of new OBSTs, or remove entry for dead OBST
C
      NEW_OBST_IF: IF (NDO.GT.0) THEN
      CALL RE_ALLOCATE_OBST(NM,NB,NDO)
      OBSTRUCTION=>M%OBSTRUCTION
      OBSTRUCTION(NB+1:NB+NDO) = TEMP_OBST(1:NDO)
      NB = NB + NDO
      ENDIF NEW_OBST_IF
C
      IF (T_CREATE.GT.0. .OR. T_REMOVE.LT.100000. .OR.
     .    HEAT_REMOVE.NE.'null' .OR. HEAT_CREATE.NE.'null') THEN
      OBSTRUCTION(NN) = TEMP_OBST(0)
      OBSTRUCTION(NN)%T_REMOVE = T_CREATE
      OBSTRUCTION(NN)%T_CREATE = T_REMOVE
      OBSTRUCTION(NN)%HEAT_REMOVE = HEAT_CREATE
      OBSTRUCTION(NN)%HEAT_CREATE = HEAT_REMOVE
      IF (HEAT_CREATE.NE.'null') OBSTRUCTION(NN)%T_REMOVE = 1000000.
      IF (HEAT_REMOVE.NE.'null') OBSTRUCTION(NN)%T_CREATE = 1000000.
      IF (RGB(1).GT.-1.) THEN 
         OBSTRUCTION(NN)%RGB(1:3) = RGB(1:3)
         OBSTRUCTION(NN)%RGB(4)   = 1.0     
         OBSTRUCTION(NN)%BCI      = -3
         ENDIF
      IF (RGB4(1).GT.-1.) THEN 
         OBSTRUCTION(NN)%RGB(1:4) = RGB4(1:4)
         OBSTRUCTION(NN)%BCI      = -3
         ENDIF
      ELSE
      OBSTRUCTION(NN) = OBSTRUCTION(NB)
      NB = NB-1
      NN = NN-1
      ENDIF
C
      ENDDO OBST_LOOP
C
      ENDDO MESH_LOOP
C
      ENDDO READ_HOLE_LOOP
C
      REWIND(LU5)
C
C
      CONTAINS
C
      SUBROUTINE RE_ALLOCATE_OBST(NM,NB,NDO)
C
      TYPE (OBSTRUCTION_TYPE), ALLOCATABLE, DIMENSION(:) :: DUMMY
      INTEGER IZERO
      INTEGER, INTENT(IN) :: NM,NDO,NB
      TYPE (MESH_TYPE), POINTER :: M
C
      M=>MESH(NM)
      ALLOCATE(DUMMY(0:NB),STAT=IZERO)
      DUMMY(0:NB) = M%OBSTRUCTION(0:NB)
C
      DEALLOCATE(M%OBSTRUCTION)
      ALLOCATE(M%OBSTRUCTION(0:NB+NDO),STAT=IZERO)
      M%OBSTRUCTION(0:NB) = DUMMY(0:NB)
C
      DEALLOCATE(DUMMY)
C
      END SUBROUTINE RE_ALLOCATE_OBST
C
      END SUBROUTINE READ_HOLE
C
C
      SUBROUTINE READ_VENT
C
      INTEGER :: NM
C
      REAL(EB) :: T_ACTIVATE,T_DEACTIVATE,T_OPEN,T_CLOSE,SPREAD_RATE,
     .            TMP_OUTSIDE
      CHARACTER(26) :: HEAT_ACTIVATE,HEAT_DEACTIVATE,
     .                 HEAT_OPEN,HEAT_CLOSE
      CHARACTER(60) :: MESH_ID
      EQUIVALENCE(T_ACTIVATE  ,T_OPEN)
      EQUIVALENCE(T_DEACTIVATE,T_CLOSE)
      EQUIVALENCE(HEAT_ACTIVATE  ,HEAT_OPEN)
      EQUIVALENCE(HEAT_DEACTIVATE,HEAT_CLOSE)
      CHARACTER(11) :: VENT_COLOR,COLOR
      EQUIVALENCE(VENT_COLOR,COLOR)
      LOGICAL REJECT_VENT,EVACUATION
      NAMELIST /VENT/ XB,IOR,CB,PBX,PBY,PBZ,SURF_ID,FYI,RGB,RGB4,
     .                VENT_COLOR,T_OPEN,T_CLOSE,COLOR,
     .                T_ACTIVATE,T_DEACTIVATE,TEXTURE_ORIGIN,
     .                OUTLINE,HEAT_ACTIVATE,HEAT_DEACTIVATE,
     .                HEAT_OPEN,HEAT_CLOSE,XYZ,EVACUATION,MESH_ID,
     .                SPREAD_RATE,TMP_OUTSIDE
C
      MESH_LOOP: DO NM=1,NMESHES
      M=>MESH(NM)
      CALL UNPACK_VAR(NM)
C
      NV = 0
      VLOOP: DO
      SURF_ID = 'null'
      CALL CHECKREAD('VENT',LU5,IOS) ; IF (IOS.EQ.1) EXIT VLOOP
      READ(LU5,NML=VENT,END=3,ERR=4,IOSTAT=IOS)
      NV = NV + 1
      IF (SURF_ID.NE.'null' .AND. SURF_ID.NE.'OPEN' .AND.
     .    SURF_ID.NE.'MIRROR') THEN
         DO N=0,NBT
         IF (SURF_ID.EQ.SURFNAME(N)) CYCLE VLOOP
         ENDDO
         NBT = NBT + 1
         SURFNAME(NBT) = SURF_ID
         ENDIF
    4 IF (IOS.GT.0) THEN
         WRITE(MESSAGE,'(A,I4)') 'ERROR: Problem with VENT number',NV+1
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO VLOOP
    3 REWIND(LU5)
C
      IF (TWO_D)                        NV = NV + 2
      IF (CYLINDRICAL .AND. M%XS.EQ.0.) NV = NV + 1
      IF (EVACUATION_ONLY(NM))          NV = NV + 2
C
      ALLOCATE(M%VENTS(NV),STAT=IZERO)
      CALL ChkMemErr('READ','VENTS',IZERO)
      VENTS=>M%VENTS
C
      VENTS(1:NV)%IBC          = 0
      VENTS(1:NV)%VTI          = 0
      VENTS(1:NV)%HEAT_INDEX_ACTIVATE   = 0
      VENTS(1:NV)%HEAT_INDEX_DEACTIVATE = 0
      VENTS(1:NV)%CLOSED       = .FALSE.
      VENTS(1:NV)%T_OPEN       = 1000000.
      VENTS(1:NV)%T_CLOSE      = 1000000.
      VENTS(1:NV)%X0           = -999.
      VENTS(1:NV)%Y0           = -999.
      VENTS(1:NV)%Z0           = -999.
      VENTS(1:NV)%FIRE_SPREAD_RATE = 0.05
C
      NVO   = NV
      N     = 0
C
      VENTLOOP: DO NN=1,NVO
C
      N       = N + 1
      IOR     = 0
      CB      = 'null'
      PBX     = -1.E6
      PBY     = -1.E6
      PBZ     = -1.E6
      SURF_ID = 'null'
      VENT_COLOR = 'null'
      MESH_ID = 'null'
      RGB        = -1.
      RGB4       = -1.
      XYZ     = -999.
      SPREAD_RATE = 0.05
      T_OPEN  = 1000000.
      T_CLOSE = 1000000.
      REJECT_VENT  = .FALSE.
      TEXTURE_ORIGIN = -999.
      OUTLINE      = .FALSE.
      HEAT_ACTIVATE   = 'null'
      HEAT_DEACTIVATE = 'null'
      TMP_OUTSIDE     = TMPO - TMPM
      IF (     EVACUATION_ONLY(NM)) EVACUATION = .TRUE.
      IF (.NOT.EVACUATION_ONLY(NM)) EVACUATION = .FALSE.
C
      IF (NN.EQ.NVO-2 .AND. CYLINDRICAL .AND. XS.EQ.0.) CB='XBAR0'
      IF (NN.EQ.NVO-1 .AND. JBAR.EQ.1)                  CB='YBAR0'
      IF (NN.EQ.NVO   .AND. JBAR.EQ.1)                  CB='YBAR'
      IF (NN.EQ.NVO-1 .AND. EVACUATION_ONLY(NM))        CB='ZBAR0'
      IF (NN.EQ.NVO   .AND. EVACUATION_ONLY(NM))        CB='ZBAR'
C
      IF (CB.EQ.'null') THEN
      CALL CHECKREAD('VENT',LU5,IOS) ; IF (IOS.EQ.1) EXIT VENTLOOP
      READ(LU5,VENT,END=37,ERR=38)    ! Read in info for VENT N
      ELSE
      SURF_ID = 'MIRROR'
      ENDIF
C
      IF (PBX.GT.-1.E5 .OR. PBY.GT.-1.E5 .OR. PBZ.GT.-1.E5) THEN
      XB(1) = XS
      XB(2) = XF
      XB(3) = YS
      XB(4) = YF
      XB(5) = ZS
      XB(6) = ZF
      IF (PBX.GT.-1.E5) XB(1:2) = PBX
      IF (PBY.GT.-1.E5) XB(3:4) = PBY
      IF (PBZ.GT.-1.E5) XB(5:6) = PBZ
      ENDIF
C
      IF (CB.NE.'null') THEN
      XB(1) = XS
      XB(2) = XF
      XB(3) = YS
      XB(4) = YF
      XB(5) = ZS
      XB(6) = ZF
      IF (CB.EQ.'XBAR0' .OR. CB.EQ.'RBAR0') XB(2) = XS
      IF (CB.EQ.'XBAR'  .OR. CB.EQ.'RBAR')  XB(1) = XF
      IF (CB.EQ.'YBAR0') XB(4) = YS
      IF (CB.EQ.'YBAR' ) XB(3) = YF
      IF (CB.EQ.'ZBAR0') XB(6) = ZS
      IF (CB.EQ.'ZBAR' ) XB(5) = ZF
      ENDIF
C
C Check that the vent is properly specified
C
      IF (MESH_ID.NE.'null' .AND. MESH_ID.NE.MESH_NAME(NM)) 
     .   REJECT_VENT = .TRUE.
C
      IF (XB(1).NE.XB(2) .AND.
     .    XB(3).NE.XB(4) .AND.
     .    XB(5).NE.XB(6)) THEN
          WRITE(MESSAGE,'(A,I4,A)') 'ERROR: VENT',NN,' must be a plane'
          CALL SHUTDOWN(MESSAGE)
          ENDIF
C
      DO I=1,5,2
      IF (XB(I).GT.XB(I+1)) THEN
         DUMMY   = XB(I)
         XB(I)   = XB(I+1)
         XB(I+1) = DUMMY
         ENDIF
      ENDDO
C
      XB(1) = MAX(XB(1),XS)
      XB(2) = MIN(XB(2),XF)
      XB(3) = MAX(XB(3),YS)
      XB(4) = MIN(XB(4),YF)
      XB(5) = MAX(XB(5),ZS)
      XB(6) = MIN(XB(6),ZF)
C
      IF (XB(1).GT.XF .OR. XB(2).LT.XS .OR.
     .    XB(3).GT.YF .OR. XB(4).LT.YS .OR.
     .    XB(5).GT.ZF .OR. XB(6).LT.ZS) REJECT_VENT = .TRUE.
C
      VT=>VENTS(N)
C
      VT%I1 = NINT( GINV(XB(1)-XS,1,NM)*RDXI   ) 
      VT%I2 = NINT( GINV(XB(2)-XS,1,NM)*RDXI   )
      VT%J1 = NINT( GINV(XB(3)-YS,2,NM)*RDETA  ) 
      VT%J2 = NINT( GINV(XB(4)-YS,2,NM)*RDETA  )
      VT%K1 = NINT( GINV(XB(5)-ZS,3,NM)*RDZETA )
      VT%K2 = NINT( GINV(XB(6)-ZS,3,NM)*RDZETA )
C
      IF (XB(1).EQ.XB(2)) THEN
         IF (VT%J1.EQ.VT%J2 .OR. VT%K1.EQ.VT%K2) REJECT_VENT=.TRUE.
         ENDIF
      IF (XB(3).EQ.XB(4)) THEN
         IF (VT%I1.EQ.VT%I2 .OR. VT%K1.EQ.VT%K2) REJECT_VENT=.TRUE.
         ENDIF
      IF (XB(5).EQ.XB(6)) THEN
         IF (VT%I1.EQ.VT%I2 .OR. VT%J1.EQ.VT%J2) REJECT_VENT=.TRUE.
         ENDIF
C
C Evacuation criteria
C
      IF (.NOT.EVACUATION .AND. EVACUATION_ONLY(NM)) REJECT_VENT=.TRUE.
      IF (EVACUATION .AND. .NOT.EVACUATION_ONLY(NM)) REJECT_VENT=.TRUE.
C
C If the VENT is to rejected
C
      IF (REJECT_VENT) THEN
         N = N-1
         NV= NV-1
         CYCLE VENTLOOP
         ENDIF
C
C Vent area
C
      VT%X1 = XB(1)
      VT%X2 = XB(2)
      VT%Y1 = XB(3)
      VT%Y2 = XB(4)
      VT%Z1 = XB(5)
      VT%Z2 = XB(6)
C
      IF (XB(1).EQ.XB(2)) VT%AREA_0 = (XB(4)-XB(3))*(XB(6)-XB(5))
      IF (XB(3).EQ.XB(4)) VT%AREA_0 = (XB(2)-XB(1))*(XB(6)-XB(5))
      IF (XB(5).EQ.XB(6)) VT%AREA_0 = (XB(2)-XB(1))*(XB(4)-XB(3))
C
      DO NNN=0,NBT
      IF (SURF_ID.EQ.SURFNAME(NNN)) VT%IBC = NNN
      ENDDO
      IF (SURF_ID.EQ.'OPEN'  ) VT%IBC = NBT + 1
      IF (SURF_ID.EQ.'MIRROR') VT%IBC = NBT + 2
      IF (SURF_ID.EQ.'OPEN')   VT%VTI =  2
      IF (SURF_ID.EQ.'MIRROR') VT%VTI = -2
      IF ((CB.NE.'null' .OR. 
     .     PBX.GT.-1.E5 .OR. PBY.GT.-1.E5 .OR. PBZ.GT.-1.E5) .AND.
     .    SURF_ID.EQ.'OPEN') VT%VTI = -2
C
      IF (VT%IBC.LE.NBT  ) VT%INDEX = 1
      IF (VT%IBC.EQ.NBT+1) VT%INDEX = 2
      IF (VT%IBC.EQ.NBT+2) VT%INDEX = 3
      VT%IOR = IOR
C
      VT%ORDINAL = NN
C
C Open and Close Logic
C
      VT%HEAT_ACTIVATE   = HEAT_ACTIVATE
      VT%HEAT_DEACTIVATE = HEAT_DEACTIVATE
C
      IF (HEAT_ACTIVATE.NE.'null'.OR.T_OPEN.LT.100000.) VT%CLOSED=.TRUE.
C
      IF (T_CLOSE.LT.0.) THEN
         VT%HEAT_INDEX_DEACTIVATE = NINT(ABS(T_CLOSE))
      ELSE
         VT%T_CLOSE    = T_CLOSE
      ENDIF
C
      IF (T_OPEN.LT.0.) THEN
         VT%HEAT_INDEX_ACTIVATE = NINT(ABS(T_OPEN))
      ELSE
         VT%T_OPEN     = T_OPEN
      ENDIF
C
      SELECT CASE(VENT_COLOR)
      CASE('WHITE')    ; VT%VCI =  0
      CASE('YELLOW')   ; VT%VCI =  1
      CASE('BLUE')     ; VT%VCI =  2
      CASE('RED')      ; VT%VCI =  3
      CASE('GREEN')    ; VT%VCI =  4
      CASE('MAGENTA')  ; VT%VCI =  5
      CASE('CYAN')     ; VT%VCI =  6
      CASE('BLACK')    ; VT%VCI =  7
      CASE('INVISIBLE'); VT%VCI =  8
      CASE DEFAULT     ; VT%VCI = 99
      END SELECT
C
      IF (VT%VCI.EQ.8) VT%VTI = -2
      IF (OUTLINE)     VT%VTI = 2
C
      IF (RGB4(1).LT.-0.5) THEN
      VT%RGB(1) = RGB(1)
      VT%RGB(2) = RGB(2)
      VT%RGB(3) = RGB(3)
      VT%RGB(4) = 1.0    
      IF (RGB(1).GE.0.) VT%VCI    = 99
      ELSE
      VT%RGB(1:4) = RGB4(1:4)
      VT%VCI      = 99
      ENDIF
C
      VT%X0 = XYZ(1)
      VT%Y0 = XYZ(2)
      VT%Z0 = XYZ(3)
      VT%FIRE_SPREAD_RATE = SPREAD_RATE
      VT%TMP_OUTSIDE      = TMP_OUTSIDE + TMPM
      TMPMIN = MIN(TMPMIN,VT%TMP_OUTSIDE)
C
      VT%TEXTURE(:) = TEXTURE_ORIGIN(:)
C
   38 ENDDO VENTLOOP
   37 REWIND(LU5)
C
C Check vents and assign orientations
C
      VENTLOOP2: DO N=1,NV
C
      VT => VENTS(N)
C
      I1 = VT%I1
      I2 = VT%I2
      J1 = VT%J1
      J2 = VT%J2
      K1 = VT%K1
      K2 = VT%K2
C
      IF (VT%IOR.EQ.0) THEN
      IF (I1.EQ.      0 .AND. I2.EQ.0) VT%IOR =  1
      IF (I1.EQ.IBAR .AND. I2.EQ.IBAR) VT%IOR = -1
      IF (J1.EQ.      0 .AND. J2.EQ.0) VT%IOR =  2
      IF (J1.EQ.JBAR .AND. J2.EQ.JBAR) VT%IOR = -2
      IF (K1.EQ.      0 .AND. K2.EQ.0) VT%IOR =  3
      IF (K1.EQ.KBAR .AND. K2.EQ.KBAR) VT%IOR = -3
      ENDIF
C
      ORIENTATION_IF: IF (VT%IOR.EQ.0) THEN
C
      IF (I1.EQ.I2) THEN
      DO K=K1+1,K2
      DO J=J1+1,J2
      IF (.NOT.SOLID(ICA(I2+1,J,K))) VT%IOR =  1
      IF (.NOT.SOLID(ICA(I2  ,J,K))) VT%IOR = -1
      ENDDO
      ENDDO
      ENDIF
C
      IF (J1.EQ.J2) THEN
      DO K=K1+1,K2
      DO I=I1+1,I2
      IF (.NOT.SOLID(ICA(I,J2+1,K))) VT%IOR =  2
      IF (.NOT.SOLID(ICA(I,J2  ,K))) VT%IOR = -2
      ENDDO
      ENDDO
      ENDIF
C
      IF (K1.EQ.K2) THEN
      DO J=J1+1,J2
      DO I=I1+1,I2
      IF (.NOT.SOLID(ICA(I,J,K2+1))) VT%IOR =  3
      IF (.NOT.SOLID(ICA(I,J,K2  ))) VT%IOR = -3
      ENDDO
      ENDDO
      ENDIF
C
      ENDIF ORIENTATION_IF
C
      IF (VT%IOR.EQ.0) THEN
         WRITE(MESSAGE,'(A,I3,A,I3)') 
     .   'ERROR: Specify orientation of VENT ',VT%ORDINAL,
     .   ', MESH NUMBER',NM
         CALL SHUTDOWN(MESSAGE)
         ENDIF
C
C Specify that an internal vent is to be of TYPE 4 or 5
C
      SELECT CASE(ABS(VT%IOR))
      CASE(1)
      IF (I1.GE.1 .AND. I1.LE.IBM1) THEN
         IF (VT%INDEX.EQ.2 .OR. VT%INDEX.EQ.3) THEN
            WRITE(MESSAGE,'(A,I3,A)') 
     .      'ERROR: OPEN or MIRROR VENT ',N,
     .      ' must be on an exterior boundary.'
            CALL SHUTDOWN(MESSAGE)
            ENDIF
         VT%INDEX = 5
         IF (.NOT.SOLID(ICA(I2+1,J2,K2)) .AND. 
     .       .NOT.SOLID(ICA(I2,J2,K2))) VT%INDEX = 4 
         ENDIF
      CASE(2)
      IF (J1.GE.1 .AND. J1.LE.JBM1) THEN
         IF (VT%INDEX.EQ.2 .OR. VT%INDEX.EQ.3) THEN
            WRITE(MESSAGE,'(A,I3,A)') 
     .      'ERROR: OPEN or MIRROR VENT ',N,
     .      ' must be on an exterior boundary.'
            CALL SHUTDOWN(MESSAGE)
            ENDIF
         VT%INDEX = 5
         IF (.NOT.SOLID(ICA(I2,J2+1,K2)) .AND. 
     .       .NOT.SOLID(ICA(I2,J2,K2))) VT%INDEX = 4 
         ENDIF
      CASE(3)
      IF (K1.GE.1 .AND. K1.LE.KBM1) THEN
         IF (VT%INDEX.EQ.2 .OR. VT%INDEX.EQ.3) THEN
            WRITE(MESSAGE,'(A,I3,A)') 
     .      'ERROR: OPEN or MIRROR VENT ',N,
     .      ' must be on an exterior boundary.'
            CALL SHUTDOWN(MESSAGE)
            ENDIF
         VT%INDEX = 5
         IF (.NOT.SOLID(ICA(I2,J2,K2+1)) .AND. 
     .       .NOT.SOLID(ICA(I2,J2,K2))) VT%INDEX = 4 
         ENDIF
      END SELECT
C
C Open up boundary cells if it is an open vent
C
      IF ( VT%INDEX.EQ.2 .AND. .NOT.VT%CLOSED) THEN
C
      SELECT CASE(VT%IOR)
      CASE( 1) ; CALL BLKCLL(NM,   0,   0,J1+1,  J2,K1+1,  K2,0)
      CASE(-1) ; CALL BLKCLL(NM,IBP1,IBP1,J1+1,  J2,K1+1,  K2,0)
      CASE( 2) ; CALL BLKCLL(NM,I1+1,  I2,   0,   0,K1+1,  K2,0)
      CASE(-2) ; CALL BLKCLL(NM,I1+1,  I2,JBP1,JBP1,K1+1,  K2,0)
      CASE( 3) ; CALL BLKCLL(NM,I1+1,  I2,J1+1,  J2,   0,   0,0)
      CASE(-3) ; CALL BLKCLL(NM,I1+1,  I2,J1+1,  J2,KBP1,KBP1,0)
      END SELECT
C
      ENDIF
C
      ENDDO VENTLOOP2
C
C Compute the volume of the domain
C
      VOL = 0
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      IF (.NOT.SOLID(ICA(I,J,K))) VOL = VOL + DX(I)*RC(I)*DY(J)*DZ(K)
      ENDDO
      ENDDO
      ENDDO
C
C Compute vent areas and check for passive openings
C
      DO N=1,NV
C
      VT => VENTS(N)
C
      VT%AREA = 0.
      I1 = VT%I1
      I2 = VT%I2
      J1 = VT%J1
      J2 = VT%J2
      K1 = VT%K1
      K2 = VT%K2
C
      SELECT CASE(ABS(VT%IOR))
      CASE(1)
      DO K=K1+1,K2
      DO J=J1+1,J2
      VT%AREA = VT%AREA + DY(J)*DZ(K)
      ENDDO
      ENDDO
      CASE(2)
      DO K=K1+1,K2
      DO I=I1+1,I2
      VT%AREA = VT%AREA + DX(I)*DZ(K)
      ENDDO
      ENDDO
      CASE(3)
      DO J=J1+1,J2
      DO I=I1+1,I2
      VT%AREA = VT%AREA + DX(I)*DY(J)
      ENDDO
      ENDDO
      END SELECT
C
      ENDDO   
C
      ENDDO MESH_LOOP
C
      END SUBROUTINE READ_VENT
C
C
      SUBROUTINE READ_INIT
C
      REAL(EB) VALUE
      NAMELIST /INIT/ XB,QUANTITY,VALUE
C
      NIB = 0
      COUNT_LOOP: DO
      CALL CHECKREAD('INIT',LU5,IOS) ; IF (IOS.EQ.1) EXIT COUNT_LOOP
      READ(LU5,NML=INIT,END=11,ERR=12,IOSTAT=IOS)
      NIB = NIB + 1
   12 IF (IOS.GT.0) THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with INIT no.',NIB+1
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO COUNT_LOOP
   11 REWIND(LU5)
C
      IF (NIB.GT.0) THEN
C
      ALLOCATE(XB1(NIB),STAT=IZERO)
      CALL ChkMemErr('READ','XB1',IZERO)
      ALLOCATE(XB2(NIB),STAT=IZERO)
      CALL ChkMemErr('READ','XB2',IZERO)
      ALLOCATE(YB1(NIB),STAT=IZERO)
      CALL ChkMemErr('READ','YB1',IZERO)
      ALLOCATE(YB2(NIB),STAT=IZERO)
      CALL ChkMemErr('READ','YB2',IZERO)
      ALLOCATE(ZB1(NIB),STAT=IZERO)
      CALL ChkMemErr('READ','ZB1',IZERO)
      ALLOCATE(ZB2(NIB),STAT=IZERO)
      CALL ChkMemErr('READ','ZB2',IZERO)
      ALLOCATE(INIT_VALUE(NIB),STAT=IZERO)
      CALL ChkMemErr('READ','INIT_VALUE',IZERO) 
      ALLOCATE(INIT_INDEX(NIB),STAT=IZERO)
      CALL ChkMemErr('READ','INIT_INDEX',IZERO) 
C
      INIT_LOOP: DO N=1,NIB
C
      VALUE = 0.
      QUANTITY = 'TEMPERATURE'
C
      CALL CHECKREAD('INIT',LU5,IOS) ; IF (IOS.EQ.1) EXIT INIT_LOOP
      READ(LU5,INIT,END=41,ERR=42) 
C
      XB1(N) = XB(1)
      XB2(N) = XB(2)
      YB1(N) = XB(3)
      YB2(N) = XB(4)
      ZB1(N) = XB(5)
      ZB2(N) = XB(6)
C
      INIT_VALUE(N) = VALUE
      IF (QUANTITY.EQ.'TEMPERATURE') THEN
         INIT_VALUE(N) = VALUE + TMPM
         INIT_INDEX(N) = -1
         TMPMIN = MIN(TMPMIN,INIT_VALUE(N))
         ENDIF
      IF (QUANTITY.EQ.'DENSITY') THEN
         INIT_INDEX(N) = -2
         RHOMAX = MAX(RHOMAX,INIT_VALUE(N))
         ENDIF
      IF (QUANTITY.EQ.SPECIES_ID(1)) INIT_INDEX(N) =  1
      IF (QUANTITY.EQ.SPECIES_ID(2)) INIT_INDEX(N) =  2
      IF (QUANTITY.EQ.SPECIES_ID(3)) INIT_INDEX(N) =  3
      IF (QUANTITY.EQ.SPECIES_ID(4)) INIT_INDEX(N) =  4
      IF (QUANTITY.EQ.SPECIES_ID(5)) INIT_INDEX(N) =  5
C
   42 ENDDO INIT_LOOP
   41 REWIND(LU5)
C
      ENDIF
C
      END SUBROUTINE READ_INIT
C
C
      SUBROUTINE READ_SPRK
C
      REAL(EB) :: T_ACTIVATE,T_DEACTIVATE,THETA_AVG,VLEN,
     .            XI,YJ,ZK
      REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: WATER_MASS
      INTEGER NM,NDTHETA,NDPHI,IPC,NNH,NNL,IMIN,IMAX
      INTEGER VELOCITY_INDEX,FLUX_INDEX,DISTRIBUTION_INDEX
      REAL(EB) :: MIN_SPRAY_ANGLE,MAX_SPRAY_ANGLE,CONSTANT_VELOCITY,
     .            SUM,SUM2,DELAY
      TYPE (SPRINKLER_TYPE), POINTER :: S
      CHARACTER(60) :: FN95
      CHARACTER(30) :: HEAT_ACTIVATE
C
      NAMELIST /SPRK/ XYZ,MAKE,T_ACTIVATE,T_DEACTIVATE,FYI,
     .                ORIENTATION,ROTATION,PART_ID,LABEL,DELAY,
     .                HEAT_ACTIVATE
      NAMELIST /PIPE/ DELAY,PRESSURE,FYI
C
C Read in pipe pressure. Number of Sprinklers (NSPR) done in READ_PART
C
      PRESSURE = -1. 
      DELAY    = 0.
      T_ACT_FIRST = 1000000.
C
      PIPE_LOOP: DO
      CALL CHECKREAD('PIPE',LU5,IOS) ; IF (IOS.EQ.1) EXIT PIPE_LOOP
      READ(LU5,PIPE,END=47,ERR=48,IOSTAT=IOS)
   48 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with PIPE line')
      ENDDO PIPE_LOOP
   47 REWIND(LU5)
C
      SYSTEM_DELAY = DELAY
C
C If no sprinklers, don't allocate or read sprinkler-specific info
C
      IF (NSPR.EQ.0) GOTO 100
C
      DROPLET_FILE = .TRUE.
C
C Allocate arrays that pertain to sprinklers 
C
      ALLOCATE(SPRINKLER_HEAD(NSPR),STAT=IZERO)
      CALL ChkMemErr('READ','SPRINKLER_HEAD',IZERO)
C
      SPRINKLER_HEAD(:)%T        = 0.
      SPRINKLER_HEAD(:)%DELAY    = 0.
      SPRINKLER_HEAD(:)%ACT_CODE = 0
      SPRINKLER_HEAD(:)%TMP_L    = TMPA
      SPRINKLER_HEAD(:)%TMP_L_S  = TMPA
      SPRINKLER_HEAD(:)%DROPLET_CLASS = NPC
      SPRINKLER_HEAD(:)%HEAT_INDEX = 0
      SPRINKLER_HEAD(:)%HEAT_ACTIVATE = 'null'
C
C Read in sprinkler locations and parameters line by line. Assign
C default parameters unless otherwise indicated.
C
      DT_SPRK = TWFIN/REAL(NFRAMES)
      MAKE    = 'null'
      NSPRO   = NSPR
      N       = 0
      NST     = 0
      ALLOCATE(SPRK_MAKE(100),STAT=IZERO)
      CALL ChkMemErr('READ','SPRK_MAKE',IZERO)
C
      SPRKLOOP: DO NN=1,NSPRO
C
      N = N+1
      T_ACTIVATE   = 1.E12
      T_DEACTIVATE = 1.E12
      PART_ID      = 'null'
      LABEL        = 'null'
      HEAT_ACTIVATE= 'null'
      DELAY        = 0.
      ORIENTATION(1:2) = 0.
      ORIENTATION(3)   =-1.
      ROTATION         = 0.
C
      CALL CHECKREAD('SPRK',LU5,IOS) ; IF (IOS.EQ.1) EXIT SPRKLOOP
      READ(LU5,SPRK) 
C
      MESH_LOOP: DO NM=1,NMESHES
      IF (XYZ(1).GT.MESH(NM)%XS .AND. XYZ(1).LT.MESH(NM)%XF .AND.
     .    XYZ(2).GT.MESH(NM)%YS .AND. XYZ(2).LT.MESH(NM)%YF .AND.
     .    XYZ(3).GT.MESH(NM)%ZS .AND. XYZ(3).LT.MESH(NM)%ZF) THEN
         SPRINKLER_HEAD(N)%MESH = NM
         EXIT MESH_LOOP
         ENDIF
      IF (NM.EQ.NMESHES) THEN
         N    = N-1
         NSPR = NSPR-1
         CYCLE SPRKLOOP
         ENDIF
      ENDDO MESH_LOOP
C
      NM = SPRINKLER_HEAD(N)%MESH
      M  =>MESH(NM)
C
      IF (MAKE.EQ.'null') 
     .   CALL SHUTDOWN('ERROR: Specify MAKE of sprinkler')
C
      SPRINKLER_HEAD(N)%T_DEACT = T_DEACTIVATE
      SPRINKLER_HEAD(N)%T_ACT   = T_ACTIVATE
      SPRINKLER_HEAD(N)%DELAY   = DELAY
      SPRINKLER_HEAD(N)%HEAT_ACTIVATE = HEAT_ACTIVATE
      SPRINKLER_HEAD(N)%X       = XYZ(1) 
      SPRINKLER_HEAD(N)%Y       = XYZ(2)
      SPRINKLER_HEAD(N)%Z       = XYZ(3)
      XI = M%CELLSI(FLOOR((XYZ(1)-M%XS)*M%RDXINT))
      YJ = M%CELLSJ(FLOOR((XYZ(2)-M%YS)*M%RDYINT))
      ZK = M%CELLSK(FLOOR((XYZ(3)-M%ZS)*M%RDZINT))
      SPRINKLER_HEAD(N)%I = FLOOR(XI+1.)
      SPRINKLER_HEAD(N)%J = FLOOR(YJ+1.)
      SPRINKLER_HEAD(N)%K = FLOOR(ZK+1.)
C
      IF (PART_ID.NE.'null') THEN
         DO IPC=1,NPC
         LP=>LAGRANGIAN(IPC)
         IF (LP%CLASS_NAME.EQ.PART_ID)
     .       SPRINKLER_HEAD(N)%DROPLET_CLASS = IPC
         ENDDO
         ENDIF
C
      IF (LABEL.NE.'null') THEN
         SPRINKLER_HEAD(N)%LABEL = LABEL
      ELSE
         WRITE(SPRINKLER_HEAD(N)%LABEL,'(A,I4.4)') 'SPRK',N
      ENDIF
C
C  Orientation=x,y,z direction vector of sprinkler pipe
C  Rotation = degrees of rotation of sprinkler head
C
      VLEN =
     . (ORIENTATION(1)**2.+ORIENTATION(2)**2.+ORIENTATION(3)**2.)**.5
      ORIENTATION = ORIENTATION/VLEN
      SPRINKLER_HEAD(N)%DIR(1:3) = ORIENTATION(1:3)
      SPRINKLER_HEAD(N)%DIR(4)   = ROTATION*TWOPI/360.
C
      DO I=1,NST
      IF (MAKE.EQ.SPRK_MAKE(I)) THEN
         SPRINKLER_HEAD(N)%INDEX  = I
         GOTO 30 
         ENDIF
      ENDDO
      NST = NST+1
      SPRK_MAKE(NST) = MAKE
      SPRINKLER_HEAD(N)%INDEX = NST
   30 CONTINUE
C
      ENDDO SPRKLOOP
C
C Allocate sprinkler type array
C
  100 CONTINUE
C
      ALLOCATE(SPRINKLER(NST))
C
C Read data from the sprinkler files
C
      READ_SPRINKLER_FILE_LOOP: DO N=1,NST
C
      S => SPRINKLER(N)
C
      FN95 = TRIM(SPRK_MAKE(N))//'.spk'
C
      LU95 = 95
      INQUIRE(FILE=FN95,EXIST=EX)
      IF (EX) THEN
      OPEN(LU95,FILE=FN95,FORM='FORMATTED',ACTION='READ',STATUS='OLD')
      ELSE 
      FN95 = TRIM(DATABASE_DIRECTORY)//TRIM(SPRK_MAKE(N))//'.spk'
      INQUIRE(FILE=FN95,EXIST=EX)
      IF (EX) THEN
      OPEN(LU95,FILE=FN95,FORM='FORMATTED',ACTION='READ',STATUS='OLD')
      ELSE
      WRITE(MESSAGE,'(4A)')
     . 'ERROR: Sprinkler file ',TRIM(SPRK_MAKE(N))//'.spk',
     . ' cannot be found.'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
      ENDIF
C
      S%RTI = 100.
      CALL SEARCH_KEYWORD('RTI',LU95,IOS)
      IF (IOS.EQ.0) READ(LU95,*) S%RTI
C
      S%K_FACTOR = 1.
      CALL SEARCH_KEYWORD('K-FACTOR',LU95,IOS)
      IF (IOS.EQ.0) READ(LU95,*) S%K_FACTOR
C
      S%C_FACTOR = 0.
      CALL SEARCH_KEYWORD('C-FACTOR',LU95,IOS)
      IF (IOS.EQ.0) READ(LU95,*) S%C_FACTOR
C
      S%TMP_ACT = 5000.
      CALL SEARCH_KEYWORD('ACTIVATION_TEMPERATURE',LU95,IOS)
      IF (IOS.EQ.0) READ(LU95,*) S%TMP_ACT
      S%TMP_ACT = S%TMP_ACT + TMPM
C
      S%P_FACTOR = 1.
      S%OPERATING_PRESSURE = 1.
      CALL SEARCH_KEYWORD('OPERATING_PRESSURE',LU95,IOS)
      IF (IOS.EQ.0) READ(LU95,*) S%OPERATING_PRESSURE
C
      S%OFFSET = 0.2
      CALL SEARCH_KEYWORD('OFFSET_DISTANCE',LU95,IOS)
      IF (IOS.EQ.0) READ(LU95,*) S%OFFSET
C
      S%WATER_TEMPERATURE = 20.
      CALL SEARCH_KEYWORD('WATER_TEMPERATURE',LU95,IOS)
      IF (IOS.EQ.0) READ(LU95,*) S%WATER_TEMPERATURE
      S%WATER_TEMPERATURE = S%WATER_TEMPERATURE + TMPM
C
C Allocate some necessary arrays for distributions
C
      ALLOCATE(S%NTHETA(3),STAT=IZERO)
      CALL ChkMemErr('READ','NTHETA',IZERO)
      ALLOCATE(S%NPHI(3),STAT=IZERO)
      CALL ChkMemErr('READ','NPHI',IZERO)
C
C Assign droplet velocities as a function of phi and theta
C
      VELOCITY_INDEX = 1
      MIN_SPRAY_ANGLE = 60.
      MAX_SPRAY_ANGLE = 75.
      CONSTANT_VELOCITY = 10.
      CALL SEARCH_KEYWORD('VELOCITY',LU95,IOS)
      IF (IOS.EQ.0) READ(LU95,*) VELOCITY_INDEX
      SELECT CASE(VELOCITY_INDEX)
      CASE(1)
      IF (IOS.EQ.0)
     .READ(LU95,*) MIN_SPRAY_ANGLE,MAX_SPRAY_ANGLE,CONSTANT_VELOCITY
      MIN_SPRAY_ANGLE = MIN_SPRAY_ANGLE*PI/180.
      MAX_SPRAY_ANGLE = MAX_SPRAY_ANGLE*PI/180.
      S%NTHETA(:) = 180
      S%NPHI(:)   = 2
      NDTHETA = S%NTHETA(1)
      NDPHI   = S%NPHI(1)
      ALLOCATE(S%DROPLET_VELOCITY(NDTHETA,NDPHI),STAT=IZERO)
      CALL ChkMemErr('READ','DROPLET_VELOCITY',IZERO)
      S%DROPLET_VELOCITY = CONSTANT_VELOCITY
      CASE(2)
      READ(LU95,*) S%NTHETA(1),S%NPHI(1)
      NDTHETA = S%NTHETA(1)
      NDPHI   = S%NPHI(1)
      ALLOCATE(S%DROPLET_VELOCITY(NDTHETA,NDPHI),STAT=IZERO)
      CALL ChkMemErr('READ','DROPLET_VELOCITY',IZERO)
      DO I=1,S%NTHETA(1)
      READ(LU95,*) (S%DROPLET_VELOCITY(I,J),J=1,S%NPHI(1))
      ENDDO
      END SELECT
C
C Assign water flux as a function of phi and theta
C
      FLUX_INDEX = 1
      CALL SEARCH_KEYWORD('FLUX',LU95,IOS)
      IF (IOS.EQ.0) READ(LU95,*) FLUX_INDEX
      SELECT CASE(FLUX_INDEX)
      CASE(1)
      NDTHETA = S%NTHETA(1)
      NDPHI   = S%NPHI(1)
      ALLOCATE(S%WATER_FLUX(NDTHETA,NDPHI),STAT=IZERO)
      CALL ChkMemErr('READ','WATER_FLUX',IZERO) 
      S%WATER_FLUX = 0.
      S%DTHETA = PI/REAL(S%NTHETA(1),EB)
      IMIN =   MIN(S%NTHETA(1),FLOOR(MIN_SPRAY_ANGLE/S%DTHETA)+1)
      IMAX = MIN(S%NTHETA(1),CEILING(MAX_SPRAY_ANGLE/S%DTHETA))
      IF (IMIN.GT.IMAX) IMAX = IMIN
      S%WATER_FLUX(IMIN:IMAX,1:S%NPHI(1)) = 1.
      CASE(2)
      READ(LU95,*) S%NTHETA(2),S%NPHI(2)
      NDTHETA = S%NTHETA(2)
      NDPHI   = S%NPHI(2)
      ALLOCATE(S%WATER_FLUX(NDTHETA,NDPHI),STAT=IZERO)
      CALL ChkMemErr('READ','WATER_FLUX',IZERO) 
      S%WATER_FLUX = 0.
      DO I=1,S%NTHETA(2)
      READ(LU95,*) (S%WATER_FLUX(I,J),J=1,S%NPHI(2))
      ENDDO
      END SELECT
C
C Assign droplet size distribution as a function of phi and theta
C
      DISTRIBUTION_INDEX = 1
      CALL SEARCH_KEYWORD('SIZE_DISTRIBUTION',LU95,IOS)
      IF (IOS.EQ.0) READ(LU95,*) DISTRIBUTION_INDEX
      SELECT CASE(DISTRIBUTION_INDEX)
      CASE(1)
      ALLOCATE(S%DROPLET_DIAMETER(1,1),STAT=IZERO)
      CALL ChkMemErr('READ','DROPLET_DIAMETER',IZERO)
      S%NTHETA(3)=1 ; S%NPHI(3)=1
      IF (IOS.EQ.0) READ(LU95,*) S%DROPLET_DIAMETER(1,1),S%GAMMA
      CASE(2)
      READ(LU95,*) S%NTHETA(3),S%NPHI(3),S%GAMMA
      ALLOCATE(S%DROPLET_DIAMETER(S%NTHETA(3),S%NPHI(3)),STAT=IZERO)
      CALL ChkMemErr('READ','DROPLET_DIAMETER',IZERO)
      S%DROPLET_DIAMETER = 0.
      DO I=1,S%NTHETA(3)
      READ(LU95,*) (S%DROPLET_DIAMETER(I,J),J=1,S%NPHI(3))
      ENDDO
      END SELECT
C
      S%DROPLET_DIAMETER = S%DROPLET_DIAMETER*1.E-6
      S%SIGMA = 1.15/S%GAMMA
C
      ALLOCATE(S%CDF_DROP(0:NDC,S%NTHETA(3),S%NPHI(3)),STAT=IZERO)
      CALL ChkMemErr('READ','CDF_DROP',IZERO)
      ALLOCATE(S%RDROP(0:NDC,S%NTHETA(3),S%NPHI(3)),STAT=IZERO)
      CALL ChkMemErr('READ','RDROP',IZERO)
      ALLOCATE(S%IL_DROP(NSTRATA,S%NTHETA(3),S%NPHI(3)),STAT=IZERO)
      CALL ChkMemErr('READ','IL_DROP',IZERO)
      ALLOCATE(S%IU_DROP(NSTRATA,S%NTHETA(3),S%NPHI(3)),STAT=IZERO)
      CALL ChkMemErr('READ','IU_DROP',IZERO)
      ALLOCATE(S%W_DROP(NSTRATA,S%NTHETA(3),S%NPHI(3)),STAT=IZERO)
      CALL ChkMemErr('READ','W_DROP',IZERO)
C
      IF (LU95.GT.0) CLOSE(LU95)
      ENDDO READ_SPRINKLER_FILE_LOOP
C
C Compute sprinkler spray cumulative distributions
C
      SPRAY_DISTRIBUTION_LOOP: DO N=1,NST
C
      S => SPRINKLER(N)
C
      NDTHETA = S%NTHETA(2)
      NDPHI   = S%NPHI(2)
C
      ALLOCATE(S%THETA(0:NDTHETA),STAT=IZERO)
      CALL ChkMemErr('READ','THETA',IZERO)
      ALLOCATE(S%PHI(0:NDPHI),STAT=IZERO)
      CALL ChkMemErr('READ','PHI',IZERO)
C
      N_TE = 1000
      ALLOCATE(S%I_TABLE(0:N_TE),STAT=IZERO)
      CALL ChkMemErr('READ','I_TABLE',IZERO)
      ALLOCATE(S%J_TABLE(0:N_TE),STAT=IZERO)
      CALL ChkMemErr('READ','J_TABLE',IZERO)
      ALLOCATE(WATER_MASS(S%NTHETA(2),S%NPHI(2)),STAT=IZERO)
      CALL ChkMemErr('READ','WATER_MASS',IZERO)
C
      S%DPHI   = TWOPI/REAL(S%NPHI(2),EB)
      S%DTHETA = PI/REAL(S%NTHETA(2),EB)
C
      DO J=0,S%NPHI(2)
      S%PHI(J) = J*S%DPHI
      ENDDO
      DO I=0,S%NTHETA(2)
      S%THETA(I) = I*S%DTHETA
      ENDDO
C
      SUM = 0.
      DO J=1,S%NPHI(2)
      DO I=1,S%NTHETA(2)
      THETA_AVG = 0.5*(S%THETA(I)+S%THETA(I-1))
      WATER_MASS(I,J) = S%WATER_FLUX(I,J)*S%DTHETA*S%DPHI*SIN(THETA_AVG)
      SUM = SUM + WATER_MASS(I,J)
      ENDDO
      ENDDO
C
      SUM2=0.
      DO J=1,S%NPHI(2)
      DO I=1,S%NTHETA(2)
      IF (WATER_MASS(I,J).GT.0.) THEN
      NNL = NINT(SUM2*N_TE)
      SUM2 = SUM2 + WATER_MASS(I,J)/SUM
      NNH = NINT(SUM2*N_TE)
      S%I_TABLE(NNL:NNH) = I
      S%J_TABLE(NNL:NNH) = J
      ENDIF
      ENDDO
      ENDDO
C
      DEALLOCATE(WATER_MASS)
C
      ENDDO SPRAY_DISTRIBUTION_LOOP
C
      REWIND(LU5)
C
      END SUBROUTINE READ_SPRK
C
C
      SUBROUTINE READ_HEAT
C
      REAL(EB) :: RTI,ACTIVATION_TEMPERATURE,XI,YJ,ZK
      INTEGER :: NHDO,NM,NNN
      NAMELIST /HEAT/ XYZ,RTI,ACTIVATION_TEMPERATURE,FYI,LABEL
      TYPE(SPRINKLER_HEAD_TYPE), POINTER :: SH
C
C Read all sprinkler lines and count the number of sprinklers
C
      NHD = 0
      COUNT_HEAT_LOOP: DO
      CALL CHECKREAD('HEAT',LU5,IOS) ; IF(IOS.EQ.1) EXIT COUNT_HEAT_LOOP
      READ(LU5,HEAT,END=7,ERR=8,IOSTAT=IOS)
      NHD = NHD + 1
    8 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with HEAT line')
      ENDDO COUNT_HEAT_LOOP
    7 REWIND(LU5)
C
C If no heat detectors, record in .smv file and get out
C
      IF (NHD.EQ.0) GOTO 100
C
C Allocate arrays for heat detectors 
C
      ALLOCATE(HEAT_DETECTOR(NHD),STAT=IZERO)
      CALL ChkMemErr('READ','HEAT_DETECTOR',IZERO)
C
      HEAT_DETECTOR(:)%TMP_L   = TMPA
      HEAT_DETECTOR(:)%TMP_L_S = TMPA
C
C Read in sprinkler locations and parameters line by line. Assign
C default parameters unless otherwise indicated.
C
      NHDO    = NHD
      N       = 0
      DT_HEAT = TWFIN/REAL(NFRAMES)
C
      HDLOOP: DO NN=1,NHDO
C
      N = N+1
C
      CALL CHECKREAD('HEAT',LU5,IOS) ; IF (IOS.EQ.1) EXIT HDLOOP
      ACTIVATION_TEMPERATURE = 74.
      RTI = 100.
      LABEL = 'null'
C
      READ(LU5,HEAT) 
C
      MESH_LOOP: DO NM=1,NMESHES
C
      IF (XYZ(1).GT.MESH(NM)%XS .AND. XYZ(1).LT.MESH(NM)%XF .AND.
     .    XYZ(2).GT.MESH(NM)%YS .AND. XYZ(2).LT.MESH(NM)%YF .AND.
     .    XYZ(3).GT.MESH(NM)%ZS .AND. XYZ(3).LT.MESH(NM)%ZF) THEN
         HEAT_DETECTOR(N)%MESH = NM
         EXIT MESH_LOOP
         ENDIF
      IF (NM.EQ.NMESHES) THEN
         N    = N-1
         NHD  = NHD-1
         CYCLE HDLOOP
         ENDIF
C
      ENDDO MESH_LOOP
C
C Check for obstructions, vents, sprinklers linked to HEAT detectors
C
      IF (LABEL.NE.'null') THEN
C
      DO NM=1,NMESHES
      M  =>MESH(NM)
      OBST_LOOP: DO NNN=1,M%NB
      OB=>M%OBSTRUCTION(NNN)
      IF (OB%HEAT_REMOVE.EQ.'ALL') OB%HEAT_INDEX_REMOVE = NHD+1
      IF (OB%HEAT_CREATE.EQ.'ALL') OB%HEAT_INDEX_CREATE = NHD+1
      IF (LABEL.EQ.OB%HEAT_REMOVE) OB%HEAT_INDEX_REMOVE = N
      IF (LABEL.EQ.OB%HEAT_CREATE) OB%HEAT_INDEX_CREATE = N
      ENDDO OBST_LOOP
C
      VENT_LOOP: DO NNN=1,M%NV
      VT=>M%VENTS(NNN)
      IF (VT%HEAT_ACTIVATE.EQ.'ALL')   VT%HEAT_INDEX_ACTIVATE = NHD+1
      IF (VT%HEAT_DEACTIVATE.EQ.'ALL') VT%HEAT_INDEX_DEACTIVATE=NHD+1
      IF (LABEL.EQ.VT%HEAT_ACTIVATE)   VT%HEAT_INDEX_ACTIVATE = N
      IF (LABEL.EQ.VT%HEAT_DEACTIVATE) VT%HEAT_INDEX_DEACTIVATE=N
      ENDDO VENT_LOOP
      ENDDO
C
      DO NNN=1,NSPR
      SH=>SPRINKLER_HEAD(NNN)
      IF (LABEL.EQ.SH%HEAT_ACTIVATE) SH%HEAT_INDEX = N
      ENDDO
C
      ENDIF
C
C Apply properties of HEAT detector
C
      NM = HEAT_DETECTOR(N)%MESH
      M  =>MESH(NM)
C
      HEAT_DETECTOR(N)%X = XYZ(1)
      HEAT_DETECTOR(N)%Y = XYZ(2) 
      HEAT_DETECTOR(N)%Z = XYZ(3)
      HEAT_DETECTOR(N)%TMP_ACT = ACTIVATION_TEMPERATURE + TMPM
      HEAT_DETECTOR(N)%RTI     = RTI
      XI = M%CELLSI(FLOOR((XYZ(1)-M%XS)*M%RDXINT))
      YJ = M%CELLSJ(FLOOR((XYZ(2)-M%YS)*M%RDYINT))
      ZK = M%CELLSK(FLOOR((XYZ(3)-M%ZS)*M%RDZINT))
      HEAT_DETECTOR(N)%I = FLOOR(XI+1.)
      HEAT_DETECTOR(N)%J = FLOOR(YJ+1.)
      HEAT_DETECTOR(N)%K = FLOOR(ZK+1.)
C
      IF (LABEL.NE.'null') THEN
         HEAT_DETECTOR(N)%LABEL = LABEL
      ELSE
         WRITE(HEAT_DETECTOR(N)%LABEL,'(A,I4.4)') 'HEAT',N
      ENDIF
C
      ENDDO HDLOOP
C
C Write out heat detector coordinates to .smv file
C
  100 REWIND(LU5)
C
      END SUBROUTINE READ_HEAT
C
C
      SUBROUTINE READ_SMOD
C
      REAL(EB) ACTIVATION_OBSCURATION,
     .         ALPHA_E,ALPHA_C,BETA_E,BETA_C,LENGTH
C
      INTEGER  NSDO,NM,NU_WORK
      NAMELIST /SMOD/XYZ,ACTIVATION_OBSCURATION,FYI,IOR,
     .   ALPHA_E,BETA_E,ALPHA_C,BETA_C,LABEL,LENGTH
      EQUIVALENCE(LENGTH,ALPHA_C)
C
C Read all smoke detector lines and count the number of detectors 
C
      NSD = 0
      COUNT_SMOK_LOOP: DO
      CALL CHECKREAD('SMOD',LU5,IOS) ; IF(IOS.EQ.1)EXIT COUNT_SMOK_LOOP
      READ(LU5,SMOD,END=7,ERR=8,IOSTAT=IOS)
      NSD = NSD + 1
    8 IF (IOS.GT.0) THEN
         WRITE(LU6,*) ' ERROR: Problem with SMOD line'
         STOP
         ENDIF
      ENDDO COUNT_SMOK_LOOP
    7 REWIND(LU5)
C
C If no smoke detectors, get out
C
      IF (NSD.EQ.0) GOTO 100
      NU_WORK = TWFIN*200
C
C Allocate arrays for smoke detectors
C
      ALLOCATE(SMOKE_DETECTOR(NSD),STAT=IZERO)
      CALL ChkMemErr('READ','SMOKE_DETECTOR',IZERO)
C
C Read in smoke detector locations and parameters line by line. Assign
C default parameters unless otherwise indicated.
C
      NSDO    = NSD
      N       = 0
      DT_SMOKE = TWFIN/REAL(NFRAMES)
C
      SDLOOP: DO NN=1,NSDO
C
      N = N+1
C
      LABEL = 'null'                       ! Reset defaults before each read
      ACTIVATION_OBSCURATION        =  1.0 ! Default setting of detector alarm 1.0%/m
C     Default setting of a typical ionization detector
      ALPHA_E  =  0.0
      BETA_E   = -1.0
      ALPHA_C  =  1.8 
      BETA_C   = -1.0
C
      CALL CHECKREAD('SMOD',LU5,IOS) ; IF (IOS.EQ.1) EXIT SDLOOP
      READ(LU5,SMOD)
C
C Determine which mesh smoke detector is in (if any)
C
      MESH_LOOP: DO NM=1,NMESHES
      IF (XYZ(1).GT.MESH(NM)%XS .AND. XYZ(1).LT.MESH(NM)%XF .AND.
     .    XYZ(2).GT.MESH(NM)%YS .AND. XYZ(2).LT.MESH(NM)%YF .AND.
     .    XYZ(3).GT.MESH(NM)%ZS .AND. XYZ(3).LT.MESH(NM)%ZF) THEN
         SMOKE_DETECTOR(N)%MESH = NM
         EXIT MESH_LOOP
         ENDIF
      IF(NM.EQ.NMESHES) THEN
         N     = N-1
         NSD   = NSD-1 
         CYCLE SDLOOP 
         ENDIF
      ENDDO MESH_LOOP
C
C Apply properties of SMOKE detector N
C
      ALLOCATE(SMOKE_DETECTOR(N)%YSMOK_E(NU_WORK),STAT=IZERO)
      CALL ChkMemErr('READ','YSMOK_E',IZERO)
      ALLOCATE(SMOKE_DETECTOR(N)%TAR_CURVE(NU_WORK),STAT=IZERO)
      CALL ChkMemErr('READ','TAR_CURVE',IZERO)
      SMOKE_DETECTOR(N)%YSMOK_E   = 0.0
      SMOKE_DETECTOR(N)%TAR_CURVE = 1.0E+30
      SMOKE_DETECTOR(N)%YSMOKE_IN = 0. 
C
      SMOKE_DETECTOR(N)%ALPHA_E = ALPHA_E
      SMOKE_DETECTOR(N)%ALPHA_C = ALPHA_C
      SMOKE_DETECTOR(N)%BETA_E  = BETA_E
      SMOKE_DETECTOR(N)%BETA_C  = BETA_C
C
      SMOKE_DETECTOR(N)%X = XYZ(1)
      SMOKE_DETECTOR(N)%Y = XYZ(2)
      SMOKE_DETECTOR(N)%Z = XYZ(3)
C
      SMOKE_DETECTOR(N)%YSMOKE_ACT = ACTIVATION_OBSCURATION
C
      IF (LABEL.NE.'null') THEN
         SMOKE_DETECTOR(N)%LABEL = LABEL
      ELSE
         WRITE(SMOKE_DETECTOR(N)%LABEL,'(A,I4.4)') 'SMOKE',N
      ENDIF
C
      ENDDO SDLOOP
C
  100 REWIND(LU5)
C
      END SUBROUTINE READ_SMOD
C
C
      SUBROUTINE READ_THCP
C
      INTEGER NM,MESH_NUMBER,K_LOW,K_HIGH
      REAL(EB) DEPTH,DIAMETER,EMISSIVITY
      TYPE (THERMOCOUPLE_TYPE), POINTER :: TC
      NAMELIST /THCP/ XYZ,DTSAM,QUANTITY,IOR,LABEL,FYI,
     .                XB,DEPTH,DIAMETER,EMISSIVITY,K_LOW,K_HIGH
C
      NTC = 0
      TCLOOP: DO
      CALL CHECKREAD('THCP',LU5,IOS) ; IF (IOS.EQ.1) EXIT TCLOOP
      READ(LU5,NML=THCP,END=11,ERR=12,IOSTAT=IOS)
      NTC = NTC + 1
   12 IF (IOS.GT.0) THEN
         WRITE(MESSAGE,'(A,I4)') 'ERROR: Problem with THCP no.',NTC+1
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO TCLOOP
   11 REWIND(LU5)
C
      DTSAM  = TWFIN/REAL(NFRAMES)
C
      IF (NTC.GT.0) THEN
C
      ALLOCATE(THERMOCOUPLE(NTC),STAT=IZERO)
      CALL ChkMemErr('READ','THERMOCOUPLE',IZERO)
C
      THERMOCOUPLE(1:NTC)%INDEX = 5
      THERMOCOUPLE(1:NTC)%VALUE = 0.
      THERMOCOUPLE(1:NTC)%COUNT = 0
      THERMOCOUPLE(1:NTC)%IOR   = 0
      THERMOCOUPLE(1:NTC)%IW    = 0
      THERMOCOUPLE(1:NTC)%I_DEPTH = 1
      THERMOCOUPLE(1:NTC)%I1    = -1
      THERMOCOUPLE(1:NTC)%I2    = -1
      THERMOCOUPLE(1:NTC)%J1    = -1
      THERMOCOUPLE(1:NTC)%J2    = -1
      THERMOCOUPLE(1:NTC)%K1    = -1
      THERMOCOUPLE(1:NTC)%K2    = -1
      THERMOCOUPLE(1:NTC)%K_LOW = -1
      THERMOCOUPLE(1:NTC)%K_HIGH= -1
C
      NTCO = NTC
      N    = 0
      K_LOW = -1
      K_HIGH= -1
C
      THCPLOOP: DO NN=1,NTCO
      N    = N+1
      IOR  = 0
      DEPTH= 0.
      DIAMETER = 0.001
      EMISSIVITY=0.85
      XB   = -1.E6
      SELECT CASE(N)
      CASE(1:9)        ; WRITE(LABEL,'(A3,I1)') 'TC ',N
      CASE(10:99)      ; WRITE(LABEL,'(A3,I2)') 'TC ',N
      CASE(100:999)    ; WRITE(LABEL,'(A3,I3)') 'TC ',N
      CASE(1000:9999)  ; WRITE(LABEL,'(A3,I4)') 'TC ',N
      END SELECT
C
      CALL CHECKREAD('THCP',LU5,IOS) ; IF (IOS.EQ.1) EXIT THCPLOOP
      READ(LU5,THCP,END=41,ERR=42) 
C
      IF (XB(1).GT.-1.E5) THEN
         XYZ(1) = 0.5*(XB(1)+XB(2))
         XYZ(2) = 0.5*(XB(3)+XB(4))
         XYZ(3) = 0.5*(XB(5)+XB(6))
      ENDIF
C
C Check for bad THCP lines
C
      BAD = .FALSE.
C
      MESH_LOOP: DO NM=1,NMESHES
      M=>MESH(NM)
      IF (XYZ(1).GE.M%XS .AND. XYZ(1).LE.M%XF .AND.
     .    XYZ(2).GE.M%YS .AND. XYZ(2).LE.M%YF .AND.
     .    XYZ(3).GE.M%ZS .AND. XYZ(3).LE.M%ZF) THEN
        MESH_NUMBER = NM
        EXIT MESH_LOOP
        ENDIF
      IF (NM.EQ.NMESHES) BAD = .TRUE.
      ENDDO MESH_LOOP
C
      SUCCESS = .FALSE.
      SEARCH1: DO ND=-NDATA,NDATA
      IF (QUANTITY.EQ.CDATA(ND)) THEN
         SUCCESS  = .TRUE.
         EXIT SEARCH1
         ENDIF
      ENDDO SEARCH1
      IF (.NOT.SUCCESS) THEN
         WRITE(MESSAGE,'(3A)') 
     .   ' WARNING: THCP quantity ',TRIM(QUANTITY),' not found'
         BAD = .TRUE.
         ENDIF
C
      IF (LDATA(ND,2) .AND. .NOT.MIXTURE_FRACTION) BAD = .TRUE.
      IF (LDATA(ND,3) .AND. .NOT.EVAPORATION)      BAD = .TRUE.
C
      IF (BAD) THEN
         N  = N-1
         NTC= NTC-1
         CYCLE THCPLOOP
         ENDIF
C
C Assign indices for the thermocouple
C
      TC => THERMOCOUPLE(N)
C
      TC%ORDINAL = NN
      TC%MESH    = MESH_NUMBER
      TC%LABEL   = LABEL
C
      IF (K_LOW.GT.-1) THEN
      TC%K_LOW = K_LOW
      ELSE
      TC%K_LOW = 1
      ENDIF
C
      IF (K_HIGH.GT.-1) THEN
      TC%K_HIGH = K_HIGH
      ELSE
      TC%K_HIGH = MESH(MESH_NUMBER)%KBAR
      ENDIF
C
      SEARCH: DO ND=-NDATA,NDATA
      IF (QUANTITY.EQ.CDATA(ND)) THEN
         TC%INDEX = IDATA(ND)
         EXIT SEARCH
         ENDIF
      ENDDO SEARCH
C
      IF (TC%INDEX.LT.0 .AND. IOR.EQ.0) THEN
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: Specify orientation of THCP '
     .   ,NN,' using the parameter IOR'
         CALL SHUTDOWN(MESSAGE)
         ENDIF
C
      TC%X          = XYZ(1)
      TC%Y          = XYZ(2)
      TC%Z          = XYZ(3)
      TC%IOR        = IOR
      TC%DEPTH      = DEPTH
      TC%DIAMETER   = DIAMETER
      TC%EMISSIVITY = EMISSIVITY
C
      IF (XB(1).GT.-1.E5) THEN
      NM = TC%MESH
      M=>MESH(NM)
      XB(1) = MAX(XB(1),M%XS)
      XB(2) = MIN(XB(2),M%XF)
      XB(3) = MAX(XB(3),M%YS)
      XB(4) = MIN(XB(4),M%YF)
      XB(5) = MAX(XB(5),M%ZS)
      XB(6) = MIN(XB(6),M%ZF)
      TC%I1 = NINT( GINV(XB(1)-M%XS,1,NM)*M%RDXI)
      TC%I2 = NINT( GINV(XB(2)-M%XS,1,NM)*M%RDXI)
      TC%J1 = NINT( GINV(XB(3)-M%YS,2,NM)*M%RDETA)
      TC%J2 = NINT( GINV(XB(4)-M%YS,2,NM)*M%RDETA)
      TC%K1 = NINT( GINV(XB(5)-M%ZS,3,NM)*M%RDZETA)
      TC%K2 = NINT( GINV(XB(6)-M%ZS,3,NM)*M%RDZETA)
      IF (TC%I1.LT.TC%I2) TC%I1 = TC%I1 + 1
      IF (TC%J1.LT.TC%J2) TC%J1 = TC%J1 + 1
      IF (TC%K1.LT.TC%K2) TC%K1 = TC%K1 + 1
      IF (XB(1).EQ.XB(2)) TC%IOR = 1
      IF (XB(3).EQ.XB(4)) TC%IOR = 2
      IF (XB(5).EQ.XB(6)) TC%IOR = 3
      ENDIF
C
   42 ENDDO THCPLOOP
   41 REWIND(LU5)
C
      ENDIF
C
      DTTC = DTSAM
C
      END SUBROUTINE READ_THCP
C
C
      SUBROUTINE READ_PL3D
C
      NAMELIST /PL3D/ DTSAM,QUANTITIES,FYI,WRITE_XYZ
C
      QUANTITIES(1) = 'TEMPERATURE'
      QUANTITIES(2) = 'U-VELOCITY'
      QUANTITIES(3) = 'V-VELOCITY'
      QUANTITIES(4) = 'W-VELOCITY'
      QUANTITIES(5) = 'HRRPUV'
      IF (ISOTHERMAL) QUANTITIES(5) = 'PRESSURE'
      DTSAM         = TWFIN/5.
      WRITE_XYZ     = .FALSE.
C
      PL3D_LOOP: DO
      CALL CHECKREAD('PL3D',LU5,IOS) ; IF (IOS.EQ.1) EXIT PL3D_LOOP
      READ(LU5,PL3D,END=27,ERR=28,IOSTAT=IOS)
   28 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with PL3D line')
      ENDDO PL3D_LOOP
   27 REWIND(LU5)
C
      WPLT = DTSAM
      PLOOP: DO N=1,5
C
      DO ND=-NDATA,NDATA
      IF (QUANTITIES(N).EQ.CDATA(ND)) THEN
         IPLOT3D(N) = IDATA(ND)
         IF (LDATA(ND,2) .AND. .NOT.MIXTURE_FRACTION) THEN
            WRITE(MESSAGE,'(3A)') 
     .      'ERROR: PLOT3D quantity ',TRIM(QUANTITIES(N)),
     .      ' not appropriate for non-fire case'
            CALL SHUTDOWN(MESSAGE)
            ENDIF 
         IF (LDATA(ND,3) .AND. .NOT.EVAPORATION) THEN
            WRITE(MESSAGE,'(3A)') 
     .      'ERROR: PLOT3D quantity ',TRIM(QUANTITIES(N)),
     .      ' not appropriate for non-sprinkler case'
            CALL SHUTDOWN(MESSAGE)
            ENDIF 
         CYCLE PLOOP
         ENDIF
      ENDDO
      WRITE(MESSAGE,'(3A)') 'ERROR: PLOT3D quantity ',
     .        TRIM(QUANTITIES(N)),' not found'
      CALL SHUTDOWN(MESSAGE)
      ENDDO PLOOP
C
      END SUBROUTINE READ_PL3D
C
C
      SUBROUTINE READ_ISOF
C
      REAL(EB) VALUE(10)
      INTEGER NIF0,REDUCE_TRIANGLES
      NAMELIST /ISOF/ DTSAM,QUANTITY,FYI,VALUE,
     .                REDUCE_TRIANGLES,COLOR_QUANTITY
C
      NIF = 0
      IF (MIXTURE_FRACTION) NIF = 2
      NIF0 = NIF+1
C
      ISOLOOP: DO
      CALL CHECKREAD('ISOF',LU5,IOS) ; IF (IOS.EQ.1) EXIT ISOLOOP
      READ(LU5,NML=ISOF,END=9,ERR=10,IOSTAT=IOS)
      NIF = NIF + 1
   10 IF (IOS.GT.0) THEN
      IF (MIXTURE_FRACTION) THEN
         WRITE(MESSAGE,'(A,I2)') 'ERROR: Problem with ISOF no.',NIF-1
         ELSE
         WRITE(MESSAGE,'(A,I2)') 'ERROR: Problem with ISOF no.',NIF+1
         ENDIF
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO ISOLOOP
    9 REWIND(LU5)
C
      DTSAM  = TWFIN/REAL(NFRAMES)
C
      ALLOCATE(INDIF(NIF),STAT=IZERO)
      CALL ChkMemErr('READ','INDIF',IZERO)
      ALLOCATE(INDIF2(NIF),STAT=IZERO)
      CALL ChkMemErr('READ','INDIF2',IZERO) ; INDIF2 = 0
      ALLOCATE(ISOLEVEL(10,NIF),STAT=IZERO)
      CALL ChkMemErr('READ','ISOLEVEL',IZERO)
      ALLOCATE(NLEVELS(NIF),STAT=IZERO)
      CALL ChkMemErr('READ','NLEVELS',IZERO)
      ALLOCATE(REDUCETRIANGLES(NIF),STAT=IZERO)
      CALL ChkMemErr('READ','REDUCETRIANGLES',IZERO)
      REDUCETRIANGLES = 1
C
      IF (MIXTURE_FRACTION) THEN
         INDIF(1) = 50 + IFUEL
         NLEVELS(1) = 1
         ISOLEVEL(1,1) = Z_F
         INDIF(2) = 11
         NLEVELS(2) = 1
         ISOLEVEL(1,2) = HRRPUA_SHEET/(7.*MESH(1)%DXMIN)
         ENDIF
C
      ISOFLOOP: DO N=NIF0,NIF
      QUANTITY = 'nulliso'
      COLOR_QUANTITY = 'nulliso'
      VALUE    = -999.
      REDUCE_TRIANGLES = 1
C
      CALL CHECKREAD('ISOF',LU5,IOS) ; IF (IOS.EQ.1) EXIT ISOFLOOP
      READ(LU5,ISOF,END=43,ERR=44) 
C
      REDUCETRIANGLES(N) = REDUCE_TRIANGLES
      SUCCESS = .FALSE.
      SEARCH: DO ND=-NDATA,NDATA
      IF (QUANTITY.EQ.CDATA(ND)) THEN
         INDIF(N) = IDATA(ND)
         VALUE_LOOP: DO I=1,10
         IF (VALUE(I).EQ.-999.) EXIT VALUE_LOOP
         NLEVELS(N) = I
         ISOLEVEL(I,N) = VALUE(I)
         ENDDO VALUE_LOOP
         SUCCESS = .TRUE.
         EXIT SEARCH
         ENDIF
      ENDDO SEARCH
C
      IF (.NOT.SUCCESS) THEN
      WRITE(MESSAGE,'(3A)') 
     .   'ERROR: ISOF quantity ',TRIM(QUANTITY),' not found'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
C
      SEARCH2: DO ND=-NDATA,NDATA
      IF (COLOR_QUANTITY.EQ.CDATA(ND)) THEN
         INDIF2(N) = IDATA(ND)
         EXIT SEARCH2
         ENDIF
      ENDDO SEARCH2
C
   44 ENDDO ISOFLOOP
   43 REWIND(LU5)
C
      DTIF = DTSAM
C
      IF (DTSAM.GT.TWFIN) NIF = 0
C
      END SUBROUTINE READ_ISOF
C
C
      SUBROUTINE READ_SLCF
C
      INTEGER :: NM,MESH_NUMBER,K_LOW,K_HIGH
C
      LOGICAL VECTOR
      NAMELIST /SLCF/ DTSAM,XB,QUANTITY,CB,FYI,
     .                PBX,PBY,PBZ,VECTOR,MESH_NUMBER,K_LOW,K_HIGH
C
      ISPDIAG = 0
C
      MESH_LOOP: DO NM=1,NMESHES
      M=>MESH(NM)
      CALL UNPACK_VAR(NM)
C
      NSF  = 0
      NSFO = 0
      SLCLOOP: DO
      VECTOR  = .FALSE.
      MESH_NUMBER=NM
      CALL CHECKREAD('SLCF',LU5,IOS) ; IF (IOS.EQ.1) EXIT SLCLOOP
      READ(LU5,NML=SLCF,END=9,ERR=10,IOSTAT=IOS)
      NSFO = NSFO + 1
      IF (MESH_NUMBER.NE.NM) CYCLE SLCLOOP
      NSF  = NSF + 1
      IF (VECTOR .AND. JBAR.EQ.1) NSF = NSF + 2
      IF (VECTOR .AND. JBAR.GT.1) NSF = NSF + 3
   10 IF (IOS.GT.0) THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with SLCF no.',NSF+1
         CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO SLCLOOP
    9 REWIND(LU5)
C
      DTSAM  = TWFIN/REAL(NFRAMES)
C
      ALLOCATE(M%ISP1(NSF),STAT=IZERO)
      CALL ChkMemErr('READ','ISP1',IZERO)
      ALLOCATE(M%ISP2(NSF),STAT=IZERO)
      CALL ChkMemErr('READ','ISP1',IZERO)
      ALLOCATE(M%JSP1(NSF),STAT=IZERO)
      CALL ChkMemErr('READ','ISP1',IZERO)
      ALLOCATE(M%JSP2(NSF),STAT=IZERO)
      CALL ChkMemErr('READ','ISP1',IZERO)
      ALLOCATE(M%KSP1(NSF),STAT=IZERO)
      CALL ChkMemErr('READ','ISP1',IZERO)
      ALLOCATE(M%KSP2(NSF),STAT=IZERO)
      CALL ChkMemErr('READ','ISP1',IZERO)
      ALLOCATE(M%INDSP(NSF),STAT=IZERO)
      CALL ChkMemErr('READ','ISP1',IZERO)
C
      CALL UNPACK_VAR(NM)
C
      N = 0
      K_LOW  = 1
      K_HIGH = KBAR
C
      SLCFLOOP: DO NN=1,NSFO
      QUANTITY = 'null'
      CB     = 'null'
      PBX    = -1.E6
      PBY    = -1.E6
      PBZ    = -1.E6
      VECTOR = .FALSE.
      MESH_NUMBER=NM
C
      CALL CHECKREAD('SLCF',LU5,IOS) ; IF (IOS.EQ.1) EXIT SLCFLOOP
      READ(LU5,SLCF,END=43,ERR=44) 
      IF (MESH_NUMBER.NE.NM) CYCLE SLCFLOOP
C
      IF (PBX.GT.-1.E5 .OR. PBY.GT.-1.E5 .OR. PBZ.GT.-1.E5) THEN
      XB(1) = XS
      XB(2) = XF
      XB(3) = YS
      XB(4) = YF
      XB(5) = ZS
      XB(6) = ZF
      IF (PBX.GT.-1.E5) XB(1:2) = PBX
      IF (PBY.GT.-1.E5) XB(3:4) = PBY
      IF (PBZ.GT.-1.E5) XB(5:6) = PBZ
      ENDIF
C
      IF (CB.NE.'null') THEN
      XB(1) = XS
      XB(2) = XF
      XB(3) = YS
      XB(4) = YF
      XB(5) = ZS
      XB(6) = ZF
      IF (CB.EQ.'XBAR0' .OR. CB.EQ.'RBAR0') XB(2) = XS
      IF (CB.EQ.'XBAR'  .OR. CB.EQ.'RBAR')  XB(1) = XF
      IF (CB.EQ.'YBAR0') XB(4) = YS
      IF (CB.EQ.'YBAR' ) XB(3) = YF
      IF (CB.EQ.'ZBAR0') XB(6) = ZS
      IF (CB.EQ.'ZBAR' ) XB(5) = ZF
      ENDIF
C
      XB(1) = MAX(XB(1),XS)
      XB(2) = MIN(XB(2),XF)
      XB(3) = MAX(XB(3),YS)
      XB(4) = MIN(XB(4),YF)
      XB(5) = MAX(XB(5),ZS)
      XB(6) = MIN(XB(6),ZF)
C
C Throw out bad slices
C
      BAD = .FALSE.
C
      IF (XB(1).GT.XF .OR. XB(2).LT.XS .OR.
     .    XB(3).GT.YF .OR. XB(4).LT.YS .OR.
     .    XB(5).GT.ZF .OR. XB(6).LT.ZS) BAD = .TRUE.
C
      SUCCESS = .FALSE.
      SEARCH1: DO ND=-NDATA,NDATA
      IF (QUANTITY.EQ.CDATA(ND)) THEN
         SUCCESS = .TRUE.
         EXIT SEARCH1
         ENDIF
      ENDDO SEARCH1
      IF (.NOT.SUCCESS) THEN
      BAD = .TRUE.
      WRITE(MESSAGE,'(3A)') 
     .    ' WARNING: SLCF quantity ',TRIM(QUANTITY),' not found'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
C
      IF (LDATA(ND,2) .AND. .NOT.MIXTURE_FRACTION)  BAD = .TRUE.
      IF (LDATA(ND,3) .AND. .NOT.EVAPORATION)       BAD = .TRUE.
C
      IF (BAD) THEN
          NSF = NSF - 1
          IF (VECTOR .AND. JBAR.EQ.1) NSF = NSF - 2
          IF (VECTOR .AND. JBAR.GT.1) NSF = NSF - 3
          CYCLE SLCFLOOP
          ENDIF
C
C Process vector quantities
C
      NITER = 1
      IF (VECTOR .AND. JBAR.EQ.1) NITER = 3
      IF (VECTOR .AND. JBAR.GT.1) NITER = 4
C
      VECTORLOOP: DO ITER=1,NITER
C
      N = N + 1
      ISP1(N) = NINT( GINV(XB(1)-XS,1,NM)*RDXI)
      ISP2(N) = NINT( GINV(XB(2)-XS,1,NM)*RDXI)
      JSP1(N) = NINT( GINV(XB(3)-YS,2,NM)*RDETA)
      JSP2(N) = NINT( GINV(XB(4)-YS,2,NM)*RDETA)
      KSP1(N) = NINT( GINV(XB(5)-ZS,3,NM)*RDZETA)
      KSP2(N) = NINT( GINV(XB(6)-ZS,3,NM)*RDZETA)
C
      IF (ITER.EQ.2) QUANTITY = 'U-VELOCITY' 
      IF (ITER.EQ.3 .AND. JBAR.GT.1) QUANTITY = 'V-VELOCITY' 
      IF (ITER.EQ.3 .AND. JBAR.EQ.1) QUANTITY = 'W-VELOCITY' 
      IF (ITER.EQ.4) QUANTITY = 'W-VELOCITY' 
C
      SEARCH: DO ND=-NDATA,NDATA
      IF (QUANTITY.EQ.CDATA(ND)) THEN
         INDSP(N) = IDATA(ND)
         EXIT SEARCH
         ENDIF
      ENDDO SEARCH
C
      ENDDO VECTORLOOP
C 
   44 ENDDO SLCFLOOP
   43 REWIND(LU5)
C
      DTSF = DTSAM
C
      LAYER_HEIGHT=.FALSE.
      K_LAYER_CALC_LO = K_LOW
      K_LAYER_CALC_HI = K_HIGH
      DO N=1,NSF
      IF (INDSP(N).GE.31 .AND. INDSP(N).LE.33) ISPDIAG=1
      IF (INDSP(N).GE.105 .AND. INDSP(N).LE.109) LAYER_HEIGHT=.TRUE.
      ENDDO
C
      ENDDO MESH_LOOP
C
      END SUBROUTINE READ_SLCF
C
C
      SUBROUTINE READ_BNDF
C
      NAMELIST /BNDF/ DTSAM,QUANTITY,FYI
C
      NBF = 0
      BFLOOP: DO
      CALL CHECKREAD('BNDF',LU5,IOS) ; IF (IOS.EQ.1) EXIT BFLOOP
      READ(LU5,NML=BNDF,END=209,ERR=210,IOSTAT=IOS)
      NBF = NBF + 1
  210 IF (IOS.GT.0) CALL SHUTDOWN('ERROR: Problem with BNDF line')
      ENDDO BFLOOP
  209 REWIND(LU5)
C
      DTSAM  = 2.*TWFIN/REAL(NFRAMES)
      ACCUMULATE_WATER = .FALSE.
C
      ALLOCATE(INDBF(NBF),STAT=IZERO)
      CALL ChkMemErr('READ','INDBG',IZERO)
C
      BNDFLOOP: DO N=1,NBF
      QUANTITY = 'WALL_TEMPERATURE'
      CALL CHECKREAD('BNDF',LU5,IOS) ; IF (IOS.EQ.1) EXIT BNDFLOOP
      READ(LU5,BNDF,END=243,ERR=244)
C
      IF (QUANTITY.EQ.'AWMPUA') ACCUMULATE_WATER = .TRUE.
C
      SUCCESS = .FALSE.
      SEARCH: DO ND=-NDATA,NDATA
      IF (QUANTITY.EQ.CDATA(ND)) THEN
         INDBF(N) = IDATA(ND)
         SUCCESS  = .TRUE.
         EXIT SEARCH
         ENDIF
      ENDDO SEARCH
      IF (.NOT.SUCCESS) THEN
      WRITE(MESSAGE,'(3A)') 
     .   'ERROR: BNDF quantity ',TRIM(QUANTITY),' not found'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
C
      IF (INDBF(N).GT.0 .OR. INDBF(N).EQ.-6) 
     .   CALL SHUTDOWN('ERROR: BNDF QUANTITY not appropriate')
C
  244 ENDDO BNDFLOOP
  243 REWIND(LU5)
C
      DTBF = DTSAM
C
      END SUBROUTINE READ_BNDF
C
C
      SUBROUTINE DEFINE_OUTPUT_QUANTITIES
C
      ALLOCATE(CDATA(-NDATA:NDATA),STAT=IZERO)
      CALL ChkMemErr('READ','CDATA',IZERO) ; CDATA = 'null'
      ALLOCATE(SDATA(-NDATA:NDATA),STAT=IZERO)
      CALL ChkMemErr('READ','SDATA',IZERO) ; SDATA = 'null'
      ALLOCATE(UDATA(-NDATA:NDATA),STAT=IZERO)
      CALL ChkMemErr('READ','UDATA',IZERO) ; UDATA = 'null'
      ALLOCATE(IDATA(-NDATA:NDATA),STAT=IZERO)
      CALL ChkMemErr('READ','IDATA',IZERO) ; IDATA = 0
      ALLOCATE(LDATA(-NDATA:NDATA,5),STAT=IZERO)
      CALL ChkMemErr('READ','IDATA',IZERO) ; LDATA = .FALSE.
C
      CDATA(0)  = 'SMOKE/WATER'             ; IDATA(0)  =  0
      UDATA(0)  = '  '                      ; SDATA(0)  = '  '
      CDATA(1)  = 'DENSITY'                 ; IDATA(1)  =  1
      UDATA(1)  = 'kg/m3'                   ; SDATA(1)  = 'density'
      LDATA(1,1)= .TRUE.
      CDATA(2)  = 'U-MOMENTUM'              ; IDATA(2)  =  2
      UDATA(2)  = 'kg/s/m2'                 ; SDATA(2)  = 'u-mom'
      LDATA(2,1)= .TRUE.
      CDATA(3)  = 'V-MOMENTUM'              ; IDATA(3)  =  3
      UDATA(3)  = 'kg/s/m2'                 ; SDATA(3)  = 'v-mom'
      LDATA(3,1)= .TRUE.
      CDATA(4)  = 'W-MOMENTUM'              ; IDATA(4)  =  4
      UDATA(4)  = 'kg/s/m2'                 ; SDATA(4)  = 'w-mom'
      LDATA(4,1)= .TRUE.
      CDATA(5)  = 'TEMPERATURE'             ; IDATA(5)  =  5
      UDATA(5)  = 'C'                       ; SDATA(5)  = 'temp'
      CDATA(6)  = 'U-VELOCITY'              ; IDATA(6)  =  6
      UDATA(6)  = 'm/s'                     ; SDATA(6)  = 'U-VEL'
      CDATA(7)  = 'V-VELOCITY'              ; IDATA(7)  =  7
      UDATA(7)  = 'm/s'                     ; SDATA(7)  = 'V-VEL'
      CDATA(8)  = 'W-VELOCITY'              ; IDATA(8)  =  8
      UDATA(8)  = 'm/s'                     ; SDATA(8)  = 'W-VEL'
      CDATA(9)  = 'PRESSURE'                ; IDATA(9)  =  9
      UDATA(9)  = 'Pa'                      ; SDATA(9)  = 'pres'
      CDATA(10) = 'VELOCITY'                ; IDATA(10) = 10
      UDATA(10) = 'm/s'                     ; SDATA(10) = 'vel'
      CDATA(11) = 'HRRPUV'                  ; IDATA(11) = 11
      UDATA(11) = 'kW/m3'                   ; SDATA(11) = 'hrrpuv'
      CDATA(12) = 'H'                       ; IDATA(12) = 12
      UDATA(12) = '(m/s)^2'                 ; SDATA(12) = 'head'
      CDATA(13) = 'KINEMATIC_VISCOSITY'     ; IDATA(13) = 13
      UDATA(13) = 'm2/s'                    ; SDATA(13) = 'visc'
      CDATA(14) = 'DIVERGENCE'              ; IDATA(14) = 14
      UDATA(14) = '1/s'                     ; SDATA(14) = 'div'
      CDATA(15) = 'WMPUV'                   ; IDATA(15) = 15
      UDATA(15) = 'kg/m3'                   ; SDATA(15) = 'wmpuv'
      LDATA(15,3) = .TRUE.
      CDATA(16) = 'ABSORPTION_COEFFICIENT'  ; IDATA(16) = 16
      UDATA(16) = '1/m'                     ; SDATA(16) = 'kappa'
      CDATA(17) = 'DYNAMIC_VISCOSITY'       ; IDATA(17) = 17
      UDATA(17) = 'kg/m/s'                  ; SDATA(17) = 'visc'
C
C Radiation
C
      CDATA(18) = 'RADIANT_INTENSITY'       ; IDATA(18) = 18
      UDATA(18) = 'kW/m^2'                  ; SDATA(18) = 'inten'
      CDATA(19) = 'RADIATION_LOSS'          ; IDATA(19) = 19
      UDATA(19) = 'kW/m3'                   ; SDATA(19) = 'loss'
      CDATA(20) = 'RADIANT_FLUX'            ; IDATA(20) = 20
      UDATA(20) = 'kW/m^2'                  ; SDATA(20) = 'flux'
      CDATA(21) = 'RADIANT_FLUX_X'          ; IDATA(21) = 21
      UDATA(21) = 'kW/m^2'                  ; SDATA(21) = 'flux_x'
      CDATA(22) = 'RADIANT_FLUX_Y'          ; IDATA(22) = 22
      UDATA(22) = 'kW/m^2'                  ; SDATA(22) = 'flux_y'
      CDATA(23) = 'RADIANT_FLUX_Z'          ; IDATA(23) = 23
      UDATA(23) = 'kW/m^2'                  ; SDATA(23) = 'flux_z'
C
C Strain and Vorticity
C
      CDATA(24) = 'STRAIN_RATE_X'           ; IDATA(24) = 24
      UDATA(24) = '1/s'                     ; SDATA(24) = 'strain_x'
      CDATA(25) = 'STRAIN_RATE_Y'           ; IDATA(25) = 25
      UDATA(25) = '1/s'                     ; SDATA(25) = 'strain_y'
      CDATA(26) = 'STRAIN_RATE_Z'           ; IDATA(26) = 26
      UDATA(26) = '1/s'                     ; SDATA(26) = 'strain_z'
      CDATA(27) = 'VORTICITY_X'             ; IDATA(27) = 27
      UDATA(27) = '1/s'                     ; SDATA(27) = 'vort_x'  
      CDATA(28) = 'VORTICITY_Y'             ; IDATA(28) = 28
      UDATA(28) = '1/s'                     ; SDATA(28) = 'vort_y'  
      CDATA(29) = 'VORTICITY_Z'             ; IDATA(29) = 29
      UDATA(29) = '1/s'                     ; SDATA(29) = 'vort_z'  
C
C Water
C
      CDATA(31) = 'DROPLET_FLUX_X'          ; IDATA(31) = 31
      UDATA(31) = 'kg/s/m2'                 ; SDATA(31) = 'xflux'
      CDATA(32) = 'DROPLET_FLUX_Y'          ; IDATA(32) = 32
      UDATA(32) = 'kg/s/m2'                 ; SDATA(32) = 'yflux'
      CDATA(33) = 'DROPLET_FLUX_Z'          ; IDATA(33) = 33
      UDATA(33) = 'kg/s/m2'                 ; SDATA(33) = 'zflux'
      CDATA(34) = 'DIAMETER'                ; IDATA(34) = 34
      UDATA(34) = 'mu-m'                    ; SDATA(34) = 'diam'
      LDATA(31:34,3) = .TRUE.
      CDATA(36) = 'DROPLET_PHASE'           ; IDATA(36) = 36
      UDATA(36) = ' '                       ; SDATA(36) = 'phase'
      CDATA(37) = 'DROPLET_TEMPERATURE'     ; IDATA(37) = 37
      UDATA(37) = 'C'                       ; SDATA(37) = 'temp'
      LDATA(36:37,3) = .TRUE. 
      CDATA(38) = 'WATER_RADIATION_LOSS'    ; IDATA(38) = 38
      UDATA(38) = 'kW/m3'                   ; SDATA(38) = 'H2O_rad'
      CDATA(39) = 'AGE'                     ; IDATA(39) = 39
      UDATA(39) = 's'                       ; SDATA(39) = 'age'
C
C Mixture Fraction related variables
C
      CDATA(41) = 'fuel'                    ; IDATA(41) = 41
      UDATA(41) = 'mol/mol'                 ; SDATA(41) = 'X_f'
      CDATA(42) = 'oxygen'                  ; IDATA(42) = 42
      UDATA(42) = 'mol/mol'                 ; SDATA(42) = 'X_O2'
      CDATA(43) = 'nitrogen'                ; IDATA(43) = 43
      UDATA(43) = 'mol/mol'                 ; SDATA(43) = 'X_N2'
      CDATA(44) = 'water vapor'             ; IDATA(44) = 44
      UDATA(44) = 'mol/mol'                 ; SDATA(44) = 'X_H2O'
      CDATA(45) = 'carbon dioxide'          ; IDATA(45) = 45
      UDATA(45) = 'mol/mol'                 ; SDATA(45) = 'X_CO2'
      CDATA(46) = 'carbon monoxide'         ; IDATA(46) = 46
      UDATA(46) = 'ppm'                     ; SDATA(46) = 'X_CO'
      CDATA(47) = 'hydrogen'                ; IDATA(47) = 47
      UDATA(47) = 'ppm'                     ; SDATA(47) = 'X_H2'
      CDATA(48) = 'soot volume fraction'    ; IDATA(48) = 48
      UDATA(48) = 'ppm'                     ; SDATA(48) = 'f_v'
      LDATA(48,1)=.TRUE.
      CDATA(49) = 'extinction coefficient'  ; IDATA(49) = 49
      UDATA(49) = '1/m'                     ; SDATA(49) = 'ext'
      LDATA(49,1)=.TRUE.
      CDATA(50) = 'visibility'              ; IDATA(50) = 50
      UDATA(50) = 'm'                       ; SDATA(50) = 'vis'
      LDATA(50,1)=.TRUE.
C
      LDATA(41:50,2) = .TRUE.
C
C Integrated Quantities
C
      CDATA(101)  = 'VOLUME FLOW'           ; IDATA(101)  =  101
      UDATA(101)  = 'm3/s'                  ; SDATA(101)  = 'vflow'
      CDATA(102)  = 'MASS FLOW'             ; IDATA(102)  =  102
      UDATA(102)  = 'kg/s'                  ; SDATA(102)  = 'mflow'
      CDATA(103)  = 'HEAT FLOW'             ; IDATA(103)  =  103
      UDATA(103)  = 'kW'                    ; SDATA(103)  = 'hflow'
      CDATA(104)  = 'HRR'                   ; IDATA(104)  =  104
      UDATA(104)  = 'kW'                    ; SDATA(104)  = 'hrr'
      CDATA(105)  = 'LAYER HEIGHT'          ; IDATA(105)  =  105
      UDATA(105)  = 'm'                     ; SDATA(105)  = 'layer'
      CDATA(106)  = 'UPPER TEMPERATURE'     ; IDATA(106)  =  106
      UDATA(106)  = 'C'                     ; SDATA(106)  = 'u-tmp'
      CDATA(107)  = 'LOWER TEMPERATURE'     ; IDATA(107)  =  107
      UDATA(107)  = 'C'                     ; SDATA(107)  = 'l-tmp'
      CDATA(108)  = 'UPPER KAPPA'           ; IDATA(108)  =  108
      UDATA(108)  = '1/m'                   ; SDATA(108)  = 'u-kap'
      CDATA(109)  = 'UPPER INTENSITY'       ; IDATA(109)  =  109
      UDATA(109)  = 'kW/m2'                 ; SDATA(109)  = 'u-rad'
C
C Thermocouple
C
      CDATA(110)  = 'THERMOCOUPLE'          ; IDATA(110)  =  110
      UDATA(110)  = 'C'                     ; SDATA(110)  = 'tc'
      LDATA(110,1)= .TRUE.
C
      CDATA(111)  = 'VOLUME FLOW +'         ; IDATA(111)  =  111
      UDATA(111)  = 'm3/s'                  ; SDATA(111)  = 'vflow+'
      CDATA(112)  = 'MASS FLOW +'           ; IDATA(112)  =  112
      UDATA(112)  = 'kg/s'                  ; SDATA(112)  = 'mflow+'
      CDATA(113)  = 'HEAT FLOW +'           ; IDATA(113)  =  113
      UDATA(113)  = 'kW'                    ; SDATA(113)  = 'hflow+'
      CDATA(116)  = 'VOLUME FLOW -'         ; IDATA(116)  =  116
      UDATA(116)  = 'm3/s'                  ; SDATA(116)  = 'vflow-'
      CDATA(117)  = 'MASS FLOW -'           ; IDATA(117)  =  117
      UDATA(117)  = 'kg/s'                  ; SDATA(117)  = 'mflow-'
      CDATA(118)  = 'HEAT FLOW -'           ; IDATA(118)  =  118
      UDATA(118)  = 'kW'                    ; SDATA(118)  = 'hflow-'
C
C Alternative Mixture Fraction quantities
C
      CDATA(141) = 'fuel mass fraction'     ; IDATA(141) = 141
      UDATA(141) = 'kg/kg'                  ; SDATA(141) = 'Y_f'
      CDATA(142) = 'oxygen mass fraction'   ; IDATA(142) = 142
      UDATA(142) = 'kg/kg'                  ; SDATA(142) = 'Y_O2'
      CDATA(143) = 'nitrogen mass fraction' ; IDATA(143) = 143
      UDATA(143) = 'kg/kg'                  ; SDATA(143) = 'Y_N2'
      CDATA(144) = 'water vapor mass fraction' ; IDATA(144) = 144
      UDATA(144) = 'kg/kg'                  ; SDATA(144) = 'Y_H2O'
      CDATA(145) = 'carbon dioxide mass fraction'  ; IDATA(145) = 145
      UDATA(145) = 'kg/kg'                  ; SDATA(145) = 'Y_CO2'
      CDATA(146) = 'carbon monoxide mass fraction' ; IDATA(146) = 146
      UDATA(146) = 'kg/kg'                  ; SDATA(146) = 'Y_CO'
      CDATA(147) = 'hydrogen mass fraction'        ; IDATA(147) = 147
      UDATA(147) = 'kg/kg'                  ; SDATA(147) = 'Y_H2'
      CDATA(148) = 'soot density'           ; IDATA(148) = 148
      UDATA(148) = 'mg/m3'                  ; SDATA(148) = 'soot'
      LDATA(148,1)=.TRUE.
C
      LDATA(141:148,2) = .TRUE.
C
C Species Mass Fractions
C
      CDATA(51) = TRIM(SPECIES_ID(1))       ; IDATA(51) = 51
      UDATA(51) = 'kg/kg'                   ; SDATA(51) = 'spec_1'
      CDATA(52) = TRIM(SPECIES_ID(2))       ; IDATA(52) = 52
      UDATA(52) = 'kg/kg'                   ; SDATA(52) = 'spec_2'
      CDATA(53) = TRIM(SPECIES_ID(3))       ; IDATA(53) = 53
      UDATA(53) = 'kg/kg'                   ; SDATA(53) = 'spec_3'
      CDATA(54) = TRIM(SPECIES_ID(4))       ; IDATA(54) = 54
      UDATA(54) = 'kg/kg'                   ; SDATA(54) = 'spec_4'
      CDATA(55) = TRIM(SPECIES_ID(5))       ; IDATA(55) = 55
      UDATA(55) = 'kg/kg'                   ; SDATA(55) = 'spec_5'
C
C Species Fluxes in the x-direction
C
      CDATA(61) = TRIM(SPECIES_ID(1))//'_FLUX_X' ; IDATA(61) = 61
      UDATA(61) = 'kg/s/m2'                 ; SDATA(61) = 'u*rho*1'
      CDATA(62) = TRIM(SPECIES_ID(2))//'_FLUX_X' ; IDATA(62) = 62
      UDATA(62) = 'kg/s/m2'                 ; SDATA(62) = 'u*rho*2'
      CDATA(63) = TRIM(SPECIES_ID(3))//'_FLUX_X' ; IDATA(63) = 63
      UDATA(63) = 'kg/s/m2'                 ; SDATA(63) = 'u*rho*3'
      CDATA(64) = TRIM(SPECIES_ID(4))//'_FLUX_X' ; IDATA(64) = 64
      UDATA(64) = 'kg/s/m2'                 ; SDATA(64) = 'u*rho*4'
      CDATA(65) = TRIM(SPECIES_ID(5))//'_FLUX_X' ; IDATA(65) = 65
      UDATA(65) = 'kg/s/m2'                 ; SDATA(65) = 'u*rho*5'
      LDATA(61:65,1) = .TRUE.
C
C Species Fluxes in the y-direction
C
      CDATA(71) = TRIM(SPECIES_ID(1))//'_FLUX_Y' ; IDATA(71) = 71
      UDATA(71) = 'kg/s/m2'                 ; SDATA(71) = 'v*rho*1'
      CDATA(72) = TRIM(SPECIES_ID(2))//'_FLUX_Y' ; IDATA(72) = 72
      UDATA(72) = 'kg/s/m2'                ; SDATA(72) = 'v*rho*2'
      CDATA(73) = TRIM(SPECIES_ID(3))//'_FLUX_Y' ; IDATA(73) = 73
      UDATA(73) = 'kg/s/m2'                ; SDATA(73) = 'v*rho*3'
      CDATA(74) = TRIM(SPECIES_ID(4))//'_FLUX_Y' ; IDATA(74) = 74
      UDATA(74) = 'kg/s/m2'                ; SDATA(74) = 'v*rho*4'
      CDATA(75) = TRIM(SPECIES_ID(5))//'_FLUX_Y' ; IDATA(75) = 75
      UDATA(75) = 'kg/s/m2'                ; SDATA(75) = 'v*rho*5'
      LDATA(71:75,1) = .TRUE.
C
C Species Fluxes in the z-direction
C
      CDATA(81) = TRIM(SPECIES_ID(1))//'_FLUX_Z' ; IDATA(81) = 81
      UDATA(81) = 'kg/s/m2'                ; SDATA(81) = 'w*rho*1'
      CDATA(82) = TRIM(SPECIES_ID(2))//'_FLUX_Z' ; IDATA(82) = 82
      UDATA(82) = 'kg/s/m2'                ; SDATA(82) = 'w*rho*2'
      CDATA(83) = TRIM(SPECIES_ID(3))//'_FLUX_Z' ; IDATA(83) = 83
      UDATA(83) = 'kg/s/m2'                ; SDATA(83) = 'w*rho*3'
      CDATA(84) = TRIM(SPECIES_ID(4))//'_FLUX_Z' ; IDATA(84) = 84
      UDATA(84) = 'kg/s/m2'                ; SDATA(84) = 'w*rho*4'
      CDATA(85) = TRIM(SPECIES_ID(5))//'_FLUX_Z' ; IDATA(85) = 85
      UDATA(85) = 'kg/s/m2'                ; SDATA(85) = 'w*rho*5'
      LDATA(81:85,1) = .TRUE.
C
C Species Volume Fractions
C
      CDATA(91) = TRIM(SPECIES_ID(1))//'_VF'     ; IDATA(91) = 91
      UDATA(91) = '  '                      ; SDATA(91) = 'spec_1_vf'
      CDATA(92) = TRIM(SPECIES_ID(2))//'_VF'     ; IDATA(92) = 92
      UDATA(92) = '  '                      ; SDATA(92) = 'spec_2_vf'
      CDATA(93) = TRIM(SPECIES_ID(3))//'_VF'     ; IDATA(93) = 93
      UDATA(93) = '  '                      ; SDATA(93) = 'spec_3_vf'
      CDATA(94) = TRIM(SPECIES_ID(4))//'_VF'     ; IDATA(94) = 94
      UDATA(94) = '  '                      ; SDATA(94) = 'spec_4_vf'
      CDATA(95) = TRIM(SPECIES_ID(5))//'_VF'     ; IDATA(95) = 95
      UDATA(95) = '  '                      ; SDATA(95) = 'spec_5_vf'
      LDATA(91:95,1) = .TRUE.
C
C Boundary Quantities
C
      CDATA(-1) = 'RADIATIVE_FLUX'          ; IDATA(-1) = -1
      UDATA(-1) = 'kW/m2'                   ; SDATA(-1) = 'rad'
      CDATA(-2) = 'CONVECTIVE_FLUX'         ; IDATA(-2) = -2
      UDATA(-2) = 'kW/m2'                   ; SDATA(-2) = 'con'
      CDATA(-3) = 'NORMAL_VELOCITY'         ; IDATA(-3) = -3
      UDATA(-3) = 'm/s'                     ; SDATA(-3) = 'vel'
      CDATA(-4) = 'GAS_TEMPERATURE'         ; IDATA(-4) = -4
      UDATA(-4) = 'C'                       ; SDATA(-4) = 'temp'
      CDATA(-5) = 'WALL_TEMPERATURE'        ; IDATA(-5) = -5
      UDATA(-5) = 'C'                       ; SDATA(-5) = 'temp'
      CDATA(-6) = 'INSIDE_WALL_TEMPERATURE' ; IDATA(-6) = -6
      UDATA(-6) = 'C'                       ; SDATA(-6) = 'inside'
      CDATA(-7) = 'BURNING_RATE'            ; IDATA(-7) = -7
      UDATA(-7) = 'kg/m2/s'                 ; SDATA(-7) = 'burn'
      CDATA(-8) = 'WMPUA'                   ; IDATA(-8) = -8
      UDATA(-8) = 'kg/m2'                   ; SDATA(-8) = 'wmpua'
      CDATA(-9) = 'WCPUA'                   ; IDATA(-9) = -9
      UDATA(-9) = 'kW/m2'                   ; SDATA(-9) = 'wcpua'
      CDATA(-10)= 'HEAT_FLUX'               ; IDATA(-10)= -10
      UDATA(-10)= 'kW/m2'                   ; SDATA(-10)= 'heat'
      CDATA(-11)= 'PRESSURE_COEFFICIENT'    ; IDATA(-11)= -11
      UDATA(-11)= ' '                       ; SDATA(-11)= 'c_p'
      CDATA(-12)= 'BACK_WALL_TEMPERATURE'   ; IDATA(-12)= -12
      UDATA(-12)= 'C'                       ; SDATA(-12)= 'back'
      CDATA(-13)= 'GAUGE_HEAT_FLUX'         ; IDATA(-13)= -13
      UDATA(-13)= 'kW/m2'                   ; SDATA(-13)= 'gauge'
      CDATA(-14)= 'CONDUCTION'              ; IDATA(-14)= -14
      UDATA(-14)= 'kW/m2'                   ; SDATA(-14)= 'con'
      CDATA(-15)= 'PYROLYSIS'               ; IDATA(-15)= -15
      UDATA(-15)= 'kW/m2'                   ; SDATA(-15)= 'pyr'
      CDATA(-16)= 'MASS_LOSS'               ; IDATA(-16)= -16
      UDATA(-16)= 'kg/m2'                   ; SDATA(-16)= 'mloss'
      CDATA(-17)= 'BURN_DEPTH'              ; IDATA(-17)= -17
      UDATA(-17)= 'm'                       ; SDATA(-17)= 'b-depth'
      CDATA(-18)= 'AWMPUA'                  ; IDATA(-18) = -18
      UDATA(-18)= 'kg/m2'                   ; SDATA(-18) = 'awmpua'
      CDATA(-19)= 'INCIDENT_HEAT_FLUX'      ; IDATA(-19)= -19
      UDATA(-19)= 'kW/m2'                   ; SDATA(-19)= 'in_flux'
      CDATA(-20)= 'CHAR_DEPTH'              ; IDATA(-20)= -20
      UDATA(-20)= 'm'                       ; SDATA(-20)= 'c-depth'
C
      CDATA(-21)= TRIM(SPECIES_ID(1))//'_FLUX'   ; IDATA(-21) = -21
      UDATA(-21)= 'kg/s/m2'                 ; SDATA(-21) = 'spec1flux'
      CDATA(-22)= TRIM(SPECIES_ID(2))//'_FLUX'   ; IDATA(-22) = -22
      UDATA(-22)= 'kg/s/m2'                 ; SDATA(-22) = 'spec2flux'
      CDATA(-23)= TRIM(SPECIES_ID(3))//'_FLUX'   ; IDATA(-23) = -23
      UDATA(-23)= 'kg/s/m2'                 ; SDATA(-23) = 'spec3flux'
      CDATA(-24)= TRIM(SPECIES_ID(4))//'_FLUX'   ; IDATA(-24) = -24
      UDATA(-24)= 'kg/s/m2'                 ; SDATA(-24) = 'spec4flux'
      CDATA(-25)= TRIM(SPECIES_ID(5))//'_FLUX'   ; IDATA(-25) = -25
      UDATA(-25)= 'kg/s/m2'                 ; SDATA(-25) = 'spec5flux'
C
      CDATA(-30)= 'HEAT_TRANSFER_COEFFICIENT'    ; IDATA(-30)= -30
      UDATA(-30)= 'W/m2/K'                  ; SDATA(-30)= 'h'
      CDATA(-31)= 'RADIOMETER'              ; IDATA(-31)= -31
      UDATA(-31)= 'kW/m2'                   ; SDATA(-31)= 'radio'
C
      END SUBROUTINE DEFINE_OUTPUT_QUANTITIES
C
C
      SUBROUTINE SEARCH_KEYWORD(NAME,LU,IOS)
C
      INTEGER, INTENT(OUT) :: IOS
      INTEGER, INTENT(IN)  :: LU
      CHARACTER(*), INTENT(IN) :: NAME
      CHARACTER(40) TEXT
C
      IF (LU.LT.0) THEN
      IOS = -1
      RETURN
      ENDIF
C
      IOS = 1
      REWIND(LU)
      READLOOP: DO
      READ(LU,'(A)',END=10) TEXT
C
      IF (TRIM(TEXT).EQ.TRIM(NAME)) THEN
         IOS = 0
         RETURN
         ELSE
         CYCLE READLOOP
         ENDIF
C
      ENDDO READLOOP
C
   10 RETURN
      END SUBROUTINE SEARCH_KEYWORD
C
C
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP,IERROR)
C
C Solve a linear system of equations with Gauss-Jordon elimination
C
      INTEGER M,MP,N,NP
      REAL(EB) A(NP,NP),B(NP,MP)
      INTEGER, PARAMETER :: NMAX=50
      INTEGER I,ICOL,IROW,J,K,L,LL,INDXC(NMAX),INDXR(NMAX),
     .           IPIV(NMAX),IERROR
      REAL(EB) BIG,DUM,PIVINV
C
      DO J=1,N
      IPIV(J) = 0
      ENDDO
C
      DO 22 I=1,N
      BIG = 0.
      DO 13 J=1,N
      IF (IPIV(J).NE.1) THEN
         DO 12 K=1,N
         IF (IPIV(K).EQ.0) THEN
            IF (ABS(A(J,K)).GE.BIG) THEN
               BIG = ABS(A(J,K))
               IROW = J
               ICOL = K
               ENDIF
            ELSE IF (IPIV(K).GT.1) THEN
               IERROR = 103   ! Singular matrix in gaussj
               RETURN
               ENDIF
   12    ENDDO
      ENDIF
   13 ENDDO 
      IPIV(ICOL) = IPIV(ICOL) + 1
      IF (IROW.NE.ICOL) THEN
         DO 14 L=1,N
         DUM = A(IROW,L)
         A(IROW,L) = A(ICOL,L)
         A(ICOL,L) = DUM
   14    ENDDO 
         DO 15 L=1,M
         DUM = B(IROW,L)
         B(IROW,L) = B(ICOL,L)
         B(ICOL,L) = DUM
   15    ENDDO 
      ENDIF
      INDXR(I) = IROW
      INDXC(I) = ICOL
      IF (A(ICOL,ICOL).EQ.0.) THEN
         IERROR = 103  ! Singular matrix in gaussj
         RETURN
         ENDIF
      PIVINV = 1./A(ICOL,ICOL)
      A(ICOL,ICOL) = 1.
      DO 16 L=1,N
      A(ICOL,L) = A(ICOL,L)*PIVINV
   16 ENDDO 
      DO 17 L=1,M
      B(ICOL,L) = B(ICOL,L)*PIVINV
   17 ENDDO 
      DO 21 LL=1,N
      IF (LL.NE.ICOL) THEN
         DUM = A(LL,ICOL)
         A(LL,ICOL) = 0.
         DO 18 L=1,N
         A(LL,L) = A(LL,L) - A(ICOL,L)*DUM
   18    ENDDO
         DO 19 L=1,M
         B(LL,L) = B(LL,L) - B(ICOL,L)*DUM
   19    ENDDO 
      ENDIF
   21 ENDDO
   22 ENDDO 
      DO 24 L=N,1,-1
      IF (INDXR(L).NE.INDXC(L)) THEN
         DO 23 K=1,N
         DUM = A(K,INDXR(L))
         A(K,INDXR(L)) = A(K,INDXC(L))
         A(K,INDXC(L)) = DUM
   23    ENDDO 
      ENDIF
   24 ENDDO 
C
      END SUBROUTINE GAUSSJ
C
C
      END MODULE READ
