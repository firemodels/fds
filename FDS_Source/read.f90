MODULE READ_INPUT
 
USE PRECISION_PARAMETERS
USE MESH_VARIABLES
USE GLOBAL_CONSTANTS
USE TRAN
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: SECOND, CHECKREAD, SHUTDOWN, COLOR2RGB
USE MEMORY_FUNCTIONS, ONLY: ChkMemErr
USE COMP_FUNCTIONS, ONLY: GET_INPUT_FILE
USE EVAC, ONLY: READ_EVAC
 
IMPLICIT NONE
PRIVATE
CHARACTER(255), PARAMETER :: readid='$Id$'
CHARACTER(255), PARAMETER :: readrev='$Revision$'
CHARACTER(255), PARAMETER :: readdate='$Date$'

PUBLIC READ_DATA, GET_REV_read, COLOR2RGB

CHARACTER(30) :: LABEL,MB
CHARACTER(100) :: MESSAGE,FYI
CHARACTER(30) :: ID,SURF_DEFAULT,BACKGROUND_SPECIES,EVAC_SURF_DEFAULT
LOGICAL :: SUCCESS,EX,THICKEN_OBSTRUCTIONS,BAD
LOGICAL :: WATER_EVAPORATION=.FALSE.,FUEL_EVAPORATION=.FALSE.
REAL(EB) :: XB(6),TEXTURE_ORIGIN(3)
REAL(EB) :: PBX,PBY,PBZ,MW_BACKGROUND,HUMIDITY
REAL(EB) :: CP_BACKGROUND,H_BACKGROUND,TREF_BACKGROUND
REAL(EB) :: MU_USER(0:MAX_SPECIES),K_USER(0:MAX_SPECIES),D_USER(0:MAX_SPECIES),EPSK(0:MAX_SPECIES),SIG(0:MAX_SPECIES),MW_MIN,MW_MAX
INTEGER  :: I,J,K,IZERO,IOS
TYPE (MESH_TYPE), POINTER :: M
TYPE(OBSTRUCTION_TYPE), POINTER :: OB
TYPE (VENTS_TYPE), POINTER :: VT
TYPE(SPECIES_TYPE), POINTER :: SS,SS0
TYPE(SURFACE_TYPE), POINTER :: SF
TYPE(MATERIAL_TYPE), POINTER :: ML
TYPE(REACTION_TYPE), POINTER :: RN
 
 
CONTAINS
 
 
SUBROUTINE READ_DATA

! Create an array of output QUANTITY names that are included in the various NAMELIST groups
 
CALL DEFINE_OUTPUT_QUANTITIES

! Get the name of the input file by reading the command line argument

CALL GET_INPUT_FILE

! If no input file is given, just print out the version number and stop

IF (FN_INPUT(1:1)==' ') THEN
   IF (MYID==0) THEN
      WRITE(LU_ERR,'(/A)') "Fire Dynamics Simulator"
      IF (.NOT.PARALLEL)   WRITE(LU_ERR,'(/A,A,A)') "Version: ",TRIM(VERSION_STRING)," Serial"
      IF (PARALLEL)        WRITE(LU_ERR,'(/A,A,A)') "Version: ",TRIM(VERSION_STRING)," Parallel"
      WRITE(LU_ERR,'(A,I4)') "SVN Revision Number: ",SVN_REVISION_NUMBER
      WRITE(LU_ERR,'(A,A)') "Compile Date: ",TRIM(COMPILE_DATE)
      WRITE(LU_ERR,'(/A)')  "Consult FDS Users Guide Chapter, Running FDS, for further instructions."
      WRITE(LU_ERR,'(/A)')  "Hit Enter to Escape..."
   ENDIF
   READ(5,*)
   STOP
ENDIF

! Stop FDS if the input file cannot be found in the current directory

INQUIRE(FILE=FN_INPUT,EXIST=EX)
IF (.NOT.EX) THEN
   IF (MYID==0) WRITE(LU_ERR,'(A,A,A)') "ERROR: The file, ", TRIM(FN_INPUT),", does not exist in the current directory"
   STOP
ENDIF

! Open the input file

OPEN(LU_INPUT,FILE=FN_INPUT,ACTION='READ')

! Read the input file, NAMELIST group by NAMELIST group

CALL READ_DEAD    ! Scan input file looking for old NAMELIST groups, and stop the run if they exist
CALL READ_HEAD
CALL READ_MULT
CALL READ_MESH
CALL READ_TRAN
CALL READ_TIME
CALL READ_MISC
CALL READ_PRES
CALL READ_RADI
CALL READ_PROP
CALL READ_PART
CALL READ_DEVC
CALL READ_CTRL
CALL READ_TREE
CALL READ_MATL
CALL READ_SURF
CALL READ_OBST
CALL READ_VENT
CALL READ_REAC
CALL READ_SPEC
CALL READ_EVAC
CALL PROC_SPEC    ! Set up various SPECies constructs
CALL PROC_PART    ! Set up various PARTicle constructs
CALL PROC_SURF_1  ! Set up SURFace constructs for species
CALL READ_RAMP    ! Read in all RAMPs, assuming they have all been identified previously
CALL READ_TABL   ! Read in all TABLs, assuming they have all been identified previously
CALL PROC_MATL    ! Set up various MATeriaL constructs
CALL PROC_SURF_2  ! Set up remaining SURFace constructs
CALL READ_DUMP
CALL READ_CLIP
CALL READ_INIT
CALL READ_ZONE
CALL PROC_WALL    ! Set up grid for 1-D heat transfer in solids
CALL PROC_CTRL    ! Set up various ConTRoL constructs
CALL PROC_PROP    ! Set up various PROPerty constructs
CALL PROC_DEVC    ! Set up various DEViCe constructs
CALL READ_PROF
CALL READ_SLCF
CALL READ_ISOF
CALL READ_BNDF

! Close the input file, and never open it again
 
CLOSE (LU_INPUT)

! Set QUANTITY ambient values

CALL SET_QUANTITIES_AMBIENT
 
END SUBROUTINE READ_DATA



SUBROUTINE READ_DEAD
 
! Look for outdated NAMELIST groups and stop the run if any are found
 
REWIND(LU_INPUT)
CALL CHECKREAD('GRID',LU_INPUT,IOS)
IF (IOS==0) CALL SHUTDOWN('ERROR: GRID is no longer a valid NAMELIST group. Read User Guide discussion on MESH.')
REWIND(LU_INPUT)
CALL CHECKREAD('HEAT',LU_INPUT,IOS)
IF (IOS==0) CALL SHUTDOWN('ERROR: HEAT is no longer a valid NAMELIST group. Read User Guide discussion on PROP and DEVC.')
REWIND(LU_INPUT)
CALL CHECKREAD('PDIM',LU_INPUT,IOS)
IF (IOS==0) CALL SHUTDOWN('ERROR: PDIM is no longer a valid NAMELIST group. Read User Guide discussion on MESH.')
REWIND(LU_INPUT)
CALL CHECKREAD('PIPE',LU_INPUT,IOS)
IF (IOS==0) CALL SHUTDOWN('ERROR: PIPE is no longer a valid NAMELIST group. Read User Guide discussion on PROP and DEVC.')
REWIND(LU_INPUT)
CALL CHECKREAD('PL3D',LU_INPUT,IOS)
IF (IOS==0) CALL SHUTDOWN('ERROR: PL3D is no longer a valid NAMELIST group. Read User Guide discussion on DUMP.')
REWIND(LU_INPUT)
CALL CHECKREAD('SMOD',LU_INPUT,IOS)
IF (IOS==0) CALL SHUTDOWN('ERROR: SMOD is no longer a valid NAMELIST group. Read User Guide discussion on DEVC.')
REWIND(LU_INPUT)
CALL CHECKREAD('SPRK',LU_INPUT,IOS)
IF (IOS==0) CALL SHUTDOWN('ERROR: SPRK is no longer a valid NAMELIST group. Read User Guide discussion on PROP and DEVC.')
REWIND(LU_INPUT)
CALL CHECKREAD('THCP',LU_INPUT,IOS)
IF (IOS==0) CALL SHUTDOWN('ERROR: THCP is no longer a valid NAMELIST group. Read User Guide discussion on DEVC.')

REWIND(LU_INPUT)
 
END SUBROUTINE READ_DEAD

 
SUBROUTINE READ_HEAD
 
NAMELIST /HEAD/ TITLE,CHID,FYI
 
CHID    = 'output'
TITLE   = '      '
 
REWIND(LU_INPUT)
HEAD_LOOP: DO
   CALL CHECKREAD('HEAD',LU_INPUT,IOS)
   IF (IOS==1) EXIT HEAD_LOOP
   READ(LU_INPUT,HEAD,END=13,ERR=14,IOSTAT=IOS)
   14 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with HEAD line')
ENDDO HEAD_LOOP
13 REWIND(LU_INPUT)
 
CLOOP: DO I=1,39
   IF (CHID(I:I)=='.') CALL SHUTDOWN('ERROR: No periods allowed in CHID')
   IF (CHID(I:I)==' ') EXIT CLOOP
ENDDO CLOOP
 
! Define and look for a stop file

FN_STOP = TRIM(CHID)//'.stop'
INQUIRE(FILE=FN_STOP,EXIST=EX)
IF (EX) THEN
   WRITE(MESSAGE,'(A,A,A)') "ERROR: Remove the file, ",TRIM(FN_STOP),", from the current directory"
   CALL SHUTDOWN(MESSAGE)
ENDIF
 
END SUBROUTINE READ_HEAD
 
 
SUBROUTINE READ_MESH

INTEGER :: IJK(3),NM,CURRENT_MPI_PROCESS,MPI_PROCESS,RGB(3),LEVEL,N_MESH_NEW,N,II,JJ,KK,NMESHES_READ,NNN,IBAR2,JBAR2,KBAR2
INTEGER, DIMENSION(0:100) :: I_MG,J_MG,K_MG
LOGICAL :: EVACUATION, EVAC_HUMANS
REAL(EB) :: EVAC_Z_OFFSET,XB(6),XB1,XB2,XB3,XB4,XB5,XB6
CHARACTER(25) :: COLOR
CHARACTER(30) :: MULT_ID
NAMELIST /MESH/ IJK,FYI,ID,SYNCHRONIZE,EVACUATION,EVAC_HUMANS,CYLINDRICAL,XB,RGB,COLOR,EVAC_Z_OFFSET, &
                MPI_PROCESS,I_MG,J_MG,K_MG,LEVEL,MULT_ID
TYPE (MESH_TYPE), POINTER :: M
TYPE (MULTIPLIER_TYPE), POINTER :: MR
 
NMESHES = 0
NMESHES_READ = 0
 
REWIND(LU_INPUT)
COUNT_MESH_LOOP: DO
   CALL CHECKREAD('MESH',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_MESH_LOOP
   MULT_ID = 'null'
   READ(LU_INPUT,MESH,END=15,ERR=16,IOSTAT=IOS)
   N_MESH_NEW = 1
   IF (MULT_ID/='null') THEN
      DO N=1,N_MULT
         MR => MULTIPLIER(N)
         IF (MULT_ID==MR%ID) N_MESH_NEW = MR%N_COPIES
      ENDDO
   ENDIF
   NMESHES      = NMESHES + N_MESH_NEW
   NMESHES_READ = NMESHES_READ + 1
   16 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with MESH line.')
ENDDO COUNT_MESH_LOOP
15 CONTINUE


! Allocate parameters associated with the mesh.
 
ALLOCATE(MESHES(NMESHES),STAT=IZERO)
CALL ChkMemErr('READ','MESHES',IZERO)
ALLOCATE(PROCESS(NMESHES),STAT=IZERO)
CALL ChkMemErr('READ','PROCESS',IZERO)
ALLOCATE(MESH_NAME(NMESHES),STAT=IZERO)
CALL ChkMemErr('READ','MESH_NAME',IZERO)
ALLOCATE(TUSED(N_TIMERS_DIM,NMESHES),STAT=IZERO)
CALL ChkMemErr('READ','TUSED',IZERO)
ALLOCATE(SYNC_TIME_STEP(NMESHES),STAT=IZERO)
CALL ChkMemErr('READ','SYNC_TIME_STEP',IZERO)
SYNC_TIME_STEP = .TRUE.
ALLOCATE(CHANGE_TIME_STEP(NMESHES),STAT=IZERO)
CALL ChkMemErr('READ','CHANGE_TIME_STEP',IZERO)
CHANGE_TIME_STEP = .FALSE.
ALLOCATE(EVACUATION_ONLY(NMESHES),STAT=IZERO)
CALL ChkMemErr('READ','EVACUATION_ONLY',IZERO)
EVACUATION_ONLY = .FALSE.
ALLOCATE(EVACUATION_GRID(NMESHES),STAT=IZERO)
CALL ChkMemErr('READ','EVACUATION_GRID',IZERO)
EVACUATION_GRID = .FALSE.
ALLOCATE(EVACUATION_Z_OFFSET(NMESHES),STAT=IZERO)
CALL ChkMemErr('READ','EVACUATION_Z_OFFSET',IZERO)
EVACUATION_Z_OFFSET = 1.0_EB

! Read in the Mesh lines from Input file

REWIND(LU_INPUT)

NM = 0
 
MESH_LOOP: DO N=1,NMESHES_READ

   ! Set MESH defaults

   IJK(1)= 10
   IJK(2)= 10
   IJK(3)= 10
   I_MG  = -1
   J_MG  = -1
   K_MG  = -1
   TWO_D = .FALSE.
   XB(1) = 0._EB
   XB(2) = 1._EB
   XB(3) = 0._EB
   XB(4) = 1._EB
   XB(5) = 0._EB
   XB(6) = 1._EB
   RGB   = -1
   COLOR = 'null'
   CYLINDRICAL = .FALSE.
   ID = 'null'
   SYNCHRONIZE = .TRUE.
   EVACUATION  = .FALSE.
   EVAC_Z_OFFSET = 1.0_EB
   EVAC_HUMANS = .FALSE.
   MPI_PROCESS = -1
   LEVEL = 0
   MULT_ID = 'null'

   ! Read the MESH line

   CALL CHECKREAD('MESH',LU_INPUT,IOS)
   IF (IOS==1) EXIT MESH_LOOP
   READ(LU_INPUT,MESH)

   ! Parameters needed by PRESSURE_CORRECTION scheme

   IBAR2 = 1
   DO I=1,100
      IF (I_MG(I)>-1) IBAR2 = IBAR2 + 1
   ENDDO
   I_MG(0)     = 0
   I_MG(IBAR2) = IJK(1)

   JBAR2 = 1
   DO J=1,100
      IF (J_MG(J)>-1) JBAR2 = JBAR2 + 1
   ENDDO
   J_MG(0)     = 0
   J_MG(JBAR2) = IJK(2)

   KBAR2 = 1
   DO K=1,100
      IF (K_MG(K)>-1) KBAR2 = KBAR2 + 1
   ENDDO
   K_MG(0)     = 0
   K_MG(KBAR2) = IJK(3)

   ! Multiply meshes if need be

   MR => MULTIPLIER(0)
   DO NNN=1,N_MULT
      IF (MULT_ID==MULTIPLIER(NNN)%ID) MR => MULTIPLIER(NNN)
   ENDDO

   K_MULT_LOOP: DO KK=MR%K_LOWER,MR%K_UPPER
      J_MULT_LOOP: DO JJ=MR%J_LOWER,MR%J_UPPER
         I_MULT_LOOP: DO II=MR%I_LOWER,MR%I_UPPER

            IF (.NOT.MR%SEQUENTIAL) THEN
               XB1 = XB(1) + MR%DX0 + II*MR%DXB(1)
               XB2 = XB(2) + MR%DX0 + II*MR%DXB(2)
               XB3 = XB(3) + MR%DY0 + JJ*MR%DXB(3)
               XB4 = XB(4) + MR%DY0 + JJ*MR%DXB(4)
               XB5 = XB(5) + MR%DZ0 + KK*MR%DXB(5)
               XB6 = XB(6) + MR%DZ0 + KK*MR%DXB(6)
            ELSE
               XB1 = XB(1) + MR%DX0 + II*MR%DXB(1)
               XB2 = XB(2) + MR%DX0 + II*MR%DXB(2)
               XB3 = XB(3) + MR%DY0 + II*MR%DXB(3)
               XB4 = XB(4) + MR%DY0 + II*MR%DXB(4)
               XB5 = XB(5) + MR%DZ0 + II*MR%DXB(5)
               XB6 = XB(6) + MR%DZ0 + II*MR%DXB(6)
            ENDIF

            ! Increase the MESH counter by 1

            NM = NM + 1

            ! Determine which PROCESS to assign the MESH to

            IF (MPI_PROCESS>-1) THEN
               CURRENT_MPI_PROCESS = MPI_PROCESS
               IF (CURRENT_MPI_PROCESS>NUMPROCS-1) THEN
                  WRITE(MESSAGE,'(A)') 'ERROR: MPI_PROCESS greater than total number of processes'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
               RESTRICT_MESH_ASSIGNMENT = .FALSE.
            ELSE
               IF (     PARALLEL) CURRENT_MPI_PROCESS = MIN(NM-1,NUMPROCS-1)
               IF (.NOT.PARALLEL) CURRENT_MPI_PROCESS = 0
            ENDIF

            ! Fill in MESH related variables

            M => MESHES(NM)
            M%MESH_LEVEL = LEVEL
            M%IBAR = IJK(1)
            M%JBAR = IJK(2)
            M%KBAR = IJK(3)
            IBAR_MAX = MAX(IBAR_MAX,M%IBAR)
            JBAR_MAX = MAX(JBAR_MAX,M%JBAR)
            KBAR_MAX = MAX(KBAR_MAX,M%KBAR)
            M%NEWC = 2*M%IBAR*M%JBAR+2*M%IBAR*M%KBAR+2*M%JBAR*M%KBAR
            IF (.NOT.SYNCHRONIZE) SYNC_TIME_STEP(NM)  = .FALSE.
            IF (EVACUATION)  EVACUATION_ONLY(NM) = .TRUE.
            IF (EVACUATION)  SYNC_TIME_STEP(NM)  = .FALSE.
            IF (EVAC_HUMANS) EVACUATION_GRID(NM) = .TRUE.
            IF (EVACUATION)  EVACUATION_Z_OFFSET(NM) = EVAC_Z_OFFSET
            IF (EVACUATION)  M%NEWC = 2*M%IBAR*M%KBAR+2*M%JBAR*M%KBAR
            IF (M%JBAR==1) TWO_D = .TRUE.
            IF (TWO_D .AND. M%JBAR/=1) THEN
               WRITE(MESSAGE,'(A)') 'ERROR: IJK(2) must be 1 for all grids in 2D Calculation'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (EVACUATION .AND. M%KBAR/=1) THEN
               WRITE(MESSAGE,'(A)') 'ERROR: IJK(3) must be 1 for all evacuation grids'
               CALL SHUTDOWN(MESSAGE)
            ENDIF

            ! Associate the MESH with the PROCESS
         
            PROCESS(NM) = CURRENT_MPI_PROCESS
            IF (MYID==0 .AND. PARALLEL) WRITE(LU_ERR,'(A,I3,A,I3)') 'Mesh ',NM,' is assigned to Process ',PROCESS(NM)
            IF (EVACUATION_ONLY(NM) .AND. PARALLEL) EVAC_PROCESS = NUMPROCS-1

            ! Mesh boundary colors
   
            IF (ANY(RGB<0) .AND. COLOR=='null') COLOR = 'BLACK'
            IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
            ALLOCATE(M%RGB(3))
            M%RGB = RGB
   
            ! Mesh Geometry and Name
   
            WRITE(MESH_NAME(NM),'(A,I3)') 'MESH',NM
            IF (ID/='null') MESH_NAME(NM) = ID
   
            ! Parameters needed by PRESSURE_CORRECTION scheme
   
            M%IBAR2 = IBAR2
            M%JBAR2 = JBAR2
            M%KBAR2 = KBAR2

            EVAC_IF_NOT: IF (.NOT.EVACUATION_ONLY(NM)) THEN
            ALLOCATE(M%I_LO(M%IBAR2))
            ALLOCATE(M%I_HI(M%IBAR2))
            ALLOCATE(M%J_LO(M%JBAR2))
            ALLOCATE(M%J_HI(M%JBAR2))
            ALLOCATE(M%K_LO(M%KBAR2))
            ALLOCATE(M%K_HI(M%KBAR2))

            DO I=1,M%IBAR2
               M%I_LO(I) = I_MG(I-1) + 1
               M%I_HI(I) = I_MG(I)
            ENDDO

            DO J=1,M%JBAR2
               M%J_LO(J) = J_MG(J-1) + 1
               M%J_HI(J) = J_MG(J)
            ENDDO

            DO K=1,M%KBAR2
               M%K_LO(K) = K_MG(K-1) + 1
               M%K_HI(K) = K_MG(K)
            ENDDO
            ENDIF EVAC_IF_NOT

            ! Process Physical Coordinates

            IF (XB1 > XB2) THEN
               WRITE(MESSAGE,'(A,I2)') 'ERROR: XMIN > XMAX on MESH ', NM
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (XB3 > XB4) THEN
               WRITE(MESSAGE,'(A,I2)') 'ERROR: YMIN > YMAX on MESH ', NM
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (XB5 > XB6) THEN
               WRITE(MESSAGE,'(A,I2)') 'ERROR: ZMIN > ZMAX on MESH ', NM
               CALL SHUTDOWN(MESSAGE)
            ENDIF

            M%XS    = XB1
            M%XF    = XB2
            M%YS    = XB3
            M%YF    = XB4
            M%ZS    = XB5
            M%ZF    = XB6
            XS_MIN  = MIN(XS_MIN,M%XS)
            XF_MAX  = MAX(XF_MAX,M%XF)
            YS_MIN  = MIN(YS_MIN,M%YS)
            YF_MAX  = MAX(YF_MAX,M%YF)
            ZS_MIN  = MIN(ZS_MIN,M%ZS)
            ZF_MAX  = MAX(ZF_MAX,M%ZF)
            M%DXI   = (M%XF-M%XS)/REAL(M%IBAR,EB)
            M%DETA  = (M%YF-M%YS)/REAL(M%JBAR,EB)
            M%DZETA = (M%ZF-M%ZS)/REAL(M%KBAR,EB)
            M%RDXI  = 1._EB/M%DXI
            M%RDETA = 1._EB/M%DETA
            M%RDZETA= 1._EB/M%DZETA
            M%IBM1  = M%IBAR-1
            M%JBM1  = M%JBAR-1
            M%KBM1  = M%KBAR-1
            M%IBP1  = M%IBAR+1
            M%JBP1  = M%JBAR+1
            M%KBP1  = M%KBAR+1

         ENDDO I_MULT_LOOP
      ENDDO J_MULT_LOOP
   ENDDO K_MULT_LOOP

ENDDO MESH_LOOP
REWIND(LU_INPUT)
 
! Start the timing arrays
 
TUSED      = 0._EB
TUSED(1,:) = SECOND()
 
END SUBROUTINE READ_MESH



SUBROUTINE READ_TRAN
USE MATH_FUNCTIONS, ONLY : GAUSSJ
 
! Compute the polynomial transform function for the vertical coordinate
 
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: A,XX
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ND
REAL(EB) :: PC,CC,COEF,XI,ETA,ZETA
INTEGER  IEXP,IC,IDERIV,N,K,IERROR,IOS,I,MESH_NUMBER, NIPX,NIPY,NIPZ,NIPXS,NIPYS,NIPZS,NIPXF,NIPYF,NIPZF,NM
TYPE (MESH_TYPE), POINTER :: M
TYPE (TRAN_TYPE), POINTER :: T
NAMELIST /TRNX/ IDERIV,CC,PC,FYI,MESH_NUMBER
NAMELIST /TRNY/ IDERIV,CC,PC,FYI,MESH_NUMBER
NAMELIST /TRNZ/ IDERIV,CC,PC,FYI,MESH_NUMBER
 
! Scan the input file, counting the number of NAMELIST entries
 
ALLOCATE(TRANS(NMESHES))
 
MESH_LOOP: DO NM=1,NMESHES

   M => MESHES(NM)
   T => TRANS(NM)
 
   DO N=1,3
      T%NOC(N) = 0
      TRNLOOP: DO
         IF (EVACUATION_ONLY(NM)) EXIT TRNLOOP
         SELECT CASE (N)
            CASE(1)
               CALL CHECKREAD('TRNX',LU_INPUT,IOS)
               IF (IOS==1) EXIT TRNLOOP
               MESH_NUMBER = 1
               READ(LU_INPUT,NML=TRNX,END=17,ERR=18,IOSTAT=IOS)
               IF (MESH_NUMBER/=NM) CYCLE TRNLOOP
            CASE(2)
               CALL CHECKREAD('TRNY',LU_INPUT,IOS)
               IF (IOS==1) EXIT TRNLOOP
               MESH_NUMBER = 1
               READ(LU_INPUT,NML=TRNY,END=17,ERR=18,IOSTAT=IOS)
               IF (MESH_NUMBER/=NM) CYCLE TRNLOOP
            CASE(3)
               CALL CHECKREAD('TRNZ',LU_INPUT,IOS)
               IF (IOS==1) EXIT TRNLOOP
               MESH_NUMBER = 1
               READ(LU_INPUT,NML=TRNZ,END=17,ERR=18,IOSTAT=IOS)
               IF (MESH_NUMBER/=NM) CYCLE TRNLOOP
         END SELECT
         T%NOC(N) = T%NOC(N) + 1
         18 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with TRN* line')
      ENDDO TRNLOOP
      17 REWIND(LU_INPUT)
   ENDDO
 
   T%NOCMAX = MAX(T%NOC(1),T%NOC(2),T%NOC(3))
   ALLOCATE(A(T%NOCMAX+1,T%NOCMAX+1))
   ALLOCATE(XX(T%NOCMAX+1,3))
   ALLOCATE(ND(T%NOCMAX+1,3))
   ALLOCATE(T%C1(0:T%NOCMAX+1,3))
   T%C1               = 0._EB
   T%C1(1,1:3)        = 1._EB
   ALLOCATE(T%C2(0:T%NOCMAX+1,3))
   ALLOCATE(T%C3(0:T%NOCMAX+1,3))
   ALLOCATE(T%CCSTORE(T%NOCMAX,3))
   ALLOCATE(T%PCSTORE(T%NOCMAX,3))
   ALLOCATE(T%IDERIVSTORE(T%NOCMAX,3))
 
   T%ITRAN  = 0
 
   DO IC=1,3
      NLOOP:  DO N=1,T%NOC(IC)
         IDERIV = -1
         IF (IC==1) THEN
            LOOP1: DO
               CALL CHECKREAD('TRNX',LU_INPUT,IOS)
               IF (IOS==1) EXIT NLOOP
               MESH_NUMBER = 1
               READ(LU_INPUT,TRNX,END=1,ERR=2)
               IF (MESH_NUMBER==NM) EXIT LOOP1
            ENDDO LOOP1
         ENDIF
         IF (IC==2) THEN
            LOOP2: DO
               CALL CHECKREAD('TRNY',LU_INPUT,IOS)
               IF (IOS==1) EXIT NLOOP
               MESH_NUMBER = 1
               READ(LU_INPUT,TRNY,END=1,ERR=2)
               IF (MESH_NUMBER==NM) EXIT LOOP2
            ENDDO LOOP2
         ENDIF
         IF (IC==3) THEN
            LOOP3: DO
               CALL CHECKREAD('TRNZ',LU_INPUT,IOS)
               IF (IOS==1) EXIT NLOOP
               MESH_NUMBER = 1
               READ(LU_INPUT,TRNZ,END=1,ERR=2)
               IF (MESH_NUMBER==NM) EXIT LOOP3
            ENDDO LOOP3
         ENDIF
         T%CCSTORE(N,IC) = CC
         T%PCSTORE(N,IC) = PC
         T%IDERIVSTORE(N,IC) = IDERIV
         IF (IDERIV>=0) T%ITRAN(IC) = 1
         IF (IDERIV<0)  T%ITRAN(IC) = 2
      2 CONTINUE
        ENDDO NLOOP
      1 REWIND(LU_INPUT)
   ENDDO 

   ICLOOP: DO IC=1,3

      SELECT CASE (T%ITRAN(IC))

         CASE (1)  ! polynomial transformation
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
               IF (IC==1) CC = T%CCSTORE(N-1,IC)-M%XS
               IF (IC==2) CC = T%CCSTORE(N-1,IC)-M%YS
               IF (IC==3) CC = T%CCSTORE(N-1,IC)-M%ZS
               IF (IC==1 .AND. IDERIV==0) PC = T%PCSTORE(N-1,IC)-M%XS
               IF (IC==2 .AND. IDERIV==0) PC = T%PCSTORE(N-1,IC)-M%YS
               IF (IC==3 .AND. IDERIV==0) PC = T%PCSTORE(N-1,IC)-M%ZS
               IF (IC==1 .AND. IDERIV>0) PC = T%PCSTORE(N-1,IC)
               IF (IC==2 .AND. IDERIV>0) PC = T%PCSTORE(N-1,IC)
               IF (IC==3 .AND. IDERIV>0) PC = T%PCSTORE(N-1,IC)
               ND(N,IC) = IDERIV
               XX(N,IC) = CC
               T%C1(N,IC) = PC
            ENDDO NNLOOP

            DO K=1,T%NOC(IC)+1
               DO N=1,T%NOC(IC)+1
                  COEF = IFAC(K,ND(N,IC))
                  IEXP = K-ND(N,IC)
                  IF (IEXP<0) A(N,K) = 0._EB
                  IF (IEXP==0) A(N,K) = COEF
                  IF (IEXP>0) A(N,K) = COEF*XX(N,IC)**IEXP
               ENDDO
            ENDDO

            IERROR = 0
            CALL GAUSSJ(A,T%NOC(IC)+1,T%NOCMAX+1,T%C1(1:T%NOCMAX+1,IC),1,1,IERROR)
            IF (IERROR/=0)CALL SHUTDOWN('ERROR: Problem with grid transformation')

         CASE (2)  ! linear transformation

            T%C1(0,IC) = 0._EB
            T%C2(0,IC) = 0._EB
            DO N=1,T%NOC(IC)
               IF (IC==1) CC = T%CCSTORE(N,IC)-M%XS
               IF (IC==2) CC = T%CCSTORE(N,IC)-M%YS
               IF (IC==3) CC = T%CCSTORE(N,IC)-M%ZS
               IF (IC==1) PC = T%PCSTORE(N,IC)-M%XS
               IF (IC==2) PC = T%PCSTORE(N,IC)-M%YS
               IF (IC==3) PC = T%PCSTORE(N,IC)-M%ZS
               T%C1(N,IC) = CC
               T%C2(N,IC) = PC
            ENDDO

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

            DO N=1,T%NOC(IC)+1
               T%C3(N,IC) = (T%C2(N,IC)-T%C2(N-1,IC))/(T%C1(N,IC)-T%C1(N-1,IC))
            ENDDO
      END SELECT
   ENDDO ICLOOP

   DEALLOCATE(A)
   DEALLOCATE(XX)
   DEALLOCATE(ND)
 
   ! Set up grid stretching arrays
 
   ALLOCATE(M%R(0:M%IBAR),STAT=IZERO)
   CALL ChkMemErr('READ','R',IZERO)
   ALLOCATE(M%RC(0:M%IBAR+1),STAT=IZERO)
   CALL ChkMemErr('READ','RC',IZERO)
   M%RC = 1._EB
   ALLOCATE(M%RRN(0:M%IBP1),STAT=IZERO)
   CALL ChkMemErr('READ','RRN',IZERO)
   M%RRN = 1._EB
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
 
   ! Define X grid stretching terms
 
   M%DXMIN = 1000._EB
   DO I=1,M%IBAR
      XI    = (REAL(I,EB)-.5)*M%DXI
      M%HX(I) = GP(XI,1,NM)
      M%DX(I) = M%HX(I)*M%DXI
      M%DXMIN = MIN(M%DXMIN,M%DX(I))
      IF (M%HX(I)<=0._EB) THEN
         WRITE(MESSAGE,'(A,I2)')  'ERROR: x transformation not monotonic, mesh ',NM
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      M%RDX(I) = 1._EB/M%DX(I)
   ENDDO
 
   M%HX(0)    = M%HX(1)
   M%HX(M%IBP1) = M%HX(M%IBAR)
   M%DX(0)    = M%DX(1)
   M%DX(M%IBP1) = M%DX(M%IBAR)
   M%RDX(0)    = 1._EB/M%DX(1)
   M%RDX(M%IBP1) = 1._EB/M%DX(M%IBAR)
 
   DO I=0,M%IBAR
      XI     = I*M%DXI
      M%X(I) = M%XS + G(XI,1,NM)
      IF (CYLINDRICAL) THEN
         M%R(I) = M%X(I)
      ELSE
         M%R(I) = 1._EB
      ENDIF
      M%DXN(I)  = 0.5_EB*(M%DX(I)+M%DX(I+1))
      M%RDXN(I) = 1._EB/M%DXN(I)
   ENDDO
   M%X(0)      = M%XS
   M%X(M%IBAR) = M%XF
 
   DO I=1,M%IBAR
      M%XC(I) = 0.5_EB*(M%X(I)+M%X(I-1))
   ENDDO
   M%XC(0)      = M%XS - 0.5_EB*M%DX(0)
   M%XC(M%IBP1) = M%XF + 0.5_EB*M%DX(M%IBP1)
 
   IF (CYLINDRICAL) THEN  
      DO I=1,M%IBAR
         M%RRN(I) = 2._EB/(M%R(I)+M%R(I-1))
         M%RC(I)  = 0.5_EB*(M%R(I)+M%R(I-1))
      ENDDO
      M%RRN(0)    = M%RRN(1)
      M%RRN(M%IBP1) = M%RRN(M%IBAR)
   ENDIF
 
   ! Define Y grid stretching terms
 
   M%DYMIN = 1000._EB
   DO J=1,M%JBAR
      ETA   = (REAL(J,EB)-.5)*M%DETA
      M%HY(J) = GP(ETA,2,NM)
      M%DY(J) = M%HY(J)*M%DETA
      M%DYMIN = MIN(M%DYMIN,M%DY(J))
      IF (M%HY(J)<=0._EB) THEN
         WRITE(MESSAGE,'(A,I2)')  'ERROR: y transformation not monotonic, mesh ',NM
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      M%RDY(J) = 1._EB/M%DY(J)
   ENDDO
 
   M%HY(0)    = M%HY(1)
   M%HY(M%JBP1) = M%HY(M%JBAR)
   M%DY(0)    = M%DY(1)
   M%DY(M%JBP1) = M%DY(M%JBAR)
   M%RDY(0)    = 1._EB/M%DY(1)
   M%RDY(M%JBP1) = 1._EB/M%DY(M%JBAR)
 
   DO J=0,M%JBAR
      ETA     = J*M%DETA
      M%Y(J)    = M%YS + G(ETA,2,NM)
      M%DYN(J)  = 0.5_EB*(M%DY(J)+M%DY(J+1))
      M%RDYN(J) = 1._EB/M%DYN(J)
   ENDDO
 
   M%Y(0)      = M%YS
   M%Y(M%JBAR) = M%YF
 
   DO J=1,M%JBAR
      M%YC(J) = 0.5_EB*(M%Y(J)+M%Y(J-1))
   ENDDO
   M%YC(0)      = M%YS - 0.5_EB*M%DY(0)
   M%YC(M%JBP1) = M%YF + 0.5_EB*M%DY(M%JBP1)
 
   ! Define Z grid stretching terms
 
   M%DZMIN = 1000._EB
   DO K=1,M%KBAR
      ZETA  = (REAL(K,EB)-.5)*M%DZETA
      M%HZ(K) = GP(ZETA,3,NM)
      M%DZ(K) = M%HZ(K)*M%DZETA
      M%DZMIN = MIN(M%DZMIN,M%DZ(K))
      IF (M%HZ(K)<=0._EB) THEN
         WRITE(MESSAGE,'(A,I2)') 'ERROR: z transformation not monotonic, mesh ',NM
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      M%RDZ(K) = 1._EB/M%DZ(K)
   ENDDO
 
   M%HZ(0)    = M%HZ(1)
   M%HZ(M%KBP1) = M%HZ(M%KBAR)
   M%DZ(0)    = M%DZ(1)
   M%DZ(M%KBP1) = M%DZ(M%KBAR)
   M%RDZ(0)    = 1._EB/M%DZ(1)
   M%RDZ(M%KBP1) = 1._EB/M%DZ(M%KBAR)
 
   DO K=0,M%KBAR
      ZETA      = K*M%DZETA
      M%Z(K)    = M%ZS + G(ZETA,3,NM)
      M%DZN(K)  = 0.5_EB*(M%DZ(K)+M%DZ(K+1))
      M%RDZN(K) = 1._EB/M%DZN(K)
   ENDDO
 
   M%Z(0)      = M%ZS
   M%Z(M%KBAR) = M%ZF
 
   DO K=1,M%KBAR
      M%ZC(K) = 0.5_EB*(M%Z(K)+M%Z(K-1))
   ENDDO
   M%ZC(0)      = M%ZS - 0.5_EB*M%DZ(0)
   M%ZC(M%KBP1) = M%ZF + 0.5_EB*M%DZ(M%KBP1)
 
   ! Set up arrays that will return coordinate positions
 
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
 
   ALLOCATE(M%CELLSI(-NIPXS:NIPX+NIPXF),STAT=IZERO)
   CALL ChkMemErr('READ','CELLSI',IZERO)
   ALLOCATE(M%CELLSJ(-NIPYS:NIPY+NIPYF),STAT=IZERO)
   CALL ChkMemErr('READ','CELLSJ',IZERO)
   ALLOCATE(M%CELLSK(-NIPZS:NIPZ+NIPZF),STAT=IZERO)
   CALL ChkMemErr('READ','CELLSK',IZERO)
 
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
 
ENDDO MESH_LOOP
 
 
CONTAINS
 
INTEGER FUNCTION IFAC(II,N)
INTEGER II,N
IFAC = 1
DO I=II-N+1,II
   IFAC = IFAC*I
ENDDO
END FUNCTION IFAC

END SUBROUTINE READ_TRAN
 
 
SUBROUTINE READ_TIME
 
REAL(EB) :: DT,VEL_CHAR,TWFIN
INTEGER :: NM
NAMELIST /TIME/ DT,T_BEGIN,T_END,TWFIN,FYI,WALL_INCREMENT,SYNCHRONIZE, &
                EVAC_DT_FLOWFIELD,EVAC_DT_STEADY_STATE,TIME_SHRINK_FACTOR, &
                LOCK_TIME_STEP,DT_SPEC,RESTRICT_TIME_STEP
TYPE (MESH_TYPE), POINTER :: M
 
DT                   = -1._EB
EVAC_DT_FLOWFIELD    = 0.01_EB
EVAC_DT_STEADY_STATE = 0.05_EB
SYNCHRONIZE          = .FALSE.
IF (ANY(SYNC_TIME_STEP)) SYNCHRONIZE = .TRUE.
TIME_SHRINK_FACTOR   = 1._EB
T_BEGIN              = 0._EB
T_END                = -9999999._EB
TWFIN                = -99999999._EB
WALL_INCREMENT       = 2
 
REWIND(LU_INPUT)
READ_TIME_LOOP: DO
   CALL CHECKREAD('TIME',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_TIME_LOOP
   READ(LU_INPUT,TIME,END=21,ERR=22,IOSTAT=IOS)
   22 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with TIME line')
   IF (TWFIN > T_END) T_END=TWFIN
ENDDO READ_TIME_LOOP
21 REWIND(LU_INPUT)
 
IF (T_END<=T_BEGIN) SET_UP = .TRUE.
T_END = T_BEGIN + (T_END-T_BEGIN)/TIME_SHRINK_FACTOR
 
IF (.NOT.SYNCHRONIZE) SYNC_TIME_STEP = .FALSE.

! Get typical mesh cell sizes and starting time step

CHARACTERISTIC_CELL_SIZE = 1.E6_EB
 
MESH_LOOP: DO NM=1,NMESHES

   M=>MESHES(NM)

   IF (TWO_D) THEN
      M%CELL_SIZE = SQRT(M%DXMIN*M%DZMIN)
   ELSE
      M%CELL_SIZE = (M%DXMIN*M%DYMIN*M%DZMIN)**ONTH
   ENDIF

   CHARACTERISTIC_CELL_SIZE = MIN( CHARACTERISTIC_CELL_SIZE , M%CELL_SIZE )

   IF (DT>0._EB) THEN
      M%DT = DT
   ELSE
      VEL_CHAR = 0.2_EB*SQRT(10._EB*(M%ZF-M%ZS))
      M%DT = M%CELL_SIZE/VEL_CHAR
   ENDIF
   IF (EVACUATION_ONLY(NM)) THEN
      M%DT = EVAC_DT_FLOWFIELD
   ENDIF

ENDDO MESH_LOOP
 
END SUBROUTINE READ_TIME
 


SUBROUTINE READ_MULT

REAL(EB) :: DX,DY,DZ,DXB(6),DX0,DY0,DZ0
CHARACTER(30) :: ID
INTEGER :: N,I_LOWER,I_UPPER,J_LOWER,J_UPPER,K_LOWER,K_UPPER,N_LOWER,N_UPPER
TYPE(MULTIPLIER_TYPE), POINTER :: MR
NAMELIST /MULT/ FYI,ID,DX,DY,DZ,DXB,I_LOWER,I_UPPER,J_LOWER,J_UPPER,K_LOWER,K_UPPER,N_LOWER,N_UPPER,DX0,DY0,DZ0

N_MULT = 0
REWIND(LU_INPUT)
COUNT_MULT_LOOP: DO
   CALL CHECKREAD('MULT',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_MULT_LOOP
   READ(LU_INPUT,NML=MULT,END=9,ERR=10,IOSTAT=IOS)
   N_MULT = N_MULT + 1
   10 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I2)') 'ERROR: Problem with MULT no.',N_MULT
      CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_MULT_LOOP
9 REWIND(LU_INPUT)

ALLOCATE(MULTIPLIER(0:N_MULT),STAT=IZERO)
CALL ChkMemErr('READ','MULTIPLIER',IZERO)

READ_MULT_LOOP: DO N=0,N_MULT

   ID      = 'null'
   IF (N==0) ID = 'MULT DEFAULT'
   DX      = 0._EB
   DY      = 0._EB
   DZ      = 0._EB
   DX0     = 0._EB
   DY0     = 0._EB
   DZ0     = 0._EB
   DXB     = 0._EB
   I_LOWER = 0
   I_UPPER = 0
   J_LOWER = 0
   J_UPPER = 0
   K_LOWER = 0
   K_UPPER = 0
   N_LOWER = 0
   N_UPPER = 0

   IF (N>0) THEN
      CALL CHECKREAD('MULT',LU_INPUT,IOS)
      IF (IOS==1) EXIT READ_MULT_LOOP
      READ(LU_INPUT,MULT)
   ENDIF

   MR => MULTIPLIER(N)
   MR%ID      = ID
   MR%DXB     = DXB
   MR%DX0     = DX0
   MR%DY0     = DY0
   MR%DZ0     = DZ0
   IF (DX/=0._EB) MR%DXB(1:2) = DX
   IF (DY/=0._EB) MR%DXB(3:4) = DY
   IF (DZ/=0._EB) MR%DXB(5:6) = DZ

   MR%I_LOWER = I_LOWER
   MR%I_UPPER = I_UPPER
   MR%J_LOWER = J_LOWER
   MR%J_UPPER = J_UPPER
   MR%K_LOWER = K_LOWER
   MR%K_UPPER = K_UPPER
   MR%N_COPIES = (I_UPPER-I_LOWER+1)*(J_UPPER-J_LOWER+1)*(K_UPPER-K_LOWER+1)

   IF (N_LOWER/=0 .OR. N_UPPER/=0) THEN
      MR%SEQUENTIAL = .TRUE.
      MR%I_LOWER  = N_LOWER
      MR%I_UPPER  = N_UPPER
      MR%J_LOWER  = 0
      MR%J_UPPER  = 0
      MR%K_LOWER  = 0
      MR%K_UPPER  = 0
      MR%N_COPIES = (N_UPPER-N_LOWER+1)
   ENDIF

ENDDO READ_MULT_LOOP
REWIND(LU_INPUT)

END SUBROUTINE READ_MULT

 
SUBROUTINE READ_MISC
USE MATH_FUNCTIONS, ONLY: GET_RAMP_INDEX
REAL(EB) :: X_H2O_TMPA,X_H2O_40_C,C_HORIZONTAL,C_VERTICAL,MW,VISCOSITY,CONDUCTIVITY, &
            FORCE_VECTOR(3)=0._EB,SPECIFIC_HEAT,SPECIFIC_ENTHALPY,REFERENCE_TEMPERATURE
CHARACTER(30) :: RAMP_GX,RAMP_GY,RAMP_GZ
INTEGER :: NCG,N,NO,IC,JC,KC,NM,COARSE_I(10000),COARSE_MESH(10000),COARSE_J(10000),COARSE_K(10000)
TYPE (MESH_TYPE), POINTER :: M2
NAMELIST /MISC/ PR,SC,TMPA,GVEC,PRESSURE_RELAX_FACTOR,RELAXATION_FACTOR,FYI, &
                CSMAG,RAMP_GX,RAMP_GY,RAMP_GZ,BAROCLINIC, &
                LAPSE_RATE,ISOTHERMAL, &
                P_INF,SURF_DEFAULT,EVAC_SURF_DEFAULT, &
                C_FORCED,C_VERTICAL,C_HORIZONTAL,H_FIXED,RESTART,ASSUMED_GAS_TEMPERATURE, &
                BACKGROUND_SPECIES,MW,LES,DNS, &
                VISCOSITY,CONDUCTIVITY,NOISE, &
                RADIATION,GAMMA,BNDF_DEFAULT, &
                U0,V0,W0,HUMIDITY, &
                ALLOW_UNDERSIDE_DROPLETS,POROUS_FLOOR,SUPPRESSION,SUPPRESSION_SEARCH, &
                CO_PRODUCTION,SOOT_DEPOSITION,TEXTURE_ORIGIN,NSTRATA, &
                THICKEN_OBSTRUCTIONS, &
                EVAC_PRESSURE_ITERATIONS,EVAC_TIME_ITERATIONS, EVACUATION_MC_MODE, &
                PRESSURE_CORRECTION,CHECK_POISSON,STRATIFICATION,RESTART_CHID, &
                CFL_MAX,CFL_MIN,VN_MAX,VN_MIN,SOLID_PHASE_ONLY,SMOKE_ALBEDO,GROUND_LEVEL, &
                AL2O3,SHARED_FILE_SYSTEM, &
                FLUX_LIMITER,FREEZE_VELOCITY,CFL_VELOCITY_NORM,PERIODIC_TEST, &
                WIND_ONLY,TERRAIN_CASE,TURB_INIT,COMPUTE_VISCOSITY_TWICE, &
                CONSTANT_PROPERTIES,FLUXMAX,DYNSMAG,DSMAG_FREQ, &
                CHECK_KINETIC_ENERGY,PROJECTION,FISHPAK_BC,FORCE_VECTOR,DEBUG_OPENMP, &
                FDS6,SECONDARY_BREAKUP,COMBUSTION2, &
                SPECIFIC_HEAT,SPECIFIC_ENTHALPY,REFERENCE_TEMPERATURE
 
! Physical constants
 
R0      = 8314.472_EB        ! Universal Gas Constant (J/kmol/K) (NIST Physics Constants)
R1      = 1.986257E-03_EB    ! Universal Gas Constant (kcal/mol/K)
TMPA    = 20._EB             ! Ambient temperature (C)
GRAV    = 9.80665_EB         ! Acceleration of gravity (m/s**2)
GAMMA   = 1.4_EB             ! Heat capacity ratio for air
P_INF   = 101325._EB         ! Ambient pressure (Pa)
P_STP   = 101325._EB         ! Standard pressure (Pa)
TMPM    = 273.15_EB          ! Melting temperature of water (K)
SIGMA   = 5.67E-8_EB         ! Stefan-Boltzmann constant (W/m**2/K**4)
C_P_W   = 4184._EB           ! Specific Heat of Water (J/kg/K)
H_V_W   = 2259._EB*1000._EB  ! Heat of Vap of Water (J/kg)
MW_AIR  = 1._EB/(0.23_EB/32._EB+0.77_EB/28._EB) ! g/mol
HUMIDITY= -1._EB             ! Relative Humidity
RHO_SOOT= 1850._EB           ! Density of soot particle (kg/m3)
SMOKE_ALBEDO= 0.3            ! Albedo of smoke
RESTART_CHID = CHID
 
! Empirical heat transfer constants
 
C_VERTICAL   = 1.31_EB  ! Vertical free convection (Holman, Table 7-2)
C_HORIZONTAL = 1.52_EB  ! Horizontal free convection 
C_FORCED     = 0.037_EB ! Forced convection coefficient
H_FIXED                 = -1.      ! Fixed heat transfer coefficient, used for diagnostics
ASSUMED_GAS_TEMPERATURE = -1000.   ! Assumed gas temperature, used for diagnostics
 
! Background parameters
 
U0 = 0._EB
V0 = 0._EB
W0 = 0._EB
BACKGROUND_SPECIES = 'AIR'
VISCOSITY = -1._EB
CONDUCTIVITY = -1._EB
SPECIFIC_HEAT = -1._EB
SPECIFIC_ENTHALPY = -1._EB
REFERENCE_TEMPERATURE = 25._EB
MU_USER = -1._EB
K_USER  = -1._EB
MW      = 0._EB   

! Logical constants

RESTART        = .FALSE.
RADIATION      = .TRUE.
CO_PRODUCTION  = .FALSE.
CHECK_POISSON  = .FALSE.
BAROCLINIC     = .FALSE.
NOISE          = .TRUE.
ISOTHERMAL     = .FALSE.  
LES            = .TRUE.
DNS            = .FALSE.
STRATIFICATION = .TRUE.
SOLID_PHASE_ONLY = .FALSE.

TEXTURE_ORIGIN(1) = 0._EB
TEXTURE_ORIGIN(2) = 0._EB
TEXTURE_ORIGIN(3) = 0._EB
 
! EVACuation parameters
 
EVAC_PRESSURE_ITERATIONS = 50
EVAC_TIME_ITERATIONS     = 50
EVACUATION_MC_MODE       = .FALSE.
 
! LES parameters
 
CSMAG                = 0.20_EB  ! Smagorinsky constant
PR                   = -1.0_EB  ! Turbulent Prandtl number
SC                   = -1.0_EB  ! Turbulent Schmidt number
 
! Misc
 
RAMP_GX              = 'null'
RAMP_GY              = 'null'
RAMP_GZ              = 'null'
SURF_DEFAULT         = 'INERT'
EVAC_SURF_DEFAULT    = 'INERT'
GVEC(1)              = 0._EB        ! x-component of gravity 
GVEC(2)              = 0._EB        ! y-component of gravity 
GVEC(3)              = -GRAV        ! z-component of gravity 
LAPSE_RATE           = 0._EB       
RELAXATION_FACTOR    = 1.00_EB      ! Relaxation factor for no-flux
NSTRATA              = 7            ! Number bins for drop dist.
FLUXMAX              = 100._EB      ! Maximum momentum flux btw. gas and particles
THICKEN_OBSTRUCTIONS = .FALSE.
CFL_MAX              = 1.0_EB       ! Stability bounds
CFL_MIN              = 0.8_EB
VN_MAX               = 1.0_EB
VN_MIN               = 0.8_EB
 
REWIND(LU_INPUT)
MISC_LOOP: DO 
   CALL CHECKREAD('MISC',LU_INPUT,IOS)
   IF (IOS==1) EXIT MISC_LOOP
   READ(LU_INPUT,MISC,END=23,ERR=24,IOSTAT=IOS)
   24 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with MISC line')
ENDDO MISC_LOOP
23 REWIND(LU_INPUT)

! FDS 6.0 defaults (also see READ_REAC)

IF (FDS6) THEN
   CFL_VELOCITY_NORM=1
   FLUX_LIMITER=2
   DYNSMAG=.TRUE.
   BAROCLINIC=.TRUE.
   CHECK_KINETIC_ENERGY=.TRUE.
ENDIF
 
! Temperature conversions

H0    = 0.5_EB*(U0**2+V0**2+W0**2)
TMPA  = TMPA + TMPM
TMPA4 = TMPA**4
 
! Humidity (40% by default, but limited for high temps)
 
IF (HUMIDITY < 0._EB) THEN
   X_H2O_TMPA = MIN( 1._EB , EXP(-(H_V_W*MW_H2O/R0)*(1._EB/TMPA     -1._EB/373.15_EB)) )
   X_H2O_40_C =              EXP(-(H_V_W*MW_H2O/R0)*(1._EB/313.15_EB-1._EB/373.15_EB))
   HUMIDITY = 40._EB*MIN( 1._EB , X_H2O_40_C/X_H2O_TMPA )
ENDIF
 
! Miscellaneous
 
MW_BACKGROUND = MW
MU_USER(0) = VISCOSITY
K_USER(0)  = CONDUCTIVITY
CP_BACKGROUND = SPECIFIC_HEAT*1000._EB
H_BACKGROUND = SPECIFIC_ENTHALPY*1000._EB
TREF_BACKGROUND = REFERENCE_TEMPERATURE+273.15_EB
HCH    = C_HORIZONTAL
HCV    = C_VERTICAL
IF (ISOTHERMAL) RADIATION = .FALSE.
C_FORCED = C_FORCED*(1012._EB)*(1.8E-5_EB)**0.2_EB / (0.7_EB)**(2._EB/3._EB)
ASSUMED_GAS_TEMPERATURE = ASSUMED_GAS_TEMPERATURE + TMPM
TEX_ORI = TEXTURE_ORIGIN
FVEC = FORCE_VECTOR
 
! Do some set-up work for PRESSURE_CORRECTION scheme

NCGC = NMESHES

PRES_CORR_SETUP: IF (PRESSURE_CORRECTION) THEN

   SYNCHRONIZE = .TRUE.
   SYNC_TIME_STEP = .TRUE.

   NCGC = 0
   MESH_LOOP_2: DO NM=1,NMESHES
      IF(EVACUATION_ONLY(NM)) CYCLE MESH_LOOP_2
      M=>MESHES(NM)
      ALLOCATE(M%CGI(M%IBAR,M%JBAR,M%KBAR))
      ALLOCATE(M%CGI2(M%IBAR2,M%JBAR2,M%KBAR2))
      DO KC=1,M%KBAR2
         DO JC=1,M%JBAR2
            DO IC=1,M%IBAR2
               NCGC = NCGC+1
               M%CGI2(IC,JC,KC) = NCGC
               M%CGI(M%I_LO(IC):M%I_HI(IC), M%J_LO(JC):M%J_HI(JC), M%K_LO(KC):M%K_HI(KC)) = NCGC
            ENDDO
         ENDDO
      ENDDO
   ENDDO MESH_LOOP_2

   NCG = 0
   MESH_LOOP_2B: DO NM=1,NMESHES
      IF(EVACUATION_ONLY(NM)) CYCLE MESH_LOOP_2B
      M=>MESHES(NM)
      DO KC=1,M%KBAR2
         DO JC=1,M%JBAR2
            DO IC=1,M%IBAR2
               NCG = NCG+1
               COARSE_MESH(NCG) = NM
               COARSE_I(NCG) = IC
               COARSE_J(NCG) = JC
               COARSE_K(NCG) = KC
            ENDDO
         ENDDO
      ENDDO
   ENDDO MESH_LOOP_2B

   ! Generate coarse grid spacings to be used in PRESSURE_CORRECTION scheme

   ALLOCATE(DX_M(NCGC,NCGC),STAT=IZERO)
   CALL ChkMemErr('READ','DX_M',IZERO)
   ALLOCATE(DY_M(NCGC,NCGC),STAT=IZERO)
   CALL ChkMemErr('READ','DY_M',IZERO)
   ALLOCATE(DZ_M(NCGC,NCGC),STAT=IZERO)
   CALL ChkMemErr('READ','DZ_M',IZERO)
   DX_M = 0._EB
   DY_M = 0._EB
   DZ_M = 0._EB
   MESH_LOOP_3: DO NM=1,NMESHES
      IF(EVACUATION_ONLY(NM)) CYCLE MESH_LOOP_3
      M => MESHES(NM)
      DO KC=1,M%KBAR2
         DO JC=1,M%JBAR2
            DO IC=1,M%IBAR2
               N = M%CGI2(IC,JC,KC)
               DO NO=1,NCGC
                  M2 => MESHES(COARSE_MESH(NO))
                  DX_M(N,NO)=0.5_EB*(M%X(M%I_HI(IC))-M%X(M%I_LO(IC)-1) + M2%X(M2%I_HI(COARSE_I(NO)))-M2%X(M2%I_LO(COARSE_I(NO))-1) )
                  DY_M(N,NO)=0.5_EB*(M%Y(M%J_HI(JC))-M%Y(M%J_LO(JC)-1) + M2%Y(M2%J_HI(COARSE_J(NO)))-M2%Y(M2%J_LO(COARSE_J(NO))-1) )
                  DZ_M(N,NO)=0.5_EB*(M%Z(M%K_HI(KC))-M%Z(M%K_LO(KC)-1) + M2%Z(M2%K_HI(COARSE_K(NO)))-M2%Z(M2%K_LO(COARSE_K(NO))-1) )
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO MESH_LOOP_3

ENDIF PRES_CORR_SETUP
 
! Gravity ramp
 
I_RAMP_GX   = 0
I_RAMP_GY   = 0
I_RAMP_GZ   = 0
N_RAMP      = 0
IF (RAMP_GX/='null') CALL GET_RAMP_INDEX(RAMP_GX,'TIME',I_RAMP_GX)
IF (RAMP_GY/='null') CALL GET_RAMP_INDEX(RAMP_GY,'TIME',I_RAMP_GY)
IF (RAMP_GZ/='null') CALL GET_RAMP_INDEX(RAMP_GZ,'TIME',I_RAMP_GZ)
 
! Prandtl and Schmidt numbers
 
IF (DNS) THEN
!   BAROCLINIC       = .TRUE.
   LES              = .FALSE.
   IF (PR<0._EB) PR = 0.7_EB
   IF (SC<0._EB) SC = 1.0_EB
ELSE
   IF (PR<0._EB) PR = 0.5_EB
   IF (SC<0._EB) SC = 0.5_EB
ENDIF
 
RSC = 1._EB/SC
RPR = 1._EB/PR
 
! Check for a restart file
 
APPEND = .FALSE.
IF (RESTART .AND. RESTART_CHID == CHID) APPEND = .TRUE.
IF (RESTART) NOISE  = .FALSE.
 
! Min and Max values of species and temperature
 
TMPMIN = TMPM
IF (LAPSE_RATE < 0._EB) TMPMIN = MIN(TMPMIN,TMPA+LAPSE_RATE*MESHES(1)%ZF)
TMPMAX = 3000._EB
YYMIN  = 0._EB
YYMAX  = 1._EB

IF (.NOT. SUPPRESSION .AND. CO_PRODUCTION) THEN
   WRITE(MESSAGE,'(A)')  'Cannot set SUPPRESSION=.FALSE. when CO_PRODUCTION=.TRUE.'
   CALL SHUTDOWN(MESSAGE)
ENDIF 

IF (SOOT_DEPOSITION .AND. CO_PRODUCTION) THEN
   WRITE(MESSAGE,'(A)')  'Cannot set CO_PRODUCTION=.TRUE. when SOOT_DEPOSITION=.TRUE.'
   CALL SHUTDOWN(MESSAGE)
ENDIF 

IF (AL2O3) RHO_SOOT = 4000._EB

END SUBROUTINE READ_MISC

SUBROUTINE READ_DUMP

! Read parameters associated with output files
 
INTEGER :: N,I_DUM
NAMELIST /DUMP/ RENDER_FILE,SMOKE3D,SMOKE3D_COMPRESSION,SMOKE3D_QUANTITY,FLUSH_FILE_BUFFERS,MASS_FILE,STATE_FILE, &
                DT_CTRL,DT_PART,DT_MASS,DT_HRR,DT_DEVC,DT_PROF,DT_SLCF,DT_PL3D,DT_ISOF,DT_BNDF,DT_FLUSH, &
                NFRAMES,DT_RESTART,DEBUG,TIMING,CHECK_VOLUME_FLOW,COLUMN_DUMP_LIMIT,MAXIMUM_DROPLETS,WRITE_XYZ, &
                PLOT3D_QUANTITY,UL_PAN_DATA,STATUS_FILES,DEVC_COLUMN_LIMIT,CTRL_COLUMN_LIMIT,PLOT3D_SPEC_ID,SMOKE3D_SPEC_ID, &
                PLOT3D_PART_ID,CENTERED_HRR_TIME,CENTERED_DEVC_TIME,HEADERS
 
RENDER_FILE          = 'null'
NFRAMES              = 1000 
SMOKE3D              = .TRUE.
SMOKE3D_COMPRESSION  = 'RLE'
DEBUG                = .FALSE.
TIMING               = .FALSE.
FLUSH_FILE_BUFFERS   = .TRUE.
SMOKE3D_QUANTITY   = 'null'
SMOKE3D_SPEC_ID    = 'null'
PLOT3D_QUANTITY(1) = 'TEMPERATURE'
PLOT3D_QUANTITY(2) = 'U-VELOCITY'
PLOT3D_QUANTITY(3) = 'V-VELOCITY'
PLOT3D_QUANTITY(4) = 'W-VELOCITY'
PLOT3D_QUANTITY(5) = 'HRRPUV'
PLOT3D_PART_ID     = 'null'
PLOT3D_SPEC_ID     = 'null'
IF (ISOTHERMAL) PLOT3D_QUANTITY(5) = 'PRESSURE'
MAXIMUM_DROPLETS  = 500000
COLUMN_DUMP_LIMIT = .TRUE.  ! Limit csv files to 255 columns
 
DT_BNDF    = -1._EB
DT_RESTART = 1000000._EB
DT_FLUSH   = -1._EB
DT_DEVC    = -1._EB
DT_HRR     = -1._EB
DT_ISOF    = -1._EB
DT_MASS    = -1._EB
DT_PART    = -1._EB
DT_PL3D    = -1._EB
DT_PROF    = -1._EB
DT_SLCF    = -1._EB
DT_CTRL    = -1._EB
CENTERED_HRR_TIME = .FALSE.
CENTERED_DEVC_TIME = .FALSE.
 
REWIND(LU_INPUT)
DUMP_LOOP: DO 
   CALL CHECKREAD('DUMP',LU_INPUT,IOS)
   IF (IOS==1) EXIT DUMP_LOOP
   READ(LU_INPUT,DUMP,END=23,ERR=24,IOSTAT=IOS)
   24 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with DUMP line')
ENDDO DUMP_LOOP
23 REWIND(LU_INPUT)
IF (DT_BNDF < 0._EB) DT_BNDF = 2._EB * (T_END - T_BEGIN)/REAL(NFRAMES,EB) 
IF (DT_DEVC < 0._EB) DT_DEVC = (T_END - T_BEGIN)/REAL(NFRAMES,EB) 
IF (DT_HRR  < 0._EB) DT_HRR  = (T_END - T_BEGIN)/REAL(NFRAMES,EB)
IF (DT_ISOF < 0._EB) DT_ISOF = (T_END - T_BEGIN)/REAL(NFRAMES,EB) 
IF (DT_MASS < 0._EB) DT_MASS = (T_END - T_BEGIN)/REAL(NFRAMES,EB) 
IF (DT_PART < 0._EB) DT_PART = (T_END - T_BEGIN)/REAL(NFRAMES,EB) 
IF (DT_PL3D < 0._EB) DT_PL3D = (T_END - T_BEGIN)/5._EB 
IF (DT_PROF < 0._EB) DT_PROF = (T_END - T_BEGIN)/REAL(NFRAMES,EB) 
IF (DT_SLCF < 0._EB) DT_SLCF = (T_END - T_BEGIN)/REAL(NFRAMES,EB)
IF (DT_CTRL < 0._EB) DT_CTRL = (T_END - T_BEGIN)/REAL(NFRAMES,EB)
IF (DT_FLUSH< 0._EB) DT_FLUSH= (T_END - T_BEGIN)/REAL(NFRAMES,EB)

! Check Plot3D QUANTITIES

PLOOP: DO N=1,5
   CALL GET_QUANTITY_INDEX(PLOT3D_SMOKEVIEW_LABEL(N),PLOT3D_SMOKEVIEW_BAR_LABEL(N),PLOT3D_QUANTITY_INDEX(N),PLOT3D_SPEC_INDEX(N), &
                           PLOT3D_PART_INDEX(N),'PLOT3D',PLOT3D_QUANTITY(N),PLOT3D_SPEC_ID(N),PLOT3D_PART_ID(N))
ENDDO PLOOP

! Check SMOKE3D viability

IF (TWO_D .OR. SOLID_PHASE_ONLY)                           SMOKE3D = .FALSE.
IF (.NOT. MIXTURE_FRACTION .AND. SMOKE3D_QUANTITY=='null') SMOKE3D = .FALSE.
IF (      MIXTURE_FRACTION .AND. SMOKE3D_QUANTITY=='null') THEN
   SMOKE3D_QUANTITY = 'MASS FRACTION'
   SMOKE3D_SPEC_ID  = 'soot'
ENDIF
IF (SMOKE3D) THEN
   CALL GET_QUANTITY_INDEX(SMOKE3D_SMOKEVIEW_LABEL,SMOKE3D_SMOKEVIEW_BAR_LABEL,SMOKE3D_QUANTITY_INDEX,SMOKE3D_SPEC_INDEX,I_DUM, &
                           'SMOKE3D',SMOKE3D_QUANTITY,SMOKE3D_SPEC_ID,'null')
ENDIF

END SUBROUTINE READ_DUMP

 
SUBROUTINE READ_SPEC
 
USE DEVICE_VARIABLES, ONLY : PROPERTY_TYPE,PROPERTY,N_PROP
REAL(EB) :: MASS_FRACTION_0,MW,SIGMALJ,EPSILONKLJ,XVAP,VISCOSITY,CONDUCTIVITY,DIFFUSIVITY,MASS_EXTINCTION_COEFFICIENT, &
            SPECIFIC_HEAT,SPECIFIC_ENTHALPY,REFERENCE_TEMPERATURE
INTEGER  :: N_SPEC_READ,N_SPEC_EXTRA,N_MIX,N,NN,I,IPC
LOGICAL  :: ABSORBING
CHARACTER(30) :: FORMULA
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC
TYPE (PROPERTY_TYPE), POINTER :: PY
NAMELIST /SPEC/ MASS_FRACTION_0,MW,FYI,ID,SIGMALJ,EPSILONKLJ,CONDUCTIVITY,VISCOSITY,DIFFUSIVITY,ABSORBING, &
                MASS_EXTINCTION_COEFFICIENT,FORMULA,SPECIFIC_HEAT,SPECIFIC_ENTHALPY,REFERENCE_TEMPERATURE
 
! Zero out indices of various species
 
I_WATER   = 0
I_FUEL    = 0
I_PROG_F  = 0
I_PROG_CO = 0
I_CO2     = 0
I_CO      = 0
I_O2      = 0
I_SOOT    = 0
I_N2      = 0
I_H2      = 0
 
! Count SPEC lines and check for errors
 
N_SPEC_READ = 0
REWIND(LU_INPUT)
COUNT_SPEC_LOOP: DO
   CALL CHECKREAD('SPEC',LU_INPUT,IOS) 
   IF (IOS==1) EXIT COUNT_SPEC_LOOP
   ID = 'null'
   READ(LU_INPUT,NML=SPEC,END=29,ERR=30,IOSTAT=IOS)
   IF (ID=='null') THEN
      WRITE(MESSAGE,'(A,I2,A)') 'ERROR: Species',N_SPECIES+1, 'needs a name (ID=...)'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   N_SPECIES = N_SPECIES + 1
   N_SPEC_READ = N_SPEC_READ + 1   
   IF (ID=='WATER VAPOR')     I_WATER = N_SPECIES
   IF (ID=='CARBON DIOXIDE')  I_CO2   = N_SPECIES
   IF (ID=='CARBON MONOXIDE') I_CO    = N_SPECIES      
   IF (ID=='OXYGEN')          I_O2    = N_SPECIES            
   IF (ID=='SOOT')            I_SOOT  = N_SPECIES            
   IF (ID=='NITROGEN')        I_N2    = N_SPECIES        
   IF (ID=='HYDROGEN')        I_H2    = N_SPECIES        
   IF (MIXTURE_FRACTION .AND. (ID=='MIXTURE_FRACTION_1' .OR. ID=='MIXTURE_FRACTION_2' .OR. &
      (CO_PRODUCTION .AND. ID=='MIXTURE_FRACTION_3'))) N_SPECIES = N_SPECIES - 1
   30 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I2)') 'ERROR: Problem with SPECies number',N_SPECIES+1
         CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_SPEC_LOOP
29 REWIND(LU_INPUT)
N_SPEC_EXTRA = N_SPECIES 

! Add other "SPECies" to the list
 
N_MIX = 0
DO N=1,N_REACTIONS
   IF (REACTION(N)%MODE==MIXTURE_FRACTION_REACTION) THEN
      N_MIX  = N_MIX  + 1
      N_SPECIES = N_SPECIES + 1
      IF (N_MIX==1) THEN
         I_FUEL = N_SPECIES
      ELSEIF (N_MIX== 2) THEN
         IF(CO_PRODUCTION) THEN 
            I_PROG_CO = N_SPECIES
         ELSE
            I_PROG_F = N_SPECIES
         ENDIF
      ELSEIF (N_MIX==3) THEN
         I_PROG_F = N_SPECIES
      ENDIF
      YYMIN(N_MIX)  = 0.0_EB
      YYMIN(I_FUEL) = 0.0_EB
   ENDIF
ENDDO

! Record the total number of mixture fraction variables that are to be tracked

N_MIX_SPECIES = N_MIX

! Add additional species for soot and water vapor if desired by user or needed by FDS

IF (SOOT_DEPOSITION .AND. MIXTURE_FRACTION) THEN
   N_SPECIES = N_SPECIES + 1
   I_PROG_SOOT = N_SPECIES
ENDIF
 
IF (WATER_EVAPORATION .AND. I_WATER==0) THEN
   N_SPECIES  = N_SPECIES + 1
   I_WATER = N_SPECIES
ENDIF
 
! Allocate species-related arrays

ALLOCATE(SPECIES(0:N_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ','SPECIES',IZERO)
SPECIES(:)%YY0 = 0._EB
SPECIES(:)%MW  = 0._EB
SPECIES(0)%MW  = MW_BACKGROUND
SPECIES(0)%ID  = BACKGROUND_SPECIES
SPECIES(0)%SPECIFIC_HEAT = CP_BACKGROUND
SPECIES(0)%SPECIFIC_ENTHALPY = H_BACKGROUND
SPECIES(0)%REFERENCE_TEMPERATURE = TREF_BACKGROUND
EPSK           = 0._EB
SIG            = 0._EB
SPECIES(:)%ABSORBING=.FALSE.
SPECIES(:)%MODE=GAS_SPECIES

N_MIX=0
NN=N_SPEC_EXTRA
I_Z_MIN = -999
I_Z_MAX = -999
!NN = 0
DO N=1,N_REACTIONS
   IF (REACTION(N)%MODE==MIXTURE_FRACTION_REACTION) THEN
      NN    = NN + 1
      N_MIX = N_MIX + 1
      SPECIES(NN)%MODE = MIXTURE_FRACTION_SPECIES
      SPECIES(NN)%ABSORBING=.TRUE.
      WRITE(SPECIES(NN)%FORMULA,'(A,I1.1)') 'Z',N_MIX
      IF (N_REACTIONS > 1) WRITE(SPECIES(NN)%ID,'(A,I1.1)') 'MIXTURE_FRACTION_',N_MIX
      IF (N_REACTIONS ==1) WRITE(SPECIES(NN)%ID,'(A)'     ) 'MIXTURE_FRACTION'
      IF (I_Z_MIN < 0) I_Z_MIN = NN
      I_Z_MAX = NN      
   ENDIF
ENDDO

IF (SOOT_DEPOSITION .AND. MIXTURE_FRACTION) THEN
   NN = NN + 1
   N_MIX = N_MIX + 1
   SPECIES(NN)%MODE = MIXTURE_FRACTION_SPECIES
   SPECIES(NN)%ABSORBING=.TRUE.
   WRITE(SPECIES(NN)%FORMULA,'(A,I1.1)') 'Z',N_MIX
   IF (N_REACTIONS > 1) WRITE(SPECIES(NN)%ID,'(A,I1.1)') 'MIXTURE_FRACTION_',N_MIX
   I_Z_MAX = NN      
ENDIF
 
! Initialize initial WATER VAPOR mass fraction, if necessary
 
IF (I_WATER/=0) THEN
   XVAP  = MIN(1._EB,EXP(2259.E3_EB*MW_H2O/R0*(1._EB/373.15_EB-1._EB/ MIN(TMPA,373.15_EB))))
   SPECIES(I_WATER)%MODE = GAS_SPECIES
   SPECIES(I_WATER)%ID   = 'WATER VAPOR'
   SPECIES(I_WATER)%YY0 =HUMIDITY*0.01_EB*XVAP/(MW_AIR/MW_H2O+(1._EB-MW_AIR/MW_H2O)*XVAP)
ENDIF
IF (I_SOOT/=0) SPECIES(I_SOOT)%MODE=AEROSOL_SPECIES
 
! Read SPEC lines from input file
 
REWIND(LU_INPUT)
SPEC_LOOP: DO N=1,N_SPEC_READ
 
   CONDUCTIVITY    = -1._EB
   DIFFUSIVITY     = -1._EB
   EPSILONKLJ      = 0._EB
   MASS_EXTINCTION_COEFFICIENT = 8700._EB  ! m2/kg
   MASS_FRACTION_0 = -1._EB
   MW              = 0._EB
   SIGMALJ         = 0._EB
   VISCOSITY       = -1._EB
   ABSORBING       = .FALSE.
   ID              = 'null'
   FORMULA         = 'null'
   SPECIFIC_HEAT   = -1._EB
   SPECIFIC_ENTHALPY        = -1._EB
   REFERENCE_TEMPERATURE  = 25._EB
   CALL CHECKREAD('SPEC',LU_INPUT,IOS)
   READ(LU_INPUT,NML=SPEC)

   IF (MIXTURE_FRACTION) THEN
      SELECT CASE(ID)
         CASE('MIXTURE_FRACTION_1')
            NN = I_FUEL
            ABSORBING  = .TRUE.
            FORMULA    = 'Z1'
         CASE('MIXTURE_FRACTION_2')
            ABSORBING  = .TRUE.
            FORMULA = 'Z2'
            IF (CO_PRODUCTION) THEN
               NN = I_PROG_CO
            ELSE
               NN = I_PROG_F
            ENDIF
         CASE('MIXTURE_FRACTION_3')
            FORMULA = 'Z3'
            IF (CO_PRODUCTION) THEN
               NN = I_PROG_CO
               ABSORBING  = .TRUE.
            ELSE
               NN = N
            ENDIF
         CASE DEFAULT
            NN = N
      END SELECT
   ELSE
      NN = N
   ENDIF
   SS => SPECIES(NN) 
   MU_USER(NN) = VISCOSITY
   K_USER(NN)  = CONDUCTIVITY
   D_USER(NN)  = DIFFUSIVITY
   EPSK(NN)    = EPSILONKLJ
   SIG(NN)     = SIGMALJ
   SS%MW        = MW
   SS%ABSORBING = ABSORBING
   SS%SPECIFIC_HEAT = SPECIFIC_HEAT*1000._EB
   SS%SPECIFIC_ENTHALPY = SPECIFIC_ENTHALPY*1000._EB
   IF ((SPECIFIC_HEAT < 0._EB .AND. SPECIFIC_ENTHALPY > 0._EB).OR.(SPECIFIC_HEAT > 0._EB .AND. SPECIFIC_ENTHALPY < 0._EB)) THEN
      WRITE(MESSAGE,'(A,I2)') 'ERROR: SPECies number',N_SPECIES+1,&
         ' if either SPECIFIC_HEAT or SPECIFIC_ENTHALPY is specified, then both must be specified'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   SS%REFERENCE_TEMPERATURE = REFERENCE_TEMPERATURE + TMPM
   IF (MASS_FRACTION_0 < 0._EB) THEN
      IF (N/=I_WATER)  SS%YY0 = 0._EB
   ELSE
      SS%YY0 = MASS_FRACTION_0
   ENDIF
   IF (ID/='null') SS%ID = ID
   SS%MASS_EXTINCTION_COEFFICIENT = MASS_EXTINCTION_COEFFICIENT
   IF (FORMULA=='null') WRITE(FORMULA,'(I2.2)') NN
   SS%FORMULA = FORMULA
 
ENDDO SPEC_LOOP

! Check to see that SPECies needed by other routines have actually been read in

PART_CLASS_LOOP: DO IPC=1,N_PART
   PC=>PARTICLE_CLASS(IPC)
   IF (PC%SPEC_ID=='null') CYCLE PART_CLASS_LOOP
   DO I=1,N_SPECIES
      IF (SPECIES(I)%ID==PC%SPEC_ID) THEN
         PC%SPEC_INDEX = I
         CYCLE PART_CLASS_LOOP
      ENDIF
   ENDDO
   WRITE(MESSAGE,'(A,A,A)') 'ERROR: SPEC_ID ',TRIM(PC%SPEC_ID),' not found'
   CALL SHUTDOWN(MESSAGE)
ENDDO PART_CLASS_LOOP

PROP_LOOP: DO IPC=1,N_PROP
   PY=>PROPERTY(IPC)
   IF (PY%SPEC_ID=='null') CYCLE PROP_LOOP
   DO I=1,N_SPECIES
      IF (SPECIES(I)%ID==PY%SPEC_ID) THEN
         PY%SPEC_INDEX = I
         CYCLE PROP_LOOP
      ENDIF
   ENDDO
   WRITE(MESSAGE,'(A,A,A)') 'ERROR: SPEC_ID ',TRIM(PY%SPEC_ID),' not found'
   CALL SHUTDOWN(MESSAGE)
ENDDO PROP_LOOP

END SUBROUTINE READ_SPEC
 
 
SUBROUTINE READ_REAC
 
CHARACTER(30) :: FUEL,OXIDIZER
LOGICAL :: IDEAL
INTEGER :: NN,N_REAC_READ
REAL(EB) :: Y_O2_INFTY,Y_F_INLET, &
            H2_YIELD,SOOT_YIELD,CO_YIELD, Y_F_LFL,X_O2_LL,EPUMO2,BOF,&
            CRITICAL_FLAME_TEMPERATURE,HEAT_OF_COMBUSTION,NU(MAX_SPECIES),E,N_S(MAX_SPECIES),C,H,N,O,OTHER, &
            FUEL_HEAT_OF_FORMATION,MAXIMUM_VISIBILITY 
NAMELIST /REAC/ E,BOF,HEAT_OF_COMBUSTION,FYI,FUEL,OXIDIZER,EPUMO2,ID, N_S,&
                Y_O2_INFTY,Y_F_INLET, &
                H2_YIELD,SOOT_YIELD,CO_YIELD,Y_F_LFL,X_O2_LL,CRITICAL_FLAME_TEMPERATURE,NU,SOOT_H_FRACTION, &
                C,H,N,O,OTHER,MW_OTHER,IDEAL,MASS_EXTINCTION_COEFFICIENT,VISIBILITY_FACTOR,MAXIMUM_VISIBILITY, &
                EDDY_DISSIPATION,C_EDC,HRRPUA_SHEET,HRRPUV_AVERAGE,HRRPUV_AVERAGE_TIME
 
N_REACTIONS = 0
REWIND(LU_INPUT)
 
COUNT_REAC_LOOP: DO
   ID   = 'null'
   FUEL = 'null'      
   CALL CHECKREAD('REAC',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_REAC_LOOP
   READ(LU_INPUT,REAC,ERR=434,IOSTAT=IOS)
   N_REACTIONS = N_REACTIONS + 1
   IF (FUEL=='null') THEN
      MIXTURE_FRACTION = .TRUE.
   ENDIF
   IF (FUEL/='null' .AND. MIXTURE_FRACTION) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: cannot use both finite rate REAC and mixture fraction REAC'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (MIXTURE_FRACTION .AND. N_REACTIONS > 1) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: can not have more than one reaction when using mixture fraction'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   434 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with REAC ',N_REACTIONS+1
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO COUNT_REAC_LOOP
 
IF (FUEL_EVAPORATION) MIXTURE_FRACTION = .TRUE.
IF (MIXTURE_FRACTION .AND. N_REACTIONS==0) N_REACTIONS = 1

IF (N_REACTIONS==0) RETURN

IF (MIXTURE_FRACTION) THEN
   N_REAC_READ = 1
ELSE
   N_REAC_READ = N_REACTIONS
ENDIF

IF (MIXTURE_FRACTION .AND.      CO_PRODUCTION) N_REACTIONS = 3
IF (MIXTURE_FRACTION .AND. .NOT.CO_PRODUCTION) N_REACTIONS = 2
ALLOCATE(REACTION(N_REACTIONS),STAT=IZERO)
CALL ChkMemErr('READ','REACTION',IZERO)

! Read the input file looking for REAC lines
REWIND(LU_INPUT)
READ_REACTION_LOOP: DO NN=1,N_REAC_READ
   RN => REACTION(NN)
   CALL CHECKREAD('REAC',LU_INPUT,IOS) 
   CALL SET_REAC_DEFAULTS
   IF (IOS==0) READ(LU_INPUT,REAC)
   RN%BOF                = BOF
   RN%E                  = E*1000._EB
   RN%EPUMO2             = EPUMO2*1000._EB
   RN%FUEL               = FUEL        
   RN%FYI                = FYI
   RN%HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION*1000._EB
   IF (FUEL=='null') THEN
      RN%MODE = MIXTURE_FRACTION_REACTION
      ALLOCATE(RN%NU(9))
      RN%NU = 0._EB
   ELSE
      RN%MODE = FINITE_RATE_REACTION
      ALLOCATE(RN%N(MAX_SPECIES))
      RN%N = N_S
      ALLOCATE(RN%NU(MAX_SPECIES))
      RN%NU = NU
   ENDIF
   RN%NAME               = ID
   RN%OXIDIZER           = OXIDIZER
ENDDO READ_REACTION_LOOP

IF (H==0._EB) SOOT_H_FRACTION = 0._EB
MW_SOOT = MW_C * (1._EB - SOOT_H_FRACTION) + MW_H * SOOT_H_FRACTION   

SET_MIXTURE_FRACTION: IF (MIXTURE_FRACTION) THEN
   DO NN = 2,N_REACTIONS
      ALLOCATE(REACTION(NN)%NU(9))
      REACTION(NN)%NU = 0._EB
   ENDDO
   !Set reaction variable constants
   REACTION%SOOT_H_FRACTION = SOOT_H_FRACTION
   REACTION%MW_OTHER        = MW_OTHER
   REACTION%MW_FUEL         = C * MW_C + H * MW_H + O * MW_O + N * MW_N + OTHER * MW_OTHER
   REACTION%Y_O2_INFTY      = Y_O2_INFTY
   REACTION%Y_N2_INFTY      = 1._EB-Y_O2_INFTY
   REACTION%Y_O2_LL         = X_O2_LL*MW_O2/ (X_O2_LL*MW_O2+(1._EB-X_O2_LL)*MW_N2)
   REACTION%Y_F_LFL         = Y_F_LFL
   REACTION%CRIT_FLAME_TMP  = CRITICAL_FLAME_TEMPERATURE + TMPM
   REACTION%Y_F_INLET       = Y_F_INLET
   REACTION%Y_N2_INLET      = 1._EB-Y_F_INLET   
   REACTION%MODE            = MIXTURE_FRACTION_REACTION
   REACTION%NAME            = ID
   MW_AVG_C = REACTION(1)%Y_O2_INFTY/MW_O2 + REACTION(1)%Y_N2_INFTY/MW_N2

   SET_THREE_PARAMETER: IF (CO_PRODUCTION) THEN !Set reaction variables for three parameter mixture fraction 
!Set reaction variables for complete reaction
      RN => REACTION(2) 
      RN%IDEAL      = IDEAL
      IF (RN%IDEAL) THEN
!Compute fuel heat of formation
         RN%NU(O2_INDEX) = C + 0.5_EB * H - 0.5_EB * O
         IF (HEAT_OF_COMBUSTION < 0._EB) HEAT_OF_COMBUSTION = EPUMO2*(MW_O2*RN%NU(O2_INDEX))/RN%MW_FUEL
         FUEL_HEAT_OF_FORMATION =  HEAT_OF_COMBUSTION * RN%MW_FUEL - &
                                   (C * CO2_HEAT_OF_FORMATION + 0.5_EB * H * H2O_HEAT_OF_FORMATION)
      ENDIF
      RN%CO_YIELD   = CO_YIELD
      RN%SOOT_YIELD = SOOT_YIELD
      RN%H2_YIELD   = H2_YIELD
      RN%NU(H2_INDEX)      = H2_YIELD * RN%MW_FUEL / MW_H2
      RN%NU(SOOT_INDEX)    = SOOT_YIELD * RN%MW_FUEL / MW_SOOT
      RN%NU(CO_INDEX)      = CO_YIELD * RN%MW_FUEL / MW_CO
      RN%NU(CO2_INDEX)     = C - RN%NU(CO_INDEX) - RN%NU(SOOT_INDEX) * (1._EB - SOOT_H_FRACTION)
      IF (RN%NU(CO2_INDEX) < 0._EB) THEN
         WRITE(MESSAGE,'(A)') 'Values for SOOT_YIELD, CO_YIELD, and SOOT_H_FRACTION result in negative CO2 yield'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      RN%NU(H2O_INDEX)     = 0.5_EB * H - RN%NU(H2_INDEX) - 0.5_EB * RN%NU(SOOT_INDEX) * SOOT_H_FRACTION
      IF (RN%NU(H2O_INDEX) < 0._EB) THEN
         WRITE(MESSAGE,'(A)') 'Values for SOOT_YIELD, H2_YIELD, and SOOT_H_FRACTION result in negative H2O yield'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      RN%NU(O2_INDEX)      = RN%NU(CO2_INDEX) + 0.5_EB * (RN%NU(CO_INDEX) + RN%NU(H2O_INDEX) - O)
      RN%NU(N2_INDEX)      = 0.5_EB * N
      RN%NU(OTHER_INDEX)   = OTHER
      RN%BOF        = 2.53E12_EB
      RN%E          = 199547._EB*1000._EB
      RN%O2_F_RATIO    = (MW_O2*RN%NU(O2_INDEX))/RN%MW_FUEL
      IF (.NOT. RN%IDEAL) THEN
         IF(HEAT_OF_COMBUSTION < 0._EB) HEAT_OF_COMBUSTION = EPUMO2*RN%O2_F_RATIO
         RN%HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION*1000._EB
      ELSE
! Correct heat of combustion for minor products of combustion
         RN%HEAT_OF_COMBUSTION = (FUEL_HEAT_OF_FORMATION + RN%NU(CO2_INDEX) * CO2_HEAT_OF_FORMATION + &
                                                         RN%NU(CO_INDEX)  * CO_HEAT_OF_FORMATION + &
                                                         RN%NU(H2O_INDEX) * H2O_HEAT_OF_FORMATION) * 1000._EB /RN%MW_FUEL
      ENDIF
      
! Set reaction variables for incomplete reaction

      RN => REACTION(1)
      RN%IDEAL              = .TRUE.
      RN%CO_YIELD           = 0._EB
      RN%H2_YIELD           = H2_YIELD
      RN%SOOT_YIELD         = SOOT_YIELD
      RN%NU(H2_INDEX)              = H2_YIELD * RN%MW_FUEL / MW_H2
      RN%NU(SOOT_INDEX)            = SOOT_YIELD * RN%MW_FUEL / MW_SOOT
      RN%NU(CO_INDEX)              = C - RN%NU(SOOT_INDEX) * (1._EB - SOOT_H_FRACTION)
      RN%NU(CO2_INDEX)             = 0._EB
      IF (RN%NU(CO2_INDEX) < 0._EB) THEN
         WRITE(MESSAGE,'(A)') 'Values for SOOT_YIELD, CO_YIELD, and SOOT_H_FRACTION result in negative CO2 yield'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      RN%NU(H2O_INDEX)             = 0.5_EB * H - RN%NU(H2_INDEX) - 0.5_EB * RN%NU(SOOT_INDEX) * SOOT_H_FRACTION
      IF (RN%NU(H2O_INDEX) < 0._EB) THEN
         WRITE(MESSAGE,'(A)') 'Values for SOOT_YIELD, H2_YIELD, and SOOT_H_FRACTION result in negative H2O yield'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      RN%NU(O2_INDEX)              = RN%NU(CO2_INDEX) + 0.5_EB * (RN%NU(CO_INDEX) + RN%NU(H2O_INDEX) - O)
      RN%NU(N2_INDEX)              = 0.5_EB * N
      RN%NU(OTHER_INDEX)           = OTHER
      RN%O2_F_RATIO    = (MW_O2*RN%NU(O2_INDEX))/RN%MW_FUEL      
      RN%HEAT_OF_COMBUSTION = REACTION(2)%HEAT_OF_COMBUSTION - 1.E3_EB * REACTION(2)%NU(CO2_INDEX) /RN%MW_FUEL * &
                             (CO2_HEAT_OF_FORMATION - CO_HEAT_OF_FORMATION)
      NN = 3

   ELSE SET_THREE_PARAMETER !Set reaction variables for two parameter mixture fraction 
!Set reaction variables for complete reaction
      RN => REACTION(1) 
      RN%IDEAL      = IDEAL
      IF (RN%IDEAL) THEN
!Compute fuel heat of formation
         RN%NU(O2_INDEX) = C + 0.5_EB * H - 0.5_EB * O
         IF (HEAT_OF_COMBUSTION < 0._EB) HEAT_OF_COMBUSTION = EPUMO2*(MW_O2*RN%NU(O2_INDEX))/RN%MW_FUEL
         FUEL_HEAT_OF_FORMATION =  HEAT_OF_COMBUSTION * RN%MW_FUEL - &
                                   (C * CO2_HEAT_OF_FORMATION + 0.5_EB * H * H2O_HEAT_OF_FORMATION)
      ENDIF
      RN%CO_YIELD   = CO_YIELD
      RN%SOOT_YIELD = SOOT_YIELD
      RN%H2_YIELD   = H2_YIELD
      RN%NU(H2_INDEX)      = H2_YIELD * RN%MW_FUEL / MW_H2
      RN%NU(SOOT_INDEX)    = SOOT_YIELD * RN%MW_FUEL / MW_SOOT
      RN%NU(CO_INDEX)      = CO_YIELD * RN%MW_FUEL / MW_CO
      RN%NU(CO2_INDEX)     = C - RN%NU(CO_INDEX) - RN%NU(SOOT_INDEX) * (1._EB - SOOT_H_FRACTION)
      RN%NU(H2O_INDEX)     = 0.5_EB * H - RN%NU(H2_INDEX) - 0.5_EB * RN%NU(SOOT_INDEX) * SOOT_H_FRACTION
      RN%NU(O2_INDEX)      = RN%NU(CO2_INDEX) + 0.5_EB * (RN%NU(CO_INDEX) + RN%NU(H2O_INDEX) - O)
      RN%NU(N2_INDEX)      = 0.5_EB * N
      RN%NU(OTHER_INDEX)   = OTHER
      RN%O2_F_RATIO    = (MW_O2*RN%NU(O2_INDEX))/RN%MW_FUEL      
      IF (.NOT. RN%IDEAL) THEN
         IF(HEAT_OF_COMBUSTION < 0._EB) HEAT_OF_COMBUSTION = EPUMO2*RN%O2_F_RATIO
         RN%HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION*1000._EB
      ELSE
!Correct heat of combustion for minor products of combustion
         RN%HEAT_OF_COMBUSTION = (FUEL_HEAT_OF_FORMATION + RN%NU(CO2_INDEX) * CO2_HEAT_OF_FORMATION + &
                                                         RN%NU(CO_INDEX)  * CO_HEAT_OF_FORMATION + &
                                                         RN%NU(H2O_INDEX) * H2O_HEAT_OF_FORMATION) * 1000._EB /RN%MW_FUEL
      ENDIF
      NN = 2

   ENDIF SET_THREE_PARAMETER
   
   !Set reaction variables for extinction reaction   

   RN => REACTION(NN)
   RN%CO_YIELD           = 0._EB
   RN%NU(CO_INDEX)              = 0._EB
   RN%NU(CO2_INDEX)             = 0._EB
   RN%NU(H2_INDEX)              = 0._EB
   RN%NU(H2O_INDEX)             = 0._EB
   RN%NU(O2_INDEX)              = 0._EB
   RN%NU(N2_INDEX)              = 0._EB   
   RN%NU(OTHER_INDEX)           = 0._EB      
   RN%NU(SOOT_INDEX)            = 0._EB
   RN%SOOT_YIELD         = 0._EB
   RN%H2_YIELD           = 0._EB   
   RN%HEAT_OF_COMBUSTION = 0._EB   
   RN%EPUMO2             = 0._EB      
   RN%IDEAL              = .TRUE.
ENDIF SET_MIXTURE_FRACTION

! Set the lower limit of the extinction coefficient

EC_LL = VISIBILITY_FACTOR/MAXIMUM_VISIBILITY

! Change units of combustion quantities

IF (FDS6) HRRPUV_AVERAGE = 1.E10_EB ! effectively remove this bound
HRRPUV_AVERAGE = HRRPUV_AVERAGE*1000._EB   ! W/m3
HRRPUA_SHEET   = HRRPUA_SHEET*  1000._EB   ! W/m2

CONTAINS
 

SUBROUTINE SET_REAC_DEFAULTS
 
BOF                         = 0._EB       ! cm**3/mol-s
CO_YIELD                    = 0._EB
CRITICAL_FLAME_TEMPERATURE  = 1427._EB    ! C
E                           = 0._EB       ! kJ/kmol
EPUMO2                      = 13100._EB   ! kJ/kg
FUEL                        = 'null'
FYI                         = 'null'
H2_YIELD                    = 0._EB
HEAT_OF_COMBUSTION          = -1._EB
ID                          = 'null'
Y_F_INLET                   = 1._EB
N_S                         = -999._EB
NU                          = 0._EB
OXIDIZER                    = 'null'
IF (LES) SOOT_YIELD         = 0.01_EB
IF (DNS) SOOT_YIELD         = 0.0_EB
SOOT_H_FRACTION             = 0.1_EB
X_O2_LL                     = 0.15_EB    ! %
Y_F_LFL                     = 0.0_EB 
Y_F_INLET                   = 1._EB
Y_O2_INFTY                  = 0.23_EB
IDEAL                       = .FALSE.
C                           = 3._EB
H                           = 8._EB
O                           = 0._EB
N                           = 0._EB
OTHER                       = 0._EB
MW_OTHER                    = 28._EB        ! Nitrogen MW
MASS_EXTINCTION_COEFFICIENT = 8700._EB     ! m2/kg
MAXIMUM_VISIBILITY          = 30._EB       ! m
VISIBILITY_FACTOR           = 3._EB
 
END SUBROUTINE SET_REAC_DEFAULTS
 
END SUBROUTINE READ_REAC
 
 
SUBROUTINE PROC_SPEC
USE PHYSICAL_FUNCTIONS, ONLY : GET_MOLECULAR_WEIGHT, GET_SPECIFIC_GAS_CONSTANT, JANAF_TABLE
REAL(EB) :: EPSIJ,SIGMA2,AMW,OMEGA,TSTAR,D_N0,EPSK_N,SIG_N,MW_N,D_TMP(N_STATE_SPECIES),MW_Z2,MW_Z3
REAL(EB), ALLOCATABLE, DIMENSION(:) :: MU_TMP,CP_TMP,K_TMP,H_TMP
REAL(EB) :: YYY(1:N_SPECIES),FUEL_FRAC = 1._EB
INTEGER :: IPC,NN,N
CHARACTER(30) :: STATE_SPECIES(N_STATE_SPECIES),FORMULA
TYPE(PARTICLE_CLASS_TYPE), POINTER :: PC
LOGICAL :: ABSORBING

! Compute state relations for mixture fraction model

IF (MIXTURE_FRACTION) THEN

   REACTION_LOOP: DO N=1,N_REACTIONS
      RN => REACTION(N)
      RN%Y_N2_INFTY = 1._EB-RN%Y_O2_INFTY  ! Assumes that air is made up of only oxygen and nitrogen
      RN%Y_N2_INLET = 1._EB-RN%Y_F_INLET   ! Assumes that the fuel is only diluted by nitrogen
      RN%Z_F_CONS   = RN%NU(O2_INDEX) * MW_O2/RN%MW_FUEL * (1._EB + RN%Y_N2_INFTY/RN%Y_O2_INFTY)
      RN%Z_F        = RN%Y_O2_INFTY/(RN%Y_O2_INFTY+RN%Y_F_INLET*RN%O2_F_RATIO)
   ENDDO REACTION_LOOP
 
   IF (CO_PRODUCTION) THEN
      YYMAX(I_PROG_CO) = REACTION(1)%Z_F
      RN => REACTION(1)
   ENDIF
   YYMAX(I_PROG_F) = REACTION(1)%Z_F

   IF (SOOT_DEPOSITION) THEN
      YYMAX(I_PROG_F)    = REACTION(1)%Z_F * (1._EB - REACTION(1)%SOOT_YIELD)
      YYMAX(I_PROG_SOOT) = REACTION(1)%Z_F * REACTION(1)%SOOT_YIELD      
   ENDIF

   IF (CO_PRODUCTION) THEN
      REACTION(3)%Z_F=1._EB   
      RN => REACTION(2)
   ELSE
      REACTION(2)%Z_F=1._EB
      RN => REACTION(1)   
   ENDIF
       
   ! For mixture fraction model, adjust the evaporation rate of fuel droplets to account for difference in HoC.

   DO IPC=1,N_PART
      PC=>PARTICLE_CLASS(IPC)
      IF (PC%HEAT_OF_COMBUSTION > 0._EB) PC%ADJUST_EVAPORATION = PC%HEAT_OF_COMBUSTION/RN%HEAT_OF_COMBUSTION
   ENDDO

ENDIF
 
! Get gas properties
 
N_Y_ARRAY = 0
IF (MIXTURE_FRACTION) N_Y_ARRAY = N_STATE_SPECIES
CALL GAS_PROPS(SPECIES(0)%ID,SIG(0),EPSK(0),SPECIES(0)%MW,ABSORBING,SPECIES(0)%FORMULA)
DO N=1,N_SPECIES
   ABSORBING = .FALSE.
   CALL GAS_PROPS(SPECIES(N)%ID,SIG(N),EPSK(N),SPECIES(N)%MW,ABSORBING,SPECIES(N)%FORMULA)
   IF (ABSORBING) SPECIES(N)%ABSORBING = .TRUE.
   IF (SPECIES(N)%MODE /= MIXTURE_FRACTION_SPECIES) THEN
      IF (.NOT. MIXTURE_FRACTION) THEN
         N_Y_ARRAY = N_Y_ARRAY + 1         
         SPECIES(N)%INDEX = N_Y_ARRAY
      ELSE
         IF (N > 0 .AND. N /= I_SOOT .AND. N/=I_N2 .AND. N/=I_O2 .AND. N/=I_WATER .AND. N/=I_H2 .AND. N/=I_CO .AND. N/=I_CO2) THEN
             N_Y_ARRAY = N_Y_ARRAY + 1         
             SPECIES(N)%INDEX = N_Y_ARRAY
         ELSE
            IF (N==I_SOOT) THEN
                  SPECIES(N)%INDEX = SOOT_INDEX
            ELSEIF (N==I_N2) THEN
                  SPECIES(N)%INDEX = N2_INDEX
            ELSEIF (N==I_O2) THEN
                  SPECIES(N)%INDEX = O2_INDEX
            ELSEIF (N==I_WATER) THEN
                  SPECIES(N)%INDEX = H2O_INDEX
            ELSEIF (N==I_H2) THEN
                  SPECIES(N)%INDEX = H2_INDEX
            ELSEIF (N==I_CO) THEN
                  SPECIES(N)%INDEX = CO_INDEX
            ELSEIF (N==I_CO2) THEN
                  SPECIES(N)%INDEX = CO2_INDEX
            ENDIF
         ENDIF
      ENDIF 
   ELSE
      SPECIES(N)%INDEX = N_STATE_SPECIES + 1     
   ENDIF
ENDDO

! If this is not a mixture fraction calculation, add an index to the species array for background gas species

IF (.NOT. MIXTURE_FRACTION) THEN 
   N_Y_ARRAY=N_Y_ARRAY+1
   SPECIES(0)%INDEX = N_Y_ARRAY
ENDIF

! Compute Y to Y Arrays
 
ALLOCATE(Y2Y(N_Y_ARRAY,N_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ','Y2Y',IZERO) 
Y2Y = 0._EB 
ALLOCATE(Y2Y_C(N_Y_ARRAY),STAT=IZERO)
CALL ChkMemErr('READ','Y2Y_C',IZERO) 
Y2Y_C = 0._EB 
ALLOCATE(Y_MF_SUM_Y(N_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ','Y_MF_SUM_Y',IZERO) 
Y_MF_SUM_Y = 0._EB 

IF (MIXTURE_FRACTION) THEN
   Y2Y_C(N2_INDEX)   = REACTION(1)%Y_N2_INFTY   
   Y2Y_C(O2_INDEX)   = REACTION(1)%Y_O2_INFTY
   Y2Y(FUEL_INDEX,I_Z_MIN) = REACTION(1)%Y_F_INLET
   Y2Y(N2_INDEX,I_Z_MIN)     = REACTION(1)%Y_N2_INLET - REACTION(1)%Y_N2_INFTY  
   Y2Y(O2_INDEX,I_Z_MIN)     =                        - REACTION(1)%Y_O2_INFTY
   IF (I_Z_MIN>=2) THEN
      Y2Y(N2_INDEX,1:I_Z_MIN-1) =                        - REACTION(1)%Y_N2_INFTY     
      Y2Y(O2_INDEX,1:I_Z_MIN-1) =                        - REACTION(1)%Y_O2_INFTY   
   ENDIF
   IF (SOOT_DEPOSITION) FUEL_FRAC = 1._EB/(1._EB - REACTION(1)%SOOT_YIELD)
   DO N = 1,N_REACTIONS-1
      Y2Y(N2_INDEX,I_Z_MIN+N)    = REACTION(N)%Y_F_INLET / REACTION(N)%MW_FUEL * MW_N2 * REACTION(N)%NU(N2_INDEX)*FUEL_FRAC &
                                    -REACTION(N)%Y_N2_INFTY
      Y2Y(H2O_INDEX,I_Z_MIN+N)   = REACTION(N)%Y_F_INLET / REACTION(N)%MW_FUEL * MW_H2O * REACTION(N)%NU(H2O_INDEX)*FUEL_FRAC
      Y2Y(H2_INDEX,I_Z_MIN+N)    = REACTION(N)%Y_F_INLET / REACTION(N)%MW_FUEL * MW_H2 * REACTION(N)%NU(H2_INDEX)*FUEL_FRAC
      Y2Y(SOOT_INDEX,I_Z_MIN+N)  = REACTION(N)%Y_F_INLET / REACTION(N)%MW_FUEL * MW_SOOT * REACTION(N)%NU(SOOT_INDEX)
      Y2Y(OTHER_INDEX,I_Z_MIN+N) = REACTION(N)%Y_F_INLET / REACTION(N)%MW_FUEL * MW_OTHER * REACTION(N)%NU(OTHER_INDEX)*FUEL_FRAC
   ENDDO    
   IF (CO_PRODUCTION) THEN
      Y2Y(O2_INDEX,I_Z_MIN+1) = -REACTION(1)%Y_F_INLET / REACTION(1)%MW_FUEL * MW_O2 * REACTION(1)%NU(O2_INDEX) &
                                -REACTION(1)%Y_O2_INFTY
      Y2Y(O2_INDEX,I_Z_MIN+2) = -REACTION(2)%Y_F_INLET / REACTION(2)%MW_FUEL * MW_O2 * REACTION(2)%NU(O2_INDEX) &
                                -REACTION(2)%Y_O2_INFTY &
                                -REACTION(2)%Y_F_INLET / REACTION(2)%MW_FUEL * MW_O2 * REACTION(2)%NU(CO_INDEX) * 0.5_EB
      Y2Y(CO_INDEX,I_Z_MIN+1)  = REACTION(1)%Y_F_INLET / REACTION(1)%MW_FUEL * MW_CO  * REACTION(1)%NU(CO_INDEX)
      Y2Y(CO2_INDEX,I_Z_MIN+2) = REACTION(2)%Y_F_INLET / REACTION(2)%MW_FUEL * MW_CO2 * &
                                 (REACTION(2)%NU(CO2_INDEX)+REACTION(2)%NU(CO_INDEX))
   ELSE
      Y2Y(O2_INDEX,I_Z_MIN+1)  = -REACTION(1)%Y_F_INLET / REACTION(1)%MW_FUEL * MW_O2 * REACTION(1)%NU(O2_INDEX)*FUEL_FRAC &
                                 -REACTION(1)%Y_O2_INFTY
      Y2Y(CO_INDEX,I_Z_MIN+1)  = REACTION(1)%Y_F_INLET / REACTION(1)%MW_FUEL * MW_CO  * REACTION(1)%NU(CO_INDEX)*FUEL_FRAC
      Y2Y(CO2_INDEX,I_Z_MIN+1) = REACTION(1)%Y_F_INLET / REACTION(1)%MW_FUEL * MW_CO2 * REACTION(1)%NU(CO2_INDEX)*FUEL_FRAC      
   ENDIF
   IF (SOOT_DEPOSITION) THEN
      Y2Y(SOOT_INDEX,:)       = 0._EB
      Y2Y(O2_INDEX,I_Z_MAX)   = -REACTION(1)%Y_O2_INFTY
      Y2Y(N2_INDEX,I_Z_MAX)   = -REACTION(1)%Y_N2_INFTY 
      Y2Y(SOOT_INDEX,I_Z_MAX) = REACTION(1)%Y_F_INLET  
   ENDIF
   DO N=1,N_SPECIES
      IF (SPECIES(N)%MODE /= MIXTURE_FRACTION_SPECIES) Y2Y(SPECIES(N)%INDEX,N) = 1._EB - Y2Y_C(SPECIES(N)%INDEX)
      Y_MF_SUM_Y(N) = SUM(Y2Y(:,N))
   ENDDO
ELSE
   Y2Y_C(N_Y_ARRAY) = 1._EB
   DO N=1,N_SPECIES
      Y2Y(SPECIES(N)%INDEX,N) = 1._EB
      Y_MF_SUM_Y(N) = 1._EB
   ENDDO
   Y2Y(N_Y_ARRAY,:) = -1._EB
ENDIF
Y_MF_SUM_Y_C = SUM(Y2Y_C)

ALLOCATE(MW_AVG_Y(N_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ','MW_AVG_Y',IZERO) 
MW_AVG_Y = 0._EB
MW_AVG_Y_C = 0._EB

IF (MIXTURE_FRACTION) THEN
   MW_AVG_Y_C = REACTION(1)%Y_O2_INFTY/MW_O2 + REACTION(1)%Y_N2_INFTY/MW_N2
   
   MW_AVG_Y(I_Z_MIN) = Y2Y(FUEL_INDEX,I_Z_MIN)/RN%MW_FUEL + Y2Y(N2_INDEX,I_Z_MIN)/MW_N2  + Y2Y(O2_INDEX,I_Z_MIN)/MW_O2 + &
                       Y2Y(CO_INDEX,I_Z_MIN)/MW_CO + Y2Y(CO2_INDEX,I_Z_MIN)/MW_CO2 + Y2Y(SOOT_INDEX,I_Z_MIN)/MW_SOOT + &
                       Y2Y(H2_INDEX,I_Z_MIN)/MW_H2 + Y2Y(H2O_INDEX,I_Z_MIN)/MW_H2O + Y2Y(OTHER_INDEX,I_Z_MIN)/MW_OTHER

   MW_AVG_Y(I_Z_MIN+1) = Y2Y(FUEL_INDEX,I_Z_MIN+1)/RN%MW_FUEL + Y2Y(N2_INDEX,I_Z_MIN+1)/MW_N2 + Y2Y(O2_INDEX,I_Z_MIN+1)/MW_O2 + &
               Y2Y(CO_INDEX,I_Z_MIN+1)/MW_CO + Y2Y(CO2_INDEX,I_Z_MIN+1)/MW_CO2 + Y2Y(SOOT_INDEX,I_Z_MIN+1)/MW_SOOT + &
               Y2Y(H2_INDEX,I_Z_MIN+1)/MW_H2 + Y2Y(H2O_INDEX,I_Z_MIN+1)/MW_H2O + Y2Y(OTHER_INDEX,I_Z_MIN+1)/MW_OTHER
   IF (CO_PRODUCTION) MW_AVG_Y(I_Z_MIN+2) = &
                  Y2Y(FUEL_INDEX,I_Z_MIN+2)/RN%MW_FUEL + Y2Y(N2_INDEX,I_Z_MIN+2)/MW_N2 +  Y2Y(O2_INDEX,I_Z_MIN+2)/MW_O2 + &
                  Y2Y(CO_INDEX,I_Z_MIN+2)/MW_CO + Y2Y(CO2_INDEX,I_Z_MIN+2)/MW_CO2 + Y2Y(SOOT_INDEX,I_Z_MIN+2)/MW_SOOT + &
                  Y2Y(H2_INDEX,I_Z_MIN+2)/MW_H2 + Y2Y(H2O_INDEX,I_Z_MIN+2)/MW_H2O + Y2Y(OTHER_INDEX,I_Z_MIN+2)/MW_OTHER
   IF (SOOT_DEPOSITION) MW_AVG_Y(I_Z_MAX) = &
                  Y2Y(FUEL_INDEX,I_Z_MAX)/RN%MW_FUEL + Y2Y(N2_INDEX,I_Z_MAX)/MW_N2 +  Y2Y(O2_INDEX,I_Z_MAX)/MW_O2 + &
                  Y2Y(CO_INDEX,I_Z_MAX)/MW_CO + Y2Y(CO2_INDEX,I_Z_MAX)/MW_CO2 + Y2Y(SOOT_INDEX,I_Z_MAX)/MW_SOOT + &
                  Y2Y(H2_INDEX,I_Z_MAX)/MW_H2 + Y2Y(H2O_INDEX,I_Z_MAX)/MW_H2O + Y2Y(OTHER_INDEX,I_Z_MAX)/MW_OTHER
   DO N=1,N_SPECIES
      IF (SPECIES(N)%MODE /= MIXTURE_FRACTION_SPECIES) MW_AVG_Y(N) = 1._EB/SPECIES(N)%MW - MW_AVG_Y_C
   ENDDO
ELSE
   MW_AVG_Y_C = 1._EB/SPECIES(0)%MW
   DO N=1,N_SPECIES
      MW_AVG_Y(N) = Y2Y(SPECIES(N)%INDEX,N)/SPECIES(N)%MW-MW_AVG_Y_C
   ENDDO
ENDIF  

! Compute the initial value of R0/MW_AVG

SS0 => SPECIES(0)
 
SS0%YY0 = 1._EB
DO N=1,N_SPECIES
   SS => SPECIES(N)
   IF (SS%MODE==GAS_SPECIES) SS0%YY0 = SS0%YY0 - SS%YY0
ENDDO
 
IF (MIXTURE_FRACTION) THEN
   RCON_MF(FUEL_INDEX) = R0/REACTION(1)%MW_FUEL
   RCON_MF(O2_INDEX) = R0/MW_O2
   RCON_MF(N2_INDEX) = R0/MW_N2
   RCON_MF(H2O_INDEX) = R0/MW_H2O
   RCON_MF(CO2_INDEX) = R0/MW_CO2
   RCON_MF(CO_INDEX) = R0/MW_CO
   RCON_MF(H2_INDEX) = R0/MW_H2
   RCON_MF(SOOT_INDEX) = R0/MW_SOOT
   RCON_MF(OTHER_INDEX) = R0/REACTION(1)%MW_OTHER
   SS0%MW = 1._EB/(REACTION(1)%Y_O2_INFTY/MW_O2+REACTION(1)%Y_N2_INFTY/MW_N2)
ENDIF

MW_MIN = SS0%MW
MW_MAX = SS0%MW
 
YYY(1:N_SPECIES) = SPECIES(1:N_SPECIES)%YY0
CALL GET_SPECIFIC_GAS_CONSTANT(YYY,RSUM0)
MW_MIN = MINVAL(SPECIES%MW)
MW_MAX = MAXVAL(SPECIES%MW)
SPECIES%RCON = R0/SPECIES%MW

IF (MIXTURE_FRACTION) THEN
   IF (CO_PRODUCTION) THEN
      IF (SOOT_DEPOSITION) THEN
         YYY(I_Z_MIN) = 0._EB
         YYY(I_Z_MIN+1) = REACTION(1)%Z_F * (1._EB - REACTION(1)%SOOT_YIELD)
         YYY(I_Z_MIN+2) = 0._EB
         YYY(I_Z_MAX)   = REACTION(1)%Z_F * REACTION(1)%SOOT_YIELD
         CALL GET_MOLECULAR_WEIGHT(YYY,MW_Z2)
         YYY(I_Z_MIN) = 0._EB
         YYY(I_Z_MIN+1) = 0._EB      
         YYY(I_Z_MIN+2) = REACTION(1)%Z_F * (1._EB - REACTION(1)%SOOT_YIELD)
         YYY(I_Z_MAX)   = REACTION(1)%Z_F * REACTION(1)%SOOT_YIELD
         CALL GET_MOLECULAR_WEIGHT(YYY,MW_Z3)
      ELSE
         YYY(I_Z_MIN) = 0._EB
         YYY(I_Z_MIN+1) = REACTION(1)%Z_F
         YYY(I_Z_MAX) = 0._EB
         CALL GET_MOLECULAR_WEIGHT(YYY,MW_Z2)
         YYY(I_Z_MIN) = 0._EB
         YYY(I_Z_MIN+1) = 0._EB      
         YYY(I_Z_MAX) = REACTION(1)%Z_F
         CALL GET_MOLECULAR_WEIGHT(YYY,MW_Z3)
      ENDIF
   ELSE
      IF (SOOT_DEPOSITION) THEN
         MW_Z2 = SS0%MW
         YYY(I_Z_MIN) = 0._EB
         YYY(I_Z_MIN+1) = REACTION(1)%Z_F  * (1._EB - REACTION(1)%SOOT_YIELD)
         YYY(I_Z_MAX)   = REACTION(1)%Z_F * REACTION(1)%SOOT_YIELD
         CALL GET_MOLECULAR_WEIGHT(YYY,MW_Z3)      
      ELSE
         MW_Z2 = SS0%MW
         YYY(I_Z_MIN) = 0._EB
         YYY(I_Z_MIN+1) = REACTION(1)%Z_F
         CALL GET_MOLECULAR_WEIGHT(YYY,MW_Z3)      
      ENDIF
   ENDIF
   MW_MIN = MIN(MW_MIN,MW_Z2,MW_Z3,REACTION(1)%MW_FUEL)
   MW_MAX = MAX(MW_MAX,MW_Z2,MW_Z3,REACTION(1)%MW_FUEL)
ENDIF
 
! Compute background density from other background quantities
 
RHOA = P_INF/(TMPA*RSUM0)
 
! Compute constant-temperature specific heats
 
CP_GAMMA = SS0%RCON*GAMMA/(GAMMA-1._EB)
CPOPR = CP_GAMMA/PR
 
! Compute variable-temperature specific heats for specific species
! Compute viscosity, thermal conductivity for species 0 to N_SPECIES. Diffusivity for species 1 to N_SPECIES. 
! These terms are used for DNS, and as the lower limits in LES calcs
! Source: Poling, Prausnitz and O'Connell. Properties of Gases and Liquids, 5th ed, 2000.

IF(MIXTURE_FRACTION) THEN
   STATE_SPECIES(1) = 'ETHYLENE'
   STATE_SPECIES(2) = 'OXYGEN'
   STATE_SPECIES(3) = 'NITROGEN'
   STATE_SPECIES(4) = 'WATER VAPOR'
   STATE_SPECIES(5) = 'CARBON DIOXIDE'
   STATE_SPECIES(6) = 'CARBON MONOXIDE'
   STATE_SPECIES(7) = 'HYDROGEN'
   STATE_SPECIES(8) = 'CARBON MONOXIDE'
   STATE_SPECIES(9) = 'NITROGEN'
ENDIF

ALLOCATE(MU_TMP(N_Y_ARRAY))
MU_TMP = 0._EB
ALLOCATE(CP_TMP(N_Y_ARRAY))
CP_TMP = 0._EB
ALLOCATE(H_TMP(N_Y_ARRAY))
H_TMP = 0._EB
ALLOCATE(K_TMP(N_Y_ARRAY))
K_TMP = 0._EB
D_TMP = 0._EB

ALLOCATE(Y2CP_C(5000))
CALL ChkMemErr('READ','Y2CP_C',IZERO)    
Y2CP_C = 0._EB   

ALLOCATE(Y2CPBAR_C(5000))   
CALL ChkMemErr('READ','Y2CPBAR_C',IZERO)    
Y2CPBAR_C = 0._EB   

ALLOCATE(Y2H_G_C(5000))
CALL ChkMemErr('READ','Y2H_G_C',IZERO)    
Y2H_G_C = 0._EB   

ALLOCATE(Y2K_C(5000))
CALL ChkMemErr('READ','Y2K_C',IZERO)    
Y2K_C = 0._EB   

ALLOCATE(Y2MU_C(5000))      
CALL ChkMemErr('READ','Y2MU_C',IZERO)    
Y2MU_C = 0._EB   

ALLOCATE(Y2CP(5000,N_SPECIES))     
CALL ChkMemErr('READ','Y2CC',IZERO)    
Y2CP = 0._EB   

ALLOCATE(Y2CPBAR(5000,N_SPECIES))     
CALL ChkMemErr('READ','Y2CPBAR',IZERO)    
Y2CPBAR = 0._EB   

ALLOCATE(Y2H_G(5000,N_SPECIES))           
CALL ChkMemErr('READ','YH_GC',IZERO)    
Y2H_G = 0._EB   

ALLOCATE(Y2K(5000,N_SPECIES))      
CALL ChkMemErr('READ','Y2K',IZERO)    
Y2K = 0._EB   

ALLOCATE(Y2MU(5000,N_SPECIES))      
CALL ChkMemErr('READ','Y2MU',IZERO)    
Y2MU = 0._EB   

IF (MIXTURE_FRACTION) THEN
   ALLOCATE(Y2D_C(5000))
   CALL ChkMemErr('READ','Y2D_C',IZERO)    
   Y2D_C = 0._EB   
   ALLOCATE(Y2D(5000,I_Z_MAX-I_Z_MIN+1))
   CALL ChkMemErr('READ','Y2D_C',IZERO)    
   Y2D = 0._EB   
ENDIF

! Loop through temperatures from 1 K to 5000 K to get temperature-specific gas properties.  Data from JANAF 4

T_LOOP_1: DO J=1,5000
   IF (MIXTURE_FRACTION) THEN
      CALL JANAF_TABLE(J,CP_TMP(FUEL_INDEX),H_TMP(FUEL_INDEX),STATE_SPECIES(FUEL_INDEX),0._EB)
      CALL JANAF_TABLE(J,CP_TMP(O2_INDEX),H_TMP(O2_INDEX),STATE_SPECIES(O2_INDEX),0._EB)
      CALL JANAF_TABLE(J,CP_TMP(N2_INDEX),H_TMP(N2_INDEX),STATE_SPECIES(N2_INDEX),0._EB)
      CALL JANAF_TABLE(J,CP_TMP(CO2_INDEX),H_TMP(CO2_INDEX),STATE_SPECIES(CO2_INDEX),0._EB)      
      CALL JANAF_TABLE(J,CP_TMP(CO_INDEX),H_TMP(CO_INDEX),STATE_SPECIES(CO_INDEX),0._EB)                  
      CALL JANAF_TABLE(J,CP_TMP(H2O_INDEX),H_TMP(H2O_INDEX),STATE_SPECIES(H2O_INDEX),0._EB)                  
      CALL JANAF_TABLE(J,CP_TMP(H2_INDEX),H_TMP(H2_INDEX),STATE_SPECIES(H2_INDEX),0._EB)                  
      CALL JANAF_TABLE(J,CP_TMP(SOOT_INDEX),H_TMP(SOOT_INDEX),STATE_SPECIES(SOOT_INDEX),0._EB)                  
      CALL JANAF_TABLE(J,CP_TMP(OTHER_INDEX),H_TMP(OTHER_INDEX),STATE_SPECIES(N2_INDEX),0._EB)   
      SUB_SPECIES_LOOP2: DO NN=1,N_STATE_SPECIES
         SIG_N = -1._EB
         EPSK_N = -1._EB 
         MW_N = -.1_EB
         CALL GAS_PROPS(STATE_SPECIES(NN),SIG_N,EPSK_N,MW_N,ABSORBING,FORMULA)
         !Conductivity
         SIGMA2 = SIG_N**2
         TSTAR = J/EPSK_N
         OMEGA = 1.16145_EB*TSTAR**(-0.14874_EB) + 0.52487_EB*EXP(-0.77320_EB*TSTAR) + 2.16178_EB*EXP(-2.43787_EB*TSTAR)
         MU_TMP(NN) = 26.69E-7_EB*(MW_N*J)**0.5_EB/(SIGMA2*OMEGA)/MW_N
         K_TMP(NN)  = MU_TMP(NN) * CP_TMP(NN) / PR
         !Diffusivity
         SIGMA2 = (0.5_EB*(SIG_N+SIG(0)))**2
         EPSIJ  = SQRT(EPSK_N*EPSK(0))
         AMW    = SQRT( (MW_N+SS0%MW)/(2._EB*MW_N*SS0%MW) )
         TSTAR = J/EPSIJ
         OMEGA = 1.06036_EB/TSTAR**(0.15610_EB) + 0.19300_EB/EXP(0.47635_EB*TSTAR) + 1.03587_EB/EXP(1.52996_EB*TSTAR) + &
                  1.76474_EB/EXP(3.89411_EB*TSTAR)
         D_TMP(NN) = 0.00266E-4_EB*AMW*(J)**1.5_EB/(SIGMA2*OMEGA)
      ENDDO SUB_SPECIES_LOOP2
      Y2CP_C(J)    = SUM(Y2Y_C(:) * CP_TMP(:))
      IF (J>1) THEN
         Y2H_G_C(J)   = Y2H_G_C(J-1) + 0.5_EB*(Y2CP_C(J)+Y2CP_C(J-1)) 
         Y2CPBAR_C(J) = (Y2CPBAR_C(J-1)*(REAL(J,EB)-1._EB)+Y2CP_C(J))/REAL(J,EB)
      ELSE
         Y2H_G_C(J)   =  SUM(Y2Y_C(:)* H_TMP(:)) + Y2CP_C(J) 
         Y2CPBAR_C(J) =   Y2CP_C(J)
      ENDIF
      Y2MU_C(J)    = SUM(Y2Y_C(:) * MU_TMP(:))
      Y2K_C(J)     = SUM(Y2Y_C(:) * K_TMP(:))      
      Y2D_C(J)     = SUM(Y2Y_C(:) * D_TMP(:))            
   ENDIF 
   ! Compute terms for non-mixture fraction species
   SPECIES_LOOP_1: DO N = 0,N_SPECIES
      SS => SPECIES(N)
      IF(SS%MODE==MIXTURE_FRACTION_SPECIES) CYCLE SPECIES_LOOP_1
      SIGMA2 = (0.5_EB*(SIG(N)+SIG(0)))**2
      EPSIJ  = SQRT(EPSK(N)*EPSK(0))
      AMW    = SQRT( (SS%MW+SS0%MW)/(2._EB*SS%MW*SS0%MW) )
      TSTAR = J/EPSIJ
      OMEGA = 1.06036_EB/TSTAR**(0.15610_EB) + 0.19300_EB/EXP(0.47635_EB*TSTAR) + 1.03587_EB/EXP(1.52996_EB*TSTAR) &
            + 1.76474_EB/EXP(3.89411_EB*TSTAR)     
      D_N0 = 0.00266E-4_EB*AMW*(J)**1.5_EB/(SIGMA2*OMEGA)
      IF (D_USER(N)>=0._EB .AND. .NOT.CONSTANT_PROPERTIES) D_N0 = D_USER(N)*(J/TMPA)**1.75_EB
      IF (D_USER(N)>=0._EB .AND. CONSTANT_PROPERTIES)      D_N0 = D_USER(N)
      SS%D(J) = D_N0
      
      IF (MIXTURE_FRACTION .AND. SS%INDEX <= N_STATE_SPECIES .OR. SS%INDEX > N_Y_ARRAY) CYCLE SPECIES_LOOP_1
      CALL JANAF_TABLE(J,CP_TMP(SS%INDEX),H_TMP(SS%INDEX),SS%ID,SS%RCON)
      IF (SS%SPECIFIC_HEAT > 0._EB) CP_TMP(SS%INDEX) = SS%SPECIFIC_HEAT
      IF (SS%SPECIFIC_ENTHALPY > 0._EB) H_TMP(SS%INDEX) = SS%SPECIFIC_ENTHALPY - SS%SPECIFIC_HEAT * SS%REFERENCE_TEMPERATURE
      SIGMA2 = SIG(N)**2
      TSTAR = J/EPSK(N)
      OMEGA = 1.16145_EB*TSTAR**(-0.14874_EB) +  0.52487_EB*EXP(-0.77320_EB*TSTAR) + 2.16178_EB*EXP(-2.43787_EB*TSTAR)
      MU_TMP(SS%INDEX) = 26.69E-7_EB*(SS%MW*J)**0.5_EB/(SIGMA2*OMEGA)/SS%MW
      IF (MU_USER(N)>=0._EB .AND. .NOT.CONSTANT_PROPERTIES) MU_TMP(SS%INDEX) = MU_USER(N)*(J/TMPA)**0.75_EB/SS%MW
      IF (MU_USER(N)>=0._EB .AND. CONSTANT_PROPERTIES)      MU_TMP(SS%INDEX) = MU_USER(N)/SS%MW
      K_TMP(SS%INDEX) = MU_TMP(SS%INDEX)*CP_TMP(SS%INDEX)/PR
      IF (K_USER(N)>=0._EB .AND. .NOT.CONSTANT_PROPERTIES)  K_TMP(SS%INDEX) = K_USER(N)*(J/TMPA)**0.75_EB/SS%MW
      IF (K_USER(N)>=0._EB .AND. CONSTANT_PROPERTIES)       K_TMP(SS%INDEX) = K_USER(N)/SS%MW
   ENDDO SPECIES_LOOP_1
   ! Constant coefficients in non-mixture fraction case
   
   IF (.NOT. MIXTURE_FRACTION) THEN
      Y2CP_C(J) = CP_TMP(SPECIES(0)%INDEX)
      IF (J>1) THEN
         Y2H_G_C(J)   = Y2H_G_C(J-1) + 0.5_EB*(Y2CP_C(J)+Y2CP_C(J-1)) 
         Y2CPBAR_C(J) = (Y2CPBAR_C(J-1)*REAL(J-1,EB)+Y2CP_C(J))/REAL(J,EB)
      ELSE
         Y2H_G_C(J)   =   H_TMP(SPECIES(0)%INDEX) + Y2CP_C(J) 
         Y2CPBAR_C(J) =   Y2CP_C(J)
      ENDIF
      Y2MU_C(J)    = MU_TMP(SPECIES(0)%INDEX)
      Y2K_C(J)     = K_TMP(SPECIES(0)%INDEX)
   ENDIF
   ! Assemble coefficients for both primitive and mixture fraction species
   SPECIES_LOOP_2: DO N=1,N_SPECIES
      Y2CP(J,N)    = SUM(Y2Y(:,N) * CP_TMP(:))
      IF (J>1) THEN
         Y2H_G(J,N)   = Y2H_G(J-1,N) +0.5_EB*(Y2CP(J,N)+Y2CP(J-1,N)) 
         Y2CPBAR(J,N) = (Y2CPBAR(J-1,N)*REAL(J-1,EB)+Y2CP(J,N))/REAL(J,EB)   
      ELSE
         Y2H_G(J,N)   = SUM(Y2Y(:,N) * H_TMP(:)) + Y2CP(J,N) 
         Y2CPBAR(J,N) = Y2CP(J,N)
      ENDIF
      Y2MU(J,N)    = SUM(Y2Y(:,N) * MU_TMP(:))
      Y2K(J,N)     = SUM(Y2Y(:,N) * K_TMP(:))
      IF (MIXTURE_FRACTION) THEN
         IF (SPECIES(N)%MODE==MIXTURE_FRACTION_SPECIES) THEN
            Y2D(J,N-I_Z_MIN+1) = SUM(Y2Y(1:N_STATE_SPECIES,N) * D_TMP(1:N_STATE_SPECIES))
         ENDIF
      ENDIF
   ENDDO SPECIES_LOOP_2
ENDDO T_LOOP_1

! Set ISFUEL property for MIXTURE_FRACTION case

IF (MIXTURE_FRACTION) SPECIES(I_FUEL)%ISFUEL = .TRUE.

! For finite rate reaction, set parameters
 
FINITE_RATE_REACTION_LOOP: DO NN=1,N_REACTIONS
   RN => REACTION(NN)
   IF (RN%MODE/=FINITE_RATE_REACTION)  CYCLE FINITE_RATE_REACTION_LOOP
   DO N=1,N_SPECIES
      IF (RN%FUEL    ==SPECIES(N)%ID) THEN
         RN%I_FUEL         = N
         SPECIES(N)%ISFUEL = .TRUE.
         IF (RN%N(N) /=-999._EB) THEN
            RN%N_F         = RN%N(N)
         ELSE
            RN%N_F         = 1._EB         
         ENDIF
      ENDIF
      IF (RN%OXIDIZER==SPECIES(N)%ID) THEN
         RN%I_OXIDIZER = N
         IF (RN%N(N) /=-999._EB) THEN
            RN%N_O        = RN%N(N)
         ELSE
            RN%N_O        = 1._EB         
         ENDIF
      ENDIF
   ENDDO
   IF (NN==1) I_FUEL = RN%I_FUEL
   RN%MW_FUEL = SPECIES(RN%I_FUEL)%MW
   IF (RN%NU(RN%I_FUEL)     == 0._EB) RN%NU(RN%I_FUEL)     = -1._EB
   IF (RN%NU(RN%I_FUEL)     >  0._EB) RN%NU(RN%I_FUEL)     = -RN%NU(RN%I_FUEL)
   IF (RN%NU(RN%I_OXIDIZER) >  0._EB) RN%NU(RN%I_OXIDIZER) = -RN%NU(RN%I_OXIDIZER)
   IF (RN%NU(RN%I_OXIDIZER) == 0._EB) THEN
      WRITE(MESSAGE,'(A)')  'ERROR: Specify a stoichiometric coefficient for oxidizer'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   RN%O2_F_RATIO = SPECIES(RN%I_OXIDIZER)%MW*RN%NU(RN%I_OXIDIZER)/(SPECIES(RN%I_FUEL)%MW*RN%NU(RN%I_FUEL))
   RN%EPUMO2     = RN%HEAT_OF_COMBUSTION/RN%O2_F_RATIO
ENDDO FINITE_RATE_REACTION_LOOP

CONTAINS
  
 
SUBROUTINE GAS_PROPS(GAS_NAME,SIGMA,EPSOK,MW,ABSORBING,FORMULA)
 
! Molecular weight (g/mol) and Lennard-Jones properties
 
REAL(EB) SIGMA,EPSOK,MW,SIGMAIN,EPSOKIN,MWIN
CHARACTER(30) GAS_NAME,FORMULA
LOGICAL ABSORBING
 
SIGMAIN = SIGMA
EPSOKIN = EPSOK
MWIN    = MW
 
SELECT CASE(GAS_NAME)
   CASE('AIR')             
      SIGMA=3.711_EB 
      EPSOK= 78.6_EB  
      MW=MW_AIR
      FORMULA='Air'
   CASE('CARBON MONOXIDE') 
      SIGMA=3.690_EB 
      EPSOK= 91.7_EB  
      MW=28._EB
      ABSORBING = .TRUE.
      FORMULA='CO'
   CASE('CARBON DIOXIDE')  
      SIGMA=3.941_EB 
      EPSOK=195.2_EB  
      MW=44._EB
      ABSORBING = .TRUE.      
      FORMULA='CO2'
   CASE('ETHYLENE')        
      SIGMA=4.163_EB 
      EPSOK=224.7_EB 
      MW=28._EB
      ABSORBING = .TRUE.
      FORMULA='C2H2'
   CASE('ARGON')
      SIGMA=3.42_EB
      EPSOK= 124.0_EB
      MW= 40._EB
      FORMULA='Ar'
   CASE('HELIUM')          
      SIGMA=2.551_EB 
      EPSOK= 10.22_EB 
      MW= 4._EB
      FORMULA='He2'
   CASE('HYDROGEN')        
      SIGMA=2.827_EB 
      EPSOK= 59.7_EB
      MW= 2._EB
      FORMULA='H2'
   CASE('METHANE')         
      SIGMA=3.758_EB 
      EPSOK=148.6_EB  
      MW=16._EB
      ABSORBING = .TRUE.
      FORMULA='CH4'
   CASE('NITROGEN')        
      SIGMA=3.798_EB 
      EPSOK= 71.4_EB  
      MW=28._EB
      FORMULA='N2'
   CASE('OXYGEN')          
      SIGMA=3.467_EB 
      EPSOK=106.7_EB  
      MW=32._EB
      FORMULA='O2'
   CASE('PROPANE')         
      SIGMA=5.118_EB 
      EPSOK=237.1_EB  
      MW=44._EB
      ABSORBING = .TRUE.
      FORMULA='C3H8'
   CASE('WATER VAPOR')     
      SIGMA=2.641_EB 
      EPSOK=809.1_EB  
      MW=18._EB
      ABSORBING = .TRUE.
      FORMULA='H2O'
   CASE DEFAULT            
      SIGMA=3.711_EB 
      EPSOK= 78.6_EB
      MW=MW_AIR
   CASE('SOOT')        
      SIGMA=3.798_EB 
      EPSOK= 71.4_EB  
      MW=12._EB
      FORMULA='C'
END SELECT  
 
IF (SIGMAIN>0._EB) SIGMA = SIGMAIN
IF (EPSOKIN>0._EB) EPSOK = EPSOKIN
IF (MWIN   >0._EB) MW    = MWIN 
 
IF (GAS_NAME=='MIXTURE_FRACTION') MW = 1._EB/(REACTION(1)%Y_F_INLET/REACTION(1)%MW_FUEL+REACTION(1)%Y_N2_INLET/MW_N2)
 
END SUBROUTINE GAS_PROPS
  
END SUBROUTINE PROC_SPEC
 
SUBROUTINE PROC_PART
USE PHYSICAL_FUNCTIONS, ONLY : JANAF_TABLE_LIQUID, GET_SPECIFIC_ENTHALPY
INTEGER :: N, J, REF_TEMP
REAL(EB) :: H_L,H_V,YY_GET(1:N_SPECIES),H_G_S,H_G_S_REF
TYPE(PARTICLE_CLASS_TYPE), POINTER :: PC
TYPE(SPECIES_TYPE), POINTER :: SS

IF (N_PART == 0) RETURN

DO N = 1, N_PART
   PC => PARTICLE_CLASS(N)
   SS => SPECIES(PC%SPEC_INDEX)
   IF (.NOT.PC%EVAPORATE) CYCLE
   DO J = 1, 5000
      IF (PC%C_P(J) > 0._EB) THEN
         PC%H_L = (REAL(J,EB)-298.15)*PC%C_P(J)
      ELSE
         CALL JANAF_TABLE_LIQUID (J,PC%C_P(J),H_V,H_L,SS%ID)
         IF (J==1) THEN
            PC%H_L(J) = H_L + PC%C_P(J)
         ELSE
            PC%H_L(J) = PC%H_L(J-1) + PC%C_P(J)
         ENDIF
      ENDIF
      IF (J==1) THEN
         PC%C_P_BAR(J) = PC%C_P(J)
      ELSE
         PC%C_P_BAR(J) = (PC%C_P_BAR(J-1) * REAL(J-1,EB) + PC%C_P(J)) / REAL(J,EB)
      ENDIF
   ENDDO
   IF (PC%SPEC_INDEX>0 .AND. PC%H_V(1)<0._EB) THEN
      REF_TEMP = NINT(PC%H_V_REFERENCE_TEMPERATURE)
      YY_GET = 0._EB
      YY_GET(PC%SPEC_INDEX) = 1._EB
      CALL GET_SPECIFIC_ENTHALPY(YY_GET,H_G_S_REF,REF_TEMP)
      DO J=1,5000
         CALL GET_SPECIFIC_ENTHALPY(YY_GET,H_G_S,J)
         PC%H_V(J) = H_V + (H_G_S-H_G_S_REF) - (PC%H_L(J)-PC%H_L(REF_TEMP))
      ENDDO
   ENDIF
ENDDO
RETURN

END SUBROUTINE PROC_PART
 
SUBROUTINE READ_PART
USE DEVICE_VARIABLES, ONLY : PROPERTY_TYPE, PROPERTY, N_PROP
USE RADCONS, ONLY : NDG
INTEGER :: NUMBER_INITIAL_DROPLETS,SAMPLING_FACTOR,DROPLETS_PER_SECOND,N,NN,IPC,RGB(3),I_DUMMY
REAL(EB) :: SPECIFIC_HEAT,VAPORIZATION_TEMPERATURE,DUMMY, MELTING_TEMPERATURE,MASS_PER_VOLUME,DIAMETER, &
            GAMMA_D,AGE,INITIAL_TEMPERATURE,HEAT_OF_COMBUSTION,HEAT_OF_VAPORIZATION,DENSITY,DT_INSERT, &
            VERTICAL_VELOCITY,HORIZONTAL_VELOCITY,MAXIMUM_DIAMETER,MINIMUM_DIAMETER,SIGMA_D,LENGTH,WIDTH,ELEVATION,AZIMUTH
REAL(EB)    VEG_SV,VEG_MOISTURE,VEG_CHAR_FRACTION,VEG_DRAG_COEFFICIENT,VEG_DENSITY,VEG_BULK_DENSITY, &
            VEG_BURNING_RATE_MAX,VEG_DEHYDRATION_RATE_MAX,VEG_INITIAL_TEMPERATURE, &
            VEG_FUEL_MPV_MIN,VEG_MOIST_MPV_MIN,H_V_REFERENCE_TEMPERATURE
CHARACTER(30) :: PART_ID,SPEC_ID,QUANTITIES(1:10),SURF_ID
CHARACTER(60) :: SMOKEVIEW_ID
CHARACTER(25) :: COLOR
LOGICAL :: MASSLESS,STATIC,FUEL,WATER,TREE,MONODISPERSE,EVAPORATE
LOGICAL :: VEG_REMOVE_CHARRED,VEG_STEM
TYPE(PARTICLE_CLASS_TYPE), POINTER :: PC
TYPE(PROPERTY_TYPE), POINTER :: PY
NAMELIST /PART/ NUMBER_INITIAL_DROPLETS,FYI,DROPLETS_PER_SECOND,MASS_PER_VOLUME, &
                SAMPLING_FACTOR,ID,STATIC,MASSLESS,FUEL,WATER,TREE, &
                DENSITY,VAPORIZATION_TEMPERATURE,SPECIFIC_HEAT,HEAT_OF_VAPORIZATION, &
                MELTING_TEMPERATURE,DIAMETER,MAXIMUM_DIAMETER,MINIMUM_DIAMETER,GAMMA_D,HEAT_OF_COMBUSTION, &
                AGE,SPEC_ID,SURF_ID,INITIAL_TEMPERATURE,XB,RGB,QUANTITIES,DT_INSERT,COLOR, &
                VERTICAL_VELOCITY,HORIZONTAL_VELOCITY,GAMMA,MONODISPERSE,EVAPORATE,SIGMA_D, &
                VEG_SV,VEG_MOISTURE,VEG_CHAR_FRACTION,VEG_DRAG_COEFFICIENT,VEG_DENSITY,VEG_BULK_DENSITY, &
                VEG_BURNING_RATE_MAX,VEG_DEHYDRATION_RATE_MAX,VEG_INITIAL_TEMPERATURE,VEG_FUEL_MPV_MIN, &
                VEG_MOIST_MPV_MIN,VEG_REMOVE_CHARRED,VEG_STEM,SMOKEVIEW_ID,LENGTH,WIDTH,ELEVATION,AZIMUTH,H_V_REFERENCE_TEMPERATURE

! Determine total number of PART lines in the input file

REWIND(LU_INPUT)
N_PART = 0
COUNT_PART_LOOP: DO
   CALL CHECKREAD('PART',LU_INPUT,IOS) 
   IF (IOS==1) EXIT COUNT_PART_LOOP
   READ(LU_INPUT,PART,END=219,ERR=220,IOSTAT=IOS)
   N_PART = N_PART + 1
   220 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with PART line')
ENDDO COUNT_PART_LOOP
219 REWIND(LU_INPUT)
 
! Allocate the derived type array to hold information about the particle classes

ALLOCATE(PARTICLE_CLASS(N_PART),STAT=IZERO)
CALL ChkMemErr('READ','PARTICLE_CLASS',IZERO) 

READ_PART_LOOP: DO N=1,N_PART
 
   PC=>PARTICLE_CLASS(N)
   DENSITY                  = 1000._EB     ! kg/m3
   DT_INSERT                = 0.01_EB      ! s
   MASS_PER_VOLUME          = 1._EB        ! kg/m3
   VAPORIZATION_TEMPERATURE = 100.0_EB     ! C
   INITIAL_TEMPERATURE      = TMPA - TMPM  ! C
   MELTING_TEMPERATURE      = TMPM - TMPM  ! C
   SPECIFIC_HEAT            = -1._EB     ! kJ/kg-K
   HEAT_OF_VAPORIZATION     = -1._EB     ! kJ/kg
   H_V_REFERENCE_TEMPERATURE = 0._EB
   HEAT_OF_COMBUSTION       = -1._EB       ! kJ/kg
   DROPLETS_PER_SECOND      = 5000
   GAMMA                    = 1.4          !Specific heat ratio
   DIAMETER                 = 500._EB      ! microns
   MAXIMUM_DIAMETER         = 1.E9_EB      ! microns, meant to be infinitely large and not used
   MINIMUM_DIAMETER         =  20._EB      ! microns, below which the droplet evaporates in one time step
   MONODISPERSE             = .FALSE.
   LENGTH                   = 0.1_EB
   WIDTH                    = 0.1_EB
   AZIMUTH                  = 0._EB
   ELEVATION                = 0._EB
   GAMMA_D                  = 2.4_EB
   SIGMA_D                  = -99999.9_EB
   AGE                      = 1.E6_EB      ! s
   ID                       = 'null'
   QUANTITIES               = 'null'
   RGB                      = -1
   SPEC_ID                  = 'null'
   SURF_ID                  = 'null'
   COLOR                    = 'null'
   SAMPLING_FACTOR          = -1      
   NUMBER_INITIAL_DROPLETS  = 0
   XB                       = 0._EB
   EVAPORATE                = .TRUE.
   FUEL                     = .FALSE.
   WATER                    = .FALSE.
   STATIC                   = .FALSE.
   MASSLESS                 = .FALSE.
   TREE                     = .FALSE.
   VERTICAL_VELOCITY        = 0.5_EB
   HORIZONTAL_VELOCITY      = 0.2_EB
   SMOKEVIEW_ID             = ' '
   VEG_SV                   = 4000. !1/m
   VEG_MOISTURE             = 10.0_EB
   VEG_CHAR_FRACTION        = 0.25_EB
   VEG_DRAG_COEFFICIENT     = 1.0_EB
   VEG_DENSITY              = 540._EB !kg/m^3 dry mass
   VEG_BULK_DENSITY         = 0.3_EB !kg/m^3
   VEG_BURNING_RATE_MAX     = 0.4_EB !kg/(m^3.s)
   VEG_DEHYDRATION_RATE_MAX = 0.4_EB !kg/(m^3.s)
   VEG_INITIAL_TEMPERATURE  = TMPA - TMPM  ! C
   VEG_FUEL_MPV_MIN         = VEG_CHAR_FRACTION*VEG_BULK_DENSITY
   VEG_MOIST_MPV_MIN        = 0.01_EB*VEG_MOISTURE
   VEG_REMOVE_CHARRED       = .TRUE.
   VEG_STEM                 = .FALSE.
 
   ! Read the PART line from the input file or set up special PARTICLE_CLASS class for water droplets or tracers
 
   CALL CHECKREAD('PART',LU_INPUT,IOS) 
   IF (IOS==1) EXIT READ_PART_LOOP
   READ(LU_INPUT,PART)

   IF (ANY(RGB<0) .AND. COLOR=='null') THEN
      COLOR = 'BLACK'
      IF (WATER) COLOR = 'BLUE'
      IF (TREE)  COLOR = 'GREEN'
      IF (FUEL)  COLOR = 'YELLOW'
   ENDIF
   IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)

   ! Miscellaneous consequences of input parameters

   IF (TREE)                                    STATIC                  = .TRUE.
   IF (SURF_ID/='null')                         STATIC                  = .TRUE.
   IF (MASSLESS)                                DIAMETER                = 0._EB
   IF (FUEL)                                    SPEC_ID                 = 'MIXTURE_FRACTION_1'
   IF (WATER)                                   SPEC_ID                 = 'WATER VAPOR'
   IF (FUEL)                                    FUEL_EVAPORATION        = .TRUE.
   IF (WATER)                                   WATER_EVAPORATION       = .TRUE.
   IF (WATER)                                   GAMMA                   = 1.32_EB
   IF (WATER)                                   WATER_PART_INDEX        = N
   IF (SAMPLING_FACTOR<=0 .AND.      MASSLESS)  SAMPLING_FACTOR         = 1
   IF (SAMPLING_FACTOR<=0 .AND. .NOT.MASSLESS)  SAMPLING_FACTOR         = 10
   IF (NUMBER_INITIAL_DROPLETS>0 .AND. RESTART) NUMBER_INITIAL_DROPLETS =  0

   IF (SPEC_ID/='null' .AND. MASSLESS) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: Cannot have MASSLESS=.TRUE. with evaporating droplets'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   PC%QUANTITIES = QUANTITIES
   ! Set up arrays in case the domain is to be seeded with droplets/particles

   IF (NUMBER_INITIAL_DROPLETS>0) DROPLET_FILE  = .TRUE.
   IF (NUMBER_INITIAL_DROPLETS>0) THEN
      DO I=1,5,2
         IF (XB(I)>XB(I+1)) THEN
            DUMMY   = XB(I)
            XB(I)   = XB(I+1)
            XB(I+1) = DUMMY
         ENDIF
      ENDDO
   ENDIF
   PC%X1 = XB(1)
   PC%X2 = XB(2)
   PC%Y1 = XB(3)
   PC%Y2 = XB(4)
   PC%Z1 = XB(5)
   PC%Z2 = XB(6)
   PC%GAMMA_VAP = GAMMA
   PC%N_INITIAL = NUMBER_INITIAL_DROPLETS
 
   ! Arrays for particle size distribution

   IF (DIAMETER > 0._EB) THEN
      ALLOCATE(PC%CDF(0:NDC),STAT=IZERO)
      CALL ChkMemErr('READ','CDF',IZERO)
      ALLOCATE(PC%R_CDF(0:NDC),STAT=IZERO)
      CALL ChkMemErr('READ','R_CDF',IZERO)
      ALLOCATE(PC%IL_CDF(NSTRATA),STAT=IZERO)
      CALL ChkMemErr('READ','IL_CDF',IZERO)
      ALLOCATE(PC%IU_CDF(NSTRATA),STAT=IZERO)
      CALL ChkMemErr('READ','IU_CDF',IZERO)
      ALLOCATE(PC%W_CDF(NSTRATA),STAT=IZERO)
      CALL ChkMemErr('READ','W_CDF',IZERO)
   ENDIF
   ALLOCATE(PC%INSERT_CLOCK(NMESHES),STAT=IZERO)
   CALL ChkMemErr('READ','INSERT_CLOCK',IZERO)
   PC%INSERT_CLOCK = T_BEGIN
   
   ! Assign property data to PARTICLE_CLASS class
 
   PC%ID                 = ID
   PC%SMOKEVIEW_ID       = SMOKEVIEW_ID
   PC%DT_INSERT          = DT_INSERT
   PC%N_INSERT           = DROPLETS_PER_SECOND*DT_INSERT
   PC%SAMPLING           = SAMPLING_FACTOR
   PC%RGB                = RGB
   PC%DIAMETER           = DIAMETER*1.E-6_EB
   PC%MAXIMUM_DIAMETER   = MAXIMUM_DIAMETER*1.E-6_EB
   PC%MINIMUM_DIAMETER   = MINIMUM_DIAMETER*1.E-6_EB
   PC%KILL_RADIUS        = MINIMUM_DIAMETER*1.E-6_EB*0.25_EB
   PC%MONODISPERSE       = MONODISPERSE
   PC%LENGTH             = LENGTH
   PC%WIDTH              = WIDTH
   PC%AZIMUTH            = AZIMUTH
   PC%ELEVATION          = ELEVATION
   PC%GAMMA              = GAMMA_D
   IF ( SIGMA_D > 0._EB ) THEN
      PC%SIGMA           = SIGMA_D
   ELSE
      PC%SIGMA           = 1.15_EB/GAMMA_D
   END IF
   PC%TMP_INITIAL        = INITIAL_TEMPERATURE + TMPM
   PC%H_V_REFERENCE_TEMPERATURE = H_V_REFERENCE_TEMPERATURE + 273.15_EB
   IF ((HEAT_OF_VAPORIZATION >  0._EB .AND. SPECIFIC_HEAT <= 0._EB) .OR. &
       (HEAT_OF_VAPORIZATION <= 0._EB .AND. SPECIFIC_HEAT >  0._EB)) THEN
      WRITE(MESSAGE,'(A,I4,A)') 'ERROR: PART ' ,N,' If one of SPECIFIC_HEAT or HEAT_OF_VAPORIZATION defined, both must be'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (HEAT_OF_VAPORIZATION > 0._EB) PC%H_V = HEAT_OF_VAPORIZATION*1000._EB
   IF (SPECIFIC_HEAT > 0._EB) PC%C_P = SPECIFIC_HEAT*1000._EB
   PC%HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION*1000._EB
   PC%TMP_V              = VAPORIZATION_TEMPERATURE + TMPM
   PC%DENSITY            = DENSITY
   PC%FTPR               = FOTH*PI*DENSITY
   PC%MASS_PER_VOLUME    = MASS_PER_VOLUME
   PC%TMP_MELT           = MELTING_TEMPERATURE + TMPM
   PC%MASSLESS           = MASSLESS
   PC%LIFETIME           = AGE
   PC%EVAPORATE          = EVAPORATE
   PC%TREE               = TREE
   PC%FUEL               = FUEL
   PC%WATER              = WATER
   PC%STATIC             = STATIC
   PC%SPEC_ID            = SPEC_ID
   PC%SPEC_INDEX         = 0       ! SPECies have not yet been read in
   PC%SURF_ID            = SURF_ID
   PC%SURF_INDEX         = -1      ! SURFs have not yet been read in
   IF (SURF_ID/='null') VIRTUAL_PARTICLES = .TRUE.
   PC%ADJUST_EVAPORATION = 1._EB   ! If H_O_C>0. this parameter will have to be reset later
   PC%VERTICAL_VELOCITY  = VERTICAL_VELOCITY
   PC%HORIZONTAL_VELOCITY= HORIZONTAL_VELOCITY
 
   ! Vegetation propeties

   PC%VEG_SV                   = VEG_SV !1/m
   PC%VEG_MOISTURE             = VEG_MOISTURE
   PC%VEG_CHAR_FRACTION        = VEG_CHAR_FRACTION
   PC%VEG_DRAG_COEFFICIENT     = VEG_DRAG_COEFFICIENT
   PC%VEG_DENSITY              = VEG_DENSITY !kg/m^3
   PC%VEG_BULK_DENSITY         = VEG_BULK_DENSITY !kg/m^3
   PC%VEG_BURNING_RATE_MAX     = VEG_BURNING_RATE_MAX !kg/m^3.s
   PC%VEG_DEHYDRATION_RATE_MAX = VEG_DEHYDRATION_RATE_MAX !kg/m^3.s
   PC%VEG_INITIAL_TEMPERATURE  = VEG_INITIAL_TEMPERATURE +TMPM ! K
   PC%VEG_FUEL_MPV_MIN         = VEG_CHAR_FRACTION*VEG_BULK_DENSITY
   PC%VEG_MOIST_MPV_MIN        = 0.01_EB*VEG_MOISTURE
   PC%VEG_REMOVE_CHARRED       = VEG_REMOVE_CHARRED
   PC%VEG_STEM                 = VEG_STEM
  
ENDDO READ_PART_LOOP

! Assign PART_INDEX to Device PROPERTY array

DO N=1,N_PROP
   PY => PROPERTY(N)
   PY%PART_INDEX = 0
   IF (PY%PART_ID/='null') THEN
      DO IPC=1,N_PART
         PC=>PARTICLE_CLASS(IPC)
         IF (PC%ID==PY%PART_ID) PY%PART_INDEX = IPC
         IF (PC%ID==PY%PART_ID .AND. PC%MASSLESS) THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: PART_ID for PROP ' ,N,' cannot refer to MASSLESS particles'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO
      IF (PY%PART_INDEX==0) THEN
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: PART_ID for PROP ' ,N,' not found'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      DROPLET_FILE=.TRUE.
   ENDIF
ENDDO

! Allocate radiation arrays

N_EVAP_INDICIES = 0
PLOOP2: DO IPC=1,N_PART
   PC=>PARTICLE_CLASS(IPC)
   IF (.NOT.PC%MASSLESS) THEN
      ALLOCATE(PC%WQABS(0:NDG,1:NUMBER_SPECTRAL_BANDS))
      CALL ChkMemErr('INIT','WQABS',IZERO)
      PC%WQABS = 0._EB
      ALLOCATE(PC%WQSCA(0:NDG,1:NUMBER_SPECTRAL_BANDS))
      CALL ChkMemErr('INIT','WQSCA',IZERO)
      PC%WQSCA = 0._EB
      ALLOCATE(PC%KWR(0:NDG))
      CALL ChkMemErr('INIT','KWR',IZERO)
      PC%KWR = 0._EB
      N_EVAP_INDICIES = N_EVAP_INDICIES + 1
      PC%EVAP_INDEX   = N_EVAP_INDICIES
   ENDIF
ENDDO PLOOP2

! Determine output quantities

DO N=1,N_PART
   PC=>PARTICLE_CLASS(N)
   PC%N_QUANTITIES = 0
   IF (ANY(PC%QUANTITIES/='null')) THEN
      PART_ID='null'
      SPEC_ID='null'
      QUANTITIES_LOOP: DO NN=1,10
         IF (PC%QUANTITIES(NN)=='null') CYCLE QUANTITIES_LOOP
         PC%N_QUANTITIES = PC%N_QUANTITIES + 1
         CALL GET_QUANTITY_INDEX(PC%SMOKEVIEW_LABEL(PC%N_QUANTITIES),PC%SMOKEVIEW_BAR_LABEL(PC%N_QUANTITIES), &
                                 PC%QUANTITIES_INDEX(PC%N_QUANTITIES),I_DUMMY,I_DUMMY,'PART',PC%QUANTITIES(NN),SPEC_ID,PART_ID)
      ENDDO QUANTITIES_LOOP
   ENDIF
ENDDO 

END SUBROUTINE READ_PART
 
 
SUBROUTINE READ_TREE

USE GLOBAL_CONSTANTS
INTEGER :: IPC,N_TREES_0,NM,NN,N
REAL(EB) :: X_TREE_MIN,X_TREE_MAX,Y_TREE_MIN,Y_TREE_MAX,Z_TREE_MIN,Z_TREE_MAX, &
            X_OVERALL_MIN,X_OVERALL_MAX,Y_OVERALL_MIN,Y_OVERALL_MAX, &
            Z_OVERALL_MIN,Z_OVERALL_MAX
REAL(EB) :: CROWN_WIDTH,CROWN_WIDTH_BOTTOM,CROWN_WIDTH_TOP,CROWN_BASE_HEIGHT,TREE_HEIGHT,XYZ(3)
REAL(EB) :: RING_THICKNESS,TON_IGNITOR_ELEMENTS,TOFF_IGNITOR_ELEMENTS, &
            T_RAMPOFF_IGNITOR_ELEMENTS,T_RAMPON_IGNITOR_ELEMENTS
LOGICAL  :: IGNITOR_ELEMENTS
CHARACTER(30) :: PART_ID, LABEL
CHARACTER(30) :: FUEL_GEOM,RAMP_IGNELEM
TYPE(PARTICLE_CLASS_TYPE), POINTER :: PC
NAMELIST /TREE/ XYZ,XB,CROWN_WIDTH,CROWN_WIDTH_BOTTOM,CROWN_WIDTH_TOP,CROWN_BASE_HEIGHT,TREE_HEIGHT,PART_ID, &
                IGNITOR_ELEMENTS,RING_THICKNESS,TON_IGNITOR_ELEMENTS,  &
                TOFF_IGNITOR_ELEMENTS,T_RAMPOFF_IGNITOR_ELEMENTS,      &
                T_RAMPON_IGNITOR_ELEMENTS,FUEL_GEOM,RAMP_IGNELEM,LABEL

X_OVERALL_MIN = 1.E12_EB ; X_OVERALL_MAX = -1.E12_EB
Y_OVERALL_MIN = 1.E12_EB ; Y_OVERALL_MAX = -1.E12_EB
Z_OVERALL_MIN = 1.E12_EB ; Z_OVERALL_MAX = -1.E12_EB

! Read the TREE lines to determine how many vegetation volumes there will be

N_TREES = 0
REWIND(LU_INPUT)
COUNT_VEG_LOOP: DO
   CALL CHECKREAD('TREE',LU_INPUT,IOS) 
   IF (IOS==1) EXIT COUNT_VEG_LOOP
   READ(LU_INPUT,NML=TREE,END=11,ERR=12,IOSTAT=IOS)
   N_TREES = N_TREES + 1
   12 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with TREE line')
ENDDO COUNT_VEG_LOOP
11 REWIND(LU_INPUT)
!
! Sequentially read the TREE namelist to get shape and size parameters for each
! vegetation volume.
!
ALLOCATE(TREE_MESH(NMESHES),STAT=IZERO)
CALL ChkMemErr('READ','TREE_MESH',IZERO)
TREE_MESH = .TRUE.

IF (N_TREES==0) THEN
 TREE_MESH = .FALSE.
 RETURN
ENDIF
!
ALLOCATE(VEG_FUEL_GEOM(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','VEG_FUEL_GEOM',IZERO)
ALLOCATE(VEG_LABELS(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','VEG_LABELS',IZERO)
ALLOCATE(MESH_FOR_TREE(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','MESH_FOR_TREE',IZERO) ; MESH_FOR_TREE = 0
ALLOCATE(CROWN_W(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','CROWN_W',IZERO)
ALLOCATE(CROWN_W_TOP(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','CROWN_W_TOP',IZERO)
ALLOCATE(CROWN_W_BOTTOM(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','CROWN_W_BOTTOM',IZERO)
ALLOCATE(CROWN_B_H(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','CROWN_B_H',IZERO)
ALLOCATE(TREE_H(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','TREE_H',IZERO)
ALLOCATE(X_TREE(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','X_TREE',IZERO)
ALLOCATE(Y_TREE(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','Y_TREE',IZERO)
ALLOCATE(Z_TREE(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','Z_TREE',IZERO)
ALLOCATE(XS_RECT_VEG(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','XS_RECT_VEG',IZERO)
ALLOCATE(XF_RECT_VEG(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','XF_RECT_VEG',IZERO)
ALLOCATE(YS_RECT_VEG(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','YS_RECT_VEG',IZERO)
ALLOCATE(YF_RECT_VEG(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','YF_RECT_VEG',IZERO)
ALLOCATE(ZS_RECT_VEG(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','ZS_RECT_VEG',IZERO)
ALLOCATE(ZF_RECT_VEG(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','ZF_RECT_VEG',IZERO)
ALLOCATE(TREE_PARTICLE_CLASS(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','TREE_PARTICLE_CLASS',IZERO)

ALLOCATE(RING_THICKNESS_VEG(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','RING_THICKNESS_VEG',IZERO)
ALLOCATE(IGN_ELEM(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','IGN_ELEM',IZERO) ; IGN_ELEM = .FALSE.
ALLOCATE(TON_IGN_ELEMS(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','TON_IGN_ELEMS',IZERO)
ALLOCATE(TOFF_IGN_ELEMS(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','TOFF_IGN_ELEMS',IZERO)
ALLOCATE(T_RAMPOFF_IGN_ELEMS(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','T_RAMPOFF_IGN_ELEMS',IZERO)
ALLOCATE(T_RAMPON_IGN_ELEMS(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','T_RAMPON_IGN_ELEMS',IZERO)
!
N_TREES_0 = N_TREES
N = 0
!
VEG_LOOP: DO NN=1,N_TREES_0
   N = N + 1
!
   PART_ID   = 'null'
   LABEL = 'null'
   FUEL_GEOM = 'null'
   IGNITOR_ELEMENTS = .FALSE.
!
   CALL CHECKREAD('TREE',LU_INPUT,IOS)
   IF (IOS==1) EXIT VEG_LOOP
   READ(LU_INPUT,TREE,END=25,IOSTAT=IOS)
   IF (PART_ID  =='null') CALL SHUTDOWN('ERROR: Specify PART_ID in TREE namelist')
   IF (FUEL_GEOM=='null') CALL SHUTDOWN('ERROR: Specify FUEL_GEOM in TREE namelist')

! Find min and max extent of each vegetation volume to help identify meshes that
! contain vegetation
   IF (FUEL_GEOM == 'CONE' .OR. FUEL_GEOM == 'CYLINDER' .OR. FUEL_GEOM == 'RING') THEN
    X_TREE_MIN = XYZ(1) - 0.5*CROWN_WIDTH
    X_TREE_MAX = XYZ(1) + 0.5*CROWN_WIDTH
    Y_TREE_MIN = XYZ(2) - 0.5*CROWN_WIDTH
    Y_TREE_MAX = XYZ(2) + 0.5*CROWN_WIDTH
    Z_TREE_MIN = XYZ(3) + CROWN_BASE_HEIGHT
    Z_TREE_MAX = XYZ(3) + TREE_HEIGHT
   ENDIF

   IF (FUEL_GEOM == 'FRUSTUM') THEN
    X_TREE_MIN = MIN(XYZ(1) - 0.5*CROWN_WIDTH_BOTTOM,XYZ(1) - 0.5*CROWN_WIDTH_TOP)
    X_TREE_MAX = MIN(XYZ(1) + 0.5*CROWN_WIDTH_BOTTOM,XYZ(1) + 0.5*CROWN_WIDTH_TOP)
    Y_TREE_MIN = MIN(XYZ(2) - 0.5*CROWN_WIDTH_BOTTOM,XYZ(2) - 0.5*CROWN_WIDTH_TOP)
    Y_TREE_MAX = MIN(XYZ(2) + 0.5*CROWN_WIDTH_BOTTOM,XYZ(2) + 0.5*CROWN_WIDTH_TOP)
    Z_TREE_MIN = XYZ(3) + CROWN_BASE_HEIGHT
    Z_TREE_MAX = XYZ(3) + TREE_HEIGHT
   ENDIF

   IF (FUEL_GEOM == 'RECTANGLE') THEN
    X_TREE_MIN = XB(1)
    X_TREE_MAX = XB(2)
    Y_TREE_MIN = XB(3)
    Y_TREE_MAX = XB(4)
    Z_TREE_MIN = XB(5)
    Z_TREE_MAX = XB(6)
   ENDIF

   IF (X_TREE_MIN < X_OVERALL_MIN) X_OVERALL_MIN = X_TREE_MIN
   IF (X_TREE_MAX > X_OVERALL_MAX) X_OVERALL_MAX = X_TREE_MAX
   IF (Y_TREE_MIN < Y_OVERALL_MIN) Y_OVERALL_MIN = Y_TREE_MIN
   IF (Y_TREE_MAX > Y_OVERALL_MAX) Y_OVERALL_MAX = Y_TREE_MAX
   IF (Z_TREE_MIN < Z_OVERALL_MIN) Z_OVERALL_MIN = Z_TREE_MIN
   IF (Z_TREE_MAX > Z_OVERALL_MAX) Z_OVERALL_MAX = Z_TREE_MAX
!
   VEG_FUEL_GEOM(N)  = FUEL_GEOM
   CROWN_W(N)        = CROWN_WIDTH
   CROWN_W_BOTTOM(N) = CROWN_WIDTH_BOTTOM
   CROWN_W_TOP(N)    = CROWN_WIDTH_TOP
   CROWN_B_H(N)      = CROWN_BASE_HEIGHT
   TREE_H(N)         = TREE_HEIGHT
   X_TREE(N) = XYZ(1)
   Y_TREE(N) = XYZ(2)
   Z_TREE(N) = XYZ(3)
!
   XS_RECT_VEG(N) = XB(1)
   XF_RECT_VEG(N) = XB(2)
   YS_RECT_VEG(N) = XB(3)
   YF_RECT_VEG(N) = XB(4)
   ZS_RECT_VEG(N) = XB(5)
   ZF_RECT_VEG(N) = XB(6)
!
   VEG_LABELS(N) = LABEL
   RING_THICKNESS_VEG(N)  = RING_THICKNESS
   IGN_ELEM(N)            = IGNITOR_ELEMENTS
   TON_IGN_ELEMS(N)       = TON_IGNITOR_ELEMENTS
   TOFF_IGN_ELEMS(N)      = TOFF_IGNITOR_ELEMENTS
   T_RAMPON_IGN_ELEMS(N)  = T_RAMPON_IGNITOR_ELEMENTS
   T_RAMPOFF_IGN_ELEMS(N) = T_RAMPOFF_IGNITOR_ELEMENTS
!
   DO IPC=1,N_PART
      PC=>PARTICLE_CLASS(IPC)
      IF (PC%ID==PART_ID) TREE_PARTICLE_CLASS(N) = IPC
   ENDDO
   DROPLET_FILE=.TRUE.
!
ENDDO VEG_LOOP
25 REWIND(LU_INPUT)
!
! Find meshes that do not have vegetation
DO NM = 1, NMESHES
 IF (MESHES(NM)%XF < X_OVERALL_MIN) TREE_MESH(NM) = .FALSE. 
 IF (MESHES(NM)%XS > X_OVERALL_MAX) TREE_MESH(NM) = .FALSE. 
 IF (MESHES(NM)%YF < Y_OVERALL_MIN) TREE_MESH(NM) = .FALSE. 
 IF (MESHES(NM)%YS > Y_OVERALL_MAX) TREE_MESH(NM) = .FALSE. 
 IF (MESHES(NM)%ZF < Z_OVERALL_MIN) TREE_MESH(NM) = .FALSE. 
 IF (MESHES(NM)%ZS > Z_OVERALL_MAX) TREE_MESH(NM) = .FALSE. 
ENDDO
!
END SUBROUTINE READ_TREE
 
 
SUBROUTINE READ_PROP
USE DEVICE_VARIABLES
USE MATH_FUNCTIONS, ONLY : GET_RAMP_INDEX,GET_TABLE_INDEX
REAL(EB) :: ACTIVATION_OBSCURATION,ACTIVATION_TEMPERATURE,ALPHA_C,ALPHA_E,BETA_C,BETA_E, &
            BEAD_DIAMETER,BEAD_EMISSIVITY, &
            CONDUIT_DIAMETER,CONDUIT_THICKNESS,CABLE_MASS_PER_LENGTH,CABLE_DIAMETER, &
            CABLE_JACKET_THICKNESS,CABLE_FAILURE_TEMPERATURE, &
            C_FACTOR,CHARACTERISTIC_VELOCITY,ORIFICE_DIAMETER, &
            DROPLET_VELOCITY,FLOW_RATE,FLOW_TAU,GAUGE_TEMPERATURE,INITIAL_TEMPERATURE,K_FACTOR,LENGTH,SPRAY_ANGLE(2), &
            OFFSET,OPERATING_PRESSURE,RTI,PDPA_START,PDPA_END,PDPA_RADIUS
INTEGER :: N, PDPA_M, PDPA_N
EQUIVALENCE(LENGTH,ALPHA_C)
CHARACTER(30) :: SMOKEVIEW_ID,QUANTITY,PART_ID,FLOW_RAMP,SPRAY_PATTERN_TABLE,SPEC_ID,PRESSURE_RAMP
TYPE (PROPERTY_TYPE), POINTER :: PY

NAMELIST /PROP/ ACTIVATION_OBSCURATION,ACTIVATION_TEMPERATURE,ALPHA_C,ALPHA_E,BETA_C,BETA_E, &
                BEAD_DIAMETER,BEAD_EMISSIVITY, &
                CONDUIT_DIAMETER,CONDUIT_THICKNESS,CABLE_MASS_PER_LENGTH,CABLE_DIAMETER, &
                CABLE_JACKET_THICKNESS,CABLE_FAILURE_TEMPERATURE, &
                C_FACTOR,CHARACTERISTIC_VELOCITY,ORIFICE_DIAMETER, &
                DROPLET_VELOCITY,FLOW_RATE,FLOW_RAMP,FLOW_TAU,ID,GAUGE_TEMPERATURE,INITIAL_TEMPERATURE,K_FACTOR,LENGTH,OFFSET, &
                OPERATING_PRESSURE,PART_ID,QUANTITY,RTI,SPRAY_ANGLE,SMOKEVIEW_ID,SPEC_ID,SPRAY_PATTERN_TABLE,PRESSURE_RAMP, &
                PDPA_START,PDPA_END,PDPA_RADIUS,PDPA_M,PDPA_N

! Count the PROP lines in the input file. Note how many of these are cables.

N_PROP=0
REWIND(LU_INPUT)
COUNT_PROP_LOOP: DO
   CALL CHECKREAD('PROP',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_PROP_LOOP
   READ(LU_INPUT,PROP,ERR=34,IOSTAT=IOS)
   N_PROP = N_PROP + 1
   IF (QUANTITY=='CABLE TEMPERATURE') N_CABL = N_CABL + 1
   34 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with PROP number', N_PROP+1
         CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_PROP_LOOP
 
! Allocate the PROPERTY derived types
 
ALLOCATE(PROPERTY(0:N_PROP),STAT=IZERO)
CALL ChkMemErr('READ','PROPERTY',IZERO) 

IF (N_CABL>0) THEN
   ALLOCATE(CABLE(N_CABL),STAT=IZERO)
   CALL ChkMemErr('READ','CABLE',IZERO) 
   N_CABL = 0
ENDIF
 
! Read the PROP lines in the order listed in the input file
 
REWIND(LU_INPUT)

READ_PROP_LOOP: DO N=0,N_PROP
 
   CALL CHECKREAD('PROP',LU_INPUT,IOS)  ! Look for PROP lines in the input file
   CALL SET_PROP_DEFAULTS          ! Reset PROP NAMELIST parameters to default values 
   IF (N > 0) READ(LU_INPUT,PROP) 

   ! Pack PROP parameters into the appropriate property derived types
   PY => PROPERTY(N)
   PY%ACTIVATION_OBSCURATION   = ACTIVATION_OBSCURATION
   PY%ACTIVATION_TEMPERATURE   = ACTIVATION_TEMPERATURE   ! NOTE: Act_Temp remains in degrees C. It is just a SETPOINT.
   PY%ALPHA_C                  = ALPHA_C
   PY%ALPHA_E                  = ALPHA_E
   PY%BETA_C                   = BETA_C
   PY%BETA_E                   = BETA_E
   PY%BEAD_DIAMETER            = BEAD_DIAMETER
   PY%BEAD_EMISSIVITY          = BEAD_EMISSIVITY
   PY%CABLE_DIAMETER           = CABLE_DIAMETER
   PY%CABLE_JACKET_THICKNESS   = CABLE_JACKET_THICKNESS
   PY%CABLE_MASS_PER_LENGTH    = CABLE_MASS_PER_LENGTH
   PY%CABLE_FAILURE_TEMPERATURE= CABLE_FAILURE_TEMPERATURE
   PY%CONDUIT_DIAMETER         = CONDUIT_DIAMETER
   PY%CONDUIT_THICKNESS        = CONDUIT_THICKNESS
   IF (CONDUIT_DIAMETER>0._EB) CONDUIT = .TRUE.
   PY%C_FACTOR                 = C_FACTOR
   PY%CHARACTERISTIC_VELOCITY  = CHARACTERISTIC_VELOCITY
   PY%GAUGE_TEMPERATURE        = GAUGE_TEMPERATURE + TMPM
   PY%ID                       = ID
   PY%INITIAL_TEMPERATURE      = INITIAL_TEMPERATURE + TMPM
   PY%OFFSET                   = OFFSET
   PY%OPERATING_PRESSURE       = OPERATING_PRESSURE
   PY%PART_ID                  = PART_ID
   PY%QUANTITY                 = QUANTITY
   IF (PY%PART_ID/='null' .AND. PY%QUANTITY == 'null' ) PY%QUANTITY = 'NOZZLE'
   PY%RTI                      = RTI
   IF (SMOKEVIEW_ID /= 'null') THEN
      PY%SMOKEVIEW_ID          = SMOKEVIEW_ID
   ELSE
      SELECT CASE(PY%QUANTITY)
         CASE DEFAULT
            PY%SMOKEVIEW_ID = 'sensor'
         CASE('SPRINKLER LINK TEMPERATURE')
            PY%SMOKEVIEW_ID = 'sprinkler_pendent'
         CASE('NOZZLE')
            PY%SMOKEVIEW_ID = 'nozzle'
         CASE('LINK TEMPERATURE')
            PY%SMOKEVIEW_ID = 'heat_detector'
         CASE('CABLE TEMPERATURE')
            PY%SMOKEVIEW_ID = 'sensor'
         CASE('spot obscuration','CHAMBER OBSCURATION')
            PY%SMOKEVIEW_ID = 'smoke_detector'
      END SELECT
   ENDIF
   PY%SPEC_ID                  = SPEC_ID
   IF (PART_ID/='null' .AND. SPRAY_PATTERN_TABLE /= 'null') THEN
      CALL GET_TABLE_INDEX(SPRAY_PATTERN_TABLE,SPRAY_PATTERN,PY%SPRAY_PATTERN_INDEX)
      PY%TABLE_ID = SPRAY_PATTERN_TABLE
   ELSE
      PY%SPRAY_PATTERN_INDEX = 0
   ENDIF 
   PY%SPRAY_ANGLE(1)           = SPRAY_ANGLE(1)*PI/180._EB
   PY%SPRAY_ANGLE(2)           = SPRAY_ANGLE(2)*PI/180._EB

   IF (QUANTITY=='CABLE TEMPERATURE') THEN
      N_CABL = N_CABL + 1
      CABLE(N_CABL)%PROP_INDEX = N
   ENDIF

   PY%PDPA_START       = PDPA_START
   PY%PDPA_END         = PDPA_END
   PY%PDPA_RADIUS      = PDPA_RADIUS
   PY%PDPA_M           = PDPA_M
   PY%PDPA_N           = PDPA_N

   ! Set flow variables

   IF (PART_ID/='null' .AND. PRESSURE_RAMP /= 'null') THEN
      CALL GET_RAMP_INDEX(PRESSURE_RAMP,'PRESSURE',PY%PRESSURE_RAMP_INDEX)
   ELSE
      PY%PRESSURE_RAMP_INDEX = 0
   ENDIF

   ! Check sufficient input

   IF (PY%PRESSURE_RAMP_INDEX == 0 .AND. FLOW_RATE > 0._EB) THEN
      IF (K_FACTOR < 0._EB) K_FACTOR = 10.0_EB
   ENDIF

   IF (PART_ID /='null') THEN
      IF ((FLOW_RATE>0._EB .AND. K_FACTOR<=0._EB .AND. OPERATING_PRESSURE<=0._EB) .OR. &
          (FLOW_RATE<0._EB .AND. K_FACTOR>=0._EB .AND. OPERATING_PRESSURE<=0._EB) .OR. &
          (FLOW_RATE<0._EB .AND. K_FACTOR<=0._EB .AND. OPERATING_PRESSURE>0._EB)) THEN
         WRITE(MESSAGE,'(A,A,A)') 'Problem with PROP ',TRIM(PY%ID),', too few flow parameters'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      IF (K_FACTOR < 0._EB)           K_FACTOR           = FLOW_RATE/SQRT(OPERATING_PRESSURE)
      IF (FLOW_RATE < 0._EB)          FLOW_RATE          = K_FACTOR*SQRT(OPERATING_PRESSURE)
      IF (OPERATING_PRESSURE < 0._EB) OPERATING_PRESSURE = (FLOW_RATE/K_FACTOR)**2
      PY%K_FACTOR           = K_FACTOR
      PY%FLOW_RATE          = FLOW_RATE
      PY%OPERATING_PRESSURE = OPERATING_PRESSURE

      IF ((DROPLET_VELOCITY == 0._EB) .AND. (ORIFICE_DIAMETER == 0.) .AND. (PRESSURE_RAMP == 'null')) THEN
         WRITE(MESSAGE,'(A,A,A)') 'Warning: PROP ',TRIM(PY%ID),': droplet velocity is not defined.'
         IF (MYID==0) WRITE(LU_ERR,'(A)') TRIM(MESSAGE)
      ENDIF
      
      IF (DROPLET_VELOCITY > 0._EB) THEN
         PY%DROPLET_VELOCITY  = DROPLET_VELOCITY
      ELSEIF ((ORIFICE_DIAMETER > 0._EB) .AND. (FLOW_RATE > 0._EB)) THEN
         PY%DROPLET_VELOCITY  = (FLOW_RATE/60._EB/1000._EB)/(PI*(ORIFICE_DIAMETER/2._EB)**2)
      ENDIF

   ENDIF

   IF (FLOW_RAMP /= 'null') THEN
      CALL GET_RAMP_INDEX(FLOW_RAMP,'TIME',PY%FLOW_RAMP_INDEX)
   ELSE
      PY%FLOW_RAMP_INDEX = 0
   ENDIF 
   IF (FLOW_TAU /= 0._EB) THEN
      PY%FLOW_TAU = FLOW_TAU 
      IF (FLOW_TAU > 0._EB) PY%FLOW_RAMP_INDEX = TANH_RAMP 
      IF (FLOW_TAU < 0._EB) PY%FLOW_RAMP_INDEX = TSQR_RAMP
   ENDIF 

ENDDO READ_PROP_LOOP
 
 
CONTAINS

 
SUBROUTINE SET_PROP_DEFAULTS
 
ACTIVATION_OBSCURATION   = 3.28_EB     ! %/m
ACTIVATION_TEMPERATURE   = -273.15_EB  ! C
ALPHA_C                  = 1.8_EB      ! m, Heskestad Length Scale
ALPHA_E                  = 0.0_EB
BETA_C                   = -1.0_EB
BETA_E                   = -1.0_EB
BEAD_DIAMETER            = 0.001       ! m
BEAD_EMISSIVITY          = 0.85_EB
CABLE_FAILURE_TEMPERATURE= 400._EB     ! C
CABLE_DIAMETER           = 0.02_EB     ! m
CABLE_JACKET_THICKNESS   = 0.002_EB    ! m
CABLE_MASS_PER_LENGTH    = 0.5_EB      ! kg/m
CONDUIT_DIAMETER         = 0.0         ! m
CONDUIT_THICKNESS        = 0.0         ! m
C_FACTOR                 = 0.0_EB
CHARACTERISTIC_VELOCITY  = 1.0_EB      ! m/s
DROPLET_VELOCITY         = 0._EB       ! m/s
FLOW_RATE                = -1._EB      ! L/min
FLOW_RAMP                = 'null'
FLOW_TAU                 = 0._EB
GAUGE_TEMPERATURE        = TMPA - TMPM
INITIAL_TEMPERATURE      = TMPA - TMPM
ID                       = 'null'
K_FACTOR                 = -1.0_EB     ! L/min/atm**0.5
OFFSET                   = 0.05_EB     ! m
OPERATING_PRESSURE       = -1.0_EB     ! atm
ORIFICE_DIAMETER         = 0.0_EB      ! m
PART_ID                  = 'null'
PDPA_START               = T_BEGIN
PDPA_END                 = T_END + 1.0_EB
PDPA_RADIUS              = 0.1_EB
PDPA_M                   = 0
PDPA_N                   = 0
PRESSURE_RAMP           = 'null'
QUANTITY                 = 'null'
RTI                      = 100._EB     ! (ms)**0.5
SMOKEVIEW_ID             = 'null' 
SPEC_ID                  = 'null'
SPRAY_ANGLE(1)           = 60._EB      ! degrees
SPRAY_ANGLE(2)           = 75._EB      ! degrees
SPRAY_PATTERN_TABLE      = 'null'
END SUBROUTINE SET_PROP_DEFAULTS
 
END SUBROUTINE READ_PROP
 
SUBROUTINE PROC_PROP
USE DEVICE_VARIABLES
REAL(EB) :: TOTAL_FLOWRATE, SUBTOTAL_FLOWRATE
INTEGER :: N,NN,N_V_FACTORS
LOGICAL :: TABLE_NORMED(1:N_TABLE)
TYPE (PROPERTY_TYPE), POINTER :: PY
TYPE (TABLES_TYPE),  POINTER :: TA

TABLE_NORMED = .FALSE.

PROP_LOOP: DO N=0,N_PROP
   PY => PROPERTY(N)

   ! Set up spinkler distributrion if needed
   IF (PY%SPRAY_PATTERN_INDEX > 0) THEN
      TA => TABLES(PY%SPRAY_PATTERN_INDEX)
      ALLOCATE(PY%TABLE_ROW(1:TA%NUMBER_ROWS))
      TOTAL_FLOWRATE=0._EB
      SUBTOTAL_FLOWRATE=0._EB
      DO NN=1,TA%NUMBER_ROWS
         IF (TA%TABLE_DATA(NN,6) <=0._EB) THEN
            WRITE(MESSAGE,'(A,A,A,I5)') 'Spray Pattern Table, ',TRIM(PY%TABLE_ID),', massflux <= 0 for line ',NN
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         TOTAL_FLOWRATE = TOTAL_FLOWRATE + TA%TABLE_DATA(NN,6)
      ENDDO
      IF (TABLE_NORMED(PY%SPRAY_PATTERN_INDEX)) THEN
         DO NN=1,TA%NUMBER_ROWS
            SUBTOTAL_FLOWRATE = SUBTOTAL_FLOWRATE + TA%TABLE_DATA(NN,6)
            PY%TABLE_ROW(NN) = SUBTOTAL_FLOWRATE/TOTAL_FLOWRATE
         ENDDO
      ELSE
         DO NN=1,TA%NUMBER_ROWS
            TA%TABLE_DATA(NN,1) = TA%TABLE_DATA(NN,1) * PI/180._EB
            TA%TABLE_DATA(NN,2) = TA%TABLE_DATA(NN,2) * PI/180._EB
            TA%TABLE_DATA(NN,3) = TA%TABLE_DATA(NN,3) * PI/180._EB
            TA%TABLE_DATA(NN,4) = TA%TABLE_DATA(NN,4) * PI/180._EB
            SUBTOTAL_FLOWRATE = SUBTOTAL_FLOWRATE + TA%TABLE_DATA(NN,6)
            PY%TABLE_ROW(NN) = SUBTOTAL_FLOWRATE/TOTAL_FLOWRATE
         ENDDO
         TABLE_NORMED(PY%SPRAY_PATTERN_INDEX) = .TRUE.
      ENDIF
      PY%TABLE_ROW(TA%NUMBER_ROWS) = 1._EB
   END IF

   ! Set up pressure dependence
   IF (PY%PRESSURE_RAMP_INDEX > 0) THEN
      IF (PY%SPRAY_PATTERN_INDEX > 0) THEN
         N_V_FACTORS = TA%NUMBER_ROWS
      ELSE
         N_V_FACTORS = 1
      ENDIF
      ALLOCATE(PY%V_FACTOR(1:N_V_FACTORS))
      IF (PY%SPRAY_PATTERN_INDEX > 0) THEN
         DO NN=1,TA%NUMBER_ROWS
            PY%V_FACTOR(NN) = TA%TABLE_DATA(NN,5)/SQRT(PY%OPERATING_PRESSURE)
         ENDDO
      ELSE
         PY%V_FACTOR = PY%DROPLET_VELOCITY/SQRT(PY%OPERATING_PRESSURE)
      ENDIF
   ENDIF
ENDDO PROP_LOOP

END SUBROUTINE PROC_PROP



SUBROUTINE READ_MATL

USE MATH_FUNCTIONS, ONLY : GET_RAMP_INDEX
USE DEVICE_VARIABLES, ONLY : CABLE,N_CABL,PROPERTY,CONDUIT
CHARACTER(30) :: CONDUCTIVITY_RAMP,SPECIFIC_HEAT_RAMP
REAL(EB) :: EMISSIVITY,CONDUCTIVITY,SPECIFIC_HEAT,DENSITY,ABSORPTION_COEFFICIENT,BOILING_TEMPERATURE, &
            INITIAL_VAPOR_FLUX
REAL(EB), DIMENSION(1:MAX_REACTIONS) :: A,E,HEAT_OF_REACTION,NU_FUEL,NU_WATER,NU_RESIDUE,N_S,N_T,&
                                        REFERENCE_RATE,REFERENCE_TEMPERATURE,THRESHOLD_TEMPERATURE
REAL(EB), DIMENSION(1:MAX_REACTIONS,1:MAX_SPECIES) :: NU_GAS,HEAT_OF_COMBUSTION
CHARACTER(30), DIMENSION(1:MAX_REACTIONS) :: RESIDUE
INTEGER :: N,NN,NNN,IOS,NR,N_REACTIONS,N_CONDUIT
NAMELIST /MATL/ ID,FYI,SPECIFIC_HEAT,CONDUCTIVITY,CONDUCTIVITY_RAMP,SPECIFIC_HEAT_RAMP, &
                REFERENCE_TEMPERATURE, REFERENCE_RATE, THRESHOLD_TEMPERATURE, &
                EMISSIVITY,HEAT_OF_REACTION,DENSITY,RESIDUE, &
                HEAT_OF_COMBUSTION,A,E,NU_FUEL,NU_WATER,NU_RESIDUE,NU_GAS,N_S,N_T,N_REACTIONS,&
                ABSORPTION_COEFFICIENT,BOILING_TEMPERATURE,INITIAL_VAPOR_FLUX

! Count the MATL lines in the input file
 
REWIND(LU_INPUT)
N_MATL = 0
COUNT_MATL_LOOP: DO
   CALL CHECKREAD('MATL',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_MATL_LOOP
   READ(LU_INPUT,MATL,ERR=34,IOSTAT=IOS)
   N_MATL = N_MATL + 1
   MATL_NAME(N_MATL) = ID
   34 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with MATL number', N_MATL+1
         CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_MATL_LOOP

N_CONDUIT = 0
IF (CONDUIT) N_CONDUIT = 2
 
! Allocate the MATERIAL derived type
 
ALLOCATE(MATERIAL(1:N_MATL+N_CABL+N_CONDUIT),STAT=IZERO)
CALL ChkMemErr('READ','MATERIAL',IZERO) 

! Add an additional MATL entry for each CABLE

DO N=1,N_CABL
   N_MATL = N_MATL+1
   MATL_NAME(N_MATL) = PROPERTY(CABLE(N)%PROP_INDEX)%ID
   MATERIAL(N_MATL)%USER_DEFINED = .FALSE.
   MATERIAL(N_MATL)%PROP_INDEX   = CABLE(N)%PROP_INDEX
   MATERIAL(N_MATL)%CABL_INDEX   = N
ENDDO

IF (CONDUIT) THEN
   N_MATL = N_MATL+1
   MATL_NAME(N_MATL) = 'STEEL CONDUIT'
   MATERIAL(N_MATL)%USER_DEFINED = .FALSE.
   MATERIAL(N_MATL)%CONDUIT      = .TRUE.
   N_MATL = N_MATL+1
   MATL_NAME(N_MATL) = 'AIR'
   MATERIAL(N_MATL)%USER_DEFINED = .FALSE.
   MATERIAL(N_MATL)%AIR          = .TRUE.
ENDIF
 
! Read the MATL lines in the order listed in the input file
 
REWIND(LU_INPUT)
 
READ_MATL_LOOP: DO N=1,N_MATL
 
   ML => MATERIAL(N)

   ! Read user defined MATL lines

   IF (ML%USER_DEFINED) THEN
      CALL CHECKREAD('MATL',LU_INPUT,IOS)
      CALL SET_MATL_DEFAULTS
      READ(LU_INPUT,MATL) 
   ENDIF

   ! Define cable properties for THIEF model

   IF (ML%CABL_INDEX>0) THEN
      CALL SET_MATL_DEFAULTS
      ID            = MATL_NAME(N)
      DENSITY       = PROPERTY(ML%PROP_INDEX)%CABLE_MASS_PER_LENGTH/(PI*(0.5_EB*PROPERTY(ML%PROP_INDEX)%CABLE_DIAMETER)**2)
      CONDUCTIVITY  = 0.2_EB
      SPECIFIC_HEAT = 1.5_EB
      EMISSIVITY    = 0.95_EB
   ENDIF

   IF (ML%CONDUIT) THEN
      CALL SET_MATL_DEFAULTS
      ID            = MATL_NAME(N)
      DENSITY       = 7850.
      CONDUCTIVITY  = 48._EB
      SPECIFIC_HEAT = 0.46_EB
      EMISSIVITY    = 0.85_EB
   ENDIF

   IF (ML%AIR) THEN
      CALL SET_MATL_DEFAULTS
      ID            = MATL_NAME(N)
      DENSITY       = 1.2_EB
      CONDUCTIVITY  = 0.026_EB
      SPECIFIC_HEAT = 1.012_EB
      ABSORPTION_COEFFICIENT = 0._EB
   ENDIF

   ! Do some error checking on the inputs

   NOT_BOILING: IF (BOILING_TEMPERATURE>4000._EB) THEN

      IF ( ( ANY(THRESHOLD_TEMPERATURE>-TMPM) .OR. ANY(REFERENCE_TEMPERATURE>-TMPM) .OR. ANY(E>=0._EB) .OR. &
             ANY(HEAT_OF_REACTION/=0._EB) ) .AND. N_REACTIONS==0) THEN
         WRITE(MESSAGE,'(A,I2,A)') 'ERROR: Problem with MATL number ',N,'. A reaction parameter is used, but N_REACTIONS=0'  
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      DO NN=1,N_REACTIONS
         IF (REFERENCE_TEMPERATURE(NN)<-TMPM  .AND. E(NN)< 0._EB) THEN
            WRITE(MESSAGE,'(A,I2,A,I2,A)') 'ERROR: Problem with MATL ',N,', REACTION ',NN,'. Set REFERENCE_TEMPERATURE or E'  
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         IF (NU_FUEL(NN)==0._EB .AND. NU_WATER(NN)==0._EB .AND. NU_RESIDUE(NN)==0._EB .AND. SUM(NU_GAS(NN,:))==0._EB) THEN
            WRITE(MESSAGE,'(A,I2,A,I2,A)') 'WARNING: MATL ',N,', REACTION ',NN,'. No product yields (NUs) set'  
            IF (MYID==0) WRITE(LU_ERR,'(A)') TRIM(MESSAGE)
         ENDIF
      ENDDO
   ELSE NOT_BOILING ! Is liquid
      IF (HEAT_OF_REACTION(1) == 0._EB) THEN
         WRITE(MESSAGE,'(A,I2,A)') 'ERROR: HEAT_OF_REACTION should be greater than zero for liquid MATL ',N,'.'  
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDIF NOT_BOILING

   ! Error checking for thermal properties
   IF ( DENSITY == 0._EB ) THEN
      WRITE(MESSAGE,'(A,I2,A)') 'ERROR: Problem with MATL number ',N,': DENSITY=0' 
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF ( CONDUCTIVITY == 0._EB .AND. CONDUCTIVITY_RAMP == 'null' ) THEN
      WRITE(MESSAGE,'(A,I2,A)') 'ERROR: Problem with MATL number ',N,': CONDUCTIVITY = 0' 
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF ( SPECIFIC_HEAT == 0._EB .AND. SPECIFIC_HEAT_RAMP == 'null' ) THEN
      WRITE(MESSAGE,'(A,I2,A)') 'ERROR: Problem with MATL number ',N,': SPECIFIC_HEAT = 0' 
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF ( SPECIFIC_HEAT > 10._EB) WRITE(LU_ERR,'(A,I2)') 'WARNING: SPECIFIC_HEAT units are kJ/kg/K check MATL number ',N
 
   ! Pack MATL parameters into the MATERIAL derived type
 
   ML%A(:)                 = A(:)
   ML%ADJUST_BURN_RATE     = 1._EB
   ML%C_S                  = 1000._EB*SPECIFIC_HEAT/TIME_SHRINK_FACTOR
   ML%E(:)                 = 1000._EB*E(:)
   ML%EMISSIVITY           = EMISSIVITY
   ML%FYI                  = FYI
   ML%HEAT_OF_COMBUSTION   = 1000._EB*HEAT_OF_COMBUSTION
   ML%H_R(:)               = 1000._EB*HEAT_OF_REACTION(:)
   ML%INIT_VAPOR_FLUX      = INITIAL_VAPOR_FLUX
   ML%KAPPA_S              = ABSORPTION_COEFFICIENT
   ML%K_S                  = CONDUCTIVITY
   ML%N_REACTIONS          = N_REACTIONS
   ML%N_S(:)               = N_S(:)
   ML%N_T(:)               = N_T(:)
   ML%NU_FUEL(:)           = NU_FUEL(:)
   ML%NU_RESIDUE(:)        = NU_RESIDUE(:)
   ML%NU_WATER(:)          = NU_WATER(:)
   ML%NU_GAS(:,:)          = NU_GAS(:,:)
   ML%RAMP_C_S             = SPECIFIC_HEAT_RAMP
   ML%RAMP_K_S             = CONDUCTIVITY_RAMP
   ML%RHO_S                = DENSITY
   ML%RESIDUE_MATL_NAME(:) = RESIDUE(:)
   ML%TMP_BOIL             = BOILING_TEMPERATURE + TMPM
   ML%TMP_THR(:)           = THRESHOLD_TEMPERATURE(:) + TMPM
   ML%TMP_REF(:)           = REFERENCE_TEMPERATURE(:) + TMPM
 
   ! Additional logic

   IF (BOILING_TEMPERATURE<5000._EB) THEN
      ML%PYROLYSIS_MODEL = PYROLYSIS_LIQUID
      ML%N_REACTIONS = 1
      IF (ML%NU_FUEL(1)==0._EB) ML%NU_FUEL(1) = 1._EB
   ELSE
      ML%PYROLYSIS_MODEL = PYROLYSIS_SOLID
      IF (N_REACTIONS==0) ML%PYROLYSIS_MODEL = PYROLYSIS_NONE
   ENDIF

   DO NN=1,ML%N_REACTIONS
      IF (NU_FUEL(NN) > 0._EB) MIXTURE_FRACTION = .TRUE.
   ENDDO

   IF (ML%RAMP_K_S/='null') THEN
      CALL GET_RAMP_INDEX(ML%RAMP_K_S,'TEMPERATURE',NR)
      ML%K_S = -NR
   ENDIF

   IF (ML%RAMP_C_S/='null') THEN
      CALL GET_RAMP_INDEX(ML%RAMP_C_S,'TEMPERATURE',NR)
      ML%C_S = -NR
   ENDIF

   DO NN=1,ML%N_REACTIONS
      IF (ML%TMP_REF(NN) > 0._EB  .AND. ML%E(NN)< 0._EB) ML%E(NN) = -LOG(REFERENCE_RATE(NN)/ML%A(NN))*R0*ML%TMP_REF(NN)
   ENDDO

ENDDO READ_MATL_LOOP
 
! Assign a material index to the RESIDUEs

DO N=1,N_MATL
   ML => MATERIAL(N)
   ML%RESIDUE_MATL_INDEX = 0
   DO NN=1,ML%N_REACTIONS
      DO NNN=1,N_MATL
         IF (MATL_NAME(NNN)==ML%RESIDUE_MATL_NAME(NN)) ML%RESIDUE_MATL_INDEX(NN) = NNN
      ENDDO
      IF (ML%RESIDUE_MATL_INDEX(NN)==0 .AND. ML%NU_RESIDUE(NN)>0._EB) THEN
         WRITE(MESSAGE,'(5A)') 'ERROR: Residue material ', TRIM(ML%RESIDUE_MATL_NAME(NN)),' of ', &
            TRIM(MATL_NAME(N)), ' is not defined.'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDDO
ENDDO

!Check for duplicate names
IF (N_MATL>1) THEN
   DO N=1,N_MATL-1
      DO NN=N+1,N_MATL
         IF(MATL_NAME(N)==MATL_NAME(NN)) THEN
            WRITE(MESSAGE,'(A,A)') 'Duplicate material name: ',TRIM(MATL_NAME(N))
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO
   ENDDO
ENDIF
CONTAINS
 
SUBROUTINE SET_MATL_DEFAULTS
 
A                      = 1.E13_EB    ! 1/s
ABSORPTION_COEFFICIENT = 5.0E4_EB    ! 1/m, corresponds to 99.3% drop within 1E-4 m distance.
BOILING_TEMPERATURE    = 5000._EB    ! C
CONDUCTIVITY           = 0.0_EB      ! W/m/K
CONDUCTIVITY_RAMP      = 'null'
DENSITY                = 0._EB       ! kg/m3
E                      = -1._EB      ! kJ/kmol
EMISSIVITY             = 0.9_EB
FYI                    = 'null'
HEAT_OF_COMBUSTION     = -1._EB      ! kJ/kg
HEAT_OF_REACTION       = 0._EB       ! kJ/kg
ID                     = 'null'
INITIAL_VAPOR_FLUX     = 5.0E-4_EB   ! m3/m2
THRESHOLD_TEMPERATURE  = -TMPM       ! 0 K
N_REACTIONS            = 0
N_S                    = 1._EB
N_T                    = 0._EB
NU_FUEL                = 0._EB
NU_GAS                 = 0._EB
NU_RESIDUE             = 0._EB
NU_WATER               = 0._EB
REFERENCE_TEMPERATURE  = -1000._EB
REFERENCE_RATE         = 0.001_EB
RESIDUE                = 'null'
SPECIFIC_HEAT          = 0.0_EB      ! kJ/kg/K
SPECIFIC_HEAT_RAMP     = 'null'
 
END SUBROUTINE SET_MATL_DEFAULTS
 
END SUBROUTINE READ_MATL



SUBROUTINE PROC_MATL

USE MATH_FUNCTIONS
INTEGER :: N,J,NS

PROC_MATL_LOOP: DO N=1,N_MATL

   ML => MATERIAL(N)

   ! Adjust burn rate if heat of combustion is different from the gas phase reaction value

   DO J = 1,ML%N_REACTIONS   

      IF (MIXTURE_FRACTION) THEN
         IF (N_REACTIONS > 0) THEN
            IF (CO_PRODUCTION) THEN
               RN => REACTION(2)         
            ELSE
               RN => REACTION(1)
            ENDIF
            IF (ML%HEAT_OF_COMBUSTION(J,I_FUEL)>0._EB .AND. RN%HEAT_OF_COMBUSTION>0._EB)  &
             ML%ADJUST_BURN_RATE(J,I_FUEL) = ML%HEAT_OF_COMBUSTION(J,I_FUEL)/(RN%Y_F_INLET*RN%HEAT_OF_COMBUSTION)
         ENDIF
      ELSEIF (N_REACTIONS>0) THEN
         RN => REACTION(1)
         DO NS = 1,N_SPECIES
            IF (ML%HEAT_OF_COMBUSTION(J,NS)>0._EB .AND. RN%HEAT_OF_COMBUSTION>0._EB)  &
                ML%ADJUST_BURN_RATE(J,NS) = ML%HEAT_OF_COMBUSTION(J,NS)/RN%HEAT_OF_COMBUSTION
         ENDDO
      ENDIF
   ENDDO

   ! Collect species yields to NU_GAS array

   DO J = 1,ML%N_REACTIONS
      IF (MIXTURE_FRACTION) THEN
         IF (I_FUEL /= 0 .AND. ML%NU_FUEL(J)  > 0._EB) ML%NU_GAS(J,I_FUEL)  = ML%NU_FUEL(J)
         IF (I_WATER /= 0.AND. ML%NU_WATER(J) > 0._EB) ML%NU_GAS(J,I_WATER) = ML%NU_WATER(J)
      ELSE
         IF (ML%NU_FUEL(J)==0._EB) THEN
            DO NS = 1,N_SPECIES
               IF (SPECIES(NS)%ISFUEL) ML%NU_FUEL(J) = ML%NU_FUEL(J) + ML%NU_GAS(J,NS)
            ENDDO
         ENDIF
      ENDIF
   ENDDO
   
   IF (ML%RAMP_C_S/='null') THEN
         IF (RAMPS(-INT(ML%C_S))%DEPENDENT_DATA(1) > 10._EB) &
            WRITE(LU_ERR,'(A,I2)') 'WARNING: SPECIFIC_HEAT units are kJ/kg/K check MATL number ',N
   ENDIF   



ENDDO PROC_MATL_LOOP

END SUBROUTINE PROC_MATL


SUBROUTINE READ_SURF
 
USE MATH_FUNCTIONS, ONLY : GET_RAMP_INDEX
USE DEVICE_VARIABLES, ONLY : CABLE,N_CABL,PROPERTY_TYPE,PROPERTY
CHARACTER(30) :: PART_ID,RAMP_MF(0:MAX_SPECIES),RAMP_Q,RAMP_V,RAMP_T,MATL_ID(MAX_LAYERS,MAX_MATERIALS),&
                 PROFILE,BACKING,GEOMETRY,NAME_LIST(MAX_MATERIALS)
LOGICAL :: ADIABATIC,BURN_AWAY,SHRINK,POROUS,FREE_SLIP,NO_SLIP
CHARACTER(60) :: TEXTURE_MAP
CHARACTER(25) :: COLOR
TYPE (PROPERTY_TYPE), POINTER :: PY
REAL(EB) :: TAU_Q,TAU_V,TAU_T,TAU_MF(0:MAX_SPECIES),HRRPUA,MLRPUA,TEXTURE_WIDTH,TEXTURE_HEIGHT,VEL_T(2), &
            E_COEFFICIENT,VOLUME_FLUX,TMP_FRONT,TMP_INNER,THICKNESS(MAX_LAYERS),VEL, &
            MASS_FLUX(0:MAX_SPECIES),MASS_FRACTION(0:MAX_SPECIES), Z0,PLE,CONVECTIVE_HEAT_FLUX,PARTICLE_MASS_FLUX, &
            TRANSPARENCY,EXTERNAL_FLUX,TMP_BACK,MASS_FLUX_TOTAL,STRETCH_FACTOR,&
            MATL_MASS_FRACTION(MAX_LAYERS,MAX_MATERIALS),EMISSIVITY,CELL_SIZE_FACTOR,MAX_PRESSURE,&
            IGNITION_TEMPERATURE,HEAT_OF_VAPORIZATION,REGRID_FACTOR,NET_HEAT_FLUX,LAYER_DIVIDE,SURFACE_DENSITY, &
            ROUGHNESS
INTEGER :: NPPC,N,IOS,NL,NN,NNN,N_LIST,N_LIST2,INDEX_LIST(MAX_MATERIALS_TOTAL),LEAK_PATH(2),DUCT_PATH(2),RGB(3),NR,IL
NAMELIST /SURF/ TMP_FRONT,TMP_INNER,THICKNESS,MASS_FRACTION,VEL,VEL_T,NPPC, &
                E_COEFFICIENT,CONVECTIVE_HEAT_FLUX,TAU_Q,TAU_V,TAU_T,RAMP_Q,RAMP_T,TAU_MF, &
                RAMP_MF,PART_ID,RAMP_V,VOLUME_FLUX, PROFILE,PLE,Z0,ID,MASS_FLUX,PARTICLE_MASS_FLUX, &
                FYI,MATL_ID,BACKING,TMP_BACK,HRRPUA,MLRPUA,SHRINK,CELL_SIZE_FACTOR, &
                TEXTURE_MAP,TEXTURE_WIDTH,TEXTURE_HEIGHT,RGB,TRANSPARENCY, BURN_AWAY,LEAK_PATH,DUCT_PATH,ADIABATIC, &
                EXTERNAL_FLUX,MASS_FLUX_TOTAL,GEOMETRY,STRETCH_FACTOR,MATL_MASS_FRACTION,EMISSIVITY,COLOR,POROUS,MAX_PRESSURE,&
                IGNITION_TEMPERATURE,HEAT_OF_VAPORIZATION,REGRID_FACTOR,NET_HEAT_FLUX,LAYER_DIVIDE,SURFACE_DENSITY, &
                ROUGHNESS,FREE_SLIP,NO_SLIP
             
! Count the SURF lines in the input file

REWIND(LU_INPUT)
N_SURF = 0
COUNT_SURF_LOOP: DO
   CALL CHECKREAD('SURF',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_SURF_LOOP
   READ(LU_INPUT,SURF,ERR=34,IOSTAT=IOS)
   N_SURF = N_SURF + 1
   34 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with SURF number', N_SURF+1
         CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_SURF_LOOP

! Allocate the SURFACE derived type, leaving space for SURF entries not defined explicitly by the user
 
N_SURF_RESERVED = 4
ALLOCATE(SURFACE(0:N_SURF+N_CABL+N_SURF_RESERVED),STAT=IZERO)
CALL ChkMemErr('READ','SURFACE',IZERO) 

! Count the SURF lines in the input file

REWIND(LU_INPUT)
NN = 0
COUNT_SURF_LOOP_AGAIN: DO
   CALL CHECKREAD('SURF',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_SURF_LOOP_AGAIN
   READ(LU_INPUT,SURF)
   NN = NN+1
   SURFACE(NN)%ID = ID
ENDDO COUNT_SURF_LOOP_AGAIN

! Add an additional SURF entry for each CABLE

DO N=1,N_CABL
   N_SURF = N_SURF+1
   SURFACE(N_SURF)%ID           = PROPERTY(CABLE(N)%PROP_INDEX)%ID
   SURFACE(N_SURF)%USER_DEFINED = .FALSE.
   SURFACE(N_SURF)%PROP_INDEX   = CABLE(N)%PROP_INDEX
   SURFACE(N_SURF)%CABL_INDEX   = N
ENDDO

! Add three extra surface types to the list that has already been compiled
 
INERT_SURF_INDEX                   = 0
OPEN_SURF_INDEX                    = N_SURF + 1
MIRROR_SURF_INDEX                  = N_SURF + 2
INTERPOLATED_SURF_INDEX            = N_SURF + 3
PERIODIC_SURF_INDEX                = N_SURF + 4

N_SURF                             = N_SURF + N_SURF_RESERVED

SURFACE(INERT_SURF_INDEX)%ID       = 'INERT'
SURFACE(OPEN_SURF_INDEX)%ID        = 'OPEN'
SURFACE(MIRROR_SURF_INDEX)%ID      = 'MIRROR'
SURFACE(INTERPOLATED_SURF_INDEX)%ID= 'INTERPOLATED'
SURFACE(PERIODIC_SURF_INDEX)%ID    = 'PERIODIC'
SURFACE(N_SURF-N_SURF_RESERVED+1:N_SURF)%USER_DEFINED = .FALSE.
SURFACE(0)%USER_DEFINED               = .FALSE.
 
! Check if SURF_DEFAULT exists

CALL CHECK_SURF_NAME(SURF_DEFAULT,EX)
IF (.NOT.EX) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: SURF_DEFAULT not found'
   CALL SHUTDOWN(MESSAGE)
ENDIF

! Add evacuation boundary type if necessary

CALL CHECK_SURF_NAME(EVAC_SURF_DEFAULT,EX)
IF (.NOT.EX) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: EVAC_SURF_DEFAULT not found'
   CALL SHUTDOWN(MESSAGE)
ENDIF

! Read the SURF lines
 
REWIND(LU_INPUT)
READ_SURF_LOOP: DO N=0,N_SURF
 
   SF => SURFACE(N)
   CALL SET_SURF_DEFAULTS
 
   ! Read the user defined SURF lines

   READ_LOOP: DO
      IF (.NOT.SF%USER_DEFINED) EXIT READ_LOOP
      CALL CHECKREAD('SURF',LU_INPUT,IOS)
      CALL SET_SURF_DEFAULTS
      READ(LU_INPUT,SURF) 
      EXIT READ_LOOP
   ENDDO READ_LOOP

   ! Provide SURF properties for CABLEs (THIEF Model)

   IF (SF%CABL_INDEX>0) THEN
      PY => PROPERTY(SF%PROP_INDEX)
      GEOMETRY       = 'CYLINDRICAL'
      STRETCH_FACTOR = 1._EB
      IF (PY%CONDUIT_THICKNESS>0._EB) THEN
         MATL_ID(1,1) = 'STEEL CONDUIT'
         MATL_ID(2,1) = 'AIR'
         MATL_ID(3,1) = PY%ID
         THICKNESS(1) = PY%CONDUIT_THICKNESS
         THICKNESS(2) = 0.5_EB*(PY%CONDUIT_DIAMETER-PY%CABLE_DIAMETER) - PY%CONDUIT_THICKNESS
         THICKNESS(3) = 0.5_EB*PY%CABLE_DIAMETER
      ELSE
         MATL_ID(1,1)   = PY%ID
         THICKNESS(1)   = 0.5_EB*PY%CABLE_DIAMETER
      ENDIF
   ENDIF

   ! Check SURF parameters for potential problems

   LAYER_LOOP: DO IL=1,MAX_LAYERS

      IF (ADIABATIC .AND. MATL_ID(IL,1)/='null') THEN
         WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' is ADIABATIC and cannot have a MATL_ID'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   
      IF (THICKNESS(IL)>0._EB .AND. MATL_ID(IL,1)=='null') THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: SURF '//TRIM(SF%ID)//' must have a MATL_ID for Layer ',IL
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   
      IF (THICKNESS(IL)<0._EB .AND. MATL_ID(IL,1)/='null') THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: SURF '//TRIM(SF%ID)// ' must have a specified THICKNESS for Layer ',IL
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   
   ENDDO LAYER_LOOP
   
   ! Identify the default SURF 

   IF (ID==SURF_DEFAULT) DEFAULT_SURF_INDEX = N

   ! Pack SURF parameters into the SURFACE derived type
 
   SF                      => SURFACE(N)
   SF%ADIABATIC            = ADIABATIC
   SELECT CASE(BACKING)
      CASE('VOID')
         SF%BACKING        = VOID
      CASE('INSULATED')
         SF%BACKING        = INSULATED
      CASE('EXPOSED')
         SF%BACKING        = EXPOSED
      CASE DEFAULT
         WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//', BACKING '//TRIM(BACKING)//' not recognized'
         CALL SHUTDOWN(MESSAGE)
   END SELECT
   SF%BURN_AWAY            = BURN_AWAY
   SF%CELL_SIZE_FACTOR     = CELL_SIZE_FACTOR
   SF%CONVECTIVE_HEAT_FLUX = 1000._EB*CONVECTIVE_HEAT_FLUX
   SF%NET_HEAT_FLUX        = 1000._EB*NET_HEAT_FLUX
   SF%DUCT_PATH            = DUCT_PATH
   SF%E_COEFFICIENT        = E_COEFFICIENT
   SF%EMISSIVITY           = EMISSIVITY   
   SF%FREE_SLIP            = FREE_SLIP
   SF%NO_SLIP              = NO_SLIP
   SF%FYI                  = FYI    
   SF%EXTERNAL_FLUX        = 1000._EB*EXTERNAL_FLUX
   SELECT CASE(GEOMETRY)
      CASE('CARTESIAN')
         SF%GEOMETRY       = SURF_CARTESIAN      
      CASE('CYLINDRICAL')
         SF%GEOMETRY       = SURF_CYLINDRICAL
         SF%BACKING        = INSULATED
      CASE('SPHERICAL')
         SF%GEOMETRY       = SURF_SPHERICAL
         SF%BACKING        = INSULATED
   END SELECT
   SF%H_V                  = 1000._EB*HEAT_OF_VAPORIZATION
   SF%HRRPUA               = 1000._EB*HRRPUA
   SF%MLRPUA               = MLRPUA
   SF%LAYER_DIVIDE         = LAYER_DIVIDE
   SF%LEAK_PATH            = LEAK_PATH
   SF%MASS_FLUX(:)         = MASS_FLUX(:)
   SF%MASS_FRACTION(:)     = MASS_FRACTION(:)
   SF%MAX_PRESSURE         = MAX_PRESSURE
   SF%REGRID_FACTOR        = REGRID_FACTOR
   SF%REGRID_FACTOR        = MIN(SF%REGRID_FACTOR,1.0_EB)
   SF%SHRINK               = .FALSE.
   SF%NPPC                 = NPPC
   SF%PARTICLE_MASS_FLUX   = PARTICLE_MASS_FLUX
   SF%PART_ID              = PART_ID
   SF%PLE                  = PLE
   SF%POROUS               = POROUS
   SELECT CASE (PROFILE)
      CASE('null')
         SF%PROFILE        = 0
      CASE('ATMOSPHERIC')
         SF%PROFILE        = ATMOSPHERIC
      CASE('PARABOLIC')
         SF%PROFILE        = PARABOLIC
      CASE('1D-PARABOLIC')
         SF%PROFILE        = ONED_PARABOLIC
   END SELECT
   SF%RAMP_MF              = RAMP_MF
   SF%RAMP_Q               = RAMP_Q 
   SF%RAMP_V               = RAMP_V  
   SF%RAMP_T               = RAMP_T  
   IF (COLOR/='null') THEN
      IF (COLOR=='INVISIBLE') THEN
         TRANSPARENCY = 0._EB
      ELSE
         CALL COLOR2RGB(RGB,COLOR)
      ENDIF
   ENDIF
   IF (ANY(RGB< 0)) THEN
      RGB(1) = 255
      RGB(2) = 204
      RGB(3) = 102
   ENDIF
   SF%RGB                  = RGB
   SF%ROUGHNESS            = ROUGHNESS
   SF%TRANSPARENCY         = TRANSPARENCY
   SF%STRETCH_FACTOR       = STRETCH_FACTOR
   SF%STRETCH_FACTOR       = MAX(1.0_EB,SF%STRETCH_FACTOR)
   SF%SURFACE_DENSITY      = SURFACE_DENSITY
   SF%TAU(0:)              = TAU_MF(0:)
   SF%TAU(TIME_HEAT)       = TAU_Q
   SF%TAU(TIME_VELO)       = TAU_V
   SF%TAU(TIME_TEMP)       = TAU_T
   SF%TEXTURE_MAP          = TEXTURE_MAP
   SF%TEXTURE_WIDTH        = TEXTURE_WIDTH
   SF%TEXTURE_HEIGHT       = TEXTURE_HEIGHT
   SF%TMP_IGN              = IGNITION_TEMPERATURE + TMPM
   SF%VEL                  = VEL
   SF%VEL_T                = VEL_T
   SF%VOLUME_FLUX          = VOLUME_FLUX
   SF%Z0                   = Z0
   SF%MASS_FLUX_TOTAL      = MASS_FLUX_TOTAL
   
   ! Set various logical parameters

   IF (SF%VEL_T(1)/=0._EB .OR. SF%VEL_T(2)/=0._EB) SF%SPECIFIED_TANGENTIAL_VELOCITY = .TRUE.

   IF (SF%HRRPUA>0._EB .OR. SF%MLRPUA>0._EB) MIXTURE_FRACTION=.TRUE.
   
   ! Count the number of layers for the surface, and compile a LIST of all material names and indices

   SF%N_LAYERS = 0
   N_LIST = 0
   NAME_LIST = 'null'
   SF%THICKNESS  = 0._EB
   SF%LAYER_MATL_INDEX = 0
   SF%LAYER_DENSITY    = 0._EB
   INDEX_LIST = -1
   COUNT_LAYERS: DO NL=1,MAX_LAYERS
      IF (THICKNESS(NL) < 0._EB) EXIT COUNT_LAYERS
      SF%N_LAYERS = SF%N_LAYERS + 1
      SF%LAYER_THICKNESS(NL) = THICKNESS(NL)
      SF%N_LAYER_MATL(NL) = 0
      IF (NL==1) SF%EMISSIVITY = 0._EB
      COUNT_LAYER_MATL: DO NN=1,MAX_MATERIALS
         IF (MATL_ID(NL,NN) == 'null') CYCLE COUNT_LAYER_MATL
         N_LIST = N_LIST + 1
         NAME_LIST(N_LIST) = MATL_ID(NL,NN)
         SF%N_LAYER_MATL(NL) = SF%N_LAYER_MATL(NL) + 1
         SF%LAYER_MATL_NAME(NL,NN) = MATL_ID(NL,NN)
         SF%LAYER_MATL_FRAC(NL,NN) = MATL_MASS_FRACTION(NL,NN)
         DO NNN=1,N_MATL
            IF (MATL_NAME(NNN)==NAME_LIST(N_LIST)) THEN
               INDEX_LIST(N_LIST) = NNN
               SF%LAYER_MATL_INDEX(NL,NN) = NNN
               SF%LAYER_DENSITY(NL) = SF%LAYER_DENSITY(NL)+SF%LAYER_MATL_FRAC(NL,NN)/MATERIAL(NNN)%RHO_S
               IF (NL==1) SF%EMISSIVITY = SF%EMISSIVITY + &
                  MATERIAL(NNN)%EMISSIVITY*SF%LAYER_MATL_FRAC(NL,NN)/MATERIAL(NNN)%RHO_S ! volume based
            ENDIF
         ENDDO
         IF (INDEX_LIST(N_LIST)<0) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'MATL: ',TRIM(NAME_LIST(N_LIST)),', on SURF: ',TRIM(SF%ID),', does not exist'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO COUNT_LAYER_MATL
      IF (SF%LAYER_DENSITY(NL) > 0._EB) SF%LAYER_DENSITY(NL) = 1./SF%LAYER_DENSITY(NL)
      IF (NL==1) SF%EMISSIVITY = SF%EMISSIVITY*SF%LAYER_DENSITY(NL)
      SF%THICKNESS = SF%THICKNESS + SF%LAYER_THICKNESS(NL)
   ENDDO COUNT_LAYERS

   ! Define mass flux division point

   IF (SF%LAYER_DIVIDE < 0._EB) THEN
      IF (SF%BACKING==EXPOSED) THEN 
         SF%LAYER_DIVIDE = 0.5_EB * REAL(SF%N_LAYERS,EB)
      ELSE
         SF%LAYER_DIVIDE = REAL(SF%N_LAYERS+1)
      ENDIF
   ENDIF

   ! Add residue materials

   DO I = 1,MAX_STEPS    ! repeat the residue loop to find chained reactions - allows MAX_STEPS steps
      N_LIST2 = N_LIST
      DO NN = 1,N_LIST2
         ML=>MATERIAL(INDEX_LIST(NN))
         ADD_REAC_MATL: DO NNN=1,ML%N_REACTIONS
            IF (ML%RESIDUE_MATL_NAME(NNN) == 'null') CYCLE ADD_REAC_MATL
            IF (ANY(NAME_LIST==ML%RESIDUE_MATL_NAME(NNN))) CYCLE ADD_REAC_MATL
            N_LIST = N_LIST + 1
            IF (N_LIST.GT.MAX_MATERIALS_TOTAL) CALL SHUTDOWN('ERROR: Too many materials in the surface.')
            NAME_LIST (N_LIST) = ML%RESIDUE_MATL_NAME(NNN)
            INDEX_LIST(N_LIST) = ML%RESIDUE_MATL_INDEX(NNN)
         ENDDO ADD_REAC_MATL
      ENDDO
   ENDDO

   ! Eliminate multiply counted materials from the list

   N_LIST2 = N_LIST
   WEED_MATL_LIST: DO NN=1,N_LIST
      DO NNN=1,NN-1
         IF (NAME_LIST(NNN)==NAME_LIST(NN)) THEN
            NAME_LIST(NN)  = 'null'
            INDEX_LIST(NN) = 0 
            N_LIST2 = N_LIST2-1
            CYCLE WEED_MATL_LIST
         ENDIF
      ENDDO
   ENDDO WEED_MATL_LIST

   ! Allocate parameters indexed by layer

   SF%N_MATL     = N_LIST2
   SF%THERMALLY_THICK = .FALSE.
   IF (SF%N_LAYERS > 0) THEN
      SF%THERMALLY_THICK = .TRUE.
      SF%TMP_INNER                    = TMP_INNER + TMPM
      SF%TMP_FRONT                    = TMP_INNER + TMPM
      SF%TMP_BACK                     = TMP_BACK  + TMPM
      ALLOCATE(SF%N_LAYER_CELLS(SF%N_LAYERS))            ! The number of cells in each layer
      ALLOCATE(SF%MIN_DIFFUSIVITY(SF%N_LAYERS))          ! The smallest diffusivity of materials in each layer
      ALLOCATE(SF%MATL_NAME(SF%N_MATL))                  ! The list of all material names associated with the surface
      ALLOCATE(SF%MATL_INDEX(SF%N_MATL))                 ! The list of all material indices associated with the surface
      ALLOCATE(SF%RESIDUE_INDEX(SF%N_MATL,MAX_REACTIONS))! Each material associated with the surface has a RESIDUE
   ELSE
      SF%TMP_FRONT                  = TMP_FRONT + TMPM
      SF%TMP_INNER                  = SF%TMP_FRONT
      SF%TMP_BACK                   = SF%TMP_FRONT
   ENDIF
   TMPMIN        = MIN(TMPMIN,SF%TMP_FRONT,SF%TMP_INNER,SF%TMP_BACK)

   ! Store the names and indices of all materials associated with the surface

   NNN = 0
   DO NN=1,N_LIST
      IF (NAME_LIST(NN)/='null') THEN
         NNN = NNN + 1
         SF%MATL_NAME(NNN)  = NAME_LIST(NN)
         SF%MATL_INDEX(NNN) = INDEX_LIST(NN)
      ENDIF
   ENDDO

   ! Store the RESIDUE indices and detect (possibly) shrinking surfaces

   DO NN=1,SF%N_MATL
      ML => MATERIAL(SF%MATL_INDEX(NN))
      IF (ML%N_REACTIONS>0 .AND. SF%TMP_IGN<5000._EB) THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: SURF '//TRIM(SF%ID)// ' cannot have a REACting MATL and IGNITION_TEMPERATURE'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      DO J=1,ML%N_REACTIONS
         DO NNN=1,SF%N_MATL
            IF (ML%RESIDUE_MATL_INDEX(J)==SF%MATL_INDEX(NNN)) SF%RESIDUE_INDEX(NN,J) = NNN
         ENDDO
         IF (ML%NU_RESIDUE(J).EQ.0._EB) SF%SHRINK = .TRUE.
      ENDDO
      IF (ML%PYROLYSIS_MODEL==PYROLYSIS_LIQUID) SF%SHRINK = .TRUE.
   ENDDO
   IF (.NOT. SHRINK ) SF%SHRINK = SHRINK

   ! Thermal boundary conditions
   IF (SF%ADIABATIC .AND. (SF%NET_HEAT_FLUX < 1.E12_EB .OR. SF%CONVECTIVE_HEAT_FLUX/=0._EB)) THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: SURF '//TRIM(SF%ID)//&
                                 ' cannot have both ADIABATIC and NET_HEAT_FLUX or CONVECTIVE_HEAT_FLUX'
         CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (SF%NET_HEAT_FLUX < 1.E12_EB .AND. SF%CONVECTIVE_HEAT_FLUX/=0._EB) THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: SURF '//TRIM(SF%ID)// ' cannot have both NET_HEAT_FLUX or CONVECTIVE_HEAT_FLUX'
         CALL SHUTDOWN(MESSAGE)
   ENDIF
   
   SF%THERMAL_BC_INDEX = SPECIFIED_TEMPERATURE
   IF (SF%ADIABATIC) THEN
                                       SF%THERMAL_BC_INDEX = NET_FLUX_BC
                                       SF%NET_HEAT_FLUX = 0._EB
   ENDIF
   IF (SF%NET_HEAT_FLUX < 1.E12_EB)    SF%THERMAL_BC_INDEX = NET_FLUX_BC
   IF (SF%CONVECTIVE_HEAT_FLUX/=0._EB) SF%THERMAL_BC_INDEX = CONVECTIVE_FLUX_BC
   IF (SF%THERMALLY_THICK)             SF%THERMAL_BC_INDEX = THERMALLY_THICK
   IF (SF%PROFILE==ATMOSPHERIC)        SF%THERMAL_BC_INDEX = INFLOW_OUTFLOW

   ! Ramps

   IF (SF%RAMP_Q/='null') THEN
      CALL GET_RAMP_INDEX(SF%RAMP_Q,'TIME',NR)
      SF%RAMP_INDEX(TIME_HEAT) = NR
   ELSE
      IF (SF%TAU(TIME_HEAT) > 0._EB) SF%RAMP_INDEX(TIME_HEAT) = TANH_RAMP
      IF (SF%TAU(TIME_HEAT) < 0._EB) SF%RAMP_INDEX(TIME_HEAT) = TSQR_RAMP
   ENDIF

   IF (SF%RAMP_V/='null') THEN
      CALL GET_RAMP_INDEX(SF%RAMP_V,'TIME',NR)
      SF%RAMP_INDEX(TIME_VELO) = NR
   ELSE
      IF (SF%TAU(TIME_VELO) > 0._EB) SF%RAMP_INDEX(TIME_VELO) = TANH_RAMP
      IF (SF%TAU(TIME_VELO) < 0._EB) SF%RAMP_INDEX(TIME_VELO) = TSQR_RAMP
   ENDIF

   IF (SF%RAMP_T/='null') THEN
      CALL GET_RAMP_INDEX(SF%RAMP_T,'TIME',NR)
      SF%RAMP_INDEX(TIME_TEMP) = NR
   ELSE
      IF (SF%TAU(TIME_TEMP) > 0._EB) SF%RAMP_INDEX(TIME_TEMP) = TANH_RAMP
      IF (SF%TAU(TIME_TEMP) < 0._EB) SF%RAMP_INDEX(TIME_TEMP) = TSQR_RAMP
   ENDIF

ENDDO READ_SURF_LOOP
 
 
CONTAINS
 
SUBROUTINE SET_SURF_DEFAULTS
 
ADIABATIC               = .FALSE.
BACKING                 = 'VOID'
BURN_AWAY               = .FALSE.
CELL_SIZE_FACTOR        = 1.0
COLOR                   = 'null'
CONVECTIVE_HEAT_FLUX    = 0._EB
NET_HEAT_FLUX           = 1.E12_EB
DUCT_PATH               = 0 
E_COEFFICIENT           = 0._EB
EMISSIVITY              = 0.9_EB
EXTERNAL_FLUX           = 0._EB
FREE_SLIP               = .FALSE.
NO_SLIP                 = .FALSE.
FYI                     = 'null'
GEOMETRY                = 'CARTESIAN'
HEAT_OF_VAPORIZATION    = 0._EB
HRRPUA                  = 0._EB
ID                      = 'null'
IGNITION_TEMPERATURE    = 5000._EB
LAYER_DIVIDE            = -1._EB
LEAK_PATH               = -1 
MASS_FLUX               = 0._EB
MASS_FLUX_TOTAL         = 0._EB
MASS_FRACTION           = -1._EB
MATL_ID                 = 'null'
MATL_MASS_FRACTION      = 0._EB
MATL_MASS_FRACTION(:,1) = 1._EB
MAX_PRESSURE            = 1.E12_EB
MLRPUA                  = 0._EB
NPPC                    = 1
PARTICLE_MASS_FLUX      = 0._EB
PART_ID                 = 'null'
PLE                     = 0.3_EB
POROUS                  = .FALSE.
PROFILE                 = 'null'
RAMP_MF                 = 'null'
RAMP_Q                  = 'null'
RAMP_V                  = 'null'
RAMP_T                  = 'null'
REGRID_FACTOR           = 0.9
RGB                     = -1
SHRINK                  = .TRUE.
IF (LES) ROUGHNESS      = 0._EB !4.5E-5_EB  ! meters, commercial steel
IF (DNS) ROUGHNESS      = 0._EB
STRETCH_FACTOR          = 2._EB
SURFACE_DENSITY         = -1._EB
TAU_MF                  =  1._EB
TAU_Q                   = 1._EB
TAU_V                   = 1._EB
TAU_T                   = 1._EB
TEXTURE_MAP             = 'null'
TEXTURE_WIDTH           = 1._EB
TEXTURE_HEIGHT          = 1._EB
THICKNESS               = -1._EB
TMP_BACK                = TMPA-TMPM
TMP_FRONT               = TMPA-TMPM 
TMP_INNER               = TMPA-TMPM
TRANSPARENCY            = 1._EB
VEL_T                   = 0._EB
VEL                     = 0._EB
VOLUME_FLUX             = 0._EB
Z0                      = 10._EB
 
END SUBROUTINE SET_SURF_DEFAULTS
 
END SUBROUTINE READ_SURF


 
SUBROUTINE PROC_SURF_1

! Go through the SURF types and process

USE MATH_FUNCTIONS, ONLY : GET_RAMP_INDEX
INTEGER :: N,NSPC,NR,IPC
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC
 
PROCESS_SURF_LOOP: DO N=0,N_SURF
 
   SF => SURFACE(N)
 
   ! Get ramps for the surface mass fraction and flux

   DO NSPC=0,N_SPECIES
      IF (SF%RAMP_MF(NSPC)/='null') THEN
         CALL GET_RAMP_INDEX(SF%RAMP_MF(NSPC),'TIME',NR)
         SF%RAMP_INDEX(NSPC) = NR
      ELSE
         IF (SF%TAU(NSPC) > 0._EB) SF%RAMP_INDEX(NSPC) = TANH_RAMP 
         IF (SF%TAU(NSPC) < 0._EB) SF%RAMP_INDEX(NSPC) = TSQR_RAMP 
      ENDIF
   ENDDO 
   
   DO IPC=1,N_PART
      PC=>PARTICLE_CLASS(IPC)
      IF (PC%SURF_ID==SF%ID) PC%SURF_INDEX = N
   ENDDO

ENDDO PROCESS_SURF_LOOP   
 
END SUBROUTINE PROC_SURF_1



SUBROUTINE PROC_SURF_2

! Go through the SURF types and process
 
INTEGER :: IPC,N,NN,NNN,NL,NSPC
REAL(EB) :: ADJUSTED_LAYER_DENSITY
LOGICAL :: BURNING,BLOWING,SUCKING
TYPE(PARTICLE_CLASS_TYPE), POINTER :: PC
 
PROCESS_SURF_LOOP: DO N=0,N_SURF
 
   SF => SURFACE(N)
   IF (SF%THERMALLY_THICK) ML => MATERIAL(SF%LAYER_MATL_INDEX(1,1))
   
   ! Particle Information
 
   SF%PART_INDEX = 0
   IF (SF%PART_ID/='null') THEN
      DO IPC=1,N_PART
         PC=>PARTICLE_CLASS(IPC)
         IF (PC%ID==SF%PART_ID)  SF%PART_INDEX = IPC
      ENDDO
      IF (SF%PART_INDEX==0) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: PART_ID '//TRIM(SF%PART_ID)//' not found'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      DROPLET_FILE=.TRUE.
   ENDIF

  ! Determine if the surface is combustible/burning
 
   SF%PYROLYSIS_MODEL = PYROLYSIS_NONE
   BURNING  = .FALSE.
   BLOWING  = .FALSE.
   SUCKING  = .FALSE.
   IF (SF%N_LAYERS > 0) THEN
      DO NN=1,SF%N_MATL
         ML => MATERIAL(SF%MATL_INDEX(NN))
         IF (ML%PYROLYSIS_MODEL/=PYROLYSIS_NONE) THEN
            SF%PYROLYSIS_MODEL = PYROLYSIS_MATERIAL
            IF (ANY(ML%NU_FUEL>0._EB))  THEN
               BURNING = .TRUE.
               SF%TAU(TIME_HEAT) = 0._EB
            ENDIF
         ENDIF
      ENDDO   
   ENDIF
   IF (SF%HRRPUA>0._EB .OR. SF%MLRPUA>0._EB) THEN
      BURNING = .TRUE.
      SF%PYROLYSIS_MODEL = PYROLYSIS_SPECIFIED
   ENDIF

   IF (SF%VEL<0._EB .OR. SF%VOLUME_FLUX<0._EB) BLOWING = .TRUE.
   IF (SF%VEL>0._EB .OR. SF%VOLUME_FLUX>0._EB) SUCKING = .TRUE.
   IF (BLOWING .OR. SUCKING) SF%SPECIFIED_NORMAL_VELOCITY = .TRUE.

   IF (BURNING .AND. (BLOWING .OR. SUCKING)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' cannot have a specified velocity or volume flux'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
 
   ! set predefined HRRPUA

   BURNING_IF: IF (BURNING) THEN
      IF (SF%HRRPUA>0._EB) THEN
         IF (CO_PRODUCTION) THEN
            RN => REACTION(2)         
         ELSE
            RN => REACTION(1)
         ENDIF
         SF%MASS_FLUX(I_FUEL) = SF%HRRPUA/ (RN%HEAT_OF_COMBUSTION*RN%Y_F_INLET)
      ENDIF
      IF (SF%MLRPUA>0._EB) SF%MASS_FLUX(I_FUEL) = SF%MLRPUA
      SF%TAU(I_FUEL)        = SF%TAU(TIME_HEAT)
      SF%RAMP_MF(I_FUEL)    = SF%RAMP_Q
      SF%RAMP_INDEX(I_FUEL) = SF%RAMP_INDEX(TIME_HEAT) 
   ENDIF BURNING_IF

   ! Compute surface density

   IF (SF%SURFACE_DENSITY < 0._EB) THEN
      SF%SURFACE_DENSITY = 0._EB
      DO NL=1,SF%N_LAYERS
         ADJUSTED_LAYER_DENSITY = 0._EB
         MATL_LOOP:DO NN=1,SF%N_LAYER_MATL(NL)
            NNN = SF%LAYER_MATL_INDEX(NL,NN)
            ML => MATERIAL(NNN)
            ADJUSTED_LAYER_DENSITY = ADJUSTED_LAYER_DENSITY + SF%LAYER_MATL_FRAC(NL,NN)/ML%RHO_S
         ENDDO MATL_LOOP
         IF (ADJUSTED_LAYER_DENSITY > 0._EB) ADJUSTED_LAYER_DENSITY = 1./ADJUSTED_LAYER_DENSITY
         SF%SURFACE_DENSITY = SF%SURFACE_DENSITY + SF%LAYER_THICKNESS(NL)*ADJUSTED_LAYER_DENSITY
      ENDDO
   ENDIF

   IF ((SF%SURFACE_DENSITY == 0._EB) .AND. SF%BURN_AWAY) THEN
      WRITE(MESSAGE,'(A)') 'WARNING: SURF '//TRIM(SF%ID)//' has BURN_AWAY set but zero combustible density'
      IF (MYID==0) WRITE(LU_ERR,'(A)') TRIM(MESSAGE)
   ENDIF   

   ! Ignition Time

   SF%T_IGN = T_BEGIN 
   IF (SF%TMP_IGN<5000._EB)                    SF%T_IGN = T_END
   IF (SF%PYROLYSIS_MODEL==PYROLYSIS_MATERIAL) SF%T_IGN = T_END

   ! Species Arrays and Method of Mass Transfer (SPECIES_BC_INDEX)
 
   SF%SPECIES_BC_INDEX = NO_MASS_FLUX

   IF (ANY(SF%MASS_FRACTION>=0._EB) .AND. (ANY(SF%MASS_FLUX/=0._EB).OR.BURNING)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' cannot specify mass fraction with mass flux and/or burning'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
   IF (ANY(SF%MASS_FRACTION>=0._EB) .AND. SUCKING) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' cannot specify both mass fraction and outflow velocity'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
   IF (ANY(SF%LEAK_PATH>=0) .AND. (BLOWING .OR. SUCKING)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' cannot leak and blow at the same time'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
   IF (ANY(SF%MASS_FLUX/=0._EB) .AND. (BLOWING .OR. SUCKING)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' cannot have both a mass flux and specified velocity'
      CALL SHUTDOWN(MESSAGE)
      ENDIF

   IF (BLOWING .OR. SUCKING)                      SF%SPECIES_BC_INDEX = SPECIFIED_MASS_FRACTION
   IF (ANY(SF%MASS_FRACTION>=0._EB))              SF%SPECIES_BC_INDEX = SPECIFIED_MASS_FRACTION
   IF (ANY(SF%MASS_FLUX    /=0._EB) .OR. &
       SF%PYROLYSIS_MODEL==PYROLYSIS_MATERIAL) SF%SPECIES_BC_INDEX = SPECIFIED_MASS_FLUX

   IF (SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FRACTION) THEN
      SPECIES_LOOP: DO NSPC=0,N_SPECIES
         IF (SF%MASS_FRACTION(NSPC)<0._EB) SF%MASS_FRACTION(NSPC) = 0._EB
      ENDDO SPECIES_LOOP
   ENDIF
 
   ! Texture map info

   SF%SURF_TYPE = 0
   IF (SF%TEXTURE_MAP/='null') SF%SURF_TYPE = 1
  
   ! Set BCs for various boundary types
 
   SF%VELOCITY_BC_INDEX = WALL_MODEL 
   IF (DNS) SF%VELOCITY_BC_INDEX = NO_SLIP_BC
   IF (SF%FREE_SLIP) SF%VELOCITY_BC_INDEX = FREE_SLIP_BC
   IF (SF%NO_SLIP)   SF%VELOCITY_BC_INDEX = NO_SLIP_BC

   IF (N==OPEN_SURF_INDEX) THEN
      SF%THERMAL_BC_INDEX = INFLOW_OUTFLOW
      SF%SPECIES_BC_INDEX = NO_MASS_FLUX
      SF%VELOCITY_BC_INDEX = FREE_SLIP_BC
      SF%SURF_TYPE = 2
   ENDIF
   IF (N==MIRROR_SURF_INDEX) THEN
      SF%THERMAL_BC_INDEX = NO_CONVECTION
      SF%SPECIES_BC_INDEX = NO_MASS_FLUX
      SF%VELOCITY_BC_INDEX = FREE_SLIP_BC
      SF%SURF_TYPE = -2
      SF%EMISSIVITY = 0._EB
   ENDIF
   IF (N==INTERPOLATED_SURF_INDEX) THEN
      SF%THERMAL_BC_INDEX = INTERPOLATED_BC
      SF%SPECIES_BC_INDEX = INTERPOLATED_BC
      SF%VELOCITY_BC_INDEX = INTERPOLATED_VELOCITY_BC
   ENDIF
   IF (N==PERIODIC_SURF_INDEX) THEN
      SF%THERMAL_BC_INDEX = INTERPOLATED_BC
      SF%SPECIES_BC_INDEX = INTERPOLATED_BC
      SF%VELOCITY_BC_INDEX = INTERPOLATED_VELOCITY_BC
   ENDIF
 
ENDDO PROCESS_SURF_LOOP
 
END SUBROUTINE PROC_SURF_2



SUBROUTINE PROC_WALL

! Set up 1-D grids and arrays for thermally-thick calcs

USE GEOMETRY_FUNCTIONS
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP

INTEGER :: IBC,N,NL,NWP_MAX
REAL(EB) :: K_S_0,C_S_0,SMALLEST_CELL_SIZE(MAX_LAYERS)

! Calculate ambient temperature thermal DIFFUSIVITY for each MATERIAL, to be used in determining number of solid cells

DO N=1,N_MATL
   ML => MATERIAL(N)
   IF (ML%K_S>0._EB) THEN
      K_S_0 = ML%K_S
   ELSE
      K_S_0 = EVALUATE_RAMP(TMPA,0._EB,-NINT(ML%K_S))
   ENDIF
   IF (ML%C_S>0._EB) THEN
      C_S_0 = ML%C_S
   ELSE
      C_S_0 = EVALUATE_RAMP(TMPA,0._EB,-NINT(ML%C_S))*1000._EB
   ENDIF
   ML%DIFFUSIVITY = K_S_0/(C_S_0*ML%RHO_S)
ENDDO

NWP_MAX = 0  ! For some utility arrays, need to know the greatest number of points of all surface types
 
! Loop through all surfaces, looking for those that are thermally-thick (have layers).
! Compute smallest cell size for each layer such that internal cells double in size.
! Each layer should have an odd number of cells.

SURF_GRID_LOOP: DO IBC=0,N_SURF

   SF => SURFACE(IBC)
   IF (SF%THERMAL_BC_INDEX /= THERMALLY_THICK) CYCLE SURF_GRID_LOOP

   ! Compute number of points per layer, and then sum up to get total points for the surface

   SF%N_CELLS = 0
   DO NL=1,SF%N_LAYERS
      SF%MIN_DIFFUSIVITY(NL) = 1000000._EB
      DO N = 1,SF%N_LAYER_MATL(NL) 
         ML => MATERIAL(SF%LAYER_MATL_INDEX(NL,N))
         SF%MIN_DIFFUSIVITY(NL) = MIN(SF%MIN_DIFFUSIVITY(NL),ML%DIFFUSIVITY)
      ENDDO
      CALL GET_N_LAYER_CELLS(SF%MIN_DIFFUSIVITY(NL),SF%LAYER_THICKNESS(NL),SF%STRETCH_FACTOR, &
                             SF%CELL_SIZE_FACTOR,SF%N_LAYER_CELLS(NL),SMALLEST_CELL_SIZE(NL))
      SF%N_CELLS = SF%N_CELLS + SF%N_LAYER_CELLS(NL)
   ENDDO

! Allocate arrays to hold x_s, 1/dx_s (center to center, RDXN), 1/dx_s (edge to edge, RDX)

   NWP_MAX = MAX(NWP_MAX,SF%N_CELLS)
   ALLOCATE(SF%DX(1:SF%N_CELLS))
   ALLOCATE(SF%RDX(0:SF%N_CELLS+1))
   ALLOCATE(SF%RDXN(0:SF%N_CELLS))
   ALLOCATE(SF%DX_WGT(0:SF%N_CELLS))
   ALLOCATE(SF%X_S(0:SF%N_CELLS))
   ALLOCATE(SF%LAYER_INDEX(0:SF%N_CELLS+1))
   ALLOCATE(SF%MF_FRAC(1:SF%N_CELLS))

! Compute node coordinates 

   CALL GET_WALL_NODE_COORDINATES(SF%N_CELLS,SF%N_LAYERS,SF%N_LAYER_CELLS, &
         SMALLEST_CELL_SIZE(1:SF%N_LAYERS),SF%STRETCH_FACTOR,SF%X_S)

   CALL GET_WALL_NODE_WEIGHTS(SF%N_CELLS,SF%N_LAYERS,SF%N_LAYER_CELLS,SF%LAYER_THICKNESS,SF%GEOMETRY, &
         SF%X_S,SF%LAYER_DIVIDE,SF%DX,SF%RDX,SF%RDXN,SF%DX_WGT,SF%DXF,SF%DXB,SF%LAYER_INDEX,SF%MF_FRAC)

! Determine if surface has internal radiation

   SF%INTERNAL_RADIATION = .FALSE.
   DO NL=1,SF%N_LAYERS
   DO N =1,SF%N_LAYER_MATL(NL)
      ML => MATERIAL(SF%LAYER_MATL_INDEX(NL,N))
      IF (ML%KAPPA_S<5.0E4_EB) SF%INTERNAL_RADIATION = .TRUE.
   ENDDO
   ENDDO
 
ENDDO SURF_GRID_LOOP
 
ALLOCATE(AAS(NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','AAS',IZERO)
ALLOCATE(CCS(NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','CCS',IZERO)
ALLOCATE(BBS(NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','BBS',IZERO)
ALLOCATE(DDS(NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','DDS',IZERO)
ALLOCATE(DDT(NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','DDT',IZERO)
ALLOCATE(K_S(0:NWP_MAX+1),STAT=IZERO)
CALL ChkMemErr('INIT','K_S',IZERO)
ALLOCATE(C_S(0:NWP_MAX+1),STAT=IZERO)
CALL ChkMemErr('INIT','C_S',IZERO)
ALLOCATE(Q_S(1:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','Q_S',IZERO)
ALLOCATE(RHO_S(0:NWP_MAX+1),STAT=IZERO)
CALL ChkMemErr('INIT','RHO_S',IZERO)
ALLOCATE(RHOCBAR(1:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','RHOCBAR',IZERO)
ALLOCATE(KAPPA_S(1:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','KAPPA_S',IZERO)
ALLOCATE(X_S_NEW(0:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','X_S_NEW',IZERO)
ALLOCATE(DX_S(1:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','DX_S',IZERO)
ALLOCATE(RDX_S(0:NWP_MAX+1),STAT=IZERO)
CALL ChkMemErr('INIT','RDX_S',IZERO)
ALLOCATE(RDXN_S(0:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','RDXN_S',IZERO)
ALLOCATE(R_S(0:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','R_S',IZERO)
ALLOCATE(DX_WGT_S(0:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','DX_WGT_S',IZERO)
ALLOCATE(LAYER_INDEX(0:NWP_MAX+1),STAT=IZERO)
CALL ChkMemErr('INIT','LAYER_INDEX',IZERO)
ALLOCATE(MF_FRAC(1:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','MF_FRAC',IZERO)

END SUBROUTINE PROC_WALL


SUBROUTINE READ_PRES

REAL(EB) :: PLACE_HOLDER
NAMELIST /PRES/ PLACE_HOLDER

REWIND(LU_INPUT)
READ_LOOP: DO
   CALL CHECKREAD('PRES',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_LOOP
   READ(LU_INPUT,PRES,END=23,ERR=24,IOSTAT=IOS)
   24 IF (IOS>0) THEN
      CALL SHUTDOWN(' ERROR: Problem with PRES line')
   ENDIF
ENDDO READ_LOOP
23 REWIND(LU_INPUT)

END SUBROUTINE READ_PRES


SUBROUTINE READ_RADI
USE RADCONS
NAMELIST /RADI/ TIME_STEP_INCREMENT,NUMBER_RADIATION_ANGLES,ANGLE_INCREMENT,KAPPA0, &
                WIDE_BAND_MODEL,CH4_BANDS,PATH_LENGTH,NMIEANG,RADTMP,RADIATIVE_FRACTION

IF (LES) RADIATIVE_FRACTION = 0.35_EB
IF (DNS) RADIATIVE_FRACTION = 0.00_EB
NUMBER_RADIATION_ANGLES = 100
TIME_STEP_INCREMENT     = 3
IF (TWO_D) THEN
   NUMBER_RADIATION_ANGLES = 50
   TIME_STEP_INCREMENT     = 2
ENDIF
 
KAPPA0          =   0._EB
RADTMP          = 900._EB
WIDE_BAND_MODEL = .FALSE.
CH4_BANDS       = .FALSE.
NMIEANG         = 15
PATH_LENGTH     = -1.0_EB ! calculate path based on the geometry
ANGLE_INCREMENT = -1
 
REWIND(LU_INPUT)
READ_LOOP: DO
   CALL CHECKREAD('RADI',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_LOOP
   READ(LU_INPUT,RADI,END=23,ERR=24,IOSTAT=IOS)
   24 IF (IOS>0) THEN
      CALL SHUTDOWN(' ERROR: Problem with RADI line')
   ENDIF
ENDDO READ_LOOP
23 REWIND(LU_INPUT)

RADTMP = RADTMP + TMPM

IF (WIDE_BAND_MODEL) THEN
   IF (CH4_BANDS) THEN
      NUMBER_SPECTRAL_BANDS = 9
   ELSE
      NUMBER_SPECTRAL_BANDS = 6
   ENDIF
   TIME_STEP_INCREMENT=MAX(1,TIME_STEP_INCREMENT)
   ANGLE_INCREMENT = 1
   UIIDIM=NUMBER_SPECTRAL_BANDS
ELSE
   NUMBER_SPECTRAL_BANDS = 1
   IF (ANGLE_INCREMENT < 0) ANGLE_INCREMENT = MAX(1,MIN(5,NUMBER_RADIATION_ANGLES/15))
   UIIDIM = ANGLE_INCREMENT
ENDIF

END SUBROUTINE READ_RADI


SUBROUTINE READ_CLIP

INTEGER  :: N 
REAL(EB) :: MINIMUM_DENSITY,MAXIMUM_DENSITY,MINIMUM_MASS_FRACTION(MAX_SPECIES),MAXIMUM_MASS_FRACTION(MAX_SPECIES), &
            MINIMUM_TEMPERATURE,MAXIMUM_TEMPERATURE
NAMELIST /CLIP/ MINIMUM_DENSITY,MAXIMUM_DENSITY,FYI,MINIMUM_MASS_FRACTION,MAXIMUM_MASS_FRACTION, &
                MINIMUM_TEMPERATURE,MAXIMUM_TEMPERATURE
 
! Check for user-defined mins and maxes.
 
MINIMUM_DENSITY       = -999._EB
MAXIMUM_DENSITY       = -999._EB
MINIMUM_TEMPERATURE   = -999._EB
MAXIMUM_TEMPERATURE   = -999._EB
MINIMUM_MASS_FRACTION = -999._EB
MAXIMUM_MASS_FRACTION = -999._EB
 
REWIND(LU_INPUT)
CLIP_LOOP: DO
   CALL CHECKREAD('CLIP',LU_INPUT,IOS) 
   IF (IOS==1) EXIT CLIP_LOOP
   READ(LU_INPUT,CLIP,END=431,ERR=432,IOSTAT=IOS)
   432 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with CLIP line')
ENDDO CLIP_LOOP
431 REWIND(LU_INPUT)
 
IF (MINIMUM_TEMPERATURE>-TMPM) TMPMIN = MINIMUM_TEMPERATURE + TMPM
IF (MAXIMUM_TEMPERATURE>-TMPM) TMPMAX = MAXIMUM_TEMPERATURE + TMPM

IF (TMPMAX > 5000._EB) CALL SHUTDOWN('MAXIMUM_TEMPERATURE cannot be greater than 4726.85 C (5000 K)')

IF (MINIMUM_DENSITY>0._EB) THEN
   RHOMIN = MINIMUM_DENSITY
ELSE
   RHOMIN = 0.1_EB*RHOA
ENDIF
IF (MAXIMUM_DENSITY>0._EB) THEN
   RHOMAX = MAXIMUM_DENSITY
ELSE
   RHOMAX = 3.0_EB*P_INF*MW_MAX/(R0*(TMPMIN+1._EB))  ! The 1 added to TMPMIN is to prevent a divide by zero error
ENDIF
DO N=1,N_SPECIES
   IF (MINIMUM_MASS_FRACTION(N)>-1._EB) YYMIN(N) = MINIMUM_MASS_FRACTION(N)
   IF (MAXIMUM_MASS_FRACTION(N)>-1._EB) YYMAX(N) = MAXIMUM_MASS_FRACTION(N)
ENDDO
 
END SUBROUTINE READ_CLIP
 
 
SUBROUTINE READ_RAMP
 
REAL(EB) :: X,T,F,TM
INTEGER  :: I,II,NN,N,NUMBER_INTERPOLATION_POINTS
CHARACTER(30) :: DEVC_ID,CTRL_ID
TYPE(RAMPS_TYPE), POINTER :: RP
NAMELIST /RAMP/ X,T,F,ID,FYI,NUMBER_INTERPOLATION_POINTS,DEVC_ID,CTRL_ID
 
IF (N_RAMP==0) RETURN
ALLOCATE(RAMPS(N_RAMP),STAT=IZERO)
CALL ChkMemErr('READ','RAMPS',IZERO)

! Count the number of points in each ramp
 
REWIND(LU_INPUT)
COUNT_RAMP_POINTS: DO N=1,N_RAMP
   RP => RAMPS(N)
   REWIND(LU_INPUT)
   RP%NUMBER_DATA_POINTS = 0
   SEARCH_LOOP: DO
      CALL CHECKREAD('RAMP',LU_INPUT,IOS)
      IF (IOS==1) EXIT SEARCH_LOOP
      READ(LU_INPUT,NML=RAMP,ERR=56,IOSTAT=IOS)
      IF (ID/=RAMP_ID(N)) CYCLE SEARCH_LOOP
      RP%NUMBER_DATA_POINTS = RP%NUMBER_DATA_POINTS + 1
      56 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with RAMP '//TRIM(RAMP_ID(N)) )
   ENDDO SEARCH_LOOP
   IF (RP%NUMBER_DATA_POINTS==0) THEN
      WRITE(MESSAGE,'(A,A,A)') 'ERROR: RAMP ',TRIM(RAMP_ID(N)), ' not found'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO COUNT_RAMP_POINTS

! Read the ramp functions
READ_RAMP_LOOP: DO N=1,N_RAMP
   RP => RAMPS(N)
   RP%DEVC_ID = 'null'
   RP%CTRL_ID = 'null'
   ALLOCATE(RP%INDEPENDENT_DATA(1:RP%NUMBER_DATA_POINTS))
   ALLOCATE(RP%DEPENDENT_DATA(1:RP%NUMBER_DATA_POINTS))
   REWIND(LU_INPUT)
   NN = 0
   NUMBER_INTERPOLATION_POINTS=5000
   SEARCH_LOOP2: DO 
      DEVC_ID = 'null'
      CTRL_ID = 'null'
      X       = -1.E6_EB
      CALL CHECKREAD('RAMP',LU_INPUT,IOS) 
      IF (IOS==1) EXIT SEARCH_LOOP2
      READ(LU_INPUT,RAMP)
      IF (ID/=RAMP_ID(N)) CYCLE SEARCH_LOOP2
      IF (DEVC_ID /='null') RP%DEVC_ID = DEVC_ID
      IF (CTRL_ID /='null') RP%CTRL_ID = CTRL_ID      
      IF (X>-1.E5_EB) THEN
         RAMP_TYPE(N) = 'X COORDINATE'
         SPATIAL_GRAVITY_VARIATION = .TRUE.
         STRATIFICATION = .FALSE.
         T = X
      ENDIF
      IF (RAMP_TYPE(N)=='TEMPERATURE') T = T + TMPM
      IF (RAMP_TYPE(N)=='TIME')        T = T_BEGIN + (T-T_BEGIN)/TIME_SHRINK_FACTOR
      NN = NN+1
      RP%INDEPENDENT_DATA(NN) = T
      IF (NN>1) THEN
         IF (T<=RP%INDEPENDENT_DATA(NN-1)) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: RAMP ',TRIM(RAMP_ID(N)), ' variable T must be monotonically increasing'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDIF
      RP%DEPENDENT_DATA(NN) = F
      RP%NUMBER_INTERPOLATION_POINTS = NUMBER_INTERPOLATION_POINTS
   ENDDO SEARCH_LOOP2
   RP%T_MIN = MINVAL(RP%INDEPENDENT_DATA)
   RP%T_MAX = MAXVAL(RP%INDEPENDENT_DATA)
   RP%SPAN = RP%T_MAX - RP%T_MIN
ENDDO READ_RAMP_LOOP
 
! Set up interpolated ramp values in INTERPOLATED_DATA and get control or device index
 
DO N=1,N_RAMP
   RP => RAMPS(N)
   RP%DT = RP%SPAN/REAL(RP%NUMBER_INTERPOLATION_POINTS,EB)   
   ALLOCATE(RAMPS(N)%INTERPOLATED_DATA(0:RP%NUMBER_INTERPOLATION_POINTS))
   RAMPS(N)%INTERPOLATED_DATA(0) = RP%DEPENDENT_DATA(1)
   DO I=1,RP%NUMBER_INTERPOLATION_POINTS-1
      TM = RP%INDEPENDENT_DATA(1) + REAL(I,EB)*RP%DT
      TLOOP: DO II=1,RP%NUMBER_DATA_POINTS-1
         IF (TM>=RP%INDEPENDENT_DATA(II) .AND. TM<RP%INDEPENDENT_DATA(II+1)) THEN
            RP%INTERPOLATED_DATA(I) = RP%DEPENDENT_DATA(II) +  (TM-RP%INDEPENDENT_DATA(II)) * &
                          (RP%DEPENDENT_DATA(II+1)-RP%DEPENDENT_DATA(II))/(RP%INDEPENDENT_DATA(II+1)-RP%INDEPENDENT_DATA(II))
            EXIT TLOOP
         ENDIF
      ENDDO TLOOP
   ENDDO
   RP%INTERPOLATED_DATA(RP%NUMBER_INTERPOLATION_POINTS) = RP%DEPENDENT_DATA(RP%NUMBER_DATA_POINTS)

   !Get Device or Control Index

   CALL SEARCH_CONTROLLER('RAMP',CTRL_ID,DEVC_ID,RP%DEVC_INDEX,RP%CTRL_INDEX,N)
       
ENDDO

END SUBROUTINE READ_RAMP
 
 
SUBROUTINE READ_TABL
 
REAL(EB) :: TABLE_DATA(6)
INTEGER  :: NN,N
TYPE(TABLES_TYPE), POINTER :: TA
NAMELIST /TABL/ ID,FYI,TABLE_DATA
 
IF (N_TABLE==0) RETURN

ALLOCATE(TABLES(N_TABLE),STAT=IZERO)
CALL ChkMemErr('READ','TABLES',IZERO)

! Count the number of points in each table
 
REWIND(LU_INPUT)
COUNT_TABLE_POINTS: DO N=1,N_TABLE
   TA => TABLES(N)
   REWIND(LU_INPUT)
   TA%NUMBER_ROWS = 0
   SELECT CASE (TABLE_TYPE(N))
      CASE (SPRAY_PATTERN)
         TA%NUMBER_COLUMNS = 6
   END SELECT
   SEARCH_LOOP: DO
      CALL CHECKREAD('TABL',LU_INPUT,IOS)
      IF (IOS==1) EXIT SEARCH_LOOP
      TABLE_DATA = -999._EB
      READ(LU_INPUT,NML=TABL,ERR=56,IOSTAT=IOS)
      IF (ID/=TABLE_ID(N)) CYCLE SEARCH_LOOP
      TA%NUMBER_ROWS = TA%NUMBER_ROWS + 1
      SELECT CASE(TABLE_TYPE(N))
         CASE (SPRAY_PATTERN)
            MESSAGE='null'
            IF (TABLE_DATA(1)<0. .OR.           TABLE_DATA(1)>180) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad 1st lattitude'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(2)<TABLE_DATA(1).OR. TABLE_DATA(2)>180) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad 2nd lattitude'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(3)<-180. .OR.        TABLE_DATA(3)>360) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad 1st longitude'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(4)<TABLE_DATA(3).OR. TABLE_DATA(4)>360) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad 2nd longitude'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(5)<0) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad velocity'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(6)<0) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad mass flow'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
      END SELECT
         
      56 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with TABLE '//TRIM(TABLE_ID(N)) )
   ENDDO SEARCH_LOOP
   IF (TA%NUMBER_ROWS==0) THEN
      WRITE(MESSAGE,'(A,A,A)') 'ERROR: TABLE ',TRIM(TABLE_ID(N)), ' not found'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO COUNT_TABLE_POINTS

READ_TABL_LOOP: DO N=1,N_TABLE
   TA => TABLES(N)
   ALLOCATE(TA%TABLE_DATA(TA%NUMBER_ROWS,TA%NUMBER_COLUMNS),STAT=IZERO)
   CALL ChkMemErr('READ','TA%TABLE_DATA',IZERO)
   REWIND(LU_INPUT)
   NN = 0
   SEARCH_LOOP2: DO 
      CALL CHECKREAD('TABL',LU_INPUT,IOS) 
      IF (IOS==1) EXIT SEARCH_LOOP2
      READ(LU_INPUT,TABL)
      IF (ID/=TABLE_ID(N)) CYCLE SEARCH_LOOP2
      NN = NN+1
      TA%TABLE_DATA(NN,:) = TABLE_DATA(1:TA%NUMBER_COLUMNS)
   ENDDO SEARCH_LOOP2
ENDDO READ_TABL_LOOP

END SUBROUTINE READ_TABL

 
SUBROUTINE READ_OBST

USE GEOMETRY_FUNCTIONS, ONLY: BLOCK_CELL
USE DEVICE_VARIABLES, ONLY : DEVICE, N_DEVC
USE CONTROL_VARIABLES, ONLY : CONTROL, N_CTRL 
TYPE(OBSTRUCTION_TYPE), POINTER :: OB2,OBT
TYPE(MULTIPLIER_TYPE), POINTER :: MR
TYPE(OBSTRUCTION_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: TEMP_OBSTRUCTION
INTEGER :: NM,NOM,N_OBST_O,NNN,IC,N,NN,NNNN,N_NEW_OBST,RGB(3),N_OBST_NEW,II,JJ,KK
CHARACTER(30) :: DEVC_ID,SURF_ID,SURF_IDS(3),SURF_ID6(6),CTRL_ID,MULT_ID
CHARACTER(60) :: MESH_ID
CHARACTER(25) :: COLOR
REAL(EB) :: TRANSPARENCY,DUMMY,XB1,XB2,XB3,XB4,XB5,XB6,BULK_DENSITY,VOL_ADJUSTED,VOL_SPECIFIED
LOGICAL :: SAWTOOTH,EMBEDDED,THICKEN,PERMIT_HOLE,ALLOW_VENT,EVACUATION, REMOVABLE,BNDF_FACE(-3:3),BNDF_OBST,OUTLINE
NAMELIST /OBST/ XB,SURF_ID,SURF_IDS,SURF_ID6,FYI,BNDF_FACE,BNDF_OBST, &
                SAWTOOTH,RGB,TRANSPARENCY,TEXTURE_ORIGIN,THICKEN, OUTLINE,DEVC_ID,CTRL_ID,COLOR, &
                PERMIT_HOLE,ALLOW_VENT,EVACUATION,MESH_ID,REMOVABLE,MULT_ID,BULK_DENSITY
 
MESH_LOOP: DO NM=1,NMESHES

   IF (PROCESS(NM)/=MYID .AND. MYID/=EVAC_PROCESS) CYCLE MESH_LOOP

   M=>MESHES(NM)
   CALL POINT_TO_MESH(NM)
 
   ! Count OBST lines
 
   REWIND(LU_INPUT)
   N_OBST = 0
   COUNT_OBST_LOOP: DO
      CALL CHECKREAD('OBST',LU_INPUT,IOS)
      IF (IOS==1) EXIT COUNT_OBST_LOOP
      MULT_ID = 'null'
      READ(LU_INPUT,NML=OBST,END=1,ERR=2,IOSTAT=IOS)
      N_OBST_NEW = 1
      IF (MULT_ID/='null') THEN
         DO N=1,N_MULT
            MR => MULTIPLIER(N)
            IF (MULT_ID==MR%ID) N_OBST_NEW = MR%N_COPIES
         ENDDO
      ENDIF
      N_OBST = N_OBST + N_OBST_NEW
      2 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I5)')  'ERROR: Problem with OBSTruction number',N_OBST+1
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDDO COUNT_OBST_LOOP
   1 REWIND(LU_INPUT)
 
   ! Allocate OBSTRUCTION array

   ALLOCATE(M%OBSTRUCTION(0:N_OBST),STAT=IZERO)
   CALL ChkMemErr('READ','OBSTRUCTION',IZERO)
   OBSTRUCTION=>M%OBSTRUCTION
 
   N        = 0
   N_OBST_O = N_OBST
 
   READ_OBST_LOOP: DO NN=1,N_OBST_O

      MULT_ID  = 'null'
      SURF_ID  = 'null'
      SURF_IDS = 'null'
      SURF_ID6 = 'null'
      COLOR    = 'null'
      MESH_ID     = 'null'
      RGB         = -1
      BULK_DENSITY= -1._EB
      TRANSPARENCY=  1._EB
      BNDF_FACE   = BNDF_DEFAULT
      BNDF_OBST   = BNDF_DEFAULT
      SAWTOOTH    = .TRUE.
      THICKEN     = THICKEN_OBSTRUCTIONS
      OUTLINE     = .FALSE.
      TEXTURE_ORIGIN = -999._EB
      DEVC_ID     = 'null'
      CTRL_ID     = 'null'
      PERMIT_HOLE = .TRUE.
      ALLOW_VENT  = .TRUE.
      REMOVABLE   = .TRUE.
      IF (.NOT.EVACUATION_ONLY(NM)) EVACUATION = .FALSE.
      IF (     EVACUATION_ONLY(NM)) EVACUATION = .TRUE.
      IF (     EVACUATION_ONLY(NM)) REMOVABLE = .FALSE.
 
      ! Read the OBST line

      CALL CHECKREAD('OBST',LU_INPUT,IOS)
      IF (IOS==1) EXIT READ_OBST_LOOP
      READ(LU_INPUT,OBST,END=35)

      ! Reorder OBST coordinates if necessary

      DO I=1,5,2
         IF (XB(I)>XB(I+1)) THEN
            DUMMY   = XB(I)
            XB(I)   = XB(I+1)
            XB(I+1) = DUMMY
         ENDIF
      ENDDO

      ! Loop over all possible multiples of the OBST

      MR => MULTIPLIER(0)
      DO NNN=1,N_MULT
         IF (MULT_ID==MULTIPLIER(NNN)%ID) MR => MULTIPLIER(NNN)
      ENDDO

      K_MULT_LOOP: DO KK=MR%K_LOWER,MR%K_UPPER
         J_MULT_LOOP: DO JJ=MR%J_LOWER,MR%J_UPPER
            I_MULT_LOOP: DO II=MR%I_LOWER,MR%I_UPPER
 
               IF (.NOT.MR%SEQUENTIAL) THEN
                  XB1 = XB(1) + MR%DX0 + II*MR%DXB(1)
                  XB2 = XB(2) + MR%DX0 + II*MR%DXB(2)
                  XB3 = XB(3) + MR%DY0 + JJ*MR%DXB(3)
                  XB4 = XB(4) + MR%DY0 + JJ*MR%DXB(4)
                  XB5 = XB(5) + MR%DZ0 + KK*MR%DXB(5)
                  XB6 = XB(6) + MR%DZ0 + KK*MR%DXB(6)
               ELSE
                  XB1 = XB(1) + MR%DX0 + II*MR%DXB(1)
                  XB2 = XB(2) + MR%DX0 + II*MR%DXB(2)
                  XB3 = XB(3) + MR%DY0 + II*MR%DXB(3)
                  XB4 = XB(4) + MR%DY0 + II*MR%DXB(4)
                  XB5 = XB(5) + MR%DZ0 + II*MR%DXB(5)
                  XB6 = XB(6) + MR%DZ0 + II*MR%DXB(6)
               ENDIF

               ! Increase the OBST counter

               N = N + 1
 
               ! Evacuation criteria
 
               IF (MESH_ID/=MESH_NAME(NM) .AND. MESH_ID/='null') THEN
                     N = N-1
                     N_OBST = N_OBST-1
                     CYCLE I_MULT_LOOP
               ENDIF
 
               IF ((.NOT.EVACUATION .AND. EVACUATION_ONLY(NM)) .OR. (EVACUATION .AND. .NOT.EVACUATION_ONLY(NM))) THEN
                     N = N-1
                     N_OBST = N_OBST-1
                     CYCLE I_MULT_LOOP
               ENDIF
 
               ! Include obstructions within half a grid cell of the computational boundary

               IF (XB2>=XS-0.5_EB*DX(0)    .AND. XB2<XS) THEN
                  XB1 = XS
                  XB2 = XS
               ENDIF
               IF (XB1< XF+0.5_EB*DX(IBP1) .AND. XB1>XF) THEN
                  XB1 = XF
                  XB2 = XF
               ENDIF
               IF (XB4>=YS-0.5_EB*DY(0)    .AND. XB4<YS) THEN
                  XB3 = YS
                  XB4 = YS
               ENDIF
               IF (XB3< YF+0.5_EB*DY(JBP1) .AND. XB3>YF) THEN
                  XB3 = YF
                  XB4 = YF
               ENDIF
               IF (XB6>=ZS-0.5_EB*DZ(0)    .AND. XB6<ZS .AND. .NOT.EVACUATION_ONLY(NM)) THEN
                  XB5 = ZS
                  XB6 = ZS
               ENDIF
               IF (XB5< ZF+0.5_EB*DZ(KBP1) .AND. XB5>ZF .AND. .NOT.EVACUATION_ONLY(NM)) THEN
                  XB5 = ZF
                  XB6 = ZF
               ENDIF
 
               ! Throw out obstructions that are not within computational domain

               XB1 = MAX(XB1,XS)
               XB2 = MIN(XB2,XF)
               XB3 = MAX(XB3,YS)
               XB4 = MIN(XB4,YF)
               XB5 = MAX(XB5,ZS)
               XB6 = MIN(XB6,ZF)
               IF (XB1>XF .OR. XB2<XS .OR. XB3>YF .OR. XB4<YS .OR. XB5>ZF .OR. XB6<ZS) THEN
                  N = N-1
                  N_OBST = N_OBST-1
                  CYCLE I_MULT_LOOP
               ENDIF
 
               ! Begin processing of OBSTruction
 
               OB=>OBSTRUCTION(N)
 
               OB%X1 = XB1
               OB%X2 = XB2
               OB%Y1 = XB3
               OB%Y2 = XB4
               OB%Z1 = XB5
               OB%Z2 = XB6
               OB%I1 = NINT( GINV(XB1-XS,1,NM)*RDXI   ) 
               OB%I2 = NINT( GINV(XB2-XS,1,NM)*RDXI   )
               OB%J1 = NINT( GINV(XB3-YS,2,NM)*RDETA  ) 
               OB%J2 = NINT( GINV(XB4-YS,2,NM)*RDETA  )
               OB%K1 = NINT( GINV(XB5-ZS,3,NM)*RDZETA ) 
               OB%K2 = NINT( GINV(XB6-ZS,3,NM)*RDZETA )
 
               ! If desired, thicken small obstructions
 
               IF (THICKEN .AND. OB%I1==OB%I2) THEN
                  OB%I1 = GINV(.5_EB*(XB1+XB2)-XS,1,NM)*RDXI
                  OB%I2 = MIN(OB%I1+1,IBAR)
               ENDIF
               IF (THICKEN .AND. OB%J1==OB%J2) THEN
                  OB%J1 = GINV(.5_EB*(XB3+XB4)-YS,2,NM)*RDETA
                  OB%J2 = MIN(OB%J1+1,JBAR)
               ENDIF
               IF (THICKEN .AND. OB%K1==OB%K2) THEN
                  OB%K1 = GINV(.5_EB*(XB5+XB6)-ZS,3,NM)*RDZETA
                  OB%K2 = MIN(OB%K1+1,KBAR)
               ENDIF

               ! Thicken evacuation mesh obstructions in the z direction

               IF (EVACUATION_ONLY(NM) .AND. EVACUATION .AND. OB%K1==OB%K2) THEN
                  OB%K1 = GINV(.5_EB*(XB5+XB6)-ZS,3,NM)*RDZETA
                  OB%K2 = MIN(OB%K1+1,KBAR)
               ENDIF
 
               ! Throw out obstructions that are too small
 
               IF ((OB%I1==OB%I2.AND.OB%J1==OB%J2) .OR. (OB%I1==OB%I2.AND.OB%K1==OB%K2) .OR. (OB%J1==OB%J2.AND.OB%K1==OB%K2)) THEN
                  N = N-1
                  N_OBST= N_OBST-1
                  CYCLE I_MULT_LOOP
               ENDIF

               IF (OB%I1==OB%I2 .OR. OB%J1==OB%J2 .OR. OB%K1==OB%K2) OB%THIN = .TRUE.
 
               ! Check to see if obstacle is completely embedded in another
 
               EMBEDDED = .FALSE.
               EMBED_LOOP: DO NNN=1,N-1
                  OB2=>OBSTRUCTION(NNN)
                  IF (OB%I1>OB2%I1 .AND. OB%I2<OB2%I2 .AND. &
                      OB%J1>OB2%J1 .AND. OB%J2<OB2%J2 .AND. &
                      OB%K1>OB2%K1 .AND. OB%K2<OB2%K2) THEN
                     EMBEDDED = .TRUE.
                     EXIT EMBED_LOOP
                  ENDIF
               ENDDO EMBED_LOOP
 
               IF (EMBEDDED  .AND. DEVC_ID=='null' .AND.  REMOVABLE .AND. CTRL_ID=='null' ) THEN
                     N = N-1
                     N_OBST= N_OBST-1
                     CYCLE I_MULT_LOOP
               ENDIF

               ! Check if the SURF IDs exist
         
               IF (EVACUATION_ONLY(NM)) SURF_ID=EVAC_SURF_DEFAULT
         
               IF (SURF_ID/='null') CALL CHECK_SURF_NAME(SURF_ID,EX)
               IF (.NOT.EX) THEN
                  WRITE(MESSAGE,'(A,A,A)')  'ERROR: SURF_ID ',TRIM(SURF_ID),' does not exist'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
         
               DO NNNN=1,3
                  IF (EVACUATION_ONLY(NM)) SURF_IDS(NNNN)=EVAC_SURF_DEFAULT
                  IF (SURF_IDS(NNNN)/='null') CALL CHECK_SURF_NAME(SURF_IDS(NNNN),EX)
                  IF (.NOT.EX) THEN
                     WRITE(MESSAGE,'(A,A,A)')  'ERROR: SURF_ID ',TRIM(SURF_IDS(NNNN)),' does not exist'
                     CALL SHUTDOWN(MESSAGE)
                  ENDIF
               ENDDO

               DO NNNN=1,6
                  IF (EVACUATION_ONLY(NM)) SURF_ID6(NNNN)=EVAC_SURF_DEFAULT
                  IF (SURF_ID6(NNNN)/='null') CALL CHECK_SURF_NAME(SURF_ID6(NNNN),EX)
                  IF (.NOT.EX) THEN
                     WRITE(MESSAGE,'(A,A,A)')  'ERROR: SURF_ID ',TRIM(SURF_ID6(NNNN)),' does not exist'
                     CALL SHUTDOWN(MESSAGE)
                  ENDIF
               ENDDO
 
               ! Save boundary condition info for obstacles
 
               OB%IBC(:) = DEFAULT_SURF_INDEX
          
               DO NNN=0,N_SURF
                  IF (SURF_ID    ==SURFACE(NNN)%ID) OB%IBC(:)    = NNN
                  IF (SURF_IDS(1)==SURFACE(NNN)%ID) OB%IBC(3)    = NNN
                  IF (SURF_IDS(2)==SURFACE(NNN)%ID) OB%IBC(-2:2) = NNN
                  IF (SURF_IDS(3)==SURFACE(NNN)%ID) OB%IBC(-3)   = NNN
                  IF (SURF_ID6(1)==SURFACE(NNN)%ID) OB%IBC(-1)   = NNN
                  IF (SURF_ID6(2)==SURFACE(NNN)%ID) OB%IBC( 1)   = NNN
                  IF (SURF_ID6(3)==SURFACE(NNN)%ID) OB%IBC(-2)   = NNN
                  IF (SURF_ID6(4)==SURFACE(NNN)%ID) OB%IBC( 2)   = NNN
                  IF (SURF_ID6(5)==SURFACE(NNN)%ID) OB%IBC(-3)   = NNN
                  IF (SURF_ID6(6)==SURFACE(NNN)%ID) OB%IBC( 3)   = NNN
               ENDDO
         
               ! Determine if the OBST is CONSUMABLE and check if POROUS inappropriately applied
         
               FACE_LOOP: DO NNN=-3,3
                  IF (NNN==0) CYCLE FACE_LOOP
                  IF (SURFACE(OB%IBC(NNN))%BURN_AWAY) THEN
                     OB%CONSUMABLE = .TRUE.
                     IF (.NOT.SAWTOOTH) THEN
                        WRITE(MESSAGE,'(A,I5,A)')  'ERROR: OBST number',N,' cannot have a BURN_AWAY SURF_ID and SAWTOOTH=.FALSE.' 
                        CALL SHUTDOWN(MESSAGE)
                     ENDIF
                  ENDIF
                  IF (SURFACE(OB%IBC(NNN))%POROUS .AND. .NOT.OB%THIN) THEN
                     WRITE(MESSAGE,'(A,I5,A)')  'ERROR: OBST number',N,' must be zero cells thick if it is to be POROUS'
                     CALL SHUTDOWN(MESSAGE)
                  ENDIF
               ENDDO FACE_LOOP

               ! Calculate the increase or decrease in the obstruction volume over the user-specified

               VOL_SPECIFIED = (OB%X2-OB%X1)*(OB%Y2-OB%Y1)*(OB%Z2-OB%Z1)
               VOL_ADJUSTED  = (X(OB%I2)-X(OB%I1))*(Y(OB%J2)-Y(OB%J1))*(Z(OB%K2)-Z(OB%K1))
               IF (VOL_SPECIFIED>0._EB .AND..NOT.EVACUATION_ONLY(NM)) THEN
                  OB%VOLUME_ADJUST = VOL_ADJUSTED/VOL_SPECIFIED
               ELSE
                  OB%VOLUME_ADJUST = 0._EB
               ENDIF
          
               ! Creation and removal logic
          
               OB%DEVC_ID = DEVC_ID
               OB%CTRL_ID = CTRL_ID
               OB%HIDDEN = .FALSE.
         
               CALL SEARCH_CONTROLLER('OBST',CTRL_ID,DEVC_ID,OB%DEVC_INDEX,OB%CTRL_INDEX,N)
               IF (DEVC_ID /='null') THEN
                  IF (.NOT.DEVICE(OB%DEVC_INDEX)%INITIAL_STATE) OB%HIDDEN = .TRUE.
                  OB%REMOVABLE = .TRUE.
               ENDIF
               IF (CTRL_ID /='null') THEN
                  IF (.NOT.CONTROL(OB%CTRL_INDEX)%INITIAL_STATE) OB%HIDDEN = .TRUE.
                  OB%REMOVABLE = .TRUE.
               ENDIF
         
               IF (OB%CONSUMABLE .AND..NOT.EVACUATION_ONLY(NM))    OB%REMOVABLE = .TRUE.      

               ! Choose obstruction color index

               SELECT CASE (COLOR)
                  CASE ('INVISIBLE')
                     OB%COLOR_INDICATOR = -3
                     RGB(1) = 255
                     RGB(2) = 204
                     RGB(3) = 102
                     TRANSPARENCY = 0._EB
                  CASE ('null')
                     IF (ANY (RGB<0)) THEN
                        OB%COLOR_INDICATOR = -1
                     ELSE
                        OB%COLOR_INDICATOR = -3
                     ENDIF
                  CASE DEFAULT
                     CALL COLOR2RGB(RGB,COLOR)
                     OB%COLOR_INDICATOR = -3
               END SELECT
               OB%RGB  = RGB
               OB%TRANSPARENCY = TRANSPARENCY

               ! Miscellaneous assignments
 
               OB%TEXTURE(:) = TEXTURE_ORIGIN(:)  ! Origin of texture map
               OB%ORDINAL = NN  ! Order of OBST in original input file
               OB%PERMIT_HOLE = PERMIT_HOLE
               OB%ALLOW_VENT  = ALLOW_VENT

               ! Only allow the use of BULK_DENSITY if the obstruction has a non-zero volume

               OB%BULK_DENSITY = BULK_DENSITY
               IF (OB%VOLUME_ADJUST==0._EB .AND. BULK_DENSITY>0._EB) THEN
                  WRITE(MESSAGE,'(A,I4,A)') 'ERROR: OBSTruction ',NN,' has no volume and thus cannot have a BULK_DENSITY'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
 
               ! Make obstruction invisible if it's within a finer mesh
 
               DO NOM=1,NM-1
                  IF (XB1>MESHES(NOM)%XS .AND. XB2<MESHES(NOM)%XF .AND. &
                      XB3>MESHES(NOM)%YS .AND. XB4<MESHES(NOM)%YF .AND. &
                      XB5>MESHES(NOM)%ZS .AND. XB6<MESHES(NOM)%ZF) OB%COLOR_INDICATOR=-2
               ENDDO
 
               ! Prevent drawing of boundary info if desired
          
               IF (BNDF_DEFAULT) THEN
                  OB%SHOW_BNDF(:) = BNDF_FACE(:)
                  IF (.NOT.BNDF_OBST) OB%SHOW_BNDF(:) = .FALSE.
               ELSE
                  OB%SHOW_BNDF(:) = BNDF_FACE(:)
                  IF (BNDF_OBST) OB%SHOW_BNDF(:) = .TRUE.
               ENDIF
 
               ! Smooth obstacles if desired
          
               IF (.NOT.SAWTOOTH .AND..NOT.EVACUATION_ONLY(NM)) THEN
                  OB%TYPE_INDICATOR = 3
                  OB%SAWTOOTH = .FALSE.
               ENDIF
 
               ! In Smokeview, draw the outline of the obstruction
 
               IF (OUTLINE) OB%TYPE_INDICATOR = 2

            ENDDO I_MULT_LOOP
         ENDDO J_MULT_LOOP
      ENDDO K_MULT_LOOP
      
   ENDDO READ_OBST_LOOP
35 REWIND(LU_INPUT)
 
ENDDO MESH_LOOP
 
! Read HOLEs and cut out blocks
 
CALL READ_HOLE
 
! Look for OBSTructions that are meant to BURN_AWAY and break them up into single cell blocks

MESH_LOOP_2: DO NM=1,NMESHES

   IF (PROCESS(NM)/=MYID .AND. MYID/=EVAC_PROCESS) CYCLE MESH_LOOP_2

   M=>MESHES(NM)
   CALL POINT_TO_MESH(NM)

   N_OBST_O = N_OBST
   DO N=1,N_OBST_O
      OB => OBSTRUCTION(N)
      IF (OB%CONSUMABLE .AND..NOT.EVACUATION_ONLY(NM)) THEN

         N_NEW_OBST = MAX(1,OB%I2-OB%I1)*MAX(1,OB%J2-OB%J1)*MAX(1,OB%K2-OB%K1)
         IF (N_NEW_OBST > 1) THEN

            ! Create a temporary array of obstructions with the same properties as the one being replaced, except coordinates

            ALLOCATE(TEMP_OBSTRUCTION(N_NEW_OBST))
            TEMP_OBSTRUCTION = OBSTRUCTION(N)
            NN = 0
            DO K=OB%K1,MAX(OB%K1,OB%K2-1)
               DO J=OB%J1,MAX(OB%J1,OB%J2-1)
                  DO I=OB%I1,MAX(OB%I1,OB%I2-1)
                     NN = NN + 1
                     OBT=>TEMP_OBSTRUCTION(NN)
                     OBT%I1 = I
                     OBT%I2 = MIN(I+1,OB%I2)
                     OBT%J1 = J
                     OBT%J2 = MIN(J+1,OB%J2)
                     OBT%K1 = K
                     OBT%K2 = MIN(K+1,OB%K2)
                     OBT%X1 = M%X(OBT%I1)
                     OBT%X2 = M%X(OBT%I2)
                     OBT%Y1 = M%Y(OBT%J1)
                     OBT%Y2 = M%Y(OBT%J2)
                     OBT%Z1 = M%Z(OBT%K1)
                     OBT%Z2 = M%Z(OBT%K2)
                  ENDDO
                ENDDO
            ENDDO

            CALL RE_ALLOCATE_OBST(NM,N_OBST,N_NEW_OBST-1)
            OBSTRUCTION=>M%OBSTRUCTION
            OBSTRUCTION(N) = TEMP_OBSTRUCTION(1)
            OBSTRUCTION(N_OBST+1:N_OBST+N_NEW_OBST-1) = TEMP_OBSTRUCTION(2:N_NEW_OBST)
            N_OBST = N_OBST + N_NEW_OBST-1
            DEALLOCATE(TEMP_OBSTRUCTION)

         ENDIF
      ENDIF
   ENDDO

ENDDO MESH_LOOP_2


! Go through all meshes, recording which cells are solid
 
MESH_LOOP_3: DO NM=1,NMESHES

   IF (PROCESS(NM)/=MYID .AND. MYID/=EVAC_PROCESS) CYCLE MESH_LOOP_3

   M=>MESHES(NM)
   CALL POINT_TO_MESH(NM)
 
   ! Compute areas of obstruction faces, both actual (AB0) and FDS approximated (AB)
 
   DO N=1,N_OBST
      OB=>OBSTRUCTION(N)
      OB%INPUT_AREA(1) = (OB%Y2-OB%Y1)*(OB%Z2-OB%Z1)
      OB%INPUT_AREA(2) = (OB%X2-OB%X1)*(OB%Z2-OB%Z1)
      OB%INPUT_AREA(3) = (OB%X2-OB%X1)*(OB%Y2-OB%Y1)
      OB%FDS_AREA(1)   = (Y(OB%J2)-Y(OB%J1))*(Z(OB%K2)-Z(OB%K1))
      OB%FDS_AREA(2)   = (X(OB%I2)-X(OB%I1))*(Z(OB%K2)-Z(OB%K1))
      OB%FDS_AREA(3)   = (X(OB%I2)-X(OB%I1))*(Y(OB%J2)-Y(OB%J1))
      OB%DIMENSIONS(1) = OB%I2 - OB%I1
      OB%DIMENSIONS(2) = OB%J2 - OB%J1
      OB%DIMENSIONS(3) = OB%K2 - OB%K1
   ENDDO
 
   ! Create main blockage index array (ICA)
 
   ALLOCATE(M%CELL_INDEX(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO) 
   CALL ChkMemErr('READ','CELL_INDEX',IZERO) 
   CELL_INDEX=>M%CELL_INDEX
 
   CELL_INDEX = 0 
   CELL_COUNT = 0
 
   DO K=0,KBP1
      IF (EVACUATION_ONLY(NM) .AND. .NOT.(K==1)) CYCLE
      DO J=0,JBP1
         DO I=0,1
            IF (CELL_INDEX(I,J,K)==0) THEN
               CELL_COUNT = CELL_COUNT + 1
               CELL_INDEX(I,J,K) = CELL_COUNT
            ENDIF
         ENDDO
         DO I=IBAR,IBP1
            IF (CELL_INDEX(I,J,K)==0) THEN
               CELL_COUNT = CELL_COUNT + 1
               CELL_INDEX(I,J,K) = CELL_COUNT
            ENDIF
         ENDDO
      ENDDO
   ENDDO
 
   DO K=0,KBP1
      IF (EVACUATION_ONLY(NM) .AND. .NOT.(K==1)) CYCLE
      DO I=0,IBP1
         DO J=0,1
            IF (CELL_INDEX(I,J,K)==0) THEN
               CELL_COUNT = CELL_COUNT + 1
               CELL_INDEX(I,J,K) = CELL_COUNT
            ENDIF
         ENDDO
         DO J=JBAR,JBP1
            IF (CELL_INDEX(I,J,K)==0) THEN
               CELL_COUNT = CELL_COUNT + 1
               CELL_INDEX(I,J,K) = CELL_COUNT
            ENDIF
         ENDDO
      ENDDO
   ENDDO
 
   DO J=0,JBP1
      DO I=0,IBP1
         DO K=0,1
            IF (EVACUATION_ONLY(NM) .AND. .NOT.(K==1)) CYCLE
            IF (CELL_INDEX(I,J,K)==0) THEN
               CELL_COUNT = CELL_COUNT + 1
               CELL_INDEX(I,J,K) = CELL_COUNT
            ENDIF
         ENDDO
         DO K=KBAR,KBP1
            IF (EVACUATION_ONLY(NM) .AND. .NOT.(K==1)) CYCLE
            IF (CELL_INDEX(I,J,K)==0) THEN
               CELL_COUNT = CELL_COUNT + 1
               CELL_INDEX(I,J,K) = CELL_COUNT
            ENDIF
         ENDDO
      ENDDO
   ENDDO
 
   DO N=1,N_OBST
      OB=>OBSTRUCTION(N)
      DO K=OB%K1,OB%K2+1
         IF (EVACUATION_ONLY(NM) .AND. .NOT.(K==1)) CYCLE
         DO J=OB%J1,OB%J2+1
            DO I=OB%I1,OB%I2+1
               IF (CELL_INDEX(I,J,K)==0) THEN
                  CELL_COUNT = CELL_COUNT + 1
                  CELL_INDEX(I,J,K) = CELL_COUNT
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
 
   ! Store in SOLID which cells are solid and which are not
 
   ALLOCATE(M%SOLID(0:CELL_COUNT),STAT=IZERO) 
   CALL ChkMemErr('READ','SOLID',IZERO) 
   M%SOLID = .FALSE.
   ALLOCATE(M%EXTERIOR(0:CELL_COUNT),STAT=IZERO) 
   CALL ChkMemErr('READ','EXTERIOR',IZERO) 
   M%EXTERIOR = .FALSE.
   IF(.NOT.EVACUATION_ONLY(NM) .AND. MYID==EVAC_PROCESS) THEN
      ! EVAC_PROCESS needs only the allocations of SOLID and CELL_INDEX
      ! of the fire meshes.
      DEALLOCATE(M%OBSTRUCTION)
      CYCLE MESH_LOOP_3
   ENDIF
   SOLID=>M%SOLID
   ALLOCATE(M%OBST_INDEX_C(0:CELL_COUNT),STAT=IZERO) 
   CALL ChkMemErr('READ','OBST_INDEX_C',IZERO) 
   M%OBST_INDEX_C = 0
   OBST_INDEX_C=>M%OBST_INDEX_C 
 
   CALL BLOCK_CELL(NM,   0,   0,   0,JBP1,   0,KBP1,1,0)
   CALL BLOCK_CELL(NM,IBP1,IBP1,   0,JBP1,   0,KBP1,1,0)
   CALL BLOCK_CELL(NM,   0,IBP1,   0,   0,   0,KBP1,1,0)
   CALL BLOCK_CELL(NM,   0,IBP1,JBP1,JBP1,   0,KBP1,1,0)
   CALL BLOCK_CELL(NM,   0,IBP1,   0,JBP1,   0,   0,1,0)
   CALL BLOCK_CELL(NM,   0,IBP1,   0,JBP1,KBP1,KBP1,1,0)
 
   DO N=1,N_OBST
      OB=>OBSTRUCTION(N)
      CALL BLOCK_CELL(NM,OB%I1+1,OB%I2,OB%J1+1,OB%J2,OB%K1+1,OB%K2,1,N)
   ENDDO
 
   ALLOCATE(M%I_CELL(CELL_COUNT),STAT=IZERO) 
   CALL ChkMemErr('READ','I_CELL',IZERO) 
   M%I_CELL = -1
   ALLOCATE(M%J_CELL(CELL_COUNT),STAT=IZERO) 
   CALL ChkMemErr('READ','J_CELL',IZERO) 
   M%J_CELL = -1
   ALLOCATE(M%K_CELL(CELL_COUNT),STAT=IZERO) 
   CALL ChkMemErr('READ','K_CELL',IZERO) 
   M%K_CELL = -1
   I_CELL=>M%I_CELL 
   J_CELL=>M%J_CELL 
   K_CELL=>M%K_CELL
 
   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
         IC = CELL_INDEX(I,J,K)
            IF (IC>0) THEN
               I_CELL(IC) = I
               J_CELL(IC) = J
               K_CELL(IC) = K
            ENDIF
         ENDDO
      ENDDO
   ENDDO
 
ENDDO MESH_LOOP_3

END SUBROUTINE READ_OBST


SUBROUTINE READ_HOLE
USE CONTROL_VARIABLES, ONLY : CONTROl, N_CTRL
USE DEVICE_VARIABLES, ONLY : DEVICE, N_DEVC
CHARACTER(30) :: DEVC_ID,CTRL_ID
CHARACTER(60) :: MESH_ID
CHARACTER(25) :: COLOR
LOGICAL :: EVACUATION
INTEGER :: NM,N_HOLE,NN,NDO,N,I1,I2,J1,J2,K1,K2,RGB(3)
REAL(EB) :: X1,X2,Y1,Y2,Z1,Z2,TRANSPARENCY, DUMMY
NAMELIST /HOLE/ XB,FYI,RGB,TRANSPARENCY,EVACUATION,MESH_ID,COLOR,DEVC_ID,CTRL_ID
TYPE(OBSTRUCTION_TYPE), ALLOCATABLE, DIMENSION(:) :: TEMP_OBST
LOGICAL, ALLOCATABLE, DIMENSION(:) :: TEMP_HOLE_EVAC

ALLOCATE(TEMP_OBST(0:6))

N_HOLE  = 0
REWIND(LU_INPUT)
COUNT_LOOP: DO
   CALL CHECKREAD('HOLE',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_LOOP
   READ(LU_INPUT,NML=HOLE,END=1,ERR=2,IOSTAT=IOS)
   N_HOLE = N_HOLE + 1
   2 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I5)')  'ERROR: Problem with HOLE number',N_HOLE+1
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO COUNT_LOOP
1 REWIND(LU_INPUT)

IF(ANY(EVACUATION_ONLY)) THEN 
   ALLOCATE(TEMP_HOLE_EVAC(1:N_HOLE))
   READ_HOLE_EVAC_LOOP: DO N=1,N_HOLE
      EVACUATION = .TRUE.
      CALL CHECKREAD('HOLE',LU_INPUT,IOS)
      IF (IOS==1) EXIT READ_HOLE_EVAC_LOOP
      READ(LU_INPUT,HOLE)
      TEMP_HOLE_EVAC(1:N_HOLE) = EVACUATION
   ENDDO READ_HOLE_EVAC_LOOP
   REWIND(LU_INPUT)
   !EVAC: Now TEMP_HOLE_EVAC(:) is FALSE,  if EVACUATION=.FALSE.
   !EVAC: otherwise it is true (no evac given or evac=true is given)
ENDIF
 
READ_HOLE_LOOP: DO N=1,N_HOLE
 
   DEVC_ID  = 'null'
   CTRL_ID  = 'null'
   MESH_ID  = 'null'
   COLOR    = 'null'
   RGB      = -1
   TRANSPARENCY  = 1._EB
   EVACUATION = .FALSE.
 
   CALL CHECKREAD('HOLE',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_HOLE_LOOP
   READ(LU_INPUT,HOLE)
 
   DO I=1,5,2
      IF (XB(I)>XB(I+1)) THEN
         DUMMY   = XB(I)
         XB(I)   = XB(I+1)
         XB(I+1) = DUMMY
      ENDIF
   ENDDO
 
   MESH_LOOP: DO NM=1,NMESHES
      M=>MESHES(NM)
      CALL POINT_TO_MESH(NM)
 
      ! Evacuation criteria
 
      IF (MESH_ID/='null' .AND. MESH_ID/=MESH_NAME(NM)) CYCLE MESH_LOOP
      IF (EVACUATION .AND. .NOT.EVACUATION_ONLY(NM)) CYCLE MESH_LOOP
      IF (EVACUATION_ONLY(NM)) THEN 
         IF (.NOT.TEMP_HOLE_EVAC(N)) CYCLE MESH_LOOP
      ENDIF
!      IF ((.NOT.EVACUATION .AND. EVACUATION_ONLY(NM)) .OR. (EVACUATION .AND. .NOT.EVACUATION_ONLY(NM))) CYCLE MESH_LOOP
 
      ! Check if hole is contained within the current mesh
 
      X1 = XB(1)
      X2 = XB(2)
      Y1 = XB(3)
      Y2 = XB(4)
      Z1 = XB(5)
      Z2 = XB(6)
 
      IF (X1>=XF .OR. X2<=XS .OR. Y1>YF .OR. Y2<=YS .OR. Z1>ZF .OR. Z2<=ZS) CYCLE MESH_LOOP
 
      X1 = MAX(X1,XS-0.001_EB*DX(0))
      X2 = MIN(X2,XF+0.001_EB*DX(IBP1))
      Y1 = MAX(Y1,YS-0.001_EB*DY(0))
      Y2 = MIN(Y2,YF+0.001_EB*DY(JBP1))
      Z1 = MAX(Z1,ZS-0.001_EB*DZ(0))
      Z2 = MIN(Z2,ZF+0.001_EB*DZ(KBP1))
 
      I1 = NINT( GINV(XB(1)-XS,1,NM)*RDXI   ) 
      I2 = NINT( GINV(XB(2)-XS,1,NM)*RDXI   )
      J1 = NINT( GINV(XB(3)-YS,2,NM)*RDETA  ) 
      J2 = NINT( GINV(XB(4)-YS,2,NM)*RDETA  )
      K1 = NINT( GINV(XB(5)-ZS,3,NM)*RDZETA ) 
      K2 = NINT( GINV(XB(6)-ZS,3,NM)*RDZETA )
 
      NN=0
      OBST_LOOP: DO
         NN=NN+1
         IF (NN>N_OBST) EXIT OBST_LOOP
         OB=>OBSTRUCTION(NN)
         IF (.NOT.OB%PERMIT_HOLE) CYCLE OBST_LOOP
 
         ! TEMP_OBST(0) is the intersection of HOLE and OBST
 
         TEMP_OBST(0)    = OBSTRUCTION(NN)
 
         TEMP_OBST(0)%I1 = MAX(I1,OB%I1)
         TEMP_OBST(0)%I2 = MIN(I2,OB%I2)
         TEMP_OBST(0)%J1 = MAX(J1,OB%J1)
         TEMP_OBST(0)%J2 = MIN(J2,OB%J2)
         TEMP_OBST(0)%K1 = MAX(K1,OB%K1)
         TEMP_OBST(0)%K2 = MIN(K2,OB%K2)
 
         TEMP_OBST(0)%X1 = MAX(X1,OB%X1)
         TEMP_OBST(0)%X2 = MIN(X2,OB%X2)
         TEMP_OBST(0)%Y1 = MAX(Y1,OB%Y1)
         TEMP_OBST(0)%Y2 = MIN(Y2,OB%Y2)
         TEMP_OBST(0)%Z1 = MAX(Z1,OB%Z1)
         TEMP_OBST(0)%Z2 = MIN(Z2,OB%Z2)
 
         ! Ignore OBSTs that do not intersect with HOLE or are merely sliced by the hole.
 
         IF (TEMP_OBST(0)%I2-TEMP_OBST(0)%I1<0 .OR. TEMP_OBST(0)%J2-TEMP_OBST(0)%J1<0 .OR. &
             TEMP_OBST(0)%K2-TEMP_OBST(0)%K1<0) CYCLE OBST_LOOP
         IF (TEMP_OBST(0)%I2-TEMP_OBST(0)%I1==0) THEN 
            IF (OB%I1<TEMP_OBST(0)%I1 .OR.  OB%I2>TEMP_OBST(0)%I2) CYCLE OBST_LOOP
         ENDIF
         IF (TEMP_OBST(0)%J2-TEMP_OBST(0)%J1==0) THEN
            IF (OB%J1<TEMP_OBST(0)%J1 .OR.  OB%J2>TEMP_OBST(0)%J2) CYCLE OBST_LOOP
         ENDIF
         IF (TEMP_OBST(0)%K2-TEMP_OBST(0)%K1==0) THEN
            IF (OB%K1<TEMP_OBST(0)%K1 .OR.  OB%K2>TEMP_OBST(0)%K2) CYCLE OBST_LOOP
         ENDIF
 
         IF (TEMP_OBST(0)%X2<=X1 .OR. TEMP_OBST(0)%X1>=X2 .OR. TEMP_OBST(0)%Y2<=Y1 .OR. TEMP_OBST(0)%Y1>=Y2 .OR. &
            TEMP_OBST(0)%Z2<=Z1 .OR. TEMP_OBST(0)%Z1>=Z2)  CYCLE OBST_LOOP
 
         ! Start counting new OBSTs that need to be created
 
         NDO=0
 
         IF ((OB%I1<I1.AND.I1<OB%I2) .OR. (XB(1)>=XS.AND.I1==0.AND.OB%I1==0)) THEN
            NDO=NDO+1
            TEMP_OBST(NDO)=OBSTRUCTION(NN)
            TEMP_OBST(NDO)%I1 = OB%I1
            TEMP_OBST(NDO)%I2 = I1
            TEMP_OBST(NDO)%X1 = OB%X1
            TEMP_OBST(NDO)%X2 = X1
         ENDIF
 
         IF ((OB%I1<I2.AND.I2<OB%I2) .OR. (XB(2)<=XF.AND.I2==IBAR.AND.OB%I2==IBAR)) THEN
            NDO=NDO+1
            TEMP_OBST(NDO)=OBSTRUCTION(NN)
            TEMP_OBST(NDO)%I1 = I2
            TEMP_OBST(NDO)%I2 = OB%I2
            TEMP_OBST(NDO)%X1 = X2 
            TEMP_OBST(NDO)%X2 = OB%X2
         ENDIF
 
         IF ((OB%J1<J1.AND.J1<OB%J2) .OR. (XB(3)>=YS.AND.J1==0.AND.OB%J1==0)) THEN
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
 
         IF ((OB%J1<J2.AND.J2<OB%J2) .OR. (XB(4)<=YF.AND.J2==JBAR.AND.OB%J2==JBAR)) THEN
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
 
         IF ((OB%K1<K1.AND.K1<OB%K2) .OR. (XB(5)>=ZS.AND.K1==0.AND.OB%K1==0)) THEN
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
 
         IF ((OB%K1<K2.AND.K2<OB%K2) .OR. (XB(6)<=ZF.AND.K2==KBAR.AND.OB%K2==KBAR)) THEN
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
 
         ! Maintain ordinal rank of original obstruction, but negate it. This will be a code for Smokeview.
 
         TEMP_OBST(:)%ORDINAL = -OB%ORDINAL
 
         ! Re-allocate space of new OBSTs, or remove entry for dead OBST
 
         NEW_OBST_IF: IF (NDO>0) THEN
               CALL RE_ALLOCATE_OBST(NM,N_OBST,NDO)
               OBSTRUCTION=>M%OBSTRUCTION
               OBSTRUCTION(N_OBST+1:N_OBST+NDO) = TEMP_OBST(1:NDO)
               N_OBST = N_OBST + NDO
         ENDIF NEW_OBST_IF
 
         ! If the HOLE is to be created or removed, save it in OBSTRUCTION(NN), the original OBST that was broken up

         DEVC_OR_CTRL: IF (DEVC_ID/='null' .OR. CTRL_ID/='null') THEN

            OBSTRUCTION(NN) = TEMP_OBST(0)
            OB => OBSTRUCTION(NN)
            OB%DEVC_ID = DEVC_ID
            OB%CTRL_ID = CTRL_ID
            CALL SEARCH_CONTROLLER('HOLE',CTRL_ID,DEVC_ID,OB%DEVC_INDEX,OB%CTRL_INDEX,N)
            IF (DEVC_ID /='null') THEN
               OB%REMOVABLE   = .TRUE.
               OB%HOLE_FILLER = .TRUE.
               OB%CTRL_INDEX = -1
               IF (DEVICE(OB%DEVC_INDEX)%INITIAL_STATE) OB%HIDDEN = .TRUE.
            ENDIF
            IF (CTRL_ID /='null') THEN
               OB%REMOVABLE   = .TRUE.
               OB%HOLE_FILLER = .TRUE.
               OB%DEVC_INDEX = -1
               IF (CONTROL(OB%CTRL_INDEX)%INITIAL_STATE) OB%HIDDEN = .TRUE.
            ENDIF
            
            IF (OB%CONSUMABLE)    OB%REMOVABLE = .TRUE.

            SELECT CASE (COLOR)
               CASE ('INVISIBLE')
                  OB%COLOR_INDICATOR = -3
                  OB%RGB(1) = 255
                  OB%RGB(2) = 204
                  OB%RGB(3) = 102
                  OB%TRANSPARENCY = 0._EB
               CASE ('null')
                  IF (ANY(RGB>0)) THEN
                     OB%COLOR_INDICATOR = -3
                     OB%RGB  = RGB
                     OB%TRANSPARENCY = TRANSPARENCY
                  ENDIF
               CASE DEFAULT
                  CALL COLOR2RGB(RGB,COLOR)
                  OB%COLOR_INDICATOR = -3
                  OB%RGB  = RGB
                  OB%TRANSPARENCY = TRANSPARENCY
            END SELECT

         ELSE DEVC_OR_CTRL

            OBSTRUCTION(NN) = OBSTRUCTION(N_OBST)
            N_OBST = N_OBST-1
            NN = NN-1

         ENDIF DEVC_OR_CTRL
 
      ENDDO OBST_LOOP
   ENDDO MESH_LOOP
ENDDO READ_HOLE_LOOP
 
REWIND(LU_INPUT)

IF(ANY(EVACUATION_ONLY)) DEALLOCATE(TEMP_HOLE_EVAC)
DEALLOCATE(TEMP_OBST)
END SUBROUTINE READ_HOLE
 
 
SUBROUTINE RE_ALLOCATE_OBST(NM,N_OBST,NDO)
TYPE (OBSTRUCTION_TYPE), ALLOCATABLE, DIMENSION(:) :: DUMMY
INTEGER, INTENT(IN) :: NM,NDO,N_OBST
TYPE (MESH_TYPE), POINTER :: M
M=>MESHES(NM)
ALLOCATE(DUMMY(0:N_OBST))
DUMMY(0:N_OBST) = M%OBSTRUCTION(0:N_OBST)
DEALLOCATE(M%OBSTRUCTION)
ALLOCATE(M%OBSTRUCTION(0:N_OBST+NDO))
M%OBSTRUCTION(0:N_OBST) = DUMMY(0:N_OBST)
DEALLOCATE(DUMMY)
END SUBROUTINE RE_ALLOCATE_OBST
 
 
SUBROUTINE READ_VENT

USE GEOMETRY_FUNCTIONS, ONLY : BLOCK_CELL
USE DEVICE_VARIABLES, ONLY : DEVICE, N_DEVC
USE CONTROL_VARIABLES, ONLY : CONTROL, N_CTRL
USE MATH_FUNCTIONS, ONLY: GET_RAMP_INDEX

INTEGER :: N,NN,NM,NNN,NVO,IOR,I1,I2,J1,J2,K1,K2,RGB(3)
REAL(EB) :: SPREAD_RATE,TRANSPARENCY,DUMMY,XYZ(3),TMP_EXTERIOR,MASS_FRACTION(MAX_SPECIES),DYNAMIC_PRESSURE
CHARACTER(30) :: DEVC_ID,CTRL_ID,SURF_ID,PRESSURE_RAMP
CHARACTER(60) :: MESH_ID
CHARACTER(25) :: COLOR
LOGICAL :: REJECT_VENT,EVACUATION,OUTLINE
NAMELIST /VENT/ XB,IOR,MB,PBX,PBY,PBZ,SURF_ID,FYI,RGB,TRANSPARENCY,COLOR, &
                TEXTURE_ORIGIN,OUTLINE,DEVC_ID,CTRL_ID, &
                XYZ,EVACUATION,MESH_ID,SPREAD_RATE,TMP_EXTERIOR,MASS_FRACTION,DYNAMIC_PRESSURE,PRESSURE_RAMP
 
MESH_LOOP_1: DO NM=1,NMESHES

!  IF (PROCESS(NM)/=MYID) CYCLE MESH_LOOP_1

   M=>MESHES(NM)
   CALL POINT_TO_MESH(NM)
 
   REWIND(LU_INPUT)
   N_VENT = 0
   COUNT_VENT_LOOP: DO
      CALL CHECKREAD('VENT',LU_INPUT,IOS) 
      IF (IOS==1) EXIT COUNT_VENT_LOOP
      READ(LU_INPUT,NML=VENT,END=3,ERR=4,IOSTAT=IOS)
      N_VENT = N_VENT + 1
      4 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I4)') 'ERROR: Problem with VENT ',N_VENT+1
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDDO COUNT_VENT_LOOP
   3 REWIND(LU_INPUT)
 
   IF (TWO_D)                         N_VENT = N_VENT + 2
   IF (CYLINDRICAL .AND. M%XS==0._EB) N_VENT = N_VENT + 1
   IF (EVACUATION_ONLY(NM))           N_VENT = N_VENT + 2
 
   ALLOCATE(M%VENTS(N_VENT),STAT=IZERO)
   CALL ChkMemErr('READ','VENTS',IZERO)
   VENTS=>M%VENTS
 
   NVO   = N_VENT
   N     = 0
 
   REWIND(LU_INPUT)
   READ_VENT_LOOP: DO NN=1,NVO
 
      N       = N + 1
      IOR     = 0
      MB      = 'null'
      PBX     = -1.E6_EB
      PBY     = -1.E6_EB
      PBZ     = -1.E6_EB
      SURF_ID = 'null'
      COLOR   = 'null'
      MESH_ID = 'null'
      RGB     =-1
      TRANSPARENCY = 1._EB
      DYNAMIC_PRESSURE = 0._EB
      PRESSURE_RAMP = 'null'
      XYZ     = -1.E6_EB
      SPREAD_RATE = 0.05_EB
      TMP_EXTERIOR = -1000.
      MASS_FRACTION = -1._EB
      REJECT_VENT  = .FALSE.
      TEXTURE_ORIGIN = -999._EB
      OUTLINE      = .FALSE.
      DEVC_ID  = 'null'
      CTRL_ID  = 'null'
      EVACUATION = .FALSE.
 
      IF (NN==NVO-2 .AND. CYLINDRICAL .AND. XS==0._EB) MB='XMIN'
      IF (NN==NVO-1 .AND. TWO_D)                       MB='YMIN'
      IF (NN==NVO   .AND. TWO_D)                       MB='YMAX'
      IF (NN==NVO-1 .AND. EVACUATION_ONLY(NM))         MB='ZMIN'
      IF (NN==NVO   .AND. EVACUATION_ONLY(NM))         MB='ZMAX'
 
      IF (MB=='null') THEN
         CALL CHECKREAD('VENT',LU_INPUT,IOS) 
         IF (IOS==1) EXIT READ_VENT_LOOP
         READ(LU_INPUT,VENT,END=37,ERR=38)    ! Read in info for VENT N
      ELSE
         SURF_ID = 'MIRROR'
      ENDIF
 
      IF (PBX>-1.E5_EB .OR. PBY>-1.E5_EB .OR. PBZ>-1.E5_EB) THEN
         XB(1) = XS
         XB(2) = XF
         XB(3) = YS
         XB(4) = YF
         XB(5) = ZS
         XB(6) = ZF
         IF (PBX>-1.E5_EB) XB(1:2) = PBX
         IF (PBY>-1.E5_EB) XB(3:4) = PBY
         IF (PBZ>-1.E5_EB) XB(5:6) = PBZ
      ENDIF
 
      IF (MB/='null') THEN
         XB(1) = XS
         XB(2) = XF
         XB(3) = YS
         XB(4) = YF
         XB(5) = ZS
         XB(6) = ZF
         SELECT CASE (MB)
            CASE('XMIN')
                XB(2) = XS
            CASE('XMAX')
                XB(1) = XF            
            CASE('YMIN')
                XB(4) = YS            
            CASE('YMAX')
                XB(3) = YF
            CASE('ZMIN')
                XB(6) = ZS                                                                
            CASE('ZMAX')      
                XB(5) = ZF                                                                      
            CASE DEFAULT
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: MB specified for VENT',NN,' is not XMIN, XMAX, YMIN, YMAX, ZMIN, or ZMAX'
               CALL SHUTDOWN(MESSAGE)
         END SELECT
      ENDIF
 
      ! Check that the vent is properly specified
 
      IF (MESH_ID/='null' .AND. MESH_ID/=MESH_NAME(NM))  REJECT_VENT = .TRUE.
 
      IF (XB(3)==XB(4) .AND. TWO_D .AND. NN<NVO-1) THEN
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: VENT',NN,' cannot be specified on a y boundary in a 2-D calculation'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
 
      IF (XB(1)/=XB(2) .AND. XB(3)/=XB(4) .AND. XB(5)/=XB(6)) THEN
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: VENT',NN,' must be a plane'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
 
      DO I=1,5,2
         IF (XB(I)>XB(I+1)) THEN
            DUMMY   = XB(I)
            XB(I)   = XB(I+1)
            XB(I+1) = DUMMY
         ENDIF
      ENDDO
      
      VT=>VENTS(N)
      
      IF (XB(1)==XB(2)) VT%TOTAL_INPUT_AREA = (XB(4)-XB(3))*(XB(6)-XB(5))
      IF (XB(3)==XB(4)) VT%TOTAL_INPUT_AREA = (XB(2)-XB(1))*(XB(6)-XB(5))
      IF (XB(5)==XB(6)) VT%TOTAL_INPUT_AREA = (XB(2)-XB(1))*(XB(4)-XB(3))

      XB(1) = MAX(XB(1),XS)
      XB(2) = MIN(XB(2),XF)
      XB(3) = MAX(XB(3),YS)
      XB(4) = MIN(XB(4),YF)
      XB(5) = MAX(XB(5),ZS)
      XB(6) = MIN(XB(6),ZF)
 
      IF (XB(1)>XF .OR. XB(2)<XS .OR. XB(3)>YF .OR. XB(4)<YS .OR. XB(5)>ZF .OR. XB(6)<ZS) REJECT_VENT = .TRUE.
 
      VT%I1 = NINT( GINV(XB(1)-XS,1,NM)*RDXI   ) 
      VT%I2 = NINT( GINV(XB(2)-XS,1,NM)*RDXI   )
      VT%J1 = NINT( GINV(XB(3)-YS,2,NM)*RDETA  ) 
      VT%J2 = NINT( GINV(XB(4)-YS,2,NM)*RDETA  )
      VT%K1 = NINT( GINV(XB(5)-ZS,3,NM)*RDZETA )
      VT%K2 = NINT( GINV(XB(6)-ZS,3,NM)*RDZETA )
 
      IF (XB(1)==XB(2)) THEN
         IF (VT%J1==VT%J2 .OR. VT%K1==VT%K2) REJECT_VENT=.TRUE.
      ENDIF
      IF (XB(3)==XB(4)) THEN
         IF (VT%I1==VT%I2 .OR. VT%K1==VT%K2) REJECT_VENT=.TRUE.
      ENDIF
      IF (XB(5)==XB(6)) THEN
         IF (VT%I1==VT%I2 .OR. VT%J1==VT%J2) REJECT_VENT=.TRUE.
      ENDIF
 
      ! Evacuation criteria
 
      IF (.NOT.EVACUATION .AND. EVACUATION_ONLY(NM)) REJECT_VENT=.TRUE.
      IF (EVACUATION .AND. .NOT.EVACUATION_ONLY(NM)) REJECT_VENT=.TRUE.
 
      ! If the VENT is to rejected
 
      IF (REJECT_VENT) THEN
         N = N-1
         N_VENT = N_VENT-1
         CYCLE READ_VENT_LOOP
      ENDIF
 
      ! Vent area
 
      VT%X1 = XB(1)
      VT%X2 = XB(2)
      VT%Y1 = XB(3)
      VT%Y2 = XB(4)
      VT%Z1 = XB(5)
      VT%Z2 = XB(6)
 
      IF (XB(1)==XB(2)) VT%INPUT_AREA = (XB(4)-XB(3))*(XB(6)-XB(5))
      IF (XB(3)==XB(4)) VT%INPUT_AREA = (XB(2)-XB(1))*(XB(6)-XB(5))
      IF (XB(5)==XB(6)) VT%INPUT_AREA = (XB(2)-XB(1))*(XB(4)-XB(3))
 
      ! Check the SURF_ID against the list of SURF's

      CALL CHECK_SURF_NAME(SURF_ID,EX)
      IF (.NOT.EX) THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: SURF_ID ',TRIM(SURF_ID),' not found'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      ! Assign IBC, Index of the Boundary Condition

      VT%IBC = DEFAULT_SURF_INDEX
      DO NNN=0,N_SURF
         IF (SURF_ID==SURFACE(NNN)%ID) VT%IBC = NNN
      ENDDO

      IF (SURF_ID=='OPEN')                            VT%TYPE_INDICATOR =  2
      IF (SURF_ID=='MIRROR' .OR. SURF_ID=='PERIODIC') VT%TYPE_INDICATOR = -2
      IF ((MB/='null' .OR.  PBX>-1.E5_EB .OR. PBY>-1.E5_EB .OR. PBZ>-1.E5_EB) .AND. SURF_ID=='OPEN') VT%TYPE_INDICATOR = -2
 
      VT%BOUNDARY_TYPE = SOLID_BOUNDARY
      IF (VT%IBC==OPEN_SURF_INDEX)     VT%BOUNDARY_TYPE = OPEN_BOUNDARY
      IF (VT%IBC==MIRROR_SURF_INDEX)   VT%BOUNDARY_TYPE = MIRROR_BOUNDARY
      IF (VT%IBC==PERIODIC_SURF_INDEX) VT%BOUNDARY_TYPE = PERIODIC_BOUNDARY
      VT%IOR = IOR
 
      VT%ORDINAL = NN
 
      ! Activate and Deactivate logic

      VT%ACTIVATED = .TRUE.
      VT%DEVC_ID   = DEVC_ID
      VT%CTRL_ID   = CTRL_ID
      CALL SEARCH_CONTROLLER('VENT',CTRL_ID,DEVC_ID,VT%DEVC_INDEX,VT%CTRL_INDEX,N)
      IF (DEVC_ID /= 'null') THEN
         IF (.NOT.DEVICE(VT%DEVC_INDEX)%INITIAL_STATE) VT%ACTIVATED = .FALSE.
      ENDIF
      IF (CTRL_ID /= 'null') THEN
         IF (.NOT.CONTROL(VT%CTRL_INDEX)%INITIAL_STATE) VT%ACTIVATED = .FALSE.
      ENDIF

      IF ( (VT%BOUNDARY_TYPE==OPEN_BOUNDARY .OR. VT%BOUNDARY_TYPE==MIRROR_BOUNDARY .OR. VT%BOUNDARY_TYPE==PERIODIC_BOUNDARY) .AND. &
           (VT%DEVC_ID /= 'null' .OR. VT%CTRL_ID /= 'null') ) THEN
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: OPEN, MIRROR, OR PERIODIC VENT',NN,' cannot be controlled by a device'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      ! Set the VENT color index

      SELECT CASE(COLOR)
         CASE('INVISIBLE')
            VT%COLOR_INDICATOR = 8
            TRANSPARENCY = 0._EB
         CASE('null')
            VT%COLOR_INDICATOR = 99
         CASE DEFAULT
            VT%COLOR_INDICATOR = 99
            CALL COLOR2RGB(RGB,COLOR)
      END SELECT
      IF (VT%COLOR_INDICATOR==8) VT%TYPE_INDICATOR = -2
      IF (OUTLINE)               VT%TYPE_INDICATOR =  2
      VT%RGB = RGB
      VT%TRANSPARENCY = TRANSPARENCY 

      ! Parameters for specified spread of a fire over a VENT
 
      VT%X0 = XYZ(1)
      VT%Y0 = XYZ(2)
      VT%Z0 = XYZ(3)
      VT%FIRE_SPREAD_RATE = SPREAD_RATE / TIME_SHRINK_FACTOR

      ! Dynamic Pressure

      VT%DYNAMIC_PRESSURE = DYNAMIC_PRESSURE
      IF (PRESSURE_RAMP/='null') CALL GET_RAMP_INDEX(PRESSURE_RAMP,'TIME',VT%PRESSURE_RAMP_INDEX)

      ! Miscellaneous
 
      VT%TMP_EXTERIOR = TMP_EXTERIOR + TMPM
      IF (VT%TMP_EXTERIOR>0._EB) TMPMIN = MIN(TMPMIN,VT%TMP_EXTERIOR) 
      VT%MASS_FRACTION = MASS_FRACTION

      VT%TEXTURE(:) = TEXTURE_ORIGIN(:)
      
38 CONTINUE
   ENDDO READ_VENT_LOOP
37 REWIND(LU_INPUT)

ENDDO MESH_LOOP_1

! Go through all the meshes again, but this time only if PROCESS(NM)==MYID

MESH_LOOP_2: DO NM=1,NMESHES

   IF (PROCESS(NM)/=MYID) CYCLE MESH_LOOP_2

   M=>MESHES(NM)
   CALL POINT_TO_MESH(NM)
 
   ! Check vents and assign orientations
 
   VENTLOOP2: DO N=1,N_VENT
 
      VT => VENTS(N)
 
      I1 = VT%I1
      I2 = VT%I2
      J1 = VT%J1
      J2 = VT%J2
      K1 = VT%K1
      K2 = VT%K2
 
      IF (VT%IOR==0) THEN
         IF (I1==      0 .AND. I2==0) VT%IOR =  1
         IF (I1==IBAR .AND. I2==IBAR) VT%IOR = -1
         IF (J1==      0 .AND. J2==0) VT%IOR =  2
         IF (J1==JBAR .AND. J2==JBAR) VT%IOR = -2
         IF (K1==      0 .AND. K2==0) VT%IOR =  3
         IF (K1==KBAR .AND. K2==KBAR) VT%IOR = -3
      ENDIF
 
      ORIENTATION_IF: IF (VT%IOR==0) THEN
         IF (I1==I2) THEN
            DO K=K1+1,K2
               DO J=J1+1,J2
                  IF (.NOT.SOLID(CELL_INDEX(I2+1,J,K))) VT%IOR =  1
                  IF (.NOT.SOLID(CELL_INDEX(I2  ,J,K))) VT%IOR = -1
               ENDDO
            ENDDO
         ENDIF
         IF (J1==J2) THEN
            DO K=K1+1,K2
               DO I=I1+1,I2
                  IF (.NOT.SOLID(CELL_INDEX(I,J2+1,K))) VT%IOR =  2
                  IF (.NOT.SOLID(CELL_INDEX(I,J2  ,K))) VT%IOR = -2
               ENDDO
            ENDDO
         ENDIF
         IF (K1==K2) THEN
            DO J=J1+1,J2
               DO I=I1+1,I2
                  IF (.NOT.SOLID(CELL_INDEX(I,J,K2+1))) VT%IOR =  3
                  IF (.NOT.SOLID(CELL_INDEX(I,J,K2  ))) VT%IOR = -3
               ENDDO
            ENDDO
         ENDIF
      ENDIF ORIENTATION_IF
 
      IF (VT%IOR==0) THEN
         WRITE(MESSAGE,'(A,I3,A,I3)')  'ERROR: Specify orientation of VENT ',VT%ORDINAL, ', MESH NUMBER',NM
         CALL SHUTDOWN(MESSAGE)
      ENDIF
 
      ! Other error messages for VENTs
 
      SELECT CASE(ABS(VT%IOR))
         CASE(1)
            IF (I1>=1 .AND. I1<=IBM1) THEN
               IF (VT%BOUNDARY_TYPE==OPEN_BOUNDARY.OR.VT%BOUNDARY_TYPE==MIRROR_BOUNDARY.OR.VT%BOUNDARY_TYPE==PERIODIC_BOUNDARY) THEN
                  WRITE(MESSAGE,'(A,I3,A)')  'ERROR: OPEN, MIRROR, OR PERIODIC VENT ',VT%ORDINAL, ' must be an exterior boundary.'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
               VT%BOUNDARY_TYPE = SOLID_BOUNDARY
               IF (.NOT.SOLID(CELL_INDEX(I2+1,J2,K2)) .AND.  .NOT.SOLID(CELL_INDEX(I2,J2,K2))) THEN
                  WRITE(MESSAGE,'(A,I3,A)')  'ERROR: VENT ',VT%ORDINAL, ' must be attached to a solid obstruction'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
            ENDIF
         CASE(2)
            IF (J1>=1 .AND. J1<=JBM1) THEN
               IF (VT%BOUNDARY_TYPE==OPEN_BOUNDARY.OR.VT%BOUNDARY_TYPE==MIRROR_BOUNDARY.OR.VT%BOUNDARY_TYPE==PERIODIC_BOUNDARY) THEN
                  WRITE(MESSAGE,'(A,I3,A)')  'ERROR: OPEN, MIRROR, OR PERIODIC VENT ',VT%ORDINAL, ' must be an exterior boundary.'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
               VT%BOUNDARY_TYPE = SOLID_BOUNDARY
               IF (.NOT.SOLID(CELL_INDEX(I2,J2+1,K2)) .AND.  .NOT.SOLID(CELL_INDEX(I2,J2,K2))) THEN
                  WRITE(MESSAGE,'(A,I3,A)')  'ERROR: VENT ',VT%ORDINAL, ' must be attached to a solid obstruction'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
            ENDIF
         CASE(3)
            IF (K1>=1 .AND. K1<=KBM1) THEN
               IF (VT%BOUNDARY_TYPE==OPEN_BOUNDARY.OR.VT%BOUNDARY_TYPE==MIRROR_BOUNDARY.OR.VT%BOUNDARY_TYPE==PERIODIC_BOUNDARY) THEN
                  WRITE(MESSAGE,'(A,I3,A)')  'ERROR: OPEN, MIRROR, OR PERIODIC VENT ',VT%ORDINAL, ' must be an exterior boundary.'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
               VT%BOUNDARY_TYPE = SOLID_BOUNDARY
               IF (.NOT.SOLID(CELL_INDEX(I2,J2,K2+1)) .AND. .NOT.SOLID(CELL_INDEX(I2,J2,K2))) THEN
                  WRITE(MESSAGE,'(A,I3,A)')  'ERROR: VENT ',VT%ORDINAL, ' must be attached to a solid obstruction'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
            ENDIF
      END SELECT
 
      ! Open up boundary cells if it is an open vent
 
      IF ( VT%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
         SELECT CASE(VT%IOR)
            CASE( 1) 
               CALL BLOCK_CELL(NM,   0,   0,J1+1,  J2,K1+1,  K2,0,0)
            CASE(-1) 
               CALL BLOCK_CELL(NM,IBP1,IBP1,J1+1,  J2,K1+1,  K2,0,0)
            CASE( 2) 
               CALL BLOCK_CELL(NM,I1+1,  I2,   0,   0,K1+1,  K2,0,0)
            CASE(-2) 
               CALL BLOCK_CELL(NM,I1+1,  I2,JBP1,JBP1,K1+1,  K2,0,0)
            CASE( 3) 
               CALL BLOCK_CELL(NM,I1+1,  I2,J1+1,  J2,   0,   0,0,0)
            CASE(-3) 
               CALL BLOCK_CELL(NM,I1+1,  I2,J1+1,  J2,KBP1,KBP1,0,0)
         END SELECT
      ENDIF
 
   ENDDO VENTLOOP2
 
   ! Compute vent areas and check for passive openings
 
   VENT_LOOP_3: DO N=1,N_VENT
 
      VT => VENTS(N)
 
      VT%FDS_AREA = 0._EB
      I1 = VT%I1
      I2 = VT%I2
      J1 = VT%J1
      J2 = VT%J2
      K1 = VT%K1
      K2 = VT%K2
 
      SELECT CASE(ABS(VT%IOR))
         CASE(1)
            DO K=K1+1,K2
               DO J=J1+1,J2
                  VT%FDS_AREA = VT%FDS_AREA + DY(J)*DZ(K)
               ENDDO
            ENDDO
         CASE(2)
            DO K=K1+1,K2
               DO I=I1+1,I2
                  VT%FDS_AREA = VT%FDS_AREA + DX(I)*DZ(K)
               ENDDO
            ENDDO
         CASE(3)
            DO J=J1+1,J2
               DO I=I1+1,I2
                  VT%FDS_AREA = VT%FDS_AREA + DX(I)*DY(J)
               ENDDO
            ENDDO
      END SELECT
 
   ENDDO  VENT_LOOP_3
 
ENDDO MESH_LOOP_2
 
END SUBROUTINE READ_VENT
 
 
SUBROUTINE READ_INIT
USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT
REAL(EB) :: TEMPERATURE,DENSITY,MASS_FRACTION(MAX_SPECIES),RR_SUM,YY_GET(1:N_SPECIES)
INTEGER  :: N,NN
TYPE(INITIALIZATION_TYPE), POINTER :: IN
NAMELIST /INIT/ XB,TEMPERATURE,DENSITY,MASS_FRACTION
 
N_INIT = 0
REWIND(LU_INPUT)
COUNT_LOOP: DO
   CALL CHECKREAD('INIT',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_LOOP
   READ(LU_INPUT,NML=INIT,END=11,ERR=12,IOSTAT=IOS)
   N_INIT = N_INIT + 1
   12 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with INIT no.',N_INIT+1
      CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_LOOP
11 REWIND(LU_INPUT)
 
! If there are no INIT lines, return

IF (N_INIT==0) RETURN
 
ALLOCATE(INITIALIZATION(N_INIT),STAT=IZERO)
CALL ChkMemErr('READ','INITIALIZATION',IZERO)
 
INIT_LOOP: DO N=1,N_INIT
   IN => INITIALIZATION(N)
   XB(1)         = -1000000._EB
   XB(2)         =  1000000._EB
   XB(3)         = -1000000._EB
   XB(4)         =  1000000._EB
   XB(5)         = -1000000._EB
   XB(6)         =  1000000._EB
   TEMPERATURE   = -1000._EB
   DENSITY       = -1000._EB
   MASS_FRACTION =  0._EB
   DO NN=1,N_SPECIES
      MASS_FRACTION(NN) = SPECIES(NN)%YY0
   ENDDO
 
   CALL CHECKREAD('INIT',LU_INPUT,IOS)
   IF (IOS==1) EXIT INIT_LOOP
   READ(LU_INPUT,INIT) 
 
   IN%X1 = XB(1)
   IN%X2 = XB(2)
   IN%Y1 = XB(3)
   IN%Y2 = XB(4)
   IN%Z1 = XB(5)
   IN%Z2 = XB(6)
   IN%TEMPERATURE   = TEMPERATURE + TMPM
   IN%DENSITY       = DENSITY
   IN%MASS_FRACTION = MASS_FRACTION

   IF (DENSITY     > 0._EB) RHOMAX = MAX(RHOMAX,IN%DENSITY)
   IF (TEMPERATURE > 0._EB) TMPMIN = MIN(TMPMIN,IN%TEMPERATURE)
   YY_GET(:) = IN%MASS_FRACTION
   CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RR_SUM)

   IF (IN%TEMPERATURE > 0._EB .AND. IN%DENSITY < 0._EB) THEN
      IN%DENSITY        = P_INF/(IN%TEMPERATURE*RR_SUM)
      IN%ADJUST_DENSITY = .TRUE.
   ENDIF
   IF (IN%TEMPERATURE < 0._EB .AND. IN%DENSITY > 0._EB) THEN
      IN%TEMPERATURE        = P_INF/(IN%DENSITY*RR_SUM)
      IN%ADJUST_TEMPERATURE = .TRUE.
   ENDIF
   IF (IN%TEMPERATURE < 0._EB .AND. IN%DENSITY < 0._EB) THEN
      IN%TEMPERATURE    = TMPA
      IN%DENSITY        = P_INF/(IN%TEMPERATURE*RR_SUM)
      IN%ADJUST_DENSITY = .TRUE.
   ENDIF

ENDDO INIT_LOOP
REWIND(LU_INPUT)

END SUBROUTINE READ_INIT


SUBROUTINE READ_ZONE
 
REAL(EB) :: LEAK_AREA(0:MAX_LEAK_PATHS)
INTEGER  :: N,NM
LOGICAL :: SEALED,READ_ZONE_LINES
CHARACTER(30) :: ID
NAMELIST /ZONE/ XB,LEAK_AREA,ID
 
N_ZONE = 0
REWIND(LU_INPUT)
COUNT_ZONE_LOOP: DO
   CALL CHECKREAD('ZONE',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_ZONE_LOOP
   READ(LU_INPUT,NML=ZONE,END=11,ERR=12,IOSTAT=IOS)
   N_ZONE = N_ZONE + 1
   12 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with ZONE no.',N_ZONE+1
      CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_ZONE_LOOP
11 REWIND(LU_INPUT)
 
! Check to see if there are any OPEN vents. If there are not, and there are no declared pressure ZONEs, stop with an error.

SEALED = .TRUE.
IF (ALL(EVACUATION_ONLY)) SEALED = .FALSE.    

DO NM=1,NMESHES
   IF (.NOT.EVACUATION_ONLY(NM)) THEN      
      M => MESHES(NM)
      DO N=1,M%N_VENT
         VT => M%VENTS(N)
         IF (VT%BOUNDARY_TYPE==OPEN_BOUNDARY)     SEALED = .FALSE.
         IF (VT%BOUNDARY_TYPE==PERIODIC_BOUNDARY) SEALED = .FALSE.
      ENDDO
   END IF
ENDDO

! If the whole domain lacks on OPEN or PERIODIC boundary, assume it to be one big pressure zone

READ_ZONE_LINES = .TRUE.
IF (SEALED .AND. N_ZONE==0) THEN
   N_ZONE = 1
   READ_ZONE_LINES = .FALSE.
ENDIF

! Stop the calculation if the user has declared ISOTHERMAL=.TRUE.

IF (ISOTHERMAL .AND. N_ZONE>0) THEN
   WRITE(MESSAGE,'(A)')  ' ERROR: ISOTHERMAL=.TRUE. must not be declared if there are sealed compartments'
   CALL SHUTDOWN(MESSAGE)
ENDIF

! Make sure that there are no leak paths to undefined pressure ZONEs

DO N=0,N_SURF
   SF => SURFACE(N)
   IF (SF%LEAK_PATH(1)>N_ZONE .OR. SF%LEAK_PATH(2)>N_ZONE) SF%LEAK_PATH = -1
ENDDO

! Allocate array to indicate if pressure ZONEs are connected

ALLOCATE(CONNECTED_ZONES(0:N_ZONE,0:N_ZONE,NMESHES),STAT=IZERO)
CALL ChkMemErr('READ','CONNECTED_ZONES',IZERO)
CONNECTED_ZONES = .FALSE.

! If there are no ZONE lines, return

IF (N_ZONE==0) RETURN
 
! Allocate ZONE arrays

ALLOCATE(P_ZONE(N_ZONE),STAT=IZERO)
CALL ChkMemErr('READ','P_ZONE',IZERO)

! Read in and process ZONE lines

READ_ZONE_LOOP: DO N=1,N_ZONE
 
   WRITE(ID,'(A,I2.2)') 'ZONE_',N
   LEAK_AREA     = 0._EB
   XB(1)         = -1000000._EB
   XB(2)         =  1000000._EB
   XB(3)         = -1000000._EB
   XB(4)         =  1000000._EB
   XB(5)         = -1000000._EB
   XB(6)         =  1000000._EB
 
   IF (READ_ZONE_LINES) THEN
      CALL CHECKREAD('ZONE',LU_INPUT,IOS)
      IF (IOS==1) EXIT READ_ZONE_LOOP
      READ(LU_INPUT,ZONE) 
   ENDIF
 
   P_ZONE(N)%ID = ID
   P_ZONE(N)%LEAK_AREA = LEAK_AREA
   P_ZONE(N)%X1 = XB(1)
   P_ZONE(N)%X2 = XB(2)
   P_ZONE(N)%Y1 = XB(3)
   P_ZONE(N)%Y2 = XB(4)
   P_ZONE(N)%Z1 = XB(5)
   P_ZONE(N)%Z2 = XB(6)
 
ENDDO READ_ZONE_LOOP
REWIND(LU_INPUT)

END SUBROUTINE READ_ZONE

 
SUBROUTINE READ_DEVC

! Just read in the DEViCes and the store the info in DEVICE()

USE DEVICE_VARIABLES, ONLY: DEVICE_TYPE, DEVICE, N_DEVC
INTEGER  :: N,NN,NM,MESH_NUMBER,N_DEVCO,IOR,TRIP_DIRECTION
REAL(EB) :: DEPTH,ORIENTATION(3),ROTATION,SETPOINT,FLOWRATE,BYPASS_FLOWRATE,DELAY,XYZ(3),CONVERSION_FACTOR
CHARACTER(30) :: QUANTITY,PROP_ID,CTRL_ID,DEVC_ID,SURF_ID,STATISTICS,PART_ID,MATL_ID,SPEC_ID,UNITS
LOGICAL :: INITIAL_STATE,LATCH,DRY,TIME_AVERAGED
TYPE (DEVICE_TYPE), POINTER :: DV
NAMELIST /DEVC/ DEPTH,FYI,IOR,ID,ORIENTATION,PROP_ID,QUANTITY,ROTATION,XB,XYZ,INITIAL_STATE,LATCH,TRIP_DIRECTION,CTRL_ID,& 
                SETPOINT,DEVC_ID,FLOWRATE,DELAY,BYPASS_FLOWRATE,SURF_ID,STATISTICS,PART_ID,MATL_ID,SPEC_ID,DRY,CONVERSION_FACTOR,&
                UNITS,TIME_AVERAGED

N_DEVC = 0
REWIND(LU_INPUT)
COUNT_DEVC_LOOP: DO
   CALL CHECKREAD('DEVC',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_DEVC_LOOP
   READ(LU_INPUT,NML=DEVC,END=11,ERR=12,IOSTAT=IOS)
   N_DEVC = N_DEVC + 1
   12 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I4)') 'ERROR: Problem with DEVC no.',N_DEVC+1
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO COUNT_DEVC_LOOP
11 REWIND(LU_INPUT)
 
IF (N_DEVC==0) RETURN

! Allocate DEVICE array and set initial values of all to 0

ALLOCATE(DEVICE(N_DEVC),STAT=IZERO)
CALL ChkMemErr('READ','DEVICE',IZERO)
 
! Read in the DEVC lines

N_DEVCO = N_DEVC
N       = 0

READ_DEVC_LOOP: DO NN=1,N_DEVCO

   N          = N+1
   CALL CHECKREAD('DEVC',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_DEVC_LOOP
   CALL SET_DEVC_DEFAULTS
   READ(LU_INPUT,DEVC) 


   IF (XB(1)>-1.E5_EB) THEN
      XYZ(1) = 0.5_EB*(XB(1)+XB(2))
      XYZ(2) = 0.5_EB*(XB(3)+XB(4))
      XYZ(3) = 0.5_EB*(XB(5)+XB(6))
   ELSE
      IF (XYZ(1) < -1.E5_EB) THEN
         WRITE(MESSAGE,'(A,I5,A)')  ' ERROR: DEVC ',NN,' must have coordinates, even if it is not a point quantity'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDIF

   ! Determine which mesh the device is in

   BAD = .FALSE.
   MESH_LOOP: DO NM=1,NMESHES
      IF (.NOT.EVACUATION_ONLY(NM)) THEN      
         M=>MESHES(NM)
         IF (XYZ(1)>=M%XS .AND. XYZ(1)<=M%XF .AND. XYZ(2)>=M%YS .AND. XYZ(2)<=M%YF .AND. XYZ(3)>=M%ZS .AND. XYZ(3)<=M%ZF) THEN
            MESH_NUMBER = NM
            EXIT MESH_LOOP
         ENDIF
      ENDIF
      IF (NM==NMESHES) BAD = .TRUE.
   ENDDO MESH_LOOP

   ! Make sure there is either a QUANTITY or PROP_ID for the DEVICE

   IF (QUANTITY=='null' .AND. PROP_ID=='null') THEN
      WRITE(MESSAGE,'(A,I5,A)')  ' ERROR: DEVC ',NN,' must have either an output QUANTITY or PROP_ID'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   IF (BAD) THEN
      IF (QUANTITY=='TIME') THEN
         XYZ(1) = MESHES(1)%XS
         XYZ(2) = MESHES(1)%YS
         XYZ(3) = MESHES(1)%ZS
         MESH_NUMBER = 1
      ELSE
         N      = N-1
         N_DEVC = N_DEVC-1
         WRITE(MESSAGE,'(A,I3,A)') 'WARNING: DEVC ',NN,', is not within any mesh.'  
         IF (MYID==0) WRITE(LU_ERR,'(A)') TRIM(MESSAGE)
         CYCLE READ_DEVC_LOOP
      ENDIF
   ENDIF

   ! Assign properties to the DEVICE array

   DV => DEVICE(N)
   M  => MESHES(MESH_NUMBER)

   DV%CONVERSION_FACTOR= CONVERSION_FACTOR
   DV%DEPTH            = DEPTH
   DV%IOR              = IOR
   DV%ID               = ID
   DV%MESH             = MESH_NUMBER
   DV%ORDINAL          = NN
   DV%ORIENTATION(1:3) = ORIENTATION(1:3)/SQRT(ORIENTATION(1)**2+ORIENTATION(2)**2+ORIENTATION(3)**2)
   DV%PROP_ID          = PROP_ID
   DV%CTRL_ID          = CTRL_ID   
   DV%DEVC_ID          = DEVC_ID   
   DV%SURF_ID          = SURF_ID            
   DV%PART_ID          = PART_ID            
   DV%MATL_ID          = MATL_ID            
   DV%SPEC_ID          = SPEC_ID            
   DV%QUANTITY         = QUANTITY
   DV%ROTATION         = ROTATION*TWOPI/360._EB
   DV%SETPOINT         = SETPOINT
   DV%LATCH            = LATCH
   DV%TRIP_DIRECTION   = TRIP_DIRECTION
   DV%INITIAL_STATE    = INITIAL_STATE
   DV%CURRENT_STATE    = INITIAL_STATE
   DV%PRIOR_STATE      = INITIAL_STATE
   DV%FLOWRATE         = FLOWRATE
   DV%BYPASS_FLOWRATE  = BYPASS_FLOWRATE
   DV%STATISTICS       = STATISTICS   
   DV%TIME_AVERAGED    = TIME_AVERAGED 
   DV%SURF_INDEX       = 0
   DV%UNITS            = UNITS
   DV%DELAY            = DELAY / TIME_SHRINK_FACTOR
   DV%X1               = XB(1)
   DV%X2               = XB(2)
   DV%Y1               = XB(3)
   DV%Y2               = XB(4)
   DV%Z1               = XB(5)
   DV%Z2               = XB(6)
   DV%X                = XYZ(1)
   DV%Y                = XYZ(2)
   DV%Z                = XYZ(3)
   DV%DRY              = DRY

   ! Coordinates for non-point devices

   IF (XB(1)>-1.E5_EB .OR. STATISTICS/='null') THEN
      NM = DV%MESH
      M=>MESHES(NM)
      XB(1) = MAX(XB(1),M%XS)
      XB(2) = MIN(XB(2),M%XF)
      XB(3) = MAX(XB(3),M%YS)
      XB(4) = MIN(XB(4),M%YF)
      XB(5) = MAX(XB(5),M%ZS)
      XB(6) = MIN(XB(6),M%ZF)
      DV%X1 = XB(1)
      DV%X2 = XB(2)
      DV%Y1 = XB(3)
      DV%Y2 = XB(4)
      DV%Z1 = XB(5)
      DV%Z2 = XB(6)
      DV%I1 = NINT( GINV(XB(1)-M%XS,1,NM)*M%RDXI)
      DV%I2 = NINT( GINV(XB(2)-M%XS,1,NM)*M%RDXI)
      DV%J1 = NINT( GINV(XB(3)-M%YS,2,NM)*M%RDETA)
      DV%J2 = NINT( GINV(XB(4)-M%YS,2,NM)*M%RDETA)
      DV%K1 = NINT( GINV(XB(5)-M%ZS,3,NM)*M%RDZETA)
      DV%K2 = NINT( GINV(XB(6)-M%ZS,3,NM)*M%RDZETA)
      IF (DV%I1<DV%I2) DV%I1 = DV%I1 + 1
      IF (DV%J1<DV%J2) DV%J1 = DV%J1 + 1
      IF (DV%K1<DV%K2) DV%K1 = DV%K1 + 1
      IF (XB(1)==XB(2) .AND. IOR==0) DV%IOR = 1
      IF (XB(3)==XB(4) .AND. IOR==0) DV%IOR = 2
      IF (XB(5)==XB(6) .AND. IOR==0) DV%IOR = 3
   ENDIF
   
ENDDO READ_DEVC_LOOP
REWIND(LU_INPUT)

CONTAINS

SUBROUTINE SET_DEVC_DEFAULTS

CONVERSION_FACTOR = 1._EB
DEPTH            = 0._EB
IOR              = 0
SELECT CASE(N)
   CASE(1:9)
      WRITE(ID,'(A7,I1)') 'Device_',N
   CASE(10:99)
      WRITE(ID,'(A7,I2)') 'Device_',N
   CASE(100:999)
      WRITE(ID,'(A7,I3)') 'Device_',N
   CASE(1000:9999)
      WRITE(ID,'(A7,I4)') 'Device_',N
END SELECT
ORIENTATION(1:3) = (/0._EB,0._EB,-1._EB/)
PROP_ID          = 'null'
CTRL_ID          = 'null'
DEVC_ID          = 'null'
SURF_ID          = 'null'
PART_ID          = 'null'
MATL_ID          = 'null'
SPEC_ID          = 'null'
FLOWRATE         = 0._EB
DELAY            = 0._EB
BYPASS_FLOWRATE  = 0._EB
QUANTITY         = 'null'
ROTATION         = 0._EB
XB(1)            = -1.E6_EB
XB(2)            =  1.E6_EB
XB(3)            = -1.E6_EB
XB(4)            =  1.E6_EB
XB(5)            = -1.E6_EB
XB(6)            =  1.E6_EB
INITIAL_STATE    = .FALSE.
LATCH            = .TRUE.
SETPOINT         = 1.E20_EB
STATISTICS       = 'null'
TRIP_DIRECTION   = 1
TIME_AVERAGED    = .TRUE.
UNITS            = 'null'
XYZ              = -1.E6_EB
DRY              = .FALSE.

END SUBROUTINE SET_DEVC_DEFAULTS

END SUBROUTINE READ_DEVC



SUBROUTINE READ_CTRL

! Just read in the ConTRoL parameters and store in the array CONTROL

USE CONTROL_VARIABLES
USE MATH_FUNCTIONS, ONLY : GET_RAMP_INDEX

LOGICAL :: INITIAL_STATE, LATCH
INTEGER :: CYCLES,N,NC
REAL(EB) :: SETPOINT(2), DELAY, CYCLE_TIME
CHARACTER(30) :: ID,FUNCTION_TYPE,INPUT_ID(40),RAMP_ID,ON_BOUND
TYPE (CONTROL_TYPE), POINTER :: CF
NAMELIST /CTRL/ ID,LATCH,INITIAL_STATE,FUNCTION_TYPE,SETPOINT,DELAY,CYCLE_TIME,INPUT_ID,RAMP_ID,CYCLES,N,ON_BOUND
 
N_CTRL = 0
REWIND(LU_INPUT)
COUNT_CTRL_LOOP: DO
   CALL CHECKREAD('CTRL',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_CTRL_LOOP
   READ(LU_INPUT,NML=CTRL,END=11,ERR=12,IOSTAT=IOS)
   N_CTRL = N_CTRL + 1
   12 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I4)') 'ERROR: Problem with CTRL no.',N_CTRL+1
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO COUNT_CTRL_LOOP
11 REWIND(LU_INPUT)
 
IF (N_CTRL==0) RETURN

! Allocate CONTROL array and set initial values of all to 0

ALLOCATE(CONTROL(N_CTRL),STAT=IZERO)
CALL ChkMemErr('READ','CONTROL',IZERO)
 
! Read in the CTRL lines

READ_CTRL_LOOP: DO NC=1,N_CTRL

   CALL CHECKREAD('CTRL',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_CTRL_LOOP
   CALL SET_CTRL_DEFAULTS
   READ(LU_INPUT,CTRL) 

   ! Make sure there is either a FUNCTION_TYPE type for the CTRL

   IF (FUNCTION_TYPE=='null') THEN
      WRITE(MESSAGE,'(A,I5,A)')  ' ERROR: CTRL ',NC,' must have a FUNCTION_TYPE'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   ! Assign properties to the CONTROL array

   CF => CONTROL(NC)
   CF%ID              = ID
   CF%LATCH           = LATCH
   CF%INITIAL_STATE   = INITIAL_STATE
   CF%CURRENT_STATE   = INITIAL_STATE   
   CF%PRIOR_STATE     = INITIAL_STATE      
   CF%SETPOINT        = SETPOINT
   CF%DELAY           = DELAY / TIME_SHRINK_FACTOR
   CF%CYCLE_TIME      = CYCLE_TIME
   CF%CYCLES          = CYCLES
   CF%RAMP_ID         = RAMP_ID
   CF%N               = N   
   CF%INPUT_ID        = INPUT_ID
   IF (ON_BOUND=='UPPER') THEN
      CF%ON_BOUND = 1
   ELSE
      CF%ON_BOUND = -1
   ENDIF   
   !Assign control index 
   SELECT CASE(FUNCTION_TYPE)
      CASE('ALL')
         CF%CONTROL_INDEX = AND_GATE
      CASE('ANY')
         CF%CONTROL_INDEX = OR_GATE
      CASE('ONLY')
         CF%CONTROL_INDEX = XOR_GATE
      CASE('AT_LEAST')
         CF%CONTROL_INDEX = X_OF_N_GATE
      CASE('TIME_DELAY')
         CF%CONTROL_INDEX = TIME_DELAY
      CASE('DEADBAND')
         CF%CONTROL_INDEX = DEADBAND
      CASE('CYCLING')
         CF%CONTROL_INDEX = CYCLING
      CASE('CUSTOM')
         CF%CONTROL_INDEX = CUSTOM
         CALL GET_RAMP_INDEX(RAMP_ID,'CONTROL',CF%RAMP_INDEX)         
         CF%LATCH = .FALSE.
      CASE('KILL')
         CF%CONTROL_INDEX = KILL
      CASE('RESTART')
         CF%CONTROL_INDEX = CORE_DUMP
      CASE DEFAULT
         WRITE(MESSAGE,'(A,I5,A)')  ' ERROR: CTRL ',NC,' FUNCTION_TYPE not recognized'
         CALL SHUTDOWN(MESSAGE)
   END SELECT
   
ENDDO READ_CTRL_LOOP
REWIND(LU_INPUT)

CONTAINS

SUBROUTINE SET_CTRL_DEFAULTS
   ID            = 'null'
   LATCH         = .TRUE.
   INITIAL_STATE = .FALSE.
   SETPOINT      = 1000000._EB
   DELAY         = 0._EB
   CYCLE_TIME    = 1000000._EB
   CYCLES        = 1
   FUNCTION_TYPE = 'null'
   RAMP_ID       = 'null'
   INPUT_ID      = 'null'
   ON_BOUND      = 'LOWER'
   N             = 1
END SUBROUTINE SET_CTRL_DEFAULTS

END SUBROUTINE READ_CTRL



SUBROUTINE PROC_CTRL

! Process the CONTROL function parameters

USE DEVICE_VARIABLES, ONLY : DEVICE, N_DEVC, GAS_CELL_RAD_FLUX, GAS_CELL_RAD_DEVC_INDEX, N_GAS_CELL_RAD_DEVC
USE CONTROL_VARIABLES
INTEGER :: NC,NN,NNN
TYPE (CONTROL_TYPE), POINTER :: CF

PROC_CTRL_LOOP: DO NC = 1, N_CTRL

   CF => CONTROL(NC)

   ! setup input array

   CF%N_INPUTS = 0
   INPUT_COUNT: DO
      IF (CF%INPUT_ID(CF%N_INPUTS+1)=='null') EXIT INPUT_COUNT
      CF%N_INPUTS = CF%N_INPUTS + 1
   END DO INPUT_COUNT
   IF (CF%N_INPUTS==0) THEN
      WRITE(MESSAGE,'(A,I5,A)')  ' ERROR: CTRL ',NC,' must have at least one input'
      CALL SHUTDOWN(MESSAGE)
   ENDIF   
   
   ALLOCATE (CF%INPUT(CF%N_INPUTS),STAT=IZERO)
   CALL ChkMemErr('READ','CF%INPUT',IZERO)
   ALLOCATE (CF%INPUT_TYPE(CF%N_INPUTS),STAT=IZERO)
   CALL ChkMemErr('READ','CF%INPUT_TYPE',IZERO)
   
   BUILD_INPUT: DO NN = 1, CF%N_INPUTS
      CTRL_LOOP: DO NNN = 1, N_CTRL
         IF(CONTROL(NNN)%ID == CF%INPUT_ID(NN)) THEN
            CF%INPUT(NN) = NNN
            CF%INPUT_TYPE(NN) = CONTROL_INPUT
            IF (CF%CONTROL_INDEX == CUSTOM) THEN
               WRITE(MESSAGE,'(A,I5,A)')  ' ERROR: CUSTOM CTRL ',NC,' cannot have another CTRL as input'
               CALL SHUTDOWN(MESSAGE)
            ENDIF   
            CYCLE BUILD_INPUT
         ENDIF
      END DO CTRL_LOOP
      DEVC_LOOP: DO NNN = 1, N_DEVC
         IF(DEVICE(NNN)%ID == CF%INPUT_ID(NN)) THEN
            CF%INPUT(NN) = NNN
            CF%INPUT_TYPE(NN) = DEVICE_INPUT
            CYCLE BUILD_INPUT
         ENDIF
      END DO DEVC_LOOP
   WRITE(MESSAGE,'(A,I5,A,A)')  ' ERROR: CTRL ',NC,' cannot locate item for input ', TRIM(CF%INPUT_ID(NN))
   CALL SHUTDOWN(MESSAGE)
   END DO BUILD_INPUT

END DO PROC_CTRL_LOOP  
   
END SUBROUTINE PROC_CTRL   


SUBROUTINE PROC_DEVC

! Process the DEViCes

USE CONTROL_VARIABLES
USE DEVICE_VARIABLES, ONLY : DEVICE_TYPE, DEVICE, N_DEVC, PROPERTY, PROPERTY_TYPE, N_PROP, &
                             GAS_CELL_RAD_FLUX, GAS_CELL_RAD_DEVC_INDEX, N_GAS_CELL_RAD_DEVC
 
INTEGER  :: N,NN,NNN,NM,QUANTITY_INDEX,MAXCELLS,I,J,K
REAL(EB) :: XX,YY,ZZ,XX1,YY1,ZZ1,DISTANCE,SCANDISTANCE,DX,DY,DZ
TYPE (DEVICE_TYPE),  POINTER :: DV
TYPE (PROPERTY_TYPE),  POINTER :: PY
 
IF (N_DEVC==0) RETURN

! Set initial values for DEViCes

DEVICE(1:N_DEVC)%VALUE = 0._EB
DEVICE(1:N_DEVC)%COUNT = 0

PROC_DEVC_LOOP: DO N=1,N_DEVC

   DV => DEVICE(N)

   ! If the Device has a SURF_ID, get the SURF_INDEX

   IF (DV%SURF_ID/='null') THEN
      DO NN=1,N_SURF
         IF (SURFACE(NN)%ID==DV%SURF_ID) DV%SURF_INDEX = NN
      ENDDO
   ENDIF

   ! Check if the device PROPERTY exists and is appropriate

   DV%PROP_INDEX = 0
   IF (DV%PROP_ID /='null') THEN
      CALL GET_PROPERTY_INDEX(DV%PROP_INDEX,'DEVC',DV%PROP_ID)
      IF (DV%QUANTITY=='null' .AND. PROPERTY(DV%PROP_INDEX)%QUANTITY=='null') THEN
         WRITE(MESSAGE,'(5A)')  ' ERROR: DEVC ',TRIM(DV%ID),' or DEVC PROPerty ',TRIM(DV%PROP_ID),' must have a QUANTITY' 
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      IF (DV%QUANTITY=='null' .AND. PROPERTY(DV%PROP_INDEX)%QUANTITY/='null') DV%QUANTITY = PROPERTY(DV%PROP_INDEX)%QUANTITY
   ENDIF

   ! Check if the output QUANTITY exists and is appropriate

   IF (DV%QUANTITY /= 'null') THEN

      CALL GET_QUANTITY_INDEX(DV%SMOKEVIEW_LABEL,DV%SMOKEVIEW_BAR_LABEL,QUANTITY_INDEX,DV%SPEC_INDEX,DV%PART_INDEX, &
                              'DEVC',DV%QUANTITY,DV%SPEC_ID,DV%PART_ID)

      IF (OUTPUT_QUANTITY(QUANTITY_INDEX)%INTEGRATED .AND. DV%X1<=-1.E6_EB) THEN
         WRITE(MESSAGE,'(3A)')  ' ERROR: DEVC QUANTITY ',TRIM(DV%QUANTITY),' requires coordinates using XB'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (.NOT.OUTPUT_QUANTITY(QUANTITY_INDEX)%INTEGRATED .AND. DV%STATISTICS=='null' .AND. DV%X1>-1.E6_EB) THEN
         WRITE(MESSAGE,'(3A)')  ' ERROR: DEVC QUANTITY ',TRIM(DV%QUANTITY),' requires coordinates using XYZ, not XB'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (QUANTITY_INDEX<0 .AND. DV%IOR==0 .AND. DV%STATISTICS=='null') THEN
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: Specify orientation of DEVC ' ,N,' using the parameter IOR'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (QUANTITY_INDEX < 0 .AND. (DV%STATISTICS=='MASS MEAN' .OR. DV%STATISTICS=='VOLUME MEAN' .OR. &
                                    DV%STATISTICS=='VOLUME INTEGRAL' .OR. DV%STATISTICS=='MASS INTEGRAL') ) THEN
         WRITE(MESSAGE,'(A,I4)') 'ERROR: Invalid STATISTICS specified for wall DEVC ',N
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (QUANTITY_INDEX > 0 .AND. DV%STATISTICS=='SURFACE INTEGRAL') THEN
         WRITE(MESSAGE,'(A,I4)') 'ERROR: Invalid STATISTICS specified for gas DEVC ',N
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (QUANTITY_INDEX > 0 .AND. DV%STATISTICS/='null' .AND. DV%STATISTICS/='TIME INTEGRAL' .AND. DV%I1<0) THEN
         WRITE(MESSAGE,'(A,I4)') 'ERROR: XB required when geometrical STATISTICS specified for gas DEVC ',N
         CALL SHUTDOWN(MESSAGE)
      ENDIF

   ENDIF
   
  ! Assign properties to the DEVICE array

   M  => MESHES(DV%MESH)

   DV%T_CHANGE         = 10000000._EB
   DV%I                = MAX( 1 , MIN( M%IBAR , FLOOR(GINV(DV%X-M%XS,1,DV%MESH)*M%RDXI)  +1 ) )
   DV%J                = MAX( 1 , MIN( M%JBAR , FLOOR(GINV(DV%Y-M%YS,2,DV%MESH)*M%RDETA) +1 ) )
   DV%K                = MAX( 1 , MIN( M%KBAR , FLOOR(GINV(DV%Z-M%ZS,3,DV%MESH)*M%RDZETA)+1 ) )
   DV%OUTPUT_INDEX     = QUANTITY_INDEX
   IF (DV%UNITS=='null') DV%UNITS = OUTPUT_QUANTITY(DV%OUTPUT_INDEX)%UNITS
   DV%CTRL_INDEX       = 0
   DV%QUANTITY         = OUTPUT_QUANTITY(QUANTITY_INDEX)%NAME
   DV%T                = 0._EB
   DV%TMP_L            = TMPA
   DV%TI_VALUE         = 0._EB
   DV%TI_T             = 0._EB
   
   ! Do initialization of special models
   
   SPECIAL_QUANTITIES: SELECT CASE (DV%QUANTITY)

      CASE ('spot obscuration','CHAMBER OBSCURATION') 

         IF (DV%PROP_INDEX<1) THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DEVC ' ,N,' is a smoke detector and must have a PROP_ID'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         IF (PROPERTY(DV%PROP_INDEX)%SPEC_INDEX==0 .AND. .NOT.MIXTURE_FRACTION) THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DEVC ' ,N,' is a smoke detector and requires a fire or smoke source'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         ALLOCATE(DV%T_E(-1:1000))
         ALLOCATE(DV%Y_E(-1:1000))
         DV%T_E      = T_BEGIN - M%DT
         DV%Y_E      = 0._EB
         DV%N_T_E    = -1
         DV%Y_C      = 0._EB
         DV%SETPOINT = PROPERTY(DV%PROP_INDEX)%ACTIVATION_OBSCURATION
   
      CASE ('LINK TEMPERATURE','SPRINKLER LINK TEMPERATURE') 

         IF (DV%PROP_INDEX<1) THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DEVC ' ,N,' must have a PROP_ID'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         IF (PROPERTY(DV%PROP_INDEX)%ACTIVATION_TEMPERATURE <= -273.15_EB) THEN
            WRITE(MESSAGE,'(A,A)') 'ERROR: ACTIVATION_TEMPERATURE needed for PROP ',TRIM(DV%PROP_ID)
            CALL SHUTDOWN(MESSAGE)
         ENDIF

         DV%SETPOINT = PROPERTY(DV%PROP_INDEX)%ACTIVATION_TEMPERATURE
         DV%TMP_L    = PROPERTY(DV%PROP_INDEX)%INITIAL_TEMPERATURE

      CASE ('RADIATIVE_FLUX_GAS','RADIATIVE HEAT FLUX GAS')

         SELECT CASE(DV%IOR)
            CASE(-1)
               DV%ORIENTATION(1:3) = (/-1._EB, 0._EB, 0._EB/)
            CASE( 1)
               DV%ORIENTATION(1:3) = (/ 1._EB, 0._EB, 0._EB/)
            CASE(-2)
               DV%ORIENTATION(1:3) = (/ 0._EB,-1._EB, 0._EB/)
            CASE( 2)
               DV%ORIENTATION(1:3) = (/ 0._EB, 1._EB, 0._EB/)
            CASE(-3)
               DV%ORIENTATION(1:3) = (/ 0._EB, 0._EB,-1._EB/)
            CASE( 3)
               DV%ORIENTATION(1:3) = (/ 0._EB, 0._EB, 1._EB/)
         END SELECT
         GAS_CELL_RAD_FLUX    = .TRUE.
         N_GAS_CELL_RAD_DEVC  = N_GAS_CELL_RAD_DEVC + 1
         DV%GAS_CELL_RAD_FLUX = .TRUE.

      CASE ('SOLID DENSITY')

         IF (DV%MATL_ID=='null') THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DEVC ' ,N,' must have a MATL_ID'
            CALL SHUTDOWN(MESSAGE)
         ENDIF

      CASE ('CABLE TEMPERATURE')

         IF (DV%PROP_INDEX<1) THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DEVC ' ,N,' must have a PROP_ID'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         PY => PROPERTY(DV%PROP_INDEX)
         DV%SETPOINT = PY%CABLE_FAILURE_TEMPERATURE
         DV%TMP_L    = PY%INITIAL_TEMPERATURE
         GAS_CELL_RAD_FLUX   = .TRUE.
         N_GAS_CELL_RAD_DEVC = N_GAS_CELL_RAD_DEVC + 1
         DV%GAS_CELL_RAD_FLUX = .TRUE.
         IF (DV%DEPTH==0._EB) THEN
            IF (PY%CONDUIT_DIAMETER>0._EB) THEN
               DV%DEPTH = 0.5_EB*(PY%CONDUIT_DIAMETER-PY%CABLE_DIAMETER) + PY%CABLE_JACKET_THICKNESS
            ELSE
               DV%DEPTH = PY%CABLE_JACKET_THICKNESS
            ENDIF
         ENDIF
         DO NNN=1,N_SURF
            IF (SURFACE(NNN)%ID==PY%ID) DV%SURF_INDEX = NNN
         ENDDO

      CASE ('LAYER HEIGHT','UPPER TEMPERATURE','LOWER TEMPERATURE') 

         DV%K1 = MAX(1     ,DV%K1)
         DV%K2 = MIN(M%KBAR,DV%K2)

      CASE ('path obscuration','PATH OBSCURATION')

         IF (DV%PROP_INDEX>0) THEN
            IF (PROPERTY(DV%PROP_INDEX)%SPEC_INDEX==0 .AND. .NOT.MIXTURE_FRACTION) THEN
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DEVC ' ,N,' is a smoke detector and requires a fire or smoke source'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
         ELSE
            IF (.NOT.MIXTURE_FRACTION) THEN
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DEVC ' ,N,' is a smoke detector and requires a fire or smoke source'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
         ENDIF
         NM = DV%MESH
         M=>MESHES(NM)
         DISTANCE = SQRT((DV%X1-DV%X2)**2 + (DV%Y1-DV%Y2)**2 + (DV%Z1-DV%Z2)**2)
         SCANDISTANCE = 0.0001_EB * DISTANCE
         DX = (DV%X2-DV%X1) * 0.0001_EB
         DY = (DV%Y2-DV%Y1) * 0.0001_EB
         DZ = (DV%Z2-DV%Z1) * 0.0001_EB
         XX = DV%X1
         YY = DV%Y1
         ZZ = DV%Z1
         MAXCELLS = 2*MAX(M%IBAR,M%JBAR,M%KBAR)
         ALLOCATE(DV%I_PATH(MAXCELLS))
         ALLOCATE(DV%J_PATH(MAXCELLS))
         ALLOCATE(DV%K_PATH(MAXCELLS))
         ALLOCATE(DV%D_PATH(MAXCELLS))
         DV%D_PATH    = 0._EB
         DV%I_PATH = INT(GINV(DV%X1-M%XS,1,NM)*M%RDXI)   + 1
         DV%J_PATH = INT(GINV(DV%Y1-M%YS,2,NM)*M%RDETA)  + 1
         DV%K_PATH = INT(GINV(DV%Z1-M%ZS,3,NM)*M%RDZETA) + 1
         DV%N_PATH    = 1
         NN = 1
         DO NNN=1,10000
            XX = XX + DX
            I = INT(GINV(XX-M%XS,1,NM)*M%RDXI)   + 1
            YY = YY + DY
            J = INT(GINV(YY-M%YS,2,NM)*M%RDETA)  + 1
            ZZ = ZZ + DZ
            K = INT(GINV(ZZ-M%ZS,3,NM)*M%RDZETA) + 1
            IF (I==DV%I_PATH(NN) .AND. J==DV%J_PATH(NN) .AND. K==DV%K_PATH(NN)) THEN
               DV%D_PATH(NN) = DV%D_PATH(NN) + SCANDISTANCE
            ELSE
               NN = NN + 1
               DV%I_PATH(NN) = I
               DV%J_PATH(NN) = J
               DV%K_PATH(NN) = K
               XX1 = DX
               YY1 = DY
               ZZ1 = DZ
               IF (I/=DV%I_PATH(NN-1)) XX1 = XX-M%X(DV%I_PATH(NN-1))
               IF (J/=DV%J_PATH(NN-1)) YY1 = YY-M%Y(DV%J_PATH(NN-1))
               IF (K/=DV%K_PATH(NN-1)) ZZ1 = ZZ-M%Z(DV%K_PATH(NN-1))
               DV%D_PATH(NN)   = SCANDISTANCE - SQRT(XX1**2+YY1**2+ZZ1**2)
               DV%D_PATH(NN-1) = DV%D_PATH(NN-1) + SCANDISTANCE - DV%D_PATH(NN)
            ENDIF
         ENDDO
         DV%N_PATH = NN
                     
      CASE ('CONTROL')

         DO NN=1,N_CTRL
            IF (CONTROL(NN)%ID==DV%CTRL_ID) DV%CTRL_INDEX = NN
         ENDDO
         IF (DV%CTRL_ID/='null' .AND. DV%CTRL_INDEX<=0) THEN
            WRITE(MESSAGE,'(A,A,A)')  'ERROR: CONTROL ',TRIM(DV%CTRL_ID),' does not exist'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         DV%SETPOINT = 0.5
         DV%TRIP_DIRECTION = 1

      CASE ('aspiration','ASPIRATION')

         ! Check either for a specified SMOKE SPECies, or if mixture fraction model is being used
         IF (DV%PROP_INDEX>0) THEN
            IF (PROPERTY(DV%PROP_INDEX)%SPEC_INDEX==0 .AND. .NOT.MIXTURE_FRACTION) THEN
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: DEVC ' ,N,' is a smoke detector and requires a fire or smoke source'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
         ENDIF
         ! Count number of inputs for detector and verify that input is DENSITY with a specified SPEC_ID for smoke
         NNN = 0
         DO NN=1,N_DEVC
            IF (DEVICE(NN)%DEVC_ID==DV%ID) THEN
               IF (DEVICE(NN)%QUANTITY/='DENSITY' .OR. DEVICE(NN)%SPEC_ID=='null') THEN
                  WRITE(MESSAGE,'(A,A,A)')  &
                     'ERROR: DEVICE ',TRIM(DEVICE(NN)%ID)," must use QUANTITY='DENSITY' and a SPEC_ID"
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
               NNN = NNN + 1
            ENDIF
         ENDDO
         ALLOCATE(DV%DEVC_INDEX(NNN),STAT=IZERO)
         CALL ChkMemErr('READ','DV%DEVC_INDEX',IZERO)
         DV%DEVC_INDEX = -1
         ALLOCATE(DV%YY_SOOT(NNN,0:100))
         CALL ChkMemErr('READ','DV%YY_SOOT',IZERO)
         DV%YY_SOOT = 0._EB
         ALLOCATE(DV%TIME_ARRAY(0:100))
         CALL ChkMemErr('READ','DV%TIME_ARRAY',IZERO)
         DV%TIME_ARRAY = 0._EB
         DV%TOTAL_FLOWRATE = DV%BYPASS_FLOWRATE
         DV%DT             = -1._EB
         DV%N_INPUTS = NNN
         NNN = 1
         DO NN=1,N_DEVC
            IF (DEVICE(NN)%DEVC_ID==DV%ID) THEN
               DV%TOTAL_FLOWRATE  = DV%TOTAL_FLOWRATE + DEVICE(NN)%FLOWRATE
               DV%DT = MAX(DV%DT,DEVICE(NN)%DELAY)
               IF (NN > N) THEN
                  WRITE(MESSAGE,'(A,A,A)')  &
                     'ERROR: ASPIRATION DEVICE ',TRIM(DV%ID),' is not listed after all its inputs'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
               DV%DEVC_INDEX(NNN)     = NN
               NNN = NNN + 1
            ENDIF
         ENDDO
         DV%DT = DV%DT * 0.01_EB

      CASE ('FED')
         IF (DV%STATISTICS /= 'null' .AND. DV%STATISTICS /= 'TIME INTEGRAL') THEN
            WRITE(MESSAGE,'(A)') 'ERROR: Only TIME INTEGRAL statistics can be used with FED devices'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         DV%STATISTICS = 'TIME INTEGRAL'

   END SELECT SPECIAL_QUANTITIES

ENDDO PROC_DEVC_LOOP

! Set up arrays to handle devices that need a radiation heat flux other than at a boundary cell

IF (GAS_CELL_RAD_FLUX) ALLOCATE (GAS_CELL_RAD_DEVC_INDEX(N_GAS_CELL_RAD_DEVC))
NNN = 0
DO NN = 1, N_DEVC
   IF (DEVICE(NN)%GAS_CELL_RAD_FLUX) THEN 
      NNN = NNN + 1
      GAS_CELL_RAD_DEVC_INDEX(NNN) = NN
   ENDIF
ENDDO

END SUBROUTINE PROC_DEVC


SUBROUTINE READ_PROF
 
INTEGER :: N,NM,MESH_NUMBER,NN,N_PROFO,IOR
REAL(EB) :: XYZ(3)
CHARACTER(30) :: QUANTITY
TYPE (PROFILE_TYPE), POINTER :: PF
NAMELIST /PROF/ XYZ,QUANTITY,IOR,ID,FYI
 
N_PROF = 0
REWIND(LU_INPUT)
COUNT_PROF_LOOP: DO
   CALL CHECKREAD('PROF',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_PROF_LOOP
   READ(LU_INPUT,NML=PROF,END=11,ERR=12,IOSTAT=IOS)
   N_PROF = N_PROF + 1
   12 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I4)') 'ERROR: Problem with PROF no.',N_PROF+1
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO COUNT_PROF_LOOP
11 REWIND(LU_INPUT)

IF (N_PROF==0) RETURN
 
ALLOCATE(PROFILE(N_PROF),STAT=IZERO)
CALL ChkMemErr('READ','PROFILE',IZERO)
 
PROFILE(1:N_PROF)%QUANTITY = 'TEMPERATURE'
PROFILE(1:N_PROF)%IOR   = 0
PROFILE(1:N_PROF)%IW    = 0
 
N_PROFO = N_PROF
N       = 0
 
PROF_LOOP: DO NN=1,N_PROFO
   N    = N+1
   IOR  = 0
   SELECT CASE(N)
      CASE(1:9)        
         WRITE(ID,'(A,I1)') 'PROFILE ',N
      CASE(10:99)      
         WRITE(ID,'(A,I2)') 'PROFILE ',N
      CASE(100:999)    
         WRITE(ID,'(A,I3)') 'PROFILE ',N
   END SELECT
 
   CALL CHECKREAD('PROF',LU_INPUT,IOS)
   IF (IOS==1) EXIT PROF_LOOP
   READ(LU_INPUT,PROF) 
 
! Check for bad PROF quantities or coordinates

   IF (IOR==0) THEN
      WRITE(MESSAGE,'(A,I4,A)') 'ERROR: Specify orientation of PROF ' ,NN,' using the parameter IOR'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   BAD = .FALSE.
 
   MESH_LOOP: DO NM=1,NMESHES
      IF (.NOT.EVACUATION_ONLY(NM)) THEN      
         M=>MESHES(NM)
         IF (XYZ(1)>=M%XS .AND. XYZ(1)<=M%XF .AND. XYZ(2)>=M%YS .AND. XYZ(2)<=M%YF .AND. XYZ(3)>=M%ZS .AND. XYZ(3)<=M%ZF) THEN
            MESH_NUMBER = NM
            EXIT MESH_LOOP
         ENDIF
      ENDIF
      IF (NM==NMESHES) BAD = .TRUE.
   ENDDO MESH_LOOP
 
   IF (BAD) THEN
      N      = N-1
      N_PROF = N_PROF-1
      CYCLE PROF_LOOP
   ENDIF
 
! Assign parameters to the PROFILE array
 
   PF => PROFILE(N)
   PF%ORDINAL = NN
   PF%MESH    = MESH_NUMBER
   PF%ID   = ID
   PF%QUANTITY = QUANTITY
   PF%X       = XYZ(1)
   PF%Y       = XYZ(2)
   PF%Z       = XYZ(3)
   PF%IOR     = IOR
 
ENDDO PROF_LOOP
REWIND(LU_INPUT)
 
END SUBROUTINE READ_PROF



SUBROUTINE READ_ISOF
 
REAL(EB) :: VALUE(10)
CHARACTER(30) :: QUANTITY,COLOR_QUANTITY,SPEC_ID,COLOR_SPEC_ID
INTEGER :: REDUCE_TRIANGLES,N,I_DUM
TYPE(ISOSURFACE_FILE_TYPE), POINTER :: IS
NAMELIST /ISOF/ QUANTITY,FYI,VALUE,REDUCE_TRIANGLES,COLOR_QUANTITY,SPEC_ID,COLOR_SPEC_ID
 
N_ISOF = 0
REWIND(LU_INPUT)
COUNT_ISOF_LOOP: DO
   CALL CHECKREAD('ISOF',LU_INPUT,IOS) 
   IF (IOS==1) EXIT COUNT_ISOF_LOOP
   READ(LU_INPUT,NML=ISOF,END=9,ERR=10,IOSTAT=IOS)
   N_ISOF = N_ISOF + 1
   10 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I2)') 'ERROR: Problem with ISOF no.',N_ISOF
      CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_ISOF_LOOP
9 REWIND(LU_INPUT)
 
ALLOCATE(ISOSURFACE_FILE(N_ISOF),STAT=IZERO)
CALL ChkMemErr('READ','ISOSURFACE_FILE',IZERO)

READ_ISOF_LOOP: DO N=1,N_ISOF
   IS => ISOSURFACE_FILE(N)
   QUANTITY         = 'null'
   COLOR_QUANTITY   = 'null'
   SPEC_ID          = 'null'
   COLOR_SPEC_ID    = 'null'
   VALUE            = -999._EB
   REDUCE_TRIANGLES = 1
 
   CALL CHECKREAD('ISOF',LU_INPUT,IOS) 
   IF (IOS==1) EXIT READ_ISOF_LOOP
   READ(LU_INPUT,ISOF) 
 
   IS%REDUCE_TRIANGLES = REDUCE_TRIANGLES

   CALL GET_QUANTITY_INDEX(IS%SMOKEVIEW_LABEL,IS%SMOKEVIEW_BAR_LABEL,IS%INDEX,IS%SPEC_INDEX,I_DUM,'ISOF',QUANTITY,SPEC_ID,'null')
   IF (COLOR_QUANTITY/='null') &
      CALL GET_QUANTITY_INDEX(IS%SMOKEVIEW_LABEL,IS%SMOKEVIEW_BAR_LABEL,IS%INDEX2,IS%SPEC_INDEX2,I_DUM,'ISOF', &
                              COLOR_QUANTITY,COLOR_SPEC_ID,'null')

   VALUE_LOOP: DO I=1,10
      IF (VALUE(I)==-999._EB) EXIT VALUE_LOOP
      IS%N_VALUES = I
      IS%VALUE(I) = VALUE(I)
   ENDDO VALUE_LOOP
 
ENDDO READ_ISOF_LOOP
REWIND(LU_INPUT)
 
END SUBROUTINE READ_ISOF
 
 
SUBROUTINE READ_SLCF

REAL(EB) :: MAXIMUM_VALUE,MINIMUM_VALUE
REAL(EB) :: AGL_SLICE
INTEGER :: N,NN,NM,MESH_NUMBER,N_SLCF_O,NITER,ITER
LOGICAL :: VECTOR,RLE,TWO_BYTE,CELL_CENTERED, FIRE_LINE, EVACUATION
CHARACTER(30) :: QUANTITY,SPEC_ID,PART_ID
TYPE (SLICE_TYPE), POINTER :: SL
NAMELIST /SLCF/ XB,QUANTITY,FYI,PBX,PBY,PBZ,VECTOR,MESH_NUMBER,RLE,MAXIMUM_VALUE,MINIMUM_VALUE,TWO_BYTE,SPEC_ID, &
                AGL_SLICE,PART_ID,CELL_CENTERED, FIRE_LINE, EVACUATION

MESH_LOOP: DO NM=1,NMESHES

   M=>MESHES(NM)
   CALL POINT_TO_MESH(NM)

   N_SLCF   = 0
   N_SLCF_O = 0
   REWIND(LU_INPUT)
   COUNT_SLCF_LOOP: DO
      VECTOR  = .FALSE.
      EVACUATION  = .FALSE.
      MESH_NUMBER=NM
      CALL CHECKREAD('SLCF',LU_INPUT,IOS)
      IF (IOS==1) EXIT COUNT_SLCF_LOOP
      READ(LU_INPUT,NML=SLCF,END=9,ERR=10,IOSTAT=IOS)
      N_SLCF_O = N_SLCF_O + 1
      IF (MESH_NUMBER/=NM) CYCLE COUNT_SLCF_LOOP
      IF (.NOT.EVACUATION_ONLY(NM) .AND.      EVACUATION) CYCLE COUNT_SLCF_LOOP
      IF (     EVACUATION_ONLY(NM) .AND. .NOT.EVACUATION) CYCLE COUNT_SLCF_LOOP
      N_SLCF  = N_SLCF + 1
      IF (VECTOR .AND. TWO_D) N_SLCF = N_SLCF + 2
      IF (VECTOR .AND. .NOT. TWO_D) N_SLCF = N_SLCF + 3
      10 IF (IOS>0) THEN
         ! WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with SLCF no.',N_SLCF+1
         WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with SLCF no.',N_SLCF_O+1
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDDO COUNT_SLCF_LOOP
   9 CONTINUE   

   ALLOCATE(M%SLICE(N_SLCF),STAT=IZERO)
   CALL ChkMemErr('READ','ISP1',IZERO)
   CALL POINT_TO_MESH(NM)  ! Reset the pointers after the allocation

   N = 0
   N_TERRAIN_SLCF = 0

   REWIND(LU_INPUT)
   SLCF_LOOP: DO NN=1,N_SLCF_O
      QUANTITY = 'null'
      PBX      = -1.E6_EB
      PBY      = -1.E6_EB
      PBZ      = -1.E6_EB
      VECTOR   = .FALSE.
      MESH_NUMBER=NM
      RLE      = .FALSE.
      MINIMUM_VALUE = 0._EB
      MAXIMUM_VALUE = 0._EB
      TWO_BYTE = .FALSE.
      AGL_SLICE = -1._EB
      SPEC_ID  = 'null'
      PART_ID  = 'null'
      CELL_CENTERED = .FALSE.
      FIRE_LINE=.FALSE.
      EVACUATION  = .FALSE.
 
      CALL CHECKREAD('SLCF',LU_INPUT,IOS)
      IF (IOS==1) EXIT SLCF_LOOP
      READ(LU_INPUT,SLCF) 
      IF (MESH_NUMBER/=NM) CYCLE SLCF_LOOP
      IF (.NOT.EVACUATION_ONLY(NM) .AND.      EVACUATION) CYCLE SLCF_LOOP
      IF (     EVACUATION_ONLY(NM) .AND. .NOT.EVACUATION) CYCLE SLCF_LOOP
 
      IF (PBX>-1.E5_EB .OR. PBY>-1.E5_EB .OR. PBZ>-1.E5_EB) THEN
         XB(1) = XS
         XB(2) = XF
         XB(3) = YS
         XB(4) = YF
         XB(5) = ZS
         XB(6) = ZF
         IF (PBX>-1.E5_EB) XB(1:2) = PBX
         IF (PBY>-1.E5_EB) XB(3:4) = PBY
         IF (PBZ>-1.E5_EB) XB(5:6) = PBZ
      ENDIF
 
      XB(1) = MAX(XB(1),XS)
      XB(2) = MIN(XB(2),XF)
      XB(3) = MAX(XB(3),YS)
      XB(4) = MIN(XB(4),YF)
      XB(5) = MAX(XB(5),ZS)
      XB(6) = MIN(XB(6),ZF)
 
      ! Reject a slice if it is beyond the bounds of the current mesh
 
      IF (XB(1)>XF .OR. XB(2)<XS .OR. XB(3)>YF .OR. XB(4)<YS .OR. XB(5)>ZF .OR. XB(6)<ZS) THEN
         N_SLCF = N_SLCF - 1
         IF (VECTOR .AND. TWO_D) N_SLCF = N_SLCF - 2
         IF (VECTOR .AND. .NOT. TWO_D) N_SLCF = N_SLCF - 3
         CYCLE SLCF_LOOP
      ENDIF
 
      ! Process vector quantities
 
      NITER = 1
      IF (VECTOR .AND. TWO_D) NITER = 3
      IF (VECTOR .AND. .NOT. TWO_D)  NITER = 4
 
      VECTORLOOP: DO ITER=1,NITER
         N = N + 1
         SL=>SLICE(N)
         SL%TWO_BYTE = TWO_BYTE
         SL%I1 = NINT( GINV(XB(1)-XS,1,NM)*RDXI)
         SL%I2 = NINT( GINV(XB(2)-XS,1,NM)*RDXI)
         SL%J1 = NINT( GINV(XB(3)-YS,2,NM)*RDETA)
         SL%J2 = NINT( GINV(XB(4)-YS,2,NM)*RDETA)
         SL%K1 = NINT( GINV(XB(5)-ZS,3,NM)*RDZETA)
         SL%K2 = NINT( GINV(XB(6)-ZS,3,NM)*RDZETA)
         SL%RLE = RLE
         SL%MINMAX(1) = MINIMUM_VALUE
         SL%MINMAX(2) = MAXIMUM_VALUE
         IF (ITER==2)                    QUANTITY = 'U-VELOCITY' 
         IF (ITER==3 .AND. .NOT. TWO_D)  QUANTITY = 'V-VELOCITY' 
         IF (ITER==3 .AND. TWO_D)        QUANTITY = 'W-VELOCITY' 
         IF (ITER==4)                    QUANTITY = 'W-VELOCITY'
         IF (ITER==1 .AND. FIRE_LINE)    QUANTITY = 'TEMPERATURE'
         IF ( ITER==1) THEN
            SL%FIRE_LINE = FIRE_LINE
           ELSE
            SL%FIRE_LINE = .FALSE.
         ENDIF
         CALL GET_QUANTITY_INDEX(SL%SMOKEVIEW_LABEL,SL%SMOKEVIEW_BAR_LABEL,SL%INDEX,SL%SPEC_INDEX,SL%PART_INDEX, &
                                 'SLCF',QUANTITY,SPEC_ID,PART_ID)
         !for terrain slices, AGL=above ground level
         !FIRE_LINE==.TRUE. ==> a terrain slice at one grid cell above ground with quantity temperature
         !                      smokeview will display only regions where temperature is above 200 C.
         !                      This is currently hard wired.
         IF (ITER == 1 .AND. (AGL_SLICE > -1._EB .OR. FIRE_LINE ) ) THEN 
            SL%TERRAIN_SLICE = .TRUE.
            IF (FIRE_LINE) THEN
               SL%SMOKEVIEW_LABEL = "Fire line"
               SL%SMOKEVIEW_BAR_LABEL = "Fire_line"
            ENDIF 
            IF (AGL_SLICE .LE. -1._EB .AND. FIRE_LINE) THEN
               AGL_SLICE = M%Z(1) - M%Z(0)
            ENDIF
            SL%SLICE_AGL     = AGL_SLICE
            N_TERRAIN_SLCF   = N_TERRAIN_SLCF + 1
         ENDIF
         IF (ITER==2 .OR. ITER==3 .OR. ITER ==4) THEN
            IF (SLICE(N-1)%TERRAIN_SLICE) THEN
               SL%TERRAIN_SLICE =  .TRUE. 
               SL%SLICE_AGL     = SLICE(N-1)%SLICE_AGL
               N_TERRAIN_SLCF   = N_TERRAIN_SLCF + 1
            ENDIF
         ENDIF
         
         ! disable cell centered for velocity
         IF (QUANTITY=='VELOCITY'   .OR. &
             QUANTITY=='U-VELOCITY' .OR. &
             QUANTITY=='V-VELOCITY' .OR. &
             QUANTITY=='W-VELOCITY') THEN
             CELL_CENTERED = .FALSE.
         ENDIF
         SL%CELL_CENTERED = CELL_CENTERED
         
      ENDDO VECTORLOOP
  
   ENDDO SLCF_LOOP

   ALLOCATE(M%K_AGL_SLICE(0:IBP1,0:JBP1,N_TERRAIN_SLCF),STAT=IZERO)
   CALL ChkMemErr('READ','K_AGL_SLICE',IZERO)
   M%K_AGL_SLICE = 0.0_EB
   N = 0
   DO NN = 1, N_SLCF
      SL=>SLICE(NN)
      IF(SL%TERRAIN_SLICE) THEN
        TERRAIN_CASE = .TRUE.
        N = N + 1 
!        M%K_AGL_SLICE(0:IBP1,0:JBP1,N) =  MAX(SL%SLICE_AGL*M%RDZ(1),1._EB)
        M%K_AGL_SLICE(0:IBP1,0:JBP1,N) =  SL%SLICE_AGL*M%RDZ(1)
!subtract one because bottom of doamin will be accounted for when cycling
!through walls cells
        M%K_AGL_SLICE(0:IBP1,0:JBP1,N) =  M%K_AGL_SLICE(0:IBP1,0:JBP1,N) - 1  
      ENDIF
   ENDDO

   N_SLCF_MAX = MAX(N_SLCF_MAX,N_SLCF) 
 
ENDDO MESH_LOOP

END SUBROUTINE READ_SLCF


SUBROUTINE READ_BNDF

USE DEVICE_VARIABLES
INTEGER :: N
LOGICAL :: CELL_CENTERED
CHARACTER(30) :: QUANTITY,PROP_ID,SPEC_ID,PART_ID
NAMELIST /BNDF/ CELL_CENTERED,QUANTITY,FYI,PROP_ID,SPEC_ID,PART_ID,RECOUNT_DRIP
TYPE(BOUNDARY_FILE_TYPE), POINTER :: BF
 
N_BNDF = 0
REWIND(LU_INPUT)
COUNT_BNDF_LOOP: DO
   CALL CHECKREAD('BNDF',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_BNDF_LOOP
   READ(LU_INPUT,NML=BNDF,END=209,ERR=210,IOSTAT=IOS)
   N_BNDF = N_BNDF + 1
   210 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with BNDF line')
ENDDO COUNT_BNDF_LOOP
209 REWIND(LU_INPUT)
 
RECOUNT_DRIP = .FALSE.
 
ALLOCATE(BOUNDARY_FILE(N_BNDF),STAT=IZERO)
CALL ChkMemErr('READ','BOUNDARY_FILE',IZERO)
 
READ_BNDF_LOOP: DO N=1,N_BNDF
   BF => BOUNDARY_FILE(N)
   CELL_CENTERED = .FALSE.
   PART_ID  = 'null'
   PROP_ID  = 'null'
   SPEC_ID  = 'null'
   QUANTITY = 'WALL_TEMPERATURE'
   CALL CHECKREAD('BNDF',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_BNDF_LOOP
   READ(LU_INPUT,BNDF)
   
   ! Look to see if output QUANTITY exists

   CALL GET_QUANTITY_INDEX(BF%SMOKEVIEW_LABEL,BF%SMOKEVIEW_BAR_LABEL,BF%INDEX,BF%SPEC_INDEX,BF%PART_INDEX, &
                           'BNDF',QUANTITY,SPEC_ID,PART_ID)

   ! Assign miscellaneous attributes to the boundary file

   BF%CELL_CENTERED = CELL_CENTERED

   ! Check to see if PROP_ID exists

   BF%PROP_INDEX = 0
   IF(PROP_ID/='null')  CALL GET_PROPERTY_INDEX(BF%PROP_INDEX,'BNDF',PROP_ID)

ENDDO READ_BNDF_LOOP
REWIND(LU_INPUT)
 
END SUBROUTINE READ_BNDF
 
 
SUBROUTINE DEFINE_OUTPUT_QUANTITIES
 
! Define all OUTPUT_QUANTITYs and their attributes

INTEGER :: N

ALLOCATE(OUTPUT_QUANTITY(-N_OUTPUT_QUANTITIES:N_OUTPUT_QUANTITIES),STAT=IZERO)
CALL ChkMemErr('READ','OUTPUT_QUANTITY',IZERO) 

DO N=-N_OUTPUT_QUANTITIES,N_OUTPUT_QUANTITIES
   IF (N < 0) THEN
      OUTPUT_QUANTITY(N)%SOLID_PHASE = .TRUE.
      OUTPUT_QUANTITY(N)%ISOF_APPROPRIATE = .FALSE.
      OUTPUT_QUANTITY(N)%SLCF_APPROPRIATE = .FALSE.
      OUTPUT_QUANTITY(N)%BNDF_APPROPRIATE = .TRUE.
   ENDIF
ENDDO
 
OUTPUT_QUANTITY(0)%NAME        = 'SMOKE/WATER'             
OUTPUT_QUANTITY(0)%UNITS       = '  '                      
OUTPUT_QUANTITY(0)%SHORT_NAME  = '  '

OUTPUT_QUANTITY(1)%NAME        = 'DENSITY'                 
OUTPUT_QUANTITY(1)%UNITS       = 'kg/m3'                   
OUTPUT_QUANTITY(1)%SHORT_NAME  = 'rho'

OUTPUT_QUANTITY(2)%NAME        = 'F_X'                     
OUTPUT_QUANTITY(2)%UNITS       = 'm/s2'                    
OUTPUT_QUANTITY(2)%SHORT_NAME  = 'f_x'
OUTPUT_QUANTITY(2)%CELL_POSITION = CELL_FACE
OUTPUT_QUANTITY(2)%IOR = 1

OUTPUT_QUANTITY(3)%NAME  = 'F_Y'                     
OUTPUT_QUANTITY(3)%UNITS  = 'm/s2'                    
OUTPUT_QUANTITY(3)%SHORT_NAME  = 'f_y'
OUTPUT_QUANTITY(3)%CELL_POSITION = CELL_FACE
OUTPUT_QUANTITY(3)%IOR = 2

OUTPUT_QUANTITY(4)%NAME  = 'F_Z'                     
OUTPUT_QUANTITY(4)%UNITS  = 'm/s2'                    
OUTPUT_QUANTITY(4)%SHORT_NAME  = 'f_z'
OUTPUT_QUANTITY(4)%CELL_POSITION = CELL_FACE
OUTPUT_QUANTITY(4)%IOR = 3

OUTPUT_QUANTITY(5)%NAME  = 'TEMPERATURE'             
OUTPUT_QUANTITY(5)%UNITS  = 'C'                       
OUTPUT_QUANTITY(5)%SHORT_NAME  = 'temp'

OUTPUT_QUANTITY(6)%NAME  = 'U-VELOCITY'              
OUTPUT_QUANTITY(6)%UNITS  = 'm/s'                    
OUTPUT_QUANTITY(6)%SHORT_NAME  = 'U-VEL'
OUTPUT_QUANTITY(6)%CELL_POSITION = CELL_FACE
OUTPUT_QUANTITY(6)%IOR = 1

OUTPUT_QUANTITY(7)%NAME  = 'V-VELOCITY'              
OUTPUT_QUANTITY(7)%UNITS = 'm/s'                     
OUTPUT_QUANTITY(7)%SHORT_NAME  = 'V-VEL'
OUTPUT_QUANTITY(7)%CELL_POSITION = CELL_FACE
OUTPUT_QUANTITY(7)%IOR = 2

OUTPUT_QUANTITY(8)%NAME  = 'W-VELOCITY'              
OUTPUT_QUANTITY(8)%UNITS  = 'm/s'                     
OUTPUT_QUANTITY(8)%SHORT_NAME  = 'W-VEL'
OUTPUT_QUANTITY(8)%CELL_POSITION = CELL_FACE
OUTPUT_QUANTITY(8)%IOR = 3

OUTPUT_QUANTITY(6:8)%PART_APPROPRIATE = .TRUE.

OUTPUT_QUANTITY(9)%NAME  = 'PRESSURE'                
OUTPUT_QUANTITY(9)%UNITS  = 'Pa'                      
OUTPUT_QUANTITY(9)%SHORT_NAME  = 'pres'

OUTPUT_QUANTITY(10)%NAME = 'VELOCITY'                
OUTPUT_QUANTITY(10)%UNITS = 'm/s'                     
OUTPUT_QUANTITY(10)%SHORT_NAME = 'vel'

OUTPUT_QUANTITY(11)%NAME = 'HRRPUV'                  
OUTPUT_QUANTITY(11)%UNITS = 'kW/m3'                   
OUTPUT_QUANTITY(11)%SHORT_NAME = 'hrrpuv'

OUTPUT_QUANTITY(12)%NAME = 'H'                       
OUTPUT_QUANTITY(12)%UNITS = '(m/s)^2'                 
OUTPUT_QUANTITY(12)%SHORT_NAME = 'head'

OUTPUT_QUANTITY(13)%NAME = 'MIXTURE FRACTION'              
OUTPUT_QUANTITY(13)%UNITS = 'kg/kg'                     
OUTPUT_QUANTITY(13)%SHORT_NAME = 'Z'
OUTPUT_QUANTITY(13)%MASS_FRACTION = .TRUE.
OUTPUT_QUANTITY(13)%MIXTURE_FRACTION_ONLY = .TRUE.

OUTPUT_QUANTITY(14)%NAME = 'DIVERGENCE'              
OUTPUT_QUANTITY(14)%UNITS = '1/s'                     
OUTPUT_QUANTITY(14)%SHORT_NAME = 'div'

OUTPUT_QUANTITY(15)%NAME = 'MIXING TIME'              
OUTPUT_QUANTITY(15)%UNITS = 's'                     
OUTPUT_QUANTITY(15)%SHORT_NAME = 'mix'

OUTPUT_QUANTITY(16)%NAME = 'ABSORPTION COEFFICIENT'  
OUTPUT_QUANTITY(16)%UNITS = '1/m'                     
OUTPUT_QUANTITY(16)%SHORT_NAME = 'kappa'

OUTPUT_QUANTITY(17)%NAME = 'VISCOSITY'               
OUTPUT_QUANTITY(17)%UNITS = 'kg/m/s'                  
OUTPUT_QUANTITY(17)%SHORT_NAME = 'visc'
 
OUTPUT_QUANTITY(18)%NAME = 'INTEGRATED INTENSITY'       
OUTPUT_QUANTITY(18)%UNITS = 'kW/m^2'                  
OUTPUT_QUANTITY(18)%SHORT_NAME = 'U'

OUTPUT_QUANTITY(19)%NAME = 'RADIATION LOSS'          
OUTPUT_QUANTITY(19)%UNITS = 'kW/m3'                   
OUTPUT_QUANTITY(19)%SHORT_NAME = 'loss'

OUTPUT_QUANTITY(20)%NAME = 'WATER RADIATION LOSS'    
OUTPUT_QUANTITY(20)%UNITS = 'kW/m3'                   
OUTPUT_QUANTITY(20)%SHORT_NAME = 'H2O_rad'

OUTPUT_QUANTITY(21)%NAME = 'RELATIVE HUMIDITY'
OUTPUT_QUANTITY(21)%UNITS = '%'
OUTPUT_QUANTITY(21)%SHORT_NAME = 'humid'

OUTPUT_QUANTITY(22)%NAME = 'HP'                       
OUTPUT_QUANTITY(22)%UNITS = '(m/s)^2'                 
OUTPUT_QUANTITY(22)%SHORT_NAME = 'hp'

OUTPUT_QUANTITY(23)%NAME = 'KINETIC ENERGY'                
OUTPUT_QUANTITY(23)%UNITS = 'm^2/s^2'                     
OUTPUT_QUANTITY(23)%SHORT_NAME = 'ke'

! Strain and Vorticity
 
OUTPUT_QUANTITY(24)%NAME = 'STRAIN RATE X'           
OUTPUT_QUANTITY(24)%UNITS = '1/s'                     
OUTPUT_QUANTITY(24)%SHORT_NAME = 'strain_x'

OUTPUT_QUANTITY(25)%NAME = 'STRAIN RATE Y'           
OUTPUT_QUANTITY(25)%UNITS = '1/s'                     
OUTPUT_QUANTITY(25)%SHORT_NAME = 'strain_y'

OUTPUT_QUANTITY(26)%NAME = 'STRAIN RATE Z'           
OUTPUT_QUANTITY(26)%UNITS = '1/s'                     
OUTPUT_QUANTITY(26)%SHORT_NAME = 'strain_z'

OUTPUT_QUANTITY(27)%NAME = 'VORTICITY X'             
OUTPUT_QUANTITY(27)%UNITS = '1/s'                     
OUTPUT_QUANTITY(27)%SHORT_NAME = 'vort_x'  

OUTPUT_QUANTITY(28)%NAME = 'VORTICITY Y'             
OUTPUT_QUANTITY(28)%UNITS = '1/s'                     
OUTPUT_QUANTITY(28)%SHORT_NAME = 'vort_y'  

OUTPUT_QUANTITY(29)%NAME = 'VORTICITY Z'             
OUTPUT_QUANTITY(29)%UNITS = '1/s'                     
OUTPUT_QUANTITY(29)%SHORT_NAME = 'vort_z'  

OUTPUT_QUANTITY(24:29)%CELL_POSITION = CELL_EDGE

OUTPUT_QUANTITY(30)%NAME = 'C_DYNSMAG'
OUTPUT_QUANTITY(30)%UNITS = '  '
OUTPUT_QUANTITY(30)%SHORT_NAME = 'c_smag'

OUTPUT_QUANTITY(31)%NAME = 'SPECIFIC HEAT'
OUTPUT_QUANTITY(31)%UNITS = 'kJ/kg/K'
OUTPUT_QUANTITY(31)%SHORT_NAME = 'c_p'

OUTPUT_QUANTITY(32)%NAME = 'HRRPUA'
OUTPUT_QUANTITY(32)%UNITS = 'kW/m2'
OUTPUT_QUANTITY(32)%SHORT_NAME = 'hrrpua'  

OUTPUT_QUANTITY(33)%NAME = 'CONDUCTIVITY'
OUTPUT_QUANTITY(33)%UNITS = 'W/m/K'
OUTPUT_QUANTITY(33)%SHORT_NAME = 'k'

! QUANTITY's the refer to droplet properties that are used to color the droplets in Smokeview

OUTPUT_QUANTITY(34)%NAME = 'DROPLET DIAMETER'                     
OUTPUT_QUANTITY(34)%UNITS = 'mu-m'                       
OUTPUT_QUANTITY(34)%SHORT_NAME = 'diam'
 
OUTPUT_QUANTITY(35)%NAME = 'DROPLET VELOCITY'                     
OUTPUT_QUANTITY(35)%UNITS = 'm/s'                       
OUTPUT_QUANTITY(35)%SHORT_NAME = 'vel'
 
OUTPUT_QUANTITY(36)%NAME = 'DROPLET PHASE'                     
OUTPUT_QUANTITY(36)%UNITS = ' '                       
OUTPUT_QUANTITY(36)%SHORT_NAME = 'ior'
 
OUTPUT_QUANTITY(37)%NAME = 'DROPLET TEMPERATURE'                     
OUTPUT_QUANTITY(37)%UNITS = 'C'                       
OUTPUT_QUANTITY(37)%SHORT_NAME = 'temp'
 
OUTPUT_QUANTITY(38)%NAME = 'DROPLET MASS'                     
OUTPUT_QUANTITY(38)%UNITS = 'mu-g'                       
OUTPUT_QUANTITY(38)%SHORT_NAME = 'mass'
 
OUTPUT_QUANTITY(39)%NAME = 'DROPLET AGE'                     
OUTPUT_QUANTITY(39)%UNITS = 's'                       
OUTPUT_QUANTITY(39)%SHORT_NAME = 'age'

OUTPUT_QUANTITY(34:39)%PART_APPROPRIATE = .TRUE.
OUTPUT_QUANTITY(34:39)%SLCF_APPROPRIATE = .FALSE.
OUTPUT_QUANTITY(34:39)%DEVC_APPROPRIATE = .FALSE.

OUTPUT_QUANTITY(40)%NAME = 'MOLECULAR WEIGHT'                     
OUTPUT_QUANTITY(40)%UNITS = 'g/mol'                       
OUTPUT_QUANTITY(40)%SHORT_NAME = 'MW'

! Time and benchmarking stats

OUTPUT_QUANTITY(41)%NAME = 'TIME'
OUTPUT_QUANTITY(41)%UNITS = 's'
OUTPUT_QUANTITY(41)%SHORT_NAME = 't'

OUTPUT_QUANTITY(42)%NAME = 'TIME STEP'                       
OUTPUT_QUANTITY(42)%UNITS = 's'                 
OUTPUT_QUANTITY(42)%SHORT_NAME = 'dt'

OUTPUT_QUANTITY(43)%NAME = 'WALL CLOCK TIME'                       
OUTPUT_QUANTITY(43)%UNITS = 's'                 
OUTPUT_QUANTITY(43)%SHORT_NAME = 't_wall'

OUTPUT_QUANTITY(44)%NAME = 'CPU TIME'                       
OUTPUT_QUANTITY(44)%UNITS = 's'                 
OUTPUT_QUANTITY(44)%SHORT_NAME = 't_cpu'

OUTPUT_QUANTITY(45)%NAME = 'ITERATION'                       
OUTPUT_QUANTITY(45)%UNITS = ' '                 
OUTPUT_QUANTITY(45)%SHORT_NAME = 'step'

OUTPUT_QUANTITY(41:45)%SLCF_APPROPRIATE = .FALSE.
OUTPUT_QUANTITY(41:45)%ISOF_APPROPRIATE = .FALSE.

! Enthalpy

OUTPUT_QUANTITY(46)%NAME = 'SPECIFIC ENTHALPY'                       
OUTPUT_QUANTITY(46)%UNITS = 'kJ/kg'                 
OUTPUT_QUANTITY(46)%SHORT_NAME = 'h'

OUTPUT_QUANTITY(47)%NAME = 'ENTHALPY'                       
OUTPUT_QUANTITY(47)%UNITS = 'kJ/m^3'                 
OUTPUT_QUANTITY(47)%SHORT_NAME = 'Q'

OUTPUT_QUANTITY(48)%NAME = 'AVERAGE SPECIFIC HEAT'                       
OUTPUT_QUANTITY(48)%UNITS = 'kJ/kg/K'                 
OUTPUT_QUANTITY(48)%SHORT_NAME = 'c_p_bar'

! Measures of turbulence resolution (diagnostics)

OUTPUT_QUANTITY(50)%NAME = 'TURBULENCE RESOLUTION'
OUTPUT_QUANTITY(50)%UNITS = '  '
OUTPUT_QUANTITY(50)%SHORT_NAME = 'mtr'

!OUTPUT_QUANTITY(51)%NAME = 'SCALAR RESOLUTION'
!OUTPUT_QUANTITY(51)%UNITS = '  '
!OUTPUT_QUANTITY(51)%SPEC_ID_REQUIRED = .TRUE.
!OUTPUT_QUANTITY(51)%SHORT_NAME = 'msr'

! Miscellaneous outputs

OUTPUT_QUANTITY(60)%NAME = 'ACTUATED SPRINKLERS'
OUTPUT_QUANTITY(60)%UNITS = ' '
OUTPUT_QUANTITY(60)%SHORT_NAME = 'act'
OUTPUT_QUANTITY(60)%TIME_AVERAGED = .FALSE.
 
! Species-specific quantities that require a SPEC_ID

OUTPUT_QUANTITY(90)%NAME='MASS FRACTION'
OUTPUT_QUANTITY(90)%UNITS='kg/kg'
OUTPUT_QUANTITY(90)%MASS_FRACTION = .TRUE.
OUTPUT_QUANTITY(90)%SHORT_NAME = 'Y'

OUTPUT_QUANTITY(91)%NAME = 'MASS FLUX X'
OUTPUT_QUANTITY(91)%UNITS = 'kg/s/m2'
OUTPUT_QUANTITY(91)%SHORT_NAME = 'u*rho*Y'
OUTPUT_QUANTITY(91)%CELL_POSITION = CELL_FACE
OUTPUT_QUANTITY(91)%IOR = 1

OUTPUT_QUANTITY(92)%NAME = 'MASS FLUX Y'
OUTPUT_QUANTITY(92)%UNITS = 'kg/s/m2'
OUTPUT_QUANTITY(92)%SHORT_NAME = 'v*rho*Y'
OUTPUT_QUANTITY(92)%CELL_POSITION = CELL_FACE
OUTPUT_QUANTITY(92)%IOR = 2

OUTPUT_QUANTITY(93)%NAME = 'MASS FLUX Z'
OUTPUT_QUANTITY(93)%UNITS = 'kg/s/m2'
OUTPUT_QUANTITY(93)%SHORT_NAME = 'w*rho*Y'
OUTPUT_QUANTITY(93)%CELL_POSITION = CELL_FACE
OUTPUT_QUANTITY(93)%IOR = 3

OUTPUT_QUANTITY(94)%NAME = 'VOLUME FRACTION'
OUTPUT_QUANTITY(94)%UNITS = 'mol/mol'
OUTPUT_QUANTITY(94)%SHORT_NAME = 'X'

OUTPUT_QUANTITY(95)%NAME = 'VISIBILITY'
OUTPUT_QUANTITY(95)%UNITS = 'm'
OUTPUT_QUANTITY(95)%SHORT_NAME = 'VIS'

OUTPUT_QUANTITY(96)%NAME = 'SOOT VOLUME FRACTION'
OUTPUT_QUANTITY(96)%UNITS = ' '
OUTPUT_QUANTITY(96)%SHORT_NAME = 'f_v'

OUTPUT_QUANTITY(97)%NAME = 'EXTINCTION COEFFICIENT'
OUTPUT_QUANTITY(97)%UNITS = '1/m'
OUTPUT_QUANTITY(97)%SHORT_NAME = 'ext_coef'

OUTPUT_QUANTITY(98)%NAME = 'OPTICAL DENSITY'
OUTPUT_QUANTITY(98)%UNITS = '1/m'
OUTPUT_QUANTITY(98)%SHORT_NAME = 'OD'

OUTPUT_QUANTITY(90:98)%SPEC_ID_REQUIRED = .TRUE.

! Miscellaneous

OUTPUT_QUANTITY(100)%NAME  = 'PRESSURE ZONE'                   
OUTPUT_QUANTITY(100)%UNITS  = ' '                    
OUTPUT_QUANTITY(100)%SHORT_NAME  = 'zone'
 
! Integrated Quantities
 
OUTPUT_QUANTITY(104)%NAME  = 'HRR'                   
OUTPUT_QUANTITY(104)%UNITS  = 'kW'                    
OUTPUT_QUANTITY(104)%SHORT_NAME  = 'hrr'
OUTPUT_QUANTITY(104)%INTEGRATED  = .TRUE.

OUTPUT_QUANTITY(105)%NAME  = 'LAYER HEIGHT'          
OUTPUT_QUANTITY(105)%UNITS  = 'm'                     
OUTPUT_QUANTITY(105)%SHORT_NAME  = 'layer'

OUTPUT_QUANTITY(106)%NAME  = 'UPPER TEMPERATURE'     
OUTPUT_QUANTITY(106)%UNITS  = 'C'                     
OUTPUT_QUANTITY(106)%SHORT_NAME  = 'u-tmp'

OUTPUT_QUANTITY(107)%NAME  = 'LOWER TEMPERATURE'     
OUTPUT_QUANTITY(107)%UNITS  = 'C'                     
OUTPUT_QUANTITY(107)%SHORT_NAME  = 'l-tmp'

OUTPUT_QUANTITY(104:107)%INTEGRATED = .TRUE.
OUTPUT_QUANTITY(104:107)%SLCF_APPROPRIATE = .FALSE.

! Fractional Effective Dose (FED)

OUTPUT_QUANTITY(109)%NAME  = 'FED'
OUTPUT_QUANTITY(109)%UNITS = ''             
OUTPUT_QUANTITY(109)%SHORT_NAME  = 'fed'
OUTPUT_QUANTITY(109)%MIXTURE_FRACTION_ONLY = .TRUE.

! Model of a TC

OUTPUT_QUANTITY(110)%NAME  = 'THERMOCOUPLE'          
OUTPUT_QUANTITY(110)%UNITS = 'C'                     
OUTPUT_QUANTITY(110)%SHORT_NAME  = 'tc'
 
! Mass and Energy Flows

OUTPUT_QUANTITY(111)%NAME  = 'VOLUME FLOW'
OUTPUT_QUANTITY(111)%UNITS  = 'm3/s'
OUTPUT_QUANTITY(111)%SHORT_NAME  = 'vflow'

OUTPUT_QUANTITY(112)%NAME  = 'MASS FLOW'
OUTPUT_QUANTITY(112)%UNITS  = 'kg/s'
OUTPUT_QUANTITY(112)%SHORT_NAME  = 'mflow'

OUTPUT_QUANTITY(113)%NAME  = 'HEAT FLOW'
OUTPUT_QUANTITY(113)%UNITS = 'kW'
OUTPUT_QUANTITY(113)%SHORT_NAME  = 'hflow'

OUTPUT_QUANTITY(114)%NAME  = 'VOLUME FLOW +'         
OUTPUT_QUANTITY(114)%UNITS  = 'm3/s'                  
OUTPUT_QUANTITY(114)%SHORT_NAME  = 'vflow+'

OUTPUT_QUANTITY(115)%NAME  = 'MASS FLOW +'           
OUTPUT_QUANTITY(115)%UNITS  = 'kg/s'                  
OUTPUT_QUANTITY(115)%SHORT_NAME  = 'mflow+'

OUTPUT_QUANTITY(116)%NAME  = 'HEAT FLOW +'           
OUTPUT_QUANTITY(116)%UNITS  = 'kW'                    
OUTPUT_QUANTITY(116)%SHORT_NAME  = 'hflow+'

OUTPUT_QUANTITY(117)%NAME  = 'VOLUME FLOW -'         
OUTPUT_QUANTITY(117)%UNITS  = 'm3/s'                  
OUTPUT_QUANTITY(117)%SHORT_NAME  = 'vflow-'

OUTPUT_QUANTITY(118)%NAME  = 'MASS FLOW -'           
OUTPUT_QUANTITY(118)%UNITS  = 'kg/s'                  
OUTPUT_QUANTITY(118)%SHORT_NAME  = 'mflow-'

OUTPUT_QUANTITY(119)%NAME  = 'HEAT FLOW -'           
OUTPUT_QUANTITY(119)%UNITS  = 'kW'                    
OUTPUT_QUANTITY(119)%SHORT_NAME  = 'hflow-'

OUTPUT_QUANTITY(111:119)%INTEGRATED = .TRUE.
OUTPUT_QUANTITY(111:119)%SLCF_APPROPRIATE = .FALSE.

! Special integration of HRRPUV 

OUTPUT_QUANTITY(120)%NAME  = 'HRRPUL'
OUTPUT_QUANTITY(120)%UNITS  = 'kW/m'
OUTPUT_QUANTITY(120)%SHORT_NAME  = 'hrrpul'

OUTPUT_QUANTITY(120)%INTEGRATED = .TRUE.
OUTPUT_QUANTITY(120)%SLCF_APPROPRIATE = .TRUE.

! Sprinklers and Detectors

OUTPUT_QUANTITY(155)%NAME = 'PATH OBSCURATION'
OUTPUT_QUANTITY(155)%UNITS = '%'
OUTPUT_QUANTITY(155)%SHORT_NAME = 'total obs'
OUTPUT_QUANTITY(155)%INTEGRATED = .TRUE.
OUTPUT_QUANTITY(155)%SLCF_APPROPRIATE = .FALSE.

OUTPUT_QUANTITY(156)%NAME = 'SPRINKLER LINK TEMPERATURE'
OUTPUT_QUANTITY(156)%UNITS = 'C'
OUTPUT_QUANTITY(156)%SHORT_NAME = 'link'
OUTPUT_QUANTITY(156)%SLCF_APPROPRIATE = .FALSE.

OUTPUT_QUANTITY(157)%NAME = 'LINK TEMPERATURE'
OUTPUT_QUANTITY(157)%UNITS = 'C'
OUTPUT_QUANTITY(157)%SHORT_NAME = 'link'
OUTPUT_QUANTITY(157)%SLCF_APPROPRIATE = .FALSE.

OUTPUT_QUANTITY(158)%NAME = 'CHAMBER OBSCURATION'
OUTPUT_QUANTITY(158)%UNITS = '%/m'
OUTPUT_QUANTITY(158)%SHORT_NAME = 'obs'
OUTPUT_QUANTITY(158)%SLCF_APPROPRIATE = .FALSE.

OUTPUT_QUANTITY(160)%NAME = 'CONTROL'
OUTPUT_QUANTITY(160)%UNITS = ' '
OUTPUT_QUANTITY(160)%SHORT_NAME = 'output'
OUTPUT_QUANTITY(160)%SLCF_APPROPRIATE = .FALSE.

OUTPUT_QUANTITY(161)%NAME = 'ASPIRATION'
OUTPUT_QUANTITY(161)%UNITS = '%/m'
OUTPUT_QUANTITY(161)%SHORT_NAME = 'obs'
OUTPUT_QUANTITY(161)%SLCF_APPROPRIATE = .FALSE.

OUTPUT_QUANTITY(163)%NAME = 'RADIATIVE HEAT FLUX GAS'
OUTPUT_QUANTITY(163)%UNITS = 'km/m^2'
OUTPUT_QUANTITY(163)%SHORT_NAME = 'rad'
OUTPUT_QUANTITY(163)%SLCF_APPROPRIATE = .FALSE.

OUTPUT_QUANTITY(164)%NAME = 'CABLE TEMPERATURE'
OUTPUT_QUANTITY(164)%UNITS = 'C'
OUTPUT_QUANTITY(164)%SHORT_NAME = 'cable'
OUTPUT_QUANTITY(164)%SLCF_APPROPRIATE = .FALSE.

! Droplet mass quantities

OUTPUT_QUANTITY(170)%NAME = 'MPUV'
OUTPUT_QUANTITY(170)%UNITS = 'kg/m3'
OUTPUT_QUANTITY(170)%SHORT_NAME = 'mpuv'

OUTPUT_QUANTITY(171)%NAME = 'ADD'
OUTPUT_QUANTITY(171)%UNITS = 'mu-m'
OUTPUT_QUANTITY(171)%SHORT_NAME = 'radii'

OUTPUT_QUANTITY(172)%NAME = 'ADT'
OUTPUT_QUANTITY(172)%UNITS = 'C'
OUTPUT_QUANTITY(172)%SHORT_NAME = 'temp'

OUTPUT_QUANTITY(170:172)%ISOF_APPROPRIATE = .FALSE.
OUTPUT_QUANTITY(170:172)%PART_ID_REQUIRED = .TRUE.

OUTPUT_QUANTITY(173)%NAME = 'DROPLET FLUX X'
OUTPUT_QUANTITY(173)%UNITS = 'kg/s/m2'
OUTPUT_QUANTITY(173)%SHORT_NAME = 'flux_x'

OUTPUT_QUANTITY(174)%NAME = 'DROPLET FLUX Y'
OUTPUT_QUANTITY(174)%UNITS = 'kg/s/m2'
OUTPUT_QUANTITY(174)%SHORT_NAME = 'flux_y'

OUTPUT_QUANTITY(175)%NAME = 'DROPLET FLUX Z'
OUTPUT_QUANTITY(175)%UNITS = 'kg/s/m2'
OUTPUT_QUANTITY(175)%SHORT_NAME = 'flux_z'

OUTPUT_QUANTITY(173:175)%ISOF_APPROPRIATE = .FALSE.
OUTPUT_QUANTITY(173:175)%INTEGRATED_DROPLETS = .TRUE.

! PDPA

OUTPUT_QUANTITY(231)%NAME = 'PDPA'
OUTPUT_QUANTITY(231)%UNITS = ' '
OUTPUT_QUANTITY(231)%SHORT_NAME = ' '

OUTPUT_QUANTITY(231)%ISOF_APPROPRIATE = .FALSE.
OUTPUT_QUANTITY(231)%PART_APPROPRIATE = .FALSE.
OUTPUT_QUANTITY(231)%SLCF_APPROPRIATE = .FALSE.

! Boundary Quantities (Negative indices)
 
OUTPUT_QUANTITY(-1)%NAME = 'RADIATIVE HEAT FLUX'          
OUTPUT_QUANTITY(-1)%UNITS = 'kW/m2'                   
OUTPUT_QUANTITY(-1)%SHORT_NAME = 'rad'

OUTPUT_QUANTITY(-2)%NAME = 'CONVECTIVE HEAT FLUX'         
OUTPUT_QUANTITY(-2)%UNITS = 'kW/m2'                   
OUTPUT_QUANTITY(-2)%SHORT_NAME = 'con'

OUTPUT_QUANTITY(-3)%NAME = 'NORMAL VELOCITY'         
OUTPUT_QUANTITY(-3)%UNITS = 'm/s'                     
OUTPUT_QUANTITY(-3)%SHORT_NAME = 'vel'

OUTPUT_QUANTITY(-4)%NAME = 'GAS TEMPERATURE'         
OUTPUT_QUANTITY(-4)%UNITS = 'C'                       
OUTPUT_QUANTITY(-4)%SHORT_NAME = 'temp'

OUTPUT_QUANTITY(-5)%NAME = 'WALL TEMPERATURE'        
OUTPUT_QUANTITY(-5)%UNITS = 'C'                       
OUTPUT_QUANTITY(-5)%SHORT_NAME = 'temp'

OUTPUT_QUANTITY(-6)%NAME = 'INSIDE WALL TEMPERATURE' 
OUTPUT_QUANTITY(-6)%UNITS = 'C'                       
OUTPUT_QUANTITY(-6)%SHORT_NAME = 'inside'
OUTPUT_QUANTITY(-6)%INSIDE_SOLID = .TRUE.
OUTPUT_QUANTITY(-6)%BNDF_APPROPRIATE = .FALSE.

OUTPUT_QUANTITY(-7)%NAME = 'BURNING RATE'            
OUTPUT_QUANTITY(-7)%UNITS = 'kg/m2/s'                 
OUTPUT_QUANTITY(-7)%SHORT_NAME = 'burn'

OUTPUT_QUANTITY(-10)%NAME= 'NET HEAT FLUX'               
OUTPUT_QUANTITY(-10)%UNITS= 'kW/m2'                   
OUTPUT_QUANTITY(-10)%SHORT_NAME= 'net'

OUTPUT_QUANTITY(-11)%NAME= 'PRESSURE COEFFICIENT'    
OUTPUT_QUANTITY(-11)%UNITS= ' '                       
OUTPUT_QUANTITY(-11)%SHORT_NAME= 'c_p'

OUTPUT_QUANTITY(-12)%NAME= 'BACK WALL TEMPERATURE'   
OUTPUT_QUANTITY(-12)%UNITS= 'C'                       
OUTPUT_QUANTITY(-12)%SHORT_NAME= 'back'

OUTPUT_QUANTITY(-13)%NAME= 'GAUGE HEAT FLUX'         
OUTPUT_QUANTITY(-13)%UNITS= 'kW/m2'                   
OUTPUT_QUANTITY(-13)%SHORT_NAME= 'gauge'

OUTPUT_QUANTITY(-19)%NAME= 'INCIDENT HEAT FLUX'      
OUTPUT_QUANTITY(-19)%UNITS= 'kW/m2'                   
OUTPUT_QUANTITY(-19)%SHORT_NAME= 'in_flux'
 
OUTPUT_QUANTITY(-20)%NAME = 'MASS FLUX'
OUTPUT_QUANTITY(-20)%UNITS= 'kg/s/m2'
OUTPUT_QUANTITY(-20)%SHORT_NAME = 'mass_flux'
OUTPUT_QUANTITY(-20)%SPEC_ID_REQUIRED = .TRUE.

OUTPUT_QUANTITY(-21)%NAME= 'HEAT TRANSFER COEFFICIENT'    
OUTPUT_QUANTITY(-21)%UNITS= 'W/m2/K'                  
OUTPUT_QUANTITY(-21)%SHORT_NAME= 'h'

OUTPUT_QUANTITY(-22)%NAME= 'RADIOMETER'    
OUTPUT_QUANTITY(-22)%UNITS= 'kW/m2'                  
OUTPUT_QUANTITY(-22)%SHORT_NAME= 'radio'

OUTPUT_QUANTITY(-23)%NAME= 'ADIABATIC SURFACE TEMPERATURE'    
OUTPUT_QUANTITY(-23)%UNITS= 'C'                  
OUTPUT_QUANTITY(-23)%SHORT_NAME= 'AST'

OUTPUT_QUANTITY(-24)%NAME= 'WALL THICKNESS'    
OUTPUT_QUANTITY(-24)%UNITS= 'm'      
OUTPUT_QUANTITY(-24)%SHORT_NAME= 'thick'

OUTPUT_QUANTITY(-25)%NAME= 'SURFACE DENSITY'
OUTPUT_QUANTITY(-25)%UNITS= 'kg/m2'
OUTPUT_QUANTITY(-25)%SHORT_NAME= 'dens'

OUTPUT_QUANTITY(-26)%NAME= 'SOLID DENSITY'
OUTPUT_QUANTITY(-26)%UNITS= 'kg/m3'
OUTPUT_QUANTITY(-26)%SHORT_NAME= 'rho_s'
OUTPUT_QUANTITY(-26)%MATL_ID_REQUIRED = .TRUE.
OUTPUT_QUANTITY(-26)%INSIDE_SOLID = .TRUE.
OUTPUT_QUANTITY(-26)%BNDF_APPROPRIATE = .FALSE.

OUTPUT_QUANTITY(-27)%NAME= 'EMISSIVITY'
OUTPUT_QUANTITY(-27)%UNITS= ''
OUTPUT_QUANTITY(-27)%SHORT_NAME= 'emiss'

OUTPUT_QUANTITY(-28)%NAME= 'SOOT SURFACE DENSITY'
OUTPUT_QUANTITY(-28)%UNITS= 'kg/m2'
OUTPUT_QUANTITY(-28)%SHORT_NAME= 'soot den'

! Solid Phase Particle Outputs

OUTPUT_QUANTITY(-30)%NAME = 'MPUA'
OUTPUT_QUANTITY(-30)%UNITS = 'kg/m2'
OUTPUT_QUANTITY(-30)%SHORT_NAME = 'mpua'

OUTPUT_QUANTITY(-31)%NAME = 'CPUA'
OUTPUT_QUANTITY(-31)%UNITS = 'kW/m2'
OUTPUT_QUANTITY(-31)%SHORT_NAME = 'cpua'

OUTPUT_QUANTITY(-32)%NAME = 'AMPUA'
OUTPUT_QUANTITY(-32)%UNITS= 'kg/m2'
OUTPUT_QUANTITY(-32)%SHORT_NAME = 'ampua'

OUTPUT_QUANTITY(-32:-30)%PART_ID_REQUIRED = .TRUE.

END SUBROUTINE DEFINE_OUTPUT_QUANTITIES


SUBROUTINE SET_QUANTITIES_AMBIENT
 
! Define OUTPUT_QUANTITYs that have fixed names

OUTPUT_QUANTITY(1)%AMBIENT_VALUE = RHOA
OUTPUT_QUANTITY(5)%AMBIENT_VALUE = TMPA-TMPM
OUTPUT_QUANTITY(110)%AMBIENT_VALUE = TMPA-TMPM
OUTPUT_QUANTITY(-6:-4)%AMBIENT_VALUE = TMPA-TMPM
OUTPUT_QUANTITY(-12)%AMBIENT_VALUE = TMPA-TMPM
OUTPUT_QUANTITY(-23)%AMBIENT_VALUE = TMPA-TMPM

END SUBROUTINE SET_QUANTITIES_AMBIENT


SUBROUTINE SEARCH_KEYWORD(NAME,LU,IOS)
 
INTEGER, INTENT(OUT) :: IOS
INTEGER, INTENT(IN)  :: LU
CHARACTER(*), INTENT(IN) :: NAME
CHARACTER(40) :: TEXT
 
IF (LU<0) THEN
   IOS = -1
   RETURN
ENDIF
 
IOS = 1
REWIND(LU)
READLOOP: DO
   READ(LU,'(A)',END=10) TEXT
   IF (TRIM(TEXT)==TRIM(NAME)) THEN
      IOS = 0
      RETURN
   ELSE
      CYCLE READLOOP
   ENDIF
ENDDO READLOOP
 
10 RETURN
END SUBROUTINE SEARCH_KEYWORD


SUBROUTINE CHECK_SURF_NAME(NAME,EXISTS)

LOGICAL, INTENT(OUT) :: EXISTS
CHARACTER(*), INTENT(IN) :: NAME
INTEGER :: NS

EXISTS = .FALSE.
DO NS=0,N_SURF
   IF (NAME==SURFACE(NS)%ID) EXISTS = .TRUE.
ENDDO

END SUBROUTINE CHECK_SURF_NAME


SUBROUTINE GET_QUANTITY_INDEX(SMOKEVIEW_LABEL,SMOKEVIEW_BAR_LABEL,OUTPUT_INDEX,SPEC_INDEX,PART_INDEX,OUTTYPE, &
                              QUANTITY,SPEC_ID,PART_ID)

CHARACTER(*), INTENT(INOUT) :: QUANTITY
CHARACTER(*), INTENT(OUT) :: SMOKEVIEW_LABEL,SMOKEVIEW_BAR_LABEL
CHARACTER(*) :: SPEC_ID,PART_ID
CHARACTER(*), INTENT(IN) :: OUTTYPE
INTEGER, INTENT(OUT) :: OUTPUT_INDEX,SPEC_INDEX,PART_INDEX
INTEGER :: ND,NS

IF (QUANTITY=='OPTICAL DENSITY'        .AND. SPEC_ID=='null' .AND. MIXTURE_FRACTION) SPEC_ID='soot'
IF (QUANTITY=='EXTINCTION COEFFICIENT' .AND. SPEC_ID=='null' .AND. MIXTURE_FRACTION) SPEC_ID='soot'
IF (QUANTITY=='SOOT VOLUME FRACTION'   .AND. SPEC_ID=='null' .AND. MIXTURE_FRACTION) SPEC_ID='soot'
IF (QUANTITY=='VISIBILITY'             .AND. SPEC_ID=='null' .AND. MIXTURE_FRACTION) SPEC_ID='soot'

SPEC_INDEX = 0
PART_INDEX = 0

! Assign SPEC_INDEX when SPEC_ID is specified
IF (SPEC_ID/='null') THEN
   DO NS=1,N_SPECIES
      IF (SPEC_ID==SPECIES(NS)%ID) THEN
         IF (MIXTURE_FRACTION) THEN
            IF (SPECIES(NS)%INDEX <= N_STATE_SPECIES) THEN
               SPEC_INDEX = -SPECIES(NS)%INDEX
            ELSE
               SPEC_INDEX = NS + N_Y_ARRAY
            ENDIF
         ELSE
            SPEC_INDEX = NS
         ENDIF
         EXIT
      ENDIF
   ENDDO
ENDIF

IF (SPEC_ID/='null' .AND. SPEC_INDEX==0) THEN
   DO NS=1,N_STATE_SPECIES
      IF (SPEC_ID==MF_SPEC_ID(NS)) THEN
         SPEC_INDEX = -NS
         EXIT
      ENDIF
   ENDDO
ENDIF

! Convert old mixture fraction based quantities

IF (SPEC_INDEX==0 .AND. (QUANTITY=='fuel'.OR.QUANTITY=='fuel mass fraction')) THEN
   SPEC_INDEX = -FUEL_INDEX
   IF (QUANTITY=='fuel')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='fuel mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (SPEC_INDEX==0 .AND. (QUANTITY=='oxygen'.OR.QUANTITY=='oxygen mass fraction')) THEN
   SPEC_INDEX = -O2_INDEX
   IF (QUANTITY=='oxygen')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='oxygen mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (SPEC_INDEX==0 .AND. (QUANTITY=='nitrogen'.OR.QUANTITY=='nitrogen mass fraction')) THEN
   SPEC_INDEX = -N2_INDEX
   IF (QUANTITY=='nitrogen')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='nitrogen mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (SPEC_INDEX==0 .AND. (QUANTITY=='water vapor'.OR.QUANTITY=='water vapor mass fraction')) THEN
   SPEC_INDEX = -H2O_INDEX
   IF (QUANTITY=='water vapor')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='water vapor mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (SPEC_INDEX==0 .AND. (QUANTITY=='carbon dioxide'.OR.QUANTITY=='carbon dioxide mass fraction')) THEN
   SPEC_INDEX = -CO2_INDEX
   IF (QUANTITY=='carbon dioxide')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='carbon dioxide mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (SPEC_INDEX==0 .AND. (QUANTITY=='carbon monoxide'.OR.QUANTITY=='carbon monoxide mass fraction')) THEN
   SPEC_INDEX = -CO_INDEX
   IF (QUANTITY=='carbon monoxide')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='carbon monoxide mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (SPEC_INDEX==0 .AND. (QUANTITY=='hydrogen'.OR.QUANTITY=='hydrogen mass fraction')) THEN
   SPEC_INDEX = -H2_INDEX
   IF (QUANTITY=='hydrogen')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='hydrogen mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (SPEC_INDEX==0 .AND. (QUANTITY=='soot'.OR.QUANTITY=='soot mass fraction')) THEN
   SPEC_INDEX = -SOOT_INDEX
   IF (QUANTITY=='soot')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='soot mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (SPEC_INDEX==0 .AND. (QUANTITY=='other'.OR.QUANTITY=='other mass fraction')) THEN
   SPEC_INDEX = -OTHER_INDEX
   IF (QUANTITY=='other')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='other mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF

! Assign PART_INDEX when PART_ID is specified

IF (PART_ID/='null') THEN
   DO NS=1,N_PART
      IF (PART_ID==PARTICLE_CLASS(NS)%ID) THEN
         PART_INDEX = NS
         EXIT
      ENDIF
   ENDDO
ENDIF

! Convert outdated QUANTITY names to new format

IF (MIXTURE_FRACTION .AND. SPEC_ID=='null') THEN
   DO NS = I_Z_MIN, I_Z_MAX
      IF (QUANTITY == SPECIES(NS)%ID) THEN
         SPEC_ID = SPECIES(NS)%ID
         SPEC_INDEX = NS+N_Y_ARRAY
         QUANTITY = 'MASS FRACTION'
         EXIT
      ENDIF
   END DO
ENDIF

IF (SPEC_ID =='null') THEN
   DO NS=1,N_SPECIES
      IF (QUANTITY==TRIM(SPECIES(NS)%ID)) THEN
         QUANTITY   = 'MASS FRACTION'
         SPEC_ID = SPECIES(NS)%ID      
         SPEC_INDEX = SPECIES(NS)%INDEX
         IF (MIXTURE_FRACTION) THEN
            IF(SPECIES(NS)%INDEX > N_STATE_SPECIES) THEN
               SPEC_INDEX = NS + N_Y_ARRAY
            ELSE
               SPEC_INDEX = -SPECIES(NS)%INDEX
               SPEC_ID = MF_SPEC_ID(SPECIES(NS)%INDEX)
            ENDIF
         ENDIF
         EXIT
      ENDIF
      IF (QUANTITY==TRIM(SPECIES(NS)%ID)//'_VF') THEN
         QUANTITY   = 'VOLUME FRACTION'
         SPEC_ID = SPECIES(NS)%ID      
         SPEC_INDEX = SPECIES(NS)%INDEX
         IF (MIXTURE_FRACTION) THEN
            IF(SPECIES(NS)%INDEX > N_STATE_SPECIES) THEN
               SPEC_INDEX = NS + N_Y_ARRAY
            ELSE
               SPEC_INDEX = -SPECIES(NS)%INDEX
               SPEC_ID = MF_SPEC_ID(SPECIES(NS)%INDEX)
            ENDIF
         ENDIF
         EXIT
      ENDIF
      IF (QUANTITY==TRIM(SPECIES(NS)%ID)//'_FLUX_X') THEN
         QUANTITY   = 'MASS FLUX X'
         SPEC_ID = SPECIES(NS)%ID      
         SPEC_INDEX = SPECIES(NS)%INDEX
         IF (MIXTURE_FRACTION) THEN
            IF(SPECIES(NS)%INDEX > N_STATE_SPECIES) THEN
               SPEC_INDEX = NS + N_Y_ARRAY
            ELSE
               SPEC_INDEX = -SPECIES(NS)%INDEX
               SPEC_ID = MF_SPEC_ID(SPECIES(NS)%INDEX)
            ENDIF
         ENDIF
         EXIT
      ENDIF
      IF (QUANTITY==TRIM(SPECIES(NS)%ID)//'_FLUX_Y') THEN
         QUANTITY   = 'MASS FLUX Y'
         SPEC_ID = SPECIES(NS)%ID      
         SPEC_INDEX = SPECIES(NS)%INDEX
         IF (MIXTURE_FRACTION) THEN
            IF(SPECIES(NS)%INDEX > N_STATE_SPECIES) THEN
               SPEC_INDEX = NS + N_Y_ARRAY
            ELSE
               SPEC_INDEX = -SPECIES(NS)%INDEX
               SPEC_ID = MF_SPEC_ID(SPECIES(NS)%INDEX)
            ENDIF
         ENDIF
         EXIT
      ENDIF
      IF (QUANTITY==TRIM(SPECIES(NS)%ID)//'_FLUX_Z') THEN
         QUANTITY   = 'MASS FLUX Z'
         SPEC_ID = SPECIES(NS)%ID      
         SPEC_INDEX = SPECIES(NS)%INDEX
         IF (MIXTURE_FRACTION) THEN
            IF(SPECIES(NS)%INDEX > N_STATE_SPECIES) THEN
               SPEC_INDEX = NS + N_Y_ARRAY
            ELSE
               SPEC_INDEX = -SPECIES(NS)%INDEX
               SPEC_ID = MF_SPEC_ID(SPECIES(NS)%INDEX)
            ENDIF
         ENDIF
         EXIT
      ENDIF
      IF (QUANTITY==TRIM(SPECIES(NS)%ID)//'_FLUX') THEN
         QUANTITY   = 'MASS FLUX'
         SPEC_ID = SPECIES(NS)%ID      
         SPEC_INDEX = SPECIES(NS)%INDEX
         IF (MIXTURE_FRACTION) THEN
            IF(SPECIES(NS)%INDEX > N_STATE_SPECIES) THEN
               SPEC_INDEX = NS + N_Y_ARRAY
            ELSE
               SPEC_INDEX = -SPECIES(NS)%INDEX
               SPEC_ID = MF_SPEC_ID(SPECIES(NS)%INDEX)
            ENDIF
         ENDIF
         EXIT
      ENDIF
   ENDDO
ENDIF
      
DO NS=1,N_PART
   IF (QUANTITY==TRIM(PARTICLE_CLASS(NS)%ID)//'_MPUV') THEN
      QUANTITY   = 'MPUV'
      PART_INDEX = NS
   ENDIF
   IF (QUANTITY==TRIM(PARTICLE_CLASS(NS)%ID)//'_ADD') THEN
      QUANTITY   = 'ADD'
      PART_INDEX = NS
   ENDIF
   IF (QUANTITY==TRIM(PARTICLE_CLASS(NS)%ID)//'_ADT') THEN
      QUANTITY   = 'ADT'
      PART_INDEX = NS
   ENDIF
   IF (QUANTITY==TRIM(PARTICLE_CLASS(NS)%ID)//'_FLUX_X') THEN
      QUANTITY   = 'DROPLET FLUX X'
      PART_INDEX = NS
   ENDIF
   IF (QUANTITY==TRIM(PARTICLE_CLASS(NS)%ID)//'_FLUX_Y') THEN
      QUANTITY   = 'DROPLET FLUX Y'
      PART_INDEX = NS
   ENDIF
   IF (QUANTITY==TRIM(PARTICLE_CLASS(NS)%ID)//'_FLUX_Z') THEN
      QUANTITY   = 'DROPLET FLUX Z'
      PART_INDEX = NS
   ENDIF
   IF (QUANTITY==TRIM(PARTICLE_CLASS(NS)%ID)//'_MPUA') THEN
      QUANTITY   = 'MPUA'
      PART_INDEX = NS
   ENDIF
   IF (QUANTITY==TRIM(PARTICLE_CLASS(NS)%ID)//'_CPUA') THEN
      QUANTITY   = 'CPUA'
      PART_INDEX = NS
   ENDIF
   IF (QUANTITY==TRIM(PARTICLE_CLASS(NS)%ID)//'_AMPUA') THEN
      QUANTITY   = 'AMPUA'
      PART_INDEX = NS
   ENDIF
ENDDO
! Convert various old mixture fraction based names

IF (QUANTITY=='visibility') THEN
   QUANTITY = 'VISIBILITY'
   SPEC_INDEX = -SOOT_INDEX
ENDIF
IF (QUANTITY=='extinction coefficient') THEN
   QUANTITY = 'EXTINCTION COEFFICIENT'
   SPEC_INDEX = -SOOT_INDEX
ENDIF
IF (QUANTITY=='optical depth') THEN
   QUANTITY = 'OPTICAL DENSITY'
   SPEC_INDEX = -SOOT_INDEX
ENDIF
IF (QUANTITY=='soot density') THEN
   QUANTITY = 'DENSITY'
   SPEC_ID  = 'soot'
   SPEC_INDEX = -SOOT_INDEX
ENDIF
IF (QUANTITY=='soot volume fraction') THEN
   QUANTITY = 'SOOT VOLUME FRACTION'
   SPEC_INDEX = -SOOT_INDEX
ENDIF

! Convert old names to new

IF (QUANTITY=='ABSORPTION_COEFFICIENT')        QUANTITY = 'ABSORPTION COEFFICIENT'
IF (QUANTITY=='ADIABATIC_SURFACE_TEMPERATURE') QUANTITY = 'ADIABATIC SURFACE TEMPERATURE'
IF (QUANTITY=='aspiration')                    QUANTITY = 'ASPIRATION'
IF (QUANTITY=='BACK_WALL_TEMPERATURE')         QUANTITY = 'BACK WALL TEMPERATURE'
IF (QUANTITY=='BURNING_RATE')                  QUANTITY = 'BURNING RATE'
IF (QUANTITY=='CONVECTIVE_FLUX')               QUANTITY = 'CONVECTIVE HEAT FLUX'
IF (QUANTITY=='DROPLET_FLUX_X')                QUANTITY = 'DROPLET FLUX X'
IF (QUANTITY=='DROPLET_FLUX_Y')                QUANTITY = 'DROPLET FLUX Y'
IF (QUANTITY=='DROPLET_FLUX_Z')                QUANTITY = 'DROPLET FLUX Z'
IF (QUANTITY=='DROPLET_AGE')                   QUANTITY = 'DROPLET AGE'
IF (QUANTITY=='DROPLET_DIAMETER')              QUANTITY = 'DROPLET DIAMETER'
IF (QUANTITY=='DROPLET_TEMPERATURE')           QUANTITY = 'DROPLET TEMPERATURE'
IF (QUANTITY=='VEG_TEMPERATURE')               QUANTITY = 'DROPLET TEMPERATURE'
IF (QUANTITY=='DROPLET_VELOCITY')              QUANTITY = 'DROPLET VELOCITY'
IF (QUANTITY=='DROPLET_MASS')                  QUANTITY = 'DROPLET MASS'
IF (QUANTITY=='GAUGE_HEAT_FLUX')               QUANTITY = 'GAUGE HEAT FLUX'
IF (QUANTITY=='HEAT_FLUX')                     QUANTITY = 'NET HEAT FLUX'
IF (QUANTITY=='INCIDENT_HEAT_FLUX')            QUANTITY = 'INCIDENT HEAT FLUX'
IF (QUANTITY=='INSIDE_WALL_TEMPERATURE')       QUANTITY = 'INSIDE WALL TEMPERATURE'
IF (QUANTITY=='MIXTURE_FRACTION')              QUANTITY = 'MIXTURE FRACTION'
IF (QUANTITY=='path obscuration')              QUANTITY = 'PATH OBSCURATION'
IF (QUANTITY=='PRESSURE_COEFFICIENT')          QUANTITY = 'PRESSURE COEFFICIENT'
IF (QUANTITY=='spot obscuration')              QUANTITY = 'CHAMBER OBSCURATION'
IF (QUANTITY=='RADIATIVE_FLUX')                QUANTITY = 'RADIATIVE HEAT FLUX'
IF (QUANTITY=='RADIATIVE_FLUX_GAS')            QUANTITY = 'RADIATIVE HEAT FLUX GAS'
IF (QUANTITY=='WALL_TEMPERATURE')              QUANTITY = 'WALL TEMPERATURE'
IF (QUANTITY=='WALL_THICKNESS')                QUANTITY = 'WALL THICKNESS'

! Loop over all possible output quantities and assign an index number to match the desired QUANTITY

DO ND=-N_OUTPUT_QUANTITIES,N_OUTPUT_QUANTITIES

  IF (QUANTITY==OUTPUT_QUANTITY(ND)%NAME) THEN

     OUTPUT_INDEX = ND

     IF (OUTPUT_QUANTITY(ND)%MIXTURE_FRACTION_ONLY .AND. .NOT.MIXTURE_FRACTION) THEN
        WRITE(MESSAGE,'(5A)')  'ERROR: ',TRIM(OUTTYPE),' QUANTITY ',TRIM(QUANTITY), &
                               ' can only be used when the MIXTURE_FRACTION model is active'
        CALL SHUTDOWN(MESSAGE)
     ENDIF

     IF (SPEC_INDEX<0 .AND. .NOT.MIXTURE_FRACTION) THEN
        WRITE(MESSAGE,'(5A)')  'ERROR: ',TRIM(OUTTYPE),' QUANTITY ',TRIM(QUANTITY), &
                               ' can only be used when the MIXTURE_FRACTION model is active'
        CALL SHUTDOWN(MESSAGE)
     ENDIF

     IF (OUTPUT_QUANTITY(ND)%SPEC_ID_REQUIRED .AND. SPEC_INDEX==0) THEN
        WRITE(MESSAGE,'(3A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),' requires a SPEC_ID'
        CALL SHUTDOWN(MESSAGE)
     ENDIF

     IF (OUTPUT_QUANTITY(ND)%PART_ID_REQUIRED .AND. PART_INDEX<1) THEN
        WRITE(MESSAGE,'(3A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),' requires a PART_ID'
        CALL SHUTDOWN(MESSAGE)
     ENDIF

     IF (QUANTITY=='RELATIVE HUMIDITY' .AND. .NOT. MIXTURE_FRACTION .AND. I_WATER==0) THEN
        WRITE(MESSAGE,'(A)')  'ERROR: RELATIVE HUMIDTY requires water droplets, SPEC=WATER VAPOR, or MIXTURE_FRACTION'
        CALL SHUTDOWN(MESSAGE)
     END IF

     SELECT CASE (TRIM(OUTTYPE))
        CASE ('SLCF')
              ! Throw out bad slices
           IF (.NOT. OUTPUT_QUANTITY(ND)%SLCF_APPROPRIATE) THEN
              WRITE(MESSAGE,'(3A)')  ' ERROR: The QUANTITY ',TRIM(QUANTITY),' is not appropriate for SLCF'
              CALL SHUTDOWN(MESSAGE)
           ENDIF
        CASE ('DEVC')
            IF (.NOT.OUTPUT_QUANTITY(ND)%DEVC_APPROPRIATE) THEN
               WRITE(MESSAGE,'(3A)')  ' ERROR: The QUANTITY ',TRIM(QUANTITY),' is not appropriate for DEVC'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (QUANTITY=='AMPUA') ACCUMULATE_WATER = .TRUE.
        CASE ('PART')
           IF (.NOT. OUTPUT_QUANTITY(ND)%PART_APPROPRIATE) THEN
              WRITE(MESSAGE,'(3A)') 'ERROR: ',TRIM(QUANTITY),' is not a particle output QUANTITY'
              CALL SHUTDOWN(MESSAGE)
           ENDIF
        CASE ('BNDF')
           IF (.NOT. OUTPUT_QUANTITY(ND)%BNDF_APPROPRIATE) THEN
              WRITE(MESSAGE,'(3A)')  ' ERROR: The QUANTITY ',TRIM(QUANTITY),' is not appropriate for BNDF'
              CALL SHUTDOWN(MESSAGE)
           ENDIF
           IF (QUANTITY=='AMPUA') ACCUMULATE_WATER = .TRUE.
        CASE('ISOF')
           IF (.NOT.OUTPUT_QUANTITY(ND)%ISOF_APPROPRIATE) THEN
              WRITE(MESSAGE,'(3A)')  'ERROR: ISOF quantity ',TRIM(QUANTITY),' not appropriate for isosurface'
              CALL SHUTDOWN(MESSAGE)
           ENDIF
        CASE ('PLOT3D')
            IF (OUTPUT_QUANTITY(ND)%SOLID_PHASE) THEN
               WRITE(MESSAGE,'(5A)') 'ERROR: ',TRIM(OUTTYPE),'_QUANTITY ',TRIM(QUANTITY), ' not appropriate for gas phase'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (.NOT.OUTPUT_QUANTITY(ND)%SLCF_APPROPRIATE) THEN
               WRITE(MESSAGE,'(5A)') 'ERROR: ',TRIM(OUTTYPE),'_QUANTITY ',TRIM(QUANTITY), ' not appropriate for Plot3D'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
        CASE ('SMOKE3D')
            IF (SMOKE3D .AND. .NOT.OUTPUT_QUANTITY(ND)%MASS_FRACTION) THEN
               WRITE(MESSAGE,'(5A)') 'ERROR: ',TRIM(OUTTYPE),'_QUANTITY ',TRIM(QUANTITY), ' must be a mass fraction'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
        CASE DEFAULT
     END SELECT

     ! Assign Smokeview Label

     IF (SPEC_INDEX>0) THEN
        IF (MIXTURE_FRACTION .AND. SPEC_INDEX <= N_STATE_SPECIES) THEN
           SMOKEVIEW_LABEL = TRIM(MF_SPEC_ID(SPEC_INDEX))//' '//TRIM(QUANTITY)
           SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)//'_'//TRIM(SPECIES(SPEC_INDEX)%FORMULA)
        ELSEIF (MIXTURE_FRACTION .AND. SPEC_INDEX > N_Y_ARRAY) THEN   
           SPEC_INDEX = SPEC_INDEX - N_Y_ARRAY        
           SMOKEVIEW_LABEL = TRIM(SPECIES(SPEC_INDEX)%ID)//' '//TRIM(QUANTITY)
           SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)//'_'//TRIM(SPECIES(SPEC_INDEX)%FORMULA)
        ELSE
           DO NS = 1, N_SPECIES
              IF (SPECIES(NS)%INDEX == SPEC_INDEX) THEN
                 SMOKEVIEW_LABEL = TRIM(SPECIES(NS)%ID)//' '//TRIM(QUANTITY)
                 SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)//'_'//TRIM(SPECIES(NS)%FORMULA)
              ENDIF
           ENDDO
        ENDIF
     ELSEIF (SPEC_INDEX<0) THEN
        SMOKEVIEW_LABEL = TRIM(MF_SPEC_ID(-SPEC_INDEX))//' '//TRIM(QUANTITY)
        SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)//'_'//TRIM(MF_SPEC_FORMULA(-SPEC_INDEX))
      ELSEIF (PART_INDEX>0) THEN
        SMOKEVIEW_LABEL = TRIM(PARTICLE_CLASS(PART_INDEX)%ID)//' '//TRIM(QUANTITY)
        SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)
     ELSE
        SMOKEVIEW_LABEL = TRIM(QUANTITY)
        SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)
     ENDIF

     RETURN
  ENDIF
ENDDO

! If no match for desired QUANTITY is found, stop the job

WRITE(MESSAGE,'(5A)') 'ERROR: ',TRIM(OUTTYPE),' QUANTITY ',TRIM(QUANTITY), ' not found'
CALL SHUTDOWN(MESSAGE)

END SUBROUTINE GET_QUANTITY_INDEX



SUBROUTINE GET_PROPERTY_INDEX(P_INDEX,OUTTYPE,PROP_ID)

USE DEVICE_VARIABLES
CHARACTER(*), INTENT(IN) :: PROP_ID
CHARACTER(*), INTENT(IN) :: OUTTYPE
INTEGER, INTENT(INOUT) :: P_INDEX
INTEGER :: NN

DO NN=1,N_PROP
  IF (PROP_ID==PROPERTY(NN)%ID) THEN
     P_INDEX = NN
     SELECT CASE (TRIM(OUTTYPE))
        CASE ('SLCF')
        CASE ('DEVC')
        CASE ('PART')
        CASE ('BNDF')
        CASE ('PLOT3D')
        CASE DEFAULT
     END SELECT
     RETURN
  ENDIF
ENDDO

WRITE(MESSAGE,'(5A)')  ' ERROR: ',TRIM(OUTTYPE),' PROP_ID ',TRIM(PROP_ID),' not found'
CALL SHUTDOWN(MESSAGE)

END SUBROUTINE GET_PROPERTY_INDEX


SUBROUTINE SEARCH_CONTROLLER(NAME,CTRL_ID,DEVC_ID,DEVICE_INDEX,CONTROL_INDEX,INPUT_INDEX)

USE DEVICE_VARIABLES, ONLY: DEVICE,N_DEVC
USE CONTROL_VARIABLES, ONLY: CONTROL,N_CTRL
CHARACTER(*), INTENT (IN) :: NAME,CTRL_ID,DEVC_ID
INTEGER :: I, DEVICE_INDEX,CONTROL_INDEX
INTEGER , INTENT(IN) :: INPUT_INDEX
 
! There cannot be both a device and controller for any given entity

IF (DEVC_ID /= 'null' .AND. CTRL_ID /='null') THEN
   WRITE(MESSAGE,'(A,A,1X,I3,A)')  'ERROR: ',TRIM(NAME),INPUT_INDEX,' has both a device (DEVC) and a control (CTRL) specified'
   CALL SHUTDOWN(MESSAGE)
ENDIF

! Search for device

IF (DEVC_ID /= 'null') THEN
   DO I=1,N_DEVC
      IF (DEVICE(I)%ID==DEVC_ID) THEN
         DEVICE_INDEX = I
         RETURN
      ENDIF
   ENDDO
   WRITE(MESSAGE,'(A,A,A)')  'ERROR: DEVICE ',TRIM(DEVC_ID),' does not exist'
   CALL SHUTDOWN(MESSAGE)
ENDIF

! Search for controller

IF (CTRL_ID /= 'null') THEN
   DO I=1,N_CTRL
      IF (CONTROL(I)%ID==CTRL_ID) THEN
         CONTROL_INDEX = I
         RETURN
      ENDIF
   ENDDO
   WRITE(MESSAGE,'(A,A,A)')  'ERROR: CONTROL ',TRIM(CTRL_ID),' does not exist'
   CALL SHUTDOWN(MESSAGE)
ENDIF

END SUBROUTINE SEARCH_CONTROLLER
 

SUBROUTINE GET_REV_read(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') readrev(INDEX(readrev,':')+1:LEN_TRIM(readrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') readdate

END SUBROUTINE GET_REV_read
 
END MODULE READ_INPUT

