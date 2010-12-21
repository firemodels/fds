MODULE READ_INPUT
 
USE PRECISION_PARAMETERS
USE MESH_VARIABLES
USE GLOBAL_CONSTANTS
USE TRAN
USE MESH_POINTERS
USE OUTPUT_DATA
USE COMP_FUNCTIONS, ONLY: SECOND, CHECKREAD, SHUTDOWN, CHECK_XB
USE MEMORY_FUNCTIONS, ONLY: ChkMemErr
USE COMP_FUNCTIONS, ONLY: GET_INPUT_FILE
USE MISC_FUNCTIONS, ONLY: SEARCH_CONTROLLER
USE EVAC, ONLY: READ_EVAC
USE HVAC_ROUTINES, ONLY: READ_HVAC
 
IMPLICIT NONE
PRIVATE
CHARACTER(255), PARAMETER :: readid='$Id$'
CHARACTER(255), PARAMETER :: readrev='$Revision$'
CHARACTER(255), PARAMETER :: readdate='$Date$'

PUBLIC READ_DATA, GET_REV_read

CHARACTER(30) :: LABEL,MB
CHARACTER(100) :: MESSAGE,FYI
CHARACTER(30) :: ID,SURF_DEFAULT,EVAC_SURF_DEFAULT,SPECIES_ID(1:MAX_SPECIES)='null'
CHARACTER(30), ALLOCATABLE, DIMENSION(:,:) :: MATL_SPEC_ID
LOGICAL :: SUCCESS,EX,THICKEN_OBSTRUCTIONS,BAD,DEFINED_BACKGROUND=.FALSE.
REAL(EB) :: XB(6),TEXTURE_ORIGIN(3)
REAL(EB) :: PBX,PBY,PBZ,HUMIDITY
REAL(EB) :: MW_MIN,MW_MAX
INTEGER  :: I,J,K,IZERO,IOS,N_MIX_READ=0
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
      IF (.NOT.USE_MPI .AND. .NOT.USE_OPENMP)   &
         WRITE(LU_ERR,'(/A,A,A)') "Version: ",TRIM(VERSION_STRING),"; MPI Disabled; OpenMP Disabled"
      IF (     USE_MPI .AND. .NOT.USE_OPENMP)   &
         WRITE(LU_ERR,'(/A,A,A)') "Version: ",TRIM(VERSION_STRING),"; MPI Enabled; OpenMP Disabled"
      IF (.NOT.USE_MPI .AND.      USE_OPENMP)   &
         WRITE(LU_ERR,'(/A,A,A)') "Version: ",TRIM(VERSION_STRING),"; MPI Disabled; OpenMP Enabled"
      IF (     USE_MPI .AND.      USE_OPENMP)   &
         WRITE(LU_ERR,'(/A,A,A)') "Version: ",TRIM(VERSION_STRING),"; MPI Enabled; OpenMP Enabled"
   ENDIF
   IF (USE_OPENMP .and. .NOT.USE_MPI) &
      WRITE(LU_ERR,'(A,I3)') 'Number of available OpenMP threads: ',OPENMP_AVAILABLE_THREADS
   IF (MYID==0) THEN
      WRITE(LU_ERR,'(A,I4)') "SVN Revision Number: ",SVN_REVISION_NUMBER
      WRITE(LU_ERR,'(A,A)') "Compile Date: ",TRIM(COMPILE_DATE)
      WRITE(LU_ERR,'(/A)')  "Consult FDS Users Guide Chapter, Running FDS, for further instructions."
      WRITE(LU_ERR,'(/A)')  "Hit Enter to Escape..."
      READ(5,*,ERR=2,END=2)
   ENDIF
 2 STOP
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
CALL READ_MISC
CALL READ_MULT
CALL READ_EVAC(1)    ! Read some evacuation input from the main evacuation meshes
CALL READ_MESH
CALL READ_TRAN
CALL READ_TIME
CALL READ_PRES
CALL READ_REAC_1
CALL READ_SPEC
CALL READ_SMIX
CALL PROC_SPEC    ! Set up various SPECies constructs
CALL READ_REAC_2
CALL READ_RADI
CALL READ_PROP
CALL READ_PART
CALL READ_DEVC
CALL READ_CTRL
CALL READ_GEOM    ! Experimental complex geometry ~RJM
CALL READ_TREE
CALL READ_MATL
CALL READ_SURF
CALL READ_OBST
CALL READ_VENT
CALL READ_ZONE
CALL READ_EVAC(2)
CALL READ_HVAC
CALL PROC_PART    ! Set up various PARTicle constructs
CALL PROC_SURF_1  ! Set up SURFace constructs for species
CALL READ_RAMP    ! Read in all RAMPs, assuming they have all been identified previously
CALL READ_TABL   ! Read in all TABLs, assuming they have all been identified previously
CALL PROC_MATL    ! Set up various MATeriaL constructs
CALL PROC_SURF_2  ! Set up remaining SURFace constructs
CALL READ_DUMP
CALL READ_CLIP
CALL READ_INIT
CALL PROC_WALL    ! Set up grid for 1-D heat transfer in solids
CALL PROC_CTRL    ! Set up various ConTRoL constructs
CALL PROC_PROP    ! Set up various PROPerty constructs
CALL PROC_DEVC    ! Set up various DEViCe constructs
CALL PROC_OBST    ! if an obst has a PROP then make sure PROP exists
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
INTEGER :: NAMELENGTH
NAMELIST /HEAD/ TITLE,CHID,FYI
 
CHID    = 'null'
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

IF (TRIM(CHID)=='null') THEN
   NAMELENGTH = LEN_TRIM(FN_INPUT)
   ROOTNAME: DO I=NAMELENGTH,2,-1      
      IF (FN_INPUT(I:I)=='.') THEN
         WRITE(CHID,'(A)') FN_INPUT(1:I-1)
         EXIT ROOTNAME
      ENDIF
   END DO ROOTNAME
ENDIF

! Define and look for a stop file

FN_STOP = TRIM(CHID)//'.stop'
INQUIRE(FILE=FN_STOP,EXIST=EX)
IF (EX) THEN
   WRITE(MESSAGE,'(A,A,A)') "ERROR: Remove the file, ",TRIM(FN_STOP),", from the current directory"
   CALL SHUTDOWN(MESSAGE)
ENDIF
 
END SUBROUTINE READ_HEAD
 
 
SUBROUTINE READ_MESH
USE EVAC, ONLY: N_DOORS, N_EXITS, N_CO_EXITS, EVAC_EMESH_EXITS_TYPE, EMESH_EXITS, EMESH_ID, EMESH_IJK, EMESH_XB, &
     EMESH_NM, N_DOOR_MESHES, EMESH_NFIELDS, EVAC_FDS6, HUMAN_SMOKE_HEIGHT, EVAC_DELTA_SEE, &
     EMESH_STAIRS, EVAC_EMESH_STAIRS_TYPE, N_STRS

INTEGER :: IJK(3),NM,CURRENT_MPI_PROCESS,MPI_PROCESS,RGB(3),LEVEL,N_MESH_NEW,N,II,JJ,KK,NMESHES_READ,NNN,NEVAC_MESHES
LOGICAL :: EVACUATION, EVAC_HUMANS
REAL(EB) :: EVAC_Z_OFFSET,XB1,XB2,XB3,XB4,XB5,XB6
CHARACTER(25) :: COLOR
CHARACTER(30) :: MULT_ID
NAMELIST /MESH/ IJK,FYI,ID,SYNCHRONIZE,EVACUATION,EVAC_HUMANS,CYLINDRICAL,XB,RGB,COLOR,EVAC_Z_OFFSET, &
                MPI_PROCESS,LEVEL,MULT_ID
TYPE (MESH_TYPE), POINTER :: M
TYPE (MULTIPLIER_TYPE), POINTER :: MR
 
NMESHES = 0
NMESHES_READ = 0
NEVAC_MESHES = 0
 
REWIND(LU_INPUT)
COUNT_MESH_LOOP: DO
   CALL CHECKREAD('MESH',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_MESH_LOOP
   MULT_ID = 'null'
   EVACUATION  = .FALSE.
   READ(LU_INPUT,MESH,END=15,ERR=16,IOSTAT=IOS)
   NMESHES_READ = NMESHES_READ + 1
   IF (NO_EVACUATION .AND. EVACUATION) CYCLE COUNT_MESH_LOOP ! skip evacuation meshes
   IF (EVACUATION_DRILL .AND. .NOT.EVACUATION) CYCLE COUNT_MESH_LOOP ! skip fire meshes
   IF (EVACUATION_MC_MODE .AND. .NOT.EVACUATION) CYCLE COUNT_MESH_LOOP ! skip fire meshes
   IF (EVACUATION .AND. EVAC_FDS6) NEVAC_MESHES = NEVAC_MESHES + 1
   N_MESH_NEW = 1
   IF (MULT_ID/='null') THEN
      DO N=1,N_MULT
         MR => MULTIPLIER(N)
         IF (MULT_ID==MR%ID) N_MESH_NEW = MR%N_COPIES
      ENDDO
   ENDIF
   NMESHES      = NMESHES + N_MESH_NEW
!   NMESHES_READ = NMESHES_READ + 1
   16 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with MESH line.')
ENDDO COUNT_MESH_LOOP
15 CONTINUE

IF (.NOT. NO_EVACUATION) NMESHES = NMESHES + N_DOOR_MESHES + NEVAC_MESHES
IF (.NOT. NO_EVACUATION .AND. EVAC_FDS6) NMESHES = NMESHES + N_STRS

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

IF (NMESHES<1) CALL SHUTDOWN('ERROR: No MESH line(s) defined.')

NM = 0
 
MESH_LOOP: DO N=1,NMESHES_READ

   ! Set MESH defaults

   IJK(1)= 10
   IJK(2)= 10
   IJK(3)= 10
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

   IF (NO_EVACUATION .AND. EVACUATION) CYCLE MESH_LOOP ! skip evacuation meshes
   IF (EVACUATION_DRILL .AND. .NOT.EVACUATION) CYCLE MESH_LOOP ! skip fire meshes
   IF (EVACUATION_MC_MODE .AND. .NOT.EVACUATION) CYCLE MESH_LOOP ! skip fire meshes

   ! Reorder XB coordinates if necessary

   CALL CHECK_XB(XB)

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
            ELSE
               IF (     USE_MPI) CURRENT_MPI_PROCESS = MIN(NM-1,NUMPROCS-1)
               IF (.NOT.USE_MPI) CURRENT_MPI_PROCESS = 0
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
            M%N_EXTERNAL_WALL_CELLS = 2*M%IBAR*M%JBAR+2*M%IBAR*M%KBAR+2*M%JBAR*M%KBAR
            IF (.NOT.SYNCHRONIZE) SYNC_TIME_STEP(NM)  = .FALSE.
            IF (EVACUATION)  EVACUATION_ONLY(NM) = .TRUE.
            IF (EVACUATION)  SYNC_TIME_STEP(NM)  = .FALSE.
            IF (EVAC_HUMANS) EVACUATION_GRID(NM) = .TRUE.
            IF (EVACUATION)  EVACUATION_Z_OFFSET(NM) = EVAC_Z_OFFSET
            IF (EVACUATION)  M%N_EXTERNAL_WALL_CELLS = 2*M%IBAR*M%KBAR+2*M%JBAR*M%KBAR
            IF (EVAC_FDS6 .AND. EVACUATION .AND. .NOT.EVAC_HUMANS) THEN
               WRITE(MESSAGE,'(A)') 'ERROR: NO DOOR FLOW EVACUATION MESHES IN FDS6'
               CALL SHUTDOWN(MESSAGE)
            ENDIF

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
            IF (MYID==0 .AND. USE_MPI) WRITE(LU_ERR,'(A,I3,A,I3)') 'Mesh ',NM,' is assigned to Process ',PROCESS(NM)
            IF (EVACUATION_ONLY(NM) .AND. USE_MPI) EVAC_PROCESS = NUMPROCS-1

            ! Mesh boundary colors
   
            IF (ANY(RGB<0) .AND. COLOR=='null') COLOR = 'BLACK'
            IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
            ALLOCATE(M%RGB(3))
            M%RGB = RGB
   
            ! Mesh Geometry and Name
   
            WRITE(MESH_NAME(NM),'(A,I3)') 'MESH',NM
            IF (ID/='null') MESH_NAME(NM) = ID
   
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
            IF (EVACUATION .AND. ABS(XB5 - XB6) <= ZERO_P) THEN
               WRITE(MESSAGE,'(A,I2)') 'ERROR: ZMIN = ZMAX on evacuation MESH ', NM
               CALL SHUTDOWN(MESSAGE)
            ENDIF

            M%XS    = XB1
            M%XF    = XB2
            M%YS    = XB3
            M%YF    = XB4
            M%ZS    = XB5
            M%ZF    = XB6
            IF (.NOT.EVACUATION ) THEN
               XS_MIN  = MIN(XS_MIN,M%XS)
               XF_MAX  = MAX(XF_MAX,M%XF)
               YS_MIN  = MIN(YS_MIN,M%YS)
               YF_MAX  = MAX(YF_MAX,M%YF)
               ZS_MIN  = MIN(ZS_MIN,M%ZS)
               ZF_MAX  = MAX(ZF_MAX,M%ZF)
            ENDIF
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

! Min and Max values of temperature
 
TMPMIN = TMPM
IF (LAPSE_RATE < 0._EB) TMPMIN = MIN(TMPMIN,TMPA+LAPSE_RATE*ZF_MAX)
TMPMAX = 3000._EB

REWIND(LU_INPUT)

! Define the additional evacuation door flow meshes

!Timo: Mesh counter NM is now fire meshes plus main evac meshes
IF (.NOT. NO_EVACUATION) CALL DEFINE_EVACUATION_MESHES(NM)

! Start the timing arrays
 
TUSED      = 0._EB
TUSED(1,:) = SECOND()

CONTAINS

  SUBROUTINE DEFINE_EVACUATION_MESHES(NM)
    IMPLICIT NONE
    ! Passed variables
    INTEGER, INTENT(INOUT) :: NM
    ! Local variables
    INTEGER :: N, N_END, I, J, NN, JMAX, NM_OLD, I_MAIN_EVAC_MESH
    REAL(EB) :: Z_MID

    N = 0
    DO I = 1, NM
       IF (EVACUATION_GRID(I) .AND. EVACUATION_ONLY(I)) THEN
          N = N + 1 ! Main evacuation mesh index for EMESH_EXITS(N) array
          EMESH_NM(N) = I
       END IF
    END DO

    NM_OLD = NM
    LOOP_EMESHES: DO N = 1, NEVAC_MESHES
       ! Additional meshes for the main evacuation meshes. These will be 
       ! at different z level than the corresponding main evacuation mesh.

       I_MAIN_EVAC_MESH = NM_OLD - NEVAC_MESHES + N

       ! Set MESH defaults

       RGB   = MESHES(I_MAIN_EVAC_MESH)%RGB
       COLOR = 'null'
       ID = TRIM(TRIM('Emesh_' // MESH_NAME(I_MAIN_EVAC_MESH)))
       MPI_PROCESS = -1
       LEVEL = 0
       EVACUATION = .TRUE.
       EVAC_HUMANS = .FALSE.

       ! Increase the MESH counter by 1

       NM = NM + 1

       ! Fill in MESH related variables
       
       M => MESHES(NM)
       M%MESH_LEVEL = LEVEL
       M%IBAR = MESHES(I_MAIN_EVAC_MESH)%IBAR
       M%JBAR = MESHES(I_MAIN_EVAC_MESH)%JBAR
       M%KBAR = MESHES(I_MAIN_EVAC_MESH)%KBAR
       IBAR_MAX = MAX(IBAR_MAX,M%IBAR)
       JBAR_MAX = MAX(JBAR_MAX,M%JBAR)
       KBAR_MAX = MAX(KBAR_MAX,M%KBAR)
       EVACUATION_ONLY(NM) = .TRUE.
       SYNC_TIME_STEP(NM)  = .FALSE.
       EVACUATION_GRID(NM) = .FALSE.
       EVACUATION_Z_OFFSET(NM) = EVAC_Z_OFFSET
       M%N_EXTERNAL_WALL_CELLS = 2*M%IBAR*M%KBAR+2*M%JBAR*M%KBAR
       IF (EVACUATION .AND. M%KBAR/=1) THEN
          WRITE(MESSAGE,'(A)') 'ERROR: IJK(3) must be 1 for all evacuation grids'
          CALL SHUTDOWN(MESSAGE)
       ENDIF

       ! Associate the MESH with the PROCESS
         
       PROCESS(NM) = CURRENT_MPI_PROCESS
       IF (MYID==0 .AND. USE_MPI) WRITE(LU_ERR,'(A,I3,A,I3)') 'Mesh ',NM,' is assigned to Process ',PROCESS(NM)
       IF (EVACUATION_ONLY(NM) .AND. USE_MPI) EVAC_PROCESS = NUMPROCS-1

       ! Mesh boundary colors
   
       IF (ANY(RGB<0) .AND. COLOR=='null') COLOR = 'BLACK'
       IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)
       ALLOCATE(M%RGB(3))
       M%RGB = RGB
   
       ! Mesh Geometry and Name
   
       WRITE(MESH_NAME(NM),'(A,I3)') 'MESH',NM
       IF (ID/='null') MESH_NAME(NM) = ID
   
       Z_MID = 0.5_EB*(MESHES(I_MAIN_EVAC_MESH)%ZS + MESHES(I_MAIN_EVAC_MESH)%ZF)
       Z_MID = Z_MID - EVACUATION_Z_OFFSET(I_MAIN_EVAC_MESH)  + HUMAN_SMOKE_HEIGHT
       M%XS    = MESHES(I_MAIN_EVAC_MESH)%XS
       M%XF    = MESHES(I_MAIN_EVAC_MESH)%XF
       M%YS    = MESHES(I_MAIN_EVAC_MESH)%YS
       M%YF    = MESHES(I_MAIN_EVAC_MESH)%YF
       M%ZS    = Z_MID - EVAC_DELTA_SEE
       M%ZF    = Z_MID + EVAC_DELTA_SEE
       M%DXI   = MESHES(I_MAIN_EVAC_MESH)%DXI
       M%DETA  = MESHES(I_MAIN_EVAC_MESH)%DETA
       M%DZETA = (M%ZF-M%ZS)/REAL(M%KBAR,EB)
       M%RDXI  = MESHES(I_MAIN_EVAC_MESH)%RDXI
       M%RDETA = MESHES(I_MAIN_EVAC_MESH)%RDETA
       M%RDZETA= 1._EB/M%DZETA
       M%IBM1  = M%IBAR-1
       M%JBM1  = M%JBAR-1
       M%KBM1  = M%KBAR-1
       M%IBP1  = M%IBAR+1
       M%JBP1  = M%JBAR+1
       M%KBP1  = M%KBAR+1
       ! WRITE (LU_ERR,FMT='(A,I5,3A)') ' EVAC: Mesh number ', NM, ' name ', TRIM(ID), ' defined for evacuation'
    END DO LOOP_EMESHES

    N_END = N_EXITS - N_CO_EXITS + N_DOORS
    LOOP_EXITS: DO N = 1, N_END
       I = EMESH_EXITS(N)%EMESH  ! The main evacuation mesh index (for EMESH_EXITS(I) array)
       IF (.NOT.EMESH_EXITS(N)%DEFINE_MESH) CYCLE LOOP_EXITS

       EMESH_EXITS(N)%MAINMESH = EMESH_NM(EMESH_EXITS(N)%EMESH)        ! The 1,...,NMESHES index
       ! Only main evacuation meshes in FDS6
       IF (EVAC_FDS6) EMESH_EXITS(N)%IMESH = EMESH_EXITS(N)%MAINMESH   ! The mesh index (all meshes included)

       ! Set MESH defaults

       IJK(1)= EMESH_IJK(1,I)
       IJK(2)= EMESH_IJK(2,I)
       IJK(3)= EMESH_IJK(3,I)

       ALLOCATE(EMESH_EXITS(N)%U_EVAC(0:IJK(1)+1,0:IJK(2)+1),STAT=IZERO)
       CALL ChkMemErr('READ','EMESH_EXITS(N)%U_EVAC',IZERO)
       ALLOCATE(EMESH_EXITS(N)%V_EVAC(0:IJK(1)+1,0:IJK(2)+1),STAT=IZERO)
       CALL ChkMemErr('READ','EMESH_EXITS(N)%V_EVAC',IZERO)

       IF (EVAC_FDS6) CYCLE LOOP_EXITS

       TWO_D = .FALSE.
       XB(1) = EMESH_XB(1,I)
       XB(2) = EMESH_XB(2,I)
       XB(3) = EMESH_XB(3,I)
       XB(4) = EMESH_XB(4,I)
       XB(5) = EMESH_XB(5,I)
       XB(6) = EMESH_XB(6,I)
       ID = TRIM(TRIM('Emesh_' // TRIM(EMESH_EXITS(N)%ID)))
       SYNCHRONIZE = .FALSE.
       EVACUATION  = .TRUE.
       EVAC_HUMANS = .FALSE.
       EVAC_Z_OFFSET = 1.0_EB
       MPI_PROCESS = -1
       LEVEL = 0
       MULT_ID = 'null'

       ! Increase the MESH counter by 1

       NM = NM + 1
       EMESH_EXITS(N)%IMESH = NM   ! The mesh index (all meshes included)

       ! Fill in MESH related variables
       
       M => MESHES(NM)
       M%MESH_LEVEL = LEVEL
       M%IBAR = IJK(1)
       M%JBAR = IJK(2)
       M%KBAR = IJK(3)
       IBAR_MAX = MAX(IBAR_MAX,M%IBAR)
       JBAR_MAX = MAX(JBAR_MAX,M%JBAR)
       KBAR_MAX = MAX(KBAR_MAX,M%KBAR)
       EVACUATION_ONLY(NM) = .TRUE.
       SYNC_TIME_STEP(NM)  = .FALSE.
       EVACUATION_GRID(NM) = .FALSE.
       EVACUATION_Z_OFFSET(NM) = EVAC_Z_OFFSET
       M%N_EXTERNAL_WALL_CELLS = 2*M%IBAR*M%KBAR+2*M%JBAR*M%KBAR
       TWO_D = .FALSE.
       IF (EVACUATION .AND. M%KBAR/=1) THEN
          WRITE(MESSAGE,'(A)') 'ERROR: IJK(3) must be 1 for all evacuation grids'
          CALL SHUTDOWN(MESSAGE)
       ENDIF

       ! Associate the MESH with the PROCESS
         
       PROCESS(NM) = CURRENT_MPI_PROCESS
       IF (MYID==0 .AND. USE_MPI) WRITE(LU_ERR,'(A,I3,A,I3)') 'Mesh ',NM,' is assigned to Process ',PROCESS(NM)
       IF (EVACUATION_ONLY(NM) .AND. USE_MPI) EVAC_PROCESS = NUMPROCS-1

       ! Mesh boundary colors
   
       ALLOCATE(M%RGB(3))
       M%RGB = EMESH_EXITS(N)%RGB
   
       ! Mesh Geometry and Name
   
       WRITE(MESH_NAME(NM),'(A,I3)') 'MESH',NM
       IF (ID/='null') MESH_NAME(NM) = ID
   
       M%XS    = XB(1)
       M%XF    = XB(2)
       M%YS    = XB(3)
       M%YF    = XB(4)
       M%ZS    = XB(5)
       M%ZF    = XB(6)
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
       WRITE (LU_ERR,FMT='(A,I5,3A)') ' EVAC: Mesh number ', NM, ' name ', TRIM(ID), ' defined for evacuation'
    END DO LOOP_EXITS

    NN = 0
    JMAX = 0
    DO I = 1, NM
       IF (EVACUATION_GRID(I) .AND. EVACUATION_ONLY(I)) THEN
          J = 0  ! Index of the flow field (for a main evacuation mesh)
          NN = NN + 1 ! Main evacuation mesh index
          ! NN = EMESH_INDEX(NM)
          EMESH_NFIELDS(NN) = 0  ! How many fields for this main evacuation mesh
          LOOP_EXITS_0: DO N = 1, N_END
             IF (.NOT.EMESH_EXITS(N)%DEFINE_MESH) CYCLE LOOP_EXITS_0
             IF (.NOT.EMESH_EXITS(N)%EMESH == NN) CYCLE LOOP_EXITS_0
             J = J + 1
             EMESH_EXITS(N)%I_DOORS_EMESH = J
             EMESH_NFIELDS(NN) = J
          END DO LOOP_EXITS_0
          WRITE(LU_ERR,FMT='(A,I5,3A,I5,A)') ' EVAC: Emesh ',NN,' ',TRIM(EMESH_ID(NN)),' has ',&
               EMESH_NFIELDS(NN),' door flow fields'
       END IF
    END DO

    IF (EVAC_FDS6) THEN
       ! Next line should be executed only once during a FDS+Evac run
       JMAX = MAXVAL(EMESH_NFIELDS)
       EVAC_TIME_ITERATIONS = EVAC_TIME_ITERATIONS*JMAX

       LOOP_STAIRS: DO N = 1, N_STRS
          ! Evacuation meshes for the stairs.

          ! Set MESH defaults

          RGB   = EMESH_STAIRS(N)%RGB
          COLOR = 'null'
          ID = TRIM('Emesh_' // TRIM(EMESH_STAIRS(N)%ID))
          MPI_PROCESS = -1
          LEVEL = 0
          EVACUATION = .TRUE.
          EVAC_HUMANS = .TRUE.

          ! Increase the MESH counter by 1

          NM = NM + 1
          EMESH_STAIRS(N)%IMESH = NM
          !EMESH_STAIRS(N)%EMESH = NM

          ! Fill in MESH related variables
       
          M => MESHES(NM)
          M%MESH_LEVEL = LEVEL
          M%IBAR = EMESH_STAIRS(N)%IBAR
          M%JBAR = EMESH_STAIRS(N)%JBAR
          M%KBAR = EMESH_STAIRS(N)%KBAR
          IBAR_MAX = MAX(IBAR_MAX,M%IBAR)
          JBAR_MAX = MAX(JBAR_MAX,M%JBAR)
          KBAR_MAX = MAX(KBAR_MAX,M%KBAR)
          EVACUATION_ONLY(NM) = .TRUE.
          SYNC_TIME_STEP(NM)  = .FALSE.
          EVACUATION_GRID(NM) = .TRUE.
          EVACUATION_Z_OFFSET(NM) = EVAC_Z_OFFSET
          M%N_EXTERNAL_WALL_CELLS = 2*M%IBAR*M%KBAR+2*M%JBAR*M%KBAR
          IF (EVACUATION .AND. M%KBAR/=1) THEN
             WRITE(MESSAGE,'(A)') 'ERROR: IJK(3) must be 1 for all evacuation grids'
             CALL SHUTDOWN(MESSAGE)
          ENDIF

          ! Associate the MESH with the PROCESS
         
          PROCESS(NM) = CURRENT_MPI_PROCESS
          IF (MYID==0 .AND. USE_MPI) WRITE(LU_ERR,'(A,I3,A,I3)') 'Mesh ',NM,' is assigned to Process ',PROCESS(NM)
          IF (EVACUATION_ONLY(NM) .AND. USE_MPI) EVAC_PROCESS = NUMPROCS-1

          ! Mesh boundary colors
   
          ALLOCATE(M%RGB(3))
          M%RGB = EMESH_STAIRS(N)%RGB
   
          ! Mesh Geometry and Name
   
          WRITE(MESH_NAME(NM),'(A,I3)') 'MESH',NM
          IF (ID/='null') MESH_NAME(NM) = ID
   
          M%XS    = EMESH_STAIRS(N)%XB(1)
          M%XF    = EMESH_STAIRS(N)%XB(2)
          M%YS    = EMESH_STAIRS(N)%XB(3)
          M%YF    = EMESH_STAIRS(N)%XB(4)
          M%ZS    = EMESH_STAIRS(N)%XB(5)
          M%ZF    = EMESH_STAIRS(N)%XB(6)
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
          WRITE (LU_ERR,FMT='(A,I5,3A)') ' EVAC: Mesh number ', NM, ' name ', TRIM(ID), ' defined for evacuation'
       END DO LOOP_STAIRS
    END IF

    IF (EVACUATION_ONLY(NM)) THEN
       DO N = 1, NMESHES
          M => MESHES(NM)
          XS_MIN  = MIN(XS_MIN,M%XS)
          XF_MAX  = MAX(XF_MAX,M%XF)
          YS_MIN  = MIN(YS_MIN,M%YS)
          YF_MAX  = MAX(YF_MAX,M%YF)
          ZS_MIN  = MIN(ZS_MIN,M%ZS)
          ZF_MAX  = MAX(ZF_MAX,M%ZF)
       END DO
    ENDIF

    RETURN
  END SUBROUTINE DEFINE_EVACUATION_MESHES

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
 
   NIPX   = 500*M%IBAR
   NIPY   = 500*M%JBAR
   NIPZ   = 500*M%KBAR
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
 
INTEGER FUNCTION IFAC(II,NN)
INTEGER, INTENT(IN) :: II,NN
INTEGER :: III
IFAC = 1
DO III=II-NN+1,II
   IFAC = IFAC*III
ENDDO
END FUNCTION IFAC

END SUBROUTINE READ_TRAN
 
 
SUBROUTINE READ_TIME
IMPLICIT NONE 
REAL(EB) :: DT,VEL_CHAR,TWFIN
INTEGER :: NM
NAMELIST /TIME/ DT,T_BEGIN,T_END,TWFIN,FYI,WALL_INCREMENT,SYNCHRONIZE, &
                EVAC_DT_FLOWFIELD,EVAC_DT_STEADY_STATE,TIME_SHRINK_FACTOR, &
                LOCK_TIME_STEP,RESTRICT_TIME_STEP
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
      VEL_CHAR = 0.2_EB*SQRT(10._EB*(ZF_MAX-ZS_MIN))
      M%DT = M%CELL_SIZE/VEL_CHAR
   ENDIF
   IF (EVACUATION_ONLY(NM)) THEN
      M%DT = EVAC_DT_FLOWFIELD
   ENDIF

   M%PART_CFL = 0.0

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
   IF (ABS(DX)>ZERO_P) MR%DXB(1:2) = DX
   IF (ABS(DY)>ZERO_P) MR%DXB(3:4) = DY
   IF (ABS(DZ)>ZERO_P) MR%DXB(5:6) = DZ

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
USE PHYSICAL_FUNCTIONS, ONLY : AMBIENT_WATER_VAPOR
USE MATH_FUNCTIONS, ONLY: GET_RAMP_INDEX
REAL(EB) :: C_HORIZONTAL,C_VERTICAL,FORCE_VECTOR(3)=0._EB,H_FIXED
CHARACTER(30) :: RAMP_GX,RAMP_GY,RAMP_GZ
LOGICAL :: PRESSURE_CORRECTION=.FALSE.
NAMELIST /MISC/ PR,SC,TMPA,GVEC,PRESSURE_RELAX_FACTOR,RELAXATION_FACTOR,FYI, &
                CSMAG,RAMP_GX,RAMP_GY,RAMP_GZ,BAROCLINIC, &
                LAPSE_RATE, &
                P_INF,SURF_DEFAULT,EVAC_SURF_DEFAULT, &
                C_FORCED,C_FORCED_CYLINDER,C_FORCED_SPHERE,C_VERTICAL,C_HORIZONTAL,H_FIXED, &
                RESTART,ASSUMED_GAS_TEMPERATURE, &
                LES,DNS,NOISE, CONVECTION,GAMMA,BNDF_DEFAULT, &
                U0,V0,W0,HUMIDITY, ALLOW_SURFACE_DROPLETS,&
                ALLOW_UNDERSIDE_DROPLETS,POROUS_FLOOR, &
                TEXTURE_ORIGIN,NSTRATA, &
                THICKEN_OBSTRUCTIONS,PARTICLE_CFL_MAX, &
                EVAC_PRESSURE_ITERATIONS,EVAC_TIME_ITERATIONS,EVACUATION_MC_MODE, &
                EVACUATION_DRILL, NO_EVACUATION, &
                PRESSURE_CORRECTION,CHECK_POISSON,STRATIFICATION,RESTART_CHID,PARTICLE_CFL, &
                CFL_MAX,CFL_MIN,VN_MAX,VN_MIN,SOLID_PHASE_ONLY,SMOKE_ALBEDO,GROUND_LEVEL, &
                AL2O3,SHARED_FILE_SYSTEM, &
                FLUX_LIMITER,FREEZE_VELOCITY,CFL_VELOCITY_NORM,PERIODIC_TEST, &
                WIND_ONLY,TERRAIN_CASE,COMPUTE_VISCOSITY_TWICE, &
                CONSTANT_PROPERTIES,DYNSMAG,DSMAG_FREQ, &
                CHECK_KINETIC_ENERGY,FISHPAK_BC,FORCE_VECTOR,DEBUG_OPENMP, &
                FDS6,CLIP_MASS_FRACTION,STORE_MU_DNS,CHECK_VN,CHECK_GR, &
                IMMERSED_BOUNDARY_METHOD, &
                H_EDDY,H_LOGLAW,H_CHILTON_COLBURN,LIMITING_DT_RATIO, &
                VAN_DRIEST,HRRPUVCUT_MAX,RUN_AVG_FAC,THERMOPHORETIC_DEPOSITION,TURBULENT_DEPOSITION, &
                VEG_LEVEL_SET,CP_FTMP,HRRPUV_MAX_SMV,TERRAIN_IMAGE,NEW_EVAP, &
                SCALAR_ENERGY_TOLERANCE,TKE_TOLERANCE,MEAN_FORCING,RFAC_FORCING,USE_MAX_FILTER_WIDTH, &
                WFDS,PATCH_VELOCITY,OVERWRITE,UVW_FILE
 
! Physical constants
 
TMPA         = 20._EB                                              ! Ambient temperature (C)
GAMMA        = 1.4_EB                                              ! Heat capacity ratio for air
P_INF        = 101325._EB                                          ! Ambient pressure (Pa)
SIGMA        = 5.67E-8_EB                                          ! Stefan-Boltzmann constant (W/m**2/K**4)
MU_AIR_0     = 1.8E-5_EB                                           ! Dynamic Viscosity of Air at 20 C (kg/m/s)
CP_AIR_0     = 1012._EB                                            ! Specific Heat of Air at 20 C (J/kg/K)
PR_AIR       = 0.7_EB
K_AIR_0      = MU_AIR_0*CP_AIR_0/PR_AIR                            ! Thermal Conductivity of Air at 20 C (W/m/K)
MW_AIR       = 1._EB/(Y_O2_INFTY/32._EB+(1._EB-Y_O2_INFTY)/28._EB) ! g/mol
RHO_SOOT     = 1850._EB                                            ! Density of soot particle (kg/m3)
SMOKE_ALBEDO = 0.3                                                 ! Albedo of smoke

! Variables that are initially set to -1. These are to be reset later.

HUMIDITY     = -1._EB             ! Relative Humidity
 
! Empirical heat transfer constants
 
C_VERTICAL        = 1.31_EB  ! Vertical free convection (Holman, Table 7-2)
C_HORIZONTAL      = 1.52_EB  ! Horizontal free convection 
C_FORCED          = 0.037_EB ! Forced convection coefficient for plates
C_FORCED_CYLINDER = 0.664_EB ! Forced convection coefficient for cylinders
C_FORCED_SPHERE   = 0.6_EB   ! Forced convection coefficient for spheres
PR_ONTH           = PR_AIR**ONTH

H_FIXED                 = -1.      ! obsolete, moved to SURF
ASSUMED_GAS_TEMPERATURE = -1000.   ! Assumed gas temperature, used for diagnostics
 
! Background parameters
 
U0 = 0._EB
V0 = 0._EB
W0 = 0._EB

! Miscellaneous constants

RESTART_CHID   = CHID
RESTART        = .FALSE.
RADIATION      = .TRUE.
CO_PRODUCTION  = .FALSE.
CHECK_POISSON  = .FALSE.
NOISE          = .TRUE.
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
NO_EVACUATION            = .FALSE.
EVACUATION_DRILL         = .FALSE.

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
RUN_AVG_FAC          = 0.5
THICKEN_OBSTRUCTIONS = .FALSE.
CFL_MAX              = 1.0_EB       ! Stability bounds
CFL_MIN              = 0.8_EB
VN_MAX               = 1.0_EB
VN_MIN               = 0.8_EB
VEG_LEVEL_SET        = .FALSE.
TERRAIN_IMAGE        = 'xxxnull'

 
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
   IF (DNS) THEN
      FLUX_LIMITER=4
      VN_MIN=0.4_EB
      VN_MAX=0.5_EB
   ELSE
      FLUX_LIMITER=2
      DYNSMAG=.TRUE.
   ENDIF
   BAROCLINIC=.TRUE.
   CHECK_VN=.TRUE.
   CHECK_GR=.TRUE.
   !!CLIP_MASS_FRACTION=.TRUE.
   !!RESTRICT_TIME_STEP=.FALSE.
   CP_FTMP=.TRUE.
   EXTINCTION2=.TRUE.
   SUPPRESSION_SEARCH=.FALSE. ! tied to EXTINCTION2=T

   ! reread the line to pick up any user-specified options

   REWIND(LU_INPUT)
   CALL CHECKREAD('MISC',LU_INPUT,IOS)
   IF (IOS==0) READ(LU_INPUT,MISC)
   REWIND(LU_INPUT)
ENDIF

! Temperature conversions

H0    = 0.5_EB*(U0**2+V0**2+W0**2)
TMPA  = TMPA + TMPM
TMPA4 = TMPA**4
 
! Humidity (40% by default, but limited for high temps)
 
IF (HUMIDITY < 0._EB) THEN
   HUMIDITY = 40._EB*MIN(1._EB,AMBIENT_WATER_VAPOR(40._EB,313.15_EB)/AMBIENT_WATER_VAPOR(40._EB,TMPA))
ENDIF

! Miscellaneous
 
HCH    = C_HORIZONTAL
HCV    = C_VERTICAL
ASSUMED_GAS_TEMPERATURE = ASSUMED_GAS_TEMPERATURE + TMPM
TEX_ORI = TEXTURE_ORIGIN
FVEC = FORCE_VECTOR
IF (IMMERSED_BOUNDARY_METHOD==2) STORE_MU_DNS=.TRUE.
GRAV = SQRT(DOT_PRODUCT(GVEC,GVEC))
 
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
 
! Min and Max values of species
 
YYMIN  = 0._EB
YYMAX  = 1._EB

IF (H_FIXED>=0._EB) THEN
   WRITE(MESSAGE,'(A)')  'ERROR: H_FIXED moved to SURF line'
   CALL SHUTDOWN(MESSAGE)
ENDIF 

IF (PRESSURE_CORRECTION) THEN
   WRITE(MESSAGE,'(A)')  'ERROR: PRESSURE_CORRECTION is a deprecated feature'
   CALL SHUTDOWN(MESSAGE)
ENDIF 

IF (FLUX_LIMITER<0 .OR. FLUX_LIMITER>4) THEN
   WRITE(MESSAGE,'(A)')  'ERROR: Permissible values for FLUX_LIMITER=0,1,2,3,4'
   CALL SHUTDOWN(MESSAGE)
ENDIF 

IF (AL2O3) RHO_SOOT = 4000._EB

! Level set based model of firespread in vegetation
IF(VEG_LEVEL_SET) WIND_ONLY = .TRUE.

IF (HRRPUV_MAX_SMV<0.0) THEN
   HRRPUV_MAX_SMV=1200.0
   USE_HRRPUV_MAX_SMV=0
ELSE
   USE_HRRPUV_MAX_SMV=1
ENDIF

END SUBROUTINE READ_MISC



SUBROUTINE READ_DUMP

! Read parameters associated with output files
 
INTEGER :: N,I_DUM
NAMELIST /DUMP/ RENDER_FILE,SMOKE3D,SMOKE3D_QUANTITY,FLUSH_FILE_BUFFERS,MASS_FILE, &
                DT_CTRL,DT_PART,DT_MASS,DT_HRR,DT_DEVC,DT_DEVC_LINE,DT_PROF,DT_SLCF,DT_PL3D,DT_ISOF,DT_BNDF,DT_FLUSH,DT_VEG, &
                NFRAMES,DT_RESTART,DEBUG,TIMING,CHECK_VOLUME_FLOW,COLUMN_DUMP_LIMIT,MAXIMUM_DROPLETS,WRITE_XYZ, &
                PLOT3D_QUANTITY,STATUS_FILES,DEVC_COLUMN_LIMIT,CTRL_COLUMN_LIMIT,PLOT3D_SPEC_ID,SMOKE3D_SPEC_ID, &
                PLOT3D_PART_ID,CENTERED_HRR_TIME,CENTERED_DEVC_TIME,PLOT3D_VELO_INDEX,VELOCITY_ERROR_FILE, &
                DT_SL3D
 
RENDER_FILE          = 'null'
NFRAMES              = 1000 
SMOKE3D              = .TRUE.
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
PLOT3D_VELO_INDEX  = 0
MAXIMUM_DROPLETS  = 500000
COLUMN_DUMP_LIMIT = .TRUE.  ! Limit csv files to 255 columns
 
DT_BNDF      = -1._EB
DT_RESTART   = 1000000._EB
DT_FLUSH     = -1._EB
DT_DEVC      = -1._EB
DT_DEVC_LINE = 0.25_EB*(T_END-T_BEGIN)
DT_HRR       = -1._EB
DT_ISOF      = -1._EB
DT_MASS      = -1._EB
DT_PART      = -1._EB
IF (FDS6) THEN
   DT_PL3D = 1.E10_EB
ELSE
   DT_PL3D = -1._EB
ENDIF
DT_PROF      = -1._EB
DT_SLCF      = -1._EB
DT_SL3D      = -1._EB
DT_CTRL      = -1._EB
DT_VEG       = -1._EB
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
IF (DT_SL3D < 0._EB) DT_SL3D = (T_END - T_BEGIN)/5._EB 
IF (DT_CTRL < 0._EB) DT_CTRL = (T_END - T_BEGIN)/REAL(NFRAMES,EB)
IF (DT_FLUSH< 0._EB) DT_FLUSH= (T_END - T_BEGIN)/REAL(NFRAMES,EB)
IF (DT_VEG  < 0._EB) DT_VEG  = (T_END - T_BEGIN)/REAL(NFRAMES,EB)

! Check Plot3D QUANTITIES

PLOOP: DO N=1,5
   CALL GET_QUANTITY_INDEX(PLOT3D_SMOKEVIEW_LABEL(N),PLOT3D_SMOKEVIEW_BAR_LABEL(N),PLOT3D_QUANTITY_INDEX(N),I_DUM, &
                           PLOT3D_SPEC_INDEX(N),PLOT3D_PART_INDEX(N),I_DUM,I_DUM,'PLOT3D', &
                           PLOT3D_QUANTITY(N),'null',PLOT3D_SPEC_ID(N),PLOT3D_PART_ID(N),'null','null')
               
ENDDO PLOOP

! Check SMOKE3D viability

IF (TWO_D .OR. SOLID_PHASE_ONLY)                           SMOKE3D = .FALSE.

IF (SMOKE3D_QUANTITY=='null') THEN
   IF (SOOT_INDEX > 0)  THEN
      SMOKE3D_QUANTITY = 'MASS FRACTION'
      SMOKE3D_SPEC_ID  = 'SOOT'
   ELSE
      IF(N_REACTIONS > 0)  THEN
         SMOKE3D_QUANTITY = 'HRRPUV'
      ELSE
         SMOKE3D = .FALSE.
      ENDIF      
   ENDIF
ENDIF

IF (SMOKE3D) THEN
   CALL GET_QUANTITY_INDEX(SMOKE3D_SMOKEVIEW_LABEL,SMOKE3D_SMOKEVIEW_BAR_LABEL,SMOKE3D_QUANTITY_INDEX,I_DUM, &
                           SMOKE3D_SPEC_INDEX,I_DUM,I_DUM,I_DUM,'SMOKE3D', &
                           SMOKE3D_QUANTITY,'null',SMOKE3D_SPEC_ID,'null','null','null')
ENDIF

END SUBROUTINE READ_DUMP

 
SUBROUTINE READ_SPEC

USE PHYSICAL_FUNCTIONS, ONLY : AMBIENT_WATER_VAPOR
USE PROPERTY_DATA, ONLY: GAS_PROPS
USE DEVICE_VARIABLES, ONLY : PROPERTY_TYPE,PROPERTY,N_PROP
REAL(EB) :: MASS_FRACTION_0,MW,SIGMALJ,EPSILONKLJ,VISCOSITY,CONDUCTIVITY,DIFFUSIVITY,MASS_EXTINCTION_COEFFICIENT, &
            SPECIFIC_HEAT,SPECIFIC_ENTHALPY,REFERENCE_TEMPERATURE
INTEGER  :: N_SPEC_READ,N,MODE,N2,I_SKIP=0
LOGICAL  :: ABSORBING,LUMPED_SPECIES_ONLY,BACKGROUND,BACKGROUND_ONLY
CHARACTER(30) :: FORMULA
NAMELIST /SPEC/ MASS_FRACTION_0,MW,FYI,ID,SIGMALJ,EPSILONKLJ,CONDUCTIVITY,VISCOSITY,DIFFUSIVITY,ABSORBING,BACKGROUND, &
                MASS_EXTINCTION_COEFFICIENT,FORMULA,SPECIFIC_HEAT,SPECIFIC_ENTHALPY,REFERENCE_TEMPERATURE,LUMPED_SPECIES_ONLY, &
                BACKGROUND_ONLY

!Init variables for no &SPEC case

ABSORBING                   = .FALSE.
CONDUCTIVITY                = -1._EB
DIFFUSIVITY                 = -1._EB
EPSILONKLJ                  =  0._EB
LUMPED_SPECIES_ONLY         = .FALSE.
MASS_EXTINCTION_COEFFICIENT = 8700._EB  ! m2/kg
MASS_FRACTION_0             = -1._EB
MODE                        = GAS_SPECIES
MW                          =  0._EB
REFERENCE_TEMPERATURE       = 25._EB
SIGMALJ                     =  0._EB
SPECIFIC_HEAT               = -1._EB
SPECIFIC_ENTHALPY           = -1._EB
VISCOSITY                   = -1._EB

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
   LUMPED_SPECIES_ONLY=.FALSE.
   BACKGROUND = .FALSE.
   BACKGROUND_ONLY = .FALSE.
   CALL CHECKREAD('SPEC',LU_INPUT,IOS) 
   IF (IOS==1) EXIT COUNT_SPEC_LOOP
   ID = 'null'
   READ(LU_INPUT,NML=SPEC,END=29,ERR=30,IOSTAT=IOS)
   IF (ID=='null') THEN
      WRITE(MESSAGE,'(A,I2,A)') 'ERROR: Species',N_SPECIES+1, 'needs a name (ID=...)'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   N_SPECIES   = N_SPECIES   + 1
   N_SPEC_READ = N_SPEC_READ + 1
   IF (BACKGROUND_ONLY) BACKGROUND = .TRUE.   
   IF (LUMPED_SPECIES_ONLY) THEN
      IF (SIMPLE_CHEMISTRY) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Cannot define additional LUMPED_SPECIES_ONLY species when using simple chemistry REAC input'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      N_LUMPED_SPECIES = N_LUMPED_SPECIES + 1
   ELSE
      IF (BACKGROUND) THEN 
         DEFINED_BACKGROUND = .TRUE.
         IF (SIMPLE_CHEMISTRY) THEN
            WRITE(MESSAGE,'(A)') 'ERROR: Cannot define both a BACKGROUND SPEC with simple chemistry REAC input'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         IF (I_BACKGROUND /= 0) THEN
            WRITE(MESSAGE,'(A)') 'ERROR: Can only define one BACKGROUND species'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         I_BACKGROUND = N_SPECIES
         IF (.NOT.BACKGROUND_ONLY) N_GAS_SPECIES = N_GAS_SPECIES + 1
      ELSE
         N_GAS_SPECIES = N_GAS_SPECIES + 1
      ENDIF
   ENDIF   
   IF (ID=='WATER VAPOR') THEN
      IF (H2O_INDEX > 0) THEN
         N_SPECIES = N_SPECIES - 1
      ELSE
         SPECIES_ID(N_SPECIES)=ID
      ENDIF
   ELSEIF (ID=='CARBON DIOXIDE')THEN
      IF (CO2_INDEX > 0) THEN
         N_SPECIES = N_SPECIES - 1
      ELSE
         SPECIES_ID(N_SPECIES)=ID
      ENDIF
   ELSEIF (ID=='CARBON MONOXIDE')THEN
      IF (CO_INDEX > 0) THEN
         N_SPECIES = N_SPECIES - 1
      ELSE
         SPECIES_ID(N_SPECIES)=ID
      ENDIF
   ELSEIF (ID=='OXYGEN')THEN
      IF (O2_INDEX > 0) THEN
         N_SPECIES = N_SPECIES - 1
      ELSE
         SPECIES_ID(N_SPECIES)=ID
      ENDIF
   ELSEIF (ID=='SOOT')THEN
      IF (SOOT_INDEX > 0) THEN
         N_SPECIES = N_SPECIES - 1
      ELSE
         SPECIES_ID(N_SPECIES)=ID
      ENDIF
   ELSEIF (ID=='NITROGEN')THEN
      IF (N2_INDEX > 0) THEN
         N_SPECIES = N_SPECIES - 1
      ELSE
         SPECIES_ID(N_SPECIES)=ID
      ENDIF
   ELSEIF (ID=='HYDROGEN')THEN
      IF (H2_INDEX > 0) THEN
         N_SPECIES = N_SPECIES - 1
      ELSE
         SPECIES_ID(N_SPECIES)=ID
      ENDIF
   ELSE
      SPECIES_ID(N_SPECIES)=ID
   ENDIF
   30 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I2)') 'ERROR: Problem with SPECies number',N_SPECIES+1
         CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_SPEC_LOOP
29 REWIND(LU_INPUT)

!Set background if undefined
IF (I_BACKGROUND == 0) THEN
   N_SPECIES = N_SPECIES + 1
   I_BACKGROUND = N_SPECIES
   I_SKIP = I_SKIP+1
   SPECIES_ID(N_SPECIES) = 'AIR'
   BACKGROUND_ONLY = .TRUE.
ENDIF

! Allocate species-related arrays

ALLOCATE(SPECIES(0:N_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ','SPECIES',IZERO)

SPECIES(1:N_SPECIES)%ID = SPECIES_ID(1:N_SPECIES)
IF (OTHER_INDEX > 0) SPECIES(OTHER_INDEX)%MW = MW_OTHER
IF (SIMPLE_CHEMISTRY) SPECIES(1:N_LUMPED_SPECIES)%LUMPED_SPECIES_ONLY=.TRUE.
! Loop through all species and assign properties
N2 = 0 
REWIND(LU_INPUT)
SPEC_LOOP: DO N=1,N_SPEC_READ
   ! Set default values for all species

   WRITE(FORMULA,'(I2.2)') N

   ABSORBING                   = .FALSE.
   CONDUCTIVITY                = -1._EB
   DIFFUSIVITY                 = -1._EB
   EPSILONKLJ                  =  0._EB
   LUMPED_SPECIES_ONLY         = .FALSE.
   MASS_EXTINCTION_COEFFICIENT = 8700._EB  ! m2/kg
   MASS_FRACTION_0             = -1._EB
   MODE                        = GAS_SPECIES
   MW                          =  0._EB
   REFERENCE_TEMPERATURE       = 25._EB
   SIGMALJ                     =  0._EB
   SPECIFIC_HEAT               = -1._EB
   SPECIFIC_ENTHALPY           = -1._EB
   VISCOSITY                   = -1._EB

   ! Either read in the SPEC line

   CALL CHECKREAD('SPEC',LU_INPUT,IOS)
   READ(LU_INPUT,NML=SPEC)
   IF (ID=='WATER VAPOR') THEN
      IF (H2O_INDEX > 0) THEN
         SPECIES(H2O_INDEX)%LUMPED_SPECIES_ONLY = LUMPED_SPECIES_ONLY
         IF (.NOT. LUMPED_SPECIES_ONLY .AND. SIMPLE_CHEMISTRY .AND. MASS_FRACTION_0 > 0) THEN
            WRITE(LU_ERR,'(A)') 'WARNING: Specified WATER VAPOR mass fraction will be in addition to the HUMDITY amount in'
            WRITE(LU_ERR,'(A)') 'the background species mixture'
            SPECIES(H2O_INDEX)%YY0 = MASS_FRACTION_0
         ENDIF         
         N2 = N2 + 1
         CYCLE
      ELSE
         H2O_INDEX = N
      ENDIF         
   ENDIF
   IF (ID=='CARBON DIOXIDE') THEN
      IF (CO2_INDEX > 0) THEN
      SPECIES(CO2_INDEX)%LUMPED_SPECIES_ONLY = LUMPED_SPECIES_ONLY
         IF (.NOT. LUMPED_SPECIES_ONLY .AND. SIMPLE_CHEMISTRY .AND. MASS_FRACTION_0 > 0) THEN
            WRITE(LU_ERR,'(A)') 'WARNING: Specified CARBON DIOXIDE mass fraction will be in addition to the HUMDITY amount in'
            WRITE(LU_ERR,'(A)') 'the background species mixture'
            SPECIES(CO2_INDEX)%YY0 = MASS_FRACTION_0
         ENDIF
         N2 = N2 + 1
         CYCLE
      ELSE
         CO2_INDEX = N
      ENDIF         
   ENDIF
   IF (ID=='CARBON MONOXIDE') THEN
      IF (CO_INDEX > 0) THEN
         SPECIES(CO_INDEX)%LUMPED_SPECIES_ONLY = LUMPED_SPECIES_ONLY
         IF (.NOT. LUMPED_SPECIES_ONLY .AND. SIMPLE_CHEMISTRY .AND. MASS_FRACTION_0 > 0) &
            SPECIES(CO_INDEX)%YY0 = MASS_FRACTION_0    
         N2 = N2 + 1
         CYCLE
      ELSE
         CO_INDEX = N
      ENDIF         
   ENDIF
   IF (ID=='NITROGEN') THEN
      IF (N2_INDEX > 0) THEN
         SPECIES(N2_INDEX)%LUMPED_SPECIES_ONLY = LUMPED_SPECIES_ONLY
         IF (.NOT. LUMPED_SPECIES_ONLY .AND. SIMPLE_CHEMISTRY .AND. MASS_FRACTION_0 > 0) THEN
            WRITE(LU_ERR,'(A)') 'WARNING: Specified NITROGEN mass fraction will be in addition to the HUMDITY amount in'
            WRITE(LU_ERR,'(A)') 'the background species mixture'
            SPECIES(N2_INDEX)%YY0 = MASS_FRACTION_0
         ENDIF
         N2 = N2 + 1
         CYCLE
      ELSE
         N2_INDEX = N
      ENDIF         
   ENDIF      
   IF (ID=='OXYGEN') THEN
      IF (O2_INDEX > 0) THEN
      SPECIES(O2_INDEX)%LUMPED_SPECIES_ONLY = LUMPED_SPECIES_ONLY
         IF (.NOT. LUMPED_SPECIES_ONLY .AND. SIMPLE_CHEMISTRY .AND. MASS_FRACTION_0 > 0) THEN
            WRITE(LU_ERR,'(A)') 'WARNING: Specified OXYGEN mass fraction will be in addition to the HUMDITY amount in'
            WRITE(LU_ERR,'(A)') 'the background species mixture and it will not participate in combustion'
            SPECIES(O2_INDEX)%YY0 = MASS_FRACTION_0
         ENDIF
         N2 = N2 + 1
         CYCLE
      ELSE
         O2_INDEX = N
      ENDIF         
   ENDIF      
   IF (ID=='HYDROGEN') THEN
      IF (H2_INDEX > 0) THEN
         SPECIES(H2_INDEX)%LUMPED_SPECIES_ONLY = LUMPED_SPECIES_ONLY
         IF (.NOT. LUMPED_SPECIES_ONLY .AND. SIMPLE_CHEMISTRY .AND. MASS_FRACTION_0 > 0) &
            SPECIES(H2_INDEX)%YY0 = MASS_FRACTION_0    
         N2 = N2 + 1
         CYCLE
      ELSE
         H2_INDEX = N
      ENDIF         
   ENDIF      
   IF (ID=='SOOT') THEN
      IF (SOOT_INDEX > 0) THEN
         SPECIES(SOOT_INDEX)%LUMPED_SPECIES_ONLY = LUMPED_SPECIES_ONLY
         IF (.NOT. LUMPED_SPECIES_ONLY .AND. SIMPLE_CHEMISTRY .AND. MASS_FRACTION_0 > 0) &
            SPECIES(SOOT_INDEX)%YY0 = MASS_FRACTION_0    
         N2 = N2 + 1
         CYCLE
      ELSE
         SOOT_INDEX = N
      ENDIF         
   ENDIF      
 
   ! If the user has not specified an initial mass fraction, do it here

   IF (MASS_FRACTION_0<0._EB) THEN
      IF (ID=='WATER VAPOR' .AND. .NOT. SIMPLE_CHEMISTRY) THEN
         MASS_FRACTION_0 = AMBIENT_WATER_VAPOR(HUMIDITY,TMPA)
      ELSE
         MASS_FRACTION_0 = 0._EB
      ENDIF
   ENDIF

   ! Special case for soot

   IF (ID=='SOOT') MODE = AEROSOL_SPECIES

   ! Fill up the SPECIES array
   SS => SPECIES(N)

   SS%MU_USER = VISCOSITY
   SS%K_USER  = CONDUCTIVITY
   SS%D_USER  = DIFFUSIVITY
   SS%EPSK    = EPSILONKLJ
   SS%SIG     = SIGMALJ

   SS%ABSORBING                   = ABSORBING
   SS%FORMULA                     = FORMULA
   SS%ID                          = ID
   SS%MASS_EXTINCTION_COEFFICIENT = MASS_EXTINCTION_COEFFICIENT
   SS%MW                          = MW
   SS%MODE                        = MODE
   SS%REFERENCE_TEMPERATURE       = REFERENCE_TEMPERATURE + TMPM
   SS%SPECIFIC_HEAT               = SPECIFIC_HEAT*1000._EB
   SS%SPECIFIC_ENTHALPY           = SPECIFIC_ENTHALPY*1000._EB
   SS%YY0                         = MASS_FRACTION_0
   SS%LUMPED_SPECIES_ONLY         = LUMPED_SPECIES_ONLY
   IF ((SPECIFIC_HEAT < 0._EB .AND. SPECIFIC_ENTHALPY > 0._EB).OR.(SPECIFIC_HEAT > 0._EB .AND. SPECIFIC_ENTHALPY < 0._EB)) THEN
      WRITE(MESSAGE,'(A,I2,A)') 'ERROR: SPECies ',N,': Specify both SPECIFIC_HEAT and SPECIFIC_ENTHALPY'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
 
ENDDO SPEC_LOOP

IF (I_BACKGROUND > N_SPEC_READ) THEN
   SS => SPECIES(I_BACKGROUND)
   SS%MU_USER = VISCOSITY
   SS%K_USER  = CONDUCTIVITY
   SS%D_USER  = DIFFUSIVITY
   SS%EPSK    = EPSILONKLJ
   SS%SIG     = SIGMALJ
   SS%ABSORBING                   = ABSORBING
   SS%FORMULA                     = 'AIR'
   SS%ID                          = SPECIES_ID(I_BACKGROUND)
   SS%MASS_EXTINCTION_COEFFICIENT = MASS_EXTINCTION_COEFFICIENT
   SS%MW                          = MW
   SS%MODE                        = MODE
   SS%REFERENCE_TEMPERATURE       = REFERENCE_TEMPERATURE + TMPM
   SS%SPECIFIC_HEAT               = SPECIFIC_HEAT*1000._EB
   SS%SPECIFIC_ENTHALPY           = SPECIFIC_ENTHALPY*1000._EB
   SS%YY0                         = MASS_FRACTION_0
   SS%LUMPED_SPECIES_ONLY         = LUMPED_SPECIES_ONLY
ENDIF

N2 = 0

DO N = 1,N_SPECIES
   SS => SPECIES(N)
   CALL GAS_PROPS(SS%ID,SS%SIG,SS%EPSK,SS%MW,SS%ABSORBING,SS%FORMULA)    
   IF (SPECIES(N)%LUMPED_SPECIES_ONLY) CYCLE
   IF (N==I_BACKGROUND .AND. BACKGROUND_ONLY) CYCLE
   N2 = N2 + 1
   SPECIES(N)%Y_INDEX = N2  
   SELECT CASE (TRIM(SPECIES(N)%ID))
      CASE('WATER VAPOR') 
         I_WATER = SPECIES(N)%Y_INDEX
      CASE('CARBON DIOXIDE')
         I_CO2 = SPECIES(N)%Y_INDEX
      CASE('CARBON MONOXIDE')
         I_CO = SPECIES(N)%Y_INDEX
      CASE('OXYGEN')
         I_O2 = SPECIES(N)%Y_INDEX
      CASE('SOOT')
         I_SOOT = SPECIES(N)%Y_INDEX
      CASE('NITROGEN')
         I_N2 = SPECIES(N)%Y_INDEX
      CASE('HYDROGEN')
         I_H2 = SPECIES(N)%Y_INDEX
   END SELECT 
ENDDO

SPECIES(0) = SPECIES(I_BACKGROUND)

END SUBROUTINE READ_SPEC

SUBROUTINE READ_SMIX
USE PROPERTY_DATA, ONLY: GAS_PROPS
CHARACTER(30) :: SPEC_ID(1:MAX_SPECIES),ID,FORMULA
REAL(EB) :: VOLUME_FRACTION(1:MAX_SPECIES),MASS_FRACTION(1:MAX_SPECIES),CONVERSION
INTEGER :: N,NS,NS2
LOGICAL :: BACKGROUND
TYPE(SPECIES_MIXTURE_TYPE), POINTER :: SM
TYPE(SPECIES_TYPE), POINTER, DIMENSION(:) :: SPEC_ALLOC
NAMELIST /SMIX/ BACKGROUND,SPEC_ID,VOLUME_FRACTION,MASS_FRACTION,ID

IF (N_MIX_SPECIES == 0) THEN
   COUNT_SMIX_LOOP: DO
      BACKGROUND = .FALSE.
      ID = 'null'
      CALL CHECKREAD('SMIX',LU_INPUT,IOS) 
      IF (IOS==1) EXIT COUNT_SMIX_LOOP
      READ(LU_INPUT,NML=SMIX,END=29,ERR=30,IOSTAT=IOS)
      N_MIX_SPECIES   = N_MIX_SPECIES + 1
      30 IF (IOS>0) THEN
            WRITE(MESSAGE,'(A,I2)') 'ERROR: Problem with Species MIXture number ',N_MIX_SPECIES+1
            CALL SHUTDOWN(MESSAGE)
         ENDIF
   ENDDO COUNT_SMIX_LOOP
   29 REWIND(LU_INPUT)

   ALLOCATE(SPECIES_MIXTURE(1:N_MIX_SPECIES))

   READ_SMIX_LOOP: DO N=1,N_MIX_SPECIES   
   
      SM => SPECIES_MIXTURE(N)
            
      MASS_FRACTION = -1._EB
      VOLUME_FRACTION = -1._EB
      BACKGROUND = .FALSE.
      CONVERSION = 0._EB
      SPEC_ID = 'null'
      READ(LU_INPUT,NML=SMIX)
  
      IF (MASS_FRACTION(1) > 0._EB .AND. VOLUME_FRACTION(1) > 0._EB) THEN
         WRITE(MESSAGE,'(A,I2,A)') 'ERROR: Problem with Species MIXture number ',N, &
                                   '. Specify one of MASS_FRACTION or VOLUME_FRACTION'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (BACKGROUND) THEN
         IF (N_MIX_SPECIES > 1) THEN
            WRITE(MESSAGE,'(A)') 'ERROR: Background mixture must be defined first'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         IF (.NOT. DEFINED_BACKGROUND) THEN
            DEFINED_BACKGROUND = .TRUE.
            I_BACKGROUND = N+N_SPECIES
         ELSE
            WRITE(MESSAGE,'(A)') 'ERROR: Can only define one background species'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         IF(TRIM(ID)=='null') ID = 'LUMPED SPECIES 0'
      ELSE   
         IF (TRIM(ID)=='null') THEN
            IF (N_MIX_SPECIES < 10) THEN
               WRITE(ID,'(A,I1)') 'LUMPED SPECIES ',N_MIX_SPECIES
            ELSE
               WRITE(ID,'(A,I2)') 'LUMPED SPECIES ',N_MIX_SPECIES
            ENDIF
         ENDIF
      ENDIF   
      
      SM%N_SUB_SPECIES = 0
      SM%ID = ID
      
      SM => SPECIES_MIXTURE(N)
      COUNT_SPEC: DO NS = 1,MAX_SPECIES
         IF (TRIM(SPEC_ID(NS)) /= 'null') THEN
            SM%N_SUB_SPECIES = SM%N_SUB_SPECIES + 1
         ELSE
            EXIT
         ENDIF
      ENDDO COUNT_SPEC

      ALLOCATE (SM%SPEC_ID(SM%N_SUB_SPECIES))
      SM%SPEC_ID(1:SM%N_SUB_SPECIES) = SPEC_ID(1:SM%N_SUB_SPECIES)
      ALLOCATE (SM%SPEC_INDEX(SM%N_SUB_SPECIES))
      SM%SPEC_INDEX = -1
      ALLOCATE (SM%VOLUME_FRACTION(SM%N_SUB_SPECIES))
      SM%VOLUME_FRACTION(1:SM%N_SUB_SPECIES) = VOLUME_FRACTION(1:SM%N_SUB_SPECIES)
      ALLOCATE (SM%MASS_FRACTION(SM%N_SUB_SPECIES))
      SM%MASS_FRACTION(1:SM%N_SUB_SPECIES) = MASS_FRACTION(1:SM%N_SUB_SPECIES)
      DO NS = 1,SM%N_SUB_SPECIES
         IF (MASS_FRACTION(NS) < 0._EB .AND. VOLUME_FRACTION(NS) < 0._EB) THEN
            WRITE(MESSAGE,'(A,I2,A,I2)') 'ERROR: Problem with Species MIXture number ',N,'--Unspecified FRACTION for subspecies ',NS
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         FIND_SPEC_ID: DO NS2 = 1,N_SPECIES-N_MIX_SPECIES
            IF (TRIM(SPECIES(NS2)%ID) == TRIM(SPEC_ID(NS))) THEN
               SM%SPEC_INDEX(NS) = NS2
               EXIT FIND_SPEC_ID
            ENDIF
            IF(NS2==N_SPECIES-N_MIX_SPECIES) THEN
               WRITE(MESSAGE,'(A,I2,A,A)') 'ERROR: Problem with Species MIXture number ',N,'--Could not find SPEC:',SPEC_ID(NS)
               CALL SHUTDOWN(MESSAGE)
            ENDIF
         ENDDO FIND_SPEC_ID
         IF (MASS_FRACTION(NS) >= 0._EB) THEN
            CONVERSION = CONVERSION + MASS_FRACTION(NS) / SPECIES(SM%SPEC_INDEX(NS))%MW
         ELSE
            CONVERSION = CONVERSION + VOLUME_FRACTION(NS) * SPECIES(SM%SPEC_INDEX(NS))%MW
         ENDIF
      ENDDO
      IF (MASS_FRACTION(1) > 0._EB) THEN
         DO NS = 1,SM%N_SUB_SPECIES
            SM%VOLUME_FRACTION(NS) = MASS_FRACTION(NS) / SPECIES(SM%SPEC_INDEX(NS))%MW / CONVERSION
         ENDDO
      ELSE
         DO NS = 1,SM%N_SUB_SPECIES
            SM%MASS_FRACTION(NS) = VOLUME_FRACTION(NS) * SPECIES(SM%SPEC_INDEX(NS))%MW / CONVERSION
         ENDDO
      ENDIF
      SM%MASS_FRACTION = SM%MASS_FRACTION / SUM(SM%MASS_FRACTION)
      SM%VOLUME_FRACTION = SM%VOLUME_FRACTION / SUM(SM%VOLUME_FRACTION)
      SM%MW = 0._EB
      DO NS = 1,SM%N_SUB_SPECIES
         SM%MW = SM%MW + SM%VOLUME_FRACTION(NS) * SPECIES(SM%SPEC_INDEX(NS))%MW
      ENDDO     
      SM%RSUM = R0/SM%MW
   ENDDO READ_SMIX_LOOP

   REWIND (LU_INPUT)   

ENDIF

PROC_SMIX_IF: IF (N_MIX_SPECIES > 0) THEN
   IF (SIMPLE_CHEMISTRY) I_BACKGROUND = N_SPECIES+1
   IF (I_BACKGROUND > N_SPECIES) THEN 
      N_SPECIES = N_SPECIES - 1
      I_BACKGROUND = I_BACKGROUND - 1
   ENDIF

   ALLOCATE(SPEC_ALLOC(0:N_SPECIES))
   SPEC_ALLOC(0:N_SPECIES) = SPECIES(0:N_SPECIES)
   DEALLOCATE(SPECIES)

   N_SPECIES = N_SPECIES + N_MIX_SPECIES

   ALLOCATE(SPECIES(0:N_SPECIES))
   SPECIES(0:N_SPECIES - N_MIX_SPECIES) = SPEC_ALLOC(0:N_SPECIES - N_MIX_SPECIES)
   DEALLOCATE(SPEC_ALLOC)

   I_Z_MIN = 999
   I_Z_MAX =-999

   DO NS=1,N_MIX_SPECIES
      SM=>SPECIES_MIXTURE(NS)
      DO NS2=1,SM%N_SUB_SPECIES
         SPECIES(SM%SPEC_INDEX(NS2))%GAS_SPECIES_ONLY=.FALSE.
      END DO
      SS => SPECIES(N_SPECIES - N_MIX_SPECIES + NS) 
      IF (I_BACKGROUND > N_SPECIES - N_MIX_SPECIES + 1) THEN
         WRITE(FORMULA,'(A,I1)') 'Z',NS-1
      ELSE
         WRITE(FORMULA,'(A,I1)') 'Z',NS   
      ENDIF
      IF (N_SPECIES - N_MIX_SPECIES + NS == I_BACKGROUND) CALL GAS_PROPS(SS%ID,SS%SIG,SS%EPSK,SS%MW,SS%ABSORBING,SS%FORMULA)    
      SS%FORMULA = FORMULA
      SS%ABSORBING = .TRUE.
      SS%MODE = LUMPED_SPECIES
      SS%MW = SPECIES_MIXTURE(NS)%MW  
      SS%ID = SPECIES_MIXTURE(NS)%ID
      IF (N_SPECIES - N_MIX_SPECIES + NS == I_BACKGROUND) SPECIES(0) = SPECIES(N_SPECIES - N_MIX_SPECIES + NS)
      IF (I_BACKGROUND > N_SPECIES - N_MIX_SPECIES .AND. NS==1) CYCLE
      N_GAS_SPECIES = N_GAS_SPECIES + 1
      I_Z_MIN = MIN(N_GAS_SPECIES,I_Z_MIN)
      I_Z_MAX = MAX(N_GAS_SPECIES,I_Z_MAX)
      SS%Y_INDEX = N_GAS_SPECIES
   ENDDO

   IF (SIMPLE_CHEMISTRY) THEN
      I_FUEL = I_Z_MIN
      IF (CO_PRODUCTION) THEN
         I_PROG_CO = I_FUEL + 1
         I_PROG_F = I_PROG_CO+1
      ELSE
         I_PROG_F = I_FUEL + 1
         IF (SOOT_DEPOSITION) I_PROG_SOOT = I_PROG_F + 1
      ENDIF
   ENDIF
ENDIF PROC_SMIX_IF

END SUBROUTINE READ_SMIX
 
SUBROUTINE READ_REAC_1
USE PHYSICAL_FUNCTIONS, ONLY : AMBIENT_WATER_VAPOR
CHARACTER(30) :: FUEL,SPEC_ID(1:MAX_SPECIES)
LOGICAL :: IDEAL
INTEGER :: NN,NNN
REAL(EB) :: Y_F_INLET, &
            H2_YIELD,SOOT_YIELD,CO_YIELD, Y_F_LFL,X_O2_LL,EPUMO2,BOF,&
            CRITICAL_FLAME_TEMPERATURE,HEAT_OF_COMBUSTION,MASS(MAX_SPECIES),NU(0:MAX_SPECIES),E,N_S(0:MAX_SPECIES),C,H,N,O,OTHER, &
            FUEL_HEAT_OF_FORMATION,MAXIMUM_VISIBILITY,AUTO_IGNITION_TEMPERATURE,SUM1,THRESHOLD_TEMPERATURE
TYPE (SPECIES_MIXTURE_TYPE), POINTER :: SM
NAMELIST /REAC/ E,BOF,HEAT_OF_COMBUSTION,FYI,FUEL,EPUMO2,ID, N_S,EXTINCTION2,&
                Y_O2_INFTY,Y_F_INLET, CO_PRODUCTION,SOOT_DEPOSITION,SUPPRESSION,SUPPRESSION_SEARCH,&
                H2_YIELD,SOOT_YIELD,CO_YIELD,Y_F_LFL,X_O2_LL,CRITICAL_FLAME_TEMPERATURE,NU,SOOT_H_FRACTION, &
                C,H,N,O,OTHER,MW_OTHER,IDEAL,VISIBILITY_FACTOR,MAXIMUM_VISIBILITY, &
                BETA_EDC,Y_P_MIN_EDC,HRRPUA_SHEET,HRRPUV_AVERAGE,FIXED_MIX_TIME, &
                TAU_CHEM,TAU_FLAME,AUTO_IGNITION_TEMPERATURE, &
                MASS,SPEC_ID,THRESHOLD_TEMPERATURE

 
N_REACTIONS = 0
REWIND(LU_INPUT)
 
COUNT_REAC_LOOP: DO
   CALL SET_REAC_DEFAULTS_1
   CALL CHECKREAD('REAC',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_REAC_LOOP
   READ(LU_INPUT,REAC,END=435,ERR=434,IOSTAT=IOS)
   N_REACTIONS = N_REACTIONS + 1
   IF (BOF < 0._EB .AND. E < 0._EB .AND. TRIM(SPEC_ID(1))=='null') SIMPLE_CHEMISTRY = .TRUE.
   IF (SIMPLE_CHEMISTRY .AND. N_REACTIONS > 1) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: can not have more than one reaction when using simple chemistry'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   434 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with REAC ',N_REACTIONS+1
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO COUNT_REAC_LOOP

435 REWIND(LU_INPUT)

IF (N_REACTIONS == 0) RETURN

CALL SET_REAC_DEFAULTS_1
READ(LU_INPUT,REAC)
REWIND(LU_INPUT)

IF ((BOF > 0._EB .OR. E > 0._EB) .AND. (C>0._EB .OR. H>0._EB)) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: cannot use both finite rate REAC and simple chemistry'
   CALL SHUTDOWN(MESSAGE)
ENDIF

IF (SIMPLE_CHEMISTRY .AND. C<0._EB .AND. H<0._EB) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: Must specify fuel chemistry using C and H when using simple chemistry'
   CALL SHUTDOWN(MESSAGE)
ENDIF

IF (.NOT. SUPPRESSION .AND. CO_PRODUCTION) THEN
   WRITE(MESSAGE,'(A)')  'ERROR: Cannot set SUPPRESSION=.FALSE. when CO_PRODUCTION=.TRUE.'
   CALL SHUTDOWN(MESSAGE)
ENDIF 

IF (SOOT_DEPOSITION .AND. CO_PRODUCTION) THEN
   WRITE(MESSAGE,'(A)')  'ERROR: Cannot set CO_PRODUCTION=.TRUE. when SOOT_DEPOSITION=.TRUE.'
   CALL SHUTDOWN(MESSAGE)
ENDIF

IF (.NOT. SIMPLE_CHEMISTRY) RETURN

IF     (     CO_PRODUCTION) THEN
   N_REACTIONS = 2
   N_MIX_SPECIES = 4
ELSEIF (.NOT.CO_PRODUCTION) THEN
   N_REACTIONS = 1
   IF (SOOT_DEPOSITION) THEN
      N_MIX_SPECIES = 4
   ELSE
      N_MIX_SPECIES = 3
   ENDIF
ENDIF

ALLOCATE(SPECIES_MIXTURE(1:N_MIX_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ','SPECIES_MIXTURE',IZERO)

ALLOCATE(REACTION(N_REACTIONS),STAT=IZERO)
CALL ChkMemErr('READ','REACTION',IZERO)

! Read the input file looking for REAC lines
RN => REACTION(1)
RN%AUTO_IGNITION_TEMPERATURE = AUTO_IGNITION_TEMPERATURE + TMPM
RN%EPUMO2             = EPUMO2*1000._EB
RN%FYI                = FYI
RN%HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION*1000._EB
IF (FUEL=='null') FUEL = 'ETHYLENE'
RN%NAME               = ID
RN%FUEL               = FUEL    
RN%Y_O2_LL            = X_O2_LL
RN%Y_F_LFL            = Y_F_LFL

!Fuel and Soot Molecular weights
MW_SOOT = MW_C * (1._EB - SOOT_H_FRACTION) + MW_H * SOOT_H_FRACTION   
MW_FUEL = C * MW_C + H * MW_H + O * MW_O + N * MW_N + OTHER * MW_OTHER

!Set species indices for lumped species lookup
IF (H <= ZERO_P) THEN
   SOOT_H_FRACTION = 0._EB
   H2_YIELD = 0._EB
ENDIF
FUEL_INDEX = 1
SPECIES_ID(1) = FUEL
O2_INDEX = 2
SPECIES_ID(2) = 'OXYGEN'
N2_INDEX = 3
SPECIES_ID(3) = 'NITROGEN'
CO2_INDEX = 4
SPECIES_ID(4) = 'CARBON DIOXIDE'
H2O_INDEX = 5
SPECIES_ID(5) = 'WATER VAPOR'

N_LUMPED_SPECIES = 5
IF (CO_PRODUCTION .OR. CO_YIELD > 0._EB)THEN
   N_LUMPED_SPECIES = N_LUMPED_SPECIES + 1
   CO_INDEX = N_LUMPED_SPECIES
   SPECIES_ID(N_LUMPED_SPECIES) = 'CARBON MONOXIDE'
ENDIF     
IF (SOOT_DEPOSITION .OR. SOOT_YIELD > 0._EB)THEN
   N_LUMPED_SPECIES = N_LUMPED_SPECIES + 1
   SOOT_INDEX = N_LUMPED_SPECIES
   SPECIES_ID(N_LUMPED_SPECIES) = 'SOOT'
ENDIF     
IF (H2_YIELD > 0._EB)THEN
   N_LUMPED_SPECIES = N_LUMPED_SPECIES + 1
   H2_INDEX = N_LUMPED_SPECIES
   SPECIES_ID(N_LUMPED_SPECIES) = 'HYDROGEN'
ENDIF     
IF (OTHER > 0._EB) THEN
   N_LUMPED_SPECIES = N_LUMPED_SPECIES + 1
   OTHER_INDEX = N_LUMPED_SPECIES
   SPECIES_ID(N_LUMPED_SPECIES) = 'OTHER'
ENDIF    
N_SPECIES = N_LUMPED_SPECIES

DO NN=1,N_MIX_SPECIES
   SM => SPECIES_MIXTURE(NN)
   SM%N_SUB_SPECIES=N_LUMPED_SPECIES
   ALLOCATE (SM%SPEC_ID(N_LUMPED_SPECIES))
   ALLOCATE (SM%SPEC_INDEX(N_LUMPED_SPECIES))
   DO NNN = 1,N_LUMPED_SPECIES
      SM%SPEC_INDEX(NNN) = NNN
   ENDDO
   SM%SPEC_ID(1:N_LUMPED_SPECIES) = SPECIES_ID(1:N_LUMPED_SPECIES)
   ALLOCATE (SM%MASS_FRACTION(N_LUMPED_SPECIES))
   SM%MASS_FRACTION = 0._EB
   ALLOCATE (SM%VOLUME_FRACTION(N_LUMPED_SPECIES))
   SM%VOLUME_FRACTION = 0._EB
   WRITE(SM%ID,'(A,I1)') 'LUMPED SPECIES ',NN-1
ENDDO

!Background Mixture uses NASA Earth Fact Sheet (updated 2007.01) for O2 and CO2 with N2 making up the balance + User Humidity
DEFINED_BACKGROUND = .TRUE.

SM => SPECIES_MIXTURE(1)

SM%MASS_FRACTION(O2_INDEX)  = Y_O2_INFTY
SM%MASS_FRACTION(CO2_INDEX)  = 0.0058_EB 
SM%MASS_FRACTION(N2_INDEX)  = 1._EB - SM%MASS_FRACTION(O2_INDEX) - SM%MASS_FRACTION(CO2_INDEX)
SM%MASS_FRACTION(H2O_INDEX)  = AMBIENT_WATER_VAPOR(HUMIDITY,TMPA)
SM%MASS_FRACTION(1:4) = SM%MASS_FRACTION(1:4) * (1._EB-SM%MASS_FRACTION(H2O_INDEX))


SM%VOLUME_FRACTION(N2_INDEX) = SM%MASS_FRACTION(N2_INDEX)/MW_N2
SM%VOLUME_FRACTION(O2_INDEX) = SM%MASS_FRACTION(O2_INDEX)/MW_O2
SM%VOLUME_FRACTION(CO2_INDEX) = SM%MASS_FRACTION(CO2_INDEX)/MW_CO2
SM%VOLUME_FRACTION(H2O_INDEX) = SM%MASS_FRACTION(H2O_INDEX)/MW_H2O

SM%MW = 1._EB/SUM(SM%VOLUME_FRACTION)
SM%RSUM = R0/SM%MW
SM%VOLUME_FRACTION = SM%VOLUME_FRACTION / SUM(SM%VOLUME_FRACTION)

SM => SPECIES_MIXTURE(2)

SM%MASS_FRACTION(FUEL_INDEX) = Y_F_INLET
SM%VOLUME_FRACTION(FUEL_INDEX) = SM%MASS_FRACTION(FUEL_INDEX) / MW_FUEL
SM%MASS_FRACTION(N2_INDEX) = 1._EB - Y_F_INLET
SM%VOLUME_FRACTION(N2_INDEX) = SM%MASS_FRACTION(N2_INDEX) / MW_N2
SM%MW = 1._EB / SUM(SM%VOLUME_FRACTION)
SM%RSUM = R0/SM%MW
SM%VOLUME_FRACTION = SM%VOLUME_FRACTION / SUM(SM%VOLUME_FRACTION)

!Set reaction variable constants
REACTION%Y_O2_LL         = X_O2_LL*MW_O2/ (X_O2_LL*MW_O2+(1._EB-X_O2_LL)*MW_N2)
REACTION%Y_F_LFL         = Y_F_LFL
REACTION%CRIT_FLAME_TMP  = CRITICAL_FLAME_TEMPERATURE + TMPM
REACTION%NAME            = ID

NU = 0._EB

SET_THREE_PARAMETER: IF (CO_PRODUCTION) THEN !Set reaction variables for three parameter mixture fraction 
!Set reaction variables for the incomplete reaction
   RN => REACTION(1)
   RN%MODE = MIXING_CONTROLLED
   RN%FUEL_INDEX = 1
   RN%IDEAL      = IDEAL
   IF (RN%IDEAL) THEN
!Compute fuel heat of formation
      NU(O2_INDEX) = C + 0.5_EB * H - 0.5_EB * O
      IF (HEAT_OF_COMBUSTION < 0._EB) HEAT_OF_COMBUSTION = EPUMO2*MW_O2*NU(O2_INDEX)/MW_FUEL
   ENDIF

   NU(N2_INDEX) = 0.5_EB * N   
   NU(CO_INDEX) = C
   NU(H2O_INDEX) = 0.5_EB * H
   IF (H2_YIELD > 0._EB) THEN
      NU(H2_INDEX)  = H2_YIELD * MW_FUEL / MW_H2
      NU(H2O_INDEX) = NU(H2O_INDEX) - NU(H2_INDEX)
   ENDIF
   IF (SOOT_YIELD > 0._EB) THEN
      NU(SOOT_INDEX) = SOOT_YIELD * MW_FUEL / MW_SOOT
      NU(CO_INDEX)  = NU(CO_INDEX) - NU(SOOT_INDEX) * (1._EB - SOOT_H_FRACTION)
      IF (H > 0._EB) NU(H2O_INDEX) = NU(H2O_INDEX)- 0.5_EB * NU(SOOT_INDEX) * SOOT_H_FRACTION
   ENDIF
   IF (OTHER > 0._EB) NU(OTHER_INDEX) = OTHER
   NU(O2_INDEX)  = 0.5_EB*(NU(CO_INDEX) - O)
   IF (H>0._EB) NU(O2_INDEX) = NU(O2_INDEX) + 0.5*NU(H2O_INDEX)

   IF (H > 0._EB) THEN
      IF (NU(H2O_INDEX) < 0._EB) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Values for SOOT_YIELD, H2_YIELD, and SOOT_H_FRACTION result in negative H2O yield'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDIF
   
   SM => SPECIES_MIXTURE(3)
   
   SM%VOLUME_FRACTION(N2_INDEX) = NU(N2_INDEX) +& !N2 from fuel molecule
                    SPECIES_MIXTURE(1)%VOLUME_FRACTION(N2_INDEX)*NU(O2_INDEX)/SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)+&!N2 air
                    SPECIES_MIXTURE(2)%VOLUME_FRACTION(N2_INDEX)/SPECIES_MIXTURE(2)%VOLUME_FRACTION(FUEL_INDEX) !N2 fuel stream
   IF (SPECIES_MIXTURE(2)%VOLUME_FRACTION(2)>0)  SM%VOLUME_FRACTION(N2_INDEX) =  SM%VOLUME_FRACTION(N2_INDEX) + &
                    SPECIES_MIXTURE(2)%VOLUME_FRACTION(1)/SPECIES_MIXTURE(2)%VOLUME_FRACTION(2) !N2 as fuel stream diluent
   SM%MASS_FRACTION(N2_INDEX) = SM%VOLUME_FRACTION(N2_INDEX) * MW_N2
   SM%VOLUME_FRACTION(CO2_INDEX) = &
                    SPECIES_MIXTURE(1)%VOLUME_FRACTION(CO2_INDEX)*NU(O2_INDEX)/SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)!CO2 air
   SM%MASS_FRACTION(CO2_INDEX) = SM%VOLUME_FRACTION(CO2_INDEX) * MW_CO2   
   SM%VOLUME_FRACTION(CO_INDEX) = NU(CO_INDEX)
   SM%MASS_FRACTION(CO_INDEX) =  SM%VOLUME_FRACTION(CO_INDEX) * MW_CO
   SM%VOLUME_FRACTION(H2O_INDEX) = NU(H2O_INDEX)  +&!H2O from fuel
                    SPECIES_MIXTURE(1)%VOLUME_FRACTION(H2O_INDEX)*NU(O2_INDEX)/SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)!H2O air
   SM%MASS_FRACTION(H2O_INDEX) = SM%VOLUME_FRACTION(H2O_INDEX) * MW_H2O
   IF (SOOT_YIELD > 0._EB)  THEN
      SM%VOLUME_FRACTION(SOOT_INDEX) = NU(SOOT_INDEX)
      SM%MASS_FRACTION(SOOT_INDEX) = SM%MASS_FRACTION(SOOT_INDEX) * MW_SOOT      
   ENDIF
   IF (H2_YIELD > 0._EB)  THEN
      SM%VOLUME_FRACTION(H2_INDEX) = NU(H2_INDEX)
      SM%MASS_FRACTION(H2_INDEX) = SM%VOLUME_FRACTION(H2_INDEX) * MW_H2
   ENDIF
   IF (OTHER > 0._EB)  THEN
      SM%VOLUME_FRACTION(OTHER_INDEX) = NU(OTHER_INDEX)
      SM%MASS_FRACTION(OTHER_INDEX) = SM%VOLUME_FRACTION(OTHER_INDEX) * MW_OTHER
   ENDIF
   SUM1 = SUM(SM%VOLUME_FRACTION)
   SM%MW = SUM(SM%MASS_FRACTION)/SUM1
   SM%RSUM = R0/SM%MW
   SM%MASS_FRACTION = SM%MASS_FRACTION / SUM(SM%MASS_FRACTION)
   SM%VOLUME_FRACTION = SM%VOLUME_FRACTION / SUM1
   
   RN%N_SPECIES = 2
   ALLOCATE(RN%SPEC_ID(RN%N_SPECIES))
   ALLOCATE(RN%SPECIES(RN%N_SPECIES))   
   ALLOCATE(RN%NU(0:RN%N_SPECIES))
   RN%NAME = 'Fuel+Air=>CO+Products'

   RN%SPEC_ID(1) = 'LUMPED SPECIES 1'   
   RN%SPEC_ID(2) = 'LUMPED SPECIES 2'
   RN%SPECIES(1) = 2
   RN%NU(1) = -1._EB
   RN%NU(0) = -NU(O2_INDEX)/SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)*SPECIES_MIXTURE(2)%VOLUME_FRACTION(FUEL_INDEX) 
   RN%SPECIES(2) = 3
   RN%NU(2) = SUM1 * SPECIES_MIXTURE(2)%VOLUME_FRACTION(FUEL_INDEX) 

   IF (.NOT. RN%IDEAL) THEN
      IF(HEAT_OF_COMBUSTION < 0._EB) HEAT_OF_COMBUSTION = EPUMO2 * (NU(O2_INDEX)+0.5_EB*(NU(CO_INDEX)-CO_YIELD*MW_FUEL/MW_CO)) *&
                                                          MW_O2/MW_FUEL

   ELSE
! Correct heat of combustion for minor products of combustion
      HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION - (C-NU(CO_INDEX)) * CO2_HEAT_OF_FORMATION * MW_CO2/MW_FUEL
      IF (H > 0._EB)  HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION - (0.5_EB*H-NU(H2O_INDEX)) * H2O_HEAT_OF_FORMATION * MW_H2O/MW_FUEL
   ENDIF
! Correct heat of combustion for incomplete reaction   
   RN%HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION*1000._EB - NU(CO_INDEX)*(CO2_HEAT_OF_FORMATION-CO_HEAT_OF_FORMATION)*1.E3_EB/MW_FUEL

! Set reaction variables for complete reaction
!*****
   RN => REACTION(2)
   RN%THRESHOLD_TEMP = THRESHOLD_TEMPERATURE + TMPM
   RN%MODE = MIXED
   RN%FUEL_INDEX = 2
   NU(CO2_INDEX) = NU(CO_INDEX) - CO_YIELD * MW_FUEL/MW_CO
   NU(CO_INDEX)  = CO_YIELD * MW_FUEL/MW_CO
   IF (NU(CO2_INDEX) < 0._EB) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: Values for SOOT_YIELD, CO_YIELD, and SOOT_H_FRACTION result in negative CO2 yield'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   NU(O2_INDEX) = NU(O2_INDEX) + 0.5_EB * NU(CO2_INDEX)

   SM => SPECIES_MIXTURE(4)

   SM%VOLUME_FRACTION(N2_INDEX) = NU(N2_INDEX) +& !N2 from fuel molecule
                    SPECIES_MIXTURE(1)%VOLUME_FRACTION(N2_INDEX)*NU(O2_INDEX)/SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)+&!N2 air
                    SPECIES_MIXTURE(2)%VOLUME_FRACTION(N2_INDEX)/SPECIES_MIXTURE(2)%VOLUME_FRACTION(FUEL_INDEX) !N2 fuel stream
   SM%MASS_FRACTION(N2_INDEX) = SM%VOLUME_FRACTION(N2_INDEX) * MW_N2
   SM%VOLUME_FRACTION(CO2_INDEX) = NU(CO2_INDEX) +&!CO2 from fuel
                    SPECIES_MIXTURE(1)%VOLUME_FRACTION(CO2_INDEX)*NU(O2_INDEX)/SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)!CO2 air
   SM%MASS_FRACTION(CO2_INDEX) =SM%VOLUME_FRACTION(CO2_INDEX) * MW_CO2
   SM%VOLUME_FRACTION(CO_INDEX) = NU(CO_INDEX)
   SM%MASS_FRACTION(CO_INDEX) = SM%VOLUME_FRACTION(CO_INDEX) * MW_CO
   SM%VOLUME_FRACTION(H2O_INDEX) = NU(H2O_INDEX) +&!H2O from fuel
                    SPECIES_MIXTURE(1)%VOLUME_FRACTION(H2O_INDEX)*NU(O2_INDEX)/SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)!H2O air
   SM%MASS_FRACTION(H2O_INDEX) = SM%VOLUME_FRACTION(H2O_INDEX) * MW_H2O
   IF (SOOT_YIELD > 0._EB)  THEN
      SM%VOLUME_FRACTION(SOOT_INDEX) = NU(SOOT_INDEX)
      SM%MASS_FRACTION(SOOT_INDEX) = SM%VOLUME_FRACTION(SOOT_INDEX) * MW_SOOT      
   ENDIF
   IF (H2_YIELD > 0._EB)  THEN
      SM%VOLUME_FRACTION(H2_INDEX) = NU(H2_INDEX)
      SM%MASS_FRACTION(H2_INDEX) = SM%VOLUME_FRACTION(H2_INDEX) * MW_H2
   ENDIF
   IF (OTHER > 0._EB)  THEN
      SM%VOLUME_FRACTION(OTHER_INDEX) = NU(OTHER_INDEX)
      SM%MASS_FRACTION(OTHER_INDEX) = SM%VOLUME_FRACTION(OTHER_INDEX) * MW_OTHER
   ENDIF
   SUM1 = SUM(SM%VOLUME_FRACTION)
   SM%MW = SUM(SM%MASS_FRACTION)/SUM1
   SM%RSUM = R0/SM%MW
   SM%MASS_FRACTION = SM%MASS_FRACTION / SUM(SM%MASS_FRACTION)
   SM%VOLUME_FRACTION = SM%VOLUME_FRACTION / SUM1
   
   RN%N_SPECIES = 2
   ALLOCATE(RN%SPEC_ID(RN%N_SPECIES))
   ALLOCATE(RN%SPECIES(RN%N_SPECIES))   
   ALLOCATE(RN%NU(0:RN%N_SPECIES))
   RN%NAME = 'CO+Products+Air=>CO2+Products'
   RN%SPEC_ID(1) = 'LUMPED SPECIES 2'
   RN%SPEC_ID(2) = 'LUMPED SPECIES 3'
   
   RN%SPECIES(1) = 3
   RN%NU(1) = -1._EB
   RN%NU(0) = -0.5_EB*SPECIES_MIXTURE(3)%VOLUME_FRACTION(CO_INDEX)*NU(CO2_INDEX)/(NU(CO2_INDEX)+NU(CO_INDEX))/&
               SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)
   RN%SPECIES(2) = 4
   RN%NU(2) = 1._EB-RN%NU(0)*(1._EB-SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX))

!*****
   RN%BOF = 2.53E12_EB
   RN%E = 199547._EB*1000._EB
   RN%HEAT_OF_COMBUSTION = (CO2_HEAT_OF_FORMATION - CO_HEAT_OF_FORMATION)*1000._EB/MW_CO

ELSE SET_THREE_PARAMETER !Set reaction variables for two parameter mixture fraction 
!Set reaction variables for complete reaction
   RN => REACTION(1)    
   RN%MODE = MIXING_CONTROLLED  
   RN%IDEAL      = IDEAL
   RN%FUEL_INDEX = 1   
   IF (RN%IDEAL) THEN
!Compute fuel heat of formation
      NU(O2_INDEX) = C + 0.5_EB * H - 0.5_EB * O
      IF (HEAT_OF_COMBUSTION < 0._EB) HEAT_OF_COMBUSTION = EPUMO2*(MW_O2*NU(O2_INDEX))/MW_FUEL
   ENDIF
   
   NU(CO2_INDEX) = C
   IF (CO_YIELD > 0._EB) THEN
      NU(CO_INDEX)  = CO_YIELD * MW_FUEL / MW_CO
      NU(CO2_INDEX) = NU(CO2_INDEX) - NU(CO_INDEX)
   ENDIF
   IF (H > 0._EB) NU(H2O_INDEX) = 0.5_EB * H
   IF (H2_YIELD > 0._EB) THEN
      NU(H2_INDEX)  = H2_YIELD * MW_FUEL / MW_H2
      NU(H2O_INDEX) = NU(H2O_INDEX) - NU(H2_INDEX)
   ENDIF
   IF (SOOT_YIELD > 0._EB) THEN
      NU(SOOT_INDEX) = SOOT_YIELD * MW_FUEL / MW_SOOT
      NU(CO2_INDEX)  = NU(CO2_INDEX) - NU(SOOT_INDEX) * (1._EB - SOOT_H_FRACTION)
      IF (H > 0._EB) NU(H2O_INDEX) = NU(H2O_INDEX)- 0.5_EB * NU(SOOT_INDEX) * SOOT_H_FRACTION
   ENDIF
   IF (OTHER > 0) NU(OTHER_INDEX) = OTHER
   IF (NU(CO2_INDEX) < 0._EB) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: Values for SOOT_YIELD, CO_YIELD, and SOOT_H_FRACTION result in negative CO2 yield'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (H > 0._EB) THEN
      IF (NU(H2O_INDEX) < 0._EB) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Values for SOOT_YIELD, H2_YIELD, and SOOT_H_FRACTION result in negative H2O yield'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDIF
   NU(O2_INDEX)  = NU(CO2_INDEX) - 0.5_EB * O
   IF (CO_YIELD > 0._EB) NU(O2_INDEX) = NU(O2_INDEX) + 0.5_EB * NU(CO_INDEX)
   IF (H>0._EB) NU(O2_INDEX) = NU(O2_INDEX) + 0.5*NU(H2O_INDEX)
   NU(N2_INDEX) = 0.5_EB * N
   IF (OTHER > 0._EB) NU(OTHER_INDEX) = OTHER
   
   SM => SPECIES_MIXTURE(3)
   SM%VOLUME_FRACTION(N2_INDEX) = NU(N2_INDEX) +& !N2 from fuel molecule
                    SPECIES_MIXTURE(1)%VOLUME_FRACTION(N2_INDEX)*NU(O2_INDEX)/SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)+&!N2 air
                    SPECIES_MIXTURE(2)%VOLUME_FRACTION(N2_INDEX)/SPECIES_MIXTURE(2)%VOLUME_FRACTION(FUEL_INDEX) !N2 fuel stream
   SM%MASS_FRACTION(N2_INDEX) =  SM%VOLUME_FRACTION(N2_INDEX) * MW_N2
   SM%VOLUME_FRACTION(CO2_INDEX) = NU(CO2_INDEX) + &!CO2 from fuel
                    SPECIES_MIXTURE(1)%VOLUME_FRACTION(CO2_INDEX)*NU(O2_INDEX)/SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)!CO2 air
   SM%MASS_FRACTION(CO2_INDEX) = SM%VOLUME_FRACTION(CO2_INDEX) * MW_CO2 
   SM%VOLUME_FRACTION(H2O_INDEX) = NU(H2O_INDEX) + & !H2O from fuel
                    SPECIES_MIXTURE(1)%VOLUME_FRACTION(H2O_INDEX)*NU(O2_INDEX)/SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)!H2O air
   SM%MASS_FRACTION(H2O_INDEX) =  SM%VOLUME_FRACTION(H2O_INDEX) * MW_H2O
   IF (CO_YIELD > 0._EB) THEN
      SM%VOLUME_FRACTION(CO_INDEX) = NU(CO_INDEX)
      SM%MASS_FRACTION(CO_INDEX) = SM%VOLUME_FRACTION(CO_INDEX) * MW_CO
   ENDIF
   IF (SOOT_YIELD > 0._EB)  THEN
      IF (SOOT_DEPOSITION) THEN
         SM%VOLUME_FRACTION(SOOT_INDEX) = 0._EB
         SM%MASS_FRACTION(SOOT_INDEX) = 0._EB
         SPECIES_MIXTURE(4)%VOLUME_FRACTION(SOOT_INDEX) = 1._EB
         SPECIES_MIXTURE(4)%MASS_FRACTION(SOOT_INDEX) = 1._EB
         SPECIES_MIXTURE(4)%MW = MW_SOOT
         SPECIES_MIXTURE(4)%RSUM = R0/MW_SOOT
      ELSE
         SM%VOLUME_FRACTION(SOOT_INDEX) = NU(SOOT_INDEX)
         SM%MASS_FRACTION(SOOT_INDEX) = SM%VOLUME_FRACTION(SOOT_INDEX) * MW_SOOT      
      ENDIF
   ENDIF
   IF (H2_YIELD > 0._EB)  THEN
      SM%VOLUME_FRACTION(H2_INDEX) = NU(H2_INDEX)
      SM%MASS_FRACTION(H2_INDEX) = SM%VOLUME_FRACTION(H2_INDEX) * MW_H2
   ENDIF
   IF (OTHER > 0._EB)  THEN
      SM%VOLUME_FRACTION(OTHER_INDEX) = NU(OTHER_INDEX)
      SM%MASS_FRACTION(OTHER_INDEX) = SM%VOLUME_FRACTION(OTHER_INDEX) * MW_OTHER
   ENDIF
   SUM1 = SUM(SM%VOLUME_FRACTION)
   SM%MW = SUM(SM%MASS_FRACTION)/SUM1
   SM%RSUM = R0/SM%MW
   SM%MASS_FRACTION = SM%MASS_FRACTION / SUM(SM%MASS_FRACTION)
   SM%VOLUME_FRACTION = SM%VOLUME_FRACTION / SUM1
   IF (SOOT_DEPOSITION) THEN
      RN%N_SPECIES=3
      ALLOCATE(RN%SPEC_ID(RN%N_SPECIES))
      ALLOCATE(RN%SPECIES(RN%N_SPECIES))   
      ALLOCATE(RN%NU(0:RN%N_SPECIES))
      RN%NAME = 'Fuel+Air=>CO2+Soot+Products'
      RN%SPEC_ID(1) = 'LUMPED SPECIES 1'
      RN%SPEC_ID(2) = 'LUMPED SPECIES2'
      RN%SPEC_ID(3) = 'LUMPED SPECIES 3'      
      
      RN%SPECIES(1) = 2
      RN%NU(1) = -1._EB
      RN%NU(0) = -NU(O2_INDEX)/SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)*SPECIES_MIXTURE(2)%VOLUME_FRACTION(FUEL_INDEX)
      RN%SPECIES(2) = 3
      RN%NU(2) = 1._EB
      RN%SPECIES(3) = 4
      RN%NU(3) = SOOT_YIELD*MW_FUEL/MW_SOOT
   ELSE
      RN%N_SPECIES=2
      ALLOCATE(RN%SPEC_ID(RN%N_SPECIES))
      ALLOCATE(RN%SPECIES(RN%N_SPECIES))   
      ALLOCATE(RN%NU(0:RN%N_SPECIES))   
      RN%NAME = 'Fuel+Air=>Products'            
      RN%SPEC_ID(1) = 'LUMPED SPECIES 1'
      RN%SPEC_ID(2) = 'LUMPED SPECIES 2'
      
      RN%SPECIES(1) = 2
      RN%NU(1) = -1._EB
      RN%NU(0) = -NU(O2_INDEX)/SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)*SPECIES_MIXTURE(2)%VOLUME_FRACTION(FUEL_INDEX)
      RN%SPECIES(2) = 3
      RN%NU(2) = SUM1
   ENDIF
   IF (.NOT. RN%IDEAL) THEN
      IF(HEAT_OF_COMBUSTION < 0._EB) HEAT_OF_COMBUSTION = &
                                     -EPUMO2*RN%NU(0)*MW_O2*SPECIES_MIXTURE(1)%VOLUME_FRACTION(O2_INDEX)/&
                                     (SPECIES_MIXTURE(2)%VOLUME_FRACTION(FUEL_INDEX)*MW_FUEL)
   ELSE
! Correct heat of combustion for minor products of combustion
      IF (SOOT_YIELD > 0._EB) HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION - NU(SOOT_INDEX)*(1._EB-SOOT_H_FRACTION)*&
                                                   CO2_HEAT_OF_FORMATION/MW_FUEL
      IF (CO_YIELD > 0._EB) HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION - NU(CO_INDEX)*&
                                                 (CO2_HEAT_OF_FORMATION - CO_HEAT_OF_FORMATION)/MW_FUEL
      IF (H > 0._EB)  HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION - (0.5_EB*H-NU(H2O_INDEX))*H2O_HEAT_OF_FORMATION/MW_FUEL
   ENDIF

   RN%HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION*1000._EB

ENDIF SET_THREE_PARAMETER 

! Set the lower limit of the extinction coefficient

EC_LL = VISIBILITY_FACTOR/MAXIMUM_VISIBILITY

! Change units of combustion quantities

HRRPUV_AVERAGE = HRRPUV_AVERAGE*1000._EB   ! W/m3
HRRPUA_SHEET   = HRRPUA_SHEET*  1000._EB   ! W/m2

CONTAINS

SUBROUTINE SET_REAC_DEFAULTS_1

AUTO_IGNITION_TEMPERATURE   = -TMPM
BOF                         = -1._EB       ! cm**3/mol-s
C                           = -1._EB
CO_YIELD                    = 0._EB
CRITICAL_FLAME_TEMPERATURE  = 1427._EB    ! C
E                           = -1._EB       ! kJ/kmol
EPUMO2                      = 13100._EB   ! kJ/kg
FUEL                        = 'null'
FYI                         = 'null'
H                           = -1._EB
H2_YIELD                    = 0._EB
HEAT_OF_COMBUSTION          = -1._EB
IF (LES) HRRPUA_SHEET       = 0._EB
IF (DNS) HRRPUA_SHEET       = 200._EB
IF (LES) HRRPUV_AVERAGE     = 2500._EB
IF (DNS) HRRPUV_AVERAGE     = 0._EB
IF (FDS6) THEN
   HRRPUV_AVERAGE       = 2500._EB
   HRRPUA_SHEET         = 200._EB
ENDIF
ID                          = 'null'
IDEAL                       = .FALSE.
MAXIMUM_VISIBILITY          = 30._EB       ! m
MW_OTHER                    = 28._EB        ! Nitrogen MW
N                           = 0._EB
NU                          = 0._EB
N_S                         = -999._EB
O                           = 0._EB
OTHER                       = 0._EB
SOOT_H_FRACTION             = 0.1_EB
IF (LES) SOOT_YIELD         = 0.0_EB
IF (DNS) SOOT_YIELD         = 0.0_EB
SPEC_ID                     = 'null'
THRESHOLD_TEMPERATURE       = 300._EB      ! For SIMPLE_CHEMSITRY w/ CO_PRODUCTION 
X_O2_LL                     = 0.15_EB    ! %
Y_F_LFL                     = 0.0_EB 
Y_F_INLET                   = 1._EB
VISIBILITY_FACTOR           = 3._EB
 
END SUBROUTINE SET_REAC_DEFAULTS_1
 
 
END SUBROUTINE READ_REAC_1


SUBROUTINE READ_REAC_2
 
CHARACTER(30) :: FUEL,SPEC_ID(1:MAX_SPECIES)
LOGICAL :: IDEAL
INTEGER :: NR,N_REAC_READ,NS,NS2
REAL(EB) :: Y_F_INLET, MASS_REACTANT, MASS_PRODUCT, &
            H2_YIELD,SOOT_YIELD,CO_YIELD, Y_F_LFL,X_O2_LL,EPUMO2,BOF,&
            CRITICAL_FLAME_TEMPERATURE,HEAT_OF_COMBUSTION,MASS(MAX_SPECIES),NU(0:MAX_SPECIES),E,N_S(0:MAX_SPECIES),C,H,N,O,OTHER, &
            FUEL_HEAT_OF_FORMATION,MAXIMUM_VISIBILITY,AUTO_IGNITION_TEMPERATURE,SUM1,THRESHOLD_TEMPERATURE
TYPE (SPECIES_MIXTURE_TYPE), POINTER :: SM
NAMELIST /REAC/ E,BOF,HEAT_OF_COMBUSTION,FYI,FUEL,EPUMO2,ID, N_S,EXTINCTION2,&
                Y_O2_INFTY,Y_F_INLET, CO_PRODUCTION,SOOT_DEPOSITION,SUPPRESSION,SUPPRESSION_SEARCH,&
                H2_YIELD,SOOT_YIELD,CO_YIELD,Y_F_LFL,X_O2_LL,CRITICAL_FLAME_TEMPERATURE,NU,SOOT_H_FRACTION, &
                C,H,N,O,OTHER,MW_OTHER,IDEAL,VISIBILITY_FACTOR,MAXIMUM_VISIBILITY, &
                BETA_EDC,Y_P_MIN_EDC,HRRPUA_SHEET,HRRPUV_AVERAGE,FIXED_MIX_TIME, &
                TAU_CHEM,TAU_FLAME,AUTO_IGNITION_TEMPERATURE, &
                MASS,SPEC_ID,THRESHOLD_TEMPERATURE
                
IF (SIMPLE_CHEMISTRY) THEN
   DO NR=1,N_REACTIONS
      RN => REACTION(NR)
      DO NS=1,RN%N_SPECIES
         RN%SPECIES(NS) = N_GAS_SPECIES - N_MIX_SPECIES + RN%SPECIES(NS)
      ENDDO
      ALLOCATE (RN%NU_SPECIES(1:N_SPECIES-N_MIX_SPECIES))
      CALL GET_RN_NU_SPECIES(NR)
   ENDDO
   RETURN
ENDIF 

N_REACTIONS = 0
REWIND(LU_INPUT)
 
COUNT_REAC_LOOP: DO
   CALL SET_REAC_DEFAULTS_2
   CALL CHECKREAD('REAC',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_REAC_LOOP
   READ(LU_INPUT,REAC,ERR=434,IOSTAT=IOS)
   N_REACTIONS = N_REACTIONS + 1
   434 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with REAC ',N_REACTIONS+1
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO COUNT_REAC_LOOP

IF (N_REACTIONS==0) RETURN
 
N_REAC_READ = N_REACTIONS

ALLOCATE(REACTION(N_REACTIONS),STAT=IZERO)
CALL ChkMemErr('READ','REACTION',IZERO)

! Read the input file looking for REAC lines
REWIND(LU_INPUT)
READ_REACTION_LOOP: DO NR=1,N_REAC_READ
   RN => REACTION(NR)
   CALL CHECKREAD('REAC',LU_INPUT,IOS) 
   CALL SET_REAC_DEFAULTS_2
   IF (IOS==0) READ(LU_INPUT,REAC)
   RN%AUTO_IGNITION_TEMPERATURE = AUTO_IGNITION_TEMPERATURE + TMPM
   RN%BOF                = BOF
   RN%E                  = E*1000._EB
   RN%FYI                = FYI
   RN%HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION*1000._EB
   RN%THRESHOLD_TEMP     = THRESHOLD_TEMPERATURE + TMPM
   RN%NAME               = ID
   IF (TRIM(FUEL)=='null') THEN
      WRITE(MESSAGE,'(A,I3,A)') 'ERROR: Problem with REAC ',NR,'. FUEL must be defined'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   RN%FUEL               = FUEL
   IF (BOF > 0._EB) THEN
      IF (THRESHOLD_TEMPERATURE > TMPM) THEN
         RN%MODE = MIXED
      ELSE
         RN%MODE = FINITE_RATE
      ENDIF
   ELSE
     RN%MODE = MIXING_CONTROLLED
   ENDIF
   
   IF (TRIM(SPEC_ID(1))=='null') THEN
      WRITE(MESSAGE,'(A,I3,A)') 'ERROR: Problem with REAC ',NR,'. SPEC_ID must be defined for finite rate reactions'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   NS2 = 0
   DO NS=1,MAX_SPECIES
      IF (TRIM(SPEC_ID(NS))/='null') THEN 
         NS2=NS2+1
      ELSE
         EXIT
      ENDIF
   ENDDO
   ALLOCATE(RN%N(0:NS2))
   RN%N = N_S(0:NS2)
   ALLOCATE(RN%NU(0:NS2))
   RN%NU=NU(0:NS2)
   ALLOCATE(RN%NU_SPECIES(1:N_SPECIES-N_MIX_SPECIES))      
   RN%NU_SPECIES=0._EB
   RN%N_SPECIES=NS2
   ALLOCATE(RN%SPECIES(RN%N_SPECIES)) 
   ALLOCATE(RN%SPEC_ID(0:N_SPECIES))
   RN%SPEC_ID(0)=SPECIES(0)%ID
   RN%SPEC_ID(1:RN%N_SPECIES)=SPEC_ID(1:RN%N_SPECIES)
   MASS_REACTANT = 0._EB
   MASS_PRODUCT = 0._EB   
   IF (RN%NU(0) < -ZERO_P) MASS_REACTANT = MASS_REACTANT + RN%NU(0)*SPECIES(Y2SPEC(RN%SPECIES(NS)))%MW
   IF (RN%NU(0) > ZERO_P)  MASS_PRODUCT = MASS_PRODUCT + RN%NU(0)*SPECIES(Y2SPEC(RN%SPECIES(NS)))%MW

   DO NS=1,RN%N_SPECIES
      IF (TRIM(SPEC_ID(NS))==TRIM(SPECIES(I_BACKGROUND)%ID)) THEN
         WRITE(MESSAGE,'(A,I3,A)') 'ERROR: Problem with REAC ',NR,&
                                   '. Background species is N_S(0) or NU(0).  Do not list in SPEC_ID' 
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      DO NS2=1,N_GAS_SPECIES
         IF (TRIM(SPEC_ID(NS))==TRIM(SPECIES(Y2SPEC(NS2))%ID)) THEN
            RN%SPECIES(NS)=NS2
            EXIT
         ENDIF
      ENDDO
      IF (RN%SPECIES(NS)==-1) THEN
         WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR: Problem with REAC ',NR,'. Species ',TRIM(SPEC_ID(NS)),&
         ' not specified or is LUMPED_SPECIES_ONLY' 
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      IF (RN%NU(NS) < -ZERO_P) MASS_REACTANT = MASS_REACTANT + RN%NU(NS)*SPECIES(Y2SPEC(RN%SPECIES(NS)))%MW
      IF (RN%NU(NS) > ZERO_P)  MASS_PRODUCT = MASS_PRODUCT + RN%NU(NS)*SPECIES(Y2SPEC(RN%SPECIES(NS)))%MW
   ENDDO  
   IF (MASS_PRODUCT > ZERO_P .AND. MASS_REACTANT < -ZERO_P) THEN
      IF ((MASS_PRODUCT+MASS_REACTANT)/MASS_PRODUCT > 1.E-6_EB) THEN
         WRITE(MESSAGE,'(A,I3,A)') 'ERROR: Problem with REAC ',NR,'. Mass of products does not equal mass of reactants.' 
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDIF
   RN%FUEL_INDEX = -1
   IF (I_FUEL == I_BACKGROUND) THEN
      RN%FUEL_INDEX = 0
   ELSE
      DO NS=1,RN%N_SPECIES
         IF (TRIM(SPEC_ID(NS))==TRIM(FUEL)) THEN
            RN%FUEL_INDEX= RN%SPECIES(NS)
            IF (RN%NU(NS)>0._EB) THEN
               WRITE(MESSAGE,'(A,I3,A,A,A)') 'ERROR: Problem with REAC ',NR,'. Fuel species, ',TRIM(FUEL),&
               ', has NU > 0.' 
               CALL SHUTDOWN(MESSAGE)
            ENDIF                      
            RN%NU = -RN%NU/RN%NU(NS)
            EXIT
         ENDIF
      ENDDO
   ENDIF
   IF (RN%FUEL_INDEX==-1) THEN
      WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR: Problem with REAC ',NR,'. Fuel species ',TRIM(FUEL),&
      ' not specified on REAC line' 
      CALL SHUTDOWN(MESSAGE)
   ENDIF      
   CALL GET_RN_NU_SPECIES(NR)     
ENDDO READ_REACTION_LOOP

! Set the lower limit of the extinction coefficient
EC_LL = VISIBILITY_FACTOR/MAXIMUM_VISIBILITY

! Change units of combustion quantities

HRRPUV_AVERAGE = HRRPUV_AVERAGE*1000._EB   ! W/m3
HRRPUA_SHEET   = HRRPUA_SHEET*  1000._EB   ! W/m2

CONTAINS
 

SUBROUTINE SET_REAC_DEFAULTS_2

AUTO_IGNITION_TEMPERATURE   = -TMPM        ! C 
BOF                         = -1._EB       ! cm**3/mol-s
CRITICAL_FLAME_TEMPERATURE  = 1427._EB    ! C
E                           = -1._EB       ! kJ/kmol
FUEL                        = 'null'
FYI                         = 'null'
HEAT_OF_COMBUSTION          = -1._EB
IF (LES) HRRPUA_SHEET       = 0._EB
IF (DNS) HRRPUA_SHEET       = 200._EB
IF (LES) HRRPUV_AVERAGE     = 2500._EB
IF (DNS) HRRPUV_AVERAGE     = 0._EB
IF (FDS6) THEN
   HRRPUV_AVERAGE       = 2500._EB
   HRRPUA_SHEET         = 200._EB
ENDIF
ID                          = 'null'
N_S                         = -999._EB
NU                          = 0._EB
SPEC_ID                     = 'null'
MAXIMUM_VISIBILITY          = 30._EB       ! m
THRESHOLD_TEMPERATURE       = -TMPM
VISIBILITY_FACTOR           = 3._EB
 
END SUBROUTINE SET_REAC_DEFAULTS_2

SUBROUTINE GET_RN_NU_SPECIES(NR)
INTEGER, INTENT(IN):: NR
INTEGER :: NS,NS2
TYPE (REACTION_TYPE),POINTER :: RN
TYPE (SPECIES_MIXTURE_TYPE), POINTER :: SM
 
RN => REACTION(NR)
RN%NU_SPECIES = 0._EB

IF (ABS(RN%NU(0))> ZERO_P) THEN
   IF (I_BACKGROUND >= N_SPECIES - N_MIX_SPECIES) THEN
      SM => SPECIES_MIXTURE(1)
      DO NS2 = 1,SM%N_SUB_SPECIES
         RN%NU_SPECIES(SM%SPEC_INDEX(NS2)) =  RN%NU_SPECIES(SM%SPEC_INDEX(NS2)) + RN%NU(0)*SM%VOLUME_FRACTION(NS2)
      ENDDO
   ELSE
      RN%NU_SPECIES(I_BACKGROUND) = RN%NU(0)
   ENDIF
ENDIF

DO NS=1,RN%N_SPECIES
   IF (Y2SPEC(RN%SPECIES(NS)) >= N_SPECIES - N_MIX_SPECIES) THEN
      SM => SPECIES_MIXTURE(Y2SPEC(RN%SPECIES(NS)) - N_SPECIES + N_MIX_SPECIES)
      DO NS2 = 1,SM%N_SUB_SPECIES
         RN%NU_SPECIES(SM%SPEC_INDEX(NS2)) =  RN%NU_SPECIES(SM%SPEC_INDEX(NS2)) + RN%NU(NS)*SM%VOLUME_FRACTION(NS2)
      ENDDO
   ELSE
      RN%NU_SPECIES(Y2SPEC(RN%SPECIES(NS))) = RN%NU_SPECIES(Y2SPEC(RN%SPECIES(NS))) + RN%NU(NS)
   ENDIF
ENDDO

END SUBROUTINE GET_RN_NU_SPECIES 

END SUBROUTINE READ_REAC_2

 
 
SUBROUTINE PROC_SPEC
USE PHYSICAL_FUNCTIONS, ONLY : GET_MOLECULAR_WEIGHT, GET_SPECIFIC_GAS_CONSTANT
USE PROPERTY_DATA, ONLY: JANAF_TABLE, CALC_GAS_PROPS, GAS_PROPS
REAL(EB) :: YYY(1:N_GAS_SPECIES)
REAL(EB), ALLOCATABLE, DIMENSION(:) :: MU_TMP,CP_TMP,K_TMP,H_TMP,D_TMP
INTEGER :: IPC,NN,N,N2
CHARACTER(30) :: FORMULA
TYPE(PARTICLE_CLASS_TYPE), POINTER :: PC=>NULL()
TYPE(SPECIES_MIXTURE_TYPE), POINTER :: SM=>NULL()

! Set lookup array
ALLOCATE(Y2SPEC(1:N_GAS_SPECIES))
NN = 0
DO N=1,N_SPECIES
   IF (SPECIES(N)%Y_INDEX > 0) THEN
      NN = NN + 1
      Y2SPEC(NN) = N
   ENDIF
ENDDO

! Get gas properties
N_Y_ARRAY = N_SPECIES

IF (SIMPLE_CHEMISTRY) N_Y_ARRAY = N_Y_ARRAY - N_MIX_SPECIES

ALLOCATE(Y2Y(N_Y_ARRAY,N_GAS_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ','Y2Y',IZERO) 
Y2Y = 0._EB 
ALLOCATE(Y2Y_C(N_Y_ARRAY),STAT=IZERO)
CALL ChkMemErr('READ','Y2Y_C',IZERO) 
Y2Y_C = 0._EB 
ALLOCATE(Y_MF_SUM_Y(N_Y_ARRAY),STAT=IZERO)
CALL ChkMemErr('READ','Y_MF_SUM_Y',IZERO) 
Y_MF_SUM_Y = 0._EB 

IF (N_MIX_SPECIES > 0 .AND. I_BACKGROUND > N_SPECIES - N_MIX_SPECIES) THEN
   SM =>SPECIES_MIXTURE(1)
   DO NN = 1,SM%N_SUB_SPECIES
      Y2Y_C(SM%SPEC_INDEX(NN))=SM%MASS_FRACTION(NN)
   ENDDO
ELSE
   Y2Y_C(I_BACKGROUND)=1._EB
ENDIF

! Compute Y to Y Arrays
DO N=1,N_SPECIES-N_MIX_SPECIES   
   IF (SPECIES(N)%Y_INDEX > 0) THEN
      Y2Y(N,SPECIES(N)%Y_INDEX) = 1._EB
   ENDIF
ENDDO

IF (N_MIX_SPECIES > 0) THEN
   N2 = N_SPECIES-N_MIX_SPECIES+1
   IF (DEFINED_BACKGROUND) N2 = N2 + 1
   DO N=N2,N_SPECIES
      SM =>SPECIES_MIXTURE(N-N_SPECIES+N_MIX_SPECIES)
      DO NN = 1,SM%N_SUB_SPECIES
         Y2Y(SM%SPEC_INDEX(NN),SPECIES(N)%Y_INDEX)=SM%MASS_FRACTION(NN)
      ENDDO
   ENDDO
ENDIF

DO N=1,N_GAS_SPECIES
   Y2Y(:,N)=Y2Y(:,N)-Y2Y_C(:)
ENDDO
Y_MF_SUM_Y_C = SUM(Y2Y_C)

ALLOCATE(MW_AVG_Y(N_GAS_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ','MW_AVG_Y',IZERO) 

MW_AVG_Y = 0._EB
MW_AVG_Y_C = 0._EB
IF (N_MIX_SPECIES > 0 .AND. DEFINED_BACKGROUND) THEN
   MW_AVG_Y_C =  1._EB/SPECIES_MIXTURE(1)%MW
ELSE
   MW_AVG_Y_C =  1._EB/SPECIES(0)%MW
ENDIF

DO N=1,N_SPECIES
   IF (SPECIES(N)%Y_INDEX > 0) MW_AVG_Y(SPECIES(N)%Y_INDEX) = 1._EB/SPECIES(N)%MW - MW_AVG_Y_C
ENDDO

! Compute the initial value of R0/MW_AVG

SS0 => SPECIES(0)

ALLOCATE(YY_AMBIENT(N_GAS_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ','YY_AMBIENT',IZERO) 
YY_AMBIENT = 0._EB
 
SS0%YY0 = 1._EB
DO N=1,N_GAS_SPECIES
   SS => SPECIES(Y2SPEC(N))
   IF (SS%MODE==GAS_SPECIES) SS0%YY0 = SS0%YY0 - SS%YY0
   YY_AMBIENT(N) = SS%YY0
ENDDO
YYY = YY_AMBIENT
CALL GET_SPECIFIC_GAS_CONSTANT(YYY,RSUM0)

IF (N_MIX_SPECIES > 0.AND. DEFINED_BACKGROUND)  SS0%MW = SPECIES_MIXTURE(1)%MW
SPECIES%RCON = R0/SPECIES%MW

MW_MIN = MINVAL(SPECIES%MW)
MW_MAX = MAXVAL(SPECIES%MW)
 
! Compute background density from other background quantities
 
RHOA = P_INF/(TMPA*RSUM0)
 
! Compute constant-temperature specific heats
 
CP_GAMMA = SS0%RCON*GAMMA/(GAMMA-1._EB)
CPOPR = CP_GAMMA/PR

! Compute variable-temperature specific heats for specific species
! Compute viscosity, thermal conductivity for species 0 to N_SPECIES. Diffusivity for species 1 to N_SPECIES. 
! These terms are used for DNS, and as the lower limits in LES calcs

ALLOCATE(D_TMP(N_Y_ARRAY))
D_TMP = 0._EB
ALLOCATE(MU_TMP(N_Y_ARRAY))
MU_TMP = 0._EB
ALLOCATE(CP_TMP(N_Y_ARRAY))
CP_TMP = 0._EB
ALLOCATE(H_TMP(N_Y_ARRAY))
H_TMP = 0._EB
ALLOCATE(K_TMP(N_Y_ARRAY))
K_TMP = 0._EB
D_TMP = 0._EB

ALLOCATE(Y2CP_C(0:5000))
CALL ChkMemErr('READ','Y2CP_C',IZERO)    
Y2CP_C = 0._EB   

ALLOCATE(Y2CPBAR_C(0:5000))   
CALL ChkMemErr('READ','Y2CPBAR_C',IZERO)    
Y2CPBAR_C = 0._EB   

ALLOCATE(Y2H_G_C(0:5000))
CALL ChkMemErr('READ','Y2H_G_C',IZERO)    
Y2H_G_C = 0._EB   

ALLOCATE(Y2K_C(0:5000))
CALL ChkMemErr('READ','Y2K_C',IZERO)    
Y2K_C = 0._EB   

ALLOCATE(Y2MU_C(0:5000))      
CALL ChkMemErr('READ','Y2MU_C',IZERO)    
Y2MU_C = 0._EB   

ALLOCATE(Y2CP(0:5000,N_GAS_SPECIES))     
CALL ChkMemErr('READ','Y2CC',IZERO)    
Y2CP = 0._EB   

ALLOCATE(Y2CPBAR(0:5000,N_GAS_SPECIES))     
CALL ChkMemErr('READ','Y2CPBAR',IZERO)    
Y2CPBAR = 0._EB   

ALLOCATE(Y2H_G(0:5000,N_GAS_SPECIES))           
CALL ChkMemErr('READ','YH_GC',IZERO)    
Y2H_G = 0._EB   

ALLOCATE(Y2K(0:5000,N_GAS_SPECIES))      
CALL ChkMemErr('READ','Y2K',IZERO)    
Y2K = 0._EB   

ALLOCATE(Y2MU(0:5000,N_GAS_SPECIES))      
CALL ChkMemErr('READ','Y2MU',IZERO)    
Y2MU = 0._EB   

IF (N_MIX_SPECIES > 0) THEN
   ALLOCATE(Y2D_C(0:5000))
   CALL ChkMemErr('READ','Y2D_C',IZERO)    
   Y2D_C = 0._EB   
   ALLOCATE(Y2D(0:5000,I_Z_MAX-I_Z_MIN+1))
   CALL ChkMemErr('READ','Y2D_C',IZERO)    
   Y2D = 0._EB   
ENDIF

! Loop through temperatures from 1 K to 5000 K to get temperature-specific gas properties.  Data from JANAF 4
DO J = 1,5000
   DO N = 1,N_SPECIES-N_MIX_SPECIES
      SS => SPECIES(N)
      CALL CALC_GAS_PROPS(J,N,D_TMP(N),MU_TMP(N),K_TMP(N),CP_TMP(N),H_TMP(N),SS%ISFUEL)
      IF (SS%ISFUEL .AND. FUEL_INDEX <=0) FUEL_INDEX = N
      SS%D(J) = D_TMP(N)
   ENDDO
   DO N = 1,N_GAS_SPECIES
      Y2CP(J,N)    = SUM(Y2Y(:,N) * CP_TMP(:))
      IF (J>1) THEN
         Y2H_G(J,N)   = Y2H_G(J-1,N) +0.5_EB*(Y2CP(J,N)+Y2CP(J-1,N)) 
         Y2CPBAR(J,N) = Y2H_G(J,N)/REAL(J,EB)   
      ELSE
         Y2H_G(J,N)   = SUM(Y2Y(:,N) * H_TMP(:)) + Y2CP(J,N) 
         Y2CPBAR(J,N) = Y2H_G(J,N)
      ENDIF
      Y2MU(J,N)    = SUM(Y2Y(:,N) * MU_TMP(:))
      Y2K(J,N)     = SUM(Y2Y(:,N) * K_TMP(:)) 
      IF (SPECIES(N)%MODE==LUMPED_SPECIES) Y2D(J,N-I_Z_MIN+1) = SUM(Y2Y(:,N) * D_TMP(:))
   END DO
   IF (N_MIX_SPECIES > 0) THEN
      Y2CP_C(J)    = SUM(Y2Y_C(:) * CP_TMP(:))
      IF (J>1) THEN
         Y2H_G_C(J)   = Y2H_G_C(J-1) + 0.5_EB*(Y2CP_C(J)+Y2CP_C(J-1)) 
         Y2CPBAR_C(J) = Y2H_G_C(J)/REAL(J,EB)
      ELSE
         Y2H_G_C(0)  = SUM(Y2Y_C(:) * H_TMP(:))
         Y2H_G_C(J)   = Y2H_G_C(0) + Y2CP_C(J) 
         Y2CPBAR_C(J) = Y2H_G_C(J)
      ENDIF
      Y2MU_C(J)    = SUM(Y2Y_C(:) * MU_TMP(:))
      Y2K_C(J)     = SUM(Y2Y_C(:) * K_TMP(:))      
      Y2D_C(J)     = SUM(Y2Y_C(:) * D_TMP(:))            
   ELSE
      SS => SPECIES(0)
      CALL CALC_GAS_PROPS(J,0,D_TMP(1),MU_TMP(1),K_TMP(1),CP_TMP(1),H_TMP(1),SS%ISFUEL)
      SS%D(J) = D_TMP(1)
      Y2CP_C(J) = CP_TMP(1)
      IF (J>1) THEN
         Y2H_G_C(J)   = Y2H_G_C(J-1) + 0.5_EB*(Y2CP_C(J)+Y2CP_C(J-1)) 
         Y2CPBAR_C(J) = (Y2CPBAR_C(J-1)*REAL(J-1,EB)+Y2CP_C(J))/REAL(J,EB)
      ELSE
         Y2H_G_C(0 )  = H_TMP(1)      
         Y2H_G_C(J)   = H_TMP(1) + Y2CP_C(J) 
         Y2CPBAR_C(J) = Y2H_G_C(J)
      ENDIF
      Y2MU_C(J)    = MU_TMP(1)
      Y2K_C(J)     = K_TMP(1)      
   ENDIF
ENDDO

DEALLOCATE(D_TMP)
DEALLOCATE(MU_TMP)
DEALLOCATE(CP_TMP)
DEALLOCATE(H_TMP)
DEALLOCATE(K_TMP)

Y2CP_C(0) = Y2CP_C(1)
Y2CP(0,:) = Y2CP(1,:)
Y2CPBAR_C(0) = Y2CP_C(1)
Y2CPBAR(0,:) = Y2CPBAR(1,:)
Y2K_C(0) = Y2K_C(1)
Y2K(0,:)  = Y2K(1,:)
Y2MU_C(0)=Y2MU_C(1)
Y2MU(0,:) = Y2MU(1,:)

IF (N_MIX_SPECIES > 0) THEN
   Y2D_C(0) = Y2D_C(1)
   Y2D(0,:) = Y2D(1,:)
ENDIF

! Set ISFUEL property for SIMPLE_CHEMISTRY case

IF (SIMPLE_CHEMISTRY) SPECIES(I_FUEL)%ISFUEL = .TRUE.

END SUBROUTINE PROC_SPEC

 
SUBROUTINE PROC_PART
USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_ENTHALPY
USE PROPERTY_DATA, ONLY: JANAF_TABLE_LIQUID
CHARACTER(30) :: SPEC_ID
INTEGER :: N, J, ITMP
REAL(EB) :: H_L,H_V,YY_GET(1:N_GAS_SPECIES),H_G_S,H_G_S_REF,H_L_REF,TMP_REF,TMP_MELT,TMP_V,TMP_WGT
TYPE(PARTICLE_CLASS_TYPE), POINTER :: PC

IF (N_PART == 0) RETURN

PART_LOOP: DO N=1,N_PART

   PC => PARTICLE_CLASS(N)   

   ! Assign device or controller

   CALL SEARCH_CONTROLLER('PART',PC%CTRL_ID,PC%DEVC_ID,PC%DEVC_INDEX,PC%CTRL_INDEX,N)

   ! Only process particles or droplets that have mass and evaporate

   IF (.NOT.PC%EVAPORATE .OR. PC%MASSLESS .OR. PC%SURF_ID/='null') CYCLE PART_LOOP
   SPEC_ID = SPECIES(PC%SPEC_INDEX)%ID
   IF (SIMPLE_CHEMISTRY) THEN
      IF(PC%Y_INDEX==I_FUEL) SPEC_ID = REACTION(1)%FUEL
   ENDIF
   TMP_REF = -1._EB
   TMP_MELT = -1._EB
   TMP_V = -1._EB
   DO J = 1, 5000
      IF (PC%C_P(J) > 0._EB) THEN
         PC%H_L(J) = (REAL(J,EB)-PC%TMP_MELT)*PC%C_P(J)
      ELSE
         CALL JANAF_TABLE_LIQUID (J,PC%C_P(J),H_V,H_L,TMP_REF,TMP_MELT,TMP_V,SPEC_ID,PC%FUEL)
         IF (J==1) THEN         
            IF (PC%C_P(J) < 0._EB .AND. .NOT.PC%TREE) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: PARTicle class ',TRIM(PC%ID),' requires CP, H_V, TMP_MELT, TMP_V, and T_REF'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (PC%H_V_REFERENCE_TEMPERATURE < 0._EB) PC%H_V_REFERENCE_TEMPERATURE=TMP_REF
            IF (PC%TMP_V < 0._EB) PC%TMP_V = TMP_V
            IF (PC%TMP_MELT < 0._EB) PC%TMP_MELT = TMP_MELT
            PC%H_L(J) = H_L + PC%C_P(J)
         ELSE
            PC%H_L(J) = PC%H_L(J-1) + 0.5_EB*(PC%C_P(J)+PC%C_P(J-1))
         ENDIF
      ENDIF
      IF (J==1) THEN
         PC%C_P_BAR(J) =  PC%H_L(J)
      ELSE
         PC%C_P_BAR(J) = PC%H_L(J) / REAL(J,EB)
      ENDIF
   ENDDO

   ! Determine the properties of the droplet
   IF (PC%SPEC_INDEX>0) THEN
      IF(PC%H_V(1) > 0._EB) H_V = PC%H_V(1)
      ITMP = INT(PC%H_V_REFERENCE_TEMPERATURE)
      TMP_WGT  = PC%H_V_REFERENCE_TEMPERATURE - REAL(ITMP,EB)
      H_L_REF = PC%H_L(ITMP)+TMP_WGT*(PC%H_L(ITMP+1)-PC%H_L(ITMP))
      YY_GET = 0._EB
      YY_GET(PC%Y_INDEX) = 1._EB
      CALL GET_SPECIFIC_ENTHALPY(YY_GET,H_G_S_REF,PC%H_V_REFERENCE_TEMPERATURE)      
      DO J=1,5000
         CALL GET_SPECIFIC_ENTHALPY(YY_GET,H_G_S,REAL(J,EB))
         PC%H_V(J) = H_V + (H_G_S-H_G_S_REF) - (PC%H_L(J)-H_L_REF)
      ENDDO
   ENDIF
ENDDO PART_LOOP

END SUBROUTINE PROC_PART


 
SUBROUTINE READ_PART
USE DEVICE_VARIABLES, ONLY : PROPERTY_TYPE, PROPERTY, N_PROP
USE RADCONS, ONLY : NDG
USE PROPERTY_DATA, ONLY: JANAF_TABLE_LIQUID
USE MATH_FUNCTIONS, ONLY: NORM2
INTEGER :: SAMPLING_FACTOR,DROPLETS_PER_SECOND,N,NN,IPC,RGB(3),I_DUMMY,NUMBER_INITIAL_DROPLETS, &
           NUMBER_RADIATION_ANGLES
REAL(EB) :: SPECIFIC_HEAT,VAPORIZATION_TEMPERATURE,MELTING_TEMPERATURE,DIAMETER, &
            GAMMA_D,AGE,INITIAL_TEMPERATURE,HEAT_OF_COMBUSTION,HEAT_OF_VAPORIZATION,DENSITY,DT_INSERT, &
            VERTICAL_VELOCITY,HORIZONTAL_VELOCITY,MAXIMUM_DIAMETER,MINIMUM_DIAMETER,SIGMA_D, &
            SURFACE_TENSION,BREAKUP_CHILD_DIAMETER,BREAKUP_CHILD_GAMMA_D,BREAKUP_CHILD_SIGMA_D,&
            DENSE_VOLUME_FRACTION
REAL(EB) :: VEG_SV,VEG_MOISTURE,VEG_CHAR_FRACTION,VEG_DRAG_COEFFICIENT,VEG_DENSITY,VEG_BULK_DENSITY, &
            VEG_BURNING_RATE_MAX,VEG_DEHYDRATION_RATE_MAX,VEG_INITIAL_TEMPERATURE, &
            VEG_FUEL_MPV_MIN,VEG_MOIST_MPV_MIN,H_V_REFERENCE_TEMPERATURE,USER_DRAG_COEFFICIENT
REAL(EB) :: C_P,H_V,H_L,TMP_REF,TMP_MELT,TMP_V !Used only to inquire if species is fuel 
REAL(EB), DIMENSION(1000,3)::ORIENTATION
CHARACTER(30) :: PART_ID,SPEC_ID,DEVC_ID,CTRL_ID,QUANTITIES(1:10),SURF_ID,DRAG_LAW,PROP_ID
CHARACTER(25) :: COLOR
CHARACTER(25) :: VEG_DEGRADATION
LOGICAL :: MASSLESS,STATIC,TREE,MONODISPERSE,EVAPORATE,BREAKUP
LOGICAL :: VEG_REMOVE_CHARRED,VEG_STEM,VEG_CHAR_OXIDATION
TYPE(PARTICLE_CLASS_TYPE), POINTER :: PC
TYPE(PROPERTY_TYPE), POINTER :: PY
NAMELIST /PART/ FYI,DROPLETS_PER_SECOND,SAMPLING_FACTOR,ID,STATIC,MASSLESS,TREE, &
                DENSITY,VAPORIZATION_TEMPERATURE,SPECIFIC_HEAT,HEAT_OF_VAPORIZATION,NUMBER_INITIAL_DROPLETS, &
                MELTING_TEMPERATURE,DIAMETER,MAXIMUM_DIAMETER,MINIMUM_DIAMETER,GAMMA_D,HEAT_OF_COMBUSTION, &
                AGE,DEVC_ID,CTRL_ID,SPEC_ID,SURF_ID,INITIAL_TEMPERATURE,RGB,QUANTITIES,DT_INSERT,COLOR, &
                VERTICAL_VELOCITY,HORIZONTAL_VELOCITY,GAMMA,MONODISPERSE,EVAPORATE,SIGMA_D, &
                VEG_SV,VEG_MOISTURE,VEG_CHAR_FRACTION,VEG_DRAG_COEFFICIENT,VEG_DENSITY,VEG_BULK_DENSITY, &
                VEG_BURNING_RATE_MAX,VEG_DEHYDRATION_RATE_MAX,VEG_INITIAL_TEMPERATURE,VEG_FUEL_MPV_MIN, &
                VEG_MOIST_MPV_MIN,VEG_REMOVE_CHARRED,VEG_STEM,VEG_DEGRADATION,VEG_CHAR_OXIDATION, &
                H_V_REFERENCE_TEMPERATURE,DRAG_LAW,USER_DRAG_COEFFICIENT, &
                SURFACE_TENSION,BREAKUP,BREAKUP_CHILD_DIAMETER,BREAKUP_CHILD_GAMMA_D,BREAKUP_CHILD_SIGMA_D,PROP_ID,&
                DENSE_VOLUME_FRACTION,NUMBER_RADIATION_ANGLES,ORIENTATION

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

IF (N_PART>0) DROPLET_FILE = .TRUE.
 
! Allocate the derived type array to hold information about the particle classes

ALLOCATE(PARTICLE_CLASS(N_PART),STAT=IZERO)
CALL ChkMemErr('READ','PARTICLE_CLASS',IZERO) 

N_EVAP_INDICES = 0
IPC = 0
READ_PART_LOOP: DO N=1,N_PART
   
   IPC = IPC + 1
   PC=>PARTICLE_CLASS(IPC)

   BREAKUP                  = .FALSE.
   BREAKUP_CHILD_DIAMETER   = 3._EB/7._EB  ! ratio of child Sauter mean to parent size in Bag breakup regime
   BREAKUP_CHILD_GAMMA_D    = 2.4_EB
   BREAKUP_CHILD_SIGMA_D    = -99999.9_EB
   CTRL_ID                  = 'null'
   DENSE_VOLUME_FRACTION    = 1.E-5_EB     ! Limiting volume fraction for drag reduction
   DENSITY                  = 1000._EB     ! kg/m3
   DEVC_ID                  = 'null'
   DROPLETS_PER_SECOND      = -1
   DT_INSERT                = -1._EB
   EVAPORATE                = .TRUE.
   VAPORIZATION_TEMPERATURE = -300._EB     ! C
   INITIAL_TEMPERATURE      = TMPA - TMPM  ! C
   MELTING_TEMPERATURE      = -300.        ! C
   SPECIFIC_HEAT            = -1._EB     ! kJ/kg-K
   HEAT_OF_VAPORIZATION     = -1._EB     ! kJ/kg
   H_V_REFERENCE_TEMPERATURE = -300._EB
   HEAT_OF_COMBUSTION       = -1._EB       ! kJ/kg
   GAMMA                    = 1.4_EB       !Specific heat ratio
   DIAMETER                 = 500._EB      ! microns
   MAXIMUM_DIAMETER         = 1.E9_EB      ! microns, meant to be infinitely large and not used
   MINIMUM_DIAMETER         =  20._EB      ! microns, below which the droplet evaporates in one time step
   MONODISPERSE             = .FALSE.
   NUMBER_INITIAL_DROPLETS  = 0
   NUMBER_RADIATION_ANGLES  = 1
   GAMMA_D                  = 2.4_EB
   SIGMA_D                  = -99999.9_EB
   AGE                      = 1.E6_EB      ! s
   ID                       = 'null'
   PROP_ID                  = 'null'
   ORIENTATION(1,:)         = (/0._EB,0._EB,1._EB/)
   QUANTITIES               = 'null'
   RGB                      = -1
   SPEC_ID                  = 'null'
   SURF_ID                  = 'null'
   SURFACE_TENSION          = 72.8E-3_EB  ! N/m, applies for water
   COLOR                    = 'null'
   SAMPLING_FACTOR          = -1      
   STATIC                   = .FALSE.
   MASSLESS                 = .FALSE.
   TREE                     = .FALSE.
   VERTICAL_VELOCITY        = 0.5_EB
   HORIZONTAL_VELOCITY      = 0.2_EB
   VEG_SV                   = 4000. ! 1/m
   VEG_MOISTURE             = 10.0_EB
   VEG_CHAR_FRACTION        = 0.25_EB
   VEG_DRAG_COEFFICIENT     = 1.0_EB
   VEG_DENSITY              = 540._EB ! kg/m3 dry mass
   VEG_BULK_DENSITY         = 0.3_EB  ! kg/m3
   VEG_BURNING_RATE_MAX     = 0.4_EB  ! kg/m3/s
   VEG_DEHYDRATION_RATE_MAX = 0.4_EB  ! kg/m3/s
   VEG_INITIAL_TEMPERATURE  = TMPA - TMPM  ! C
   VEG_FUEL_MPV_MIN         = VEG_CHAR_FRACTION*VEG_BULK_DENSITY
   VEG_MOIST_MPV_MIN        = 0.01_EB*VEG_MOISTURE
   VEG_REMOVE_CHARRED       = .TRUE.
   VEG_DEGRADATION          = 'LINEAR'
   VEG_CHAR_OXIDATION       = .FALSE.
   VEG_STEM                 = .FALSE.
   DRAG_LAW                 = 'SPHERE'
   USER_DRAG_COEFFICIENT    = -1._EB
 
   ! Read the PART line from the input file or set up special PARTICLE_CLASS class for water droplets or tracers
 
   CALL CHECKREAD('PART',LU_INPUT,IOS) 
   IF (IOS==1) EXIT READ_PART_LOOP
   READ(LU_INPUT,PART)

   ! Error messages

   IF (NUMBER_INITIAL_DROPLETS>0) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: NUMBER_INITIAL_DROPLETS is now an INIT parameter. Consult the FDS Users Guide'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (DT_INSERT>0._EB) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: DT_INSERT is no longer specified on the PART line. Consult the FDS Users Guide'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (DROPLETS_PER_SECOND>0) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: DROPLETS_PER_SECOND is no longer specified on the PART line. Consult the FDS Users Guide'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   ! Miscellaneous consequences of input parameters

   IF (TREE)                                    STATIC                  = .TRUE.
   IF (TREE)                                    WFDS                    = .TRUE.
   IF (SURF_ID/='null')                         STATIC                  = .TRUE.
   IF (MASSLESS)                                DIAMETER                = 0._EB
   IF (TRIM(SPEC_ID)=='WATER VAPOR') THEN
      WATER_PART_INDEX = N
      GAMMA = 1.32_EB
      PC%WATER = .TRUE.
      IF (ANY(RGB<0) .AND. COLOR=='null') COLOR='BLUE'
   ENDIF
   IF (SIMPLE_CHEMISTRY) THEN
      IF (TRIM(SPEC_ID)==TRIM(REACTION(1)%FUEL)) THEN
         PC%FUEL = .TRUE.
         SPEC_ID= 'LUMPED SPECIES 1'
         IF (ANY(RGB<0) .AND. COLOR=='null') COLOR='YELLOW'
      ENDIF
   ELSE
     CALL JANAF_TABLE_LIQUID (1,C_P,H_V,H_L,TMP_REF,TMP_MELT,TMP_V,SPEC_ID,PC%FUEL) 
     IF (PC%FUEL .AND. ANY(RGB<0) .AND. COLOR=='null') COLOR = 'YELLOW'
   ENDIF

   IF (ANY(RGB<0) .AND. COLOR=='null') THEN
      COLOR = 'BLACK'
      IF (TREE)  COLOR = 'GREEN'
   ENDIF

   ! Set default colors for Smokeview

   IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)

   IF (SAMPLING_FACTOR<=0 .AND.      MASSLESS)  SAMPLING_FACTOR         = 1
   IF (SAMPLING_FACTOR<=0 .AND. .NOT.MASSLESS)  SAMPLING_FACTOR         = 10
    PC%SPEC_ID = SPEC_ID
   IF (PC%SPEC_ID/='null') THEN
      IF (MASSLESS) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Cannot have MASSLESS=.TRUE. with evaporating droplets'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      DO NN=1,N_SPECIES
         IF (TRIM(SPECIES(NN)%ID)==TRIM(PC%SPEC_ID)) THEN
            IF (SPECIES(NN)%Y_INDEX <=0) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: PART SPEC_ID ',TRIM(PC%SPEC_ID),' not a tracked gas species'
               CALL SHUTDOWN(MESSAGE)
            ENDIF         
            PC%SPEC_INDEX = NN           
            PC%Y_INDEX = SPECIES(NN)%Y_INDEX 
         ENDIF
      ENDDO
      IF(PC%SPEC_INDEX < 0) WRITE(MESSAGE,'(A,A,A)') 'ERROR: PART SPEC_ID ',TRIM(PC%SPEC_ID),' not found'
   ENDIF
   IF ((MELTING_TEMPERATURE > -300._EB .OR. SPECIFIC_HEAT > 0._EB .OR. HEAT_OF_VAPORIZATION > 0._EB .OR. &
        H_V_REFERENCE_TEMPERATURE > -300._EB .OR. HEAT_OF_COMBUSTION > 0._EB) .AND. TRIM(SPEC_ID)=='null') THEN
      WRITE(MESSAGE,'(A,I3,A)') 'ERROR: PART line ',N,' requires a SPEC_ID.'
      CALL SHUTDOWN(MESSAGE)
   ENDIF  

   PC%QUANTITIES = QUANTITIES
   PC%GAMMA_VAP  = GAMMA
 
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
   IF (BREAKUP) THEN
      ALLOCATE(PC%CHILD_CDF(0:NDC),STAT=IZERO)
      CALL ChkMemErr('READ','CHILD_CDF',IZERO)
      ALLOCATE(PC%CHILD_R_CDF(0:NDC),STAT=IZERO)
      CALL ChkMemErr('READ','CHILD_R_CDF',IZERO)
   ENDIF      
   
   ! Assign property data to PARTICLE_CLASS class
 
   PC%ID                 = ID
   PC%BREAKUP            = BREAKUP
   PC%BREAKUP_CHILD_DIAMETER = BREAKUP_CHILD_DIAMETER
   PC%BREAKUP_CHILD_GAMMA    = BREAKUP_CHILD_GAMMA_D
   IF ( BREAKUP_CHILD_SIGMA_D > 0._EB ) THEN
      PC%BREAKUP_CHILD_SIGMA = BREAKUP_CHILD_SIGMA_D
   ELSE
      PC%BREAKUP_CHILD_SIGMA = 1.15_EB/BREAKUP_CHILD_GAMMA_D
   END IF
   PC%CTRL_ID            = CTRL_ID
   PC%DENSE_VOLUME_FRACTION = DENSE_VOLUME_FRACTION
   PC%DEVC_ID            = DEVC_ID
   PC%TMP_INITIAL        = INITIAL_TEMPERATURE + TMPM
   PC%SAMPLING           = SAMPLING_FACTOR
   PC%RGB                = RGB
   PC%DIAMETER           = DIAMETER*1.E-6_EB
   PC%MAXIMUM_DIAMETER   = MAXIMUM_DIAMETER*1.E-6_EB
   PC%MINIMUM_DIAMETER   = MINIMUM_DIAMETER*1.E-6_EB
   PC%KILL_RADIUS        = MINIMUM_DIAMETER*1.E-6_EB*0.25_EB
   PC%MONODISPERSE       = MONODISPERSE
   PC%PROP_ID            = PROP_ID
   PC%GAMMA              = GAMMA_D
   IF ( SIGMA_D > 0._EB ) THEN
      PC%SIGMA           = SIGMA_D
   ELSE
      PC%SIGMA           = 1.15_EB/GAMMA_D
   END IF
   PC%TMP_INITIAL        = INITIAL_TEMPERATURE + TMPM
   IF ((HEAT_OF_VAPORIZATION >  0._EB .AND. SPECIFIC_HEAT <= 0._EB) .OR. &
       (HEAT_OF_VAPORIZATION <= 0._EB .AND. SPECIFIC_HEAT >  0._EB)) THEN
      WRITE(MESSAGE,'(A,I4,A)') 'ERROR: PART ' ,N,' If one of SPECIFIC_HEAT or HEAT_OF_VAPORIZATION defined, both must be'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (HEAT_OF_VAPORIZATION > 0._EB) PC%H_V = HEAT_OF_VAPORIZATION*1000._EB
   IF (SPECIFIC_HEAT > 0._EB) THEN
      PC%C_P = SPECIFIC_HEAT*1000._EB   
      PC%H_V = HEAT_OF_VAPORIZATION*1000._EB
      IF (MELTING_TEMPERATURE < -TMPM) THEN
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: PART ' ,N,' MELTING_TEMPERATURE not set'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      PC%TMP_MELT = MELTING_TEMPERATURE + TMPM
      IF (H_V_REFERENCE_TEMPERATURE < -TMPM) H_V_REFERENCE_TEMPERATURE = MELTING_TEMPERATURE
      PC%H_V_REFERENCE_TEMPERATURE = H_V_REFERENCE_TEMPERATURE + 273.15_EB
      IF (VAPORIZATION_TEMPERATURE< -TMPM) THEN
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: PART ' ,N,' VAPORIZATION_TEMPERATURE not set'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      PC%TMP_V = VAPORIZATION_TEMPERATURE + TMPM
   ENDIF
   PC%HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION*1000._EB
   PC%DENSITY            = DENSITY
   PC%FTPR               = FOTH*PI*DENSITY
   PC%MASSLESS           = MASSLESS
   PC%LIFETIME           = AGE
   PC%EVAPORATE          = EVAPORATE
   PC%TREE               = TREE
   PC%STATIC             = STATIC
   PC%SPEC_ID            = SPEC_ID
   PC%SURF_ID            = SURF_ID
   PC%SURF_INDEX         = -1      ! SURFs have not yet been read in
   PC%SURFACE_TENSION    = SURFACE_TENSION
   IF (SURF_ID/='null') VIRTUAL_PARTICLES = .TRUE.
   PC%ADJUST_EVAPORATION  = 1._EB   ! If H_O_C>0. this parameter will have to be reset later
   PC%VERTICAL_VELOCITY   = VERTICAL_VELOCITY
   PC%HORIZONTAL_VELOCITY = HORIZONTAL_VELOCITY
   PC%USER_DRAG_COEFFICIENT = USER_DRAG_COEFFICIENT
   IF (USER_DRAG_COEFFICIENT>=0._EB) DRAG_LAW = 'USER'
   SELECT CASE(DRAG_LAW)
      CASE('SPHERE')
         PC%DRAG_LAW = SPHERE_DRAG
      CASE('CYLINDER')
         PC%DRAG_LAW = CYLINDER_DRAG
      CASE('USER')
         PC%DRAG_LAW = USER_DRAG
      CASE DEFAULT
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: unrecognized drag law on PART line'
         CALL SHUTDOWN(MESSAGE)
   END SELECT
 
   ! Vegetation propeties

   PC%VEG_SV                   = VEG_SV !1/m
   PC%VEG_MOISTURE             = VEG_MOISTURE
   PC%VEG_CHAR_FRACTION        = VEG_CHAR_FRACTION
   PC%VEG_DRAG_COEFFICIENT     = VEG_DRAG_COEFFICIENT
   PC%VEG_DENSITY              = VEG_DENSITY !kg/m3
   PC%VEG_BULK_DENSITY         = VEG_BULK_DENSITY !kg/m3
   PC%VEG_BURNING_RATE_MAX     = VEG_BURNING_RATE_MAX !kg/m3.s
   PC%VEG_DEHYDRATION_RATE_MAX = VEG_DEHYDRATION_RATE_MAX !kg/m3.s
   PC%VEG_INITIAL_TEMPERATURE  = VEG_INITIAL_TEMPERATURE +TMPM ! K
   PC%VEG_FUEL_MPV_MIN         = VEG_CHAR_FRACTION*VEG_BULK_DENSITY
   PC%VEG_MOIST_MPV_MIN        = 0.01_EB*VEG_MOISTURE
   PC%VEG_REMOVE_CHARRED       = VEG_REMOVE_CHARRED
   PC%VEG_STEM                 = VEG_STEM
   PC%VEG_DEGRADATION          = VEG_DEGRADATION
   PC%VEG_CHAR_OXIDATION       = VEG_CHAR_OXIDATION

   ! Set evaporation index
   IF (.NOT.PC%MASSLESS) THEN
      N_EVAP_INDICES = N_EVAP_INDICES + 1
      PC%EVAP_INDEX   = N_EVAP_INDICES
   ENDIF

   ! Split particles

   PC%N_SPLIT                  = NUMBER_RADIATION_ANGLES
   IF (PC%N_SPLIT>1) THEN   
      ALLOCATE(PC%ORIENTATION(1:PC%N_SPLIT,1:3),STAT=IZERO)
      CALL ChkMemErr('READ','ORIENTATION',IZERO)
      DO NN = 1,PC%N_SPLIT
         PC%ORIENTATION(NN,1:3)   = ORIENTATION(NN,1:3)/ NORM2(ORIENTATION(NN,:))
      ENDDO
   ENDIF
   
ENDDO READ_PART_LOOP

! Assign PART_INDEX to Device PROPERTY array

DO N=1,N_PROP
   PY => PROPERTY(N)
   IF (PY%PART_ID/='null') THEN
      DO IPC=1,N_PART
         PC=>PARTICLE_CLASS(IPC)
         IF (PC%ID==PY%PART_ID) PY%PART_INDEX = IPC
         IF (PC%ID==PY%PART_ID .AND. PC%MASSLESS) THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: PART_ID for PROP ' ,N,' cannot refer to MASSLESS particles'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO
      IF (PY%PART_INDEX<0) THEN
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: PART_ID for PROP ' ,N,' not found'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      DROPLET_FILE=.TRUE.
   ENDIF
ENDDO

! Allocate radiation arrays

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
                                 PC%QUANTITIES_INDEX(PC%N_QUANTITIES),I_DUMMY, &
                                 I_DUMMY,I_DUMMY,I_DUMMY,I_DUMMY,'PART', &
                                 PC%QUANTITIES(NN),'null',SPEC_ID,PART_ID,'null','null')   
      ENDDO QUANTITIES_LOOP
   ENDIF
ENDDO 

IF (SIMPLE_CHEMISTRY) THEN

   ! For mixture fraction model, adjust the evaporation rate of fuel droplets to account for difference in HoC.

   DO IPC=1,N_PART
      PC=>PARTICLE_CLASS(IPC)
      IF (PC%HEAT_OF_COMBUSTION > 0._EB) PC%ADJUST_EVAPORATION = PC%HEAT_OF_COMBUSTION/REACTION(1)%HEAT_OF_COMBUSTION
   ENDDO

ENDIF


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
LOGICAL  :: IGNITOR_ELEMENTS,OUTPUT_TREE
CHARACTER(30) :: PART_ID, LABEL
CHARACTER(30) :: FUEL_GEOM,RAMP_IGNELEM
TYPE(PARTICLE_CLASS_TYPE), POINTER :: PC
NAMELIST /TREE/ XYZ,XB,CROWN_WIDTH,CROWN_WIDTH_BOTTOM,CROWN_WIDTH_TOP,CROWN_BASE_HEIGHT,TREE_HEIGHT,PART_ID, &
                IGNITOR_ELEMENTS,RING_THICKNESS,TON_IGNITOR_ELEMENTS,  &
                TOFF_IGNITOR_ELEMENTS,T_RAMPOFF_IGNITOR_ELEMENTS,      &
                T_RAMPON_IGNITOR_ELEMENTS,FUEL_GEOM,RAMP_IGNELEM,LABEL,OUTPUT_TREE

X_OVERALL_MIN = 1.E12_EB ; X_OVERALL_MAX = -1.E12_EB
Y_OVERALL_MIN = 1.E12_EB ; Y_OVERALL_MAX = -1.E12_EB
Z_OVERALL_MIN = 1.E12_EB ; Z_OVERALL_MAX = -1.E12_EB

! Read the TREE lines to determine how many vegetation volumes there will be
! and how many trees volumes will have data ouputs

N_TREES = 0
N_TREES_OUT = 0
REWIND(LU_INPUT)
COUNT_VEG_LOOP: DO
   CALL CHECKREAD('TREE',LU_INPUT,IOS) 
   IF (IOS==1) EXIT COUNT_VEG_LOOP
   OUTPUT_TREE=.FALSE.
   READ(LU_INPUT,NML=TREE,END=11,ERR=12,IOSTAT=IOS)
   N_TREES = N_TREES + 1
   IF (OUTPUT_TREE) N_TREES_OUT = N_TREES_OUT + 1
   12 IF (IOS>0) THEN 
         WRITE(MESSAGE,'(A,I6)') 'ERROR: Problem with TREE line ',N_TREES+1
         CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_VEG_LOOP
11 REWIND(LU_INPUT)
!
! Sequentially read the TREE namelist to get shape and size parameters for each
! vegetation volume.
!
ALLOCATE(TREE_MESH(NMESHES),STAT=IZERO)
CALL ChkMemErr('READ','TREE_MESH',IZERO); TREE_MESH = .TRUE.

IF (N_TREES==0) THEN
 TREE_MESH = .FALSE.
 RETURN
ENDIF
!
TREE_CASE = .TRUE.
!
IF (N_TREES_OUT /= 0) THEN
 ALLOCATE(TREE_MESH_OUT(NMESHES),STAT=IZERO)
 CALL ChkMemErr('READ','TREE_MESH_OUT',IZERO) ; TREE_MESH_OUT=.FALSE.
 ALLOCATE(TREE_OUTPUT_DATA(N_TREES_OUT,4,NMESHES),STAT=IZERO)
 CALL ChkMemErr('READ','TREE_OUTPUT_DATA',IZERO)
 ALLOCATE(TREE_OUTPUT_DATA_TOTAL(N_TREES_OUT,4),STAT=IZERO)
 CALL ChkMemErr('READ','TREE_OUTPUT_DATA_TOTAL',IZERO)
ENDIF

ALLOCATE(N_TREE_OUT(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','N_TREE_OUT',IZERO) ; N_TREE_OUT = 0
ALLOCATE(VEG_FUEL_GEOM(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','VEG_FUEL_GEOM',IZERO)
ALLOCATE(VEG_LABELS(N_TREES),STAT=IZERO)
CALL ChkMemErr('READ','VEG_LABELS',IZERO)
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
N_TREES_OUT = 0
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
   OUTPUT_TREE = .FALSE.
   READ(LU_INPUT,TREE,END=25,IOSTAT=IOS)
   IF (PART_ID  =='null') CALL SHUTDOWN('ERROR: Specify PART_ID in TREE namelist')
   IF (FUEL_GEOM=='null') CALL SHUTDOWN('ERROR: Specify FUEL_GEOM in TREE namelist')
!
! Identify trees that user requested diagnostics for
   IF (OUTPUT_TREE) THEN
    N_TREES_OUT   = N_TREES_OUT + 1
    N_TREE_OUT(N) = N_TREES_OUT
   ENDIF

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
   IF (LABEL == 'null') THEN
     VEG_LABELS(N) = 'no_veg_data_ouput'
   ELSE
     VEG_LABELS(N) = LABEL
   ENDIF
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
! Find meshes that require user specified ouput of vegetation quantities
   FIND_TREE_MESH_OUT: DO NM = 1, NMESHES
     IF (.NOT. OUTPUT_TREE) CYCLE FIND_TREE_MESH_OUT
     IF (TREE_MESH_OUT(NM)) CYCLE FIND_TREE_MESH_OUT
     IF (X_TREE_MIN <= MESHES(NM)%XF .AND. X_TREE_MIN >= MESHES(NM)%XS) &
         TREE_MESH_OUT(NM) = .TRUE.
     IF (X_TREE_MAX <= MESHES(NM)%XF .AND. X_TREE_MAX >= MESHES(NM)%XS) &
         TREE_MESH_OUT(NM) = .TRUE.
     IF (MESHES(NM)%XS <= X_TREE_MAX .AND. MESHES(NM)%XS >= X_TREE_MIN) &
         TREE_MESH_OUT(NM) = .TRUE.
     IF (MESHES(NM)%XF <= X_TREE_MAX .AND. MESHES(NM)%XF >= X_TREE_MIN) &
         TREE_MESH_OUT(NM) = .TRUE.
   ENDDO FIND_TREE_MESH_OUT
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
            BEAD_DIAMETER,BEAD_EMISSIVITY,BEAD_SPECIFIC_HEAT,BEAD_DENSITY,BEAD_H_FIXED, &
            CONDUIT_DIAMETER,CONDUIT_THICKNESS,CABLE_MASS_PER_LENGTH,CABLE_DIAMETER, &
            CABLE_JACKET_THICKNESS,CABLE_FAILURE_TEMPERATURE, &
            C_FACTOR,CHARACTERISTIC_VELOCITY,DT_INSERT,ORIFICE_DIAMETER, &
            DROPLET_VELOCITY,FLOW_RATE,FLOW_TAU,GAUGE_TEMPERATURE,INITIAL_TEMPERATURE,K_FACTOR,LENGTH,SPRAY_ANGLE(2), &
            OFFSET,OPERATING_PRESSURE,RTI,PDPA_START,PDPA_END,PDPA_RADIUS
INTEGER :: I,N,NN,PDPA_M,PDPA_N,DROPLETS_PER_SECOND
LOGICAL :: PDPA_INTEGRATE,PDPA_NORMALIZE
EQUIVALENCE(LENGTH,ALPHA_C)
CHARACTER(30) :: SMOKEVIEW_ID(SMOKEVIEW_OBJECTS_DIMENSION),QUANTITY='null',PART_ID='null',FLOW_RAMP='null', &
                 SPRAY_PATTERN_TABLE='null',SPEC_ID='null',&
                 PRESSURE_RAMP='null',SMOKEVIEW_PARAMETERS(SMOKEVIEW_OBJECTS_DIMENSION)
TYPE (PROPERTY_TYPE), POINTER :: PY

NAMELIST /PROP/ ACTIVATION_OBSCURATION,ACTIVATION_TEMPERATURE,ALPHA_C,ALPHA_E,BETA_C,BETA_E, &
                BEAD_DIAMETER,BEAD_EMISSIVITY,BEAD_H_FIXED,BEAD_SPECIFIC_HEAT,BEAD_DENSITY, &
                CONDUIT_DIAMETER,CONDUIT_THICKNESS,CABLE_MASS_PER_LENGTH,CABLE_DIAMETER, &
                CABLE_JACKET_THICKNESS,CABLE_FAILURE_TEMPERATURE, &
                C_FACTOR,CHARACTERISTIC_VELOCITY,DROPLETS_PER_SECOND,DT_INSERT,ORIFICE_DIAMETER, &
                DROPLET_VELOCITY,FLOW_RATE,FLOW_RAMP,FLOW_TAU,ID,GAUGE_TEMPERATURE,INITIAL_TEMPERATURE,K_FACTOR,LENGTH,OFFSET, &
                OPERATING_PRESSURE,PART_ID,QUANTITY,RTI,SPRAY_ANGLE,SMOKEVIEW_ID,SPEC_ID,SPRAY_PATTERN_TABLE,PRESSURE_RAMP, &
                PDPA_START,PDPA_END,PDPA_RADIUS,PDPA_M,PDPA_N,PDPA_INTEGRATE,PDPA_NORMALIZE,SMOKEVIEW_PARAMETERS

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
   PY%BEAD_DENSITY             = BEAD_DENSITY
   PY%BEAD_DIAMETER            = BEAD_DIAMETER
   PY%BEAD_EMISSIVITY          = BEAD_EMISSIVITY
   PY%BEAD_H_FIXED             = BEAD_H_FIXED    
   PY%BEAD_SPECIFIC_HEAT       = BEAD_SPECIFIC_HEAT*1000._EB
   PY%CABLE_DIAMETER           = CABLE_DIAMETER
   PY%CABLE_JACKET_THICKNESS   = CABLE_JACKET_THICKNESS
   PY%CABLE_MASS_PER_LENGTH    = CABLE_MASS_PER_LENGTH
   PY%CABLE_FAILURE_TEMPERATURE= CABLE_FAILURE_TEMPERATURE
   PY%CONDUIT_DIAMETER         = CONDUIT_DIAMETER
   PY%CONDUIT_THICKNESS        = CONDUIT_THICKNESS
   IF (CONDUIT_DIAMETER>0._EB) CONDUIT = .TRUE.
   PY%C_FACTOR                 = C_FACTOR
   PY%CHARACTERISTIC_VELOCITY  = CHARACTERISTIC_VELOCITY
   PY%DT_INSERT                = DT_INSERT
   PY%GAUGE_TEMPERATURE        = GAUGE_TEMPERATURE + TMPM
   PY%ID                       = ID
   PY%INITIAL_TEMPERATURE      = INITIAL_TEMPERATURE + TMPM
   PY%N_INSERT                 = DROPLETS_PER_SECOND*DT_INSERT
   PY%OFFSET                   = OFFSET
   PY%OPERATING_PRESSURE       = OPERATING_PRESSURE
   ALLOCATE(PY%PARTICLE_INSERT_CLOCK(NMESHES),STAT=IZERO)
   CALL ChkMemErr('READ','PARTICLE_INSERT_CLOCK',IZERO)
   PY%PARTICLE_INSERT_CLOCK    = T_BEGIN
   PY%PART_ID                  = PART_ID
   PY%QUANTITY                 = QUANTITY
   IF (PY%PART_ID/='null' .AND. PY%QUANTITY == 'null' ) PY%QUANTITY = 'NOZZLE'
   PY%RTI                      = RTI
   IF (SMOKEVIEW_ID(1)/='null') THEN
      PY%SMOKEVIEW_ID          = SMOKEVIEW_ID
      PY%N_SMOKEVIEW_IDS = 0
      DO NN=1,SMOKEVIEW_OBJECTS_DIMENSION
         IF (SMOKEVIEW_ID(NN)/='null') PY%N_SMOKEVIEW_IDS = PY%N_SMOKEVIEW_IDS + 1
      ENDDO
   ELSE
      PY%N_SMOKEVIEW_IDS = 1
      SELECT CASE(PY%QUANTITY)
         CASE DEFAULT
            PY%SMOKEVIEW_ID(1) = 'sensor'
         CASE('SPRINKLER LINK TEMPERATURE')
            PY%SMOKEVIEW_ID(1) = 'sprinkler_pendent'
         CASE('NOZZLE')
            PY%SMOKEVIEW_ID(1) = 'nozzle'
         CASE('LINK TEMPERATURE')
            PY%SMOKEVIEW_ID(1) = 'heat_detector'
         CASE('CABLE TEMPERATURE')
            PY%SMOKEVIEW_ID(1) = 'sensor'
         CASE('spot obscuration','CHAMBER OBSCURATION')
            PY%SMOKEVIEW_ID(1) = 'smoke_detector'
         CASE('THERMOCOUPLE')
            PY%SMOKEVIEW_ID(1) = 'thermocouple'
      END SELECT
   ENDIF
   PY%SMOKEVIEW_PARAMETERS = SMOKEVIEW_PARAMETERS
   PY%N_SMOKEVIEW_PARAMETERS = 0
   DO I=1,SMOKEVIEW_OBJECTS_DIMENSION
      IF (PY%SMOKEVIEW_PARAMETERS(I)/='null') PY%N_SMOKEVIEW_PARAMETERS = PY%N_SMOKEVIEW_PARAMETERS + 1
   ENDDO
   PY%SPEC_ID              = SPEC_ID
   IF (PART_ID/='null' .AND. SPRAY_PATTERN_TABLE /= 'null') THEN
      CALL GET_TABLE_INDEX(SPRAY_PATTERN_TABLE,SPRAY_PATTERN,PY%SPRAY_PATTERN_INDEX)
      PY%TABLE_ID = SPRAY_PATTERN_TABLE
   ELSE
      PY%SPRAY_PATTERN_INDEX = 0
   ENDIF 
   PY%SPRAY_ANGLE(1) = SPRAY_ANGLE(1)*PI/180._EB
   PY%SPRAY_ANGLE(2) = SPRAY_ANGLE(2)*PI/180._EB

   ! THIEF model

   IF (QUANTITY=='CABLE TEMPERATURE') THEN
      N_CABL = N_CABL + 1
      CABLE(N_CABL)%PROP_INDEX = N
   ENDIF

   ! PDPA model

   PY%PDPA_START       = PDPA_START
   PY%PDPA_END         = PDPA_END
   PY%PDPA_RADIUS      = PDPA_RADIUS
   PY%PDPA_M           = PDPA_M
   PY%PDPA_N           = PDPA_N
   PY%PDPA_INTEGRATE   = PDPA_INTEGRATE
   PY%PDPA_NORMALIZE   = PDPA_NORMALIZE
   IF (PY%QUANTITY == 'NUMBER CONCENTRATION') THEN 
      PY%PDPA_M        = 0
      PY%PDPA_N        = 0
   ENDIF
   IF ((PY%QUANTITY == 'MASS CONCENTRATION').OR. &
       (PY%QUANTITY == 'ENTHALPY')          .OR. &
       (PY%QUANTITY == 'DROPLET FLUX X')    .OR. &
       (PY%QUANTITY == 'DROPLET FLUX Y')    .OR. &
       (PY%QUANTITY == 'DROPLET FLUX Z')) THEN
      PY%PDPA_M        = 3
      PY%PDPA_N        = 0
   ENDIF

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

   IF (PART_ID /='null' .AND. ABS(PDPA_RADIUS) <= ZERO_P) THEN
      IF ((FLOW_RATE>0._EB .AND. K_FACTOR<=0._EB .AND. OPERATING_PRESSURE<=0._EB) .OR. &
          (FLOW_RATE<0._EB .AND. K_FACTOR>=0._EB .AND. OPERATING_PRESSURE<=0._EB) .OR. &
          (FLOW_RATE<0._EB .AND. K_FACTOR<=0._EB .AND. OPERATING_PRESSURE>0._EB)) THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: Problem with PROP ',TRIM(PY%ID),', too few flow parameters'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      IF (K_FACTOR < 0._EB .AND. OPERATING_PRESSURE > 0._EB)  K_FACTOR           = FLOW_RATE/SQRT(OPERATING_PRESSURE)
      IF (FLOW_RATE < 0._EB .AND. OPERATING_PRESSURE > 0._EB) FLOW_RATE          = K_FACTOR*SQRT(OPERATING_PRESSURE)
      IF (OPERATING_PRESSURE < 0._EB .AND. K_FACTOR > 0._EB)  OPERATING_PRESSURE = (FLOW_RATE/K_FACTOR)**2
      PY%K_FACTOR           = K_FACTOR
      PY%FLOW_RATE          = FLOW_RATE
      PY%OPERATING_PRESSURE = OPERATING_PRESSURE

      IF (DROPLET_VELOCITY<=ZERO_P .AND. ORIFICE_DIAMETER<=ZERO_P .AND. &
          PRESSURE_RAMP=='null' .AND. SPRAY_PATTERN_TABLE=='null') THEN
         WRITE(MESSAGE,'(A,A,A)') 'WARNING: PROP ',TRIM(PY%ID),' droplet velocity is not defined.'
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
   IF (ABS(FLOW_TAU) > ZERO_P) THEN
      PY%FLOW_TAU = FLOW_TAU 
      IF (FLOW_TAU > 0._EB) PY%FLOW_RAMP_INDEX = TANH_RAMP 
      IF (FLOW_TAU < 0._EB) PY%FLOW_RAMP_INDEX = TSQR_RAMP
   ENDIF 

   ! Check for SPEC_ID
   IF (PY%SPEC_ID/='null') THEN
      DO I=1,N_SPECIES
         IF (TRIM(SPECIES(I)%ID)==TRIM(PY%SPEC_ID)) THEN
            PY%SPEC_INDEX = I
            PY%Y_INDEX = SPECIES(I)%Y_INDEX
            EXIT
         ENDIF
      ENDDO
      IF (PY%SPEC_INDEX < 0) THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: PROP SPEC_ID ',TRIM(PY%SPEC_ID),' not found'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
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
BEAD_DENSITY             = 8908._EB    ! kg/m3 (Nickel)
BEAD_DIAMETER            = 0.001       ! m
BEAD_EMISSIVITY          = 0.85_EB
BEAD_H_FIXED             = -1._EB      ! W/m2/K
BEAD_SPECIFIC_HEAT       = 0.44_EB     ! kJ/kg/K (Nickel)
CABLE_FAILURE_TEMPERATURE= 400._EB     ! C
CABLE_DIAMETER           = 0.02_EB     ! m
CABLE_JACKET_THICKNESS   = 0.002_EB    ! m
CABLE_MASS_PER_LENGTH    = 0.5_EB      ! kg/m
CONDUIT_DIAMETER         = 0.0         ! m
CONDUIT_THICKNESS        = 0.0         ! m
C_FACTOR                 = 0.0_EB
CHARACTERISTIC_VELOCITY  = 1.0_EB      ! m/s
DROPLET_VELOCITY         = 0._EB       ! m/s
DROPLETS_PER_SECOND      = 5000
DT_INSERT                = 0.01        ! s
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
PDPA_RADIUS              = 0.0_EB
PDPA_M                   = 0
PDPA_N                   = 0
PDPA_INTEGRATE           = .TRUE.
PDPA_NORMALIZE           = .TRUE.
PRESSURE_RAMP            = 'null'
QUANTITY                 = 'null'
RTI                      = 100._EB     ! (ms)**0.5
SMOKEVIEW_ID             = 'null' 
SMOKEVIEW_PARAMETERS     = 'null' 
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
            WRITE(MESSAGE,'(A,A,A,I5)') 'ERROR: Spray Pattern Table, ',TRIM(PY%TABLE_ID),', massflux <= 0 for line ',NN
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
CHARACTER(30) :: CONDUCTIVITY_RAMP,SPECIFIC_HEAT_RAMP,SPEC_ID(MAX_SPECIES)
REAL(EB) :: EMISSIVITY,CONDUCTIVITY,SPECIFIC_HEAT,DENSITY,ABSORPTION_COEFFICIENT,BOILING_TEMPERATURE, &
            INITIAL_VAPOR_FLUX,PEAK_REACTION_RATE
REAL(EB), DIMENSION(1:MAX_REACTIONS) :: A,E,HEATING_RATE,PYROLYSIS_RANGE,HEAT_OF_REACTION,NU_FUEL,NU_WATER,NU_RESIDUE, &
                                        N_S,N_T,REFERENCE_RATE,REFERENCE_TEMPERATURE,THRESHOLD_TEMPERATURE,HEAT_OF_COMBUSTION, &
                                        THRESHOLD_SIGN
REAL(EB), DIMENSION(1:MAX_REACTIONS,1:MAX_SPECIES) :: NU_GAS
LOGICAL, DIMENSION(1:MAX_REACTIONS) :: PCR
CHARACTER(30), DIMENSION(1:MAX_REACTIONS) :: RESIDUE
INTEGER :: N,NN,NNN,IOS,NR,N_REACTIONS,N_CONDUIT
NAMELIST /MATL/ ID,FYI,SPECIFIC_HEAT,CONDUCTIVITY,CONDUCTIVITY_RAMP,SPECIFIC_HEAT_RAMP, &
                REFERENCE_TEMPERATURE, REFERENCE_RATE, THRESHOLD_TEMPERATURE, THRESHOLD_SIGN, &
                EMISSIVITY,HEATING_RATE,PYROLYSIS_RANGE,HEAT_OF_REACTION,DENSITY,RESIDUE, &
                HEAT_OF_COMBUSTION,A,E,NU_FUEL,NU_WATER,NU_RESIDUE,NU_GAS,N_S,N_T,N_REACTIONS,&
                ABSORPTION_COEFFICIENT,BOILING_TEMPERATURE,INITIAL_VAPOR_FLUX,PCR,SPEC_ID

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

ALLOCATE(MATL_SPEC_ID(1:N_MATL,1:MAX_SPECIES))
MATL_SPEC_ID(1:N_MATL,1:MAX_SPECIES)='null'

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

      IF ( ( ANY(THRESHOLD_TEMPERATURE>-TMPM) .OR. ANY(REFERENCE_TEMPERATURE>-TMPM) .OR. ANY(A>=0._EB) .OR. ANY(E>=0._EB) .OR. &
             ANY(ABS(HEAT_OF_REACTION)>ZERO_P) ) .AND. N_REACTIONS==0) THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: Problem with MATL number ',TRIM(ID),'. A reaction parameter is used, but N_REACTIONS=0'  
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      DO NN=1,N_REACTIONS
         IF (REFERENCE_TEMPERATURE(NN)<-TMPM  .AND. (E(NN)< 0._EB .OR. A(NN)<0._EB)) THEN
            WRITE(MESSAGE,'(A,A,A,I2,A)') 'ERROR: Problem with MATL ',TRIM(ID),', REAC ',NN,'. Set REFERENCE_TEMPERATURE or E, A'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         IF (ABS(NU_FUEL(NN))<=ZERO_P .AND. ABS(NU_WATER(NN))<=ZERO_P .AND. &
             ABS(NU_RESIDUE(NN))<=ZERO_P .AND. ABS(SUM(NU_GAS(NN,:)))<=ZERO_P) THEN
            WRITE(MESSAGE,'(A,A,A,I2,A)') 'WARNING: MATL ',TRIM(ID),', REAC ',NN,'. No product yields (NUs) set'  
            IF (MYID==0) WRITE(LU_ERR,'(A)') TRIM(MESSAGE)
         ENDIF
      ENDDO
   ELSE NOT_BOILING ! Is liquid
      IF (ABS(HEAT_OF_REACTION(1))<=ZERO_P) THEN
         WRITE(MESSAGE,'(A,A)') 'ERROR: HEAT_OF_REACTION should be greater than zero for liquid MATL ',TRIM(ID)  
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDIF NOT_BOILING

   ! Error checking for thermal properties

   IF (ABS(DENSITY) <=ZERO_P ) THEN
      WRITE(MESSAGE,'(A,A,A)') 'ERROR: Problem with MATL ',TRIM(ID),': DENSITY=0' 
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (ABS(CONDUCTIVITY) <=ZERO_P .AND. CONDUCTIVITY_RAMP == 'null' ) THEN
      WRITE(MESSAGE,'(A,A,A)') 'ERROR: Problem with MATL ',TRIM(ID),': CONDUCTIVITY = 0' 
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (ABS(SPECIFIC_HEAT)<=ZERO_P .AND. SPECIFIC_HEAT_RAMP == 'null' ) THEN
      WRITE(MESSAGE,'(A,A,A)') 'ERROR: Problem with MATL ',TRIM(ID),': SPECIFIC_HEAT = 0' 
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (SPECIFIC_HEAT > 10._EB) WRITE(LU_ERR,'(A,A)') 'WARNING: SPECIFIC_HEAT units are kJ/kg/K check MATL ',TRIM(ID)
 
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
   ML%NU_GAS               = NU_GAS
   MATL_SPEC_ID(N,:)       = SPEC_ID
   ML%RAMP_C_S             = SPECIFIC_HEAT_RAMP
   ML%RAMP_K_S             = CONDUCTIVITY_RAMP
   ML%RHO_S                = DENSITY
   ML%RESIDUE_MATL_NAME(:) = RESIDUE(:)
   ML%HEATING_RATE(:)      = HEATING_RATE(:)/60._EB
   ML%PYROLYSIS_RANGE(:)   = PYROLYSIS_RANGE(:)
   ML%PCR(:)               = PCR(:)
   ML%TMP_BOIL             = BOILING_TEMPERATURE + TMPM
   ML%TMP_THR(:)           = THRESHOLD_TEMPERATURE(:) + TMPM
   ML%TMP_REF(:)           = REFERENCE_TEMPERATURE(:) + TMPM
   ML%THR_SIGN(:)          = THRESHOLD_SIGN
   ML%RATE_REF(:)          = REFERENCE_RATE(:) 
 
   ! Additional logic

   IF (BOILING_TEMPERATURE<5000._EB) THEN
      ML%PYROLYSIS_MODEL = PYROLYSIS_LIQUID
      ML%N_REACTIONS = 1
      IF (ABS(ML%NU_FUEL(1))<=ZERO_P) ML%NU_FUEL(1) = 1._EB
   ELSE
      ML%PYROLYSIS_MODEL = PYROLYSIS_SOLID
      IF (N_REACTIONS==0) ML%PYROLYSIS_MODEL = PYROLYSIS_NONE
   ENDIF

   DO NN=1,ML%N_REACTIONS
      IF (NU_FUEL(NN) > 0._EB) SIMPLE_CHEMISTRY = .TRUE.
   ENDDO

   IF (ML%RAMP_K_S/='null') THEN
      CALL GET_RAMP_INDEX(ML%RAMP_K_S,'TEMPERATURE',NR)
      ML%K_S = -NR
   ENDIF

   IF (ML%RAMP_C_S/='null') THEN
      CALL GET_RAMP_INDEX(ML%RAMP_C_S,'TEMPERATURE',NR)
      ML%C_S = -NR
   ENDIF

   ! Determine A and E if REFERENCE_TEMPERATURE is specified

   DO NN=1,ML%N_REACTIONS
      IF (ML%TMP_REF(NN) > 0._EB) THEN
         IF (ML%RATE_REF(NN) > 0._EB) THEN
            PEAK_REACTION_RATE = ML%RATE_REF(NN)
         ELSE
            PEAK_REACTION_RATE = 2._EB*ML%HEATING_RATE(NN)*(1._EB-ML%NU_RESIDUE(NN))/ML%PYROLYSIS_RANGE(NN)
         ENDIF
         ML%E(NN) = EXP(1._EB)*PEAK_REACTION_RATE*R0*ML%TMP_REF(NN)**2/ML%HEATING_RATE(NN)
         ML%A(NN) = EXP(1._EB)*PEAK_REACTION_RATE*EXP(ML%E(NN)/(R0*ML%TMP_REF(NN)))
      ENDIF
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
         WRITE(MESSAGE,'(5A)') 'ERROR: Residue ', TRIM(ML%RESIDUE_MATL_NAME(NN)),' of ',TRIM(MATL_NAME(N)),' is not defined.'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDDO
ENDDO

! Check for duplicate names

IF (N_MATL>1) THEN
   DO N=1,N_MATL-1
      DO NN=N+1,N_MATL
         IF(MATL_NAME(N)==MATL_NAME(NN)) THEN
            WRITE(MESSAGE,'(A,A)') 'ERROR: Duplicate material name: ',TRIM(MATL_NAME(N))
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO
   ENDDO
ENDIF

CONTAINS
 
SUBROUTINE SET_MATL_DEFAULTS
 
A                      = -1._EB      ! 1/s
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
THRESHOLD_SIGN         = 1.0
N_REACTIONS            = 0
N_S                    = 1._EB
N_T                    = 0._EB
NU_FUEL                = 0._EB
NU_GAS                 = 0._EB
NU_RESIDUE             = 0._EB
NU_WATER               = 0._EB
PCR                    = .FALSE.
REFERENCE_RATE         = -1._EB
REFERENCE_TEMPERATURE  = -1000._EB
RESIDUE                = 'null'
SPECIFIC_HEAT          = 0.0_EB      ! kJ/kg/K
SPECIFIC_HEAT_RAMP     = 'null'
SPEC_ID                = 'null'
HEATING_RATE           = 5._EB       ! K/min
PYROLYSIS_RANGE        = 80._EB      ! K or C
 
END SUBROUTINE SET_MATL_DEFAULTS
 
END SUBROUTINE READ_MATL



SUBROUTINE PROC_MATL

USE MATH_FUNCTIONS
INTEGER :: N,J,NS,NS2,NR,RN_SPEC_INDEX(MAX_SPECIES)
REAL(EB) :: NU_GAS(1:MAX_REACTIONS,1:MAX_SPECIES)

PROC_MATL_LOOP: DO N=1,N_MATL
   
   RN_SPEC_INDEX = -1
   ML => MATERIAL(N)

   !Get species for NU_GAS and HEAT_OF COMBUSTION
   FIND_SPECIES:IF (ANY(ML%NU_GAS > 0._EB)) THEN
      NU_GAS = ML%NU_GAS
      IF(TRIM(MATL_SPEC_ID(N,1))=='null') THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: MATL ',TRIM(MATL_NAME(N)),' requires SPEC_ID'
         CALL SHUTDOWN(MESSAGE)                  
      ENDIF
      
      !Build lookup array for MATL SPEC_ID input
      DO NS=1,N_SPECIES
         IF (TRIM(MATL_SPEC_ID(N,NS))=='null') EXIT
         DO NS2=1,N_GAS_SPECIES
            IF(TRIM(MATL_SPEC_ID(N,NS))==TRIM(SPECIES(Y2SPEC(NS2))%ID)) THEN
               RN_SPEC_INDEX(NS)=NS2
               EXIT
            ENDIF
         ENDDO
         IF (RN_SPEC_INDEX(NS)==-1) THEN
            WRITE(MESSAGE,'(A,A,A,A)') 'ERROR: SPEC:',TRIM(MATL_SPEC_ID(N,NS)),&
                                       ' is not a tracked species on MATL ',TRIM(MATL_NAME(N))
            CALL SHUTDOWN(MESSAGE)                  
         ENDIF      
      ENDDO
      
      ML%NU_GAS=0._EB
      
      !Build NU array for MATL
      REACTION_SEARCH: DO NR=1,ML%N_REACTIONS
         IF(ALL(ABS(NU_GAS(NR,:))<=ZERO_P)) CYCLE REACTION_SEARCH
         SPEC_LOOP_1: DO NS=1,N_SPECIES
            IF(ABS(NU_GAS(NR,NS)) <=ZERO_P) CYCLE SPEC_LOOP_1
            IF(RN_SPEC_INDEX(NS)==-1) THEN
               WRITE(MESSAGE,'(A,I2,A)') 'ERROR: Reacation ',NR,' has more species than are defined.'
               CALL SHUTDOWN(MESSAGE)                  
            ENDIF
            ML%NU_GAS(NR,RN_SPEC_INDEX(NS))=NU_GAS(NR,NS)
          ENDDO SPEC_LOOP_1
      ENDDO REACTION_SEARCH
   ENDIF FIND_SPECIES

   ! Adjust burn rate if heat of combustion is different from the gas phase reaction value

   DO J = 1,MAX(1,ML%N_REACTIONS)

      IF (SIMPLE_CHEMISTRY) THEN
         IF (N_REACTIONS > 0) THEN
            IF (CO_PRODUCTION) THEN
               RN => REACTION(2)         
            ELSE
               RN => REACTION(1)
            ENDIF
            IF (ML%HEAT_OF_COMBUSTION(J)>0._EB .AND. RN%HEAT_OF_COMBUSTION>0._EB)  &
             ML%ADJUST_BURN_RATE(J,I_FUEL) = ML%HEAT_OF_COMBUSTION(J)/(RN%Y_F_INLET*RN%HEAT_OF_COMBUSTION)
         ENDIF
      ELSEIF (N_REACTIONS>0) THEN
         RN => REACTION(1)
         DO NS = 1,N_GAS_SPECIES
            IF (ML%HEAT_OF_COMBUSTION(J)>0._EB .AND. RN%HEAT_OF_COMBUSTION>0._EB)  &
                ML%ADJUST_BURN_RATE(J,NS) = ML%HEAT_OF_COMBUSTION(J)/RN%HEAT_OF_COMBUSTION
         ENDDO
      ENDIF
   ENDDO

   ! Collect species yields to NU_GAS array

   DO J = 1,ML%N_REACTIONS
      IF (SIMPLE_CHEMISTRY) THEN
         IF (I_FUEL /= 0 .AND. ML%NU_FUEL(J)  > 0._EB) ML%NU_GAS(J,I_FUEL)  = ML%NU_FUEL(J)
         IF (I_WATER /= 0.AND. ML%NU_WATER(J) > 0._EB) ML%NU_GAS(J,I_WATER) = ML%NU_WATER(J)
      ELSE
         IF (ABS(ML%NU_FUEL(J))<=ZERO_P) THEN
            DO NS = 1,N_GAS_SPECIES
               IF (SPECIES(NS)%ISFUEL) ML%NU_FUEL(J) = ML%NU_FUEL(J) + ML%NU_GAS(J,NS)
            ENDDO
         ENDIF
      ENDIF
   ENDDO
   
   IF (ML%RAMP_C_S/='null') THEN
         IF (RAMPS(-INT(ML%C_S))%DEPENDENT_DATA(1) > 10._EB) &
            WRITE(LU_ERR,'(A,A)') 'WARNING: SPECIFIC_HEAT units are kJ/kg/K check MATL ',TRIM(ID)
   ENDIF   

ENDDO PROC_MATL_LOOP

IF (ALLOCATED(MATL_SPEC_ID)) DEALLOCATE(MATL_SPEC_ID)

END SUBROUTINE PROC_MATL


SUBROUTINE READ_SURF
 
USE MATH_FUNCTIONS, ONLY : GET_RAMP_INDEX
USE DEVICE_VARIABLES, ONLY : CABLE,N_CABL,PROPERTY_TYPE,PROPERTY
CHARACTER(30) :: PART_ID,RAMP_MF(0:MAX_SPECIES),RAMP_Q,RAMP_V,RAMP_T,MATL_ID(MAX_LAYERS,MAX_MATERIALS),&
                 PROFILE,BACKING,GEOMETRY,NAME_LIST(MAX_MATERIALS),EXTERNAL_FLUX_RAMP,RAMP_EF,RAMP_PART,SPEC_ID(MAX_SPECIES)
EQUIVALENCE(EXTERNAL_FLUX_RAMP,RAMP_EF)
LOGICAL :: ADIABATIC,BURN_AWAY,SHRINK,POROUS,FREE_SLIP,NO_SLIP
CHARACTER(60) :: TEXTURE_MAP
CHARACTER(25) :: COLOR
TYPE (PROPERTY_TYPE), POINTER :: PY
REAL(EB) :: TAU_Q,TAU_V,TAU_T,TAU_MF(0:MAX_SPECIES),HRRPUA,MLRPUA,TEXTURE_WIDTH,TEXTURE_HEIGHT,VEL_T(2), &
            TAU_EXTERNAL_FLUX,TAU_EF,E_COEFFICIENT,VOLUME_FLUX,TMP_FRONT,TMP_INNER(MAX_LAYERS),THICKNESS(MAX_LAYERS),VEL, &
            MASS_FLUX(0:MAX_SPECIES),MASS_FRACTION(0:MAX_SPECIES), Z0,PLE,CONVECTIVE_HEAT_FLUX,PARTICLE_MASS_FLUX, &
            TRANSPARENCY,EXTERNAL_FLUX,TMP_BACK,MASS_FLUX_TOTAL,STRETCH_FACTOR,CONVECTION_LENGTH_SCALE, &
            MATL_MASS_FRACTION(MAX_LAYERS,MAX_MATERIALS),CELL_SIZE_FACTOR,MAX_PRESSURE,&
            IGNITION_TEMPERATURE,HEAT_OF_VAPORIZATION,REGRID_FACTOR,NET_HEAT_FLUX,LAYER_DIVIDE,SURFACE_DENSITY, &
            ROUGHNESS,RADIUS,LENGTH,WIDTH,DT_INSERT,H_FIXED,TAU_PART,EMISSIVITY,EMISSIVITY_BACK,EMISSIVITY_DEFAULT, &
            SPREAD_RATE,XYZ(3)
EQUIVALENCE(TAU_EXTERNAL_FLUX,TAU_EF)
INTEGER :: NPPC,N,IOS,NL,NN,NNN,N_LIST,N_LIST2,INDEX_LIST(MAX_MATERIALS_TOTAL),LEAK_PATH(2),DUCT_PATH(2),RGB(3),NR,IL
INTEGER ::  VEGETATION_LAYERS
REAL(EB) :: VEGETATION_CDRAG,VEGETATION_CHAR_FRACTION,VEGETATION_ELEMENT_DENSITY,VEGETATION_HEIGHT, &
            VEGETATION_INITIAL_TEMP,VEGETATION_LOAD,VEGETATION_LSET_IGNITE_TIME,VEGETATION_MOISTURE,VEGETATION_SVRATIO, &
            FIRELINE_MLR_MAX,VEGETATION_GROUND_TEMP,VEG_LSET_ROS_HEAD,VEG_LSET_ROS_FLANK,VEG_LSET_ROS_BACK, &
            VEG_LSET_WIND_EXP
LOGICAL :: VEGETATION,VEGETATION_NO_BURN,VEGETATION_LINEAR_DEGRAD,VEGETATION_ARRHENIUS_DEGRAD,VEG_LEVEL_SET_SPREAD
NAMELIST /SURF/ TMP_FRONT,TMP_INNER,THICKNESS,MASS_FRACTION,VEL,VEL_T,NPPC, &
                E_COEFFICIENT,CONVECTIVE_HEAT_FLUX,TAU_Q,TAU_V,TAU_T,RAMP_Q,RAMP_T,TAU_MF,TAU_EXTERNAL_FLUX,TAU_EF, &
                RAMP_MF,PART_ID,RAMP_V,VOLUME_FLUX, PROFILE,PLE,Z0,ID,MASS_FLUX,PARTICLE_MASS_FLUX, &
                FYI,MATL_ID,BACKING,TMP_BACK,HRRPUA,MLRPUA,SHRINK,CELL_SIZE_FACTOR, &
                TEXTURE_MAP,TEXTURE_WIDTH,TEXTURE_HEIGHT,RGB,TRANSPARENCY, BURN_AWAY,LEAK_PATH,DUCT_PATH,ADIABATIC, &
                MASS_FLUX_TOTAL,GEOMETRY,STRETCH_FACTOR,MATL_MASS_FRACTION,EMISSIVITY,COLOR,POROUS,MAX_PRESSURE,&
                IGNITION_TEMPERATURE,HEAT_OF_VAPORIZATION,REGRID_FACTOR,NET_HEAT_FLUX,LAYER_DIVIDE,SURFACE_DENSITY, &
                ROUGHNESS,FREE_SLIP,NO_SLIP,EXTERNAL_FLUX,EXTERNAL_FLUX_RAMP,RAMP_EF,RADIUS,LENGTH,WIDTH,DT_INSERT,H_FIXED, &
                RAMP_PART,TAU_PART,EMISSIVITY_BACK,SPREAD_RATE,XYZ, &
                VEGETATION,VEGETATION_CDRAG,VEGETATION_CHAR_FRACTION,VEGETATION_ELEMENT_DENSITY,VEGETATION_HEIGHT, &
                VEGETATION_INITIAL_TEMP,VEGETATION_LSET_IGNITE_TIME,VEGETATION_LOAD,VEGETATION_MOISTURE, &
                VEGETATION_NO_BURN,VEGETATION_SVRATIO,VEGETATION_LAYERS,CONVECTION_LENGTH_SCALE, &
                FIRELINE_MLR_MAX,VEGETATION_GROUND_TEMP,VEGETATION_LINEAR_DEGRAD,VEGETATION_ARRHENIUS_DEGRAD, &
                VEG_LSET_ROS_HEAD,VEG_LSET_ROS_FLANK,VEG_LSET_ROS_BACK,VEG_LEVEL_SET_SPREAD,VEG_LSET_WIND_EXP,SPEC_ID

             
! Count the SURF lines in the input file

REWIND(LU_INPUT)
N_SURF = 0
COUNT_SURF_LOOP: DO
   HRRPUA = 0._EB
   MLRPUA = 0._EB
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
 
N_SURF_RESERVED = 5
IF (.NOT. NO_EVACUATION) N_SURF_RESERVED = 6
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
HVAC_SURF_INDEX                    = N_SURF + 5
IF (.NOT. NO_EVACUATION) EVACUATION_SURF_INDEX = N_SURF + 6

N_SURF                             = N_SURF + N_SURF_RESERVED

SURFACE(INERT_SURF_INDEX)%ID       = 'INERT'
SURFACE(OPEN_SURF_INDEX)%ID        = 'OPEN'
SURFACE(MIRROR_SURF_INDEX)%ID      = 'MIRROR'
SURFACE(INTERPOLATED_SURF_INDEX)%ID= 'INTERPOLATED'
SURFACE(PERIODIC_SURF_INDEX)%ID    = 'PERIODIC'
SURFACE(HVAC_SURF_INDEX)%ID        = 'HVAC'
SURFACE(N_SURF-N_SURF_RESERVED+1:N_SURF)%USER_DEFINED = .FALSE.
SURFACE(0)%USER_DEFINED               = .FALSE.
IF (.NOT. NO_EVACUATION) SURFACE(EVACUATION_SURF_INDEX)%ID = 'EVACUATION_OUTFLOW'
 
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
      CONVECTION_LENGTH_SCALE = PI*SUM(THICKNESS)
   ENDIF

   ! Vegetation parameters
   IF(VEGETATION) WFDS = .TRUE.
   
   ! Level set vegetation fire spread specific
   SF%VEG_LSET_SPREAD    = VEG_LEVEL_SET_SPREAD
   SF%VEG_LSET_ROS_HEAD  = VEG_LSET_ROS_HEAD !head fire rate of spread m/s
   SF%VEG_LSET_ROS_FLANK = VEG_LSET_ROS_FLANK !flank fire
   SF%VEG_LSET_ROS_BACK  = VEG_LSET_ROS_BACK !flank fire
   SF%VEG_LSET_WIND_EXP  = VEG_LSET_WIND_EXP !exponent on wind cosine in ROS formula

   ! Boundary Vegetation specific

   SF%VEGETATION = VEGETATION !T or F
   SF%VEG_NO_BURN  = VEGETATION_NO_BURN
!  IF(SF%VEGETATION) ADIABATIC = .TRUE.
   SF%VEG_CHARFRAC = VEGETATION_CHAR_FRACTION
   SF%VEG_MOISTURE = VEGETATION_MOISTURE
   SF%VEG_HEIGHT   = VEGETATION_HEIGHT
   SF%VEG_INITIAL_TEMP = VEGETATION_INITIAL_TEMP
   SF%VEG_GROUND_TEMP  = VEGETATION_GROUND_TEMP
   IF (VEGETATION_GROUND_TEMP /= -99._EB) SF%VEG_GROUND_ZERO_RAD = .FALSE.
   SF%VEG_LOAD     = VEGETATION_LOAD
   SF%FIRELINE_MLR_MAX = FIRELINE_MLR_MAX
!  SF%VEG_DEHYDRATION_RATE_MAX = SRF_VEG_DEHYDRATION_RATE_MAX
   SF%VEG_PACKING  = VEGETATION_LOAD/VEGETATION_HEIGHT/VEGETATION_ELEMENT_DENSITY
   SF%VEG_SVRATIO  = VEGETATION_SVRATIO
   SF%VEG_KAPPA    = 0.25_EB*VEGETATION_SVRATIO*SF%VEG_PACKING
   SF%NVEG_L       = INT(1. + VEGETATION_HEIGHT*3._EB*SF%VEG_KAPPA)
   IF(VEGETATION_LAYERS > 0) SF%NVEG_L = VEGETATION_LAYERS
   SF%VEG_DRAG_INI = VEGETATION_CDRAG*0.375_EB*SF%VEG_PACKING*SF%VEG_SVRATIO
   SF%VEG_LSET_IGNITE_T = VEGETATION_LSET_IGNITE_TIME
   SF%VEG_LINEAR_DEGRAD    = VEGETATION_LINEAR_DEGRAD
   SF%VEG_ARRHENIUS_DEGRAD = VEGETATION_ARRHENIUS_DEGRAD
   IF(VEGETATION_ARRHENIUS_DEGRAD) SF%VEG_LINEAR_DEGRAD = .FALSE.

   ALLOCATE(SF%VEG_FUEL_FLUX_L(SF%NVEG_L),STAT=IZERO)
   CALL ChkMemErr('READ_SURF','VEG_FUEL_FLUX_L',IZERO)
   ALLOCATE(SF%VEG_MOIST_FLUX_L(SF%NVEG_L),STAT=IZERO)
   CALL ChkMemErr('READ_SURF','VEG_MOIST_FLUX_L',IZERO)
   ALLOCATE(SF%VEG_DIVQNET_L(SF%NVEG_L),STAT=IZERO)
   CALL ChkMemErr('READ_SURF','VEG_DIVQNET',IZERO)

   ALLOCATE(SF%VEG_FINCM_RADFCT_L(0:SF%NVEG_L),STAT=IZERO) !add index for mult veg
   CALL ChkMemErr('READ_SURF','VEG_FINCM_RADFCT_L',IZERO)
   ALLOCATE(SF%VEG_FINCP_RADFCT_L(0:SF%NVEG_L),STAT=IZERO)
   CALL ChkMemErr('READ','VEG_FINCP_RADFCT_L',IZERO) 

   ALLOCATE(SF%VEG_SEMISSP_RADFCT_L(0:SF%NVEG_L,0:SF%NVEG_L),STAT=IZERO) !add index for mult veg
   CALL ChkMemErr('READ','VEG_SEMISSP_RADFCT_L',IZERO) 
   ALLOCATE(SF%VEG_SEMISSM_RADFCT_L(0:SF%NVEG_L,0:SF%NVEG_L),STAT=IZERO)
   CALL ChkMemErr('READ','VEG_SEMISSM_RADFCT_L',IZERO) 
   

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
   SF%CONV_LENGTH          = CONVECTION_LENGTH_SCALE
   SF%NET_HEAT_FLUX        = 1000._EB*NET_HEAT_FLUX
   SF%DUCT_PATH            = DUCT_PATH
   SF%DT_INSERT            = DT_INSERT
   SF%E_COEFFICIENT        = E_COEFFICIENT
   SF%EMISSIVITY           = EMISSIVITY   
   SF%EMISSIVITY_BACK      = EMISSIVITY_BACK
   SF%FIRE_SPREAD_RATE     = SPREAD_RATE
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
      CASE DEFAULT
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: SURF ',TRIM(SF%ID),' GEOMETRY not recognized'
         CALL SHUTDOWN(MESSAGE)
   END SELECT
   SF%H_V                  = 1000._EB*HEAT_OF_VAPORIZATION
   SF%HRRPUA               = 1000._EB*HRRPUA
   SF%MLRPUA               = MLRPUA
   SF%LAYER_DIVIDE         = LAYER_DIVIDE
   SF%LEAK_PATH            = LEAK_PATH
   SF%LENGTH               = LENGTH
   SF%MASS_FLUX            = 0._EB
   SF%MASS_FLUX(0)         = MASS_FLUX(0)
   SF%MASS_FRACTION        = -1._EB
   SF%MASS_FRACTION(0)     = MASS_FRACTION(0)
   SF%MAX_PRESSURE         = MAX_PRESSURE
   ALLOCATE(SF%PARTICLE_INSERT_CLOCK(NMESHES),STAT=IZERO)
   CALL ChkMemErr('READ','PARTICLE_INSERT_CLOCK',IZERO)
   SF%PARTICLE_INSERT_CLOCK = T_BEGIN
   SF%RADIUS               = RADIUS
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
   SF%RAMP_EF              = EXTERNAL_FLUX_RAMP
   SF%RAMP_MF              = 'null'
   SF%RAMP_MF(0)           = RAMP_MF(0)
   SF%RAMP_Q               = RAMP_Q 
   SF%RAMP_V               = RAMP_V  
   SF%RAMP_T               = RAMP_T  
   SF%RAMP_PART            = RAMP_PART
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
   SF%TAU(0)               = TAU_MF(0)
   SF%TAU(TIME_HEAT)       = TAU_Q
   SF%TAU(TIME_VELO)       = TAU_V
   SF%TAU(TIME_TEMP)       = TAU_T
   SF%TAU(TIME_EFLUX)      = TAU_EXTERNAL_FLUX
   SF%TAU(TIME_PART)       = TAU_PART
   SF%TEXTURE_MAP          = TEXTURE_MAP
   SF%TEXTURE_WIDTH        = TEXTURE_WIDTH
   SF%TEXTURE_HEIGHT       = TEXTURE_HEIGHT
   SF%TMP_IGN              = IGNITION_TEMPERATURE + TMPM
   SF%VEL                  = VEL
   SF%VEL_T                = VEL_T
   SF%VOLUME_FLUX          = VOLUME_FLUX
   SF%WIDTH                = WIDTH   
   SF%Z0                   = Z0
   SF%MASS_FLUX_TOTAL      = MASS_FLUX_TOTAL
   SF%H_FIXED              = H_FIXED
   SF%XYZ                  = XYZ

   IF (ANY(MASS_FLUX(1:N_GAS_SPECIES)>0._EB) .OR. ANY(MASS_FRACTION(1:N_GAS_SPECIES)>0._EB)) THEN
      IF(SPEC_ID(1)=='null') THEN
         WRITE (MESSAGE,'(A,A,A)') 'ERROR: Problem with SURF:',TRIM(SF%ID),&
                                   '. Must define SPEC_ID when using MASS_FLUX or MASS_FRACTION'
         CALL SHUTDOWN(MESSAGE)
      ELSE 
         DO NN=1,N_GAS_SPECIES
            IF (TRIM(SPEC_ID(NN))=='null') EXIT
            DO NNN=0,N_SPECIES
               IF (TRIM(SPECIES(NNN)%ID)==TRIM(SPEC_ID(NN))) THEN
                  IF (SPECIES(NNN)%Y_INDEX<=0) THEN
                     WRITE(MESSAGE,'(A,A,A,A)') 'ERROR: SPEC:',TRIM(SPEC_ID(NN)),'is not a tracked species on SURF:',TRIM(SF%ID)
                     CALL SHUTDOWN(MESSAGE)
                  ENDIF
                  SF%MASS_FLUX(SPECIES(NNN)%Y_INDEX)    = MASS_FLUX(NN)
                  SF%MASS_FRACTION(SPECIES(NNN)%Y_INDEX)= MASS_FRACTION(NN)
                  SF%TAU(SPECIES(NNN)%Y_INDEX)          = TAU_MF(NN)
                  SF%RAMP_MF(SPECIES(NNN)%Y_INDEX)      = RAMP_MF(NN)
                  EXIT
               ENDIF
               IF (NNN==N_SPECIES) THEN
                  WRITE(MESSAGE,'(A,A,A,A)') 'ERROR: Problem with SURF:',TRIM(SF%ID),' SPEC ',&
                     TRIM(SPEC_ID(NN)),' not found'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
   ENDIF
   
   ! Set various logical parameters

   IF (ABS(SF%VEL_T(1))>ZERO_P .OR. ABS(SF%VEL_T(2))>ZERO_P) SF%SPECIFIED_TANGENTIAL_VELOCITY = .TRUE.

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
      IF (NL==1) EMISSIVITY = 0._EB
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
               IF (NL==1) EMISSIVITY = EMISSIVITY + &
                  MATERIAL(NNN)%EMISSIVITY*SF%LAYER_MATL_FRAC(NL,NN)/MATERIAL(NNN)%RHO_S ! volume based
            ENDIF
         ENDDO
         IF (INDEX_LIST(N_LIST)<0) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: MATL_ID ',TRIM(NAME_LIST(N_LIST)),', on SURF: ',TRIM(SF%ID),', does not exist'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO COUNT_LAYER_MATL
      IF (SF%LAYER_DENSITY(NL) > 0._EB) SF%LAYER_DENSITY(NL) = 1./SF%LAYER_DENSITY(NL)
      IF (NL==1) EMISSIVITY = EMISSIVITY*SF%LAYER_DENSITY(NL)
      SF%THICKNESS = SF%THICKNESS + SF%LAYER_THICKNESS(NL)
   ENDDO COUNT_LAYERS

   ! Set front side emissivity
   IF (SF%EMISSIVITY < 0._EB) THEN
      IF (SF%N_LAYERS > 0) THEN
         SF%EMISSIVITY = EMISSIVITY
      ELSE
         SF%EMISSIVITY = EMISSIVITY_DEFAULT
      ENDIF
   ENDIF
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
            IF (N_LIST>MAX_MATERIALS_TOTAL) CALL SHUTDOWN('ERROR: Too many materials in the surface.')
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
      IF (SF%TMP_INNER(1)>0._EB)           SF%TMP_FRONT = SF%TMP_INNER(1)
      IF (SF%TMP_INNER(SF%N_LAYERS)>0._EB) SF%TMP_BACK  = SF%TMP_INNER(SF%N_LAYERS)
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
         IF (ABS(ML%NU_RESIDUE(J))<=ZERO_P) SF%SHRINK = .TRUE.
      ENDDO
      IF (ML%PYROLYSIS_MODEL==PYROLYSIS_LIQUID) SF%SHRINK = .TRUE.
   ENDDO
   IF (.NOT. SHRINK .AND. SF%SHRINK) THEN
      SF%SHRINK = .FALSE.
      WRITE(MESSAGE,'(A,A,A)') 'WARNING: SURF ',TRIM(SF%ID),' has SHRINK set .FALSE. while reactions suggest .TRUE.'
      IF (MYID==0) WRITE(LU_ERR,'(A)') TRIM(MESSAGE)
   ENDIF

   ! Thermal boundary conditions
   IF (SF%ADIABATIC .AND. (SF%NET_HEAT_FLUX < 1.E12_EB .OR. ABS(SF%CONVECTIVE_HEAT_FLUX)>ZERO_P)) THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: SURF '//TRIM(SF%ID)//&
                                 ' cannot have both ADIABATIC and NET_HEAT_FLUX or CONVECTIVE_HEAT_FLUX'
         CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (SF%NET_HEAT_FLUX < 1.E12_EB .AND. ABS(SF%CONVECTIVE_HEAT_FLUX)>ZERO_P) THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: SURF '//TRIM(SF%ID)// ' cannot have both NET_HEAT_FLUX or CONVECTIVE_HEAT_FLUX'
         CALL SHUTDOWN(MESSAGE)
   ENDIF
   
   SF%THERMAL_BC_INDEX = SPECIFIED_TEMPERATURE
   IF (SF%ADIABATIC) THEN
                                       SF%THERMAL_BC_INDEX = NET_FLUX_BC
                                       SF%NET_HEAT_FLUX = 0._EB
   ENDIF
   IF (SF%NET_HEAT_FLUX < 1.E12_EB)    SF%THERMAL_BC_INDEX = NET_FLUX_BC
   IF (ABS(SF%CONVECTIVE_HEAT_FLUX)>ZERO_P) SF%THERMAL_BC_INDEX = CONVECTIVE_FLUX_BC
   IF (SF%THERMALLY_THICK)             SF%THERMAL_BC_INDEX = THERMALLY_THICK
   IF (SF%PROFILE==ATMOSPHERIC)        SF%THERMAL_BC_INDEX = INFLOW_OUTFLOW
   IF (SF%VEGETATION)                  SF%THERMAL_BC_INDEX = VEG_BNDRY_FUEL

   ! Set convection length scale according to THICKNESS or RADIUS

   IF (SF%CONV_LENGTH < 0._EB) THEN
      IF (SF%THERMALLY_THICK) THEN
         RADIUS = SF%THICKNESS
      ELSE
         IF ((SF%GEOMETRY/=SURF_CARTESIAN) .AND. (RADIUS < 0._EB)) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: SURF ',TRIM(SF%ID),' needs a RADIUS'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDIF

      SELECT CASE(SF%GEOMETRY)
         CASE(SURF_CARTESIAN)
            SF%CONV_LENGTH = 0._EB
            IF ((LENGTH<0._EB).AND.(WIDTH<0._EB)) SF%CONV_LENGTH = 1.0_EB
            IF (LENGTH > 0._EB) SF%CONV_LENGTH = LENGTH
            IF (WIDTH > 0._EB) SF%CONV_LENGTH = SQRT(SF%CONV_LENGTH**2 + WIDTH**2)
         CASE(SURF_CYLINDRICAL)
            SF%CONV_LENGTH = PI*RADIUS
         CASE(SURF_SPHERICAL)
            SF%CONV_LENGTH = 2._EB*RADIUS
      END SELECT
   ENDIF

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

   IF (SF%RAMP_EF/='null') THEN
      CALL GET_RAMP_INDEX(SF%RAMP_EF,'TIME',NR)
      SF%RAMP_INDEX(TIME_EFLUX) = NR
   ELSE
      IF (SF%TAU(TIME_EFLUX) > 0._EB) SF%RAMP_INDEX(TIME_EFLUX) = TANH_RAMP
      IF (SF%TAU(TIME_EFLUX) < 0._EB) SF%RAMP_INDEX(TIME_EFLUX) = TSQR_RAMP
   ENDIF 

   IF (SF%RAMP_PART/='null') THEN
      CALL GET_RAMP_INDEX(SF%RAMP_PART,'TIME',NR)
      SF%RAMP_INDEX(TIME_PART) = NR
   ELSE
      IF (SF%TAU(TIME_PART) > 0._EB) SF%RAMP_INDEX(TIME_PART) = TANH_RAMP
      IF (SF%TAU(TIME_PART) < 0._EB) SF%RAMP_INDEX(TIME_PART) = TSQR_RAMP
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
CONVECTION_LENGTH_SCALE = -1._EB
NET_HEAT_FLUX           = 1.E12_EB
DT_INSERT               = 0.01_EB
DUCT_PATH               = 0 
E_COEFFICIENT           = 0._EB
EMISSIVITY              = -1._EB
EMISSIVITY_DEFAULT      = 0.9_EB
EMISSIVITY_BACK         = -1._EB
EXTERNAL_FLUX           = 0._EB
EXTERNAL_FLUX_RAMP      = 'null'
FREE_SLIP               = .FALSE.
NO_SLIP                 = .FALSE.
FYI                     = 'null'
GEOMETRY                = 'CARTESIAN'
HEAT_OF_VAPORIZATION    = 0._EB
H_FIXED                 = -1._EB
HRRPUA                  = 0._EB
ID                      = 'null'
IGNITION_TEMPERATURE    = 5000._EB
LAYER_DIVIDE            = -1._EB
LEAK_PATH               = -1 
LENGTH                  = -1._EB
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
RADIUS                  = -1._EB
RAMP_MF                 = 'null'
RAMP_Q                  = 'null'
RAMP_V                  = 'null'
RAMP_T                  = 'null'
RAMP_PART               = 'null'
REGRID_FACTOR           = 0.9
RGB                     = -1
SHRINK                  = .TRUE.
IF (LES) ROUGHNESS      = 0._EB !4.5E-5_EB  ! meters, commercial steel
IF (DNS) ROUGHNESS      = 0._EB
SPEC_ID                 = 'null'
SPREAD_RATE             = 0.05_EB
STRETCH_FACTOR          = 2._EB
SURFACE_DENSITY         = -1._EB
TAU_MF                  =  1._EB
TAU_Q                   = 1._EB
TAU_V                   = 1._EB
TAU_T                   = 1._EB
TAU_PART                = 0._EB
TAU_EXTERNAL_FLUX       = 0.001_EB
TEXTURE_MAP             = 'null'
TEXTURE_WIDTH           = 1._EB
TEXTURE_HEIGHT          = 1._EB
THICKNESS               = -1._EB
TMP_BACK                = -TMPM-1._EB
TMP_FRONT               = -TMPM-1._EB
TMP_INNER               = -TMPM-1._EB
TRANSPARENCY            = 1._EB
VEL_T                   = 0._EB
VEL                     = 0._EB
VOLUME_FLUX             = 0._EB
WIDTH                   = -1._EB
XYZ                     = -1.E6_EB
Z0                      = 10._EB

VEGETATION                   = .FALSE.
VEGETATION_NO_BURN           = .FALSE.
VEGETATION_CDRAG             = 1.0_EB
VEGETATION_CHAR_FRACTION     = 0.20_EB
VEGETATION_ELEMENT_DENSITY   = 512._EB !kg/m^3
VEGETATION_HEIGHT            = 0.50_EB !m
VEGETATION_INITIAL_TEMP      = TMPA-TMPM
VEGETATION_GROUND_TEMP       = -99._EB
VEGETATION_LOAD              = 0.30_EB !kg/m^2
FIRELINE_MLR_MAX             = 999. !kg/m/s w*R*(1-ChiChar) 
!SRF_VEG_DEHYDRATION_RATE_MAX = 999. !kg/m^2/s   
VEGETATION_LAYERS            = 0
VEGETATION_MOISTURE          = 0.06_EB
VEGETATION_SVRATIO           = 12000_EB !1/m
VEGETATION_LSET_IGNITE_TIME  = 1.0E20_EB
VEGETATION_LINEAR_DEGRAD     = .TRUE.
VEGETATION_ARRHENIUS_DEGRAD  = .FALSE.
VEG_LSET_ROS_HEAD            = 0.0_EB
VEG_LSET_ROS_FLANK           = 0.0_EB
VEG_LSET_ROS_BACK            = 0.0_EB
VEG_LEVEL_SET_SPREAD         = .FALSE.
VEG_LSET_WIND_EXP            = 1.0_EB
 
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

   DO NSPC=0,N_GAS_SPECIES
      IF (TRIM(SF%RAMP_MF(NSPC))/='null') THEN
         CALL GET_RAMP_INDEX(SF%RAMP_MF(NSPC),'TIME',NR)
         SF%RAMP_INDEX(NSPC) = NR
      ELSE
         IF (SF%TAU(NSPC) > 0._EB) SF%RAMP_INDEX(NSPC) = TANH_RAMP 
         IF (SF%TAU(NSPC) < 0._EB) SF%RAMP_INDEX(NSPC) = TSQR_RAMP 
      ENDIF
   ENDDO 
   
   DO IPC=1,N_PART
      PC=>PARTICLE_CLASS(IPC)
      IF (PC%SURF_ID==SF%ID) THEN
         PC%SURF_INDEX = N
         IF (.NOT.SF%THERMALLY_THICK .AND. SF%RADIUS<0._EB) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: SURF ',TRIM(SF%ID),' needs a RADIUS'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDIF
   ENDDO

ENDDO PROCESS_SURF_LOOP   
 
! If a particle class uses a SURF line, make sure the SURF ID exists

DO IPC=1,N_PART
   PC=>PARTICLE_CLASS(IPC)
   IF (PC%SURF_ID/='null' .AND. PC%SURF_INDEX<1) THEN
      WRITE(MESSAGE,'(A,A,A)') 'ERROR: SURF ',TRIM(PC%SURF_ID),' not found'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO

END SUBROUTINE PROC_SURF_1



SUBROUTINE PROC_SURF_2

! Go through the SURF types and process
 
INTEGER :: IPC,N,NN,NNN,NL,NSPC
REAL(EB) :: ADJUSTED_LAYER_DENSITY
INTEGER  :: IVEG_L,IIVEG_L
REAL(EB) :: DETA_VEG,DZVEG_L,ETA_H,ETAFM_VEG,ETAFP_VEG 
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

   ! Make decisions based on whether there is forced ventilation at the surface

   BLOWING  = .FALSE.
   SUCKING  = .FALSE.
   IF (SF%VEL<0._EB .OR. SF%VOLUME_FLUX<0._EB .OR. SF%MASS_FLUX_TOTAL < 0._EB) BLOWING = .TRUE.
   IF (SF%VEL>0._EB .OR. SF%VOLUME_FLUX>0._EB .OR. SF%MASS_FLUX_TOTAL > 0._EB) SUCKING = .TRUE.
   IF (BLOWING .OR. SUCKING) SF%SPECIFIED_NORMAL_VELOCITY = .TRUE.
   IF (SUCKING) SF%FREE_SLIP = .TRUE.

   IF (BURNING .AND. (BLOWING .OR. SUCKING)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' cannot have a specified velocity or volume flux'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
 
   ! set predefined HRRPUA

   BURNING_IF: IF (BURNING .AND. .NOT.ALL(EVACUATION_ONLY)) THEN
      IF (SF%HRRPUA>0._EB) THEN
         IF (CO_PRODUCTION) THEN
            RN => REACTION(2)         
         ELSE
            RN => REACTION(1)
         ENDIF
         SF%MASS_FLUX(I_FUEL) = SF%HRRPUA/ (RN%HEAT_OF_COMBUSTION*SPECIES_MIXTURE(2)%MASS_FRACTION(FUEL_INDEX))
      ENDIF
      IF (SF%MLRPUA>0._EB) SF%MASS_FLUX(I_FUEL) = SF%MLRPUA
      ! Adjust burning rate according to the difference of heats of combustion
      IF (SF%N_LAYERS > 0) THEN
         ML => MATERIAL(SF%MATL_INDEX(1))
         SF%ADJUST_BURN_RATE(I_FUEL) = ML%ADJUST_BURN_RATE(1,I_FUEL)
         SF%MASS_FLUX(I_FUEL)        = SF%MASS_FLUX(I_FUEL)/SF%ADJUST_BURN_RATE(I_FUEL)
      ENDIF
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

   IF ((ABS(SF%SURFACE_DENSITY) <=ZERO_P) .AND. SF%BURN_AWAY) THEN
      WRITE(MESSAGE,'(A,A,A)') 'WARNING: SURF ',TRIM(SF%ID),' has BURN_AWAY set but zero combustible density'
      IF (MYID==0) WRITE(LU_ERR,'(A)') TRIM(MESSAGE)
   ENDIF   

   ! Ignition Time

   SF%T_IGN = T_BEGIN 
   IF (SF%TMP_IGN<5000._EB)                    SF%T_IGN = HUGE(T_END)
   IF (SF%PYROLYSIS_MODEL==PYROLYSIS_MATERIAL) SF%T_IGN = HUGE(T_END)

   ! Species Arrays and Method of Mass Transfer (SPECIES_BC_INDEX)

   SF%SPECIES_BC_INDEX = NO_MASS_FLUX

   IF (ANY(SF%MASS_FRACTION>=0._EB) .AND. (ANY(ABS(SF%MASS_FLUX)>ZERO_P) .OR. SF%PYROLYSIS_MODEL/= PYROLYSIS_NONE)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' cannot specify mass fraction with mass flux and/or pyrolysis'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
   IF (ANY(SF%MASS_FRACTION>=0._EB) .AND. SUCKING) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' cannot specify both mass fraction and outflow velocity'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
   IF (ANY(SF%LEAK_PATH>=0) .AND. (BLOWING .OR. SUCKING .OR. SF%PYROLYSIS_MODEL/= PYROLYSIS_NONE)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' cannot leak and specify flow or pyrolysis at the same time'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
   IF (ANY(ABS(SF%MASS_FLUX)>ZERO_P) .AND. (BLOWING .OR. SUCKING)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' cannot have both a mass flux and specified velocity'
      CALL SHUTDOWN(MESSAGE)
      ENDIF

   IF (BLOWING .OR. SUCKING)         SF%SPECIES_BC_INDEX = SPECIFIED_MASS_FRACTION
   IF (ANY(SF%MASS_FRACTION>=0._EB)) SF%SPECIES_BC_INDEX = SPECIFIED_MASS_FRACTION
   IF (ANY(ABS(SF%MASS_FLUX)>ZERO_P) .OR. &
       SF%PYROLYSIS_MODEL==PYROLYSIS_MATERIAL) SF%SPECIES_BC_INDEX = SPECIFIED_MASS_FLUX

   IF (SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FRACTION) THEN
      IF (ALL(SF%MASS_FRACTION< 0._EB)) THEN
         SF%MASS_FRACTION(1:N_GAS_SPECIES) = YY_AMBIENT(1:N_GAS_SPECIES)
      ELSE
         SPECIES_LOOP: DO NSPC=1,N_GAS_SPECIES
            IF (SF%MASS_FRACTION(NSPC)<0._EB) SF%MASS_FRACTION(NSPC) = 0._EB
         ENDDO SPECIES_LOOP
         IF (SUM(SF%MASS_FRACTION(1:N_GAS_SPECIES)) > 1._EB) THEN
               WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' sum of specified MASS_FRACTION > 1'
               CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDIF
      SF%MASS_FRACTION(0) = 1._EB - SUM(SF%MASS_FRACTION(1:N_GAS_SPECIES))
   ENDIF

   ! Boundary fuel model for vegetation

   IF (SF%VEGETATION) SF%SPECIES_BC_INDEX = SPECIFIED_MASS_FLUX
 
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
      SF%EMISSIVITY = 1._EB      
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
   IF (N==HVAC_SURF_INDEX) THEN
      SF%THERMAL_BC_INDEX = HVAC_BOUNDARY
      SF%SPECIES_BC_INDEX = HVAC_BOUNDARY
   ENDIF
   IF (N==EVACUATION_SURF_INDEX) THEN
      SF%THERMAL_BC_INDEX = INFLOW_OUTFLOW
      SF%SPECIES_BC_INDEX = SPECIFIED_MASS_FRACTION
      SF%SPECIFIED_NORMAL_VELOCITY = .TRUE.
      SF%FREE_SLIP = .TRUE.
      SF%VELOCITY_BC_INDEX = FREE_SLIP_BC
      SF%VEL                   = +0.000001_EB ! VEL
      SF%TAU(TIME_VELO)        = 0.1_EB! TAU_V
      SF%RAMP_INDEX(TIME_VELO) = TANH_RAMP
   ENDIF

! Boundary surface vegetation

   DZVEG_L  = SF%VEG_HEIGHT/REAL(SF%NVEG_L,EB)
   DETA_VEG = SF%VEG_KAPPA*DZVEG_L 

! Factors for computing decay of +/- incident fluxes
   SF%VEG_FINCM_RADFCT_L(:) =  0.0_EB
   SF%VEG_FINCP_RADFCT_L(:) =  0.0_EB
   ETA_H = SF%VEG_KAPPA*SF%VEG_HEIGHT
   DO IVEG_L = 0,SF%NVEG_L
    ETAFM_VEG = IVEG_L*DETA_VEG
    ETAFP_VEG = ETA_H - ETAFM_VEG
    SF%VEG_FINCM_RADFCT_L(IVEG_L) = EXP(-ETAFM_VEG)
    SF%VEG_FINCP_RADFCT_L(IVEG_L) = EXP(-ETAFP_VEG)
   ENDDO

!  Integrand for computing +/- self emission fluxes
   SF%VEG_SEMISSP_RADFCT_L(:,:) = 0.0_EB
   SF%VEG_SEMISSM_RADFCT_L(:,:) = 0.0_EB
! q+
   DO IIVEG_L = 0,SF%NVEG_L !grid coordinate
    DO IVEG_L = IIVEG_L,SF%NVEG_L !integrand index
     ETAFM_VEG = (IVEG_L-IIVEG_L)*DETA_VEG
     ETAFP_VEG = ETAFM_VEG + DETA_VEG
     SF%VEG_SEMISSP_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-ETAFM_VEG) - EXP(-ETAFP_VEG)
    ENDDO
   ENDDO
! q-
   DO IIVEG_L = 0,SF%NVEG_L
    DO IVEG_L = 1,IIVEG_L
     ETAFM_VEG = (IIVEG_L-IVEG_L)*DETA_VEG
     ETAFP_VEG = ETAFM_VEG + DETA_VEG
     SF%VEG_SEMISSM_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-ETAFM_VEG) - EXP(-ETAFP_VEG)
    ENDDO
   ENDDO

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

USE SCARC_SOLVER, ONLY: SCARC_METHOD   , SCARC_DEBUG     , SCARC_CASE , &
                        SCARC_EPS_REL  , SCARC_EPS_DIVG  , SCARC_BREL , &
                        SCARC_MG_NIT   , SCARC_CG_NIT    , SCARC_BICG_NIT    , SCARC_SM_NIT    , SCARC_CO_NIT    , &
                        SCARC_MG_EPS   , SCARC_CG_EPS    , SCARC_BICG_EPS    , SCARC_SM_EPS    , SCARC_CO_EPS    , &
                        SCARC_MG_PRECON, SCARC_CG_PRECON , SCARC_BICG_PRECON , SCARC_SM_PRECON , SCARC_CO_PRECON , &
                        SCARC_MG_OMEGA , SCARC_CG_OMEGA  , SCARC_BICG_OMEGA  , SCARC_SM_OMEGA  , SCARC_CO_OMEGA  , &
                        SCARC_MG_NLDIFF

NAMELIST /PRES/ VELOCITY_TOLERANCE,MAX_PRESSURE_ITERATIONS,PRESSIT_ACCELERATOR,PRESSIT_SCALE_FACTOR, &
                SCARC_METHOD   , SCARC_DEBUG     , SCARC_CASE , &
                SCARC_EPS_REL  , SCARC_EPS_DIVG  , SCARC_BREL , &
                SCARC_MG_NIT   , SCARC_CG_NIT    , SCARC_BICG_NIT    , SCARC_SM_NIT    , SCARC_CO_NIT    , &
                SCARC_MG_EPS   , SCARC_CG_EPS    , SCARC_BICG_EPS    , SCARC_SM_EPS    , SCARC_CO_EPS    , &
                SCARC_MG_PRECON, SCARC_CG_PRECON , SCARC_BICG_PRECON , SCARC_SM_PRECON , SCARC_CO_PRECON , &
                SCARC_MG_OMEGA , SCARC_CG_OMEGA  , SCARC_BICG_OMEGA  , SCARC_SM_OMEGA  , SCARC_CO_OMEGA  , &
                SCARC_MG_NLDIFF

REWIND(LU_INPUT)
READ_LOOP: DO
   CALL CHECKREAD('PRES',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_LOOP
   READ(LU_INPUT,PRES,END=23,ERR=24,IOSTAT=IOS)
   24 IF (IOS>0) THEN
      CALL SHUTDOWN('ERROR: Problem with PRES line')
   ENDIF
ENDDO READ_LOOP
23 REWIND(LU_INPUT)

IF (VELOCITY_TOLERANCE<1.E4) THEN
   ITERATE_PRESSURE = .TRUE.
ENDIF

IF (SCARC_METHOD /= 'null') PRES_METHOD = 'SCARC'

END SUBROUTINE READ_PRES


SUBROUTINE READ_RADI
USE RADCONS
NAMELIST /RADI/ TIME_STEP_INCREMENT,NUMBER_RADIATION_ANGLES,ANGLE_INCREMENT,KAPPA0, RADIATION,&
                WIDE_BAND_MODEL,CH4_BANDS,PATH_LENGTH,NMIEANG,RADTMP,RADIATIVE_FRACTION,RADCAL_FUEL

IF (LES) RADIATIVE_FRACTION = 0.35_EB
IF (DNS) RADIATIVE_FRACTION = 0.00_EB
NUMBER_RADIATION_ANGLES = 100
TIME_STEP_INCREMENT     = 3
IF (TWO_D) THEN
   NUMBER_RADIATION_ANGLES = 60
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
      CALL SHUTDOWN('ERROR: Problem with RADI line')
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
      IF (RP%DEVC_ID =='null') RP%DEVC_ID = DEVC_ID
      IF (RP%CTRL_ID =='null') RP%CTRL_ID = CTRL_ID      
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
   ALLOCATE(RAMPS(N)%INTERPOLATED_DATA(0:RP%NUMBER_INTERPOLATION_POINTS+1))
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
   RP%INTERPOLATED_DATA(RP%NUMBER_INTERPOLATION_POINTS)   = RP%DEPENDENT_DATA(RP%NUMBER_DATA_POINTS)
   RP%INTERPOLATED_DATA(RP%NUMBER_INTERPOLATION_POINTS+1) = RP%DEPENDENT_DATA(RP%NUMBER_DATA_POINTS)

   ! Get Device or Control Index

   CALL SEARCH_CONTROLLER('RAMP',RP%CTRL_ID,RP%DEVC_ID,RP%DEVC_INDEX,RP%CTRL_INDEX,N)
   
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
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad 1st lattitude'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(2)<TABLE_DATA(1).OR. TABLE_DATA(2)>180) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad 2nd lattitude'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(3)<-180. .OR.        TABLE_DATA(3)>360) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad 1st longitude'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(4)<TABLE_DATA(3).OR. TABLE_DATA(4)>360) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad 2nd longitude'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(5)<0) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad velocity'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(6)<0) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad mass flow'
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
CHARACTER(30) :: ID,DEVC_ID,PROP_ID,SURF_ID,SURF_IDS(3),SURF_ID6(6),CTRL_ID,MULT_ID
CHARACTER(60) :: MESH_ID
CHARACTER(25) :: COLOR
LOGICAL :: EVACUATION_OBST
REAL(EB) :: TRANSPARENCY,XB1,XB2,XB3,XB4,XB5,XB6,BULK_DENSITY,VOL_ADJUSTED,VOL_SPECIFIED
LOGICAL :: SAWTOOTH,EMBEDDED,THICKEN,PERMIT_HOLE,ALLOW_VENT,EVACUATION, REMOVABLE,BNDF_FACE(-3:3),BNDF_OBST,OUTLINE
NAMELIST /OBST/ XB,ID,PROP_ID,SURF_ID,SURF_IDS,SURF_ID6,FYI,BNDF_FACE,BNDF_OBST, &
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
 
   IF (EVACUATION_ONLY(NM)) CALL DEFINE_EVACUATION_OBSTS(NM,1)

   ! Allocate OBSTRUCTION array

   ALLOCATE(M%OBSTRUCTION(0:N_OBST),STAT=IZERO)
   CALL ChkMemErr('READ','OBSTRUCTION',IZERO)
   OBSTRUCTION=>M%OBSTRUCTION
 
   N        = 0
   N_OBST_O = N_OBST
 
   READ_OBST_LOOP: DO NN=1,N_OBST_O

      ID       = 'null'
      MULT_ID  = 'null'
      PROP_ID  = 'null'
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
      IF (     EVACUATION_ONLY(NM)) REMOVABLE  = .FALSE.
 
      ! Read the OBST line

      EVACUATION_OBST = .FALSE.
      IF (EVACUATION_ONLY(NM)) CALL DEFINE_EVACUATION_OBSTS(NM,2)
      EVACUATION_OBSTS: IF (.NOT. EVACUATION_OBST) THEN
         CALL CHECKREAD('OBST',LU_INPUT,IOS)
         IF (IOS==1) EXIT READ_OBST_LOOP
         READ(LU_INPUT,OBST,END=35)
      END IF EVACUATION_OBSTS

      ! Reorder OBST coordinates if necessary

      CALL CHECK_XB(XB)

      IF (EVACUATION_ONLY(NM)) THEN
         DEVC_ID    = 'null'
         CTRL_ID    = 'null'
         PROP_ID    = 'null'
      END IF

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

               ! Thicken evacuation mesh obstructions in the z direction

               IF (EVACUATION_ONLY(NM) .AND. EVACUATION) THEN
                  OB%Z1 = ZS
                  OB%Z2 = ZF
                  XB5   = ZS
                  XB6   = ZF
               ENDIF

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
          
               NNNN = 0
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
                  IF (TRIM(SURFACE(NNN)%ID)==TRIM(EVAC_SURF_DEFAULT)) NNNN = NNN
               ENDDO

               ! Fire + evacuation calculation: draw obsts as outlines by default

               IF (.NOT.OUTLINE .AND. EVACUATION_ONLY(NM) .AND. .NOT.ALL(EVACUATION_ONLY)) THEN
                  IF (SURFACE(NNNN)%TRANSPARENCY < 0.99999_EB .AND. .NOT.OUTLINE) THEN
                     OUTLINE = .FALSE.
                  ELSE
                     OUTLINE = .TRUE.
                  ENDIF
               ENDIF
         
               ! Determine if the OBST is CONSUMABLE and check if POROUS inappropriately applied
         
               FACE_LOOP: DO NNN=-3,3
                  IF (NNN==0) CYCLE FACE_LOOP
                  IF (SURFACE(OB%IBC(NNN))%BURN_AWAY) THEN
                     OB%CONSUMABLE = .TRUE.
                     IF (.NOT.SAWTOOTH) THEN
                        IF (ID=='null')WRITE(MESSAGE,'(A,I5,A)')'ERROR: OBST ',N,       ' cannot be BURN_AWAY and SAWTOOTH=.FALSE.'
                        IF (ID/='null')WRITE(MESSAGE,'(A,A,A)') 'ERROR: OBST ',TRIM(ID),' cannot be BURN_AWAY and SAWTOOTH=.FALSE.'
                        CALL SHUTDOWN(MESSAGE)
                     ENDIF
                  ENDIF
                  IF (SURFACE(OB%IBC(NNN))%POROUS .AND. .NOT.OB%THIN) THEN
                     IF (ID=='null') WRITE(MESSAGE,'(A,I5,A)') 'ERROR: OBST ',N,       ' must be zero cells thick to be POROUS'
                     IF (ID/='null') WRITE(MESSAGE,'(A,A,A)')  'ERROR: OBST ',TRIM(ID),' must be zero cells thick to be POROUS'
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
          
               ! Property ID
          
               OB%PROP_ID = PROP_ID
         
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

               IF (EVACUATION_ONLY(NM)) BULK_DENSITY = -1._EB
               OB%BULK_DENSITY = BULK_DENSITY
               IF (OB%VOLUME_ADJUST==0._EB .AND. OB%BULK_DENSITY>0._EB) OB%BULK_DENSITY = -1._EB
             !    IF (ID=='null') WRITE(MESSAGE,'(A,I4,A)') 'ERROR: OBST ',NN,      ' has no volume and cannot have a BULK_DENSITY'
             !    IF (ID/='null') WRITE(MESSAGE,'(A,A,A)')  'ERROR: OBST ',TRIM(ID),' has no volume and cannot have a BULK_DENSITY'
             !    CALL SHUTDOWN(MESSAGE)
             ! ENDIF
 
               ! Make obstruction invisible if it's within a finer mesh
 
               DO NOM=1,NM-1
                  IF (EVACUATION_ONLY(NOM)) CYCLE
                  IF (EVACUATION_ONLY(NM)) CYCLE
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
               IF (EVACUATION_ONLY(NM)) OB%SHOW_BNDF(:) = .FALSE.
 
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
   SOLID=>M%SOLID
   ALLOCATE(M%OBST_INDEX_C(0:CELL_COUNT),STAT=IZERO) 
   CALL ChkMemErr('READ','OBST_INDEX_C',IZERO) 
   M%OBST_INDEX_C = 0
   OBST_INDEX_C=>M%OBST_INDEX_C 

   ! Make all exterior cells solid
 
   CALL BLOCK_CELL(NM,   0,   0,   0,JBP1,   0,KBP1,1,0)
   CALL BLOCK_CELL(NM,IBP1,IBP1,   0,JBP1,   0,KBP1,1,0)
   CALL BLOCK_CELL(NM,   0,IBP1,   0,   0,   0,KBP1,1,0)
   CALL BLOCK_CELL(NM,   0,IBP1,JBP1,JBP1,   0,KBP1,1,0)
   CALL BLOCK_CELL(NM,   0,IBP1,   0,JBP1,   0,   0,1,0)
   CALL BLOCK_CELL(NM,   0,IBP1,   0,JBP1,KBP1,KBP1,1,0)

   ! Block off cells filled by obstructions

   DO N=1,N_OBST
      OB=>OBSTRUCTION(N)
      CALL BLOCK_CELL(NM,OB%I1+1,OB%I2,OB%J1+1,OB%J2,OB%K1+1,OB%K2,1,N)
   ENDDO
 
   ! EVAC_PROCESS needs only the allocations of SOLID, CELL_INDEX, OBST_INDEX_C and OBSTRUCTION

   IF (.NOT.EVACUATION_ONLY(NM) .AND. MYID==EVAC_PROCESS) THEN
      DEALLOCATE(M%EXTERIOR)
      CYCLE MESH_LOOP_3
   ENDIF

   ! Create arrays to hold cell indices

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

CONTAINS

  SUBROUTINE DEFINE_EVACUATION_OBSTS(NM,IMODE)
    !
    ! Define the evacuation OBSTs for the doors/exits, if needed.  A VENT should always
    ! be defined on an OBST that is at least one grid cell thick or the VENT should be
    ! on the outer boundary of the evacuation mesh, which is by default solid.
    ! The core of the STRS meshes are also defined.
    !
    USE EVAC, ONLY: N_DOORS, N_EXITS, N_CO_EXITS, EVAC_EMESH_EXITS_TYPE, EMESH_EXITS, N_DOOR_MESHES, EVAC_FDS6, &
         N_STRS, EMESH_STAIRS, EVAC_EMESH_STAIRS_TYPE
    IMPLICIT NONE
    ! Passed variables
    INTEGER, INTENT(IN) :: NM, IMODE
    ! Local variables
    INTEGER :: N, N_END, I, I1, I2, J1, J2
    REAL(EB) :: TINY

    TINY = 0.1_EB*MIN(MESHES(NM)%DXI, MESHES(NM)%DETA)
    N_END = N_EXITS - N_CO_EXITS + N_DOORS
    IMODE_1_IF: IF (IMODE==1) THEN
       NEND_LOOP_1: DO N = 1, N_END

          IF (.NOT.EMESH_EXITS(N)%DEFINE_MESH) CYCLE NEND_LOOP_1
          IF (EMESH_EXITS(N)%IMESH==NM .OR. EMESH_EXITS(N)%MAINMESH==NM) THEN
             EMESH_EXITS(N)%I_OBST = 0

             ! Move EMESH_EXITS(N)%XB to mesh cell boundaries
             EMESH_EXITS(N)%XB(1) = MAX(EMESH_EXITS(N)%XB(1),MESHES(NM)%XS)
             EMESH_EXITS(N)%XB(2) = MIN(EMESH_EXITS(N)%XB(2),MESHES(NM)%XF)
             EMESH_EXITS(N)%XB(3) = MAX(EMESH_EXITS(N)%XB(3),MESHES(NM)%YS)
             EMESH_EXITS(N)%XB(4) = MIN(EMESH_EXITS(N)%XB(4),MESHES(NM)%YF)
             
             I1 = NINT(GINV(EMESH_EXITS(N)%XB(1)-MESHES(NM)%XS,1,NM)*MESHES(NM)%RDXI ) 
             I2 = NINT(GINV(EMESH_EXITS(N)%XB(2)-MESHES(NM)%XS,1,NM)*MESHES(NM)%RDXI )
             J1 = NINT(GINV(EMESH_EXITS(N)%XB(3)-MESHES(NM)%YS,2,NM)*MESHES(NM)%RDETA) 
             J2 = NINT(GINV(EMESH_EXITS(N)%XB(4)-MESHES(NM)%YS,2,NM)*MESHES(NM)%RDETA)
             
             EMESH_EXITS(N)%XB(1) = MESHES(NM)%X(I1)
             EMESH_EXITS(N)%XB(2) = MESHES(NM)%X(I2)
             EMESH_EXITS(N)%XB(3) = MESHES(NM)%Y(J1)
             EMESH_EXITS(N)%XB(4) = MESHES(NM)%Y(J2)

             ! Check if the exit/door is at the mesh boundary, then no OBST is needed.
             SELECT CASE (EMESH_EXITS(N)%IOR)
             CASE (-1)
                IF (MESHES(NM)%XS >= EMESH_EXITS(N)%XB(1)-TINY) CYCLE NEND_LOOP_1
             CASE (+1)
                IF (MESHES(NM)%XF <= EMESH_EXITS(N)%XB(2)+TINY) CYCLE NEND_LOOP_1
             CASE (-2)
                IF (MESHES(NM)%YS >= EMESH_EXITS(N)%XB(3)-TINY) CYCLE NEND_LOOP_1
             CASE (+2)
                IF (MESHES(NM)%YF <= EMESH_EXITS(N)%XB(4)+TINY) CYCLE NEND_LOOP_1
             END SELECT
             N_OBST = N_OBST + 1
             EMESH_EXITS(N)%I_OBST = N_OBST
             EVACUATION_OBST = .TRUE.
             IF (.NOT.EVAC_FDS6 .AND. EMESH_EXITS(N)%IMESH==NM) EXIT NEND_LOOP_1  ! One VENT per door flow mesh
          END IF
       END DO NEND_LOOP_1

       NSTRS_LOOP_1: DO N = 1, N_STRS

          IF (.NOT.EMESH_STAIRS(N)%DEFINE_MESH) CYCLE NSTRS_LOOP_1
          IF (.NOT.EVAC_FDS6) EXIT  NSTRS_LOOP_1  ! Automatic STRS meshes are for FDS6
          IF (EMESH_STAIRS(N)%IMESH==NM) THEN

             ! Move EMESH_STAIRS(N)%XB_CORE to mesh cell boundaries
             EMESH_STAIRS(N)%XB_CORE(1) = MAX(EMESH_STAIRS(N)%XB_CORE(1),MESHES(NM)%XS)
             EMESH_STAIRS(N)%XB_CORE(2) = MIN(EMESH_STAIRS(N)%XB_CORE(2),MESHES(NM)%XF)
             EMESH_STAIRS(N)%XB_CORE(3) = MAX(EMESH_STAIRS(N)%XB_CORE(3),MESHES(NM)%YS)
             EMESH_STAIRS(N)%XB_CORE(4) = MIN(EMESH_STAIRS(N)%XB_CORE(4),MESHES(NM)%YF)
             
             I1 = NINT(GINV(EMESH_STAIRS(N)%XB_CORE(1)-MESHES(NM)%XS,1,NM)*MESHES(NM)%RDXI ) 
             I2 = NINT(GINV(EMESH_STAIRS(N)%XB_CORE(2)-MESHES(NM)%XS,1,NM)*MESHES(NM)%RDXI )
             J1 = NINT(GINV(EMESH_STAIRS(N)%XB_CORE(3)-MESHES(NM)%YS,2,NM)*MESHES(NM)%RDETA) 
             J2 = NINT(GINV(EMESH_STAIRS(N)%XB_CORE(4)-MESHES(NM)%YS,2,NM)*MESHES(NM)%RDETA)
             
             EMESH_STAIRS(N)%XB_CORE(1) = MESHES(NM)%X(I1)
             EMESH_STAIRS(N)%XB_CORE(2) = MESHES(NM)%X(I2)
             EMESH_STAIRS(N)%XB_CORE(3) = MESHES(NM)%Y(J1)
             EMESH_STAIRS(N)%XB_CORE(4) = MESHES(NM)%Y(J2)

             N_OBST = N_OBST + 1
             EMESH_STAIRS(N)%I_OBST = N_OBST
             EVACUATION_OBST = .TRUE.
          END IF
       END DO NSTRS_LOOP_1

    END IF IMODE_1_IF

    IMODE_2_IF: IF (IMODE==2) THEN
       NEND_LOOP_2: DO N = 1, N_END
          IF (.NOT.EMESH_EXITS(N)%DEFINE_MESH) CYCLE NEND_LOOP_2
          IF (EMESH_EXITS(N)%I_OBST==NN .AND. (EMESH_EXITS(N)%IMESH==NM .OR. EMESH_EXITS(N)%MAINMESH==NM)) THEN
             EVACUATION_OBST = .TRUE.
             EVACUATION = .TRUE.
             REMOVABLE = .FALSE.
             THICKEN = .TRUE.
             PERMIT_HOLE = .FALSE.
             ALLOW_VENT = .TRUE.
             MESH_ID = TRIM(MESH_NAME(NM))
             XB(1) = EMESH_EXITS(N)%XB(1)
             XB(2) = EMESH_EXITS(N)%XB(2)
             XB(3) = EMESH_EXITS(N)%XB(3)
             XB(4) = EMESH_EXITS(N)%XB(4)
             XB(5) = EMESH_EXITS(N)%XB(5)
             XB(6) = EMESH_EXITS(N)%XB(6)
             SELECT CASE (EMESH_EXITS(N)%IOR)
             CASE (-1)
                XB(1) = MAX(MESHES(NM)%XS, XB(1) - 0.49_EB*MESHES(NM)%DXI)
             CASE (+1)
                XB(2) = MIN(MESHES(NM)%XF, XB(2) + 0.49_EB*MESHES(NM)%DXI)
             CASE (-2)
                XB(3) = MAX(MESHES(NM)%YS, XB(3) - 0.49_EB*MESHES(NM)%DETA)
             CASE (+2)
                XB(4) = MIN(MESHES(NM)%YF, XB(4) + 0.49_EB*MESHES(NM)%DETA)
             END SELECT
             RGB(:) = EMESH_EXITS(N)%RGB(:)
             ID = TRIM('Eobst_' // TRIM(MESH_NAME(NM)))
             IF (.NOT.EVAC_FDS6 .AND. EMESH_EXITS(N)%IMESH==NM) EXIT NEND_LOOP_2  ! One VENT per door flow mesh
          END IF
       END DO NEND_LOOP_2

       NSTRS_LOOP_2: DO N = 1, N_STRS
          IF (.NOT.EMESH_STAIRS(N)%DEFINE_MESH) CYCLE NSTRS_LOOP_2
          IF (EMESH_STAIRS(N)%I_OBST==NN .AND. EMESH_STAIRS(N)%IMESH==NM) THEN
             EVACUATION_OBST = .TRUE.
             EVACUATION = .TRUE.
             REMOVABLE = .FALSE.
             ! THICKEN = .TRUE.
             PERMIT_HOLE = .FALSE.
             ALLOW_VENT = .FALSE.
             MESH_ID = TRIM(MESH_NAME(NM))
             XB(1) = EMESH_STAIRS(N)%XB_CORE(1)
             XB(2) = EMESH_STAIRS(N)%XB_CORE(2)
             XB(3) = EMESH_STAIRS(N)%XB_CORE(3)
             XB(4) = EMESH_STAIRS(N)%XB_CORE(4)
             XB(5) = EMESH_STAIRS(N)%XB(5)
             XB(6) = EMESH_STAIRS(N)%XB(6)
             RGB(:) = EMESH_STAIRS(N)%RGB(:)
             ID = TRIM('Eobst_' // TRIM(MESH_NAME(NM)))
          END IF
       END DO NSTRS_LOOP_2

    END IF IMODE_2_IF

    RETURN
  END SUBROUTINE DEFINE_EVACUATION_OBSTS

END SUBROUTINE READ_OBST



SUBROUTINE READ_HOLE

USE CONTROL_VARIABLES, ONLY : CONTROl, N_CTRL
USE DEVICE_VARIABLES, ONLY : DEVICE, N_DEVC
CHARACTER(30) :: DEVC_ID,CTRL_ID,MULT_ID
CHARACTER(60) :: MESH_ID
CHARACTER(25) :: COLOR
LOGICAL :: EVACUATION_HOLE
LOGICAL :: EVACUATION
INTEGER :: NM,N_HOLE,NN,NDO,N,I1,I2,J1,J2,K1,K2,RGB(3),N_HOLE_NEW,N_HOLE_O,II,JJ,KK,NNN
REAL(EB) :: X1,X2,Y1,Y2,Z1,Z2,TRANSPARENCY
NAMELIST /HOLE/ XB,FYI,RGB,TRANSPARENCY,EVACUATION,MESH_ID,COLOR,DEVC_ID,CTRL_ID,MULT_ID
TYPE(OBSTRUCTION_TYPE), ALLOCATABLE, DIMENSION(:) :: TEMP_OBST
TYPE(MULTIPLIER_TYPE), POINTER :: MR
LOGICAL, ALLOCATABLE, DIMENSION(:) :: TEMP_HOLE_EVAC

ALLOCATE(TEMP_OBST(0:6))

N_HOLE    = 0
N_HOLE_O  = 0
REWIND(LU_INPUT)

COUNT_LOOP: DO
   CALL CHECKREAD('HOLE',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_LOOP
   MULT_ID = 'null'
   READ(LU_INPUT,NML=HOLE,END=1,ERR=2,IOSTAT=IOS)
   N_HOLE_NEW = 1
   IF (MULT_ID/='null') THEN
      DO N=1,N_MULT
         MR => MULTIPLIER(N)
         IF (MULT_ID==MR%ID) N_HOLE_NEW = MR%N_COPIES
      ENDDO
   ENDIF
   N_HOLE_O = N_HOLE_O + 1
   N_HOLE   = N_HOLE   + N_HOLE_NEW
   2 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I5)')  'ERROR: Problem with HOLE number',N_HOLE_O+1
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO COUNT_LOOP
1 REWIND(LU_INPUT)

CALL DEFINE_EVACUATION_HOLES(1)
! TEMP_HOLE_EVAC(:) indicates if the given HOLE is to be used in the EVACUATION routine

IF (ANY(EVACUATION_ONLY)) THEN 
   ALLOCATE(TEMP_HOLE_EVAC(1:N_HOLE_O))
   READ_HOLE_EVAC_LOOP: DO N=1,N_HOLE_O
      EVACUATION_HOLE = .FALSE.
      CALL DEFINE_EVACUATION_HOLES(2)
      EVACUATION = .TRUE.
      IF (.NOT.EVACUATION_HOLE) THEN
         CALL CHECKREAD('HOLE',LU_INPUT,IOS)
         IF (IOS==1) EXIT READ_HOLE_EVAC_LOOP
         READ(LU_INPUT,HOLE)
      END IF
      TEMP_HOLE_EVAC(N) = EVACUATION
   ENDDO READ_HOLE_EVAC_LOOP
   REWIND(LU_INPUT)
ENDIF
 
READ_HOLE_LOOP: DO N=1,N_HOLE_O
 
   ! Set default values for the HOLE namelist parameters

   DEVC_ID  = 'null'
   CTRL_ID  = 'null'
   MESH_ID  = 'null'
   MULT_ID  = 'null'
   COLOR    = 'null'
   RGB      = -1
   TRANSPARENCY  = 1._EB
   EVACUATION = .FALSE.

   ! Read the HOLE line
 
   EVACUATION_HOLE = .FALSE.
   IF (ANY(EVACUATION_ONLY)) CALL DEFINE_EVACUATION_HOLES(3)
   EVACUATION_HOLES: IF (.NOT. EVACUATION_HOLE) THEN
      CALL CHECKREAD('HOLE',LU_INPUT,IOS)
      IF (IOS==1) EXIT READ_HOLE_LOOP
      READ(LU_INPUT,HOLE)
   END IF EVACUATION_HOLES
 
   ! Re-order coordinates, if necessary

   CALL CHECK_XB(XB)

   IF (ALL(EVACUATION_ONLY)) THEN
      DEVC_ID    = 'null'
      CTRL_ID    = 'null'
   END IF

   ! Loop over all the meshes to determine where the HOLE is
 
   MESH_LOOP: DO NM=1,NMESHES

      M=>MESHES(NM)
      CALL POINT_TO_MESH(NM)
 
      ! Evacuation criteria
 
      IF (MESH_ID/='null' .AND. MESH_ID/=MESH_NAME(NM)) CYCLE MESH_LOOP
      IF (EVACUATION .AND. .NOT.EVACUATION_ONLY(NM)) CYCLE MESH_LOOP
      IF (EVACUATION_ONLY(NM)) THEN 
         IF (.NOT.TEMP_HOLE_EVAC(N)) CYCLE MESH_LOOP
      ENDIF
 
      ! Loop over all possible multiples of the HOLE

      MR => MULTIPLIER(0)
      DO NNN=1,N_MULT
         IF (MULT_ID==MULTIPLIER(NNN)%ID) MR => MULTIPLIER(NNN)
      ENDDO

      K_MULT_LOOP: DO KK=MR%K_LOWER,MR%K_UPPER
         J_MULT_LOOP: DO JJ=MR%J_LOWER,MR%J_UPPER
            I_MULT_LOOP: DO II=MR%I_LOWER,MR%I_UPPER

               IF (.NOT.MR%SEQUENTIAL) THEN
                  X1 = XB(1) + MR%DX0 + II*MR%DXB(1)
                  X2 = XB(2) + MR%DX0 + II*MR%DXB(2)
                  Y1 = XB(3) + MR%DY0 + JJ*MR%DXB(3)
                  Y2 = XB(4) + MR%DY0 + JJ*MR%DXB(4)
                  Z1 = XB(5) + MR%DZ0 + KK*MR%DXB(5)
                  Z2 = XB(6) + MR%DZ0 + KK*MR%DXB(6)
               ELSE
                  X1 = XB(1) + MR%DX0 + II*MR%DXB(1)
                  X2 = XB(2) + MR%DX0 + II*MR%DXB(2)
                  Y1 = XB(3) + MR%DY0 + II*MR%DXB(3)
                  Y2 = XB(4) + MR%DY0 + II*MR%DXB(4)
                  Z1 = XB(5) + MR%DZ0 + II*MR%DXB(5)
                  Z2 = XB(6) + MR%DZ0 + II*MR%DXB(6)
               ENDIF

               ! Check if hole is contained within the current mesh
 
               IF (X1>=XF .OR. X2<=XS .OR. Y1>YF .OR. Y2<=YS .OR. Z1>ZF .OR. Z2<=ZS) CYCLE I_MULT_LOOP
 
               ! Assign mesh-limited bounds

               X1 = MAX(X1,XS-0.001_EB*DX(0))
               X2 = MIN(X2,XF+0.001_EB*DX(IBP1))
               Y1 = MAX(Y1,YS-0.001_EB*DY(0))
               Y2 = MIN(Y2,YF+0.001_EB*DY(JBP1))
               Z1 = MAX(Z1,ZS-0.001_EB*DZ(0))
               Z2 = MIN(Z2,ZF+0.001_EB*DZ(KBP1))
 
               I1 = NINT( GINV(X1-XS,1,NM)*RDXI   ) 
               I2 = NINT( GINV(X2-XS,1,NM)*RDXI   )
               J1 = NINT( GINV(Y1-YS,2,NM)*RDETA  ) 
               J2 = NINT( GINV(Y2-YS,2,NM)*RDETA  )
               K1 = NINT( GINV(Z1-ZS,3,NM)*RDZETA ) 
               K2 = NINT( GINV(Z2-ZS,3,NM)*RDZETA )
 
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
         
                  DEVC_OR_CTRL: IF ((DEVC_ID/='null' .OR. CTRL_ID/='null') .AND. .NOT.EVACUATION_ONLY(NM)) THEN
         
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
            ENDDO I_MULT_LOOP
         ENDDO J_MULT_LOOP
      ENDDO K_MULT_LOOP
   ENDDO MESH_LOOP
ENDDO READ_HOLE_LOOP
 
REWIND(LU_INPUT)

IF(ANY(EVACUATION_ONLY)) DEALLOCATE(TEMP_HOLE_EVAC)
DEALLOCATE(TEMP_OBST)

CONTAINS

  SUBROUTINE DEFINE_EVACUATION_HOLES(IMODE)
    !
    ! Clear the STRS meshes by a hole with size of XB of the stairs.
    ! The core is put there applying permit_hole=false.
    USE EVAC, ONLY: EVAC_FDS6, N_STRS, EMESH_STAIRS, EVAC_EMESH_STAIRS_TYPE
    IMPLICIT NONE
    ! Passed variables
    INTEGER, INTENT(IN) :: IMODE
    ! Local variables
    INTEGER :: I

    IF (.NOT.ANY(EVACUATION_ONLY)) RETURN
    IMODE_1_IF: IF (IMODE==1) THEN
       NSTRS_LOOP_1: DO I = 1, N_STRS
          IF (.NOT.EMESH_STAIRS(I)%DEFINE_MESH) CYCLE NSTRS_LOOP_1
          N_HOLE_O = N_HOLE_O + 1
          N_HOLE   = N_HOLE   + 1 ! No mult for evacuation strs meshes
          EMESH_STAIRS(I)%I_HOLE = N_HOLE_O
          EVACUATION_HOLE = .TRUE.
       END DO NSTRS_LOOP_1
    END IF IMODE_1_IF

    IMODE_2_IF: IF (IMODE==2) THEN
       EVACUATION_HOLE = .FALSE.
       NSTRS_LOOP_2: DO I = 1, N_STRS
          IF (.NOT.EMESH_STAIRS(I)%DEFINE_MESH) CYCLE NSTRS_LOOP_2
          IF (.NOT.EMESH_STAIRS(I)%I_HOLE==N) CYCLE NSTRS_LOOP_2
          EVACUATION_HOLE = .TRUE.
          EXIT NSTRS_LOOP_2
       END DO NSTRS_LOOP_2
    END IF IMODE_2_IF

    IMODE_3_IF: IF (IMODE==3) THEN
       EVACUATION_HOLE = .FALSE.
       NSTRS_LOOP_3: DO I = 1, N_STRS
          IF (.NOT.EMESH_STAIRS(I)%DEFINE_MESH) CYCLE NSTRS_LOOP_3
          IF (.NOT.EMESH_STAIRS(I)%I_HOLE==N) CYCLE NSTRS_LOOP_3
          EVACUATION_HOLE = .TRUE.
          RGB = EMESH_STAIRS(I)%RGB
          XB = EMESH_STAIRS(I)%XB
          EVACUATION = .TRUE.
          MESH_ID = TRIM(MESH_NAME(EMESH_STAIRS(I)%IMESH))
          EXIT NSTRS_LOOP_3
       END DO NSTRS_LOOP_3
    END IF IMODE_3_IF

  END SUBROUTINE DEFINE_EVACUATION_HOLES

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

INTEGER :: N,NN,NM,NNN,NVO,IOR,I1,I2,J1,J2,K1,K2,RGB(3),N_EDDY
REAL(EB) :: SPREAD_RATE,TRANSPARENCY,XYZ(3),TMP_EXTERIOR,MASS_FRACTION(MAX_SPECIES),DYNAMIC_PRESSURE, &
            MEAN_VELOCITY,RELAXATION_FACTOR,REYNOLDS_STRESS(3,3)=0._EB,L_EDDY,VEL_RMS
CHARACTER(30) :: ID,DEVC_ID,CTRL_ID,SURF_ID,PRESSURE_RAMP
CHARACTER(60) :: MESH_ID
CHARACTER(25) :: COLOR
LOGICAL :: REJECT_VENT,EVACUATION,OUTLINE,EVACUATION_VENT
NAMELIST /VENT/ XB,ID,IOR,MB,PBX,PBY,PBZ,SURF_ID,FYI,RGB,TRANSPARENCY,COLOR, &
                TEXTURE_ORIGIN,OUTLINE,DEVC_ID,CTRL_ID, &
                XYZ,EVACUATION,MESH_ID,SPREAD_RATE,TMP_EXTERIOR,MASS_FRACTION,DYNAMIC_PRESSURE,PRESSURE_RAMP, &
                MEAN_VELOCITY,RELAXATION_FACTOR,REYNOLDS_STRESS,L_EDDY,N_EDDY,VEL_RMS
 
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

   IF (EVACUATION_ONLY(NM)) CALL DEFINE_EVACUATION_VENTS(NM,1)
 
   IF (TWO_D)                         N_VENT = N_VENT + 2
   IF (CYLINDRICAL .AND. M%XS<=ZERO_P) N_VENT = N_VENT + 1
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
      ID      = 'null'
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
      MEAN_VELOCITY = -999._EB
      RELAXATION_FACTOR = 0.5_EB
      N_EDDY=0
      L_EDDY=0._EB
      VEL_RMS=0._EB
 
      IF (NN==NVO-2 .AND. CYLINDRICAL .AND. XS<=ZERO_P) MB='XMIN'
      IF (NN==NVO-1 .AND. TWO_D)                        MB='YMIN'
      IF (NN==NVO   .AND. TWO_D)                        MB='YMAX'
      IF (NN==NVO-1 .AND. EVACUATION_ONLY(NM))          MB='ZMIN'
      IF (NN==NVO   .AND. EVACUATION_ONLY(NM))          MB='ZMAX'
 
      IF (MB=='null') THEN
         EVACUATION_VENT = .FALSE.
         IF (EVACUATION_ONLY(NM)) CALL DEFINE_EVACUATION_VENTS(NM,2)
         EVACUATION_VENTS: IF (.NOT. EVACUATION_VENT) THEN
            CALL CHECKREAD('VENT',LU_INPUT,IOS) 
            IF (IOS==1) EXIT READ_VENT_LOOP
            READ(LU_INPUT,VENT,END=37,ERR=38)    ! Read in info for VENT N
         END IF EVACUATION_VENTS
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
 
      IF (ABS(XB(3)-XB(4))<=ZERO_P  .AND. TWO_D .AND. NN<NVO-1) THEN
         IF (ID=='null')WRITE(MESSAGE,'(A,I4,A)')'ERROR: VENT ',NN,      ' cannot be specified on a y boundary in a 2D calculation'
         IF (ID/='null')WRITE(MESSAGE,'(A,A,A)') 'ERROR: VENT ',TRIM(ID),' cannot be specified on a y boundary in a 2D calculation'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
 
      IF (ABS(XB(1)-XB(2))>ZERO_P  .AND. ABS(XB(3)-XB(4))>ZERO_P  .AND.ABS(XB(5)-XB(6))>ZERO_P ) THEN
         IF (ID=='null') WRITE(MESSAGE,'(A,I4,A)') 'ERROR: VENT ',NN,      ' must be a plane'
         IF (ID/='null') WRITE(MESSAGE,'(A,A,A)')  'ERROR: VENT ',TRIM(ID),' must be a plane'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
 
      CALL CHECK_XB(XB)

      IF (ALL(EVACUATION_ONLY)) THEN
         DEVC_ID    = 'null'
         CTRL_ID    = 'null'
      END IF
      
      VT=>VENTS(N)
      
      IF (ABS(XB(1)-XB(2))<=ZERO_P ) VT%TOTAL_INPUT_AREA = (XB(4)-XB(3))*(XB(6)-XB(5))
      IF (ABS(XB(3)-XB(4))<=ZERO_P ) VT%TOTAL_INPUT_AREA = (XB(2)-XB(1))*(XB(6)-XB(5))
      IF (ABS(XB(5)-XB(6))<=ZERO_P ) VT%TOTAL_INPUT_AREA = (XB(2)-XB(1))*(XB(4)-XB(3))

      XB(1) = MAX(XB(1),XS-DX(0))
      XB(2) = MIN(XB(2),XF+DX(IBP1))
      XB(3) = MAX(XB(3),YS-DY(0))
      XB(4) = MIN(XB(4),YF+DY(JBP1))
      XB(5) = MAX(XB(5),ZS-DZ(0))
      XB(6) = MIN(XB(6),ZF+DZ(KBP1))
 
      IF (XB(1)>XF+DX(IBP1) .OR. XB(2)<XS-DX(0) .OR. &
          XB(3)>YF+DY(JBP1) .OR. XB(4)<YS-DY(0) .OR. &
          XB(5)>ZF+DZ(KBP1) .OR. XB(6)<ZS-DZ(0)) REJECT_VENT = .TRUE.
 
      VT%I1 = NINT( GINV(XB(1)-XS,1,NM)*RDXI   ) 
      VT%I2 = NINT( GINV(XB(2)-XS,1,NM)*RDXI   )
      VT%J1 = NINT( GINV(XB(3)-YS,2,NM)*RDETA  ) 
      VT%J2 = NINT( GINV(XB(4)-YS,2,NM)*RDETA  )
      VT%K1 = NINT( GINV(XB(5)-ZS,3,NM)*RDZETA )
      VT%K2 = NINT( GINV(XB(6)-ZS,3,NM)*RDZETA )
 
      ! Thicken evacuation mesh vents in the z direction

      IF (EVACUATION_ONLY(NM) .AND. EVACUATION .AND. VT%K1==VT%K2 .AND. .NOT.REJECT_VENT) THEN
         VT%K1 = GINV(.5_EB*(XB(5)+XB(6))-ZS,3,NM)*RDZETA
         VT%K2 = KBAR
         XB(5) = ZS
         XB(6) = ZF
         IF (ABS(XB(1)-XB(2))>ZERO_P  .AND. ABS(XB(3)-XB(4))>ZERO_P ) THEN
            IF (ID=='null') WRITE(MESSAGE,'(A,I4,A)') 'ERROR: Evacuation VENT ',NN,      ' must be a vertical plane'
            IF (ID/='null') WRITE(MESSAGE,'(A,A,A)')  'ERROR: Evacuation VENT ',TRIM(ID),' must be a vertical plane'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDIF

      IF (ABS(XB(1)-XB(2))<=ZERO_P ) THEN
         IF (VT%J1==VT%J2 .OR. VT%K1==VT%K2) REJECT_VENT=.TRUE.
         IF (VT%I1>IBAR .OR. VT%I2<0)        REJECT_VENT=.TRUE.
      ENDIF
      IF (ABS(XB(3)-XB(4))<=ZERO_P ) THEN
         IF (VT%I1==VT%I2 .OR. VT%K1==VT%K2) REJECT_VENT=.TRUE.
         IF (VT%J1>JBAR .OR. VT%J2<0)        REJECT_VENT=.TRUE.
      ENDIF
      IF (ABS(XB(5)-XB(6))<=ZERO_P ) THEN
         IF (VT%I1==VT%I2 .OR. VT%J1==VT%J2) REJECT_VENT=.TRUE.
         IF (VT%K1>KBAR .OR. VT%K2<0)        REJECT_VENT=.TRUE.
      ENDIF

      ! Evacuation criteria
 
      IF (.NOT.EVACUATION .AND. EVACUATION_ONLY(NM)) REJECT_VENT=.TRUE.
      IF (EVACUATION .AND. .NOT.EVACUATION_ONLY(NM)) REJECT_VENT=.TRUE.
 
      IF (ALL(EVACUATION_ONLY)) THEN
         DEVC_ID    = 'null'
         CTRL_ID    = 'null'
      END IF

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
 
      IF (ABS(XB(1)-XB(2))<=ZERO_P ) VT%INPUT_AREA = (XB(4)-XB(3))*(XB(6)-XB(5))
      IF (ABS(XB(3)-XB(4))<=ZERO_P ) VT%INPUT_AREA = (XB(2)-XB(1))*(XB(6)-XB(5))
      IF (ABS(XB(5)-XB(6))<=ZERO_P ) VT%INPUT_AREA = (XB(2)-XB(1))*(XB(4)-XB(3))
 
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
      IF (VT%IBC==HVAC_SURF_INDEX)     VT%BOUNDARY_TYPE = HVAC_BOUNDARY
      VT%IOR = IOR
      VT%ORDINAL = NN
 
      ! Activate and Deactivate logic

      VT%ACTIVATED = .TRUE.
      VT%DEVC_ID   = DEVC_ID
      VT%CTRL_ID   = CTRL_ID
      VT%ID        = ID      
      CALL SEARCH_CONTROLLER('VENT',CTRL_ID,DEVC_ID,VT%DEVC_INDEX,VT%CTRL_INDEX,N)
      IF (DEVC_ID /= 'null') THEN
         IF (.NOT.DEVICE(VT%DEVC_INDEX)%INITIAL_STATE) VT%ACTIVATED = .FALSE.
      ENDIF
      IF (CTRL_ID /= 'null') THEN
         IF (.NOT.CONTROL(VT%CTRL_INDEX)%INITIAL_STATE) VT%ACTIVATED = .FALSE.
      ENDIF

      IF ( (VT%BOUNDARY_TYPE==OPEN_BOUNDARY .OR. VT%BOUNDARY_TYPE==MIRROR_BOUNDARY .OR. VT%BOUNDARY_TYPE==PERIODIC_BOUNDARY) .AND. &
           (VT%DEVC_ID /= 'null' .OR. VT%CTRL_ID /= 'null') ) THEN
         IF (ID=='null') WRITE(MESSAGE,'(A,I4,A)') 'ERROR: VENT ',NN,      ' cannot be controlled by a device'
         IF (ID/='null') WRITE(MESSAGE,'(A,A,A)')  'ERROR: VENT ',TRIM(ID),' cannot be controlled by a device'
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
      
      ! Boundary Forcing
      
      VT%FVEL = MEAN_VELOCITY
      VT%RFAC = RELAXATION_FACTOR
      IF (VT%FVEL>-100._EB) VT%BOUNDARY_FORCING=.TRUE.
      
      ! Synthetic Eddy Method
      
      VT%N_EDDY = N_EDDY
      VT%L_EDDY = L_EDDY
      IF (VEL_RMS>0._EB) THEN
         VT%R_IJ=0._EB
         VT%R_IJ(1,1)=VEL_RMS**2
         VT%R_IJ(2,2)=VEL_RMS**2
         VT%R_IJ(3,3)=VEL_RMS**2
      ELSE
         VT%R_IJ = REYNOLDS_STRESS
      ENDIF
      
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
 
   VENT_LOOP_2: DO N=1,N_VENT
 
      VT => VENTS(N)
 
      I1 = MAX(0,VT%I1)
      I2 = MIN(IBAR,VT%I2)
      J1 = MAX(0,VT%J1)
      J2 = MIN(JBAR,VT%J2)
      K1 = MAX(0,VT%K1)
      K2 = MIN(KBAR,VT%K2)
 
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
 
   ENDDO VENT_LOOP_2
 
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
 
      VT%GHOST_CELLS_ONLY = .TRUE.
 
      SELECT CASE(ABS(VT%IOR))
         CASE(1)
            DO K=K1+1,K2
               DO J=J1+1,J2
                  IF (J>=1 .AND. J<=JBAR .AND. K>=1 .AND. K<=KBAR) VT%GHOST_CELLS_ONLY = .FALSE.
                  VT%FDS_AREA = VT%FDS_AREA + DY(J)*DZ(K)
               ENDDO
            ENDDO
         CASE(2)
            DO K=K1+1,K2
               DO I=I1+1,I2
                  IF (I>=1 .AND. I<=IBAR .AND. K>=1 .AND. K<=KBAR) VT%GHOST_CELLS_ONLY = .FALSE.
                  VT%FDS_AREA = VT%FDS_AREA + DX(I)*DZ(K)
               ENDDO
            ENDDO
         CASE(3)
            DO J=J1+1,J2
               DO I=I1+1,I2
                  IF (I>=1 .AND. I<=IBAR .AND. J>=1 .AND. J<=JBAR) VT%GHOST_CELLS_ONLY = .FALSE.
                  VT%FDS_AREA = VT%FDS_AREA + DX(I)*DY(J)
               ENDDO
            ENDDO
      END SELECT
 
   ENDDO  VENT_LOOP_3
   
   ! Allocate arrays for turbulent inflow boundary conditions
 
   VENT_LOOP_4: DO N=1,N_VENT
      VT => VENTS(N)    
      EDDY_IF: IF (VT%N_EDDY>0) THEN
         SELECT CASE(ABS(VT%IOR))
            CASE(1)
               ALLOCATE(VT%U_EDDY(VT%J1+1:VT%J2,VT%K1+1:VT%K2),STAT=IZERO)
               CALL ChkMemErr('READ_VENT','U_PRIME',IZERO)
               ALLOCATE(VT%V_EDDY(VT%J1+1:VT%J2,VT%K1+1:VT%K2),STAT=IZERO)
               CALL ChkMemErr('READ_VENT','V_PRIME',IZERO)
               ALLOCATE(VT%W_EDDY(VT%J1+1:VT%J2,VT%K1+1:VT%K2),STAT=IZERO)
               CALL ChkMemErr('READ_VENT','W_PRIME',IZERO)
            CASE(2)
               ALLOCATE(VT%U_EDDY(VT%I1+1:VT%I2,VT%K1+1:VT%K2),STAT=IZERO)
               CALL ChkMemErr('READ_VENT','U_PRIME',IZERO)
               ALLOCATE(VT%V_EDDY(VT%I1+1:VT%I2,VT%K1+1:VT%K2),STAT=IZERO)
               CALL ChkMemErr('READ_VENT','V_PRIME',IZERO)
               ALLOCATE(VT%W_EDDY(VT%I1+1:VT%I2,VT%K1+1:VT%K2),STAT=IZERO)
               CALL ChkMemErr('READ_VENT','W_PRIME',IZERO)
            CASE(3)
               ALLOCATE(VT%U_EDDY(VT%I1+1:VT%I2,VT%J1+1:VT%J2),STAT=IZERO)
               CALL ChkMemErr('READ_VENT','U_PRIME',IZERO)
               ALLOCATE(VT%V_EDDY(VT%I1+1:VT%I2,VT%J1+1:VT%J2),STAT=IZERO)
               CALL ChkMemErr('READ_VENT','V_PRIME',IZERO)
               ALLOCATE(VT%W_EDDY(VT%I1+1:VT%I2,VT%J1+1:VT%J2),STAT=IZERO)
               CALL ChkMemErr('READ_VENT','W_PRIME',IZERO)
         END SELECT
         ALLOCATE(VT%X_EDDY(VT%N_EDDY),STAT=IZERO)
         CALL ChkMemErr('READ_VENT','X_EDDY',IZERO)
         ALLOCATE(VT%Y_EDDY(VT%N_EDDY),STAT=IZERO)
         CALL ChkMemErr('READ_VENT','Y_EDDY',IZERO)
         ALLOCATE(VT%Z_EDDY(VT%N_EDDY),STAT=IZERO)
         CALL ChkMemErr('READ_VENT','Z_EDDY',IZERO)
         ALLOCATE(VT%CU_EDDY(VT%N_EDDY),STAT=IZERO)
         CALL ChkMemErr('READ_VENT','CU_EDDY',IZERO)
         ALLOCATE(VT%CV_EDDY(VT%N_EDDY),STAT=IZERO)
         CALL ChkMemErr('READ_VENT','CV_EDDY',IZERO)
         ALLOCATE(VT%CW_EDDY(VT%N_EDDY),STAT=IZERO)
         CALL ChkMemErr('READ_VENT','CW_EDDY',IZERO)
      ENDIF EDDY_IF
   ENDDO VENT_LOOP_4
 
ENDDO MESH_LOOP_2
 

CONTAINS

  SUBROUTINE DEFINE_EVACUATION_VENTS(NM,IMODE)
    !
    ! Define the evacuation outflow VENTs for the doors/exits.
    !
    USE EVAC, ONLY: N_DOORS, N_EXITS, N_CO_EXITS, EVAC_EMESH_EXITS_TYPE, EMESH_EXITS, N_DOOR_MESHES, EVAC_FDS6
    IMPLICIT NONE
    ! Passed variables
    INTEGER, INTENT(IN) :: NM, IMODE
    ! Local variables
    INTEGER :: N, N_END

    N_END = N_EXITS - N_CO_EXITS + N_DOORS
    IMODE_1_IF: IF (IMODE==1) THEN
       NEND_LOOP_1: DO N = 1, N_END
          IF (.NOT.EMESH_EXITS(N)%DEFINE_MESH) CYCLE NEND_LOOP_1
          IF (EMESH_EXITS(N)%IMESH==NM .OR. EMESH_EXITS(N)%MAINMESH==NM) THEN
             N_VENT = N_VENT + 1
             EMESH_EXITS(N)%I_VENT = N_VENT
             EVACUATION_VENT = .TRUE.
             EVACUATION = .TRUE.
             IF (.NOT.EVAC_FDS6 .AND. EMESH_EXITS(N)%IMESH==NM) EXIT NEND_LOOP_1
          END IF
       END DO NEND_LOOP_1
    END IF IMODE_1_IF

    IMODE_2_IF: IF (IMODE==2) THEN
       ! Evacuation VENTs (for the outflow vents) need: XB, EVACUATION, RGB, MESH_ID, SURF_ID, IOR
       NEND_LOOP_2: DO N = 1, N_END
          IF (.NOT.EMESH_EXITS(N)%DEFINE_MESH) CYCLE NEND_LOOP_2
          IF (EMESH_EXITS(N)%I_VENT==NN .AND. (EMESH_EXITS(N)%IMESH==NM .OR. EMESH_EXITS(N)%MAINMESH==NM)) THEN
             EVACUATION_VENT = .TRUE.
             EVACUATION = .TRUE.
             SURF_ID = 'EVACUATION_OUTFLOW'
             MESH_ID = TRIM(MESH_NAME(NM))
             XB(1) = EMESH_EXITS(N)%XB(1)
             XB(2) = EMESH_EXITS(N)%XB(2)
             XB(3) = EMESH_EXITS(N)%XB(3)
             XB(4) = EMESH_EXITS(N)%XB(4)
             XB(5) = EMESH_EXITS(N)%XB(5)
             XB(6) = EMESH_EXITS(N)%XB(6)
             RGB(:) = EMESH_EXITS(N)%RGB(:)
             ID = TRIM('Event_' // TRIM(MESH_NAME(NM)))
             IF (.NOT.EVAC_FDS6 .AND. EMESH_EXITS(N)%IMESH==NM) EXIT NEND_LOOP_2   ! One VENT per door flow mesh
          END IF
       END DO NEND_LOOP_2
    END IF IMODE_2_IF

    RETURN
  END SUBROUTINE DEFINE_EVACUATION_VENTS

END SUBROUTINE READ_VENT
 
 

SUBROUTINE READ_INIT

USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT
USE COMP_FUNCTIONS, ONLY: GET_FILE_NUMBER
REAL(EB) :: TEMPERATURE,DENSITY,MASS_FRACTION(1:MAX_SPECIES),RR_SUM,YY_GET(1:N_GAS_SPECIES),MASS_PER_VOLUME, &
            MASS_PER_TIME,DT_INSERT,XYZ(3),UVW(3),DIAMETER
INTEGER  :: N,NN,NNN,II,JJ,KK,NS,NS2,NUMBER_INITIAL_DROPLETS,N_INIT_NEW,N_INIT_READ
CHARACTER(30) :: PART_ID,SHAPE,MULT_ID,ORIGIN='CENTER',PROF_ID,SPEC_ID(1:MAX_SPECIES)
TYPE(INITIALIZATION_TYPE), POINTER :: IN
TYPE(MULTIPLIER_TYPE), POINTER :: MR
NAMELIST /INIT/ XB,TEMPERATURE,DENSITY,MASS_FRACTION,MASS_PER_TIME,MASS_PER_VOLUME,NUMBER_INITIAL_DROPLETS, &
                PART_ID,SHAPE,DT_INSERT,MULT_ID,PROF_ID,ORIGIN,XYZ,UVW,DIAMETER,SPEC_ID
 
N_INIT = 0
N_INIT_READ = 0
REWIND(LU_INPUT)

COUNT_LOOP: DO
   CALL CHECKREAD('INIT',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_LOOP
   MULT_ID = 'null'
   READ(LU_INPUT,NML=INIT,END=11,ERR=12,IOSTAT=IOS)
   N_INIT_READ = N_INIT_READ + 1
   12 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with INIT number ',N_INIT_READ+1
      CALL SHUTDOWN(MESSAGE)
      ENDIF
   N_INIT_NEW = 1
   IF (MULT_ID/='null') THEN
      DO N=1,N_MULT
         MR => MULTIPLIER(N)
         IF (MULT_ID==MR%ID) N_INIT_NEW = MR%N_COPIES
      ENDDO
   ENDIF
   N_INIT = N_INIT + N_INIT_NEW
ENDDO COUNT_LOOP
11 REWIND(LU_INPUT)
 
! If there are no INIT lines, return

IF (N_INIT==0) RETURN
 
ALLOCATE(INITIALIZATION(N_INIT),STAT=IZERO)
CALL ChkMemErr('READ','INITIALIZATION',IZERO)
 
NN = 0

INIT_LOOP: DO N=1,N_INIT_READ

   DENSITY                    = -1000._EB
   DIAMETER                   =  0._EB
   DT_INSERT                  = -1._EB
   MASS_FRACTION              =  0._EB
   MASS_PER_TIME              = -1._EB
   MASS_PER_VOLUME            = -1._EB
   MULT_ID                    = 'null'
   NUMBER_INITIAL_DROPLETS    =  0
   PART_ID                    = 'null'
   PROF_ID                    = 'null'
   SHAPE                      = 'BLOCK'
   SPEC_ID = 'null'
   TEMPERATURE                = -1000._EB
   UVW                        =  0._EB
   XYZ                        = -1.E6_EB
   XB(1)                      = -1000000._EB
   XB(2)                      =  1000000._EB
   XB(3)                      = -1000000._EB
   XB(4)                      =  1000000._EB
   XB(5)                      = -1000000._EB
   XB(6)                      =  1000000._EB
 
   CALL CHECKREAD('INIT',LU_INPUT,IOS)
   IF (IOS==1) EXIT INIT_LOOP
   READ(LU_INPUT,INIT) 
 
   ! Reorder XB coordinates if necessary

   CALL CHECK_XB(XB)

   ! Loop over all possible multiples of the INIT

   MR => MULTIPLIER(0)
   DO NNN=1,N_MULT
      IF (MULT_ID==MULTIPLIER(NNN)%ID) MR => MULTIPLIER(NNN)
   ENDDO

   K_MULT_LOOP: DO KK=MR%K_LOWER,MR%K_UPPER
      J_MULT_LOOP: DO JJ=MR%J_LOWER,MR%J_UPPER
         I_MULT_LOOP: DO II=MR%I_LOWER,MR%I_UPPER

            NN = NN + 1
            IN => INITIALIZATION(NN)

            ! Store the input parameters

            IF (.NOT.MR%SEQUENTIAL) THEN
               IN%X1 = XB(1) + MR%DX0 + II*MR%DXB(1)
               IN%X2 = XB(2) + MR%DX0 + II*MR%DXB(2)
               IN%Y1 = XB(3) + MR%DY0 + JJ*MR%DXB(3)
               IN%Y2 = XB(4) + MR%DY0 + JJ*MR%DXB(4)
               IN%Z1 = XB(5) + MR%DZ0 + KK*MR%DXB(5)
               IN%Z2 = XB(6) + MR%DZ0 + KK*MR%DXB(6)
            ELSE
               IN%X1 = XB(1) + MR%DX0 + II*MR%DXB(1)
               IN%X2 = XB(2) + MR%DX0 + II*MR%DXB(2)
               IN%Y1 = XB(3) + MR%DY0 + II*MR%DXB(3)
               IN%Y2 = XB(4) + MR%DY0 + II*MR%DXB(4)
               IN%Z1 = XB(5) + MR%DZ0 + II*MR%DXB(5)
               IN%Z2 = XB(6) + MR%DZ0 + II*MR%DXB(6)
            ENDIF

            IN%VOLUME        = (IN%X2-IN%X1)*(IN%Y2-IN%Y1)*(IN%Z2-IN%Z1)
            IN%TEMPERATURE   = TEMPERATURE + TMPM
            IN%DENSITY       = DENSITY
            IN%SHAPE         = SHAPE
         
            IF (DENSITY     > 0._EB) RHOMAX = MAX(RHOMAX,IN%DENSITY)

            SPEC_INIT_IF:IF (N_GAS_SPECIES > 0._EB .AND. ANY(MASS_FRACTION > ZERO_P)) THEN
               IF (SPEC_ID(1)=='null') THEN
                  WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR: Problem with INIT number ',N,'. SPEC_ID must be used with MASS_FRACTION'
                  CALL SHUTDOWN(MESSAGE)
               ELSE 
                  DO NS=1,N_GAS_SPECIES
                     IF (SPEC_ID(NS)=='null') EXIT
                     DO NS2=1,N_SPECIES
                        IF (TRIM(SPEC_ID(NS))==TRIM(SPECIES(NS2)%ID)) THEN
                           IF(SPECIES(NS2)%Y_INDEX<=0) THEN
                              WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR: Problem with INIT number ',N,' SPEC ',&
                              TRIM(SPECIES(NS2)%ID),' is not a tracked species'
                              CALL SHUTDOWN(MESSAGE)
                           ENDIF
                           IN%MASS_FRACTION(SPECIES(NS2)%Y_INDEX)=MASS_FRACTION(NS)
                           EXIT
                        ENDIF
                        IF (NS2==N_SPECIES)  THEN
                           WRITE(MESSAGE,'(A,I3,A,A,A)') 'ERROR: Problem with INIT number ',N,' SPEC ',&
                              TRIM(SPEC_ID(NS)),' not found'
                              CALL SHUTDOWN(MESSAGE)
                        ENDIF
                     ENDDO   
                  ENDDO  

                  IF (SUM(IN%MASS_FRACTION) > 1._EB) THEN
                     WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR: Problem with INIT number ',N,'. Sum of specified mass fractions > 1'
                     CALL SHUTDOWN(MESSAGE)
                  ENDIF

                  YY_GET = IN%MASS_FRACTION(1:N_GAS_SPECIES)
                  CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RR_SUM)
               ENDIF            
            ELSE SPEC_INIT_IF
                IF (N_GAS_SPECIES > 0) IN%MASS_FRACTION(1:N_GAS_SPECIES) = YY_AMBIENT(1:N_GAS_SPECIES)            
                RR_SUM = RSUM0
            ENDIF SPEC_INIT_IF
            
            IF (TEMPERATURE > 0._EB) TMPMIN = MIN(TMPMIN,IN%TEMPERATURE)
         
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
         
            ! Special case where INIT is used to introduce a block of particles
         
            IN%MASS_PER_TIME           = MASS_PER_TIME
            IN%MASS_PER_VOLUME         = MASS_PER_VOLUME
            IN%NUMBER_INITIAL_DROPLETS = NUMBER_INITIAL_DROPLETS

            IF ( (IN%MASS_PER_TIME>0._EB .OR. IN%MASS_PER_VOLUME>0._EB) .AND. IN%VOLUME<=ZERO_P) THEN
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: INIT ',N,' XB has no volume'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            
            IF (DT_INSERT<0._EB) THEN
               IN%DT_INSERT = T_END-T_BEGIN+10.
            ELSE
               IN%DT_INSERT = DT_INSERT
            ENDIF

            ! Set up a clock to keep track of particle insertions

            ALLOCATE(IN%PARTICLE_INSERT_CLOCK(NMESHES),STAT=IZERO)
            CALL ChkMemErr('READ','PARTICLE_INSERT_CLOCK',IZERO)
            IN%PARTICLE_INSERT_CLOCK = T_BEGIN
         
            ! Assign an index to identify the particle class

            IF (PART_ID/='null') THEN
               DO NS=1,N_PART
                  IF (PART_ID==PARTICLE_CLASS(NS)%ID) THEN
                     IN%PART_INDEX = NS
                     DROPLET_FILE = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
               IF (IN%PART_INDEX<1) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: PART_ID ',TRIM(PART_ID),' does not exist'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
               IN%NUMBER_INITIAL_DROPLETS = NUMBER_INITIAL_DROPLETS*PARTICLE_CLASS(IN%PART_INDEX)%N_SPLIT
            ENDIF
            
            ! Special case: INIT is used to patch a velocity profile

            PATCH_VELOCITY_IF: IF (PATCH_VELOCITY) THEN
               SELECT CASE(ORIGIN)
                  CASE('CENTER')
                     IN%X0 = 0.5_EB*(IN%X1+IN%X2)
                     IN%Y0 = 0.5_EB*(IN%Y1+IN%Y2)
                     IN%Z0 = 0.5_EB*(IN%Z1+IN%Z2)
                  CASE('CORNER')
                     IN%X0 = IN%X1
                     IN%Y0 = IN%Y1
                     IN%Z0 = IN%Z1
                  CASE DEFAULT
                     WRITE(MESSAGE,'(A)') 'ERROR: ORIGIN specified incorrectly on INIT'
                     CALL SHUTDOWN(MESSAGE)
               END SELECT
               IN%PROF_ID = PROF_ID
            ENDIF PATCH_VELOCITY_IF
            
            ! Special case: POINTWISE_DROPLET_INIT

            POINTWISE_DROPLET_IF: IF (ANY(XYZ>-1000._EB)) THEN
               IF (NUMBER_INITIAL_DROPLETS>1) THEN
                  WRITE(MESSAGE,'(A)') 'ERROR: NUMBER_INITIAL_DROPLETS=1 for pointwise insert'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
               IN%SHAPE = 'POINT'
               IN%NUMBER_INITIAL_DROPLETS=1     ! one droplet per INIT line
               IN%POINTWISE_DROPLET_INIT=.TRUE. ! flag for insert loop
               IN%X0 = XYZ(1) + MR%DX0
               IN%Y0 = XYZ(2) + MR%DY0
               IN%Z0 = XYZ(3) + MR%DZ0
               IN%U0 = UVW(1)
               IN%V0 = UVW(2)
               IN%W0 = UVW(3)
               IN%RADIUS = DIAMETER*0.5E-6_EB        ! convert diameter (microns) to radius (meters)
            ENDIF POINTWISE_DROPLET_IF
         
         ENDDO I_MULT_LOOP
      ENDDO J_MULT_LOOP
   ENDDO K_MULT_LOOP

ENDDO INIT_LOOP

REWIND(LU_INPUT)

END SUBROUTINE READ_INIT


SUBROUTINE READ_ZONE
 
INTEGER, PARAMETER :: MAX_LEAK_PATHS=200
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
      WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with ZONE number ',N_ZONE+1
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

   ALLOCATE(P_ZONE(N)%LEAK_AREA(0:N_ZONE),STAT=IZERO)
   CALL ChkMemErr('READ','LEAK_AREA',IZERO)
 
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
   P_ZONE(N)%LEAK_AREA(0:N_ZONE) = LEAK_AREA(0:N_ZONE)
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

USE DEVICE_VARIABLES, ONLY: DEVICE_TYPE, DEVICE, N_DEVC, N_DEVC_TIME, N_DEVC_LINE,MAX_DEVC_LINE_POINTS
INTEGER  :: NN,NM,MESH_NUMBER,N_DEVC_READ,IOR,TRIP_DIRECTION,VELO_INDEX,POINTS,I_POINT
REAL(EB) :: DEPTH,ORIENTATION(3),ROTATION,SETPOINT,FLOWRATE,BYPASS_FLOWRATE,DELAY,XYZ(3),CONVERSION_FACTOR
CHARACTER(30) :: QUANTITY,PROP_ID,CTRL_ID,DEVC_ID,SURF_ID,STATISTICS,PART_ID,MATL_ID,SPEC_ID,UNITS,DUCT_ID,NODE_ID,X_ID,Y_ID,Z_ID
LOGICAL :: INITIAL_STATE,LATCH,DRY,TIME_AVERAGED,EVACUATION,HIDE_COORDINATES
TYPE (DEVICE_TYPE), POINTER :: DV
NAMELIST /DEVC/ DEPTH,FYI,IOR,ID,ORIENTATION,PROP_ID,QUANTITY,ROTATION,XB,XYZ,INITIAL_STATE,LATCH,TRIP_DIRECTION,CTRL_ID,& 
                SETPOINT,DEVC_ID,FLOWRATE,DELAY,BYPASS_FLOWRATE,SURF_ID,STATISTICS,PART_ID,MATL_ID,SPEC_ID,DRY,CONVERSION_FACTOR,&
                UNITS,TIME_AVERAGED,VELO_INDEX,DUCT_ID,NODE_ID,POINTS,EVACUATION,X_ID,Y_ID,Z_ID,HIDE_COORDINATES

! Read the input file and count the number of DEVC lines

N_DEVC = 0
N_DEVC_READ = 0
MAX_DEVC_LINE_POINTS = 0
N_DEVC_TIME = 0
N_DEVC_LINE = 0

REWIND(LU_INPUT)
COUNT_DEVC_LOOP: DO
   CALL CHECKREAD('DEVC',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_DEVC_LOOP
   POINTS = 1
   READ(LU_INPUT,NML=DEVC,END=11,ERR=12,IOSTAT=IOS)
   N_DEVC      = N_DEVC      + POINTS
   N_DEVC_READ = N_DEVC_READ + 1
   MAX_DEVC_LINE_POINTS = MAX(MAX_DEVC_LINE_POINTS,POINTS)
   12 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I4)') 'ERROR: Problem with DEVC number ',N_DEVC_READ+1
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO COUNT_DEVC_LOOP
11 REWIND(LU_INPUT)

IF (N_DEVC==0) RETURN

! Allocate DEVICE array to hold all information for each device

ALLOCATE(DEVICE(N_DEVC),STAT=IZERO)
CALL ChkMemErr('READ','DEVICE',IZERO)
 
! Read in the DEVC lines, keeping track of TIME-history devices, and LINE array devices

N_DEVC      = 0

READ_DEVC_LOOP: DO NN=1,N_DEVC_READ

   CALL CHECKREAD('DEVC',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_DEVC_LOOP
   CALL SET_DEVC_DEFAULTS
   READ(LU_INPUT,DEVC) 

   ! Error statement of POINTS>1 and XB not specified

   IF (POINTS>1 .AND. ANY(XB<-1.E5_EB)) THEN
      WRITE(MESSAGE,'(A,A,A)')  'ERROR: DEVC ',TRIM(ID),' must have coordinates given in terms of XB'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (POINTS>1 .AND. STATISTICS/='null') THEN
      WRITE(MESSAGE,'(A,A,A)')  'ERROR: DEVC ',TRIM(ID),' cannot use POINTS>1 and STATISTICS at the same time'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   ! Reorder XB coordinates if necessary

   IF (POINTS==1) CALL CHECK_XB(XB)

   ! Process the point devices along a line, if necessary

   POINTS_LOOP: DO I_POINT=1,POINTS

      IF (XB(1)>-1.E5_EB) THEN
         IF (POINTS > 1) THEN
            XYZ(1) = XB(1) + (XB(2)-XB(1))*REAL(I_POINT-1,EB)/REAL(MAX(POINTS-1,1),EB)
            XYZ(2) = XB(3) + (XB(4)-XB(3))*REAL(I_POINT-1,EB)/REAL(MAX(POINTS-1,1),EB)
            XYZ(3) = XB(5) + (XB(6)-XB(5))*REAL(I_POINT-1,EB)/REAL(MAX(POINTS-1,1),EB)
         ELSE
            XYZ(1) = XB(1) + (XB(2)-XB(1))/2._EB
            XYZ(2) = XB(3) + (XB(4)-XB(3))/2._EB
            XYZ(3) = XB(5) + (XB(6)-XB(5))/2._EB
         ENDIF
      ELSE
         IF (XYZ(1) < -1.E5_EB .AND. DUCT_ID=='null' .AND. NODE_ID=='null') THEN
            WRITE(MESSAGE,'(A,A,A)')  'ERROR: DEVC ',TRIM(ID),' must have coordinates, even if it is not a point quantity'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDIF

      ! Determine which mesh the device is in
   
      BAD = .TRUE.
      MESH_LOOP: DO NM=1,NMESHES
         IF (EVACUATION_ONLY(NM)) CYCLE MESH_LOOP
         M=>MESHES(NM)
         IF (XYZ(1)>=M%XS .AND. XYZ(1)<=M%XF .AND. XYZ(2)>=M%YS .AND. XYZ(2)<=M%YF .AND. XYZ(3)>=M%ZS .AND. XYZ(3)<=M%ZF) THEN
            MESH_NUMBER = NM
            BAD = .FALSE.
            EXIT MESH_LOOP
         ENDIF
      ENDDO MESH_LOOP
   
      ! Make sure there is either a QUANTITY or PROP_ID for the DEVICE
   
      IF (QUANTITY=='null' .AND. PROP_ID=='null') THEN
         WRITE(MESSAGE,'(A,A,A)')  'ERROR: DEVC ',TRIM(ID),' must have either an output QUANTITY or PROP_ID'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   
      IF (BAD) THEN
         IF (DUCT_ID/='null' .OR. NODE_ID/='null') THEN
            XYZ(1) = MESHES(1)%XS
            XYZ(2) = MESHES(1)%YS
            XYZ(3) = MESHES(1)%ZS
            MESH_NUMBER = 1
         ELSE
            IF (ALL(EVACUATION_ONLY)) CYCLE READ_DEVC_LOOP
            WRITE(MESSAGE,'(A,A,A)') 'WARNING: DEVC ',TRIM(ID),' is not within any mesh.'  
            IF (MYID==0) WRITE(LU_ERR,'(A)') TRIM(MESSAGE)
            CYCLE READ_DEVC_LOOP
         ENDIF
      ENDIF

      ! Determine if the DEVC is a TIME or LINE device

      IF (POINTS==1)                   N_DEVC_TIME = N_DEVC_TIME + 1
      IF (POINTS>1 .AND. I_POINT==1)   N_DEVC_LINE = N_DEVC_LINE + 1
      
      ! Assign properties to the DEVICE array
   
      N_DEVC = N_DEVC + 1

      DV => DEVICE(N_DEVC)
      M  => MESHES(MESH_NUMBER)
   
      DV%CONVERSION_FACTOR= CONVERSION_FACTOR
      DV%DEPTH            = DEPTH
      DV%IOR              = IOR
      IF (ID=='null') WRITE(ID,'(A7,I4.4)') 'Device_',N_DEVC
      DV%ID               = ID
      IF (POINTS>1)DV%LINE= N_DEVC_LINE
      DV%POINT            = I_POINT
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
      DV%DUCT_ID          = DUCT_ID 
      DV%NODE_ID          = NODE_ID   
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
      IF (X_ID=='null') X_ID = TRIM(ID)//'-x'
      IF (Y_ID=='null') Y_ID = TRIM(ID)//'-y'
      IF (Z_ID=='null') Z_ID = TRIM(ID)//'-z'
      DV%X_ID             = X_ID
      DV%Y_ID             = Y_ID
      DV%Z_ID             = Z_ID
      DV%DRY              = DRY
      DV%EVACUATION       = EVACUATION
      DV%VELO_INDEX       = VELO_INDEX

      IF (POINTS > 1) THEN
         IF (.NOT.HIDE_COORDINATES) THEN
            IF (ABS(XB(1)-XB(2))> ZERO_P .AND. ABS(XB(3)-XB(4))<=ZERO_P .AND. ABS(XB(5)-XB(6))<=ZERO_P) DV%LINE_COORD_CODE = 1
            IF (ABS(XB(1)-XB(2))<=ZERO_P .AND. ABS(XB(3)-XB(4))> ZERO_P .AND. ABS(XB(5)-XB(6))<=ZERO_P) DV%LINE_COORD_CODE = 2
            IF (ABS(XB(1)-XB(2))<=ZERO_P .AND. ABS(XB(3)-XB(4))<=ZERO_P .AND. ABS(XB(5)-XB(6))> ZERO_P) DV%LINE_COORD_CODE = 3
            IF (ABS(XB(1)-XB(2))> ZERO_P .AND. ABS(XB(3)-XB(4))> ZERO_P .AND. ABS(XB(5)-XB(6))<=ZERO_P) DV%LINE_COORD_CODE = 12
            IF (ABS(XB(1)-XB(2))> ZERO_P .AND. ABS(XB(3)-XB(4))<=ZERO_P .AND. ABS(XB(5)-XB(6))> ZERO_P) DV%LINE_COORD_CODE = 13
            IF (ABS(XB(1)-XB(2))<=ZERO_P .AND. ABS(XB(3)-XB(4))> ZERO_P .AND. ABS(XB(5)-XB(6))> ZERO_P) DV%LINE_COORD_CODE = 23
         ELSE
            DV%LINE_COORD_CODE = 0
         ENDIF
      ENDIF
    
   ENDDO POINTS_LOOP

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
      IF (ABS(XB(1)-XB(2))<=ZERO_P .AND. IOR==0) DV%IOR = 1
      IF (ABS(XB(3)-XB(4))<=ZERO_P .AND. IOR==0) DV%IOR = 2
      IF (ABS(XB(5)-XB(6))<=ZERO_P .AND. IOR==0) DV%IOR = 3
   ENDIF
   
ENDDO READ_DEVC_LOOP

REWIND(LU_INPUT)

CONTAINS

SUBROUTINE SET_DEVC_DEFAULTS

CONVERSION_FACTOR = 1._EB
DEPTH            = 0._EB
IOR              = 0
ID               = 'null'
ORIENTATION(1:3) = (/0._EB,0._EB,-1._EB/)
PROP_ID          = 'null'
CTRL_ID          = 'null'
DEVC_ID          = 'null'
SURF_ID          = 'null'
PART_ID          = 'null'
MATL_ID          = 'null'
SPEC_ID          = 'null'
DUCT_ID          = 'null'
NODE_ID          = 'null'
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
POINTS           = 1
SETPOINT         = 1.E20_EB
STATISTICS       = 'null'
TRIP_DIRECTION   = 1
TIME_AVERAGED    = .TRUE.
UNITS            = 'null'
VELO_INDEX       = 0
XYZ              = -1.E6_EB
X_ID             = 'null'
Y_ID             = 'null'
Z_ID             = 'null'
HIDE_COORDINATES = .FALSE.
DRY              = .FALSE.
EVACUATION       = .FALSE.

END SUBROUTINE SET_DEVC_DEFAULTS

END SUBROUTINE READ_DEVC



SUBROUTINE READ_CTRL

! Just read in the ConTRoL parameters and store in the array CONTROL

USE CONTROL_VARIABLES
USE MATH_FUNCTIONS, ONLY : GET_RAMP_INDEX

LOGICAL :: INITIAL_STATE, LATCH, EVACUATION
INTEGER :: CYCLES,N,NC
REAL(EB) :: SETPOINT(2), DELAY, CYCLE_TIME
CHARACTER(30) :: ID,FUNCTION_TYPE,INPUT_ID(40),RAMP_ID,ON_BOUND
TYPE (CONTROL_TYPE), POINTER :: CF
NAMELIST /CTRL/ ID,LATCH,INITIAL_STATE,FUNCTION_TYPE,SETPOINT,DELAY,CYCLE_TIME,INPUT_ID,RAMP_ID,CYCLES,N,ON_BOUND, &
                EVACUATION
 
N_CTRL = 0
REWIND(LU_INPUT)
COUNT_CTRL_LOOP: DO
   CALL CHECKREAD('CTRL',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_CTRL_LOOP
   READ(LU_INPUT,NML=CTRL,END=11,ERR=12,IOSTAT=IOS)
   N_CTRL = N_CTRL + 1
   12 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I4)') 'ERROR: Problem with CTRL number ',N_CTRL+1
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
      WRITE(MESSAGE,'(A,I5,A)')  'ERROR: CTRL ',NC,' must have a FUNCTION_TYPE'
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
   CF%EVACUATION      = EVACUATION
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
         WRITE(MESSAGE,'(A,I5,A)')  'ERROR: CTRL ',NC,' FUNCTION_TYPE not recognized'
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
   EVACUATION    = .FALSE.
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
      WRITE(MESSAGE,'(A,I5,A)')  'ERROR: CTRL ',NC,' must have at least one input'
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
               WRITE(MESSAGE,'(A,I5,A)')  'ERROR: CUSTOM CTRL ',NC,' cannot have another CTRL as input'
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
      IF (ALL(EVACUATION_ONLY)) CYCLE BUILD_INPUT
   WRITE(MESSAGE,'(A,I5,A,A)')  'ERROR: CTRL ',NC,' cannot locate item for input ', TRIM(CF%INPUT_ID(NN))
   CALL SHUTDOWN(MESSAGE)
   END DO BUILD_INPUT

END DO PROC_CTRL_LOOP  
   
END SUBROUTINE PROC_CTRL   

SUBROUTINE PROC_OBST

INTEGER :: NM, N

MESH_LOOP: DO NM=1,NMESHES

   IF(EVACUATION_ONLY(NM)) CYCLE MESH_LOOP

   M=>MESHES(NM)
   DO N=1,M%N_OBST
      OB=>M%OBSTRUCTION(N)
      IF (OB%PROP_ID /='null') THEN
         CALL GET_PROPERTY_INDEX(OB%PROP_INDEX,'OBST',OB%PROP_ID)
      ENDIF
   END DO
END DO MESH_LOOP
  
END SUBROUTINE PROC_OBST


SUBROUTINE PROC_DEVC

! Process the DEViCes
USE COMP_FUNCTIONS, ONLY : CHANGE_UNITS
USE CONTROL_VARIABLES
USE DEVICE_VARIABLES, ONLY : DEVICE_TYPE, DEVICE, N_DEVC, PROPERTY, PROPERTY_TYPE, N_PROP, &
                             GAS_CELL_RAD_FLUX, GAS_CELL_RAD_DEVC_INDEX, N_GAS_CELL_RAD_DEVC
 
INTEGER  :: N,NN,NNN,NM,QUANTITY_INDEX,MAXCELLS,I,J,K,I_DUM
REAL(EB) :: XX,YY,ZZ,XX1,YY1,ZZ1,DISTANCE,SCANDISTANCE,DX,DY,DZ
TYPE (DEVICE_TYPE),  POINTER :: DV
TYPE (PROPERTY_TYPE),  POINTER :: PY
 
IF (N_DEVC==0) RETURN

! Set initial values for DEViCes

DEVICE(1:N_DEVC)%VALUE = 0._EB
DEVICE(1:N_DEVC)%TIME_INTERVAL = 0._EB

PROC_DEVC_LOOP: DO N=1,N_DEVC

   DV => DEVICE(N)
   
   ! Check for HVAC outputs with no HVAC inputs
   IF ((DV%DUCT_ID/='null' .OR. DV%NODE_ID/='null') .AND. .NOT. HVAC_SOLVE) THEN
      WRITE(MESSAGE,'(A)')  'ERROR: HVAC outputs specified with no HVAC inputs'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

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
         WRITE(MESSAGE,'(5A)')  'ERROR: DEVC ',TRIM(DV%ID),' or DEVC PROPerty ',TRIM(DV%PROP_ID),' must have a QUANTITY' 
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      IF (DV%QUANTITY=='null' .AND. PROPERTY(DV%PROP_INDEX)%QUANTITY/='null') DV%QUANTITY = PROPERTY(DV%PROP_INDEX)%QUANTITY
   ENDIF

   ! Check if the output QUANTITY exists and is appropriate

   IF (DV%QUANTITY /= 'null') THEN
      CALL GET_QUANTITY_INDEX(DV%SMOKEVIEW_LABEL,DV%SMOKEVIEW_BAR_LABEL,QUANTITY_INDEX,I_DUM, &
                              DV%SPEC_INDEX,DV%PART_INDEX,DV%DUCT_INDEX,DV%NODE_INDEX,'DEVC', &
                              DV%QUANTITY,'null',DV%SPEC_ID,DV%PART_ID,DV%DUCT_ID,DV%NODE_ID)
                              
      IF (OUTPUT_QUANTITY(QUANTITY_INDEX)%INTEGRATED .AND. DV%X1<=-1.E6_EB) THEN
         WRITE(MESSAGE,'(3A)')  'ERROR: DEVC QUANTITY ',TRIM(DV%QUANTITY),' requires coordinates using XB'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (.NOT.OUTPUT_QUANTITY(QUANTITY_INDEX)%INTEGRATED .AND. DV%STATISTICS=='null' .AND. DV%X1>-1.E6_EB .AND. DV%LINE==0) THEN
         WRITE(MESSAGE,'(3A)')  'ERROR: DEVC QUANTITY ',TRIM(DV%QUANTITY),' requires coordinates using XYZ, not XB'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (QUANTITY_INDEX<0 .AND. DV%IOR==0 .AND. DV%STATISTICS=='null') THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: Specify orientation of DEVC ',TRIM(DV%ID),' using the parameter IOR'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (QUANTITY_INDEX < 0 .AND. (DV%STATISTICS=='MASS MEAN' .OR. DV%STATISTICS=='VOLUME MEAN' .OR. &
                                    DV%STATISTICS=='VOLUME INTEGRAL' .OR. DV%STATISTICS=='MASS INTEGRAL') ) THEN
         WRITE(MESSAGE,'(A,A)') 'ERROR: Invalid STATISTICS specified for wall DEVC ',TRIM(DV%ID)
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (QUANTITY_INDEX > 0 .AND. DV%STATISTICS=='SURFACE INTEGRAL') THEN
         WRITE(MESSAGE,'(A,A)') 'ERROR: Invalid STATISTICS specified for gas DEVC ',TRIM(DV%ID)
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (QUANTITY_INDEX > 0 .AND. DV%STATISTICS/='null' .AND. DV%STATISTICS/='TIME INTEGRAL' .AND. DV%I1<0) THEN
         WRITE(MESSAGE,'(A,A)') 'ERROR: XB required when geometrical STATISTICS specified for gas DEVC ',TRIM(DV%ID)
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
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' is a smoke detector and must have a PROP_ID'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         IF (PROPERTY(DV%PROP_INDEX)%SPEC_INDEX==0 .AND. .NOT.SIMPLE_CHEMISTRY) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' is a smoke detector and requires a fire or smoke source'
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
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' must have a PROP_ID'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         IF (PROPERTY(DV%PROP_INDEX)%ACTIVATION_TEMPERATURE <= -273.15_EB) THEN
            WRITE(MESSAGE,'(A,A)') 'ERROR: ACTIVATION_TEMPERATURE needed for PROP ',TRIM(DV%PROP_ID)
            CALL SHUTDOWN(MESSAGE)
         ENDIF

         DV%SETPOINT = PROPERTY(DV%PROP_INDEX)%ACTIVATION_TEMPERATURE
         DV%TMP_L    = PROPERTY(DV%PROP_INDEX)%INITIAL_TEMPERATURE

      CASE ('THERMOCOUPLE')

         DV%TMP_L = PROPERTY(DV%PROP_INDEX)%INITIAL_TEMPERATURE

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
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' must have a MATL_ID'
            CALL SHUTDOWN(MESSAGE)
         ENDIF

      CASE ('CABLE TEMPERATURE')

         IF (DV%PROP_INDEX<1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' must have a PROP_ID'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         PY => PROPERTY(DV%PROP_INDEX)
         DV%SETPOINT = PY%CABLE_FAILURE_TEMPERATURE
         DV%TMP_L    = PY%INITIAL_TEMPERATURE
         GAS_CELL_RAD_FLUX   = .TRUE.
         N_GAS_CELL_RAD_DEVC = N_GAS_CELL_RAD_DEVC + 1
         DV%GAS_CELL_RAD_FLUX = .TRUE.
         IF (ABS(DV%DEPTH)<=ZERO_P) THEN
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
            IF (PROPERTY(DV%PROP_INDEX)%SPEC_INDEX==0 .AND. .NOT.SIMPLE_CHEMISTRY) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' is a smoke detector and requires a fire or smoke source'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
         ELSE
            IF (.NOT.SIMPLE_CHEMISTRY) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' is a smoke detector and requires a fire or smoke source'
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

         ! Check either for a specified SMOKE SPECies, or if simple chemistry model is being used
         IF (DV%PROP_INDEX>0) THEN
            IF (PROPERTY(DV%PROP_INDEX)%SPEC_INDEX==0 .AND. .NOT.SIMPLE_CHEMISTRY) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' is a smoke detector and requires a fire or smoke source'
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
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: ASPIRATION DEVICE ',TRIM(DV%ID),' is not listed after all its inputs'
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

   IF (DV%STATISTICS/='null') CALL CHANGE_UNITS(DV%QUANTITY,DV%UNITS,DV%STATISTICS,MYID,LU_ERR)

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
REAL(EB) :: XYZ(3),P0=0._EB,PX(3)=0._EB,PXX(3,3)=0._EB
CHARACTER(30) :: QUANTITY
TYPE (PROFILE_TYPE), POINTER :: PF
NAMELIST /PROF/ XYZ,QUANTITY,IOR,ID,FYI,P0,PX,PXX
 
N_PROF = 0
REWIND(LU_INPUT)
COUNT_PROF_LOOP: DO
   CALL CHECKREAD('PROF',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_PROF_LOOP
   READ(LU_INPUT,NML=PROF,END=11,ERR=12,IOSTAT=IOS)
   N_PROF = N_PROF + 1
   12 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I4)') 'ERROR: Problem with PROF number ',N_PROF+1
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
   
   NOT_PATCH_PROFILE_IF: IF (.NOT.PATCH_VELOCITY) THEN
 
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
   
   ELSE NOT_PATCH_PROFILE_IF
   
      ! Experimental patch velocity profile
      
      PF => PROFILE(N)
      PF%P0 = P0   ! value at origin of Taylor expansion
      DO J=1,3
         PF%PX(J) = PX(J)          ! first derivative of P evaluated at origin
         DO I=1,3
            PF%PXX(I,J) = PXX(I,J) ! second derivative of P evaluated at origin
         ENDDO
      ENDDO
      PF%ID  = ID
      PF%QUANTITY = QUANTITY
      
   ENDIF NOT_PATCH_PROFILE_IF
 
ENDDO PROF_LOOP
REWIND(LU_INPUT)
 
END SUBROUTINE READ_PROF



SUBROUTINE READ_ISOF
 
LOGICAL :: AUTO_HRRPUV_VALUE
REAL(EB) :: VALUE(10)
CHARACTER(30) :: QUANTITY,COLOR_QUANTITY,SPEC_ID,COLOR_SPEC_ID
INTEGER :: REDUCE_TRIANGLES,N,I_DUM,VELO_INDEX
TYPE(ISOSURFACE_FILE_TYPE), POINTER :: IS
NAMELIST /ISOF/ QUANTITY,FYI,VALUE,REDUCE_TRIANGLES,COLOR_QUANTITY,SPEC_ID,COLOR_SPEC_ID,VELO_INDEX,AUTO_HRRPUV_VALUE
 
N_ISOF = 0
REWIND(LU_INPUT)
COUNT_ISOF_LOOP: DO
   CALL CHECKREAD('ISOF',LU_INPUT,IOS) 
   IF (IOS==1) EXIT COUNT_ISOF_LOOP
   READ(LU_INPUT,NML=ISOF,END=9,ERR=10,IOSTAT=IOS)
   N_ISOF = N_ISOF + 1
   10 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I2)') 'ERROR: Problem with ISOF number ',N_ISOF
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
   VELO_INDEX       = 0
   AUTO_HRRPUV_VALUE      = .FALSE.
 
   CALL CHECKREAD('ISOF',LU_INPUT,IOS) 
   IF (IOS==1) EXIT READ_ISOF_LOOP
   READ(LU_INPUT,ISOF) 
 
   IS%REDUCE_TRIANGLES = REDUCE_TRIANGLES
   IS%VELO_INDEX       = VELO_INDEX

   CALL GET_QUANTITY_INDEX(IS%SMOKEVIEW_LABEL,IS%SMOKEVIEW_BAR_LABEL,IS%INDEX,I_DUM, &
                           IS%SPEC_INDEX,I_DUM,I_DUM,I_DUM,'ISOF', &
                           QUANTITY,'null',SPEC_ID,'null','null','null')
                                              
   IF (COLOR_QUANTITY/='null') &
      CALL GET_QUANTITY_INDEX(IS%SMOKEVIEW_LABEL2,IS%SMOKEVIEW_BAR_LABEL2,IS%INDEX2,I_DUM, &
                              IS%SPEC_INDEX2,I_DUM,I_DUM,I_DUM,'ISOF', &
                              COLOR_QUANTITY,'null',COLOR_SPEC_ID,'null','null','null')
                              
   IF(QUANTITY=='HRRPUV'.AND.AUTO_HRRPUV_VALUE)THEN
      IS%N_VALUES=1
      IS%VALUE(1)=MIN(HRRPUVCUT_MAX,20._EB/CHARACTERISTIC_CELL_SIZE)
   ELSE
      VALUE_LOOP: DO I=1,10
         IF (VALUE(I)<=-998._EB) EXIT VALUE_LOOP
         IS%N_VALUES = I
         IS%VALUE(I) = VALUE(I)
      ENDDO VALUE_LOOP
  ENDIF
 
ENDDO READ_ISOF_LOOP
REWIND(LU_INPUT)
 
END SUBROUTINE READ_ISOF
 
 
SUBROUTINE READ_SLCF

USE EVAC, ONLY: EVAC_FDS6
REAL(EB) :: MAXIMUM_VALUE,MINIMUM_VALUE
REAL(EB) :: AGL_SLICE
INTEGER :: N,NN,NM,MESH_NUMBER,N_SLCF_O,NITER,ITER,VELO_INDEX,I_DUM
LOGICAL :: VECTOR,RLE,TWO_BYTE,CELL_CENTERED, FIRE_LINE, EVACUATION,LEVEL_SET_FIRE_LINE
CHARACTER(30) :: QUANTITY,SPEC_ID,PART_ID,QUANTITY2
TYPE (SLICE_TYPE), POINTER :: SL
NAMELIST /SLCF/ XB,QUANTITY,FYI,PBX,PBY,PBZ,VECTOR,MESH_NUMBER,RLE,MAXIMUM_VALUE,MINIMUM_VALUE,TWO_BYTE,SPEC_ID, &
                AGL_SLICE,PART_ID,CELL_CENTERED,FIRE_LINE,EVACUATION,VELO_INDEX,LEVEL_SET_FIRE_LINE,QUANTITY2

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
      !IF (.NOT.EVAC_FDS6 .AND. EVACUATION_ONLY(NM) .AND. .NOT.EVACUATION_GRID(NM)) CYCLE COUNT_SLCF_LOOP
      IF (EVAC_FDS6 .AND. EVACUATION_ONLY(NM) .AND. .NOT.EVACUATION_GRID(NM)) CYCLE COUNT_SLCF_LOOP
      N_SLCF  = N_SLCF + 1
      IF (VECTOR .AND. TWO_D) N_SLCF = N_SLCF + 2
      IF (VECTOR .AND. .NOT. TWO_D) N_SLCF = N_SLCF + 3
      10 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with SLCF number ',N_SLCF_O+1
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
      QUANTITY  = 'null'
      QUANTITY2 = 'null'
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
      VELO_INDEX = 0
      LEVEL_SET_FIRE_LINE = .FALSE.
 
      CALL CHECKREAD('SLCF',LU_INPUT,IOS)
      IF (IOS==1) EXIT SLCF_LOOP
      READ(LU_INPUT,SLCF) 
      IF (MESH_NUMBER/=NM) CYCLE SLCF_LOOP
      IF (.NOT.EVACUATION_ONLY(NM) .AND.      EVACUATION) CYCLE SLCF_LOOP
      IF (     EVACUATION_ONLY(NM) .AND. .NOT.EVACUATION) CYCLE SLCF_LOOP
      IF (     EVAC_FDS6 .AND. EVACUATION_ONLY(NM) .AND. .NOT.EVACUATION_GRID(NM)) CYCLE SLCF_LOOP
      !IF (.NOT.EVAC_FDS6 .AND. EVACUATION_ONLY(NM) .AND. .NOT.EVACUATION_GRID(NM)) CYCLE SLCF_LOOP
 
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
         IF (ITER==1 .AND. LEVEL_SET_FIRE_LINE) QUANTITY = 'TEMPERATURE'
         IF (ITER==1) THEN
            SL%FIRE_LINE = FIRE_LINE
            SL%LEVEL_SET_FIRE_LINE = LEVEL_SET_FIRE_LINE
            IF (LEVEL_SET_FIRE_LINE .AND. .NOT. VEG_LEVEL_SET) THEN
               WRITE(MESSAGE,'(A)') "ERROR: VEG_LEVEL_SET must be TRUE on MISC line to run the LS model"
               CALL SHUTDOWN(MESSAGE)
            ENDIF
         ELSE
            SL%FIRE_LINE = .FALSE.
            SL%LEVEL_SET_FIRE_LINE = .FALSE.
            SPEC_ID      = 'null'
         ENDIF
         SL%VELO_INDEX = VELO_INDEX
         CALL GET_QUANTITY_INDEX(SL%SMOKEVIEW_LABEL,SL%SMOKEVIEW_BAR_LABEL,SL%INDEX,SL%INDEX2, &
                                 SL%SPEC_INDEX,SL%PART_INDEX,I_DUM,I_DUM,'SLCF', &
                                 QUANTITY,QUANTITY2,SPEC_ID,PART_ID,'null','null')
                                 
         ! For terrain slices, AGL=above ground level
         ! FIRE_LINE==.TRUE. means a terrain slice at one grid cell above ground with quantity temperature.
         ! Smokeview will display only regions where temperature is above 200 C. This is currently hard wired.
         ! LEVEL_SET_FIRE_LINE = .TRUE. will create a slice file  !!! this is not fully functional !!!


         IF (ITER == 1 .AND. (AGL_SLICE > -1._EB .OR. FIRE_LINE .OR.  LEVEL_SET_FIRE_LINE ) ) THEN 
            SL%TERRAIN_SLICE = .TRUE.
            IF (FIRE_LINE) THEN
               SL%SMOKEVIEW_LABEL = "Fire line"
               SL%SMOKEVIEW_BAR_LABEL = "Fire_line"
            ENDIF 
            IF (LEVEL_SET_FIRE_LINE) THEN
               SL%SMOKEVIEW_LABEL = "Level Set Fire line"
               SL%SMOKEVIEW_BAR_LABEL = "LS_Fire_line"
            ENDIF 
            IF (AGL_SLICE <= -1._EB .AND. FIRE_LINE)           AGL_SLICE = M%Z(1) - M%Z(0)
            IF (AGL_SLICE <= -1._EB .AND. LEVEL_SET_FIRE_LINE) AGL_SLICE = 0._EB
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
         
         ! Disable cell centered for velocity

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
   M%K_AGL_SLICE = 0
   N = 0
   DO NN = 1, N_SLCF
      SL=>SLICE(NN)
      IF(SL%TERRAIN_SLICE) THEN
        TERRAIN_CASE = .TRUE.
        N = N + 1 
        M%K_AGL_SLICE(0:IBP1,0:JBP1,N) =  SL%SLICE_AGL*M%RDZ(1)
        ! Subtract one because bottom of domain will be accounted for when cycling through walls cells
        M%K_AGL_SLICE(0:IBP1,0:JBP1,N) =  MAX(0,M%K_AGL_SLICE(0:IBP1,0:JBP1,N)-1)
      ENDIF
   ENDDO

   N_SLCF_MAX = MAX(N_SLCF_MAX,N_SLCF) 

   IF(VEG_LEVEL_SET) THEN
     ALLOCATE(M%LS_Z_TERRAIN(0:IBP1,0:JBP1),STAT=IZERO) ; CALL ChkMemErr('READ','LS_Z_TERRAIN',IZERO)
   ENDIF
 
ENDDO MESH_LOOP

END SUBROUTINE READ_SLCF


SUBROUTINE READ_BNDF

USE DEVICE_VARIABLES
INTEGER :: N,I_DUM
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
   210 IF (IOS>0) THEN 
         WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with BNDF number ',N_BNDF+1
         CALL SHUTDOWN(MESSAGE)
       ENDIF
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

   CALL GET_QUANTITY_INDEX(BF%SMOKEVIEW_LABEL,BF%SMOKEVIEW_BAR_LABEL,BF%INDEX,I_DUM, &
                           BF%SPEC_INDEX,BF%PART_INDEX,I_DUM,I_DUM,'BNDF', &
                           QUANTITY,'null',SPEC_ID,PART_ID,'null','null')
                           
   ! Assign miscellaneous attributes to the boundary file

   BF%CELL_CENTERED = CELL_CENTERED

   ! Check to see if PROP_ID exists

   BF%PROP_INDEX = 0
   IF(PROP_ID/='null')  CALL GET_PROPERTY_INDEX(BF%PROP_INDEX,'BNDF',PROP_ID)

ENDDO READ_BNDF_LOOP
REWIND(LU_INPUT)
 
END SUBROUTINE READ_BNDF
 

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


SUBROUTINE GET_QUANTITY_INDEX(SMOKEVIEW_LABEL,SMOKEVIEW_BAR_LABEL,OUTPUT_INDEX,OUTPUT2_INDEX, &
                              SPEC_INDEX,PART_INDEX,DUCT_INDEX,NODE_INDEX,OUTTYPE, &
                              QUANTITY,QUANTITY2,SPEC_ID_IN,PART_ID,DUCT_ID,NODE_ID)

CHARACTER(*), INTENT(INOUT) :: QUANTITY
CHARACTER(*), INTENT(OUT) :: SMOKEVIEW_LABEL,SMOKEVIEW_BAR_LABEL
CHARACTER(*) :: SPEC_ID_IN,PART_ID,DUCT_ID,NODE_ID
CHARACTER(30) :: SPEC_ID
CHARACTER(*), INTENT(IN) :: OUTTYPE,QUANTITY2
INTEGER, INTENT(OUT) :: OUTPUT_INDEX,SPEC_INDEX,PART_INDEX,DUCT_INDEX,NODE_INDEX,OUTPUT2_INDEX
INTEGER :: ND,NS,NN

SPEC_ID = SPEC_ID_IN

IF (QUANTITY=='OPTICAL DENSITY'        .AND. SPEC_ID=='null' .AND. SIMPLE_CHEMISTRY) SPEC_ID='SOOT'
IF (QUANTITY=='EXTINCTION COEFFICIENT' .AND. SPEC_ID=='null' .AND. SIMPLE_CHEMISTRY) SPEC_ID='SOOT'
IF (QUANTITY=='SOOT VOLUME FRACTION'   .AND. SPEC_ID=='null' .AND. SIMPLE_CHEMISTRY) SPEC_ID='SOOT'
IF (QUANTITY=='VISIBILITY'             .AND. SPEC_ID=='null' .AND. SIMPLE_CHEMISTRY) SPEC_ID='SOOT'

IF (SPEC_ID=='fuel')            SPEC_ID='FUEL'
IF (SPEC_ID=='nitrogen')        SPEC_ID='NITROGEN'
IF (SPEC_ID=='oxygen')          SPEC_ID='OXYGEN'
IF (SPEC_ID=='water vapor')     SPEC_ID='WATER VAPOR'
IF (SPEC_ID=='carbon dioxide')  SPEC_ID='CARBON DIOXIDE'
IF (SPEC_ID=='carbon monoxide') SPEC_ID='CARBON MONOXIDE'
IF (SPEC_ID=='hydrogen')        SPEC_ID='HYDROGEN'
IF (SPEC_ID=='other')           SPEC_ID='OTHER'
IF (SPEC_ID=='soot')            SPEC_ID='SOOT'

SPEC_INDEX = 0
PART_INDEX = 0
DUCT_INDEX = 0
NODE_INDEX = 0
OUTPUT2_INDEX = 0

! Convert old mixture fraction based quantities

IF (QUANTITY=='fuel'.OR.QUANTITY=='fuel mass fraction') THEN
   SPEC_ID = 'FUEL'
   IF (QUANTITY=='fuel')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='fuel mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (QUANTITY=='oxygen'.OR.QUANTITY=='oxygen mass fraction') THEN
   SPEC_ID = 'OXYGEN'
   IF (QUANTITY=='oxygen')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='oxygen mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (QUANTITY=='nitrogen'.OR.QUANTITY=='nitrogen mass fraction') THEN
   SPEC_ID = 'NITROGEN'
   IF (QUANTITY=='nitrogen')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='nitrogen mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (QUANTITY=='water vapor'.OR.QUANTITY=='water vapor mass fraction') THEN
   SPEC_ID = 'WATER VAPOR'
   IF (QUANTITY=='water vapor')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='water vapor mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (QUANTITY=='carbon dioxide'.OR.QUANTITY=='carbon dioxide mass fraction') THEN
   SPEC_ID = 'CARBON DIOXIDE'
   IF (QUANTITY=='carbon dioxide')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='carbon dioxide mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (QUANTITY=='carbon monoxide'.OR.QUANTITY=='carbon monoxide mass fraction') THEN
   SPEC_ID = 'CARBON MONOXIDE'
   IF (QUANTITY=='carbon monoxide')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='carbon monoxide mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (QUANTITY=='hydrogen'.OR.QUANTITY=='hydrogen mass fraction') THEN
   SPEC_ID = 'HYDROGEN'
   IF (QUANTITY=='hydrogen')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='hydrogen mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (QUANTITY=='soot'.OR.QUANTITY=='soot mass fraction') THEN
   SPEC_ID = 'SOOT'
   IF (QUANTITY=='soot')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='soot mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF
IF (QUANTITY=='other'.OR.QUANTITY=='other mass fraction') THEN
   SPEC_ID = 'OTHER'
   IF (QUANTITY=='other')               QUANTITY   = 'VOLUME FRACTION'
   IF (QUANTITY=='other mass fraction') QUANTITY   = 'MASS FRACTION'
ENDIF

! Assign SPEC_INDEX when SPEC_ID is specified

IF (SPEC_ID/='null') THEN
   DO NS=1,N_SPECIES
      IF (SPEC_ID==SPECIES(NS)%ID) THEN
         SPEC_INDEX = NS
         EXIT
      ENDIF
   ENDDO
ENDIF

IF (SPEC_ID/='null' .AND. SPEC_INDEX==0) THEN
   IF (SPEC_ID=='FUEL' .AND. SIMPLE_CHEMISTRY) THEN
      SPEC_INDEX=1
   ELSE
      WRITE(MESSAGE,'(A,A,A)')  'ERROR: ',TRIM(SPEC_ID),' is neither explicitly specified nor part of the mixture fraction model'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDIF 

! Assign HVAC indexes

IF (DUCT_ID/='null') THEN
   DO ND = 1, N_DUCTS
      IF (DUCT_ID==DUCT(ND)%ID) THEN
         DUCT_INDEX = ND
         EXIT
      ENDIF
   ENDDO
ENDIF

IF (NODE_ID/='null') THEN
   DO NN = 1, N_DUCTNODES
      IF (NODE_ID==DUCTNODE(NN)%ID) THEN
         NODE_INDEX = NN
         EXIT
      ENDIF
   ENDDO
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
IF (SPEC_ID =='null') THEN
   DO NS=1,N_SPECIES
      IF (QUANTITY==TRIM(SPECIES(NS)%ID)) THEN
         QUANTITY   = 'MASS FRACTION'
         SPEC_ID = SPECIES(NS)%ID    
         SPEC_INDEX = NS
         EXIT
      ENDIF
      IF (QUANTITY==TRIM(SPECIES(NS)%ID)//'_VF') THEN
         QUANTITY   = 'VOLUME FRACTION'
         SPEC_ID = SPECIES(NS)%ID      
         SPEC_INDEX = NS
         EXIT
      ENDIF
      IF (QUANTITY==TRIM(SPECIES(NS)%ID)//'_FLUX_X') THEN
         QUANTITY   = 'MASS FLUX X'
         SPEC_ID = SPECIES(NS)%ID      
         SPEC_INDEX = NS
         EXIT
      ENDIF
      IF (QUANTITY==TRIM(SPECIES(NS)%ID)//'_FLUX_Y') THEN
         QUANTITY   = 'MASS FLUX Y'
         SPEC_ID = SPECIES(NS)%ID      
         SPEC_INDEX = NS
         EXIT
      ENDIF
      IF (QUANTITY==TRIM(SPECIES(NS)%ID)//'_FLUX_Z') THEN
         QUANTITY   = 'MASS FLUX Z'
         SPEC_ID = SPECIES(NS)%ID      
         SPEC_INDEX = NS
         EXIT
      ENDIF
      IF (QUANTITY==TRIM(SPECIES(NS)%ID)//'_FLUX') THEN
         QUANTITY   = 'MASS FLUX'
         SPEC_ID = SPECIES(NS)%ID      
         SPEC_INDEX = NS
         EXIT
      ENDIF
   ENDDO
ENDIF
 
! Convert old particle quantity names

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
   SPEC_ID = 'SOOT'
   SPEC_INDEX = SOOT_INDEX
ENDIF
IF (QUANTITY=='extinction coefficient') THEN
   QUANTITY = 'EXTINCTION COEFFICIENT'
   SPEC_ID = 'SOOT'
   SPEC_INDEX = SOOT_INDEX
ENDIF
IF (QUANTITY=='optical depth') THEN
   QUANTITY = 'OPTICAL DENSITY'
   SPEC_INDEX = SOOT_INDEX
ENDIF
IF (QUANTITY=='soot density') THEN
   QUANTITY = 'DENSITY'
   SPEC_ID = 'SOOT'
   SPEC_INDEX = SOOT_INDEX
ENDIF
IF (QUANTITY=='soot volume fraction') THEN
   QUANTITY = 'SOOT VOLUME FRACTION'
   SPEC_ID = 'SOOT'
   SPEC_INDEX = SOOT_INDEX
ENDIF

IF (TRIM(SPEC_ID)=='SOOT' .AND. SPEC_INDEX == 0) THEN
   WRITE(MESSAGE,'(A,A)')  'ERROR: No SOOT species specified for QUANTITY=',TRIM(QUANTITY)
   CALL SHUTDOWN(MESSAGE)
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
IF (QUANTITY=='path obscuration')              QUANTITY = 'PATH OBSCURATION'
IF (QUANTITY=='PRESSURE_COEFFICIENT')          QUANTITY = 'PRESSURE COEFFICIENT'
IF (QUANTITY=='spot obscuration')              QUANTITY = 'CHAMBER OBSCURATION'
IF (QUANTITY=='RADIATIVE_FLUX')                QUANTITY = 'RADIATIVE HEAT FLUX'
IF (QUANTITY=='RADIATIVE_FLUX_GAS')            QUANTITY = 'RADIATIVE HEAT FLUX GAS'
IF (QUANTITY=='WALL_TEMPERATURE')              QUANTITY = 'WALL TEMPERATURE'
IF (QUANTITY=='WALL_THICKNESS')                QUANTITY = 'WALL THICKNESS'
IF (QUANTITY=='WATER RADIATION LOSS')          QUANTITY = 'PARTICLE RADIATION LOSS'

! Loop over all possible output quantities and assign an index number to match the desired QUANTITY
DO ND=-N_OUTPUT_QUANTITIES,N_OUTPUT_QUANTITIES
   IF (OUTPUT_QUANTITY(ND)%NAME=='null') CYCLE
   IF (QUANTITY2==OUTPUT_QUANTITY(ND)%NAME) THEN
      
      OUTPUT2_INDEX=ND
      
      IF (OUTPUT_QUANTITY(ND)%SPEC_ID_REQUIRED .AND. SPEC_INDEX==0) THEN
         WRITE(MESSAGE,'(3A)')  'ERROR: Output QUANTITY2 ',TRIM(QUANTITY2),' requires a SPEC_ID'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      
      ! QUANTITY2 only works with SLCF at the moment
      IF (.NOT.OUTPUT_QUANTITY(ND)%SLCF_APPROPRIATE) THEN
          WRITE(MESSAGE,'(3A)')  'ERROR: The QUANTITY2 ',TRIM(QUANTITY2),' is not appropriate for SLCF'
          CALL SHUTDOWN(MESSAGE)
      ENDIF

   ENDIF
ENDDO

QUANTITY_INDEX_LOOP: DO ND=-N_OUTPUT_QUANTITIES,N_OUTPUT_QUANTITIES

  QUANTITY_IF:IF (QUANTITY==OUTPUT_QUANTITY(ND)%NAME) THEN

     OUTPUT_INDEX = ND

      IF (OUTPUT_QUANTITY(ND)%QUANTITY2_REQUIRED .AND. OUTPUT2_INDEX==0) THEN
         WRITE(MESSAGE,'(3A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),' requires a QUANTITY2'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

     IF (OUTPUT_QUANTITY(ND)%LUMPED_SPECIES_ONLY .AND. .NOT.SIMPLE_CHEMISTRY) THEN
        WRITE(MESSAGE,'(5A)')  'ERROR: ',TRIM(OUTTYPE),' QUANTITY ',TRIM(QUANTITY), &
                               ' can only be used when the SIMPLE_CHEMISTRY model is active'
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
     IF (OUTPUT_QUANTITY(ND)%DUCT_ID_REQUIRED .AND. DUCT_INDEX<1) THEN
        WRITE(MESSAGE,'(3A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),' requires a DUCT_ID'
        CALL SHUTDOWN(MESSAGE)
     ENDIF

     IF (OUTPUT_QUANTITY(ND)%NODE_ID_REQUIRED .AND. NODE_INDEX<1) THEN
        WRITE(MESSAGE,'(3A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),' requires a NODE_ID'
        CALL SHUTDOWN(MESSAGE)
     ENDIF


     IF (( QUANTITY=='RELATIVE HUMIDITY' .OR. QUANTITY=='HUMIDITY').AND. H2O_INDEX==0) THEN
        WRITE(MESSAGE,'(A)')  'ERROR: RELATIVE HUMIDTY and HUMIDTY require SPEC=WATER VAPOR, or SIMPLE_CHEMISTRY'
        CALL SHUTDOWN(MESSAGE)
     END IF

      SELECT CASE (TRIM(OUTTYPE))
         CASE ('SLCF')
            ! Throw out bad slices
            IF (.NOT. OUTPUT_QUANTITY(ND)%SLCF_APPROPRIATE) THEN
               WRITE(MESSAGE,'(3A)')  'ERROR: The QUANTITY ',TRIM(QUANTITY),' is not appropriate for SLCF'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
         CASE ('DEVC')
            IF (.NOT.OUTPUT_QUANTITY(ND)%DEVC_APPROPRIATE) THEN
               WRITE(MESSAGE,'(3A)')  'ERROR: The QUANTITY ',TRIM(QUANTITY),' is not appropriate for DEVC'
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
               WRITE(MESSAGE,'(3A)')  'ERROR: The QUANTITY ',TRIM(QUANTITY),' is not appropriate for BNDF'
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
             IF (SMOKE3D .AND. (.NOT.OUTPUT_QUANTITY(ND)%MASS_FRACTION .AND. ND/=11)) THEN
                WRITE(MESSAGE,'(5A)') 'ERROR: ',TRIM(OUTTYPE),'_QUANTITY ',TRIM(QUANTITY), ' must be a mass fraction'
                CALL SHUTDOWN(MESSAGE)
             ENDIF
         CASE DEFAULT
      END SELECT

     ! Assign Smokeview Label

     IF (SPEC_INDEX>0) THEN
        SMOKEVIEW_LABEL = TRIM(SPECIES(SPEC_INDEX)%ID)//' '//TRIM(QUANTITY)
        SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)//'_'//TRIM(SPECIES(SPEC_INDEX)%FORMULA)
     ELSEIF (PART_INDEX>0) THEN
        SMOKEVIEW_LABEL = TRIM(PARTICLE_CLASS(PART_INDEX)%ID)//' '//TRIM(QUANTITY)
        SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)
     ELSEIF (OUTPUT2_INDEX/=0) THEN
        SMOKEVIEW_LABEL = TRIM(QUANTITY)//' '//TRIM(QUANTITY2)
        SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)//' '//TRIM(OUTPUT_QUANTITY(OUTPUT2_INDEX)%SHORT_NAME)
     ELSE
        SMOKEVIEW_LABEL = TRIM(QUANTITY)
        SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)
     ENDIF

      RETURN
   ENDIF QUANTITY_IF
ENDDO QUANTITY_INDEX_LOOP

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
        CASE ('OBST')
        CASE ('BNDF')
        CASE ('PLOT3D')
        CASE DEFAULT
     END SELECT
     RETURN
  ENDIF
ENDDO

WRITE(MESSAGE,'(5A)')  'ERROR: ',TRIM(OUTTYPE),' PROP_ID ',TRIM(PROP_ID),' not found'
CALL SHUTDOWN(MESSAGE)

END SUBROUTINE GET_PROPERTY_INDEX


!! The following are experimental routines for specification of complex geometry. ~RJM
SUBROUTINE READ_GEOM

USE GLOBAL_CONSTANTS, ONLY: LU_INPUT,N_GEOM
USE MATH_FUNCTIONS, ONLY: NORM2
CHARACTER(30) :: ID,SHAPE,TFILE
REAL(EB) :: XYZ(3),ORIENTATION(3),XB(6),RADIUS,VELOCITY(3),OMEGA,RGB(3),PIXELS,ROUGHNESS
LOGICAL :: HOLE,TWO_SIDED
INTEGER :: IOS,IZERO,NG
TYPE(GEOMETRY_TYPE), POINTER :: G
NAMELIST /GEOM/ ID,SHAPE,XB,XYZ,ORIENTATION,RADIUS,VELOCITY,OMEGA,HOLE,RGB,PIXELS,TWO_SIDED,TFILE,ROUGHNESS

N_GEOM=0
REWIND(LU_INPUT)
COUNT_GEOM_LOOP: DO
   CALL CHECKREAD('GEOM',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_GEOM_LOOP
   READ(LU_INPUT,NML=GEOM,END=11,ERR=12,IOSTAT=IOS)
   N_GEOM=N_GEOM+1
   12 IF (IOS>0) CALL SHUTDOWN('ERROR: problem with GEOM line')
ENDDO COUNT_GEOM_LOOP
11 REWIND(LU_INPUT)

IF (N_GEOM==0) RETURN

! Allocate GEOMETRY array

ALLOCATE(GEOMETRY(N_GEOM),STAT=IZERO)
CALL ChkMemErr('READ','GEOMETRY',IZERO)

READ_GEOM_LOOP: DO NG=1,N_GEOM
   G=>GEOMETRY(NG)
   
   CALL CHECKREAD('GEOM',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_GEOM_LOOP
   
   ! Set defaults
   
   ID = 'geom'
   TFILE = 'null'
   SHAPE = 'SPHERE'
   HOLE = .FALSE.
   XYZ = 0._EB
   ORIENTATION = (/0._EB,0._EB,1._EB/)
   XB = 0._EB
   RADIUS = 1._EB
   VELOCITY = 0._EB
   OMEGA = 0._EB
   RGB = (/192,192,192/)
   PIXELS = 1.0
   TWO_SIDED = .FALSE.
   ROUGHNESS = 0._EB
   
   ! Read the GEOM line
   
   READ(LU_INPUT,GEOM,END=35)
   
   G%X1 = XB(1)
   G%X2 = XB(2)
   G%Y1 = XB(3)
   G%Y2 = XB(4)
   G%Z1 = XB(5)
   G%Z2 = XB(6)
   G%X0 = XYZ(1)
   G%Y0 = XYZ(2)
   G%Z0 = XYZ(3)
   G%X  = XYZ(1)
   G%Y  = XYZ(2)
   G%Z  = XYZ(3)
   G%XOR = ORIENTATION(1)/NORM2(ORIENTATION)
   G%YOR = ORIENTATION(2)/NORM2(ORIENTATION)
   G%ZOR = ORIENTATION(3)/NORM2(ORIENTATION)
   G%U0 = VELOCITY(1)
   G%V0 = VELOCITY(2)
   G%W0 = VELOCITY(3)
   G%OMEGA = OMEGA
   G%OMEGA_X = G%OMEGA*G%XOR
   G%OMEGA_Y = G%OMEGA*G%YOR
   G%OMEGA_Z = G%OMEGA*G%ZOR
   G%RADIUS = RADIUS
   G%SHAPE = SHAPE
   G%HOLE = HOLE
   G%RGB = RGB
   G%PIXELS = PIXELS
   G%TWO_SIDED = TWO_SIDED
   G%ID = ID
   G%TFILE = TFILE
   G%ROUGHNESS = ROUGHNESS
   
   IF (ABS(G%U0)>0._EB .OR. ABS(G%V0)>0._EB .OR. ABS(G%W0)>0._EB) G%TRANSLATE = .TRUE.
   IF (ABS(G%OMEGA)>0._EB) G%ROTATE = .TRUE.
   
   SELECT CASE(G%SHAPE)
      CASE('BOX')
         G%ISHAPE = IBOX
         G%X0 = 0.5_EB*(G%X1+G%X2)
         G%Y0 = 0.5_EB*(G%Y1+G%Y2)
         G%Z0 = 0.5_EB*(G%Z1+G%Z2)
         G%SMVOBJECT = 'movingbox'
      CASE('SPHERE')
         G%ISHAPE = ISPHERE
         G%SMVOBJECT = 'movingsphere'
      CASE('CYLINDER')
         G%ISHAPE = ICYLINDER
         G%SMVOBJECT = 'tube'
      CASE('PLANE')
         G%ISHAPE = IPLANE
         G%SMVOBJECT = 'plane'
      CASE DEFAULT
         CALL SHUTDOWN('ERROR: unrecognized SHAPE on GEOM line')
   END SELECT
   
   ! Allocate bounding box arrays
   
   ALLOCATE(G%MIN_I(NMESHES),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','MIN_I',IZERO)
   ALLOCATE(G%MAX_I(NMESHES),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','MAX_I',IZERO)
   ALLOCATE(G%MIN_J(NMESHES),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','MIN_J',IZERO)
   ALLOCATE(G%MAX_J(NMESHES),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','MAX_J',IZERO)
   ALLOCATE(G%MIN_K(NMESHES),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','MIN_K',IZERO)
   ALLOCATE(G%MAX_K(NMESHES),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','MAX_K',IZERO)
   
ENDDO READ_GEOM_LOOP
35 REWIND(LU_INPUT)

END SUBROUTINE READ_GEOM


SUBROUTINE GET_REV_read(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') readrev(INDEX(readrev,':')+1:LEN_TRIM(readrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') readdate

END SUBROUTINE GET_REV_read
 
END MODULE READ_INPUT

