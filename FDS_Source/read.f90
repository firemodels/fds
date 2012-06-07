MODULE READ_INPUT
 
USE PRECISION_PARAMETERS
USE MESH_VARIABLES
USE GLOBAL_CONSTANTS
USE TRAN
USE MESH_POINTERS
USE OUTPUT_DATA
USE COMP_FUNCTIONS, ONLY: SECOND, CHECKREAD, SHUTDOWN, CHECK_XB, SEARCH_INPUT_FILE
USE MEMORY_FUNCTIONS, ONLY: ChkMemErr
USE COMP_FUNCTIONS, ONLY: GET_INPUT_FILE
USE MISC_FUNCTIONS, ONLY: SEARCH_CONTROLLER
USE EVAC, ONLY: READ_EVAC
USE HVAC_ROUTINES, ONLY: READ_HVAC,PROC_HVAC
USE COMPLEX_GEOMETRY, ONLY: READ_GEOM,READ_VERT,READ_FACE,READ_VOLU
 
IMPLICIT NONE
PRIVATE
CHARACTER(255), PARAMETER :: readid='$Id$'
CHARACTER(255), PARAMETER :: readrev='$Revision$'
CHARACTER(255), PARAMETER :: readdate='$Date$'

PUBLIC READ_DATA, GET_REV_read

CHARACTER(30) :: LABEL,MB,ODE_SOLVER
CHARACTER(100) :: MESSAGE,FYI
CHARACTER(30) :: ID,SURF_DEFAULT='INERT',EVAC_SURF_DEFAULT='INERT'
LOGICAL :: SUCCESS,EX,THICKEN_OBSTRUCTIONS,BAD,SIMPLE_CHEMISTRY=.FALSE.,BACKGROUND_DECLARED=.FALSE.,IDEAL=.FALSE.
REAL(EB) :: XB(6),TEXTURE_ORIGIN(3)
REAL(EB) :: PBX,PBY,PBZ
REAL(EB) :: MW_MIN,MW_MAX
REAL(EB) :: REAC_ATOM_ERROR,REAC_MASS_ERROR
INTEGER  :: I,J,K,IZERO,IOS,BACKGROUND_SPEC_INDEX=-1,BACKGROUND_SMIX_INDEX=-1
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE(OBSTRUCTION_TYPE), POINTER :: OB=>NULL()
TYPE (VENTS_TYPE), POINTER :: VT=>NULL()
TYPE(SURFACE_TYPE), POINTER :: SF=>NULL()
TYPE(MATERIAL_TYPE), POINTER :: ML=>NULL()
TYPE(REACTION_TYPE), POINTER :: RN=>NULL()
 
 
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
      WRITE(LU_ERR,'(A,I5)') "SVN Revision Number: ",SVN_REVISION_NUMBER
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
CALL READ_MESH(1)    ! Read and initialize some evacuation mesh parameters
CALL READ_EVAC(1)    ! Read mesh structure related evacuation input
CALL READ_MESH(2)    ! Read the fire (and evacuation) meshes
CALL READ_TRAN
CALL READ_TIME
CALL READ_PRES
CALL READ_REAC
CALL READ_SPEC
CALL READ_SMIX
CALL PROC_REAC
CALL READ_RADI
CALL READ_PROP
CALL READ_PART
CALL READ_DEVC
CALL READ_CTRL
CALL READ_TREE
CALL READ_MATL
CALL READ_SURF
CALL READ_GEOM    ! Experimental complex geometry ~RJM
CALL READ_VERT
CALL READ_FACE
CALL READ_CSVF    ! read in csv file to pass on to smokeview ~GPF
CALL READ_VOLU
CALL READ_OBST
CALL READ_VENT
CALL READ_ZONE
CALL READ_EVAC(2) ! Read the main evacuation inputs (no mesh related)
CALL READ_HVAC
CALL PROC_SURF_1  ! Set up SURFace constructs for species
CALL READ_RAMP    ! Read in all RAMPs, assuming they have all been identified previously
CALL PROC_SMIX
CALL PROC_HVAC
CALL PROC_MATL    ! Set up various MATeriaL constructs
CALL PROC_SURF_2  ! Set up remaining SURFace constructs
CALL READ_DUMP
CALL READ_CLIP
CALL PROC_WALL    ! Set up grid for 1-D heat transfer in solids
CALL PROC_PART    ! Set up various PARTicle constructs
CALL READ_INIT
CALL READ_TABL    ! Read in all TABLs, assuming they have all been identified previously
CALL PROC_CTRL    ! Set up various ConTRoL constructs
CALL PROC_PROP    ! Set up various PROPerty constructs
CALL PROC_DEVC    ! Set up various DEViCe constructs
CALL PROC_OBST    ! if an obst has a PROP then make sure PROP exists
CALL READ_PROF
CALL READ_SLCF
CALL READ_ISOF
CALL READ_BNDF
CALL READ_BNDE

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
NAMELIST /HEAD/ CHID,FYI,TITLE
 
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
 
 
SUBROUTINE READ_MESH(IMODE)
USE EVAC, ONLY: N_DOORS, N_EXITS, N_CO_EXITS, EVAC_EMESH_EXITS_TYPE, EMESH_EXITS, EMESH_ID, EMESH_IJK, EMESH_XB, &
     EMESH_NM, N_DOOR_MESHES, EMESH_NFIELDS, EVAC_FDS6, HUMAN_SMOKE_HEIGHT, EVAC_DELTA_SEE, &
     EMESH_STAIRS, EVAC_EMESH_STAIRS_TYPE, N_STRS, INPUT_EVAC_GRIDS, NO_EVAC_MESHES

INTEGER, INTENT(IN) :: IMODE
INTEGER :: IJK(3),NM,CURRENT_MPI_PROCESS,MPI_PROCESS,RGB(3),LEVEL,N_MESH_NEW,N,II,JJ,KK,NMESHES_READ,NNN,NEVAC_MESHES
LOGICAL :: EVACUATION, EVAC_HUMANS
REAL(EB) :: EVAC_Z_OFFSET,XB1,XB2,XB3,XB4,XB5,XB6
CHARACTER(25) :: COLOR
CHARACTER(30) :: MULT_ID
NAMELIST /MESH/ COLOR,CYLINDRICAL,EVACUATION,EVAC_HUMANS,EVAC_Z_OFFSET, FYI,ID,IJK,LEVEL,MPI_PROCESS,MULT_ID,RGB,&
                SYNCHRONIZE,XB
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (MULTIPLIER_TYPE), POINTER :: MR=>NULL()
 
NMESHES = 0
NMESHES_READ = 0
NEVAC_MESHES = 0
IF (IMODE==1) THEN
   NO_EVAC_MESHES = .TRUE.
   INPUT_EVAC_GRIDS   = 0
   IF(NO_EVACUATION) THEN
      N_EVAC = 0
      RETURN
   END IF
END IF

REWIND(LU_INPUT)
COUNT_MESH_LOOP: DO
   CALL CHECKREAD('MESH',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_MESH_LOOP
   MULT_ID = 'null'
   EVACUATION  = .FALSE.
   EVAC_HUMANS  = .FALSE.
   READ(LU_INPUT,MESH,END=15,ERR=16,IOSTAT=IOS)
   NMESHES_READ = NMESHES_READ + 1
   IF (NO_EVACUATION .AND. EVACUATION) CYCLE COUNT_MESH_LOOP ! skip evacuation meshes
   IF (EVACUATION_DRILL .AND. .NOT.EVACUATION) CYCLE COUNT_MESH_LOOP ! skip fire meshes
   IF (EVACUATION_MC_MODE .AND. .NOT.EVACUATION) CYCLE COUNT_MESH_LOOP ! skip fire meshes
   IF (EVACUATION .AND. EVAC_FDS6) NEVAC_MESHES = NEVAC_MESHES + 1
   IF (IMODE==1 .AND. EVAC_HUMANS) NO_EVAC_MESHES = .FALSE.
   IF (IMODE==1 .AND. EVAC_HUMANS) INPUT_EVAC_GRIDS = INPUT_EVAC_GRIDS + 1
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

EVAC_MODE_IF: IF (IMODE==1) THEN
   REWIND(LU_INPUT)
   IF (NO_EVAC_MESHES) THEN
      NO_EVACUATION      = .TRUE.
      EVACUATION_DRILL   = .FALSE.
      EVACUATION_MC_MODE = .FALSE.
      N_EVAC             = 0
      RETURN
   END IF
   ALLOCATE(EMESH_ID(     MAX(1,INPUT_EVAC_GRIDS)), STAT=IZERO)
   CALL ChkMemErr('READ_EVAC','EMESH_ID',IZERO) 
   ALLOCATE(EMESH_XB(6,   MAX(1,INPUT_EVAC_GRIDS)), STAT=IZERO)
   CALL ChkMemErr('READ_EVAC','EMESH_XB',IZERO) 
   ALLOCATE(EMESH_IJK(3,  MAX(1,INPUT_EVAC_GRIDS)), STAT=IZERO)
   CALL ChkMemErr('READ_EVAC','EMESH_IJK',IZERO) 

   NM = 0
   EVAC_MESH_LOOP: DO N = 1, NMESHES_READ
      ! Set evacuation MESH defaults
      IJK(1)= 10
      IJK(2)= 10
      IJK(3)= 1
      XB(1) = 0._EB
      XB(2) = 1._EB
      XB(3) = 0._EB
      XB(4) = 1._EB
      XB(5) = 0._EB
      XB(6) = 1._EB
      RGB   = -1
      COLOR   = 'null'
      ID      = 'null'
      EVACUATION  = .FALSE.
      EVAC_HUMANS = .FALSE.
      ! Read the MESH line
      CALL CHECKREAD('MESH', LU_INPUT, IOS)
      IF (IOS==1) EXIT EVAC_MESH_LOOP
      READ(LU_INPUT, MESH)
      IF (.NOT.EVACUATION)                   CYCLE EVAC_MESH_LOOP ! skip fire meshes
      IF (.NOT.EVAC_HUMANS .AND. EVACUATION) CYCLE EVAC_MESH_LOOP ! skip additional evac meshes
      NM = NM + 1
      ! Reorder XB coordinates if necessary
      CALL CHECK_XB(XB)
      EMESH_ID(NM)    = TRIM(ID)
      EMESH_IJK(1,NM) = IJK(1)
      EMESH_IJK(2,NM) = IJK(2)
      EMESH_IJK(3,NM) = IJK(3)
      EMESH_XB(1,NM)  = XB(1)
      EMESH_XB(2,NM)  = XB(2)
      EMESH_XB(3,NM)  = XB(3)
      EMESH_XB(4,NM)  = XB(4)
      EMESH_XB(5,NM)  = XB(5)
      EMESH_XB(6,NM)  = XB(6)
   END DO EVAC_MESH_LOOP
   REWIND(LU_INPUT)
   RETURN
END IF EVAC_MODE_IF

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
            IF (EVACUATION .AND. ABS(XB5 - XB6) <= SPACING(XB(6))) THEN
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
          EVAC_Z_OFFSET = EMESH_STAIRS(N)%EVAC_Z_OFFSET
          
          ! Increase the MESH counter by 1

          NM = NM + 1
          EMESH_STAIRS(N)%IMESH = NM

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

    IF (ALL(EVACUATION_ONLY)) THEN
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
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (TRAN_TYPE), POINTER :: T=>NULL()
NAMELIST /TRNX/ CC,FYI,IDERIV,MESH_NUMBER,PC
NAMELIST /TRNY/ CC,FYI,IDERIV,MESH_NUMBER,PC
NAMELIST /TRNZ/ CC,FYI,IDERIV,MESH_NUMBER,PC
 
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

REAL(EB) :: DT,VEL_CHAR
INTEGER :: NM
NAMELIST /TIME/ DT,EVAC_DT_FLOWFIELD,EVAC_DT_STEADY_STATE,FYI,LIMITING_DT_RATIO,LOCK_TIME_STEP,RESTRICT_TIME_STEP, &
                SYNCHRONIZE,T_BEGIN,T_END,TIME_SHRINK_FACTOR,WALL_INCREMENT, &
                TWFIN ! Backward compatibility
TYPE (MESH_TYPE), POINTER :: M=>NULL()
 
DT                   = -1._EB
EVAC_DT_FLOWFIELD    = 0.01_EB
EVAC_DT_STEADY_STATE = 0.05_EB
TIME_SHRINK_FACTOR   = 1._EB
T_BEGIN              = 0._EB
T_END                = 1._EB
SYNCHRONIZE          = .FALSE.
IF (ANY(SYNC_TIME_STEP)) SYNCHRONIZE = .TRUE.
WALL_INCREMENT       = 2
 
REWIND(LU_INPUT)
READ_TIME_LOOP: DO
   CALL CHECKREAD('TIME',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_TIME_LOOP
   READ(LU_INPUT,TIME,END=21,ERR=22,IOSTAT=IOS)
   22 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with TIME line')
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

   IF (.NOT.EVACUATION_ONLY(NM)) CHARACTERISTIC_CELL_SIZE = MIN( CHARACTERISTIC_CELL_SIZE , M%CELL_SIZE )

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
TYPE(MULTIPLIER_TYPE), POINTER :: MR=>NULL()
NAMELIST /MULT/ DX,DXB,DX0,DY,DY0,DZ,DZ0,FYI,ID,I_LOWER,I_UPPER,J_LOWER,J_UPPER,K_LOWER,K_UPPER,N_LOWER,N_UPPER

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

USE MATH_FUNCTIONS, ONLY: GET_RAMP_INDEX
REAL(EB) :: C_HORIZONTAL,C_VERTICAL,FORCE_VECTOR(3)=0._EB,MAXIMUM_VISIBILITY
CHARACTER(30) :: RAMP_GX,RAMP_GY,RAMP_GZ,TURBULENCE_MODEL
NAMELIST /MISC/ AL2O3,ALLOW_SURFACE_PARTICLES,ALLOW_UNDERSIDE_PARTICLES,ASSUMED_GAS_TEMPERATURE,BAROCLINIC,BNDF_DEFAULT,&
                CDF_CUTOFF,CFL_MAX,CFL_MIN,CFL_VELOCITY_NORM,CHECK_GR,CHECK_HT,CHECK_VN,CLIP_MASS_FRACTION,&
                CONSTANT_SPECIFIC_HEAT,C_DEARDORFF,C_SMAGORINSKY,C_VREMAN,C_FORCED,&
                C_FORCED_CYLINDER,C_FORCED_SPHERE,C_G,C_HORIZONTAL,C_VERTICAL,DNS,ENTHALPY_TRANSPORT,&
                EVACUATION_DRILL,EVACUATION_MC_MODE,EVAC_PRESSURE_ITERATIONS,EVAC_TIME_ITERATIONS,FLUX_LIMITER,FORCE_VECTOR,&
                FREEZE_VELOCITY,FYI,GAMMA,GRAVITATIONAL_DEPOSITION,GROUND_LEVEL,GVEC,HRRPUVCUT_MAX,HRRPUV_MAX_SMV, &
                IMMERSED_BOUNDARY_METHOD,INITIAL_UNMIXED_FRACTION,LAPSE_RATE,&
                LES,MAXIMUM_VISIBILITY,MEAN_FORCING,MIXING_LAYER_H0,MIXING_LAYER_THETA0,MIXING_LAYER_U0,&
                NEW_EVAP,NOISE,NOISE_VELOCITY,NO_EVACUATION,&
                OVERWRITE,PERIODIC_TEST,POROUS_FLOOR,PR,PROJECTION,P_INF,RAMP_GX,RAMP_GY,RAMP_GZ,RICHARDSON_ERROR_TOLERANCE,&
                RESTART,RESTART_CHID,RUN_AVG_FAC,SC,SECOND_ORDER_PARTICLE_TRANSPORT,SHARED_FILE_SYSTEM,&
                SMOKE_ALBEDO,SOLID_PHASE_ONLY,&
                STRATIFICATION,TENSOR_DIFFUSIVITY,TERRAIN_CASE,TERRAIN_IMAGE,TEXTURE_ORIGIN,THERMOPHORETIC_DEPOSITION,&
                THICKEN_OBSTRUCTIONS,TMPA,TURBULENT_DEPOSITION,TURBULENCE_MODEL,U0,UVW_FILE,V0,VEG_LEVEL_SET,VISIBILITY_FACTOR,&
                VN_MAX,VN_MIN,W0,WFDS,WIND_ONLY, &
                EVAC_SURF_DEFAULT,SUPPRESSION,SURF_DEFAULT,RADIATION ! Backward compatibility


! Physical constants
 
TMPA         = 20._EB                                              ! Ambient temperature (C)
GAMMA        = 1.4_EB                                              ! Heat capacity ratio for air
P_INF        = 101325._EB                                          ! Ambient pressure (Pa)
SIGMA        = 5.67E-8_EB                                          ! Stefan-Boltzmann constant (W/m**2/K**4)
MU_AIR_0     = 1.8E-5_EB                                           ! Dynamic Viscosity of Air at 20 C (kg/m/s)
CP_AIR_0     = 1012._EB                                            ! Specific Heat of Air at 20 C (J/kg/K)
PR_AIR       = 0.7_EB
K_AIR_0      = MU_AIR_0*CP_AIR_0/PR_AIR                            ! Thermal Conductivity of Air at 20 C (W/m/K)
RHO_SOOT     = 1850._EB                                            ! Density of soot particle (kg/m3)

! Empirical heat transfer constants
 
C_VERTICAL        = 1.31_EB  ! Vertical free convection (Holman, Table 7-2)
C_HORIZONTAL      = 1.52_EB  ! Horizontal free convection 
C_FORCED          = 0.037_EB ! Forced convection coefficient for plates
C_FORCED_CYLINDER = 0.664_EB ! Forced convection coefficient for cylinders
C_FORCED_SPHERE   = 0.6_EB   ! Forced convection coefficient for spheres
PR_ONTH           = PR_AIR**ONTH
ASSUMED_GAS_TEMPERATURE = -1000.   ! Assumed gas temperature, used for diagnostics
 
! Background parameters
 
U0 = 0._EB
V0 = 0._EB
W0 = 0._EB

! Miscellaneous constants

RESTART_CHID   = CHID
RESTART        = .FALSE.
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
 
PR                   = -1.0_EB  ! Turbulent Prandtl number
SC                   = -1.0_EB  ! Turbulent Schmidt number
 
! Misc
 
RAMP_GX              = 'null'
RAMP_GY              = 'null'
RAMP_GZ              = 'null'
GVEC(1)              = 0._EB        ! x-component of gravity 
GVEC(2)              = 0._EB        ! y-component of gravity 
GVEC(3)              = -GRAV        ! z-component of gravity 
LAPSE_RATE           = 0._EB       
THICKEN_OBSTRUCTIONS = .FALSE.
CFL_MAX              = 1.0_EB       ! Stability bounds
CFL_MIN              = 0.8_EB
VN_MAX               = 1.0_EB
VN_MIN               = 0.8_EB
VEG_LEVEL_SET        = .FALSE.
TERRAIN_IMAGE        = 'xxxnull'
MAXIMUM_VISIBILITY   = 30._EB    !m
VISIBILITY_FACTOR    = 3._EB
TURBULENCE_MODEL     = 'null'
 
REWIND(LU_INPUT)
MISC_LOOP: DO 
   CALL CHECKREAD('MISC',LU_INPUT,IOS)
   IF (IOS==1) EXIT MISC_LOOP
   READ(LU_INPUT,MISC,END=23,ERR=24,IOSTAT=IOS)
   24 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with MISC line')
ENDDO MISC_LOOP
23 REWIND(LU_INPUT)

! FDS 6 options

IF (DNS) THEN
   FLUX_LIMITER  = 4
   CHECK_VN = .TRUE.
   VN_MIN = 0.4_EB
   VN_MAX = 0.5_EB
   IF (TURBULENCE_MODEL/='null') THEN
      WRITE(MESSAGE,'(A,A,A)')  'ERROR: TURBULENCE_MODEL, ',TRIM(TURBULENCE_MODEL),', is not appropriate for DNS.'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ELSE
   FLUX_LIMITER  = 2
   TURBULENCE_MODEL='DEARDORFF'
ENDIF

! Re-read the line to pick up any user-specified options

REWIND(LU_INPUT)
CALL CHECKREAD('MISC',LU_INPUT,IOS)
IF (IOS==0) READ(LU_INPUT,MISC)
REWIND(LU_INPUT)

! Temperature conversions

H0    = 0.5_EB*(U0**2+V0**2+W0**2)
TMPA  = TMPA + TMPM
TMPA4 = TMPA**4
 
! Miscellaneous
 
HCH    = C_HORIZONTAL
HCV    = C_VERTICAL
ASSUMED_GAS_TEMPERATURE = ASSUMED_GAS_TEMPERATURE + TMPM
TEX_ORI = TEXTURE_ORIGIN
FVEC = FORCE_VECTOR
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
 
IF (FLUX_LIMITER<0 .OR. FLUX_LIMITER>5) THEN
   WRITE(MESSAGE,'(A)')  'ERROR on MISC: Permissible values for FLUX_LIMITER=0:5'
   CALL SHUTDOWN(MESSAGE)
ENDIF

! Turbulence model

SELECT CASE (TRIM(TURBULENCE_MODEL))
   CASE ('CONSTANT SMAGORINSKY')
      TURB_MODEL=CONSMAG
   CASE ('DYNAMIC SMAGORINSKY')
      TURB_MODEL=DYNSMAG
   CASE ('DEARDORFF')
      TURB_MODEL=DEARDORFF
   CASE ('VREMAN')
      TURB_MODEL=VREMAN
   CASE ('null')
      TURB_MODEL=0
   CASE DEFAULT
      WRITE(MESSAGE,'(A,A,A)')  'ERROR: TURBULENCE_MODEL, ',TRIM(TURBULENCE_MODEL),', is not recognized.'
      CALL SHUTDOWN(MESSAGE)
END SELECT

IF (AL2O3) RHO_SOOT = 4000._EB

! Level set based model of firespread in vegetation

IF(VEG_LEVEL_SET) WIND_ONLY = .TRUE.

IF (HRRPUV_MAX_SMV<0.0) THEN
   HRRPUV_MAX_SMV=1200.0
   USE_HRRPUV_MAX_SMV=0
ELSE
   USE_HRRPUV_MAX_SMV=1
ENDIF

! Set the lower limit of the extinction coefficient

EC_LL = VISIBILITY_FACTOR/MAXIMUM_VISIBILITY

END SUBROUTINE READ_MISC



SUBROUTINE READ_DUMP

! Read parameters associated with output files
 
REAL(EB) :: DT_DEFAULT
INTEGER :: N,I_DUM,SIG_FIGS,SIG_FIGS_EXP
NAMELIST /DUMP/ COLUMN_DUMP_LIMIT,CTRL_COLUMN_LIMIT,CUTCELL_DATA_FILE,&
                DEBUG,DEVC_COLUMN_LIMIT,DT_BNDE,DT_BNDF,DT_CTRL,DT_DEVC,DT_DEVC_LINE,DT_FLUSH,&
                DT_GEOM,DT_HRR,DT_ISOF,DT_MASS,DT_PART,DT_PL3D,DT_PROF,DT_RESTART,DT_SL3D,DT_SLCF,DT_VEG,&
                EB_PART_FILE,FLUSH_FILE_BUFFERS,&
                MASS_FILE,MAXIMUM_PARTICLES,NFRAMES,PLOT3D_PART_ID,PLOT3D_QUANTITY,PLOT3D_SPEC_ID,PLOT3D_VELO_INDEX,&
                RENDER_FILE,SIG_FIGS,SIG_FIGS_EXP,SMOKE3D,SMOKE3D_QUANTITY,SMOKE3D_SPEC_ID,STATUS_FILES,TIMING,&
                UVW_TIMER,VELOCITY_ERROR_FILE,WRITE_XYZ
 
! Set defaults

DEBUG              = .FALSE.
FLUSH_FILE_BUFFERS = .TRUE.
MAXIMUM_PARTICLES   = 500000
NFRAMES            = 1000 
PLOT3D_QUANTITY(1) = 'TEMPERATURE'
PLOT3D_QUANTITY(2) = 'U-VELOCITY'
PLOT3D_QUANTITY(3) = 'V-VELOCITY'
PLOT3D_QUANTITY(4) = 'W-VELOCITY'
PLOT3D_QUANTITY(5) = 'HRRPUV'
PLOT3D_PART_ID     = 'null'
PLOT3D_SPEC_ID     = 'null'
PLOT3D_VELO_INDEX  = 0
RENDER_FILE        = 'null'
SMOKE3D            = .TRUE.
SMOKE3D_QUANTITY   = 'null'
SMOKE3D_SPEC_ID    = 'null'
SIG_FIGS           = 8
SIG_FIGS_EXP       = 3
TIMING             = .FALSE.
UVW_TIMER          = 1.E10_EB
 
DT_GEOM      = 1.E10_EB
DT_BNDE      = -1._EB
DT_BNDF      = -1._EB
DT_RESTART   = 1000000._EB
DT_FLUSH     = -1._EB
DT_DEVC      = -1._EB
DT_DEVC_LINE = 0.25_EB*(T_END-T_BEGIN)
DT_HRR       = -1._EB
DT_ISOF      = -1._EB
DT_MASS      = -1._EB
DT_PART      = -1._EB
DT_PL3D      =  1.E10_EB
DT_PROF      = -1._EB
DT_SLCF      = -1._EB
DT_SL3D      = 0.2_EB*(T_END-T_BEGIN)
DT_CTRL      = -1._EB
DT_VEG       = -1._EB
 
! Read the DUMP line

REWIND(LU_INPUT)
DUMP_LOOP: DO 
   CALL CHECKREAD('DUMP',LU_INPUT,IOS)
   IF (IOS==1) EXIT DUMP_LOOP
   READ(LU_INPUT,DUMP,END=23,ERR=24,IOSTAT=IOS)
   24 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with DUMP line')
ENDDO DUMP_LOOP
23 REWIND(LU_INPUT)

! Set output time intervals

DT_DEFAULT = (T_END - T_BEGIN)/REAL(NFRAMES,EB)

IF (DT_BNDE < 0._EB) THEN ; DT_BNDE = 2._EB*DT_DEFAULT ; ELSE ; DT_BNDE = DT_BNDE /TIME_SHRINK_FACTOR ; ENDIF
IF (DT_BNDF < 0._EB) THEN ; DT_BNDF = 2._EB*DT_DEFAULT ; ELSE ; DT_BNDF = DT_BNDF /TIME_SHRINK_FACTOR ; ENDIF
IF (DT_DEVC < 0._EB) THEN ; DT_DEVC = DT_DEFAULT       ; ELSE ; DT_DEVC = DT_DEVC /TIME_SHRINK_FACTOR ; ENDIF
IF (DT_HRR  < 0._EB) THEN ; DT_HRR  = DT_DEFAULT       ; ELSE ; DT_HRR  = DT_HRR  /TIME_SHRINK_FACTOR ; ENDIF
IF (DT_ISOF < 0._EB) THEN ; DT_ISOF = DT_DEFAULT       ; ELSE ; DT_ISOF = DT_ISOF /TIME_SHRINK_FACTOR ; ENDIF
IF (DT_MASS < 0._EB) THEN ; DT_MASS = DT_DEFAULT       ; ELSE ; DT_MASS = DT_MASS /TIME_SHRINK_FACTOR ; ENDIF
IF (DT_PART < 0._EB) THEN ; DT_PART = DT_DEFAULT       ; ELSE ; DT_PART = DT_PART /TIME_SHRINK_FACTOR ; ENDIF
IF (DT_PROF < 0._EB) THEN ; DT_PROF = DT_DEFAULT       ; ELSE ; DT_PROF = DT_PROF /TIME_SHRINK_FACTOR ; ENDIF
IF (DT_SLCF < 0._EB) THEN ; DT_SLCF = DT_DEFAULT       ; ELSE ; DT_SLCF = DT_SLCF /TIME_SHRINK_FACTOR ; ENDIF
IF (DT_CTRL < 0._EB) THEN ; DT_CTRL = DT_DEFAULT       ; ELSE ; DT_CTRL = DT_CTRL /TIME_SHRINK_FACTOR ; ENDIF
IF (DT_FLUSH< 0._EB) THEN ; DT_FLUSH= DT_DEFAULT       ; ELSE ; DT_FLUSH= DT_FLUSH/TIME_SHRINK_FACTOR ; ENDIF
IF (DT_VEG  < 0._EB) THEN ; DT_VEG  = DT_DEFAULT       ; ELSE ; DT_VEG  = DT_VEG  /TIME_SHRINK_FACTOR ; ENDIF

! Check Plot3D QUANTITIES

PLOOP: DO N=1,5
   CALL GET_QUANTITY_INDEX(PLOT3D_SMOKEVIEW_LABEL(N),PLOT3D_SMOKEVIEW_BAR_LABEL(N),PLOT3D_QUANTITY_INDEX(N),I_DUM, &
                           PLOT3D_Y_INDEX(N),PLOT3D_Z_INDEX(N),PLOT3D_PART_INDEX(N),I_DUM,I_DUM,'PLOT3D', &
                           PLOT3D_QUANTITY(N),'null',PLOT3D_SPEC_ID(N),PLOT3D_PART_ID(N),'null','null')
ENDDO PLOOP

! Check SMOKE3D viability

IF (TWO_D .OR. SOLID_PHASE_ONLY) SMOKE3D = .FALSE.

IF (SMOKE3D_QUANTITY=='null') THEN
   IF (SOOT_INDEX > 0)  THEN
      SMOKE3D_QUANTITY = 'MASS FRACTION'
      SMOKE3D_SPEC_ID  = 'SOOT'
   ELSE
      IF (N_REACTIONS > 0)  THEN
         SMOKE3D_QUANTITY = 'HRRPUV'
      ELSE
         SMOKE3D = .FALSE.
      ENDIF      
   ENDIF
ENDIF

IF (SMOKE3D) THEN
   CALL GET_QUANTITY_INDEX(SMOKE3D_SMOKEVIEW_LABEL,SMOKE3D_SMOKEVIEW_BAR_LABEL,SMOKE3D_QUANTITY_INDEX,I_DUM, &
                           SMOKE3D_Y_INDEX,SMOKE3D_Z_INDEX,I_DUM,I_DUM,I_DUM,'SMOKE3D', &
                           SMOKE3D_QUANTITY,'null',SMOKE3D_SPEC_ID,'null','null','null')
ENDIF

! Set format of real number output

WRITE(FMT_R,'(A,I2.2,A,I2.2,A,I1.1)') 'ES',SIG_FIGS+SIG_FIGS_EXP+4,'.',SIG_FIGS-1,'E',SIG_FIGS_EXP

END SUBROUTINE READ_DUMP

 

SUBROUTINE READ_SPEC

USE PHYSICAL_FUNCTIONS, ONLY : AMBIENT_WATER_VAPOR
USE MATH_FUNCTIONS, ONLY : GET_RAMP_INDEX
USE PROPERTY_DATA, ONLY: GAS_PROPS,FED_PROPS
REAL(EB) :: MASS_FRACTION_0,MW,SIGMALJ,EPSILONKLJ,VISCOSITY,CONDUCTIVITY,DIFFUSIVITY,MASS_EXTINCTION_COEFFICIENT, &
            SPECIFIC_HEAT,REFERENCE_ENTHALPY,REFERENCE_TEMPERATURE,FIC_CONCENTRATION,FLD_LETHAL_DOSE,HUMIDITY,&
            SPECIFIC_HEAT_LIQUID,DENSITY_LIQUID,VAPORIZATION_TEMPERATURE,HEAT_OF_VAPORIZATION,MELTING_TEMPERATURE,&
            H_V_REFERENCE_TEMPERATURE,MEAN_DIAMETER,CONDUCTIVITY_SOLID,DENSITY_SOLID
INTEGER  :: N_SPEC_READ,N,NN,NR
LOGICAL  :: ABSORBING,SMIX_COMPONENT_ONLY,BACKGROUND,PREDEFINED(MAX_SPECIES),AEROSOL
CHARACTER(30) :: PREDEFINED_SPEC_ID(MAX_SPECIES),RAMP_CP,RAMP_CP_L,RAMP_K,RAMP_MU,RAMP_D
CHARACTER(100) :: FORMULA
TYPE(SPECIES_TYPE), POINTER :: SS=>NULL()
NAMELIST /SPEC/ ABSORBING,AEROSOL,BACKGROUND,CONDUCTIVITY,CONDUCTIVITY_SOLID,DENSITY_LIQUID,DENSITY_SOLID,DIFFUSIVITY,EPSILONKLJ,&
                FIC_CONCENTRATION,FLD_LETHAL_DOSE, &
                FORMULA,FYI,HEAT_OF_VAPORIZATION,HUMIDITY,H_V_REFERENCE_TEMPERATURE,ID,MASS_EXTINCTION_COEFFICIENT,&
                MASS_FRACTION_0,MEAN_DIAMETER,MELTING_TEMPERATURE,MW,RAMP_CP,RAMP_CP_L,RAMP_D,RAMP_K,RAMP_MU,REFERENCE_ENTHALPY,&
                REFERENCE_TEMPERATURE,SIGMALJ,SMIX_COMPONENT_ONLY,SPECIFIC_HEAT,SPECIFIC_HEAT_LIQUID,VAPORIZATION_TEMPERATURE,&
                VISCOSITY

! Look through the input file to see if the word 'BACKGROUND' appears anywhere

CALL SEARCH_INPUT_FILE('BACKGROUND',LU_INPUT,BACKGROUND_DECLARED)

! Create predefined inputs related to simple chemistry mode

PREDEFINED = .FALSE.

IF (SIMPLE_CHEMISTRY) THEN
   N_SPECIES = 7
   PREDEFINED(1:7)       = .TRUE.
   PREDEFINED_SPEC_ID(1) = REACTION(1)%FUEL
   PREDEFINED_SPEC_ID(2) = 'NITROGEN'
   PREDEFINED_SPEC_ID(3) = 'OXYGEN'
   PREDEFINED_SPEC_ID(4) = 'CARBON DIOXIDE'
   PREDEFINED_SPEC_ID(5) = 'CARBON MONOXIDE'
   PREDEFINED_SPEC_ID(6) = 'WATER VAPOR'
   PREDEFINED_SPEC_ID(7) = 'SOOT'
ELSE
   N_SPECIES = 0
ENDIF

N_SPEC_READ = 0

! Count SPEC lines and check for errors

REWIND(LU_INPUT)
COUNT_SPEC_LOOP: DO
   CALL SET_SPEC_DEFAULT
   CALL CHECKREAD('SPEC',LU_INPUT,IOS) 
   IF (IOS==1) EXIT COUNT_SPEC_LOOP
   READ(LU_INPUT,NML=SPEC,END=29,ERR=30,IOSTAT=IOS)
   N_SPEC_READ   = N_SPEC_READ   + 1   
   IF (ID=='null') THEN
      WRITE(MESSAGE,'(A,I2,A)') 'ERROR: Species ',N_SPEC_READ+1, ' needs a name (ID=...)'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   
   N_SPECIES = N_SPECIES + 1
   
   ! If simple chemistry mode, check to see if the species has already been defined.

   IF (SIMPLE_CHEMISTRY) THEN
      DO NN=1,N_SPECIES-1
         IF (TRIM(PREDEFINED_SPEC_ID(NN))==TRIM(ID)) THEN
            PREDEFINED(NN) = .FALSE.
            N_SPECIES = N_SPECIES - 1
            EXIT
         ENDIF
      ENDDO
   ENDIF
 
   IF (BACKGROUND) THEN 
      IF (N_SPEC_READ > 1) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: BACKGROUND species must be defined as the first SPEC input'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      IF (SIMPLE_CHEMISTRY) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Cannot define a BACKGROUND species if using the simple chemistry'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      IF (SMIX_COMPONENT_ONLY) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Cannot define a SMIX_COMPONENT_ONLY species as the BACKGROUND species'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      BACKGROUND_SPEC_INDEX = N_SPEC_READ
   ENDIF
   IF (SMIX_COMPONENT_ONLY .AND. (MASS_FRACTION_0>0._EB.OR.HUMIDITY>0._EB)) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Cannot define MASS_FRACTION_0 or HUMIDITY for a SMIX_COMPONENT_ONLY species'
         CALL SHUTDOWN(MESSAGE)
   ENDIF
      
   30 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I2)') 'ERROR: Problem with SPECies number',N_SPECIES+1
         CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_SPEC_LOOP
29 REWIND(LU_INPUT)      

! If no background species has been declared, assume it is AIR.
IF (.NOT. SIMPLE_CHEMISTRY .AND. BACKGROUND_SPEC_INDEX<1 .AND. .NOT.BACKGROUND_DECLARED) THEN
   N_SPECIES = N_SPECIES + 1
   PREDEFINED(1)         = .TRUE.
   PREDEFINED_SPEC_ID(1) = 'AIR'
   BACKGROUND_SPEC_INDEX = 1
ENDIF

! Allocate the primitive species array.

ALLOCATE(SPECIES(N_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ','SPECIES',IZERO) 

! Read and process species

SPEC_READ_LOOP: DO N=1,N_SPECIES

   CALL SET_SPEC_DEFAULT

   IF (PREDEFINED(N)) THEN
      ID = PREDEFINED_SPEC_ID(N)
      SMIX_COMPONENT_ONLY = .TRUE.
      IF (N==BACKGROUND_SPEC_INDEX) BACKGROUND = .TRUE.
   ELSE
      READ(LU_INPUT,NML=SPEC)
   ENDIF

   SS => SPECIES(N)

   SS%ABSORBING                   = ABSORBING
   SS%BACKGROUND                  = BACKGROUND
   SS%K_USER                      = CONDUCTIVITY
   SS%CONDUCTIVITY_SOLID          = CONDUCTIVITY_SOLID 
   SS%D_USER                      = DIFFUSIVITY
   SS%DENSITY_SOLID               = DENSITY_SOLID
   SS%EPSK                        = EPSILONKLJ
   SS%FIC_CONCENTRATION           = FIC_CONCENTRATION
   SS%FLD_LETHAL_DOSE             = FLD_LETHAL_DOSE
   SS%FORMULA                     = FORMULA
   SS%ID                          = ID
   SS%MASS_EXTINCTION_COEFFICIENT = MAX(0._EB,MASS_EXTINCTION_COEFFICIENT)
   SS%MEAN_DIAMETER               = MEAN_DIAMETER
   SS%MU_USER                     = VISCOSITY
   SS%MW                          = MW
   SS%RAMP_CP                     = RAMP_CP
   SS%RAMP_CP_L                   = RAMP_CP_L
   SS%RAMP_D                      = RAMP_D
   SS%RAMP_K                      = RAMP_K
   SS%RAMP_MU                     = RAMP_MU
   SS%REFERENCE_TEMPERATURE       = REFERENCE_TEMPERATURE + TMPM
   SS%SIG                         = SIGMALJ
   SS%SMIX_COMPONENT_ONLY         = SMIX_COMPONENT_ONLY   
   SS%SPECIFIC_HEAT               = SPECIFIC_HEAT*1000._EB
   SS%REFERENCE_ENTHALPY          = REFERENCE_ENTHALPY*1000._EB
   SS%YY0                         = MAX(0._EB,MASS_FRACTION_0)

   SS%DENSITY_LIQUID              = DENSITY_LIQUID

   IF ((HEAT_OF_VAPORIZATION >  0._EB .AND. SPECIFIC_HEAT_LIQUID <= 0._EB) .OR. &
       (HEAT_OF_VAPORIZATION <= 0._EB .AND. SPECIFIC_HEAT_LIQUID >  0._EB)) THEN
      WRITE(MESSAGE,'(A,I4,A)') 'ERROR: SPEC ' ,N,' If one of SPECIFIC_HEAT_LIQUID or HEAT_OF_VAPORIZATION defined, both must be'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (SPECIFIC_HEAT_LIQUID > 0._EB) THEN
      IF (MELTING_TEMPERATURE < -TMPM) THEN
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: PART ' ,N,' MELTING_TEMPERATURE not set'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      SS%SPECIFIC_HEAT_LIQUID        = SPECIFIC_HEAT_LIQUID*1000._EB
      SS%HEAT_OF_VAPORIZATION        = HEAT_OF_VAPORIZATION*1000._EB
      SS%TMP_MELT = MELTING_TEMPERATURE + TMPM
      IF (H_V_REFERENCE_TEMPERATURE < -TMPM) H_V_REFERENCE_TEMPERATURE = MELTING_TEMPERATURE
      SS%H_V_REFERENCE_TEMPERATURE = H_V_REFERENCE_TEMPERATURE + 273.15_EB
      IF (VAPORIZATION_TEMPERATURE< -TMPM) THEN
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: PART ' ,N,' VAPORIZATION_TEMPERATURE not set'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      SS%TMP_V = VAPORIZATION_TEMPERATURE + TMPM
   ENDIF

   IF (SS%SPECIFIC_HEAT < 0._EB .AND. SS%REFERENCE_ENTHALPY > -1.E20_EB .AND. RAMP_CP=='null') THEN
      WRITE(MESSAGE,'(A,I2,A)') 'ERROR: SPECies ',N,': Specify both SPECIFIC_HEAT or RAMP_CP and REFERENCE_ENTHALPY'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (SS%SPECIFIC_HEAT > 0._EB .AND. SS%REFERENCE_ENTHALPY < -1.E20_EB) THEN
      SS%REFERENCE_ENTHALPY = SS%SPECIFIC_HEAT * SS%REFERENCE_TEMPERATURE
   ENDIF
   IF (SS%RAMP_CP/='null' .AND. SS%REFERENCE_ENTHALPY < -1.E20_EB) THEN
      WRITE(MESSAGE,'(A,I2,A)') 'ERROR: SPECies ',N,': Specify both RAMP_CP and REFERENCE_ENTHALPY'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   CALL GAS_PROPS(SS%ID,SS%SIG,SS%EPSK,SS%MW,SS%ABSORBING,SS%FORMULA,SS%LISTED,SS%ATOMS)  
   CALL FED_PROPS(SS%ID,SS%FLD_LETHAL_DOSE,SS%FIC_CONCENTRATION)         

   IF (TRIM(SS%FORMULA)=='null') WRITE(SS%FORMULA,'(A,I0)') 'SPEC_',N

   ! For simple chemistry Determine if the species is the one specified on the REAC line(s)

   IF (SIMPLE_CHEMISTRY) THEN
      IF (TRIM(ID)==TRIM(REACTION(1)%FUEL)) THEN
         FUEL_INDEX = N
         SS%MW = REACTION(1)%MW_FUEL
      ENDIF
      IF (TRIM(ID)=='SOOT') SS%MW = REACTION(1)%MW_SOOT
   ENDIF

   SS%RCON = R0/SS%MW
   SS%MODE = GAS_SPECIES

   ! Special processing of certain species

   SELECT CASE (ID)
      CASE('WATER VAPOR')
         H2O_INDEX = N
         IF (MASS_FRACTION_0<0._EB) THEN           
            IF (HUMIDITY < 0._EB) HUMIDITY = 40._EB
            SS%YY0 = AMBIENT_WATER_VAPOR(HUMIDITY,TMPA)
         ENDIF
      CASE('CARBON DIOXIDE')
         CO2_INDEX = N
      CASE('CARBON MONOXIDE')
         CO_INDEX = N
      CASE('OXYGEN')
         O2_INDEX = N
      CASE('NITROGEN')
         N2_INDEX = N
      CASE('HYDROGEN')
         H2_INDEX = N
      CASE('HYDROGEN CYANIDE')
         HCN_INDEX = N
      CASE('NITRIC OXIDE')
         NO_INDEX = N
      CASE('NITROGEN DIOXIDE')
         NO2_INDEX = N
      CASE('SOOT')
         SOOT_INDEX = N
         IF (MASS_EXTINCTION_COEFFICIENT < 0._EB) SS%MASS_EXTINCTION_COEFFICIENT = 8700._EB
         IF (.NOT. SIMPLE_CHEMISTRY .AND. TRIM(SS%FORMULA)=='Soot') SS%ATOMS(6) = 1._EB
   END SELECT

   IF (AEROSOL) SS%MODE = AEROSOL_SPECIES

   ! Get ramps
   IF (SS%RAMP_CP/='null') THEN
      CALL GET_RAMP_INDEX(SS%RAMP_CP,'TEMPERATURE',NR)
      SS%RAMP_CP_INDEX = NR
   ENDIF
   IF (SS%RAMP_CP_L/='null') THEN
      CALL GET_RAMP_INDEX(SS%RAMP_CP_L,'TEMPERATURE',NR)
      SS%RAMP_CP_L_INDEX = NR
   ENDIF
   IF (SS%RAMP_D/='null') THEN
      CALL GET_RAMP_INDEX(SS%RAMP_D,'TEMPERATURE',NR)
      SS%RAMP_D_INDEX = NR
   ENDIF
   IF (SS%RAMP_K/='null') THEN
      CALL GET_RAMP_INDEX(SS%RAMP_K,'TEMPERATURE',NR)
      SS%RAMP_K_INDEX = NR
   ENDIF
   IF (SS%RAMP_MU/='null') THEN
      CALL GET_RAMP_INDEX(SS%RAMP_MU,'TEMPERATURE',NR)
      SS%RAMP_MU_INDEX = NR
   ENDIF
   

ENDDO SPEC_READ_LOOP      

CONTAINS


SUBROUTINE SET_SPEC_DEFAULT

ABSORBING                   = .FALSE.
AEROSOL                    = .FALSE.
BACKGROUND                  = .FALSE.
CONDUCTIVITY                = -1._EB
CONDUCTIVITY_SOLID          = 0.26 ! Ben-Dor, et al. 2002. (~10 x air)
DENSITY_SOLID               = 1800._EB !Slowik, et al. 2004
DIFFUSIVITY                 = -1._EB
EPSILONKLJ                  =  0._EB
FIC_CONCENTRATION           =  0._EB
FLD_LETHAL_DOSE             =  0._EB
FORMULA                     = 'null'
FYI                         = 'null'
HUMIDITY                    = -1._EB
ID                          = 'null'
SMIX_COMPONENT_ONLY         = .FALSE.
MASS_EXTINCTION_COEFFICIENT = -1._EB  ! m2/kg
MASS_FRACTION_0             = -1._EB
MEAN_DIAMETER               =  1.E-6_EB
MW                          =  0._EB
REFERENCE_TEMPERATURE       = 25._EB
SIGMALJ                     =  0._EB
SPECIFIC_HEAT               = -1._EB
REFERENCE_ENTHALPY           = -2.E20_EB
VISCOSITY                   = -1._EB

DENSITY_LIQUID              = -1._EB
HEAT_OF_VAPORIZATION        = -1._EB     ! kJ/kg
H_V_REFERENCE_TEMPERATURE   = -300._EB
MELTING_TEMPERATURE         = -300.        ! C
SPECIFIC_HEAT_LIQUID        = -1._EB     ! kJ/kg-K
VAPORIZATION_TEMPERATURE    = -300._EB     ! C

RAMP_CP                     = 'null'
RAMP_CP_L                   = 'null'
RAMP_D                      = 'null'
RAMP_K                      = 'null'
RAMP_MU                     = 'null'

END SUBROUTINE SET_SPEC_DEFAULT

END SUBROUTINE READ_SPEC



SUBROUTINE READ_SMIX

CHARACTER(30) :: SPEC_ID(1:MAX_SPECIES),ID
CHARACTER(100) :: FORMULA
REAL(EB) :: VOLUME_FRACTION(MAX_SPECIES),MASS_FRACTION(MAX_SPECIES),CONVERSION,MASS_FRACTION_0
INTEGER :: N,NS,NS2,N_SUB_SPECIES,N_SMIX_READ,SPEC_INDEX(0:MAX_SPECIES)=-1,Y_INDEX(N_SPECIES)
LOGICAL :: BACKGROUND,PREDEFINED(0:MAX_SPECIES)
TYPE(SPECIES_MIXTURE_TYPE), POINTER :: SM=>NULL()
NAMELIST /SMIX/ BACKGROUND,ID,MASS_FRACTION,MASS_FRACTION_0,SPEC_ID,VOLUME_FRACTION

! Set up an array to hold information about explicitly and implicitly defined species mixtures

PREDEFINED = .FALSE.
N_TRACKED_SPECIES = 0
N_SMIX_READ       = 0

IF (SIMPLE_CHEMISTRY) THEN
    PREDEFINED(0:2) = .TRUE.
    N_TRACKED_SPECIES = 2
ENDIF

! Count the explicitly defined species mixtures

REWIND(LU_INPUT)

COUNT_SMIX_LOOP: DO
   BACKGROUND = .FALSE.
   ID = 'null'
   CALL CHECKREAD('SMIX',LU_INPUT,IOS) 
   IF (IOS==1) EXIT COUNT_SMIX_LOOP
   READ(LU_INPUT,NML=SMIX,END=29,ERR=30,IOSTAT=IOS)
   N_SMIX_READ       = N_SMIX_READ + 1
   N_TRACKED_SPECIES = N_TRACKED_SPECIES + 1
   IF (BACKGROUND) THEN
      IF (N_SMIX_READ > 1) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Background mixture must be defined first'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      IF (SIMPLE_CHEMISTRY) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Cannot define a background SMIX with simple chemistry'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      BACKGROUND_SMIX_INDEX = N_SMIX_READ
      N_TRACKED_SPECIES = N_TRACKED_SPECIES - 1
   ENDIF
   
   30 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I2)') 'ERROR: Problem with Species MIXture number ',N_SMIX_READ+1
         CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_SMIX_LOOP
29 REWIND(LU_INPUT)

! If there is no BACKGROUND species or species mixture, stop.

IF (.NOT.SIMPLE_CHEMISTRY .AND. BACKGROUND_SPEC_INDEX<0 .AND. BACKGROUND_SMIX_INDEX<0) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: Declare one SPEC or SMIX to be the BACKGROUND'
   CALL SHUTDOWN(MESSAGE)
ENDIF

! Look through the primitive species and create SMIX entries either for a background or a non-component_only.

DO N=1,N_SPECIES
   IF (.NOT.SPECIES(N)%SMIX_COMPONENT_ONLY .OR. SPECIES(N)%BACKGROUND) THEN
      IF (SPECIES(N)%BACKGROUND) THEN
         IF (BACKGROUND_SMIX_INDEX>=0) CYCLE
         IF (SIMPLE_CHEMISTRY)         CYCLE
         PREDEFINED(0) = .TRUE.
         SPEC_INDEX(0) = N
      ELSE
         N_TRACKED_SPECIES = N_TRACKED_SPECIES + 1
         PREDEFINED(N_TRACKED_SPECIES) = .TRUE.
         SPEC_INDEX(N_TRACKED_SPECIES) = N
      ENDIF
   ENDIF
ENDDO

! Allocate the SPECIES_MIXTURE array

ALLOCATE(SPECIES_MIXTURE(0:N_TRACKED_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ_SMIX','SPECIES_MIXTURE',IZERO)

! Read and/or process all of the species mixtures

READ_SMIX_LOOP: DO N=0,N_TRACKED_SPECIES  

   SM => SPECIES_MIXTURE(N)

   WRITE(SM%FORMULA,'(A,I2)') 'Z',N   
   MASS_FRACTION   = 0._EB
   MASS_FRACTION_0 = 0._EB
   VOLUME_FRACTION = 0._EB
   BACKGROUND      = .FALSE.
   CONVERSION      = 0._EB
   SPEC_ID         = 'null'

   IF (PREDEFINED(N)) THEN
      CALL SETUP_PREDEFINED_SMIX(N)
   ELSE
      READ(LU_INPUT,NML=SMIX)
   ENDIF

   IF (TRIM(ID)=='null') THEN
      WRITE(MESSAGE,'(A,I2,A)') 'ERROR: SMIX ',N,' needs an ID'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   IF (ANY(MASS_FRACTION>0._EB) .AND. ANY(VOLUME_FRACTION>0._EB)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: Do not combine MASS_FRACTION and VOLUME_FRACTION on the same SMIX line'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   SM%ID = ID
   SM%ZZ0 = MAX(0._EB,MASS_FRACTION_0)
   
   ! Count the number of species included in the mixture

   N_SUB_SPECIES = 0
   COUNT_SPEC: DO NS=1,N_SPECIES
      IF (TRIM(SPEC_ID(NS)) /= 'null') THEN
         N_SUB_SPECIES = N_SUB_SPECIES + 1
      ELSE
         EXIT
      ENDIF
   ENDDO COUNT_SPEC

   ! Allocate arrays to store the species id, mass, volume fractions

   ALLOCATE (SM%SPEC_ID(N_SPECIES),STAT=IZERO)
   ALLOCATE (SM%VOLUME_FRACTION(N_SPECIES),STAT=IZERO)
   ALLOCATE (SM%MASS_FRACTION(N_SPECIES),STAT=IZERO)
   SM%SPEC_ID         = 'null'
   SM%VOLUME_FRACTION = 0._EB
   SM%MASS_FRACTION   = 0._EB

   Y_INDEX = -1
   DO NS = 1,N_SUB_SPECIES
      FIND_SPEC_ID: DO NS2 = 1,N_SPECIES
         IF (TRIM(SPECIES(NS2)%ID) == TRIM(SPEC_ID(NS))) THEN
            SM%SPEC_ID(NS2) = SPEC_ID(NS)
            Y_INDEX(NS)  = NS2
            IF (N_SUB_SPECIES==1) THEN
               SM%FORMULA = SPECIES(NS2)%FORMULA
               IF (SPECIES(NS2)%MODE == AEROSOL_SPECIES) THEN
                  SM%DEPOSITING = .TRUE.
                  SM%MEAN_DIAMETER = SPECIES(NS2)%MEAN_DIAMETER
                  SM%DENSITY_SOLID = SPECIES(NS2)%DENSITY_SOLID
                  SM%CONDUCTIVITY_SOLID=SPECIES(NS2)%CONDUCTIVITY_SOLID
               ENDIF
            ENDIF
            EXIT FIND_SPEC_ID
         ENDIF
      ENDDO FIND_SPEC_ID
      IF (MASS_FRACTION(NS)>0._EB)     CONVERSION = CONVERSION + MASS_FRACTION(NS)   / SPECIES(Y_INDEX(NS))%MW
      IF (VOLUME_FRACTION(NS)>0._EB)   CONVERSION = CONVERSION + VOLUME_FRACTION(NS) * SPECIES(Y_INDEX(NS))%MW
   ENDDO

   IF (ANY(MASS_FRACTION>0._EB)) THEN
      DO NS = 1,N_SUB_SPECIES
         SM%VOLUME_FRACTION(Y_INDEX(NS)) = MASS_FRACTION(NS) / SPECIES(Y_INDEX(NS))%MW / CONVERSION
         SM%MASS_FRACTION(Y_INDEX(NS))   = MASS_FRACTION(NS)
      ENDDO
   ENDIF

   IF (ANY(VOLUME_FRACTION>0._EB)) THEN
      DO NS = 1,N_SUB_SPECIES
         SM%MASS_FRACTION(Y_INDEX(NS))   = VOLUME_FRACTION(NS) * SPECIES(Y_INDEX(NS))%MW / CONVERSION
         SM%VOLUME_FRACTION(Y_INDEX(NS)) = VOLUME_FRACTION(NS)
      ENDDO
   ENDIF

   ! Normalize mass and volume fractions, plus stoichiometric coefficient

   SM%MASS_FRACTION = SM%MASS_FRACTION / SUM(SM%MASS_FRACTION)
   SM%ADJUST_NU = SUM(SM%VOLUME_FRACTION)
   SM%VOLUME_FRACTION = SM%VOLUME_FRACTION / SUM(SM%VOLUME_FRACTION)

   ! Calculate the molecular weight and extinction coefficient for the mixture

   SM%MW = 0._EB
   DO NS = 1,N_SPECIES
      IF (SM%MASS_FRACTION(NS) <ZERO_P) CYCLE
      SM%MASS_EXTINCTION_COEFFICIENT = SM%MASS_FRACTION(NS) * SPECIES(NS)%MASS_EXTINCTION_COEFFICIENT      
      SM%MW = SM%MW + SM%VOLUME_FRACTION(NS) * SPECIES(NS)%MW !! *SM%ADJUST_NU ::term for potential non-normalized inputs
      IF (SPECIES(NS)%FORMULA(1:5)=='SPEC_') SM%VALID_ATOMS=.FALSE.
      SM%ATOMS = SM%ATOMS + SM%VOLUME_FRACTION(NS)*SPECIES(NS)%ATOMS !! *SM%ADJUST_NU ::term for potential non-normalized inputs
   ENDDO     
   SM%RCON = R0/SM%MW

ENDDO READ_SMIX_LOOP

IF (SPECIES_MIXTURE(0)%DEPOSITING) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: Cannot define the background as a depositing species'
   CALL SHUTDOWN(MESSAGE)
ENDIF

REWIND (LU_INPUT)   

! Normalize the initial mass fractions of the lumped species if necessary

IF (SUM(SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0) > 1._EB) &
  SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0 = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0/ &
                                         SUM(SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0)

SPECIES_MIXTURE(0)%ZZ0 = 1._EB - SUM(SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0)

CONTAINS


SUBROUTINE SETUP_PREDEFINED_SMIX(N)

! Set up the SMIX line either for the SIMPLE_CHEMISTRY mode or for a primitive species

INTEGER :: N
TYPE(REACTION_TYPE), POINTER :: RN

SIMPLE_CHEMISTRY_IF:IF (SIMPLE_CHEMISTRY) THEN
   MASS_FRACTION = 0._EB
   RN => REACTION(1)
   SELECT CASE(N)
      CASE(0)
         ID               = 'AIR'
         FORMULA          = 'Z0'
         SPEC_ID(1)       = 'WATER VAPOR'
         SPEC_ID(2)       = 'OXYGEN'
         SPEC_ID(3)       = 'CARBON DIOXIDE'
         SPEC_ID(4)       = 'NITROGEN'
         MASS_FRACTION(1) = Y_H2O_INFTY
         MASS_FRACTION(2) = Y_O2_INFTY*(1._EB-Y_H2O_INFTY)
         MASS_FRACTION(3) = Y_CO2_INFTY*(1._EB-Y_H2O_INFTY)
         MASS_FRACTION(4) = 1._EB-SUM(MASS_FRACTION)
      CASE(1)
         ID               = RN%FUEL
         FORMULA          = 'Z1'
         SPEC_ID(1)       = RN%FUEL
         MASS_FRACTION(1) = 1._EB  
      CASE(2)
         ID                 = 'PRODUCTS'
         FORMULA            = 'Z2'
         SPEC_ID(1)         = 'CARBON MONOXIDE'
         SPEC_ID(2)         = 'SOOT'
         SPEC_ID(3)         = 'WATER VAPOR'
         SPEC_ID(4)         = 'CARBON DIOXIDE'
         SPEC_ID(5)         = 'NITROGEN'
         RN%NU_CO           = (SPECIES(FUEL_INDEX)%MW/MW_CO)     *RN%CO_YIELD
         RN%NU_SOOT         = (SPECIES(FUEL_INDEX)%MW/RN%MW_SOOT)*RN%SOOT_YIELD
         RN%NU_H2O          = 0.5_EB*RN%H - 0.5_EB*RN%NU_SOOT*RN%SOOT_H_FRACTION
         RN%NU_CO2          = RN%C - RN%NU_CO - RN%NU_SOOT*(1._EB-RN%SOOT_H_FRACTION)
         RN%NU_O2           = RN%NU_CO2 + 0.5_EB*(RN%NU_CO+RN%NU_H2O-RN%O)
         RN%NU_N2           = RN%N*0.5_EB
         VOLUME_FRACTION(1) = RN%NU_CO
         VOLUME_FRACTION(2) = RN%NU_SOOT
         VOLUME_FRACTION(3) = RN%NU_H2O + SPECIES_MIXTURE(0)%VOLUME_FRACTION(H2O_INDEX)*RN%NU_O2 / &
                                          SPECIES_MIXTURE(0)%VOLUME_FRACTION(O2_INDEX)
         VOLUME_FRACTION(4) = RN%NU_CO2 + SPECIES_MIXTURE(0)%VOLUME_FRACTION(CO2_INDEX)*RN%NU_O2 / &
                                          SPECIES_MIXTURE(0)%VOLUME_FRACTION(O2_INDEX)
         VOLUME_FRACTION(5) = RN%NU_N2  + SPECIES_MIXTURE(0)%VOLUME_FRACTION(N2_INDEX)*RN%NU_O2 / &
                                          SPECIES_MIXTURE(0)%VOLUME_FRACTION(O2_INDEX)
         VOLUME_FRACTION    = VOLUME_FRACTION/SUM(VOLUME_FRACTION)
         SPECIES(SOOT_INDEX)%ATOMS=0._EB
         SPECIES(SOOT_INDEX)%ATOMS(1)=RN%SOOT_H_FRACTION
         SPECIES(SOOT_INDEX)%ATOMS(6)=1._EB-RN%SOOT_H_FRACTION         
      CASE DEFAULT
         VOLUME_FRACTION(1) = 1._EB
         MASS_FRACTION_0    = SPECIES(SPEC_INDEX(N))%YY0
         SPEC_ID(1)         = SPECIES(SPEC_INDEX(N))%ID
         ID                 = SPEC_ID(1)
   END SELECT
ELSE SIMPLE_CHEMISTRY_IF
   VOLUME_FRACTION(1) = 1._EB
   MASS_FRACTION_0    = SPECIES(SPEC_INDEX(N))%YY0
   SPEC_ID(1)         = SPECIES(SPEC_INDEX(N))%ID
   ID                 = SPEC_ID(1)  
ENDIF SIMPLE_CHEMISTRY_IF

END SUBROUTINE SETUP_PREDEFINED_SMIX

END SUBROUTINE READ_SMIX



SUBROUTINE PROC_SMIX

! Create the Z to Y transformation matrix and fill up the gas property tables

USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE PROPERTY_DATA, ONLY: JANAF_TABLE, CALC_GAS_PROPS, GAS_PROPS
REAL(EB), ALLOCATABLE, DIMENSION(:) :: MU_TMP,CP_TMP,K_TMP,H_TMP,D_TMP,ZZ_GET
REAL(EB) :: CP1,CP2,H1
INTEGER :: NN,N
TYPE(SPECIES_TYPE), POINTER :: SS=>NULL()
TYPE(SPECIES_MIXTURE_TYPE), POINTER :: SM=>NULL()

! Setup the array to convert the tracked species array to array of all primitive species

ALLOCATE(Z2Y(N_SPECIES,0:N_TRACKED_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ','Z2Y',IZERO) 
Z2Y = 0._EB 

DO N=0,N_TRACKED_SPECIES
   SM => SPECIES_MIXTURE(N)
   DO NN=1,N_SPECIES
      Z2Y(NN,N) = SM%MASS_FRACTION(NN)
   ENDDO
ENDDO

! Set up the arrays of molecular weights

ALLOCATE(MWR_Z(0:N_TRACKED_SPECIES),STAT=IZERO)
CALL ChkMemErr('READ','MW_AVG_Y',IZERO) 

MWR_Z(0:N_TRACKED_SPECIES) = 1._EB/SPECIES_MIXTURE(0:N_TRACKED_SPECIES)%MW

ALLOCATE(ZZ_GET(0:N_TRACKED_SPECIES))
ZZ_GET(0:N_TRACKED_SPECIES) = SPECIES_MIXTURE(0:N_TRACKED_SPECIES)%ZZ0
CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM0)
DEALLOCATE(ZZ_GET)

MW_MIN = MINVAL(SPECIES_MIXTURE%MW)
MW_MAX = MAXVAL(SPECIES_MIXTURE%MW)
 
! Compute background density from other background quantities
 
RHOA = P_INF/(TMPA*RSUM0)
 
! Compute constant-temperature specific heats
 
CP_GAMMA = SPECIES_MIXTURE(0)%RCON*GAMMA/(GAMMA-1._EB)
CPOPR = CP_GAMMA/PR

! Compute gas properties for primitive species 1 to N_SPECIES. 

ALLOCATE(D_TMP(N_SPECIES))
D_TMP = 0._EB
ALLOCATE(MU_TMP(N_SPECIES))
MU_TMP = 0._EB
ALLOCATE(CP_TMP(N_SPECIES))
CP_TMP = 0._EB
ALLOCATE(H_TMP(N_SPECIES))
H_TMP = 0._EB
ALLOCATE(K_TMP(N_SPECIES))
K_TMP = 0._EB

ALLOCATE(CPBAR_Z(0:5000,0:N_TRACKED_SPECIES))   
CALL ChkMemErr('READ','CPBAR_Z',IZERO)    
CPBAR_Z = 0._EB   

ALLOCATE(CP_AVG_Z(0:5000,0:N_TRACKED_SPECIES))   
CALL ChkMemErr('READ','CP_AVG_Z',IZERO)    
CP_AVG_Z = 0._EB      

ALLOCATE(K_Z(0:5000,0:N_TRACKED_SPECIES))
CALL ChkMemErr('READ','K_Z',IZERO)    
K_Z = 0._EB   

ALLOCATE(MU_Z(0:5000,0:N_TRACKED_SPECIES))      
CALL ChkMemErr('READ','MU_Z',IZERO)    
MU_Z = 0._EB   

ALLOCATE(CP_Z(0:5000,0:N_TRACKED_SPECIES))     
CALL ChkMemErr('READ','CP_Z',IZERO)    
CP_Z = 0._EB   

ALLOCATE(D_Z(0:5000,0:N_TRACKED_SPECIES))      
CALL ChkMemErr('READ','D_Z',IZERO)    
D_Z = 0._EB   

! Adjust reference enthalpy to 0 K if a RAMP_CP is given

DO N=1,N_SPECIES
   SS => SPECIES(N)
   IF(SS%RAMP_CP_INDEX > 0) THEN
      IF (NINT(SS%REFERENCE_TEMPERATURE)==0) CYCLE
      CP2 = EVALUATE_RAMP(1._EB,1._EB,SS%RAMP_CP_INDEX)*1000._EB
      H1 = 0._EB
      DO J=1,NINT(SS%REFERENCE_TEMPERATURE)
         CP1 = CP2
         CP2 = EVALUATE_RAMP(REAL(J,EB),1._EB,SS%RAMP_CP_INDEX)*1000._EB
         H1 = H1 + 0.5_EB*(CP1+CP2)
      ENDDO
      SS%REFERENCE_ENTHALPY = SS%REFERENCE_ENTHALPY - H1
   ENDIF
END DO


! Loop through temperatures from 1 K to 5000 K to get temperature-specific gas properties.  Data from JANAF 4

TABLE_LOOP: DO J=1,5000

   ! For each primitive species, get its property values at temperature J

   DO N=1,N_SPECIES
      SS => SPECIES(N)
      CALL CALC_GAS_PROPS(J,N,D_TMP(N),MU_TMP(N),K_TMP(N),CP_TMP(N),H_TMP(N),SS%ISFUEL)
      IF (SS%RAMP_CP_INDEX>0) THEN
         CP_TMP(N) = EVALUATE_RAMP(REAL(J,EB),0._EB,SS%RAMP_CP_INDEX)*1000._EB
         H_TMP(N) = SS%REFERENCE_ENTHALPY
      ENDIF
      IF (SS%RAMP_D_INDEX>0)  D_TMP(N) = EVALUATE_RAMP(REAL(J,EB),1._EB,SS%RAMP_D_INDEX)
      IF (SS%RAMP_K_INDEX>0)  K_TMP(N) = EVALUATE_RAMP(REAL(J,EB),1._EB,SS%RAMP_K_INDEX)/SS%MW
      IF (SS%RAMP_MU_INDEX>0) MU_TMP(N) = EVALUATE_RAMP(REAL(J,EB),1._EB,SS%RAMP_MU_INDEX)/SS%MW
   ENDDO

   ! For each tracked species, store the mass-weighted property values

   DO N=0,N_TRACKED_SPECIES
      D_Z(J,N)  = SUM(Z2Y(:,N) * D_TMP(:))
      CP_Z(J,N) = SUM(Z2Y(:,N) * CP_TMP(:))
      MU_Z(J,N) = SUM(Z2Y(:,N) * MU_TMP(:))
      K_Z(J,N)  = SUM(Z2Y(:,N) * K_TMP(:))
      IF (J==1) CP_Z(0,N) = CP_Z(1,N)
      CP_AVG_Z(J,N) = (CP_AVG_Z(J-1,N)*REAL(J-1,EB) + 0.5_EB*(CP_Z(J,N)+CP_Z(J-1,N)))/REAL(J,EB)
      IF (J>1) THEN      
         CPBAR_Z(J,N) = (CPBAR_Z(J-1,N)*(REAL(J,EB)-1._EB)+0.5_EB*(CP_Z(J,N)+CP_Z(J-1,N)))/REAL(J,EB)
      ELSE
         CPBAR_Z(0,N) = SUM(Z2Y(:,N) * H_TMP(:))
         CPBAR_Z(J,N) = CPBAR_Z(0,N) + CP_Z(J,N)
      ENDIF
   ENDDO
   
ENDDO TABLE_LOOP

DEALLOCATE(D_TMP)
DEALLOCATE(MU_TMP)
DEALLOCATE(CP_TMP)
DEALLOCATE(H_TMP)
DEALLOCATE(K_TMP)

END SUBROUTINE PROC_SMIX



SUBROUTINE READ_REAC

USE PHYSICAL_FUNCTIONS, ONLY : AMBIENT_WATER_VAPOR
USE PROPERTY_DATA, ONLY : ELEMENT,GET_FORMULA_WEIGHT,MAKE_PERIODIC_TABLE,SIMPLE_SPECIES_MW,GAS_PROPS
CHARACTER(30) :: FUEL
CHARACTER(100) :: FORMULA
CHARACTER(30), TARGET  :: SMIX_ID(MAX_SPECIES),SPEC_ID(MAX_SPECIES)
CHARACTER(255) :: EQUATION
CHARACTER(100) :: FWD_ID
INTEGER :: NR,NS,NS2
REAL(EB) :: SOOT_YIELD,CO_YIELD,EPUMO2,A,A_LEAN, &
            CRITICAL_FLAME_TEMPERATURE,HEAT_OF_COMBUSTION,NU(MAX_SPECIES),E,E_LEAN,N_S(MAX_SPECIES),C,H,N,O, &
            AUTO_IGNITION_TEMPERATURE,THRESHOLD_TEMPERATURE,HUMIDITY,SOOT_H_FRACTION,N_T
REAL(EB) :: E_TMP=0._EB,S_TMP=0._EB,ATOM_COUNTS(118),MW_FUEL=0._EB
LOGICAL :: A_TMP,L_TMP,CHECK_ATOM_BALANCE,EDCM,FAST_CHEMISTRY,PHI_DEPENDENCE,REVERSIBLE
NAMELIST /REAC/ A,A_LEAN,AUTO_IGNITION_TEMPERATURE,BETA_EDC,C,CHECK_ATOM_BALANCE,CO_PRODUCTION,CO_YIELD,CRITICAL_FLAME_TEMPERATURE,&
                E,E_LEAN,EDCM,EPUMO2,EQUATION,FIXED_MIX_TIME,FORMULA,FUEL,FWD_ID,FYI,H,HEAT_OF_COMBUSTION,HRRPUA_SHEET,&
                HRRPUV_AVERAGE,HUMIDITY,ID,IDEAL,N,NU,N_S,N_T,O,ODE_SOLVER,PHI_DEPENDENCE,REAC_ATOM_ERROR,REAC_MASS_ERROR,&
                REVERSIBLE,SMIX_ID,SOOT_H_FRACTION,SOOT_YIELD,SPEC_ID,SUPPRESSION,TAU_CHEM,TAU_FLAME,THRESHOLD_TEMPERATURE,&
                Y_CO2_INFTY,Y_O2_INFTY,Y_P_MIN_EDC
                
CALL MAKE_PERIODIC_TABLE
CALL SIMPLE_SPECIES_MW
ATOM_COUNTS = 0._EB
N_REACTIONS = 0
REWIND(LU_INPUT)
 
COUNT_REAC_LOOP: DO
   CALL CHECKREAD('REAC',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_REAC_LOOP
   CALL SET_REAC_DEFAULTS
   READ(LU_INPUT,REAC,END=435,ERR=434,IOSTAT=IOS)
   N_REACTIONS = N_REACTIONS + 1
   IF (A < 0._EB .AND. E < 0._EB .AND. TRIM(SMIX_ID(1))=='null' .AND. TRIM(EQUATION)=='null') THEN
      SIMPLE_CHEMISTRY = .TRUE.
      A = 1.E16_EB
      E = 0._EB
   ENDIF
   IF (.NOT.SIMPLE_CHEMISTRY .AND. TRIM(SMIX_ID(1))=='null' .AND. TRIM(EQUATION)=='null') THEN
      WRITE(MESSAGE,'(A,I3,A)') 'ERROR: Problem with REAC ',N_REACTIONS,'. SMIX_ID and NU arrays or EQUATION must be defined'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   434 IF (IOS>0) THEN
      WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with REAC ',N_REACTIONS+1
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO COUNT_REAC_LOOP

435 REWIND(LU_INPUT)

ALLOCATE(REACTION(N_REACTIONS),STAT=IZERO)

! Read and store the reaction parameters

REAC_READ_LOOP: DO NR=1,N_REACTIONS

   CALL SET_REAC_DEFAULTS
   READ(LU_INPUT,REAC)

   IF (FUEL=='null' .AND. ID/='null') FUEL = ID ! Backward compatibility

   IF (FUEL=='null' .AND. ID=='null') THEN
      WRITE(MESSAGE,'(A,I1,A)') 'ERROR: REAC ',NR,' requires a FUEL'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   RN => REACTION(NR)

   IF (SIMPLE_CHEMISTRY) THEN
      IF(C<=ZERO_P .AND. H<=ZERO_P) THEN
         IF (TRIM(FORMULA)=='null') THEN
            CALL GAS_PROPS(FUEL,S_TMP,E_TMP,MW_FUEL,A_TMP,FORMULA,L_TMP,ATOM_COUNTS)   
         ELSE         
            CALL GET_FORMULA_WEIGHT(FORMULA,MW_FUEL,ATOM_COUNTS)   
         ENDIF
         IF (ATOM_COUNTS(1)+ATOM_COUNTS(6)+ATOM_COUNTS(7)+ATOM_COUNTS(8) - SUM(ATOM_COUNTS) < 0._EB) THEN
            WRITE(MESSAGE,'(A)') 'ERROR: Fuel FORMULA for SIMPLE_CHEMISTRY can only contain C,H,O, and N'
            CALL SHUTDOWN(MESSAGE)
         ELSE
            C = ATOM_COUNTS(6)
            H = ATOM_COUNTS(1)
            O = ATOM_COUNTS(8)
            N = ATOM_COUNTS(7)
         ENDIF
         IF (C<=ZERO_P .AND. H<=ZERO_P) THEN
            WRITE(MESSAGE,'(A)') 'ERROR: Must specify fuel chemistry using C and/or H when using simple chemistry'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ELSE
         MW_FUEL = ELEMENT(6)%MASS*C+ELEMENT(1)%MASS*H+ELEMENT(8)%MASS*O+ELEMENT(7)%MASS*N   
      ENDIF
   ENDIF
   RN%A                         = A
   RN%A_IN                      = A
   RN%A_LEAN                    = A_LEAN
   RN%AUTO_IGNITION_TEMPERATURE = AUTO_IGNITION_TEMPERATURE + TMPM
   RN%C                         = C
   RN%CHECK_ATOM_BALANCE        = CHECK_ATOM_BALANCE
   RN%CO_YIELD                  = CO_YIELD
   RN%CRIT_FLAME_TMP            = CRITICAL_FLAME_TEMPERATURE + TMPM
   RN%E                         = E*1000._EB
   RN%E_IN                      = E
   RN%E_LEAN                    = E_LEAN*1000._EB
   RN%EDCM                      = EDCM
   RN%EQUATION                  = EQUATION
   RN%EPUMO2                    = EPUMO2*1000._EB
   RN%FAST_CHEMISTRY            = FAST_CHEMISTRY
   RN%FUEL                      = FUEL
   RN%FWD_ID                    = FWD_ID
   RN%FYI                       = FYI
   RN%H                         = H
   RN%HEAT_OF_COMBUSTION        = HEAT_OF_COMBUSTION*1000._EB
   RN%ID                        = ID
   RN%MW_FUEL                   = MW_FUEL
   RN%MW_SOOT                   = ELEMENT(6)%MASS * (1._EB-SOOT_H_FRACTION) + ELEMENT(1)%MASS*SOOT_H_FRACTION
   RN%N                         = N
   RN%N_T                       = N_T
   RN%O                         = O
   RN%PHI_DEPENDENCE            = PHI_DEPENDENCE
   RN%REVERSIBLE                = REVERSIBLE 
   RN%SOOT_H_FRACTION           = SOOT_H_FRACTION
   RN%SOOT_YIELD                = SOOT_YIELD
   RN%THRESHOLD_TEMP            = THRESHOLD_TEMPERATURE + TMPM

   ! Background water vapor. Only used for simple chemistry model.

   Y_H2O_INFTY = AMBIENT_WATER_VAPOR(HUMIDITY,TMPA)

   ! Determine the type of reaction

   IF (A > 0._EB) THEN 
     IF (.NOT. EDCM) THEN
        RN%MODE = FINITE_RATE
     ELSE
        RN%MODE = EDDY_DISSIPATION_CONCEPT
     ENDIF
   ELSE
     RN%MODE = EDDY_DISSIPATION
   ENDIF

   IF (A > 1.E15_EB .AND. E == 0._EB) RN%FAST_CHEMISTRY=.TRUE.
      
   ! Determine the number of stoichiometric coefficients for this reaction

   IF (.NOT.SIMPLE_CHEMISTRY) THEN
      NS2 = 0
      DO NS=1,MAX_SPECIES
         IF (TRIM(SMIX_ID(NS))/='null') THEN
            NS2=NS2+1
         ELSE
            EXIT
         ENDIF
      ENDDO
      RN%N_SMIX = NS2
      NS2 = 0
      IF(TRIM(RN%EQUATION)/='null') RN%N_SMIX = MAX_SPECIES
      DO NS=1,MAX_SPECIES
         IF (TRIM(SPEC_ID(NS))/='null') THEN
            NS2=NS2+1
         ELSE
            EXIT
         ENDIF
      ENDDO
      RN%N_SPEC = NS2
   ELSE
      RN%N_SMIX = 3
      RN%N_SPEC = 0
   ENDIF

   ! Store the "read in" values of N_S, NU, and SMIX_ID for use in PROC_REAC.
   IF (RN%N_SPEC > 0) THEN
      ALLOCATE(RN%N_S_READ(RN%N_SPEC))
      RN%N_S_READ(1:RN%N_SPEC) = N_S(1:RN%N_SPEC)
      ALLOCATE(RN%SPEC_ID_READ(RN%N_SPEC))
      RN%SPEC_ID_READ = 'null'
      RN%SPEC_ID_READ(1:RN%N_SPEC)=SPEC_ID(1:RN%N_SPEC)
   ENDIF

   ALLOCATE(RN%NU_READ(RN%N_SMIX))
   RN%NU_READ(1:RN%N_SMIX) = NU(1:RN%N_SMIX)

   ALLOCATE(RN%SMIX_ID_READ(RN%N_SMIX))
   RN%SMIX_ID_READ = 'null'   
   RN%SMIX_ID_READ(1:RN%N_SMIX)=SMIX_ID(1:RN%N_SMIX)   

END DO REAC_READ_LOOP

SIMP_CHEM: DO NR =1,N_REACTIONS
   RN => REACTION(NR)
   IF (.NOT. RN%FAST_CHEMISTRY) THEN
      SIMPLE_CHEM = .FALSE.
      EXIT SIMP_CHEM
   ENDIF
ENDDO SIMP_CHEM

REWIND(LU_INPUT)

CONTAINS

SUBROUTINE SET_REAC_DEFAULTS

AUTO_IGNITION_TEMPERATURE   = -TMPM
A                           = -1._EB     ! cm**3/mol-s
A_LEAN                      = -1._EB
C                           = 0._EB
CHECK_ATOM_BALANCE          = .TRUE.
CO_YIELD                    = 0._EB
CRITICAL_FLAME_TEMPERATURE  = 1427._EB   ! C
E                           = -1._EB     ! kJ/kmol
E_LEAN                      = -1._EB
EDCM                        = .FALSE.
EPUMO2                      = 13100._EB  ! kJ/kg
EQUATION                    = 'null'
FAST_CHEMISTRY              = .FALSE.
FORMULA                     = 'null'
FUEL                        = 'null'
FWD_ID                      = 'null'
FYI                         = 'null'
H                           = 0._EB
HEAT_OF_COMBUSTION          = -2.E20_EB
HRRPUV_AVERAGE              = 2500._EB
HRRPUA_SHEET                = 200._EB
HUMIDITY                    = 40._EB
ID                          = 'null'
N                           = 0._EB
NU                          = 0._EB
N_S                         = -999._EB
N_T                         = 0._EB
O                           = 0._EB
ODE_SOLVER                  = 'null'
PHI_DEPENDENCE              = .FALSE.
REAC_ATOM_ERROR             = 1.E-4_EB
REAC_MASS_ERROR             = 1.E-4_EB
REVERSIBLE                  = .FALSE.
SOOT_H_FRACTION             = 0.1_EB
SOOT_YIELD                  = 0.0_EB
SMIX_ID                     = 'null'
SPEC_ID                     = 'null'
THRESHOLD_TEMPERATURE       = 300._EB   ! For SIMPLE_CHEMSITRY w/ CO_PRODUCTION
 
END SUBROUTINE SET_REAC_DEFAULTS

END SUBROUTINE READ_REAC



SUBROUTINE PROC_REAC
USE PROPERTY_DATA, ONLY : PARSE_EQUATION, SHUTDOWN_ATOM
REAL(EB) :: MASS_PRODUCT,MASS_REACTANT,REACTION_BALANCE(118)
INTEGER :: NS,NS2,NR,NSPEC,NRR,NFR,R_COUNT,F_COUNT
LOGICAL :: NAME_FOUND,SKIP_ATOM_BALANCE
TYPE (SPECIES_MIXTURE_TYPE), POINTER :: SM
TYPE(REACTION_TYPE), POINTER :: RR=>NULL(),FR=>NULL()

IF (N_REACTIONS <=0) RETURN

!Basic input error checking
IF (SIMPLE_CHEMISTRY .AND. N_REACTIONS > 1) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: can not have more than one reaction when using simple chemistry'
   CALL SHUTDOWN(MESSAGE)
ENDIF

! The following information is what the user would have entered into the input file in the more general case

IF (SIMPLE_CHEMISTRY) THEN
   RN => REACTION(1)
   RN%SMIX_ID_READ(1) = RN%FUEL
   RN%SMIX_ID_READ(2) = 'AIR'
   RN%SMIX_ID_READ(3) = 'PRODUCTS'
   RN%NU_READ(1)      = -1._EB
   RN%NU_READ(2)      = -RN%NU_O2*SPECIES_MIXTURE(1)%VOLUME_FRACTION(FUEL_INDEX)/SPECIES_MIXTURE(0)%VOLUME_FRACTION(O2_INDEX)
   RN%N_SMIX          = 3
   RN%NU_READ(3)      = -(RN%NU_READ(1)*SPECIES_MIXTURE(1)%MW+RN%NU_READ(2)*SPECIES_MIXTURE(0)%MW)/SPECIES_MIXTURE(2)%MW
ENDIF

REAC_LOOP: DO NR=1,N_REACTIONS

   RN => REACTION(NR)
   
   IF ((RN%A > 0._EB .OR. RN%E > 0._EB) .AND. (RN%C>ZERO_P .OR. RN%H>ZERO_P)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: cannot use both finite rate REAC and simple chemistry'
      CALL SHUTDOWN(MESSAGE)
   ENDIF 
   
   IF (.NOT. SIMPLE_CHEMISTRY .AND. RN%HEAT_OF_COMBUSTION <-1.E10_EB) THEN
      WRITE(MESSAGE,'(A,I3,A)') 'ERROR: Problem with REAC ',NR,'. HEAT_OF_COMBUSTION not set.'
      CALL SHUTDOWN(MESSAGE)
   ENDIF 
  
   IF (TRIM(RN%EQUATION)/='null') THEN
      IF(ANY(ABS(RN%NU_READ)>ZERO_P)) THEN
         WRITE(MESSAGE,'(A,I3,A)') 'ERROR: Problem with REAC ',NR,'. Cannot set NUs if an EQUATION is specified.'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      CALL PARSE_EQUATION(NR)
      RN%N_SMIX = 0
      DO NS=1,N_TRACKED_SPECIES+1
         IF(ABS(RN%NU_READ(NS))>ZERO_P) THEN
            RN%N_SMIX = RN%N_SMIX+1
         ENDIF
      ENDDO
   ENDIF

   IF (TRIM(RN%FUEL)=='null') THEN
      WRITE(MESSAGE,'(A,I3,A)') 'ERROR: Problem with REAC ',N_REACTIONS,'. FUEL must be defined'
      CALL SHUTDOWN(MESSAGE)
   ENDIF  
   ! Allocate the arrays that are going to carry the mixture stoichiometry to the rest of the code

   ALLOCATE(RN%SMIX_ID(0:N_TRACKED_SPECIES))
   ALLOCATE(RN%NU(0:N_TRACKED_SPECIES))
   ALLOCATE(RN%SPEC_ID(1:N_SPECIES))
   ALLOCATE(RN%N_S(1:N_SPECIES))
   RN%SMIX_ID = 'null'
   RN%SPEC_ID = 'null'
   RN%NU      = 0._EB
   RN%N_S     = -999._EB

   ! Transfer SMIX_ID, SPEC_ID, NU, and N_S that were indexed by the order they were read in
   ! to now be indexed by the SMIX or SPEC index
   DO NS=1,RN%N_SMIX
      IF (TRIM(RN%SMIX_ID_READ(NS))=='null') CYCLE   
      NAME_FOUND = .FALSE.
      DO NS2=0,N_TRACKED_SPECIES
         IF (TRIM(RN%SMIX_ID_READ(NS))==TRIM(SPECIES_MIXTURE(NS2)%ID)) THEN
            RN%SMIX_ID(NS2) = RN%SMIX_ID_READ(NS)
            RN%NU(NS2)      = RN%NU_READ(NS)
            NAME_FOUND = .TRUE.
            EXIT
         ENDIF
         IF (TRIM(RN%EQUATION)/='null') THEN
            IF (TRIM(RN%SMIX_ID_READ(NS))==TRIM(SPECIES_MIXTURE(NS2)%FORMULA)) THEN
               RN%SMIX_ID(NS2) = SPECIES_MIXTURE(NS2)%ID
               RN%NU(NS2)      = RN%NU_READ(NS)
               NAME_FOUND = .TRUE.
               EXIT
            ENDIF
         ENDIF
      ENDDO
      IF (.NOT. NAME_FOUND) THEN
         WRITE(MESSAGE,'(A,I3,A,A,A)') 'ERROR: Problem with REAC ',NR,'. Tracked species ',TRIM(RN%SMIX_ID_READ(NS)),' not found.'
         CALL SHUTDOWN(MESSAGE)      
      ENDIF
   ENDDO

   ! Look for indices of fuels, oxidizers, and products. Normalize the stoichiometric coefficients by that of the fuel.
   DO NS2=0,N_TRACKED_SPECIES
      IF (RN%NU(NS2)>ZERO_P) I_PRODUCTS = NS2
      DO NSPEC=1,N_SPECIES
         IF (SPECIES_MIXTURE(NS2)%SPEC_ID(NSPEC)==RN%FUEL .OR. SPECIES_MIXTURE(NS2)%ID==RN%FUEL) THEN
            RN%FUEL_SMIX_INDEX = NS2
            RN%NU = -RN%NU/RN%NU(NS2)
            EXIT
         ENDIF
      ENDDO
   ENDDO


   ! Adjust mol/cm^3/s based rate to kg/m^3/s rate
   RN%RHO_EXPONENT = 0._EB
   DO NS=1,RN%N_SPEC
      IF (TRIM(RN%SPEC_ID_READ(NS))=='null') CYCLE
      IF (RN%A < 0.0_EB) CYCLE
      NAME_FOUND = .FALSE.
      DO NS2=1,N_SPECIES
         IF (TRIM(RN%SPEC_ID_READ(NS))==TRIM(SPECIES(NS2)%ID)) THEN
            RN%SPEC_ID(NS2) = RN%SPEC_ID_READ(NS)
            RN%N_S(NS2)     = RN%N_S_READ(NS)
            RN%A            = RN%A * (1000._EB*SPECIES(NS2)%MW)**(-RN%N_S(NS2))
            RN%RHO_EXPONENT = RN%RHO_EXPONENT + RN%N_S(NS2)
            NAME_FOUND = .TRUE.
            EXIT
         ENDIF
      ENDDO
      IF (.NOT. NAME_FOUND) THEN
         WRITE(MESSAGE,'(A,I3,A,A,A)') 'ERROR: Problem with REAC ',NR,'. Primitive species ', &
                                       TRIM(RN%SPEC_ID_READ(NS)),' not found.'
         CALL SHUTDOWN(MESSAGE)      
      ENDIF     
   ENDDO
   
   RN%RHO_EXPONENT = RN%RHO_EXPONENT - 1._EB
   RN%A = RN%A *1000._EB*SPECIES_MIXTURE(RN%FUEL_SMIX_INDEX)%MW

   IF (RN%FUEL/='null' .AND. RN%FUEL_SMIX_INDEX<1) THEN
      WRITE(MESSAGE,'(A,I3,A,A,A)') 'ERROR: Problem with REAC ',NR,'. Fuel ',TRIM(RN%FUEL),' not found.'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   ! Compute the primitive species reaction coefficients

   ALLOCATE(RN%NU_SPECIES(N_SPECIES))
   RN%NU_SPECIES = 0._EB
   DO NS=0,N_TRACKED_SPECIES
      SM => SPECIES_MIXTURE(NS)
      RN%NU(NS) = RN%NU(NS)*SM%ADJUST_NU
      DO NS2 = 1,N_SPECIES
         RN%NU_SPECIES(NS2) =  RN%NU_SPECIES(NS2) + RN%NU(NS)*SM%VOLUME_FRACTION(NS2)
      ENDDO
      IF (SM%ID=='WATER VAPOR')  I_WATER = NS
      IF (SM%ID=='CARBON DIOXIDE') I_CO2 = NS
   ENDDO

   ! Check atom balance of the reaction
   IF (.NOT. SIMPLE_CHEMISTRY .AND. RN%CHECK_ATOM_BALANCE) THEN
      SKIP_ATOM_BALANCE = .FALSE.
      REACTION_BALANCE = 0._EB
      DO NS=0,N_TRACKED_SPECIES
         IF (ABS(RN%NU(NS))>ZERO_P .AND. .NOT. SPECIES_MIXTURE(NS)%VALID_ATOMS) SKIP_ATOM_BALANCE = .TRUE.
         REACTION_BALANCE = REACTION_BALANCE + RN%NU(NS)*SPECIES_MIXTURE(NS)%ATOMS
      ENDDO
      IF (ANY(ABS(REACTION_BALANCE)>REAC_ATOM_ERROR) .AND. .NOT. SKIP_ATOM_BALANCE) &
         CALL SHUTDOWN_ATOM(REACTION_BALANCE,NR,REAC_ATOM_ERROR)
   ENDIF

   ! Check the mass balance of the reaction

   MASS_REACTANT = 0._EB
   MASS_PRODUCT  = 0._EB   
 
   DO NS=0,N_TRACKED_SPECIES
      IF (RN%NU(NS) < -ZERO_P) MASS_REACTANT = MASS_REACTANT + RN%NU(NS)*SPECIES_MIXTURE(NS)%MW
      IF (RN%NU(NS) >  ZERO_P) MASS_PRODUCT  = MASS_PRODUCT  + RN%NU(NS)*SPECIES_MIXTURE(NS)%MW
   ENDDO  
   IF (ABS(MASS_PRODUCT) < ZERO_P .OR. ABS(MASS_REACTANT) < ZERO_P) THEN
      IF (ABS(MASS_PRODUCT) <ZERO_P) WRITE(MESSAGE,'(A,I3,A)') 'ERROR: Problem with REAC ',NR,'. Products not specified.'
      IF (ABS(MASS_REACTANT)<ZERO_P) WRITE(MESSAGE,'(A,I3,A)') 'ERROR: Problem with REAC ',NR,'. Reactants not specified.'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (ABS(MASS_PRODUCT+MASS_REACTANT)/ABS(MASS_PRODUCT) > REAC_MASS_ERROR) THEN
      WRITE(MESSAGE,'(A,I3,A,F8.3,A,F8.3)') 'ERROR: Problem with REAC ',NR,'. Mass of products, ',MASS_PRODUCT, &
         ', does not equal mass of reactants,',-MASS_REACTANT
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   ! Heat of Combustion calculation for SIMPLE_CHEMISTRY

   IF (SIMPLE_CHEMISTRY) THEN
      IF (RN%HEAT_OF_COMBUSTION<0._EB) THEN
         RN%HEAT_OF_COMBUSTION = -RN%EPUMO2*RN%NU_SPECIES(O2_INDEX)*SPECIES(O2_INDEX)%MW/SPECIES(FUEL_INDEX)%MW
      ELSE
         IF (IDEAL) THEN
            RN%HEAT_OF_COMBUSTION = RN%HEAT_OF_COMBUSTION*SPECIES(FUEL_INDEX)%MW*0.001 !J/kg -> J/mol
            RN%HEAT_OF_COMBUSTION = RN%HEAT_OF_COMBUSTION - RN%NU_CO*(CO2_HEAT_OF_FORMATION - CO_HEAT_OF_FORMATION) &
                                                          - RN%NU_SOOT*CO2_HEAT_OF_FORMATION*(1._EB-RN%SOOT_H_FRACTION) &
                                                          - RN%NU_SOOT*H2O_HEAT_OF_FORMATION*RN%SOOT_H_FRACTION*0.5_EB
            RN%HEAT_OF_COMBUSTION = RN%HEAT_OF_COMBUSTION/SPECIES(FUEL_INDEX)%MW*1000._EB !J/mol->J/kg
         ENDIF
      ENDIF   
   ENDIF

   IF (NR==1) REACTION%HOC_COMPLETE = RN%HEAT_OF_COMBUSTION 


   IF (TRIM(ODE_SOLVER)/='null') THEN
      SELECT CASE (TRIM(ODE_SOLVER))
         CASE ('SINGLE EXACT')
            COMBUSTION_ODE = SINGLE_EXACT      
         CASE ('EXPLICIT EULER')
            COMBUSTION_ODE = EXPLICIT_EULER
         CASE ('RUNGE-KUTTA 2')
            COMBUSTION_ODE = RUNGE_KUTTA_2
         !CASE ('RUNGE-KUTTA 4')
         !   COMBUSTION_ODE = RUNGE_KUTTA_4
         CASE ('RK2 RICHARDSON')
            COMBUSTION_ODE = RK2_RICHARDSON
         CASE ('EDCM RK2')
            COMBUSTION_ODE = EDCM_RK2   
         CASE DEFAULT
            WRITE(MESSAGE,'(A)') 'ERROR: Problem with REAC. Name of ODE_SOLVER is not recognized.'
            CALL SHUTDOWN(MESSAGE)
      END SELECT
   ELSE
      IF (N_REACTIONS ==1 .AND. REACTION(1)%MODE==EDDY_DISSIPATION) THEN
         COMBUSTION_ODE = SINGLE_EXACT
      ELSE
         COMBUSTION_ODE = RK2_RICHARDSON
      ENDIF
   ENDIF
   
ENDDO REAC_LOOP

! Calculates Heat of Combustion for reverse reactions

REVERSE_LOOP: DO NRR=1,N_REACTIONS
   RR => REACTION(NRR)
   R_COUNT=0
   F_COUNT=0
    IF (.NOT. RR%REVERSIBLE) CYCLE REVERSE_LOOP
   R_COUNT = 1
   FORWARD_LOOP: DO NFR=1,N_REACTIONS
      FR => REACTION(NFR)     
      IF(.NOT. FR%ID == RR%FWD_ID) CYCLE FORWARD_LOOP
      F_COUNT = 1
      RR%HEAT_OF_COMBUSTION=-FR%HEAT_OF_COMBUSTION*(SPECIES_MIXTURE(FR%FUEL_SMIX_INDEX)%MW/SPECIES_MIXTURE(RR%FUEL_SMIX_INDEX)%MW)
      ENDDO FORWARD_LOOP
   
   IF(.NOT. R_COUNT /= F_COUNT) CYCLE REVERSE_LOOP
      WRITE(MESSAGE,'(A)') 'ERROR: Problem with REAC. Forward reaction is not recognized.'
      CALL SHUTDOWN(MESSAGE)
ENDDO REVERSE_LOOP       
 
! Change units of combustion quantities

HRRPUV_AVERAGE = HRRPUV_AVERAGE*1000._EB   ! W/m3
HRRPUA_SHEET   = HRRPUA_SHEET*  1000._EB   ! W/m2

END SUBROUTINE PROC_REAC



SUBROUTINE READ_PART

USE MATH_FUNCTIONS, ONLY : GET_RAMP_INDEX
USE DEVICE_VARIABLES, ONLY : PROPERTY_TYPE, PROPERTY, N_PROP
USE RADCONS, ONLY : NDG
USE MATH_FUNCTIONS, ONLY: NORM2,GET_TABLE_INDEX
INTEGER :: SAMPLING_FACTOR,N,NN,ILPC,IPC,RGB(3),I_DUMMY(6),N_STRATA
REAL(EB) :: DIAMETER, DENSITY,&
            GAMMA_D,AGE,INITIAL_TEMPERATURE,HEAT_OF_COMBUSTION, &
            VERTICAL_VELOCITY,HORIZONTAL_VELOCITY,MAXIMUM_DIAMETER,MINIMUM_DIAMETER,SIGMA_D, &
            SURFACE_TENSION,BREAKUP_RATIO,BREAKUP_GAMMA_D,BREAKUP_SIGMA_D,&
            DENSE_VOLUME_FRACTION,REAL_REFRACTIVE_INDEX,COMPLEX_REFRACTIVE_INDEX
REAL(EB) :: VEG_SV,VEG_MOISTURE,VEG_CHAR_FRACTION,VEG_DRAG_COEFFICIENT,VEG_DENSITY,VEG_BULK_DENSITY, &
            VEG_BURNING_RATE_MAX,VEG_DEHYDRATION_RATE_MAX,VEG_INITIAL_TEMPERATURE, &
            VEG_FUEL_MPV_MIN,VEG_MOIST_MPV_MIN,USER_DRAG_COEFFICIENT,FREE_AREA_FRACTION
REAL(EB), DIMENSION(1000,3) :: ORIENTATION
REAL(EB), DIMENSION(3) :: OR_TEMP
CHARACTER(30) :: SPEC_ID,DEVC_ID,CTRL_ID,QUANTITIES(1:10),SURF_ID,DRAG_LAW,PROP_ID,RADIATIVE_PROPERTY_TABLE='null',&
                 CNF_RAMP_ID='null',BREAKUP_CNF_RAMP_ID='null',DISTRIBUTION,BREAKUP_DISTRIBUTION
CHARACTER(25) :: COLOR
CHARACTER(25) :: VEG_DEGRADATION
LOGICAL :: MASSLESS,STATIC,TREE,MONODISPERSE,BREAKUP,CHECK_DISTRIBUTION,FUEL,WATER,&
           TURBULENT_DISPERSION
LOGICAL :: VEG_REMOVE_CHARRED,VEG_STEM,VEG_CHAR_OXIDATION
TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
NAMELIST /PART/ AGE,BREAKUP,BREAKUP_RATIO,BREAKUP_GAMMA_D,BREAKUP_SIGMA_D,CHECK_DISTRIBUTION,&
                BREAKUP_CNF_RAMP_ID,BREAKUP_DISTRIBUTION,CNF_RAMP_ID,COLOR,COMPLEX_REFRACTIVE_INDEX,&
                CTRL_ID,DENSE_VOLUME_FRACTION,&
                DENSITY,DEVC_ID,DIAMETER,DISTRIBUTION,DRAG_LAW,FREE_AREA_FRACTION,FYI,GAMMA_D,HEAT_OF_COMBUSTION,&
                HORIZONTAL_VELOCITY,ID,INITIAL_TEMPERATURE,MASSLESS,MAXIMUM_DIAMETER,MINIMUM_DIAMETER,MONODISPERSE,&
                N_STRATA,ORIENTATION,PROP_ID,QUANTITIES,RADIATIVE_PROPERTY_TABLE,REAL_REFRACTIVE_INDEX,RGB,&
                SAMPLING_FACTOR,SIGMA_D,SPEC_ID,STATIC,SURFACE_TENSION,SURF_ID,TREE,&
                TURBULENT_DISPERSION,USER_DRAG_COEFFICIENT,VEG_BULK_DENSITY,&
                VEG_BURNING_RATE_MAX,VEG_CHAR_FRACTION,VEG_CHAR_OXIDATION,VEG_DEGRADATION,VEG_DEHYDRATION_RATE_MAX,VEG_DENSITY,&
                VEG_DRAG_COEFFICIENT,VEG_FUEL_MPV_MIN,VEG_INITIAL_TEMPERATURE,VEG_MOISTURE,VEG_MOIST_MPV_MIN,VEG_REMOVE_CHARRED,&
                VEG_STEM,VEG_SV,VERTICAL_VELOCITY,&
                FUEL,WATER ! Backward compatibility

! Determine total number of PART lines in the input file

REWIND(LU_INPUT)
N_LAGRANGIAN_CLASSES = 0
COUNT_PART_LOOP: DO
   CALL CHECKREAD('PART',LU_INPUT,IOS) 
   IF (IOS==1) EXIT COUNT_PART_LOOP
   READ(LU_INPUT,PART,END=219,ERR=220,IOSTAT=IOS)
   N_LAGRANGIAN_CLASSES = N_LAGRANGIAN_CLASSES + 1
   220 IF (IOS>0) CALL SHUTDOWN('ERROR: Problem with PART line')
ENDDO COUNT_PART_LOOP
219 REWIND(LU_INPUT)

IF (N_LAGRANGIAN_CLASSES>0) PARTICLE_FILE = .TRUE.
 
! Allocate the derived type array to hold information about the particle classes

ALLOCATE(LAGRANGIAN_PARTICLE_CLASS(N_LAGRANGIAN_CLASSES),STAT=IZERO)
CALL ChkMemErr('READ','N_LAGRANGIAN_CLASSES',IZERO) 

N_LP_ARRAY_INDICES = 0
IPC = 0
ILPC = 0

READ_PART_LOOP: DO N=1,N_LAGRANGIAN_CLASSES
   
   BREAKUP                  = .FALSE.
   BREAKUP_RATIO            = 3._EB/7._EB  ! ratio of child Sauter mean to parent size in Bag breakup regime
   BREAKUP_GAMMA_D          = 2.4_EB
   BREAKUP_SIGMA_D          = -99999.9_EB
   CTRL_ID                  = 'null'
   DENSE_VOLUME_FRACTION    = 1.E-5_EB     ! Limiting volume fraction for drag reduction
   DENSITY                  = -1._EB
   DEVC_ID                  = 'null'
   INITIAL_TEMPERATURE      = TMPA - TMPM  ! C
   HEAT_OF_COMBUSTION       = -1._EB       ! kJ/kg
   DIAMETER                 = 500._EB      ! microns
   MAXIMUM_DIAMETER         = 1.E9_EB      ! microns, meant to be infinitely large and not used
   MINIMUM_DIAMETER         = -1._EB       ! microns, below which the PARTICLE evaporates in one time step
   MONODISPERSE             = .FALSE.
   N_STRATA                 = 7
   GAMMA_D                  = 2.4_EB
   SIGMA_D                  = -99999.9_EB
   AGE                      = 1.E6_EB      ! s
   ID                       = 'null'
   PROP_ID                  = 'null'
   ORIENTATION              = 0._EB
   QUANTITIES               = 'null'
   RADIATIVE_PROPERTY_TABLE = 'null'
   RGB                      = -1
   SPEC_ID                  = 'null'
   SURF_ID                  = 'null'
   SURFACE_TENSION          = 72.8E-3_EB  ! N/m, applies for water
   COLOR                    = 'null'
   SAMPLING_FACTOR          = -1      
   STATIC                   = .FALSE.
   MASSLESS                 = .FALSE.
   TREE                     = .FALSE.
   TURBULENT_DISPERSION     = .FALSE.
   REAL_REFRACTIVE_INDEX    = 1.33_EB
   COMPLEX_REFRACTIVE_INDEX = 0.01_EB
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
   DISTRIBUTION             = 'ROSIN-RAMMLER-LOGNORMAL'
   CNF_RAMP_ID              = 'null'
   CHECK_DISTRIBUTION       = .FALSE.
   BREAKUP_DISTRIBUTION     = 'ROSIN-RAMMLER-LOGNORMAL'
   BREAKUP_CNF_RAMP_ID      = 'null'
   FREE_AREA_FRACTION       = 0.5_EB
   WATER                    = .FALSE. ! Backward compatibility
   FUEL                     = .FALSE. ! Backward compatibility

   ! Read the PART line from the input file or set up special PARTICLE_CLASS class for water PARTICLEs or tracers
 
   CALL CHECKREAD('PART',LU_INPUT,IOS) 
   IF (IOS==1) EXIT READ_PART_LOOP
   READ(LU_INPUT,PART)
   
   LPC => LAGRANGIAN_PARTICLE_CLASS(N)
   
   ! Backward compatibility

   IF (FUEL) THEN
      IF (N_REACTIONS>0) THEN
         SPEC_ID = REACTION(1)%FUEL
      ELSE
         WRITE(MESSAGE,'(A)') 'ERROR: Cannot have FUEL particles without a REAC line.'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDIF
         
   IF (WATER) SPEC_ID = 'WATER VAPOR'

   ! If the user specifically specifies SURF_ID
   
   IF (SURF_ID/='null') VIRTUAL_PARTICLES = .TRUE.

   ! Miscellaneous consequences of input parameters

   IF (TREE) THEN 
      STATIC = .TRUE.
      WFDS_FE = .TRUE.
      SURF_ID = 'VEGETATION'
      IF (SAMPLING_FACTOR<=0) SAMPLING_FACTOR = 1
      IF (SPEC_ID/='null') THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Cannot have TREE=.TRUE. with a SPEC_ID'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDIF

   IF (SPEC_ID/='null') THEN
      SURF_ID = 'DROPLET'
      IF (SAMPLING_FACTOR<=0) SAMPLING_FACTOR = 10
   ENDIF

   IF (MASSLESS) THEN 
      DIAMETER = 0._EB
      SURF_ID  = 'MASSLESS PARTICLE'
      IF (SAMPLING_FACTOR<=0) SAMPLING_FACTOR = 1
   ENDIF

   ! If particle class has no ID at this point, stop.

   IF (SURF_ID=='null') THEN
      WRITE(MESSAGE,'(A,I2,A)') 'ERROR: PART ',N,' needs a SURF_ID.'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   ! Set default colors for Smokeview

   IF (TRIM(SPEC_ID)=='WATER VAPOR') THEN
      IF (ANY(RGB<0) .AND. COLOR=='null') COLOR='BLUE'
   ENDIF

   IF (ANY(RGB<0) .AND. COLOR=='null') THEN
      COLOR = 'BLACK'
      IF (TREE)  COLOR = 'GREEN'
   ENDIF

   IF (COLOR /= 'null') CALL COLOR2RGB(RGB,COLOR)

   ! Determine if the SPEC_ID is OK

   LPC%SPEC_ID = SPEC_ID
   IF (LPC%SPEC_ID/='null') THEN
      IF (MASSLESS) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: Cannot have MASSLESS=.TRUE. with evaporating PARTICLEs'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      DO NN=0,N_TRACKED_SPECIES
         IF (TRIM(SPECIES_MIXTURE(NN)%ID)==TRIM(LPC%SPEC_ID)) THEN
            LPC%Z_INDEX = NN
            SPECIES_MIXTURE(NN)%EVAPORATING = .TRUE.
            EXIT
         ENDIF
      ENDDO
      IF(LPC%Z_INDEX < 0) THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: PART SPEC_ID ',TRIM(LPC%SPEC_ID),' not found'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      DO NN=1,N_SPECIES
         IF (SPECIES_MIXTURE(LPC%Z_INDEX)%MASS_FRACTION(NN)>0._EB) THEN
            IF (LPC%Y_INDEX > 0) THEN
               WRITE(MESSAGE,'(A,I3,A)') 'ERROR: PART line ',N,'.  Particles cannot evaporate to a lumped species.'
               CALL SHUTDOWN(MESSAGE)            
            ELSE
               LPC%Y_INDEX = NN        
            ENDIF
         ENDIF
      ENDDO
      IF (SPECIES(LPC%Y_INDEX)%DENSITY_LIQUID > 0._EB .AND. DENSITY < 0._EB) DENSITY=SPECIES(LPC%Y_INDEX)%DENSITY_LIQUID
   ENDIF

   ! Arrays for particle size distribution

   LPC%N_STRATA = N_STRATA

   IF (DIAMETER > 0._EB) THEN
      ALLOCATE(LPC%CDF(0:NDC),STAT=IZERO)
      CALL ChkMemErr('READ','CDF',IZERO)
      ALLOCATE(LPC%R_CDF(0:NDC),STAT=IZERO)
      CALL ChkMemErr('READ','R_CDF',IZERO)
      ALLOCATE(LPC%IL_CDF(LPC%N_STRATA),STAT=IZERO)
      CALL ChkMemErr('READ','IL_CDF',IZERO)
      ALLOCATE(LPC%IU_CDF(LPC%N_STRATA),STAT=IZERO)
      CALL ChkMemErr('READ','IU_CDF',IZERO)
      ALLOCATE(LPC%W_CDF(LPC%N_STRATA),STAT=IZERO)
      CALL ChkMemErr('READ','W_CDF',IZERO)
   ENDIF

   ! Arrays related to particle break-up model

   IF (BREAKUP) THEN
      ALLOCATE(LPC%BREAKUP_CDF(0:NDC),STAT=IZERO)
      CALL ChkMemErr('READ','BREAKUP_CDF',IZERO)
      ALLOCATE(LPC%BREAKUP_R_CDF(0:NDC),STAT=IZERO)
      CALL ChkMemErr('READ','BREAKUP_R_CDF',IZERO)
   ENDIF      

   ! Radiative property table

   IF (RADIATIVE_PROPERTY_TABLE /= 'null') THEN
      CALL GET_TABLE_INDEX(RADIATIVE_PROPERTY_TABLE,PART_RADIATIVE_PROPERTY,LPC%RADIATIVE_PROPERTY_INDEX)
      LPC%RADIATIVE_PROPERTY_TABLE_ID = RADIATIVE_PROPERTY_TABLE
   ELSE
      LPC%RADIATIVE_PROPERTY_INDEX = 0
   ENDIF

   ! Assign property data to LAGRANGIAN_PARTICLE_CLASS class
 
   LPC%ID                     = ID
   LPC%BREAKUP                = BREAKUP
   LPC%BREAKUP_RATIO          = BREAKUP_RATIO   
   LPC%BREAKUP_GAMMA          = BREAKUP_GAMMA_D
   IF ( BREAKUP_SIGMA_D > 0._EB ) THEN
      LPC%BREAKUP_SIGMA = BREAKUP_SIGMA_D
   ELSE
      LPC%BREAKUP_SIGMA = 1.15_EB/BREAKUP_GAMMA_D
   END IF
   LPC%CTRL_ID                = CTRL_ID
   LPC%DENSE_VOLUME_FRACTION  = DENSE_VOLUME_FRACTION
   LPC%DENSITY                = DENSITY
   LPC%DEVC_ID            = DEVC_ID
   LPC%TMP_INITIAL        = INITIAL_TEMPERATURE + TMPM
   LPC%SAMPLING           = SAMPLING_FACTOR
   LPC%RGB                = RGB
   LPC%DIAMETER           = DIAMETER*1.E-6_EB
   LPC%MAXIMUM_DIAMETER   = MAXIMUM_DIAMETER*1.E-6_EB
   IF (MINIMUM_DIAMETER<0._EB) MINIMUM_DIAMETER=0.005_EB*DIAMETER
   LPC%MINIMUM_DIAMETER   = MINIMUM_DIAMETER*1.E-6_EB
   LPC%KILL_RADIUS        = MINIMUM_DIAMETER*1.E-6_EB*0.25_EB
   LPC%MONODISPERSE       = MONODISPERSE
   LPC%PROP_ID            = PROP_ID
   LPC%QUANTITIES         = QUANTITIES
   LPC%GAMMA              = GAMMA_D
   IF ( SIGMA_D > 0._EB ) THEN
      LPC%SIGMA           = SIGMA_D
   ELSE
      LPC%SIGMA           = 1.15_EB/GAMMA_D
   END IF
   LPC%DISTRIBUTION       = DISTRIBUTION
   LPC%CHECK_DISTRIBUTION = CHECK_DISTRIBUTION
   LPC%BREAKUP_DISTRIBUTION = BREAKUP_DISTRIBUTION
   LPC%CNF_RAMP_ID        = CNF_RAMP_ID
   LPC%BREAKUP_CNF_RAMP_ID  = BREAKUP_CNF_RAMP_ID
   
   IF(LPC%CNF_RAMP_ID/='null') THEN
        CALL GET_RAMP_INDEX(LPC%CNF_RAMP_ID,'DIAMETER',LPC%CNF_RAMP_INDEX)
   ENDIF
   IF(LPC%BREAKUP_CNF_RAMP_ID/='null') THEN
        CALL GET_RAMP_INDEX(LPC%BREAKUP_CNF_RAMP_ID,'DIAMETER',LPC%BREAKUP_CNF_RAMP_INDEX)
   ENDIF
   
   LPC%TMP_INITIAL        = INITIAL_TEMPERATURE + TMPM
   LPC%REAL_REFRACTIVE_INDEX = REAL_REFRACTIVE_INDEX
   LPC%COMPLEX_REFRACTIVE_INDEX = COMPLEX_REFRACTIVE_INDEX
   IF (LPC%REAL_REFRACTIVE_INDEX <= 0._EB .OR. LPC%COMPLEX_REFRACTIVE_INDEX < 0._EB) THEN
      WRITE(MESSAGE,'(A,A)') 'Bad refractive index on PART line ',LPC%ID
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   LPC%HEAT_OF_COMBUSTION = HEAT_OF_COMBUSTION*1000._EB
   LPC%FTPR               = FOTH*PI*DENSITY
   LPC%MASSLESS           = MASSLESS
   LPC%LIFETIME           = AGE
   LPC%TREE               = TREE
   LPC%TURBULENT_DISPERSION = TURBULENT_DISPERSION
   LPC%STATIC             = STATIC
   LPC%SPEC_ID            = SPEC_ID
   LPC%SURF_ID            = SURF_ID
   LPC%SURF_INDEX         = -1   
   LPC%SURFACE_TENSION    = SURFACE_TENSION
   LPC%ADJUST_EVAPORATION  = 1._EB   ! If H_O_C>0. this parameter will have to be reset later
   LPC%VERTICAL_VELOCITY   = VERTICAL_VELOCITY
   LPC%HORIZONTAL_VELOCITY = HORIZONTAL_VELOCITY
   LPC%USER_DRAG_COEFFICIENT = USER_DRAG_COEFFICIENT
   IF (USER_DRAG_COEFFICIENT>=0._EB) DRAG_LAW = 'USER'
   IF (LPC%TREE) DRAG_LAW='TREE'

   SELECT CASE(DRAG_LAW)
      CASE('SPHERE')
         LPC%DRAG_LAW = SPHERE_DRAG
      CASE('CYLINDER')
         LPC%DRAG_LAW = CYLINDER_DRAG
      CASE('USER')
         LPC%DRAG_LAW = USER_DRAG
      CASE('SCREEN')
         LPC%DRAG_LAW = SCREEN_DRAG
      CASE('TREE')
         LPC%DRAG_LAW = TREE_DRAG
      CASE DEFAULT
         WRITE(MESSAGE,'(A)') 'ERROR: unrecognized drag law on PART line'
         CALL SHUTDOWN(MESSAGE)
   END SELECT
    
   LPC%FREE_AREA_FRACTION = FREE_AREA_FRACTION

   ! Vegetation properties

   LPC%VEG_SV                   = VEG_SV !1/m
   LPC%VEG_MOISTURE             = VEG_MOISTURE
   LPC%VEG_CHAR_FRACTION        = VEG_CHAR_FRACTION
   LPC%VEG_DRAG_COEFFICIENT     = VEG_DRAG_COEFFICIENT
   LPC%VEG_DENSITY              = VEG_DENSITY !kg/m3
   LPC%VEG_BULK_DENSITY         = VEG_BULK_DENSITY !kg/m3
   LPC%VEG_BURNING_RATE_MAX     = VEG_BURNING_RATE_MAX !kg/m3.s
   LPC%VEG_DEHYDRATION_RATE_MAX = VEG_DEHYDRATION_RATE_MAX !kg/m3.s
   LPC%VEG_INITIAL_TEMPERATURE  = VEG_INITIAL_TEMPERATURE +TMPM ! K
   LPC%VEG_FUEL_MPV_MIN         = VEG_CHAR_FRACTION*VEG_BULK_DENSITY
   LPC%VEG_MOIST_MPV_MIN        = 0.01_EB*VEG_MOISTURE
   LPC%VEG_REMOVE_CHARRED       = VEG_REMOVE_CHARRED
   LPC%VEG_STEM                 = VEG_STEM
   LPC%VEG_DEGRADATION          = VEG_DEGRADATION
   LPC%VEG_CHAR_OXIDATION       = VEG_CHAR_OXIDATION

   ! Count and process the number of orientations for the particle

   LPC%N_ORIENTATION = 0

   DO NN=1,10
      IF (ANY(ABS(ORIENTATION(NN,1:3))>ZERO_P)) LPC%N_ORIENTATION = LPC%N_ORIENTATION + 1
   ENDDO

   IF (LPC%N_ORIENTATION>0) THEN   
      ALLOCATE(LPC%ORIENTATION(1:LPC%N_ORIENTATION,1:3),STAT=IZERO)
      CALL ChkMemErr('READ','ORIENTATION',IZERO)
      DO NN=1,LPC%N_ORIENTATION
         OR_TEMP(1:3) = ORIENTATION(NN,1:3)
         LPC%ORIENTATION(NN,1:3)   = ORIENTATION(NN,1:3)/ NORM2(OR_TEMP)
      ENDDO
   ENDIF

   ! Determine the number of slots to create in the particle evaporation and radiation arrays

   IF (.NOT. LPC%MASSLESS) THEN
      N_LP_ARRAY_INDICES = N_LP_ARRAY_INDICES + 1
      LPC%ARRAY_INDEX =  N_LP_ARRAY_INDICES
   ENDIF
   
ENDDO READ_PART_LOOP

! Allocate radiation arrays

PLOOP2: DO ILPC=1,N_LAGRANGIAN_CLASSES
   LPC=>LAGRANGIAN_PARTICLE_CLASS(ILPC)
   IF (.NOT. LPC%MASSLESS) THEN
      ALLOCATE(LPC%WQABS(0:NDG,1:NUMBER_SPECTRAL_BANDS))
      CALL ChkMemErr('INIT','WQABS',IZERO)
      LPC%WQABS = 0._EB
      ALLOCATE(LPC%WQSCA(0:NDG,1:NUMBER_SPECTRAL_BANDS))
      CALL ChkMemErr('INIT','WQSCA',IZERO)
      LPC%WQSCA = 0._EB
      ALLOCATE(LPC%R50(0:NDG))
      CALL ChkMemErr('INIT','R50',IZERO)
      LPC%R50 = 0._EB
   ENDIF
ENDDO PLOOP2

! Determine output quantities

DO ILPC=1,N_LAGRANGIAN_CLASSES
   LPC=>LAGRANGIAN_PARTICLE_CLASS(ILPC)
   LPC%N_QUANTITIES = 0
   IF (ANY(LPC%QUANTITIES/='null')) THEN
      QUANTITIES_LOOP: DO N=1,10
         IF (LPC%QUANTITIES(N)=='null') CYCLE QUANTITIES_LOOP
         LPC%N_QUANTITIES = LPC%N_QUANTITIES + 1
         CALL GET_QUANTITY_INDEX(LPC%SMOKEVIEW_LABEL(LPC%N_QUANTITIES),LPC%SMOKEVIEW_BAR_LABEL(LPC%N_QUANTITIES), &
                                 LPC%QUANTITIES_INDEX(LPC%N_QUANTITIES),I_DUMMY(1), &
                                 I_DUMMY(2),I_DUMMY(3),I_DUMMY(4),I_DUMMY(5),I_DUMMY(6),'PART', &
                                 LPC%QUANTITIES(N),'null','null','null','null','null')   
      ENDDO QUANTITIES_LOOP
   ENDIF
ENDDO 

! Adjust the evaporation rate of fuel PARTICLEs to account for difference in HoC.

DO ILPC=1,N_LAGRANGIAN_CLASSES
   LPC=>LAGRANGIAN_PARTICLE_CLASS(ILPC)
   IF (LPC%HEAT_OF_COMBUSTION > 0._EB) LPC%ADJUST_EVAPORATION = LPC%HEAT_OF_COMBUSTION/REACTION(1)%HEAT_OF_COMBUSTION
ENDDO

END SUBROUTINE READ_PART
 

 
SUBROUTINE PROC_PART

USE PROPERTY_DATA, ONLY: JANAF_TABLE_LIQUID
INTEGER :: N,NN,J,ITMP
REAL(EB) :: H_L,H_V,CPBAR,H_G_S,H_G_S_REF,H_L_REF,TMP_REF,TMP_MELT,TMP_V,TMP_WGT,DENSITY,MASS,VOLUME
TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
TYPE(SPECIES_TYPE),POINTER:: SS=>NULL()
TYPE(SURFACE_TYPE),POINTER:: SF=>NULL()

IF (N_LAGRANGIAN_CLASSES == 0) RETURN

PART_LOOP: DO N=1,N_LAGRANGIAN_CLASSES

   LPC => LAGRANGIAN_PARTICLE_CLASS(N)   
   SF  => SURFACE(LPC%SURF_INDEX)

   ! Assign device or controller

   CALL SEARCH_CONTROLLER('PART',LPC%CTRL_ID,LPC%DEVC_ID,LPC%DEVC_INDEX,LPC%CTRL_INDEX,N)

   ! Get density if the particles are liquid droplets or have mass

   IF (LPC%SURF_INDEX==DROPLET_SURF_INDEX) THEN
      CALL JANAF_TABLE_LIQUID (1,CPBAR,H_V,H_L,TMP_REF,TMP_MELT,TMP_V,SPECIES(LPC%Y_INDEX)%ID,LPC%FUEL,DENSITY)
      IF (LPC%DENSITY < 0._EB) LPC%DENSITY = DENSITY
   ENDIF

   IF (SF%THERMALLY_THICK) THEN
      MASS = 0._EB
      VOLUME = 0._EB
      DO NN=1,SF%N_CELLS
         SELECT CASE (SF%GEOMETRY)
            CASE (SURF_CARTESIAN)
               MASS = MASS + SF%LENGTH*SF%WIDTH*(SF%X_S(NN)-SF%X_S(NN-1))*SUM(SF%RHO_0(NN,1:SF%N_MATL))
               VOLUME = VOLUME + SF%LENGTH*SF%WIDTH*(SF%X_S(NN)-SF%X_S(NN-1))
            CASE (SURF_CYLINDRICAL)
               MASS = MASS + SF%LENGTH*PI*(SF%X_S(NN)**2-SF%X_S(NN-1)**2)*SUM(SF%RHO_0(NN,1:SF%N_MATL))
               VOLUME = VOLUME + SF%LENGTH*PI*(SF%X_S(NN)**2-SF%X_S(NN-1)**2)
            CASE (SURF_SPHERICAL)
               MASS = MASS + FOTHPI*(SF%X_S(NN)**3-SF%X_S(NN-1)**3)*SUM(SF%RHO_0(NN,1:SF%N_MATL))
               VOLUME = VOLUME + FOTHPI*(SF%X_S(NN)**3-SF%X_S(NN-1)**3)
         END SELECT
      ENDDO
      LPC%DENSITY = MASS/VOLUME
      LPC%FTPR = FOTH*PI*LPC%DENSITY               
   ENDIF

   ! Set the flag to do particle exchanges between meshes

   OMESH_PARTICLES=.TRUE.
   
   ! Only process DROPLETs

   SURF_OR_SPEC: IF (LPC%SURF_INDEX==DROPLET_SURF_INDEX) THEN
    
      SS => SPECIES(LPC%Y_INDEX)
      ALLOCATE(SS%C_P_L(0:5000),STAT=IZERO)
      CALL ChkMemErr('PROC_PART','SS%C_P_L',IZERO)
      SS%C_P_L=SS%SPECIFIC_HEAT_LIQUID
      ALLOCATE(SS%C_P_L_BAR(0:5000),STAT=IZERO)
      CALL ChkMemErr('PROC_PART','SS%C_P_L_BAR',IZERO)
      ALLOCATE(SS%H_L(0:5000),STAT=IZERO)
      CALL ChkMemErr('PROC_PART','SS%H_L',IZERO)
      ALLOCATE(SS%H_V(0:5000),STAT=IZERO)
      CALL ChkMemErr('PROC_PART','SS%H_V',IZERO)

      TMP_REF = -1._EB
      TMP_MELT = -1._EB
      TMP_V = -1._EB
      DO J = 1, 5000
         IF (SS%C_P_L(J) > 0._EB) THEN
            SS%H_L(J) = (REAL(J,EB)-SS%TMP_MELT)*SS%C_P_L(J)
            IF (J==1) THEN
               CALL JANAF_TABLE_LIQUID (J,CPBAR,H_V,H_L,TMP_REF,TMP_MELT,TMP_V,SS%ID,LPC%FUEL,DENSITY)
               IF (SS%H_V_REFERENCE_TEMPERATURE < 0._EB) SS%H_V_REFERENCE_TEMPERATURE=TMP_REF
               IF (SS%TMP_V < 0._EB) SS%TMP_V = TMP_V
               IF (LPC%DENSITY < 0._EB) LPC%DENSITY = DENSITY
               IF (LPC%DENSITY < 0._EB) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: PARTicle class ',TRIM(SS%ID),' requires a density'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
               LPC%FTPR = FOTH*PI*LPC%DENSITY               
               IF (SS%TMP_MELT < 0._EB) SS%TMP_MELT = TMP_MELT
            ENDIF
         ELSE
            CALL JANAF_TABLE_LIQUID (J,SS%C_P_L(J),H_V,H_L,TMP_REF,TMP_MELT,TMP_V,SS%ID,LPC%FUEL,DENSITY)
            IF (J==1) THEN         
               IF (LPC%DENSITY < 0._EB) LPC%DENSITY = DENSITY
               IF (LPC%DENSITY < 0._EB) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: PARTicle class ',TRIM(SS%ID),' requires a density'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF   
               LPC%FTPR = FOTH*PI*LPC%DENSITY                                 
               IF (SS%C_P_L(J) < 0._EB .AND. .NOT. LPC%TREE) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: PARTicle class ',TRIM(SS%ID),' requires CP, H_V, TMP_MELT, TMP_V, and T_REF'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
               IF (SS%H_V_REFERENCE_TEMPERATURE < 0._EB) SS%H_V_REFERENCE_TEMPERATURE=TMP_REF
               IF (SS%TMP_V < 0._EB) SS%TMP_V = TMP_V
               IF (SS%TMP_MELT < 0._EB) SS%TMP_MELT = TMP_MELT
               SS%H_L(J) = H_L + SS%C_P_L(J)
            ELSE
               SS%H_L(J) = SS%H_L(J-1) + 0.5_EB*(SS%C_P_L(J)+SS%C_P_L(J-1))
            ENDIF
         ENDIF
      END DO

      SS%C_P_L(0) = SS%C_P_L(1)
      SS%H_L(0) = SS%H_L(1)
      
      ! Adjust liquid H_L to force H_V at H_V_REFERENCE_TEMPERATURE

      IF(SS%HEAT_OF_VAPORIZATION > 0._EB) H_V = SS%HEAT_OF_VAPORIZATION
      ITMP = INT(SS%H_V_REFERENCE_TEMPERATURE)
      TMP_WGT  = SS%H_V_REFERENCE_TEMPERATURE - REAL(ITMP,EB)
      H_L_REF = SS%H_L(ITMP)+TMP_WGT*(SS%H_L(ITMP+1)-SS%H_L(ITMP))
      H_G_S_REF=(CPBAR_Z(ITMP,LPC%Z_INDEX)+TMP_WGT*(CPBAR_Z(ITMP+1,LPC%Z_INDEX)-CPBAR_Z(ITMP,LPC%Z_INDEX)))*&
               SS%H_V_REFERENCE_TEMPERATURE
      SS%H_L = SS%H_L + (H_G_S_REF - H_L_REF) - H_V

      ! Determine the properties of the PARTICLE

      DO J=1,5000
         H_G_S = CPBAR_Z(J,LPC%Z_INDEX)*REAL(J,EB)
         SS%H_V(J) = H_G_S - SS%H_L(J)         
         IF (J==1) THEN
            SS%C_P_L_BAR(J) =  SS%H_L(J)
         ELSE
            SS%C_P_L_BAR(J) = SS%H_L(J) / REAL(J,EB)
         ENDIF
      ENDDO         
      
      SS%H_V(0) = SS%H_V(1)  
      SS%C_P_L_BAR(0) = SS%H_L(1)
      
   ENDIF SURF_OR_SPEC
      
ENDDO PART_LOOP

END SUBROUTINE PROC_PART
 
 
SUBROUTINE READ_TREE

USE GLOBAL_CONSTANTS
INTEGER :: ILPC,N_TREES_0,NM,NN,N
REAL(EB) :: X_TREE_MIN,X_TREE_MAX,Y_TREE_MIN,Y_TREE_MAX,Z_TREE_MIN,Z_TREE_MAX, &
            X_OVERALL_MIN,X_OVERALL_MAX,Y_OVERALL_MIN,Y_OVERALL_MAX, &
            Z_OVERALL_MIN,Z_OVERALL_MAX
REAL(EB) :: CROWN_WIDTH,CROWN_WIDTH_BOTTOM,CROWN_WIDTH_TOP,CROWN_BASE_HEIGHT,TREE_HEIGHT,XYZ(3)
REAL(EB) :: RING_THICKNESS,TON_IGNITOR_ELEMENTS,TOFF_IGNITOR_ELEMENTS, &
            T_RAMPOFF_IGNITOR_ELEMENTS,T_RAMPON_IGNITOR_ELEMENTS
LOGICAL  :: IGNITOR_ELEMENTS,OUTPUT_TREE
CHARACTER(30) :: PART_ID, LABEL
CHARACTER(30) :: FUEL_GEOM,RAMP_IGNELEM
TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
NAMELIST /TREE/ CROWN_BASE_HEIGHT,CROWN_WIDTH,CROWN_WIDTH_BOTTOM,CROWN_WIDTH_TOP,FUEL_GEOM,IGNITOR_ELEMENTS,LABEL,OUTPUT_TREE,&
                PART_ID,RAMP_IGNELEM,RING_THICKNESS,TOFF_IGNITOR_ELEMENTS,TON_IGNITOR_ELEMENTS,TREE_HEIGHT,&
                T_RAMPOFF_IGNITOR_ELEMENTS,T_RAMPON_IGNITOR_ELEMENTS,XB,XYZ

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
   DO ILPC=1,N_LAGRANGIAN_CLASSES
      LPC=>LAGRANGIAN_PARTICLE_CLASS(ILPC)
      IF (LPC%ID==PART_ID) TREE_PARTICLE_CLASS(N) = ILPC
   ENDDO
   PARTICLE_FILE=.TRUE.
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
USE PHYSICAL_FUNCTIONS, ONLY : SPRAY_ANGLE_DISTRIBUTION
REAL(EB) :: ACTIVATION_OBSCURATION,ACTIVATION_TEMPERATURE,ALPHA_C,ALPHA_E,BETA_C,BETA_E, &
            BEAD_DIAMETER,BEAD_EMISSIVITY,BEAD_SPECIFIC_HEAT,BEAD_DENSITY,BEAD_H_FIXED, &
            C_FACTOR,CHARACTERISTIC_VELOCITY,ORIFICE_DIAMETER,DROPLET_VELOCITY, &
            PARTICLE_VELOCITY,FLOW_RATE,FLOW_TAU,GAUGE_TEMPERATURE,INITIAL_TEMPERATURE,K_FACTOR,LENGTH,SPRAY_ANGLE(2,2), &
            OFFSET,OPERATING_PRESSURE,RTI,PDPA_START,PDPA_END,PDPA_RADIUS,MASS_FLOW_RATE,&
            SPRAY_PATTERN_MU=-1._EB,SPRAY_PATTERN_BETA=5._EB, &
            PDPA_HISTOGRAM_LIMITS(2)=0._EB, &
            P0=0._EB,PX(3)=0._EB,PXX(3,3)=0._EB
EQUIVALENCE(PARTICLE_VELOCITY,DROPLET_VELOCITY)
INTEGER ::I,N,NN,PDPA_M,PDPA_N,PARTICLES_PER_SECOND,VELOCITY_COMPONENT=0,PDPA_HISTOGRAM_NBINS=-1
LOGICAL :: PDPA_INTEGRATE,PDPA_NORMALIZE,PDPA_HISTOGRAM=.FALSE.
EQUIVALENCE(LENGTH,ALPHA_C)
CHARACTER(30) :: SMOKEVIEW_ID(SMOKEVIEW_OBJECTS_DIMENSION),QUANTITY='null',PART_ID='null',FLOW_RAMP='null', &
                 SPRAY_PATTERN_TABLE='null',SPEC_ID='null',&
                 PRESSURE_RAMP='null',SMOKEVIEW_PARAMETERS(SMOKEVIEW_OBJECTS_DIMENSION), &
                 SPRAY_PATTERN_SHAPE='GAUSSIAN'
TYPE (PROPERTY_TYPE), POINTER :: PY=>NULL()

NAMELIST /PROP/ ACTIVATION_OBSCURATION,ACTIVATION_TEMPERATURE,ALPHA_C,ALPHA_E,BEAD_DENSITY,BEAD_DIAMETER,BEAD_EMISSIVITY,&
                BEAD_H_FIXED,BEAD_SPECIFIC_HEAT,BETA_C,BETA_E,&
                CHARACTERISTIC_VELOCITY,C_FACTOR,PARTICLES_PER_SECOND,&
                PARTICLE_VELOCITY,FLOW_RAMP,FLOW_RATE,FLOW_TAU,GAUGE_TEMPERATURE,ID,INITIAL_TEMPERATURE,K_FACTOR,LENGTH, &
                MASS_FLOW_RATE,OFFSET,OPERATING_PRESSURE,ORIFICE_DIAMETER,P0,PART_ID,PDPA_END,PDPA_HISTOGRAM, &
                PDPA_HISTOGRAM_LIMITS,PDPA_HISTOGRAM_NBINS,PDPA_INTEGRATE,PDPA_M,PDPA_N,PDPA_NORMALIZE,PDPA_RADIUS,&
                PDPA_START,PRESSURE_RAMP,PX,PXX,QUANTITY,RTI,SMOKEVIEW_ID,SMOKEVIEW_PARAMETERS,SPEC_ID,SPRAY_ANGLE,&
                SPRAY_PATTERN_BETA,SPRAY_PATTERN_MU,SPRAY_PATTERN_SHAPE,SPRAY_PATTERN_TABLE,VELOCITY_COMPONENT,&
                DROPLET_VELOCITY !Backwards compatability

! Count the PROP lines in the input file. Note how many of these are cables.

N_PROP=0
REWIND(LU_INPUT)
COUNT_PROP_LOOP: DO
   CALL CHECKREAD('PROP',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_PROP_LOOP
   READ(LU_INPUT,PROP,ERR=34,IOSTAT=IOS)
   N_PROP = N_PROP + 1
   34 IF (IOS>0) THEN
         WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with PROP number', N_PROP+1
         CALL SHUTDOWN(MESSAGE)
      ENDIF
ENDDO COUNT_PROP_LOOP
 
! Allocate the PROPERTY derived types
 
ALLOCATE(PROPERTY(0:N_PROP),STAT=IZERO)
CALL ChkMemErr('READ','PROPERTY',IZERO) 

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
   PY%BEAD_SPECIFIC_HEAT       = BEAD_SPECIFIC_HEAT*1000._EB/TIME_SHRINK_FACTOR
   PY%C_FACTOR                 = C_FACTOR
   PY%CHARACTERISTIC_VELOCITY  = CHARACTERISTIC_VELOCITY
   PY%GAUGE_TEMPERATURE        = GAUGE_TEMPERATURE + TMPM
   PY%ID                       = ID
   PY%INITIAL_TEMPERATURE      = INITIAL_TEMPERATURE + TMPM
   PY%PARTICLES_PER_SECOND      = PARTICLES_PER_SECOND
   PY%OFFSET                   = OFFSET
   PY%OPERATING_PRESSURE       = OPERATING_PRESSURE
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
   PY%SPRAY_ANGLE = SPRAY_ANGLE*PI/180._EB
   IF(ANY(PY%SPRAY_ANGLE(1:2,2)<0)) PY%SPRAY_ANGLE(1:2,2)=PY%SPRAY_ANGLE(1:2,1)
   SPRAY_PATTERN_MU=SPRAY_PATTERN_MU*PI/180._EB
   IF (PART_ID/='null' .AND. SPRAY_PATTERN_TABLE == 'null' ) THEN
      ALLOCATE(PY%SPRAY_LON_CDF(0:NDC2),PY%SPRAY_LON(0:NDC2),PY%SPRAY_LAT(0:NDC2),PY%SPRAY_LAT_CDF(0:NDC2,0:NDC2))
      IF(SPRAY_PATTERN_MU<0._EB) THEN
         IF(SPRAY_ANGLE(1,1)>0._EB) THEN
            SPRAY_PATTERN_MU=0.5_EB*SUM(PY%SPRAY_ANGLE(1:2,1))
         ELSE
            SPRAY_PATTERN_MU=0._EB
         ENDIF
      ENDIF
      CALL SPRAY_ANGLE_DISTRIBUTION(PY%SPRAY_LON,PY%SPRAY_LAT,PY%SPRAY_LON_CDF,PY%SPRAY_LAT_CDF, &
                                      SPRAY_PATTERN_BETA,SPRAY_PATTERN_MU,PY%SPRAY_ANGLE &
                                      ,SPRAY_PATTERN_SHAPE,NDC2)
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
       (PY%QUANTITY == 'PARTICLE FLUX X')    .OR. &
       (PY%QUANTITY == 'PARTICLE FLUX Y')    .OR. &
       (PY%QUANTITY == 'PARTICLE FLUX Z')) THEN
      PY%PDPA_M        = 3
      PY%PDPA_N        = 0
   ENDIF

   ! Histograms of PDPA data
   PY%PDPA_HISTOGRAM             = PDPA_HISTOGRAM
   PY%PDPA_HISTOGRAM_NBINS       = PDPA_HISTOGRAM_NBINS
   PY%PDPA_HISTOGRAM_LIMITS      = PDPA_HISTOGRAM_LIMITS
   IF(PDPA_HISTOGRAM) THEN
      HISTOGRAM_FILE=.TRUE.
      IF(PDPA_HISTOGRAM_NBINS<2) THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: Problem with PROP ',TRIM(PY%ID),', PDPA_HISTOGRAM needs PDPA_NBINS>2'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
       IF(ABS(PDPA_HISTOGRAM_LIMITS(1)-PDPA_HISTOGRAM_LIMITS(2)) < ZERO_P) THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: Problem with PROP ',TRIM(PY%ID),', PDPA_HISTOGRAM needs PDPA_LIMITS'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      
   ENDIF
   
   PATCH_VELOCITY_IF: IF (VELOCITY_COMPONENT>0) THEN
      PY%I_VEL = VELOCITY_COMPONENT
      PY%P0 = P0  ! value at origin of Taylor expansion
      DO J=1,3
         PY%PX(J) = PX(J)  ! first derivative of P evaluated at origin
         DO I=1,3
            IF (I>J) PXX(I,J)=PXX(J,I) ! make symmetric
            PY%PXX(I,J) = PXX(I,J) ! second derivative of P evaluated at origin
         ENDDO
      ENDDO
   ENDIF PATCH_VELOCITY_IF

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
      IF (MASS_FLOW_RATE > 0._EB) THEN
         PY%MASS_FLOW_RATE = MASS_FLOW_RATE
         IF (ABS(PARTICLE_VELOCITY) <= ZERO_P) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: Problem with PROP ',TRIM(PY%ID),', must specify PARTICLE_VELOCITY with MASS_FLOW_RATE'
            CALL SHUTDOWN(MESSAGE)         
         ELSE
            PY%PARTICLE_VELOCITY  = PARTICLE_VELOCITY
         ENDIF
      ELSE
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

         IF (PARTICLE_VELOCITY<=ZERO_P .AND. ORIFICE_DIAMETER<=ZERO_P .AND. &
            PRESSURE_RAMP=='null' .AND. SPRAY_PATTERN_TABLE=='null') THEN
            WRITE(MESSAGE,'(A,A,A)') 'WARNING: PROP ',TRIM(PY%ID),' PARTICLE velocity is not defined.'
            IF (MYID==0) WRITE(LU_ERR,'(A)') TRIM(MESSAGE)
         ENDIF
         
         IF (PARTICLE_VELOCITY > 0._EB) THEN
            PY%PARTICLE_VELOCITY  = PARTICLE_VELOCITY
         ELSEIF ((ORIFICE_DIAMETER > 0._EB) .AND. (FLOW_RATE > 0._EB)) THEN
            PY%PARTICLE_VELOCITY  = (FLOW_RATE/60._EB/1000._EB)/(PI*(ORIFICE_DIAMETER/2._EB)**2)
         ENDIF
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
      CALL GET_SPEC_OR_SMIX_INDEX(PY%SPEC_ID,PY%Y_INDEX,PY%Z_INDEX)
      IF (PY%Y_INDEX<1 .AND. PY%Z_INDEX<0) THEN
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
C_FACTOR                 = 0.0_EB
CHARACTERISTIC_VELOCITY  = 1.0_EB      ! m/s
PARTICLE_VELOCITY         = 0._EB       ! m/s
PARTICLES_PER_SECOND      = 5000
!DT_INSERT                = 0.01        ! s
FLOW_RATE                = -1._EB      ! L/min
FLOW_RAMP                = 'null'
FLOW_TAU                 = 0._EB
GAUGE_TEMPERATURE        = TMPA - TMPM
INITIAL_TEMPERATURE      = TMPA - TMPM
ID                       = 'null'
K_FACTOR                 = -1.0_EB     ! L/min/bar**0.5
MASS_FLOW_RATE           = -1._EB      ! kg/s
OFFSET                   = 0.05_EB     ! m
OPERATING_PRESSURE       = -1.0_EB     ! bar
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
SPRAY_ANGLE(1,1)           = 60._EB      ! degrees
SPRAY_ANGLE(2,1)           = 75._EB      ! degrees
SPRAY_ANGLE(1,2)           = -999._EB      ! degrees
SPRAY_ANGLE(2,2)           = -999._EB      ! degrees
SPRAY_PATTERN_TABLE      = 'null'
SPRAY_PATTERN_SHAPE      = 'GAUSSIAN'
SPRAY_PATTERN_MU         = -1._EB
SPRAY_PATTERN_BETA       = 5.0_EB
PDPA_HISTOGRAM           =.FALSE.
END SUBROUTINE SET_PROP_DEFAULTS
 
END SUBROUTINE READ_PROP
 


SUBROUTINE PROC_PROP
USE DEVICE_VARIABLES
REAL(EB) :: TOTAL_FLOWRATE, SUBTOTAL_FLOWRATE
INTEGER :: N,NN,N_V_FACTORS,ILPC
LOGICAL :: TABLE_NORMED(1:N_TABLE)
TYPE (PROPERTY_TYPE), POINTER :: PY=>NULL()
TYPE (TABLES_TYPE),  POINTER :: TA=>NULL()
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE),POINTER :: LPC=>NULL()

TABLE_NORMED = .FALSE.

PROP_LOOP: DO N=0,N_PROP
   PY => PROPERTY(N)

   ! Check to see if density is defined for sprinkler / nozzle
   IF (PY%FLOW_RATE > 0._EB .OR. PY%OPERATING_PRESSURE > 0._EB) THEN
      DO NN=1,N_LAGRANGIAN_CLASSES
         IF (TRIM(LAGRANGIAN_PARTICLE_CLASS(NN)%ID)==TRIM(PY%PART_ID)) THEN
            IF(LAGRANGIAN_PARTICLE_CLASS(NN)%DENSITY <0._EB) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: Problem with PROP ',TRIM(PY%ID),', referenced PART_ID has no specified density.'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
         ENDIF
      ENDDO
   ENDIF

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
         PY%V_FACTOR = PY%PARTICLE_VELOCITY/SQRT(PY%OPERATING_PRESSURE)
      ENDIF
   ENDIF

   ! Assign PART_INDEX to Device PROPERTY array

   IF (PY%PART_ID/='null') THEN
      DO ILPC=1,N_LAGRANGIAN_CLASSES
         LPC => LAGRANGIAN_PARTICLE_CLASS(ILPC)
         IF (LPC%ID==PY%PART_ID) PY%PART_INDEX = ILPC
         IF (LPC%ID==PY%PART_ID .AND. LPC%MASSLESS) THEN
            WRITE(MESSAGE,'(A,I4,A)') 'ERROR: PART_ID for PROP ' ,N,' cannot refer to MASSLESS particles'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO
      IF (PY%PART_INDEX<0) THEN
         WRITE(MESSAGE,'(A,I4,A)') 'ERROR: PART_ID for PROP ' ,N,' not found'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      PARTICLE_FILE=.TRUE.
      IF (PY%FLOW_RATE > 0._EB .AND. LPC%DENSITY < 0._EB) THEN
         WRITE(MESSAGE,'(A,A,A)') 'PROP ERROR: PARTicle class ',TRIM(LPC%ID),' requires a density'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   ENDIF

ENDDO PROP_LOOP

END SUBROUTINE PROC_PROP



SUBROUTINE READ_MATL

USE MATH_FUNCTIONS, ONLY : GET_RAMP_INDEX
USE DEVICE_VARIABLES, ONLY : PROPERTY
CHARACTER(30) :: CONDUCTIVITY_RAMP,SPECIFIC_HEAT_RAMP,SPEC_ID(MAX_SPECIES,MAX_REACTIONS)
REAL(EB) :: EMISSIVITY,CONDUCTIVITY,SPECIFIC_HEAT,DENSITY,ABSORPTION_COEFFICIENT,BOILING_TEMPERATURE, &
            INITIAL_VAPOR_FLUX,PEAK_REACTION_RATE
REAL(EB), DIMENSION(MAX_MATERIALS,MAX_REACTIONS) :: NU_MATL
REAL(EB), DIMENSION(MAX_REACTIONS) :: A,E,HEATING_RATE,PYROLYSIS_RANGE,HEAT_OF_REACTION, &
                                        N_S,N_T,REFERENCE_RATE,REFERENCE_TEMPERATURE,THRESHOLD_TEMPERATURE,HEAT_OF_COMBUSTION, &
                                        THRESHOLD_SIGN
REAL(EB), DIMENSION(MAX_SPECIES,MAX_REACTIONS) :: NU_SPEC
LOGICAL, DIMENSION(MAX_REACTIONS) :: PCR
CHARACTER(30), DIMENSION(MAX_MATERIALS,MAX_REACTIONS) :: MATL_ID
INTEGER :: N,NN,NNN,IOS,NR,N_REACTIONS
NAMELIST /MATL/ A,ABSORPTION_COEFFICIENT,BOILING_TEMPERATURE,CONDUCTIVITY,CONDUCTIVITY_RAMP,DENSITY,E,EMISSIVITY,FYI,&
                HEATING_RATE,HEAT_OF_COMBUSTION,HEAT_OF_REACTION,ID,INITIAL_VAPOR_FLUX,MATL_ID,NU_MATL,NU_SPEC,N_REACTIONS,&
                N_S,N_T,PCR,PYROLYSIS_RANGE,REFERENCE_RATE,REFERENCE_TEMPERATURE,SPECIFIC_HEAT,SPECIFIC_HEAT_RAMP,SPEC_ID,&
                THRESHOLD_SIGN,THRESHOLD_TEMPERATURE

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

! Allocate the MATERIAL derived type
 
ALLOCATE(MATERIAL(1:N_MATL),STAT=IZERO)
CALL ChkMemErr('READ','MATERIAL',IZERO) 

! Read the MATL lines in the order listed in the input file
 
REWIND(LU_INPUT)

READ_MATL_LOOP: DO N=1,N_MATL
 
   ML => MATERIAL(N)

   ! Read user defined MATL lines

   CALL CHECKREAD('MATL',LU_INPUT,IOS)
   CALL SET_MATL_DEFAULTS
   READ(LU_INPUT,MATL) 

   ! Do some error checking on the inputs

   NOT_BOILING: IF (BOILING_TEMPERATURE>4000._EB) THEN

      IF ( ( ANY(THRESHOLD_TEMPERATURE>-TMPM) .OR. ANY(REFERENCE_TEMPERATURE>-TMPM) .OR. ANY(A>=0._EB) .OR. ANY(E>=0._EB) .OR. &
             ANY(ABS(HEAT_OF_REACTION)>ZERO_P) ) .AND. N_REACTIONS==0) THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: Problem with MATL number ',TRIM(ID),'. A reaction parameter is used, but N_REACTIONS=0'  
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      DO NR=1,N_REACTIONS
         IF (REFERENCE_TEMPERATURE(NR)<-TMPM  .AND. (E(NR)< 0._EB .OR. A(NR)<0._EB)) THEN
            WRITE(MESSAGE,'(A,A,A,I2,A)') 'ERROR: Problem with MATL ',TRIM(ID),', REAC ',NR,'. Set REFERENCE_TEMPERATURE or E, A'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         IF (ABS(SUM(NU_MATL(:,NR)))<=ZERO_P .AND. ABS(SUM(NU_SPEC(:,NR)))<=ZERO_P) THEN
            WRITE(MESSAGE,'(A,A,A,I2,A)') 'WARNING: MATL ',TRIM(ID),', REAC ',NR,'. No product yields (NUs) set'  
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
   ML%NU_RESIDUE           = NU_MATL
   ML%NU_SPEC              = NU_SPEC
   ML%SPEC_ID              = SPEC_ID
   ML%RAMP_C_S             = SPECIFIC_HEAT_RAMP
   ML%RAMP_K_S             = CONDUCTIVITY_RAMP
   ML%RHO_S                = DENSITY
   ML%RESIDUE_MATL_NAME    = MATL_ID
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
   ELSE
      ML%PYROLYSIS_MODEL = PYROLYSIS_SOLID
      IF (N_REACTIONS==0) ML%PYROLYSIS_MODEL = PYROLYSIS_NONE
   ENDIF

   IF (ML%RAMP_K_S/='null') THEN
      CALL GET_RAMP_INDEX(ML%RAMP_K_S,'TEMPERATURE',NR)
      ML%K_S = -NR
   ENDIF

   IF (ML%RAMP_C_S/='null') THEN
      CALL GET_RAMP_INDEX(ML%RAMP_C_S,'TEMPERATURE',NR)
      ML%C_S = -NR
   ENDIF

   ! Determine A and E if REFERENCE_TEMPERATURE is specified

   DO NR=1,ML%N_REACTIONS
      IF (ML%TMP_REF(NR) > 0._EB) THEN
         IF (ML%RATE_REF(NR) > 0._EB) THEN
            PEAK_REACTION_RATE = ML%RATE_REF(NR)
         ELSE
            PEAK_REACTION_RATE = 2._EB*ML%HEATING_RATE(NR)*(1._EB-SUM(ML%NU_RESIDUE(:,NR)))/ML%PYROLYSIS_RANGE(NR)
         ENDIF
         ML%E(NR) = EXP(1._EB)*PEAK_REACTION_RATE*R0*ML%TMP_REF(NR)**2/ML%HEATING_RATE(NR)
         ML%A(NR) = EXP(1._EB)*PEAK_REACTION_RATE*EXP(ML%E(NR)/(R0*ML%TMP_REF(NR)))
      ENDIF

      ML%N_RESIDUE(NR) = 0
      DO NN=1,MAX_MATERIALS
         IF (ML%RESIDUE_MATL_NAME(NN,NR)/='null') ML%N_RESIDUE(NR) = ML%N_RESIDUE(NR) + 1
      ENDDO
   ENDDO

ENDDO READ_MATL_LOOP
 
! Assign a material index to the RESIDUEs

DO N=1,N_MATL
   ML => MATERIAL(N)
   ML%RESIDUE_MATL_INDEX = 0
   DO NR=1,ML%N_REACTIONS
      DO NN=1,ML%N_RESIDUE(NR)
         DO NNN=1,N_MATL
            IF (MATL_NAME(NNN)==ML%RESIDUE_MATL_NAME(NN,NR)) ML%RESIDUE_MATL_INDEX(NN,NR) = NNN
         ENDDO
         IF (ML%RESIDUE_MATL_INDEX(NN,NR)==0 .AND. ML%NU_RESIDUE(NN,NR)>0._EB) THEN
            WRITE(MESSAGE,'(5A)') 'ERROR: Residue ', TRIM(ML%RESIDUE_MATL_NAME(NN,NR)),' of ',TRIM(MATL_NAME(N)),' is not defined.'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDDO
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
NU_SPEC                = 0._EB
NU_MATL                = 0._EB
PCR                    = .FALSE.
REFERENCE_RATE         = -1._EB
REFERENCE_TEMPERATURE  = -1000._EB
MATL_ID                = 'null'
SPECIFIC_HEAT          = 0.0_EB      ! kJ/kg/K
SPECIFIC_HEAT_RAMP     = 'null'
SPEC_ID                = 'null'
HEATING_RATE           = 5._EB       ! K/min
PYROLYSIS_RANGE        = 80._EB      ! K or C
 
END SUBROUTINE SET_MATL_DEFAULTS
 
END SUBROUTINE READ_MATL



SUBROUTINE PROC_MATL

! Process Materials -- do some additional set-up work with materials

INTEGER :: N,J,NS,NS2,NR,Z_INDEX(N_TRACKED_SPECIES,MAX_REACTIONS)

PROC_MATL_LOOP: DO N=1,N_MATL
   
   ML => MATERIAL(N)

   ! Convert ML%NU_SPEC(I_ORDINAL,I_REACTION) and ML%SPEC_ID(I_ORDINAL,I_REACTION) to ML%NU_GAS(I_SPECIES,I_REACTION)

   Z_INDEX = -1
   DO NR=1,ML%N_REACTIONS
      DO NS=1,MAX_SPECIES
         
         IF (TRIM(ML%SPEC_ID(NS,NR))/='null' .NEQV. ML%NU_SPEC(NS,NR)>ZERO_P) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: MATL ',TRIM(MATL_NAME(N)),' requires both a SPEC_ID and NU_SPEC'
            CALL SHUTDOWN(MESSAGE)                  
         ENDIF
         IF (TRIM(ML%SPEC_ID(NS,NR))=='null') EXIT
         IF (NS==2 .AND. ML%PYROLYSIS_MODEL==PYROLYSIS_LIQUID) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: MATL ',TRIM(MATL_NAME(N)),' can only specify one SPEC_ID for a liquid'
            CALL SHUTDOWN(MESSAGE)                  
         ENDIF
         DO NS2=0,N_TRACKED_SPECIES
            IF (TRIM(ML%SPEC_ID(NS,NR))==TRIM(SPECIES_MIXTURE(NS2)%ID)) THEN
               Z_INDEX(NS,NR) = NS2
               ML%NU_GAS(Z_INDEX(NS,NR),NR) = ML%NU_SPEC(NS,NR)
               EXIT
            ENDIF
         ENDDO
         IF (Z_INDEX(NS,NR)==-1) THEN
            WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: SPECies ',TRIM(ML%SPEC_ID(NS,NR)),&
                                         ' corresponding to MATL ',TRIM(MATL_NAME(N)),' is not a tracked species'
            CALL SHUTDOWN(MESSAGE)                  
         ENDIF      

      ENDDO
   ENDDO
      
   ! Adjust burn rate if heat of combustion is different from the gas phase reaction value

   DO J=1,MAX(1,ML%N_REACTIONS)
      IF (N_REACTIONS>0) THEN
         RN => REACTION(1)
         DO NS = 1,N_TRACKED_SPECIES
            IF (ML%HEAT_OF_COMBUSTION(J)>0._EB .AND. RN%HEAT_OF_COMBUSTION>0._EB)  &
                ML%ADJUST_BURN_RATE(NS,J) = ML%HEAT_OF_COMBUSTION(J)/RN%HEAT_OF_COMBUSTION
         ENDDO
      ENDIF
   ENDDO

   ! Check units of specific heat

   IF (ML%RAMP_C_S/='null') THEN
         IF (RAMPS(-INT(ML%C_S))%DEPENDENT_DATA(1) > 10._EB) &
            WRITE(LU_ERR,'(A,A)') 'WARNING: SPECIFIC_HEAT units are kJ/kg/K check MATL ',TRIM(ID)
   ENDIF   

ENDDO PROC_MATL_LOOP

END SUBROUTINE PROC_MATL


SUBROUTINE READ_SURF
 
USE MATH_FUNCTIONS, ONLY : GET_RAMP_INDEX
USE DEVICE_VARIABLES, ONLY : PROPERTY_TYPE,PROPERTY
CHARACTER(30) :: PART_ID,RAMP_MF(MAX_SPECIES),RAMP_Q,RAMP_V,RAMP_T,MATL_ID(MAX_LAYERS,MAX_MATERIALS),&
                 PROFILE,BACKING,GEOMETRY,NAME_LIST(MAX_MATERIALS),EXTERNAL_FLUX_RAMP,RAMP_EF,RAMP_PART,SPEC_ID(MAX_SPECIES)
EQUIVALENCE(EXTERNAL_FLUX_RAMP,RAMP_EF)
LOGICAL :: ADIABATIC,BURN_AWAY,SHRINK,FREE_SLIP,NO_SLIP
CHARACTER(60) :: TEXTURE_MAP,BOUNDARY_CONDITION_FILENAME
CHARACTER(25) :: COLOR
REAL(EB) :: TAU_Q,TAU_V,TAU_T,TAU_MF(MAX_SPECIES),HRRPUA,MLRPUA,TEXTURE_WIDTH,TEXTURE_HEIGHT,VEL_T(2), &
            TAU_EXTERNAL_FLUX,TAU_EF,E_COEFFICIENT,VOLUME_FLUX,TMP_FRONT,TMP_INNER(MAX_LAYERS),THICKNESS(MAX_LAYERS),VEL, &
            MASS_FLUX(MAX_SPECIES),MASS_FRACTION(MAX_SPECIES), Z0,PLE,CONVECTIVE_HEAT_FLUX,PARTICLE_MASS_FLUX, &
            TRANSPARENCY,EXTERNAL_FLUX,TMP_BACK,MASS_FLUX_TOTAL,STRETCH_FACTOR(MAX_LAYERS),CONVECTION_LENGTH_SCALE, &
            MATL_MASS_FRACTION(MAX_LAYERS,MAX_MATERIALS),CELL_SIZE_FACTOR,MAX_PRESSURE,&
            IGNITION_TEMPERATURE,HEAT_OF_VAPORIZATION,REGRID_FACTOR,NET_HEAT_FLUX,LAYER_DIVIDE,SURFACE_DENSITY, &
            ROUGHNESS,RADIUS,LENGTH,WIDTH,DT_INSERT,H_FIXED,TAU_PART,EMISSIVITY,EMISSIVITY_BACK,EMISSIVITY_DEFAULT, &
            SPREAD_RATE,XYZ(3),MINIMUM_LAYER_THICKNESS
EQUIVALENCE(TAU_EXTERNAL_FLUX,TAU_EF)
INTEGER :: NPPC,N,IOS,NL,NN,NNN,NRM,N_LIST,N_LIST2,INDEX_LIST(MAX_MATERIALS_TOTAL),LEAK_PATH(2),DUCT_PATH(2),RGB(3),NR,IL
INTEGER ::  VEGETATION_LAYERS
REAL(EB) :: VEGETATION_CDRAG,VEGETATION_CHAR_FRACTION,VEGETATION_ELEMENT_DENSITY,VEGETATION_HEIGHT, &
            VEGETATION_INITIAL_TEMP,VEGETATION_LOAD,VEGETATION_LSET_IGNITE_TIME,VEGETATION_MOISTURE,VEGETATION_SVRATIO, &
            FIRELINE_MLR_MAX,VEGETATION_GROUND_TEMP,VEG_LSET_ROS_HEAD,VEG_LSET_ROS_FLANK,VEG_LSET_ROS_BACK, &
            VEG_LSET_WIND_EXP
LOGICAL :: VEGETATION,VEGETATION_NO_BURN,VEGETATION_LINEAR_DEGRAD,VEGETATION_ARRHENIUS_DEGRAD,VEG_LEVEL_SET_SPREAD,&
           DEFAULT,EVAC_DEFAULT
NAMELIST /SURF/ ADIABATIC,BACKING,BOUNDARY_CONDITION_FILENAME,BURN_AWAY,CELL_SIZE_FACTOR,COLOR,CONVECTION_LENGTH_SCALE, &
                CONVECTIVE_HEAT_FLUX,DEFAULT,&
                DT_INSERT,EMISSIVITY,EMISSIVITY_BACK,EVAC_DEFAULT,EXTERNAL_FLUX,E_COEFFICIENT,&
                FIRELINE_MLR_MAX,FREE_SLIP,FYI,GEOMETRY,HEAT_OF_VAPORIZATION,HRRPUA,H_FIXED,ID,IGNITION_TEMPERATURE,LAYER_DIVIDE,&
                LEAK_PATH,LENGTH,MASS_FLUX,MASS_FLUX_TOTAL,MASS_FRACTION,MATL_ID,MATL_MASS_FRACTION,&
                MINIMUM_LAYER_THICKNESS,MLRPUA,NET_HEAT_FLUX,&
                NO_SLIP,NPPC,PARTICLE_MASS_FLUX,PART_ID,PLE,PROFILE,RADIUS,RAMP_EF,RAMP_MF,RAMP_PART,RAMP_Q,RAMP_T,RAMP_V,&
                REGRID_FACTOR,RGB,ROUGHNESS,SHRINK,SPEC_ID,SPREAD_RATE,STRETCH_FACTOR,SURFACE_DENSITY,&
                TAU_EF,TAU_MF,TAU_PART,TAU_Q,TAU_T,TAU_V,TEXTURE_HEIGHT,TEXTURE_MAP,TEXTURE_WIDTH,THICKNESS,&
                TMP_BACK,TMP_FRONT,TMP_INNER,TRANSPARENCY,VEGETATION,VEGETATION_ARRHENIUS_DEGRAD,VEGETATION_CDRAG,&
                VEGETATION_CHAR_FRACTION,VEGETATION_ELEMENT_DENSITY,VEGETATION_GROUND_TEMP,VEGETATION_HEIGHT,&
                VEGETATION_INITIAL_TEMP,VEGETATION_LAYERS,VEGETATION_LINEAR_DEGRAD,VEGETATION_LOAD,VEGETATION_LSET_IGNITE_TIME,&
                VEGETATION_MOISTURE,VEGETATION_NO_BURN,VEGETATION_SVRATIO,VEG_LEVEL_SET_SPREAD,VEG_LSET_ROS_BACK,&
                VEG_LSET_ROS_FLANK,VEG_LSET_ROS_HEAD,VEG_LSET_WIND_EXP,VEL,VEL_T,VOLUME_FLUX,WIDTH,XYZ,Z0,&
                EXTERNAL_FLUX_RAMP,TAU_EXTERNAL_FLUX ! Backwards compatability??
            
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
 
N_SURF_RESERVED = 9
ALLOCATE(SURFACE(0:N_SURF+N_SURF_RESERVED),STAT=IZERO)
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

! Add extra surface types to the list that has already been compiled
 
INERT_SURF_INDEX             = 0
OPEN_SURF_INDEX              = N_SURF + 1
MIRROR_SURF_INDEX            = N_SURF + 2
INTERPOLATED_SURF_INDEX      = N_SURF + 3
PERIODIC_SURF_INDEX          = N_SURF + 4
HVAC_SURF_INDEX              = N_SURF + 5
MASSLESS_PARTICLE_SURF_INDEX = N_SURF + 6
DROPLET_SURF_INDEX           = N_SURF + 7
VEGETATION_SURF_INDEX        = N_SURF + 8
EVACUATION_SURF_INDEX        = N_SURF + 9

N_SURF = N_SURF + N_SURF_RESERVED

SURFACE(INERT_SURF_INDEX)%ID             = 'INERT'
SURFACE(OPEN_SURF_INDEX)%ID              = 'OPEN'
SURFACE(MIRROR_SURF_INDEX)%ID            = 'MIRROR'
SURFACE(INTERPOLATED_SURF_INDEX)%ID      = 'INTERPOLATED'
SURFACE(PERIODIC_SURF_INDEX)%ID          = 'PERIODIC'
SURFACE(HVAC_SURF_INDEX)%ID              = 'HVAC'
SURFACE(MASSLESS_PARTICLE_SURF_INDEX)%ID = 'MASSLESS PARTICLE'
SURFACE(DROPLET_SURF_INDEX)%ID           = 'DROPLET'
SURFACE(VEGETATION_SURF_INDEX)%ID        = 'VEGETATION'
SURFACE(EVACUATION_SURF_INDEX)%ID        = 'EVACUATION_OUTFLOW'

SURFACE(0)%USER_DEFINED                               = .FALSE.
SURFACE(N_SURF-N_SURF_RESERVED+1:N_SURF)%USER_DEFINED = .FALSE.
 
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

   IF (DEFAULT) SURF_DEFAULT = TRIM(ID)
   IF (EVAC_DEFAULT) EVAC_SURF_DEFAULT = TRIM(ID)

   ! Vegetation parameters

   IF(VEGETATION) WFDS_BNDRYFUEL = .TRUE.
   
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
   IF (ABS(VEGETATION_GROUND_TEMP+99._EB)>ZERO_P) SF%VEG_GROUND_ZERO_RAD = .FALSE.
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
   SF%BC_FILENAME          = BOUNDARY_CONDITION_FILENAME
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
   SF%MASS_FRACTION        = 0._EB
   SF%MAX_PRESSURE         = MAX_PRESSURE
   SF%MINIMUM_LAYER_THICKNESS = MINIMUM_LAYER_THICKNESS
   SF%NRA                  = NUMBER_RADIATION_ANGLES
   SF%NSB                  = NUMBER_SPECTRAL_BANDS
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

   ! Error checking
   
   IF (ANY(MASS_FLUX>0._EB) .AND. ANY(MASS_FRACTION>0._EB))  THEN
      WRITE (MESSAGE,'(A,A,A)') 'ERROR: Problem with SURF:',TRIM(SF%ID),&
                                '. Cannot use both MASS_FLUX and MASS_FRACTION'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   IF (ANY(MASS_FLUX>0._EB) .AND. ABS(VEL)>ZERO_P)  THEN
      WRITE (MESSAGE,'(A,A,A)') 'ERROR: Problem with SURF:',TRIM(SF%ID),&
                                '. Cannot use both MASS_FLUX and VEL'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   
   IF (ANY(MASS_FRACTION<0._EB))  THEN
      WRITE (MESSAGE,'(A,A,A)') 'ERROR: Problem with SURF:',TRIM(SF%ID),&
                                '. Cannot use a negative MASS_FRACTION'
      CALL SHUTDOWN(MESSAGE)
   ENDIF   
   
   IF (ANY(MASS_FLUX>0._EB) .OR. ANY(MASS_FRACTION>0._EB)) THEN
      IF (SPEC_ID(1)=='null') THEN
         WRITE (MESSAGE,'(A,A,A)') 'ERROR: Problem with SURF:',TRIM(SF%ID),&
                                   '. Must define SPEC_ID when using MASS_FLUX or MASS_FRACTION'
         CALL SHUTDOWN(MESSAGE)
      ELSE 
         DO NN=1,MAX_SPECIES
            IF (TRIM(SPEC_ID(NN))=='null') EXIT
            DO NNN=0,N_TRACKED_SPECIES
               IF (TRIM(SPECIES_MIXTURE(NNN)%ID)==TRIM(SPEC_ID(NN))) THEN
                  SF%MASS_FLUX(NNN)    = MASS_FLUX(NN)
                  SF%MASS_FRACTION(NNN)= MASS_FRACTION(NN)
                  SF%TAU(NNN)          = TAU_MF(NN)
                  SF%RAMP_MF(NNN)      = RAMP_MF(NN)
                  EXIT
               ENDIF
               IF (NNN==N_TRACKED_SPECIES) THEN
                  WRITE(MESSAGE,'(A,A,A,A,A)') 'ERROR: Problem with SURF:',TRIM(SF%ID),' SPEC ',TRIM(SPEC_ID(NN)),' not found'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      IF (SUM(SF%MASS_FRACTION) > 1._EB) THEN
         WRITE (MESSAGE,'(A,A,A)') 'ERROR: Problem with SURF:',TRIM(SF%ID),'. SUM(MASS_FRACTION) > 1'
            CALL SHUTDOWN(MESSAGE)
      ELSE
         IF (SF%MASS_FRACTION(0)<ZERO_P .AND. SUM(SF%MASS_FRACTION(1:N_TRACKED_SPECIES)) > 0._EB) &
            SF%MASS_FRACTION(0) = 1._EB - SUM(SF%MASS_FRACTION(1:N_TRACKED_SPECIES))
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
   ALLOCATE(SF%LAYER_THICKNESS(MAX_LAYERS))
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
         DO NR=1,ML%N_REACTIONS
            DO NNN=1,ML%N_RESIDUE(NR)
               IF (ML%RESIDUE_MATL_NAME(NNN,NR) == 'null') CYCLE
               IF (ANY(NAME_LIST==ML%RESIDUE_MATL_NAME(NNN,NR))) CYCLE 
               N_LIST = N_LIST + 1
               IF (N_LIST>MAX_MATERIALS_TOTAL) CALL SHUTDOWN('ERROR: Too many materials in the surface.')
               NAME_LIST (N_LIST) = ML%RESIDUE_MATL_NAME(NNN,NR)
               INDEX_LIST(N_LIST) = ML%RESIDUE_MATL_INDEX(NNN,NR)
            ENDDO
         ENDDO 
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
      SF%TMP_INNER                             = TMP_INNER + TMPM
      IF (SF%TMP_INNER(1)>=0._EB) SF%TMP_FRONT = SF%TMP_INNER(1)
      SF%TMP_BACK                              = TMP_BACK + TMPM
      ALLOCATE(SF%N_LAYER_CELLS(SF%N_LAYERS))            ! The number of cells in each layer
      ALLOCATE(SF%MIN_DIFFUSIVITY(SF%N_LAYERS))          ! The smallest diffusivity of materials in each layer
      ALLOCATE(SF%MATL_NAME(SF%N_MATL))                  ! The list of all material names associated with the surface
      ALLOCATE(SF%MATL_INDEX(SF%N_MATL))                 ! The list of all material indices associated with the surface
      ALLOCATE(SF%RESIDUE_INDEX(SF%N_MATL,MAX_MATERIALS,MAX_REACTIONS))! Each material associated with the surface has a RESIDUE
   ELSE
      SF%TMP_FRONT                  = TMP_FRONT + TMPM
      SF%TMP_INNER                  = SF%TMP_FRONT
      SF%TMP_BACK                   = SF%TMP_FRONT
   ENDIF
   DO NN = 1,SF%N_LAYERS
      IF (TMP_INNER(NN)>= -TMPM) TMPMIN = MIN(TMPMIN,TMP_INNER(NN)+TMPM)
   ENDDO
   IF (TMP_FRONT >= -TMPM) TMPMIN = MIN(TMPMIN,TMP_FRONT+TMPM)
   IF (TMP_BACK >= -TMPM) TMPMIN = MIN(TMPMIN,TMP_BACK+TMPM)
   IF (ASSUMED_GAS_TEMPERATURE >= 0._EB) TMPMIN = MIN(TMPMIN,ASSUMED_GAS_TEMPERATURE)

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
      DO NR=1,ML%N_REACTIONS
         DO NRM=1,ML%N_RESIDUE(NR)
            DO NNN=1,SF%N_MATL
               IF (ML%RESIDUE_MATL_INDEX(NRM,NR)==SF%MATL_INDEX(NNN)) SF%RESIDUE_INDEX(NN,NRM,NR) = NNN
            ENDDO
         ENDDO
         IF (SUM(ML%NU_RESIDUE(:,NR))<=ZERO_P) SF%SHRINK = .TRUE.
      ENDDO
      IF (ML%PYROLYSIS_MODEL==PYROLYSIS_LIQUID) SF%SHRINK = .TRUE.
   ENDDO

   ! Check for SHRINKage

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
   IF (TRIM(SF%BC_FILENAME)/='null')   SF%THERMAL_BC_INDEX = SPECIFIED_TEMPERATURE_FROM_FILE

   ! Set convection length scale according to THICKNESS or RADIUS

   IF (SF%CONV_LENGTH < 0._EB) THEN
      IF (SF%THERMALLY_THICK) THEN
         RADIUS = SF%THICKNESS
      ELSE
         IF ((SF%GEOMETRY/=SURF_CARTESIAN)) THEN
            IF (SF%THERMALLY_THICK) THEN
               RADIUS = SF%THICKNESS
            ELSE
               IF (RADIUS < 0._EB) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: SURF ',TRIM(SF%ID),' needs a RADIUS'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
            ENDIF
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
BOUNDARY_CONDITION_FILENAME = 'null'
BURN_AWAY               = .FALSE.
CELL_SIZE_FACTOR        = 1.0
COLOR                   = 'null'
CONVECTIVE_HEAT_FLUX    = 0._EB
CONVECTION_LENGTH_SCALE = -1._EB
NET_HEAT_FLUX           = 1.E12_EB
DEFAULT                 = .FALSE.
DT_INSERT               = 0.01_EB
DUCT_PATH               = 0 
E_COEFFICIENT           = 0._EB
EMISSIVITY              = -1._EB
EMISSIVITY_DEFAULT      = 0.9_EB
EMISSIVITY_BACK         = -1._EB
EVAC_DEFAULT            = .FALSE.
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
MASS_FRACTION           = 0._EB
MATL_ID                 = 'null'
MATL_MASS_FRACTION      = 0._EB
MATL_MASS_FRACTION(:,1) = 1._EB
MAX_PRESSURE            = 1.E12_EB
MINIMUM_LAYER_THICKNESS = 1.E-6_EB
MLRPUA                  = 0._EB
NPPC                    = 1
PARTICLE_MASS_FLUX      = 0._EB
PART_ID                 = 'null'
PLE                     = 0.3_EB
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
INTEGER :: N,NSPC,NR,ILPC
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
 
PROCESS_SURF_LOOP: DO N=0,N_SURF
 
   SF => SURFACE(N)
 
   ! Get ramps for the surface mass fraction and flux

   DO NSPC=0,N_TRACKED_SPECIES
      IF (TRIM(SF%RAMP_MF(NSPC))/='null') THEN
         CALL GET_RAMP_INDEX(SF%RAMP_MF(NSPC),'TIME',NR)
         SF%RAMP_INDEX(NSPC) = NR
      ELSE
         IF (SF%TAU(NSPC) > 0._EB) SF%RAMP_INDEX(NSPC) = TANH_RAMP 
         IF (SF%TAU(NSPC) < 0._EB) SF%RAMP_INDEX(NSPC) = TSQR_RAMP 
      ENDIF
   ENDDO 
   
   ! Look for particle classes that use SURF for property info

   DO ILPC=1,N_LAGRANGIAN_CLASSES

      LPC=>LAGRANGIAN_PARTICLE_CLASS(ILPC)
      IF (LPC%SURF_ID==SF%ID) THEN
         LPC%SURF_INDEX = N

         IF (LPC%SURF_INDEX==DROPLET_SURF_INDEX           .OR. &
             LPC%SURF_INDEX==MASSLESS_PARTICLE_SURF_INDEX .OR. &
             LPC%SURF_INDEX==VEGETATION_SURF_INDEX) CYCLE

         SELECT CASE (SF%GEOMETRY)
            CASE(SURF_CARTESIAN)
               IF (SF%LENGTH <0._EB) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: SURF ',TRIM(SF%ID),' needs a LENGTH'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
               IF (SF%WIDTH <0._EB) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: SURF ',TRIM(SF%ID),' needs a WIDTH'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
            CASE(SURF_CYLINDRICAL)
               IF (SF%LENGTH <0._EB) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: SURF ',TRIM(SF%ID),' needs a LENGTH'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF               
         END SELECT      
      ENDIF
   ENDDO

ENDDO PROCESS_SURF_LOOP   
 
! If a particle class uses a SURF line, make sure the SURF ID exists

DO ILPC=1,N_LAGRANGIAN_CLASSES
   LPC=>LAGRANGIAN_PARTICLE_CLASS(ILPC)
   IF (LPC%SURF_INDEX<0) THEN
      WRITE(MESSAGE,'(A,A,A)') 'ERROR: SURF ',TRIM(LPC%SURF_ID),' not found'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
ENDDO

END SUBROUTINE PROC_SURF_1



SUBROUTINE PROC_SURF_2

! Go through the SURF types and process
 
INTEGER :: ILPC,N,NN,NNN,NL
REAL(EB) :: ADJUSTED_LAYER_DENSITY,R_L(0:MAX_LAYERS)
INTEGER  :: IVEG_L,IIVEG_L,I_FUEL,I_GRAD
REAL(EB) :: DETA_VEG,DZVEG_L,ETA_H,ETAFM_VEG,ETAFP_VEG 
LOGICAL :: BURNING,BLOWING,SUCKING

TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
 
PROCESS_SURF_LOOP: DO N=0,N_SURF
 
   SF => SURFACE(N)
   IF (SF%THERMALLY_THICK) ML => MATERIAL(SF%LAYER_MATL_INDEX(1,1))
   
   SELECT CASE(SF%GEOMETRY)
      CASE(SURF_CARTESIAN)    ; I_GRAD = 1
      CASE(SURF_CYLINDRICAL)  ; I_GRAD = 2
      CASE(SURF_SPHERICAL)    ; I_GRAD = 3
   END SELECT

   ! Particle Information
 
   SF%PART_INDEX = 0
   IF (SF%PART_ID/='null') THEN
      DO ILPC=1,N_LAGRANGIAN_CLASSES
         LPC=>LAGRANGIAN_PARTICLE_CLASS(ILPC)
         IF (LPC%ID==SF%PART_ID)  SF%PART_INDEX = ILPC
      ENDDO
      IF (SF%PART_INDEX==0) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: PART_ID '//TRIM(SF%PART_ID)//' not found'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      PARTICLE_FILE=.TRUE.
   ENDIF

  ! Determine if the surface is combustible/burning
 
   SF%PYROLYSIS_MODEL = PYROLYSIS_NONE
   BURNING  = .FALSE.
   DO NL=1,SF%N_LAYERS
      DO NN=1,SF%N_LAYER_MATL(NL)
         NNN = SF%LAYER_MATL_INDEX(NL,NN)
         ML => MATERIAL(NNN)   
         IF (ML%PYROLYSIS_MODEL/=PYROLYSIS_NONE) THEN
            SF%PYROLYSIS_MODEL = PYROLYSIS_MATERIAL
            SF%STRETCH_FACTOR(NL) = 1._EB
            IF (N_REACTIONS>0) THEN
               IF (REACTION(1)%FUEL_SMIX_INDEX>=0) THEN
                  IF (ANY(ML%NU_SPEC(REACTION(1)%FUEL_SMIX_INDEX,:)>0._EB))  THEN
                     BURNING = .TRUE.
                     SF%TAU(TIME_HEAT) = 0._EB
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDDO

   IF (SF%HRRPUA>0._EB .OR. SF%MLRPUA>0._EB) THEN
      IF (N_REACTIONS > 1) THEN
         WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' has HRRPUA or MLRPUA set and there is more than one reaction'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      BURNING = .TRUE.
      SF%PYROLYSIS_MODEL = PYROLYSIS_SPECIFIED
   ENDIF

   IF (BURNING .AND. N_REACTIONS==0) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' indicates burning, but there is no REAC line'
      CALL SHUTDOWN(MESSAGE)
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
         RN => REACTION(1)
         SF%MASS_FLUX(RN%FUEL_SMIX_INDEX) = SF%HRRPUA/RN%HOC_COMPLETE
      ENDIF
      IF (SF%MLRPUA>0._EB) THEN
         RN => REACTION(1)
         SF%MASS_FLUX(RN%FUEL_SMIX_INDEX) = SF%MLRPUA
      ENDIF
      ! Adjust burning rate according to the difference of heats of combustion
      I_FUEL = REACTION(1)%FUEL_SMIX_INDEX
      IF (SF%N_LAYERS > 0) THEN
         ML => MATERIAL(SF%MATL_INDEX(1))
         SF%ADJUST_BURN_RATE(I_FUEL) = ML%ADJUST_BURN_RATE(I_FUEL,1)
         SF%MASS_FLUX(I_FUEL)        = SF%MASS_FLUX(I_FUEL)/SF%ADJUST_BURN_RATE(I_FUEL)
      ENDIF
      SF%TAU(I_FUEL)        = SF%TAU(TIME_HEAT)
      SF%RAMP_MF(I_FUEL)    = SF%RAMP_Q
      SF%RAMP_INDEX(I_FUEL) = SF%RAMP_INDEX(TIME_HEAT) 
   ENDIF BURNING_IF

   ! Compute surface density

   IF (SF%SURFACE_DENSITY < 0._EB) THEN
      SF%SURFACE_DENSITY = 0._EB
      R_L(0) = SF%THICKNESS
      DO NL=1,SF%N_LAYERS
         ADJUSTED_LAYER_DENSITY = 0._EB
         MATL_LOOP:DO NN=1,SF%N_LAYER_MATL(NL)
            NNN = SF%LAYER_MATL_INDEX(NL,NN)
            ML => MATERIAL(NNN)
            ADJUSTED_LAYER_DENSITY = ADJUSTED_LAYER_DENSITY + SF%LAYER_MATL_FRAC(NL,NN)/ML%RHO_S
         ENDDO MATL_LOOP
         IF (ADJUSTED_LAYER_DENSITY > 0._EB) ADJUSTED_LAYER_DENSITY = 1./ADJUSTED_LAYER_DENSITY
         R_L(NL) = R_L(NL-1)-SF%LAYER_THICKNESS(NL)
         SF%SURFACE_DENSITY = SF%SURFACE_DENSITY + ADJUSTED_LAYER_DENSITY * &
            (R_L(NL-1)**I_GRAD-R_L(NL)**I_GRAD)/(REAL(I_GRAD,EB)*SF%THICKNESS**(I_GRAD-1))
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

   IF (ANY(SF%MASS_FRACTION>0._EB) .AND. (ANY(ABS(SF%MASS_FLUX)>ZERO_P) .OR. SF%PYROLYSIS_MODEL/= PYROLYSIS_NONE)) THEN
      WRITE(MESSAGE,'(A)') 'ERROR: SURF '//TRIM(SF%ID)//' cannot specify mass fraction with mass flux and/or pyrolysis'
      CALL SHUTDOWN(MESSAGE)
      ENDIF
   IF (ANY(SF%MASS_FRACTION>0._EB) .AND. SUCKING) THEN
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

   IF (BLOWING .OR. SUCKING)        SF%SPECIES_BC_INDEX = SPECIFIED_MASS_FRACTION
   IF (ANY(SF%MASS_FRACTION>0._EB)) SF%SPECIES_BC_INDEX = SPECIFIED_MASS_FRACTION
   IF (ANY(ABS(SF%MASS_FLUX)>ZERO_P) .OR. &
       SF%PYROLYSIS_MODEL==PYROLYSIS_MATERIAL) SF%SPECIES_BC_INDEX = SPECIFIED_MASS_FLUX

   IF (SF%SPECIES_BC_INDEX==SPECIFIED_MASS_FRACTION) THEN
      IF (ALL(ABS(SF%MASS_FRACTION)< ZERO_P)) SF%MASS_FRACTION(0:N_TRACKED_SPECIES) = SPECIES_MIXTURE(0:N_TRACKED_SPECIES)%ZZ0
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
      SF%SPECIES_BC_INDEX = INFLOW_OUTFLOW_MASS_FLUX
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
   IF (N==MASSLESS_PARTICLE_SURF_INDEX) THEN
      SF%NRA = 1
      SF%NSB = 1
   ENDIF
   IF (N==DROPLET_SURF_INDEX) THEN
      SF%NRA = 1
      SF%NSB = 1
   ENDIF
   IF (N==VEGETATION_SURF_INDEX) THEN
      SF%NRA = 1
      SF%NSB = 1
   ENDIF
   IF (N==EVACUATION_SURF_INDEX) THEN
      SF%THERMAL_BC_INDEX = INFLOW_OUTFLOW
      SF%SPECIES_BC_INDEX = SPECIFIED_MASS_FRACTION
      SF%SPECIFIED_NORMAL_VELOCITY = .TRUE.
      SF%FREE_SLIP = .TRUE.
      SF%VELOCITY_BC_INDEX = FREE_SLIP_BC
      SF%VEL                   = +0.000001_EB ! VEL
      SF%TAU(TIME_VELO)        = 0.1_EB ! TAU_V
      SF%RAMP_INDEX(TIME_VELO) = TANH_RAMP
   ENDIF

   ! Do not allow N_LAYERS or N_CELLS to be zero

   IF (.NOT.SF%THERMALLY_THICK) THEN
      SF%N_LAYERS = 1
      SF%N_CELLS  = 1
      SF%N_MATL   = 1
      ALLOCATE(SF%N_LAYER_CELLS(SF%N_LAYERS))
      ALLOCATE(SF%X_S(0:SF%N_CELLS))
      SF%X_S = 0._EB
      ALLOCATE(SF%RHO_0(0:SF%N_CELLS+1,SF%N_MATL))
      SF%RHO_0 = 0._EB
      SF%TMP_INNER(:) = TMPA
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

INTEGER :: SURF_INDEX,N,NL,II,IL,NN
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

SURF_GRID_LOOP: DO SURF_INDEX=0,N_SURF

   SF => SURFACE(SURF_INDEX)
   IF (SF%THERMAL_BC_INDEX /= THERMALLY_THICK) CYCLE SURF_GRID_LOOP

   ! Compute number of points per layer, and then sum up to get total points for the surface

   SF%N_CELLS = 0
   DO NL=1,SF%N_LAYERS
      SF%MIN_DIFFUSIVITY(NL) = 1000000._EB
      DO N = 1,SF%N_LAYER_MATL(NL) 
         ML => MATERIAL(SF%LAYER_MATL_INDEX(NL,N))
         SF%MIN_DIFFUSIVITY(NL) = MIN(SF%MIN_DIFFUSIVITY(NL),ML%DIFFUSIVITY)
      ENDDO
      CALL GET_N_LAYER_CELLS(SF%MIN_DIFFUSIVITY(NL),SF%LAYER_THICKNESS(NL),SF%STRETCH_FACTOR(NL), &
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
   ALLOCATE(SF%RHO_0(0:SF%N_CELLS+1,SF%N_MATL))

! Compute node coordinates 

   CALL GET_WALL_NODE_COORDINATES(SF%N_CELLS,SF%N_LAYERS,SF%N_LAYER_CELLS, &
         SMALLEST_CELL_SIZE(1:SF%N_LAYERS),SF%STRETCH_FACTOR(1:SF%N_LAYERS),SF%X_S)

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

! Initialize the material densities of the solid

   SF%RHO_0 = 0._EB

   DO II=0,SF%N_CELLS+1
      IL = SF%LAYER_INDEX(II)
      IF (SF%TMP_INNER(IL)<=0._EB) SF%TMP_INNER(IL) = TMPA
      DO NN=1,SF%N_LAYER_MATL(IL)
         DO N=1,SF%N_MATL
            IF (SF%LAYER_MATL_INDEX(IL,NN)==SF%MATL_INDEX(N)) &
               SF%RHO_0(II,N) = SF%LAYER_MATL_FRAC(IL,NN)*SF%LAYER_DENSITY(IL)
         ENDDO
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
ALLOCATE(R_S_NEW(0:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','R_S_NEW',IZERO)
ALLOCATE(DX_WGT_S(0:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','DX_WGT_S',IZERO)
ALLOCATE(LAYER_INDEX(0:NWP_MAX+1),STAT=IZERO)
CALL ChkMemErr('INIT','LAYER_INDEX',IZERO)
ALLOCATE(MF_FRAC(1:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','MF_FRAC',IZERO)
ALLOCATE(SHRINK_FACTOR(1:NWP_MAX),STAT=IZERO)
CALL ChkMemErr('INIT','SHRINK_FACTOR',IZERO)

END SUBROUTINE PROC_WALL


SUBROUTINE READ_PRES

USE SCRC, ONLY: SCARC_METHOD , SCARC_KRYLOV , SCARC_MULTIGRID, SCARC_SMOOTH  , SCARC_PRECON, &
                SCARC_COARSE , SCARC_INITIAL, SCARC_SYSTEM   , SCARC_ACCURACY, SCARC_DEBUG , &
                SCARC_MULTIGRID_CYCLE, SCARC_MULTIGRID_LEVEL , SCARC_MULTIGRID_COARSENING  , &
                SCARC_MULTIGRID_ITERATIONS, SCARC_MULTIGRID_ACCURACY, SCARC_MULTIGRID_INTERPOL, &
                SCARC_KRYLOV_ITERATIONS, SCARC_KRYLOV_ACCURACY, &
                SCARC_SMOOTH_ITERATIONS, SCARC_SMOOTH_ACCURACY, SCARC_SMOOTH_OMEGA, &
                SCARC_PRECON_ITERATIONS, SCARC_PRECON_ACCURACY, SCARC_PRECON_OMEGA, &
                SCARC_COARSE_ITERATIONS, SCARC_COARSE_ACCURACY

NAMELIST /PRES/ CHECK_POISSON,MAX_PRESSURE_ITERATIONS, &
                PRESSURE_RELAX_TIME,RELAXATION_FACTOR, &
                SCARC_METHOD , SCARC_KRYLOV , SCARC_MULTIGRID, SCARC_SMOOTH  , SCARC_PRECON, &
                SCARC_COARSE , SCARC_INITIAL, SCARC_SYSTEM   , SCARC_ACCURACY, SCARC_DEBUG , &
                SCARC_MULTIGRID_CYCLE, SCARC_MULTIGRID_LEVEL , SCARC_MULTIGRID_COARSENING  , &
                SCARC_MULTIGRID_ITERATIONS, SCARC_MULTIGRID_ACCURACY, SCARC_MULTIGRID_INTERPOL, &
                SCARC_KRYLOV_ITERATIONS, SCARC_KRYLOV_ACCURACY, &
                SCARC_SMOOTH_ITERATIONS, SCARC_SMOOTH_ACCURACY, SCARC_SMOOTH_OMEGA, &
                SCARC_PRECON_ITERATIONS, SCARC_PRECON_ACCURACY, SCARC_PRECON_OMEGA, &
                SCARC_COARSE_ITERATIONS, SCARC_COARSE_ACCURACY, VELOCITY_TOLERANCE

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

IF (SCARC_METHOD /= 'null') PRES_METHOD = 'SCARC'

!IF (VELOCITY_TOLERANCE>100._EB.OR.(VELOCITY_TOLERANCE==0._EB.AND.PRES_METHOD=='SCARC')) THEN
IF (VELOCITY_TOLERANCE>100._EB) THEN
   ITERATE_PRESSURE = .FALSE.
ELSE
   ITERATE_PRESSURE = .TRUE.
   IF (VELOCITY_TOLERANCE<ZERO_P) VELOCITY_TOLERANCE = 0.5_EB*MESHES(1)%CELL_SIZE
ENDIF

END SUBROUTINE READ_PRES


SUBROUTINE READ_RADI

USE RADCONS
NAMELIST /RADI/ ANGLE_INCREMENT,CH4_BANDS,KAPPA0,NMIEANG,NUMBER_RADIATION_ANGLES,PATH_LENGTH,RADCAL_FUEL,&
                RADIATION,RADIATIVE_FRACTION,RADTMP,TIME_STEP_INCREMENT,WIDE_BAND_MODEL
REAL(EB) THETALOW,THETAUP
INTEGER  NRA,N

! Set default values

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
 
! Read radiation parameters

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

! Calculate actual number of radiation angles and determine the angular discretization

IF (.NOT.RADIATION) THEN
   NUMBER_RADIATION_ANGLES = 1
ELSE

   NRA = NUMBER_RADIATION_ANGLES

! Determine the number of polar angles (theta)

   IF (CYLINDRICAL) THEN
      NRT = NINT(SQRT(REAL(NRA)))
   ELSEIF (TWO_D) THEN
      NRT = 1
   ELSE
      NRT = 2*NINT(0.5_EB*1.17*REAL(NRA)**(1._EB/2.26))
   ENDIF      
 
! Determine number of azimuthal angles (phi)
 
   ALLOCATE(NRP(1:NRT),STAT=IZERO)
   CALL ChkMemErr('INIT','NRP',IZERO)

   N = 0
   DO I=1,NRT
      IF (CYLINDRICAL) THEN
         NRP(I) = NINT(REAL(NRA)/(REAL(NRT)))
      ELSEIF (TWO_D) THEN
         NRP(I) = 4*NINT(0.25_EB*REAL(NRA))
      ELSE
         THETALOW = PI*REAL(I-1)/REAL(NRT)
         THETAUP  = PI*REAL(I)/REAL(NRT)
         NRP(I) = NINT(0.5_EB*REAL(NRA)*(COS(THETALOW)-COS(THETAUP)))
         NRP(I) = MAX(4,NRP(I))
         NRP(I) = 4*NINT(0.25_EB*REAL(NRP(I)))
      ENDIF
      N = N + NRP(I)
   ENDDO
   NUMBER_RADIATION_ANGLES = N 

ENDIF

END SUBROUTINE READ_RADI


SUBROUTINE READ_CLIP

REAL(EB) :: MAXIMUM_DENSITY,MINIMUM_DENSITY,MINIMUM_TEMPERATURE,MAXIMUM_TEMPERATURE
NAMELIST /CLIP/ FYI,MAXIMUM_DENSITY,MAXIMUM_TEMPERATURE,MINIMUM_DENSITY,MINIMUM_TEMPERATURE
 
! Check for user-defined mins and maxes.
 
MINIMUM_DENSITY       = -999._EB
MAXIMUM_DENSITY       = -999._EB
MINIMUM_TEMPERATURE   = -999._EB
MAXIMUM_TEMPERATURE   = -999._EB
 
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
 
END SUBROUTINE READ_CLIP
 
 
SUBROUTINE READ_RAMP
 
REAL(EB) :: X,T,F,TM
INTEGER  :: I,II,NN,N,NUMBER_INTERPOLATION_POINTS
CHARACTER(30) :: DEVC_ID,CTRL_ID
TYPE(RAMPS_TYPE), POINTER :: RP=>NULL()
NAMELIST /RAMP/ CTRL_ID,DEVC_ID,F,FYI,ID,NUMBER_INTERPOLATION_POINTS,T,X

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
   IF (RP%NUMBER_DATA_POINTS<2) THEN
      IF (RP%NUMBER_DATA_POINTS==0) WRITE(MESSAGE,'(A,A,A)') 'ERROR: RAMP ',TRIM(RAMP_ID(N)), ' not found'
      IF (RP%NUMBER_DATA_POINTS==1) WRITE(MESSAGE,'(A,A,A)') 'ERROR: RAMP ',TRIM(RAMP_ID(N)), ' has only one point'
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
 
REAL(EB) :: TABLE_DATA(9)
INTEGER  :: NN,N
TYPE(TABLES_TYPE), POINTER :: TA=>NULL()
NAMELIST /TABL/ FYI,ID,TABLE_DATA

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
      CASE (PART_RADIATIVE_PROPERTY)
         TA%NUMBER_COLUMNS = 3
   END SELECT
   SEARCH_LOOP: DO
      CALL CHECKREAD('TABL',LU_INPUT,IOS)
      IF (IOS==1) EXIT SEARCH_LOOP
      TABLE_DATA = -999._EB
      READ(LU_INPUT,NML=TABL,ERR=56,IOSTAT=IOS)
      IF (ID/=TABLE_ID(N)) CYCLE SEARCH_LOOP
      TA%NUMBER_ROWS = TA%NUMBER_ROWS + 1
      MESSAGE='null'
      SELECT CASE(TABLE_TYPE(N))
         CASE (SPRAY_PATTERN)
            IF (TABLE_DATA(1)<0._EB .OR.           TABLE_DATA(1)>180._EB) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad 1st lattitude'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(2)<TABLE_DATA(1).OR. TABLE_DATA(2)>180._EB) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad 2nd lattitude'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(3)<-180._EB .OR.        TABLE_DATA(3)>360._EB) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad 1st longitude'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(4)<TABLE_DATA(3).OR. TABLE_DATA(4)>360._EB) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad 2nd longitude'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(5)<0._EB) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad velocity'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(6)<0._EB) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad mass flow'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
         CASE (PART_RADIATIVE_PROPERTY)
            IF (TABLE_DATA(1)<0._EB) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad wave length'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(2)<=0._EB) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad real index'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            IF (TABLE_DATA(3)< 0._EB) THEN
               WRITE(MESSAGE,'(A,I5,A,A,A)') 'ERROR: Row ',TA%NUMBER_ROWS,' of ',TRIM(TABLE_ID(N)),' has a bad complex index'
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
TYPE(OBSTRUCTION_TYPE), POINTER :: OB2=>NULL(),OBT=>NULL()
TYPE(MULTIPLIER_TYPE), POINTER :: MR=>NULL()
TYPE(OBSTRUCTION_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: TEMP_OBSTRUCTION
INTEGER :: NM,NOM,N_OBST_O,NNN,IC,N,NN,NNNN,N_NEW_OBST,RGB(3),N_OBST_NEW,II,JJ,KK,EVAC_N
CHARACTER(30) :: ID,DEVC_ID,PROP_ID,SURF_ID,SURF_IDS(3),SURF_ID6(6),CTRL_ID,MULT_ID
CHARACTER(60) :: MESH_ID
CHARACTER(25) :: COLOR
LOGICAL :: EVACUATION_OBST
REAL(EB) :: TRANSPARENCY,XB1,XB2,XB3,XB4,XB5,XB6,BULK_DENSITY,VOL_ADJUSTED,VOL_SPECIFIED
LOGICAL :: SAWTOOTH,EMBEDDED,THICKEN,PERMIT_HOLE,ALLOW_VENT,EVACUATION, REMOVABLE,BNDF_FACE(-3:3),BNDF_OBST,OUTLINE,NOTERRAIN
NAMELIST /OBST/ ALLOW_VENT,BNDF_FACE,BNDF_OBST,BULK_DENSITY,COLOR,CTRL_ID,DEVC_ID,EVACUATION,FYI,ID,MESH_ID,MULT_ID,NOTERRAIN,&
                OUTLINE,PERMIT_HOLE,PROP_ID,REMOVABLE,RGB,SAWTOOTH,SURF_ID,SURF_ID6,SURF_IDS,TEXTURE_ORIGIN,THICKEN,&
                TRANSPARENCY,XB

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
 
   IF (EVACUATION_ONLY(NM)) CALL DEFINE_EVACUATION_OBSTS(NM,1,0)

   ! Allocate OBSTRUCTION array

   ALLOCATE(M%OBSTRUCTION(0:N_OBST),STAT=IZERO)
   CALL ChkMemErr('READ','OBSTRUCTION',IZERO)
   OBSTRUCTION=>M%OBSTRUCTION
 
   N        = 0
   N_OBST_O = N_OBST
   EVAC_N   = 1
 
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
      NOTERRAIN      = .FALSE.
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
      IF (EVACUATION_ONLY(NM)) CALL DEFINE_EVACUATION_OBSTS(NM,2,EVAC_N)
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
 
               EVAC_N = EVAC_N + 1
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
               
               OB%NOTERRAIN = NOTERRAIN

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
 
               OB%SURF_INDEX(:) = DEFAULT_SURF_INDEX
          
               NNNN = 0
               DO NNN=0,N_SURF
                  IF (SURF_ID    ==SURFACE(NNN)%ID) OB%SURF_INDEX(:)    = NNN
                  IF (SURF_IDS(1)==SURFACE(NNN)%ID) OB%SURF_INDEX(3)    = NNN
                  IF (SURF_IDS(2)==SURFACE(NNN)%ID) OB%SURF_INDEX(-2:2) = NNN
                  IF (SURF_IDS(3)==SURFACE(NNN)%ID) OB%SURF_INDEX(-3)   = NNN
                  IF (SURF_ID6(1)==SURFACE(NNN)%ID) OB%SURF_INDEX(-1)   = NNN
                  IF (SURF_ID6(2)==SURFACE(NNN)%ID) OB%SURF_INDEX( 1)   = NNN
                  IF (SURF_ID6(3)==SURFACE(NNN)%ID) OB%SURF_INDEX(-2)   = NNN
                  IF (SURF_ID6(4)==SURFACE(NNN)%ID) OB%SURF_INDEX( 2)   = NNN
                  IF (SURF_ID6(5)==SURFACE(NNN)%ID) OB%SURF_INDEX(-3)   = NNN
                  IF (SURF_ID6(6)==SURFACE(NNN)%ID) OB%SURF_INDEX( 3)   = NNN
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
         
               ! Determine if the OBST is CONSUMABLE
         
               FACE_LOOP: DO NNN=-3,3
                  IF (NNN==0) CYCLE FACE_LOOP
                  IF (SURFACE(OB%SURF_INDEX(NNN))%BURN_AWAY) THEN
                     OB%CONSUMABLE = .TRUE.
                     IF (.NOT.SAWTOOTH) THEN
                        IF (ID=='null')WRITE(MESSAGE,'(A,I5,A)')'ERROR: OBST ',N,       ' cannot be BURN_AWAY and SAWTOOTH=.FALSE.'
                        IF (ID/='null')WRITE(MESSAGE,'(A,A,A)') 'ERROR: OBST ',TRIM(ID),' cannot be BURN_AWAY and SAWTOOTH=.FALSE.'
                        CALL SHUTDOWN(MESSAGE)
                     ENDIF
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
               IF (ABS(OB%VOLUME_ADJUST)<ZERO_P .AND. OB%BULK_DENSITY>0._EB) OB%BULK_DENSITY = -1._EB
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

   ! EVAC_PROCESS needs only the allocations of SOLID, CELL_INDEX, OBST_INDEX_C, OBSTRUCTION, and I,J,K_CELL

   IF (.NOT.EVACUATION_ONLY(NM) .AND. MYID==EVAC_PROCESS) THEN
      DEALLOCATE(M%EXTERIOR)
      CYCLE MESH_LOOP_3
   ENDIF
 
ENDDO MESH_LOOP_3

CONTAINS

  SUBROUTINE DEFINE_EVACUATION_OBSTS(NM,IMODE,EVAC_N)
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
    INTEGER, INTENT(IN) :: NM, IMODE, EVAC_N
    ! Local variables
    INTEGER :: N, N_END, I1, I2, J1, J2
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
          IF (EMESH_EXITS(N)%I_OBST==EVAC_N .AND. (EMESH_EXITS(N)%IMESH==NM .OR. EMESH_EXITS(N)%MAINMESH==NM)) THEN
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
          IF (EMESH_STAIRS(N)%I_OBST==EVAC_N .AND. EMESH_STAIRS(N)%IMESH==NM) THEN
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
NAMELIST /HOLE/ COLOR,CTRL_ID,DEVC_ID,FYI,EVACUATION,MESH_ID,MULT_ID,RGB,TRANSPARENCY,XB
TYPE(OBSTRUCTION_TYPE), ALLOCATABLE, DIMENSION(:) :: TEMP_OBST
TYPE(MULTIPLIER_TYPE), POINTER :: MR=>NULL()
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
      ! DEVC_ID    = 'null'
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
    REAL(EB) :: TINY_X, TINY_Y, TINY_Z

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
          TINY_X = 0.01_EB*(EMESH_STAIRS(I)%XB(2)-EMESH_STAIRS(I)%XB(1))/EMESH_STAIRS(I)%IBAR
          TINY_Y = 0.01_EB*(EMESH_STAIRS(I)%XB(4)-EMESH_STAIRS(I)%XB(3))/EMESH_STAIRS(I)%JBAR
          TINY_Z = 0.01_EB
          XB(1) = XB(1)-TINY_X ; XB(2) = XB(2)+TINY_X
          XB(3) = XB(3)-TINY_Y ; XB(4) = XB(4)+TINY_Y
          XB(5) = XB(5)-TINY_Z ; XB(6) = XB(6)+TINY_Z
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
TYPE (MESH_TYPE), POINTER :: M=>NULL()
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
REAL(EB) :: SPREAD_RATE,TRANSPARENCY,XYZ(3),TMP_EXTERIOR,DYNAMIC_PRESSURE, &
            REYNOLDS_STRESS(3,3),L_EDDY,VEL_RMS,L_EDDY_IJ(3,3),UVW(3)
CHARACTER(30) :: ID,DEVC_ID,CTRL_ID,SURF_ID,PRESSURE_RAMP,TMP_EXTERIOR_RAMP
CHARACTER(60) :: MESH_ID
CHARACTER(25) :: COLOR
LOGICAL :: REJECT_VENT,EVACUATION,OUTLINE,EVACUATION_VENT
NAMELIST /VENT/ COLOR,CTRL_ID,DEVC_ID,DYNAMIC_PRESSURE,EVACUATION,FYI,ID,IOR,L_EDDY,L_EDDY_IJ,MB,MESH_ID,N_EDDY,OUTLINE,&
                PBX,PBY,PBZ,PRESSURE_RAMP,REYNOLDS_STRESS,RGB,SPREAD_RATE,SURF_ID,TEXTURE_ORIGIN,TMP_EXTERIOR,TMP_EXTERIOR_RAMP,&
                TRANSPARENCY,UVW,VEL_RMS,XB,XYZ
 
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
      TMP_EXTERIOR_RAMP = 'null'
      REJECT_VENT  = .FALSE.
      TEXTURE_ORIGIN = -999._EB
      OUTLINE      = .FALSE.
      DEVC_ID  = 'null'
      CTRL_ID  = 'null'
      EVACUATION = .FALSE.
      N_EDDY=0
      L_EDDY=0._EB
      L_EDDY_IJ=0._EB
      VEL_RMS=0._EB
      REYNOLDS_STRESS=0._EB
      UVW = -1.E12_EB
 
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
 
      IF (ABS(XB(3)-XB(4))<=SPACING(XB(4))  .AND. TWO_D .AND. NN<NVO-1) THEN
         IF (ID=='null')WRITE(MESSAGE,'(A,I4,A)')'ERROR: VENT ',NN,      ' cannot be specified on a y boundary in a 2D calculation'
         IF (ID/='null')WRITE(MESSAGE,'(A,A,A)') 'ERROR: VENT ',TRIM(ID),' cannot be specified on a y boundary in a 2D calculation'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
 
      IF (ABS(XB(1)-XB(2))>SPACING(XB(2))  .AND. ABS(XB(3)-XB(4))>SPACING(XB(4))  .AND.ABS(XB(5)-XB(6))>SPACING(XB(6)) ) THEN
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
      
      IF (ABS(XB(1)-XB(2))<=SPACING(XB(2)) ) VT%TOTAL_INPUT_AREA = (XB(4)-XB(3))*(XB(6)-XB(5))
      IF (ABS(XB(3)-XB(4))<=SPACING(XB(4)) ) VT%TOTAL_INPUT_AREA = (XB(2)-XB(1))*(XB(6)-XB(5))
      IF (ABS(XB(5)-XB(6))<=SPACING(XB(6)) ) VT%TOTAL_INPUT_AREA = (XB(2)-XB(1))*(XB(4)-XB(3))

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
         IF (ABS(XB(1)-XB(2))>SPACING(XB(2))  .AND. ABS(XB(3)-XB(4))>SPACING(XB(4)) ) THEN
            IF (ID=='null') WRITE(MESSAGE,'(A,I4,A)') 'ERROR: Evacuation VENT ',NN,      ' must be a vertical plane'
            IF (ID/='null') WRITE(MESSAGE,'(A,A,A)')  'ERROR: Evacuation VENT ',TRIM(ID),' must be a vertical plane'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      ENDIF

      IF (ABS(XB(1)-XB(2))<=SPACING(XB(2)) ) THEN
         IF (VT%J1==VT%J2 .OR. VT%K1==VT%K2) REJECT_VENT=.TRUE.
         IF (VT%I1>IBAR .OR. VT%I2<0)        REJECT_VENT=.TRUE.
      ENDIF
      IF (ABS(XB(3)-XB(4))<=SPACING(XB(4)) ) THEN
         IF (VT%I1==VT%I2 .OR. VT%K1==VT%K2) REJECT_VENT=.TRUE.
         IF (VT%J1>JBAR .OR. VT%J2<0)        REJECT_VENT=.TRUE.
      ENDIF
      IF (ABS(XB(5)-XB(6))<=SPACING(XB(6)) ) THEN
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
 
      IF (ABS(XB(1)-XB(2))<=SPACING(XB(2)) ) VT%INPUT_AREA = (XB(4)-XB(3))*(XB(6)-XB(5))
      IF (ABS(XB(3)-XB(4))<=SPACING(XB(4)) ) VT%INPUT_AREA = (XB(2)-XB(1))*(XB(6)-XB(5))
      IF (ABS(XB(5)-XB(6))<=SPACING(XB(6)) ) VT%INPUT_AREA = (XB(2)-XB(1))*(XB(4)-XB(3))
 
      ! Check the SURF_ID against the list of SURF's

      CALL CHECK_SURF_NAME(SURF_ID,EX)
      IF (.NOT.EX) THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: SURF_ID ',TRIM(SURF_ID),' not found'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      ! Assign SURF_INDEX, Index of the Boundary Condition

      VT%SURF_INDEX = DEFAULT_SURF_INDEX
      DO NNN=0,N_SURF
         IF (SURF_ID==SURFACE(NNN)%ID) VT%SURF_INDEX = NNN
      ENDDO

      IF (SURF_ID=='OPEN')                            VT%TYPE_INDICATOR =  2
      IF (SURF_ID=='MIRROR' .OR. SURF_ID=='PERIODIC') VT%TYPE_INDICATOR = -2
      IF ((MB/='null' .OR.  PBX>-1.E5_EB .OR. PBY>-1.E5_EB .OR. PBZ>-1.E5_EB) .AND. SURF_ID=='OPEN') VT%TYPE_INDICATOR = -2
 
      VT%BOUNDARY_TYPE = SOLID_BOUNDARY
      IF (VT%SURF_INDEX==OPEN_SURF_INDEX)     VT%BOUNDARY_TYPE = OPEN_BOUNDARY
      IF (VT%SURF_INDEX==MIRROR_SURF_INDEX)   VT%BOUNDARY_TYPE = MIRROR_BOUNDARY
      IF (VT%SURF_INDEX==PERIODIC_SURF_INDEX) VT%BOUNDARY_TYPE = PERIODIC_BOUNDARY
      IF (VT%SURF_INDEX==HVAC_SURF_INDEX)     VT%BOUNDARY_TYPE = HVAC_BOUNDARY
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
      
      ! Synthetic Eddy Method
      
      VT%N_EDDY = N_EDDY
      IF (L_EDDY>0._EB) THEN
         VT%SIGMA_IJ = L_EDDY
      ELSE
         VT%SIGMA_IJ = L_EDDY_IJ ! Modified SEM (Jarrin, Ch. 7)
      ENDIF
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
      IF (TMP_EXTERIOR_RAMP/='null') CALL GET_RAMP_INDEX(TMP_EXTERIOR_RAMP,'TIME',VT%TMP_EXTERIOR_RAMP_INDEX)

      VT%TEXTURE(:) = TEXTURE_ORIGIN(:)
      
      VT%UVW = UVW
      IF (ALL(VT%UVW > -1.E12_EB)) THEN
         VT%UVW = VT%UVW/SQRT(VT%UVW(1)**2+VT%UVW(2)**2+VT%UVW(3)**2)
      ENDIF

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
USE DEVICE_VARIABLES, ONLY: DEVICE_TYPE, DEVICE, N_DEVC
REAL(EB) :: TEMPERATURE,DENSITY,MASS_FRACTION(1:MAX_SPECIES),RR_SUM,ZZ_GET(0:N_TRACKED_SPECIES),MASS_PER_VOLUME, &
            MASS_PER_TIME,DT_INSERT,UVW(3),HRRPUV,XYZ(3),DX,DY,DZ
INTEGER  :: N,NN,NNN,II,JJ,KK,NS,NS2,NUMBER_INITIAL_PARTICLES,N_PARTICLES,N_INIT_NEW,N_INIT_READ,N_PARTICLES_PER_CELL
LOGICAL  :: CELL_CENTERED
EQUIVALENCE(NUMBER_INITIAL_PARTICLES,N_PARTICLES)
CHARACTER(30) :: ID,CTRL_ID,DEVC_ID,PART_ID,SHAPE,MULT_ID,SPEC_ID(1:MAX_SPECIES)
TYPE(INITIALIZATION_TYPE), POINTER :: IN=>NULL()
TYPE(MULTIPLIER_TYPE), POINTER :: MR=>NULL()
TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
TYPE(DEVICE_TYPE), POINTER :: DV
NAMELIST /INIT/ CELL_CENTERED,CTRL_ID,DENSITY,DEVC_ID,DT_INSERT,DX,DY,DZ,HRRPUV,ID,MASS_FRACTION,MASS_PER_TIME,MASS_PER_VOLUME,&
                MULT_ID,N_PARTICLES,N_PARTICLES_PER_CELL,PART_ID,SHAPE,SPEC_ID,TEMPERATURE,UVW,XB,XYZ,&
                NUMBER_INITIAL_PARTICLES !Backwards compatability

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

   ! Set default values

   CELL_CENTERED             = .FALSE.
   CTRL_ID                   = 'null'
   DENSITY                   = -1000._EB
   DEVC_ID                   = 'null'
   DT_INSERT                 = -1._EB
   DX                        =  0._EB
   DY                        =  0._EB
   DZ                        =  0._EB
   HRRPUV                    =  0._EB
   ID                        = 'null'
   MASS_FRACTION             =  0._EB
   MASS_PER_TIME             = -1._EB
   MASS_PER_VOLUME           = -1._EB
   MULT_ID                   = 'null'
   N_PARTICLES               =  0
   N_PARTICLES_PER_CELL     = 0
   PART_ID                   = 'null'
   SHAPE                     = 'BLOCK'
   SPEC_ID                   = 'null'
   TEMPERATURE               = -1000._EB
   UVW                       = 0._EB
   XB(1)                     = -1000000._EB
   XB(2)                     =  1000000._EB
   XB(3)                     = -1000000._EB
   XB(4)                     =  1000000._EB
   XB(5)                     = -1000000._EB
   XB(6)                     =  1000000._EB
   XYZ                       = -1000000._EB
 
   ! Read in the INIT lines

   CALL CHECKREAD('INIT',LU_INPUT,IOS)
   IF (IOS==1) EXIT INIT_LOOP
   READ(LU_INPUT,INIT) 

   ! Transform XYZ into XB if necessary

   IF (ANY(XYZ>-100000._EB)) THEN
      XB(1:2) = XYZ(1)
      XB(3:4) = XYZ(2)
      XB(5:6) = XYZ(3)
   ENDIF
 
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

            IN%CELL_CENTERED = CELL_CENTERED
            IN%DX            = DX
            IN%DY            = DY
            IN%DZ            = DZ
            IN%ID            = ID
            IN%CTRL_ID       = CTRL_ID
            IN%DEVC_ID       = DEVC_ID
            CALL SEARCH_CONTROLLER('INIT',IN%CTRL_ID,IN%DEVC_ID,IN%DEVC_INDEX,IN%CTRL_INDEX,N)
            IN%VOLUME        = (IN%X2-IN%X1)*(IN%Y2-IN%Y1)*(IN%Z2-IN%Z1)
            IN%TEMPERATURE   = TEMPERATURE + TMPM
            IN%DENSITY       = DENSITY
            IN%SHAPE         = SHAPE
            IN%HRRPUV        = HRRPUV*1000._EB
            IF (HRRPUV > ZERO_P) INIT_HRRPUV = .TRUE.
            IF (DENSITY > 0._EB) RHOMAX = MAX(RHOMAX,IN%DENSITY)

            SPEC_INIT_IF:IF (N_TRACKED_SPECIES > 0._EB .AND. ANY(MASS_FRACTION > ZERO_P)) THEN
               IF (SPEC_ID(1)=='null') THEN
                  WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR: Problem with INIT number ',N,'. SPEC_ID must be used with MASS_FRACTION'
                  CALL SHUTDOWN(MESSAGE)
               ELSE 
                  DO NS=1,MAX_SPECIES
                     IF (SPEC_ID(NS)=='null') EXIT
                     DO NS2=0,N_TRACKED_SPECIES
                        IF (NS2>0 .AND. TRIM(SPEC_ID(NS))==TRIM(SPECIES_MIXTURE(NS2)%ID)) THEN
                           IN%MASS_FRACTION(NS2)=MASS_FRACTION(NS)
                           EXIT
                        ENDIF
                        IF (NS2==N_TRACKED_SPECIES)  THEN
                           WRITE(MESSAGE,'(A,I3,A,A,A)') 'ERROR: Problem with INIT number ',N,' tracked species ',&
                              TRIM(SPEC_ID(NS)),' not found'
                              CALL SHUTDOWN(MESSAGE)
                        ENDIF
                     ENDDO   
                  ENDDO  

                  IF (SUM(IN%MASS_FRACTION) > 1._EB) THEN
                     WRITE(MESSAGE,'(A,I3,A,A)') 'ERROR: Problem with INIT number ',N,'. Sum of specified mass fractions > 1'
                     CALL SHUTDOWN(MESSAGE)
                  ENDIF

                  ZZ_GET(1:N_TRACKED_SPECIES) = IN%MASS_FRACTION(1:N_TRACKED_SPECIES)
                  CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RR_SUM)
               ENDIF            
            ELSE SPEC_INIT_IF
                IF (N_TRACKED_SPECIES > 0) IN%MASS_FRACTION(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
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
         
            IN%MASS_PER_TIME   = MASS_PER_TIME
            IN%MASS_PER_VOLUME = MASS_PER_VOLUME
            
            IF(N_PARTICLES_PER_CELL>0 .AND. N_PARTICLES>0) THEN
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: INIT ',N,' Cannot use both N_PARTICLES and N_PARTICLES_PER_CELL'
               CALL SHUTDOWN(MESSAGE)               
            ENDIF
            
            IN%N_PARTICLES     = N_PARTICLES
            IN%N_PARTICLES_PER_CELL = N_PARTICLES_PER_CELL

            IF ( IN%MASS_PER_VOLUME>0._EB .AND. IN%VOLUME<=ZERO_P) THEN
               WRITE(MESSAGE,'(A,I4,A)') 'ERROR: INIT ',N,' XB has no volume'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
            
            IN%DT_INSERT = DT_INSERT
            IF (DT_INSERT>0._EB) IN%SINGLE_INSERTION = .FALSE.

            ! Set up a clock to keep track of particle insertions

            ALLOCATE(IN%PARTICLE_INSERT_CLOCK(NMESHES),STAT=IZERO)
            CALL ChkMemErr('READ','PARTICLE_INSERT_CLOCK',IZERO)
            IN%PARTICLE_INSERT_CLOCK = T_BEGIN
         
            ! Assign an index to identify the particle class

            IF (PART_ID/='null') THEN
               DO NS=1,N_LAGRANGIAN_CLASSES
                  IF (PART_ID==LAGRANGIAN_PARTICLE_CLASS(NS)%ID) THEN
                     IN%PART_INDEX = NS
                     PARTICLE_FILE = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
               IF (IN%PART_INDEX<1) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: PART_ID ',TRIM(PART_ID),' does not exist'
                  CALL SHUTDOWN(MESSAGE)
               ENDIF
               LPC => LAGRANGIAN_PARTICLE_CLASS(IN%PART_INDEX)
               IN%N_PARTICLES = N_PARTICLES*MAX(1,LPC%N_ORIENTATION)
               IF (IN%MASS_PER_TIME>0._EB .OR. IN%MASS_PER_VOLUME>0._EB) THEN
                  IF (LPC%DENSITY < 0._EB) THEN
                     WRITE(MESSAGE,'(A,A,A)') 'INIT ERROR: PARTicle class ',TRIM(LPC%ID),' requires a density'
                     CALL SHUTDOWN(MESSAGE)
                  ENDIF
               ENDIF
            ENDIF
            
            ! Initial velocity components

            IN%U0 = UVW(1)
            IN%V0 = UVW(2)
            IN%W0 = UVW(3)
         
         ENDDO I_MULT_LOOP
      ENDDO J_MULT_LOOP
   ENDDO K_MULT_LOOP

ENDDO INIT_LOOP

! Check if there are any devices that refer to INIT lines

DEVICE_LOOP: DO NN=1,N_DEVC
   DV => DEVICE(NN)
   IF (DV%INIT_ID=='null') CYCLE
   DO I=1,N_INIT
      IN => INITIALIZATION(I)
      IF (IN%ID==DV%INIT_ID) CYCLE DEVICE_LOOP
   ENDDO 
   WRITE(MESSAGE,'(A,A,A)') 'ERROR: The INIT_ID for DEVC ',TRIM(DV%ID),' cannot be found.'
   CALL SHUTDOWN(MESSAGE)
ENDDO DEVICE_LOOP

! Rewind the input file and return

REWIND(LU_INPUT)

END SUBROUTINE READ_INIT


SUBROUTINE READ_ZONE
 
INTEGER, PARAMETER :: MAX_LEAK_PATHS=200
REAL(EB) :: LEAK_AREA(0:MAX_LEAK_PATHS)
INTEGER  :: N,NM,NN
LOGICAL :: SEALED,READ_ZONE_LINES
CHARACTER(30) :: ID
NAMELIST /ZONE/ ID,LEAK_AREA,XB
 
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
   IF (N > 1) THEN
      DO NN = 1,N-1
         IF(P_ZONE(NN)%LEAK_AREA(N)>0._EB) THEN
            IF(P_ZONE(N)%LEAK_AREA(NN) > 0._EB) THEN
               WRITE(MESSAGE,'(A,I3,A,I3)')  'ERROR: LEAK_AREA specified twice for ZONE ',N,' and ',NN
               CALL SHUTDOWN(MESSAGE)           
            ELSE
               P_ZONE(N)%LEAK_AREA(NN) = P_ZONE(NN)%LEAK_AREA(N)
            ENDIF
         ENDIF
      ENDDO
   ENDIF
 
ENDDO READ_ZONE_LOOP
REWIND(LU_INPUT)

END SUBROUTINE READ_ZONE

 
SUBROUTINE READ_DEVC

! Just read in the DEViCes and the store the info in DEVICE()

USE DEVICE_VARIABLES, ONLY: DEVICE_TYPE, DEVICE, N_DEVC, N_DEVC_TIME, N_DEVC_LINE,MAX_DEVC_LINE_POINTS, DEVC_PIPE_OPERATING
INTEGER  :: NN,NM,MESH_NUMBER,N_DEVC_READ,IOR,TRIP_DIRECTION,VELO_INDEX,POINTS,I_POINT,PIPE_INDEX
REAL(EB) :: DEPTH,ORIENTATION(3),ROTATION,SETPOINT,FLOWRATE,BYPASS_FLOWRATE,DELAY,XYZ(3),CONVERSION_FACTOR,SMOOTHING_FACTOR
CHARACTER(30) :: QUANTITY,PROP_ID,CTRL_ID,DEVC_ID,INIT_ID,SURF_ID,STATISTICS,PART_ID,MATL_ID,SPEC_ID,UNITS, &
                 DUCT_ID,NODE_ID(2),X_ID,Y_ID,Z_ID,NO_UPDATE_DEVC_ID,NO_UPDATE_CTRL_ID
LOGICAL :: INITIAL_STATE,LATCH,DRY,TIME_AVERAGED,EVACUATION,HIDE_COORDINATES,RELATIVE,OUTPUT
TYPE (DEVICE_TYPE), POINTER :: DV=>NULL()
NAMELIST /DEVC/ BYPASS_FLOWRATE,CONVERSION_FACTOR,CTRL_ID,DELAY,DEPTH,DEVC_ID,DRY,DUCT_ID,EVACUATION,FLOWRATE,FYI,&
                HIDE_COORDINATES,ID,INITIAL_STATE,INIT_ID,IOR,LATCH,MATL_ID,NODE_ID, &
                NO_UPDATE_DEVC_ID,NO_UPDATE_CTRL_ID,ORIENTATION,OUTPUT,PART_ID,PIPE_INDEX,POINTS,&
                PROP_ID,QUANTITY,&
                RELATIVE,ROTATION,SETPOINT,SMOOTHING_FACTOR,SPEC_ID,STATISTICS,SURF_ID,TIME_AVERAGED,TRIP_DIRECTION,UNITS,&
                VELO_INDEX,XB,XYZ,X_ID,Y_ID,Z_ID

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

   ! Error statement involving POINTS

   IF (POINTS>1 .AND. ANY(XB<-1.E5_EB) .AND. INIT_ID=='null') THEN
      WRITE(MESSAGE,'(A,A,A)')  'ERROR: DEVC ',TRIM(ID),' must have coordinates given in terms of XB'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (POINTS>1 .AND. STATISTICS/='null' .AND. STATISTICS/='RMS') THEN
      WRITE(MESSAGE,'(A,A,A)')  'ERROR: DEVC ',TRIM(ID),' cannot use POINTS>1 and STATISTICS at the same time'
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   IF (POINTS<2 .AND. STATISTICS=='RMS') THEN
      WRITE(MESSAGE,'(A,A,A)')  'ERROR: DEVC ',TRIM(ID),' cannot compute RMS STATISTICS. Set POINTS>1.'
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
         IF (XYZ(1) < -1.E5_EB .AND. DUCT_ID=='null' .AND. NODE_ID(1)=='null' .AND. INIT_ID=='null') THEN
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
            IF (ABS(XYZ(1)-M%XS)<ZERO_P) XYZ(1) = XYZ(1) + 0.01_EB*M%DX(1)
            IF (ABS(XYZ(1)-M%XF)<ZERO_P) XYZ(1) = XYZ(1) - 0.01_EB*M%DX(M%IBAR)
            IF (ABS(XYZ(2)-M%YS)<ZERO_P) XYZ(2) = XYZ(2) + 0.01_EB*M%DY(1)
            IF (ABS(XYZ(2)-M%YF)<ZERO_P) XYZ(2) = XYZ(2) - 0.01_EB*M%DY(M%JBAR)
            IF (ABS(XYZ(3)-M%ZS)<ZERO_P) XYZ(3) = XYZ(3) + 0.01_EB*M%DZ(1)
            IF (ABS(XYZ(3)-M%ZF)<ZERO_P) XYZ(3) = XYZ(3) - 0.01_EB*M%DZ(M%KBAR)
            MESH_NUMBER = NM
            BAD = .FALSE.
            EXIT MESH_LOOP
         ENDIF
      ENDDO MESH_LOOP

      ! Process EVAC meshes

      EVACUATION_MESH_LOOP: DO NM=1,NMESHES
         IF (.NOT.EVACUATION_ONLY(NM)) CYCLE EVACUATION_MESH_LOOP
         M=>MESHES(NM)
         IF (XYZ(1)>=M%XS .AND. XYZ(1)<=M%XF .AND. XYZ(2)>=M%YS .AND. XYZ(2)<=M%YF .AND. XYZ(3)>=M%ZS .AND. XYZ(3)<=M%ZF) THEN
            IF (BAD) MESH_NUMBER = NM
            IF (.NOT.BAD .AND. EVACUATION .AND. QUANTITY=='TIME' .AND. SETPOINT<=T_BEGIN) MESH_NUMBER = NM
            BAD = .FALSE.
            EXIT EVACUATION_MESH_LOOP
         ENDIF
      ENDDO EVACUATION_MESH_LOOP
   
      ! Make sure there is either a QUANTITY or PROP_ID for the DEVICE
   
      IF (QUANTITY=='null' .AND. PROP_ID=='null') THEN
         WRITE(MESSAGE,'(A,A,A)')  'ERROR: DEVC ',TRIM(ID),' must have either an output QUANTITY or PROP_ID'
         CALL SHUTDOWN(MESSAGE)
      ENDIF
   
      IF (BAD) THEN
         IF (DUCT_ID/='null' .OR. NODE_ID(1)/='null' .OR. INIT_ID/='null') THEN
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

      IF (QUANTITY=='TIME') OUTPUT = .FALSE. ! Don't print out clocks

      IF (POINTS==1 .AND. OUTPUT)      N_DEVC_TIME = N_DEVC_TIME + 1
      IF (POINTS>1 .AND. I_POINT==1)   N_DEVC_LINE = N_DEVC_LINE + 1
      
      ! Assign properties to the DEVICE array
   
      N_DEVC = N_DEVC + 1

      DV => DEVICE(N_DEVC)
      M  => MESHES(MESH_NUMBER)
   
      DV%RELATIVE         = RELATIVE
      DV%CONVERSION_FACTOR= CONVERSION_FACTOR
      DV%DEPTH            = DEPTH
      DV%IOR              = IOR
      DV%ID               = ID
      IF (POINTS>1)DV%LINE= N_DEVC_LINE
      DV%POINT            = I_POINT
      DV%MESH             = MESH_NUMBER
      DV%ORDINAL          = NN
      DV%ORIENTATION(1:3) = ORIENTATION(1:3)/SQRT(ORIENTATION(1)**2+ORIENTATION(2)**2+ORIENTATION(3)**2)
      DV%PROP_ID          = PROP_ID
      DV%CTRL_ID          = CTRL_ID   
      DV%DEVC_ID          = DEVC_ID   
      DV%CTRL_ID          = CTRL_ID   
      DV%SURF_ID          = SURF_ID            
      DV%PART_ID          = PART_ID            
      DV%MATL_ID          = MATL_ID            
      DV%SPEC_ID          = SPEC_ID            
      DV%DUCT_ID          = DUCT_ID 
      DV%INIT_ID          = INIT_ID 
      DV%NODE_ID          = NODE_ID   
      DV%QUANTITY         = QUANTITY
      DV%ROTATION         = ROTATION*TWOPI/360._EB
      DV%SETPOINT         = SETPOINT
      DV%LATCH            = LATCH
      DV%OUTPUT           = OUTPUT
      DV%TRIP_DIRECTION   = TRIP_DIRECTION
      DV%INITIAL_STATE    = INITIAL_STATE
      DV%CURRENT_STATE    = INITIAL_STATE
      DV%PRIOR_STATE      = INITIAL_STATE
      DV%FLOWRATE         = FLOWRATE
      DV%BYPASS_FLOWRATE  = BYPASS_FLOWRATE
      DV%SMOOTHING_FACTOR = SMOOTHING_FACTOR
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
      DV%PIPE_INDEX       = PIPE_INDEX
      DV%NO_UPDATE_DEVC_ID = NO_UPDATE_DEVC_ID
      DV%NO_UPDATE_CTRL_ID = NO_UPDATE_CTRL_ID

      IF (POINTS > 1) THEN
         IF (.NOT.HIDE_COORDINATES) THEN
            IF (ABS(XB(1)-XB(2))> SPACING(XB(2)) .AND. ABS(XB(3)-XB(4))<=SPACING(XB(4)) .AND. &
               ABS(XB(5)-XB(6))<=SPACING(XB(6))) DV%LINE_COORD_CODE = 1
            IF (ABS(XB(1)-XB(2))<=SPACING(XB(2)) .AND. ABS(XB(3)-XB(4))> SPACING(XB(4)) .AND. &
               ABS(XB(5)-XB(6))<=SPACING(XB(6))) DV%LINE_COORD_CODE = 2
            IF (ABS(XB(1)-XB(2))<=SPACING(XB(2)) .AND. ABS(XB(3)-XB(4))<=SPACING(XB(4)) .AND. &
               ABS(XB(5)-XB(6))> SPACING(XB(6))) DV%LINE_COORD_CODE = 3
            IF (ABS(XB(1)-XB(2))> SPACING(XB(2)) .AND. ABS(XB(3)-XB(4))> SPACING(XB(4)) .AND. &
               ABS(XB(5)-XB(6))<=SPACING(XB(6))) DV%LINE_COORD_CODE = 12
            IF (ABS(XB(1)-XB(2))> SPACING(XB(2)) .AND. ABS(XB(3)-XB(4))<=SPACING(XB(4)) .AND. &
               ABS(XB(5)-XB(6))> SPACING(XB(6))) DV%LINE_COORD_CODE = 13
            IF (ABS(XB(1)-XB(2))<=SPACING(XB(2)) .AND. ABS(XB(3)-XB(4))> SPACING(XB(4)) .AND. &
               ABS(XB(5)-XB(6))> SPACING(XB(6))) DV%LINE_COORD_CODE = 23
         ELSE
            DV%LINE_COORD_CODE = 0
         ENDIF
      ENDIF
    
   ENDDO POINTS_LOOP

   ! Coordinates for non-point devices

   IF ((XB(1)>-1.E5_EB.OR.STATISTICS/='null') .AND. POINTS==1) THEN
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
      IF (ABS(XB(1)-XB(2))<=SPACING(XB(2)) .AND. IOR==0) DV%IOR = 1
      IF (ABS(XB(3)-XB(4))<=SPACING(XB(4)) .AND. IOR==0) DV%IOR = 2
      IF (ABS(XB(5)-XB(6))<=SPACING(XB(6)) .AND. IOR==0) DV%IOR = 3
   ENDIF
   
ENDDO READ_DEVC_LOOP

ALLOCATE (DEVC_PIPE_OPERATING(MAXVAL(DEVICE%PIPE_INDEX)))
DEVC_PIPE_OPERATING = 0

REWIND(LU_INPUT)

CONTAINS

SUBROUTINE SET_DEVC_DEFAULTS

RELATIVE         = .FALSE.
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
INIT_ID          = 'null'
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
OUTPUT           = .TRUE.
POINTS           = 1
SETPOINT         = 1.E20_EB
SMOOTHING_FACTOR = 0._EB
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
PIPE_INDEX       = 1
NO_UPDATE_DEVC_ID= 'null'
NO_UPDATE_CTRL_ID= 'null'

END SUBROUTINE SET_DEVC_DEFAULTS

END SUBROUTINE READ_DEVC



SUBROUTINE READ_CTRL

! Just read in the ConTRoL parameters and store in the array CONTROL

USE CONTROL_VARIABLES
USE MATH_FUNCTIONS, ONLY : GET_RAMP_INDEX

LOGICAL :: INITIAL_STATE, LATCH, EVACUATION
INTEGER :: CYCLES,N,NC,TRIP_DIRECTION
REAL(EB) :: SETPOINT(2), DELAY, CYCLE_TIME,CONSTANT,PROPORTIONAL_GAIN,INTEGRAL_GAIN,DIFFERENTIAL_GAIN,TARGET_VALUE
CHARACTER(30) :: ID,FUNCTION_TYPE,INPUT_ID(40),RAMP_ID,ON_BOUND
TYPE (CONTROL_TYPE), POINTER :: CF=>NULL()
NAMELIST /CTRL/  CONSTANT,CYCLES,CYCLE_TIME,DELAY,DIFFERENTIAL_GAIN,EVACUATION,FUNCTION_TYPE,ID,INITIAL_STATE,INTEGRAL_GAIN,&
                 INPUT_ID,LATCH,N,ON_BOUND,PROPORTIONAL_GAIN,RAMP_ID,&
                 SETPOINT,TARGET_VALUE,TRIP_DIRECTION
 
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
   CF%CONSTANT        = CONSTANT
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
   CF%TRIP_DIRECTION  = TRIP_DIRECTION
   CF%PROPORTIONAL_GAIN =PROPORTIONAL_GAIN
   CF%INTEGRAL_GAIN = INTEGRAL_GAIN
   CF%DIFFERENTIAL_GAIN = DIFFERENTIAL_GAIN
   CF%TARGET_VALUE = TARGET_VALUE
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
      CASE('SUM')
         CF%CONTROL_INDEX = CF_SUM
      CASE('SUBTRACT')
         CF%CONTROL_INDEX = CF_SUBTRACT
      CASE('MULTIPLY')
         CF%CONTROL_INDEX = CF_MULTIPLY
      CASE('DIVIDE')
         CF%CONTROL_INDEX = CF_DIVIDE
      CASE('POWER')
         CF%CONTROL_INDEX = CF_POWER      
      CASE('PID')
         CF%CONTROL_INDEX = CF_PID
         IF (CF%TARGET_VALUE<-1.E30_EB) THEN
            WRITE(MESSAGE,'(A,I5,A)')  'ERROR: CTRL ',NC,', PID controller must be given a TARGET_VALUE'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      CASE DEFAULT
         WRITE(MESSAGE,'(A,I5,A)')  'ERROR: CTRL ',NC,' FUNCTION_TYPE not recognized'
         CALL SHUTDOWN(MESSAGE)
   END SELECT
   
ENDDO READ_CTRL_LOOP
REWIND(LU_INPUT)

CONTAINS

SUBROUTINE SET_CTRL_DEFAULTS
   CONSTANT      = -9.E30_EB
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
   TRIP_DIRECTION = 1
   PROPORTIONAL_GAIN = 1._EB
   INTEGRAL_GAIN     = 0._EB
   DIFFERENTIAL_GAIN = 0._EB
   TARGET_VALUE      = 0._EB
   
END SUBROUTINE SET_CTRL_DEFAULTS

END SUBROUTINE READ_CTRL



SUBROUTINE PROC_CTRL

! Process the CONTROL function parameters

USE CONTROL_VARIABLES
USE DEVICE_VARIABLES, ONLY: N_DEVC,DEVICE
LOGICAL :: CONSTANT_SPECIFIED
INTEGER :: NC,NN,NNN
TYPE (CONTROL_TYPE), POINTER :: CF=>NULL()

PROC_CTRL_LOOP: DO NC = 1, N_CTRL

   CF => CONTROL(NC)
   CONSTANT_SPECIFIED = .FALSE.
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
   SELECT CASE (CF%CONTROL_INDEX)
      CASE (CF_SUBTRACT,CF_DIVIDE,CF_POWER)
         IF (CF%N_INPUTS /= 2) THEN
            WRITE(MESSAGE,'(A,I5,A)')  'ERROR: CTRL ',NC,' must have at only two inputs'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      CASE (CF_SUM,CF_MULTIPLY)
      CASE DEFAULT
         IF (ANY(CF%INPUT_ID=='CONSTANT')) THEN
            WRITE(MESSAGE,'(A,I5,A)')  'ERROR: CTRL ',NC,' the INTPUT_ID of CONSTANT cannot be used'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
      END SELECT
   ALLOCATE (CF%INPUT(CF%N_INPUTS),STAT=IZERO)
   CALL ChkMemErr('READ','CF%INPUT',IZERO)
   ALLOCATE (CF%INPUT_TYPE(CF%N_INPUTS),STAT=IZERO)
   CALL ChkMemErr('READ','CF%INPUT_TYPE',IZERO)
   
   BUILD_INPUT: DO NN = 1, CF%N_INPUTS
      IF (CF%INPUT_ID(NN)=='CONSTANT') THEN
         IF (CONSTANT_SPECIFIED) THEN
            WRITE(MESSAGE,'(A,I5,A)')  'ERROR: CTRL ',NC,' can only specify one input as a constant value'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         IF (CF%CONSTANT < -8.E30_EB) THEN
            WRITE(MESSAGE,'(A,I5,A)')  'ERROR: CTRL ',NC,' has the INPUT_ID CONSTANT but no constant value was specified'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         CF%INPUT_TYPE(NN) = CONSTANT_INPUT
         CONSTANT_SPECIFIED = .TRUE.
         CYCLE BUILD_INPUT
      ENDIF
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
USE DEVICE_VARIABLES, ONLY : DEVICE_TYPE, DEVICE, N_DEVC, PROPERTY, PROPERTY_TYPE
 
INTEGER  :: N,NN,NNN,NM,QUANTITY_INDEX,MAXCELLS,I,J,K,I_DUM
REAL(EB) :: XX,YY,ZZ,XX1,YY1,ZZ1,DISTANCE,SCANDISTANCE,DX,DY,DZ
TYPE (DEVICE_TYPE),  POINTER :: DV=>NULL()
 
IF (N_DEVC==0) RETURN

! Set initial values for DEViCes

DEVICE(1:N_DEVC)%VALUE = 0._EB
DEVICE(1:N_DEVC)%TIME_INTERVAL = 0._EB

PROC_DEVC_LOOP: DO N=1,N_DEVC

   DV => DEVICE(N)
   
   ! Check for HVAC outputs with no HVAC inputs

   IF ((DV%DUCT_ID/='null' .OR. DV%NODE_ID(1)/='null') .AND. .NOT. HVAC_SOLVE) THEN
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
                              DV%Y_INDEX,DV%Z_INDEX,DV%PART_INDEX,DV%DUCT_INDEX,DV%NODE_INDEX(1),'DEVC', &
                              DV%QUANTITY,'null',DV%SPEC_ID,DV%PART_ID,DV%DUCT_ID,DV%NODE_ID(1))
                              
      IF (OUTPUT_QUANTITY(QUANTITY_INDEX)%INTEGRATED .AND. DV%X1<=-1.E6_EB) THEN
         WRITE(MESSAGE,'(3A)')  'ERROR: DEVC QUANTITY ',TRIM(DV%QUANTITY),' requires coordinates using XB'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (.NOT.OUTPUT_QUANTITY(QUANTITY_INDEX)%INTEGRATED .AND. DV%STATISTICS=='null' .AND. DV%X1>-1.E6_EB .AND. DV%LINE==0) THEN
         WRITE(MESSAGE,'(3A)')  'ERROR: DEVC QUANTITY ',TRIM(DV%QUANTITY),' requires coordinates using XYZ, not XB'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (QUANTITY_INDEX<0 .AND. DV%IOR==0 .AND. DV%STATISTICS=='null' .AND. DV%INIT_ID=='null') THEN
         WRITE(MESSAGE,'(A,A,A)') 'ERROR: Specify orientation of DEVC ',TRIM(DV%ID),' using the parameter IOR'
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (QUANTITY_INDEX < 0 .AND. (DV%STATISTICS=='MASS MEAN' .OR. DV%STATISTICS=='VOLUME MEAN' .OR. &
                                    DV%STATISTICS=='VOLUME INTEGRAL' .OR. DV%STATISTICS=='MASS INTEGRAL' .OR. &
                                    DV%STATISTICS=='AREA INTEGRAL') ) THEN
         WRITE(MESSAGE,'(A,A)') 'ERROR: Invalid STATISTICS specified for wall DEVC ',TRIM(DV%ID)
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (QUANTITY_INDEX > 0 .AND. DV%STATISTICS=='SURFACE INTEGRAL') THEN
         WRITE(MESSAGE,'(A,A)') 'ERROR: Invalid STATISTICS specified for gas DEVC ',TRIM(DV%ID)
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (QUANTITY_INDEX > 0 .AND. DV%STATISTICS/='null' .AND. DV%STATISTICS/='TIME INTEGRAL' .AND. &
          DV%STATISTICS/='RMS' .AND. DV%I1<0) THEN
         WRITE(MESSAGE,'(A,A)') 'ERROR: XB required when geometrical STATISTICS specified for gas DEVC ',TRIM(DV%ID)
         CALL SHUTDOWN(MESSAGE)
      ENDIF
      
      IF (TRIM(DV%QUANTITY)=='NODE PRESSURE DIFFERENCE') THEN
         CALL GET_QUANTITY_INDEX(DV%SMOKEVIEW_LABEL,DV%SMOKEVIEW_BAR_LABEL,QUANTITY_INDEX,I_DUM, &
                                 DV%Y_INDEX,DV%Z_INDEX,DV%PART_INDEX,DV%DUCT_INDEX,DV%NODE_INDEX(2),'DEVC', &
                                 DV%QUANTITY,'null',DV%SPEC_ID,DV%PART_ID,DV%DUCT_ID,DV%NODE_ID(2))
         IF (DV%NODE_INDEX(1)==DV%NODE_INDEX(2)) THEN
            WRITE(MESSAGE,'(A,A)') 'ERROR: NODE PRESSURE DIFFERENCE node 1 = node 2 ',TRIM(DV%ID)
            CALL SHUTDOWN(MESSAGE)
         ENDIF     
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
   
   ! Initialize histogram
   IF(PROPERTY(DV%PROP_INDEX)%PDPA_HISTOGRAM) THEN
      ALLOCATE(DV%PDPA_HISTOGRAM_COUNTS(PROPERTY(DV%PROP_INDEX)%PDPA_HISTOGRAM_NBINS))
      DV%PDPA_HISTOGRAM_COUNTS(:)=0._EB
   ENDIF
   ! Do initialization of special models
   
   SPECIAL_QUANTITIES: SELECT CASE (DV%QUANTITY)

      CASE ('CHAMBER OBSCURATION') 

         IF (DV%PROP_INDEX<1) THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' is a smoke detector and must have a PROP_ID'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         IF (PROPERTY(DV%PROP_INDEX)%Y_INDEX<0 .AND. PROPERTY(DV%PROP_INDEX)%Z_INDEX<0) THEN
            IF (SOOT_INDEX<1) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' is a smoke detector and requires a smoke source'
               CALL SHUTDOWN(MESSAGE)
            ELSE
               PROPERTY(DV%PROP_INDEX)%Y_INDEX = SOOT_INDEX
            ENDIF
         ENDIF
         ALLOCATE(DV%T_E(-1:1000))
         ALLOCATE(DV%Y_E(-1:1000))
         DV%T_E      = T_BEGIN - M%DT
         DV%Y_E      = 0._EB
         DV%N_T_E    = -1
         DV%Y_C      = 0._EB
         DV%SETPOINT = PROPERTY(DV%PROP_INDEX)%ACTIVATION_OBSCURATION
         IF (PROPERTY(DV%PROP_INDEX)%Y_INDEX>0) DV%Y_INDEX = PROPERTY(DV%PROP_INDEX)%Y_INDEX
         IF (PROPERTY(DV%PROP_INDEX)%Z_INDEX>0) DV%Z_INDEX = PROPERTY(DV%PROP_INDEX)%Z_INDEX
   
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

      CASE ('SOLID DENSITY')

         IF (DV%MATL_ID=='null') THEN
            WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' must have a MATL_ID'
            CALL SHUTDOWN(MESSAGE)
         ENDIF

      CASE ('LAYER HEIGHT','UPPER TEMPERATURE','LOWER TEMPERATURE') 

         DV%K1 = MAX(1     ,DV%K1)
         DV%K2 = MIN(M%KBAR,DV%K2)

      CASE ('PATH OBSCURATION')

         IF (DV%PROP_INDEX>0) THEN
            IF (PROPERTY(DV%PROP_INDEX)%Y_INDEX<1 .AND. PROPERTY(DV%PROP_INDEX)%Z_INDEX<1) THEN
               IF (SOOT_INDEX<1) THEN
                  WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' is a smoke detector and requires a smoke source'
                  CALL SHUTDOWN(MESSAGE)
               ELSE
                  PROPERTY(DV%PROP_INDEX)%Y_INDEX = SOOT_INDEX
               ENDIF
            ENDIF
         ELSE
            IF (SOOT_INDEX <=0) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' is a smoke detector and requires a smoke source'
               CALL SHUTDOWN(MESSAGE)
            ENDIF
         ENDIF
         IF (PROPERTY(DV%PROP_INDEX)%Y_INDEX>0) DV%Y_INDEX = PROPERTY(DV%PROP_INDEX)%Y_INDEX
         IF (PROPERTY(DV%PROP_INDEX)%Z_INDEX>0) DV%Z_INDEX = PROPERTY(DV%PROP_INDEX)%Z_INDEX
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

      CASE ('ASPIRATION')

         ! Check either for a specified SMOKE SPECies, or if simple chemistry model is being used
         IF (DV%PROP_INDEX>0) THEN
            IF (PROPERTY(DV%PROP_INDEX)%Y_INDEX<1 .AND. PROPERTY(DV%PROP_INDEX)%Z_INDEX<1 .AND. SOOT_INDEX<1) THEN
               WRITE(MESSAGE,'(A,A,A)') 'ERROR: DEVC ',TRIM(DV%ID),' is a smoke detector and requires a smoke source'
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
         
      CASE ('VELOCITY PATCH') ! error statements to come
         PATCH_VELOCITY = .TRUE.
         ALLOCATE(DV%DEVC_INDEX(1),STAT=IZERO)
         DV%DEVC_INDEX(1) = 0
         DO NN=1,N_DEVC
            IF (DEVICE(NN)%ID==DV%DEVC_ID) DV%DEVC_INDEX(1) = NN
         ENDDO
         IF (DV%DEVC_INDEX(1)==0) THEN
            WRITE(MESSAGE,'(A)') 'ERROR: A VELOCITY PATCH DEVC line needs a DEVC_ID to control it'
            CALL SHUTDOWN(MESSAGE)
         ENDIF
         
   END SELECT SPECIAL_QUANTITIES

   IF (DV%STATISTICS/='null') CALL CHANGE_UNITS(DV%QUANTITY,DV%UNITS,DV%STATISTICS,MYID,LU_ERR)

   IF (DV%NO_UPDATE_DEVC_ID/='null' .OR. DV%NO_UPDATE_CTRL_ID/='null') &
      CALL SEARCH_CONTROLLER('DEVC',DV%NO_UPDATE_CTRL_ID,DV%NO_UPDATE_DEVC_ID,DV%NO_UPDATE_DEVC_INDEX,DV%NO_UPDATE_CTRL_INDEX,N)

ENDDO PROC_DEVC_LOOP

END SUBROUTINE PROC_DEVC


SUBROUTINE READ_PROF
 
INTEGER :: N,NM,MESH_NUMBER,NN,N_PROFO,IOR
REAL(EB) :: XYZ(3)
CHARACTER(30) :: QUANTITY
TYPE (PROFILE_TYPE), POINTER :: PF=>NULL()
NAMELIST /PROF/ FYI,ID,IOR,QUANTITY,XYZ
 
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
CHARACTER(30) :: QUANTITY,SPEC_ID,COLOR_SPEC_ID
INTEGER :: REDUCE_TRIANGLES,N,I_DUM,VELO_INDEX
TYPE(ISOSURFACE_FILE_TYPE), POINTER :: IS=>NULL()
NAMELIST /ISOF/ COLOR_SPEC_ID,FYI,QUANTITY,REDUCE_TRIANGLES,SPEC_ID,VALUE,VELO_INDEX
 
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
   SPEC_ID          = 'null'
   COLOR_SPEC_ID    = 'null'
   VALUE            = -999._EB
   REDUCE_TRIANGLES = 1
   VELO_INDEX       = 0
 
   CALL CHECKREAD('ISOF',LU_INPUT,IOS) 
   IF (IOS==1) EXIT READ_ISOF_LOOP
   READ(LU_INPUT,ISOF) 
 
   IS%REDUCE_TRIANGLES = REDUCE_TRIANGLES
   IS%VELO_INDEX       = VELO_INDEX

   CALL GET_QUANTITY_INDEX(IS%SMOKEVIEW_LABEL,IS%SMOKEVIEW_BAR_LABEL,IS%INDEX,I_DUM, &
                           IS%Y_INDEX,IS%Z_INDEX,I_DUM,I_DUM,I_DUM,'ISOF', &
                           QUANTITY,'null',SPEC_ID,'null','null','null')
                                              
   VALUE_LOOP: DO I=1,10
      IF (VALUE(I)<=-998._EB) EXIT VALUE_LOOP
      IS%N_VALUES = I
      IS%VALUE(I) = VALUE(I)
   ENDDO VALUE_LOOP
 
ENDDO READ_ISOF_LOOP

REWIND(LU_INPUT)
 
END SUBROUTINE READ_ISOF
 
 
SUBROUTINE READ_SLCF

USE EVAC, ONLY: EVAC_FDS6
REAL(EB) :: MAXIMUM_VALUE,MINIMUM_VALUE
REAL(EB) :: AGL_SLICE
INTEGER :: N,NN,NM,MESH_NUMBER,N_SLCF_O,NITER,ITER,VELO_INDEX,I_DUM
LOGICAL :: VECTOR,CELL_CENTERED, FIRE_LINE, EVACUATION,LEVEL_SET_FIRE_LINE
CHARACTER(30) :: QUANTITY,SPEC_ID,PART_ID,QUANTITY2
TYPE (SLICE_TYPE), POINTER :: SL=>NULL()
NAMELIST /SLCF/ AGL_SLICE,CELL_CENTERED,EVACUATION,FIRE_LINE,FYI,ID,LEVEL_SET_FIRE_LINE,MAXIMUM_VALUE,MESH_NUMBER,MINIMUM_VALUE,&
                PART_ID,PBX,PBY,PBZ,QUANTITY,QUANTITY2,SPEC_ID,VECTOR,VELO_INDEX,XB

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
      ID       = 'null'
      MESH_NUMBER=NM
      MINIMUM_VALUE = 0._EB
      MAXIMUM_VALUE = 0._EB
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
         SL%ID = ID
         SL%I1 = NINT( GINV(XB(1)-XS,1,NM)*RDXI)
         SL%I2 = NINT( GINV(XB(2)-XS,1,NM)*RDXI)
         SL%J1 = NINT( GINV(XB(3)-YS,2,NM)*RDETA)
         SL%J2 = NINT( GINV(XB(4)-YS,2,NM)*RDETA)
         SL%K1 = NINT( GINV(XB(5)-ZS,3,NM)*RDZETA)
         SL%K2 = NINT( GINV(XB(6)-ZS,3,NM)*RDZETA)
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
                                 SL%Y_INDEX,SL%Z_INDEX,SL%PART_INDEX,I_DUM,I_DUM,'SLCF', &
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

         ! IF (QUANTITY=='VELOCITY'   .OR. &
         !    QUANTITY=='U-VELOCITY' .OR. &
         !    QUANTITY=='V-VELOCITY' .OR. &
         !    QUANTITY=='W-VELOCITY') THEN
         !    CELL_CENTERED = .FALSE.
         ! ENDIF
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
NAMELIST /BNDF/ CELL_CENTERED,FYI,PART_ID,PROP_ID,QUANTITY,SPEC_ID
TYPE(BOUNDARY_FILE_TYPE), POINTER :: BF=>NULL()
 
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
                           BF%Y_INDEX,BF%Z_INDEX,BF%PART_INDEX,I_DUM,I_DUM,'BNDF', &
                           QUANTITY,'null',SPEC_ID,PART_ID,'null','null')
                           
   ! Assign miscellaneous attributes to the boundary file

   BF%CELL_CENTERED = CELL_CENTERED

   ! Check to see if PROP_ID exists

   BF%PROP_INDEX = 0
   IF(PROP_ID/='null')  CALL GET_PROPERTY_INDEX(BF%PROP_INDEX,'BNDF',PROP_ID)

ENDDO READ_BNDF_LOOP
REWIND(LU_INPUT)
 
END SUBROUTINE READ_BNDF


SUBROUTINE READ_BNDE

USE DEVICE_VARIABLES
INTEGER :: N,I_DUM
LOGICAL :: CELL_CENTERED
CHARACTER(30) :: QUANTITY,PROP_ID,SPEC_ID,PART_ID
NAMELIST /BNDE/ CELL_CENTERED,FYI,PART_ID,PROP_ID,QUANTITY,SPEC_ID
TYPE(BOUNDARY_ELEMENT_FILE_TYPE), POINTER :: BE=>NULL()
 
N_BNDE = 0
REWIND(LU_INPUT)
COUNT_BNDE_LOOP: DO
   CALL CHECKREAD('BNDE',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_BNDE_LOOP
   READ(LU_INPUT,NML=BNDE,END=309,ERR=310,IOSTAT=IOS)
   N_BNDE = N_BNDE + 1
   310 IF (IOS>0) THEN 
         WRITE(MESSAGE,'(A,I3)') 'ERROR: Problem with BNDE number ',N_BNDE+1
         CALL SHUTDOWN(MESSAGE)
       ENDIF
ENDDO COUNT_BNDE_LOOP
309 REWIND(LU_INPUT)
 
ALLOCATE(BOUNDARY_ELEMENT_FILE(N_BNDE),STAT=IZERO)
CALL ChkMemErr('READ','BOUNDARY_ELEMENT_FILE',IZERO)
 
READ_BNDE_LOOP: DO N=1,N_BNDE
   BE => BOUNDARY_ELEMENT_FILE(N)
   CELL_CENTERED = .TRUE.
   PART_ID  = 'null'
   PROP_ID  = 'null'
   SPEC_ID  = 'null'
   QUANTITY = 'WALL_TEMPERATURE'
   CALL CHECKREAD('BNDE',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_BNDE_LOOP
   READ(LU_INPUT,BNDE)
   
   ! Look to see if output QUANTITY exists

   CALL GET_QUANTITY_INDEX(BE%SMOKEVIEW_LABEL,BE%SMOKEVIEW_BAR_LABEL,BE%INDEX,I_DUM, &
                           BE%Y_INDEX,BE%Z_INDEX,BE%PART_INDEX,I_DUM,I_DUM,'BNDE', &
                           QUANTITY,'null',SPEC_ID,PART_ID,'null','null')
                           
   ! Assign miscellaneous attributes to the boundary file

   BE%CELL_CENTERED = CELL_CENTERED

   ! Check to see if PROP_ID exists

   BE%PROP_INDEX = 0
   IF(PROP_ID/='null')  CALL GET_PROPERTY_INDEX(BE%PROP_INDEX,'BNDE',PROP_ID)

ENDDO READ_BNDE_LOOP
REWIND(LU_INPUT)
 
END SUBROUTINE READ_BNDE
 

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
                              Y_INDEX,Z_INDEX,PART_INDEX,DUCT_INDEX,NODE_INDEX,OUTTYPE, &
                              QUANTITY,QUANTITY2,SPEC_ID_IN,PART_ID,DUCT_ID,NODE_ID)
CHARACTER(*), INTENT(INOUT) :: QUANTITY
CHARACTER(*), INTENT(OUT) :: SMOKEVIEW_LABEL,SMOKEVIEW_BAR_LABEL
CHARACTER(*) :: SPEC_ID_IN,PART_ID,DUCT_ID,NODE_ID
CHARACTER(30) :: SPEC_ID
CHARACTER(*), INTENT(IN) :: OUTTYPE,QUANTITY2
INTEGER, INTENT(OUT) :: OUTPUT_INDEX,Y_INDEX,Z_INDEX,PART_INDEX,DUCT_INDEX,NODE_INDEX,OUTPUT2_INDEX
INTEGER :: ND,NS,NN

! Backward compatibility

IF (QUANTITY=='oxygen') THEN
   QUANTITY    = 'VOLUME FRACTION'
   SPEC_ID_IN  = 'OXYGEN'
ENDIF
IF (QUANTITY=='carbon monoxide') THEN
   QUANTITY    = 'VOLUME FRACTION'
   SPEC_ID_IN  = 'CARBON MONOXIDE'
ENDIF
IF (QUANTITY=='carbon dioxide') THEN
   QUANTITY    = 'VOLUME FRACTION'
   SPEC_ID_IN  = 'CARBON DIOXIDE'
ENDIF
IF (QUANTITY=='soot') THEN
   QUANTITY    = 'VOLUME FRACTION'
   SPEC_ID_IN  = 'SOOT'
ENDIF
IF (QUANTITY=='soot density') THEN
   QUANTITY    = 'DENSITY'
   SPEC_ID_IN  = 'SOOT'
ENDIF
IF (QUANTITY=='fuel') THEN
   QUANTITY    = 'VOLUME FRACTION'
   WRITE(SPEC_ID_IN,'(A)') REACTION(1)%FUEL
ENDIF

DO ND=-N_OUTPUT_QUANTITIES,N_OUTPUT_QUANTITIES
   IF (QUANTITY==OUTPUT_QUANTITY(ND)%OLD_NAME) QUANTITY = OUTPUT_QUANTITY(ND)%NAME
ENDDO

! Initialize indices

Y_INDEX = -1
Z_INDEX = -1

SPEC_ID = SPEC_ID_IN

IF (QUANTITY=='OPTICAL DENSITY'        .AND. SPEC_ID=='null') SPEC_ID='SOOT'
IF (QUANTITY=='EXTINCTION COEFFICIENT' .AND. SPEC_ID=='null') SPEC_ID='SOOT'
IF (QUANTITY=='SOOT VOLUME FRACTION'   .AND. SPEC_ID=='null') SPEC_ID='SOOT'
IF (QUANTITY=='VISIBILITY'             .AND. SPEC_ID=='null') SPEC_ID='SOOT'

PART_INDEX = 0
DUCT_INDEX = 0
NODE_INDEX = 0
OUTPUT2_INDEX = 0

! Look for the appropriate SPEC or SMIX index

IF (SPEC_ID/='null') THEN
   CALL GET_SPEC_OR_SMIX_INDEX(SPEC_ID,Y_INDEX,Z_INDEX)
   IF (Z_INDEX<0 .AND. Y_INDEX<1) THEN
      WRITE(MESSAGE,'(A,A,A)')  'ERROR: ',TRIM(SPEC_ID),' is not explicitly specified'
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

IF (TRIM(QUANTITY)=='FILTER LOADING') THEN
   Y_INDEX = -999
   DO NS = 1,N_TRACKED_SPECIES
      IF (TRIM(SPECIES_MIXTURE(NS)%ID)==TRIM(SPEC_ID)) THEN
         Z_INDEX = NS
         EXIT
      ENDIF
   ENDDO
   IF (Z_INDEX<0) THEN
      WRITE(MESSAGE,'(A,A,A)')  'ERROR: FILTER LOADING. ',TRIM(SPEC_ID),' is not a tracked species'
      CALL SHUTDOWN(MESSAGE)
   ENDIF  
ENDIF

! Assign PART_INDEX when PART_ID is specified

IF (PART_ID/='null') THEN
   DO NS=1,N_LAGRANGIAN_CLASSES
      IF (PART_ID==LAGRANGIAN_PARTICLE_CLASS(NS)%ID) THEN
         PART_INDEX = NS
         EXIT
      ENDIF
   ENDDO
ENDIF

! Loop over all possible output quantities and assign an index number to match the desired QUANTITY

DO ND=-N_OUTPUT_QUANTITIES,N_OUTPUT_QUANTITIES
   IF (OUTPUT_QUANTITY(ND)%NAME=='null') CYCLE
   IF (QUANTITY2==OUTPUT_QUANTITY(ND)%NAME) THEN
      
      OUTPUT2_INDEX=ND
      
      IF (OUTPUT_QUANTITY(ND)%SPEC_ID_REQUIRED .AND. (Y_INDEX<1 .AND. Z_INDEX<0)) THEN
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

      IF (OUTPUT_QUANTITY(ND)%SPEC_ID_REQUIRED .AND. (Y_INDEX<1 .AND. Z_INDEX<0)) THEN
         IF (SPEC_ID=='null') THEN
            WRITE(MESSAGE,'(3A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),' requires a SPEC_ID'
         ELSE
            WRITE(MESSAGE,'(5A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),'. SPEC_ID ',TRIM(SPEC_ID),' not found.'
         ENDIF
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (OUTPUT_QUANTITY(ND)%PART_ID_REQUIRED .AND. PART_INDEX<1) THEN
         IF (PART_ID=='null') THEN
            WRITE(MESSAGE,'(3A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),' requires a PART_ID'
         ELSE
            WRITE(MESSAGE,'(5A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),'. PART_ID ',TRIM(PART_ID),' not found.'
         ENDIF
         CALL SHUTDOWN(MESSAGE)
      ENDIF
        
      IF (OUTPUT_QUANTITY(ND)%DUCT_ID_REQUIRED .AND. DUCT_INDEX<1) THEN
         IF (DUCT_ID=='null') THEN
            WRITE(MESSAGE,'(3A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),' requires a DUCT_ID'
         ELSE
            WRITE(MESSAGE,'(5A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),'. DUCT_ID ',TRIM(DUCT_ID),' not found.'
         ENDIF
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (OUTPUT_QUANTITY(ND)%NODE_ID_REQUIRED .AND. NODE_INDEX<1) THEN
         IF (NODE_ID=='null') THEN
            WRITE(MESSAGE,'(3A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),' requires a NODE_ID'
         ELSE
            WRITE(MESSAGE,'(5A)')  'ERROR: Output QUANTITY ',TRIM(QUANTITY),'. NODE_ID ',TRIM(NODE_ID),' not found.'
         ENDIF
         CALL SHUTDOWN(MESSAGE)
      ENDIF

      IF (( QUANTITY=='RELATIVE HUMIDITY' .OR. QUANTITY=='HUMIDITY').AND. H2O_INDEX==0) THEN
         WRITE(MESSAGE,'(A)')  'ERROR: RELATIVE HUMIDITY and HUMIDITY require SPEC=WATER VAPOR'
         CALL SHUTDOWN(MESSAGE)
      END IF

      IF (TRIM(QUANTITY)=='SURFACE DEPOSITION') THEN
         Y_INDEX = -999
         DO NS=0,N_TRACKED_SPECIES
            IF (TRIM(SPEC_ID)==TRIM(SPECIES_MIXTURE(NS)%ID)) THEN
               Z_INDEX = NS
               EXIT
            ENDIF
         ENDDO
         IF (Z_INDEX < 0) THEN
            WRITE(MESSAGE,'(A,A,A,A)')'ERROR: SURFACE DEPOSITION for ',TRIM(SPEC_ID),' is invalid as species', &
                                    ' is not a tracked species'
            CALL SHUTDOWN(MESSAGE)         
         ELSEIF (Z_INDEX==0) THEN
            WRITE(MESSAGE,'(A)')  'ERROR: Cannot select background species for deposition'
            CALL SHUTDOWN(MESSAGE)         
         ENDIF
         IF(.NOT. SPECIES_MIXTURE(Z_INDEX)%DEPOSITING) THEN
            WRITE(MESSAGE,'(A,A,A)')'ERROR: SURFACE DEPOSITION for ',TRIM(SPEC_ID),' is not an aerosol tracked species'
            CALL SHUTDOWN(MESSAGE)         
         ENDIF
         IF (SPECIES_MIXTURE(Z_INDEX)%AWM_INDEX < 0) THEN
            N_SURFACE_DENSITY_SPECIES = N_SURFACE_DENSITY_SPECIES + 1
            SPECIES_MIXTURE(Z_INDEX)%AWM_INDEX = N_SURFACE_DENSITY_SPECIES
         ENDIF
      ENDIF

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
         CASE ('BNDE')
            IF (.NOT. OUTPUT_QUANTITY(ND)%BNDE_APPROPRIATE) THEN
               WRITE(MESSAGE,'(3A)')  'ERROR: The QUANTITY ',TRIM(QUANTITY),' is not appropriate for BNDE'
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

      IF (Y_INDEX>0) THEN
         SMOKEVIEW_LABEL = TRIM(SPECIES(Y_INDEX)%ID)//' '//TRIM(QUANTITY)
         SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)//'_'//TRIM(SPECIES(Y_INDEX)%FORMULA)
      ELSEIF(Z_INDEX>=0) THEN
         SMOKEVIEW_LABEL = TRIM(SPECIES_MIXTURE(Z_INDEX)%ID)//' '//TRIM(QUANTITY)
         SMOKEVIEW_BAR_LABEL = TRIM(OUTPUT_QUANTITY(ND)%SHORT_NAME)//'_'//TRIM(SPECIES_MIXTURE(Z_INDEX)%FORMULA)
      ELSEIF (PART_INDEX>0) THEN
         SMOKEVIEW_LABEL = TRIM(LAGRANGIAN_PARTICLE_CLASS(PART_INDEX)%ID)//' '//TRIM(QUANTITY)
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



SUBROUTINE GET_SPEC_OR_SMIX_INDEX(SPEC_ID,Y_INDX,Z_INDX)

! Find the appropriate SPEC or SMIX index for the given SPEC_ID

CHARACTER(*), INTENT(IN) :: SPEC_ID
INTEGER, INTENT(OUT) :: Y_INDX,Z_INDX
INTEGER :: NS

Y_INDX = -999
Z_INDX = -999

DO NS=0,N_TRACKED_SPECIES
   IF (TRIM(SPEC_ID)==TRIM(SPECIES_MIXTURE(NS)%ID)) THEN
      Z_INDX = NS
      RETURN
   ENDIF
ENDDO

DO NS=1,N_SPECIES
   IF (TRIM(SPEC_ID)==TRIM(SPECIES(NS)%ID)) THEN
      Y_INDX = NS
      RETURN
    ENDIF
ENDDO

END SUBROUTINE GET_SPEC_OR_SMIX_INDEX



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

SUBROUTINE GET_REV_read(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') readrev(INDEX(readrev,':')+2:LEN_TRIM(readrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') readdate

END SUBROUTINE GET_REV_read

SUBROUTINE READ_CSVF
USE OUTPUT_DATA

CHARACTER(256) :: CSVFILE,UVWFILE='null'
NAMELIST /CSVF/ CSVFILE,UVWFILE

N_CSVF=0
REWIND(LU_INPUT)
COUNT_CSVF_LOOP: DO
   CALL CHECKREAD('CSVF',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_CSVF_LOOP
   READ(LU_INPUT,NML=CSVF,END=16,ERR=17,IOSTAT=IOS)
   N_CSVF=N_CSVF+1
   16 IF (IOS>0) CALL SHUTDOWN('ERROR: problem with CSVF line')
ENDDO COUNT_CSVF_LOOP
17 REWIND(LU_INPUT)

IF (N_CSVF==0) RETURN

! Allocate CSVFINFO array

ALLOCATE(CSVFINFO(N_CSVF),STAT=IZERO)
CALL ChkMemErr('READ','CSVF',IZERO)

READ_CSVF_LOOP: DO I=1,N_CSVF
   
   CALL CHECKREAD('CSVF',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_CSVF_LOOP
   
   ! Read the CSVF line
   
   READ(LU_INPUT,CSVF,END=37)
   
   CSVFINFO(I)%CSVFILE = TRIM(CSVFILE)
   IF (TRIM(UVWFILE)/='null' .AND. I>NMESHES) THEN
      CALL SHUTDOWN('Problem with CSVF line: UVWFILE must be in order with MESH.')
   ELSE
      CSVFINFO(I)%UVWFILE = UVWFILE
      UVW_RESTART = .TRUE.
   ENDIF

ENDDO READ_CSVF_LOOP
37 REWIND(LU_INPUT)

END SUBROUTINE READ_CSVF

 
END MODULE READ_INPUT

