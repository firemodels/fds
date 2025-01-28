#ifndef GITHASH_PP
#define GITHASH_PP "unknown"
#endif
#ifndef GITDATE_PP
#define GITDATE_PP "unknown"
#endif
#ifndef BUILDDATE_PP
#define BUILDDATE_PP "unknown"
#endif

!> \brief A collection of utility routines used throughout FDS.

MODULE COMP_FUNCTIONS

USE PRECISION_PARAMETERS
IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS


!> \brief Returns the wall clock time in seconds.

REAL(EB) FUNCTION CURRENT_TIME()
#ifdef WITHOUT_MPIF08
USE MPI
#else
USE MPI_F08
#endif
CURRENT_TIME = MPI_WTIME()
END FUNCTION CURRENT_TIME

!> \brief Return the current time and date in ISO_8601 format.
!> \param DATE A character string containing the date and time of day.

SUBROUTINE GET_DATE_ISO_8601(DATE)

INTEGER :: VALUES(8)
CHARACTER(LABEL_LENGTH), INTENT(OUT) :: DATE
INTEGER :: HOURS_TZ
INTEGER :: MINUTES_TZ

CALL DATE_AND_TIME(VALUES=VALUES)

HOURS_TZ = VALUES(4)/60
MINUTES_TZ = MOD(VALUES(4),60)

WRITE(DATE,'(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,".",I3.3,SP,I3.2,":",SS,I2.2)') &
        VALUES(1), &
        VALUES(2), VALUES(3), VALUES(5), VALUES(6), VALUES(7), &
        VALUES(8), HOURS_TZ, ABS(MINUTES_TZ)

END SUBROUTINE GET_DATE_ISO_8601


!> \brief Return the current time and date
!> \param DATE A character string containing the date and time of day

SUBROUTINE GET_DATE(DATE)

INTEGER :: DATE_TIME(8)
CHARACTER(10) :: BIG_BEN(3),MONTH
CHARACTER(LABEL_LENGTH), INTENT(OUT) :: DATE

CALL DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)

SELECT CASE(DATE_TIME(2))
   CASE(1)
      MONTH='January'
   CASE(2)
      MONTH='February'
   CASE(3)
      MONTH='March'
   CASE(4)
      MONTH='April'
   CASE(5)
      MONTH='May'
   CASE(6)
      MONTH='June'
   CASE(7)
      MONTH='July'
   CASE(8)
      MONTH='August'
   CASE(9)
      MONTH='September'
   CASE(10)
      MONTH='October'
   CASE(11)
      MONTH='November'
   CASE(12)
      MONTH='December'
END SELECT

WRITE(DATE,'(A,I3,A,I4,2X,I2.2,A,I2.2,A,I2.2)') TRIM(MONTH),DATE_TIME(3),', ',DATE_TIME(1), &
                                                DATE_TIME(5),':',DATE_TIME(6),':',DATE_TIME(7)
END SUBROUTINE GET_DATE


!> \brief Pause the code for DELTA_T s
!> \param DELTA_T Time to pause (s)

SUBROUTINE FDS_SLEEP(DELTA_T)

INTEGER, INTENT(IN) :: DELTA_T
REAL(EB) :: S1,S2

S1 = CURRENT_TIME()

DO
  S2 = CURRENT_TIME()
  IF (S2-S1>REAL(DELTA_T,EB)) EXIT
ENDDO

END SUBROUTINE FDS_SLEEP


!> \brief Stop the code gracefully after writing a message.
!> \param MESSAGE Character string containing an explanation for shutting down.
!> \param PROCESS_0_ONLY If .TRUE., only MPI process 0 should write the message.

SUBROUTINE SHUTDOWN(MESSAGE,PROCESS_0_ONLY)

USE GLOBAL_CONSTANTS, ONLY: MY_RANK,LU_ERR,CHID,SETUP_STOP,STOP_STATUS,N_MPI_PROCESSES
CHARACTER(*) MESSAGE
LOGICAL, INTENT(IN), OPTIONAL :: PROCESS_0_ONLY
CHARACTER(20) :: FMT

IF (STOP_STATUS /= SETUP_STOP) THEN
   IF (PRESENT(PROCESS_0_ONLY)) THEN
      WRITE(FMT,'("(/A,A,I",I0,",A,A,A)")') CEILING(LOG10(REAL(N_MPI_PROCESSES+1,EB))) + 1
      IF (.NOT.PROCESS_0_ONLY) WRITE(LU_ERR,FMT) TRIM(MESSAGE),' (MPI Process:',MY_RANK,', CHID: ',TRIM(CHID),')'
      IF (PROCESS_0_ONLY .AND. MY_RANK==0) WRITE(LU_ERR,'(/A,A,A,A)') TRIM(MESSAGE),' (CHID: ',TRIM(CHID),')'
   ELSE
      IF (MY_RANK==0) WRITE(LU_ERR,'(/A,A,A,A)') TRIM(MESSAGE),' (CHID: ',TRIM(CHID),')'
   ENDIF
ENDIF

STOP_STATUS = SETUP_STOP

END SUBROUTINE SHUTDOWN


!> \brief Return current system memory usage.
!> \param VALUE_RSS Resident stack size
!> Return the memory used by the given process at the time of the call. This only works under Linux because
!> the system file called '/proc/PID/status' is queried for the VmRSS value (kB).
!> The DEVC output QUANTITY 'RAM' outputs this value in units of MB.
!> For non-linux or non-Intel builds, this routine just returns 0 because it is non-standard Fortran and it makes use
!> of linux system files.

SUBROUTINE SYSTEM_MEM_USAGE(VALUE_RSS)

#ifdef USE_IFPORT
   USE IFPORT  ! Intel Fortran extension library. This is needed for GETPID.
#endif
INTEGER, INTENT(OUT) :: VALUE_RSS
CHARACTER(200):: FILENAME=' '
CHARACTER(80) :: LINE
CHARACTER(8)  :: PID_CHAR=' '
INTEGER :: PID
LOGICAL :: IFXST

VALUE_RSS=0  ! return zero if the memory file is not found

PID = 0
#ifdef USE_IFPORT
   PID = GETPID() ! GETPID is non-standard Fortran.
#endif

IF (PID==0) RETURN

WRITE(PID_CHAR,'(I8)') PID
FILENAME = '/proc/'//TRIM(ADJUSTL(PID_CHAR))//'/status'

INQUIRE(FILE=FILENAME,EXIST=IFXST)
IF (.NOT.IFXST) then
   WRITE (*,*) 'WARNING: Memory system file does not exist'
  RETURN
ENDIF

OPEN(1000,FILE=FILENAME,ACTION='READ')
DO
   READ (1000,'(A)',END=120) LINE
   IF (LINE(1:6)=='VmRSS:') THEN
      READ (LINE(7:),*) VALUE_RSS
      EXIT
   ENDIF
ENDDO
120 CONTINUE
CLOSE(1000)

END SUBROUTINE SYSTEM_MEM_USAGE


!> \brief Read the name of the FDS input file, which is the first argument after the fds command itself.

SUBROUTINE GET_INPUT_FILE
USE GLOBAL_CONSTANTS, ONLY: FN_INPUT
IF (FN_INPUT=='null') CALL GET_COMMAND_ARGUMENT(1,FN_INPUT)
END SUBROUTINE GET_INPUT_FILE


!> \brief Assign a unique integer to be used as a file logical unit.

INTEGER FUNCTION GET_FILE_NUMBER()
USE GLOBAL_CONSTANTS, ONLY: MY_RANK,FILE_COUNTER
FILE_COUNTER(MY_RANK) = FILE_COUNTER(MY_RANK) + 1
GET_FILE_NUMBER = FILE_COUNTER(MY_RANK)
END FUNCTION GET_FILE_NUMBER


!> \brief Look through the FDS input file for the given namelist variable and then stop at that line.
!> \param NAME The four character namelist variable
!> \param LU Logical unit of the FDS input file
!> \param IOS Error code

SUBROUTINE CHECKREAD(NAME,LU,IOS)

USE GLOBAL_CONSTANTS, ONLY: INPUT_FILE_LINE_NUMBER,STOP_STATUS,SETUP_STOP,MY_RANK,LU_ERR
INTEGER :: II
INTEGER, INTENT(OUT) :: IOS
INTEGER, INTENT(IN) :: LU
CHARACTER(4), INTENT(IN) :: NAME
CHARACTER(80) TEXT
IOS = 1

READLOOP: DO
   READ(LU,'(A)',END=10) TEXT
   INPUT_FILE_LINE_NUMBER = INPUT_FILE_LINE_NUMBER + 1
   TLOOP: DO II=1,72
      IF (TEXT(II:II)/='&' .AND. TEXT(II:II)/=' ') EXIT TLOOP
      IF (TEXT(II:II)=='&') THEN
         IF (TEXT(II+1:II+4)==NAME) THEN
            IF (TEXT(II+5:II+5)==' ') THEN
               BACKSPACE(LU)
               IOS = 0
               EXIT READLOOP
            ELSE
               IF (MY_RANK==0) WRITE(LU_ERR,'(/A,I0,A,A)') 'ERROR: Input line ',INPUT_FILE_LINE_NUMBER,&
                                                      ' is not formatted properly; NAMELIST: ',NAME
               STOP_STATUS = SETUP_STOP
               IOS = 1
               EXIT READLOOP
            ENDIF
         ELSE
            CYCLE READLOOP
         ENDIF
      ENDIF
   ENDDO TLOOP
ENDDO READLOOP

10 RETURN
END SUBROUTINE CHECKREAD


!> \brief Look for odd or illegal characters in the FDS input file.
!> \param LU Logical Unit of FDS input file.
!> \param IOS Error code, 0 means the bad text was found, 1 means it was not.
!> \param TEXT Illegal character found in FDS input file.

SUBROUTINE SCAN_INPUT_FILE(LU,IOS,TEXT)

INTEGER, INTENT(OUT) :: IOS
INTEGER, INTENT(IN) :: LU
CHARACTER(80), INTENT(OUT) :: TEXT

IOS = 1

READLOOP: DO
   READ(LU,'(A)',END=10) TEXT
   IF (IACHAR(TEXT(1:1))==13 .AND. TEXT(2:2)=='&') THEN
      IOS = 0
      EXIT READLOOP
   ENDIF
ENDDO READLOOP

10 RETURN
END SUBROUTINE SCAN_INPUT_FILE


!> \brief Look for a certain TEXT string in the input file.
!> \param LU Logical Unit of the FDS input file
!> \param TEXT Character string to search for
!> \param FOUND T if the TEXT is found; F if it is not

SUBROUTINE SEARCH_INPUT_FILE(LU,TEXT,FOUND)

! Look for TEXT in the input file.

INTEGER, INTENT(IN) :: LU
LOGICAL, INTENT(OUT) :: FOUND
CHARACTER(*), INTENT(IN) :: TEXT
CHARACTER(NAMELIST_LENGTH) :: LINE
INTEGER :: IND

FOUND = .FALSE.
REWIND(LU)

READLOOP: DO
   READ(LU,'(A)',END=10) LINE
   IND = INDEX(LINE,TEXT)
   IF (IND>0) THEN
      FOUND = .TRUE.
      EXIT READLOOP
   ENDIF
ENDDO READLOOP

10 RETURN
END SUBROUTINE SEARCH_INPUT_FILE


!> \brief Read through a comma-delimited output file and stop when the time, T, exceeds the last time read.
!> \param LU Logical Unit of the file to be appended
!> \param N_TEXT_LINES Number of header lines in the file
!> \param T The time at which to append the data

SUBROUTINE APPEND_FILE(LU,N_TEXT_LINES,T)

USE GLOBAL_CONSTANTS, ONLY: CLIP_RESTART_FILES
INTEGER, INTENT(IN) :: LU,N_TEXT_LINES
INTEGER :: ITER
REAL(EB) :: FILE_TIME
REAL(EB), INTENT(IN) :: T

READ_TEXT_LINES: DO ITER=1,N_TEXT_LINES
   READ(LU,*,END=10)
ENDDO READ_TEXT_LINES

READ_LOOP: DO
   READ(LU,*,END=10) FILE_TIME
   IF (FILE_TIME>T .AND. CLIP_RESTART_FILES) THEN
      BACKSPACE(LU)
      EXIT READ_LOOP
   ENDIF
ENDDO READ_LOOP

RETURN
10 BACKSPACE(LU)

END SUBROUTINE APPEND_FILE


!> \brief Reorder an input sextuple XB if needed.
!> \param XB User-specified real sextuplet.

SUBROUTINE CHECK_XB(XB)
REAL(EB) :: DUMMY,XB(6)
INTEGER  :: I
DO I=1,5,2
   IF (XB(I)>XB(I+1)) THEN
      DUMMY   = XB(I)
      XB(I)   = XB(I+1)
      XB(I+1) = DUMMY
   ENDIF
ENDDO
END SUBROUTINE CHECK_XB


!> \brief Change the units of the output quantity if it is an integrated quantity.
!> \param QUANTITY Name of output quantity.
!> \param UNITS Quantity units to be changed.
!> \param SPATIAL_STATISTIC Type of spatial integration.
!> \param TEMPORAL_STATISTIC Type of time integration.
!> \param LU_ERR Logical Unit of error output file.

SUBROUTINE CHANGE_UNITS(QUANTITY,UNITS,SPATIAL_STATISTIC,TEMPORAL_STATISTIC,LU_ERR)

USE GLOBAL_CONSTANTS, ONLY: MY_RANK
CHARACTER(LABEL_LENGTH), INTENT(IN) :: QUANTITY,SPATIAL_STATISTIC,TEMPORAL_STATISTIC
INTEGER, INTENT(IN) :: LU_ERR
CHARACTER(LABEL_LENGTH), INTENT(INOUT) :: UNITS
CHARACTER(LABEL_LENGTH) :: NEW_UNITS
CHARACTER(MESSAGE_LENGTH) :: MESSAGE
INTEGER :: I,UNIT_L

UNIT_L = LEN(TRIM(UNITS))
NEW_UNITS = UNITS
I = 1
SELECT CASE (SPATIAL_STATISTIC)
   CASE('VOLUME INTEGRAL')
      I = INDEX(UNITS,'/m3')
      IF (I/=0) WRITE(NEW_UNITS,'(A,A)') UNITS(1:I-1),UNITS(I+3:UNIT_L)
      IF (TRIM(UNITS)=='1/s') THEN
         NEW_UNITS = 'm3/s'
         I = 1
      ENDIF
      IF (TRIM(UNITS)=='Pa') THEN
         NEW_UNITS = 'J'
         I = 1
      ENDIF
   CASE('MASS INTEGRAL')
      I = INDEX(UNITS,'/kg')
      IF (I/=0) WRITE(NEW_UNITS,'(A,A)') UNITS(1:I-1),UNITS(I+3:UNIT_L)
   CASE('AREA INTEGRAL','SURFACE INTEGRAL')
      I = INDEX(UNITS,'/m2')
      IF (I/=0) WRITE(NEW_UNITS,'(A,A)') UNITS(1:I-1),UNITS(I+3:UNIT_L)
      IF (I==0 .AND. UNITS=='m/s') THEN
         NEW_UNITS = 'm3/s'
         I = 1
      ENDIF
   CASE ('VOLUME')
      NEW_UNITS = 'm3'
      I=1
   CASE ('AREA','SURFACE AREA')
      NEW_UNITS = 'm2'
      I=1
   CASE ('MINLOC X','MINLOC Y','MINLOC Z','MAXLOC X','MAXLOC Y','MAXLOC Z','CENTROID X','CENTROID Y','CENTROID Z')
      NEW_UNITS = 'm'
      I=1
END SELECT

UNITS = NEW_UNITS

IF (I==0) THEN
   WRITE(MESSAGE,'(A,A,A,A)')  'WARNING: Problem with units compatibility of SPATIAL_STATISITIC ',TRIM(SPATIAL_STATISTIC), &
                               ' with the QUANTITY ',TRIM(QUANTITY)
   IF (MY_RANK==0) WRITE(LU_ERR,'(A)') TRIM(MESSAGE)
ENDIF

I = 1

SELECT CASE (TEMPORAL_STATISTIC)
   CASE('TIME INTEGRAL')
      I = INDEX(UNITS,'/s')
      IF (I/=0) THEN
         WRITE(NEW_UNITS,'(A,A)') UNITS(1:I-1),UNITS(I+2:UNIT_L)
      ELSE
         I = INDEX(UNITS,'kW')
         IF (I==1) WRITE(NEW_UNITS,'(A,A)') 'kJ',UNITS(3:UNIT_L)
         IF (I>1)  WRITE(NEW_UNITS,'(A,A,A)') UNITS(1:I-1),'kJ',UNITS(I+2:UNIT_L)
      ENDIF
      IF (QUANTITY=='FED') I=1
   CASE('MAX TIME','MIN TIME')
      NEW_UNITS = 's'
END SELECT

UNITS = NEW_UNITS

IF (I==0) THEN
   WRITE(MESSAGE,'(A,A,A,A)')  'WARNING: Problem with units compatibility of TEMPORAL_STATISTIC ',TRIM(TEMPORAL_STATISTIC), &
                               ' with the QUANTITY ',TRIM(QUANTITY)
   IF (MY_RANK==0) WRITE(LU_ERR,'(A)') TRIM(MESSAGE)
ENDIF

END SUBROUTINE CHANGE_UNITS


!> \brief Process the DUMP parameters associated with output files
!> \param T Current time (s)

SUBROUTINE INITIALIZE_OUTPUT_CLOCKS(T)

USE MESH_VARIABLES, ONLY: MESHES
USE TYPES, ONLY: RAMPS, RAMPS_TYPE
USE GLOBAL_CONSTANTS, ONLY: T_BEGIN,T_END,NFRAMES,TIME_SHRINK_FACTOR
USE OUTPUT_CLOCKS
USE GLOBAL_CONSTANTS, ONLY: LOWER_MESH_INDEX,UPPER_MESH_INDEX
REAL(EB), INTENT(IN) :: T
REAL(EB) :: DT_DEFAULT
INTEGER :: N,NM,IS
TYPE(RAMPS_TYPE), POINTER :: RP
LOGICAL :: TWO_D_SLICES,THREE_D_SLICES

! Set output time intervals

DT_DEFAULT = (T_END-T_BEGIN)/REAL(NFRAMES,EB)

CALL SET_OUTPUT_CLOCK(DT_BNDF   ,RAMP_BNDF_INDEX,BNDF_CLOCK,BNDF_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_CPU    ,RAMP_CPU_INDEX , CPU_CLOCK, CPU_COUNTER,.FALSE.,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_CTRL   ,RAMP_CTRL_INDEX,CTRL_CLOCK,CTRL_COUNTER,.FALSE.,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_DEVC   ,RAMP_DEVC_INDEX,DEVC_CLOCK,DEVC_COUNTER,.FALSE.,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_FLUSH  ,RAMP_FLSH_INDEX,FLSH_CLOCK,FLSH_COUNTER,.FALSE.,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_GEOM   ,RAMP_GEOM_INDEX,GEOM_CLOCK,GEOM_COUNTER,.FALSE.,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_HRR    ,RAMP_HRR_INDEX , HRR_CLOCK, HRR_COUNTER,.FALSE.,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_HVAC   ,RAMP_HVAC_INDEX,HVAC_CLOCK,HVAC_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_ISOF   ,RAMP_ISOF_INDEX,ISOF_CLOCK,ISOF_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_MASS   ,RAMP_MASS_INDEX,MASS_CLOCK,MASS_COUNTER,.FALSE.,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_PART   ,RAMP_PART_INDEX,PART_CLOCK,PART_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_PL3D   ,RAMP_PL3D_INDEX,PL3D_CLOCK,PL3D_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_PROF   ,RAMP_PROF_INDEX,PROF_CLOCK,PROF_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_RADF   ,RAMP_RADF_INDEX,RADF_CLOCK,RADF_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_RESTART,RAMP_RSRT_INDEX,RSRT_CLOCK,RSRT_COUNTER,.FALSE.,.FALSE.)
CALL SET_OUTPUT_CLOCK(DT_SLCF   ,RAMP_SLCF_INDEX,SLCF_CLOCK,SLCF_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_SL3D   ,RAMP_SL3D_INDEX,SL3D_CLOCK,SL3D_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_SMOKE3D,RAMP_SM3D_INDEX,SM3D_CLOCK,SM3D_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_UVW    ,RAMP_UVW_INDEX , UVW_CLOCK, UVW_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_TMP    ,RAMP_TMP_INDEX , TMP_CLOCK, TMP_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_SPEC   ,RAMP_SPEC_INDEX,SPEC_CLOCK,SPEC_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_SLCF_VTK    ,RAMP_SLCF_INDEX,  SLCF_VTK_CLOCK, SLCF_VTK_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_SL3D_VTK    ,RAMP_SL3D_INDEX,  SL3D_VTK_CLOCK, SL3D_VTK_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_SM3D_VTK    ,RAMP_SM3D_INDEX,  SM3D_VTK_CLOCK, SM3D_VTK_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_BNDF_VTK    ,RAMP_BNDF_INDEX,  BNDF_VTK_CLOCK, BNDF_VTK_COUNTER,.TRUE. ,.TRUE. )
CALL SET_OUTPUT_CLOCK(DT_PART_VTK    ,RAMP_PART_INDEX,  PART_VTK_CLOCK, PART_VTK_COUNTER,.TRUE. ,.TRUE. )

! Special cases

GEOM_CLOCK(0) = T_BEGIN

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (PL3D_CLOCK(0)==T_BEGIN .AND. PL3D_COUNTER(NM)==0) PL3D_COUNTER(NM) = 1
   THREE_D_SLICES = .FALSE.
   TWO_D_SLICES   = .FALSE.
   DO IS=1,MESHES(NM)%N_SLCF
      IF (MESHES(NM)%SLICE(IS)%THREE_D) THEN
         THREE_D_SLICES = .TRUE.
      ELSE
         TWO_D_SLICES = .TRUE.
      ENDIF
   ENDDO
   IF (.NOT.THREE_D_SLICES) SL3D_COUNTER(NM) = UBOUND(SL3D_CLOCK,DIM=1)
   IF (.NOT.TWO_D_SLICES)   SLCF_COUNTER(NM) = UBOUND(SLCF_CLOCK,DIM=1)
ENDDO

CONTAINS

SUBROUTINE SET_OUTPUT_CLOCK(DT_OUT,RAMP_INDEX,CLOCK,COUNTER,MESH_SPECIFIC_COUNTER,TRUNCATE_LAST_INTERVAL)

INTEGER, INTENT(IN) :: RAMP_INDEX
REAL(EB), ALLOCATABLE, DIMENSION(:) :: CLOCK
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTER
REAL(EB) :: DT_OUT
INTEGER :: NF
CHARACTER(MESSAGE_LENGTH) :: MESSAGE
LOGICAL, INTENT(IN) :: MESH_SPECIFIC_COUNTER,TRUNCATE_LAST_INTERVAL

IF (RAMP_INDEX>0) THEN
   RP => RAMPS(RAMP_INDEX)
   IF (RP%EXTERNAL_FILE) THEN
      WRITE(MESSAGE,'(2A)') 'ERROR(ZZZ): Cannot set EXTERNAL_FILE=T for a RAMP used for an output clock. RAMP: ',TRIM(RP%ID)
      CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.); RETURN
   ENDIF
   NF = RP%NUMBER_DATA_POINTS
   ALLOCATE(CLOCK(0:NF))
   CLOCK(0:NF-1) = RP%INDEPENDENT_DATA(1:NF)
   CLOCK(NF)     = HUGE(EB)
ELSE
   IF (DT_OUT<0._EB) THEN
      DT_OUT = DT_DEFAULT
   ELSEIF (DT_OUT==0._EB) THEN
      DT_OUT = 1._EB
   ELSE
      DT_OUT = MAX(1.E-6_EB,DT_OUT/TIME_SHRINK_FACTOR)
   ENDIF
   NF = CEILING((T_END-T_BEGIN)/DT_OUT)
   ALLOCATE(CLOCK(0:NF))
   DO N=0,NF
      CLOCK(N) = T_BEGIN + N*DT_OUT
   ENDDO
   IF (TRUNCATE_LAST_INTERVAL) CLOCK(NF) = T_END
   IF (DT_OUT>(T_END-T_BEGIN)) CLOCK = HUGE(EB)
ENDIF

DO N=0,NF
   IF (T==T_BEGIN .OR. T<CLOCK(N)) THEN
      IF (MESH_SPECIFIC_COUNTER) THEN
         DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
            IF (NM==LOWER_MESH_INDEX) ALLOCATE(COUNTER(LOWER_MESH_INDEX:UPPER_MESH_INDEX))
            COUNTER(NM) = N
         ENDDO
      ELSE
         ALLOCATE(COUNTER(1))
         COUNTER(1) = N
      ENDIF
      EXIT
   ENDIF
ENDDO

END SUBROUTINE SET_OUTPUT_CLOCK

END SUBROUTINE INITIALIZE_OUTPUT_CLOCKS

!> \brief Read in values for RAMP or CTRL controlled by an external file.

SUBROUTINE READ_EXTERNAL_FILE(FAILED)
USE GLOBAL_CONSTANTS, ONLY: N_RAMP,EXTERNAL_RAMP,EXTERNAL_CTRL,EXTERNAL_FILENAME,EXTERNAL_HEARTBEAT_FILENAME,&
                            LU_EXTERNAL,LU_EXTERNAL_HEARTBEAT,HEARTBEAT_STOP,HEARTBEAT_FAIL,DT_EXTERNAL_HEARTBEAT,STOP_STATUS
USE TYPES, ONLY: RAMPS
USE CONTROL_VARIABLES, ONLY: N_CTRL,CONTROL
INTEGER:: NC, NR
CHARACTER(LABEL_LENGTH) :: INPUT_LABEL,INPUT_TYPE
CHARACTER(250) :: TEST_FILENAME
REAL(EB) :: INPUT_REAL
INTEGER :: T_HEARTBEAT,DT_PAUSE
LOGICAL :: HEARTBEAT_FOUND
LOGICAL,INTENT(OUT) :: FAILED

FAILED = .FALSE.
EXTERNAL_RAMP = RAMPS%LAST
EXTERNAL_CTRL = CONTROL%CURRENT_STATE

IF (DT_EXTERNAL_HEARTBEAT > 0) THEN
   T_HEARTBEAT = 0
   HEARTBEAT_LOOP: DO
      INQUIRE(FILE=EXTERNAL_HEARTBEAT_FILENAME,EXIST=HEARTBEAT_FOUND)
      IF (HEARTBEAT_FOUND) THEN
         EXIT HEARTBEAT_LOOP
      ELSE
         DT_PAUSE = MIN(1,DT_EXTERNAL_HEARTBEAT - T_HEARTBEAT)
         IF (DT_PAUSE <= 0._EB) THEN
            IF (HEARTBEAT_FAIL) STOP_STATUS = HEARTBEAT_STOP
            FAILED = .TRUE.
            RETURN
         ELSE
            CALL FDS_SLEEP(DT_PAUSE)
            T_HEARTBEAT = T_HEARTBEAT + DT_PAUSE
         ENDIF
      ENDIF
   ENDDO HEARTBEAT_LOOP
   OPEN(LU_EXTERNAL_HEARTBEAT,FILE=EXTERNAL_HEARTBEAT_FILENAME,STATUS='OLD',ERR=100)
   READ(UNIT=LU_EXTERNAL_HEARTBEAT,END=300,ERR=300,FMT=*) TEST_FILENAME
   EXTERNAL_FILENAME = TEST_FILENAME
300 CONTINUE
   CLOSE(LU_EXTERNAL_HEARTBEAT,STATUS='DELETE')
ENDIF

OPEN(LU_EXTERNAL,FILE=EXTERNAL_FILENAME,ERR=100,STATUS='OLD')

DO
   READ(UNIT=LU_EXTERNAL,END=200,ERR=110,FMT=*) INPUT_TYPE,INPUT_LABEL,INPUT_REAL
   IF (INPUT_TYPE=='RAMP') THEN
      DO NR=1,N_RAMP
         IF (TRIM(RAMPS(NR)%ID) == TRIM(INPUT_LABEL)) THEN
            EXTERNAL_RAMP(NR) = INPUT_REAL
            EXIT
         ENDIF
      ENDDO
   ELSEIF (INPUT_TYPE=='CTRL')  THEN
      DO NC=1,N_CTRL
         IF (TRIM(CONTROL(NC)%ID) == TRIM(INPUT_LABEL)) THEN
            IF (INPUT_REAL < 0) THEN
               EXTERNAL_CTRL(NC) = .FALSE.
            ELSE
               EXTERNAL_CTRL(NC) = .TRUE.
            ENDIF
            EXIT
         ENDIF
      ENDDO
   ENDIF

ENDDO

RETURN

110 CLOSE(LU_EXTERNAL)
100 FAILED = .TRUE.
RETURN

200 CLOSE(LU_EXTERNAL)
RETURN

END SUBROUTINE READ_EXTERNAL_FILE


END MODULE COMP_FUNCTIONS


!> \brief This module contains a routine that allows one to switch the left and right sides of an equality

MODULE COMP_OPERATORS

USE PRECISION_PARAMETERS
IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE
PUBLIC EQUATE

INTERFACE EQUATE
   MODULE PROCEDURE EQUATE_REALS
   MODULE PROCEDURE EQUATE_REAL_VECTORS
   MODULE PROCEDURE EQUATE_INTEGERS
   MODULE PROCEDURE EQUATE_INTEGER_VECTORS
   MODULE PROCEDURE EQUATE_LOGICALS
   MODULE PROCEDURE EQUATE_LOGICAL_VECTORS
END INTERFACE

CONTAINS

SUBROUTINE EQUATE_REALS(A,B,SWAP)
REAL(EB), INTENT(INOUT) :: A,B
LOGICAL, INTENT(IN) :: SWAP
IF (SWAP) THEN
   B = A
ELSE
   A = B
ENDIF
END SUBROUTINE EQUATE_REALS

SUBROUTINE EQUATE_REAL_VECTORS(A,B,SWAP)
REAL(EB), INTENT(INOUT), DIMENSION(:) :: A,B
LOGICAL, INTENT(IN) :: SWAP
IF (SIZE(A)==0) RETURN
IF (SWAP) THEN
   B = A
ELSE
   A = B
ENDIF
END SUBROUTINE EQUATE_REAL_VECTORS

SUBROUTINE EQUATE_INTEGERS(A,B,SWAP)
INTEGER, INTENT(INOUT) :: A,B
LOGICAL, INTENT(IN) :: SWAP
IF (SWAP) THEN
   B = A
ELSE
   A = B
ENDIF
END SUBROUTINE EQUATE_INTEGERS

SUBROUTINE EQUATE_INTEGER_VECTORS(A,B,SWAP)
INTEGER, INTENT(INOUT), DIMENSION(:) :: A,B
LOGICAL, INTENT(IN) :: SWAP
IF (SIZE(A)==0) RETURN
IF (SWAP) THEN
   B = A
ELSE
   A = B
ENDIF
END SUBROUTINE EQUATE_INTEGER_VECTORS

SUBROUTINE EQUATE_LOGICALS(A,B,SWAP)
LOGICAL, INTENT(INOUT) :: A,B
LOGICAL, INTENT(IN) :: SWAP
IF (SWAP) THEN
   B = A
ELSE
   A = B
ENDIF
END SUBROUTINE EQUATE_LOGICALS

SUBROUTINE EQUATE_LOGICAL_VECTORS(A,B,SWAP)
LOGICAL, INTENT(INOUT), DIMENSION(:) :: A,B
LOGICAL, INTENT(IN) :: SWAP
IF (SWAP) THEN
   B = A
ELSE
   A = B
ENDIF
END SUBROUTINE EQUATE_LOGICAL_VECTORS

END MODULE COMP_OPERATORS


!> \brief Routines that do various mathematical manipulations

MODULE MATH_FUNCTIONS

USE PRECISION_PARAMETERS
IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS


!> \brief Add new value to a growing histogram
!> \param NBINS Number of bins in the histogram
!> \param LIMITS The range of values underlying the histogram
!> \param COUNTS Array of weight currently assigned to each bin
!> \param VAL New value to be added
!> \param WEIGHT Relative weight associated with the new value

SUBROUTINE UPDATE_HISTOGRAM(NBINS,LIMITS,COUNTS,VAL,WEIGHT)
INTEGER,INTENT(IN)::NBINS
REAL(EB), INTENT(IN)::LIMITS(2),VAL,WEIGHT
REAL(EB), INTENT(INOUT) :: COUNTS(NBINS)
INTEGER::IND=0
IND=MIN(NBINS,MAX(CEILING((VAL-LIMITS(1))/(LIMITS(2)-LIMITS(1))*NBINS),1))
COUNTS(IND)=COUNTS(IND)+WEIGHT
END SUBROUTINE UPDATE_HISTOGRAM

!> \brief Calculate the value of polynomial function.
!> \param N Number of coefficients in the polynomial
!> \param TEMP The independent variable
!> \param COEF The array of coefficients

REAL(EB) FUNCTION POLYVAL(N,TEMP,COEF)
INTEGER, INTENT(IN) :: N
REAL(EB), INTENT(IN) :: TEMP,COEF(N)
INTEGER :: I
POLYVAL = 0._EB
DO I=1,N
   POLYVAL  = POLYVAL  + COEF(I)*TEMP**(I-1)
ENDDO
END FUNCTION POLYVAL


!> \brief Determine the index of a table with a given name
!> \param ID Name of the table
!> \param TYPE Kind of table
!> \param TABLE_INDEX Index of the table

SUBROUTINE GET_TABLE_INDEX(ID,TYPE,TABLE_INDEX)

USE GLOBAL_CONSTANTS, ONLY: N_TABLE,TABLE_ID,TABLE_TYPE
CHARACTER(*), INTENT(IN) :: ID
INTEGER, INTENT(IN) :: TYPE
INTEGER, INTENT(OUT) :: TABLE_INDEX
INTEGER :: NT

IF (ID=='null') THEN
   TABLE_INDEX = 0
   RETURN
ENDIF

SEARCH: DO NT=1,N_TABLE
   IF (ID==TABLE_ID(NT)) THEN
      TABLE_INDEX = NT
      RETURN
   ENDIF
ENDDO SEARCH

N_TABLE                = N_TABLE + 1
TABLE_INDEX            = N_TABLE
TABLE_ID(TABLE_INDEX)   = ID
TABLE_TYPE(TABLE_INDEX) = TYPE

END SUBROUTINE GET_TABLE_INDEX


!> \brief Return an interpolated value from a given RAMP function
!> \param RAMP_INPUT Independent variable of the RAMP function
!> \param RAMP_INDEX Index of the RAMP function
!> \param TAU Time scale used for reserved t-squared or tanh ramps

REAL(EB) FUNCTION EVALUATE_RAMP(RAMP_INPUT,RAMP_INDEX,TAU)

USE GLOBAL_CONSTANTS, ONLY : EXTERNAL_RAMP
USE TYPES, ONLY: RAMPS
USE DEVICE_VARIABLES, ONLY: DEVICE
USE CONTROL_VARIABLES, ONLY: CONTROL
REAL(EB), INTENT(IN) :: RAMP_INPUT
REAL(EB), INTENT(IN), OPTIONAL :: TAU
REAL(EB):: RAMP_POSITION
INTEGER:: I
INTEGER, INTENT(IN)  :: RAMP_INDEX

SELECT CASE(RAMP_INDEX)
   CASE(-2)
      EVALUATE_RAMP = MAX( TANH(RAMP_INPUT/TAU), 0._EB )
   CASE(-1)
      EVALUATE_RAMP = MIN( (RAMP_INPUT/TAU)**2 , 1.0_EB )
   CASE( 0)
      EVALUATE_RAMP = 1._EB
   CASE(1:)
      IF (RAMPS(RAMP_INDEX)%DEVC_INDEX > 0) THEN
         RAMP_POSITION = &
            MAX(0._EB,MIN(RAMPS(RAMP_INDEX)%SPAN,DEVICE(RAMPS(RAMP_INDEX)%DEVC_INDEX)%SMOOTHED_VALUE - RAMPS(RAMP_INDEX)%T_MIN))
      ELSEIF (RAMPS(RAMP_INDEX)%CTRL_INDEX > 0) THEN
         RAMP_POSITION = &
            MAX(0._EB,MIN(RAMPS(RAMP_INDEX)%SPAN,CONTROL(RAMPS(RAMP_INDEX)%CTRL_INDEX)%INSTANT_VALUE - RAMPS(RAMP_INDEX)%T_MIN))
      ELSEIF (RAMPS(RAMP_INDEX)%DEVC_DEP_INDEX > 0) THEN
         EVALUATE_RAMP = DEVICE(RAMPS(RAMP_INDEX)%DEVC_DEP_INDEX)%SMOOTHED_VALUE
         RETURN
      ELSEIF (RAMPS(RAMP_INDEX)%CTRL_DEP_INDEX > 0) THEN
         EVALUATE_RAMP = CONTROL(RAMPS(RAMP_INDEX)%CTRL_DEP_INDEX)%INSTANT_VALUE
         RAMPS(RAMP_INDEX)%LAST = EVALUATE_RAMP
         RETURN
      ELSEIF (RAMPS(RAMP_INDEX)%EXTERNAL_FILE) THEN
         EVALUATE_RAMP = EXTERNAL_RAMP(RAMP_INDEX)
         RETURN
      ELSE
         RAMP_POSITION = &
            MAX(0._EB,MIN(RAMPS(RAMP_INDEX)%SPAN,RAMP_INPUT - RAMPS(RAMP_INDEX)%T_MIN))
      ENDIF
      I = MIN(UBOUND(RAMPS(RAMP_INDEX)%INTERPOLATED_DATA,1),&
          MAX(LBOUND(RAMPS(RAMP_INDEX)%INTERPOLATED_DATA,1),NINT(RAMP_POSITION*RAMPS(RAMP_INDEX)%RDT)))
      EVALUATE_RAMP = RAMPS(RAMP_INDEX)%INTERPOLATED_DATA(I)
END SELECT

END FUNCTION EVALUATE_RAMP


!> \brief Compute inverse erfc function
!> \param Y Y=ERFC(X), where X is returned by IERFC(Y)

REAL(EB) FUNCTION IERFC(Y)

REAL(EB), INTENT(IN) :: Y
REAL(EB) :: QA,QB,QC,QD,Q0,Q1,Q2,Q3,Q4,PA,PB,P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10, &
            P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,X,S,T,U,W,Z
PARAMETER ( &
QA = 9.16461398268964-01_EB, &
QB = 2.31729200323405-01_EB, &
QC = 4.88826640273108-01_EB, &
QD = 1.24610454613712-01_EB, &
Q0 = 4.99999303439796-01_EB, &
Q1 = 1.16065025341614-01_EB, &
Q2 = 1.50689047360223-01_EB, &
Q3 = 2.69999308670029-01_EB, &
Q4 = -7.28846765585675-02_EB)
PARAMETER ( &
PA = 3.97886080735226000+00_EB, &
PB = 1.20782237635245222-01_EB, &
P0 = 2.44044510593190935-01_EB, &
P1 = 4.34397492331430115-01_EB, &
P2 = 6.86265948274097816-01_EB, &
P3 = 9.56464974744799006-01_EB, &
P4 = 1.16374581931560831+00_EB, &
P5 = 1.21448730779995237+00_EB, &
P6 = 1.05375024970847138+00_EB, &
P7 = 7.13657635868730364-01_EB, &
P8 = 3.16847638520135944-01_EB, &
P9 = 1.47297938331485121-02_EB, &
P10 = -1.05872177941595488-01_EB, &
P11 = -7.43424357241784861-02_EB)
PARAMETER ( &
P12 = 2.20995927012179067-03_EB, &
P13 = 3.46494207789099922-02_EB, &
P14 = 1.42961988697898018-02_EB, &
P15 = -1.18598117047771104-02_EB, &
P16 = -1.12749169332504870-02_EB, &
P17 = 3.39721910367775861-03_EB, &
P18 = 6.85649426074558612-03_EB, &
P19 = -7.71708358954120939-04_EB, &
P20 = -3.51287146129100025-03_EB, &
P21 = 1.05739299623423047-04_EB, &
P22 = 1.12648096188977922-03_EB)

Z = Y
IF (Y  > 1._EB) Z = 2._EB - Y
W = QA - LOG(Z)
U = SQRT(W)
S = (QC + LOG(U)) / W
T = 1._EB / (U + QB)

X = U * (1._EB - S * (0.5_EB + S * QD)) - ((((Q4 * T + Q3) * T + Q2) * T + Q1) * T + Q0) * T
T = PA / (PA + X)
U = T - 0.5_EB

S = (((((((((P22 * U + P21) * U + P20) * U + P19) * U + P18) * U + P17) * U + P16) * U + P15) * U + P14) * U + P13) * U + P12

S = ((((((((((((S * U + P11) * U + P10) * U +  P9) * U + P8) * U + P7) * U + P6) * U + P5) * U + P4) * U + P3) &
    * U + P2) * U + P1) * U + P0) * T - Z * EXP(X * X - PB)

X = X + S * (1._EB + X * S)

IF (Y > 1._EB) X = -X

IERFC = X

END FUNCTION IERFC


!> \brief Solve a linear system of equations with Gauss-Jordon elimination (Press, Numerical Recipes)
!> \param A Primary matrix of A*x=b
!> \param N Dimension of A
!> \param NP Dimension of array containing A
!> \param B Array of vectors B of the linear system
!> \param M Number of columns of B
!> \param MP Number of columns of the array holding the B vectors
!> \param IERROR Error code

SUBROUTINE GAUSSJ(A,N,NP,B,M,MP,IERROR)

INTEGER, INTENT(IN) :: M,MP,N,NP
REAL(EB), INTENT(INOUT) :: A(NP,NP),B(NP,MP)
INTEGER, INTENT(OUT) :: IERROR
REAL(EB) :: BIG,DUM,PIVINV
INTEGER :: I,ICOL=0,IROW=0,J,K,L,LL,INDXC(NP),INDXR(NP),IPIV(NP)

IERROR = 0
IPIV(1:N) = 0

DO I=1,N
   BIG = 0._EB
   DO J=1,N
      IF (IPIV(J)/=1) THEN
         DO K=1,N
            IF (IPIV(K)==0) THEN
               IF (ABS(A(J,K))>=BIG) THEN
                  BIG = ABS(A(J,K))
                  IROW = J
                  ICOL = K
               ENDIF
            ELSE IF (IPIV(K)>1) THEN
               IERROR = 103   ! Singular matrix in gaussj
               RETURN
            ENDIF
         ENDDO
      ENDIF
   ENDDO
   IPIV(ICOL) = IPIV(ICOL) + 1
   IF (IROW/=ICOL) THEN
      DO L=1,N
         DUM = A(IROW,L)
         A(IROW,L) = A(ICOL,L)
         A(ICOL,L) = DUM
      ENDDO
      DO L=1,M
         DUM = B(IROW,L)
         B(IROW,L) = B(ICOL,L)
         B(ICOL,L) = DUM
      ENDDO
   ENDIF
   INDXR(I) = IROW
   INDXC(I) = ICOL
   IF (ABS(A(ICOL,ICOL))<=TWO_EPSILON_EB) THEN
      IERROR = 103  ! Singular matrix in gaussj
      RETURN
      ENDIF
   PIVINV = 1._EB/A(ICOL,ICOL)
   A(ICOL,ICOL) = 1._EB
   A(ICOL,1:N) = A(ICOL,1:N) * PIVINV
   B(ICOL,1:M) = B(ICOL,1:M) * PIVINV
   DO LL=1,N
      IF (LL/=ICOL) THEN
         DUM = A(LL,ICOL)
         A(LL,ICOL) = 0._EB
         A(LL,1:N) = A(LL,1:N) - A(ICOL,1:N)*DUM
         B(LL,1:M) = B(LL,1:M) - B(ICOL,1:M)*DUM
      ENDIF
   ENDDO
ENDDO
DO L=N,1,-1
   IF (INDXR(L)/=INDXC(L)) THEN
      DO K=1,N
         DUM = A(K,INDXR(L))
         A(K,INDXR(L)) = A(K,INDXC(L))
         A(K,INDXC(L)) = DUM
      ENDDO
   ENDIF
ENDDO

END SUBROUTINE GAUSSJ


!> \brief Solve a linear system of equations for m=n and m/=n
!> \param A Primary matrix of A*x=b
!> \param X Solution vector of A*x=b
!> \param M Column dimension of A and dimension of x
!> \param N Row dimension of A and dimension of b
!> \param B Constant vector of b A*x=b
!> \param IERR Error code

SUBROUTINE LINEAR_SYSTEM_SOLVE(M,N,A,B,X,IERR)
INTEGER, INTENT(IN) :: M,N
REAL(EB), INTENT(INOUT) :: A(N,M),B(N),X(M)
REAL(EB) :: AT(M,N),AAT(N,N),ATA(M,M),ATB(M)
INTEGER, INTENT(OUT) :: IERR

IERR = 0
! System is underdetermined - find a minimal solution
! Solution is given by x = A^T t, solve t = (A A^T)**-1 b, get x as A^T t.
IF (M > N) THEN
   AT = TRANSPOSE(A)
   AAT = MATMUL(A,AT)
   CALL GAUSSJ(AAT,N,N,B,1,1,IERR)
   IF (IERR > 0) THEN
      X = 0._EB
   ELSE
      X = MATMUL(AT,B)
   ENDIF
! System is overdetermined - find least squares solution
! Solution is x = (A^T A)**-1 A^T b
ELSEIF (N > M) THEN
   AT = TRANSPOSE(A)
   ATA = MATMUL(AT,A)
   ATB = MATMUL(AT,B)
   CALL GAUSSJ(ATA,M,M,ATB,1,1,IERR)
   IF (IERR > 0) THEN
      X = 0._EB
   ELSE
      IERR = 200
      X = ATB
   ENDIF
! Solution is x = A**-1 b
ELSE
   CALL GAUSSJ(A,M,M,B,1,1,IERR)
   IF (IERR > 0) THEN
      X = 0._EB
   ELSE
      X = B
   ENDIF
ENDIF

END SUBROUTINE LINEAR_SYSTEM_SOLVE


!> \brief Linearly interpolate the value of a given function at a given point
!> \param X Independent variable
!> \param Y Dependent variable
!> \param XI Point to interpolate
!> \param ANS Value of the function at the point XI

SUBROUTINE INTERPOLATE1D(X,Y,XI,ANS)

REAL(EB), INTENT(IN), DIMENSION(:) :: X, Y
REAL(EB), INTENT(IN) :: XI
REAL(EB), INTENT(OUT) :: ANS
INTEGER I, UX,LX

UX = UBOUND(X,1)
LX = LBOUND(X,1)

IF (XI <= X(LX)) THEN
   ANS = Y(LX)
ELSEIF (XI >= X(UX)) THEN
   ANS = Y(UX)
ELSE
   L1: DO I=LX,UX-1
      IF (ABS(XI -X(I)) <= SPACING(X(I))) THEN
         ANS = Y(I)
         EXIT L1
      ELSEIF (X(I+1)>XI) THEN
         ANS = Y(I)+(XI-X(I))/(X(I+1)-X(I)) * (Y(I+1)-Y(I))
         EXIT L1
      ENDIF
   ENDDO L1
ENDIF

END SUBROUTINE INTERPOLATE1D


!> \brief Interpolate a 1D array of numbers
!> \param LOWER Lower index of the array X
!> \param X Array of numbers
!> \param XI Real number representing a fractional array index
!> \param ANS Interpolated value at XI

SUBROUTINE INTERPOLATE1D_UNIFORM(LOWER,X,XI,ANS)

INTEGER, INTENT(IN) :: LOWER
REAL(EB), INTENT(IN), DIMENSION(LOWER:) :: X
REAL(EB), INTENT(IN) :: XI
REAL(EB), INTENT(OUT) :: ANS
INTEGER I, UX,LX
REAL(EB) :: FRAC

UX = UBOUND(X,1)
LX = LBOUND(X,1)

IF (XI <= LX) THEN
   ANS = X(LX)
ELSEIF (XI >= UX) THEN
   ANS = X(UX)
ELSE
   I = INT(XI)
   FRAC = XI - REAL(I,EB)
   ANS = X(I) + FRAC*(X(I+1)-X(I))
ENDIF

END SUBROUTINE INTERPOLATE1D_UNIFORM


!> \brief Randomly choose a point from a normal distribution
!> \param MEAN Mean of the normal distribution
!> \param SIGMA Standard deviation

REAL(EB) FUNCTION NORMAL(MEAN,SIGMA)

REAL(EB), INTENT(IN) :: MEAN,SIGMA
REAL(EB) :: TMP,FAC,GSAVE,RSQ,R1,R2
REAL     :: RN
INTEGER :: FLAG
SAVE FLAG,GSAVE
DATA FLAG /0/

IF (FLAG==0) THEN
   RSQ=2.0_EB
   DO WHILE(RSQ>=1.0_EB.OR.RSQ==0.0_EB)
      CALL RANDOM_NUMBER(RN)
      R1=2.0_EB*REAL(RN,EB)-1.0_EB
      CALL RANDOM_NUMBER(RN)
      R2=2.0_EB*REAL(RN,EB)-1.0_EB
      RSQ=R1*R1+R2*R2
   ENDDO
   FAC=SQRT(-2.0_EB*LOG(RSQ)/RSQ)
   GSAVE=R1*FAC
   TMP=R2*FAC
   FLAG=1
ELSE
   TMP=GSAVE
   FLAG=0
ENDIF
NORMAL=TMP*SIGMA+MEAN

END FUNCTION NORMAL


!> \brief Compute the cross product of two triplets, A x B = C
!> \param C The resulting vector
!> \param A First vector
!> \param B Second vector

SUBROUTINE CROSS_PRODUCT(C,A,B)

REAL(EB), INTENT(IN) :: A(3),B(3)
REAL(EB), INTENT(OUT) :: C(3)

C(1) = A(2)*B(3)-A(3)*B(2)
C(2) = A(3)*B(1)-A(1)*B(3)
C(3) = A(1)*B(2)-A(2)*B(1)

END SUBROUTINE CROSS_PRODUCT


!> \brief Randomly choose a value from the distribution with given CDF
!> \param CDF Cumulative Distribution Function
!> \param VAR Independent variable
!> \param NPTS Number of points in the CDF
!> \param CHOICE Randomly chosen value

SUBROUTINE RANDOM_CHOICE(CDF,VAR,NPTS,CHOICE)

INTEGER,  INTENT(IN)  :: NPTS
REAL(EB), INTENT(IN)  :: CDF(0:NPTS),VAR(0:NPTS)
REAL(EB), INTENT(OUT) :: CHOICE
INTEGER  :: IT
REAL(EB) :: CFRAC,A,B
REAL(EB) :: RN
REAL     :: RN2

CALL RANDOM_NUMBER(RN2)
RN = REAL(RN2,EB)
A = MINVAL(CDF)
B = MAXVAL(CDF)
RN = A + (B-A)*RN

CDF_LOOP: DO IT=1,NPTS
   IF (CDF(IT) > RN) THEN
      CFRAC  = (RN-CDF(IT-1))/(CDF(IT)-CDF(IT-1))
      CHOICE = VAR(IT-1) + (VAR(IT)-VAR(IT-1))*CFRAC
      EXIT CDF_LOOP
   ENDIF
ENDDO CDF_LOOP

END SUBROUTINE RANDOM_CHOICE


REAL(EB) FUNCTION MINMOD2(X,Y)
REAL(EB), INTENT(IN) :: X,Y
MINMOD2 = 0.5_EB*(SIGN(1._EB,X)+SIGN(1._EB,Y))*MIN(ABS(X),ABS(Y))
END FUNCTION MINMOD2


REAL(EB) FUNCTION MINMOD4(W,X,Y,Z)
REAL(EB), INTENT(IN) :: W,X,Y,Z
MINMOD4 = 0.125_EB*(SIGN(1._EB,W)+SIGN(1._EB,X))* &
          ABS( (SIGN(1._EB,W)+SIGN(1._EB,Y))*(SIGN(1._EB,W)+SIGN(1._EB,Z)) )* &
          MIN(ABS(W),ABS(X),ABS(Y),ABS(Z))
END FUNCTION MINMOD4


!> \brief Generate pairs of normally distributed pseudo-random numbers with zero mean and unit variance based on the
!> Box-Muller transformation.
!> \param Z0 Output value 1
!> \param Z1 Output value 2

SUBROUTINE BOX_MULLER(Z0,Z1)

REAL(EB), INTENT(OUT) :: Z0,Z1
REAL(EB) :: U1,U2,A

CALL RANDOM_NUMBER(U1)
CALL RANDOM_NUMBER(U2)
A = SQRT(-2._EB*LOG(U1))
Z0 = A*COS(TWOPI*U2)
Z1 = A*SIN(TWOPI*U2)

END SUBROUTINE BOX_MULLER


!> \brief Flux-limiting function related to mass transfer B-number

REAL(EB) FUNCTION F_B(B)
REAL(EB), INTENT(IN) :: B

IF (B<=0._EB) THEN
   F_B = 1._EB
ELSE
   F_B = (1._EB+B)**0.7_EB*LOG(1._EB+B)/B
ENDIF
END FUNCTION F_B


!> \brief This subroutine computes the flux-limited scalar value on a face.
!> \param A Array of velocity components (m/s)
!> \param U Array of scalars
!> \param F Array of flux-limited scalars
!> \param I1 Lower I index
!> \param I2 Upper I index
!> \param J1 Lower J index
!> \param J2 Upper J index
!> \param K1 Lower K index
!> \param K2 Upper K index
!> \param IOR Orientation index (1, 2, or 3)
!> \param LIMITER Indicator of the flux limiting scheme

!> There are 6 options for flux LIMITER:
!>
!> CENTRAL_LIMITER  = 0
!> GODUNOV_LIMITER  = 1
!> SUPERBEE_LIMITER = 2
!> MINMOD_LIMITER   = 3
!> CHARM_LIMITER    = 4
!> MP5_LIMITER      = 5
!>
!> Example: x-direction (IOR=1)
!>
!>                   location of face
!>
!>                        F(I,J,K)
!>   |     o     |     o     |     o     |     o     |
!>                        A(I,J,K)
!>     U(I-1,J,K)   U(I,J,K)   U(I+1,J,K)  U(I+2,J,K)

SUBROUTINE GET_SCALAR_FACE_VALUE(A,U,F,I1,I2,J1,J2,K1,K2,IOR,LIMITER)

REAL(EB), INTENT(IN) :: A(0:,0:,0:),U(0:,0:,0:)
REAL(EB), INTENT(OUT) :: F(0:,0:,0:)
INTEGER, INTENT(IN) :: LIMITER,I1,I2,J1,J2,K1,K2,IOR
REAL(EB) :: R,B,DU_UP,DU_LOC,V(-2:2)
INTEGER :: I,J,K,IM1,JM1,KM1,IP1,JP1,KP1,IP2,JP2,KP2

SELECT CASE(IOR)
   CASE(1) ; IM1=-1 ; JM1= 0 ; KM1= 0 ; IP1=1 ; JP1=0 ; KP1=0 ; IP2=2 ; JP2=0 ; KP2=0
   CASE(2) ; IM1= 0 ; JM1=-1 ; KM1= 0 ; IP1=0 ; JP1=1 ; KP1=0 ; IP2=0 ; JP2=2 ; KP2=0
   CASE(3) ; IM1= 0 ; JM1= 0 ; KM1=-1 ; IP1=0 ; JP1=0 ; KP1=1 ; IP2=0 ; JP2=0 ; KP2=2
END SELECT

!$OMP PARALLEL IF(I2>1)
SELECT CASE(LIMITER)
   CASE(0) ! central differencing
      !$OMP DO SCHEDULE(STATIC)
      DO K=K1,K2
         DO J=J1,J2
            DO I=I1,I2
                  F(I,J,K) = 0.5_EB*(U(I,J,K) + U(I+IP1,J+JP1,K+KP1))
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(1) ! first-order upwinding
      !$OMP DO SCHEDULE(STATIC)
      DO K=K1,K2
         DO J=J1,J2
            DO I=I1,I2
               IF (A(I,J,K)>0._EB) THEN
                  F(I,J,K) = U(I,J,K)
               ELSE
                  F(I,J,K) = U(I+IP1,J+JP1,K+KP1)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(2) ! SUPERBEE, Roe (1986)
      !$OMP DO SCHEDULE(STATIC) PRIVATE(DU_UP,DU_LOC,B,R)
      DO K=K1,K2
         DO J=J1,J2
            DO I=I1,I2
               DU_LOC = U(I+IP1,J+JP1,K+KP1)-U(I,J,K)
               IF (A(I,J,K)>0._EB) THEN
                  DU_UP  = U(I,J,K)-U(I+IM1,J+JM1,K+KM1)
                  R = DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC))
                  B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
                  F(I,J,K) = U(I,J,K) + 0.5_EB*B*DU_LOC
               ELSE
                  DU_UP  = U(I+IP2,J+JP2,K+KP2)-U(I+IP1,J+JP1,K+KP1)
                  R = DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC))
                  B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
                  F(I,J,K) = U(I+IP1,J+JP1,K+KP1) - 0.5_EB*B*DU_LOC
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(3) ! MINMOD
      !$OMP DO SCHEDULE(STATIC) PRIVATE(DU_UP,DU_LOC,B,R)
      DO K=K1,K2
         DO J=J1,J2
            DO I=I1,I2
               DU_LOC = U(I+IP1,J+JP1,K+KP1)-U(I,J,K)
               IF (A(I,J,K)>0._EB) THEN
                  DU_UP  = U(I,J,K)-U(I+IM1,J+JM1,K+KM1)
                  R = DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC))
                  B = MAX(0._EB,MIN(R,1._EB))
                  F(I,J,K) = U(I,J,K) + 0.5_EB*B*DU_LOC
               ELSE
                  DU_UP  = U(I+IP2,J+JP2,K+KP2)-U(I+IP1,J+JP1,K+KP1)
                  R = DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC))
                  B = MAX(0._EB,MIN(R,1._EB))
                  F(I,J,K) = U(I+IP1,J+JP1,K+KP1) - 0.5_EB*B*DU_LOC
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(4) ! CHARM
      !$OMP DO SCHEDULE(STATIC) PRIVATE(DU_UP,DU_LOC,B,R)
      DO K=K1,K2
         DO J=J1,J2
            DO I=I1,I2
               DU_LOC = U(I+IP1,J+JP1,K+KP1)-U(I,J,K)
               IF (A(I,J,K)>0._EB) THEN
                  DU_UP  = U(I,J,K)-U(I+IM1,J+JM1,K+KM1)
                  R = 0._EB
                  B = 0._EB
                  IF (ABS(DU_UP)>TWO_EPSILON_EB) R = DU_LOC/DU_UP
                  IF (R>0._EB) B = R*(3._EB*R+1._EB)/((R+1._EB)**2)
                  F(I,J,K) = U(I,J,K) + 0.5_EB*B*DU_UP
               ELSE
                  DU_UP  = U(I+IP2,J+JP2,K+KP2)-U(I+IP1,J+JP1,K+KP1)
                  R = 0._EB
                  B = 0._EB
                  IF (ABS(DU_UP)>TWO_EPSILON_EB) R = DU_LOC/DU_UP
                  IF (R>0._EB) B = R*(3._EB*R+1._EB)/((R+1._EB)**2)
                  F(I,J,K) = U(I+IP1,J+JP1,K+KP1) - 0.5_EB*B*DU_UP
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(5) ! MP5, Suresh and Huynh (1997)
      DO K=K1,K2
         DO J=J1,J2
            DO I=I1,I2
               IF (A(I,J,K)>0._EB) THEN
                  V = (/2._EB*U(I+IM1,J+JM1,K+KM1)-U(I,J,K),&
                       U(I+IM1,J+JM1,K+KM1),U(I,J,K),U(I+IP1,J+JP1,K+KP1),U(I+IP2,J+JP2,K+KP2)/)
                  F(I,J,K) = MP5()
               ELSE
                  V = (/2._EB*U(I+IP2,J+JP2,K+KP2)-U(I+IP1,J+JP1,K+KP1),&
                       U(I+IP2,J+JP2,K+KP2),U(I+IP1,J+JP1,K+KP1),U(I,J,K),U(I+IM1,J+JM1,K+KM1)/)
                  F(I,J,K) = MP5()
               ENDIF
            ENDDO
         ENDDO
      ENDDO
END SELECT
!$OMP END PARALLEL

CONTAINS

REAL(EB) FUNCTION MP5()
REAL(EB), PARAMETER :: B1 = 0.016666666666667_EB, B2 = 1.333333333333_EB, ALPHA=4._EB, EPSM=1.E-10_EB
REAL(EB) :: VOR,VMP,DJM1,DJ,DJP1,DM4JPH,DM4JMH,VUL,VAV,VMD,VLC,VMIN,VMAX

! Monotonicity preserving 5th-order scheme (MP5) of Suresh and Huynh, JCP 136, 83-99 (1997)

VOR = B1*(2._EB*V(-2)-13._EB*V(-1)+47._EB*V(0)+27._EB*V(1)-3._EB*V(2))
VMP = V(0) + MINMOD2(V(1)-V(0),ALPHA*(V(0)-V(-1)))
IF ((VOR-V(0))*(VOR-VMP)<EPSM) THEN
   MP5=VOR
ELSE
   DJM1 = V(-2)-2._EB*V(-1)+V(0)
   DJ   = V(-1)-2._EB*V(0) +V(1)
   DJP1 = V(0) -2._EB*V(1) +V(2)
   DM4JPH = MINMOD4(4._EB*DJ-DJP1,4._EB*DJP1-DJ,DJ,DJP1)
   DM4JMH = MINMOD4(4._EB*DJ-DJM1,4._EB*DJM1-DJ,DJ,DJM1)
   VUL = V(0) + ALPHA*(V(0)-V(-1))
   VAV = 0.5_EB*(V(0)+V(1))
   VMD = VAV - 0.5_EB*DM4JPH
   VLC = V(0) + 0.5_EB*(V(0)-V(-1)) + B2*DM4JMH
   VMIN = MAX(MIN(V(0),V(1),VMD),MIN(V(0),VUL,VLC))
   VMAX = MIN(MAX(V(0),V(1),VMD),MAX(V(0),VUL,VLC))
   MP5 = VOR + MINMOD2(VMIN-VOR,VMAX-VOR)
ENDIF

END FUNCTION MP5

END SUBROUTINE GET_SCALAR_FACE_VALUE


!> \brief Random fluctuation, Theta'(t+dt) = R^2*Theta'(t) + Normal(0,sqrt(1-R^2)*SIGMA) ; R = exp(-dt/TAU)
!> \param SIGMA Standard deviation of time series
!> \param TAU Time scale or period of the time function (s)

SUBROUTINE RANDOM_WIND_FLUCTUATIONS(SIGMA,TAU)

USE TYPES, ONLY: RESERVED_RAMPS_TYPE,RESERVED_RAMPS,N_RESERVED_RAMPS
USE GLOBAL_CONSTANTS, ONLY: T_END,T_BEGIN
TYPE(RESERVED_RAMPS_TYPE), POINTER :: RRP
REAL(EB), INTENT(IN) :: SIGMA,TAU
REAL(EB) :: LCC,DT_THETA
INTEGER :: I

N_RESERVED_RAMPS = N_RESERVED_RAMPS + 1
RRP => RESERVED_RAMPS(N_RESERVED_RAMPS)
RRP%NUMBER_DATA_POINTS = 1001
ALLOCATE(RRP%INDEPENDENT_DATA(RRP%NUMBER_DATA_POINTS))
ALLOCATE(RRP%DEPENDENT_DATA(RRP%NUMBER_DATA_POINTS))
DT_THETA = (T_END-T_BEGIN)/REAL(RRP%NUMBER_DATA_POINTS-1,EB)
RRP%INDEPENDENT_DATA(1) = T_BEGIN
RRP%DEPENDENT_DATA(1)   = 0._EB
LCC = EXP(-DT_THETA/TAU)  ! Lagrangian Correlation Coefficient, R
DO I=2,RRP%NUMBER_DATA_POINTS
   RRP%INDEPENDENT_DATA(I) = RRP%INDEPENDENT_DATA(I-1) + DT_THETA
   RRP%DEPENDENT_DATA(I)   = LCC**2*RRP%DEPENDENT_DATA(I-1) + NORMAL(0._EB,SQRT(1._EB-LCC**2)*SIGMA)
ENDDO

END SUBROUTINE RANDOM_WIND_FLUCTUATIONS

END MODULE MATH_FUNCTIONS


!> \brief Functions for physical quantities

MODULE PHYSICAL_FUNCTIONS

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
IMPLICIT NONE (TYPE,EXTERNAL)

LOGICAL :: DEF_QREF_ARRAYS = .FALSE.
REAL(EB) :: QFLAME(1:5,0:20,0:30,1:6),ABSF(1:5,0:20,0:30,1:6)

PRIVATE :: DEF_QREF_ARRAYS,QFLAME,ABSF

CONTAINS


!> \brief Check if the species mass fractions are in bounds
!> \param ZZ_IN Array of mass fractions

LOGICAL FUNCTION IS_REALIZABLE(ZZ_IN)

REAL(EB), INTENT(IN) :: ZZ_IN(1:N_TRACKED_SPECIES)
REAL(EB), PARAMETER:: ZERO_MINUS=-EPSILON(0._FB),ONE_PLUS=1._EB+EPSILON(1._FB)

IF (ANY(ZZ_IN<ZERO_MINUS) .OR. SUM(ZZ_IN)>ONE_PLUS) THEN
   IS_REALIZABLE=.FALSE.
ELSE
   IS_REALIZABLE=.TRUE.
ENDIF

END FUNCTION IS_REALIZABLE


!> \brief Clip mass fractions between zero and one and redistribute clipped mass to most abundant species
!> \param ZZ_GET Array of mass fractions

SUBROUTINE GET_REALIZABLE_MF(ZZ_GET)

REAL(EB), INTENT(INOUT) :: ZZ_GET(1:N_TRACKED_SPECIES)
REAL(EB) :: SUM_OTHER_SPECIES
INTEGER :: N_ZZ_MAX

! clip mass fractions
ZZ_GET=MAX(0._EB,MIN(1._EB,ZZ_GET))

! absorb all error in most abundant species
N_ZZ_MAX = MAXLOC(ZZ_GET,1)

SUM_OTHER_SPECIES = SUM(ZZ_GET) - ZZ_GET(N_ZZ_MAX)
ZZ_GET(N_ZZ_MAX) = 1._EB - SUM_OTHER_SPECIES

END SUBROUTINE GET_REALIZABLE_MF


!> \brief Determine the mass fraction of a primitive species
!> \param Z_IN Array of lumped mass fractions
!> \param INDEX Index of desired primitive species
!> \param Y_OUT Mass fraction of desired primitive species

SUBROUTINE GET_MASS_FRACTION(Z_IN,INDEX,Y_OUT)

INTEGER, INTENT(IN) :: INDEX
REAL(EB) :: Z_IN(1:N_TRACKED_SPECIES)
REAL(EB), INTENT(OUT) :: Y_OUT

Y_OUT = DOT_PRODUCT(Z2Y(INDEX,1:N_TRACKED_SPECIES),Z_IN)
Y_OUT = MIN(1._EB,MAX(0._EB,Y_OUT))

END SUBROUTINE GET_MASS_FRACTION


!> \brief Determine the mass fractions of all primitive species
!> \param Z_IN Array of lumped species mass fractions
!> \param Y_OUT Array of primitive species mass fractions

SUBROUTINE GET_MASS_FRACTION_ALL(Z_IN,Y_OUT)

REAL(EB), INTENT(IN)  :: Z_IN(1:N_TRACKED_SPECIES)
REAL(EB), INTENT(OUT) :: Y_OUT(1:N_SPECIES)
INTEGER :: I

DO I=1,N_SPECIES
   Y_OUT(I) = DOT_PRODUCT(Z2Y(I,1:N_TRACKED_SPECIES),Z_IN)
ENDDO

Y_OUT = MIN(1._EB,MAX(0._EB,Y_OUT))

END SUBROUTINE GET_MASS_FRACTION_ALL


!> \brief Compute ratio of molecular weight of all gas components except the one specified to the mol. wgt. of the one specified
!> \param INDEX_IN Index of the gas species
!> \param MW_RATIO W_all/W_in
!> \param Y_IN Optional mass fraction of primitive species
!> \param Z_IN Optional mass fraction of lumped species

SUBROUTINE GET_MW_RATIO(INDEX_IN,MW_RATIO,Y_IN,Z_IN)

REAL(EB), INTENT(IN), OPTIONAL :: Y_IN(1:N_SPECIES), Z_IN(1:N_TRACKED_SPECIES)
INTEGER, INTENT(IN):: INDEX_IN
REAL(EB), INTENT(OUT) :: MW_RATIO
INTEGER:: NS

MW_RATIO = 0._EB
IF (PRESENT(Y_IN)) THEN
   IF (ABS(Y_IN(INDEX_IN)-1._EB) > TWO_EPSILON_EB) THEN
      DO NS=1,N_SPECIES
         IF (NS==INDEX_IN) CYCLE
         MW_RATIO = MW_RATIO + Y_IN(NS)/SPECIES(NS)%MW
      ENDDO
      IF (MW_RATIO<=TWO_EPSILON_EB) THEN
         MW_RATIO=SPECIES_MIXTURE(1)%MW
      ELSE
         MW_RATIO = (1._EB-Y_IN(INDEX_IN))/MW_RATIO
      ENDIF
   ELSE
      MW_RATIO=SPECIES_MIXTURE(1)%MW
   ENDIF
   MW_RATIO = MW_RATIO/SPECIES(INDEX_IN)%MW
ELSE
   IF (ABS(Z_IN(INDEX_IN)-1._EB) > TWO_EPSILON_EB) THEN
      DO NS=1,N_TRACKED_SPECIES
         IF (NS==INDEX_IN) CYCLE
         MW_RATIO = MW_RATIO + Z_IN(NS)/SPECIES_MIXTURE(NS)%MW
      ENDDO
      IF (MW_RATIO<=TWO_EPSILON_EB) THEN
         MW_RATIO=SPECIES_MIXTURE(1)%MW
      ELSE
         MW_RATIO = (1._EB-Z_IN(INDEX_IN))/MW_RATIO
      ENDIF
   ELSE
      MW_RATIO=SPECIES_MIXTURE(1)%MW
   ENDIF
   MW_RATIO = MW_RATIO/SPECIES_MIXTURE(INDEX_IN)%MW
ENDIF

END SUBROUTINE GET_MW_RATIO


SUBROUTINE GET_EQUIL_DATA(MW,TMP_L,PRES_IN,H_V,H_V_EFF,T_BOIL_EFF,X_EQ,H_V_DATA)

REAL(EB), INTENT(IN):: MW,TMP_L,PRES_IN
REAL(EB), INTENT(IN) :: H_V_DATA(:)
REAL(EB), INTENT(INOUT):: T_BOIL_EFF
REAL(EB), INTENT(OUT):: H_V,H_V_EFF,X_EQ
REAL(EB):: DHOR

H_V=H_V_DATA(MIN(I_MAX_TEMP,NINT(TMP_L)))
H_V_EFF=H_V_DATA(MIN(I_MAX_TEMP,NINT(T_BOIL_EFF)))
DHOR = H_V_EFF*MW/R0
T_BOIL_EFF = MAX(0._EB,DHOR*T_BOIL_EFF/(DHOR-T_BOIL_EFF*LOG(PRES_IN/P_STP)+TWO_EPSILON_EB))
H_V_EFF=H_V_DATA(MIN(I_MAX_TEMP,NINT(T_BOIL_EFF)))
H_V_EFF = 0.5_EB*(H_V+H_V_EFF)
X_EQ = MIN(1._EB,EXP(H_V_EFF*MW/R0*(1._EB/T_BOIL_EFF-1._EB/TMP_L)))

END SUBROUTINE GET_EQUIL_DATA


!> \brief Compute volume fraction of a primitive species
!> \param Y_INDEX Index of primitive species
!> \param ZZ_ARRAY Array of lumped species
!> \param R_MIX R/W of species mixture

REAL(EB) FUNCTION GET_VOLUME_FRACTION(Y_INDEX,ZZ_ARRAY,R_MIX)

INTEGER, INTENT(IN) :: Y_INDEX
REAL(EB), INTENT(IN) :: R_MIX,ZZ_ARRAY(:)
REAL(EB) :: MASS_FRACTION,RCON

CALL GET_MASS_FRACTION(ZZ_ARRAY,Y_INDEX,MASS_FRACTION)
RCON = SPECIES(Y_INDEX)%RCON
GET_VOLUME_FRACTION = RCON*MASS_FRACTION/R_MIX

END FUNCTION GET_VOLUME_FRACTION


!> \brief Determine molecular weight of the gas mixture
!> \param Z_IN Array of lumped species mass fractions
!> \param MW_OUT Average molecular weight (g/mol)

SUBROUTINE GET_MOLECULAR_WEIGHT(Z_IN,MW_OUT)

REAL(EB), INTENT(IN)  :: Z_IN(1:N_TRACKED_SPECIES)
REAL(EB), INTENT(OUT) :: MW_OUT

MW_OUT =  1._EB/DOT_PRODUCT(MWR_Z,Z_IN)

END SUBROUTINE GET_MOLECULAR_WEIGHT


!> \brief Compute R/W for a gas mixture
!> \param Z_IN Array of lumped species mass fractions
!> \param RSUM_OUT R/W_avg

SUBROUTINE GET_SPECIFIC_GAS_CONSTANT(Z_IN,RSUM_OUT)

REAL(EB) :: Z_IN(1:N_TRACKED_SPECIES)
REAL(EB), INTENT(OUT) :: RSUM_OUT

RSUM_OUT =  R0 * DOT_PRODUCT(MWR_Z,Z_IN)

END SUBROUTINE GET_SPECIFIC_GAS_CONSTANT


!> \brief Get specific heat of the gas mixture at a specified temperature
!> \param Z_IN Array of lumped species mass fractions
!> \param CP_OUT Specific heat of the mixture (J/kg/K)
!> \param TMPG Gas mixture temperature (K)

SUBROUTINE GET_SPECIFIC_HEAT(Z_IN,CP_OUT,TMPG)

INTEGER :: ITMP
REAL(EB), INTENT(IN) :: TMPG
REAL(EB) :: Z_IN(1:N_TRACKED_SPECIES)
REAL(EB), INTENT(OUT) :: CP_OUT

ITMP = MIN(I_MAX_TEMP,NINT(TMPG))
CP_OUT  = DOT_PRODUCT(CP_Z(ITMP,1:N_TRACKED_SPECIES),Z_IN)

END SUBROUTINE GET_SPECIFIC_HEAT

SUBROUTINE GET_SPECIFIC_HEAT_INTERP(Z_IN,CP_OUT,TMPG)

INTEGER :: ITMP
REAL(EB), INTENT(IN) :: TMPG
REAL(EB) :: Z_IN(1:N_TRACKED_SPECIES)
REAL(EB), INTENT(OUT) :: CP_OUT

ITMP = MIN(I_MAX_TEMP-1,INT(TMPG))
CP_OUT = DOT_PRODUCT(CP_Z(ITMP,1:N_TRACKED_SPECIES),Z_IN)+(TMPG-REAL(ITMP,EB))* &
          DOT_PRODUCT(CP_Z(ITMP+1,1:N_TRACKED_SPECIES)-CP_Z(ITMP,1:N_TRACKED_SPECIES),Z_IN)

END SUBROUTINE GET_SPECIFIC_HEAT_INTERP

SUBROUTINE GET_SPECIFIC_HEAT_Z(N,CP_OUT,TMPG)

INTEGER :: ITMP
REAL(EB), INTENT(IN) :: TMPG
INTEGER, INTENT(IN) :: N
REAL(EB), INTENT(OUT) :: CP_OUT

ITMP = MIN(I_MAX_TEMP-1,INT(TMPG))
CP_OUT  = CP_Z(ITMP,N) + (TMPG-REAL(ITMP,EB))*(CP_Z(ITMP+1,N)-CP_Z(ITMP,N))

END SUBROUTINE  GET_SPECIFIC_HEAT_Z

SUBROUTINE GET_SPECIFIC_HEAT_TMP_DERIVATIVE(Z_IN,DCPDT,TMPG)

INTEGER :: ITMP
REAL(EB), INTENT(IN) :: TMPG
REAL(EB) :: Z_IN(1:N_TRACKED_SPECIES)
REAL(EB), INTENT(OUT) :: DCPDT
REAL(EB) :: CP_OUT_1, CP_OUT_2

ITMP = MIN(I_MAX_TEMP-1,INT(TMPG))
CP_OUT_1  = DOT_PRODUCT(CP_Z(ITMP,1:N_TRACKED_SPECIES),Z_IN)
CP_OUT_2  = DOT_PRODUCT(CP_Z(ITMP+1,1:N_TRACKED_SPECIES),Z_IN)
DCPDT = (CP_OUT_2-CP_OUT_1)

END SUBROUTINE GET_SPECIFIC_HEAT_TMP_DERIVATIVE

!> \brief Get sensible enthalpy of the gas mixture at a specified temperature
!> \param Z_IN Array of lumped species mass fractions
!> \param H_S_OUT Specific heat of the mixture (J/kg)
!> \param TMPG Gas mixture temperature (K)

SUBROUTINE GET_SENSIBLE_ENTHALPY(Z_IN,H_S_OUT,TMPG)

INTEGER :: ITMP
REAL(EB), INTENT(IN) :: TMPG
REAL(EB) :: Z_IN(1:N_TRACKED_SPECIES)
REAL(EB), INTENT(OUT) :: H_S_OUT

ITMP = MIN(I_MAX_TEMP-1,INT(TMPG))

H_S_OUT = DOT_PRODUCT(H_SENS_Z(ITMP,1:N_TRACKED_SPECIES),Z_IN)+(TMPG-REAL(ITMP,EB))* &
          DOT_PRODUCT(H_SENS_Z(ITMP+1,1:N_TRACKED_SPECIES)-H_SENS_Z(ITMP,1:N_TRACKED_SPECIES),Z_IN)

END SUBROUTINE GET_SENSIBLE_ENTHALPY

!> \brief Get average specific heat of the gas mixture up to a specified temperature
!> \param Z_IN Array of lumped species mass fractions
!> \param CPBAR_OUT Average specific heat of the mixture (J/kg/K)
!> \param TMPG Gas mixture temperature (K)

SUBROUTINE GET_AVERAGE_SPECIFIC_HEAT(Z_IN,CPBAR_OUT,TMPG)

INTEGER :: ITMP
REAL(EB), INTENT(IN) :: TMPG
REAL(EB) :: Z_IN(1:N_TRACKED_SPECIES)
REAL(EB), INTENT(OUT) :: CPBAR_OUT

ITMP = MIN(I_MAX_TEMP,NINT(TMPG))
CPBAR_OUT = DOT_PRODUCT(CPBAR_Z(ITMP,1:N_TRACKED_SPECIES),Z_IN)

END SUBROUTINE GET_AVERAGE_SPECIFIC_HEAT


!> \brief Get sensible enthalpy of lumped species N at a specified temperature
!> \param N Index of lumped species
!> \param TMPG Gas mixture temperature (K)
!> \param H_S Sensible enthalpy of lumped species N (J/kg)

SUBROUTINE GET_SENSIBLE_ENTHALPY_Z(N,TMPG,H_S)

INTEGER :: ITMP
REAL(EB), INTENT(IN) :: TMPG
INTEGER, INTENT(IN) :: N
REAL(EB), INTENT(OUT) :: H_S

ITMP = MIN(I_MAX_TEMP-1,INT(TMPG))
H_S = H_SENS_Z(ITMP,N)+(TMPG-REAL(ITMP,EB))*(H_SENS_Z(ITMP+1,N)-H_SENS_Z(ITMP,N))

END SUBROUTINE GET_SENSIBLE_ENTHALPY_Z


!> \brief Get average specific heat of lumped species N up to a specified temperature
!> \param N Index of lumped species
!> \param TMPG Gas mixture temperature (K)
!> \param CPBAR_OUT Average specific heat of lumped species N (J/kg/K)

SUBROUTINE GET_AVERAGE_SPECIFIC_HEAT_Z(N,TMPG,CPBAR_OUT)

INTEGER :: ITMP
REAL(EB), INTENT(IN) :: TMPG
INTEGER, INTENT(IN) :: N
REAL(EB), INTENT(OUT) :: CPBAR_OUT

ITMP = MIN(I_MAX_TEMP,NINT(TMPG))
CPBAR_OUT = CPBAR_Z(ITMP,N)

END SUBROUTINE GET_AVERAGE_SPECIFIC_HEAT_Z


!> \brief Get thermal conductivity of gas mixture at a specified temperature
!> \param Z_IN Array of lumped species mass fractions
!> \param K_OUT Conductivity of gas mixture (W/m/K)
!> \param TMPG Gas mixture temperature (K)

SUBROUTINE GET_CONDUCTIVITY(Z_IN,K_OUT,TMPG)

REAL(EB), INTENT(IN) :: Z_IN(1:N_TRACKED_SPECIES), TMPG
REAL(EB), INTENT(OUT) :: K_OUT
INTEGER :: ITMP

ITMP = MIN(I_MAX_TEMP,NINT(TMPG))
K_OUT = DOT_PRODUCT(K_RSQMW_Z(ITMP,1:N_TRACKED_SPECIES),Z_IN)/DOT_PRODUCT(Z_IN,RSQ_MW_Z)

END SUBROUTINE GET_CONDUCTIVITY


!> \brief Get enthalpy (J/kg) of a particle at a specified uniform temperature
!> \param I_LPC Index of particle class
!> \param TMP_S Particle temperature (K)

REAL(EB) FUNCTION GET_PARTICLE_ENTHALPY(I_LPC,TMP_S)
USE MATH_FUNCTIONS, ONLY: INTERPOLATE1D_UNIFORM
REAL(EB), INTENT(IN) :: TMP_S
REAL(EB) :: RHO,RHO_H,VOL,DTMP,H_S,THICKNESS
INTEGER :: I,N,ITMP,I_GRAD
INTEGER, INTENT(IN) :: I_LPC
TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
TYPE(SURFACE_TYPE), POINTER :: SF=>NULL()
TYPE(MATERIAL_TYPE), POINTER :: ML=>NULL()

LPC=>LAGRANGIAN_PARTICLE_CLASS(I_LPC)

IF (LPC%LIQUID_DROPLET) THEN
   CALL INTERPOLATE1D_UNIFORM(LBOUND(SPECIES(LPC%Y_INDEX)%H_L,1),SPECIES(LPC%Y_INDEX)%H_L,TMP_S,GET_PARTICLE_ENTHALPY)
ELSE
   SF=>SURFACE(LPC%SURF_INDEX)
   SELECT CASE(SF%GEOMETRY)
      CASE(SURF_CARTESIAN)                          ; I_GRAD = 1
      CASE(SURF_CYLINDRICAL,SURF_INNER_CYLINDRICAL) ; I_GRAD = 2
      CASE(SURF_SPHERICAL)                          ; I_GRAD = 3
   END SELECT
   RHO_H = 0._EB
   RHO = 0._EB
   ITMP = MIN(I_MAX_TEMP-1,INT(TMP_S))
   DTMP = TMP_S-REAL(ITMP,EB)
   THICKNESS = SUM(SF%LAYER_THICKNESS)
   DO I=1,SUM(SF%N_LAYER_CELLS)
      IF (SF%GEOMETRY==SURF_INNER_CYLINDRICAL) THEN
         VOL = (SF%INNER_RADIUS+SF%X_S(I))**I_GRAD - (SF%INNER_RADIUS+SF%X_S(I-1))**I_GRAD
      ELSE
         VOL = (THICKNESS+SF%INNER_RADIUS-SF%X_S(I-1))**I_GRAD - (THICKNESS+SF%INNER_RADIUS-SF%X_S(I))**I_GRAD
      ENDIF
      MATL_REMESH: DO N=1,SF%N_MATL
         IF (SF%RHO_0(I,N)<=TWO_EPSILON_EB) CYCLE MATL_REMESH
         ML  => MATERIAL(SF%MATL_INDEX(N))
         H_S = ML%H(ITMP)+DTMP*(ML%H(ITMP+1)-ML%H(ITMP))
         RHO_H = RHO_H + SF%RHO_0(I,N) * H_S * VOL
         RHO = RHO + SF%RHO_0(I,N) * VOL
      ENDDO MATL_REMESH
   ENDDO
   GET_PARTICLE_ENTHALPY = RHO_H/RHO
ENDIF

END FUNCTION GET_PARTICLE_ENTHALPY


!> \brief Get viscosity of gas mixture at a specified temperature
!> \param Z_IN Array of lumped species mass fractions
!> \param MU_OUT Viscosity of gas mixture (kg/m/s)
!> \param TMPG Gas mixture temperature (K)

SUBROUTINE GET_VISCOSITY(Z_IN,MU_OUT,TMPG)

INTEGER :: ITMP
REAL(EB), INTENT(IN)  :: TMPG,Z_IN(1:N_TRACKED_SPECIES)
REAL(EB), INTENT(OUT) :: MU_OUT

ITMP = MIN(I_MAX_TEMP,NINT(TMPG))

MU_OUT = DOT_PRODUCT(MU_RSQMW_Z(ITMP,1:N_TRACKED_SPECIES),Z_IN)/DOT_PRODUCT(Z_IN,RSQ_MW_Z)

END SUBROUTINE GET_VISCOSITY


SUBROUTINE GET_ENTHALPY(Z_IN,H_OUT,TMPG)

INTEGER :: ITMP
REAL(EB), INTENT(IN) :: TMPG
REAL(EB) :: Z_IN(1:N_TRACKED_SPECIES),DTMP
REAL(EB), INTENT(OUT) :: H_OUT

IF (TMPG>=REAL(I_MAX_TEMP,EB)) THEN
   H_OUT = DOT_PRODUCT(CPBAR_Z(I_MAX_TEMP,1:N_TRACKED_SPECIES),Z_IN)*REAL(I_MAX_TEMP,EB) + &
           DOT_PRODUCT(CP_Z(I_MAX_TEMP,1:N_TRACKED_SPECIES),Z_IN)*(TMPG-REAL(I_MAX_TEMP,EB))
ELSE
   ITMP = INT(TMPG)
   DTMP = TMPG-REAL(ITMP)
   H_OUT = DOT_PRODUCT(CPBAR_Z(ITMP,1:N_TRACKED_SPECIES),Z_IN)
   H_OUT = H_OUT+DTMP*(DOT_PRODUCT(CPBAR_Z(ITMP+1,1:N_TRACKED_SPECIES),Z_IN)-H_OUT)
   H_OUT = H_OUT*TMPG
ENDIF

END SUBROUTINE GET_ENTHALPY


SUBROUTINE GET_ENTHALPY_Z(N,TMPG,H_OUT)

INTEGER :: ITMP
REAL(EB), INTENT(IN) :: TMPG
INTEGER, INTENT(IN) :: N
REAL(EB), INTENT(OUT) :: H_OUT

IF (TMPG>=REAL(I_MAX_TEMP,EB)) THEN
   H_OUT = CPBAR_Z(I_MAX_TEMP,N)*REAL(I_MAX_TEMP,EB) + &
           CP_Z(I_MAX_TEMP,N)*(TMPG-REAL(I_MAX_TEMP,EB))
ELSE
   ITMP = MIN(I_MAX_TEMP-1,INT(TMPG))
   H_OUT = CPBAR_Z(ITMP,N)+(TMPG-REAL(ITMP,EB))*(CPBAR_Z(ITMP+1,N)-CPBAR_Z(ITMP,N))
   H_OUT = H_OUT*TMPG
ENDIF

END SUBROUTINE GET_ENTHALPY_Z


!> \brief Compute surface emissivity
!> \param ONE_D Pointer to BOUNDARY_ONE_D derived type variable
!> \param NODE_INDEX Interior node index

SUBROUTINE GET_EMISSIVITY(ONE_D,NODE_INDEX,EMISSIVITY)

INTEGER, INTENT(IN) :: NODE_INDEX
INTEGER :: N
REAL(EB), INTENT(OUT) :: EMISSIVITY
REAL(EB) :: VOLSUM
TYPE (BOUNDARY_ONE_D_TYPE), POINTER :: ONE_D
TYPE (MATERIAL_TYPE), POINTER :: ML

EMISSIVITY = 0._EB
VOLSUM     = 0._EB

DO N=1,ONE_D%N_MATL
   IF (ONE_D%MATL_COMP(N)%RHO(NODE_INDEX)<=TWO_EPSILON_EB) CYCLE
   ML => MATERIAL(ONE_D%MATL_INDEX(N))
   VOLSUM  = VOLSUM  + ONE_D%MATL_COMP(N)%RHO(NODE_INDEX)/ML%RHO_S
   EMISSIVITY = EMISSIVITY + ONE_D%MATL_COMP(N)%RHO(NODE_INDEX)*ML%EMISSIVITY/ML%RHO_S
ENDDO
IF (VOLSUM > 0._EB) EMISSIVITY = EMISSIVITY/VOLSUM

END SUBROUTINE GET_EMISSIVITY


!> \brief Extract temperature given species mass fractions and enthalpy

SUBROUTINE GET_TEMPERATURE(TMP,ETOT,ZZ_GET)
REAL(EB), INTENT(IN) :: ETOT, ZZ_GET(1:N_TRACKED_SPECIES)
REAL(EB), INTENT(INOUT) :: TMP
INTEGER :: ITCOUNT
REAL(EB) :: CP, CP2, DCPDT, HGAS, TGUESS

TGUESS = TMP
ITCOUNT = 0
CP_LOOP: DO ! Newton method to find solution of T (and hence cpbar) from enthalpy
   ITCOUNT = ITCOUNT + 1
   CALL GET_ENTHALPY(ZZ_GET,HGAS,TGUESS)
   CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP,TGUESS)
   IF (TGUESS>1._EB) THEN
      CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP2,TGUESS-1._EB)
      DCPDT = CP - CP2
   ELSE
      CALL GET_AVERAGE_SPECIFIC_HEAT(ZZ_GET,CP2,TGUESS+1._EB)
      DCPDT = CP2- CP
   ENDIF
   CP = HGAS / TGUESS
   TMP =TGUESS+(ETOT-HGAS)/(CP+TGUESS*DCPDT)
   IF (ABS(TMP - TGUESS) < TWO_EPSILON_EB) EXIT CP_LOOP
   IF (ABS(TMP - TGUESS)/(TMP+TWO_EPSILON_EB) < 0.0005_EB) EXIT CP_LOOP
   IF (ITCOUNT > 10) THEN
      TMP = 0.5_EB*(TMP+TGUESS)
      EXIT CP_LOOP
   ENDIF
   TGUESS = MAX(0._EB,TMP)
ENDDO CP_LOOP

END SUBROUTINE GET_TEMPERATURE

SUBROUTINE MOLAR_CONC_TO_MASS_FRAC(CC_IN,ZZ_OUT)
   REAL(EB) :: CC_IN(1:N_TRACKED_SPECIES),ZZ_OUT(1:N_TRACKED_SPECIES)
   ZZ_OUT(1:N_TRACKED_SPECIES) = CC_IN(1:N_TRACKED_SPECIES)*SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%MW
   ! Check for negative mass fraction, and rescale to accomodate negative values
   WHERE(ZZ_OUT<0._EB) ZZ_OUT=0._EB
   ZZ_OUT = ZZ_OUT / SUM(ZZ_OUT)
END SUBROUTINE MOLAR_CONC_TO_MASS_FRAC

SUBROUTINE CALC_EQUIV_RATIO (ZZ, EQUIV)
REAL(EB), INTENT(IN) :: ZZ(N_TRACKED_SPECIES)
REAL(EB), INTENT(OUT) :: EQUIV
REAL(EB) :: NUMER, DENOM
INTEGER :: NS

NUMER = 0._EB
DENOM = 0._EB

DO NS=1, N_TRACKED_SPECIES
   NUMER = NUMER + ZZ(NS)*SPECIES_MIXTURE(NS)%OXR
   DENOM = DENOM + ZZ(NS)*SPECIES_MIXTURE(NS)%OXA
ENDDO

IF (DENOM < TWO_EPSILON_EB) THEN
   EQUIV= 0._EB;
ELSE
   EQUIV = NUMER/ DENOM;
ENDIF
END SUBROUTINE CALC_EQUIV_RATIO


REAL(EB) FUNCTION DRAG(RE,DRAG_LAW,KN)

! drag coefficient

INTEGER, INTENT(IN) :: DRAG_LAW
REAL(EB), INTENT(IN) :: RE
REAL(EB), OPTIONAL, INTENT(IN) :: KN

IF (RE<TWO_EPSILON_EB) THEN
   DRAG = 0._EB
   RETURN
ENDIF

SELECT CASE(DRAG_LAW)

   ! see J.P. Holman 7th Ed. Fig. 6-10
   CASE(SPHERE_DRAG)
      IF (RE<1._EB) THEN
         IF (PRESENT(KN)) THEN
            DRAG = 24._EB/RE/CUNNINGHAM(KN)
         ELSE
            DRAG = MIN(240000._EB,24._EB/RE)
         ENDIF
      ELSEIF (RE<1000._EB) THEN
         !!DRAG = 24._EB*(1._EB+0.15_EB*RE**0.687_EB)/RE ! see Crowe, Sommerfeld, Tsuji, 1998, Eq. (4.51)
         DRAG = 24._EB*(0.85_EB+0.15_EB*RE**0.687_EB)/RE ! matches Stokes drag at RE=1 (RJM)
      ELSEIF (RE>=1000._EB) THEN
         DRAG = 0.44_EB
      ENDIF

   ! see J.P. Holman 7th Ed. Fig. 6-9
   CASE(CYLINDER_DRAG)
      IF (RE<=1._EB) THEN
         DRAG = 10._EB/(RE**0.8_EB)
      ELSEIF (RE>1._EB .AND. RE<1000._EB) THEN
         DRAG = 10._EB*(0.6_EB+0.4_EB*RE**0.8_EB)/RE
      ELSEIF (RE>=1000._EB) THEN
         DRAG = 1._EB
      ENDIF

   ! see Hoerner 1965, Fig. 3-26
   CASE(DISK_DRAG)
      DRAG=20.37_EB/RE+1.17_EB/(1._EB+1._EB/RE)

   CASE(USER_DRAG)
      DRAG = 1._EB ! PC%DRAG_COEFFICIENT set elsewhere

   CASE DEFAULT
      DRAG = 0._EB

END SELECT

END FUNCTION DRAG


REAL(EB) FUNCTION CUNNINGHAM(KN)
REAL(EB), INTENT(IN) :: KN
REAL(EB), PARAMETER :: K1=1.257_EB,K2=0.4_EB,K3=1.1_EB

CUNNINGHAM = 1._EB+K1*KN+K2*KN*EXP(-K3/KN)

END FUNCTION CUNNINGHAM


!> \brief Compute the surface (mass or energy) density of a wall cell
!> \param MODE Indicator of returned quantity (0) kg/m2, (1) kg/m3, (2) kJ/m2, (3) kJ/m3
!> \param SF Pointer to SURFACE
!> \param ONE_D Pointer to BOUNDARY_ONE_D
!> \param MATL_INDEX (Optional) Index of the material component

REAL(EB) FUNCTION SURFACE_DENSITY(MODE,SF,ONE_D,MATL_INDEX)

INTEGER, INTENT(IN) :: MODE
INTEGER, INTENT(IN), OPTIONAL :: MATL_INDEX
INTEGER :: I_GRAD,NWP,II2,N,ITMP
REAL(EB) :: WGT,R_S(0:NWP_MAX),EPUM,DTMP
TYPE(BOUNDARY_ONE_D_TYPE), INTENT(IN), POINTER :: ONE_D
TYPE(SURFACE_TYPE), INTENT(IN), POINTER :: SF
TYPE(MATERIAL_TYPE), POINTER :: ML

THERMALLY_THICK_IF: IF (SF%THERMAL_BC_INDEX/=THERMALLY_THICK) THEN

   SURFACE_DENSITY = 0._EB

ELSE THERMALLY_THICK_IF

   SELECT CASE(SF%GEOMETRY)
      CASE(SURF_CARTESIAN)                           ; I_GRAD = 1
      CASE(SURF_CYLINDRICAL,SURF_INNER_CYLINDRICAL)  ; I_GRAD = 2
      CASE(SURF_SPHERICAL)                           ; I_GRAD = 3
   END SELECT

   NWP = SUM(ONE_D%N_LAYER_CELLS)
   IF (SF%GEOMETRY==SURF_INNER_CYLINDRICAL) THEN
      R_S(0:NWP) = SF%INNER_RADIUS + ONE_D%X(0:NWP)
   ELSE
      R_S(0:NWP) = SF%INNER_RADIUS + ONE_D%X(NWP) - ONE_D%X(0:NWP)
   ENDIF

   SURFACE_DENSITY = 0._EB
   NUMBER_WALL_POINTS_LOOP: DO II2=1,NWP
      AREA_VOLUME_SELECT: SELECT CASE(MODE)
         CASE(0,2); WGT = ABS(R_S(II2-1)**I_GRAD-R_S(II2)**I_GRAD)/(REAL(I_GRAD,EB)*(SF%INNER_RADIUS+SF%THICKNESS)**(I_GRAD-1))
         CASE(1,3); WGT = ABS(R_S(II2-1)**I_GRAD-R_S(II2)**I_GRAD)/(SF%INNER_RADIUS+SF%THICKNESS)**I_GRAD
      END SELECT AREA_VOLUME_SELECT

      EPUM = 1._EB ! energy per unit mass
      ITMP = MIN(I_MAX_TEMP-1,INT(ONE_D%TMP(II2)))
      DTMP = ONE_D%TMP(II2) - REAL(ITMP,EB)
      IF (PRESENT(MATL_INDEX)) THEN
         ENERGY_SELECT_1: SELECT CASE(MODE)
            CASE(2,3)
               ML => MATERIAL(MATL_INDEX)
               EPUM = ML%H(ITMP)+DTMP*(ML%H(ITMP+1)-ML%H(ITMP))
         END SELECT ENERGY_SELECT_1
         SURFACE_DENSITY = SURFACE_DENSITY + ONE_D%MATL_COMP(MATL_INDEX)%RHO(II2)*WGT*EPUM
      ELSE
         DO N=1,SF%N_MATL
            ENERGY_SELECT_2: SELECT CASE(MODE)
               CASE(2,3)
                  ML => MATERIAL(N)
                  EPUM = ML%H(ITMP)+DTMP*(ML%H(ITMP+1)-ML%H(ITMP))
            END SELECT ENERGY_SELECT_2
            SURFACE_DENSITY = SURFACE_DENSITY + ONE_D%MATL_COMP(N)%RHO(II2)*WGT*EPUM
         ENDDO
      ENDIF

   ENDDO NUMBER_WALL_POINTS_LOOP

ENDIF THERMALLY_THICK_IF

END FUNCTION SURFACE_DENSITY


!> \brief Compute particle Cumulative Number Fraction (CNF) and Cumulative Volume Fraction (CVF)
!> \param DM Median particle diameter (m)
!> \param RR Array of particle radii (m)
!> \param CNF Cumulative Number Fraction
!> \param CVF Cumulative Volume Fraction
!> \param NPT Number of points in the distribution
!> \param GAMMA Parameter in the distribution function
!> \param SIGMA Parameter in the distribution function
!> \param DISTRIBUTION Type of distribution

SUBROUTINE PARTICLE_SIZE_DISTRIBUTION(DM,RR,CNF,CVF,NPT,GAMMA,SIGMA,DISTRIBUTION)

USE MATH_FUNCTIONS, ONLY: IERFC
CHARACTER(LABEL_LENGTH), INTENT(IN) :: DISTRIBUTION
REAL(EB), INTENT(IN) :: DM,GAMMA,SIGMA
INTEGER, INTENT(IN) :: NPT
REAL(EB) :: SUM1,SUM2,DD1,DI,ETRM,GFAC,SFAC,DMIN,X1
INTEGER  :: J
REAL(EB), INTENT(OUT) :: RR(0:NPT),CNF(0:NPT),CVF(0:NPT)

RR(0)  = 0._EB
CNF(0) = 0._EB
SUM1   = 0._EB
SUM2   = 0._EB

X1     = IERFC(2._EB*CNF_CUTOFF)
DMIN   = MAX(DM*EXP(-X1*SQRT(2._EB)*SIGMA),0._EB)
DD1    = (-LOG(CNF_CUTOFF)/LOG(2._EB))**(1._EB/GAMMA)*DM
DD1    = (DD1-DMIN)/REAL(NPT,EB)
GFAC   = LOG(2._EB)*GAMMA*DD1/(DM**GAMMA)
SFAC   = DD1/(SQRT(TWOPI)*SIGMA)

INTLOOP: DO J=1,NPT
   DI = DMIN + (J-0.5_EB)*DD1
   RR(J) = 0.5_EB*DI
   IF ((DI<=DM .OR. DISTRIBUTION=='LOGNORMAL') .AND. DISTRIBUTION/='ROSIN-RAMMLER') THEN
      ETRM = EXP(-(LOG(DI/DM))**2/(2._EB*SIGMA**2))
      SUM1 = SUM1 + (SFAC/DI**4)*ETRM
      SUM2 = SUM2 + (SFAC/DI)*ETRM
   ELSE
      ETRM = EXP(-LOG(2._EB)*(DI/DM)**GAMMA)
      SUM1 = SUM1 + GFAC*DI**(GAMMA-4._EB)*ETRM
      SUM2 = 1._EB - ETRM
   ENDIF
   CNF(J) = SUM1
   CVF(J) = SUM2
ENDDO INTLOOP

CNF = CNF/SUM1
CVF = CVF/SUM2

END SUBROUTINE PARTICLE_SIZE_DISTRIBUTION


SUBROUTINE SPRAY_ANGLE_DISTRIBUTION(LON,LAT,LON_CDF,LAT_CDF,BETA,MU,SPRAY_ANGLE,DISTRIBUTION_TYPE,NPT)

INTEGER,INTENT(IN) :: NPT
REAL(EB),INTENT(OUT) :: LON_CDF(0:NPT),LON(0:NPT),LAT(0:NPT),LAT_CDF(0:NPT,0:NPT)
REAL(EB),INTENT(IN) :: BETA,MU,SPRAY_ANGLE(2,2)
CHARACTER(LABEL_LENGTH),INTENT(IN) :: DISTRIBUTION_TYPE
INTEGER :: I,J
REAL(EB) :: DLON,DLAT,PDF(0:NPT,0:NPT),THETA_MAX,THETA_MIN

THETA_MAX=MAXVAL(SPRAY_ANGLE)
THETA_MIN=MINVAL(SPRAY_ANGLE)

DLAT=(THETA_MAX-THETA_MIN)/REAL(NPT,EB)
DLON=2._EB*PI/REAL(NPT,EB)
!Discretize latitude and longtitude
LAT=(/ (THETA_MIN+I*DLAT,I=0,NPT) /)
LON=(/ (0._EB+I*DLON, I=0,NPT) /)

! SPRAY_ANGLE May be different in X and Y directions
! i.e spray angle is dependent on longtitude
! SPRAY_AGLE_MIN and MAX form ellipses with semi-axes defined by
! SPRAY_ANGLE(1,1:2) and SPRAY_ANGLE(2,1:2) respectively

DO I=0,NPT
   THETA_MIN=SPRAY_ANGLE(1,1)*SPRAY_ANGLE(1,2)
   IF (THETA_MIN>0._EB) THEN
      THETA_MIN=THETA_MIN/SQRT((SPRAY_ANGLE(1,2)*COS(LON(I)))**2+(SPRAY_ANGLE(1,1)*SIN(LON(I)))**2)
   ENDIF
   THETA_MAX=SPRAY_ANGLE(2,1)*SPRAY_ANGLE(2,2)
   THETA_MAX=THETA_MAX/SQRT((SPRAY_ANGLE(2,2)*COS(LON(I)))**2+(SPRAY_ANGLE(2,1)*SIN(LON(I)))**2)
   SELECT CASE(DISTRIBUTION_TYPE)
   CASE("TRIANGLE")
      DO J=0,NPT
         IF(LAT(J)<THETA_MIN .OR. LAT(J)>THETA_MAX) THEN
           PDF(J,I)=0._EB
         ELSE
           IF(LON(I)<MU) THEN
              PDF(J,I)=2*(LAT(J)-THETA_MIN)/((THETA_MAX-THETA_MIN)*(THETA_MAX-MU))
           ELSE
              PDF(J,I)=2*(THETA_MAX-LAT(J))/((THETA_MAX-THETA_MIN)*(THETA_MAX-MU))
           ENDIF
         ENDIF
      ENDDO
   CASE("GAUSSIAN")
      DO J=0,NPT
        IF(LAT(J)<THETA_MIN .OR. LAT(J)>THETA_MAX) THEN
           PDF(J,I)=0._EB
        ELSE
           PDF(J,I)=EXP(-BETA*((LAT(J)-MU)/(THETA_MAX-THETA_MIN))**2)
        ENDIF
      ENDDO
   CASE DEFAULT ! "UNIFORM"
       DO J=0,NPT
        IF(LAT(J)<THETA_MIN .OR. LAT(J)>THETA_MAX) THEN
           PDF(J,I)=0._EB
        ELSE
           PDF(J,I)=1._EB
        ENDIF
      ENDDO
   END SELECT
ENDDO


!
DO I=0,NPT
PDF(I,:)=PDF(I,:)*SIN(LAT(I))
ENDDO

LAT_CDF=0._EB
! Latitude distribution conditional on Longtitude
DO I=1,NPT
   LAT_CDF(I,:)=LAT_CDF(I-1,:)+0.5*(PDF(I,:)+PDF(I-1,:))*DLAT
ENDDO

! Marginal longtitude distribution
LON_CDF=0._EB
DO I=1,NPT
   LON_CDF(I)=LON_CDF(I-1)+0.5*(LAT_CDF(NPT,I-1)+LAT_CDF(NPT,I))*DLON
ENDDO

! Normalize marginal longtitude distribution
LON_CDF=LON_CDF/LON_CDF(NPT)
!Normalize conditional latitude distributions
DO I=1,NPT
   LAT_CDF(:,I)=LAT_CDF(:,I)/LAT_CDF(NPT,I)
ENDDO

END SUBROUTINE SPRAY_ANGLE_DISTRIBUTION


REAL(EB) FUNCTION FED(Y_IN,RSUM,FED_ACTIVITY)

! Returns the integrand of FED (Fractional Effective Dose) calculation.

REAL(EB), INTENT(IN) :: Y_IN(1:N_TRACKED_SPECIES),RSUM
INTEGER, INTENT(IN) :: FED_ACTIVITY
INTEGER  :: N
REAL(EB) :: Y_MF_INT, TMP_1
REAL(EB), DIMENSION(3) :: CO_FED_FAC
!                at rest           light work(default) heavy work
DATA CO_FED_FAC /0.70486250E-5_EB, 2.7641667E-5_EB,    8.2925E-5_EB/

! All equations from D.A. Purser, Sec. 2, Chap. 6, SFPE Handbook, 4th Ed.
! Note: Purser uses minutes, here dt is in seconds. Conversion at the end of the function.
! Total FED dose:
! FED_dose = (FED_LCO + FED_LCN + FED_LNOx + FLD_irr)*FED_VCO2 + FED_LO2;

FED = 0._EB

! Carbon monoxide (CO)
!          at rest    light work heavy work
! RMV_FED /8.5_EB,    25.0_EB,   50.0_EB /
! D_FED   /40.0_EB,   30.0_EB,   20.0_EB /
! RMV/D   /0.2125_EB, 0.8333_EB, 2.5_EB/
!
! FED_LCO = (3.317E-5 * (C_CO)^1.036 * RMV * (dt/60)) / D;
!   with RMV=25 [l/min], D=30 [%] COHb concentration at incapacitation and C_CO in ppm
!
IF (CO_INDEX > 0) THEN
   Call GET_MASS_FRACTION(Y_IN,CO_INDEX,Y_MF_INT)
   TMP_1 = SPECIES(CO_INDEX)%RCON*Y_MF_INT*1.E6_EB/RSUM
   ! FED   = 2.764E-5_EB*TMP_1**(1.036_EB)
   FED   = CO_FED_FAC(FED_ACTIVITY)*TMP_1**(1.036_EB)
ENDIF

! Nitrogen oxides (NOx, here NO + NO2)
! FED_LNOx = C_NOx/1500 * (dt/60);
!   with C_NOx = C_NO + C_NO2, all in ppm
TMP_1 = 0._EB
IF (NO_INDEX > 0) THEN
   Call GET_MASS_FRACTION(Y_IN,NO_INDEX,Y_MF_INT)
   TMP_1 = SPECIES(NO_INDEX)%RCON*Y_MF_INT/RSUM
ENDIF
IF (NO2_INDEX > 0) THEN
   Call GET_MASS_FRACTION(Y_IN,NO2_INDEX,Y_MF_INT)
   TMP_1 = TMP_1 + SPECIES(NO2_INDEX)%RCON*Y_MF_INT/RSUM
ENDIF
IF (TMP_1 > 0._EB) FED = FED + TMP_1/0.001500_EB

! Cyanide
! FED_LCN = (exp(C_CN/43)/220 - 0.0045) * (dt/60);
!   with C_CN = C_HCN - C_NOx, all in ppm
IF (HCN_INDEX > 0) THEN
   Call GET_MASS_FRACTION(Y_IN,HCN_INDEX,Y_MF_INT)
   TMP_1 = SPECIES(HCN_INDEX)%RCON*Y_MF_INT/RSUM - TMP_1
   IF (TMP_1 > 0._EB) FED = FED + (Exp(TMP_1/0.000043_EB)/220.0_EB-0.00454545_EB)
ENDIF

! Irritants
! FLD_irr = (C_HCl/F_HCl + C_HBr/F_HBr + C_HF/F_HF + C_SO2/F_SO2 + C_NO2/F_NO2 + C_C3H4O/F_C3H4O + C_CH2O/F_CH2O) * (dt/60);
!   all in ppm
TMP_1 = 0._EB
DO N=1,N_SPECIES
   IF (SPECIES(N)%FLD_LETHAL_DOSE > 0._EB) THEN
      Call GET_MASS_FRACTION(Y_IN,N,Y_MF_INT)
      TMP_1 = TMP_1 + SPECIES(N)%RCON*Y_MF_INT/RSUM / SPECIES(N)%FLD_LETHAL_DOSE
   ENDIF
ENDDO
FED = FED + TMP_1

! Carbon dioxide (CO2) induced hyperventilation:
! FED_VCO2 = exp(0.1903*C_CO2/1E4 + 2.0004)/7.1;
!   C_CO2 in ppm
IF (CO2_INDEX > 0) THEN
   Call GET_MASS_FRACTION(Y_IN,CO2_INDEX,Y_MF_INT)
   TMP_1 = SPECIES(CO2_INDEX)%RCON*Y_MF_INT/RSUM
   If ( TMP_1 > 0.0_EB ) FED = FED * Exp( 0.1903_EB*TMP_1*100.0_EB + 2.0004_EB )/7.1_EB
ENDIF

! Low oxygen (O2)
! FED_LO2 = 1/exp(8.13 - 0.54*(0.209 - C_O2/1E6)) * (dt/60);
!   C_O2 in ppm
IF (O2_INDEX > 0) THEN
   Call GET_MASS_FRACTION(Y_IN,O2_INDEX,Y_MF_INT)
   TMP_1 = SPECIES(O2_INDEX)%RCON*Y_MF_INT/RSUM
   IF ( TMP_1 < 0.20_EB ) FED = FED + 1.0_EB  / Exp(8.13_EB-0.54_EB*(20.9_EB-100.0_EB*TMP_1))
ENDIF

! Convert the FED integrand for minutes.
FED = FED / 60._EB

END FUNCTION FED


REAL(EB) FUNCTION FIC(Y_IN,RSUM)
! Returns FIC (Fractional Incapacitating Concentration)

REAL(EB), INTENT(IN) :: Y_IN(1:N_TRACKED_SPECIES),RSUM
REAL(EB) :: Y_MF_INT
INTEGER  :: N

FIC = 0._EB
DO N=1,N_SPECIES
   IF (SPECIES(N)%FIC_CONCENTRATION > 0._EB) THEN
      Call GET_MASS_FRACTION(Y_IN,N,Y_MF_INT)
      FIC = FIC + SPECIES(N)%RCON*Y_MF_INT/RSUM / SPECIES(N)%FIC_CONCENTRATION
   ENDIF
ENDDO

END FUNCTION FIC


REAL(EB) FUNCTION WATER_VAPOR_MASS_FRACTION(HUMIDITY,TEMP,PZONE)

! Compute the water vapor mass fraction given the relative humidity and temperature

REAL(EB), INTENT(IN) :: TEMP,HUMIDITY,PZONE
REAL(EB) :: X_SAT,DHOR_T,DHOR,P_RATIO,T_BOIL_EFF
REAL(EB),PARAMETER :: T_BOIL=373.15_EB,DHOR_T_B=4916.346083_EB

DHOR_T = (H_V_H2O(INT(TEMP))+(TEMP-REAL(INT(TEMP,EB)))*(H_V_H2O(INT(TEMP)+1)-H_V_H2O(INT(TEMP))))*MW_H2O/R0
P_RATIO = PZONE/P_STP
T_BOIL_EFF = MAX(0._EB,DHOR_T_B*T_BOIL/(DHOR_T_B-T_BOIL*LOG(P_RATIO)+TWO_EPSILON_EB))
DHOR = 0.5_EB*((H_V_H2O(INT(T_BOIL_EFF))+(T_BOIL_EFF-REAL(INT(T_BOIL_EFF,EB)))*&
       (H_V_H2O(INT(T_BOIL_EFF)+1)-H_V_H2O(INT(T_BOIL_EFF))))*MW_H2O/R0+DHOR_T)
X_SAT  = MIN(1._EB,EXP(DHOR*(1._EB/T_BOIL_EFF-1._EB/TEMP)))
WATER_VAPOR_MASS_FRACTION = HUMIDITY*0.01_EB*X_SAT/(MW_AIR/MW_H2O+(1._EB-MW_AIR/MW_H2O)*HUMIDITY*0.01_EB*X_SAT)

END FUNCTION WATER_VAPOR_MASS_FRACTION


REAL(EB) FUNCTION RELATIVE_HUMIDITY(Y_H2O,TEMP,PZONE)

! Compute the relative humidity given the water vapor mass fraction and temperature

REAL (EB), INTENT(IN) :: TEMP,Y_H2O,PZONE
REAL (EB) :: X_SAT,X_H2O,DHOR,DHOR_T,T_BOIL_EFF,P_RATIO
REAL(EB),PARAMETER :: T_BOIL=373.15_EB,DHOR_T_B=4916.346083_EB

P_RATIO = PZONE/P_STP
T_BOIL_EFF = MAX(0._EB,DHOR_T_B*T_BOIL/(DHOR_T_B-T_BOIL*LOG(P_RATIO)+TWO_EPSILON_EB))
IF (TEMP >= T_BOIL_EFF) THEN
   X_SAT = 1._EB
ELSE
   DHOR_T = (H_V_H2O(INT(TEMP))+(TEMP-REAL(INT(TEMP,EB)))*(H_V_H2O(INT(TEMP)+1)-H_V_H2O(INT(TEMP))))*MW_H2O/R0
   DHOR = 0.5_EB*((H_V_H2O(INT(T_BOIL_EFF))+(T_BOIL_EFF-REAL(INT(T_BOIL_EFF,EB)))*&
          (H_V_H2O(INT(T_BOIL_EFF)+1)-H_V_H2O(INT(T_BOIL_EFF))))*MW_H2O/R0+DHOR_T)
   X_SAT  = MIN(1._EB,EXP(DHOR*(1._EB/T_BOIL_EFF-1._EB/TEMP)))
ENDIF
X_H2O = Y_H2O*MW_AIR/(MW_H2O-Y_H2O*(MW_H2O-MW_AIR))
RELATIVE_HUMIDITY = 100._EB * X_H2O / X_SAT

END FUNCTION RELATIVE_HUMIDITY


REAL(EB) FUNCTION LES_FILTER_WIDTH_FUNCTION(DX,DY,DZ)
USE GLOBAL_CONSTANTS, ONLY : LES_FILTER_WIDTH_TYPE
REAL(EB), INTENT(IN):: DX,DY,DZ

SELECT CASE(LES_FILTER_WIDTH_TYPE)
   CASE(MEAN_LES_FILTER)
      IF (TWO_D) THEN
         LES_FILTER_WIDTH_FUNCTION = SQRT(DX*DZ)
      ELSE
         LES_FILTER_WIDTH_FUNCTION = (DX*DY*DZ)**ONTH
      ENDIF
   CASE(MAX_LES_FILTER)
      IF (TWO_D) THEN
         LES_FILTER_WIDTH_FUNCTION = MAX(DX,DZ)
      ELSE
         LES_FILTER_WIDTH_FUNCTION = MAX(DX,DY,DZ)
      ENDIF
   CASE(FIXED_LES_FILTER)
         LES_FILTER_WIDTH_FUNCTION = FIXED_LES_FILTER_WIDTH
END SELECT

END FUNCTION LES_FILTER_WIDTH_FUNCTION


REAL(EB) FUNCTION GET_POTENTIAL_TEMPERATURE(TMP_IN,Z_IN)

USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP
REAL(EB), INTENT(IN) :: TMP_IN,Z_IN
REAL(EB) :: PP

IF (STRATIFICATION) THEN
   PP = EVALUATE_RAMP(Z_IN,I_RAMP_P0_Z)
ELSE
   PP = P_INF
ENDIF
GET_POTENTIAL_TEMPERATURE = TMP_IN*(1.E5_EB/PP)**GM1OG  ! GM1OG = (GAMMA-1)/GAMMA = R/CP

END FUNCTION GET_POTENTIAL_TEMPERATURE


SUBROUTINE MONIN_OBUKHOV_SIMILARITY(Z,Z_0,L,U_STAR,THETA_STAR,THETA_0,U,TMP)

REAL(EB), INTENT(IN) :: Z,Z_0,L,U_STAR,THETA_STAR,THETA_0
REAL(EB), INTENT(OUT) :: U,TMP
REAL(EB), PARAMETER :: P_REF=100000._EB,RHO_REF=1.2_EB
REAL(EB) :: PSI_M,PSI_H,THETA,KAPPA

KAPPA = VON_KARMAN_CONSTANT
CALL MONIN_OBUKHOV_STABILITY_CORRECTIONS(PSI_M,PSI_H,Z,L)
U = (U_STAR/KAPPA)*(LOG(Z/Z_0)-PSI_M)
THETA = THETA_0 + (THETA_STAR/KAPPA)*(LOG(Z/Z_0)-PSI_H)
TMP = THETA*(P_REF/(P_REF-RHO_REF*GRAV*Z))**(-0.285_EB)

END SUBROUTINE MONIN_OBUKHOV_SIMILARITY


SUBROUTINE MONIN_OBUKHOV_STABILITY_CORRECTIONS(PSI_M,PSI_H,Z,L)

! Reference: A.J. Dyer. A review of flux profile relationships. Boundary-Layer Meteorology, 7:363-372, 1974.

REAL(EB), INTENT(IN) :: Z,L
REAL(EB), INTENT(OUT) :: PSI_M,PSI_H
REAL(EB) :: ZETA

IF (L>=0._EB) THEN
   ! stable boundary layer
   PSI_M = -5._EB*Z/L
   PSI_H = PSI_M
ELSE
   ! unstable boundary layer
   ZETA = (1._EB-16._EB*Z/L)**0.25_EB
   PSI_M = 2._EB*LOG(0.5_EB*(1._EB+ZETA)) + LOG(0.5_EB*(1._EB+ZETA**2)) - 2._EB*ATAN(ZETA) + 0.5_EB*PI
   PSI_H = 2._EB*LOG(0.5_EB*(1._EB+ZETA**2))
ENDIF

END SUBROUTINE MONIN_OBUKHOV_STABILITY_CORRECTIONS


!> \brief Computes the mass and heat transfer coeffiicents for a liquid droplet in FILM based on the selected EVAP_MODEL on MISC
!> \param H_MASS The mass transfer coefficient (m2/s)
!> \param H_HEAT The dropelt heat transfer coefficient (W/m2/K)
!> \param D_FILM Diffusivity in the film (m2/s)
!> \param K_FILM Conductivity in the film (W/m/k)
!> \param CP_FILM Specific heat in the film (J/kg/K)
!> \param RHO_FILM Density in in the film (kg/m3)
!> \param LENGTH Length scale (m)
!> \param Y_DROP Equilibrium vapor fraction for the current droplet temperature
!> \param Y_GAS Mass fraction of vapor in the gas
!> \param B_NUMBER B number for the droplet
!> \param NU_FAC_GAS Constant factor used in computing  the NUsselt number
!> \param SH_FAC_GAS Constant factor used in computing  the Sherwood number
!> \param RE_L Renyolds number
!> \param TMP_FILM Film temperature for the droplet (K)
!> \param ZZ_GET Tracked species mass fractions in the gas cell with the droplet
!> \param Z_INDEX Droplet species index in ZZ
!> \param EVAP_MODEL Indicator of evaporation model

SUBROUTINE DROPLET_H_MASS_H_HEAT_GAS(H_MASS,H_HEAT,D_FILM,K_FILM,CP_FILM,RHO_FILM,LENGTH,Y_DROP,Y_GAS,B_NUMBER,NU_FAC_GAS, &
                                     SH_FAC_GAS,RE_L,TMP_FILM,ZZ_GET,Z_INDEX,EVAP_MODEL)
USE MATH_FUNCTIONS, ONLY: F_B
REAL(EB), INTENT(IN) :: D_FILM,CP_FILM,K_FILM,RHO_FILM,LENGTH,Y_DROP,Y_GAS,NU_FAC_GAS,SH_FAC_GAS,RE_L,TMP_FILM, &
                        ZZ_GET(1:N_TRACKED_SPECIES)
INTEGER, INTENT(IN) :: Z_INDEX,EVAP_MODEL
REAL(EB), INTENT(INOUT) :: B_NUMBER
REAL(EB), INTENT(OUT) :: H_MASS,H_HEAT
REAL(EB) :: NUSSELT,SHERWOOD,LEWIS,THETA,C_GAS_DROP,C_GAS_FILM,ZZ_GET2(1:N_TRACKED_SPECIES)

SELECT CASE (EVAP_MODEL)
   CASE(RM_NO_B) ! Ranz Marshall
      NUSSELT  = 2._EB + NU_FAC_GAS*SQRT(RE_L)
      H_HEAT   = NUSSELT*K_FILM/LENGTH
      IF (Y_DROP <= Y_GAS) THEN
         H_MASS   = 0._EB
      ELSE
         SHERWOOD = 2._EB + SH_FAC_GAS*SQRT(RE_L)
         H_MASS   = SHERWOOD*D_FILM/LENGTH
      ENDIF
   CASE(RM_B) ! Sazhin M0, Eq 106 + 109 with B_T=B_M. This is the default model.
      IF (Y_DROP <= Y_GAS) THEN
         NUSSELT  = 2._EB + NU_FAC_GAS*SQRT(RE_L)
         H_HEAT   = NUSSELT*K_FILM/LENGTH
         H_MASS   = 0._EB
      ELSE
         NUSSELT  = ( 2._EB + NU_FAC_GAS*SQRT(RE_L) )*LOG(1._EB+B_NUMBER)/B_NUMBER
         H_HEAT   = NUSSELT*K_FILM/LENGTH
         SHERWOOD = ( 2._EB + SH_FAC_GAS*SQRT(RE_L) )*LOG(1._EB+B_NUMBER)/(Y_DROP-Y_GAS)
         H_MASS   = SHERWOOD*D_FILM/LENGTH
         ! above we save a divide and multiply of B_NUMBER
         ! the full model corresponding to Sazhin (108) and (109) would be
         ! SH = SH_0 * LOG(1+B_M)/B_M
         ! H_MASS = SH * D/L * B_M/(Y_D-Y_G)
      ENDIF
   CASE(RM_LEWIS_B) ! Sazhin M1, Eq 106 + 109 with eq 102.
      IF (Y_DROP <= Y_GAS) THEN
         NUSSELT  = 2._EB + NU_FAC_GAS*SQRT(RE_L)
         H_HEAT   = NUSSELT*K_FILM/LENGTH
         H_MASS   = 0._EB
      ELSE
         SHERWOOD = ( 2._EB + SH_FAC_GAS*SQRT(RE_L) )*LOG(1._EB+B_NUMBER)/(Y_DROP-Y_GAS)
         H_MASS   = SHERWOOD*D_FILM/LENGTH
         LEWIS    = K_FILM / (RHO_FILM * D_FILM * CP_FILM)
         ZZ_GET2(1:N_TRACKED_SPECIES) = 0._EB
         ZZ_GET2(Z_INDEX) = 1._EB
         CALL GET_SPECIFIC_HEAT(ZZ_GET2,C_GAS_DROP,TMP_FILM)
         CALL GET_SPECIFIC_HEAT(ZZ_GET,C_GAS_FILM,TMP_FILM)
         THETA = C_GAS_DROP/C_GAS_FILM/LEWIS
         B_NUMBER = (1._EB+B_NUMBER)**THETA-1._EB
         NUSSELT  = ( 2._EB + NU_FAC_GAS*SQRT(RE_L) )*LOG(1._EB+B_NUMBER)/B_NUMBER
         H_HEAT   = NUSSELT*K_FILM/LENGTH
      ENDIF
   CASE(RM_FL_LEWIS_B) ! Sazhin M2, Eq 116 and 117 with eq 106, 109, and 102.
      IF (Y_DROP <= Y_GAS) THEN
         NUSSELT  = 2._EB + NU_FAC_GAS*SQRT(RE_L)
         H_HEAT   = NUSSELT*K_FILM/LENGTH
         H_MASS   = 0._EB
      ELSE
         SHERWOOD = ( 2._EB + SH_FAC_GAS*SQRT(RE_L) )*LOG(1._EB+B_NUMBER)/((Y_DROP-Y_GAS)*F_B(B_NUMBER))
         H_MASS   = SHERWOOD*D_FILM/LENGTH
         LEWIS    = K_FILM / (RHO_FILM * D_FILM * CP_FILM)
         ZZ_GET2(1:N_TRACKED_SPECIES) = 0._EB
         ZZ_GET2(Z_INDEX) = 1._EB
         CALL GET_SPECIFIC_HEAT(ZZ_GET2,C_GAS_DROP,TMP_FILM)
         CALL GET_SPECIFIC_HEAT(ZZ_GET,C_GAS_FILM,TMP_FILM)
         THETA = C_GAS_DROP/C_GAS_FILM/LEWIS
         B_NUMBER = (1._EB+B_NUMBER)**THETA-1._EB
         NUSSELT  = ( 2._EB + NU_FAC_GAS*SQRT(RE_L) )*LOG(1._EB+B_NUMBER)/(B_NUMBER*F_B(B_NUMBER))
         H_HEAT   = NUSSELT*K_FILM/LENGTH
      ENDIF
END SELECT

END SUBROUTINE DROPLET_H_MASS_H_HEAT_GAS


!> \brief Compute the components of the prevailing wind
!> \param T Current time (s)
!> \param NM Current mesh

SUBROUTINE COMPUTE_WIND_COMPONENTS(T,NM)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
REAL(EB), INTENT(IN) :: T
REAL(EB) :: THETA, SIN_THETA, COS_THETA
INTEGER, INTENT(IN) :: NM
INTEGER :: K
TYPE (MESH_TYPE), POINTER :: M

M => MESHES(NM)

DO K=0,M%KBP1
   IF (I_RAMP_DIRECTION_T/=0 .OR. I_RAMP_DIRECTION_Z/=0) THEN
      IF (I_RAMP_DIRECTION_T==0) THEN
         THETA=EVALUATE_RAMP(M%ZC(K),I_RAMP_DIRECTION_Z)*DEG2RAD
      ELSEIF (I_RAMP_DIRECTION_Z==0) THEN
         THETA=EVALUATE_RAMP(T,I_RAMP_DIRECTION_T)*DEG2RAD
      ELSE
         THETA=(EVALUATE_RAMP(M%ZC(K),I_RAMP_DIRECTION_Z)+EVALUATE_RAMP(T,I_RAMP_DIRECTION_T))*DEG2RAD
      ENDIF
      SIN_THETA = -SIN(THETA)
      COS_THETA = -COS(THETA)
   ELSE
      SIN_THETA = 1._EB
      COS_THETA = 1._EB
   ENDIF
   M%U_WIND(K) = U0*EVALUATE_RAMP(M%ZC(K),I_RAMP_SPEED_Z)*EVALUATE_RAMP(T,I_RAMP_SPEED_T)*SIN_THETA
   M%V_WIND(K) = V0*EVALUATE_RAMP(M%ZC(K),I_RAMP_SPEED_Z)*EVALUATE_RAMP(T,I_RAMP_SPEED_T)*COS_THETA
   M%W_WIND(K) = W0
ENDDO

END SUBROUTINE COMPUTE_WIND_COMPONENTS


!> \brief Compute properties of the gas film for a liquid surface
!> \param N_VAP Number of evaporating fluids
!> \param FILM_FAC Linear factor for determining the filn conditions
!> \param Y_VAP Mass fraciton of the liquid vapors at the surface
!> \param Y_GAS Mass fraction of the liquid vapors in the gas cell
!> \param ZZ_INDEX Array of tracked species indicies for the evaporating liquids
!> \param TMP_S Temperature of the surface (K)
!> \param TMP_GAS Temperature of the gas cell (K)
!> \param ZZ_GAS Trakced species mass fractions in the gas cell
!> \param PB Film pressure (Pa)
!> \param TMP_FILM Film temperature (K)
!> \param MU_FILM Film viscosity (kg/m/s)
!> \param K_FILM Film conductivity (W/m/K)
!> \param CP_FILM Film specific heat (kJ/kg/K)
!> \param D_FILM Film diffusivity (m2/)
!> \param RHO_FILM Film density (kg/m3)
!> \param PR_FILM Film Prandtl number
!> \param SC_FILM Film Schmidt number

SUBROUTINE GET_FILM_PROPERTIES(N_VAP,FILM_FAC,Y_VAP,Y_GAS,ZZ_INDEX,TMP_S,TMP_GAS,ZZ_GAS,PB,TMP_FILM,MU_FILM,K_FILM,CP_FILM,D_FILM,&
                               RHO_FILM,PR_FILM,SC_FILM)

INTEGER, INTENT(IN) :: N_VAP,ZZ_INDEX(N_VAP)
REAL(EB), INTENT(IN) :: FILM_FAC,Y_VAP(N_VAP),Y_GAS(N_VAP),TMP_S,TMP_GAS,ZZ_GAS(1:N_TRACKED_SPECIES),PB
REAL(EB), INTENT(OUT) :: TMP_FILM,MU_FILM,K_FILM,CP_FILM,D_FILM,RHO_FILM,PR_FILM,SC_FILM
REAL(EB) :: X_SUM,R_FILM,ZZ_FILM(1:N_TRACKED_SPECIES),Y_FILM(N_VAP),SUM_FILM,SUM_GAS,OM_SUM_FILM
INTEGER :: I

! Take liquid surface Y and gas cell Y and compoute film Y
Y_FILM = Y_VAP + FILM_FAC*(Y_GAS-Y_VAP)
SUM_FILM = SUM(Y_FILM)
SUM_GAS = SUM(Y_GAS)
OM_SUM_FILM = 1._EB-SUM_FILM

IF (OM_SUM_FILM<TWO_EPSILON_EB) THEN
   ! If film is all vapor, just set the film Z to the vapor mass fractions.
   ZZ_FILM = 0._EB
   LOOP1: DO I=1,N_VAP
      IF (ZZ_INDEX(I)==0) CYCLE LOOP1
      ZZ_FILM(ZZ_INDEX(I)) = Y_FILM(I)
   ENDDO LOOP1
ELSE
   ! Determine the additional mass fraction of tracked species for each vapor species present
   IF (ABS(SUM_GAS-1._EB) < TWO_EPSILON_EB) THEN
      ZZ_FILM = 0._EB
      DO I=1,N_VAP
         ZZ_FILM(ZZ_INDEX(I)) = Y_FILM(I)
      ENDDO
      ZZ_FILM(1) = 1._EB - SUM_FILM
   ELSE
      ZZ_FILM = ZZ_GAS
      LOOP2: DO I=1,N_VAP
         IF (ZZ_INDEX(I)==0 .OR. Y_FILM(I)<TWO_EPSILON_EB) CYCLE LOOP2
         ZZ_FILM(ZZ_INDEX(I)) = ZZ_GAS(ZZ_INDEX(I)) + &
                                (Y_FILM(I)*(1._EB-SUM_GAS+Y_GAS(I))+Y_GAS(I)*(SUM_FILM-Y_FILM(I)-1._EB))/OM_SUM_FILM
      ENDDO LOOP2
      ZZ_FILM = ZZ_FILM/SUM(ZZ_FILM)
   ENDIF
ENDIF

! Use film mass fractions to get properties and non-dimensional parameters. D is weighed by mole fraction.

TMP_FILM = TMP_S + FILM_FAC*(TMP_GAS - TMP_S) ! LC Eq.(18)
CALL GET_VISCOSITY(ZZ_FILM,MU_FILM,TMP_FILM)
CALL GET_CONDUCTIVITY(ZZ_FILM,K_FILM,TMP_FILM)
CALL GET_SPECIFIC_HEAT(ZZ_FILM,CP_FILM,TMP_FILM)
CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_FILM,R_FILM)

X_SUM = 0._EB
LOOP3: DO I=1,N_VAP
   IF (ZZ_INDEX(I)==0 .OR. Y_FILM(I)<TWO_EPSILON_EB) CYCLE LOOP3
   D_FILM = D_Z(MIN(I_MAX_TEMP,NINT(TMP_FILM)),ZZ_INDEX(I))*ZZ_FILM(ZZ_INDEX(I))/SPECIES_MIXTURE(ZZ_INDEX(I))%MW
   X_SUM=X_SUM+ZZ_FILM(ZZ_INDEX(I))/SPECIES_MIXTURE(ZZ_INDEX(I))%MW
ENDDO LOOP3

IF (X_SUM > TWO_EPSILON_EB) THEN
   D_FILM = D_FILM/X_SUM
ELSE
   D_FILM = D_Z(NINT(TMPA),1)
ENDIF
PR_FILM = MU_FILM*CP_FILM/K_FILM
RHO_FILM = PB/(R_FILM*TMP_FILM)
PR_FILM = MU_FILM*CP_FILM/K_FILM
SC_FILM = MU_FILM/(RHO_FILM*D_FILM)

END SUBROUTINE GET_FILM_PROPERTIES

!> \brief Converts an array of liquid vapor mole fractions into mass fractions
!> \param N_MATS Number of liquids
!> \param ZZ_GET Gas cell tracked species mass factions
!> \param X_SV Liquid surface liquid vapor mole fraction
!> \param Y_SV Liquid surface liquid vapor mass fraction
!> \param MW Liquid vapor molecular weights
!> \param SMIX_INDEX Lookup for tracked species the liquids evaporate to

SUBROUTINE GET_Y_SURF(N_MATS,ZZ_GET,X_SV,Y_SV,MW,SMIX_INDEX)
INTEGER, INTENT(IN) :: N_MATS,SMIX_INDEX(N_MATS)
REAL(EB), INTENT(IN) :: ZZ_GET(1:N_TRACKED_SPECIES),X_SV(N_MATS),MW(N_MATS)
REAL(EB), INTENT(OUT) :: Y_SV(N_MATS)
REAL(EB) :: ZZ_2(1:N_TRACKED_SPECIES),MW_DRY,MASS
INTEGER :: I

ZZ_2 = ZZ_GET
Y_SV = 0._EB

!Zero out tracked species that are liquid vapors and get the MW of what is left

Y_ALL_LOOP: DO I=1,N_MATS
   IF (X_SV(I) < TWO_EPSILON_EB) CYCLE
   ZZ_2(SMIX_INDEX(I))=0._EB
ENDDO Y_ALL_LOOP

IF (SUM(ZZ_2) > TWO_EPSILON_EB) THEN
   ZZ_2 = ZZ_2 / SUM(ZZ_2)
   MW_DRY = 0._EB
   DO I=1,N_TRACKED_SPECIES
      MW_DRY = MW_DRY + ZZ_2(I)/SPECIES_MIXTURE(I)%MW
   ENDDO
   MW_DRY = 1._EB/MW_DRY
ENDIF

! Get mass based on mole fractions and MWs and convert X to Y
MASS = (1._EB-SUM(X_SV))*MW_DRY
DO I=1,N_MATS
   Y_SV(I) = X_SV(I)*MW(I)
   MASS = MASS + Y_SV(I)
ENDDO

Y_SV = Y_SV/MASS

END SUBROUTINE GET_Y_SURF


!> \brief Estimates the peak reaction temperature for a material reaction
!> \param N_MATL MATL index
!> \param NR Reaction index

SUBROUTINE GET_TMP_REF(N_MATL,NR)
INTEGER, INTENT(IN) :: N_MATL,NR
REAL(EB) :: HEATING_RATE,DT=0.01_EB,DTDT,RR_MAX,REACTION_RATE,TMP,RHO_S
TYPE(MATERIAL_TYPE), POINTER :: ML=>NULL()

ML=> MATERIAL(N_MATL)

IF (ML%RATE_REF(NR) > 0._EB) THEN
   HEATING_RATE = ML%RATE_REF(NR)
ELSE
   HEATING_RATE = TGA_HEATING_RATE
ENDIF

TMP = 0._EB
ML%TMP_REF(NR) = 0._EB
DTDT = HEATING_RATE/60._EB
RHO_S = ML%RHO_S
IF (ABS(ML%E(NR)) < TWO_EPSILON_EB) THEN
   RR_MAX = ML%A(NR)*RHO_S**ML%N_S(NR)
ELSE
   RR_MAX = 0._EB
ENDIF

DO WHILE (INT(TMP)<I_MAX_TEMP)
   TMP = TMP + DTDT * DT
   REACTION_RATE = ML%A(NR)*RHO_S**ML%N_S(NR)*EXP(-ML%E(NR)/(R0*TMP))*TMP**ML%N_T(NR)*ML%X_O2_PYRO**ML%N_O2(NR)
   IF (REACTION_RATE > RR_MAX) THEN
      ML%TMP_REF(NR) = TMP
      RR_MAX = REACTION_RATE
   ENDIF
   RHO_S = RHO_S - REACTION_RATE * DT
   IF (RHO_S<TWO_EPSILON_EB) EXIT
ENDDO

END SUBROUTINE GET_TMP_REF


!> \brief Determine if an ember ignites a substrate
!> \param NM Mesh number
!> \param IP Particle index

SUBROUTINE EMBER_IGNITION_MODEL(NM,IP)

INTEGER, INTENT(IN) :: NM,IP
TYPE(MESH_TYPE), POINTER :: M
TYPE(LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE(SURFACE_TYPE), POINTER :: SF
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
REAL(EB) :: F_N,POWER,X

M  => MESHES(NM)
LP => M%LAGRANGIAN_PARTICLE(IP)
BC => M%BOUNDARY_COORD(LP%BC_INDEX)
SF => SURFACE(M%LS_SURF_INDEX(BC%IIG,BC%JJG))
IF (SF%EMBER_IGNITION_POWER_MEAN<0._EB) RETURN

! Compute the energy generation rate of the ember

POWER = (M%BOUNDARY_PROP1(LP%B1_INDEX)%Q_RAD_OUT-M%BOUNDARY_PROP1(LP%B1_INDEX)%Q_RAD_IN - M%BOUNDARY_PROP1(LP%B1_INDEX)%Q_CON_F)*&
         0.001_EB*M%BOUNDARY_PROP1(LP%B1_INDEX)%AREA

! Compute the CDF of the POWER
F_N = 0.5_EB*(1._EB + ERF((POWER-SF%EMBER_IGNITION_POWER_MEAN)/(SR2*SF%EMBER_IGNITION_POWER_SIGMA)))
! Adjust for weighting factor
F_N = 1._EB - (1._EB - F_N)**LP%PWT

CALL RANDOM_NUMBER(X)
IF (X<F_N) THEN
   M%PHI_LS(BC%IIG,BC%JJG) = 1._EB
   BC%Z = -1.E6_EB
ENDIF

END SUBROUTINE EMBER_IGNITION_MODEL

!> \brief Computes refrence heat flux for S-pyro model using empirical data
!> \param HRRPUA Current wall cell HRRPUA (W/m2)
!> \param HOC Effective heat of combustion for all gases emitted by wall cell (J/kg)
!> \param Y_S Effective soot yield for all gases emitted by wall cell (kg/kg)
!> \param CHI_R Effective radiative fraction for all gases emitted by wall cell (kg/kg)
!> \param Q_INC current incident radiation (W/m2)

REAL(EB) FUNCTION Q_REF_FIT(HRRPUA,HOC,Y_S,CHI_R,Q_INC)
USE MATH_FUNCTIONS, ONLY: INTERPOLATE1D_UNIFORM
REAL(EB), INTENT(IN) :: HRRPUA,HOC,Y_S,Q_INC,CHI_R
INTEGER :: I_HOC,I_Y_S,I_HRRPUA,I_CHI_R
REAL(EB) :: F_HOC,F_Y_S,F_HRRPUA,F_CHI_R
REAL(EB) :: QFLA_1(2,2,2),ABSA_1(2,2,2)
REAL(EB) :: QFLA_2(2,2),ABSA_2(2,2)
REAL(EB) :: QFLA_3(2),ABSA_3(2)
REAL(EB) :: QFLA_4,ABSA_4

IF (.NOT. DEF_QREF_ARRAYS) THEN
   CALL DEFINE_QREF_ARRAYS
   DEF_QREF_ARRAYS = .TRUE.
ENDIF

I_HRRPUA = MIN(29,INT(HRRPUA*0.00001_EB))
F_HRRPUA = MIN(30._EB,HRRPUA*0.00001_EB) - REAL(I_HRRPUA,EB)
IF (F_HRRPUA < TWO_EPSILON_EB) THEN
   Q_REF_FIT = Q_INC
   RETURN
ENDIF
I_HOC = MIN(4,MAX(1,INT(HOC*1.E-7_EB)))
F_HOC = MIN(5._EB,MAX(1._EB,HOC*1.E-7_EB)) - REAL(I_HOC,EB)
I_Y_S = MIN(19,INT(Y_S*100._EB))
F_Y_S = MIN(20._EB,Y_S*100._EB) - REAL(I_Y_S,EB)
I_CHI_R = MIN(5,INT(CHI_R*10_EB))
F_CHI_R = MIN(6._EB,CHI_R*10_EB) - REAL(I_CHI_R,EB)

QFLA_1 = 0._EB
ABSA_1 = 0._EB

QFLA_2 = 0._EB
ABSA_2 = 0._EB
QFLA_1(1:2,1:2,1:2) = QFLAME(I_HOC:I_HOC+1,I_Y_S:I_Y_S+1,I_HRRPUA,I_CHI_R:I_CHI_R+1) +  F_HRRPUA * &
  (QFLAME(I_HOC:I_HOC+1,I_Y_S:I_Y_S+1,I_HRRPUA+1,I_CHI_R:I_CHI_R+1)-QFLAME(I_HOC:I_HOC+1,I_Y_S:I_Y_S+1,I_HRRPUA,I_CHI_R:I_CHI_R+1))
ABSA_1(1:2,1:2,1:2) = ABSF(I_HOC:I_HOC+1,I_Y_S:I_Y_S+1,I_HRRPUA,I_CHI_R:I_CHI_R+1) +  F_HRRPUA * &
  (ABSF(I_HOC:I_HOC+1,I_Y_S:I_Y_S+1,I_HRRPUA+1,I_CHI_R:I_CHI_R+1)-ABSF(I_HOC:I_HOC+1,I_Y_S:I_Y_S+1,I_HRRPUA,I_CHI_R:I_CHI_R+1))

QFLA_2(1:2,1:2) = QFLA_1(1:2,1:2,1)+ F_CHI_R * (QFLA_1(1:2,1:2,2) - QFLA_1(1:2,1:2,1))
ABSA_2(1:2,1:2) = ABSA_1(1:2,1:2,1)+ F_CHI_R * (ABSA_1(1:2,1:2,2) - ABSA_1(1:2,1:2,1))

QFLA_3(1:2) = QFLA_2(1:2,1)+ F_Y_S * (QFLA_2(1:2,2) - QFLA_2(1:2,1))
ABSA_3(1:2) = ABSA_2(1:2,1)+ F_Y_S * (ABSA_2(1:2,2) - ABSA_2(1:2,1))

QFLA_4 = QFLA_3(1)+ F_HOC * (QFLA_3(2) - QFLA_3(1))
ABSA_4 = ABSA_3(1)+ F_HOC * (ABSA_3(2) - ABSA_3(1))

Q_REF_FIT = QFLA_4*1000._EB + Q_INC * (1._EB - ABSA_4)

END FUNCTION Q_REF_FIT


SUBROUTINE DEFINE_QREF_ARRAYS

QFLAME(1,0,0:10,1) = (/0.00_EB,15.24_EB,17.38_EB,16.77_EB,16.17_EB,15.60_EB,15.04_EB,14.47_EB,13.90_EB,13.18_EB,12.45_EB/)
QFLAME(1,0,11:20,1) = (/11.73_EB,11.00_EB,10.42_EB,9.85_EB,9.27_EB,8.69_EB,8.31_EB,7.92_EB,7.54_EB,7.15_EB/)
QFLAME(1,0,21:30,1) = (/6.75_EB,6.36_EB,5.96_EB,5.56_EB,5.16_EB,5.07_EB,4.97_EB,4.88_EB,4.79_EB,4.69_EB/)
QFLAME(1,0,0:10,2) = (/0.00_EB,17.10_EB,19.83_EB,19.71_EB,19.60_EB,19.13_EB,18.66_EB,18.20_EB,17.73_EB,17.28_EB,16.83_EB/)
QFLAME(1,0,11:20,2) = (/16.38_EB,15.93_EB,15.37_EB,14.80_EB,14.23_EB,13.66_EB,13.20_EB,12.74_EB,12.28_EB,11.82_EB/)
QFLAME(1,0,21:30,2) = (/11.33_EB,10.83_EB,10.34_EB,9.85_EB,9.35_EB,9.18_EB,9.00_EB,8.82_EB,8.65_EB,8.47_EB/)
QFLAME(1,0,0:10,3) = (/0.00_EB,19.11_EB,23.31_EB,24.07_EB,24.83_EB,24.61_EB,24.39_EB,24.16_EB,23.94_EB,23.76_EB,23.57_EB/)
QFLAME(1,0,11:20,3) = (/23.39_EB,23.21_EB,22.55_EB,21.89_EB,21.23_EB,20.57_EB,20.01_EB,19.45_EB,18.88_EB,18.32_EB/)
QFLAME(1,0,21:30,3) = (/17.64_EB,16.95_EB,16.26_EB,15.58_EB,14.89_EB,14.59_EB,14.28_EB,13.98_EB,13.68_EB,13.37_EB/)
QFLAME(1,0,0:10,4) = (/0.00_EB,21.17_EB,26.97_EB,28.80_EB,30.64_EB,30.78_EB,30.91_EB,31.05_EB,31.18_EB,31.13_EB,31.08_EB/)
QFLAME(1,0,11:20,4) = (/31.03_EB,30.98_EB,30.44_EB,29.89_EB,29.35_EB,28.80_EB,28.11_EB,27.42_EB,26.72_EB,26.03_EB/)
QFLAME(1,0,21:30,4) = (/25.19_EB,24.35_EB,23.52_EB,22.68_EB,21.84_EB,21.25_EB,20.67_EB,20.08_EB,19.49_EB,18.90_EB/)
QFLAME(1,0,0:10,5) = (/0.00_EB,23.29_EB,30.85_EB,33.88_EB,36.91_EB,37.52_EB,38.12_EB,38.73_EB,39.33_EB,39.64_EB,39.95_EB/)
QFLAME(1,0,11:20,5) = (/40.26_EB,40.57_EB,39.95_EB,39.33_EB,38.71_EB,38.09_EB,36.96_EB,35.84_EB,34.71_EB,33.58_EB/)
QFLAME(1,0,21:30,5) = (/32.82_EB,32.07_EB,31.31_EB,30.55_EB,29.80_EB,29.01_EB,28.23_EB,27.44_EB,26.66_EB,25.87_EB/)
QFLAME(1,0,0:10,6) = (/0.00_EB,25.27_EB,34.75_EB,38.99_EB,43.23_EB,44.43_EB,45.64_EB,46.84_EB,48.05_EB,48.62_EB,49.19_EB/)
QFLAME(1,0,11:20,6) = (/49.76_EB,50.34_EB,49.86_EB,49.38_EB,48.91_EB,48.43_EB,47.37_EB,46.31_EB,45.25_EB,44.20_EB/)
QFLAME(1,0,21:30,6) = (/43.23_EB,42.27_EB,41.31_EB,40.35_EB,39.38_EB,38.21_EB,37.03_EB,35.86_EB,34.68_EB,33.50_EB/)
QFLAME(1,1,0:10,1) = (/0.00_EB,16.02_EB,18.83_EB,18.61_EB,18.38_EB,17.90_EB,17.42_EB,16.94_EB,16.45_EB,15.71_EB,14.97_EB/)
QFLAME(1,1,11:20,1) = (/14.22_EB,13.48_EB,12.85_EB,12.22_EB,11.59_EB,10.96_EB,10.54_EB,10.11_EB,9.69_EB,9.26_EB/)
QFLAME(1,1,21:30,1) = (/8.74_EB,8.21_EB,7.69_EB,7.17_EB,6.64_EB,6.51_EB,6.37_EB,6.24_EB,6.10_EB,5.96_EB/)
QFLAME(1,1,0:10,2) = (/0.00_EB,17.10_EB,19.79_EB,19.58_EB,19.37_EB,18.89_EB,18.42_EB,17.94_EB,17.47_EB,17.02_EB,16.57_EB/)
QFLAME(1,1,11:20,2) = (/16.13_EB,15.68_EB,15.14_EB,14.60_EB,14.07_EB,13.53_EB,13.05_EB,12.57_EB,12.09_EB,11.61_EB/)
QFLAME(1,1,21:30,2) = (/11.11_EB,10.61_EB,10.11_EB,9.62_EB,9.12_EB,8.96_EB,8.80_EB,8.65_EB,8.49_EB,8.33_EB/)
QFLAME(1,1,0:10,3) = (/0.00_EB,19.13_EB,23.27_EB,23.95_EB,24.63_EB,24.37_EB,24.12_EB,23.86_EB,23.60_EB,23.40_EB,23.19_EB/)
QFLAME(1,1,11:20,3) = (/22.99_EB,22.79_EB,22.21_EB,21.63_EB,21.06_EB,20.48_EB,19.91_EB,19.34_EB,18.78_EB,18.21_EB/)
QFLAME(1,1,21:30,3) = (/17.51_EB,16.80_EB,16.10_EB,15.40_EB,14.69_EB,14.43_EB,14.16_EB,13.90_EB,13.63_EB,13.37_EB/)
QFLAME(1,1,0:10,4) = (/0.00_EB,21.17_EB,26.96_EB,28.69_EB,30.43_EB,30.51_EB,30.58_EB,30.66_EB,30.74_EB,30.71_EB,30.69_EB/)
QFLAME(1,1,11:20,4) = (/30.67_EB,30.65_EB,29.99_EB,29.33_EB,28.67_EB,28.01_EB,27.42_EB,26.82_EB,26.23_EB,25.63_EB/)
QFLAME(1,1,21:30,4) = (/24.84_EB,24.05_EB,23.26_EB,22.47_EB,21.67_EB,21.01_EB,20.35_EB,19.68_EB,19.02_EB,18.36_EB/)
QFLAME(1,1,0:10,5) = (/0.00_EB,23.30_EB,30.88_EB,33.78_EB,36.68_EB,37.20_EB,37.72_EB,38.25_EB,38.77_EB,39.12_EB,39.47_EB/)
QFLAME(1,1,11:20,5) = (/39.82_EB,40.17_EB,39.58_EB,38.99_EB,38.40_EB,37.81_EB,36.76_EB,35.71_EB,34.66_EB,33.61_EB/)
QFLAME(1,1,21:30,5) = (/32.76_EB,31.90_EB,31.05_EB,30.19_EB,29.34_EB,28.43_EB,27.53_EB,26.62_EB,25.71_EB,24.81_EB/)
QFLAME(1,1,0:10,6) = (/0.00_EB,25.29_EB,34.72_EB,38.83_EB,42.94_EB,44.07_EB,45.19_EB,46.32_EB,47.45_EB,47.79_EB,48.12_EB/)
QFLAME(1,1,11:20,6) = (/48.46_EB,48.80_EB,48.37_EB,47.93_EB,47.50_EB,47.06_EB,46.22_EB,45.38_EB,44.54_EB,43.70_EB/)
QFLAME(1,1,21:30,6) = (/42.58_EB,41.45_EB,40.32_EB,39.19_EB,38.07_EB,36.97_EB,35.86_EB,34.76_EB,33.66_EB,32.56_EB/)
QFLAME(1,2,0:10,1) = (/0.00_EB,16.66_EB,19.96_EB,20.01_EB,20.06_EB,19.64_EB,19.22_EB,18.80_EB,18.37_EB,17.65_EB,16.92_EB/)
QFLAME(1,2,11:20,1) = (/16.20_EB,15.47_EB,14.81_EB,14.15_EB,13.48_EB,12.82_EB,12.36_EB,11.90_EB,11.44_EB,10.97_EB/)
QFLAME(1,2,21:30,1) = (/10.37_EB,9.77_EB,9.17_EB,8.57_EB,7.97_EB,7.78_EB,7.60_EB,7.41_EB,7.22_EB,7.03_EB/)
QFLAME(1,2,0:10,2) = (/0.00_EB,17.09_EB,20.07_EB,20.08_EB,20.09_EB,19.68_EB,19.26_EB,18.85_EB,18.43_EB,17.77_EB,17.10_EB/)
QFLAME(1,2,11:20,2) = (/16.43_EB,15.76_EB,15.16_EB,14.56_EB,13.96_EB,13.36_EB,12.93_EB,12.50_EB,12.08_EB,11.65_EB/)
QFLAME(1,2,21:30,2) = (/11.13_EB,10.62_EB,10.10_EB,9.58_EB,9.06_EB,8.94_EB,8.82_EB,8.69_EB,8.57_EB,8.44_EB/)
QFLAME(1,2,0:10,3) = (/0.00_EB,19.11_EB,23.24_EB,23.87_EB,24.51_EB,24.22_EB,23.94_EB,23.65_EB,23.37_EB,23.14_EB,22.91_EB/)
QFLAME(1,2,11:20,3) = (/22.69_EB,22.46_EB,21.92_EB,21.39_EB,20.85_EB,20.31_EB,19.76_EB,19.22_EB,18.67_EB,18.12_EB/)
QFLAME(1,2,21:30,3) = (/17.39_EB,16.65_EB,15.92_EB,15.18_EB,14.45_EB,14.18_EB,13.92_EB,13.65_EB,13.38_EB,13.12_EB/)
QFLAME(1,2,0:10,4) = (/0.00_EB,21.17_EB,26.92_EB,28.60_EB,30.28_EB,30.35_EB,30.42_EB,30.48_EB,30.55_EB,30.55_EB,30.54_EB/)
QFLAME(1,2,11:20,4) = (/30.53_EB,30.53_EB,29.88_EB,29.24_EB,28.60_EB,27.96_EB,27.31_EB,26.66_EB,26.00_EB,25.35_EB/)
QFLAME(1,2,21:30,4) = (/24.61_EB,23.86_EB,23.12_EB,22.37_EB,21.63_EB,20.84_EB,20.04_EB,19.25_EB,18.46_EB,17.67_EB/)
QFLAME(1,2,0:10,5) = (/0.00_EB,23.48_EB,30.85_EB,33.67_EB,36.48_EB,36.94_EB,37.39_EB,37.84_EB,38.30_EB,38.38_EB,38.46_EB/)
QFLAME(1,2,11:20,5) = (/38.54_EB,38.63_EB,38.21_EB,37.80_EB,37.39_EB,36.98_EB,36.00_EB,35.01_EB,34.03_EB,33.05_EB/)
QFLAME(1,2,21:30,5) = (/32.17_EB,31.29_EB,30.41_EB,29.53_EB,28.65_EB,27.77_EB,26.89_EB,26.01_EB,25.13_EB,24.26_EB/)
QFLAME(1,2,0:10,6) = (/0.00_EB,25.34_EB,34.69_EB,38.66_EB,42.63_EB,43.67_EB,44.70_EB,45.74_EB,46.77_EB,47.22_EB,47.66_EB/)
QFLAME(1,2,11:20,6) = (/48.11_EB,48.55_EB,48.11_EB,47.66_EB,47.22_EB,46.77_EB,45.78_EB,44.79_EB,43.80_EB,42.81_EB/)
QFLAME(1,2,21:30,6) = (/41.78_EB,40.75_EB,39.72_EB,38.70_EB,37.67_EB,36.61_EB,35.54_EB,34.48_EB,33.42_EB,32.35_EB/)
QFLAME(1,3,0:10,1) = (/0.00_EB,17.10_EB,20.75_EB,20.98_EB,21.22_EB,20.86_EB,20.50_EB,20.14_EB,19.78_EB,19.05_EB,18.32_EB/)
QFLAME(1,3,11:20,1) = (/17.58_EB,16.85_EB,16.17_EB,15.50_EB,14.82_EB,14.14_EB,13.64_EB,13.15_EB,12.65_EB,12.16_EB/)
QFLAME(1,3,21:30,1) = (/11.53_EB,10.91_EB,10.29_EB,9.67_EB,9.05_EB,8.80_EB,8.55_EB,8.30_EB,8.06_EB,7.81_EB/)
QFLAME(1,3,0:10,2) = (/0.00_EB,17.38_EB,20.83_EB,21.03_EB,21.24_EB,20.88_EB,20.53_EB,20.17_EB,19.82_EB,19.12_EB,18.43_EB/)
QFLAME(1,3,11:20,2) = (/17.74_EB,17.05_EB,16.41_EB,15.78_EB,15.14_EB,14.50_EB,14.03_EB,13.55_EB,13.08_EB,12.61_EB/)
QFLAME(1,3,21:30,2) = (/12.04_EB,11.47_EB,10.91_EB,10.34_EB,9.77_EB,9.57_EB,9.37_EB,9.17_EB,8.97_EB,8.77_EB/)
QFLAME(1,3,0:10,3) = (/0.00_EB,19.10_EB,23.20_EB,23.80_EB,24.39_EB,24.14_EB,23.88_EB,23.63_EB,23.37_EB,23.12_EB,22.86_EB/)
QFLAME(1,3,11:20,3) = (/22.61_EB,22.35_EB,21.79_EB,21.23_EB,20.67_EB,20.11_EB,19.60_EB,19.09_EB,18.57_EB,18.06_EB/)
QFLAME(1,3,21:30,3) = (/17.32_EB,16.59_EB,15.86_EB,15.12_EB,14.39_EB,14.11_EB,13.84_EB,13.57_EB,13.29_EB,13.02_EB/)
QFLAME(1,3,0:10,4) = (/0.00_EB,21.19_EB,26.88_EB,28.51_EB,30.14_EB,30.20_EB,30.25_EB,30.30_EB,30.36_EB,30.35_EB,30.34_EB/)
QFLAME(1,3,11:20,4) = (/30.33_EB,30.32_EB,29.68_EB,29.03_EB,28.39_EB,27.74_EB,27.12_EB,26.50_EB,25.88_EB,25.26_EB/)
QFLAME(1,3,21:30,4) = (/24.49_EB,23.73_EB,22.96_EB,22.20_EB,21.43_EB,20.68_EB,19.94_EB,19.19_EB,18.44_EB,17.69_EB/)
QFLAME(1,3,0:10,5) = (/0.00_EB,23.48_EB,30.81_EB,33.57_EB,36.33_EB,36.76_EB,37.19_EB,37.62_EB,38.05_EB,38.21_EB,38.37_EB/)
QFLAME(1,3,11:20,5) = (/38.53_EB,38.68_EB,38.20_EB,37.71_EB,37.22_EB,36.74_EB,35.74_EB,34.75_EB,33.76_EB,32.77_EB/)
QFLAME(1,3,21:30,5) = (/31.90_EB,31.03_EB,30.16_EB,29.29_EB,28.42_EB,27.55_EB,26.67_EB,25.79_EB,24.92_EB,24.04_EB/)
QFLAME(1,3,0:10,6) = (/0.00_EB,25.41_EB,34.66_EB,38.55_EB,42.45_EB,43.46_EB,44.47_EB,45.48_EB,46.50_EB,46.92_EB,47.34_EB/)
QFLAME(1,3,11:20,6) = (/47.76_EB,48.18_EB,47.78_EB,47.39_EB,47.00_EB,46.61_EB,45.64_EB,44.66_EB,43.68_EB,42.70_EB/)
QFLAME(1,3,21:30,6) = (/41.62_EB,40.54_EB,39.45_EB,38.37_EB,37.29_EB,36.24_EB,35.18_EB,34.13_EB,33.08_EB,32.03_EB/)
QFLAME(1,4,0:10,1) = (/0.00_EB,17.54_EB,21.54_EB,21.96_EB,22.37_EB,22.07_EB,21.78_EB,21.48_EB,21.19_EB,20.45_EB,19.71_EB/)
QFLAME(1,4,11:20,1) = (/18.97_EB,18.24_EB,17.54_EB,16.84_EB,16.15_EB,15.45_EB,14.92_EB,14.39_EB,13.87_EB,13.34_EB/)
QFLAME(1,4,21:30,1) = (/12.69_EB,12.05_EB,11.41_EB,10.77_EB,10.12_EB,9.82_EB,9.51_EB,9.20_EB,8.89_EB,8.58_EB/)
QFLAME(1,4,0:10,2) = (/0.00_EB,17.68_EB,21.58_EB,21.98_EB,22.38_EB,22.09_EB,21.79_EB,21.49_EB,21.20_EB,20.48_EB,19.77_EB/)
QFLAME(1,4,11:20,2) = (/19.06_EB,18.34_EB,17.67_EB,16.99_EB,16.32_EB,15.64_EB,15.12_EB,14.60_EB,14.08_EB,13.57_EB/)
QFLAME(1,4,21:30,2) = (/12.95_EB,12.33_EB,11.72_EB,11.10_EB,10.48_EB,10.21_EB,9.93_EB,9.66_EB,9.38_EB,9.10_EB/)
QFLAME(1,4,0:10,3) = (/0.00_EB,19.08_EB,23.16_EB,23.72_EB,24.27_EB,24.05_EB,23.82_EB,23.60_EB,23.38_EB,23.09_EB,22.81_EB/)
QFLAME(1,4,11:20,3) = (/22.53_EB,22.25_EB,21.66_EB,21.08_EB,20.50_EB,19.92_EB,19.43_EB,18.95_EB,18.47_EB,17.99_EB/)
QFLAME(1,4,21:30,3) = (/17.26_EB,16.53_EB,15.79_EB,15.06_EB,14.33_EB,14.05_EB,13.77_EB,13.49_EB,13.21_EB,12.93_EB/)
QFLAME(1,4,0:10,4) = (/0.00_EB,21.20_EB,26.83_EB,28.42_EB,30.00_EB,30.04_EB,30.08_EB,30.12_EB,30.16_EB,30.15_EB,30.14_EB/)
QFLAME(1,4,11:20,4) = (/30.13_EB,30.12_EB,29.47_EB,28.82_EB,28.17_EB,27.52_EB,26.93_EB,26.34_EB,25.75_EB,25.16_EB/)
QFLAME(1,4,21:30,4) = (/24.38_EB,23.59_EB,22.81_EB,22.02_EB,21.24_EB,20.53_EB,19.83_EB,19.12_EB,18.42_EB,17.71_EB/)
QFLAME(1,4,0:10,5) = (/0.00_EB,23.47_EB,30.77_EB,33.48_EB,36.18_EB,36.59_EB,37.00_EB,37.40_EB,37.81_EB,38.04_EB,38.27_EB/)
QFLAME(1,4,11:20,5) = (/38.51_EB,38.74_EB,38.18_EB,37.62_EB,37.05_EB,36.49_EB,35.49_EB,34.49_EB,33.49_EB,32.49_EB/)
QFLAME(1,4,21:30,5) = (/31.63_EB,30.77_EB,29.91_EB,29.05_EB,28.20_EB,27.32_EB,26.45_EB,25.57_EB,24.70_EB,23.82_EB/)
QFLAME(1,4,0:10,6) = (/0.00_EB,25.49_EB,34.63_EB,38.45_EB,42.26_EB,43.25_EB,44.24_EB,45.23_EB,46.22_EB,46.62_EB,47.01_EB/)
QFLAME(1,4,11:20,6) = (/47.40_EB,47.80_EB,47.46_EB,47.13_EB,46.79_EB,46.46_EB,45.49_EB,44.53_EB,43.56_EB,42.60_EB/)
QFLAME(1,4,21:30,6) = (/41.46_EB,40.32_EB,39.18_EB,38.05_EB,36.91_EB,35.87_EB,34.83_EB,33.78_EB,32.74_EB,31.70_EB/)
QFLAME(1,5,0:10,1) = (/0.00_EB,17.97_EB,22.33_EB,22.93_EB,23.53_EB,23.29_EB,23.06_EB,22.83_EB,22.59_EB,21.85_EB,21.11_EB/)
QFLAME(1,5,11:20,1) = (/20.36_EB,19.62_EB,18.91_EB,18.19_EB,17.48_EB,16.77_EB,16.21_EB,15.64_EB,15.08_EB,14.52_EB/)
QFLAME(1,5,21:30,1) = (/13.86_EB,13.19_EB,12.53_EB,11.86_EB,11.20_EB,10.83_EB,10.46_EB,10.09_EB,9.72_EB,9.35_EB/)
QFLAME(1,5,0:10,2) = (/0.00_EB,17.97_EB,22.33_EB,22.93_EB,23.53_EB,23.29_EB,23.05_EB,22.82_EB,22.58_EB,21.84_EB,21.11_EB/)
QFLAME(1,5,11:20,2) = (/20.37_EB,19.64_EB,18.92_EB,18.21_EB,17.49_EB,16.78_EB,16.22_EB,15.65_EB,15.09_EB,14.52_EB/)
QFLAME(1,5,21:30,2) = (/13.86_EB,13.19_EB,12.53_EB,11.86_EB,11.19_EB,10.84_EB,10.49_EB,10.14_EB,9.78_EB,9.43_EB/)
QFLAME(1,5,0:10,3) = (/0.00_EB,19.07_EB,23.12_EB,23.64_EB,24.15_EB,23.96_EB,23.77_EB,23.57_EB,23.38_EB,23.07_EB,22.76_EB/)
QFLAME(1,5,11:20,3) = (/22.45_EB,22.14_EB,21.53_EB,20.93_EB,20.32_EB,19.72_EB,19.27_EB,18.82_EB,18.38_EB,17.93_EB/)
QFLAME(1,5,21:30,3) = (/17.19_EB,16.46_EB,15.73_EB,15.00_EB,14.26_EB,13.98_EB,13.69_EB,13.40_EB,13.12_EB,12.83_EB/)
QFLAME(1,5,0:10,4) = (/0.00_EB,21.21_EB,26.79_EB,28.33_EB,29.86_EB,29.89_EB,29.91_EB,29.94_EB,29.96_EB,29.95_EB,29.94_EB/)
QFLAME(1,5,11:20,4) = (/29.93_EB,29.92_EB,29.26_EB,28.61_EB,27.95_EB,27.30_EB,26.74_EB,26.18_EB,25.63_EB,25.07_EB/)
QFLAME(1,5,21:30,4) = (/24.26_EB,23.46_EB,22.65_EB,21.85_EB,21.04_EB,20.38_EB,19.72_EB,19.06_EB,18.40_EB,17.74_EB/)
QFLAME(1,5,0:10,5) = (/0.00_EB,23.47_EB,30.73_EB,33.38_EB,36.03_EB,36.42_EB,36.80_EB,37.18_EB,37.56_EB,37.87_EB,38.18_EB/)
QFLAME(1,5,11:20,5) = (/38.49_EB,38.80_EB,38.16_EB,37.52_EB,36.88_EB,36.25_EB,35.24_EB,34.23_EB,33.22_EB,32.22_EB/)
QFLAME(1,5,21:30,5) = (/31.37_EB,30.52_EB,29.67_EB,28.82_EB,27.97_EB,27.10_EB,26.22_EB,25.35_EB,24.48_EB,23.61_EB/)
QFLAME(1,5,0:10,6) = (/0.00_EB,25.57_EB,34.60_EB,38.34_EB,42.08_EB,43.04_EB,44.01_EB,44.98_EB,45.95_EB,46.31_EB,46.68_EB/)
QFLAME(1,5,11:20,6) = (/47.05_EB,47.42_EB,47.14_EB,46.86_EB,46.58_EB,46.30_EB,45.35_EB,44.40_EB,43.44_EB,42.49_EB/)
QFLAME(1,5,21:30,6) = (/41.30_EB,40.10_EB,38.91_EB,37.72_EB,36.53_EB,35.50_EB,34.47_EB,33.44_EB,32.41_EB,31.38_EB/)
QFLAME(1,6,0:10,1) = (/0.00_EB,18.23_EB,22.78_EB,23.48_EB,24.18_EB,23.98_EB,23.78_EB,23.58_EB,23.38_EB,22.64_EB,21.90_EB/)
QFLAME(1,6,11:20,1) = (/21.16_EB,20.42_EB,19.70_EB,18.98_EB,18.26_EB,17.54_EB,16.97_EB,16.39_EB,15.81_EB,15.23_EB/)
QFLAME(1,6,21:30,1) = (/14.56_EB,13.89_EB,13.22_EB,12.55_EB,11.87_EB,11.46_EB,11.04_EB,10.63_EB,10.21_EB,9.80_EB/)
QFLAME(1,6,0:10,2) = (/0.00_EB,18.24_EB,22.79_EB,23.48_EB,24.18_EB,23.97_EB,23.77_EB,23.57_EB,23.37_EB,22.64_EB,21.90_EB/)
QFLAME(1,6,11:20,2) = (/21.17_EB,20.44_EB,19.72_EB,18.99_EB,18.27_EB,17.55_EB,16.97_EB,16.40_EB,15.82_EB,15.24_EB/)
QFLAME(1,6,21:30,2) = (/14.57_EB,13.89_EB,13.22_EB,12.54_EB,11.87_EB,11.47_EB,11.07_EB,10.67_EB,10.27_EB,9.87_EB/)
QFLAME(1,6,0:10,3) = (/0.00_EB,19.12_EB,23.42_EB,24.05_EB,24.68_EB,24.51_EB,24.35_EB,24.18_EB,24.01_EB,23.62_EB,23.23_EB/)
QFLAME(1,6,11:20,3) = (/22.83_EB,22.44_EB,21.81_EB,21.18_EB,20.55_EB,19.91_EB,19.44_EB,18.96_EB,18.48_EB,18.00_EB/)
QFLAME(1,6,21:30,3) = (/17.28_EB,16.55_EB,15.82_EB,15.10_EB,14.37_EB,14.06_EB,13.76_EB,13.45_EB,13.14_EB,12.84_EB/)
QFLAME(1,6,0:10,4) = (/0.00_EB,21.19_EB,26.74_EB,28.25_EB,29.77_EB,29.80_EB,29.83_EB,29.86_EB,29.89_EB,29.87_EB,29.86_EB/)
QFLAME(1,6,11:20,4) = (/29.84_EB,29.83_EB,29.18_EB,28.53_EB,27.87_EB,27.22_EB,26.66_EB,26.09_EB,25.53_EB,24.96_EB/)
QFLAME(1,6,21:30,4) = (/24.15_EB,23.34_EB,22.53_EB,21.72_EB,20.91_EB,20.25_EB,19.60_EB,18.95_EB,18.30_EB,17.65_EB/)
QFLAME(1,6,0:10,5) = (/0.00_EB,23.47_EB,30.68_EB,33.32_EB,35.95_EB,36.31_EB,36.66_EB,37.02_EB,37.37_EB,37.69_EB,38.00_EB/)
QFLAME(1,6,11:20,5) = (/38.32_EB,38.63_EB,38.00_EB,37.36_EB,36.73_EB,36.09_EB,35.11_EB,34.13_EB,33.14_EB,32.16_EB/)
QFLAME(1,6,21:30,5) = (/31.28_EB,30.40_EB,29.52_EB,28.64_EB,27.76_EB,26.92_EB,26.09_EB,25.25_EB,24.41_EB,23.57_EB/)
QFLAME(1,6,0:10,6) = (/0.00_EB,25.56_EB,34.56_EB,38.25_EB,41.94_EB,42.97_EB,44.01_EB,45.05_EB,46.09_EB,46.40_EB,46.72_EB/)
QFLAME(1,6,11:20,6) = (/47.03_EB,47.35_EB,47.00_EB,46.66_EB,46.32_EB,45.97_EB,45.06_EB,44.15_EB,43.24_EB,42.33_EB/)
QFLAME(1,6,21:30,6) = (/41.13_EB,39.93_EB,38.73_EB,37.52_EB,36.32_EB,35.33_EB,34.33_EB,33.33_EB,32.33_EB,31.34_EB/)
QFLAME(1,7,0:10,1) = (/0.00_EB,18.50_EB,23.24_EB,24.03_EB,24.83_EB,24.67_EB,24.50_EB,24.34_EB,24.17_EB,23.44_EB,22.70_EB/)
QFLAME(1,7,11:20,1) = (/21.96_EB,21.23_EB,20.50_EB,19.77_EB,19.05_EB,18.32_EB,17.73_EB,17.13_EB,16.54_EB,15.95_EB/)
QFLAME(1,7,21:30,1) = (/15.27_EB,14.59_EB,13.91_EB,13.23_EB,12.55_EB,12.08_EB,11.62_EB,11.16_EB,10.70_EB,10.24_EB/)
QFLAME(1,7,0:10,2) = (/0.00_EB,18.50_EB,23.24_EB,24.03_EB,24.83_EB,24.66_EB,24.49_EB,24.33_EB,24.16_EB,23.43_EB,22.70_EB/)
QFLAME(1,7,11:20,2) = (/21.97_EB,21.24_EB,20.51_EB,19.78_EB,19.05_EB,18.32_EB,17.73_EB,17.14_EB,16.55_EB,15.96_EB/)
QFLAME(1,7,21:30,2) = (/15.28_EB,14.59_EB,13.91_EB,13.22_EB,12.54_EB,12.09_EB,11.65_EB,11.20_EB,10.75_EB,10.30_EB/)
QFLAME(1,7,0:10,3) = (/0.00_EB,19.16_EB,23.72_EB,24.46_EB,25.20_EB,25.06_EB,24.92_EB,24.78_EB,24.64_EB,24.17_EB,23.69_EB/)
QFLAME(1,7,11:20,3) = (/23.21_EB,22.74_EB,22.08_EB,21.42_EB,20.77_EB,20.11_EB,19.60_EB,19.10_EB,18.59_EB,18.08_EB/)
QFLAME(1,7,21:30,3) = (/17.36_EB,16.64_EB,15.92_EB,15.20_EB,14.48_EB,14.15_EB,13.82_EB,13.50_EB,13.17_EB,12.84_EB/)
QFLAME(1,7,0:10,4) = (/0.00_EB,21.17_EB,26.69_EB,28.18_EB,29.67_EB,29.71_EB,29.74_EB,29.78_EB,29.81_EB,29.79_EB,29.78_EB/)
QFLAME(1,7,11:20,4) = (/29.76_EB,29.74_EB,29.09_EB,28.44_EB,27.80_EB,27.15_EB,26.58_EB,26.00_EB,25.43_EB,24.85_EB/)
QFLAME(1,7,21:30,4) = (/24.04_EB,23.22_EB,22.40_EB,21.59_EB,20.77_EB,20.13_EB,19.49_EB,18.85_EB,18.21_EB,17.57_EB/)
QFLAME(1,7,0:10,5) = (/0.00_EB,23.47_EB,30.64_EB,33.26_EB,35.87_EB,36.20_EB,36.53_EB,36.86_EB,37.18_EB,37.51_EB,37.83_EB/)
QFLAME(1,7,11:20,5) = (/38.15_EB,38.47_EB,37.84_EB,37.20_EB,36.57_EB,35.94_EB,34.98_EB,34.02_EB,33.06_EB,32.10_EB/)
QFLAME(1,7,21:30,5) = (/31.19_EB,30.28_EB,29.38_EB,28.47_EB,27.56_EB,26.75_EB,25.95_EB,25.14_EB,24.34_EB,23.53_EB/)
QFLAME(1,7,0:10,6) = (/0.00_EB,25.54_EB,34.53_EB,38.16_EB,41.80_EB,42.90_EB,44.01_EB,45.12_EB,46.22_EB,46.49_EB,46.75_EB/)
QFLAME(1,7,11:20,6) = (/47.01_EB,47.28_EB,46.87_EB,46.46_EB,46.05_EB,45.65_EB,44.78_EB,43.91_EB,43.04_EB,42.17_EB/)
QFLAME(1,7,21:30,6) = (/40.96_EB,39.75_EB,38.54_EB,37.33_EB,36.12_EB,35.16_EB,34.19_EB,33.23_EB,32.26_EB,31.29_EB/)
QFLAME(1,8,0:10,1) = (/0.00_EB,18.76_EB,23.69_EB,24.59_EB,25.48_EB,25.35_EB,25.22_EB,25.09_EB,24.96_EB,24.23_EB,23.50_EB/)
QFLAME(1,8,11:20,1) = (/22.76_EB,22.03_EB,21.30_EB,20.56_EB,19.83_EB,19.10_EB,18.49_EB,17.88_EB,17.27_EB,16.66_EB/)
QFLAME(1,8,21:30,1) = (/15.98_EB,15.29_EB,14.60_EB,13.91_EB,13.22_EB,12.71_EB,12.20_EB,11.70_EB,11.19_EB,10.68_EB/)
QFLAME(1,8,0:10,2) = (/0.00_EB,18.76_EB,23.70_EB,24.59_EB,25.48_EB,25.35_EB,25.21_EB,25.08_EB,24.95_EB,24.22_EB,23.50_EB/)
QFLAME(1,8,11:20,2) = (/22.77_EB,22.04_EB,21.30_EB,20.56_EB,19.82_EB,19.08_EB,18.48_EB,17.88_EB,17.28_EB,16.68_EB/)
QFLAME(1,8,21:30,2) = (/15.99_EB,15.29_EB,14.60_EB,13.91_EB,13.21_EB,12.72_EB,12.22_EB,11.73_EB,11.23_EB,10.74_EB/)
QFLAME(1,8,0:10,3) = (/0.00_EB,19.21_EB,24.01_EB,24.87_EB,25.73_EB,25.61_EB,25.50_EB,25.39_EB,25.27_EB,24.71_EB,24.16_EB/)
QFLAME(1,8,11:20,3) = (/23.60_EB,23.04_EB,22.36_EB,21.67_EB,20.99_EB,20.31_EB,19.77_EB,19.23_EB,18.70_EB,18.16_EB/)
QFLAME(1,8,21:30,3) = (/17.44_EB,16.73_EB,16.01_EB,15.30_EB,14.58_EB,14.24_EB,13.89_EB,13.54_EB,13.19_EB,12.85_EB/)
QFLAME(1,8,0:10,4) = (/0.00_EB,21.16_EB,26.64_EB,28.11_EB,29.58_EB,29.62_EB,29.66_EB,29.70_EB,29.74_EB,29.72_EB,29.69_EB/)
QFLAME(1,8,11:20,4) = (/29.67_EB,29.65_EB,29.01_EB,28.36_EB,27.72_EB,27.08_EB,26.50_EB,25.91_EB,25.33_EB,24.74_EB/)
QFLAME(1,8,21:30,4) = (/23.92_EB,23.10_EB,22.28_EB,21.46_EB,20.64_EB,20.01_EB,19.38_EB,18.75_EB,18.12_EB,17.49_EB/)
QFLAME(1,8,0:10,5) = (/0.00_EB,23.47_EB,30.60_EB,33.19_EB,35.79_EB,36.09_EB,36.39_EB,36.69_EB,37.00_EB,37.32_EB,37.65_EB/)
QFLAME(1,8,11:20,5) = (/37.98_EB,38.31_EB,37.67_EB,37.04_EB,36.41_EB,35.78_EB,34.85_EB,33.91_EB,32.98_EB,32.05_EB/)
QFLAME(1,8,21:30,5) = (/31.11_EB,30.17_EB,29.23_EB,28.29_EB,27.35_EB,26.58_EB,25.81_EB,25.04_EB,24.27_EB,23.50_EB/)
QFLAME(1,8,0:10,6) = (/0.00_EB,25.53_EB,34.49_EB,38.08_EB,41.66_EB,42.83_EB,44.01_EB,45.19_EB,46.36_EB,46.57_EB,46.78_EB/)
QFLAME(1,8,11:20,6) = (/46.99_EB,47.20_EB,46.73_EB,46.26_EB,45.79_EB,45.32_EB,44.49_EB,43.66_EB,42.83_EB,42.01_EB/)
QFLAME(1,8,21:30,6) = (/40.79_EB,39.57_EB,38.35_EB,37.14_EB,35.92_EB,34.99_EB,34.05_EB,33.12_EB,32.19_EB,31.25_EB/)
QFLAME(1,9,0:10,1) = (/0.00_EB,19.02_EB,24.15_EB,25.14_EB,26.13_EB,26.04_EB,25.94_EB,25.85_EB,25.75_EB,25.02_EB,24.29_EB/)
QFLAME(1,9,11:20,1) = (/23.56_EB,22.83_EB,22.09_EB,21.35_EB,20.61_EB,19.87_EB,19.25_EB,18.63_EB,18.00_EB,17.38_EB/)
QFLAME(1,9,21:30,1) = (/16.68_EB,15.98_EB,15.29_EB,14.59_EB,13.89_EB,13.34_EB,12.78_EB,12.23_EB,11.68_EB,11.12_EB/)
QFLAME(1,9,0:10,2) = (/0.00_EB,19.02_EB,24.15_EB,25.14_EB,26.13_EB,26.03_EB,25.93_EB,25.83_EB,25.74_EB,25.01_EB,24.29_EB/)
QFLAME(1,9,11:20,2) = (/23.57_EB,22.85_EB,22.10_EB,21.35_EB,20.60_EB,19.85_EB,19.24_EB,18.63_EB,18.01_EB,17.40_EB/)
QFLAME(1,9,21:30,2) = (/16.70_EB,15.99_EB,15.29_EB,14.59_EB,13.88_EB,13.34_EB,12.80_EB,12.26_EB,11.72_EB,11.18_EB/)
QFLAME(1,9,0:10,3) = (/0.00_EB,19.25_EB,24.31_EB,25.28_EB,26.25_EB,26.17_EB,26.08_EB,25.99_EB,25.90_EB,25.26_EB,24.62_EB/)
QFLAME(1,9,11:20,3) = (/23.98_EB,23.34_EB,22.63_EB,21.92_EB,21.21_EB,20.50_EB,19.94_EB,19.37_EB,18.80_EB,18.24_EB/)
QFLAME(1,9,21:30,3) = (/17.53_EB,16.82_EB,16.11_EB,15.40_EB,14.69_EB,14.32_EB,13.96_EB,13.59_EB,13.22_EB,12.85_EB/)
QFLAME(1,9,0:10,4) = (/0.00_EB,21.14_EB,26.59_EB,28.04_EB,29.49_EB,29.53_EB,29.57_EB,29.62_EB,29.66_EB,29.64_EB,29.61_EB/)
QFLAME(1,9,11:20,4) = (/29.58_EB,29.56_EB,28.92_EB,28.28_EB,27.65_EB,27.01_EB,26.42_EB,25.82_EB,25.23_EB,24.64_EB/)
QFLAME(1,9,21:30,4) = (/23.81_EB,22.98_EB,22.16_EB,21.33_EB,20.50_EB,19.88_EB,19.26_EB,18.64_EB,18.02_EB,17.40_EB/)
QFLAME(1,9,0:10,5) = (/0.00_EB,23.47_EB,30.55_EB,33.13_EB,35.71_EB,35.98_EB,36.26_EB,36.53_EB,36.81_EB,37.14_EB,37.47_EB/)
QFLAME(1,9,11:20,5) = (/37.81_EB,38.14_EB,37.51_EB,36.88_EB,36.25_EB,35.62_EB,34.72_EB,33.81_EB,32.90_EB,31.99_EB/)
QFLAME(1,9,21:30,5) = (/31.02_EB,30.05_EB,29.08_EB,28.12_EB,27.15_EB,26.41_EB,25.67_EB,24.93_EB,24.20_EB,23.46_EB/)
QFLAME(1,9,0:10,6) = (/0.00_EB,25.52_EB,34.46_EB,37.99_EB,41.52_EB,42.76_EB,44.01_EB,45.26_EB,46.50_EB,46.66_EB,46.82_EB/)
QFLAME(1,9,11:20,6) = (/46.98_EB,47.13_EB,46.60_EB,46.06_EB,45.53_EB,45.00_EB,44.21_EB,43.42_EB,42.63_EB,41.84_EB/)
QFLAME(1,9,21:30,6) = (/40.62_EB,39.39_EB,38.17_EB,36.94_EB,35.72_EB,34.82_EB,33.91_EB,33.01_EB,32.11_EB,31.21_EB/)
QFLAME(1,10,0:10,1) = (/0.00_EB,19.28_EB,24.61_EB,25.70_EB,26.79_EB,26.72_EB,26.66_EB,26.60_EB,26.54_EB,25.81_EB,25.09_EB/)
QFLAME(1,10,11:20,1) = (/24.36_EB,23.63_EB,22.89_EB,22.14_EB,21.39_EB,20.65_EB,20.01_EB,19.37_EB,18.73_EB,18.09_EB/)
QFLAME(1,10,21:30,1) = (/17.39_EB,16.68_EB,15.98_EB,15.27_EB,14.57_EB,13.97_EB,13.37_EB,12.76_EB,12.16_EB,11.56_EB/)
QFLAME(1,10,0:10,2) = (/0.00_EB,19.28_EB,24.61_EB,25.69_EB,26.78_EB,26.72_EB,26.65_EB,26.59_EB,26.53_EB,25.81_EB,25.09_EB/)
QFLAME(1,10,11:20,2) = (/24.37_EB,23.65_EB,22.89_EB,22.14_EB,21.38_EB,20.62_EB,20.00_EB,19.37_EB,18.74_EB,18.12_EB/)
QFLAME(1,10,21:30,2) = (/17.41_EB,16.69_EB,15.98_EB,15.27_EB,14.56_EB,13.97_EB,13.38_EB,12.79_EB,12.20_EB,11.61_EB/)
QFLAME(1,10,0:10,3) = (/0.00_EB,19.30_EB,24.61_EB,25.69_EB,26.78_EB,26.72_EB,26.66_EB,26.60_EB,26.53_EB,25.81_EB,25.09_EB/)
QFLAME(1,10,11:20,3) = (/24.36_EB,23.64_EB,22.90_EB,22.17_EB,21.43_EB,20.70_EB,20.10_EB,19.51_EB,18.91_EB,18.31_EB/)
QFLAME(1,10,21:30,3) = (/17.61_EB,16.91_EB,16.20_EB,15.50_EB,14.80_EB,14.41_EB,14.02_EB,13.63_EB,13.24_EB,12.86_EB/)
QFLAME(1,10,0:10,4) = (/0.00_EB,21.12_EB,26.54_EB,27.97_EB,29.39_EB,29.44_EB,29.49_EB,29.54_EB,29.59_EB,29.56_EB,29.53_EB/)
QFLAME(1,10,11:20,4) = (/29.50_EB,29.47_EB,28.84_EB,28.20_EB,27.57_EB,26.94_EB,26.33_EB,25.73_EB,25.13_EB,24.53_EB/)
QFLAME(1,10,21:30,4) = (/23.70_EB,22.86_EB,22.03_EB,21.20_EB,20.37_EB,19.76_EB,19.15_EB,18.54_EB,17.93_EB,17.32_EB/)
QFLAME(1,10,0:10,5) = (/0.00_EB,23.48_EB,30.51_EB,33.07_EB,35.62_EB,35.87_EB,36.12_EB,36.37_EB,36.62_EB,36.96_EB,37.30_EB/)
QFLAME(1,10,11:20,5) = (/37.64_EB,37.98_EB,37.35_EB,36.72_EB,36.10_EB,35.47_EB,34.59_EB,33.70_EB,32.82_EB,31.93_EB/)
QFLAME(1,10,21:30,5) = (/30.94_EB,29.94_EB,28.94_EB,27.94_EB,26.94_EB,26.24_EB,25.53_EB,24.83_EB,24.13_EB,23.42_EB/)
QFLAME(1,10,0:10,6) = (/0.00_EB,25.51_EB,34.42_EB,37.90_EB,41.38_EB,42.69_EB,44.01_EB,45.32_EB,46.64_EB,46.75_EB,46.85_EB/)
QFLAME(1,10,11:20,6) = (/46.96_EB,47.06_EB,46.46_EB,45.87_EB,45.27_EB,44.67_EB,43.92_EB,43.18_EB,42.43_EB,41.68_EB/)
QFLAME(1,10,21:30,6) = (/40.45_EB,39.22_EB,37.98_EB,36.75_EB,35.52_EB,34.65_EB,33.78_EB,32.91_EB,32.04_EB,31.17_EB/)
QFLAME(1,11,0:10,1) = (/0.00_EB,19.42_EB,24.83_EB,25.96_EB,27.09_EB,27.04_EB,26.99_EB,26.95_EB,26.90_EB,26.18_EB,25.45_EB/)
QFLAME(1,11,11:20,1) = (/24.73_EB,24.01_EB,23.26_EB,22.51_EB,21.76_EB,21.01_EB,20.37_EB,19.72_EB,19.08_EB,18.44_EB/)
QFLAME(1,11,21:30,1) = (/17.72_EB,17.01_EB,16.30_EB,15.59_EB,14.88_EB,14.25_EB,13.63_EB,13.00_EB,12.37_EB,11.75_EB/)
QFLAME(1,11,0:10,2) = (/0.00_EB,19.42_EB,24.83_EB,25.95_EB,27.08_EB,27.03_EB,26.98_EB,26.94_EB,26.89_EB,26.17_EB,25.45_EB/)
QFLAME(1,11,11:20,2) = (/24.74_EB,24.02_EB,23.26_EB,22.50_EB,21.75_EB,20.99_EB,20.36_EB,19.72_EB,19.09_EB,18.46_EB/)
QFLAME(1,11,21:30,2) = (/17.74_EB,17.02_EB,16.31_EB,15.59_EB,14.87_EB,14.26_EB,13.65_EB,13.03_EB,12.42_EB,11.81_EB/)
QFLAME(1,11,0:10,3) = (/0.00_EB,19.43_EB,24.83_EB,25.96_EB,27.08_EB,27.03_EB,26.99_EB,26.94_EB,26.89_EB,26.17_EB,25.45_EB/)
QFLAME(1,11,11:20,3) = (/24.73_EB,24.01_EB,23.27_EB,22.53_EB,21.80_EB,21.06_EB,20.45_EB,19.85_EB,19.24_EB,18.63_EB/)
QFLAME(1,11,21:30,3) = (/17.93_EB,17.22_EB,16.51_EB,15.80_EB,15.09_EB,14.66_EB,14.22_EB,13.79_EB,13.35_EB,12.92_EB/)
QFLAME(1,11,0:10,4) = (/0.00_EB,21.11_EB,26.58_EB,28.02_EB,29.45_EB,29.50_EB,29.55_EB,29.61_EB,29.66_EB,29.59_EB,29.52_EB/)
QFLAME(1,11,11:20,4) = (/29.45_EB,29.38_EB,28.76_EB,28.13_EB,27.51_EB,26.89_EB,26.28_EB,25.66_EB,25.05_EB,24.44_EB/)
QFLAME(1,11,21:30,4) = (/23.61_EB,22.79_EB,21.96_EB,21.14_EB,20.32_EB,19.70_EB,19.08_EB,18.46_EB,17.84_EB,17.22_EB/)
QFLAME(1,11,0:10,5) = (/0.00_EB,23.47_EB,30.46_EB,32.99_EB,35.52_EB,35.78_EB,36.04_EB,36.30_EB,36.56_EB,36.89_EB,37.21_EB/)
QFLAME(1,11,11:20,5) = (/37.54_EB,37.87_EB,37.25_EB,36.63_EB,36.01_EB,35.39_EB,34.51_EB,33.62_EB,32.74_EB,31.85_EB/)
QFLAME(1,11,21:30,5) = (/30.86_EB,29.87_EB,28.87_EB,27.88_EB,26.89_EB,26.17_EB,25.45_EB,24.73_EB,24.01_EB,23.28_EB/)
QFLAME(1,11,0:10,6) = (/0.00_EB,25.50_EB,34.38_EB,37.84_EB,41.29_EB,42.62_EB,43.95_EB,45.28_EB,46.61_EB,46.69_EB,46.78_EB/)
QFLAME(1,11,11:20,6) = (/46.86_EB,46.94_EB,46.35_EB,45.76_EB,45.16_EB,44.57_EB,43.83_EB,43.09_EB,42.35_EB,41.61_EB/)
QFLAME(1,11,21:30,6) = (/40.36_EB,39.11_EB,37.86_EB,36.61_EB,35.36_EB,34.51_EB,33.65_EB,32.79_EB,31.94_EB,31.08_EB/)
QFLAME(1,12,0:10,1) = (/0.00_EB,19.55_EB,25.05_EB,26.22_EB,27.39_EB,27.36_EB,27.32_EB,27.29_EB,27.26_EB,26.54_EB,25.82_EB/)
QFLAME(1,12,11:20,1) = (/25.10_EB,24.38_EB,23.63_EB,22.88_EB,22.13_EB,21.37_EB,20.72_EB,20.08_EB,19.43_EB,18.78_EB/)
QFLAME(1,12,21:30,1) = (/18.06_EB,17.34_EB,16.62_EB,15.91_EB,15.19_EB,14.54_EB,13.89_EB,13.24_EB,12.59_EB,11.93_EB/)
QFLAME(1,12,0:10,2) = (/0.00_EB,19.56_EB,25.05_EB,26.22_EB,27.38_EB,27.35_EB,27.32_EB,27.28_EB,27.25_EB,26.53_EB,25.82_EB/)
QFLAME(1,12,11:20,2) = (/25.10_EB,24.39_EB,23.63_EB,22.87_EB,22.11_EB,21.35_EB,20.72_EB,20.08_EB,19.44_EB,18.80_EB/)
QFLAME(1,12,21:30,2) = (/18.07_EB,17.35_EB,16.63_EB,15.91_EB,15.19_EB,14.55_EB,13.92_EB,13.28_EB,12.64_EB,12.01_EB/)
QFLAME(1,12,0:10,3) = (/0.00_EB,19.57_EB,25.05_EB,26.22_EB,27.38_EB,27.35_EB,27.32_EB,27.29_EB,27.25_EB,26.54_EB,25.82_EB/)
QFLAME(1,12,11:20,3) = (/25.10_EB,24.38_EB,23.64_EB,22.90_EB,22.16_EB,21.42_EB,20.80_EB,20.19_EB,19.57_EB,18.95_EB/)
QFLAME(1,12,21:30,3) = (/18.24_EB,17.53_EB,16.81_EB,16.10_EB,15.39_EB,14.91_EB,14.43_EB,13.94_EB,13.46_EB,12.98_EB/)
QFLAME(1,12,0:10,4) = (/0.00_EB,21.09_EB,26.62_EB,28.07_EB,29.51_EB,29.56_EB,29.62_EB,29.67_EB,29.73_EB,29.62_EB,29.51_EB/)
QFLAME(1,12,11:20,4) = (/29.40_EB,29.29_EB,28.68_EB,28.07_EB,27.45_EB,26.84_EB,26.22_EB,25.59_EB,24.97_EB,24.35_EB/)
QFLAME(1,12,21:30,4) = (/23.53_EB,22.71_EB,21.90_EB,21.08_EB,20.26_EB,19.63_EB,19.01_EB,18.38_EB,17.75_EB,17.12_EB/)
QFLAME(1,12,0:10,5) = (/0.00_EB,23.47_EB,30.42_EB,32.92_EB,35.42_EB,35.69_EB,35.96_EB,36.23_EB,36.50_EB,36.81_EB,37.13_EB/)
QFLAME(1,12,11:20,5) = (/37.44_EB,37.76_EB,37.15_EB,36.54_EB,35.93_EB,35.32_EB,34.43_EB,33.54_EB,32.66_EB,31.77_EB/)
QFLAME(1,12,21:30,5) = (/30.78_EB,29.80_EB,28.81_EB,27.82_EB,26.84_EB,26.10_EB,25.36_EB,24.62_EB,23.88_EB,23.15_EB/)
QFLAME(1,12,0:10,6) = (/0.00_EB,25.50_EB,34.33_EB,37.77_EB,41.21_EB,42.55_EB,43.89_EB,45.24_EB,46.58_EB,46.64_EB,46.70_EB/)
QFLAME(1,12,11:20,6) = (/46.76_EB,46.82_EB,46.23_EB,45.65_EB,45.06_EB,44.47_EB,43.74_EB,43.01_EB,42.27_EB,41.54_EB/)
QFLAME(1,12,21:30,6) = (/40.28_EB,39.01_EB,37.74_EB,36.47_EB,35.21_EB,34.36_EB,33.52_EB,32.68_EB,31.83_EB,30.99_EB/)
QFLAME(1,13,0:10,1) = (/0.00_EB,19.69_EB,25.27_EB,26.48_EB,27.69_EB,27.67_EB,27.65_EB,27.64_EB,27.62_EB,26.90_EB,26.19_EB/)
QFLAME(1,13,11:20,1) = (/25.47_EB,24.75_EB,24.00_EB,23.25_EB,22.49_EB,21.74_EB,21.08_EB,20.43_EB,19.77_EB,19.12_EB/)
QFLAME(1,13,21:30,1) = (/18.40_EB,17.67_EB,16.95_EB,16.22_EB,15.50_EB,14.82_EB,14.15_EB,13.47_EB,12.80_EB,12.12_EB/)
QFLAME(1,13,0:10,2) = (/0.00_EB,19.70_EB,25.27_EB,26.48_EB,27.68_EB,27.67_EB,27.65_EB,27.63_EB,27.61_EB,26.90_EB,26.18_EB/)
QFLAME(1,13,11:20,2) = (/25.47_EB,24.76_EB,24.00_EB,23.24_EB,22.48_EB,21.72_EB,21.07_EB,20.43_EB,19.78_EB,19.13_EB/)
QFLAME(1,13,21:30,2) = (/18.41_EB,17.68_EB,16.96_EB,16.23_EB,15.50_EB,14.84_EB,14.18_EB,13.53_EB,12.87_EB,12.21_EB/)
QFLAME(1,13,0:10,3) = (/0.00_EB,19.71_EB,25.27_EB,26.48_EB,27.68_EB,27.67_EB,27.65_EB,27.63_EB,27.61_EB,26.90_EB,26.18_EB/)
QFLAME(1,13,11:20,3) = (/25.47_EB,24.75_EB,24.01_EB,23.26_EB,22.52_EB,21.77_EB,21.15_EB,20.52_EB,19.90_EB,19.27_EB/)
QFLAME(1,13,21:30,3) = (/18.56_EB,17.84_EB,17.12_EB,16.40_EB,15.68_EB,15.16_EB,14.63_EB,14.10_EB,13.57_EB,13.04_EB/)
QFLAME(1,13,0:10,4) = (/0.00_EB,21.08_EB,26.67_EB,28.12_EB,29.56_EB,29.62_EB,29.68_EB,29.74_EB,29.80_EB,29.65_EB,29.50_EB/)
QFLAME(1,13,11:20,4) = (/29.35_EB,29.20_EB,28.60_EB,28.00_EB,27.40_EB,26.79_EB,26.16_EB,25.53_EB,24.89_EB,24.26_EB/)
QFLAME(1,13,21:30,4) = (/23.45_EB,22.64_EB,21.83_EB,21.02_EB,20.21_EB,19.57_EB,18.94_EB,18.30_EB,17.66_EB,17.03_EB/)
QFLAME(1,13,0:10,5) = (/0.00_EB,23.46_EB,30.37_EB,32.84_EB,35.31_EB,35.59_EB,35.87_EB,36.15_EB,36.43_EB,36.74_EB,37.04_EB/)
QFLAME(1,13,11:20,5) = (/37.35_EB,37.65_EB,37.05_EB,36.45_EB,35.84_EB,35.24_EB,34.35_EB,33.46_EB,32.57_EB,31.68_EB/)
QFLAME(1,13,21:30,5) = (/30.70_EB,29.72_EB,28.74_EB,27.76_EB,26.78_EB,26.03_EB,25.27_EB,24.52_EB,23.76_EB,23.01_EB/)
QFLAME(1,13,0:10,6) = (/0.00_EB,25.50_EB,34.29_EB,37.70_EB,41.12_EB,42.48_EB,43.84_EB,45.19_EB,46.55_EB,46.59_EB,46.63_EB/)
QFLAME(1,13,11:20,6) = (/46.66_EB,46.70_EB,46.12_EB,45.54_EB,44.95_EB,44.37_EB,43.65_EB,42.92_EB,42.20_EB,41.47_EB/)
QFLAME(1,13,21:30,6) = (/40.19_EB,38.90_EB,37.62_EB,36.34_EB,35.05_EB,34.22_EB,33.39_EB,32.56_EB,31.73_EB,30.90_EB/)
QFLAME(1,14,0:10,1) = (/0.00_EB,19.83_EB,25.49_EB,26.74_EB,27.99_EB,27.99_EB,27.98_EB,27.98_EB,27.98_EB,27.27_EB,26.55_EB/)
QFLAME(1,14,11:20,1) = (/25.84_EB,25.13_EB,24.37_EB,23.61_EB,22.86_EB,22.10_EB,21.44_EB,20.78_EB,20.12_EB,19.46_EB/)
QFLAME(1,14,21:30,1) = (/18.73_EB,18.00_EB,17.27_EB,16.54_EB,15.81_EB,15.11_EB,14.41_EB,13.71_EB,13.01_EB,12.30_EB/)
QFLAME(1,14,0:10,2) = (/0.00_EB,19.83_EB,25.49_EB,26.74_EB,27.99_EB,27.98_EB,27.98_EB,27.97_EB,27.97_EB,27.26_EB,26.55_EB/)
QFLAME(1,14,11:20,2) = (/25.84_EB,25.13_EB,24.37_EB,23.61_EB,22.85_EB,22.09_EB,21.43_EB,20.78_EB,20.13_EB,19.47_EB/)
QFLAME(1,14,21:30,2) = (/18.74_EB,18.01_EB,17.28_EB,16.55_EB,15.82_EB,15.14_EB,14.45_EB,13.77_EB,13.09_EB,12.40_EB/)
QFLAME(1,14,0:10,3) = (/0.00_EB,19.85_EB,25.49_EB,26.74_EB,27.98_EB,27.98_EB,27.98_EB,27.98_EB,27.97_EB,27.26_EB,26.55_EB/)
QFLAME(1,14,11:20,3) = (/25.84_EB,25.13_EB,24.38_EB,23.63_EB,22.88_EB,22.13_EB,21.50_EB,20.86_EB,20.23_EB,19.59_EB/)
QFLAME(1,14,21:30,3) = (/18.87_EB,18.15_EB,17.43_EB,16.70_EB,15.98_EB,15.40_EB,14.83_EB,14.26_EB,13.68_EB,13.11_EB/)
QFLAME(1,14,0:10,4) = (/0.00_EB,21.07_EB,26.71_EB,28.16_EB,29.62_EB,29.68_EB,29.74_EB,29.81_EB,29.87_EB,29.68_EB,29.49_EB/)
QFLAME(1,14,11:20,4) = (/29.30_EB,29.11_EB,28.52_EB,27.93_EB,27.34_EB,26.75_EB,26.10_EB,25.46_EB,24.81_EB,24.17_EB/)
QFLAME(1,14,21:30,4) = (/23.36_EB,22.56_EB,21.76_EB,20.96_EB,20.15_EB,19.51_EB,18.86_EB,18.22_EB,17.58_EB,16.93_EB/)
QFLAME(1,14,0:10,5) = (/0.00_EB,23.46_EB,30.33_EB,32.77_EB,35.21_EB,35.50_EB,35.79_EB,36.08_EB,36.37_EB,36.67_EB,36.96_EB/)
QFLAME(1,14,11:20,5) = (/37.25_EB,37.54_EB,36.95_EB,36.35_EB,35.76_EB,35.16_EB,34.27_EB,33.38_EB,32.49_EB,31.60_EB/)
QFLAME(1,14,21:30,5) = (/30.63_EB,29.65_EB,28.68_EB,27.71_EB,26.73_EB,25.96_EB,25.19_EB,24.42_EB,23.64_EB,22.87_EB/)
QFLAME(1,14,0:10,6) = (/0.00_EB,25.49_EB,34.24_EB,37.64_EB,41.03_EB,42.41_EB,43.78_EB,45.15_EB,46.52_EB,46.54_EB,46.55_EB/)
QFLAME(1,14,11:20,6) = (/46.57_EB,46.58_EB,46.00_EB,45.43_EB,44.85_EB,44.27_EB,43.56_EB,42.84_EB,42.12_EB,41.40_EB/)
QFLAME(1,14,21:30,6) = (/40.10_EB,38.80_EB,37.50_EB,36.20_EB,34.90_EB,34.08_EB,33.26_EB,32.45_EB,31.63_EB,30.81_EB/)
QFLAME(1,15,0:10,1) = (/0.00_EB,19.97_EB,25.71_EB,27.00_EB,28.29_EB,28.30_EB,28.31_EB,28.33_EB,28.34_EB,27.63_EB,26.92_EB/)
QFLAME(1,15,11:20,1) = (/26.21_EB,25.50_EB,24.74_EB,23.98_EB,23.22_EB,22.46_EB,21.80_EB,21.13_EB,20.47_EB,19.80_EB/)
QFLAME(1,15,21:30,1) = (/19.07_EB,18.33_EB,17.60_EB,16.86_EB,16.12_EB,15.40_EB,14.67_EB,13.94_EB,13.22_EB,12.49_EB/)
QFLAME(1,15,0:10,2) = (/0.00_EB,19.97_EB,25.71_EB,27.00_EB,28.29_EB,28.30_EB,28.31_EB,28.32_EB,28.33_EB,27.62_EB,26.92_EB/)
QFLAME(1,15,11:20,2) = (/26.21_EB,25.50_EB,24.74_EB,23.98_EB,23.22_EB,22.45_EB,21.79_EB,21.13_EB,20.47_EB,19.81_EB/)
QFLAME(1,15,21:30,2) = (/19.08_EB,18.34_EB,17.61_EB,16.87_EB,16.14_EB,15.43_EB,14.72_EB,14.02_EB,13.31_EB,12.60_EB/)
QFLAME(1,15,0:10,3) = (/0.00_EB,19.98_EB,25.71_EB,27.00_EB,28.29_EB,28.30_EB,28.31_EB,28.32_EB,28.33_EB,27.63_EB,26.92_EB/)
QFLAME(1,15,11:20,3) = (/26.21_EB,25.50_EB,24.75_EB,23.99_EB,23.24_EB,22.49_EB,21.85_EB,21.20_EB,20.56_EB,19.91_EB/)
QFLAME(1,15,21:30,3) = (/19.19_EB,18.46_EB,17.73_EB,17.00_EB,16.27_EB,15.65_EB,15.03_EB,14.41_EB,13.79_EB,13.17_EB/)
QFLAME(1,15,0:10,4) = (/0.00_EB,21.06_EB,26.75_EB,28.21_EB,29.68_EB,29.74_EB,29.81_EB,29.87_EB,29.94_EB,29.71_EB,29.48_EB/)
QFLAME(1,15,11:20,4) = (/29.25_EB,29.03_EB,28.44_EB,27.86_EB,27.28_EB,26.70_EB,26.04_EB,25.39_EB,24.73_EB,24.08_EB/)
QFLAME(1,15,21:30,4) = (/23.28_EB,22.48_EB,21.69_EB,20.89_EB,20.10_EB,19.45_EB,18.79_EB,18.14_EB,17.49_EB,16.83_EB/)
QFLAME(1,15,0:10,5) = (/0.00_EB,23.46_EB,30.28_EB,32.69_EB,35.10_EB,35.41_EB,35.71_EB,36.01_EB,36.31_EB,36.59_EB,36.87_EB/)
QFLAME(1,15,11:20,5) = (/37.15_EB,37.43_EB,36.85_EB,36.26_EB,35.67_EB,35.09_EB,34.20_EB,33.30_EB,32.41_EB,31.52_EB/)
QFLAME(1,15,21:30,5) = (/30.55_EB,29.58_EB,28.61_EB,27.65_EB,26.68_EB,25.89_EB,25.10_EB,24.31_EB,23.52_EB,22.73_EB/)
QFLAME(1,15,0:10,6) = (/0.00_EB,25.49_EB,34.19_EB,37.57_EB,40.95_EB,42.33_EB,43.72_EB,45.11_EB,46.49_EB,46.49_EB,46.48_EB/)
QFLAME(1,15,11:20,6) = (/46.47_EB,46.46_EB,45.89_EB,45.32_EB,44.75_EB,44.17_EB,43.46_EB,42.75_EB,42.04_EB,41.33_EB/)
QFLAME(1,15,21:30,6) = (/40.01_EB,38.70_EB,37.38_EB,36.06_EB,34.74_EB,33.94_EB,33.14_EB,32.33_EB,31.53_EB,30.73_EB/)
QFLAME(1,16,0:10,1) = (/0.00_EB,20.11_EB,25.93_EB,27.26_EB,28.59_EB,28.62_EB,28.64_EB,28.67_EB,28.70_EB,27.99_EB,27.29_EB/)
QFLAME(1,16,11:20,1) = (/26.58_EB,25.87_EB,25.11_EB,24.35_EB,23.59_EB,22.82_EB,22.15_EB,21.48_EB,20.81_EB,20.14_EB/)
QFLAME(1,16,21:30,1) = (/19.40_EB,18.66_EB,17.92_EB,17.18_EB,16.44_EB,15.68_EB,14.93_EB,14.18_EB,13.43_EB,12.68_EB/)
QFLAME(1,16,0:10,2) = (/0.00_EB,20.11_EB,25.94_EB,27.26_EB,28.59_EB,28.62_EB,28.64_EB,28.67_EB,28.69_EB,27.99_EB,27.28_EB/)
QFLAME(1,16,11:20,2) = (/26.58_EB,25.87_EB,25.11_EB,24.35_EB,23.58_EB,22.82_EB,22.15_EB,21.48_EB,20.82_EB,20.15_EB/)
QFLAME(1,16,21:30,2) = (/19.41_EB,18.67_EB,17.93_EB,17.19_EB,16.45_EB,15.72_EB,14.99_EB,14.26_EB,13.53_EB,12.80_EB/)
QFLAME(1,16,0:10,3) = (/0.00_EB,20.12_EB,25.94_EB,27.26_EB,28.59_EB,28.61_EB,28.64_EB,28.67_EB,28.69_EB,27.99_EB,27.28_EB/)
QFLAME(1,16,11:20,3) = (/26.58_EB,25.87_EB,25.11_EB,24.36_EB,23.60_EB,22.85_EB,22.20_EB,21.54_EB,20.89_EB,20.23_EB/)
QFLAME(1,16,21:30,3) = (/19.50_EB,18.77_EB,18.04_EB,17.30_EB,16.57_EB,15.90_EB,15.23_EB,14.57_EB,13.90_EB,13.23_EB/)
QFLAME(1,16,0:10,4) = (/0.00_EB,21.04_EB,26.79_EB,28.26_EB,29.74_EB,29.81_EB,29.87_EB,29.94_EB,30.01_EB,29.74_EB,29.47_EB/)
QFLAME(1,16,11:20,4) = (/29.20_EB,28.94_EB,28.37_EB,27.80_EB,27.22_EB,26.65_EB,25.99_EB,25.32_EB,24.65_EB,23.99_EB/)
QFLAME(1,16,21:30,4) = (/23.20_EB,22.41_EB,21.62_EB,20.83_EB,20.04_EB,19.38_EB,18.72_EB,18.06_EB,17.40_EB,16.74_EB/)
QFLAME(1,16,0:10,5) = (/0.00_EB,23.45_EB,30.23_EB,32.62_EB,35.00_EB,35.31_EB,35.63_EB,35.94_EB,36.25_EB,36.52_EB,36.79_EB/)
QFLAME(1,16,11:20,5) = (/37.05_EB,37.32_EB,36.74_EB,36.17_EB,35.59_EB,35.01_EB,34.12_EB,33.22_EB,32.33_EB,31.44_EB/)
QFLAME(1,16,21:30,5) = (/30.47_EB,29.51_EB,28.55_EB,27.59_EB,26.63_EB,25.82_EB,25.01_EB,24.21_EB,23.40_EB,22.59_EB/)
QFLAME(1,16,0:10,6) = (/0.00_EB,25.48_EB,34.15_EB,37.50_EB,40.86_EB,42.26_EB,43.66_EB,45.06_EB,46.46_EB,46.43_EB,46.40_EB/)
QFLAME(1,16,11:20,6) = (/46.37_EB,46.34_EB,45.77_EB,45.21_EB,44.64_EB,44.07_EB,43.37_EB,42.67_EB,41.97_EB,41.26_EB/)
QFLAME(1,16,21:30,6) = (/39.93_EB,38.59_EB,37.26_EB,35.92_EB,34.59_EB,33.80_EB,33.01_EB,32.22_EB,31.43_EB,30.64_EB/)
QFLAME(1,17,0:10,1) = (/0.00_EB,20.25_EB,26.15_EB,27.52_EB,28.89_EB,28.93_EB,28.98_EB,29.02_EB,29.06_EB,28.36_EB,27.65_EB/)
QFLAME(1,17,11:20,1) = (/26.95_EB,26.25_EB,25.48_EB,24.72_EB,23.95_EB,23.19_EB,22.51_EB,21.84_EB,21.16_EB,20.49_EB/)
QFLAME(1,17,21:30,1) = (/19.74_EB,18.99_EB,18.24_EB,17.49_EB,16.75_EB,15.97_EB,15.19_EB,14.42_EB,13.64_EB,12.86_EB/)
QFLAME(1,17,0:10,2) = (/0.00_EB,20.25_EB,26.16_EB,27.53_EB,28.89_EB,28.93_EB,28.97_EB,29.01_EB,29.05_EB,28.35_EB,27.65_EB/)
QFLAME(1,17,11:20,2) = (/26.94_EB,26.24_EB,25.48_EB,24.71_EB,23.95_EB,23.19_EB,22.51_EB,21.84_EB,21.16_EB,20.49_EB/)
QFLAME(1,17,21:30,2) = (/19.74_EB,19.00_EB,18.25_EB,17.51_EB,16.77_EB,16.01_EB,15.26_EB,14.51_EB,13.75_EB,13.00_EB/)
QFLAME(1,17,0:10,3) = (/0.00_EB,20.26_EB,26.16_EB,27.52_EB,28.89_EB,28.93_EB,28.97_EB,29.01_EB,29.05_EB,28.35_EB,27.65_EB/)
QFLAME(1,17,11:20,3) = (/26.94_EB,26.24_EB,25.48_EB,24.72_EB,23.97_EB,23.21_EB,22.54_EB,21.88_EB,21.22_EB,20.55_EB/)
QFLAME(1,17,21:30,3) = (/19.82_EB,19.08_EB,18.34_EB,17.60_EB,16.86_EB,16.15_EB,15.44_EB,14.72_EB,14.01_EB,13.29_EB/)
QFLAME(1,17,0:10,4) = (/0.00_EB,21.03_EB,26.83_EB,28.31_EB,29.80_EB,29.87_EB,29.94_EB,30.01_EB,30.07_EB,29.77_EB,29.46_EB/)
QFLAME(1,17,11:20,4) = (/29.16_EB,28.85_EB,28.29_EB,27.73_EB,27.17_EB,26.61_EB,25.93_EB,25.25_EB,24.57_EB,23.90_EB/)
QFLAME(1,17,21:30,4) = (/23.11_EB,22.33_EB,21.55_EB,20.77_EB,19.99_EB,19.32_EB,18.65_EB,17.98_EB,17.31_EB,16.64_EB/)
QFLAME(1,17,0:10,5) = (/0.00_EB,23.45_EB,30.19_EB,32.54_EB,34.89_EB,35.22_EB,35.54_EB,35.87_EB,36.19_EB,36.45_EB,36.70_EB/)
QFLAME(1,17,11:20,5) = (/36.96_EB,37.21_EB,36.64_EB,36.07_EB,35.50_EB,34.94_EB,34.04_EB,33.14_EB,32.25_EB,31.35_EB/)
QFLAME(1,17,21:30,5) = (/30.40_EB,29.44_EB,28.49_EB,27.53_EB,26.57_EB,25.75_EB,24.93_EB,24.10_EB,23.28_EB,22.46_EB/)
QFLAME(1,17,0:10,6) = (/0.00_EB,25.48_EB,34.10_EB,37.44_EB,40.78_EB,42.19_EB,43.61_EB,45.02_EB,46.44_EB,46.38_EB,46.33_EB/)
QFLAME(1,17,11:20,6) = (/46.27_EB,46.22_EB,45.66_EB,45.10_EB,44.54_EB,43.98_EB,43.28_EB,42.58_EB,41.89_EB,41.19_EB/)
QFLAME(1,17,21:30,6) = (/39.84_EB,38.49_EB,37.14_EB,35.79_EB,34.43_EB,33.66_EB,32.88_EB,32.10_EB,31.33_EB,30.55_EB/)
QFLAME(1,18,0:10,1) = (/0.00_EB,20.39_EB,26.37_EB,27.78_EB,29.19_EB,29.25_EB,29.31_EB,29.36_EB,29.42_EB,28.72_EB,28.02_EB/)
QFLAME(1,18,11:20,1) = (/27.32_EB,26.62_EB,25.85_EB,25.08_EB,24.32_EB,23.55_EB,22.87_EB,22.19_EB,21.51_EB,20.83_EB/)
QFLAME(1,18,21:30,1) = (/20.07_EB,19.32_EB,18.57_EB,17.81_EB,17.06_EB,16.26_EB,15.45_EB,14.65_EB,13.85_EB,13.05_EB/)
QFLAME(1,18,0:10,2) = (/0.00_EB,20.39_EB,26.38_EB,27.79_EB,29.19_EB,29.25_EB,29.30_EB,29.36_EB,29.42_EB,28.71_EB,28.01_EB/)
QFLAME(1,18,11:20,2) = (/27.31_EB,26.61_EB,25.85_EB,25.08_EB,24.32_EB,23.55_EB,22.87_EB,22.19_EB,21.51_EB,20.83_EB/)
QFLAME(1,18,21:30,2) = (/20.08_EB,19.33_EB,18.58_EB,17.83_EB,17.08_EB,16.31_EB,15.53_EB,14.75_EB,13.98_EB,13.20_EB/)
QFLAME(1,18,0:10,3) = (/0.00_EB,20.39_EB,26.38_EB,27.78_EB,29.19_EB,29.25_EB,29.30_EB,29.36_EB,29.41_EB,28.71_EB,28.01_EB/)
QFLAME(1,18,11:20,3) = (/27.31_EB,26.61_EB,25.85_EB,25.09_EB,24.33_EB,23.57_EB,22.89_EB,22.22_EB,21.55_EB,20.88_EB/)
QFLAME(1,18,21:30,3) = (/20.13_EB,19.39_EB,18.65_EB,17.90_EB,17.16_EB,16.40_EB,15.64_EB,14.88_EB,14.12_EB,13.36_EB/)
QFLAME(1,18,0:10,4) = (/0.00_EB,21.02_EB,26.87_EB,28.36_EB,29.85_EB,29.93_EB,30.00_EB,30.07_EB,30.14_EB,29.80_EB,29.45_EB/)
QFLAME(1,18,11:20,4) = (/29.11_EB,28.76_EB,28.21_EB,27.66_EB,27.11_EB,26.56_EB,25.87_EB,25.18_EB,24.49_EB,23.81_EB/)
QFLAME(1,18,21:30,4) = (/23.03_EB,22.26_EB,21.48_EB,20.71_EB,19.93_EB,19.26_EB,18.58_EB,17.90_EB,17.22_EB,16.55_EB/)
QFLAME(1,18,0:10,5) = (/0.00_EB,23.45_EB,30.14_EB,32.47_EB,34.79_EB,35.13_EB,35.46_EB,35.79_EB,36.13_EB,36.37_EB,36.62_EB/)
QFLAME(1,18,11:20,5) = (/36.86_EB,37.10_EB,36.54_EB,35.98_EB,35.42_EB,34.86_EB,33.96_EB,33.06_EB,32.17_EB,31.27_EB/)
QFLAME(1,18,21:30,5) = (/30.32_EB,29.37_EB,28.42_EB,27.47_EB,26.52_EB,25.68_EB,24.84_EB,24.00_EB,23.16_EB,22.32_EB/)
QFLAME(1,18,0:10,6) = (/0.00_EB,25.48_EB,34.05_EB,37.37_EB,40.69_EB,42.12_EB,43.55_EB,44.98_EB,46.41_EB,46.33_EB,46.25_EB/)
QFLAME(1,18,11:20,6) = (/46.18_EB,46.10_EB,45.54_EB,44.99_EB,44.43_EB,43.88_EB,43.19_EB,42.50_EB,41.81_EB,41.12_EB/)
QFLAME(1,18,21:30,6) = (/39.75_EB,38.39_EB,37.02_EB,35.65_EB,34.28_EB,33.52_EB,32.75_EB,31.99_EB,31.23_EB,30.46_EB/)
QFLAME(1,19,0:10,1) = (/0.00_EB,20.52_EB,26.59_EB,28.04_EB,29.49_EB,29.56_EB,29.64_EB,29.71_EB,29.78_EB,29.08_EB,28.39_EB/)
QFLAME(1,19,11:20,1) = (/27.69_EB,26.99_EB,26.22_EB,25.45_EB,24.68_EB,23.91_EB,23.23_EB,22.54_EB,21.86_EB,21.17_EB/)
QFLAME(1,19,21:30,1) = (/20.41_EB,19.65_EB,18.89_EB,18.13_EB,17.37_EB,16.54_EB,15.71_EB,14.89_EB,14.06_EB,13.23_EB/)
QFLAME(1,19,0:10,2) = (/0.00_EB,20.53_EB,26.60_EB,28.05_EB,29.50_EB,29.57_EB,29.64_EB,29.71_EB,29.78_EB,29.08_EB,28.38_EB/)
QFLAME(1,19,11:20,2) = (/27.68_EB,26.98_EB,26.22_EB,25.45_EB,24.69_EB,23.92_EB,23.23_EB,22.54_EB,21.85_EB,21.16_EB/)
QFLAME(1,19,21:30,2) = (/20.41_EB,19.66_EB,18.90_EB,18.15_EB,17.40_EB,16.60_EB,15.80_EB,15.00_EB,14.20_EB,13.40_EB/)
QFLAME(1,19,0:10,3) = (/0.00_EB,20.53_EB,26.60_EB,28.05_EB,29.49_EB,29.56_EB,29.63_EB,29.70_EB,29.77_EB,29.08_EB,28.38_EB/)
QFLAME(1,19,11:20,3) = (/27.68_EB,26.99_EB,26.22_EB,25.45_EB,24.69_EB,23.92_EB,23.24_EB,22.56_EB,21.88_EB,21.20_EB/)
QFLAME(1,19,21:30,3) = (/20.45_EB,19.70_EB,18.95_EB,18.20_EB,17.46_EB,16.65_EB,15.84_EB,15.03_EB,14.23_EB,13.42_EB/)
QFLAME(1,19,0:10,4) = (/0.00_EB,21.00_EB,26.91_EB,28.41_EB,29.91_EB,29.99_EB,30.06_EB,30.14_EB,30.21_EB,29.83_EB,29.44_EB/)
QFLAME(1,19,11:20,4) = (/29.06_EB,28.67_EB,28.13_EB,27.59_EB,27.05_EB,26.51_EB,25.81_EB,25.11_EB,24.41_EB,23.72_EB/)
QFLAME(1,19,21:30,4) = (/22.95_EB,22.18_EB,21.41_EB,20.65_EB,19.88_EB,19.19_EB,18.51_EB,17.82_EB,17.14_EB,16.45_EB/)
QFLAME(1,19,0:10,5) = (/0.00_EB,23.44_EB,30.10_EB,32.39_EB,34.69_EB,35.03_EB,35.38_EB,35.72_EB,36.07_EB,36.30_EB,36.53_EB/)
QFLAME(1,19,11:20,5) = (/36.76_EB,37.00_EB,36.44_EB,35.89_EB,35.34_EB,34.78_EB,33.88_EB,32.98_EB,32.09_EB,31.19_EB/)
QFLAME(1,19,21:30,5) = (/30.24_EB,29.30_EB,28.36_EB,27.41_EB,26.47_EB,25.61_EB,24.75_EB,23.90_EB,23.04_EB,22.18_EB/)
QFLAME(1,19,0:10,6) = (/0.00_EB,25.47_EB,34.01_EB,37.31_EB,40.60_EB,42.05_EB,43.49_EB,44.93_EB,46.38_EB,46.28_EB,46.18_EB/)
QFLAME(1,19,11:20,6) = (/46.08_EB,45.98_EB,45.43_EB,44.88_EB,44.33_EB,43.78_EB,43.10_EB,42.41_EB,41.73_EB,41.05_EB/)
QFLAME(1,19,21:30,6) = (/39.67_EB,38.28_EB,36.90_EB,35.51_EB,34.12_EB,33.37_EB,32.62_EB,31.87_EB,31.12_EB,30.37_EB/)
QFLAME(1,20,0:10,1) = (/0.00_EB,20.66_EB,26.82_EB,28.30_EB,29.79_EB,29.88_EB,29.97_EB,30.05_EB,30.14_EB,29.45_EB,28.75_EB/)
QFLAME(1,20,11:20,1) = (/28.06_EB,27.37_EB,26.59_EB,25.82_EB,25.05_EB,24.27_EB,23.58_EB,22.89_EB,22.20_EB,21.51_EB/)
QFLAME(1,20,21:30,1) = (/20.75_EB,19.98_EB,19.21_EB,18.45_EB,17.68_EB,16.83_EB,15.98_EB,15.12_EB,14.27_EB,13.42_EB/)
QFLAME(1,20,0:10,2) = (/0.00_EB,20.67_EB,26.82_EB,28.31_EB,29.80_EB,29.88_EB,29.97_EB,30.05_EB,30.14_EB,29.44_EB,28.74_EB/)
QFLAME(1,20,11:20,2) = (/28.05_EB,27.35_EB,26.59_EB,25.82_EB,25.05_EB,24.29_EB,23.59_EB,22.89_EB,22.20_EB,21.50_EB/)
QFLAME(1,20,21:30,2) = (/20.74_EB,19.99_EB,19.23_EB,18.47_EB,17.71_EB,16.89_EB,16.07_EB,15.24_EB,14.42_EB,13.60_EB/)
QFLAME(1,20,0:10,3) = (/0.00_EB,20.67_EB,26.82_EB,28.31_EB,29.79_EB,29.88_EB,29.96_EB,30.05_EB,30.13_EB,29.44_EB,28.75_EB/)
QFLAME(1,20,11:20,3) = (/28.05_EB,27.36_EB,26.59_EB,25.82_EB,25.05_EB,24.28_EB,23.59_EB,22.90_EB,22.21_EB,21.52_EB/)
QFLAME(1,20,21:30,3) = (/20.76_EB,20.01_EB,19.26_EB,18.50_EB,17.75_EB,16.90_EB,16.04_EB,15.19_EB,14.34_EB,13.48_EB/)
QFLAME(1,20,0:10,4) = (/0.00_EB,20.99_EB,26.95_EB,28.46_EB,29.97_EB,30.05_EB,30.13_EB,30.20_EB,30.28_EB,29.86_EB,29.43_EB/)
QFLAME(1,20,11:20,4) = (/29.01_EB,28.58_EB,28.05_EB,27.52_EB,26.99_EB,26.46_EB,25.75_EB,25.04_EB,24.34_EB,23.63_EB/)
QFLAME(1,20,21:30,4) = (/22.87_EB,22.11_EB,21.35_EB,20.59_EB,19.83_EB,19.13_EB,18.44_EB,17.74_EB,17.05_EB,16.35_EB/)
QFLAME(1,20,0:10,5) = (/0.00_EB,23.44_EB,30.05_EB,32.32_EB,34.58_EB,34.94_EB,35.29_EB,35.65_EB,36.01_EB,36.23_EB,36.45_EB/)
QFLAME(1,20,11:20,5) = (/36.67_EB,36.89_EB,36.34_EB,35.80_EB,35.25_EB,34.71_EB,33.81_EB,32.90_EB,32.00_EB,31.10_EB/)
QFLAME(1,20,21:30,5) = (/30.17_EB,29.23_EB,28.29_EB,27.35_EB,26.42_EB,25.54_EB,24.67_EB,23.79_EB,22.92_EB,22.04_EB/)
QFLAME(1,20,0:10,6) = (/0.00_EB,25.47_EB,33.96_EB,37.24_EB,40.52_EB,41.97_EB,43.43_EB,44.89_EB,46.35_EB,46.23_EB,46.10_EB/)
QFLAME(1,20,11:20,6) = (/45.98_EB,45.86_EB,45.31_EB,44.77_EB,44.22_EB,43.68_EB,43.00_EB,42.33_EB,41.66_EB,40.98_EB/)
QFLAME(1,20,21:30,6) = (/39.58_EB,38.18_EB,36.77_EB,35.37_EB,33.97_EB,33.23_EB,32.50_EB,31.76_EB,31.02_EB,30.29_EB/)
QFLAME(2,0,0:10,1) = (/0.00_EB,17.57_EB,21.57_EB,21.74_EB,21.91_EB,21.39_EB,20.88_EB,20.37_EB,19.86_EB,19.25_EB,18.63_EB/)
QFLAME(2,0,11:20,1) = (/18.02_EB,17.41_EB,16.80_EB,16.19_EB,15.58_EB,14.97_EB,14.55_EB,14.12_EB,13.69_EB,13.26_EB/)
QFLAME(2,0,21:30,1) = (/12.88_EB,12.51_EB,12.13_EB,11.76_EB,11.38_EB,11.09_EB,10.79_EB,10.50_EB,10.21_EB,9.92_EB/)
QFLAME(2,0,0:10,2) = (/0.00_EB,19.54_EB,24.15_EB,24.75_EB,25.35_EB,24.95_EB,24.55_EB,24.15_EB,23.75_EB,23.14_EB,22.53_EB/)
QFLAME(2,0,11:20,2) = (/21.91_EB,21.30_EB,20.88_EB,20.47_EB,20.05_EB,19.64_EB,19.40_EB,19.16_EB,18.92_EB,18.68_EB/)
QFLAME(2,0,21:30,2) = (/18.31_EB,17.94_EB,17.57_EB,17.20_EB,16.82_EB,16.60_EB,16.38_EB,16.15_EB,15.93_EB,15.70_EB/)
QFLAME(2,0,0:10,3) = (/0.00_EB,21.57_EB,27.67_EB,29.16_EB,30.66_EB,30.61_EB,30.56_EB,30.50_EB,30.45_EB,30.11_EB,29.76_EB/)
QFLAME(2,0,11:20,3) = (/29.41_EB,29.07_EB,28.71_EB,28.34_EB,27.98_EB,27.62_EB,27.64_EB,27.66_EB,27.68_EB,27.70_EB/)
QFLAME(2,0,21:30,3) = (/27.28_EB,26.86_EB,26.44_EB,26.02_EB,25.60_EB,25.28_EB,24.96_EB,24.64_EB,24.32_EB,24.00_EB/)
QFLAME(2,0,0:10,4) = (/0.00_EB,23.51_EB,31.27_EB,33.92_EB,36.56_EB,36.95_EB,37.34_EB,37.73_EB,38.12_EB,38.13_EB,38.14_EB/)
QFLAME(2,0,11:20,4) = (/38.15_EB,38.15_EB,37.51_EB,36.86_EB,36.21_EB,35.57_EB,36.16_EB,36.76_EB,37.35_EB,37.95_EB/)
QFLAME(2,0,21:30,4) = (/37.46_EB,36.98_EB,36.50_EB,36.02_EB,35.53_EB,35.10_EB,34.67_EB,34.24_EB,33.80_EB,33.37_EB/)
QFLAME(2,0,0:10,5) = (/0.00_EB,25.40_EB,34.98_EB,38.80_EB,42.63_EB,43.72_EB,44.81_EB,45.91_EB,47.00_EB,47.10_EB,47.20_EB/)
QFLAME(2,0,11:20,5) = (/47.30_EB,47.39_EB,47.05_EB,46.71_EB,46.36_EB,46.02_EB,46.72_EB,47.43_EB,48.13_EB,48.84_EB/)
QFLAME(2,0,21:30,5) = (/48.47_EB,48.10_EB,47.73_EB,47.36_EB,46.99_EB,46.52_EB,46.05_EB,45.58_EB,45.11_EB,44.63_EB/)
QFLAME(2,0,0:10,6) = (/0.00_EB,27.24_EB,38.83_EB,43.81_EB,48.78_EB,50.71_EB,52.64_EB,54.56_EB,56.49_EB,57.04_EB,57.58_EB/)
QFLAME(2,0,11:20,6) = (/58.13_EB,58.68_EB,58.47_EB,58.26_EB,58.05_EB,57.84_EB,57.63_EB,57.42_EB,57.21_EB,57.00_EB/)
QFLAME(2,0,21:30,6) = (/57.43_EB,57.86_EB,58.28_EB,58.71_EB,59.14_EB,58.74_EB,58.34_EB,57.94_EB,57.54_EB,57.13_EB/)
QFLAME(2,1,0:10,1) = (/0.00_EB,18.11_EB,22.66_EB,23.13_EB,23.61_EB,23.25_EB,22.88_EB,22.51_EB,22.15_EB,21.57_EB,20.98_EB/)
QFLAME(2,1,11:20,1) = (/20.40_EB,19.82_EB,19.23_EB,18.63_EB,18.04_EB,17.44_EB,16.98_EB,16.53_EB,16.07_EB,15.62_EB/)
QFLAME(2,1,21:30,1) = (/15.21_EB,14.81_EB,14.40_EB,14.00_EB,13.59_EB,13.27_EB,12.95_EB,12.63_EB,12.31_EB,11.99_EB/)
QFLAME(2,1,0:10,2) = (/0.00_EB,19.52_EB,24.17_EB,24.74_EB,25.31_EB,24.88_EB,24.45_EB,24.02_EB,23.58_EB,22.98_EB,22.37_EB/)
QFLAME(2,1,11:20,2) = (/21.77_EB,21.16_EB,20.72_EB,20.29_EB,19.85_EB,19.42_EB,19.18_EB,18.94_EB,18.70_EB,18.46_EB/)
QFLAME(2,1,21:30,2) = (/18.14_EB,17.83_EB,17.51_EB,17.19_EB,16.88_EB,16.59_EB,16.30_EB,16.02_EB,15.73_EB,15.44_EB/)
QFLAME(2,1,0:10,3) = (/0.00_EB,21.58_EB,27.66_EB,29.14_EB,30.61_EB,30.52_EB,30.42_EB,30.33_EB,30.23_EB,29.85_EB,29.47_EB/)
QFLAME(2,1,11:20,3) = (/29.09_EB,28.71_EB,28.36_EB,28.00_EB,27.65_EB,27.29_EB,27.30_EB,27.31_EB,27.32_EB,27.33_EB/)
QFLAME(2,1,21:30,3) = (/26.83_EB,26.32_EB,25.81_EB,25.30_EB,24.80_EB,24.55_EB,24.30_EB,24.06_EB,23.81_EB,23.56_EB/)
QFLAME(2,1,0:10,4) = (/0.00_EB,23.50_EB,31.28_EB,33.88_EB,36.48_EB,36.81_EB,37.14_EB,37.47_EB,37.80_EB,37.81_EB,37.82_EB/)
QFLAME(2,1,11:20,4) = (/37.83_EB,37.84_EB,37.22_EB,36.61_EB,35.99_EB,35.38_EB,36.00_EB,36.62_EB,37.24_EB,37.86_EB/)
QFLAME(2,1,21:30,4) = (/37.23_EB,36.61_EB,35.99_EB,35.36_EB,34.74_EB,34.41_EB,34.08_EB,33.75_EB,33.42_EB,33.09_EB/)
QFLAME(2,1,0:10,5) = (/0.00_EB,25.38_EB,34.98_EB,38.78_EB,42.57_EB,43.57_EB,44.57_EB,45.57_EB,46.57_EB,46.70_EB,46.84_EB/)
QFLAME(2,1,11:20,5) = (/46.97_EB,47.10_EB,46.75_EB,46.40_EB,46.05_EB,45.69_EB,46.31_EB,46.92_EB,47.54_EB,48.15_EB/)
QFLAME(2,1,21:30,5) = (/47.84_EB,47.54_EB,47.23_EB,46.92_EB,46.61_EB,46.08_EB,45.55_EB,45.02_EB,44.49_EB,43.96_EB/)
QFLAME(2,1,0:10,6) = (/0.00_EB,27.22_EB,38.79_EB,43.76_EB,48.73_EB,50.65_EB,52.56_EB,54.48_EB,56.40_EB,56.84_EB,57.29_EB/)
QFLAME(2,1,11:20,6) = (/57.73_EB,58.17_EB,57.92_EB,57.66_EB,57.41_EB,57.16_EB,56.92_EB,56.68_EB,56.44_EB,56.20_EB/)
QFLAME(2,1,21:30,6) = (/56.74_EB,57.28_EB,57.82_EB,58.35_EB,58.89_EB,58.51_EB,58.12_EB,57.74_EB,57.35_EB,56.96_EB/)
QFLAME(2,2,0:10,1) = (/0.00_EB,18.59_EB,23.55_EB,24.29_EB,25.02_EB,24.77_EB,24.52_EB,24.28_EB,24.03_EB,23.47_EB,22.92_EB/)
QFLAME(2,2,11:20,1) = (/22.36_EB,21.81_EB,21.18_EB,20.56_EB,19.93_EB,19.31_EB,18.86_EB,18.41_EB,17.97_EB,17.52_EB/)
QFLAME(2,2,21:30,1) = (/17.10_EB,16.68_EB,16.25_EB,15.83_EB,15.41_EB,15.08_EB,14.74_EB,14.41_EB,14.08_EB,13.75_EB/)
QFLAME(2,2,0:10,2) = (/0.00_EB,19.50_EB,24.18_EB,24.86_EB,25.55_EB,25.25_EB,24.96_EB,24.66_EB,24.37_EB,23.77_EB,23.18_EB/)
QFLAME(2,2,11:20,2) = (/22.58_EB,21.98_EB,21.45_EB,20.91_EB,20.38_EB,19.84_EB,19.49_EB,19.13_EB,18.78_EB,18.42_EB/)
QFLAME(2,2,21:30,2) = (/18.10_EB,17.78_EB,17.46_EB,17.14_EB,16.81_EB,16.52_EB,16.22_EB,15.93_EB,15.64_EB,15.34_EB/)
QFLAME(2,2,0:10,3) = (/0.00_EB,21.52_EB,27.64_EB,29.09_EB,30.54_EB,30.43_EB,30.32_EB,30.20_EB,30.09_EB,29.72_EB,29.34_EB/)
QFLAME(2,2,11:20,3) = (/28.96_EB,28.58_EB,28.16_EB,27.73_EB,27.31_EB,26.88_EB,26.89_EB,26.90_EB,26.90_EB,26.91_EB/)
QFLAME(2,2,21:30,3) = (/26.54_EB,26.17_EB,25.80_EB,25.42_EB,25.05_EB,24.69_EB,24.33_EB,23.97_EB,23.62_EB,23.26_EB/)
QFLAME(2,2,0:10,4) = (/0.00_EB,23.44_EB,31.27_EB,33.84_EB,36.41_EB,36.73_EB,37.06_EB,37.38_EB,37.70_EB,37.68_EB,37.66_EB/)
QFLAME(2,2,11:20,4) = (/37.65_EB,37.63_EB,36.99_EB,36.34_EB,35.69_EB,35.04_EB,35.54_EB,36.05_EB,36.55_EB,37.06_EB/)
QFLAME(2,2,21:30,4) = (/36.60_EB,36.14_EB,35.68_EB,35.22_EB,34.76_EB,34.35_EB,33.94_EB,33.53_EB,33.13_EB,32.72_EB/)
QFLAME(2,2,0:10,5) = (/0.00_EB,25.31_EB,34.96_EB,38.73_EB,42.50_EB,43.56_EB,44.61_EB,45.67_EB,46.72_EB,46.72_EB,46.72_EB/)
QFLAME(2,2,11:20,5) = (/46.72_EB,46.72_EB,46.38_EB,46.04_EB,45.69_EB,45.35_EB,46.04_EB,46.72_EB,47.41_EB,48.10_EB/)
QFLAME(2,2,21:30,5) = (/47.69_EB,47.28_EB,46.88_EB,46.47_EB,46.06_EB,45.56_EB,45.05_EB,44.54_EB,44.03_EB,43.53_EB/)
QFLAME(2,2,0:10,6) = (/0.00_EB,27.27_EB,38.74_EB,43.72_EB,48.70_EB,50.58_EB,52.46_EB,54.35_EB,56.23_EB,56.69_EB,57.15_EB/)
QFLAME(2,2,11:20,6) = (/57.60_EB,58.06_EB,57.79_EB,57.51_EB,57.24_EB,56.97_EB,56.73_EB,56.50_EB,56.27_EB,56.03_EB/)
QFLAME(2,2,21:30,6) = (/56.46_EB,56.88_EB,57.31_EB,57.73_EB,58.16_EB,57.79_EB,57.41_EB,57.04_EB,56.66_EB,56.29_EB/)
QFLAME(2,3,0:10,1) = (/0.00_EB,18.96_EB,24.22_EB,25.14_EB,26.07_EB,25.90_EB,25.72_EB,25.55_EB,25.37_EB,24.85_EB,24.34_EB/)
QFLAME(2,3,11:20,1) = (/23.82_EB,23.31_EB,22.67_EB,22.04_EB,21.41_EB,20.78_EB,20.33_EB,19.89_EB,19.44_EB,19.00_EB/)
QFLAME(2,3,21:30,1) = (/18.56_EB,18.13_EB,17.69_EB,17.25_EB,16.82_EB,16.47_EB,16.12_EB,15.78_EB,15.43_EB,15.08_EB/)
QFLAME(2,3,0:10,2) = (/0.00_EB,19.58_EB,24.63_EB,25.53_EB,26.42_EB,26.22_EB,26.01_EB,25.81_EB,25.61_EB,25.06_EB,24.51_EB/)
QFLAME(2,3,11:20,2) = (/23.96_EB,23.42_EB,22.85_EB,22.27_EB,21.70_EB,21.13_EB,20.75_EB,20.37_EB,19.98_EB,19.60_EB/)
QFLAME(2,3,21:30,2) = (/19.23_EB,18.87_EB,18.50_EB,18.14_EB,17.77_EB,17.44_EB,17.11_EB,16.79_EB,16.46_EB,16.13_EB/)
QFLAME(2,3,0:10,3) = (/0.00_EB,21.52_EB,27.62_EB,29.05_EB,30.49_EB,30.37_EB,30.25_EB,30.13_EB,30.01_EB,29.63_EB,29.25_EB/)
QFLAME(2,3,11:20,3) = (/28.86_EB,28.48_EB,28.05_EB,27.62_EB,27.20_EB,26.77_EB,26.78_EB,26.78_EB,26.79_EB,26.80_EB/)
QFLAME(2,3,21:30,3) = (/26.41_EB,26.01_EB,25.61_EB,25.21_EB,24.81_EB,24.47_EB,24.13_EB,23.80_EB,23.46_EB,23.12_EB/)
QFLAME(2,3,0:10,4) = (/0.00_EB,23.45_EB,31.24_EB,33.79_EB,36.34_EB,36.66_EB,36.98_EB,37.29_EB,37.61_EB,37.58_EB,37.55_EB/)
QFLAME(2,3,11:20,4) = (/37.52_EB,37.49_EB,36.86_EB,36.23_EB,35.60_EB,34.97_EB,35.47_EB,35.97_EB,36.47_EB,36.97_EB/)
QFLAME(2,3,21:30,4) = (/36.51_EB,36.04_EB,35.58_EB,35.12_EB,34.66_EB,34.24_EB,33.83_EB,33.41_EB,32.99_EB,32.58_EB/)
QFLAME(2,3,0:10,5) = (/0.00_EB,25.32_EB,34.93_EB,38.68_EB,42.43_EB,43.50_EB,44.56_EB,45.62_EB,46.69_EB,46.66_EB,46.63_EB/)
QFLAME(2,3,11:20,5) = (/46.60_EB,46.57_EB,46.23_EB,45.89_EB,45.55_EB,45.21_EB,45.88_EB,46.55_EB,47.22_EB,47.89_EB/)
QFLAME(2,3,21:30,5) = (/47.48_EB,47.07_EB,46.66_EB,46.26_EB,45.85_EB,45.32_EB,44.80_EB,44.28_EB,43.75_EB,43.23_EB/)
QFLAME(2,3,0:10,6) = (/0.00_EB,27.24_EB,38.71_EB,43.68_EB,48.64_EB,50.50_EB,52.36_EB,54.22_EB,56.08_EB,56.49_EB,56.91_EB/)
QFLAME(2,3,11:20,6) = (/57.32_EB,57.73_EB,57.45_EB,57.16_EB,56.88_EB,56.59_EB,56.41_EB,56.22_EB,56.03_EB,55.84_EB/)
QFLAME(2,3,21:30,6) = (/56.27_EB,56.70_EB,57.13_EB,57.56_EB,57.99_EB,57.55_EB,57.10_EB,56.66_EB,56.22_EB,55.78_EB/)
QFLAME(2,4,0:10,1) = (/0.00_EB,19.33_EB,24.88_EB,26.00_EB,27.12_EB,27.02_EB,26.92_EB,26.81_EB,26.71_EB,26.24_EB,25.76_EB/)
QFLAME(2,4,11:20,1) = (/25.28_EB,24.80_EB,24.16_EB,23.53_EB,22.89_EB,22.25_EB,21.81_EB,21.36_EB,20.92_EB,20.47_EB/)
QFLAME(2,4,21:30,1) = (/20.02_EB,19.58_EB,19.13_EB,18.68_EB,18.23_EB,17.86_EB,17.50_EB,17.14_EB,16.77_EB,16.41_EB/)
QFLAME(2,4,0:10,2) = (/0.00_EB,19.65_EB,25.09_EB,26.19_EB,27.29_EB,27.18_EB,27.07_EB,26.96_EB,26.85_EB,26.35_EB,25.85_EB/)
QFLAME(2,4,11:20,2) = (/25.35_EB,24.85_EB,24.24_EB,23.64_EB,23.03_EB,22.43_EB,22.01_EB,21.60_EB,21.19_EB,20.78_EB/)
QFLAME(2,4,21:30,2) = (/20.37_EB,19.96_EB,19.55_EB,19.14_EB,18.73_EB,18.37_EB,18.00_EB,17.64_EB,17.28_EB,16.91_EB/)
QFLAME(2,4,0:10,3) = (/0.00_EB,21.51_EB,27.60_EB,29.02_EB,30.43_EB,30.31_EB,30.18_EB,30.05_EB,29.92_EB,29.54_EB,29.16_EB/)
QFLAME(2,4,11:20,3) = (/28.77_EB,28.39_EB,27.95_EB,27.52_EB,27.08_EB,26.65_EB,26.66_EB,26.67_EB,26.68_EB,26.70_EB/)
QFLAME(2,4,21:30,3) = (/26.27_EB,25.85_EB,25.42_EB,25.00_EB,24.57_EB,24.25_EB,23.94_EB,23.62_EB,23.30_EB,22.98_EB/)
QFLAME(2,4,0:10,4) = (/0.00_EB,23.45_EB,31.21_EB,33.74_EB,36.26_EB,36.58_EB,36.90_EB,37.21_EB,37.53_EB,37.49_EB,37.44_EB/)
QFLAME(2,4,11:20,4) = (/37.39_EB,37.35_EB,36.73_EB,36.12_EB,35.51_EB,34.89_EB,35.39_EB,35.89_EB,36.38_EB,36.88_EB/)
QFLAME(2,4,21:30,4) = (/36.41_EB,35.95_EB,35.49_EB,35.02_EB,34.56_EB,34.13_EB,33.71_EB,33.29_EB,32.86_EB,32.44_EB/)
QFLAME(2,4,0:10,5) = (/0.00_EB,25.33_EB,34.89_EB,38.63_EB,42.36_EB,43.44_EB,44.51_EB,45.58_EB,46.65_EB,46.59_EB,46.54_EB/)
QFLAME(2,4,11:20,5) = (/46.48_EB,46.42_EB,46.08_EB,45.75_EB,45.41_EB,45.08_EB,45.73_EB,46.38_EB,47.03_EB,47.68_EB/)
QFLAME(2,4,21:30,5) = (/47.27_EB,46.86_EB,46.45_EB,46.04_EB,45.63_EB,45.09_EB,44.55_EB,44.01_EB,43.47_EB,42.93_EB/)
QFLAME(2,4,0:10,6) = (/0.00_EB,27.21_EB,38.69_EB,43.64_EB,48.58_EB,50.42_EB,52.26_EB,54.10_EB,55.93_EB,56.30_EB,56.67_EB/)
QFLAME(2,4,11:20,6) = (/57.03_EB,57.40_EB,57.11_EB,56.81_EB,56.52_EB,56.22_EB,56.08_EB,55.93_EB,55.79_EB,55.64_EB/)
QFLAME(2,4,21:30,6) = (/56.08_EB,56.51_EB,56.95_EB,57.38_EB,57.82_EB,57.31_EB,56.80_EB,56.29_EB,55.78_EB,55.27_EB/)
QFLAME(2,5,0:10,1) = (/0.00_EB,19.70_EB,25.55_EB,26.86_EB,28.17_EB,28.14_EB,28.11_EB,28.08_EB,28.06_EB,27.62_EB,27.18_EB/)
QFLAME(2,5,11:20,1) = (/26.74_EB,26.30_EB,25.65_EB,25.01_EB,24.37_EB,23.72_EB,23.28_EB,22.84_EB,22.39_EB,21.95_EB/)
QFLAME(2,5,21:30,1) = (/21.49_EB,21.03_EB,20.56_EB,20.10_EB,19.64_EB,19.26_EB,18.88_EB,18.50_EB,18.12_EB,17.74_EB/)
QFLAME(2,5,0:10,2) = (/0.00_EB,19.73_EB,25.55_EB,26.86_EB,28.17_EB,28.15_EB,28.13_EB,28.11_EB,28.09_EB,27.64_EB,27.18_EB/)
QFLAME(2,5,11:20,2) = (/26.73_EB,26.28_EB,25.64_EB,25.00_EB,24.36_EB,23.72_EB,23.28_EB,22.84_EB,22.40_EB,21.96_EB/)
QFLAME(2,5,21:30,2) = (/21.50_EB,21.05_EB,20.60_EB,20.15_EB,19.69_EB,19.29_EB,18.89_EB,18.50_EB,18.10_EB,17.70_EB/)
QFLAME(2,5,0:10,3) = (/0.00_EB,21.50_EB,27.58_EB,28.98_EB,30.38_EB,30.25_EB,30.11_EB,29.98_EB,29.84_EB,29.45_EB,29.06_EB/)
QFLAME(2,5,11:20,3) = (/28.68_EB,28.29_EB,27.85_EB,27.41_EB,26.97_EB,26.53_EB,26.54_EB,26.56_EB,26.57_EB,26.59_EB/)
QFLAME(2,5,21:30,3) = (/26.14_EB,25.69_EB,25.23_EB,24.78_EB,24.33_EB,24.03_EB,23.74_EB,23.44_EB,23.14_EB,22.85_EB/)
QFLAME(2,5,0:10,4) = (/0.00_EB,23.46_EB,31.18_EB,33.68_EB,36.19_EB,36.50_EB,36.82_EB,37.13_EB,37.45_EB,37.39_EB,37.33_EB/)
QFLAME(2,5,11:20,4) = (/37.27_EB,37.20_EB,36.61_EB,36.01_EB,35.42_EB,34.82_EB,35.31_EB,35.81_EB,36.30_EB,36.79_EB/)
QFLAME(2,5,21:30,4) = (/36.32_EB,35.86_EB,35.39_EB,34.92_EB,34.46_EB,34.03_EB,33.59_EB,33.16_EB,32.73_EB,32.30_EB/)
QFLAME(2,5,0:10,5) = (/0.00_EB,25.34_EB,34.86_EB,38.57_EB,42.29_EB,43.37_EB,44.46_EB,45.54_EB,46.62_EB,46.53_EB,46.44_EB/)
QFLAME(2,5,11:20,5) = (/46.35_EB,46.27_EB,45.93_EB,45.60_EB,45.27_EB,44.94_EB,45.58_EB,46.21_EB,46.84_EB,47.48_EB/)
QFLAME(2,5,21:30,5) = (/47.07_EB,46.65_EB,46.24_EB,45.83_EB,45.42_EB,44.86_EB,44.30_EB,43.75_EB,43.19_EB,42.64_EB/)
QFLAME(2,5,0:10,6) = (/0.00_EB,27.18_EB,38.67_EB,43.60_EB,48.53_EB,50.34_EB,52.16_EB,53.97_EB,55.79_EB,56.11_EB,56.43_EB/)
QFLAME(2,5,11:20,6) = (/56.75_EB,57.07_EB,56.76_EB,56.46_EB,56.15_EB,55.85_EB,55.75_EB,55.65_EB,55.55_EB,55.45_EB/)
QFLAME(2,5,21:30,6) = (/55.89_EB,56.33_EB,56.77_EB,57.21_EB,57.65_EB,57.07_EB,56.49_EB,55.91_EB,55.33_EB,54.75_EB/)
QFLAME(2,6,0:10,1) = (/0.00_EB,19.92_EB,25.97_EB,27.40_EB,28.82_EB,28.85_EB,28.87_EB,28.89_EB,28.91_EB,28.49_EB,28.07_EB/)
QFLAME(2,6,11:20,1) = (/27.64_EB,27.22_EB,26.58_EB,25.94_EB,25.30_EB,24.66_EB,24.22_EB,23.77_EB,23.33_EB,22.89_EB/)
QFLAME(2,6,21:30,1) = (/22.42_EB,21.95_EB,21.48_EB,21.02_EB,20.55_EB,20.16_EB,19.77_EB,19.38_EB,18.99_EB,18.60_EB/)
QFLAME(2,6,0:10,2) = (/0.00_EB,19.94_EB,25.97_EB,27.39_EB,28.82_EB,28.85_EB,28.88_EB,28.91_EB,28.94_EB,28.51_EB,28.08_EB/)
QFLAME(2,6,11:20,2) = (/27.65_EB,27.21_EB,26.57_EB,25.93_EB,25.29_EB,24.66_EB,24.22_EB,23.78_EB,23.34_EB,22.90_EB/)
QFLAME(2,6,21:30,2) = (/22.44_EB,21.98_EB,21.52_EB,21.06_EB,20.60_EB,20.19_EB,19.78_EB,19.38_EB,18.97_EB,18.56_EB/)
QFLAME(2,6,0:10,3) = (/0.00_EB,21.49_EB,27.65_EB,29.14_EB,30.64_EB,30.57_EB,30.51_EB,30.44_EB,30.37_EB,29.99_EB,29.60_EB/)
QFLAME(2,6,11:20,3) = (/29.21_EB,28.82_EB,28.35_EB,27.87_EB,27.40_EB,26.92_EB,26.86_EB,26.80_EB,26.73_EB,26.67_EB/)
QFLAME(2,6,21:30,3) = (/26.22_EB,25.77_EB,25.32_EB,24.87_EB,24.42_EB,24.10_EB,23.79_EB,23.47_EB,23.16_EB,22.85_EB/)
QFLAME(2,6,0:10,4) = (/0.00_EB,23.44_EB,31.15_EB,33.64_EB,36.12_EB,36.43_EB,36.74_EB,37.04_EB,37.35_EB,37.28_EB,37.21_EB/)
QFLAME(2,6,11:20,4) = (/37.15_EB,37.08_EB,36.52_EB,35.97_EB,35.41_EB,34.85_EB,35.30_EB,35.75_EB,36.19_EB,36.64_EB/)
QFLAME(2,6,21:30,4) = (/36.17_EB,35.70_EB,35.23_EB,34.76_EB,34.30_EB,33.86_EB,33.43_EB,33.00_EB,32.56_EB,32.13_EB/)
QFLAME(2,6,0:10,5) = (/0.00_EB,25.33_EB,34.82_EB,38.53_EB,42.23_EB,43.28_EB,44.32_EB,45.37_EB,46.42_EB,46.35_EB,46.28_EB/)
QFLAME(2,6,11:20,5) = (/46.21_EB,46.14_EB,45.79_EB,45.44_EB,45.09_EB,44.74_EB,45.39_EB,46.05_EB,46.71_EB,47.37_EB/)
QFLAME(2,6,21:30,5) = (/46.95_EB,46.52_EB,46.10_EB,45.68_EB,45.25_EB,44.72_EB,44.18_EB,43.65_EB,43.11_EB,42.58_EB/)
QFLAME(2,6,0:10,6) = (/0.00_EB,27.18_EB,38.63_EB,43.55_EB,48.47_EB,50.28_EB,52.09_EB,53.90_EB,55.71_EB,56.00_EB,56.29_EB/)
QFLAME(2,6,11:20,6) = (/56.58_EB,56.87_EB,56.54_EB,56.22_EB,55.89_EB,55.57_EB,55.56_EB,55.55_EB,55.54_EB,55.53_EB/)
QFLAME(2,6,21:30,6) = (/55.90_EB,56.27_EB,56.65_EB,57.02_EB,57.40_EB,56.82_EB,56.25_EB,55.68_EB,55.11_EB,54.54_EB/)
QFLAME(2,7,0:10,1) = (/0.00_EB,20.14_EB,26.39_EB,27.94_EB,29.48_EB,29.55_EB,29.62_EB,29.70_EB,29.77_EB,29.36_EB,28.96_EB/)
QFLAME(2,7,11:20,1) = (/28.55_EB,28.15_EB,27.51_EB,26.87_EB,26.23_EB,25.59_EB,25.15_EB,24.71_EB,24.27_EB,23.83_EB/)
QFLAME(2,7,21:30,1) = (/23.35_EB,22.88_EB,22.41_EB,21.93_EB,21.46_EB,21.06_EB,20.66_EB,20.26_EB,19.86_EB,19.47_EB/)
QFLAME(2,7,0:10,2) = (/0.00_EB,20.16_EB,26.39_EB,27.93_EB,29.47_EB,29.55_EB,29.64_EB,29.72_EB,29.80_EB,29.39_EB,28.97_EB/)
QFLAME(2,7,11:20,2) = (/28.56_EB,28.14_EB,27.50_EB,26.87_EB,26.23_EB,25.59_EB,25.15_EB,24.71_EB,24.27_EB,23.83_EB/)
QFLAME(2,7,21:30,2) = (/23.37_EB,22.90_EB,22.43_EB,21.97_EB,21.50_EB,21.09_EB,20.67_EB,20.26_EB,19.84_EB,19.43_EB/)
QFLAME(2,7,0:10,3) = (/0.00_EB,21.48_EB,27.71_EB,29.31_EB,30.90_EB,30.90_EB,30.90_EB,30.91_EB,30.91_EB,30.52_EB,30.13_EB/)
QFLAME(2,7,11:20,3) = (/29.75_EB,29.36_EB,28.85_EB,28.34_EB,27.82_EB,27.31_EB,27.17_EB,27.03_EB,26.89_EB,26.75_EB/)
QFLAME(2,7,21:30,3) = (/26.30_EB,25.85_EB,25.40_EB,24.95_EB,24.51_EB,24.17_EB,23.84_EB,23.51_EB,23.18_EB,22.85_EB/)
QFLAME(2,7,0:10,4) = (/0.00_EB,23.42_EB,31.12_EB,33.59_EB,36.06_EB,36.36_EB,36.65_EB,36.95_EB,37.25_EB,37.18_EB,37.10_EB/)
QFLAME(2,7,11:20,4) = (/37.03_EB,36.95_EB,36.43_EB,35.92_EB,35.40_EB,34.88_EB,35.29_EB,35.69_EB,36.09_EB,36.49_EB/)
QFLAME(2,7,21:30,4) = (/36.02_EB,35.55_EB,35.08_EB,34.60_EB,34.13_EB,33.70_EB,33.26_EB,32.83_EB,32.39_EB,31.96_EB/)
QFLAME(2,7,0:10,5) = (/0.00_EB,25.31_EB,34.78_EB,38.48_EB,42.17_EB,43.18_EB,44.19_EB,45.20_EB,46.22_EB,46.17_EB,46.12_EB/)
QFLAME(2,7,11:20,5) = (/46.06_EB,46.01_EB,45.64_EB,45.27_EB,44.90_EB,44.53_EB,45.21_EB,45.90_EB,46.58_EB,47.26_EB/)
QFLAME(2,7,21:30,5) = (/46.83_EB,46.39_EB,45.96_EB,45.53_EB,45.09_EB,44.58_EB,44.06_EB,43.55_EB,43.03_EB,42.52_EB/)
QFLAME(2,7,0:10,6) = (/0.00_EB,27.18_EB,38.59_EB,43.50_EB,48.42_EB,50.22_EB,52.03_EB,53.84_EB,55.64_EB,55.90_EB,56.15_EB/)
QFLAME(2,7,11:20,6) = (/56.41_EB,56.67_EB,56.32_EB,55.98_EB,55.64_EB,55.29_EB,55.37_EB,55.45_EB,55.52_EB,55.60_EB/)
QFLAME(2,7,21:30,6) = (/55.91_EB,56.22_EB,56.53_EB,56.84_EB,57.14_EB,56.58_EB,56.02_EB,55.46_EB,54.89_EB,54.33_EB/)
QFLAME(2,8,0:10,1) = (/0.00_EB,20.36_EB,26.81_EB,28.47_EB,30.14_EB,30.26_EB,30.38_EB,30.50_EB,30.62_EB,30.24_EB,29.85_EB/)
QFLAME(2,8,11:20,1) = (/29.46_EB,29.07_EB,28.44_EB,27.80_EB,27.16_EB,26.53_EB,26.09_EB,25.65_EB,25.21_EB,24.77_EB/)
QFLAME(2,8,21:30,1) = (/24.29_EB,23.81_EB,23.33_EB,22.85_EB,22.37_EB,21.96_EB,21.55_EB,21.14_EB,20.74_EB,20.33_EB/)
QFLAME(2,8,0:10,2) = (/0.00_EB,20.37_EB,26.80_EB,28.46_EB,30.12_EB,30.26_EB,30.39_EB,30.52_EB,30.66_EB,30.26_EB,29.87_EB/)
QFLAME(2,8,11:20,2) = (/29.47_EB,29.07_EB,28.44_EB,27.80_EB,27.16_EB,26.53_EB,26.09_EB,25.65_EB,25.21_EB,24.77_EB/)
QFLAME(2,8,21:30,2) = (/24.30_EB,23.83_EB,23.35_EB,22.88_EB,22.41_EB,21.98_EB,21.56_EB,21.14_EB,20.72_EB,20.29_EB/)
QFLAME(2,8,0:10,3) = (/0.00_EB,21.47_EB,27.78_EB,29.47_EB,31.16_EB,31.23_EB,31.30_EB,31.37_EB,31.44_EB,31.05_EB,30.67_EB/)
QFLAME(2,8,11:20,3) = (/30.28_EB,29.90_EB,29.35_EB,28.80_EB,28.25_EB,27.71_EB,27.49_EB,27.27_EB,27.05_EB,26.84_EB/)
QFLAME(2,8,21:30,3) = (/26.39_EB,25.94_EB,25.49_EB,25.04_EB,24.59_EB,24.24_EB,23.89_EB,23.54_EB,23.19_EB,22.85_EB/)
QFLAME(2,8,0:10,4) = (/0.00_EB,23.41_EB,31.08_EB,33.54_EB,35.99_EB,36.28_EB,36.57_EB,36.86_EB,37.15_EB,37.07_EB,36.99_EB/)
QFLAME(2,8,11:20,4) = (/36.91_EB,36.83_EB,36.35_EB,35.87_EB,35.39_EB,34.91_EB,35.27_EB,35.63_EB,35.99_EB,36.34_EB/)
QFLAME(2,8,21:30,4) = (/35.87_EB,35.39_EB,34.92_EB,34.45_EB,33.97_EB,33.53_EB,33.10_EB,32.66_EB,32.23_EB,31.79_EB/)
QFLAME(2,8,0:10,5) = (/0.00_EB,25.30_EB,34.75_EB,38.43_EB,42.10_EB,43.08_EB,44.06_EB,45.04_EB,46.02_EB,45.98_EB,45.95_EB/)
QFLAME(2,8,11:20,5) = (/45.92_EB,45.89_EB,45.50_EB,45.10_EB,44.71_EB,44.32_EB,45.03_EB,45.74_EB,46.45_EB,47.16_EB/)
QFLAME(2,8,21:30,5) = (/46.71_EB,46.27_EB,45.82_EB,45.38_EB,44.93_EB,44.44_EB,43.94_EB,43.45_EB,42.95_EB,42.46_EB/)
QFLAME(2,8,0:10,6) = (/0.00_EB,27.18_EB,38.55_EB,43.46_EB,48.36_EB,50.16_EB,51.97_EB,53.77_EB,55.57_EB,55.79_EB,56.02_EB/)
QFLAME(2,8,11:20,6) = (/56.24_EB,56.46_EB,56.10_EB,55.74_EB,55.38_EB,55.01_EB,55.18_EB,55.35_EB,55.51_EB,55.68_EB/)
QFLAME(2,8,21:30,6) = (/55.92_EB,56.16_EB,56.41_EB,56.65_EB,56.89_EB,56.34_EB,55.78_EB,55.23_EB,54.67_EB,54.12_EB/)
QFLAME(2,9,0:10,1) = (/0.00_EB,20.58_EB,27.23_EB,29.01_EB,30.79_EB,30.96_EB,31.14_EB,31.31_EB,31.48_EB,31.11_EB,30.74_EB/)
QFLAME(2,9,11:20,1) = (/30.37_EB,30.00_EB,29.36_EB,28.73_EB,28.10_EB,27.46_EB,27.03_EB,26.59_EB,26.15_EB,25.71_EB/)
QFLAME(2,9,21:30,1) = (/25.22_EB,24.73_EB,24.25_EB,23.76_EB,23.28_EB,22.86_EB,22.44_EB,22.02_EB,21.61_EB,21.19_EB/)
QFLAME(2,9,0:10,2) = (/0.00_EB,20.58_EB,27.22_EB,29.00_EB,30.78_EB,30.96_EB,31.15_EB,31.33_EB,31.52_EB,31.14_EB,30.76_EB/)
QFLAME(2,9,11:20,2) = (/30.38_EB,30.00_EB,29.37_EB,28.73_EB,28.10_EB,27.46_EB,27.02_EB,26.59_EB,26.15_EB,25.71_EB/)
QFLAME(2,9,21:30,2) = (/25.23_EB,24.75_EB,24.27_EB,23.79_EB,23.31_EB,22.88_EB,22.45_EB,22.02_EB,21.59_EB,21.16_EB/)
QFLAME(2,9,0:10,3) = (/0.00_EB,21.46_EB,27.84_EB,29.63_EB,31.42_EB,31.56_EB,31.70_EB,31.84_EB,31.97_EB,31.59_EB,31.20_EB/)
QFLAME(2,9,11:20,3) = (/30.82_EB,30.43_EB,29.85_EB,29.26_EB,28.68_EB,28.10_EB,27.80_EB,27.51_EB,27.21_EB,26.92_EB/)
QFLAME(2,9,21:30,3) = (/26.47_EB,26.02_EB,25.57_EB,25.13_EB,24.68_EB,24.31_EB,23.95_EB,23.58_EB,23.21_EB,22.85_EB/)
QFLAME(2,9,0:10,4) = (/0.00_EB,23.39_EB,31.05_EB,33.49_EB,35.93_EB,36.21_EB,36.49_EB,36.77_EB,37.05_EB,36.96_EB,36.88_EB/)
QFLAME(2,9,11:20,4) = (/36.79_EB,36.70_EB,36.26_EB,35.82_EB,35.38_EB,34.94_EB,35.26_EB,35.57_EB,35.88_EB,36.19_EB/)
QFLAME(2,9,21:30,4) = (/35.72_EB,35.24_EB,34.76_EB,34.29_EB,33.81_EB,33.37_EB,32.93_EB,32.49_EB,32.06_EB,31.62_EB/)
QFLAME(2,9,0:10,5) = (/0.00_EB,25.29_EB,34.71_EB,38.38_EB,42.04_EB,42.99_EB,43.93_EB,44.87_EB,45.81_EB,45.80_EB,45.79_EB/)
QFLAME(2,9,11:20,5) = (/45.78_EB,45.76_EB,45.35_EB,44.94_EB,44.53_EB,44.11_EB,44.85_EB,45.58_EB,46.31_EB,47.05_EB/)
QFLAME(2,9,21:30,5) = (/46.59_EB,46.14_EB,45.68_EB,45.22_EB,44.77_EB,44.29_EB,43.82_EB,43.34_EB,42.87_EB,42.40_EB/)
QFLAME(2,9,0:10,6) = (/0.00_EB,27.18_EB,38.51_EB,43.41_EB,48.31_EB,50.10_EB,51.90_EB,53.70_EB,55.50_EB,55.69_EB,55.88_EB/)
QFLAME(2,9,11:20,6) = (/56.07_EB,56.26_EB,55.88_EB,55.50_EB,55.12_EB,54.74_EB,54.99_EB,55.24_EB,55.50_EB,55.75_EB/)
QFLAME(2,9,21:30,6) = (/55.93_EB,56.11_EB,56.29_EB,56.46_EB,56.64_EB,56.09_EB,55.55_EB,55.00_EB,54.46_EB,53.91_EB/)
QFLAME(2,10,0:10,1) = (/0.00_EB,20.80_EB,27.65_EB,29.55_EB,31.45_EB,31.67_EB,31.89_EB,32.11_EB,32.34_EB,31.98_EB,31.63_EB/)
QFLAME(2,10,11:20,1) = (/31.27_EB,30.92_EB,30.29_EB,29.66_EB,29.03_EB,28.40_EB,27.96_EB,27.52_EB,27.09_EB,26.65_EB/)
QFLAME(2,10,21:30,1) = (/26.15_EB,25.66_EB,25.17_EB,24.68_EB,24.19_EB,23.76_EB,23.33_EB,22.90_EB,22.48_EB,22.05_EB/)
QFLAME(2,10,0:10,2) = (/0.00_EB,20.80_EB,27.64_EB,29.54_EB,31.43_EB,31.66_EB,31.90_EB,32.14_EB,32.37_EB,32.01_EB,31.65_EB/)
QFLAME(2,10,11:20,2) = (/31.29_EB,30.94_EB,30.30_EB,29.67_EB,29.03_EB,28.40_EB,27.96_EB,27.52_EB,27.08_EB,26.65_EB/)
QFLAME(2,10,21:30,2) = (/26.16_EB,25.67_EB,25.19_EB,24.70_EB,24.22_EB,23.78_EB,23.34_EB,22.90_EB,22.46_EB,22.02_EB/)
QFLAME(2,10,0:10,3) = (/0.00_EB,21.46_EB,27.90_EB,29.79_EB,31.68_EB,31.89_EB,32.10_EB,32.30_EB,32.51_EB,32.12_EB,31.74_EB/)
QFLAME(2,10,11:20,3) = (/31.35_EB,30.97_EB,30.35_EB,29.73_EB,29.11_EB,28.49_EB,28.12_EB,27.75_EB,27.37_EB,27.00_EB/)
QFLAME(2,10,21:30,3) = (/26.55_EB,26.11_EB,25.66_EB,25.21_EB,24.77_EB,24.38_EB,24.00_EB,23.61_EB,23.23_EB,22.85_EB/)
QFLAME(2,10,0:10,4) = (/0.00_EB,23.37_EB,31.02_EB,33.44_EB,35.87_EB,36.14_EB,36.41_EB,36.68_EB,36.95_EB,36.86_EB,36.76_EB/)
QFLAME(2,10,11:20,4) = (/36.67_EB,36.57_EB,36.17_EB,35.77_EB,35.37_EB,34.98_EB,35.24_EB,35.51_EB,35.78_EB,36.05_EB/)
QFLAME(2,10,21:30,4) = (/35.57_EB,35.09_EB,34.61_EB,34.13_EB,33.65_EB,33.21_EB,32.77_EB,32.33_EB,31.89_EB,31.45_EB/)
QFLAME(2,10,0:10,5) = (/0.00_EB,25.28_EB,34.68_EB,38.33_EB,41.98_EB,42.89_EB,43.80_EB,44.71_EB,45.61_EB,45.62_EB,45.63_EB/)
QFLAME(2,10,11:20,5) = (/45.63_EB,45.64_EB,45.20_EB,44.77_EB,44.34_EB,43.91_EB,44.66_EB,45.42_EB,46.18_EB,46.94_EB/)
QFLAME(2,10,21:30,5) = (/46.47_EB,46.01_EB,45.54_EB,45.07_EB,44.61_EB,44.15_EB,43.70_EB,43.24_EB,42.79_EB,42.34_EB/)
QFLAME(2,10,0:10,6) = (/0.00_EB,27.18_EB,38.47_EB,43.36_EB,48.25_EB,50.04_EB,51.84_EB,53.63_EB,55.43_EB,55.59_EB,55.74_EB/)
QFLAME(2,10,11:20,6) = (/55.90_EB,56.06_EB,55.66_EB,55.26_EB,54.86_EB,54.46_EB,54.80_EB,55.14_EB,55.49_EB,55.83_EB/)
QFLAME(2,10,21:30,6) = (/55.94_EB,56.05_EB,56.17_EB,56.28_EB,56.39_EB,55.85_EB,55.31_EB,54.77_EB,54.24_EB,53.70_EB/)
QFLAME(2,11,0:10,1) = (/0.00_EB,20.92_EB,27.88_EB,29.84_EB,31.80_EB,32.05_EB,32.29_EB,32.54_EB,32.79_EB,32.44_EB,32.10_EB/)
QFLAME(2,11,11:20,1) = (/31.75_EB,31.40_EB,30.77_EB,30.15_EB,29.52_EB,28.90_EB,28.46_EB,28.01_EB,27.57_EB,27.13_EB/)
QFLAME(2,11,21:30,1) = (/26.64_EB,26.15_EB,25.65_EB,25.16_EB,24.67_EB,24.24_EB,23.81_EB,23.37_EB,22.94_EB,22.51_EB/)
QFLAME(2,11,0:10,2) = (/0.00_EB,20.92_EB,27.88_EB,29.83_EB,31.78_EB,32.04_EB,32.30_EB,32.56_EB,32.82_EB,32.47_EB,32.12_EB/)
QFLAME(2,11,11:20,2) = (/31.76_EB,31.41_EB,30.78_EB,30.15_EB,29.52_EB,28.90_EB,28.45_EB,28.01_EB,27.57_EB,27.13_EB/)
QFLAME(2,11,21:30,2) = (/26.64_EB,26.16_EB,25.67_EB,25.18_EB,24.69_EB,24.25_EB,23.81_EB,23.36_EB,22.92_EB,22.48_EB/)
QFLAME(2,11,0:10,3) = (/0.00_EB,21.52_EB,28.11_EB,30.06_EB,32.01_EB,32.24_EB,32.48_EB,32.71_EB,32.94_EB,32.57_EB,32.19_EB/)
QFLAME(2,11,11:20,3) = (/31.82_EB,31.44_EB,30.83_EB,30.21_EB,29.60_EB,28.98_EB,28.60_EB,28.21_EB,27.83_EB,27.45_EB/)
QFLAME(2,11,21:30,3) = (/27.00_EB,26.55_EB,26.09_EB,25.64_EB,25.19_EB,24.80_EB,24.40_EB,24.01_EB,23.62_EB,23.22_EB/)
QFLAME(2,11,0:10,4) = (/0.00_EB,23.37_EB,30.99_EB,33.43_EB,35.86_EB,36.15_EB,36.44_EB,36.73_EB,37.02_EB,36.91_EB,36.80_EB/)
QFLAME(2,11,11:20,4) = (/36.69_EB,36.57_EB,36.16_EB,35.75_EB,35.34_EB,34.93_EB,35.17_EB,35.41_EB,35.66_EB,35.90_EB/)
QFLAME(2,11,21:30,4) = (/35.43_EB,34.95_EB,34.48_EB,34.01_EB,33.54_EB,33.11_EB,32.67_EB,32.23_EB,31.80_EB,31.36_EB/)
QFLAME(2,11,0:10,5) = (/0.00_EB,25.27_EB,34.64_EB,38.28_EB,41.92_EB,42.83_EB,43.73_EB,44.64_EB,45.54_EB,45.54_EB,45.53_EB/)
QFLAME(2,11,11:20,5) = (/45.53_EB,45.52_EB,45.10_EB,44.68_EB,44.26_EB,43.84_EB,44.58_EB,45.31_EB,46.05_EB,46.78_EB/)
QFLAME(2,11,21:30,5) = (/46.32_EB,45.86_EB,45.39_EB,44.93_EB,44.46_EB,44.01_EB,43.55_EB,43.10_EB,42.64_EB,42.19_EB/)
QFLAME(2,11,0:10,6) = (/0.00_EB,27.18_EB,38.39_EB,43.29_EB,48.19_EB,49.97_EB,51.75_EB,53.53_EB,55.30_EB,55.46_EB,55.62_EB/)
QFLAME(2,11,11:20,6) = (/55.78_EB,55.94_EB,55.54_EB,55.15_EB,54.75_EB,54.35_EB,54.74_EB,55.12_EB,55.50_EB,55.88_EB/)
QFLAME(2,11,21:30,6) = (/55.95_EB,56.02_EB,56.08_EB,56.15_EB,56.22_EB,55.68_EB,55.14_EB,54.60_EB,54.06_EB,53.52_EB/)
QFLAME(2,12,0:10,1) = (/0.00_EB,21.05_EB,28.12_EB,30.13_EB,32.15_EB,32.42_EB,32.70_EB,32.97_EB,33.25_EB,32.90_EB,32.56_EB/)
QFLAME(2,12,11:20,1) = (/32.22_EB,31.88_EB,31.26_EB,30.64_EB,30.01_EB,29.39_EB,28.95_EB,28.51_EB,28.06_EB,27.62_EB/)
QFLAME(2,12,21:30,1) = (/27.13_EB,26.63_EB,26.14_EB,25.64_EB,25.15_EB,24.71_EB,24.28_EB,23.84_EB,23.41_EB,22.97_EB/)
QFLAME(2,12,0:10,2) = (/0.00_EB,21.05_EB,28.11_EB,30.12_EB,32.13_EB,32.42_EB,32.70_EB,32.99_EB,33.27_EB,32.92_EB,32.58_EB/)
QFLAME(2,12,11:20,2) = (/32.23_EB,31.89_EB,31.26_EB,30.64_EB,30.02_EB,29.40_EB,28.95_EB,28.51_EB,28.06_EB,27.62_EB/)
QFLAME(2,12,21:30,2) = (/27.13_EB,26.64_EB,26.15_EB,25.65_EB,25.16_EB,24.72_EB,24.27_EB,23.83_EB,23.38_EB,22.94_EB/)
QFLAME(2,12,0:10,3) = (/0.00_EB,21.58_EB,28.32_EB,30.33_EB,32.34_EB,32.60_EB,32.86_EB,33.12_EB,33.38_EB,33.01_EB,32.65_EB/)
QFLAME(2,12,11:20,3) = (/32.28_EB,31.92_EB,31.30_EB,30.69_EB,30.08_EB,29.47_EB,29.08_EB,28.68_EB,28.29_EB,27.89_EB/)
QFLAME(2,12,21:30,3) = (/27.44_EB,26.98_EB,26.53_EB,26.07_EB,25.62_EB,25.22_EB,24.81_EB,24.41_EB,24.00_EB,23.60_EB/)
QFLAME(2,12,0:10,4) = (/0.00_EB,23.37_EB,30.96_EB,33.41_EB,35.86_EB,36.17_EB,36.48_EB,36.79_EB,37.10_EB,36.97_EB,36.83_EB/)
QFLAME(2,12,11:20,4) = (/36.70_EB,36.57_EB,36.15_EB,35.73_EB,35.31_EB,34.89_EB,35.10_EB,35.32_EB,35.53_EB,35.75_EB/)
QFLAME(2,12,21:30,4) = (/35.29_EB,34.82_EB,34.36_EB,33.90_EB,33.44_EB,33.00_EB,32.57_EB,32.14_EB,31.70_EB,31.27_EB/)
QFLAME(2,12,0:10,5) = (/0.00_EB,25.27_EB,34.59_EB,38.23_EB,41.86_EB,42.77_EB,43.67_EB,44.57_EB,45.47_EB,45.46_EB,45.44_EB/)
QFLAME(2,12,11:20,5) = (/45.42_EB,45.41_EB,45.00_EB,44.59_EB,44.18_EB,43.78_EB,44.49_EB,45.20_EB,45.92_EB,46.63_EB/)
QFLAME(2,12,21:30,5) = (/46.17_EB,45.70_EB,45.24_EB,44.78_EB,44.32_EB,43.86_EB,43.41_EB,42.95_EB,42.50_EB,42.04_EB/)
QFLAME(2,12,0:10,6) = (/0.00_EB,27.19_EB,38.30_EB,43.22_EB,48.14_EB,49.90_EB,51.66_EB,53.42_EB,55.18_EB,55.34_EB,55.50_EB/)
QFLAME(2,12,11:20,6) = (/55.65_EB,55.81_EB,55.42_EB,55.03_EB,54.64_EB,54.25_EB,54.67_EB,55.09_EB,55.51_EB,55.94_EB/)
QFLAME(2,12,21:30,6) = (/55.96_EB,55.98_EB,56.00_EB,56.03_EB,56.05_EB,55.51_EB,54.97_EB,54.43_EB,53.88_EB,53.34_EB/)
QFLAME(2,13,0:10,1) = (/0.00_EB,21.18_EB,28.35_EB,30.42_EB,32.50_EB,32.80_EB,33.10_EB,33.40_EB,33.70_EB,33.37_EB,33.03_EB/)
QFLAME(2,13,11:20,1) = (/32.70_EB,32.36_EB,31.74_EB,31.13_EB,30.51_EB,29.89_EB,29.44_EB,29.00_EB,28.55_EB,28.11_EB/)
QFLAME(2,13,21:30,1) = (/27.61_EB,27.12_EB,26.62_EB,26.13_EB,25.63_EB,25.19_EB,24.75_EB,24.31_EB,23.87_EB,23.43_EB/)
QFLAME(2,13,0:10,2) = (/0.00_EB,21.17_EB,28.34_EB,30.42_EB,32.49_EB,32.79_EB,33.10_EB,33.41_EB,33.72_EB,33.38_EB,33.04_EB/)
QFLAME(2,13,11:20,2) = (/32.70_EB,32.36_EB,31.75_EB,31.13_EB,30.51_EB,29.89_EB,29.45_EB,29.00_EB,28.55_EB,28.10_EB/)
QFLAME(2,13,21:30,2) = (/27.61_EB,27.12_EB,26.62_EB,26.13_EB,25.64_EB,25.19_EB,24.74_EB,24.29_EB,23.85_EB,23.40_EB/)
QFLAME(2,13,0:10,3) = (/0.00_EB,21.64_EB,28.53_EB,30.60_EB,32.67_EB,32.95_EB,33.24_EB,33.52_EB,33.81_EB,33.46_EB,33.10_EB/)
QFLAME(2,13,11:20,3) = (/32.75_EB,32.39_EB,31.78_EB,31.18_EB,30.57_EB,29.96_EB,29.56_EB,29.15_EB,28.75_EB,28.34_EB/)
QFLAME(2,13,21:30,3) = (/27.88_EB,27.42_EB,26.96_EB,26.50_EB,26.05_EB,25.63_EB,25.22_EB,24.81_EB,24.39_EB,23.98_EB/)
QFLAME(2,13,0:10,4) = (/0.00_EB,23.37_EB,30.93_EB,33.39_EB,35.85_EB,36.18_EB,36.51_EB,36.84_EB,37.17_EB,37.02_EB,36.87_EB/)
QFLAME(2,13,11:20,4) = (/36.72_EB,36.57_EB,36.14_EB,35.71_EB,35.27_EB,34.84_EB,35.03_EB,35.22_EB,35.41_EB,35.60_EB/)
QFLAME(2,13,21:30,4) = (/35.15_EB,34.69_EB,34.24_EB,33.78_EB,33.33_EB,32.90_EB,32.47_EB,32.04_EB,31.61_EB,31.18_EB/)
QFLAME(2,13,0:10,5) = (/0.00_EB,25.27_EB,34.55_EB,38.18_EB,41.81_EB,42.70_EB,43.60_EB,44.50_EB,45.40_EB,45.37_EB,45.35_EB/)
QFLAME(2,13,11:20,5) = (/45.32_EB,45.30_EB,44.90_EB,44.50_EB,44.11_EB,43.71_EB,44.40_EB,45.09_EB,45.78_EB,46.47_EB/)
QFLAME(2,13,21:30,5) = (/46.01_EB,45.55_EB,45.09_EB,44.63_EB,44.17_EB,43.72_EB,43.26_EB,42.80_EB,42.35_EB,41.89_EB/)
QFLAME(2,13,0:10,6) = (/0.00_EB,27.19_EB,38.22_EB,43.15_EB,48.08_EB,49.83_EB,51.57_EB,53.31_EB,55.05_EB,55.21_EB,55.37_EB/)
QFLAME(2,13,11:20,6) = (/55.53_EB,55.69_EB,55.31_EB,54.92_EB,54.53_EB,54.15_EB,54.61_EB,55.07_EB,55.53_EB,55.99_EB/)
QFLAME(2,13,21:30,6) = (/55.97_EB,55.95_EB,55.92_EB,55.90_EB,55.88_EB,55.34_EB,54.79_EB,54.25_EB,53.71_EB,53.16_EB/)
QFLAME(2,14,0:10,1) = (/0.00_EB,21.30_EB,28.58_EB,30.71_EB,32.85_EB,33.18_EB,33.50_EB,33.83_EB,34.16_EB,33.83_EB,33.50_EB/)
QFLAME(2,14,11:20,1) = (/33.17_EB,32.84_EB,32.23_EB,31.61_EB,31.00_EB,30.39_EB,29.94_EB,29.49_EB,29.04_EB,28.59_EB/)
QFLAME(2,14,21:30,1) = (/28.10_EB,27.60_EB,27.11_EB,26.61_EB,26.12_EB,25.67_EB,25.23_EB,24.78_EB,24.33_EB,23.89_EB/)
QFLAME(2,14,0:10,2) = (/0.00_EB,21.30_EB,28.58_EB,30.71_EB,32.84_EB,33.17_EB,33.50_EB,33.83_EB,34.17_EB,33.83_EB,33.50_EB/)
QFLAME(2,14,11:20,2) = (/33.17_EB,32.84_EB,32.23_EB,31.62_EB,31.01_EB,30.39_EB,29.94_EB,29.49_EB,29.04_EB,28.59_EB/)
QFLAME(2,14,21:30,2) = (/28.09_EB,27.60_EB,27.10_EB,26.61_EB,26.11_EB,25.66_EB,25.21_EB,24.76_EB,24.31_EB,23.86_EB/)
QFLAME(2,14,0:10,3) = (/0.00_EB,21.70_EB,28.73_EB,30.86_EB,32.99_EB,33.31_EB,33.62_EB,33.93_EB,34.25_EB,33.90_EB,33.55_EB/)
QFLAME(2,14,11:20,3) = (/33.21_EB,32.86_EB,32.26_EB,31.66_EB,31.06_EB,30.45_EB,30.04_EB,29.62_EB,29.20_EB,28.79_EB/)
QFLAME(2,14,21:30,3) = (/28.32_EB,27.86_EB,27.40_EB,26.93_EB,26.47_EB,26.05_EB,25.63_EB,25.20_EB,24.78_EB,24.36_EB/)
QFLAME(2,14,0:10,4) = (/0.00_EB,23.37_EB,30.90_EB,33.38_EB,35.85_EB,36.20_EB,36.54_EB,36.89_EB,37.24_EB,37.07_EB,36.91_EB/)
QFLAME(2,14,11:20,4) = (/36.74_EB,36.57_EB,36.13_EB,35.68_EB,35.24_EB,34.80_EB,34.96_EB,35.12_EB,35.29_EB,35.45_EB/)
QFLAME(2,14,21:30,4) = (/35.01_EB,34.56_EB,34.12_EB,33.67_EB,33.23_EB,32.80_EB,32.37_EB,31.95_EB,31.52_EB,31.09_EB/)
QFLAME(2,14,0:10,5) = (/0.00_EB,25.26_EB,34.51_EB,38.13_EB,41.75_EB,42.64_EB,43.54_EB,44.43_EB,45.33_EB,45.29_EB,45.26_EB/)
QFLAME(2,14,11:20,5) = (/45.22_EB,45.18_EB,44.80_EB,44.41_EB,44.03_EB,43.65_EB,44.31_EB,44.98_EB,45.65_EB,46.32_EB/)
QFLAME(2,14,21:30,5) = (/45.86_EB,45.40_EB,44.94_EB,44.48_EB,44.03_EB,43.57_EB,43.11_EB,42.66_EB,42.20_EB,41.75_EB/)
QFLAME(2,14,0:10,6) = (/0.00_EB,27.20_EB,38.13_EB,43.08_EB,48.03_EB,49.75_EB,51.48_EB,53.20_EB,54.93_EB,55.09_EB,55.25_EB/)
QFLAME(2,14,11:20,6) = (/55.41_EB,55.57_EB,55.19_EB,54.81_EB,54.43_EB,54.05_EB,54.55_EB,55.04_EB,55.54_EB,56.04_EB/)
QFLAME(2,14,21:30,6) = (/55.98_EB,55.91_EB,55.84_EB,55.78_EB,55.71_EB,55.17_EB,54.62_EB,54.08_EB,53.53_EB,52.99_EB/)
QFLAME(2,15,0:10,1) = (/0.00_EB,21.43_EB,28.81_EB,31.01_EB,33.20_EB,33.55_EB,33.91_EB,34.26_EB,34.61_EB,34.29_EB,33.97_EB/)
QFLAME(2,15,11:20,1) = (/33.64_EB,33.32_EB,32.71_EB,32.10_EB,31.49_EB,30.88_EB,30.43_EB,29.98_EB,29.53_EB,29.08_EB/)
QFLAME(2,15,21:30,1) = (/28.58_EB,28.09_EB,27.59_EB,27.09_EB,26.60_EB,26.15_EB,25.70_EB,25.25_EB,24.80_EB,24.35_EB/)
QFLAME(2,15,0:10,2) = (/0.00_EB,21.43_EB,28.81_EB,31.00_EB,33.19_EB,33.55_EB,33.90_EB,34.26_EB,34.61_EB,34.29_EB,33.96_EB/)
QFLAME(2,15,11:20,2) = (/33.64_EB,33.32_EB,32.71_EB,32.11_EB,31.50_EB,30.89_EB,30.44_EB,29.98_EB,29.53_EB,29.07_EB/)
QFLAME(2,15,21:30,2) = (/28.58_EB,28.08_EB,27.58_EB,27.08_EB,26.59_EB,26.13_EB,25.68_EB,25.22_EB,24.77_EB,24.32_EB/)
QFLAME(2,15,0:10,3) = (/0.00_EB,21.76_EB,28.94_EB,31.13_EB,33.32_EB,33.66_EB,34.00_EB,34.34_EB,34.68_EB,34.34_EB,34.01_EB/)
QFLAME(2,15,11:20,3) = (/33.67_EB,33.34_EB,32.74_EB,32.14_EB,31.54_EB,30.94_EB,30.52_EB,30.09_EB,29.66_EB,29.23_EB/)
QFLAME(2,15,21:30,3) = (/28.77_EB,28.30_EB,27.83_EB,27.36_EB,26.90_EB,26.47_EB,26.03_EB,25.60_EB,25.17_EB,24.73_EB/)
QFLAME(2,15,0:10,4) = (/0.00_EB,23.38_EB,30.87_EB,33.36_EB,35.84_EB,36.21_EB,36.58_EB,36.95_EB,37.31_EB,37.13_EB,36.94_EB/)
QFLAME(2,15,11:20,4) = (/36.76_EB,36.57_EB,36.12_EB,35.66_EB,35.21_EB,34.75_EB,34.89_EB,35.03_EB,35.16_EB,35.30_EB/)
QFLAME(2,15,21:30,4) = (/34.87_EB,34.43_EB,33.99_EB,33.56_EB,33.12_EB,32.70_EB,32.27_EB,31.85_EB,31.43_EB,31.00_EB/)
QFLAME(2,15,0:10,5) = (/0.00_EB,25.26_EB,34.46_EB,38.08_EB,41.69_EB,42.58_EB,43.47_EB,44.37_EB,45.26_EB,45.21_EB,45.16_EB/)
QFLAME(2,15,11:20,5) = (/45.12_EB,45.07_EB,44.70_EB,44.32_EB,43.95_EB,43.58_EB,44.23_EB,44.87_EB,45.52_EB,46.16_EB/)
QFLAME(2,15,21:30,5) = (/45.71_EB,45.25_EB,44.79_EB,44.34_EB,43.88_EB,43.42_EB,42.97_EB,42.51_EB,42.05_EB,41.60_EB/)
QFLAME(2,15,0:10,6) = (/0.00_EB,27.21_EB,38.05_EB,43.01_EB,47.97_EB,49.68_EB,51.39_EB,53.10_EB,54.80_EB,54.96_EB,55.13_EB/)
QFLAME(2,15,11:20,6) = (/55.29_EB,55.45_EB,55.07_EB,54.70_EB,54.32_EB,53.94_EB,54.48_EB,55.02_EB,55.56_EB,56.10_EB/)
QFLAME(2,15,21:30,6) = (/55.98_EB,55.87_EB,55.76_EB,55.65_EB,55.54_EB,54.99_EB,54.45_EB,53.90_EB,53.35_EB,52.81_EB/)
QFLAME(2,16,0:10,1) = (/0.00_EB,21.56_EB,29.05_EB,31.30_EB,33.55_EB,33.93_EB,34.31_EB,34.69_EB,35.07_EB,34.75_EB,34.43_EB/)
QFLAME(2,16,11:20,1) = (/34.12_EB,33.80_EB,33.20_EB,32.59_EB,31.99_EB,31.38_EB,30.93_EB,30.47_EB,30.02_EB,29.57_EB/)
QFLAME(2,16,21:30,1) = (/29.07_EB,28.57_EB,28.07_EB,27.58_EB,27.08_EB,26.63_EB,26.17_EB,25.72_EB,25.26_EB,24.81_EB/)
QFLAME(2,16,0:10,2) = (/0.00_EB,21.55_EB,29.05_EB,31.30_EB,33.55_EB,33.93_EB,34.30_EB,34.68_EB,35.06_EB,34.74_EB,34.43_EB/)
QFLAME(2,16,11:20,2) = (/34.11_EB,33.79_EB,33.19_EB,32.59_EB,31.99_EB,31.39_EB,30.94_EB,30.48_EB,30.02_EB,29.56_EB/)
QFLAME(2,16,21:30,2) = (/29.06_EB,28.56_EB,28.06_EB,27.56_EB,27.06_EB,26.60_EB,26.15_EB,25.69_EB,25.23_EB,24.78_EB/)
QFLAME(2,16,0:10,3) = (/0.00_EB,21.82_EB,29.15_EB,31.40_EB,33.65_EB,34.01_EB,34.38_EB,34.75_EB,35.11_EB,34.79_EB,34.46_EB/)
QFLAME(2,16,11:20,3) = (/34.14_EB,33.81_EB,33.22_EB,32.62_EB,32.03_EB,31.43_EB,30.99_EB,30.56_EB,30.12_EB,29.68_EB/)
QFLAME(2,16,21:30,3) = (/29.21_EB,28.74_EB,28.27_EB,27.80_EB,27.32_EB,26.88_EB,26.44_EB,26.00_EB,25.55_EB,25.11_EB/)
QFLAME(2,16,0:10,4) = (/0.00_EB,23.38_EB,30.85_EB,33.34_EB,35.84_EB,36.22_EB,36.61_EB,37.00_EB,37.38_EB,37.18_EB,36.98_EB/)
QFLAME(2,16,11:20,4) = (/36.77_EB,36.57_EB,36.10_EB,35.64_EB,35.17_EB,34.71_EB,34.82_EB,34.93_EB,35.04_EB,35.15_EB/)
QFLAME(2,16,21:30,4) = (/34.73_EB,34.30_EB,33.87_EB,33.44_EB,33.01_EB,32.59_EB,32.17_EB,31.75_EB,31.33_EB,30.91_EB/)
QFLAME(2,16,0:10,5) = (/0.00_EB,25.26_EB,34.42_EB,38.03_EB,41.63_EB,42.52_EB,43.41_EB,44.30_EB,45.19_EB,45.13_EB,45.07_EB/)
QFLAME(2,16,11:20,5) = (/45.01_EB,44.95_EB,44.59_EB,44.24_EB,43.88_EB,43.52_EB,44.14_EB,44.76_EB,45.38_EB,46.01_EB/)
QFLAME(2,16,21:30,5) = (/45.55_EB,45.10_EB,44.64_EB,44.19_EB,43.74_EB,43.28_EB,42.82_EB,42.36_EB,41.91_EB,41.45_EB/)
QFLAME(2,16,0:10,6) = (/0.00_EB,27.21_EB,37.96_EB,42.94_EB,47.92_EB,49.61_EB,51.30_EB,52.99_EB,54.68_EB,54.84_EB,55.00_EB/)
QFLAME(2,16,11:20,6) = (/55.16_EB,55.32_EB,54.95_EB,54.58_EB,54.21_EB,53.84_EB,54.42_EB,54.99_EB,55.57_EB,56.15_EB/)
QFLAME(2,16,21:30,6) = (/55.99_EB,55.84_EB,55.68_EB,55.53_EB,55.37_EB,54.82_EB,54.27_EB,53.73_EB,53.18_EB,52.63_EB/)
QFLAME(2,17,0:10,1) = (/0.00_EB,21.68_EB,29.28_EB,31.59_EB,33.90_EB,34.30_EB,34.71_EB,35.12_EB,35.52_EB,35.21_EB,34.90_EB/)
QFLAME(2,17,11:20,1) = (/34.59_EB,34.28_EB,33.68_EB,33.08_EB,32.48_EB,31.88_EB,31.42_EB,30.96_EB,30.51_EB,30.05_EB/)
QFLAME(2,17,21:30,1) = (/29.55_EB,29.06_EB,28.56_EB,28.06_EB,27.56_EB,27.10_EB,26.65_EB,26.19_EB,25.73_EB,25.27_EB/)
QFLAME(2,17,0:10,2) = (/0.00_EB,21.68_EB,29.28_EB,31.59_EB,33.90_EB,34.30_EB,34.70_EB,35.11_EB,35.51_EB,35.20_EB,34.89_EB/)
QFLAME(2,17,11:20,2) = (/34.58_EB,34.27_EB,33.67_EB,33.08_EB,32.49_EB,31.89_EB,31.43_EB,30.97_EB,30.51_EB,30.05_EB/)
QFLAME(2,17,21:30,2) = (/29.54_EB,29.04_EB,28.54_EB,28.04_EB,27.54_EB,27.08_EB,26.61_EB,26.15_EB,25.69_EB,25.23_EB/)
QFLAME(2,17,0:10,3) = (/0.00_EB,21.88_EB,29.35_EB,31.66_EB,33.97_EB,34.37_EB,34.76_EB,35.15_EB,35.55_EB,35.23_EB,34.92_EB/)
QFLAME(2,17,11:20,3) = (/34.60_EB,34.29_EB,33.70_EB,33.11_EB,32.52_EB,31.92_EB,31.47_EB,31.02_EB,30.57_EB,30.12_EB/)
QFLAME(2,17,21:30,3) = (/29.65_EB,29.17_EB,28.70_EB,28.23_EB,27.75_EB,27.30_EB,26.85_EB,26.39_EB,25.94_EB,25.49_EB/)
QFLAME(2,17,0:10,4) = (/0.00_EB,23.38_EB,30.82_EB,33.33_EB,35.83_EB,36.24_EB,36.64_EB,37.05_EB,37.46_EB,37.23_EB,37.01_EB/)
QFLAME(2,17,11:20,4) = (/36.79_EB,36.57_EB,36.09_EB,35.62_EB,35.14_EB,34.66_EB,34.75_EB,34.83_EB,34.92_EB,35.00_EB/)
QFLAME(2,17,21:30,4) = (/34.58_EB,34.17_EB,33.75_EB,33.33_EB,32.91_EB,32.49_EB,32.08_EB,31.66_EB,31.24_EB,30.82_EB/)
QFLAME(2,17,0:10,5) = (/0.00_EB,25.25_EB,34.38_EB,37.98_EB,41.58_EB,42.46_EB,43.35_EB,44.23_EB,45.12_EB,45.05_EB,44.98_EB/)
QFLAME(2,17,11:20,5) = (/44.91_EB,44.84_EB,44.49_EB,44.15_EB,43.80_EB,43.45_EB,44.05_EB,44.65_EB,45.25_EB,45.85_EB/)
QFLAME(2,17,21:30,5) = (/45.40_EB,44.95_EB,44.49_EB,44.04_EB,43.59_EB,43.13_EB,42.68_EB,42.22_EB,41.76_EB,41.30_EB/)
QFLAME(2,17,0:10,6) = (/0.00_EB,27.22_EB,37.88_EB,42.87_EB,47.86_EB,49.54_EB,51.21_EB,52.88_EB,54.55_EB,54.72_EB,54.88_EB/)
QFLAME(2,17,11:20,6) = (/55.04_EB,55.20_EB,54.84_EB,54.47_EB,54.10_EB,53.74_EB,54.35_EB,54.97_EB,55.59_EB,56.20_EB/)
QFLAME(2,17,21:30,6) = (/56.00_EB,55.80_EB,55.60_EB,55.40_EB,55.20_EB,54.65_EB,54.10_EB,53.55_EB,53.00_EB,52.45_EB/)
QFLAME(2,18,0:10,1) = (/0.00_EB,21.81_EB,29.51_EB,31.88_EB,34.25_EB,34.68_EB,35.11_EB,35.55_EB,35.98_EB,35.67_EB,35.37_EB/)
QFLAME(2,18,11:20,1) = (/35.07_EB,34.76_EB,34.16_EB,33.57_EB,32.97_EB,32.37_EB,31.91_EB,31.46_EB,31.00_EB,30.54_EB/)
QFLAME(2,18,21:30,1) = (/30.04_EB,29.54_EB,29.04_EB,28.54_EB,28.05_EB,27.58_EB,27.12_EB,26.65_EB,26.19_EB,25.73_EB/)
QFLAME(2,18,0:10,2) = (/0.00_EB,21.81_EB,29.51_EB,31.88_EB,34.25_EB,34.68_EB,35.11_EB,35.53_EB,35.96_EB,35.65_EB,35.35_EB/)
QFLAME(2,18,11:20,2) = (/35.05_EB,34.74_EB,34.16_EB,33.57_EB,32.98_EB,32.39_EB,31.93_EB,31.46_EB,31.00_EB,30.53_EB/)
QFLAME(2,18,21:30,2) = (/30.03_EB,29.52_EB,29.02_EB,28.51_EB,28.01_EB,27.55_EB,27.08_EB,26.62_EB,26.16_EB,25.69_EB/)
QFLAME(2,18,0:10,3) = (/0.00_EB,21.95_EB,29.56_EB,31.93_EB,34.30_EB,34.72_EB,35.14_EB,35.56_EB,35.98_EB,35.68_EB,35.37_EB/)
QFLAME(2,18,11:20,3) = (/35.07_EB,34.76_EB,34.18_EB,33.59_EB,33.00_EB,32.41_EB,31.95_EB,31.49_EB,31.03_EB,30.57_EB/)
QFLAME(2,18,21:30,3) = (/30.09_EB,29.61_EB,29.13_EB,28.66_EB,28.18_EB,27.72_EB,27.25_EB,26.79_EB,26.33_EB,25.87_EB/)
QFLAME(2,18,0:10,4) = (/0.00_EB,23.38_EB,30.79_EB,33.31_EB,35.83_EB,36.25_EB,36.68_EB,37.10_EB,37.53_EB,37.29_EB,37.05_EB/)
QFLAME(2,18,11:20,4) = (/36.81_EB,36.57_EB,36.08_EB,35.59_EB,35.11_EB,34.62_EB,34.68_EB,34.74_EB,34.80_EB,34.86_EB/)
QFLAME(2,18,21:30,4) = (/34.44_EB,34.03_EB,33.62_EB,33.21_EB,32.80_EB,32.39_EB,31.98_EB,31.56_EB,31.15_EB,30.74_EB/)
QFLAME(2,18,0:10,5) = (/0.00_EB,25.25_EB,34.34_EB,37.93_EB,41.52_EB,42.40_EB,43.28_EB,44.16_EB,45.05_EB,44.97_EB,44.89_EB/)
QFLAME(2,18,11:20,5) = (/44.81_EB,44.73_EB,44.39_EB,44.06_EB,43.72_EB,43.39_EB,43.96_EB,44.54_EB,45.12_EB,45.69_EB/)
QFLAME(2,18,21:30,5) = (/45.24_EB,44.79_EB,44.34_EB,43.89_EB,43.45_EB,42.99_EB,42.53_EB,42.07_EB,41.61_EB,41.16_EB/)
QFLAME(2,18,0:10,6) = (/0.00_EB,27.23_EB,37.79_EB,42.80_EB,47.81_EB,49.46_EB,51.12_EB,52.77_EB,54.43_EB,54.59_EB,54.75_EB/)
QFLAME(2,18,11:20,6) = (/54.92_EB,55.08_EB,54.72_EB,54.36_EB,54.00_EB,53.64_EB,54.29_EB,54.95_EB,55.60_EB,56.25_EB/)
QFLAME(2,18,21:30,6) = (/56.01_EB,55.76_EB,55.52_EB,55.28_EB,55.03_EB,54.48_EB,53.93_EB,53.38_EB,52.83_EB,52.27_EB/)
QFLAME(2,19,0:10,1) = (/0.00_EB,21.93_EB,29.74_EB,32.17_EB,34.60_EB,35.06_EB,35.52_EB,35.97_EB,36.43_EB,36.14_EB,35.84_EB/)
QFLAME(2,19,11:20,1) = (/35.54_EB,35.24_EB,34.65_EB,34.06_EB,33.46_EB,32.87_EB,32.41_EB,31.95_EB,31.49_EB,31.03_EB/)
QFLAME(2,19,21:30,1) = (/30.53_EB,30.03_EB,29.53_EB,29.03_EB,28.53_EB,28.06_EB,27.59_EB,27.12_EB,26.66_EB,26.19_EB/)
QFLAME(2,19,0:10,2) = (/0.00_EB,21.93_EB,29.75_EB,32.18_EB,34.61_EB,35.06_EB,35.51_EB,35.96_EB,36.41_EB,36.11_EB,35.81_EB/)
QFLAME(2,19,11:20,2) = (/35.52_EB,35.22_EB,34.64_EB,34.06_EB,33.47_EB,32.89_EB,32.42_EB,31.95_EB,31.49_EB,31.02_EB/)
QFLAME(2,19,21:30,2) = (/30.51_EB,30.00_EB,29.50_EB,28.99_EB,28.48_EB,28.02_EB,27.55_EB,27.08_EB,26.62_EB,26.15_EB/)
QFLAME(2,19,0:10,3) = (/0.00_EB,22.01_EB,29.77_EB,32.20_EB,34.63_EB,35.08_EB,35.52_EB,35.97_EB,36.42_EB,36.12_EB,35.83_EB/)
QFLAME(2,19,11:20,3) = (/35.53_EB,35.24_EB,34.65_EB,34.07_EB,33.49_EB,32.90_EB,32.43_EB,31.96_EB,31.49_EB,31.02_EB/)
QFLAME(2,19,21:30,3) = (/30.53_EB,30.05_EB,29.57_EB,29.09_EB,28.60_EB,28.13_EB,27.66_EB,27.19_EB,26.72_EB,26.25_EB/)
QFLAME(2,19,0:10,4) = (/0.00_EB,23.38_EB,30.76_EB,33.29_EB,35.82_EB,36.27_EB,36.71_EB,37.16_EB,37.60_EB,37.34_EB,37.08_EB/)
QFLAME(2,19,11:20,4) = (/36.83_EB,36.57_EB,36.07_EB,35.57_EB,35.07_EB,34.57_EB,34.61_EB,34.64_EB,34.67_EB,34.71_EB/)
QFLAME(2,19,21:30,4) = (/34.30_EB,33.90_EB,33.50_EB,33.10_EB,32.70_EB,32.29_EB,31.88_EB,31.47_EB,31.06_EB,30.65_EB/)
QFLAME(2,19,0:10,5) = (/0.00_EB,25.25_EB,34.29_EB,37.88_EB,41.46_EB,42.34_EB,43.22_EB,44.10_EB,44.97_EB,44.88_EB,44.79_EB/)
QFLAME(2,19,11:20,5) = (/44.70_EB,44.61_EB,44.29_EB,43.97_EB,43.64_EB,43.32_EB,43.87_EB,44.43_EB,44.98_EB,45.54_EB/)
QFLAME(2,19,21:30,5) = (/45.09_EB,44.64_EB,44.20_EB,43.75_EB,43.30_EB,42.84_EB,42.38_EB,41.93_EB,41.47_EB,41.01_EB/)
QFLAME(2,19,0:10,6) = (/0.00_EB,27.23_EB,37.71_EB,42.73_EB,47.75_EB,49.39_EB,51.03_EB,52.67_EB,54.31_EB,54.47_EB,54.63_EB/)
QFLAME(2,19,11:20,6) = (/54.79_EB,54.96_EB,54.60_EB,54.24_EB,53.89_EB,53.53_EB,54.23_EB,54.92_EB,55.61_EB,56.31_EB/)
QFLAME(2,19,21:30,6) = (/56.02_EB,55.73_EB,55.44_EB,55.15_EB,54.86_EB,54.31_EB,53.76_EB,53.20_EB,52.65_EB,52.10_EB/)
QFLAME(2,20,0:10,1) = (/0.00_EB,22.06_EB,29.98_EB,32.46_EB,34.95_EB,35.43_EB,35.92_EB,36.40_EB,36.89_EB,36.60_EB,36.31_EB/)
QFLAME(2,20,11:20,1) = (/36.01_EB,35.72_EB,35.13_EB,34.54_EB,33.96_EB,33.37_EB,32.90_EB,32.44_EB,31.98_EB,31.51_EB/)
QFLAME(2,20,21:30,1) = (/31.01_EB,30.51_EB,30.01_EB,29.51_EB,29.01_EB,28.54_EB,28.07_EB,27.59_EB,27.12_EB,26.65_EB/)
QFLAME(2,20,0:10,2) = (/0.00_EB,22.06_EB,29.98_EB,32.47_EB,34.96_EB,35.43_EB,35.91_EB,36.38_EB,36.85_EB,36.56_EB,36.28_EB/)
QFLAME(2,20,11:20,2) = (/35.99_EB,35.70_EB,35.12_EB,34.54_EB,33.97_EB,33.39_EB,32.92_EB,32.45_EB,31.97_EB,31.50_EB/)
QFLAME(2,20,21:30,2) = (/30.99_EB,30.48_EB,29.98_EB,29.47_EB,28.96_EB,28.49_EB,28.02_EB,27.55_EB,27.08_EB,26.61_EB/)
QFLAME(2,20,0:10,3) = (/0.00_EB,22.07_EB,29.97_EB,32.46_EB,34.96_EB,35.43_EB,35.90_EB,36.38_EB,36.85_EB,36.57_EB,36.28_EB/)
QFLAME(2,20,11:20,3) = (/36.00_EB,35.71_EB,35.13_EB,34.55_EB,33.97_EB,33.40_EB,32.91_EB,32.43_EB,31.95_EB,31.46_EB/)
QFLAME(2,20,21:30,3) = (/30.98_EB,30.49_EB,30.00_EB,29.52_EB,29.03_EB,28.55_EB,28.07_EB,27.59_EB,27.10_EB,26.62_EB/)
QFLAME(2,20,0:10,4) = (/0.00_EB,23.38_EB,30.73_EB,33.27_EB,35.82_EB,36.28_EB,36.75_EB,37.21_EB,37.67_EB,37.40_EB,37.12_EB/)
QFLAME(2,20,11:20,4) = (/36.84_EB,36.57_EB,36.06_EB,35.55_EB,35.04_EB,34.53_EB,34.54_EB,34.54_EB,34.55_EB,34.56_EB/)
QFLAME(2,20,21:30,4) = (/34.16_EB,33.77_EB,33.38_EB,32.99_EB,32.59_EB,32.19_EB,31.78_EB,31.37_EB,30.96_EB,30.56_EB/)
QFLAME(2,20,0:10,5) = (/0.00_EB,25.24_EB,34.25_EB,37.83_EB,41.40_EB,42.28_EB,43.15_EB,44.03_EB,44.90_EB,44.80_EB,44.70_EB/)
QFLAME(2,20,11:20,5) = (/44.60_EB,44.50_EB,44.19_EB,43.88_EB,43.57_EB,43.26_EB,43.79_EB,44.32_EB,44.85_EB,45.38_EB/)
QFLAME(2,20,21:30,5) = (/44.94_EB,44.49_EB,44.05_EB,43.60_EB,43.16_EB,42.70_EB,42.24_EB,41.78_EB,41.32_EB,40.86_EB/)
QFLAME(2,20,0:10,6) = (/0.00_EB,27.24_EB,37.62_EB,42.66_EB,47.70_EB,49.32_EB,50.94_EB,52.56_EB,54.18_EB,54.34_EB,54.51_EB/)
QFLAME(2,20,11:20,6) = (/54.67_EB,54.83_EB,54.48_EB,54.13_EB,53.78_EB,53.43_EB,54.16_EB,54.90_EB,55.63_EB,56.36_EB/)
QFLAME(2,20,21:30,6) = (/56.03_EB,55.69_EB,55.36_EB,55.03_EB,54.69_EB,54.14_EB,53.58_EB,53.03_EB,52.47_EB,51.92_EB/)
QFLAME(3,0,0:10,1) = (/0.00_EB,18.59_EB,23.26_EB,23.77_EB,24.28_EB,23.96_EB,23.64_EB,23.32_EB,23.01_EB,22.33_EB,21.66_EB/)
QFLAME(3,0,11:20,1) = (/20.98_EB,20.31_EB,19.80_EB,19.30_EB,18.79_EB,18.28_EB,17.79_EB,17.31_EB,16.82_EB,16.33_EB/)
QFLAME(3,0,21:30,1) = (/15.96_EB,15.58_EB,15.21_EB,14.83_EB,14.46_EB,14.16_EB,13.86_EB,13.56_EB,13.26_EB,12.96_EB/)
QFLAME(3,0,0:10,2) = (/0.00_EB,20.48_EB,25.93_EB,26.83_EB,27.72_EB,27.55_EB,27.38_EB,27.21_EB,27.04_EB,26.45_EB,25.86_EB/)
QFLAME(3,0,11:20,2) = (/25.28_EB,24.69_EB,24.26_EB,23.84_EB,23.41_EB,22.99_EB,22.54_EB,22.10_EB,21.65_EB,21.21_EB/)
QFLAME(3,0,21:30,2) = (/20.99_EB,20.78_EB,20.56_EB,20.35_EB,20.14_EB,19.94_EB,19.74_EB,19.54_EB,19.34_EB,19.14_EB/)
QFLAME(3,0,0:10,3) = (/0.00_EB,22.33_EB,29.33_EB,31.03_EB,32.72_EB,33.10_EB,33.47_EB,33.85_EB,34.23_EB,33.73_EB,33.24_EB/)
QFLAME(3,0,11:20,3) = (/32.75_EB,32.25_EB,32.05_EB,31.85_EB,31.65_EB,31.45_EB,30.90_EB,30.35_EB,29.79_EB,29.24_EB/)
QFLAME(3,0,21:30,3) = (/29.46_EB,29.68_EB,29.91_EB,30.13_EB,30.35_EB,29.92_EB,29.49_EB,29.06_EB,28.63_EB,28.20_EB/)
QFLAME(3,0,0:10,4) = (/0.00_EB,24.30_EB,32.84_EB,35.63_EB,38.42_EB,39.38_EB,40.33_EB,41.29_EB,42.24_EB,42.07_EB,41.91_EB/)
QFLAME(3,0,11:20,4) = (/41.74_EB,41.58_EB,41.48_EB,41.39_EB,41.29_EB,41.20_EB,40.67_EB,40.14_EB,39.61_EB,39.08_EB/)
QFLAME(3,0,21:30,4) = (/39.54_EB,40.01_EB,40.48_EB,40.94_EB,41.41_EB,40.99_EB,40.57_EB,40.15_EB,39.73_EB,39.31_EB/)
QFLAME(3,0,0:10,5) = (/0.00_EB,26.06_EB,36.55_EB,40.65_EB,44.75_EB,46.38_EB,48.02_EB,49.65_EB,51.29_EB,51.46_EB,51.63_EB/)
QFLAME(3,0,11:20,5) = (/51.80_EB,51.97_EB,51.77_EB,51.56_EB,51.35_EB,51.15_EB,51.06_EB,50.96_EB,50.87_EB,50.78_EB/)
QFLAME(3,0,21:30,5) = (/50.66_EB,50.54_EB,50.42_EB,50.30_EB,50.17_EB,50.42_EB,50.67_EB,50.92_EB,51.17_EB,51.42_EB/)
QFLAME(3,0,0:10,6) = (/0.00_EB,27.95_EB,40.30_EB,45.79_EB,51.27_EB,53.50_EB,55.74_EB,57.97_EB,60.21_EB,60.82_EB,61.44_EB/)
QFLAME(3,0,11:20,6) = (/62.06_EB,62.68_EB,62.96_EB,63.24_EB,63.52_EB,63.80_EB,63.69_EB,63.57_EB,63.46_EB,63.35_EB/)
QFLAME(3,0,21:30,6) = (/62.92_EB,62.49_EB,62.06_EB,61.63_EB,61.20_EB,61.97_EB,62.74_EB,63.51_EB,64.28_EB,65.05_EB/)
QFLAME(3,1,0:10,1) = (/0.00_EB,18.96_EB,24.10_EB,24.89_EB,25.68_EB,25.49_EB,25.31_EB,25.12_EB,24.94_EB,24.29_EB,23.64_EB/)
QFLAME(3,1,11:20,1) = (/22.99_EB,22.34_EB,21.85_EB,21.37_EB,20.89_EB,20.40_EB,19.92_EB,19.43_EB,18.94_EB,18.46_EB/)
QFLAME(3,1,21:30,1) = (/18.07_EB,17.68_EB,17.30_EB,16.91_EB,16.52_EB,16.21_EB,15.90_EB,15.59_EB,15.29_EB,14.98_EB/)
QFLAME(3,1,0:10,2) = (/0.00_EB,20.46_EB,25.95_EB,26.84_EB,27.72_EB,27.55_EB,27.37_EB,27.20_EB,27.02_EB,26.37_EB,25.72_EB/)
QFLAME(3,1,11:20,2) = (/25.07_EB,24.42_EB,24.04_EB,23.66_EB,23.27_EB,22.89_EB,22.34_EB,21.79_EB,21.24_EB,20.69_EB/)
QFLAME(3,1,21:30,2) = (/20.50_EB,20.32_EB,20.13_EB,19.94_EB,19.75_EB,19.51_EB,19.27_EB,19.03_EB,18.79_EB,18.55_EB/)
QFLAME(3,1,0:10,3) = (/0.00_EB,22.33_EB,29.33_EB,31.04_EB,32.74_EB,33.09_EB,33.44_EB,33.80_EB,34.15_EB,33.66_EB,33.17_EB/)
QFLAME(3,1,11:20,3) = (/32.67_EB,32.18_EB,31.99_EB,31.81_EB,31.62_EB,31.43_EB,30.83_EB,30.23_EB,29.63_EB,29.03_EB/)
QFLAME(3,1,21:30,3) = (/28.90_EB,28.77_EB,28.64_EB,28.51_EB,28.38_EB,28.27_EB,28.15_EB,28.04_EB,27.92_EB,27.81_EB/)
QFLAME(3,1,0:10,4) = (/0.00_EB,24.24_EB,32.84_EB,35.65_EB,38.47_EB,39.41_EB,40.35_EB,41.29_EB,42.23_EB,42.02_EB,41.82_EB/)
QFLAME(3,1,11:20,4) = (/41.61_EB,41.41_EB,41.28_EB,41.15_EB,41.03_EB,40.90_EB,40.44_EB,39.99_EB,39.53_EB,39.07_EB/)
QFLAME(3,1,21:30,4) = (/39.54_EB,40.01_EB,40.48_EB,40.95_EB,41.42_EB,40.97_EB,40.52_EB,40.07_EB,39.62_EB,39.16_EB/)
QFLAME(3,1,0:10,5) = (/0.00_EB,26.10_EB,36.53_EB,40.65_EB,44.78_EB,46.38_EB,47.99_EB,49.59_EB,51.20_EB,51.37_EB,51.55_EB/)
QFLAME(3,1,11:20,5) = (/51.72_EB,51.90_EB,51.61_EB,51.33_EB,51.05_EB,50.77_EB,50.74_EB,50.71_EB,50.68_EB,50.65_EB/)
QFLAME(3,1,21:30,5) = (/51.18_EB,51.71_EB,52.24_EB,52.77_EB,53.30_EB,52.97_EB,52.64_EB,52.31_EB,51.98_EB,51.65_EB/)
QFLAME(3,1,0:10,6) = (/0.00_EB,27.86_EB,40.23_EB,45.69_EB,51.15_EB,53.41_EB,55.67_EB,57.93_EB,60.19_EB,60.72_EB,61.24_EB/)
QFLAME(3,1,11:20,6) = (/61.76_EB,62.28_EB,62.61_EB,62.95_EB,63.28_EB,63.62_EB,63.46_EB,63.31_EB,63.16_EB,63.01_EB/)
QFLAME(3,1,21:30,6) = (/62.44_EB,61.88_EB,61.31_EB,60.75_EB,60.18_EB,61.02_EB,61.85_EB,62.69_EB,63.52_EB,64.35_EB/)
QFLAME(3,2,0:10,1) = (/0.00_EB,19.34_EB,24.83_EB,25.84_EB,26.85_EB,26.78_EB,26.70_EB,26.63_EB,26.56_EB,25.95_EB,25.33_EB/)
QFLAME(3,2,11:20,1) = (/24.72_EB,24.10_EB,23.63_EB,23.16_EB,22.69_EB,22.22_EB,21.75_EB,21.28_EB,20.80_EB,20.33_EB/)
QFLAME(3,2,21:30,1) = (/19.96_EB,19.58_EB,19.21_EB,18.83_EB,18.46_EB,18.10_EB,17.75_EB,17.40_EB,17.05_EB,16.70_EB/)
QFLAME(3,2,0:10,2) = (/0.00_EB,20.47_EB,25.94_EB,26.88_EB,27.81_EB,27.70_EB,27.59_EB,27.48_EB,27.36_EB,26.69_EB,26.02_EB/)
QFLAME(3,2,11:20,2) = (/25.35_EB,24.68_EB,24.27_EB,23.87_EB,23.47_EB,23.07_EB,22.61_EB,22.15_EB,21.69_EB,21.23_EB/)
QFLAME(3,2,21:30,2) = (/20.98_EB,20.74_EB,20.50_EB,20.26_EB,20.01_EB,19.72_EB,19.43_EB,19.14_EB,18.85_EB,18.57_EB/)
QFLAME(3,2,0:10,3) = (/0.00_EB,22.38_EB,29.34_EB,31.04_EB,32.74_EB,33.08_EB,33.43_EB,33.77_EB,34.11_EB,33.60_EB,33.09_EB/)
QFLAME(3,2,11:20,3) = (/32.57_EB,32.06_EB,31.83_EB,31.61_EB,31.39_EB,31.16_EB,30.65_EB,30.15_EB,29.64_EB,29.13_EB/)
QFLAME(3,2,21:30,3) = (/29.21_EB,29.29_EB,29.37_EB,29.46_EB,29.54_EB,29.18_EB,28.82_EB,28.46_EB,28.10_EB,27.74_EB/)
QFLAME(3,2,0:10,4) = (/0.00_EB,24.23_EB,32.83_EB,35.65_EB,38.46_EB,39.39_EB,40.32_EB,41.25_EB,42.18_EB,41.95_EB,41.71_EB/)
QFLAME(3,2,11:20,4) = (/41.47_EB,41.24_EB,41.12_EB,41.00_EB,40.88_EB,40.76_EB,40.27_EB,39.79_EB,39.31_EB,38.83_EB/)
QFLAME(3,2,21:30,4) = (/39.15_EB,39.47_EB,39.79_EB,40.12_EB,40.44_EB,40.16_EB,39.88_EB,39.60_EB,39.32_EB,39.04_EB/)
QFLAME(3,2,0:10,5) = (/0.00_EB,26.04_EB,36.50_EB,40.66_EB,44.82_EB,46.39_EB,47.95_EB,49.51_EB,51.08_EB,51.24_EB,51.39_EB/)
QFLAME(3,2,11:20,5) = (/51.55_EB,51.71_EB,51.48_EB,51.24_EB,51.01_EB,50.77_EB,50.74_EB,50.70_EB,50.66_EB,50.62_EB/)
QFLAME(3,2,21:30,5) = (/50.25_EB,49.88_EB,49.52_EB,49.15_EB,48.78_EB,49.19_EB,49.59_EB,49.99_EB,50.39_EB,50.80_EB/)
QFLAME(3,2,0:10,6) = (/0.00_EB,27.87_EB,40.11_EB,45.59_EB,51.06_EB,53.32_EB,55.58_EB,57.84_EB,60.09_EB,60.55_EB,61.01_EB/)
QFLAME(3,2,11:20,6) = (/61.47_EB,61.93_EB,62.24_EB,62.56_EB,62.87_EB,63.19_EB,63.06_EB,62.92_EB,62.79_EB,62.66_EB/)
QFLAME(3,2,21:30,6) = (/62.10_EB,61.54_EB,60.98_EB,60.42_EB,59.86_EB,60.70_EB,61.53_EB,62.36_EB,63.20_EB,64.03_EB/)
QFLAME(3,3,0:10,1) = (/0.00_EB,19.62_EB,25.39_EB,26.58_EB,27.76_EB,27.78_EB,27.80_EB,27.81_EB,27.83_EB,27.24_EB,26.64_EB/)
QFLAME(3,3,11:20,1) = (/26.05_EB,25.46_EB,25.00_EB,24.55_EB,24.09_EB,23.64_EB,23.17_EB,22.69_EB,22.22_EB,21.75_EB/)
QFLAME(3,3,21:30,1) = (/21.36_EB,20.96_EB,20.57_EB,20.18_EB,19.79_EB,19.43_EB,19.07_EB,18.71_EB,18.36_EB,18.00_EB/)
QFLAME(3,3,0:10,2) = (/0.00_EB,20.46_EB,26.15_EB,27.28_EB,28.41_EB,28.41_EB,28.40_EB,28.39_EB,28.39_EB,27.75_EB,27.11_EB/)
QFLAME(3,3,11:20,2) = (/26.48_EB,25.84_EB,25.43_EB,25.02_EB,24.62_EB,24.21_EB,23.74_EB,23.28_EB,22.82_EB,22.35_EB/)
QFLAME(3,3,21:30,2) = (/22.06_EB,21.76_EB,21.46_EB,21.16_EB,20.86_EB,20.55_EB,20.24_EB,19.92_EB,19.61_EB,19.29_EB/)
QFLAME(3,3,0:10,3) = (/0.00_EB,22.37_EB,29.32_EB,31.02_EB,32.71_EB,33.05_EB,33.38_EB,33.72_EB,34.05_EB,33.54_EB,33.02_EB/)
QFLAME(3,3,11:20,3) = (/32.51_EB,31.99_EB,31.76_EB,31.54_EB,31.31_EB,31.08_EB,30.54_EB,30.00_EB,29.46_EB,28.92_EB/)
QFLAME(3,3,21:30,3) = (/29.01_EB,29.10_EB,29.19_EB,29.28_EB,29.38_EB,29.05_EB,28.73_EB,28.41_EB,28.09_EB,27.76_EB/)
QFLAME(3,3,0:10,4) = (/0.00_EB,24.23_EB,32.81_EB,35.63_EB,38.45_EB,39.36_EB,40.28_EB,41.19_EB,42.11_EB,41.86_EB,41.62_EB/)
QFLAME(3,3,11:20,4) = (/41.37_EB,41.13_EB,41.02_EB,40.91_EB,40.81_EB,40.70_EB,40.17_EB,39.64_EB,39.11_EB,38.58_EB/)
QFLAME(3,3,21:30,4) = (/38.93_EB,39.28_EB,39.62_EB,39.97_EB,40.31_EB,39.98_EB,39.65_EB,39.32_EB,38.98_EB,38.65_EB/)
QFLAME(3,3,0:10,5) = (/0.00_EB,26.05_EB,36.48_EB,40.64_EB,44.81_EB,46.35_EB,47.89_EB,49.43_EB,50.98_EB,51.13_EB,51.28_EB/)
QFLAME(3,3,11:20,5) = (/51.43_EB,51.58_EB,51.32_EB,51.06_EB,50.80_EB,50.54_EB,50.47_EB,50.39_EB,50.32_EB,50.24_EB/)
QFLAME(3,3,21:30,5) = (/50.20_EB,50.17_EB,50.13_EB,50.09_EB,50.05_EB,50.18_EB,50.31_EB,50.44_EB,50.57_EB,50.70_EB/)
QFLAME(3,3,0:10,6) = (/0.00_EB,27.86_EB,40.12_EB,45.56_EB,51.01_EB,53.26_EB,55.52_EB,57.78_EB,60.04_EB,60.54_EB,61.04_EB/)
QFLAME(3,3,11:20,6) = (/61.54_EB,62.04_EB,62.22_EB,62.40_EB,62.59_EB,62.77_EB,62.70_EB,62.64_EB,62.57_EB,62.51_EB/)
QFLAME(3,3,21:30,6) = (/61.99_EB,61.47_EB,60.96_EB,60.44_EB,59.93_EB,60.64_EB,61.35_EB,62.07_EB,62.78_EB,63.50_EB/)
QFLAME(3,4,0:10,1) = (/0.00_EB,19.90_EB,25.96_EB,27.31_EB,28.67_EB,28.78_EB,28.89_EB,28.99_EB,29.10_EB,28.53_EB,27.96_EB/)
QFLAME(3,4,11:20,1) = (/27.39_EB,26.82_EB,26.38_EB,25.94_EB,25.50_EB,25.06_EB,24.59_EB,24.11_EB,23.64_EB,23.16_EB/)
QFLAME(3,4,21:30,1) = (/22.76_EB,22.35_EB,21.94_EB,21.53_EB,21.12_EB,20.75_EB,20.39_EB,20.02_EB,19.66_EB,19.30_EB/)
QFLAME(3,4,0:10,2) = (/0.00_EB,20.45_EB,26.36_EB,27.69_EB,29.01_EB,29.11_EB,29.21_EB,29.31_EB,29.41_EB,28.81_EB,28.21_EB/)
QFLAME(3,4,11:20,2) = (/27.60_EB,27.00_EB,26.59_EB,26.18_EB,25.76_EB,25.35_EB,24.88_EB,24.42_EB,23.95_EB,23.48_EB/)
QFLAME(3,4,21:30,2) = (/23.13_EB,22.77_EB,22.42_EB,22.07_EB,21.71_EB,21.37_EB,21.04_EB,20.70_EB,20.36_EB,20.02_EB/)
QFLAME(3,4,0:10,3) = (/0.00_EB,22.36_EB,29.31_EB,31.00_EB,32.69_EB,33.02_EB,33.34_EB,33.67_EB,34.00_EB,33.48_EB,32.96_EB/)
QFLAME(3,4,11:20,3) = (/32.44_EB,31.93_EB,31.70_EB,31.47_EB,31.23_EB,31.00_EB,30.43_EB,29.86_EB,29.28_EB,28.71_EB/)
QFLAME(3,4,21:30,3) = (/28.81_EB,28.91_EB,29.01_EB,29.11_EB,29.22_EB,28.93_EB,28.65_EB,28.36_EB,28.07_EB,27.79_EB/)
QFLAME(3,4,0:10,4) = (/0.00_EB,24.23_EB,32.79_EB,35.61_EB,38.44_EB,39.33_EB,40.23_EB,41.13_EB,42.03_EB,41.78_EB,41.52_EB/)
QFLAME(3,4,11:20,4) = (/41.27_EB,41.01_EB,40.92_EB,40.83_EB,40.74_EB,40.65_EB,40.07_EB,39.49_EB,38.92_EB,38.34_EB/)
QFLAME(3,4,21:30,4) = (/38.71_EB,39.08_EB,39.45_EB,39.82_EB,40.19_EB,39.81_EB,39.42_EB,39.03_EB,38.65_EB,38.26_EB/)
QFLAME(3,4,0:10,5) = (/0.00_EB,26.06_EB,36.45_EB,40.62_EB,44.79_EB,46.31_EB,47.84_EB,49.36_EB,50.88_EB,51.02_EB,51.17_EB/)
QFLAME(3,4,11:20,5) = (/51.31_EB,51.45_EB,51.17_EB,50.88_EB,50.60_EB,50.31_EB,50.20_EB,50.09_EB,49.97_EB,49.86_EB/)
QFLAME(3,4,21:30,5) = (/50.15_EB,50.45_EB,50.74_EB,51.03_EB,51.32_EB,51.18_EB,51.04_EB,50.89_EB,50.75_EB,50.61_EB/)
QFLAME(3,4,0:10,6) = (/0.00_EB,27.85_EB,40.13_EB,45.54_EB,50.96_EB,53.21_EB,55.47_EB,57.72_EB,59.98_EB,60.52_EB,61.07_EB/)
QFLAME(3,4,11:20,6) = (/61.61_EB,62.16_EB,62.20_EB,62.25_EB,62.30_EB,62.34_EB,62.34_EB,62.35_EB,62.35_EB,62.35_EB/)
QFLAME(3,4,21:30,6) = (/61.88_EB,61.41_EB,60.94_EB,60.46_EB,59.99_EB,60.58_EB,61.18_EB,61.77_EB,62.37_EB,62.96_EB/)
QFLAME(3,5,0:10,1) = (/0.00_EB,20.18_EB,26.52_EB,28.05_EB,29.58_EB,29.78_EB,29.98_EB,30.17_EB,30.37_EB,29.82_EB,29.27_EB/)
QFLAME(3,5,11:20,1) = (/28.72_EB,28.17_EB,27.75_EB,27.33_EB,26.90_EB,26.48_EB,26.00_EB,25.53_EB,25.06_EB,24.58_EB/)
QFLAME(3,5,21:30,1) = (/24.16_EB,23.73_EB,23.30_EB,22.87_EB,22.45_EB,22.08_EB,21.71_EB,21.33_EB,20.96_EB,20.59_EB/)
QFLAME(3,5,0:10,2) = (/0.00_EB,20.44_EB,26.56_EB,28.09_EB,29.62_EB,29.82_EB,30.02_EB,30.23_EB,30.43_EB,29.86_EB,29.30_EB/)
QFLAME(3,5,11:20,2) = (/28.73_EB,28.17_EB,27.75_EB,27.33_EB,26.91_EB,26.49_EB,26.02_EB,25.55_EB,25.08_EB,24.61_EB/)
QFLAME(3,5,21:30,2) = (/24.20_EB,23.79_EB,23.38_EB,22.97_EB,22.56_EB,22.20_EB,21.84_EB,21.48_EB,21.11_EB,20.75_EB/)
QFLAME(3,5,0:10,3) = (/0.00_EB,22.35_EB,29.30_EB,30.98_EB,32.67_EB,32.98_EB,33.30_EB,33.62_EB,33.94_EB,33.42_EB,32.90_EB/)
QFLAME(3,5,11:20,3) = (/32.38_EB,31.86_EB,31.63_EB,31.39_EB,31.16_EB,30.92_EB,30.32_EB,29.71_EB,29.10_EB,28.50_EB/)
QFLAME(3,5,21:30,3) = (/28.61_EB,28.72_EB,28.83_EB,28.94_EB,29.06_EB,28.81_EB,28.56_EB,28.31_EB,28.06_EB,27.81_EB/)
QFLAME(3,5,0:10,4) = (/0.00_EB,24.23_EB,32.77_EB,35.60_EB,38.42_EB,39.31_EB,40.19_EB,41.07_EB,41.96_EB,41.69_EB,41.43_EB/)
QFLAME(3,5,11:20,4) = (/41.17_EB,40.90_EB,40.82_EB,40.75_EB,40.67_EB,40.59_EB,39.97_EB,39.34_EB,38.72_EB,38.10_EB/)
QFLAME(3,5,21:30,4) = (/38.49_EB,38.89_EB,39.28_EB,39.67_EB,40.07_EB,39.63_EB,39.19_EB,38.75_EB,38.31_EB,37.87_EB/)
QFLAME(3,5,0:10,5) = (/0.00_EB,26.08_EB,36.42_EB,40.60_EB,44.78_EB,46.28_EB,47.78_EB,49.28_EB,50.78_EB,50.92_EB,51.05_EB/)
QFLAME(3,5,11:20,5) = (/51.18_EB,51.32_EB,51.01_EB,50.70_EB,50.39_EB,50.09_EB,49.93_EB,49.78_EB,49.63_EB,49.48_EB/)
QFLAME(3,5,21:30,5) = (/50.10_EB,50.73_EB,51.35_EB,51.97_EB,52.59_EB,52.18_EB,51.76_EB,51.34_EB,50.93_EB,50.51_EB/)
QFLAME(3,5,0:10,6) = (/0.00_EB,27.84_EB,40.14_EB,45.52_EB,50.90_EB,53.16_EB,55.41_EB,57.66_EB,59.92_EB,60.51_EB,61.10_EB/)
QFLAME(3,5,11:20,6) = (/61.69_EB,62.27_EB,62.19_EB,62.10_EB,62.01_EB,61.92_EB,61.99_EB,62.06_EB,62.13_EB,62.20_EB/)
QFLAME(3,5,21:30,6) = (/61.77_EB,61.34_EB,60.91_EB,60.48_EB,60.06_EB,60.53_EB,61.00_EB,61.48_EB,61.95_EB,62.42_EB/)
QFLAME(3,6,0:10,1) = (/0.00_EB,20.38_EB,26.90_EB,28.55_EB,30.20_EB,30.44_EB,30.69_EB,30.93_EB,31.18_EB,30.65_EB,30.12_EB/)
QFLAME(3,6,11:20,1) = (/29.60_EB,29.07_EB,28.66_EB,28.25_EB,27.84_EB,27.43_EB,26.95_EB,26.47_EB,26.00_EB,25.52_EB/)
QFLAME(3,6,21:30,1) = (/25.10_EB,24.68_EB,24.26_EB,23.84_EB,23.42_EB,23.05_EB,22.67_EB,22.30_EB,21.93_EB,21.56_EB/)
QFLAME(3,6,0:10,2) = (/0.00_EB,20.59_EB,26.94_EB,28.58_EB,30.22_EB,30.48_EB,30.73_EB,30.98_EB,31.24_EB,30.69_EB,30.15_EB/)
QFLAME(3,6,11:20,2) = (/29.61_EB,29.06_EB,28.65_EB,28.25_EB,27.84_EB,27.43_EB,26.96_EB,26.49_EB,26.02_EB,25.54_EB/)
QFLAME(3,6,21:30,2) = (/25.14_EB,24.73_EB,24.32_EB,23.91_EB,23.51_EB,23.14_EB,22.78_EB,22.41_EB,22.05_EB,21.69_EB/)
QFLAME(3,6,0:10,3) = (/0.00_EB,22.34_EB,29.29_EB,31.03_EB,32.77_EB,33.12_EB,33.46_EB,33.81_EB,34.16_EB,33.64_EB,33.12_EB/)
QFLAME(3,6,11:20,3) = (/32.61_EB,32.09_EB,31.84_EB,31.59_EB,31.34_EB,31.09_EB,30.51_EB,29.92_EB,29.34_EB,28.75_EB/)
QFLAME(3,6,21:30,3) = (/28.76_EB,28.77_EB,28.79_EB,28.80_EB,28.81_EB,28.55_EB,28.29_EB,28.03_EB,27.77_EB,27.51_EB/)
QFLAME(3,6,0:10,4) = (/0.00_EB,24.23_EB,32.75_EB,35.57_EB,38.39_EB,39.26_EB,40.13_EB,41.00_EB,41.87_EB,41.60_EB,41.34_EB/)
QFLAME(3,6,11:20,4) = (/41.08_EB,40.82_EB,40.73_EB,40.64_EB,40.54_EB,40.45_EB,39.82_EB,39.20_EB,38.57_EB,37.95_EB/)
QFLAME(3,6,21:30,4) = (/38.35_EB,38.76_EB,39.17_EB,39.58_EB,39.99_EB,39.56_EB,39.14_EB,38.71_EB,38.29_EB,37.86_EB/)
QFLAME(3,6,0:10,5) = (/0.00_EB,26.07_EB,36.39_EB,40.56_EB,44.73_EB,46.23_EB,47.72_EB,49.22_EB,50.71_EB,50.84_EB,50.97_EB/)
QFLAME(3,6,11:20,5) = (/51.09_EB,51.22_EB,50.98_EB,50.73_EB,50.48_EB,50.24_EB,50.06_EB,49.89_EB,49.71_EB,49.53_EB/)
QFLAME(3,6,21:30,5) = (/50.11_EB,50.68_EB,51.25_EB,51.82_EB,52.40_EB,52.00_EB,51.60_EB,51.20_EB,50.81_EB,50.41_EB/)
QFLAME(3,6,0:10,6) = (/0.00_EB,27.84_EB,40.08_EB,45.47_EB,50.86_EB,53.10_EB,55.34_EB,57.58_EB,59.82_EB,60.40_EB,60.98_EB/)
QFLAME(3,6,11:20,6) = (/61.56_EB,62.14_EB,62.01_EB,61.88_EB,61.75_EB,61.62_EB,61.71_EB,61.80_EB,61.90_EB,61.99_EB/)
QFLAME(3,6,21:30,6) = (/61.53_EB,61.07_EB,60.62_EB,60.16_EB,59.70_EB,60.19_EB,60.68_EB,61.17_EB,61.66_EB,62.15_EB/)
QFLAME(3,7,0:10,1) = (/0.00_EB,20.57_EB,27.29_EB,29.05_EB,30.81_EB,31.10_EB,31.40_EB,31.69_EB,31.98_EB,31.48_EB,30.98_EB/)
QFLAME(3,7,11:20,1) = (/30.47_EB,29.97_EB,29.57_EB,29.17_EB,28.77_EB,28.37_EB,27.89_EB,27.41_EB,26.94_EB,26.46_EB/)
QFLAME(3,7,21:30,1) = (/26.04_EB,25.63_EB,25.21_EB,24.80_EB,24.39_EB,24.02_EB,23.64_EB,23.27_EB,22.90_EB,22.53_EB/)
QFLAME(3,7,0:10,2) = (/0.00_EB,20.74_EB,27.31_EB,29.07_EB,30.83_EB,31.13_EB,31.44_EB,31.74_EB,32.04_EB,31.52_EB,31.00_EB/)
QFLAME(3,7,11:20,2) = (/30.48_EB,29.96_EB,29.56_EB,29.16_EB,28.77_EB,28.37_EB,27.90_EB,27.42_EB,26.95_EB,26.48_EB/)
QFLAME(3,7,21:30,2) = (/26.07_EB,25.67_EB,25.26_EB,24.86_EB,24.45_EB,24.08_EB,23.72_EB,23.35_EB,22.98_EB,22.62_EB/)
QFLAME(3,7,0:10,3) = (/0.00_EB,22.34_EB,29.28_EB,31.08_EB,32.87_EB,33.25_EB,33.62_EB,34.00_EB,34.38_EB,33.86_EB,33.35_EB/)
QFLAME(3,7,11:20,3) = (/32.83_EB,32.32_EB,32.05_EB,31.79_EB,31.52_EB,31.26_EB,30.69_EB,30.13_EB,29.57_EB,29.01_EB/)
QFLAME(3,7,21:30,3) = (/28.92_EB,28.83_EB,28.74_EB,28.65_EB,28.56_EB,28.29_EB,28.01_EB,27.74_EB,27.47_EB,27.20_EB/)
QFLAME(3,7,0:10,4) = (/0.00_EB,24.22_EB,32.72_EB,35.54_EB,38.36_EB,39.21_EB,40.07_EB,40.92_EB,41.78_EB,41.52_EB,41.26_EB/)
QFLAME(3,7,11:20,4) = (/41.00_EB,40.74_EB,40.63_EB,40.52_EB,40.42_EB,40.31_EB,39.68_EB,39.05_EB,38.42_EB,37.79_EB/)
QFLAME(3,7,21:30,4) = (/38.22_EB,38.64_EB,39.06_EB,39.49_EB,39.91_EB,39.50_EB,39.09_EB,38.68_EB,38.27_EB,37.86_EB/)
QFLAME(3,7,0:10,5) = (/0.00_EB,26.07_EB,36.36_EB,40.52_EB,44.69_EB,46.18_EB,47.66_EB,49.15_EB,50.64_EB,50.76_EB,50.88_EB/)
QFLAME(3,7,11:20,5) = (/51.00_EB,51.12_EB,50.94_EB,50.76_EB,50.57_EB,50.39_EB,50.19_EB,49.99_EB,49.79_EB,49.59_EB/)
QFLAME(3,7,21:30,5) = (/50.11_EB,50.63_EB,51.15_EB,51.68_EB,52.20_EB,51.82_EB,51.44_EB,51.06_EB,50.68_EB,50.30_EB/)
QFLAME(3,7,0:10,6) = (/0.00_EB,27.83_EB,40.01_EB,45.41_EB,50.81_EB,53.04_EB,55.27_EB,57.50_EB,59.73_EB,60.30_EB,60.86_EB/)
QFLAME(3,7,11:20,6) = (/61.43_EB,62.00_EB,61.83_EB,61.66_EB,61.49_EB,61.32_EB,61.43_EB,61.55_EB,61.66_EB,61.78_EB/)
QFLAME(3,7,21:30,6) = (/61.29_EB,60.80_EB,60.32_EB,59.83_EB,59.34_EB,59.85_EB,60.36_EB,60.86_EB,61.37_EB,61.88_EB/)
QFLAME(3,8,0:10,1) = (/0.00_EB,20.77_EB,27.67_EB,29.55_EB,31.42_EB,31.76_EB,32.11_EB,32.45_EB,32.79_EB,32.31_EB,31.83_EB/)
QFLAME(3,8,11:20,1) = (/31.34_EB,30.86_EB,30.48_EB,30.09_EB,29.70_EB,29.32_EB,28.84_EB,28.36_EB,27.87_EB,27.39_EB/)
QFLAME(3,8,21:30,1) = (/26.99_EB,26.58_EB,26.17_EB,25.76_EB,25.35_EB,24.98_EB,24.61_EB,24.24_EB,23.87_EB,23.50_EB/)
QFLAME(3,8,0:10,2) = (/0.00_EB,20.88_EB,27.68_EB,29.56_EB,31.43_EB,31.79_EB,32.14_EB,32.49_EB,32.85_EB,32.35_EB,31.85_EB/)
QFLAME(3,8,11:20,2) = (/31.35_EB,30.85_EB,30.47_EB,30.08_EB,29.70_EB,29.31_EB,28.84_EB,28.36_EB,27.89_EB,27.41_EB/)
QFLAME(3,8,21:30,2) = (/27.01_EB,26.61_EB,26.20_EB,25.80_EB,25.40_EB,25.03_EB,24.66_EB,24.29_EB,23.92_EB,23.55_EB/)
QFLAME(3,8,0:10,3) = (/0.00_EB,22.33_EB,29.27_EB,31.12_EB,32.97_EB,33.38_EB,33.78_EB,34.19_EB,34.60_EB,34.08_EB,33.57_EB/)
QFLAME(3,8,11:20,3) = (/33.06_EB,32.54_EB,32.26_EB,31.98_EB,31.70_EB,31.42_EB,30.88_EB,30.34_EB,29.81_EB,29.27_EB/)
QFLAME(3,8,21:30,3) = (/29.08_EB,28.89_EB,28.69_EB,28.50_EB,28.31_EB,28.03_EB,27.74_EB,27.46_EB,27.17_EB,26.89_EB/)
QFLAME(3,8,0:10,4) = (/0.00_EB,24.21_EB,32.70_EB,35.51_EB,38.32_EB,39.17_EB,40.01_EB,40.85_EB,41.69_EB,41.43_EB,41.17_EB/)
QFLAME(3,8,11:20,4) = (/40.92_EB,40.66_EB,40.54_EB,40.41_EB,40.29_EB,40.17_EB,39.54_EB,38.90_EB,38.27_EB,37.64_EB/)
QFLAME(3,8,21:30,4) = (/38.08_EB,38.52_EB,38.95_EB,39.39_EB,39.83_EB,39.43_EB,39.04_EB,38.64_EB,38.25_EB,37.85_EB/)
QFLAME(3,8,0:10,5) = (/0.00_EB,26.07_EB,36.32_EB,40.48_EB,44.64_EB,46.12_EB,47.61_EB,49.09_EB,50.57_EB,50.68_EB,50.80_EB/)
QFLAME(3,8,11:20,5) = (/50.91_EB,51.03_EB,50.90_EB,50.78_EB,50.66_EB,50.54_EB,50.31_EB,50.09_EB,49.86_EB,49.64_EB/)
QFLAME(3,8,21:30,5) = (/50.11_EB,50.59_EB,51.06_EB,51.53_EB,52.01_EB,51.64_EB,51.28_EB,50.92_EB,50.56_EB,50.20_EB/)
QFLAME(3,8,0:10,6) = (/0.00_EB,27.83_EB,39.95_EB,45.35_EB,50.76_EB,52.98_EB,55.19_EB,57.41_EB,59.63_EB,60.19_EB,60.75_EB/)
QFLAME(3,8,11:20,6) = (/61.30_EB,61.86_EB,61.65_EB,61.44_EB,61.23_EB,61.02_EB,61.16_EB,61.29_EB,61.43_EB,61.57_EB/)
QFLAME(3,8,21:30,6) = (/61.05_EB,60.53_EB,60.02_EB,59.50_EB,58.98_EB,59.51_EB,60.03_EB,60.56_EB,61.08_EB,61.61_EB/)
QFLAME(3,9,0:10,1) = (/0.00_EB,20.97_EB,28.06_EB,30.05_EB,32.03_EB,32.42_EB,32.82_EB,33.21_EB,33.60_EB,33.14_EB,32.68_EB/)
QFLAME(3,9,11:20,1) = (/32.22_EB,31.76_EB,31.39_EB,31.01_EB,30.64_EB,30.26_EB,29.78_EB,29.30_EB,28.81_EB,28.33_EB/)
QFLAME(3,9,21:30,1) = (/27.93_EB,27.53_EB,27.13_EB,26.73_EB,26.32_EB,25.95_EB,25.58_EB,25.21_EB,24.84_EB,24.47_EB/)
QFLAME(3,9,0:10,2) = (/0.00_EB,21.03_EB,28.06_EB,30.05_EB,32.04_EB,32.44_EB,32.85_EB,33.25_EB,33.65_EB,33.18_EB,32.70_EB/)
QFLAME(3,9,11:20,2) = (/32.23_EB,31.75_EB,31.38_EB,31.00_EB,30.63_EB,30.25_EB,29.78_EB,29.30_EB,28.82_EB,28.35_EB/)
QFLAME(3,9,21:30,2) = (/27.94_EB,27.54_EB,27.14_EB,26.74_EB,26.34_EB,25.97_EB,25.60_EB,25.23_EB,24.86_EB,24.48_EB/)
QFLAME(3,9,0:10,3) = (/0.00_EB,22.32_EB,29.26_EB,31.17_EB,33.07_EB,33.51_EB,33.94_EB,34.38_EB,34.82_EB,34.31_EB,33.79_EB/)
QFLAME(3,9,11:20,3) = (/33.28_EB,32.77_EB,32.47_EB,32.18_EB,31.88_EB,31.59_EB,31.07_EB,30.56_EB,30.04_EB,29.52_EB/)
QFLAME(3,9,21:30,3) = (/29.23_EB,28.94_EB,28.65_EB,28.36_EB,28.07_EB,27.77_EB,27.47_EB,27.17_EB,26.87_EB,26.58_EB/)
QFLAME(3,9,0:10,4) = (/0.00_EB,24.20_EB,32.67_EB,35.48_EB,38.29_EB,39.12_EB,39.94_EB,40.77_EB,41.60_EB,41.34_EB,41.09_EB/)
QFLAME(3,9,11:20,4) = (/40.83_EB,40.58_EB,40.44_EB,40.30_EB,40.17_EB,40.03_EB,39.39_EB,38.76_EB,38.12_EB,37.49_EB/)
QFLAME(3,9,21:30,4) = (/37.94_EB,38.39_EB,38.84_EB,39.30_EB,39.75_EB,39.37_EB,38.99_EB,38.61_EB,38.23_EB,37.85_EB/)
QFLAME(3,9,0:10,5) = (/0.00_EB,26.06_EB,36.29_EB,40.44_EB,44.60_EB,46.07_EB,47.55_EB,49.02_EB,50.50_EB,50.61_EB,50.71_EB/)
QFLAME(3,9,11:20,5) = (/50.82_EB,50.93_EB,50.87_EB,50.81_EB,50.75_EB,50.69_EB,50.44_EB,50.19_EB,49.94_EB,49.69_EB/)
QFLAME(3,9,21:30,5) = (/50.11_EB,50.54_EB,50.96_EB,51.39_EB,51.81_EB,51.47_EB,51.12_EB,50.78_EB,50.44_EB,50.09_EB/)
QFLAME(3,9,0:10,6) = (/0.00_EB,27.83_EB,39.88_EB,45.29_EB,50.71_EB,52.92_EB,55.12_EB,57.33_EB,59.54_EB,60.08_EB,60.63_EB/)
QFLAME(3,9,11:20,6) = (/61.18_EB,61.72_EB,61.47_EB,61.22_EB,60.97_EB,60.72_EB,60.88_EB,61.04_EB,61.20_EB,61.36_EB/)
QFLAME(3,9,21:30,6) = (/60.81_EB,60.26_EB,59.72_EB,59.17_EB,58.62_EB,59.17_EB,59.71_EB,60.25_EB,60.79_EB,61.33_EB/)
QFLAME(3,10,0:10,1) = (/0.00_EB,21.16_EB,28.44_EB,30.54_EB,32.64_EB,33.08_EB,33.52_EB,33.96_EB,34.40_EB,33.97_EB,33.53_EB/)
QFLAME(3,10,11:20,1) = (/33.09_EB,32.66_EB,32.29_EB,31.93_EB,31.57_EB,31.21_EB,30.72_EB,30.24_EB,29.75_EB,29.27_EB/)
QFLAME(3,10,21:30,1) = (/28.87_EB,28.48_EB,28.08_EB,27.69_EB,27.29_EB,26.92_EB,26.55_EB,26.18_EB,25.81_EB,25.44_EB/)
QFLAME(3,10,0:10,2) = (/0.00_EB,21.18_EB,28.43_EB,30.54_EB,32.65_EB,33.10_EB,33.55_EB,34.01_EB,34.46_EB,34.01_EB,33.55_EB/)
QFLAME(3,10,11:20,2) = (/33.10_EB,32.65_EB,32.28_EB,31.92_EB,31.56_EB,31.19_EB,30.71_EB,30.24_EB,29.76_EB,29.28_EB/)
QFLAME(3,10,21:30,2) = (/28.88_EB,28.48_EB,28.08_EB,27.68_EB,27.28_EB,26.91_EB,26.54_EB,26.16_EB,25.79_EB,25.42_EB/)
QFLAME(3,10,0:10,3) = (/0.00_EB,22.32_EB,29.25_EB,31.21_EB,33.17_EB,33.64_EB,34.11_EB,34.57_EB,35.04_EB,34.53_EB,34.02_EB/)
QFLAME(3,10,11:20,3) = (/33.51_EB,33.00_EB,32.69_EB,32.37_EB,32.06_EB,31.75_EB,31.26_EB,30.77_EB,30.27_EB,29.78_EB/)
QFLAME(3,10,21:30,3) = (/29.39_EB,29.00_EB,28.60_EB,28.21_EB,27.82_EB,27.51_EB,27.20_EB,26.89_EB,26.58_EB,26.27_EB/)
QFLAME(3,10,0:10,4) = (/0.00_EB,24.19_EB,32.65_EB,35.45_EB,38.26_EB,39.07_EB,39.88_EB,40.69_EB,41.51_EB,41.25_EB,41.00_EB/)
QFLAME(3,10,11:20,4) = (/40.75_EB,40.50_EB,40.35_EB,40.19_EB,40.04_EB,39.89_EB,39.25_EB,38.61_EB,37.97_EB,37.34_EB/)
QFLAME(3,10,21:30,4) = (/37.80_EB,38.27_EB,38.74_EB,39.20_EB,39.67_EB,39.30_EB,38.94_EB,38.57_EB,38.21_EB,37.84_EB/)
QFLAME(3,10,0:10,5) = (/0.00_EB,26.06_EB,36.25_EB,40.40_EB,44.55_EB,46.02_EB,47.49_EB,48.96_EB,50.43_EB,50.53_EB,50.63_EB/)
QFLAME(3,10,11:20,5) = (/50.73_EB,50.83_EB,50.83_EB,50.84_EB,50.84_EB,50.84_EB,50.57_EB,50.29_EB,50.02_EB,49.74_EB/)
QFLAME(3,10,21:30,5) = (/50.12_EB,50.49_EB,50.87_EB,51.24_EB,51.61_EB,51.29_EB,50.96_EB,50.64_EB,50.31_EB,49.99_EB/)
QFLAME(3,10,0:10,6) = (/0.00_EB,27.82_EB,39.82_EB,45.24_EB,50.66_EB,52.85_EB,55.05_EB,57.25_EB,59.44_EB,59.98_EB,60.51_EB/)
QFLAME(3,10,11:20,6) = (/61.05_EB,61.59_EB,61.29_EB,61.00_EB,60.71_EB,60.42_EB,60.60_EB,60.78_EB,60.96_EB,61.15_EB/)
QFLAME(3,10,21:30,6) = (/60.57_EB,59.99_EB,59.42_EB,58.84_EB,58.27_EB,58.83_EB,59.39_EB,59.94_EB,60.50_EB,61.06_EB/)
QFLAME(3,11,0:10,1) = (/0.00_EB,21.28_EB,28.67_EB,30.83_EB,33.00_EB,33.46_EB,33.93_EB,34.40_EB,34.86_EB,34.44_EB,34.01_EB/)
QFLAME(3,11,11:20,1) = (/33.59_EB,33.16_EB,32.81_EB,32.45_EB,32.10_EB,31.74_EB,31.26_EB,30.78_EB,30.30_EB,29.81_EB/)
QFLAME(3,11,21:30,1) = (/29.41_EB,29.02_EB,28.62_EB,28.22_EB,27.82_EB,27.44_EB,27.07_EB,26.70_EB,26.33_EB,25.96_EB/)
QFLAME(3,11,0:10,2) = (/0.00_EB,21.30_EB,28.66_EB,30.83_EB,33.00_EB,33.48_EB,33.96_EB,34.43_EB,34.91_EB,34.47_EB,34.04_EB/)
QFLAME(3,11,11:20,2) = (/33.60_EB,33.16_EB,32.80_EB,32.44_EB,32.08_EB,31.72_EB,31.25_EB,30.77_EB,30.30_EB,29.83_EB/)
QFLAME(3,11,21:30,2) = (/29.42_EB,29.02_EB,28.62_EB,28.22_EB,27.82_EB,27.44_EB,27.07_EB,26.69_EB,26.32_EB,25.94_EB/)
QFLAME(3,11,0:10,3) = (/0.00_EB,22.32_EB,29.40_EB,31.43_EB,33.47_EB,33.96_EB,34.45_EB,34.94_EB,35.43_EB,34.94_EB,34.45_EB/)
QFLAME(3,11,11:20,3) = (/33.96_EB,33.47_EB,33.16_EB,32.85_EB,32.54_EB,32.23_EB,31.74_EB,31.26_EB,30.77_EB,30.28_EB/)
QFLAME(3,11,21:30,3) = (/29.89_EB,29.49_EB,29.09_EB,28.70_EB,28.30_EB,27.99_EB,27.67_EB,27.35_EB,27.04_EB,26.72_EB/)
QFLAME(3,11,0:10,4) = (/0.00_EB,24.19_EB,32.62_EB,35.42_EB,38.22_EB,39.03_EB,39.84_EB,40.64_EB,41.45_EB,41.20_EB,40.96_EB/)
QFLAME(3,11,11:20,4) = (/40.71_EB,40.46_EB,40.30_EB,40.15_EB,39.99_EB,39.83_EB,39.21_EB,38.59_EB,37.97_EB,37.35_EB/)
QFLAME(3,11,21:30,4) = (/37.79_EB,38.23_EB,38.67_EB,39.11_EB,39.55_EB,39.19_EB,38.83_EB,38.46_EB,38.10_EB,37.73_EB/)
QFLAME(3,11,0:10,5) = (/0.00_EB,26.05_EB,36.22_EB,40.36_EB,44.50_EB,45.96_EB,47.43_EB,48.89_EB,50.36_EB,50.46_EB,50.55_EB/)
QFLAME(3,11,11:20,5) = (/50.65_EB,50.75_EB,50.75_EB,50.76_EB,50.76_EB,50.77_EB,50.49_EB,50.22_EB,49.94_EB,49.67_EB/)
QFLAME(3,11,21:30,5) = (/50.02_EB,50.38_EB,50.73_EB,51.09_EB,51.44_EB,51.12_EB,50.81_EB,50.49_EB,50.17_EB,49.85_EB/)
QFLAME(3,11,0:10,6) = (/0.00_EB,27.82_EB,39.79_EB,45.18_EB,50.56_EB,52.76_EB,54.96_EB,57.15_EB,59.35_EB,59.89_EB,60.43_EB/)
QFLAME(3,11,11:20,6) = (/60.97_EB,61.50_EB,61.20_EB,60.89_EB,60.58_EB,60.28_EB,60.47_EB,60.66_EB,60.86_EB,61.05_EB/)
QFLAME(3,11,21:30,6) = (/60.59_EB,60.12_EB,59.65_EB,59.19_EB,58.72_EB,59.19_EB,59.66_EB,60.13_EB,60.59_EB,61.06_EB/)
QFLAME(3,12,0:10,1) = (/0.00_EB,21.40_EB,28.90_EB,31.12_EB,33.35_EB,33.84_EB,34.33_EB,34.83_EB,35.32_EB,34.91_EB,34.49_EB/)
QFLAME(3,12,11:20,1) = (/34.08_EB,33.67_EB,33.32_EB,32.97_EB,32.62_EB,32.27_EB,31.79_EB,31.31_EB,30.84_EB,30.36_EB/)
QFLAME(3,12,21:30,1) = (/29.96_EB,29.55_EB,29.15_EB,28.74_EB,28.34_EB,27.97_EB,27.59_EB,27.22_EB,26.85_EB,26.48_EB/)
QFLAME(3,12,0:10,2) = (/0.00_EB,21.41_EB,28.89_EB,31.12_EB,33.36_EB,33.86_EB,34.36_EB,34.86_EB,35.36_EB,34.94_EB,34.52_EB/)
QFLAME(3,12,11:20,2) = (/34.10_EB,33.67_EB,33.32_EB,32.96_EB,32.61_EB,32.25_EB,31.78_EB,31.31_EB,30.84_EB,30.37_EB/)
QFLAME(3,12,21:30,2) = (/29.97_EB,29.56_EB,29.16_EB,28.76_EB,28.35_EB,27.97_EB,27.60_EB,27.22_EB,26.85_EB,26.47_EB/)
QFLAME(3,12,0:10,3) = (/0.00_EB,22.33_EB,29.54_EB,31.66_EB,33.77_EB,34.29_EB,34.80_EB,35.31_EB,35.83_EB,35.36_EB,34.89_EB/)
QFLAME(3,12,11:20,3) = (/34.42_EB,33.95_EB,33.64_EB,33.33_EB,33.02_EB,32.71_EB,32.23_EB,31.75_EB,31.26_EB,30.78_EB/)
QFLAME(3,12,21:30,3) = (/30.38_EB,29.98_EB,29.58_EB,29.19_EB,28.79_EB,28.46_EB,28.14_EB,27.82_EB,27.50_EB,27.18_EB/)
QFLAME(3,12,0:10,4) = (/0.00_EB,24.18_EB,32.59_EB,35.39_EB,38.19_EB,38.99_EB,39.79_EB,40.59_EB,41.40_EB,41.15_EB,40.91_EB/)
QFLAME(3,12,11:20,4) = (/40.67_EB,40.42_EB,40.26_EB,40.10_EB,39.94_EB,39.77_EB,39.17_EB,38.57_EB,37.96_EB,37.36_EB/)
QFLAME(3,12,21:30,4) = (/37.78_EB,38.19_EB,38.61_EB,39.02_EB,39.44_EB,39.08_EB,38.71_EB,38.35_EB,37.99_EB,37.62_EB/)
QFLAME(3,12,0:10,5) = (/0.00_EB,26.05_EB,36.20_EB,40.32_EB,44.45_EB,45.91_EB,47.37_EB,48.83_EB,50.29_EB,50.38_EB,50.48_EB/)
QFLAME(3,12,11:20,5) = (/50.57_EB,50.67_EB,50.67_EB,50.68_EB,50.69_EB,50.69_EB,50.42_EB,50.14_EB,49.87_EB,49.59_EB/)
QFLAME(3,12,21:30,5) = (/49.93_EB,50.26_EB,50.60_EB,50.93_EB,51.27_EB,50.96_EB,50.65_EB,50.34_EB,50.03_EB,49.72_EB/)
QFLAME(3,12,0:10,6) = (/0.00_EB,27.83_EB,39.77_EB,45.12_EB,50.47_EB,52.66_EB,54.86_EB,57.06_EB,59.25_EB,59.80_EB,60.34_EB/)
QFLAME(3,12,11:20,6) = (/60.88_EB,61.42_EB,61.10_EB,60.78_EB,60.45_EB,60.13_EB,60.34_EB,60.55_EB,60.75_EB,60.96_EB/)
QFLAME(3,12,21:30,6) = (/60.60_EB,60.25_EB,59.89_EB,59.53_EB,59.18_EB,59.55_EB,59.93_EB,60.31_EB,60.68_EB,61.06_EB/)
QFLAME(3,13,0:10,1) = (/0.00_EB,21.52_EB,29.12_EB,31.42_EB,33.71_EB,34.22_EB,34.74_EB,35.26_EB,35.77_EB,35.38_EB,34.98_EB/)
QFLAME(3,13,11:20,1) = (/34.58_EB,34.18_EB,33.84_EB,33.49_EB,33.14_EB,32.80_EB,32.33_EB,31.85_EB,31.38_EB,30.91_EB/)
QFLAME(3,13,21:30,1) = (/30.50_EB,30.09_EB,29.68_EB,29.27_EB,28.86_EB,28.49_EB,28.11_EB,27.74_EB,27.37_EB,26.99_EB/)
QFLAME(3,13,0:10,2) = (/0.00_EB,21.53_EB,29.11_EB,31.41_EB,33.71_EB,34.24_EB,34.76_EB,35.29_EB,35.81_EB,35.40_EB,35.00_EB/)
QFLAME(3,13,11:20,2) = (/34.59_EB,34.19_EB,33.83_EB,33.48_EB,33.13_EB,32.78_EB,32.31_EB,31.85_EB,31.38_EB,30.92_EB/)
QFLAME(3,13,21:30,2) = (/30.51_EB,30.10_EB,29.70_EB,29.29_EB,28.88_EB,28.51_EB,28.13_EB,27.75_EB,27.37_EB,26.99_EB/)
QFLAME(3,13,0:10,3) = (/0.00_EB,22.33_EB,29.69_EB,31.88_EB,34.07_EB,34.61_EB,35.15_EB,35.68_EB,36.22_EB,35.77_EB,35.32_EB/)
QFLAME(3,13,11:20,3) = (/34.87_EB,34.42_EB,34.11_EB,33.80_EB,33.50_EB,33.19_EB,32.71_EB,32.24_EB,31.76_EB,31.28_EB/)
QFLAME(3,13,21:30,3) = (/30.88_EB,30.48_EB,30.08_EB,29.67_EB,29.27_EB,28.94_EB,28.62_EB,28.29_EB,27.96_EB,27.63_EB/)
QFLAME(3,13,0:10,4) = (/0.00_EB,24.18_EB,32.57_EB,35.36_EB,38.15_EB,38.95_EB,39.74_EB,40.54_EB,41.34_EB,41.10_EB,40.86_EB/)
QFLAME(3,13,11:20,4) = (/40.62_EB,40.39_EB,40.22_EB,40.05_EB,39.89_EB,39.72_EB,39.13_EB,38.55_EB,37.96_EB,37.37_EB/)
QFLAME(3,13,21:30,4) = (/37.76_EB,38.15_EB,38.54_EB,38.94_EB,39.33_EB,38.96_EB,38.60_EB,38.24_EB,37.88_EB,37.51_EB/)
QFLAME(3,13,0:10,5) = (/0.00_EB,26.04_EB,36.17_EB,40.28_EB,44.39_EB,45.85_EB,47.31_EB,48.76_EB,50.22_EB,50.31_EB,50.40_EB/)
QFLAME(3,13,11:20,5) = (/50.49_EB,50.59_EB,50.60_EB,50.60_EB,50.61_EB,50.62_EB,50.34_EB,50.07_EB,49.79_EB,49.51_EB/)
QFLAME(3,13,21:30,5) = (/49.83_EB,50.15_EB,50.46_EB,50.78_EB,51.10_EB,50.79_EB,50.49_EB,50.19_EB,49.88_EB,49.58_EB/)
QFLAME(3,13,0:10,6) = (/0.00_EB,27.83_EB,39.75_EB,45.06_EB,50.37_EB,52.57_EB,54.76_EB,56.96_EB,59.16_EB,59.71_EB,60.25_EB/)
QFLAME(3,13,11:20,6) = (/60.80_EB,61.34_EB,61.00_EB,60.66_EB,60.32_EB,59.98_EB,60.20_EB,60.43_EB,60.65_EB,60.87_EB/)
QFLAME(3,13,21:30,6) = (/60.62_EB,60.37_EB,60.13_EB,59.88_EB,59.63_EB,59.92_EB,60.20_EB,60.49_EB,60.77_EB,61.06_EB/)
QFLAME(3,14,0:10,1) = (/0.00_EB,21.64_EB,29.35_EB,31.71_EB,34.06_EB,34.60_EB,35.15_EB,35.69_EB,36.23_EB,35.84_EB,35.46_EB/)
QFLAME(3,14,11:20,1) = (/35.07_EB,34.69_EB,34.35_EB,34.01_EB,33.67_EB,33.33_EB,32.86_EB,32.39_EB,31.92_EB,31.45_EB/)
QFLAME(3,14,21:30,1) = (/31.04_EB,30.63_EB,30.21_EB,29.80_EB,29.39_EB,29.01_EB,28.64_EB,28.26_EB,27.89_EB,27.51_EB/)
QFLAME(3,14,0:10,2) = (/0.00_EB,21.64_EB,29.34_EB,31.70_EB,34.07_EB,34.61_EB,35.16_EB,35.71_EB,36.26_EB,35.87_EB,35.48_EB/)
QFLAME(3,14,11:20,2) = (/35.09_EB,34.70_EB,34.35_EB,34.00_EB,33.65_EB,33.31_EB,32.84_EB,32.38_EB,31.92_EB,31.46_EB/)
QFLAME(3,14,21:30,2) = (/31.05_EB,30.64_EB,30.24_EB,29.83_EB,29.42_EB,29.04_EB,28.66_EB,28.28_EB,27.90_EB,27.52_EB/)
QFLAME(3,14,0:10,3) = (/0.00_EB,22.34_EB,29.83_EB,32.10_EB,34.38_EB,34.94_EB,35.50_EB,36.06_EB,36.62_EB,36.19_EB,35.76_EB/)
QFLAME(3,14,11:20,3) = (/35.33_EB,34.90_EB,34.59_EB,34.28_EB,33.97_EB,33.67_EB,33.20_EB,32.73_EB,32.25_EB,31.78_EB/)
QFLAME(3,14,21:30,3) = (/31.38_EB,30.97_EB,30.57_EB,30.16_EB,29.76_EB,29.42_EB,29.09_EB,28.76_EB,28.42_EB,28.09_EB/)
QFLAME(3,14,0:10,4) = (/0.00_EB,24.18_EB,32.54_EB,35.33_EB,38.11_EB,38.90_EB,39.70_EB,40.49_EB,41.29_EB,41.05_EB,40.82_EB/)
QFLAME(3,14,11:20,4) = (/40.58_EB,40.35_EB,40.18_EB,40.00_EB,39.83_EB,39.66_EB,39.09_EB,38.52_EB,37.95_EB,37.38_EB/)
QFLAME(3,14,21:30,4) = (/37.75_EB,38.12_EB,38.48_EB,38.85_EB,39.21_EB,38.85_EB,38.49_EB,38.13_EB,37.77_EB,37.41_EB/)
QFLAME(3,14,0:10,5) = (/0.00_EB,26.04_EB,36.14_EB,40.24_EB,44.34_EB,45.79_EB,47.24_EB,48.70_EB,50.15_EB,50.24_EB,50.33_EB/)
QFLAME(3,14,11:20,5) = (/50.42_EB,50.51_EB,50.52_EB,50.53_EB,50.54_EB,50.55_EB,50.27_EB,49.99_EB,49.71_EB,49.44_EB/)
QFLAME(3,14,21:30,5) = (/49.73_EB,50.03_EB,50.33_EB,50.63_EB,50.93_EB,50.63_EB,50.33_EB,50.04_EB,49.74_EB,49.44_EB/)
QFLAME(3,14,0:10,6) = (/0.00_EB,27.83_EB,39.73_EB,45.00_EB,50.27_EB,52.47_EB,54.67_EB,56.87_EB,59.06_EB,59.61_EB,60.16_EB/)
QFLAME(3,14,11:20,6) = (/60.71_EB,61.26_EB,60.91_EB,60.55_EB,60.19_EB,59.84_EB,60.07_EB,60.31_EB,60.54_EB,60.78_EB/)
QFLAME(3,14,21:30,6) = (/60.64_EB,60.50_EB,60.36_EB,60.22_EB,60.08_EB,60.28_EB,60.47_EB,60.67_EB,60.87_EB,61.06_EB/)
QFLAME(3,15,0:10,1) = (/0.00_EB,21.75_EB,29.58_EB,32.00_EB,34.41_EB,34.98_EB,35.55_EB,36.12_EB,36.69_EB,36.31_EB,35.94_EB/)
QFLAME(3,15,11:20,1) = (/35.57_EB,35.20_EB,34.86_EB,34.53_EB,34.19_EB,33.86_EB,33.39_EB,32.93_EB,32.46_EB,32.00_EB/)
QFLAME(3,15,21:30,1) = (/31.58_EB,31.16_EB,30.74_EB,30.33_EB,29.91_EB,29.53_EB,29.16_EB,28.78_EB,28.40_EB,28.03_EB/)
QFLAME(3,15,0:10,2) = (/0.00_EB,21.76_EB,29.57_EB,31.99_EB,34.42_EB,34.99_EB,35.57_EB,36.14_EB,36.71_EB,36.34_EB,35.96_EB/)
QFLAME(3,15,11:20,2) = (/35.59_EB,35.21_EB,34.87_EB,34.52_EB,34.18_EB,33.83_EB,33.38_EB,32.92_EB,32.46_EB,32.01_EB/)
QFLAME(3,15,21:30,2) = (/31.60_EB,31.19_EB,30.77_EB,30.36_EB,29.95_EB,29.57_EB,29.19_EB,28.81_EB,28.43_EB,28.05_EB/)
QFLAME(3,15,0:10,3) = (/0.00_EB,22.35_EB,29.98_EB,32.33_EB,34.68_EB,35.26_EB,35.84_EB,36.43_EB,37.01_EB,36.60_EB,36.19_EB/)
QFLAME(3,15,11:20,3) = (/35.78_EB,35.37_EB,35.07_EB,34.76_EB,34.45_EB,34.14_EB,33.68_EB,33.21_EB,32.75_EB,32.28_EB/)
QFLAME(3,15,21:30,3) = (/31.88_EB,31.47_EB,31.06_EB,30.65_EB,30.24_EB,29.90_EB,29.56_EB,29.22_EB,28.88_EB,28.54_EB/)
QFLAME(3,15,0:10,4) = (/0.00_EB,24.17_EB,32.51_EB,35.29_EB,38.07_EB,38.86_EB,39.65_EB,40.44_EB,41.23_EB,41.00_EB,40.77_EB/)
QFLAME(3,15,11:20,4) = (/40.54_EB,40.31_EB,40.13_EB,39.96_EB,39.78_EB,39.61_EB,39.05_EB,38.50_EB,37.95_EB,37.40_EB/)
QFLAME(3,15,21:30,4) = (/37.74_EB,38.08_EB,38.42_EB,38.76_EB,39.10_EB,38.74_EB,38.38_EB,38.02_EB,37.66_EB,37.30_EB/)
QFLAME(3,15,0:10,5) = (/0.00_EB,26.04_EB,36.11_EB,40.20_EB,44.29_EB,45.73_EB,47.18_EB,48.63_EB,50.08_EB,50.16_EB,50.25_EB/)
QFLAME(3,15,11:20,5) = (/50.34_EB,50.43_EB,50.44_EB,50.45_EB,50.46_EB,50.47_EB,50.20_EB,49.92_EB,49.64_EB,49.36_EB/)
QFLAME(3,15,21:30,5) = (/49.64_EB,49.92_EB,50.20_EB,50.47_EB,50.75_EB,50.46_EB,50.18_EB,49.89_EB,49.60_EB,49.31_EB/)
QFLAME(3,15,0:10,6) = (/0.00_EB,27.83_EB,39.70_EB,44.94_EB,50.18_EB,52.38_EB,54.57_EB,56.77_EB,58.97_EB,59.52_EB,60.08_EB/)
QFLAME(3,15,11:20,6) = (/60.63_EB,61.18_EB,60.81_EB,60.44_EB,60.06_EB,59.69_EB,59.94_EB,60.19_EB,60.44_EB,60.69_EB/)
QFLAME(3,15,21:30,6) = (/60.66_EB,60.63_EB,60.60_EB,60.57_EB,60.54_EB,60.64_EB,60.75_EB,60.85_EB,60.96_EB,61.06_EB/)
QFLAME(3,16,0:10,1) = (/0.00_EB,21.87_EB,29.80_EB,32.29_EB,34.77_EB,35.36_EB,35.96_EB,36.55_EB,37.14_EB,36.78_EB,36.42_EB/)
QFLAME(3,16,11:20,1) = (/36.07_EB,35.71_EB,35.38_EB,35.05_EB,34.72_EB,34.39_EB,33.93_EB,33.46_EB,33.00_EB,32.54_EB/)
QFLAME(3,16,21:30,1) = (/32.12_EB,31.70_EB,31.28_EB,30.85_EB,30.43_EB,30.05_EB,29.68_EB,29.30_EB,28.92_EB,28.54_EB/)
QFLAME(3,16,0:10,2) = (/0.00_EB,21.88_EB,29.80_EB,32.29_EB,34.78_EB,35.37_EB,35.97_EB,36.56_EB,37.16_EB,36.80_EB,36.44_EB/)
QFLAME(3,16,11:20,2) = (/36.08_EB,35.73_EB,35.39_EB,35.04_EB,34.70_EB,34.36_EB,33.91_EB,33.46_EB,33.01_EB,32.55_EB/)
QFLAME(3,16,21:30,2) = (/32.14_EB,31.73_EB,31.31_EB,30.90_EB,30.48_EB,30.10_EB,29.72_EB,29.34_EB,28.95_EB,28.57_EB/)
QFLAME(3,16,0:10,3) = (/0.00_EB,22.35_EB,30.13_EB,32.55_EB,34.98_EB,35.58_EB,36.19_EB,36.80_EB,37.40_EB,37.01_EB,36.63_EB/)
QFLAME(3,16,11:20,3) = (/36.24_EB,35.85_EB,35.54_EB,35.24_EB,34.93_EB,34.62_EB,34.16_EB,33.70_EB,33.24_EB,32.79_EB/)
QFLAME(3,16,21:30,3) = (/32.37_EB,31.96_EB,31.55_EB,31.14_EB,30.72_EB,30.38_EB,30.03_EB,29.69_EB,29.34_EB,29.00_EB/)
QFLAME(3,16,0:10,4) = (/0.00_EB,24.17_EB,32.49_EB,35.26_EB,38.04_EB,38.82_EB,39.61_EB,40.39_EB,41.17_EB,40.95_EB,40.72_EB/)
QFLAME(3,16,11:20,4) = (/40.50_EB,40.27_EB,40.09_EB,39.91_EB,39.73_EB,39.55_EB,39.02_EB,38.48_EB,37.94_EB,37.41_EB/)
QFLAME(3,16,21:30,4) = (/37.72_EB,38.04_EB,38.35_EB,38.67_EB,38.99_EB,38.63_EB,38.27_EB,37.91_EB,37.55_EB,37.19_EB/)
QFLAME(3,16,0:10,5) = (/0.00_EB,26.03_EB,36.08_EB,40.16_EB,44.23_EB,45.68_EB,47.12_EB,48.56_EB,50.01_EB,50.09_EB,50.18_EB/)
QFLAME(3,16,11:20,5) = (/50.26_EB,50.34_EB,50.36_EB,50.37_EB,50.39_EB,50.40_EB,50.12_EB,49.84_EB,49.56_EB,49.28_EB/)
QFLAME(3,16,21:30,5) = (/49.54_EB,49.80_EB,50.06_EB,50.32_EB,50.58_EB,50.30_EB,50.02_EB,49.74_EB,49.45_EB,49.17_EB/)
QFLAME(3,16,0:10,6) = (/0.00_EB,27.84_EB,39.68_EB,44.88_EB,50.08_EB,52.28_EB,54.48_EB,56.68_EB,58.88_EB,59.43_EB,59.99_EB/)
QFLAME(3,16,11:20,6) = (/60.54_EB,61.10_EB,60.71_EB,60.32_EB,59.93_EB,59.54_EB,59.81_EB,60.07_EB,60.33_EB,60.60_EB/)
QFLAME(3,16,21:30,6) = (/60.67_EB,60.75_EB,60.83_EB,60.91_EB,60.99_EB,61.01_EB,61.02_EB,61.03_EB,61.05_EB,61.06_EB/)
QFLAME(3,17,0:10,1) = (/0.00_EB,21.99_EB,30.03_EB,32.58_EB,35.12_EB,35.74_EB,36.36_EB,36.98_EB,37.60_EB,37.25_EB,36.91_EB/)
QFLAME(3,17,11:20,1) = (/36.56_EB,36.22_EB,35.89_EB,35.57_EB,35.24_EB,34.92_EB,34.46_EB,34.00_EB,33.55_EB,33.09_EB/)
QFLAME(3,17,21:30,1) = (/32.66_EB,32.24_EB,31.81_EB,31.38_EB,30.96_EB,30.58_EB,30.20_EB,29.82_EB,29.44_EB,29.06_EB/)
QFLAME(3,17,0:10,2) = (/0.00_EB,21.99_EB,30.02_EB,32.58_EB,35.13_EB,35.75_EB,36.37_EB,36.99_EB,37.61_EB,37.27_EB,36.93_EB/)
QFLAME(3,17,11:20,2) = (/36.58_EB,36.24_EB,35.90_EB,35.56_EB,35.23_EB,34.89_EB,34.44_EB,34.00_EB,33.55_EB,33.10_EB/)
QFLAME(3,17,21:30,2) = (/32.68_EB,32.27_EB,31.85_EB,31.43_EB,31.02_EB,30.63_EB,30.25_EB,29.87_EB,29.48_EB,29.10_EB/)
QFLAME(3,17,0:10,3) = (/0.00_EB,22.36_EB,30.27_EB,32.77_EB,35.28_EB,35.91_EB,36.54_EB,37.17_EB,37.80_EB,37.43_EB,37.06_EB/)
QFLAME(3,17,11:20,3) = (/36.69_EB,36.32_EB,36.02_EB,35.71_EB,35.41_EB,35.10_EB,34.65_EB,34.19_EB,33.74_EB,33.29_EB/)
QFLAME(3,17,21:30,3) = (/32.87_EB,32.46_EB,32.04_EB,31.62_EB,31.21_EB,30.86_EB,30.51_EB,30.16_EB,29.80_EB,29.45_EB/)
QFLAME(3,17,0:10,4) = (/0.00_EB,24.16_EB,32.46_EB,35.23_EB,38.00_EB,38.78_EB,39.56_EB,40.34_EB,41.12_EB,40.90_EB,40.68_EB/)
QFLAME(3,17,11:20,4) = (/40.45_EB,40.23_EB,40.05_EB,39.86_EB,39.68_EB,39.50_EB,38.98_EB,38.46_EB,37.94_EB,37.42_EB/)
QFLAME(3,17,21:30,4) = (/37.71_EB,38.00_EB,38.29_EB,38.58_EB,38.87_EB,38.51_EB,38.15_EB,37.80_EB,37.44_EB,37.08_EB/)
QFLAME(3,17,0:10,5) = (/0.00_EB,26.03_EB,36.05_EB,40.12_EB,44.18_EB,45.62_EB,47.06_EB,48.50_EB,49.94_EB,50.02_EB,50.10_EB/)
QFLAME(3,17,11:20,5) = (/50.18_EB,50.26_EB,50.28_EB,50.30_EB,50.31_EB,50.33_EB,50.05_EB,49.77_EB,49.49_EB,49.20_EB/)
QFLAME(3,17,21:30,5) = (/49.45_EB,49.69_EB,49.93_EB,50.17_EB,50.41_EB,50.13_EB,49.86_EB,49.59_EB,49.31_EB,49.04_EB/)
QFLAME(3,17,0:10,6) = (/0.00_EB,27.84_EB,39.66_EB,44.82_EB,49.98_EB,52.18_EB,54.38_EB,56.58_EB,58.78_EB,59.34_EB,59.90_EB/)
QFLAME(3,17,11:20,6) = (/60.46_EB,61.02_EB,60.61_EB,60.21_EB,59.80_EB,59.40_EB,59.67_EB,59.95_EB,60.23_EB,60.50_EB/)
QFLAME(3,17,21:30,6) = (/60.69_EB,60.88_EB,61.07_EB,61.26_EB,61.45_EB,61.37_EB,61.29_EB,61.21_EB,61.14_EB,61.06_EB/)
QFLAME(3,18,0:10,1) = (/0.00_EB,22.11_EB,30.26_EB,32.87_EB,35.48_EB,36.12_EB,36.77_EB,37.41_EB,38.05_EB,37.72_EB,37.39_EB/)
QFLAME(3,18,11:20,1) = (/37.06_EB,36.72_EB,36.40_EB,36.08_EB,35.76_EB,35.44_EB,34.99_EB,34.54_EB,34.09_EB,33.63_EB/)
QFLAME(3,18,21:30,1) = (/33.20_EB,32.77_EB,32.34_EB,31.91_EB,31.48_EB,31.10_EB,30.72_EB,30.34_EB,29.96_EB,29.58_EB/)
QFLAME(3,18,0:10,2) = (/0.00_EB,22.11_EB,30.25_EB,32.87_EB,35.49_EB,36.13_EB,36.77_EB,37.42_EB,38.06_EB,37.73_EB,37.41_EB/)
QFLAME(3,18,11:20,2) = (/37.08_EB,36.75_EB,36.42_EB,36.09_EB,35.75_EB,35.42_EB,34.98_EB,34.53_EB,34.09_EB,33.65_EB/)
QFLAME(3,18,21:30,2) = (/33.23_EB,32.81_EB,32.39_EB,31.97_EB,31.55_EB,31.17_EB,30.78_EB,30.39_EB,30.01_EB,29.62_EB/)
QFLAME(3,18,0:10,3) = (/0.00_EB,22.36_EB,30.42_EB,33.00_EB,35.58_EB,36.23_EB,36.89_EB,37.54_EB,38.19_EB,37.84_EB,37.50_EB/)
QFLAME(3,18,11:20,3) = (/37.15_EB,36.80_EB,36.49_EB,36.19_EB,35.88_EB,35.58_EB,35.13_EB,34.68_EB,34.23_EB,33.79_EB/)
QFLAME(3,18,21:30,3) = (/33.37_EB,32.95_EB,32.53_EB,32.11_EB,31.69_EB,31.34_EB,30.98_EB,30.62_EB,30.27_EB,29.91_EB/)
QFLAME(3,18,0:10,4) = (/0.00_EB,24.16_EB,32.43_EB,35.20_EB,37.96_EB,38.74_EB,39.51_EB,40.29_EB,41.06_EB,40.85_EB,40.63_EB/)
QFLAME(3,18,11:20,4) = (/40.41_EB,40.19_EB,40.00_EB,39.82_EB,39.63_EB,39.44_EB,38.94_EB,38.44_EB,37.93_EB,37.43_EB/)
QFLAME(3,18,21:30,4) = (/37.70_EB,37.96_EB,38.23_EB,38.49_EB,38.76_EB,38.40_EB,38.04_EB,37.68_EB,37.33_EB,36.97_EB/)
QFLAME(3,18,0:10,5) = (/0.00_EB,26.02_EB,36.02_EB,40.07_EB,44.13_EB,45.56_EB,47.00_EB,48.43_EB,49.87_EB,49.95_EB,50.02_EB/)
QFLAME(3,18,11:20,5) = (/50.10_EB,50.18_EB,50.20_EB,50.22_EB,50.24_EB,50.25_EB,49.97_EB,49.69_EB,49.41_EB,49.13_EB/)
QFLAME(3,18,21:30,5) = (/49.35_EB,49.57_EB,49.79_EB,50.02_EB,50.24_EB,49.97_EB,49.70_EB,49.43_EB,49.17_EB,48.90_EB/)
QFLAME(3,18,0:10,6) = (/0.00_EB,27.84_EB,39.64_EB,44.76_EB,49.89_EB,52.09_EB,54.29_EB,56.49_EB,58.69_EB,59.25_EB,59.81_EB/)
QFLAME(3,18,11:20,6) = (/60.38_EB,60.94_EB,60.52_EB,60.10_EB,59.67_EB,59.25_EB,59.54_EB,59.83_EB,60.12_EB,60.41_EB/)
QFLAME(3,18,21:30,6) = (/60.71_EB,61.01_EB,61.30_EB,61.60_EB,61.90_EB,61.73_EB,61.56_EB,61.40_EB,61.23_EB,61.06_EB/)
QFLAME(3,19,0:10,1) = (/0.00_EB,22.23_EB,30.48_EB,33.16_EB,35.83_EB,36.50_EB,37.17_EB,37.84_EB,38.51_EB,38.19_EB,37.87_EB/)
QFLAME(3,19,11:20,1) = (/37.55_EB,37.23_EB,36.92_EB,36.60_EB,36.29_EB,35.97_EB,35.53_EB,35.08_EB,34.63_EB,34.18_EB/)
QFLAME(3,19,21:30,1) = (/33.74_EB,33.31_EB,32.87_EB,32.44_EB,32.00_EB,31.62_EB,31.24_EB,30.86_EB,30.48_EB,30.10_EB/)
QFLAME(3,19,0:10,2) = (/0.00_EB,22.22_EB,30.48_EB,33.16_EB,35.84_EB,36.51_EB,37.18_EB,37.84_EB,38.51_EB,38.20_EB,37.89_EB/)
QFLAME(3,19,11:20,2) = (/37.58_EB,37.27_EB,36.94_EB,36.61_EB,36.28_EB,35.95_EB,35.51_EB,35.07_EB,34.63_EB,34.19_EB/)
QFLAME(3,19,21:30,2) = (/33.77_EB,33.35_EB,32.93_EB,32.51_EB,32.08_EB,31.70_EB,31.31_EB,30.92_EB,30.54_EB,30.15_EB/)
QFLAME(3,19,0:10,3) = (/0.00_EB,22.37_EB,30.56_EB,33.22_EB,35.88_EB,36.56_EB,37.23_EB,37.91_EB,38.59_EB,38.26_EB,37.93_EB/)
QFLAME(3,19,11:20,3) = (/37.60_EB,37.27_EB,36.97_EB,36.67_EB,36.36_EB,36.06_EB,35.61_EB,35.17_EB,34.73_EB,34.29_EB/)
QFLAME(3,19,21:30,3) = (/33.87_EB,33.44_EB,33.02_EB,32.60_EB,32.18_EB,31.82_EB,31.45_EB,31.09_EB,30.73_EB,30.36_EB/)
QFLAME(3,19,0:10,4) = (/0.00_EB,24.15_EB,32.41_EB,35.17_EB,37.92_EB,38.70_EB,39.47_EB,40.24_EB,41.01_EB,40.80_EB,40.58_EB/)
QFLAME(3,19,11:20,4) = (/40.37_EB,40.15_EB,39.96_EB,39.77_EB,39.58_EB,39.38_EB,38.90_EB,38.41_EB,37.93_EB,37.44_EB/)
QFLAME(3,19,21:30,4) = (/37.68_EB,37.92_EB,38.16_EB,38.40_EB,38.64_EB,38.29_EB,37.93_EB,37.57_EB,37.22_EB,36.86_EB/)
QFLAME(3,19,0:10,5) = (/0.00_EB,26.02_EB,35.99_EB,40.03_EB,44.07_EB,45.50_EB,46.93_EB,48.37_EB,49.80_EB,49.87_EB,49.95_EB/)
QFLAME(3,19,11:20,5) = (/50.03_EB,50.10_EB,50.12_EB,50.14_EB,50.16_EB,50.18_EB,49.90_EB,49.62_EB,49.33_EB,49.05_EB/)
QFLAME(3,19,21:30,5) = (/49.25_EB,49.46_EB,49.66_EB,49.86_EB,50.07_EB,49.81_EB,49.54_EB,49.28_EB,49.02_EB,48.76_EB/)
QFLAME(3,19,0:10,6) = (/0.00_EB,27.84_EB,39.61_EB,44.70_EB,49.79_EB,51.99_EB,54.19_EB,56.39_EB,58.59_EB,59.16_EB,59.73_EB/)
QFLAME(3,19,11:20,6) = (/60.29_EB,60.86_EB,60.42_EB,59.98_EB,59.54_EB,59.11_EB,59.41_EB,59.71_EB,60.02_EB,60.32_EB/)
QFLAME(3,19,21:30,6) = (/60.73_EB,61.13_EB,61.54_EB,61.95_EB,62.35_EB,62.09_EB,61.84_EB,61.58_EB,61.32_EB,61.06_EB/)
QFLAME(3,20,0:10,1) = (/0.00_EB,22.34_EB,30.71_EB,33.45_EB,36.18_EB,36.88_EB,37.58_EB,38.27_EB,38.97_EB,38.66_EB,38.35_EB/)
QFLAME(3,20,11:20,1) = (/38.05_EB,37.74_EB,37.43_EB,37.12_EB,36.81_EB,36.50_EB,36.06_EB,35.61_EB,35.17_EB,34.73_EB/)
QFLAME(3,20,21:30,1) = (/34.29_EB,33.85_EB,33.41_EB,32.96_EB,32.52_EB,32.14_EB,31.76_EB,31.38_EB,31.00_EB,30.61_EB/)
QFLAME(3,20,0:10,2) = (/0.00_EB,22.34_EB,30.71_EB,33.45_EB,36.19_EB,36.89_EB,37.58_EB,38.27_EB,38.96_EB,38.67_EB,38.37_EB/)
QFLAME(3,20,11:20,2) = (/38.07_EB,37.78_EB,37.45_EB,37.13_EB,36.80_EB,36.48_EB,36.04_EB,35.61_EB,35.17_EB,34.74_EB/)
QFLAME(3,20,21:30,2) = (/34.31_EB,33.89_EB,33.47_EB,33.04_EB,32.62_EB,32.23_EB,31.84_EB,31.45_EB,31.06_EB,30.67_EB/)
QFLAME(3,20,0:10,3) = (/0.00_EB,22.37_EB,30.71_EB,33.44_EB,36.18_EB,36.88_EB,37.58_EB,38.28_EB,38.98_EB,38.67_EB,38.36_EB/)
QFLAME(3,20,11:20,3) = (/38.06_EB,37.75_EB,37.45_EB,37.14_EB,36.84_EB,36.53_EB,36.10_EB,35.66_EB,35.22_EB,34.79_EB/)
QFLAME(3,20,21:30,3) = (/34.36_EB,33.94_EB,33.51_EB,33.09_EB,32.66_EB,32.29_EB,31.93_EB,31.56_EB,31.19_EB,30.82_EB/)
QFLAME(3,20,0:10,4) = (/0.00_EB,24.15_EB,32.38_EB,35.13_EB,37.89_EB,38.65_EB,39.42_EB,40.19_EB,40.95_EB,40.74_EB,40.53_EB/)
QFLAME(3,20,11:20,4) = (/40.32_EB,40.11_EB,39.92_EB,39.72_EB,39.53_EB,39.33_EB,38.86_EB,38.39_EB,37.92_EB,37.45_EB/)
QFLAME(3,20,21:30,4) = (/37.67_EB,37.88_EB,38.10_EB,38.31_EB,38.53_EB,38.17_EB,37.82_EB,37.46_EB,37.11_EB,36.75_EB/)
QFLAME(3,20,0:10,5) = (/0.00_EB,26.01_EB,35.96_EB,39.99_EB,44.02_EB,45.45_EB,46.87_EB,48.30_EB,49.73_EB,49.80_EB,49.87_EB/)
QFLAME(3,20,11:20,5) = (/49.95_EB,50.02_EB,50.04_EB,50.06_EB,50.09_EB,50.11_EB,49.82_EB,49.54_EB,49.26_EB,48.97_EB/)
QFLAME(3,20,21:30,5) = (/49.16_EB,49.34_EB,49.53_EB,49.71_EB,49.89_EB,49.64_EB,49.39_EB,49.13_EB,48.88_EB,48.63_EB/)
QFLAME(3,20,0:10,6) = (/0.00_EB,27.85_EB,39.59_EB,44.64_EB,49.70_EB,51.90_EB,54.10_EB,56.30_EB,58.50_EB,59.07_EB,59.64_EB/)
QFLAME(3,20,11:20,6) = (/60.21_EB,60.78_EB,60.32_EB,59.87_EB,59.41_EB,58.96_EB,59.28_EB,59.59_EB,59.91_EB,60.23_EB/)
QFLAME(3,20,21:30,6) = (/60.74_EB,61.26_EB,61.78_EB,62.29_EB,62.81_EB,62.46_EB,62.11_EB,61.76_EB,61.41_EB,61.06_EB/)
QFLAME(4,0,0:10,1) = (/0.00_EB,19.40_EB,24.30_EB,24.89_EB,25.49_EB,25.17_EB,24.86_EB,24.55_EB,24.23_EB,23.55_EB,22.87_EB/)
QFLAME(4,0,11:20,1) = (/22.18_EB,21.50_EB,20.94_EB,20.39_EB,19.83_EB,19.27_EB,18.80_EB,18.33_EB,17.86_EB,17.39_EB/)
QFLAME(4,0,21:30,1) = (/16.98_EB,16.57_EB,16.15_EB,15.74_EB,15.33_EB,15.05_EB,14.77_EB,14.49_EB,14.21_EB,13.93_EB/)
QFLAME(4,0,0:10,2) = (/0.00_EB,21.27_EB,27.15_EB,28.10_EB,29.06_EB,28.92_EB,28.79_EB,28.65_EB,28.52_EB,27.87_EB,27.22_EB/)
QFLAME(4,0,11:20,2) = (/26.57_EB,25.91_EB,25.54_EB,25.17_EB,24.79_EB,24.42_EB,23.96_EB,23.49_EB,23.03_EB,22.57_EB/)
QFLAME(4,0,21:30,2) = (/22.11_EB,21.64_EB,21.18_EB,20.72_EB,20.26_EB,20.15_EB,20.04_EB,19.93_EB,19.81_EB,19.70_EB/)
QFLAME(4,0,0:10,3) = (/0.00_EB,23.19_EB,30.50_EB,32.42_EB,34.34_EB,34.69_EB,35.03_EB,35.38_EB,35.73_EB,35.22_EB,34.70_EB/)
QFLAME(4,0,11:20,3) = (/34.19_EB,33.68_EB,33.87_EB,34.06_EB,34.25_EB,34.44_EB,33.73_EB,33.02_EB,32.32_EB,31.61_EB/)
QFLAME(4,0,21:30,3) = (/31.07_EB,30.52_EB,29.98_EB,29.44_EB,28.89_EB,29.06_EB,29.23_EB,29.39_EB,29.56_EB,29.73_EB/)
QFLAME(4,0,0:10,4) = (/0.00_EB,25.01_EB,33.97_EB,36.96_EB,39.94_EB,40.91_EB,41.88_EB,42.86_EB,43.83_EB,43.61_EB,43.39_EB/)
QFLAME(4,0,11:20,4) = (/43.18_EB,42.96_EB,43.41_EB,43.86_EB,44.31_EB,44.75_EB,44.71_EB,44.66_EB,44.61_EB,44.57_EB/)
QFLAME(4,0,21:30,4) = (/43.38_EB,42.19_EB,41.00_EB,39.81_EB,38.62_EB,39.00_EB,39.38_EB,39.76_EB,40.13_EB,40.51_EB/)
QFLAME(4,0,0:10,5) = (/0.00_EB,26.72_EB,37.67_EB,41.96_EB,46.24_EB,47.97_EB,49.69_EB,51.41_EB,53.13_EB,53.34_EB,53.55_EB/)
QFLAME(4,0,11:20,5) = (/53.75_EB,53.96_EB,53.91_EB,53.86_EB,53.81_EB,53.76_EB,54.47_EB,55.18_EB,55.89_EB,56.60_EB/)
QFLAME(4,0,21:30,5) = (/56.39_EB,56.18_EB,55.98_EB,55.77_EB,55.56_EB,55.63_EB,55.70_EB,55.76_EB,55.83_EB,55.90_EB/)
QFLAME(4,0,0:10,6) = (/0.00_EB,28.50_EB,41.25_EB,47.09_EB,52.92_EB,55.38_EB,57.84_EB,60.31_EB,62.77_EB,63.64_EB,64.50_EB/)
QFLAME(4,0,11:20,6) = (/65.37_EB,66.24_EB,66.39_EB,66.54_EB,66.69_EB,66.84_EB,67.60_EB,68.36_EB,69.12_EB,69.88_EB/)
QFLAME(4,0,21:30,6) = (/69.80_EB,69.72_EB,69.64_EB,69.56_EB,69.48_EB,69.47_EB,69.47_EB,69.47_EB,69.47_EB,69.46_EB/)
QFLAME(4,1,0:10,1) = (/0.00_EB,19.58_EB,25.02_EB,25.90_EB,26.78_EB,26.59_EB,26.41_EB,26.22_EB,26.03_EB,25.40_EB,24.76_EB/)
QFLAME(4,1,11:20,1) = (/24.13_EB,23.50_EB,22.95_EB,22.40_EB,21.85_EB,21.31_EB,20.84_EB,20.38_EB,19.91_EB,19.45_EB/)
QFLAME(4,1,21:30,1) = (/19.05_EB,18.66_EB,18.26_EB,17.87_EB,17.48_EB,17.15_EB,16.81_EB,16.48_EB,16.15_EB,15.82_EB/)
QFLAME(4,1,0:10,2) = (/0.00_EB,21.24_EB,27.15_EB,28.13_EB,29.11_EB,28.94_EB,28.78_EB,28.62_EB,28.45_EB,27.85_EB,27.24_EB/)
QFLAME(4,1,11:20,2) = (/26.63_EB,26.03_EB,25.60_EB,25.18_EB,24.75_EB,24.33_EB,23.85_EB,23.37_EB,22.88_EB,22.40_EB/)
QFLAME(4,1,21:30,2) = (/21.99_EB,21.58_EB,21.16_EB,20.75_EB,20.34_EB,20.38_EB,20.43_EB,20.47_EB,20.51_EB,20.56_EB/)
QFLAME(4,1,0:10,3) = (/0.00_EB,23.08_EB,30.49_EB,32.43_EB,34.37_EB,34.70_EB,35.03_EB,35.36_EB,35.68_EB,35.18_EB,34.67_EB/)
QFLAME(4,1,11:20,3) = (/34.16_EB,33.65_EB,33.82_EB,33.99_EB,34.16_EB,34.33_EB,33.67_EB,33.01_EB,32.35_EB,31.69_EB/)
QFLAME(4,1,21:30,3) = (/31.12_EB,30.55_EB,29.98_EB,29.41_EB,28.84_EB,29.14_EB,29.44_EB,29.74_EB,30.04_EB,30.34_EB/)
QFLAME(4,1,0:10,4) = (/0.00_EB,24.94_EB,33.95_EB,36.95_EB,39.96_EB,40.93_EB,41.91_EB,42.88_EB,43.86_EB,43.61_EB,43.36_EB/)
QFLAME(4,1,11:20,4) = (/43.10_EB,42.85_EB,43.34_EB,43.83_EB,44.31_EB,44.80_EB,44.02_EB,43.23_EB,42.45_EB,41.67_EB/)
QFLAME(4,1,21:30,4) = (/41.10_EB,40.54_EB,39.98_EB,39.41_EB,38.85_EB,39.52_EB,40.19_EB,40.85_EB,41.52_EB,42.19_EB/)
QFLAME(4,1,0:10,5) = (/0.00_EB,26.68_EB,37.63_EB,41.93_EB,46.23_EB,47.93_EB,49.63_EB,51.33_EB,53.03_EB,53.22_EB,53.41_EB/)
QFLAME(4,1,11:20,5) = (/53.60_EB,53.78_EB,53.84_EB,53.89_EB,53.95_EB,54.00_EB,54.63_EB,55.26_EB,55.89_EB,56.52_EB/)
QFLAME(4,1,21:30,5) = (/56.11_EB,55.71_EB,55.31_EB,54.90_EB,54.50_EB,54.74_EB,54.98_EB,55.22_EB,55.45_EB,55.69_EB/)
QFLAME(4,1,0:10,6) = (/0.00_EB,28.45_EB,40.97_EB,46.90_EB,52.83_EB,55.31_EB,57.79_EB,60.27_EB,62.75_EB,63.66_EB,64.56_EB/)
QFLAME(4,1,11:20,6) = (/65.46_EB,66.37_EB,67.12_EB,67.86_EB,68.61_EB,69.35_EB,69.29_EB,69.22_EB,69.15_EB,69.08_EB/)
QFLAME(4,1,21:30,6) = (/68.87_EB,68.67_EB,68.46_EB,68.25_EB,68.05_EB,67.92_EB,67.80_EB,67.68_EB,67.55_EB,67.43_EB/)
QFLAME(4,2,0:10,1) = (/0.00_EB,19.90_EB,25.66_EB,26.78_EB,27.89_EB,27.82_EB,27.75_EB,27.68_EB,27.61_EB,27.02_EB,26.44_EB/)
QFLAME(4,2,11:20,1) = (/25.85_EB,25.27_EB,24.75_EB,24.23_EB,23.71_EB,23.20_EB,22.73_EB,22.26_EB,21.79_EB,21.32_EB/)
QFLAME(4,2,21:30,1) = (/20.90_EB,20.49_EB,20.08_EB,19.67_EB,19.26_EB,18.94_EB,18.63_EB,18.31_EB,18.00_EB,17.68_EB/)
QFLAME(4,2,0:10,2) = (/0.00_EB,21.18_EB,27.15_EB,28.17_EB,29.20_EB,29.06_EB,28.92_EB,28.78_EB,28.64_EB,28.02_EB,27.41_EB/)
QFLAME(4,2,11:20,2) = (/26.80_EB,26.18_EB,25.77_EB,25.37_EB,24.96_EB,24.56_EB,24.09_EB,23.63_EB,23.17_EB,22.71_EB/)
QFLAME(4,2,21:30,2) = (/22.28_EB,21.86_EB,21.44_EB,21.01_EB,20.59_EB,20.59_EB,20.58_EB,20.58_EB,20.58_EB,20.58_EB/)
QFLAME(4,2,0:10,3) = (/0.00_EB,23.16_EB,30.46_EB,32.43_EB,34.39_EB,34.70_EB,35.01_EB,35.32_EB,35.62_EB,35.10_EB,34.58_EB/)
QFLAME(4,2,11:20,3) = (/34.06_EB,33.53_EB,33.74_EB,33.95_EB,34.16_EB,34.37_EB,33.68_EB,32.98_EB,32.29_EB,31.59_EB/)
QFLAME(4,2,21:30,3) = (/30.95_EB,30.31_EB,29.68_EB,29.04_EB,28.40_EB,28.79_EB,29.18_EB,29.57_EB,29.96_EB,30.35_EB/)
QFLAME(4,2,0:10,4) = (/0.00_EB,24.91_EB,33.93_EB,36.95_EB,39.97_EB,40.94_EB,41.91_EB,42.87_EB,43.84_EB,43.58_EB,43.32_EB/)
QFLAME(4,2,11:20,4) = (/43.06_EB,42.80_EB,43.26_EB,43.73_EB,44.20_EB,44.67_EB,44.50_EB,44.34_EB,44.18_EB,44.02_EB/)
QFLAME(4,2,21:30,4) = (/43.52_EB,43.03_EB,42.53_EB,42.03_EB,41.54_EB,41.57_EB,41.60_EB,41.63_EB,41.66_EB,41.69_EB/)
QFLAME(4,2,0:10,5) = (/0.00_EB,26.63_EB,37.58_EB,41.89_EB,46.20_EB,47.91_EB,49.62_EB,51.33_EB,53.04_EB,53.18_EB,53.32_EB/)
QFLAME(4,2,11:20,5) = (/53.46_EB,53.60_EB,53.80_EB,53.99_EB,54.19_EB,54.38_EB,54.52_EB,54.65_EB,54.78_EB,54.92_EB/)
QFLAME(4,2,21:30,5) = (/55.16_EB,55.40_EB,55.64_EB,55.88_EB,56.12_EB,55.75_EB,55.37_EB,55.00_EB,54.63_EB,54.25_EB/)
QFLAME(4,2,0:10,6) = (/0.00_EB,28.44_EB,41.00_EB,46.86_EB,52.73_EB,55.23_EB,57.73_EB,60.24_EB,62.74_EB,63.59_EB,64.44_EB/)
QFLAME(4,2,11:20,6) = (/65.29_EB,66.14_EB,66.16_EB,66.18_EB,66.20_EB,66.22_EB,66.88_EB,67.55_EB,68.21_EB,68.88_EB/)
QFLAME(4,2,21:30,6) = (/68.43_EB,67.97_EB,67.52_EB,67.07_EB,66.62_EB,67.06_EB,67.50_EB,67.94_EB,68.38_EB,68.82_EB/)
QFLAME(4,3,0:10,1) = (/0.00_EB,20.12_EB,26.15_EB,27.45_EB,28.74_EB,28.76_EB,28.77_EB,28.79_EB,28.80_EB,28.25_EB,27.69_EB/)
QFLAME(4,3,11:20,1) = (/27.13_EB,26.58_EB,26.09_EB,25.60_EB,25.11_EB,24.62_EB,24.15_EB,23.68_EB,23.22_EB,22.75_EB/)
QFLAME(4,3,21:30,1) = (/22.34_EB,21.93_EB,21.51_EB,21.10_EB,20.69_EB,20.35_EB,20.01_EB,19.68_EB,19.34_EB,19.00_EB/)
QFLAME(4,3,0:10,2) = (/0.00_EB,21.16_EB,27.21_EB,28.42_EB,29.63_EB,29.60_EB,29.57_EB,29.54_EB,29.51_EB,28.93_EB,28.35_EB/)
QFLAME(4,3,11:20,2) = (/27.77_EB,27.20_EB,26.78_EB,26.36_EB,25.94_EB,25.52_EB,25.06_EB,24.60_EB,24.14_EB,23.68_EB/)
QFLAME(4,3,21:30,2) = (/23.26_EB,22.84_EB,22.42_EB,22.00_EB,21.58_EB,21.46_EB,21.35_EB,21.23_EB,21.11_EB,21.00_EB/)
QFLAME(4,3,0:10,3) = (/0.00_EB,23.12_EB,30.43_EB,32.41_EB,34.38_EB,34.69_EB,35.00_EB,35.31_EB,35.62_EB,35.09_EB,34.56_EB/)
QFLAME(4,3,11:20,3) = (/34.04_EB,33.51_EB,33.71_EB,33.91_EB,34.11_EB,34.30_EB,33.62_EB,32.94_EB,32.26_EB,31.58_EB/)
QFLAME(4,3,21:30,3) = (/30.95_EB,30.33_EB,29.70_EB,29.08_EB,28.45_EB,28.81_EB,29.18_EB,29.54_EB,29.91_EB,30.27_EB/)
QFLAME(4,3,0:10,4) = (/0.00_EB,24.88_EB,33.89_EB,36.92_EB,39.96_EB,40.92_EB,41.88_EB,42.85_EB,43.81_EB,43.55_EB,43.29_EB/)
QFLAME(4,3,11:20,4) = (/43.03_EB,42.77_EB,43.23_EB,43.69_EB,44.16_EB,44.62_EB,44.45_EB,44.28_EB,44.10_EB,43.93_EB/)
QFLAME(4,3,21:30,4) = (/43.40_EB,42.87_EB,42.34_EB,41.81_EB,41.28_EB,41.37_EB,41.46_EB,41.55_EB,41.64_EB,41.73_EB/)
QFLAME(4,3,0:10,5) = (/0.00_EB,26.61_EB,37.53_EB,41.85_EB,46.17_EB,47.89_EB,49.60_EB,51.32_EB,53.04_EB,53.17_EB,53.30_EB/)
QFLAME(4,3,11:20,5) = (/53.44_EB,53.57_EB,53.74_EB,53.90_EB,54.07_EB,54.24_EB,54.41_EB,54.58_EB,54.75_EB,54.92_EB/)
QFLAME(4,3,21:30,5) = (/55.07_EB,55.23_EB,55.39_EB,55.54_EB,55.70_EB,55.41_EB,55.12_EB,54.83_EB,54.54_EB,54.25_EB/)
QFLAME(4,3,0:10,6) = (/0.00_EB,28.44_EB,41.02_EB,46.83_EB,52.65_EB,55.16_EB,57.67_EB,60.18_EB,62.69_EB,63.53_EB,64.37_EB/)
QFLAME(4,3,11:20,6) = (/65.21_EB,66.05_EB,66.12_EB,66.19_EB,66.27_EB,66.34_EB,66.72_EB,67.10_EB,67.48_EB,67.86_EB/)
QFLAME(4,3,21:30,6) = (/67.70_EB,67.55_EB,67.40_EB,67.24_EB,67.09_EB,67.41_EB,67.73_EB,68.05_EB,68.37_EB,68.69_EB/)
QFLAME(4,4,0:10,1) = (/0.00_EB,20.33_EB,26.64_EB,28.12_EB,29.59_EB,29.69_EB,29.80_EB,29.90_EB,30.00_EB,29.47_EB,28.94_EB/)
QFLAME(4,4,11:20,1) = (/28.41_EB,27.89_EB,27.42_EB,26.96_EB,26.50_EB,26.03_EB,25.57_EB,25.11_EB,24.65_EB,24.19_EB/)
QFLAME(4,4,21:30,1) = (/23.78_EB,23.36_EB,22.94_EB,22.53_EB,22.11_EB,21.76_EB,21.40_EB,21.04_EB,20.68_EB,20.32_EB/)
QFLAME(4,4,0:10,2) = (/0.00_EB,21.13_EB,27.27_EB,28.66_EB,30.05_EB,30.14_EB,30.22_EB,30.30_EB,30.39_EB,29.84_EB,29.30_EB/)
QFLAME(4,4,11:20,2) = (/28.75_EB,28.21_EB,27.78_EB,27.35_EB,26.92_EB,26.49_EB,26.03_EB,25.57_EB,25.11_EB,24.64_EB/)
QFLAME(4,4,21:30,2) = (/24.23_EB,23.81_EB,23.40_EB,22.98_EB,22.57_EB,22.34_EB,22.11_EB,21.88_EB,21.65_EB,21.42_EB/)
QFLAME(4,4,0:10,3) = (/0.00_EB,23.07_EB,30.40_EB,32.39_EB,34.37_EB,34.68_EB,34.99_EB,35.30_EB,35.61_EB,35.08_EB,34.55_EB/)
QFLAME(4,4,11:20,3) = (/34.02_EB,33.49_EB,33.68_EB,33.86_EB,34.05_EB,34.24_EB,33.57_EB,32.90_EB,32.24_EB,31.57_EB/)
QFLAME(4,4,21:30,3) = (/30.96_EB,30.34_EB,29.73_EB,29.11_EB,28.50_EB,28.84_EB,29.18_EB,29.52_EB,29.86_EB,30.20_EB/)
QFLAME(4,4,0:10,4) = (/0.00_EB,24.84_EB,33.85_EB,36.90_EB,39.94_EB,40.90_EB,41.86_EB,42.82_EB,43.78_EB,43.52_EB,43.26_EB/)
QFLAME(4,4,11:20,4) = (/43.00_EB,42.74_EB,43.20_EB,43.66_EB,44.11_EB,44.57_EB,44.39_EB,44.21_EB,44.03_EB,43.85_EB/)
QFLAME(4,4,21:30,4) = (/43.29_EB,42.72_EB,42.15_EB,41.59_EB,41.02_EB,41.18_EB,41.33_EB,41.48_EB,41.63_EB,41.78_EB/)
QFLAME(4,4,0:10,5) = (/0.00_EB,26.60_EB,37.49_EB,41.81_EB,46.14_EB,47.86_EB,49.58_EB,51.31_EB,53.03_EB,53.16_EB,53.29_EB/)
QFLAME(4,4,11:20,5) = (/53.41_EB,53.54_EB,53.68_EB,53.82_EB,53.95_EB,54.09_EB,54.30_EB,54.51_EB,54.71_EB,54.92_EB/)
QFLAME(4,4,21:30,5) = (/54.99_EB,55.06_EB,55.13_EB,55.20_EB,55.27_EB,55.06_EB,54.86_EB,54.65_EB,54.45_EB,54.24_EB/)
QFLAME(4,4,0:10,6) = (/0.00_EB,28.43_EB,41.04_EB,46.81_EB,52.57_EB,55.09_EB,57.61_EB,60.13_EB,62.65_EB,63.48_EB,64.30_EB/)
QFLAME(4,4,11:20,6) = (/65.13_EB,65.96_EB,66.08_EB,66.21_EB,66.33_EB,66.46_EB,66.55_EB,66.65_EB,66.74_EB,66.84_EB/)
QFLAME(4,4,21:30,6) = (/66.98_EB,67.12_EB,67.27_EB,67.41_EB,67.56_EB,67.76_EB,67.96_EB,68.16_EB,68.36_EB,68.56_EB/)
QFLAME(4,5,0:10,1) = (/0.00_EB,20.55_EB,27.14_EB,28.79_EB,30.45_EB,30.63_EB,30.82_EB,31.01_EB,31.19_EB,30.69_EB,30.19_EB/)
QFLAME(4,5,11:20,1) = (/29.70_EB,29.20_EB,28.76_EB,28.32_EB,27.89_EB,27.45_EB,27.00_EB,26.54_EB,26.08_EB,25.63_EB/)
QFLAME(4,5,21:30,1) = (/25.21_EB,24.79_EB,24.38_EB,23.96_EB,23.54_EB,23.16_EB,22.78_EB,22.40_EB,22.02_EB,21.64_EB/)
QFLAME(4,5,0:10,2) = (/0.00_EB,21.11_EB,27.33_EB,28.91_EB,30.48_EB,30.68_EB,30.87_EB,31.07_EB,31.26_EB,30.75_EB,30.24_EB/)
QFLAME(4,5,11:20,2) = (/29.73_EB,29.23_EB,28.79_EB,28.34_EB,27.90_EB,27.46_EB,27.00_EB,26.54_EB,26.07_EB,25.61_EB/)
QFLAME(4,5,21:30,2) = (/25.20_EB,24.79_EB,24.38_EB,23.97_EB,23.56_EB,23.21_EB,22.87_EB,22.53_EB,22.18_EB,21.84_EB/)
QFLAME(4,5,0:10,3) = (/0.00_EB,23.03_EB,30.37_EB,32.37_EB,34.36_EB,34.67_EB,34.98_EB,35.30_EB,35.61_EB,35.07_EB,34.54_EB/)
QFLAME(4,5,11:20,3) = (/34.00_EB,33.47_EB,33.64_EB,33.82_EB,33.99_EB,34.17_EB,33.52_EB,32.87_EB,32.21_EB,31.56_EB/)
QFLAME(4,5,21:30,3) = (/30.96_EB,30.36_EB,29.76_EB,29.15_EB,28.55_EB,28.87_EB,29.18_EB,29.50_EB,29.81_EB,30.12_EB/)
QFLAME(4,5,0:10,4) = (/0.00_EB,24.81_EB,33.82_EB,36.87_EB,39.93_EB,40.89_EB,41.84_EB,42.80_EB,43.76_EB,43.50_EB,43.24_EB/)
QFLAME(4,5,11:20,4) = (/42.98_EB,42.72_EB,43.17_EB,43.62_EB,44.07_EB,44.52_EB,44.33_EB,44.14_EB,43.96_EB,43.77_EB/)
QFLAME(4,5,21:30,4) = (/43.17_EB,42.57_EB,41.97_EB,41.37_EB,40.77_EB,40.98_EB,41.19_EB,41.40_EB,41.62_EB,41.83_EB/)
QFLAME(4,5,0:10,5) = (/0.00_EB,26.58_EB,37.45_EB,41.78_EB,46.11_EB,47.84_EB,49.57_EB,51.30_EB,53.03_EB,53.15_EB,53.27_EB/)
QFLAME(4,5,11:20,5) = (/53.39_EB,53.51_EB,53.62_EB,53.73_EB,53.84_EB,53.95_EB,54.19_EB,54.43_EB,54.68_EB,54.92_EB/)
QFLAME(4,5,21:30,5) = (/54.91_EB,54.89_EB,54.87_EB,54.86_EB,54.84_EB,54.72_EB,54.60_EB,54.48_EB,54.36_EB,54.24_EB/)
QFLAME(4,5,0:10,6) = (/0.00_EB,28.42_EB,41.06_EB,46.78_EB,52.49_EB,55.02_EB,57.55_EB,60.08_EB,62.60_EB,63.42_EB,64.24_EB/)
QFLAME(4,5,11:20,6) = (/65.05_EB,65.87_EB,66.04_EB,66.22_EB,66.40_EB,66.57_EB,66.38_EB,66.19_EB,66.00_EB,65.82_EB/)
QFLAME(4,5,21:30,6) = (/66.26_EB,66.70_EB,67.14_EB,67.58_EB,68.03_EB,68.11_EB,68.19_EB,68.27_EB,68.35_EB,68.43_EB/)
QFLAME(4,6,0:10,1) = (/0.00_EB,20.70_EB,27.46_EB,29.23_EB,31.00_EB,31.25_EB,31.50_EB,31.74_EB,31.99_EB,31.52_EB,31.05_EB/)
QFLAME(4,6,11:20,1) = (/30.57_EB,30.10_EB,29.68_EB,29.26_EB,28.84_EB,28.42_EB,27.96_EB,27.51_EB,27.06_EB,26.60_EB/)
QFLAME(4,6,21:30,1) = (/26.19_EB,25.78_EB,25.36_EB,24.95_EB,24.54_EB,24.15_EB,23.77_EB,23.39_EB,23.00_EB,22.62_EB/)
QFLAME(4,6,0:10,2) = (/0.00_EB,21.14_EB,27.62_EB,29.33_EB,31.03_EB,31.29_EB,31.54_EB,31.79_EB,32.05_EB,31.57_EB,31.08_EB/)
QFLAME(4,6,11:20,2) = (/30.60_EB,30.12_EB,29.70_EB,29.27_EB,28.85_EB,28.42_EB,27.97_EB,27.51_EB,27.05_EB,26.59_EB/)
QFLAME(4,6,21:30,2) = (/26.18_EB,25.76_EB,25.35_EB,24.94_EB,24.52_EB,24.17_EB,23.82_EB,23.47_EB,23.12_EB,22.76_EB/)
QFLAME(4,6,0:10,3) = (/0.00_EB,22.98_EB,30.34_EB,32.35_EB,34.36_EB,34.70_EB,35.04_EB,35.38_EB,35.73_EB,35.23_EB,34.74_EB/)
QFLAME(4,6,11:20,3) = (/34.25_EB,33.76_EB,33.85_EB,33.93_EB,34.02_EB,34.10_EB,33.49_EB,32.87_EB,32.26_EB,31.64_EB/)
QFLAME(4,6,21:30,3) = (/31.07_EB,30.50_EB,29.94_EB,29.37_EB,28.80_EB,29.06_EB,29.32_EB,29.58_EB,29.84_EB,30.10_EB/)
QFLAME(4,6,0:10,4) = (/0.00_EB,24.78_EB,33.78_EB,36.84_EB,39.90_EB,40.86_EB,41.82_EB,42.78_EB,43.74_EB,43.48_EB,43.22_EB/)
QFLAME(4,6,11:20,4) = (/42.96_EB,42.70_EB,43.15_EB,43.59_EB,44.04_EB,44.49_EB,44.24_EB,43.98_EB,43.73_EB,43.48_EB/)
QFLAME(4,6,21:30,4) = (/42.86_EB,42.23_EB,41.61_EB,40.99_EB,40.36_EB,40.64_EB,40.92_EB,41.20_EB,41.47_EB,41.75_EB/)
QFLAME(4,6,0:10,5) = (/0.00_EB,26.56_EB,37.40_EB,41.73_EB,46.06_EB,47.79_EB,49.53_EB,51.26_EB,52.99_EB,53.12_EB,53.24_EB/)
QFLAME(4,6,11:20,5) = (/53.36_EB,53.48_EB,53.53_EB,53.59_EB,53.65_EB,53.71_EB,53.92_EB,54.14_EB,54.35_EB,54.57_EB/)
QFLAME(4,6,21:30,5) = (/54.54_EB,54.51_EB,54.48_EB,54.46_EB,54.43_EB,54.37_EB,54.30_EB,54.24_EB,54.18_EB,54.12_EB/)
QFLAME(4,6,0:10,6) = (/0.00_EB,28.39_EB,41.01_EB,46.72_EB,52.42_EB,54.96_EB,57.49_EB,60.02_EB,62.55_EB,63.37_EB,64.18_EB/)
QFLAME(4,6,11:20,6) = (/64.99_EB,65.81_EB,66.04_EB,66.27_EB,66.50_EB,66.73_EB,66.59_EB,66.45_EB,66.31_EB,66.17_EB/)
QFLAME(4,6,21:30,6) = (/66.51_EB,66.84_EB,67.18_EB,67.51_EB,67.84_EB,67.89_EB,67.93_EB,67.97_EB,68.02_EB,68.06_EB/)
QFLAME(4,7,0:10,1) = (/0.00_EB,20.85_EB,27.78_EB,29.67_EB,31.56_EB,31.87_EB,32.18_EB,32.48_EB,32.79_EB,32.34_EB,31.90_EB/)
QFLAME(4,7,11:20,1) = (/31.45_EB,31.00_EB,30.60_EB,30.19_EB,29.79_EB,29.38_EB,28.93_EB,28.48_EB,28.03_EB,27.58_EB/)
QFLAME(4,7,21:30,1) = (/27.17_EB,26.76_EB,26.35_EB,25.94_EB,25.54_EB,25.15_EB,24.76_EB,24.37_EB,23.99_EB,23.60_EB/)
QFLAME(4,7,0:10,2) = (/0.00_EB,21.18_EB,27.90_EB,29.74_EB,31.58_EB,31.90_EB,32.21_EB,32.52_EB,32.83_EB,32.38_EB,31.93_EB/)
QFLAME(4,7,11:20,2) = (/31.47_EB,31.02_EB,30.61_EB,30.20_EB,29.79_EB,29.39_EB,28.93_EB,28.48_EB,28.03_EB,27.57_EB/)
QFLAME(4,7,21:30,2) = (/27.16_EB,26.74_EB,26.32_EB,25.91_EB,25.49_EB,25.13_EB,24.77_EB,24.41_EB,24.05_EB,23.69_EB/)
QFLAME(4,7,0:10,3) = (/0.00_EB,22.93_EB,30.30_EB,32.33_EB,34.36_EB,34.73_EB,35.10_EB,35.47_EB,35.84_EB,35.40_EB,34.95_EB/)
QFLAME(4,7,11:20,3) = (/34.50_EB,34.06_EB,34.05_EB,34.05_EB,34.04_EB,34.04_EB,33.46_EB,32.88_EB,32.30_EB,31.72_EB/)
QFLAME(4,7,21:30,3) = (/31.19_EB,30.65_EB,30.12_EB,29.58_EB,29.05_EB,29.25_EB,29.46_EB,29.66_EB,29.87_EB,30.07_EB/)
QFLAME(4,7,0:10,4) = (/0.00_EB,24.75_EB,33.73_EB,36.80_EB,39.87_EB,40.83_EB,41.79_EB,42.76_EB,43.72_EB,43.46_EB,43.20_EB/)
QFLAME(4,7,11:20,4) = (/42.94_EB,42.68_EB,43.13_EB,43.57_EB,44.01_EB,44.45_EB,44.14_EB,43.82_EB,43.51_EB,43.19_EB/)
QFLAME(4,7,21:30,4) = (/42.55_EB,41.90_EB,41.25_EB,40.61_EB,39.96_EB,40.30_EB,40.65_EB,40.99_EB,41.33_EB,41.67_EB/)
QFLAME(4,7,0:10,5) = (/0.00_EB,26.53_EB,37.35_EB,41.68_EB,46.01_EB,47.75_EB,49.49_EB,51.22_EB,52.96_EB,53.08_EB,53.20_EB/)
QFLAME(4,7,11:20,5) = (/53.32_EB,53.44_EB,53.45_EB,53.46_EB,53.46_EB,53.47_EB,53.66_EB,53.84_EB,54.03_EB,54.21_EB/)
QFLAME(4,7,21:30,5) = (/54.17_EB,54.13_EB,54.09_EB,54.05_EB,54.01_EB,54.01_EB,54.01_EB,54.00_EB,54.00_EB,54.00_EB/)
QFLAME(4,7,0:10,6) = (/0.00_EB,28.36_EB,40.97_EB,46.66_EB,52.35_EB,54.89_EB,57.43_EB,59.97_EB,62.50_EB,63.31_EB,64.13_EB/)
QFLAME(4,7,11:20,6) = (/64.94_EB,65.75_EB,66.03_EB,66.32_EB,66.60_EB,66.89_EB,66.80_EB,66.71_EB,66.62_EB,66.53_EB/)
QFLAME(4,7,21:30,6) = (/66.75_EB,66.98_EB,67.21_EB,67.44_EB,67.66_EB,67.67_EB,67.67_EB,67.68_EB,67.68_EB,67.69_EB/)
QFLAME(4,8,0:10,1) = (/0.00_EB,21.00_EB,28.11_EB,30.11_EB,32.12_EB,32.49_EB,32.86_EB,33.22_EB,33.59_EB,33.17_EB,32.75_EB/)
QFLAME(4,8,11:20,1) = (/32.33_EB,31.91_EB,31.52_EB,31.13_EB,30.74_EB,30.35_EB,29.90_EB,29.45_EB,29.01_EB,28.56_EB/)
QFLAME(4,8,21:30,1) = (/28.15_EB,27.75_EB,27.34_EB,26.94_EB,26.53_EB,26.14_EB,25.75_EB,25.36_EB,24.97_EB,24.58_EB/)
QFLAME(4,8,0:10,2) = (/0.00_EB,21.22_EB,28.18_EB,30.16_EB,32.13_EB,32.51_EB,32.88_EB,33.25_EB,33.62_EB,33.19_EB,32.77_EB/)
QFLAME(4,8,11:20,2) = (/32.34_EB,31.91_EB,31.52_EB,31.13_EB,30.74_EB,30.35_EB,29.90_EB,29.45_EB,29.00_EB,28.56_EB/)
QFLAME(4,8,21:30,2) = (/28.14_EB,27.71_EB,27.29_EB,26.87_EB,26.45_EB,26.09_EB,25.72_EB,25.35_EB,24.99_EB,24.62_EB/)
QFLAME(4,8,0:10,3) = (/0.00_EB,22.89_EB,30.26_EB,32.31_EB,34.35_EB,34.76_EB,35.16_EB,35.56_EB,35.96_EB,35.56_EB,35.16_EB/)
QFLAME(4,8,11:20,3) = (/34.75_EB,34.35_EB,34.26_EB,34.16_EB,34.07_EB,33.97_EB,33.43_EB,32.89_EB,32.35_EB,31.80_EB/)
QFLAME(4,8,21:30,3) = (/31.30_EB,30.80_EB,30.30_EB,29.80_EB,29.30_EB,29.45_EB,29.59_EB,29.74_EB,29.89_EB,30.04_EB/)
QFLAME(4,8,0:10,4) = (/0.00_EB,24.72_EB,33.69_EB,36.76_EB,39.83_EB,40.80_EB,41.77_EB,42.74_EB,43.70_EB,43.44_EB,43.19_EB/)
QFLAME(4,8,11:20,4) = (/42.93_EB,42.67_EB,43.11_EB,43.54_EB,43.98_EB,44.42_EB,44.04_EB,43.66_EB,43.28_EB,42.90_EB/)
QFLAME(4,8,21:30,4) = (/42.23_EB,41.57_EB,40.90_EB,40.23_EB,39.56_EB,39.97_EB,40.37_EB,40.78_EB,41.19_EB,41.59_EB/)
QFLAME(4,8,0:10,5) = (/0.00_EB,26.51_EB,37.30_EB,41.63_EB,45.96_EB,47.70_EB,49.44_EB,51.19_EB,52.93_EB,53.05_EB,53.17_EB/)
QFLAME(4,8,11:20,5) = (/53.29_EB,53.41_EB,53.36_EB,53.32_EB,53.28_EB,53.23_EB,53.39_EB,53.55_EB,53.70_EB,53.86_EB/)
QFLAME(4,8,21:30,5) = (/53.81_EB,53.76_EB,53.70_EB,53.65_EB,53.60_EB,53.65_EB,53.71_EB,53.77_EB,53.82_EB,53.88_EB/)
QFLAME(4,8,0:10,6) = (/0.00_EB,28.33_EB,40.92_EB,46.61_EB,52.29_EB,54.83_EB,57.37_EB,59.91_EB,62.45_EB,63.26_EB,64.07_EB/)
QFLAME(4,8,11:20,6) = (/64.88_EB,65.69_EB,66.03_EB,66.37_EB,66.70_EB,67.04_EB,67.00_EB,66.96_EB,66.92_EB,66.88_EB/)
QFLAME(4,8,21:30,6) = (/67.00_EB,67.12_EB,67.24_EB,67.36_EB,67.48_EB,67.45_EB,67.42_EB,67.38_EB,67.35_EB,67.32_EB/)
QFLAME(4,9,0:10,1) = (/0.00_EB,21.16_EB,28.43_EB,30.55_EB,32.68_EB,33.11_EB,33.53_EB,33.96_EB,34.39_EB,33.99_EB,33.60_EB/)
QFLAME(4,9,11:20,1) = (/33.20_EB,32.81_EB,32.44_EB,32.06_EB,31.69_EB,31.31_EB,30.87_EB,30.42_EB,29.98_EB,29.53_EB/)
QFLAME(4,9,21:30,1) = (/29.13_EB,28.73_EB,28.33_EB,27.93_EB,27.53_EB,27.14_EB,26.74_EB,26.35_EB,25.96_EB,25.56_EB/)
QFLAME(4,9,0:10,2) = (/0.00_EB,21.25_EB,28.46_EB,30.57_EB,32.68_EB,33.11_EB,33.54_EB,33.97_EB,34.40_EB,34.00_EB,33.61_EB/)
QFLAME(4,9,11:20,2) = (/33.21_EB,32.81_EB,32.43_EB,32.06_EB,31.68_EB,31.31_EB,30.86_EB,30.42_EB,29.98_EB,29.54_EB/)
QFLAME(4,9,21:30,2) = (/29.11_EB,28.69_EB,28.27_EB,27.84_EB,27.42_EB,27.04_EB,26.67_EB,26.30_EB,25.92_EB,25.55_EB/)
QFLAME(4,9,0:10,3) = (/0.00_EB,22.84_EB,30.23_EB,32.29_EB,34.35_EB,34.78_EB,35.21_EB,35.64_EB,36.07_EB,35.72_EB,35.36_EB/)
QFLAME(4,9,11:20,3) = (/35.01_EB,34.65_EB,34.46_EB,34.28_EB,34.09_EB,33.91_EB,33.40_EB,32.90_EB,32.39_EB,31.88_EB/)
QFLAME(4,9,21:30,3) = (/31.42_EB,30.95_EB,30.48_EB,30.01_EB,29.54_EB,29.64_EB,29.73_EB,29.83_EB,29.92_EB,30.02_EB/)
QFLAME(4,9,0:10,4) = (/0.00_EB,24.69_EB,33.65_EB,36.73_EB,39.80_EB,40.77_EB,41.74_EB,42.71_EB,43.68_EB,43.43_EB,43.17_EB/)
QFLAME(4,9,11:20,4) = (/42.91_EB,42.65_EB,43.09_EB,43.52_EB,43.95_EB,44.38_EB,43.94_EB,43.50_EB,43.06_EB,42.61_EB/)
QFLAME(4,9,21:30,4) = (/41.92_EB,41.23_EB,40.54_EB,39.85_EB,39.16_EB,39.63_EB,40.10_EB,40.57_EB,41.04_EB,41.51_EB/)
QFLAME(4,9,0:10,5) = (/0.00_EB,26.48_EB,37.25_EB,41.58_EB,45.91_EB,47.66_EB,49.40_EB,51.15_EB,52.90_EB,53.02_EB,53.14_EB/)
QFLAME(4,9,11:20,5) = (/53.25_EB,53.37_EB,53.28_EB,53.18_EB,53.09_EB,52.99_EB,53.12_EB,53.25_EB,53.38_EB,53.51_EB/)
QFLAME(4,9,21:30,5) = (/53.44_EB,53.38_EB,53.31_EB,53.25_EB,53.18_EB,53.30_EB,53.41_EB,53.53_EB,53.64_EB,53.76_EB/)
QFLAME(4,9,0:10,6) = (/0.00_EB,28.30_EB,40.88_EB,46.55_EB,52.22_EB,54.76_EB,57.31_EB,59.86_EB,62.40_EB,63.21_EB,64.02_EB/)
QFLAME(4,9,11:20,6) = (/64.83_EB,65.63_EB,66.02_EB,66.42_EB,66.81_EB,67.20_EB,67.21_EB,67.22_EB,67.23_EB,67.24_EB/)
QFLAME(4,9,21:30,6) = (/67.25_EB,67.26_EB,67.28_EB,67.29_EB,67.30_EB,67.23_EB,67.16_EB,67.09_EB,67.02_EB,66.95_EB/)
QFLAME(4,10,0:10,1) = (/0.00_EB,21.31_EB,28.75_EB,31.00_EB,33.24_EB,33.73_EB,34.21_EB,34.70_EB,35.19_EB,34.82_EB,34.45_EB/)
QFLAME(4,10,11:20,1) = (/34.08_EB,33.71_EB,33.35_EB,33.00_EB,32.64_EB,32.28_EB,31.84_EB,31.39_EB,30.95_EB,30.51_EB/)
QFLAME(4,10,21:30,1) = (/30.11_EB,29.72_EB,29.32_EB,28.92_EB,28.53_EB,28.13_EB,27.73_EB,27.34_EB,26.94_EB,26.54_EB/)
QFLAME(4,10,0:10,2) = (/0.00_EB,21.29_EB,28.74_EB,30.99_EB,33.23_EB,33.72_EB,34.21_EB,34.70_EB,35.19_EB,34.82_EB,34.45_EB/)
QFLAME(4,10,11:20,2) = (/34.08_EB,33.71_EB,33.35_EB,32.99_EB,32.63_EB,32.27_EB,31.83_EB,31.39_EB,30.96_EB,30.52_EB/)
QFLAME(4,10,21:30,2) = (/30.09_EB,29.66_EB,29.24_EB,28.81_EB,28.38_EB,28.00_EB,27.62_EB,27.24_EB,26.86_EB,26.47_EB/)
QFLAME(4,10,0:10,3) = (/0.00_EB,22.79_EB,30.19_EB,32.27_EB,34.35_EB,34.81_EB,35.27_EB,35.73_EB,36.19_EB,35.88_EB,35.57_EB/)
QFLAME(4,10,11:20,3) = (/35.26_EB,34.95_EB,34.67_EB,34.39_EB,34.12_EB,33.84_EB,33.37_EB,32.90_EB,32.43_EB,31.97_EB/)
QFLAME(4,10,21:30,3) = (/31.53_EB,31.10_EB,30.66_EB,30.23_EB,29.79_EB,29.83_EB,29.87_EB,29.91_EB,29.95_EB,29.99_EB/)
QFLAME(4,10,0:10,4) = (/0.00_EB,24.66_EB,33.60_EB,36.69_EB,39.77_EB,40.75_EB,41.72_EB,42.69_EB,43.67_EB,43.41_EB,43.15_EB/)
QFLAME(4,10,11:20,4) = (/42.89_EB,42.64_EB,43.07_EB,43.49_EB,43.92_EB,44.35_EB,43.84_EB,43.34_EB,42.83_EB,42.33_EB/)
QFLAME(4,10,21:30,4) = (/41.61_EB,40.90_EB,40.18_EB,39.47_EB,38.75_EB,39.29_EB,39.83_EB,40.36_EB,40.90_EB,41.44_EB/)
QFLAME(4,10,0:10,5) = (/0.00_EB,26.46_EB,37.20_EB,41.53_EB,45.86_EB,47.61_EB,49.36_EB,51.11_EB,52.87_EB,52.98_EB,53.10_EB/)
QFLAME(4,10,11:20,5) = (/53.22_EB,53.34_EB,53.19_EB,53.05_EB,52.90_EB,52.75_EB,52.85_EB,52.95_EB,53.05_EB,53.15_EB/)
QFLAME(4,10,21:30,5) = (/53.08_EB,53.00_EB,52.92_EB,52.85_EB,52.77_EB,52.94_EB,53.12_EB,53.29_EB,53.46_EB,53.64_EB/)
QFLAME(4,10,0:10,6) = (/0.00_EB,28.27_EB,40.83_EB,46.49_EB,52.15_EB,54.70_EB,57.25_EB,59.80_EB,62.35_EB,63.16_EB,63.96_EB/)
QFLAME(4,10,11:20,6) = (/64.77_EB,65.58_EB,66.02_EB,66.47_EB,66.91_EB,67.36_EB,67.42_EB,67.47_EB,67.53_EB,67.59_EB/)
QFLAME(4,10,21:30,6) = (/67.50_EB,67.40_EB,67.31_EB,67.21_EB,67.12_EB,67.01_EB,66.90_EB,66.79_EB,66.69_EB,66.58_EB/)
QFLAME(4,11,0:10,1) = (/0.00_EB,21.38_EB,28.93_EB,31.24_EB,33.54_EB,34.07_EB,34.59_EB,35.11_EB,35.63_EB,35.28_EB,34.94_EB/)
QFLAME(4,11,11:20,1) = (/34.59_EB,34.24_EB,33.89_EB,33.54_EB,33.18_EB,32.83_EB,32.40_EB,31.96_EB,31.52_EB,31.09_EB/)
QFLAME(4,11,21:30,1) = (/30.69_EB,30.30_EB,29.90_EB,29.50_EB,29.11_EB,28.72_EB,28.33_EB,27.93_EB,27.54_EB,27.15_EB/)
QFLAME(4,11,0:10,2) = (/0.00_EB,21.37_EB,28.92_EB,31.23_EB,33.54_EB,34.06_EB,34.59_EB,35.11_EB,35.63_EB,35.28_EB,34.93_EB/)
QFLAME(4,11,11:20,2) = (/34.58_EB,34.23_EB,33.88_EB,33.53_EB,33.18_EB,32.82_EB,32.39_EB,31.96_EB,31.53_EB,31.10_EB/)
QFLAME(4,11,21:30,2) = (/30.67_EB,30.25_EB,29.83_EB,29.41_EB,28.98_EB,28.61_EB,28.23_EB,27.85_EB,27.47_EB,27.09_EB/)
QFLAME(4,11,0:10,3) = (/0.00_EB,22.76_EB,30.23_EB,32.39_EB,34.55_EB,35.05_EB,35.55_EB,36.04_EB,36.54_EB,36.25_EB,35.96_EB/)
QFLAME(4,11,11:20,3) = (/35.66_EB,35.37_EB,35.09_EB,34.80_EB,34.52_EB,34.24_EB,33.78_EB,33.32_EB,32.86_EB,32.40_EB/)
QFLAME(4,11,21:30,3) = (/31.97_EB,31.54_EB,31.12_EB,30.69_EB,30.26_EB,30.26_EB,30.26_EB,30.26_EB,30.25_EB,30.25_EB/)
QFLAME(4,11,0:10,4) = (/0.00_EB,24.62_EB,33.56_EB,36.65_EB,39.74_EB,40.71_EB,41.69_EB,42.67_EB,43.64_EB,43.39_EB,43.14_EB/)
QFLAME(4,11,11:20,4) = (/42.89_EB,42.64_EB,43.05_EB,43.45_EB,43.85_EB,44.25_EB,43.76_EB,43.27_EB,42.77_EB,42.28_EB/)
QFLAME(4,11,21:30,4) = (/41.58_EB,40.88_EB,40.19_EB,39.49_EB,38.79_EB,39.25_EB,39.70_EB,40.16_EB,40.61_EB,41.06_EB/)
QFLAME(4,11,0:10,5) = (/0.00_EB,26.43_EB,37.15_EB,41.48_EB,45.82_EB,47.57_EB,49.33_EB,51.08_EB,52.83_EB,52.95_EB,53.07_EB/)
QFLAME(4,11,11:20,5) = (/53.20_EB,53.32_EB,53.14_EB,52.96_EB,52.78_EB,52.60_EB,52.70_EB,52.80_EB,52.90_EB,53.00_EB/)
QFLAME(4,11,21:30,5) = (/52.90_EB,52.81_EB,52.71_EB,52.61_EB,52.52_EB,52.74_EB,52.96_EB,53.19_EB,53.41_EB,53.63_EB/)
QFLAME(4,11,0:10,6) = (/0.00_EB,28.25_EB,40.73_EB,46.41_EB,52.09_EB,54.63_EB,57.18_EB,59.72_EB,62.27_EB,63.07_EB,63.88_EB/)
QFLAME(4,11,11:20,6) = (/64.69_EB,65.49_EB,65.87_EB,66.24_EB,66.61_EB,66.99_EB,67.02_EB,67.05_EB,67.09_EB,67.12_EB/)
QFLAME(4,11,21:30,6) = (/67.08_EB,67.03_EB,66.98_EB,66.93_EB,66.88_EB,66.71_EB,66.54_EB,66.37_EB,66.20_EB,66.03_EB/)
QFLAME(4,12,0:10,1) = (/0.00_EB,21.46_EB,29.10_EB,31.48_EB,33.85_EB,34.41_EB,34.96_EB,35.52_EB,36.08_EB,35.75_EB,35.42_EB/)
QFLAME(4,12,11:20,1) = (/35.10_EB,34.77_EB,34.42_EB,34.08_EB,33.73_EB,33.39_EB,32.96_EB,32.53_EB,32.10_EB,31.67_EB/)
QFLAME(4,12,21:30,1) = (/31.27_EB,30.87_EB,30.48_EB,30.08_EB,29.69_EB,29.30_EB,28.92_EB,28.53_EB,28.15_EB,27.76_EB/)
QFLAME(4,12,0:10,2) = (/0.00_EB,21.44_EB,29.09_EB,31.47_EB,33.84_EB,34.40_EB,34.96_EB,35.52_EB,36.07_EB,35.75_EB,35.42_EB/)
QFLAME(4,12,11:20,2) = (/35.09_EB,34.76_EB,34.42_EB,34.07_EB,33.73_EB,33.38_EB,32.96_EB,32.53_EB,32.10_EB,31.67_EB/)
QFLAME(4,12,21:30,2) = (/31.26_EB,30.84_EB,30.42_EB,30.00_EB,29.59_EB,29.21_EB,28.84_EB,28.46_EB,28.09_EB,27.71_EB/)
QFLAME(4,12,0:10,3) = (/0.00_EB,22.72_EB,30.28_EB,32.51_EB,34.75_EB,35.28_EB,35.82_EB,36.36_EB,36.90_EB,36.62_EB,36.34_EB/)
QFLAME(4,12,11:20,3) = (/36.07_EB,35.79_EB,35.50_EB,35.22_EB,34.93_EB,34.65_EB,34.19_EB,33.74_EB,33.28_EB,32.83_EB/)
QFLAME(4,12,21:30,3) = (/32.41_EB,31.99_EB,31.57_EB,31.15_EB,30.73_EB,30.69_EB,30.65_EB,30.60_EB,30.56_EB,30.52_EB/)
QFLAME(4,12,0:10,4) = (/0.00_EB,24.59_EB,33.51_EB,36.61_EB,39.70_EB,40.68_EB,41.66_EB,42.64_EB,43.62_EB,43.38_EB,43.14_EB/)
QFLAME(4,12,11:20,4) = (/42.89_EB,42.65_EB,43.03_EB,43.41_EB,43.78_EB,44.16_EB,43.68_EB,43.19_EB,42.71_EB,42.23_EB/)
QFLAME(4,12,21:30,4) = (/41.55_EB,40.87_EB,40.19_EB,39.51_EB,38.83_EB,39.21_EB,39.58_EB,39.95_EB,40.32_EB,40.69_EB/)
QFLAME(4,12,0:10,5) = (/0.00_EB,26.40_EB,37.10_EB,41.44_EB,45.77_EB,47.53_EB,49.29_EB,51.04_EB,52.80_EB,52.93_EB,53.05_EB/)
QFLAME(4,12,11:20,5) = (/53.17_EB,53.29_EB,53.08_EB,52.86_EB,52.65_EB,52.44_EB,52.54_EB,52.64_EB,52.74_EB,52.85_EB/)
QFLAME(4,12,21:30,5) = (/52.73_EB,52.62_EB,52.50_EB,52.38_EB,52.27_EB,52.54_EB,52.81_EB,53.08_EB,53.36_EB,53.63_EB/)
QFLAME(4,12,0:10,6) = (/0.00_EB,28.24_EB,40.62_EB,46.32_EB,52.03_EB,54.57_EB,57.10_EB,59.64_EB,62.18_EB,62.99_EB,63.80_EB/)
QFLAME(4,12,11:20,6) = (/64.61_EB,65.41_EB,65.71_EB,66.01_EB,66.31_EB,66.62_EB,66.62_EB,66.63_EB,66.64_EB,66.65_EB/)
QFLAME(4,12,21:30,6) = (/66.65_EB,66.65_EB,66.65_EB,66.64_EB,66.64_EB,66.41_EB,66.18_EB,65.95_EB,65.72_EB,65.49_EB/)
QFLAME(4,13,0:10,1) = (/0.00_EB,21.54_EB,29.28_EB,31.72_EB,34.15_EB,34.74_EB,35.34_EB,35.93_EB,36.52_EB,36.21_EB,35.91_EB/)
QFLAME(4,13,11:20,1) = (/35.60_EB,35.30_EB,34.96_EB,34.62_EB,34.28_EB,33.94_EB,33.52_EB,33.09_EB,32.67_EB,32.24_EB/)
QFLAME(4,13,21:30,1) = (/31.85_EB,31.45_EB,31.06_EB,30.66_EB,30.27_EB,29.89_EB,29.51_EB,29.13_EB,28.75_EB,28.37_EB/)
QFLAME(4,13,0:10,2) = (/0.00_EB,21.52_EB,29.27_EB,31.71_EB,34.15_EB,34.74_EB,35.33_EB,35.92_EB,36.52_EB,36.21_EB,35.90_EB/)
QFLAME(4,13,11:20,2) = (/35.60_EB,35.29_EB,34.96_EB,34.62_EB,34.28_EB,33.94_EB,33.52_EB,33.09_EB,32.67_EB,32.25_EB/)
QFLAME(4,13,21:30,2) = (/31.84_EB,31.42_EB,31.01_EB,30.60_EB,30.19_EB,29.82_EB,29.45_EB,29.08_EB,28.71_EB,28.33_EB/)
QFLAME(4,13,0:10,3) = (/0.00_EB,22.68_EB,30.32_EB,32.63_EB,34.94_EB,35.52_EB,36.10_EB,36.67_EB,37.25_EB,36.99_EB,36.73_EB/)
QFLAME(4,13,11:20,3) = (/36.47_EB,36.21_EB,35.92_EB,35.63_EB,35.34_EB,35.05_EB,34.60_EB,34.15_EB,33.70_EB,33.26_EB/)
QFLAME(4,13,21:30,3) = (/32.84_EB,32.43_EB,32.02_EB,31.61_EB,31.20_EB,31.12_EB,31.03_EB,30.95_EB,30.87_EB,30.78_EB/)
QFLAME(4,13,0:10,4) = (/0.00_EB,24.56_EB,33.47_EB,36.57_EB,39.67_EB,40.65_EB,41.63_EB,42.62_EB,43.60_EB,43.36_EB,43.13_EB/)
QFLAME(4,13,11:20,4) = (/42.89_EB,42.66_EB,43.01_EB,43.36_EB,43.71_EB,44.07_EB,43.59_EB,43.12_EB,42.65_EB,42.18_EB/)
QFLAME(4,13,21:30,4) = (/41.52_EB,40.86_EB,40.20_EB,39.54_EB,38.88_EB,39.16_EB,39.45_EB,39.74_EB,40.03_EB,40.32_EB/)
QFLAME(4,13,0:10,5) = (/0.00_EB,26.38_EB,37.06_EB,41.39_EB,45.73_EB,47.49_EB,49.25_EB,51.01_EB,52.77_EB,52.90_EB,53.02_EB/)
QFLAME(4,13,11:20,5) = (/53.14_EB,53.27_EB,53.02_EB,52.77_EB,52.53_EB,52.28_EB,52.38_EB,52.49_EB,52.59_EB,52.69_EB/)
QFLAME(4,13,21:30,5) = (/52.56_EB,52.42_EB,52.29_EB,52.15_EB,52.02_EB,52.34_EB,52.66_EB,52.98_EB,53.30_EB,53.62_EB/)
QFLAME(4,13,0:10,6) = (/0.00_EB,28.22_EB,40.52_EB,46.24_EB,51.96_EB,54.50_EB,57.03_EB,59.57_EB,62.10_EB,62.91_EB,63.72_EB/)
QFLAME(4,13,11:20,6) = (/64.52_EB,65.33_EB,65.56_EB,65.79_EB,66.02_EB,66.25_EB,66.23_EB,66.21_EB,66.20_EB,66.18_EB/)
QFLAME(4,13,21:30,6) = (/66.23_EB,66.27_EB,66.32_EB,66.36_EB,66.40_EB,66.11_EB,65.82_EB,65.53_EB,65.24_EB,64.95_EB/)
QFLAME(4,14,0:10,1) = (/0.00_EB,21.61_EB,29.45_EB,31.96_EB,34.46_EB,35.08_EB,35.71_EB,36.34_EB,36.96_EB,36.68_EB,36.39_EB/)
QFLAME(4,14,11:20,1) = (/36.11_EB,35.83_EB,35.49_EB,35.16_EB,34.83_EB,34.50_EB,34.08_EB,33.66_EB,33.24_EB,32.82_EB/)
QFLAME(4,14,21:30,1) = (/32.43_EB,32.03_EB,31.64_EB,31.24_EB,30.85_EB,30.47_EB,30.10_EB,29.73_EB,29.36_EB,28.98_EB/)
QFLAME(4,14,0:10,2) = (/0.00_EB,21.59_EB,29.44_EB,31.95_EB,34.45_EB,35.08_EB,35.71_EB,36.33_EB,36.96_EB,36.67_EB,36.39_EB/)
QFLAME(4,14,11:20,2) = (/36.11_EB,35.82_EB,35.49_EB,35.16_EB,34.83_EB,34.50_EB,34.08_EB,33.66_EB,33.24_EB,32.83_EB/)
QFLAME(4,14,21:30,2) = (/32.42_EB,32.01_EB,31.60_EB,31.20_EB,30.79_EB,30.42_EB,30.05_EB,29.69_EB,29.32_EB,28.95_EB/)
QFLAME(4,14,0:10,3) = (/0.00_EB,22.65_EB,30.36_EB,32.75_EB,35.14_EB,35.76_EB,36.37_EB,36.99_EB,37.61_EB,37.36_EB,37.12_EB/)
QFLAME(4,14,11:20,3) = (/36.87_EB,36.63_EB,36.34_EB,36.04_EB,35.75_EB,35.45_EB,35.01_EB,34.57_EB,34.13_EB,33.69_EB/)
QFLAME(4,14,21:30,3) = (/33.28_EB,32.88_EB,32.48_EB,32.07_EB,31.67_EB,31.55_EB,31.42_EB,31.30_EB,31.17_EB,31.05_EB/)
QFLAME(4,14,0:10,4) = (/0.00_EB,24.52_EB,33.42_EB,36.53_EB,39.63_EB,40.62_EB,41.60_EB,42.59_EB,43.57_EB,43.35_EB,43.12_EB/)
QFLAME(4,14,11:20,4) = (/42.89_EB,42.66_EB,42.99_EB,43.32_EB,43.65_EB,43.97_EB,43.51_EB,43.05_EB,42.59_EB,42.13_EB/)
QFLAME(4,14,21:30,4) = (/41.49_EB,40.84_EB,40.20_EB,39.56_EB,38.92_EB,39.12_EB,39.33_EB,39.53_EB,39.74_EB,39.95_EB/)
QFLAME(4,14,0:10,5) = (/0.00_EB,26.35_EB,37.01_EB,41.35_EB,45.68_EB,47.44_EB,49.21_EB,50.97_EB,52.74_EB,52.87_EB,52.99_EB/)
QFLAME(4,14,11:20,5) = (/53.12_EB,53.25_EB,52.96_EB,52.68_EB,52.40_EB,52.12_EB,52.22_EB,52.33_EB,52.44_EB,52.54_EB/)
QFLAME(4,14,21:30,5) = (/52.39_EB,52.23_EB,52.08_EB,51.92_EB,51.77_EB,52.14_EB,52.51_EB,52.88_EB,53.25_EB,53.62_EB/)
QFLAME(4,14,0:10,6) = (/0.00_EB,28.20_EB,40.41_EB,46.16_EB,51.90_EB,54.43_EB,56.96_EB,59.49_EB,62.02_EB,62.83_EB,63.63_EB/)
QFLAME(4,14,11:20,6) = (/64.44_EB,65.25_EB,65.41_EB,65.56_EB,65.72_EB,65.88_EB,65.83_EB,65.79_EB,65.75_EB,65.71_EB/)
QFLAME(4,14,21:30,6) = (/65.80_EB,65.89_EB,65.98_EB,66.07_EB,66.16_EB,65.81_EB,65.46_EB,65.11_EB,64.76_EB,64.40_EB/)
QFLAME(4,15,0:10,1) = (/0.00_EB,21.69_EB,29.63_EB,32.20_EB,34.76_EB,35.42_EB,36.08_EB,36.74_EB,37.41_EB,37.14_EB,36.88_EB/)
QFLAME(4,15,11:20,1) = (/36.62_EB,36.36_EB,36.03_EB,35.70_EB,35.38_EB,35.05_EB,34.64_EB,34.23_EB,33.81_EB,33.40_EB/)
QFLAME(4,15,21:30,1) = (/33.01_EB,32.61_EB,32.22_EB,31.82_EB,31.43_EB,31.06_EB,30.69_EB,30.33_EB,29.96_EB,29.59_EB/)
QFLAME(4,15,0:10,2) = (/0.00_EB,21.67_EB,29.62_EB,32.19_EB,34.76_EB,35.42_EB,36.08_EB,36.74_EB,37.40_EB,37.14_EB,36.88_EB/)
QFLAME(4,15,11:20,2) = (/36.61_EB,36.35_EB,36.03_EB,35.70_EB,35.38_EB,35.06_EB,34.64_EB,34.23_EB,33.82_EB,33.40_EB/)
QFLAME(4,15,21:30,2) = (/33.00_EB,32.60_EB,32.19_EB,31.79_EB,31.39_EB,31.03_EB,30.66_EB,30.30_EB,29.94_EB,29.57_EB/)
QFLAME(4,15,0:10,3) = (/0.00_EB,22.61_EB,30.41_EB,32.87_EB,35.34_EB,35.99_EB,36.65_EB,37.30_EB,37.96_EB,37.73_EB,37.51_EB/)
QFLAME(4,15,11:20,3) = (/37.28_EB,37.05_EB,36.75_EB,36.45_EB,36.15_EB,35.85_EB,35.42_EB,34.99_EB,34.55_EB,34.12_EB/)
QFLAME(4,15,21:30,3) = (/33.72_EB,33.33_EB,32.93_EB,32.54_EB,32.14_EB,31.97_EB,31.81_EB,31.64_EB,31.48_EB,31.31_EB/)
QFLAME(4,15,0:10,4) = (/0.00_EB,24.49_EB,33.38_EB,36.49_EB,39.60_EB,40.59_EB,41.58_EB,42.56_EB,43.55_EB,43.33_EB,43.11_EB/)
QFLAME(4,15,11:20,4) = (/42.89_EB,42.67_EB,42.97_EB,43.27_EB,43.58_EB,43.88_EB,43.43_EB,42.98_EB,42.53_EB,42.08_EB/)
QFLAME(4,15,21:30,4) = (/41.45_EB,40.83_EB,40.20_EB,39.58_EB,38.96_EB,39.08_EB,39.20_EB,39.33_EB,39.45_EB,39.57_EB/)
QFLAME(4,15,0:10,5) = (/0.00_EB,26.32_EB,36.97_EB,41.30_EB,45.63_EB,47.40_EB,49.17_EB,50.94_EB,52.71_EB,52.84_EB,52.97_EB/)
QFLAME(4,15,11:20,5) = (/53.09_EB,53.22_EB,52.91_EB,52.59_EB,52.28_EB,51.96_EB,52.07_EB,52.17_EB,52.28_EB,52.39_EB/)
QFLAME(4,15,21:30,5) = (/52.21_EB,52.04_EB,51.87_EB,51.69_EB,51.52_EB,51.94_EB,52.36_EB,52.77_EB,53.19_EB,53.61_EB/)
QFLAME(4,15,0:10,6) = (/0.00_EB,28.18_EB,40.30_EB,46.07_EB,51.84_EB,54.36_EB,56.89_EB,59.41_EB,61.93_EB,62.74_EB,63.55_EB/)
QFLAME(4,15,11:20,6) = (/64.36_EB,65.17_EB,65.25_EB,65.34_EB,65.42_EB,65.50_EB,65.44_EB,65.37_EB,65.31_EB,65.24_EB/)
QFLAME(4,15,21:30,6) = (/65.38_EB,65.52_EB,65.65_EB,65.79_EB,65.93_EB,65.51_EB,65.10_EB,64.69_EB,64.27_EB,63.86_EB/)
QFLAME(4,16,0:10,1) = (/0.00_EB,21.77_EB,29.80_EB,32.44_EB,35.07_EB,35.76_EB,36.46_EB,37.15_EB,37.85_EB,37.61_EB,37.37_EB/)
QFLAME(4,16,11:20,1) = (/37.12_EB,36.88_EB,36.56_EB,36.24_EB,35.93_EB,35.61_EB,35.20_EB,34.79_EB,34.39_EB,33.98_EB/)
QFLAME(4,16,21:30,1) = (/33.58_EB,33.19_EB,32.79_EB,32.40_EB,32.00_EB,31.64_EB,31.28_EB,30.92_EB,30.56_EB,30.20_EB/)
QFLAME(4,16,0:10,2) = (/0.00_EB,21.74_EB,29.79_EB,32.43_EB,35.06_EB,35.76_EB,36.45_EB,37.15_EB,37.84_EB,37.60_EB,37.36_EB/)
QFLAME(4,16,11:20,2) = (/37.12_EB,36.88_EB,36.56_EB,36.25_EB,35.93_EB,35.61_EB,35.21_EB,34.80_EB,34.39_EB,33.98_EB/)
QFLAME(4,16,21:30,2) = (/33.58_EB,33.18_EB,32.79_EB,32.39_EB,31.99_EB,31.63_EB,31.27_EB,30.91_EB,30.55_EB,30.19_EB/)
QFLAME(4,16,0:10,3) = (/0.00_EB,22.57_EB,30.45_EB,32.99_EB,35.54_EB,36.23_EB,36.92_EB,37.62_EB,38.31_EB,38.10_EB,37.89_EB/)
QFLAME(4,16,11:20,3) = (/37.68_EB,37.47_EB,37.17_EB,36.87_EB,36.56_EB,36.26_EB,35.83_EB,35.40_EB,34.97_EB,34.55_EB/)
QFLAME(4,16,21:30,3) = (/34.16_EB,33.77_EB,33.38_EB,33.00_EB,32.61_EB,32.40_EB,32.20_EB,31.99_EB,31.78_EB,31.58_EB/)
QFLAME(4,16,0:10,4) = (/0.00_EB,24.46_EB,33.33_EB,36.45_EB,39.57_EB,40.56_EB,41.55_EB,42.54_EB,43.53_EB,43.32_EB,43.10_EB/)
QFLAME(4,16,11:20,4) = (/42.89_EB,42.68_EB,42.95_EB,43.23_EB,43.51_EB,43.78_EB,43.35_EB,42.91_EB,42.47_EB,42.03_EB/)
QFLAME(4,16,21:30,4) = (/41.42_EB,40.82_EB,40.21_EB,39.60_EB,39.00_EB,39.04_EB,39.08_EB,39.12_EB,39.16_EB,39.20_EB/)
QFLAME(4,16,0:10,5) = (/0.00_EB,26.29_EB,36.92_EB,41.25_EB,45.59_EB,47.36_EB,49.13_EB,50.91_EB,52.68_EB,52.81_EB,52.94_EB/)
QFLAME(4,16,11:20,5) = (/53.07_EB,53.20_EB,52.85_EB,52.50_EB,52.15_EB,51.80_EB,51.91_EB,52.02_EB,52.13_EB,52.23_EB/)
QFLAME(4,16,21:30,5) = (/52.04_EB,51.85_EB,51.65_EB,51.46_EB,51.27_EB,51.74_EB,52.20_EB,52.67_EB,53.14_EB,53.61_EB/)
QFLAME(4,16,0:10,6) = (/0.00_EB,28.16_EB,40.20_EB,45.99_EB,51.78_EB,54.30_EB,56.81_EB,59.33_EB,61.85_EB,62.66_EB,63.47_EB/)
QFLAME(4,16,11:20,6) = (/64.28_EB,65.09_EB,65.10_EB,65.11_EB,65.12_EB,65.13_EB,65.04_EB,64.95_EB,64.86_EB,64.77_EB/)
QFLAME(4,16,21:30,6) = (/64.96_EB,65.14_EB,65.32_EB,65.50_EB,65.69_EB,65.21_EB,64.74_EB,64.27_EB,63.79_EB,63.32_EB/)
QFLAME(4,17,0:10,1) = (/0.00_EB,21.84_EB,29.98_EB,32.68_EB,35.37_EB,36.10_EB,36.83_EB,37.56_EB,38.29_EB,38.07_EB,37.85_EB/)
QFLAME(4,17,11:20,1) = (/37.63_EB,37.41_EB,37.10_EB,36.79_EB,36.47_EB,36.16_EB,35.76_EB,35.36_EB,34.96_EB,34.56_EB/)
QFLAME(4,17,21:30,1) = (/34.16_EB,33.77_EB,33.37_EB,32.98_EB,32.58_EB,32.23_EB,31.88_EB,31.52_EB,31.17_EB,30.81_EB/)
QFLAME(4,17,0:10,2) = (/0.00_EB,21.82_EB,29.97_EB,32.67_EB,35.36_EB,36.09_EB,36.82_EB,37.56_EB,38.29_EB,38.07_EB,37.85_EB/)
QFLAME(4,17,11:20,2) = (/37.63_EB,37.41_EB,37.10_EB,36.79_EB,36.48_EB,36.17_EB,35.77_EB,35.36_EB,34.96_EB,34.56_EB/)
QFLAME(4,17,21:30,2) = (/34.16_EB,33.77_EB,33.38_EB,32.99_EB,32.59_EB,32.24_EB,31.88_EB,31.53_EB,31.17_EB,30.81_EB/)
QFLAME(4,17,0:10,3) = (/0.00_EB,22.54_EB,30.49_EB,33.11_EB,35.73_EB,36.47_EB,37.20_EB,37.93_EB,38.67_EB,38.47_EB,38.28_EB/)
QFLAME(4,17,11:20,3) = (/38.09_EB,37.89_EB,37.59_EB,37.28_EB,36.97_EB,36.66_EB,36.24_EB,35.82_EB,35.40_EB,34.98_EB/)
QFLAME(4,17,21:30,3) = (/34.60_EB,34.22_EB,33.84_EB,33.46_EB,33.08_EB,32.83_EB,32.58_EB,32.34_EB,32.09_EB,31.84_EB/)
QFLAME(4,17,0:10,4) = (/0.00_EB,24.43_EB,33.29_EB,36.41_EB,39.53_EB,40.52_EB,41.52_EB,42.51_EB,43.51_EB,43.30_EB,43.10_EB/)
QFLAME(4,17,11:20,4) = (/42.89_EB,42.68_EB,42.94_EB,43.19_EB,43.44_EB,43.69_EB,43.26_EB,42.83_EB,42.41_EB,41.98_EB/)
QFLAME(4,17,21:30,4) = (/41.39_EB,40.80_EB,40.21_EB,39.63_EB,39.04_EB,38.99_EB,38.95_EB,38.91_EB,38.87_EB,38.83_EB/)
QFLAME(4,17,0:10,5) = (/0.00_EB,26.26_EB,36.87_EB,41.21_EB,45.54_EB,47.32_EB,49.09_EB,50.87_EB,52.65_EB,52.78_EB,52.91_EB/)
QFLAME(4,17,11:20,5) = (/53.04_EB,53.18_EB,52.79_EB,52.41_EB,52.03_EB,51.64_EB,51.75_EB,51.86_EB,51.97_EB,52.08_EB/)
QFLAME(4,17,21:30,5) = (/51.87_EB,51.66_EB,51.44_EB,51.23_EB,51.02_EB,51.53_EB,52.05_EB,52.57_EB,53.08_EB,53.60_EB/)
QFLAME(4,17,0:10,6) = (/0.00_EB,28.14_EB,40.09_EB,45.90_EB,51.72_EB,54.23_EB,56.74_EB,59.25_EB,61.77_EB,62.58_EB,63.39_EB/)
QFLAME(4,17,11:20,6) = (/64.20_EB,65.01_EB,64.95_EB,64.89_EB,64.83_EB,64.76_EB,64.65_EB,64.53_EB,64.42_EB,64.30_EB/)
QFLAME(4,17,21:30,6) = (/64.53_EB,64.76_EB,64.99_EB,65.22_EB,65.45_EB,64.91_EB,64.38_EB,63.84_EB,63.31_EB,62.78_EB/)
QFLAME(4,18,0:10,1) = (/0.00_EB,21.92_EB,30.15_EB,32.91_EB,35.68_EB,36.44_EB,37.21_EB,37.97_EB,38.74_EB,38.54_EB,38.34_EB/)
QFLAME(4,18,11:20,1) = (/38.14_EB,37.94_EB,37.63_EB,37.33_EB,37.02_EB,36.71_EB,36.32_EB,35.92_EB,35.53_EB,35.13_EB/)
QFLAME(4,18,21:30,1) = (/34.74_EB,34.35_EB,33.95_EB,33.56_EB,33.16_EB,32.82_EB,32.47_EB,32.12_EB,31.77_EB,31.42_EB/)
QFLAME(4,18,0:10,2) = (/0.00_EB,21.89_EB,30.15_EB,32.91_EB,35.67_EB,36.43_EB,37.20_EB,37.96_EB,38.73_EB,38.53_EB,38.33_EB/)
QFLAME(4,18,11:20,2) = (/38.14_EB,37.94_EB,37.64_EB,37.33_EB,37.03_EB,36.73_EB,36.33_EB,35.93_EB,35.53_EB,35.13_EB/)
QFLAME(4,18,21:30,2) = (/34.74_EB,34.36_EB,33.97_EB,33.58_EB,33.19_EB,32.84_EB,32.49_EB,32.14_EB,31.79_EB,31.43_EB/)
QFLAME(4,18,0:10,3) = (/0.00_EB,22.50_EB,30.53_EB,33.23_EB,35.93_EB,36.70_EB,37.48_EB,38.25_EB,39.02_EB,38.85_EB,38.67_EB/)
QFLAME(4,18,11:20,3) = (/38.49_EB,38.32_EB,38.00_EB,37.69_EB,37.38_EB,37.06_EB,36.65_EB,36.23_EB,35.82_EB,35.41_EB/)
QFLAME(4,18,21:30,3) = (/35.03_EB,34.66_EB,34.29_EB,33.92_EB,33.55_EB,33.26_EB,32.97_EB,32.68_EB,32.39_EB,32.10_EB/)
QFLAME(4,18,0:10,4) = (/0.00_EB,24.39_EB,33.24_EB,36.37_EB,39.50_EB,40.49_EB,41.49_EB,42.49_EB,43.48_EB,43.29_EB,43.09_EB/)
QFLAME(4,18,11:20,4) = (/42.89_EB,42.69_EB,42.92_EB,43.14_EB,43.37_EB,43.60_EB,43.18_EB,42.76_EB,42.35_EB,41.93_EB/)
QFLAME(4,18,21:30,4) = (/41.36_EB,40.79_EB,40.22_EB,39.65_EB,39.08_EB,38.95_EB,38.83_EB,38.70_EB,38.58_EB,38.46_EB/)
QFLAME(4,18,0:10,5) = (/0.00_EB,26.23_EB,36.83_EB,41.16_EB,45.50_EB,47.28_EB,49.06_EB,50.84_EB,52.61_EB,52.75_EB,52.88_EB/)
QFLAME(4,18,11:20,5) = (/53.02_EB,53.15_EB,52.74_EB,52.32_EB,51.90_EB,51.48_EB,51.59_EB,51.71_EB,51.82_EB,51.93_EB/)
QFLAME(4,18,21:30,5) = (/51.70_EB,51.46_EB,51.23_EB,51.00_EB,50.77_EB,51.33_EB,51.90_EB,52.46_EB,53.03_EB,53.60_EB/)
QFLAME(4,18,0:10,6) = (/0.00_EB,28.12_EB,39.99_EB,45.82_EB,51.66_EB,54.16_EB,56.67_EB,59.18_EB,61.68_EB,62.49_EB,63.30_EB/)
QFLAME(4,18,11:20,6) = (/64.11_EB,64.93_EB,64.79_EB,64.66_EB,64.53_EB,64.39_EB,64.25_EB,64.11_EB,63.97_EB,63.83_EB/)
QFLAME(4,18,21:30,6) = (/64.11_EB,64.38_EB,64.66_EB,64.93_EB,65.21_EB,64.61_EB,64.02_EB,63.42_EB,62.83_EB,62.23_EB/)
QFLAME(4,19,0:10,1) = (/0.00_EB,21.99_EB,30.33_EB,33.15_EB,35.98_EB,36.78_EB,37.58_EB,38.38_EB,39.18_EB,39.00_EB,38.82_EB/)
QFLAME(4,19,11:20,1) = (/38.65_EB,38.47_EB,38.17_EB,37.87_EB,37.57_EB,37.27_EB,36.88_EB,36.49_EB,36.10_EB,35.71_EB/)
QFLAME(4,19,21:30,1) = (/35.32_EB,34.92_EB,34.53_EB,34.14_EB,33.74_EB,33.40_EB,33.06_EB,32.72_EB,32.38_EB,32.03_EB/)
QFLAME(4,19,0:10,2) = (/0.00_EB,21.97_EB,30.32_EB,33.15_EB,35.97_EB,36.77_EB,37.57_EB,38.37_EB,39.17_EB,39.00_EB,38.82_EB/)
QFLAME(4,19,11:20,2) = (/38.65_EB,38.47_EB,38.17_EB,37.88_EB,37.58_EB,37.29_EB,36.89_EB,36.50_EB,36.10_EB,35.71_EB/)
QFLAME(4,19,21:30,2) = (/35.33_EB,34.94_EB,34.56_EB,34.18_EB,33.80_EB,33.45_EB,33.10_EB,32.75_EB,32.40_EB,32.05_EB/)
QFLAME(4,19,0:10,3) = (/0.00_EB,22.46_EB,30.58_EB,33.35_EB,36.13_EB,36.94_EB,37.75_EB,38.56_EB,39.38_EB,39.22_EB,39.06_EB/)
QFLAME(4,19,11:20,3) = (/38.90_EB,38.74_EB,38.42_EB,38.10_EB,37.78_EB,37.47_EB,37.06_EB,36.65_EB,36.24_EB,35.84_EB/)
QFLAME(4,19,21:30,3) = (/35.47_EB,35.11_EB,34.75_EB,34.38_EB,34.02_EB,33.69_EB,33.36_EB,33.03_EB,32.70_EB,32.37_EB/)
QFLAME(4,19,0:10,4) = (/0.00_EB,24.36_EB,33.19_EB,36.33_EB,39.46_EB,40.46_EB,41.46_EB,42.46_EB,43.46_EB,43.27_EB,43.08_EB/)
QFLAME(4,19,11:20,4) = (/42.89_EB,42.70_EB,42.90_EB,43.10_EB,43.30_EB,43.50_EB,43.10_EB,42.69_EB,42.29_EB,41.88_EB/)
QFLAME(4,19,21:30,4) = (/41.33_EB,40.78_EB,40.22_EB,39.67_EB,39.12_EB,38.91_EB,38.70_EB,38.50_EB,38.29_EB,38.08_EB/)
QFLAME(4,19,0:10,5) = (/0.00_EB,26.21_EB,36.78_EB,41.12_EB,45.45_EB,47.24_EB,49.02_EB,50.80_EB,52.58_EB,52.72_EB,52.86_EB/)
QFLAME(4,19,11:20,5) = (/52.99_EB,53.13_EB,52.68_EB,52.23_EB,51.78_EB,51.32_EB,51.44_EB,51.55_EB,51.66_EB,51.77_EB/)
QFLAME(4,19,21:30,5) = (/51.52_EB,51.27_EB,51.02_EB,50.77_EB,50.52_EB,51.13_EB,51.75_EB,52.36_EB,52.98_EB,53.59_EB/)
QFLAME(4,19,0:10,6) = (/0.00_EB,28.10_EB,39.88_EB,45.74_EB,51.59_EB,54.10_EB,56.60_EB,59.10_EB,61.60_EB,62.41_EB,63.22_EB/)
QFLAME(4,19,11:20,6) = (/64.03_EB,64.84_EB,64.64_EB,64.43_EB,64.23_EB,64.02_EB,63.86_EB,63.69_EB,63.53_EB,63.36_EB/)
QFLAME(4,19,21:30,6) = (/63.68_EB,64.00_EB,64.33_EB,64.65_EB,64.97_EB,64.31_EB,63.66_EB,63.00_EB,62.35_EB,61.69_EB/)
QFLAME(4,20,0:10,1) = (/0.00_EB,22.07_EB,30.50_EB,33.39_EB,36.29_EB,37.12_EB,37.96_EB,38.79_EB,39.62_EB,39.47_EB,39.31_EB/)
QFLAME(4,20,11:20,1) = (/39.15_EB,39.00_EB,38.70_EB,38.41_EB,38.12_EB,37.82_EB,37.44_EB,37.06_EB,36.67_EB,36.29_EB/)
QFLAME(4,20,21:30,1) = (/35.90_EB,35.50_EB,35.11_EB,34.72_EB,34.32_EB,33.99_EB,33.65_EB,33.32_EB,32.98_EB,32.64_EB/)
QFLAME(4,20,0:10,2) = (/0.00_EB,22.04_EB,30.50_EB,33.39_EB,36.28_EB,37.11_EB,37.94_EB,38.78_EB,39.61_EB,39.46_EB,39.31_EB/)
QFLAME(4,20,11:20,2) = (/39.15_EB,39.00_EB,38.71_EB,38.42_EB,38.13_EB,37.85_EB,37.46_EB,37.07_EB,36.68_EB,36.29_EB/)
QFLAME(4,20,21:30,2) = (/35.91_EB,35.53_EB,35.15_EB,34.77_EB,34.40_EB,34.05_EB,33.71_EB,33.36_EB,33.02_EB,32.67_EB/)
QFLAME(4,20,0:10,3) = (/0.00_EB,22.43_EB,30.62_EB,33.47_EB,36.33_EB,37.18_EB,38.03_EB,38.88_EB,39.73_EB,39.59_EB,39.44_EB/)
QFLAME(4,20,11:20,3) = (/39.30_EB,39.16_EB,38.84_EB,38.51_EB,38.19_EB,37.87_EB,37.47_EB,37.07_EB,36.67_EB,36.27_EB/)
QFLAME(4,20,21:30,3) = (/35.91_EB,35.55_EB,35.20_EB,34.84_EB,34.49_EB,34.12_EB,33.75_EB,33.38_EB,33.00_EB,32.63_EB/)
QFLAME(4,20,0:10,4) = (/0.00_EB,24.33_EB,33.15_EB,36.29_EB,39.43_EB,40.43_EB,41.43_EB,42.43_EB,43.44_EB,43.25_EB,43.07_EB/)
QFLAME(4,20,11:20,4) = (/42.89_EB,42.70_EB,42.88_EB,43.06_EB,43.23_EB,43.41_EB,43.01_EB,42.62_EB,42.23_EB,41.83_EB/)
QFLAME(4,20,21:30,4) = (/41.30_EB,40.76_EB,40.23_EB,39.69_EB,39.16_EB,38.87_EB,38.58_EB,38.29_EB,38.00_EB,37.71_EB/)
QFLAME(4,20,0:10,5) = (/0.00_EB,26.18_EB,36.74_EB,41.07_EB,45.41_EB,47.19_EB,48.98_EB,50.77_EB,52.55_EB,52.69_EB,52.83_EB/)
QFLAME(4,20,11:20,5) = (/52.97_EB,53.11_EB,52.62_EB,52.14_EB,51.65_EB,51.17_EB,51.28_EB,51.39_EB,51.51_EB,51.62_EB/)
QFLAME(4,20,21:30,5) = (/51.35_EB,51.08_EB,50.81_EB,50.54_EB,50.27_EB,50.93_EB,51.59_EB,52.26_EB,52.92_EB,53.59_EB/)
QFLAME(4,20,0:10,6) = (/0.00_EB,28.08_EB,39.77_EB,45.65_EB,51.53_EB,54.03_EB,56.52_EB,59.02_EB,61.52_EB,62.33_EB,63.14_EB/)
QFLAME(4,20,11:20,6) = (/63.95_EB,64.76_EB,64.49_EB,64.21_EB,63.93_EB,63.65_EB,63.46_EB,63.27_EB,63.08_EB,62.89_EB/)
QFLAME(4,20,21:30,6) = (/63.26_EB,63.63_EB,64.00_EB,64.36_EB,64.73_EB,64.01_EB,63.30_EB,62.58_EB,61.86_EB,61.15_EB/)
QFLAME(5,0,0:10,1) = (/0.00_EB,19.08_EB,24.58_EB,25.64_EB,26.71_EB,26.70_EB,26.70_EB,26.69_EB,26.69_EB,26.19_EB,25.69_EB/)
QFLAME(5,0,11:20,1) = (/25.18_EB,24.68_EB,24.18_EB,23.67_EB,23.17_EB,22.66_EB,22.23_EB,21.79_EB,21.36_EB,20.93_EB/)
QFLAME(5,0,21:30,1) = (/20.55_EB,20.17_EB,19.79_EB,19.42_EB,19.04_EB,18.71_EB,18.37_EB,18.04_EB,17.71_EB,17.38_EB/)
QFLAME(5,0,0:10,2) = (/0.00_EB,20.98_EB,27.24_EB,28.64_EB,30.04_EB,30.16_EB,30.29_EB,30.42_EB,30.54_EB,30.08_EB,29.61_EB/)
QFLAME(5,0,11:20,2) = (/29.15_EB,28.68_EB,28.26_EB,27.85_EB,27.43_EB,27.01_EB,26.78_EB,26.54_EB,26.31_EB,26.08_EB/)
QFLAME(5,0,21:30,2) = (/25.60_EB,25.13_EB,24.66_EB,24.18_EB,23.71_EB,23.36_EB,23.01_EB,22.66_EB,22.31_EB,21.95_EB/)
QFLAME(5,0,0:10,3) = (/0.00_EB,22.92_EB,30.54_EB,32.82_EB,35.11_EB,35.71_EB,36.32_EB,36.93_EB,37.53_EB,37.31_EB,37.08_EB/)
QFLAME(5,0,11:20,3) = (/36.86_EB,36.63_EB,36.28_EB,35.92_EB,35.57_EB,35.22_EB,35.07_EB,34.93_EB,34.78_EB,34.64_EB/)
QFLAME(5,0,21:30,3) = (/34.33_EB,34.02_EB,33.70_EB,33.39_EB,33.08_EB,32.93_EB,32.77_EB,32.62_EB,32.47_EB,32.31_EB/)
QFLAME(5,0,0:10,4) = (/0.00_EB,24.71_EB,33.97_EB,37.29_EB,40.61_EB,41.83_EB,43.05_EB,44.26_EB,45.48_EB,45.60_EB,45.72_EB/)
QFLAME(5,0,11:20,4) = (/45.84_EB,45.96_EB,45.81_EB,45.66_EB,45.51_EB,45.36_EB,45.09_EB,44.81_EB,44.53_EB,44.26_EB/)
QFLAME(5,0,21:30,4) = (/44.05_EB,43.84_EB,43.63_EB,43.42_EB,43.21_EB,42.90_EB,42.60_EB,42.30_EB,42.00_EB,41.70_EB/)
QFLAME(5,0,0:10,5) = (/0.00_EB,26.45_EB,37.64_EB,42.40_EB,47.16_EB,49.02_EB,50.87_EB,52.72_EB,54.58_EB,55.13_EB,55.69_EB/)
QFLAME(5,0,11:20,5) = (/56.24_EB,56.80_EB,56.85_EB,56.91_EB,56.96_EB,57.02_EB,56.71_EB,56.39_EB,56.08_EB,55.77_EB/)
QFLAME(5,0,21:30,5) = (/55.74_EB,55.72_EB,55.70_EB,55.68_EB,55.65_EB,55.31_EB,54.96_EB,54.62_EB,54.27_EB,53.93_EB/)
QFLAME(5,0,0:10,6) = (/0.00_EB,28.27_EB,41.21_EB,47.40_EB,53.60_EB,56.00_EB,58.40_EB,60.81_EB,63.21_EB,64.53_EB,65.86_EB/)
QFLAME(5,0,11:20,6) = (/67.18_EB,68.51_EB,68.93_EB,69.34_EB,69.76_EB,70.17_EB,69.64_EB,69.11_EB,68.58_EB,68.05_EB/)
QFLAME(5,0,21:30,6) = (/68.04_EB,68.03_EB,68.02_EB,68.00_EB,67.99_EB,67.95_EB,67.90_EB,67.85_EB,67.81_EB,67.76_EB/)
QFLAME(5,1,0:10,1) = (/0.00_EB,19.26_EB,25.12_EB,26.38_EB,27.64_EB,27.74_EB,27.85_EB,27.95_EB,28.05_EB,27.58_EB,27.12_EB/)
QFLAME(5,1,11:20,1) = (/26.65_EB,26.19_EB,25.71_EB,25.24_EB,24.76_EB,24.29_EB,23.86_EB,23.44_EB,23.02_EB,22.59_EB/)
QFLAME(5,1,21:30,1) = (/22.21_EB,21.83_EB,21.45_EB,21.06_EB,20.68_EB,20.35_EB,20.01_EB,19.68_EB,19.35_EB,19.01_EB/)
QFLAME(5,1,0:10,2) = (/0.00_EB,20.93_EB,27.23_EB,28.64_EB,30.06_EB,30.20_EB,30.34_EB,30.49_EB,30.63_EB,30.16_EB,29.68_EB/)
QFLAME(5,1,11:20,2) = (/29.21_EB,28.73_EB,28.32_EB,27.90_EB,27.49_EB,27.07_EB,26.73_EB,26.38_EB,26.04_EB,25.69_EB/)
QFLAME(5,1,21:30,2) = (/25.31_EB,24.93_EB,24.55_EB,24.17_EB,23.79_EB,23.44_EB,23.10_EB,22.75_EB,22.41_EB,22.06_EB/)
QFLAME(5,1,0:10,3) = (/0.00_EB,22.88_EB,30.53_EB,32.84_EB,35.14_EB,35.75_EB,36.36_EB,36.97_EB,37.58_EB,37.32_EB,37.06_EB/)
QFLAME(5,1,11:20,3) = (/36.80_EB,36.53_EB,36.26_EB,35.99_EB,35.72_EB,35.44_EB,35.24_EB,35.04_EB,34.83_EB,34.63_EB/)
QFLAME(5,1,21:30,3) = (/34.37_EB,34.10_EB,33.84_EB,33.58_EB,33.32_EB,32.83_EB,32.34_EB,31.86_EB,31.37_EB,30.88_EB/)
QFLAME(5,1,0:10,4) = (/0.00_EB,24.66_EB,33.95_EB,37.29_EB,40.63_EB,41.85_EB,43.07_EB,44.29_EB,45.51_EB,45.63_EB,45.76_EB/)
QFLAME(5,1,11:20,4) = (/45.88_EB,46.00_EB,45.87_EB,45.74_EB,45.61_EB,45.48_EB,45.52_EB,45.56_EB,45.59_EB,45.63_EB/)
QFLAME(5,1,21:30,4) = (/45.15_EB,44.67_EB,44.20_EB,43.72_EB,43.24_EB,42.84_EB,42.45_EB,42.05_EB,41.65_EB,41.26_EB/)
QFLAME(5,1,0:10,5) = (/0.00_EB,26.41_EB,37.61_EB,42.39_EB,47.17_EB,49.04_EB,50.91_EB,52.78_EB,54.65_EB,55.21_EB,55.78_EB/)
QFLAME(5,1,11:20,5) = (/56.34_EB,56.90_EB,56.92_EB,56.93_EB,56.94_EB,56.95_EB,56.56_EB,56.16_EB,55.77_EB,55.38_EB/)
QFLAME(5,1,21:30,5) = (/55.44_EB,55.50_EB,55.56_EB,55.62_EB,55.68_EB,55.38_EB,55.07_EB,54.77_EB,54.47_EB,54.17_EB/)
QFLAME(5,1,0:10,6) = (/0.00_EB,28.28_EB,41.20_EB,47.36_EB,53.53_EB,55.95_EB,58.38_EB,60.80_EB,63.22_EB,64.54_EB,65.87_EB/)
QFLAME(5,1,11:20,6) = (/67.19_EB,68.51_EB,68.90_EB,69.28_EB,69.67_EB,70.06_EB,69.52_EB,68.99_EB,68.46_EB,67.92_EB/)
QFLAME(5,1,21:30,6) = (/67.90_EB,67.87_EB,67.85_EB,67.83_EB,67.80_EB,67.63_EB,67.47_EB,67.30_EB,67.14_EB,66.97_EB/)
QFLAME(5,2,0:10,1) = (/0.00_EB,19.53_EB,25.58_EB,27.03_EB,28.47_EB,28.66_EB,28.86_EB,29.05_EB,29.24_EB,28.82_EB,28.40_EB/)
QFLAME(5,2,11:20,1) = (/27.98_EB,27.57_EB,27.10_EB,26.64_EB,26.17_EB,25.71_EB,25.30_EB,24.89_EB,24.47_EB,24.06_EB/)
QFLAME(5,2,21:30,1) = (/23.68_EB,23.31_EB,22.93_EB,22.55_EB,22.17_EB,21.84_EB,21.50_EB,21.17_EB,20.83_EB,20.50_EB/)
QFLAME(5,2,0:10,2) = (/0.00_EB,20.92_EB,27.22_EB,28.67_EB,30.13_EB,30.27_EB,30.41_EB,30.55_EB,30.69_EB,30.26_EB,29.84_EB/)
QFLAME(5,2,11:20,2) = (/29.41_EB,28.98_EB,28.54_EB,28.10_EB,27.66_EB,27.21_EB,26.88_EB,26.54_EB,26.20_EB,25.86_EB/)
QFLAME(5,2,21:30,2) = (/25.47_EB,25.09_EB,24.70_EB,24.31_EB,23.92_EB,23.58_EB,23.25_EB,22.91_EB,22.57_EB,22.23_EB/)
QFLAME(5,2,0:10,3) = (/0.00_EB,22.77_EB,30.52_EB,32.84_EB,35.16_EB,35.78_EB,36.41_EB,37.03_EB,37.65_EB,37.43_EB,37.21_EB/)
QFLAME(5,2,11:20,3) = (/36.98_EB,36.76_EB,36.36_EB,35.97_EB,35.57_EB,35.18_EB,34.96_EB,34.75_EB,34.53_EB,34.32_EB/)
QFLAME(5,2,21:30,3) = (/34.05_EB,33.79_EB,33.53_EB,33.26_EB,33.00_EB,32.70_EB,32.40_EB,32.10_EB,31.80_EB,31.50_EB/)
QFLAME(5,2,0:10,4) = (/0.00_EB,24.63_EB,33.92_EB,37.27_EB,40.62_EB,41.85_EB,43.08_EB,44.30_EB,45.53_EB,45.66_EB,45.80_EB/)
QFLAME(5,2,11:20,4) = (/45.93_EB,46.07_EB,45.97_EB,45.88_EB,45.78_EB,45.69_EB,45.29_EB,44.89_EB,44.49_EB,44.09_EB/)
QFLAME(5,2,21:30,4) = (/43.93_EB,43.78_EB,43.62_EB,43.47_EB,43.31_EB,42.83_EB,42.36_EB,41.88_EB,41.41_EB,40.93_EB/)
QFLAME(5,2,0:10,5) = (/0.00_EB,26.37_EB,37.58_EB,42.37_EB,47.17_EB,49.05_EB,50.93_EB,52.81_EB,54.69_EB,55.23_EB,55.77_EB/)
QFLAME(5,2,11:20,5) = (/56.31_EB,56.84_EB,56.88_EB,56.92_EB,56.95_EB,56.99_EB,56.60_EB,56.21_EB,55.82_EB,55.43_EB/)
QFLAME(5,2,21:30,5) = (/55.42_EB,55.42_EB,55.41_EB,55.41_EB,55.40_EB,55.12_EB,54.84_EB,54.56_EB,54.27_EB,53.99_EB/)
QFLAME(5,2,0:10,6) = (/0.00_EB,28.25_EB,41.11_EB,47.31_EB,53.51_EB,55.93_EB,58.35_EB,60.77_EB,63.18_EB,64.50_EB,65.82_EB/)
QFLAME(5,2,11:20,6) = (/67.14_EB,68.45_EB,68.86_EB,69.27_EB,69.67_EB,70.08_EB,69.57_EB,69.05_EB,68.54_EB,68.03_EB/)
QFLAME(5,2,21:30,6) = (/68.01_EB,67.99_EB,67.98_EB,67.96_EB,67.94_EB,67.80_EB,67.66_EB,67.52_EB,67.38_EB,67.24_EB/)
QFLAME(5,3,0:10,1) = (/0.00_EB,19.70_EB,25.96_EB,27.56_EB,29.15_EB,29.41_EB,29.67_EB,29.94_EB,30.20_EB,29.82_EB,29.44_EB/)
QFLAME(5,3,11:20,1) = (/29.07_EB,28.69_EB,28.24_EB,27.78_EB,27.33_EB,26.87_EB,26.47_EB,26.07_EB,25.68_EB,25.28_EB/)
QFLAME(5,3,21:30,1) = (/24.90_EB,24.52_EB,24.15_EB,23.77_EB,23.40_EB,23.06_EB,22.73_EB,22.40_EB,22.07_EB,21.74_EB/)
QFLAME(5,3,0:10,2) = (/0.00_EB,20.89_EB,27.22_EB,28.77_EB,30.31_EB,30.55_EB,30.78_EB,31.02_EB,31.26_EB,30.86_EB,30.47_EB/)
QFLAME(5,3,11:20,2) = (/30.08_EB,29.69_EB,29.24_EB,28.79_EB,28.35_EB,27.90_EB,27.55_EB,27.21_EB,26.86_EB,26.52_EB/)
QFLAME(5,3,21:30,2) = (/26.14_EB,25.75_EB,25.37_EB,24.98_EB,24.60_EB,24.26_EB,23.92_EB,23.58_EB,23.24_EB,22.90_EB/)
QFLAME(5,3,0:10,3) = (/0.00_EB,22.77_EB,30.49_EB,32.83_EB,35.17_EB,35.79_EB,36.41_EB,37.04_EB,37.66_EB,37.43_EB,37.21_EB/)
QFLAME(5,3,11:20,3) = (/36.98_EB,36.76_EB,36.36_EB,35.97_EB,35.57_EB,35.18_EB,34.96_EB,34.73_EB,34.51_EB,34.29_EB/)
QFLAME(5,3,21:30,3) = (/34.05_EB,33.82_EB,33.58_EB,33.35_EB,33.11_EB,32.78_EB,32.46_EB,32.13_EB,31.80_EB,31.48_EB/)
QFLAME(5,3,0:10,4) = (/0.00_EB,24.59_EB,33.89_EB,37.26_EB,40.62_EB,41.85_EB,43.08_EB,44.31_EB,45.54_EB,45.68_EB,45.81_EB/)
QFLAME(5,3,11:20,4) = (/45.95_EB,46.09_EB,45.98_EB,45.87_EB,45.77_EB,45.66_EB,45.39_EB,45.13_EB,44.86_EB,44.60_EB/)
QFLAME(5,3,21:30,4) = (/44.36_EB,44.12_EB,43.89_EB,43.65_EB,43.42_EB,42.94_EB,42.47_EB,42.00_EB,41.53_EB,41.06_EB/)
QFLAME(5,3,0:10,5) = (/0.00_EB,26.34_EB,37.54_EB,42.34_EB,47.14_EB,49.03_EB,50.93_EB,52.82_EB,54.72_EB,55.23_EB,55.75_EB/)
QFLAME(5,3,11:20,5) = (/56.26_EB,56.78_EB,56.85_EB,56.93_EB,57.01_EB,57.08_EB,56.65_EB,56.22_EB,55.79_EB,55.36_EB/)
QFLAME(5,3,21:30,5) = (/55.38_EB,55.40_EB,55.42_EB,55.44_EB,55.46_EB,55.18_EB,54.89_EB,54.61_EB,54.33_EB,54.05_EB/)
QFLAME(5,3,0:10,6) = (/0.00_EB,28.20_EB,41.06_EB,47.26_EB,53.46_EB,55.88_EB,58.30_EB,60.72_EB,63.14_EB,64.45_EB,65.76_EB/)
QFLAME(5,3,11:20,6) = (/67.07_EB,68.37_EB,68.78_EB,69.19_EB,69.60_EB,70.01_EB,69.47_EB,68.94_EB,68.41_EB,67.87_EB/)
QFLAME(5,3,21:30,6) = (/67.87_EB,67.86_EB,67.86_EB,67.86_EB,67.85_EB,67.71_EB,67.57_EB,67.43_EB,67.29_EB,67.15_EB/)
QFLAME(5,4,0:10,1) = (/0.00_EB,19.87_EB,26.35_EB,28.09_EB,29.82_EB,30.15_EB,30.49_EB,30.82_EB,31.16_EB,30.82_EB,30.49_EB/)
QFLAME(5,4,11:20,1) = (/30.15_EB,29.81_EB,29.37_EB,28.93_EB,28.48_EB,28.04_EB,27.65_EB,27.26_EB,26.88_EB,26.49_EB/)
QFLAME(5,4,21:30,1) = (/26.12_EB,25.74_EB,25.37_EB,25.00_EB,24.62_EB,24.29_EB,23.96_EB,23.63_EB,23.30_EB,22.97_EB/)
QFLAME(5,4,0:10,2) = (/0.00_EB,20.86_EB,27.23_EB,28.86_EB,30.50_EB,30.83_EB,31.16_EB,31.49_EB,31.82_EB,31.47_EB,31.11_EB/)
QFLAME(5,4,11:20,2) = (/30.75_EB,30.40_EB,29.94_EB,29.49_EB,29.03_EB,28.58_EB,28.23_EB,27.88_EB,27.53_EB,27.18_EB/)
QFLAME(5,4,21:30,2) = (/26.80_EB,26.42_EB,26.04_EB,25.66_EB,25.27_EB,24.93_EB,24.59_EB,24.25_EB,23.90_EB,23.56_EB/)
QFLAME(5,4,0:10,3) = (/0.00_EB,22.77_EB,30.46_EB,32.82_EB,35.18_EB,35.80_EB,36.42_EB,37.04_EB,37.66_EB,37.44_EB,37.21_EB/)
QFLAME(5,4,11:20,3) = (/36.98_EB,36.76_EB,36.36_EB,35.97_EB,35.58_EB,35.18_EB,34.95_EB,34.72_EB,34.49_EB,34.25_EB/)
QFLAME(5,4,21:30,3) = (/34.05_EB,33.84_EB,33.64_EB,33.43_EB,33.23_EB,32.87_EB,32.52_EB,32.16_EB,31.80_EB,31.45_EB/)
QFLAME(5,4,0:10,4) = (/0.00_EB,24.55_EB,33.86_EB,37.24_EB,40.61_EB,41.85_EB,43.08_EB,44.31_EB,45.55_EB,45.69_EB,45.83_EB/)
QFLAME(5,4,11:20,4) = (/45.97_EB,46.11_EB,45.99_EB,45.87_EB,45.75_EB,45.63_EB,45.50_EB,45.37_EB,45.23_EB,45.10_EB/)
QFLAME(5,4,21:30,4) = (/44.79_EB,44.47_EB,44.15_EB,43.84_EB,43.52_EB,43.05_EB,42.58_EB,42.12_EB,41.65_EB,41.18_EB/)
QFLAME(5,4,0:10,5) = (/0.00_EB,26.32_EB,37.51_EB,42.31_EB,47.11_EB,49.02_EB,50.93_EB,52.84_EB,54.75_EB,55.24_EB,55.73_EB/)
QFLAME(5,4,11:20,5) = (/56.22_EB,56.71_EB,56.83_EB,56.94_EB,57.06_EB,57.18_EB,56.71_EB,56.24_EB,55.77_EB,55.30_EB/)
QFLAME(5,4,21:30,5) = (/55.34_EB,55.38_EB,55.43_EB,55.47_EB,55.51_EB,55.23_EB,54.95_EB,54.67_EB,54.38_EB,54.10_EB/)
QFLAME(5,4,0:10,6) = (/0.00_EB,28.16_EB,41.00_EB,47.21_EB,53.42_EB,55.84_EB,58.26_EB,60.68_EB,63.10_EB,64.40_EB,65.70_EB/)
QFLAME(5,4,11:20,6) = (/67.00_EB,68.30_EB,68.70_EB,69.11_EB,69.52_EB,69.93_EB,69.38_EB,68.82_EB,68.27_EB,67.72_EB/)
QFLAME(5,4,21:30,6) = (/67.73_EB,67.74_EB,67.75_EB,67.75_EB,67.76_EB,67.62_EB,67.48_EB,67.34_EB,67.20_EB,67.06_EB/)
QFLAME(5,5,0:10,1) = (/0.00_EB,20.03_EB,26.74_EB,28.61_EB,30.49_EB,30.90_EB,31.31_EB,31.71_EB,32.12_EB,31.83_EB,31.53_EB/)
QFLAME(5,5,11:20,1) = (/31.23_EB,30.94_EB,30.50_EB,30.07_EB,29.64_EB,29.21_EB,28.83_EB,28.45_EB,28.08_EB,27.70_EB/)
QFLAME(5,5,21:30,1) = (/27.33_EB,26.96_EB,26.59_EB,26.22_EB,25.85_EB,25.52_EB,25.19_EB,24.87_EB,24.54_EB,24.21_EB/)
QFLAME(5,5,0:10,2) = (/0.00_EB,20.83_EB,27.24_EB,28.96_EB,30.68_EB,31.11_EB,31.54_EB,31.96_EB,32.39_EB,32.07_EB,31.75_EB/)
QFLAME(5,5,11:20,2) = (/31.42_EB,31.10_EB,30.64_EB,30.18_EB,29.72_EB,29.26_EB,28.91_EB,28.55_EB,28.19_EB,27.84_EB/)
QFLAME(5,5,21:30,2) = (/27.46_EB,27.08_EB,26.71_EB,26.33_EB,25.95_EB,25.61_EB,25.26_EB,24.92_EB,24.57_EB,24.23_EB/)
QFLAME(5,5,0:10,3) = (/0.00_EB,22.76_EB,30.44_EB,32.81_EB,35.19_EB,35.81_EB,36.43_EB,37.05_EB,37.67_EB,37.44_EB,37.21_EB/)
QFLAME(5,5,11:20,3) = (/36.98_EB,36.76_EB,36.36_EB,35.97_EB,35.58_EB,35.19_EB,34.95_EB,34.70_EB,34.46_EB,34.22_EB/)
QFLAME(5,5,21:30,3) = (/34.05_EB,33.87_EB,33.69_EB,33.52_EB,33.34_EB,32.96_EB,32.57_EB,32.19_EB,31.81_EB,31.42_EB/)
QFLAME(5,5,0:10,4) = (/0.00_EB,24.51_EB,33.83_EB,37.22_EB,40.61_EB,41.84_EB,43.08_EB,44.32_EB,45.56_EB,45.70_EB,45.84_EB/)
QFLAME(5,5,11:20,4) = (/45.99_EB,46.13_EB,46.00_EB,45.87_EB,45.73_EB,45.60_EB,45.60_EB,45.60_EB,45.61_EB,45.61_EB/)
QFLAME(5,5,21:30,4) = (/45.21_EB,44.82_EB,44.42_EB,44.02_EB,43.63_EB,43.16_EB,42.70_EB,42.23_EB,41.77_EB,41.31_EB/)
QFLAME(5,5,0:10,5) = (/0.00_EB,26.30_EB,37.47_EB,42.28_EB,47.08_EB,49.00_EB,50.93_EB,52.85_EB,54.77_EB,55.24_EB,55.71_EB/)
QFLAME(5,5,11:20,5) = (/56.17_EB,56.64_EB,56.80_EB,56.96_EB,57.12_EB,57.27_EB,56.76_EB,56.25_EB,55.74_EB,55.23_EB/)
QFLAME(5,5,21:30,5) = (/55.30_EB,55.37_EB,55.43_EB,55.50_EB,55.57_EB,55.28_EB,55.00_EB,54.72_EB,54.44_EB,54.16_EB/)
QFLAME(5,5,0:10,6) = (/0.00_EB,28.11_EB,40.94_EB,47.15_EB,53.37_EB,55.79_EB,58.22_EB,60.64_EB,63.06_EB,64.35_EB,65.64_EB/)
QFLAME(5,5,11:20,6) = (/66.93_EB,68.22_EB,68.63_EB,69.04_EB,69.45_EB,69.86_EB,69.28_EB,68.71_EB,68.14_EB,67.56_EB/)
QFLAME(5,5,21:30,6) = (/67.59_EB,67.61_EB,67.63_EB,67.65_EB,67.68_EB,67.53_EB,67.39_EB,67.25_EB,67.11_EB,66.96_EB/)
QFLAME(5,6,0:10,1) = (/0.00_EB,20.16_EB,27.01_EB,28.99_EB,30.96_EB,31.42_EB,31.89_EB,32.35_EB,32.81_EB,32.55_EB,32.28_EB/)
QFLAME(5,6,11:20,1) = (/32.02_EB,31.75_EB,31.32_EB,30.90_EB,30.47_EB,30.05_EB,29.68_EB,29.32_EB,28.95_EB,28.59_EB/)
QFLAME(5,6,21:30,1) = (/28.22_EB,27.85_EB,27.48_EB,27.12_EB,26.75_EB,26.42_EB,26.10_EB,25.77_EB,25.45_EB,25.12_EB/)
QFLAME(5,6,0:10,2) = (/0.00_EB,20.81_EB,27.41_EB,29.26_EB,31.11_EB,31.59_EB,32.07_EB,32.55_EB,33.03_EB,32.74_EB,32.45_EB/)
QFLAME(5,6,11:20,2) = (/32.17_EB,31.88_EB,31.43_EB,30.99_EB,30.54_EB,30.09_EB,29.74_EB,29.40_EB,29.05_EB,28.70_EB/)
QFLAME(5,6,21:30,2) = (/28.33_EB,27.95_EB,27.58_EB,27.21_EB,26.83_EB,26.50_EB,26.16_EB,25.82_EB,25.49_EB,25.15_EB/)
QFLAME(5,6,0:10,3) = (/0.00_EB,22.72_EB,30.41_EB,32.80_EB,35.20_EB,35.82_EB,36.44_EB,37.06_EB,37.68_EB,37.47_EB,37.25_EB/)
QFLAME(5,6,11:20,3) = (/37.04_EB,36.82_EB,36.44_EB,36.06_EB,35.68_EB,35.30_EB,35.06_EB,34.83_EB,34.59_EB,34.36_EB/)
QFLAME(5,6,21:30,3) = (/34.15_EB,33.94_EB,33.73_EB,33.53_EB,33.32_EB,32.93_EB,32.55_EB,32.16_EB,31.78_EB,31.39_EB/)
QFLAME(5,6,0:10,4) = (/0.00_EB,24.49_EB,33.80_EB,37.19_EB,40.58_EB,41.82_EB,43.06_EB,44.31_EB,45.55_EB,45.71_EB,45.87_EB/)
QFLAME(5,6,11:20,4) = (/46.03_EB,46.19_EB,46.04_EB,45.90_EB,45.75_EB,45.61_EB,45.55_EB,45.49_EB,45.44_EB,45.38_EB/)
QFLAME(5,6,21:30,4) = (/45.02_EB,44.67_EB,44.31_EB,43.95_EB,43.59_EB,43.17_EB,42.74_EB,42.32_EB,41.89_EB,41.46_EB/)
QFLAME(5,6,0:10,5) = (/0.00_EB,26.28_EB,37.44_EB,42.24_EB,47.04_EB,48.97_EB,50.90_EB,52.83_EB,54.76_EB,55.23_EB,55.70_EB/)
QFLAME(5,6,11:20,5) = (/56.17_EB,56.63_EB,56.79_EB,56.94_EB,57.10_EB,57.25_EB,56.74_EB,56.22_EB,55.71_EB,55.19_EB/)
QFLAME(5,6,21:30,5) = (/55.26_EB,55.33_EB,55.40_EB,55.46_EB,55.53_EB,55.25_EB,54.97_EB,54.68_EB,54.40_EB,54.12_EB/)
QFLAME(5,6,0:10,6) = (/0.00_EB,28.10_EB,40.90_EB,47.11_EB,53.32_EB,55.74_EB,58.17_EB,60.59_EB,63.02_EB,64.31_EB,65.60_EB/)
QFLAME(5,6,11:20,6) = (/66.89_EB,68.18_EB,68.59_EB,68.99_EB,69.40_EB,69.80_EB,69.24_EB,68.67_EB,68.11_EB,67.55_EB/)
QFLAME(5,6,21:30,6) = (/67.56_EB,67.58_EB,67.60_EB,67.62_EB,67.64_EB,67.50_EB,67.36_EB,67.21_EB,67.07_EB,66.93_EB/)
QFLAME(5,7,0:10,1) = (/0.00_EB,20.28_EB,27.28_EB,29.36_EB,31.43_EB,31.95_EB,32.47_EB,32.98_EB,33.50_EB,33.26_EB,33.03_EB/)
QFLAME(5,7,11:20,1) = (/32.80_EB,32.56_EB,32.15_EB,31.73_EB,31.31_EB,30.89_EB,30.54_EB,30.18_EB,29.83_EB,29.47_EB/)
QFLAME(5,7,21:30,1) = (/29.11_EB,28.74_EB,28.38_EB,28.01_EB,27.65_EB,27.33_EB,27.01_EB,26.68_EB,26.36_EB,26.04_EB/)
QFLAME(5,7,0:10,2) = (/0.00_EB,20.79_EB,27.58_EB,29.56_EB,31.54_EB,32.07_EB,32.60_EB,33.14_EB,33.67_EB,33.42_EB,33.16_EB/)
QFLAME(5,7,11:20,2) = (/32.91_EB,32.66_EB,32.22_EB,31.79_EB,31.35_EB,30.92_EB,30.58_EB,30.24_EB,29.90_EB,29.56_EB/)
QFLAME(5,7,21:30,2) = (/29.19_EB,28.82_EB,28.45_EB,28.08_EB,27.71_EB,27.38_EB,27.06_EB,26.73_EB,26.40_EB,26.08_EB/)
QFLAME(5,7,0:10,3) = (/0.00_EB,22.68_EB,30.38_EB,32.79_EB,35.20_EB,35.83_EB,36.45_EB,37.07_EB,37.70_EB,37.49_EB,37.29_EB/)
QFLAME(5,7,11:20,3) = (/37.09_EB,36.89_EB,36.52_EB,36.15_EB,35.78_EB,35.41_EB,35.18_EB,34.95_EB,34.72_EB,34.50_EB/)
QFLAME(5,7,21:30,3) = (/34.25_EB,34.01_EB,33.77_EB,33.53_EB,33.29_EB,32.91_EB,32.52_EB,32.13_EB,31.75_EB,31.36_EB/)
QFLAME(5,7,0:10,4) = (/0.00_EB,24.48_EB,33.77_EB,37.16_EB,40.56_EB,41.80_EB,43.05_EB,44.29_EB,45.54_EB,45.71_EB,45.89_EB/)
QFLAME(5,7,11:20,4) = (/46.07_EB,46.24_EB,46.09_EB,45.93_EB,45.77_EB,45.62_EB,45.50_EB,45.38_EB,45.27_EB,45.15_EB/)
QFLAME(5,7,21:30,4) = (/44.83_EB,44.51_EB,44.20_EB,43.88_EB,43.56_EB,43.17_EB,42.78_EB,42.40_EB,42.01_EB,41.62_EB/)
QFLAME(5,7,0:10,5) = (/0.00_EB,26.27_EB,37.40_EB,42.20_EB,47.00_EB,48.94_EB,50.87_EB,52.81_EB,54.74_EB,55.22_EB,55.69_EB/)
QFLAME(5,7,11:20,5) = (/56.16_EB,56.63_EB,56.78_EB,56.93_EB,57.08_EB,57.23_EB,56.71_EB,56.19_EB,55.67_EB,55.15_EB/)
QFLAME(5,7,21:30,5) = (/55.22_EB,55.29_EB,55.36_EB,55.43_EB,55.50_EB,55.21_EB,54.93_EB,54.64_EB,54.36_EB,54.08_EB/)
QFLAME(5,7,0:10,6) = (/0.00_EB,28.08_EB,40.86_EB,47.07_EB,53.27_EB,55.70_EB,58.12_EB,60.55_EB,62.97_EB,64.27_EB,65.56_EB/)
QFLAME(5,7,11:20,6) = (/66.85_EB,68.15_EB,68.55_EB,68.95_EB,69.35_EB,69.75_EB,69.19_EB,68.64_EB,68.08_EB,67.53_EB/)
QFLAME(5,7,21:30,6) = (/67.54_EB,67.56_EB,67.58_EB,67.59_EB,67.61_EB,67.46_EB,67.32_EB,67.18_EB,67.03_EB,66.89_EB/)
QFLAME(5,8,0:10,1) = (/0.00_EB,20.41_EB,27.55_EB,29.73_EB,31.91_EB,32.48_EB,33.05_EB,33.62_EB,34.19_EB,33.98_EB,33.78_EB/)
QFLAME(5,8,11:20,1) = (/33.58_EB,33.38_EB,32.97_EB,32.55_EB,32.14_EB,31.73_EB,31.39_EB,31.04_EB,30.70_EB,30.35_EB/)
QFLAME(5,8,21:30,1) = (/29.99_EB,29.63_EB,29.27_EB,28.91_EB,28.55_EB,28.23_EB,27.91_EB,27.59_EB,27.28_EB,26.96_EB/)
QFLAME(5,8,0:10,2) = (/0.00_EB,20.78_EB,27.75_EB,29.86_EB,31.97_EB,32.55_EB,33.14_EB,33.72_EB,34.31_EB,34.09_EB,33.87_EB/)
QFLAME(5,8,11:20,2) = (/33.65_EB,33.43_EB,33.01_EB,32.59_EB,32.17_EB,31.75_EB,31.42_EB,31.08_EB,30.75_EB,30.42_EB/)
QFLAME(5,8,21:30,2) = (/30.06_EB,29.69_EB,29.32_EB,28.96_EB,28.59_EB,28.27_EB,27.96_EB,27.64_EB,27.32_EB,27.00_EB/)
QFLAME(5,8,0:10,3) = (/0.00_EB,22.63_EB,30.36_EB,32.78_EB,35.21_EB,35.83_EB,36.46_EB,37.08_EB,37.71_EB,37.52_EB,37.33_EB/)
QFLAME(5,8,11:20,3) = (/37.14_EB,36.95_EB,36.59_EB,36.24_EB,35.88_EB,35.53_EB,35.30_EB,35.08_EB,34.86_EB,34.63_EB/)
QFLAME(5,8,21:30,3) = (/34.36_EB,34.09_EB,33.82_EB,33.54_EB,33.27_EB,32.88_EB,32.49_EB,32.11_EB,31.72_EB,31.33_EB/)
QFLAME(5,8,0:10,4) = (/0.00_EB,24.46_EB,33.74_EB,37.14_EB,40.53_EB,41.78_EB,43.03_EB,44.28_EB,45.53_EB,45.72_EB,45.91_EB/)
QFLAME(5,8,11:20,4) = (/46.11_EB,46.30_EB,46.13_EB,45.96_EB,45.79_EB,45.63_EB,45.45_EB,45.28_EB,45.10_EB,44.92_EB/)
QFLAME(5,8,21:30,4) = (/44.64_EB,44.36_EB,44.08_EB,43.81_EB,43.53_EB,43.18_EB,42.83_EB,42.48_EB,42.13_EB,41.78_EB/)
QFLAME(5,8,0:10,5) = (/0.00_EB,26.25_EB,37.36_EB,42.16_EB,46.97_EB,48.91_EB,50.85_EB,52.79_EB,54.73_EB,55.20_EB,55.68_EB/)
QFLAME(5,8,11:20,5) = (/56.15_EB,56.62_EB,56.77_EB,56.92_EB,57.06_EB,57.21_EB,56.69_EB,56.16_EB,55.64_EB,55.11_EB/)
QFLAME(5,8,21:30,5) = (/55.18_EB,55.25_EB,55.32_EB,55.39_EB,55.46_EB,55.18_EB,54.89_EB,54.61_EB,54.32_EB,54.03_EB/)
QFLAME(5,8,0:10,6) = (/0.00_EB,28.07_EB,40.83_EB,47.02_EB,53.22_EB,55.65_EB,58.07_EB,60.50_EB,62.93_EB,64.22_EB,65.52_EB/)
QFLAME(5,8,11:20,6) = (/66.81_EB,68.11_EB,68.51_EB,68.90_EB,69.30_EB,69.69_EB,69.15_EB,68.60_EB,68.05_EB,67.51_EB/)
QFLAME(5,8,21:30,6) = (/67.52_EB,67.54_EB,67.55_EB,67.56_EB,67.58_EB,67.43_EB,67.29_EB,67.14_EB,67.00_EB,66.85_EB/)
QFLAME(5,9,0:10,1) = (/0.00_EB,20.53_EB,27.82_EB,30.10_EB,32.38_EB,33.00_EB,33.63_EB,34.25_EB,34.87_EB,34.70_EB,34.53_EB/)
QFLAME(5,9,11:20,1) = (/34.36_EB,34.19_EB,33.79_EB,33.38_EB,32.98_EB,32.58_EB,32.24_EB,31.91_EB,31.57_EB,31.24_EB/)
QFLAME(5,9,21:30,1) = (/30.88_EB,30.52_EB,30.16_EB,29.81_EB,29.45_EB,29.13_EB,28.82_EB,28.50_EB,28.19_EB,27.87_EB/)
QFLAME(5,9,0:10,2) = (/0.00_EB,20.76_EB,27.92_EB,30.16_EB,32.40_EB,33.04_EB,33.67_EB,34.31_EB,34.95_EB,34.76_EB,34.58_EB/)
QFLAME(5,9,11:20,2) = (/34.39_EB,34.21_EB,33.80_EB,33.39_EB,32.98_EB,32.57_EB,32.25_EB,31.93_EB,31.61_EB,31.29_EB/)
QFLAME(5,9,21:30,2) = (/30.92_EB,30.56_EB,30.20_EB,29.83_EB,29.47_EB,29.16_EB,28.85_EB,28.55_EB,28.24_EB,27.93_EB/)
QFLAME(5,9,0:10,3) = (/0.00_EB,22.59_EB,30.33_EB,32.77_EB,35.21_EB,35.84_EB,36.47_EB,37.10_EB,37.72_EB,37.55_EB,37.37_EB/)
QFLAME(5,9,11:20,3) = (/37.19_EB,37.02_EB,36.67_EB,36.33_EB,35.98_EB,35.64_EB,35.42_EB,35.20_EB,34.99_EB,34.77_EB/)
QFLAME(5,9,21:30,3) = (/34.46_EB,34.16_EB,33.86_EB,33.55_EB,33.25_EB,32.86_EB,32.47_EB,32.08_EB,31.69_EB,31.30_EB/)
QFLAME(5,9,0:10,4) = (/0.00_EB,24.44_EB,33.70_EB,37.11_EB,40.51_EB,41.76_EB,43.02_EB,44.27_EB,45.52_EB,45.73_EB,45.94_EB/)
QFLAME(5,9,11:20,4) = (/46.14_EB,46.35_EB,46.17_EB,45.99_EB,45.82_EB,45.64_EB,45.40_EB,45.17_EB,44.93_EB,44.70_EB/)
QFLAME(5,9,21:30,4) = (/44.45_EB,44.21_EB,43.97_EB,43.73_EB,43.49_EB,43.18_EB,42.87_EB,42.56_EB,42.25_EB,41.94_EB/)
QFLAME(5,9,0:10,5) = (/0.00_EB,26.23_EB,37.32_EB,42.13_EB,46.93_EB,48.88_EB,50.82_EB,52.77_EB,54.72_EB,55.19_EB,55.67_EB/)
QFLAME(5,9,11:20,5) = (/56.14_EB,56.62_EB,56.76_EB,56.90_EB,57.04_EB,57.19_EB,56.66_EB,56.13_EB,55.60_EB,55.07_EB/)
QFLAME(5,9,21:30,5) = (/55.14_EB,55.21_EB,55.29_EB,55.36_EB,55.43_EB,55.14_EB,54.85_EB,54.57_EB,54.28_EB,53.99_EB/)
QFLAME(5,9,0:10,6) = (/0.00_EB,28.06_EB,40.79_EB,46.98_EB,53.17_EB,55.60_EB,58.03_EB,60.45_EB,62.88_EB,64.18_EB,65.48_EB/)
QFLAME(5,9,11:20,6) = (/66.78_EB,68.07_EB,68.47_EB,68.86_EB,69.25_EB,69.64_EB,69.10_EB,68.57_EB,68.03_EB,67.49_EB/)
QFLAME(5,9,21:30,6) = (/67.50_EB,67.51_EB,67.52_EB,67.53_EB,67.54_EB,67.40_EB,67.25_EB,67.11_EB,66.96_EB,66.82_EB/)
QFLAME(5,10,0:10,1) = (/0.00_EB,20.66_EB,28.08_EB,30.47_EB,32.85_EB,33.53_EB,34.21_EB,34.88_EB,35.56_EB,35.42_EB,35.28_EB/)
QFLAME(5,10,11:20,1) = (/35.14_EB,35.00_EB,34.61_EB,34.21_EB,33.81_EB,33.42_EB,33.09_EB,32.77_EB,32.45_EB,32.12_EB/)
QFLAME(5,10,21:30,1) = (/31.77_EB,31.41_EB,31.06_EB,30.70_EB,30.35_EB,30.04_EB,29.72_EB,29.41_EB,29.10_EB,28.79_EB/)
QFLAME(5,10,0:10,2) = (/0.00_EB,20.74_EB,28.09_EB,30.46_EB,32.83_EB,33.52_EB,34.21_EB,34.90_EB,35.59_EB,35.44_EB,35.29_EB/)
QFLAME(5,10,11:20,2) = (/35.14_EB,34.99_EB,34.59_EB,34.19_EB,33.80_EB,33.40_EB,33.09_EB,32.77_EB,32.46_EB,32.15_EB/)
QFLAME(5,10,21:30,2) = (/31.79_EB,31.43_EB,31.07_EB,30.71_EB,30.35_EB,30.05_EB,29.75_EB,29.45_EB,29.15_EB,28.85_EB/)
QFLAME(5,10,0:10,3) = (/0.00_EB,22.55_EB,30.30_EB,32.76_EB,35.22_EB,35.85_EB,36.48_EB,37.11_EB,37.74_EB,37.57_EB,37.41_EB/)
QFLAME(5,10,11:20,3) = (/37.24_EB,37.08_EB,36.75_EB,36.42_EB,36.09_EB,35.75_EB,35.54_EB,35.33_EB,35.12_EB,34.90_EB/)
QFLAME(5,10,21:30,3) = (/34.57_EB,34.23_EB,33.90_EB,33.56_EB,33.22_EB,32.83_EB,32.44_EB,32.05_EB,31.66_EB,31.27_EB/)
QFLAME(5,10,0:10,4) = (/0.00_EB,24.42_EB,33.67_EB,37.08_EB,40.49_EB,41.74_EB,43.00_EB,44.26_EB,45.51_EB,45.74_EB,45.96_EB/)
QFLAME(5,10,11:20,4) = (/46.18_EB,46.41_EB,46.22_EB,46.03_EB,45.84_EB,45.65_EB,45.35_EB,45.06_EB,44.76_EB,44.47_EB/)
QFLAME(5,10,21:30,4) = (/44.27_EB,44.06_EB,43.86_EB,43.66_EB,43.46_EB,43.19_EB,42.91_EB,42.64_EB,42.37_EB,42.09_EB/)
QFLAME(5,10,0:10,5) = (/0.00_EB,26.22_EB,37.29_EB,42.09_EB,46.89_EB,48.85_EB,50.80_EB,52.75_EB,54.70_EB,55.18_EB,55.66_EB/)
QFLAME(5,10,11:20,5) = (/56.13_EB,56.61_EB,56.75_EB,56.89_EB,57.03_EB,57.17_EB,56.63_EB,56.10_EB,55.57_EB,55.03_EB/)
QFLAME(5,10,21:30,5) = (/55.10_EB,55.18_EB,55.25_EB,55.32_EB,55.39_EB,55.10_EB,54.82_EB,54.53_EB,54.24_EB,53.95_EB/)
QFLAME(5,10,0:10,6) = (/0.00_EB,28.04_EB,40.75_EB,46.94_EB,53.12_EB,55.55_EB,57.98_EB,60.41_EB,62.84_EB,64.14_EB,65.44_EB/)
QFLAME(5,10,11:20,6) = (/66.74_EB,68.04_EB,68.43_EB,68.81_EB,69.20_EB,69.59_EB,69.06_EB,68.53_EB,68.00_EB,67.47_EB/)
QFLAME(5,10,21:30,6) = (/67.48_EB,67.49_EB,67.49_EB,67.50_EB,67.51_EB,67.36_EB,67.22_EB,67.07_EB,66.93_EB,66.78_EB/)
QFLAME(5,11,0:10,1) = (/0.00_EB,20.73_EB,28.24_EB,30.69_EB,33.13_EB,33.84_EB,34.55_EB,35.26_EB,35.98_EB,35.85_EB,35.73_EB/)
QFLAME(5,11,11:20,1) = (/35.61_EB,35.49_EB,35.11_EB,34.73_EB,34.34_EB,33.96_EB,33.64_EB,33.32_EB,33.00_EB,32.68_EB/)
QFLAME(5,11,21:30,1) = (/32.33_EB,31.98_EB,31.63_EB,31.27_EB,30.92_EB,30.62_EB,30.31_EB,30.01_EB,29.70_EB,29.40_EB/)
QFLAME(5,11,0:10,2) = (/0.00_EB,20.80_EB,28.25_EB,30.68_EB,33.11_EB,33.83_EB,34.55_EB,35.27_EB,35.98_EB,35.86_EB,35.73_EB/)
QFLAME(5,11,11:20,2) = (/35.61_EB,35.48_EB,35.10_EB,34.71_EB,34.33_EB,33.94_EB,33.63_EB,33.32_EB,33.01_EB,32.70_EB/)
QFLAME(5,11,21:30,2) = (/32.35_EB,31.99_EB,31.64_EB,31.28_EB,30.93_EB,30.63_EB,30.34_EB,30.05_EB,29.75_EB,29.46_EB/)
QFLAME(5,11,0:10,3) = (/0.00_EB,22.52_EB,30.29_EB,32.79_EB,35.29_EB,35.96_EB,36.63_EB,37.30_EB,37.97_EB,37.83_EB,37.69_EB/)
QFLAME(5,11,11:20,3) = (/37.55_EB,37.41_EB,37.08_EB,36.74_EB,36.41_EB,36.07_EB,35.85_EB,35.64_EB,35.42_EB,35.20_EB/)
QFLAME(5,11,21:30,3) = (/34.87_EB,34.53_EB,34.20_EB,33.86_EB,33.53_EB,33.16_EB,32.78_EB,32.41_EB,32.03_EB,31.66_EB/)
QFLAME(5,11,0:10,4) = (/0.00_EB,24.39_EB,33.64_EB,37.05_EB,40.46_EB,41.73_EB,42.99_EB,44.26_EB,45.53_EB,45.74_EB,45.96_EB/)
QFLAME(5,11,11:20,4) = (/46.18_EB,46.40_EB,46.21_EB,46.02_EB,45.83_EB,45.64_EB,45.35_EB,45.06_EB,44.77_EB,44.48_EB/)
QFLAME(5,11,21:30,4) = (/44.28_EB,44.08_EB,43.89_EB,43.69_EB,43.49_EB,43.23_EB,42.96_EB,42.70_EB,42.43_EB,42.17_EB/)
QFLAME(5,11,0:10,5) = (/0.00_EB,26.20_EB,37.25_EB,42.05_EB,46.86_EB,48.81_EB,50.76_EB,52.71_EB,54.66_EB,55.15_EB,55.63_EB/)
QFLAME(5,11,11:20,5) = (/56.11_EB,56.59_EB,56.73_EB,56.87_EB,57.01_EB,57.15_EB,56.63_EB,56.11_EB,55.58_EB,55.06_EB/)
QFLAME(5,11,21:30,5) = (/55.14_EB,55.21_EB,55.28_EB,55.36_EB,55.43_EB,55.14_EB,54.84_EB,54.55_EB,54.25_EB,53.96_EB/)
QFLAME(5,11,0:10,6) = (/0.00_EB,28.03_EB,40.71_EB,46.90_EB,53.08_EB,55.51_EB,57.94_EB,60.37_EB,62.80_EB,64.09_EB,65.39_EB/)
QFLAME(5,11,11:20,6) = (/66.68_EB,67.98_EB,68.34_EB,68.71_EB,69.07_EB,69.44_EB,68.94_EB,68.45_EB,67.95_EB,67.45_EB/)
QFLAME(5,11,21:30,6) = (/67.46_EB,67.47_EB,67.48_EB,67.49_EB,67.50_EB,67.36_EB,67.23_EB,67.09_EB,66.95_EB,66.82_EB/)
QFLAME(5,12,0:10,1) = (/0.00_EB,20.80_EB,28.40_EB,30.91_EB,33.41_EB,34.15_EB,34.90_EB,35.64_EB,36.39_EB,36.29_EB,36.19_EB/)
QFLAME(5,12,11:20,1) = (/36.09_EB,35.99_EB,35.61_EB,35.24_EB,34.87_EB,34.49_EB,34.18_EB,33.87_EB,33.56_EB,33.25_EB/)
QFLAME(5,12,21:30,1) = (/32.90_EB,32.55_EB,32.20_EB,31.85_EB,31.49_EB,31.20_EB,30.90_EB,30.60_EB,30.30_EB,30.00_EB/)
QFLAME(5,12,0:10,2) = (/0.00_EB,20.86_EB,28.41_EB,30.90_EB,33.39_EB,34.14_EB,34.88_EB,35.63_EB,36.38_EB,36.28_EB,36.18_EB/)
QFLAME(5,12,11:20,2) = (/36.08_EB,35.98_EB,35.61_EB,35.23_EB,34.86_EB,34.48_EB,34.18_EB,33.87_EB,33.56_EB,33.26_EB/)
QFLAME(5,12,21:30,2) = (/32.91_EB,32.56_EB,32.20_EB,31.85_EB,31.50_EB,31.21_EB,30.92_EB,30.64_EB,30.35_EB,30.06_EB/)
QFLAME(5,12,0:10,3) = (/0.00_EB,22.49_EB,30.28_EB,32.82_EB,35.35_EB,36.06_EB,36.77_EB,37.48_EB,38.19_EB,38.08_EB,37.97_EB/)
QFLAME(5,12,11:20,3) = (/37.85_EB,37.74_EB,37.40_EB,37.06_EB,36.73_EB,36.39_EB,36.17_EB,35.95_EB,35.72_EB,35.50_EB/)
QFLAME(5,12,21:30,3) = (/35.17_EB,34.83_EB,34.50_EB,34.17_EB,33.84_EB,33.48_EB,33.12_EB,32.77_EB,32.41_EB,32.05_EB/)
QFLAME(5,12,0:10,4) = (/0.00_EB,24.37_EB,33.60_EB,37.02_EB,40.43_EB,41.71_EB,42.99_EB,44.26_EB,45.54_EB,45.75_EB,45.97_EB/)
QFLAME(5,12,11:20,4) = (/46.18_EB,46.39_EB,46.20_EB,46.01_EB,45.82_EB,45.63_EB,45.34_EB,45.06_EB,44.78_EB,44.49_EB/)
QFLAME(5,12,21:30,4) = (/44.30_EB,44.11_EB,43.91_EB,43.72_EB,43.52_EB,43.27_EB,43.01_EB,42.76_EB,42.50_EB,42.25_EB/)
QFLAME(5,12,0:10,5) = (/0.00_EB,26.18_EB,37.21_EB,42.02_EB,46.82_EB,48.77_EB,50.72_EB,52.67_EB,54.62_EB,55.11_EB,55.60_EB/)
QFLAME(5,12,11:20,5) = (/56.09_EB,56.57_EB,56.71_EB,56.85_EB,56.99_EB,57.14_EB,56.62_EB,56.11_EB,55.60_EB,55.09_EB/)
QFLAME(5,12,21:30,5) = (/55.17_EB,55.24_EB,55.32_EB,55.40_EB,55.47_EB,55.17_EB,54.87_EB,54.57_EB,54.27_EB,53.97_EB/)
QFLAME(5,12,0:10,6) = (/0.00_EB,28.01_EB,40.68_EB,46.86_EB,53.04_EB,55.47_EB,57.90_EB,60.33_EB,62.76_EB,64.05_EB,65.34_EB/)
QFLAME(5,12,11:20,6) = (/66.63_EB,67.92_EB,68.26_EB,68.61_EB,68.95_EB,69.29_EB,68.83_EB,68.36_EB,67.90_EB,67.43_EB/)
QFLAME(5,12,21:30,6) = (/67.44_EB,67.46_EB,67.47_EB,67.48_EB,67.49_EB,67.36_EB,67.24_EB,67.11_EB,66.98_EB,66.85_EB/)
QFLAME(5,13,0:10,1) = (/0.00_EB,20.87_EB,28.56_EB,31.12_EB,33.68_EB,34.46_EB,35.24_EB,36.02_EB,36.80_EB,36.72_EB,36.64_EB/)
QFLAME(5,13,11:20,1) = (/36.56_EB,36.48_EB,36.12_EB,35.75_EB,35.39_EB,35.03_EB,34.73_EB,34.42_EB,34.12_EB,33.81_EB/)
QFLAME(5,13,21:30,1) = (/33.46_EB,33.11_EB,32.76_EB,32.42_EB,32.07_EB,31.78_EB,31.48_EB,31.19_EB,30.90_EB,30.61_EB/)
QFLAME(5,13,0:10,2) = (/0.00_EB,20.93_EB,28.57_EB,31.12_EB,33.67_EB,34.45_EB,35.22_EB,36.00_EB,36.78_EB,36.70_EB,36.63_EB/)
QFLAME(5,13,11:20,2) = (/36.55_EB,36.48_EB,36.11_EB,35.75_EB,35.39_EB,35.02_EB,34.72_EB,34.42_EB,34.12_EB,33.82_EB/)
QFLAME(5,13,21:30,2) = (/33.47_EB,33.12_EB,32.77_EB,32.42_EB,32.07_EB,31.79_EB,31.51_EB,31.23_EB,30.95_EB,30.67_EB/)
QFLAME(5,13,0:10,3) = (/0.00_EB,22.47_EB,30.27_EB,32.84_EB,35.41_EB,36.17_EB,36.92_EB,37.67_EB,38.42_EB,38.33_EB,38.25_EB/)
QFLAME(5,13,11:20,3) = (/38.16_EB,38.07_EB,37.73_EB,37.39_EB,37.05_EB,36.71_EB,36.48_EB,36.25_EB,36.03_EB,35.80_EB/)
QFLAME(5,13,21:30,3) = (/35.47_EB,35.14_EB,34.80_EB,34.47_EB,34.14_EB,33.80_EB,33.46_EB,33.12_EB,32.79_EB,32.45_EB/)
QFLAME(5,13,0:10,4) = (/0.00_EB,24.34_EB,33.56_EB,36.99_EB,40.41_EB,41.69_EB,42.98_EB,44.27_EB,45.55_EB,45.76_EB,45.97_EB/)
QFLAME(5,13,11:20,4) = (/46.18_EB,46.39_EB,46.19_EB,46.00_EB,45.81_EB,45.62_EB,45.34_EB,45.06_EB,44.79_EB,44.51_EB/)
QFLAME(5,13,21:30,4) = (/44.32_EB,44.13_EB,43.94_EB,43.75_EB,43.56_EB,43.31_EB,43.06_EB,42.82_EB,42.57_EB,42.32_EB/)
QFLAME(5,13,0:10,5) = (/0.00_EB,26.15_EB,37.18_EB,41.98_EB,46.78_EB,48.73_EB,50.68_EB,52.64_EB,54.59_EB,55.08_EB,55.57_EB/)
QFLAME(5,13,11:20,5) = (/56.06_EB,56.55_EB,56.69_EB,56.84_EB,56.98_EB,57.12_EB,56.62_EB,56.12_EB,55.62_EB,55.12_EB/)
QFLAME(5,13,21:30,5) = (/55.20_EB,55.28_EB,55.35_EB,55.43_EB,55.51_EB,55.20_EB,54.90_EB,54.59_EB,54.29_EB,53.98_EB/)
QFLAME(5,13,0:10,6) = (/0.00_EB,28.00_EB,40.64_EB,46.82_EB,53.01_EB,55.43_EB,57.86_EB,60.29_EB,62.72_EB,64.00_EB,65.29_EB/)
QFLAME(5,13,11:20,6) = (/66.57_EB,67.86_EB,68.18_EB,68.50_EB,68.83_EB,69.15_EB,68.71_EB,68.28_EB,67.85_EB,67.41_EB/)
QFLAME(5,13,21:30,6) = (/67.43_EB,67.44_EB,67.46_EB,67.47_EB,67.48_EB,67.36_EB,67.25_EB,67.13_EB,67.01_EB,66.89_EB/)
QFLAME(5,14,0:10,1) = (/0.00_EB,20.94_EB,28.72_EB,31.34_EB,33.96_EB,34.77_EB,35.59_EB,36.40_EB,37.21_EB,37.15_EB,37.09_EB/)
QFLAME(5,14,11:20,1) = (/37.03_EB,36.97_EB,36.62_EB,36.27_EB,35.92_EB,35.57_EB,35.27_EB,34.97_EB,34.67_EB,34.37_EB/)
QFLAME(5,14,21:30,1) = (/34.03_EB,33.68_EB,33.33_EB,32.99_EB,32.64_EB,32.36_EB,32.07_EB,31.78_EB,31.50_EB,31.21_EB/)
QFLAME(5,14,0:10,2) = (/0.00_EB,20.99_EB,28.73_EB,31.34_EB,33.95_EB,34.75_EB,35.56_EB,36.37_EB,37.18_EB,37.13_EB,37.08_EB/)
QFLAME(5,14,11:20,2) = (/37.02_EB,36.97_EB,36.62_EB,36.27_EB,35.91_EB,35.56_EB,35.26_EB,34.97_EB,34.67_EB,34.37_EB/)
QFLAME(5,14,21:30,2) = (/34.03_EB,33.68_EB,33.34_EB,32.99_EB,32.65_EB,32.37_EB,32.10_EB,31.82_EB,31.55_EB,31.27_EB/)
QFLAME(5,14,0:10,3) = (/0.00_EB,22.44_EB,30.26_EB,32.87_EB,35.48_EB,36.27_EB,37.06_EB,37.86_EB,38.65_EB,38.59_EB,38.53_EB/)
QFLAME(5,14,11:20,3) = (/38.46_EB,38.40_EB,38.06_EB,37.71_EB,37.37_EB,37.03_EB,36.79_EB,36.56_EB,36.33_EB,36.10_EB/)
QFLAME(5,14,21:30,3) = (/35.77_EB,35.44_EB,35.11_EB,34.78_EB,34.45_EB,34.13_EB,33.80_EB,33.48_EB,33.16_EB,32.84_EB/)
QFLAME(5,14,0:10,4) = (/0.00_EB,24.32_EB,33.53_EB,36.95_EB,40.38_EB,41.68_EB,42.97_EB,44.27_EB,45.57_EB,45.77_EB,45.97_EB/)
QFLAME(5,14,11:20,4) = (/46.17_EB,46.38_EB,46.19_EB,45.99_EB,45.80_EB,45.61_EB,45.34_EB,45.07_EB,44.79_EB,44.52_EB/)
QFLAME(5,14,21:30,4) = (/44.33_EB,44.15_EB,43.96_EB,43.77_EB,43.59_EB,43.35_EB,43.11_EB,42.87_EB,42.64_EB,42.40_EB/)
QFLAME(5,14,0:10,5) = (/0.00_EB,26.13_EB,37.14_EB,41.94_EB,46.74_EB,48.69_EB,50.65_EB,52.60_EB,54.55_EB,55.04_EB,55.54_EB/)
QFLAME(5,14,11:20,5) = (/56.04_EB,56.53_EB,56.68_EB,56.82_EB,56.96_EB,57.11_EB,56.62_EB,56.13_EB,55.64_EB,55.15_EB/)
QFLAME(5,14,21:30,5) = (/55.23_EB,55.31_EB,55.39_EB,55.47_EB,55.55_EB,55.24_EB,54.93_EB,54.62_EB,54.30_EB,53.99_EB/)
QFLAME(5,14,0:10,6) = (/0.00_EB,27.98_EB,40.60_EB,46.79_EB,52.97_EB,55.39_EB,57.82_EB,60.25_EB,62.68_EB,63.96_EB,65.24_EB/)
QFLAME(5,14,11:20,6) = (/66.52_EB,67.80_EB,68.10_EB,68.40_EB,68.70_EB,69.00_EB,68.60_EB,68.20_EB,67.80_EB,67.39_EB/)
QFLAME(5,14,21:30,6) = (/67.41_EB,67.43_EB,67.44_EB,67.46_EB,67.48_EB,67.37_EB,67.26_EB,67.15_EB,67.04_EB,66.93_EB/)
QFLAME(5,15,0:10,1) = (/0.00_EB,21.01_EB,28.89_EB,31.56_EB,34.24_EB,35.09_EB,35.93_EB,36.78_EB,37.62_EB,37.58_EB,37.54_EB/)
QFLAME(5,15,11:20,1) = (/37.50_EB,37.46_EB,37.12_EB,36.78_EB,36.44_EB,36.11_EB,35.81_EB,35.52_EB,35.23_EB,34.94_EB/)
QFLAME(5,15,21:30,1) = (/34.59_EB,34.25_EB,33.90_EB,33.56_EB,33.21_EB,32.93_EB,32.66_EB,32.38_EB,32.10_EB,31.82_EB/)
QFLAME(5,15,0:10,2) = (/0.00_EB,21.05_EB,28.89_EB,31.56_EB,34.23_EB,35.06_EB,35.90_EB,36.74_EB,37.57_EB,37.55_EB,37.52_EB/)
QFLAME(5,15,11:20,2) = (/37.50_EB,37.47_EB,37.13_EB,36.79_EB,36.44_EB,36.10_EB,35.81_EB,35.51_EB,35.22_EB,34.93_EB/)
QFLAME(5,15,21:30,2) = (/34.59_EB,34.24_EB,33.90_EB,33.56_EB,33.22_EB,32.95_EB,32.68_EB,32.41_EB,32.14_EB,31.87_EB/)
QFLAME(5,15,0:10,3) = (/0.00_EB,22.41_EB,30.25_EB,32.90_EB,35.54_EB,36.38_EB,37.21_EB,38.05_EB,38.88_EB,38.84_EB,38.80_EB/)
QFLAME(5,15,11:20,3) = (/38.77_EB,38.73_EB,38.38_EB,38.04_EB,37.69_EB,37.34_EB,37.11_EB,36.87_EB,36.63_EB,36.39_EB/)
QFLAME(5,15,21:30,3) = (/36.07_EB,35.74_EB,35.41_EB,35.08_EB,34.75_EB,34.45_EB,34.14_EB,33.84_EB,33.54_EB,33.23_EB/)
QFLAME(5,15,0:10,4) = (/0.00_EB,24.30_EB,33.49_EB,36.92_EB,40.36_EB,41.66_EB,42.97_EB,44.27_EB,45.58_EB,45.78_EB,45.97_EB/)
QFLAME(5,15,11:20,4) = (/46.17_EB,46.37_EB,46.18_EB,45.99_EB,45.79_EB,45.60_EB,45.33_EB,45.07_EB,44.80_EB,44.54_EB/)
QFLAME(5,15,21:30,4) = (/44.35_EB,44.17_EB,43.99_EB,43.80_EB,43.62_EB,43.39_EB,43.16_EB,42.93_EB,42.71_EB,42.48_EB/)
QFLAME(5,15,0:10,5) = (/0.00_EB,26.11_EB,37.10_EB,41.91_EB,46.71_EB,48.66_EB,50.61_EB,52.56_EB,54.51_EB,55.01_EB,55.51_EB/)
QFLAME(5,15,11:20,5) = (/56.01_EB,56.51_EB,56.66_EB,56.80_EB,56.95_EB,57.09_EB,56.61_EB,56.13_EB,55.66_EB,55.18_EB/)
QFLAME(5,15,21:30,5) = (/55.26_EB,55.34_EB,55.42_EB,55.51_EB,55.59_EB,55.27_EB,54.95_EB,54.64_EB,54.32_EB,54.00_EB/)
QFLAME(5,15,0:10,6) = (/0.00_EB,27.96_EB,40.57_EB,46.75_EB,52.93_EB,55.36_EB,57.78_EB,60.21_EB,62.63_EB,63.91_EB,65.19_EB/)
QFLAME(5,15,11:20,6) = (/66.47_EB,67.74_EB,68.02_EB,68.30_EB,68.58_EB,68.86_EB,68.49_EB,68.11_EB,67.74_EB,67.37_EB/)
QFLAME(5,15,21:30,6) = (/67.39_EB,67.41_EB,67.43_EB,67.45_EB,67.47_EB,67.37_EB,67.27_EB,67.16_EB,67.06_EB,66.96_EB/)
QFLAME(5,16,0:10,1) = (/0.00_EB,21.08_EB,29.05_EB,31.78_EB,34.52_EB,35.40_EB,36.28_EB,37.16_EB,38.04_EB,38.02_EB,38.00_EB/)
QFLAME(5,16,11:20,1) = (/37.97_EB,37.95_EB,37.63_EB,37.30_EB,36.97_EB,36.64_EB,36.36_EB,36.07_EB,35.79_EB,35.50_EB/)
QFLAME(5,16,21:30,1) = (/35.16_EB,34.82_EB,34.47_EB,34.13_EB,33.79_EB,33.51_EB,33.24_EB,32.97_EB,32.70_EB,32.43_EB/)
QFLAME(5,16,0:10,2) = (/0.00_EB,21.11_EB,29.05_EB,31.78_EB,34.51_EB,35.37_EB,36.24_EB,37.11_EB,37.97_EB,37.97_EB,37.97_EB/)
QFLAME(5,16,11:20,2) = (/37.97_EB,37.97_EB,37.63_EB,37.30_EB,36.97_EB,36.64_EB,36.35_EB,36.06_EB,35.77_EB,35.48_EB/)
QFLAME(5,16,21:30,2) = (/35.15_EB,34.81_EB,34.47_EB,34.13_EB,33.79_EB,33.53_EB,33.27_EB,33.00_EB,32.74_EB,32.48_EB/)
QFLAME(5,16,0:10,3) = (/0.00_EB,22.39_EB,30.24_EB,32.92_EB,35.61_EB,36.48_EB,37.36_EB,38.23_EB,39.11_EB,39.10_EB,39.08_EB/)
QFLAME(5,16,11:20,3) = (/39.07_EB,39.06_EB,38.71_EB,38.36_EB,38.01_EB,37.66_EB,37.42_EB,37.18_EB,36.93_EB,36.69_EB/)
QFLAME(5,16,21:30,3) = (/36.37_EB,36.04_EB,35.71_EB,35.39_EB,35.06_EB,34.77_EB,34.49_EB,34.20_EB,33.91_EB,33.62_EB/)
QFLAME(5,16,0:10,4) = (/0.00_EB,24.27_EB,33.45_EB,36.89_EB,40.33_EB,41.64_EB,42.96_EB,44.28_EB,45.59_EB,45.79_EB,45.98_EB/)
QFLAME(5,16,11:20,4) = (/46.17_EB,46.36_EB,46.17_EB,45.98_EB,45.78_EB,45.59_EB,45.33_EB,45.07_EB,44.81_EB,44.55_EB/)
QFLAME(5,16,21:30,4) = (/44.37_EB,44.19_EB,44.01_EB,43.83_EB,43.65_EB,43.43_EB,43.21_EB,42.99_EB,42.77_EB,42.55_EB/)
QFLAME(5,16,0:10,5) = (/0.00_EB,26.09_EB,37.07_EB,41.87_EB,46.67_EB,48.62_EB,50.57_EB,52.52_EB,54.47_EB,54.98_EB,55.48_EB/)
QFLAME(5,16,11:20,5) = (/55.99_EB,56.49_EB,56.64_EB,56.78_EB,56.93_EB,57.07_EB,56.61_EB,56.14_EB,55.67_EB,55.21_EB/)
QFLAME(5,16,21:30,5) = (/55.29_EB,55.37_EB,55.46_EB,55.54_EB,55.63_EB,55.30_EB,54.98_EB,54.66_EB,54.34_EB,54.01_EB/)
QFLAME(5,16,0:10,6) = (/0.00_EB,27.95_EB,40.53_EB,46.71_EB,52.89_EB,55.32_EB,57.74_EB,60.17_EB,62.59_EB,63.87_EB,65.14_EB/)
QFLAME(5,16,11:20,6) = (/66.41_EB,67.68_EB,67.94_EB,68.20_EB,68.45_EB,68.71_EB,68.37_EB,68.03_EB,67.69_EB,67.36_EB/)
QFLAME(5,16,21:30,6) = (/67.38_EB,67.40_EB,67.42_EB,67.44_EB,67.46_EB,67.37_EB,67.27_EB,67.18_EB,67.09_EB,67.00_EB/)
QFLAME(5,17,0:10,1) = (/0.00_EB,21.16_EB,29.21_EB,32.00_EB,34.80_EB,35.71_EB,36.62_EB,37.54_EB,38.45_EB,38.45_EB,38.45_EB/)
QFLAME(5,17,11:20,1) = (/38.45_EB,38.45_EB,38.13_EB,37.81_EB,37.50_EB,37.18_EB,36.90_EB,36.62_EB,36.34_EB,36.06_EB/)
QFLAME(5,17,21:30,1) = (/35.72_EB,35.38_EB,35.04_EB,34.70_EB,34.36_EB,34.09_EB,33.83_EB,33.56_EB,33.30_EB,33.03_EB/)
QFLAME(5,17,0:10,2) = (/0.00_EB,21.18_EB,29.21_EB,32.00_EB,34.79_EB,35.68_EB,36.58_EB,37.47_EB,38.37_EB,38.39_EB,38.42_EB/)
QFLAME(5,17,11:20,2) = (/38.44_EB,38.46_EB,38.14_EB,37.82_EB,37.50_EB,37.18_EB,36.90_EB,36.61_EB,36.33_EB,36.04_EB/)
QFLAME(5,17,21:30,2) = (/35.71_EB,35.37_EB,35.04_EB,34.70_EB,34.37_EB,34.11_EB,33.85_EB,33.60_EB,33.34_EB,33.08_EB/)
QFLAME(5,17,0:10,3) = (/0.00_EB,22.36_EB,30.23_EB,32.95_EB,35.67_EB,36.59_EB,37.50_EB,38.42_EB,39.34_EB,39.35_EB,39.36_EB/)
QFLAME(5,17,11:20,3) = (/39.38_EB,39.39_EB,39.04_EB,38.68_EB,38.33_EB,37.98_EB,37.73_EB,37.49_EB,37.24_EB,36.99_EB/)
QFLAME(5,17,21:30,3) = (/36.67_EB,36.34_EB,36.01_EB,35.69_EB,35.36_EB,35.09_EB,34.83_EB,34.56_EB,34.29_EB,34.02_EB/)
QFLAME(5,17,0:10,4) = (/0.00_EB,24.25_EB,33.42_EB,36.86_EB,40.30_EB,41.63_EB,42.95_EB,44.28_EB,45.61_EB,45.79_EB,45.98_EB/)
QFLAME(5,17,11:20,4) = (/46.17_EB,46.36_EB,46.16_EB,45.97_EB,45.78_EB,45.58_EB,45.33_EB,45.07_EB,44.82_EB,44.56_EB/)
QFLAME(5,17,21:30,4) = (/44.39_EB,44.21_EB,44.04_EB,43.86_EB,43.68_EB,43.47_EB,43.26_EB,43.05_EB,42.84_EB,42.63_EB/)
QFLAME(5,17,0:10,5) = (/0.00_EB,26.07_EB,37.03_EB,41.83_EB,46.63_EB,48.58_EB,50.53_EB,52.48_EB,54.43_EB,54.94_EB,55.45_EB/)
QFLAME(5,17,11:20,5) = (/55.96_EB,56.48_EB,56.62_EB,56.77_EB,56.91_EB,57.06_EB,56.60_EB,56.15_EB,55.69_EB,55.24_EB/)
QFLAME(5,17,21:30,5) = (/55.32_EB,55.41_EB,55.49_EB,55.58_EB,55.67_EB,55.34_EB,55.01_EB,54.68_EB,54.35_EB,54.02_EB/)
QFLAME(5,17,0:10,6) = (/0.00_EB,27.93_EB,40.49_EB,46.67_EB,52.85_EB,55.28_EB,57.70_EB,60.13_EB,62.55_EB,63.82_EB,65.09_EB/)
QFLAME(5,17,11:20,6) = (/66.36_EB,67.62_EB,67.86_EB,68.09_EB,68.33_EB,68.56_EB,68.26_EB,67.95_EB,67.64_EB,67.34_EB/)
QFLAME(5,17,21:30,6) = (/67.36_EB,67.38_EB,67.40_EB,67.43_EB,67.45_EB,67.37_EB,67.28_EB,67.20_EB,67.12_EB,67.04_EB/)
QFLAME(5,18,0:10,1) = (/0.00_EB,21.23_EB,29.37_EB,32.22_EB,35.08_EB,36.02_EB,36.97_EB,37.92_EB,38.86_EB,38.88_EB,38.90_EB/)
QFLAME(5,18,11:20,1) = (/38.92_EB,38.94_EB,38.63_EB,38.33_EB,38.02_EB,37.72_EB,37.44_EB,37.17_EB,36.90_EB,36.63_EB/)
QFLAME(5,18,21:30,1) = (/36.29_EB,35.95_EB,35.61_EB,35.27_EB,34.93_EB,34.67_EB,34.42_EB,34.16_EB,33.90_EB,33.64_EB/)
QFLAME(5,18,0:10,2) = (/0.00_EB,21.24_EB,29.37_EB,32.22_EB,35.07_EB,35.99_EB,36.92_EB,37.84_EB,38.77_EB,38.82_EB,38.86_EB/)
QFLAME(5,18,11:20,2) = (/38.91_EB,38.96_EB,38.65_EB,38.34_EB,38.03_EB,37.72_EB,37.44_EB,37.16_EB,36.88_EB,36.60_EB/)
QFLAME(5,18,21:30,2) = (/36.27_EB,35.93_EB,35.60_EB,35.27_EB,34.94_EB,34.69_EB,34.44_EB,34.19_EB,33.94_EB,33.69_EB/)
QFLAME(5,18,0:10,3) = (/0.00_EB,22.33_EB,30.22_EB,32.98_EB,35.74_EB,36.69_EB,37.65_EB,38.61_EB,39.57_EB,39.60_EB,39.64_EB/)
QFLAME(5,18,11:20,3) = (/39.68_EB,39.72_EB,39.36_EB,39.01_EB,38.65_EB,38.30_EB,38.05_EB,37.79_EB,37.54_EB,37.29_EB/)
QFLAME(5,18,21:30,3) = (/36.96_EB,36.64_EB,36.32_EB,35.99_EB,35.67_EB,35.42_EB,35.17_EB,34.91_EB,34.66_EB,34.41_EB/)
QFLAME(5,18,0:10,4) = (/0.00_EB,24.22_EB,33.38_EB,36.83_EB,40.28_EB,41.61_EB,42.95_EB,44.28_EB,45.62_EB,45.80_EB,45.98_EB/)
QFLAME(5,18,11:20,4) = (/46.17_EB,46.35_EB,46.16_EB,45.96_EB,45.77_EB,45.57_EB,45.32_EB,45.07_EB,44.83_EB,44.58_EB/)
QFLAME(5,18,21:30,4) = (/44.40_EB,44.23_EB,44.06_EB,43.89_EB,43.72_EB,43.51_EB,43.31_EB,43.11_EB,42.91_EB,42.71_EB/)
QFLAME(5,18,0:10,5) = (/0.00_EB,26.05_EB,36.99_EB,41.79_EB,46.59_EB,48.54_EB,50.49_EB,52.45_EB,54.40_EB,54.91_EB,55.43_EB/)
QFLAME(5,18,11:20,5) = (/55.94_EB,56.46_EB,56.60_EB,56.75_EB,56.90_EB,57.04_EB,56.60_EB,56.15_EB,55.71_EB,55.26_EB/)
QFLAME(5,18,21:30,5) = (/55.35_EB,55.44_EB,55.53_EB,55.62_EB,55.70_EB,55.37_EB,55.04_EB,54.70_EB,54.37_EB,54.04_EB/)
QFLAME(5,18,0:10,6) = (/0.00_EB,27.92_EB,40.46_EB,46.63_EB,52.81_EB,55.24_EB,57.66_EB,60.09_EB,62.51_EB,63.78_EB,65.04_EB/)
QFLAME(5,18,11:20,6) = (/66.30_EB,67.57_EB,67.78_EB,67.99_EB,68.20_EB,68.42_EB,68.14_EB,67.87_EB,67.59_EB,67.32_EB/)
QFLAME(5,18,21:30,6) = (/67.34_EB,67.37_EB,67.39_EB,67.42_EB,67.44_EB,67.37_EB,67.29_EB,67.22_EB,67.15_EB,67.07_EB/)
QFLAME(5,19,0:10,1) = (/0.00_EB,21.30_EB,29.53_EB,32.44_EB,35.36_EB,36.33_EB,37.31_EB,38.29_EB,39.27_EB,39.31_EB,39.35_EB/)
QFLAME(5,19,11:20,1) = (/39.39_EB,39.43_EB,39.14_EB,38.84_EB,38.55_EB,38.25_EB,37.99_EB,37.72_EB,37.46_EB,37.19_EB/)
QFLAME(5,19,21:30,1) = (/36.85_EB,36.52_EB,36.18_EB,35.84_EB,35.51_EB,35.25_EB,35.00_EB,34.75_EB,34.50_EB,34.25_EB/)
QFLAME(5,19,0:10,2) = (/0.00_EB,21.30_EB,29.53_EB,32.44_EB,35.35_EB,36.30_EB,37.26_EB,38.21_EB,39.16_EB,39.24_EB,39.31_EB/)
QFLAME(5,19,11:20,2) = (/39.38_EB,39.46_EB,39.16_EB,38.86_EB,38.56_EB,38.26_EB,37.98_EB,37.71_EB,37.43_EB,37.15_EB/)
QFLAME(5,19,21:30,2) = (/36.82_EB,36.50_EB,36.17_EB,35.84_EB,35.51_EB,35.27_EB,35.02_EB,34.78_EB,34.53_EB,34.29_EB/)
QFLAME(5,19,0:10,3) = (/0.00_EB,22.30_EB,30.21_EB,33.00_EB,35.80_EB,36.80_EB,37.80_EB,38.80_EB,39.80_EB,39.86_EB,39.92_EB/)
QFLAME(5,19,11:20,3) = (/39.98_EB,40.05_EB,39.69_EB,39.33_EB,38.97_EB,38.62_EB,38.36_EB,38.10_EB,37.84_EB,37.59_EB/)
QFLAME(5,19,21:30,3) = (/37.26_EB,36.94_EB,36.62_EB,36.30_EB,35.98_EB,35.74_EB,35.51_EB,35.27_EB,35.04_EB,34.80_EB/)
QFLAME(5,19,0:10,4) = (/0.00_EB,24.20_EB,33.34_EB,36.80_EB,40.25_EB,41.60_EB,42.94_EB,44.29_EB,45.63_EB,45.81_EB,45.99_EB/)
QFLAME(5,19,11:20,4) = (/46.16_EB,46.34_EB,46.15_EB,45.95_EB,45.76_EB,45.56_EB,45.32_EB,45.08_EB,44.83_EB,44.59_EB/)
QFLAME(5,19,21:30,4) = (/44.42_EB,44.25_EB,44.08_EB,43.92_EB,43.75_EB,43.55_EB,43.36_EB,43.17_EB,42.98_EB,42.78_EB/)
QFLAME(5,19,0:10,5) = (/0.00_EB,26.02_EB,36.96_EB,41.76_EB,46.55_EB,48.51_EB,50.46_EB,52.41_EB,54.36_EB,54.88_EB,55.40_EB/)
QFLAME(5,19,11:20,5) = (/55.92_EB,56.44_EB,56.58_EB,56.73_EB,56.88_EB,57.03_EB,56.60_EB,56.16_EB,55.73_EB,55.29_EB/)
QFLAME(5,19,21:30,5) = (/55.38_EB,55.47_EB,55.56_EB,55.65_EB,55.74_EB,55.40_EB,55.06_EB,54.73_EB,54.39_EB,54.05_EB/)
QFLAME(5,19,0:10,6) = (/0.00_EB,27.90_EB,40.42_EB,46.60_EB,52.77_EB,55.20_EB,57.62_EB,60.05_EB,62.47_EB,63.73_EB,64.99_EB/)
QFLAME(5,19,11:20,6) = (/66.25_EB,67.51_EB,67.70_EB,67.89_EB,68.08_EB,68.27_EB,68.03_EB,67.78_EB,67.54_EB,67.30_EB/)
QFLAME(5,19,21:30,6) = (/67.32_EB,67.35_EB,67.38_EB,67.41_EB,67.43_EB,67.37_EB,67.30_EB,67.24_EB,67.17_EB,67.11_EB/)
QFLAME(5,20,0:10,1) = (/0.00_EB,21.37_EB,29.69_EB,32.66_EB,35.63_EB,36.65_EB,37.66_EB,38.67_EB,39.69_EB,39.74_EB,39.80_EB/)
QFLAME(5,20,11:20,1) = (/39.86_EB,39.92_EB,39.64_EB,39.36_EB,39.07_EB,38.79_EB,38.53_EB,38.27_EB,38.01_EB,37.75_EB/)
QFLAME(5,20,21:30,1) = (/37.42_EB,37.08_EB,36.75_EB,36.41_EB,36.08_EB,35.83_EB,35.59_EB,35.34_EB,35.10_EB,34.85_EB/)
QFLAME(5,20,0:10,2) = (/0.00_EB,21.36_EB,29.69_EB,32.66_EB,35.63_EB,36.61_EB,37.59_EB,38.58_EB,39.56_EB,39.66_EB,39.76_EB/)
QFLAME(5,20,11:20,2) = (/39.85_EB,39.95_EB,39.66_EB,39.38_EB,39.09_EB,38.80_EB,38.53_EB,38.26_EB,37.98_EB,37.71_EB/)
QFLAME(5,20,21:30,2) = (/37.38_EB,37.06_EB,36.74_EB,36.41_EB,36.09_EB,35.85_EB,35.61_EB,35.37_EB,35.13_EB,34.89_EB/)
QFLAME(5,20,0:10,3) = (/0.00_EB,22.28_EB,30.20_EB,33.03_EB,35.86_EB,36.90_EB,37.94_EB,38.98_EB,40.02_EB,40.11_EB,40.20_EB/)
QFLAME(5,20,11:20,3) = (/40.29_EB,40.38_EB,40.02_EB,39.65_EB,39.29_EB,38.93_EB,38.67_EB,38.41_EB,38.15_EB,37.88_EB/)
QFLAME(5,20,21:30,3) = (/37.56_EB,37.24_EB,36.92_EB,36.60_EB,36.28_EB,36.06_EB,35.85_EB,35.63_EB,35.41_EB,35.20_EB/)
QFLAME(5,20,0:10,4) = (/0.00_EB,24.17_EB,33.31_EB,36.77_EB,40.22_EB,41.58_EB,42.94_EB,44.29_EB,45.65_EB,45.82_EB,45.99_EB/)
QFLAME(5,20,11:20,4) = (/46.16_EB,46.33_EB,46.14_EB,45.94_EB,45.75_EB,45.55_EB,45.32_EB,45.08_EB,44.84_EB,44.60_EB/)
QFLAME(5,20,21:30,4) = (/44.44_EB,44.27_EB,44.11_EB,43.94_EB,43.78_EB,43.60_EB,43.41_EB,43.23_EB,43.04_EB,42.86_EB/)
QFLAME(5,20,0:10,5) = (/0.00_EB,26.00_EB,36.92_EB,41.72_EB,46.52_EB,48.47_EB,50.42_EB,52.37_EB,54.32_EB,54.84_EB,55.37_EB/)
QFLAME(5,20,11:20,5) = (/55.89_EB,56.42_EB,56.57_EB,56.72_EB,56.86_EB,57.01_EB,56.59_EB,56.17_EB,55.75_EB,55.32_EB/)
QFLAME(5,20,21:30,5) = (/55.41_EB,55.51_EB,55.60_EB,55.69_EB,55.78_EB,55.44_EB,55.09_EB,54.75_EB,54.40_EB,54.06_EB/)
QFLAME(5,20,0:10,6) = (/0.00_EB,27.88_EB,40.38_EB,46.56_EB,52.74_EB,55.16_EB,57.58_EB,60.01_EB,62.43_EB,63.69_EB,64.94_EB/)
QFLAME(5,20,11:20,6) = (/66.19_EB,67.45_EB,67.62_EB,67.79_EB,67.96_EB,68.12_EB,67.91_EB,67.70_EB,67.49_EB,67.28_EB/)
QFLAME(5,20,21:30,6) = (/67.31_EB,67.34_EB,67.37_EB,67.40_EB,67.43_EB,67.37_EB,67.31_EB,67.26_EB,67.20_EB,67.14_EB/)
ABSF(1,0,0:10,1) = (/0.000_EB,0.013_EB,0.025_EB,0.025_EB,0.025_EB,0.029_EB,0.033_EB,0.037_EB,0.041_EB,0.042_EB,0.043_EB/)
ABSF(1,0,11:20,1) = (/0.044_EB,0.046_EB,0.046_EB,0.047_EB,0.048_EB,0.048_EB,0.050_EB,0.051_EB,0.052_EB,0.054_EB/)
ABSF(1,0,21:30,1) = (/0.054_EB,0.054_EB,0.054_EB,0.054_EB,0.054_EB,0.057_EB,0.059_EB,0.062_EB,0.064_EB,0.067_EB/)
ABSF(1,0,0:10,2) = (/0.000_EB,0.008_EB,0.021_EB,0.028_EB,0.034_EB,0.038_EB,0.043_EB,0.047_EB,0.051_EB,0.055_EB,0.059_EB/)
ABSF(1,0,11:20,2) = (/0.063_EB,0.067_EB,0.064_EB,0.061_EB,0.059_EB,0.056_EB,0.055_EB,0.053_EB,0.052_EB,0.051_EB/)
ABSF(1,0,21:30,2) = (/0.052_EB,0.054_EB,0.056_EB,0.057_EB,0.059_EB,0.066_EB,0.073_EB,0.080_EB,0.087_EB,0.094_EB/)
ABSF(1,0,0:10,3) = (/0.000_EB,0.006_EB,0.021_EB,0.031_EB,0.041_EB,0.045_EB,0.049_EB,0.053_EB,0.057_EB,0.065_EB,0.074_EB/)
ABSF(1,0,11:20,3) = (/0.082_EB,0.091_EB,0.084_EB,0.078_EB,0.072_EB,0.066_EB,0.064_EB,0.063_EB,0.061_EB,0.060_EB/)
ABSF(1,0,21:30,3) = (/0.061_EB,0.062_EB,0.064_EB,0.065_EB,0.066_EB,0.078_EB,0.089_EB,0.101_EB,0.112_EB,0.124_EB/)
ABSF(1,0,0:10,4) = (/0.000_EB,0.010_EB,0.023_EB,0.035_EB,0.047_EB,0.047_EB,0.046_EB,0.046_EB,0.045_EB,0.060_EB,0.074_EB/)
ABSF(1,0,11:20,4) = (/0.089_EB,0.104_EB,0.103_EB,0.103_EB,0.103_EB,0.103_EB,0.095_EB,0.088_EB,0.081_EB,0.074_EB/)
ABSF(1,0,21:30,4) = (/0.076_EB,0.077_EB,0.079_EB,0.081_EB,0.083_EB,0.099_EB,0.116_EB,0.132_EB,0.148_EB,0.164_EB/)
ABSF(1,0,0:10,5) = (/0.000_EB,0.010_EB,0.028_EB,0.043_EB,0.059_EB,0.061_EB,0.064_EB,0.066_EB,0.069_EB,0.088_EB,0.107_EB/)
ABSF(1,0,11:20,5) = (/0.126_EB,0.145_EB,0.145_EB,0.146_EB,0.146_EB,0.146_EB,0.128_EB,0.110_EB,0.092_EB,0.074_EB/)
ABSF(1,0,21:30,5) = (/0.077_EB,0.080_EB,0.083_EB,0.086_EB,0.089_EB,0.110_EB,0.131_EB,0.152_EB,0.174_EB,0.195_EB/)
ABSF(1,0,0:10,6) = (/0.000_EB,0.011_EB,0.029_EB,0.037_EB,0.045_EB,0.048_EB,0.052_EB,0.056_EB,0.059_EB,0.084_EB,0.108_EB/)
ABSF(1,0,11:20,6) = (/0.133_EB,0.157_EB,0.171_EB,0.184_EB,0.198_EB,0.212_EB,0.197_EB,0.181_EB,0.166_EB,0.151_EB/)
ABSF(1,0,21:30,6) = (/0.143_EB,0.136_EB,0.128_EB,0.120_EB,0.113_EB,0.138_EB,0.163_EB,0.188_EB,0.213_EB,0.239_EB/)
ABSF(1,1,0:10,1) = (/0.000_EB,0.027_EB,0.030_EB,0.031_EB,0.032_EB,0.035_EB,0.039_EB,0.042_EB,0.046_EB,0.047_EB,0.048_EB/)
ABSF(1,1,11:20,1) = (/0.050_EB,0.051_EB,0.051_EB,0.051_EB,0.052_EB,0.052_EB,0.053_EB,0.054_EB,0.056_EB,0.057_EB/)
ABSF(1,1,21:30,1) = (/0.058_EB,0.058_EB,0.058_EB,0.059_EB,0.059_EB,0.062_EB,0.064_EB,0.067_EB,0.069_EB,0.071_EB/)
ABSF(1,1,0:10,2) = (/0.000_EB,0.012_EB,0.024_EB,0.032_EB,0.040_EB,0.044_EB,0.047_EB,0.050_EB,0.053_EB,0.059_EB,0.064_EB/)
ABSF(1,1,11:20,2) = (/0.070_EB,0.075_EB,0.073_EB,0.072_EB,0.070_EB,0.068_EB,0.066_EB,0.064_EB,0.062_EB,0.060_EB/)
ABSF(1,1,21:30,2) = (/0.061_EB,0.063_EB,0.064_EB,0.066_EB,0.067_EB,0.072_EB,0.078_EB,0.083_EB,0.088_EB,0.093_EB/)
ABSF(1,1,0:10,3) = (/0.000_EB,0.015_EB,0.026_EB,0.038_EB,0.049_EB,0.051_EB,0.054_EB,0.056_EB,0.059_EB,0.067_EB,0.076_EB/)
ABSF(1,1,11:20,3) = (/0.084_EB,0.092_EB,0.088_EB,0.083_EB,0.079_EB,0.074_EB,0.072_EB,0.071_EB,0.069_EB,0.067_EB/)
ABSF(1,1,21:30,3) = (/0.070_EB,0.072_EB,0.075_EB,0.077_EB,0.079_EB,0.089_EB,0.099_EB,0.109_EB,0.118_EB,0.128_EB/)
ABSF(1,1,0:10,4) = (/0.000_EB,0.012_EB,0.028_EB,0.042_EB,0.057_EB,0.057_EB,0.058_EB,0.058_EB,0.059_EB,0.072_EB,0.085_EB/)
ABSF(1,1,11:20,4) = (/0.098_EB,0.111_EB,0.107_EB,0.103_EB,0.099_EB,0.095_EB,0.094_EB,0.092_EB,0.091_EB,0.089_EB/)
ABSF(1,1,21:30,4) = (/0.091_EB,0.094_EB,0.096_EB,0.099_EB,0.101_EB,0.112_EB,0.122_EB,0.133_EB,0.143_EB,0.153_EB/)
ABSF(1,1,0:10,5) = (/0.000_EB,0.014_EB,0.032_EB,0.043_EB,0.053_EB,0.057_EB,0.060_EB,0.063_EB,0.067_EB,0.086_EB,0.105_EB/)
ABSF(1,1,11:20,5) = (/0.124_EB,0.143_EB,0.147_EB,0.150_EB,0.154_EB,0.157_EB,0.143_EB,0.129_EB,0.115_EB,0.101_EB/)
ABSF(1,1,21:30,5) = (/0.100_EB,0.098_EB,0.097_EB,0.095_EB,0.094_EB,0.113_EB,0.132_EB,0.151_EB,0.170_EB,0.189_EB/)
ABSF(1,1,0:10,6) = (/0.000_EB,0.013_EB,0.033_EB,0.044_EB,0.054_EB,0.059_EB,0.065_EB,0.070_EB,0.075_EB,0.091_EB,0.108_EB/)
ABSF(1,1,11:20,6) = (/0.124_EB,0.140_EB,0.152_EB,0.164_EB,0.177_EB,0.189_EB,0.177_EB,0.165_EB,0.153_EB,0.141_EB/)
ABSF(1,1,21:30,6) = (/0.132_EB,0.122_EB,0.113_EB,0.103_EB,0.093_EB,0.113_EB,0.133_EB,0.154_EB,0.174_EB,0.194_EB/)
ABSF(1,2,0:10,1) = (/0.000_EB,0.034_EB,0.038_EB,0.038_EB,0.038_EB,0.041_EB,0.044_EB,0.047_EB,0.050_EB,0.052_EB,0.053_EB/)
ABSF(1,2,11:20,1) = (/0.055_EB,0.057_EB,0.056_EB,0.056_EB,0.056_EB,0.056_EB,0.057_EB,0.058_EB,0.060_EB,0.061_EB/)
ABSF(1,2,21:30,1) = (/0.062_EB,0.063_EB,0.065_EB,0.066_EB,0.067_EB,0.069_EB,0.071_EB,0.072_EB,0.074_EB,0.076_EB/)
ABSF(1,2,0:10,2) = (/0.000_EB,0.015_EB,0.033_EB,0.036_EB,0.039_EB,0.042_EB,0.046_EB,0.049_EB,0.053_EB,0.055_EB,0.057_EB/)
ABSF(1,2,11:20,2) = (/0.059_EB,0.061_EB,0.061_EB,0.061_EB,0.062_EB,0.062_EB,0.062_EB,0.062_EB,0.063_EB,0.063_EB/)
ABSF(1,2,21:30,2) = (/0.063_EB,0.063_EB,0.064_EB,0.064_EB,0.064_EB,0.071_EB,0.078_EB,0.085_EB,0.091_EB,0.098_EB/)
ABSF(1,2,0:10,3) = (/0.000_EB,0.018_EB,0.032_EB,0.043_EB,0.054_EB,0.056_EB,0.057_EB,0.058_EB,0.059_EB,0.069_EB,0.078_EB/)
ABSF(1,2,11:20,3) = (/0.087_EB,0.096_EB,0.092_EB,0.088_EB,0.085_EB,0.081_EB,0.080_EB,0.079_EB,0.078_EB,0.078_EB/)
ABSF(1,2,21:30,3) = (/0.079_EB,0.079_EB,0.080_EB,0.081_EB,0.082_EB,0.092_EB,0.101_EB,0.110_EB,0.119_EB,0.129_EB/)
ABSF(1,2,0:10,4) = (/0.000_EB,0.015_EB,0.032_EB,0.045_EB,0.059_EB,0.062_EB,0.066_EB,0.069_EB,0.073_EB,0.084_EB,0.094_EB/)
ABSF(1,2,11:20,4) = (/0.105_EB,0.116_EB,0.111_EB,0.106_EB,0.102_EB,0.097_EB,0.093_EB,0.089_EB,0.086_EB,0.082_EB/)
ABSF(1,2,21:30,4) = (/0.084_EB,0.087_EB,0.089_EB,0.092_EB,0.094_EB,0.103_EB,0.113_EB,0.122_EB,0.132_EB,0.141_EB/)
ABSF(1,2,0:10,5) = (/0.000_EB,0.020_EB,0.034_EB,0.044_EB,0.053_EB,0.057_EB,0.060_EB,0.064_EB,0.067_EB,0.081_EB,0.096_EB/)
ABSF(1,2,11:20,5) = (/0.110_EB,0.124_EB,0.124_EB,0.124_EB,0.124_EB,0.124_EB,0.116_EB,0.108_EB,0.100_EB,0.092_EB/)
ABSF(1,2,21:30,5) = (/0.096_EB,0.099_EB,0.103_EB,0.107_EB,0.111_EB,0.130_EB,0.148_EB,0.167_EB,0.185_EB,0.204_EB/)
ABSF(1,2,0:10,6) = (/0.000_EB,0.022_EB,0.037_EB,0.039_EB,0.041_EB,0.049_EB,0.057_EB,0.065_EB,0.073_EB,0.092_EB,0.112_EB/)
ABSF(1,2,11:20,6) = (/0.131_EB,0.151_EB,0.161_EB,0.171_EB,0.181_EB,0.191_EB,0.173_EB,0.155_EB,0.137_EB,0.119_EB/)
ABSF(1,2,21:30,6) = (/0.111_EB,0.102_EB,0.093_EB,0.084_EB,0.075_EB,0.106_EB,0.137_EB,0.168_EB,0.199_EB,0.229_EB/)
ABSF(1,3,0:10,1) = (/0.000_EB,0.038_EB,0.044_EB,0.044_EB,0.044_EB,0.047_EB,0.050_EB,0.053_EB,0.056_EB,0.057_EB,0.059_EB/)
ABSF(1,3,11:20,1) = (/0.060_EB,0.062_EB,0.061_EB,0.061_EB,0.061_EB,0.060_EB,0.061_EB,0.062_EB,0.064_EB,0.065_EB/)
ABSF(1,3,21:30,1) = (/0.067_EB,0.069_EB,0.071_EB,0.073_EB,0.075_EB,0.076_EB,0.077_EB,0.078_EB,0.079_EB,0.080_EB/)
ABSF(1,3,0:10,2) = (/0.000_EB,0.024_EB,0.040_EB,0.042_EB,0.045_EB,0.048_EB,0.051_EB,0.054_EB,0.058_EB,0.060_EB,0.061_EB/)
ABSF(1,3,11:20,2) = (/0.063_EB,0.065_EB,0.065_EB,0.065_EB,0.064_EB,0.064_EB,0.065_EB,0.065_EB,0.066_EB,0.067_EB/)
ABSF(1,3,21:30,2) = (/0.068_EB,0.069_EB,0.071_EB,0.072_EB,0.073_EB,0.078_EB,0.082_EB,0.087_EB,0.091_EB,0.096_EB/)
ABSF(1,3,0:10,3) = (/0.000_EB,0.021_EB,0.036_EB,0.046_EB,0.056_EB,0.058_EB,0.061_EB,0.063_EB,0.065_EB,0.074_EB,0.083_EB/)
ABSF(1,3,11:20,3) = (/0.092_EB,0.101_EB,0.097_EB,0.094_EB,0.090_EB,0.086_EB,0.085_EB,0.084_EB,0.083_EB,0.082_EB/)
ABSF(1,3,21:30,3) = (/0.083_EB,0.084_EB,0.085_EB,0.086_EB,0.087_EB,0.095_EB,0.103_EB,0.111_EB,0.120_EB,0.128_EB/)
ABSF(1,3,0:10,4) = (/0.000_EB,0.019_EB,0.036_EB,0.050_EB,0.064_EB,0.067_EB,0.069_EB,0.072_EB,0.075_EB,0.087_EB,0.099_EB/)
ABSF(1,3,11:20,4) = (/0.111_EB,0.122_EB,0.117_EB,0.111_EB,0.106_EB,0.100_EB,0.098_EB,0.096_EB,0.094_EB,0.091_EB/)
ABSF(1,3,21:30,4) = (/0.094_EB,0.096_EB,0.098_EB,0.100_EB,0.103_EB,0.111_EB,0.120_EB,0.129_EB,0.138_EB,0.146_EB/)
ABSF(1,3,0:10,5) = (/0.000_EB,0.025_EB,0.039_EB,0.049_EB,0.060_EB,0.062_EB,0.065_EB,0.067_EB,0.070_EB,0.085_EB,0.101_EB/)
ABSF(1,3,11:20,5) = (/0.117_EB,0.133_EB,0.132_EB,0.132_EB,0.132_EB,0.131_EB,0.123_EB,0.114_EB,0.105_EB,0.096_EB/)
ABSF(1,3,21:30,5) = (/0.100_EB,0.103_EB,0.107_EB,0.111_EB,0.114_EB,0.131_EB,0.149_EB,0.166_EB,0.184_EB,0.201_EB/)
ABSF(1,3,0:10,6) = (/0.000_EB,0.026_EB,0.041_EB,0.046_EB,0.050_EB,0.055_EB,0.061_EB,0.067_EB,0.073_EB,0.093_EB,0.113_EB/)
ABSF(1,3,11:20,6) = (/0.133_EB,0.153_EB,0.162_EB,0.172_EB,0.182_EB,0.192_EB,0.174_EB,0.156_EB,0.139_EB,0.121_EB/)
ABSF(1,3,21:30,6) = (/0.113_EB,0.105_EB,0.097_EB,0.089_EB,0.082_EB,0.111_EB,0.141_EB,0.171_EB,0.201_EB,0.231_EB/)
ABSF(1,4,0:10,1) = (/0.000_EB,0.043_EB,0.049_EB,0.049_EB,0.049_EB,0.052_EB,0.055_EB,0.058_EB,0.061_EB,0.062_EB,0.064_EB/)
ABSF(1,4,11:20,1) = (/0.065_EB,0.067_EB,0.066_EB,0.066_EB,0.065_EB,0.065_EB,0.066_EB,0.067_EB,0.067_EB,0.068_EB/)
ABSF(1,4,21:30,1) = (/0.071_EB,0.074_EB,0.077_EB,0.080_EB,0.083_EB,0.083_EB,0.084_EB,0.084_EB,0.084_EB,0.084_EB/)
ABSF(1,4,0:10,2) = (/0.000_EB,0.033_EB,0.048_EB,0.049_EB,0.050_EB,0.053_EB,0.056_EB,0.060_EB,0.063_EB,0.064_EB,0.066_EB/)
ABSF(1,4,11:20,2) = (/0.067_EB,0.069_EB,0.068_EB,0.068_EB,0.067_EB,0.067_EB,0.068_EB,0.069_EB,0.070_EB,0.071_EB/)
ABSF(1,4,21:30,2) = (/0.073_EB,0.075_EB,0.078_EB,0.080_EB,0.082_EB,0.084_EB,0.086_EB,0.089_EB,0.091_EB,0.093_EB/)
ABSF(1,4,0:10,3) = (/0.000_EB,0.024_EB,0.041_EB,0.049_EB,0.058_EB,0.061_EB,0.064_EB,0.067_EB,0.070_EB,0.079_EB,0.088_EB/)
ABSF(1,4,11:20,3) = (/0.097_EB,0.106_EB,0.102_EB,0.099_EB,0.095_EB,0.092_EB,0.090_EB,0.089_EB,0.087_EB,0.085_EB/)
ABSF(1,4,21:30,3) = (/0.087_EB,0.088_EB,0.089_EB,0.090_EB,0.091_EB,0.099_EB,0.106_EB,0.113_EB,0.120_EB,0.127_EB/)
ABSF(1,4,0:10,4) = (/0.000_EB,0.022_EB,0.041_EB,0.055_EB,0.069_EB,0.071_EB,0.073_EB,0.076_EB,0.078_EB,0.090_EB,0.103_EB/)
ABSF(1,4,11:20,4) = (/0.116_EB,0.129_EB,0.122_EB,0.116_EB,0.109_EB,0.103_EB,0.102_EB,0.102_EB,0.101_EB,0.101_EB/)
ABSF(1,4,21:30,4) = (/0.103_EB,0.105_EB,0.107_EB,0.109_EB,0.111_EB,0.119_EB,0.127_EB,0.136_EB,0.144_EB,0.152_EB/)
ABSF(1,4,0:10,5) = (/0.000_EB,0.029_EB,0.043_EB,0.054_EB,0.066_EB,0.068_EB,0.069_EB,0.071_EB,0.072_EB,0.089_EB,0.107_EB/)
ABSF(1,4,11:20,5) = (/0.124_EB,0.141_EB,0.141_EB,0.140_EB,0.139_EB,0.138_EB,0.129_EB,0.120_EB,0.110_EB,0.101_EB/)
ABSF(1,4,21:30,5) = (/0.104_EB,0.107_EB,0.111_EB,0.114_EB,0.117_EB,0.133_EB,0.149_EB,0.165_EB,0.182_EB,0.198_EB/)
ABSF(1,4,0:10,6) = (/0.000_EB,0.030_EB,0.046_EB,0.052_EB,0.059_EB,0.062_EB,0.066_EB,0.069_EB,0.073_EB,0.093_EB,0.113_EB/)
ABSF(1,4,11:20,6) = (/0.134_EB,0.154_EB,0.164_EB,0.174_EB,0.183_EB,0.193_EB,0.175_EB,0.158_EB,0.140_EB,0.122_EB/)
ABSF(1,4,21:30,6) = (/0.116_EB,0.109_EB,0.102_EB,0.095_EB,0.089_EB,0.117_EB,0.146_EB,0.175_EB,0.204_EB,0.232_EB/)
ABSF(1,5,0:10,1) = (/0.000_EB,0.047_EB,0.055_EB,0.055_EB,0.054_EB,0.057_EB,0.060_EB,0.063_EB,0.066_EB,0.068_EB,0.069_EB/)
ABSF(1,5,11:20,1) = (/0.071_EB,0.072_EB,0.071_EB,0.071_EB,0.070_EB,0.069_EB,0.070_EB,0.071_EB,0.071_EB,0.072_EB/)
ABSF(1,5,21:30,1) = (/0.076_EB,0.080_EB,0.084_EB,0.087_EB,0.091_EB,0.091_EB,0.090_EB,0.089_EB,0.089_EB,0.088_EB/)
ABSF(1,5,0:10,2) = (/0.000_EB,0.042_EB,0.056_EB,0.056_EB,0.056_EB,0.059_EB,0.062_EB,0.065_EB,0.068_EB,0.069_EB,0.070_EB/)
ABSF(1,5,11:20,2) = (/0.072_EB,0.073_EB,0.072_EB,0.071_EB,0.070_EB,0.069_EB,0.070_EB,0.072_EB,0.073_EB,0.074_EB/)
ABSF(1,5,21:30,2) = (/0.078_EB,0.081_EB,0.084_EB,0.088_EB,0.091_EB,0.091_EB,0.091_EB,0.091_EB,0.090_EB,0.090_EB/)
ABSF(1,5,0:10,3) = (/0.000_EB,0.027_EB,0.045_EB,0.053_EB,0.060_EB,0.064_EB,0.068_EB,0.072_EB,0.075_EB,0.084_EB,0.093_EB/)
ABSF(1,5,11:20,3) = (/0.102_EB,0.110_EB,0.107_EB,0.104_EB,0.101_EB,0.098_EB,0.095_EB,0.093_EB,0.091_EB,0.089_EB/)
ABSF(1,5,21:30,3) = (/0.091_EB,0.092_EB,0.093_EB,0.095_EB,0.096_EB,0.102_EB,0.108_EB,0.114_EB,0.120_EB,0.126_EB/)
ABSF(1,5,0:10,4) = (/0.000_EB,0.026_EB,0.045_EB,0.059_EB,0.074_EB,0.075_EB,0.077_EB,0.079_EB,0.080_EB,0.094_EB,0.108_EB/)
ABSF(1,5,11:20,4) = (/0.121_EB,0.135_EB,0.128_EB,0.120_EB,0.113_EB,0.106_EB,0.107_EB,0.108_EB,0.109_EB,0.110_EB/)
ABSF(1,5,21:30,4) = (/0.112_EB,0.114_EB,0.116_EB,0.118_EB,0.120_EB,0.127_EB,0.135_EB,0.142_EB,0.150_EB,0.157_EB/)
ABSF(1,5,0:10,5) = (/0.000_EB,0.033_EB,0.047_EB,0.060_EB,0.072_EB,0.073_EB,0.074_EB,0.074_EB,0.075_EB,0.093_EB,0.112_EB/)
ABSF(1,5,11:20,5) = (/0.131_EB,0.150_EB,0.149_EB,0.147_EB,0.146_EB,0.145_EB,0.135_EB,0.125_EB,0.116_EB,0.106_EB/)
ABSF(1,5,21:30,5) = (/0.109_EB,0.112_EB,0.114_EB,0.117_EB,0.120_EB,0.135_EB,0.150_EB,0.165_EB,0.180_EB,0.194_EB/)
ABSF(1,5,0:10,6) = (/0.000_EB,0.034_EB,0.050_EB,0.059_EB,0.067_EB,0.069_EB,0.070_EB,0.071_EB,0.072_EB,0.093_EB,0.114_EB/)
ABSF(1,5,11:20,6) = (/0.135_EB,0.156_EB,0.166_EB,0.175_EB,0.184_EB,0.194_EB,0.176_EB,0.159_EB,0.141_EB,0.124_EB/)
ABSF(1,5,21:30,6) = (/0.118_EB,0.112_EB,0.107_EB,0.101_EB,0.096_EB,0.123_EB,0.151_EB,0.179_EB,0.206_EB,0.234_EB/)
ABSF(1,6,0:10,1) = (/0.000_EB,0.052_EB,0.061_EB,0.060_EB,0.059_EB,0.062_EB,0.065_EB,0.068_EB,0.071_EB,0.073_EB,0.074_EB/)
ABSF(1,6,11:20,1) = (/0.076_EB,0.077_EB,0.076_EB,0.076_EB,0.075_EB,0.074_EB,0.075_EB,0.075_EB,0.076_EB,0.076_EB/)
ABSF(1,6,21:30,1) = (/0.080_EB,0.085_EB,0.089_EB,0.093_EB,0.097_EB,0.096_EB,0.095_EB,0.094_EB,0.093_EB,0.092_EB/)
ABSF(1,6,0:10,2) = (/0.000_EB,0.048_EB,0.060_EB,0.061_EB,0.061_EB,0.064_EB,0.067_EB,0.070_EB,0.072_EB,0.074_EB,0.075_EB/)
ABSF(1,6,11:20,2) = (/0.077_EB,0.078_EB,0.077_EB,0.076_EB,0.075_EB,0.074_EB,0.075_EB,0.076_EB,0.077_EB,0.079_EB/)
ABSF(1,6,21:30,2) = (/0.082_EB,0.086_EB,0.090_EB,0.094_EB,0.098_EB,0.097_EB,0.096_EB,0.096_EB,0.095_EB,0.094_EB/)
ABSF(1,6,0:10,3) = (/0.000_EB,0.030_EB,0.052_EB,0.058_EB,0.065_EB,0.068_EB,0.072_EB,0.075_EB,0.078_EB,0.086_EB,0.093_EB/)
ABSF(1,6,11:20,3) = (/0.101_EB,0.108_EB,0.105_EB,0.102_EB,0.100_EB,0.097_EB,0.095_EB,0.094_EB,0.092_EB,0.091_EB/)
ABSF(1,6,21:30,3) = (/0.093_EB,0.095_EB,0.097_EB,0.100_EB,0.102_EB,0.107_EB,0.112_EB,0.117_EB,0.122_EB,0.128_EB/)
ABSF(1,6,0:10,4) = (/0.000_EB,0.029_EB,0.049_EB,0.063_EB,0.077_EB,0.079_EB,0.081_EB,0.083_EB,0.085_EB,0.099_EB,0.112_EB/)
ABSF(1,6,11:20,4) = (/0.125_EB,0.138_EB,0.132_EB,0.125_EB,0.119_EB,0.112_EB,0.112_EB,0.113_EB,0.113_EB,0.113_EB/)
ABSF(1,6,21:30,4) = (/0.116_EB,0.118_EB,0.120_EB,0.122_EB,0.125_EB,0.132_EB,0.139_EB,0.146_EB,0.153_EB,0.160_EB/)
ABSF(1,6,0:10,5) = (/0.000_EB,0.036_EB,0.051_EB,0.065_EB,0.078_EB,0.078_EB,0.078_EB,0.077_EB,0.077_EB,0.096_EB,0.115_EB/)
ABSF(1,6,11:20,5) = (/0.134_EB,0.153_EB,0.152_EB,0.152_EB,0.151_EB,0.150_EB,0.141_EB,0.132_EB,0.123_EB,0.114_EB/)
ABSF(1,6,21:30,5) = (/0.116_EB,0.118_EB,0.120_EB,0.121_EB,0.123_EB,0.138_EB,0.153_EB,0.168_EB,0.182_EB,0.197_EB/)
ABSF(1,6,0:10,6) = (/0.000_EB,0.036_EB,0.055_EB,0.062_EB,0.070_EB,0.073_EB,0.076_EB,0.080_EB,0.083_EB,0.102_EB,0.121_EB/)
ABSF(1,6,11:20,6) = (/0.140_EB,0.159_EB,0.169_EB,0.178_EB,0.188_EB,0.198_EB,0.183_EB,0.168_EB,0.152_EB,0.137_EB/)
ABSF(1,6,21:30,6) = (/0.130_EB,0.122_EB,0.114_EB,0.107_EB,0.099_EB,0.127_EB,0.154_EB,0.182_EB,0.209_EB,0.236_EB/)
ABSF(1,7,0:10,1) = (/0.000_EB,0.056_EB,0.066_EB,0.066_EB,0.065_EB,0.068_EB,0.070_EB,0.073_EB,0.076_EB,0.077_EB,0.079_EB/)
ABSF(1,7,11:20,1) = (/0.081_EB,0.082_EB,0.081_EB,0.081_EB,0.080_EB,0.079_EB,0.079_EB,0.080_EB,0.080_EB,0.080_EB/)
ABSF(1,7,21:30,1) = (/0.085_EB,0.090_EB,0.094_EB,0.099_EB,0.104_EB,0.102_EB,0.101_EB,0.099_EB,0.098_EB,0.096_EB/)
ABSF(1,7,0:10,2) = (/0.000_EB,0.053_EB,0.065_EB,0.066_EB,0.067_EB,0.069_EB,0.072_EB,0.074_EB,0.077_EB,0.079_EB,0.080_EB/)
ABSF(1,7,11:20,2) = (/0.082_EB,0.083_EB,0.082_EB,0.081_EB,0.080_EB,0.079_EB,0.080_EB,0.081_EB,0.082_EB,0.083_EB/)
ABSF(1,7,21:30,2) = (/0.087_EB,0.091_EB,0.096_EB,0.100_EB,0.104_EB,0.103_EB,0.102_EB,0.101_EB,0.099_EB,0.098_EB/)
ABSF(1,7,0:10,3) = (/0.000_EB,0.034_EB,0.059_EB,0.064_EB,0.069_EB,0.072_EB,0.075_EB,0.078_EB,0.081_EB,0.087_EB,0.093_EB/)
ABSF(1,7,11:20,3) = (/0.099_EB,0.106_EB,0.103_EB,0.101_EB,0.098_EB,0.096_EB,0.095_EB,0.094_EB,0.093_EB,0.093_EB/)
ABSF(1,7,21:30,3) = (/0.096_EB,0.099_EB,0.102_EB,0.105_EB,0.108_EB,0.112_EB,0.116_EB,0.121_EB,0.125_EB,0.129_EB/)
ABSF(1,7,0:10,4) = (/0.000_EB,0.033_EB,0.053_EB,0.067_EB,0.081_EB,0.083_EB,0.086_EB,0.088_EB,0.090_EB,0.103_EB,0.116_EB/)
ABSF(1,7,11:20,4) = (/0.129_EB,0.142_EB,0.136_EB,0.130_EB,0.124_EB,0.118_EB,0.118_EB,0.117_EB,0.117_EB,0.117_EB/)
ABSF(1,7,21:30,4) = (/0.119_EB,0.122_EB,0.124_EB,0.127_EB,0.129_EB,0.136_EB,0.143_EB,0.149_EB,0.156_EB,0.162_EB/)
ABSF(1,7,0:10,5) = (/0.000_EB,0.038_EB,0.056_EB,0.070_EB,0.083_EB,0.083_EB,0.082_EB,0.081_EB,0.080_EB,0.099_EB,0.118_EB/)
ABSF(1,7,11:20,5) = (/0.138_EB,0.157_EB,0.156_EB,0.156_EB,0.155_EB,0.154_EB,0.146_EB,0.138_EB,0.130_EB,0.122_EB/)
ABSF(1,7,21:30,5) = (/0.123_EB,0.124_EB,0.125_EB,0.125_EB,0.126_EB,0.141_EB,0.156_EB,0.171_EB,0.185_EB,0.200_EB/)
ABSF(1,7,0:10,6) = (/0.000_EB,0.038_EB,0.059_EB,0.066_EB,0.073_EB,0.078_EB,0.083_EB,0.088_EB,0.093_EB,0.110_EB,0.128_EB/)
ABSF(1,7,11:20,6) = (/0.145_EB,0.162_EB,0.172_EB,0.182_EB,0.192_EB,0.202_EB,0.189_EB,0.176_EB,0.163_EB,0.151_EB/)
ABSF(1,7,21:30,6) = (/0.141_EB,0.132_EB,0.122_EB,0.113_EB,0.103_EB,0.130_EB,0.157_EB,0.185_EB,0.212_EB,0.239_EB/)
ABSF(1,8,0:10,1) = (/0.000_EB,0.060_EB,0.072_EB,0.071_EB,0.071_EB,0.073_EB,0.076_EB,0.078_EB,0.081_EB,0.082_EB,0.084_EB/)
ABSF(1,8,11:20,1) = (/0.086_EB,0.087_EB,0.087_EB,0.086_EB,0.085_EB,0.084_EB,0.084_EB,0.084_EB,0.084_EB,0.084_EB/)
ABSF(1,8,21:30,1) = (/0.089_EB,0.094_EB,0.100_EB,0.105_EB,0.110_EB,0.108_EB,0.106_EB,0.104_EB,0.102_EB,0.101_EB/)
ABSF(1,8,0:10,2) = (/0.000_EB,0.058_EB,0.070_EB,0.071_EB,0.072_EB,0.074_EB,0.077_EB,0.079_EB,0.081_EB,0.083_EB,0.085_EB/)
ABSF(1,8,11:20,2) = (/0.087_EB,0.089_EB,0.087_EB,0.086_EB,0.085_EB,0.083_EB,0.084_EB,0.085_EB,0.086_EB,0.087_EB/)
ABSF(1,8,21:30,2) = (/0.092_EB,0.096_EB,0.101_EB,0.106_EB,0.111_EB,0.109_EB,0.107_EB,0.106_EB,0.104_EB,0.102_EB/)
ABSF(1,8,0:10,3) = (/0.000_EB,0.038_EB,0.066_EB,0.069_EB,0.073_EB,0.076_EB,0.079_EB,0.081_EB,0.084_EB,0.089_EB,0.094_EB/)
ABSF(1,8,11:20,3) = (/0.098_EB,0.103_EB,0.101_EB,0.099_EB,0.097_EB,0.095_EB,0.095_EB,0.095_EB,0.094_EB,0.094_EB/)
ABSF(1,8,21:30,3) = (/0.098_EB,0.102_EB,0.106_EB,0.110_EB,0.113_EB,0.117_EB,0.121_EB,0.124_EB,0.128_EB,0.131_EB/)
ABSF(1,8,0:10,4) = (/0.000_EB,0.036_EB,0.057_EB,0.071_EB,0.084_EB,0.087_EB,0.090_EB,0.092_EB,0.095_EB,0.108_EB,0.120_EB/)
ABSF(1,8,11:20,4) = (/0.133_EB,0.145_EB,0.140_EB,0.135_EB,0.129_EB,0.124_EB,0.123_EB,0.122_EB,0.121_EB,0.120_EB/)
ABSF(1,8,21:30,4) = (/0.123_EB,0.126_EB,0.129_EB,0.131_EB,0.134_EB,0.140_EB,0.146_EB,0.153_EB,0.159_EB,0.165_EB/)
ABSF(1,8,0:10,5) = (/0.000_EB,0.041_EB,0.061_EB,0.075_EB,0.089_EB,0.087_EB,0.086_EB,0.084_EB,0.082_EB,0.102_EB,0.121_EB/)
ABSF(1,8,11:20,5) = (/0.141_EB,0.161_EB,0.160_EB,0.160_EB,0.159_EB,0.159_EB,0.152_EB,0.145_EB,0.138_EB,0.130_EB/)
ABSF(1,8,21:30,5) = (/0.130_EB,0.130_EB,0.130_EB,0.130_EB,0.129_EB,0.144_EB,0.159_EB,0.174_EB,0.188_EB,0.203_EB/)
ABSF(1,8,0:10,6) = (/0.000_EB,0.040_EB,0.064_EB,0.070_EB,0.075_EB,0.082_EB,0.090_EB,0.097_EB,0.104_EB,0.119_EB,0.134_EB/)
ABSF(1,8,11:20,6) = (/0.150_EB,0.165_EB,0.175_EB,0.185_EB,0.196_EB,0.206_EB,0.196_EB,0.185_EB,0.175_EB,0.164_EB/)
ABSF(1,8,21:30,6) = (/0.153_EB,0.141_EB,0.130_EB,0.118_EB,0.107_EB,0.134_EB,0.161_EB,0.188_EB,0.214_EB,0.241_EB/)
ABSF(1,9,0:10,1) = (/0.000_EB,0.065_EB,0.077_EB,0.077_EB,0.076_EB,0.079_EB,0.081_EB,0.083_EB,0.086_EB,0.087_EB,0.089_EB/)
ABSF(1,9,11:20,1) = (/0.091_EB,0.093_EB,0.092_EB,0.091_EB,0.090_EB,0.089_EB,0.088_EB,0.088_EB,0.088_EB,0.088_EB/)
ABSF(1,9,21:30,1) = (/0.094_EB,0.099_EB,0.105_EB,0.111_EB,0.116_EB,0.114_EB,0.112_EB,0.109_EB,0.107_EB,0.105_EB/)
ABSF(1,9,0:10,2) = (/0.000_EB,0.064_EB,0.074_EB,0.076_EB,0.077_EB,0.080_EB,0.082_EB,0.084_EB,0.086_EB,0.088_EB,0.090_EB/)
ABSF(1,9,11:20,2) = (/0.092_EB,0.094_EB,0.092_EB,0.091_EB,0.090_EB,0.088_EB,0.089_EB,0.090_EB,0.090_EB,0.091_EB/)
ABSF(1,9,21:30,2) = (/0.096_EB,0.102_EB,0.107_EB,0.112_EB,0.117_EB,0.115_EB,0.113_EB,0.111_EB,0.108_EB,0.106_EB/)
ABSF(1,9,0:10,3) = (/0.000_EB,0.042_EB,0.073_EB,0.075_EB,0.078_EB,0.080_EB,0.082_EB,0.085_EB,0.087_EB,0.090_EB,0.094_EB/)
ABSF(1,9,11:20,3) = (/0.097_EB,0.101_EB,0.099_EB,0.098_EB,0.096_EB,0.094_EB,0.095_EB,0.095_EB,0.095_EB,0.096_EB/)
ABSF(1,9,21:30,3) = (/0.100_EB,0.105_EB,0.110_EB,0.114_EB,0.119_EB,0.122_EB,0.125_EB,0.127_EB,0.130_EB,0.133_EB/)
ABSF(1,9,0:10,4) = (/0.000_EB,0.040_EB,0.060_EB,0.074_EB,0.088_EB,0.091_EB,0.094_EB,0.097_EB,0.100_EB,0.112_EB,0.124_EB/)
ABSF(1,9,11:20,4) = (/0.137_EB,0.149_EB,0.144_EB,0.140_EB,0.135_EB,0.130_EB,0.129_EB,0.127_EB,0.125_EB,0.124_EB/)
ABSF(1,9,21:30,4) = (/0.127_EB,0.130_EB,0.133_EB,0.136_EB,0.139_EB,0.145_EB,0.150_EB,0.156_EB,0.162_EB,0.168_EB/)
ABSF(1,9,0:10,5) = (/0.000_EB,0.044_EB,0.065_EB,0.080_EB,0.094_EB,0.092_EB,0.090_EB,0.087_EB,0.085_EB,0.105_EB,0.124_EB/)
ABSF(1,9,11:20,5) = (/0.144_EB,0.164_EB,0.164_EB,0.164_EB,0.164_EB,0.164_EB,0.157_EB,0.151_EB,0.145_EB,0.139_EB/)
ABSF(1,9,21:30,5) = (/0.137_EB,0.136_EB,0.135_EB,0.134_EB,0.132_EB,0.147_EB,0.162_EB,0.177_EB,0.191_EB,0.206_EB/)
ABSF(1,9,0:10,6) = (/0.000_EB,0.042_EB,0.068_EB,0.073_EB,0.078_EB,0.087_EB,0.096_EB,0.105_EB,0.114_EB,0.128_EB,0.141_EB/)
ABSF(1,9,11:20,6) = (/0.154_EB,0.168_EB,0.178_EB,0.189_EB,0.200_EB,0.210_EB,0.202_EB,0.194_EB,0.186_EB,0.177_EB/)
ABSF(1,9,21:30,6) = (/0.164_EB,0.151_EB,0.137_EB,0.124_EB,0.110_EB,0.137_EB,0.164_EB,0.191_EB,0.217_EB,0.244_EB/)
ABSF(1,10,0:10,1) = (/0.000_EB,0.069_EB,0.083_EB,0.082_EB,0.082_EB,0.084_EB,0.086_EB,0.088_EB,0.091_EB,0.092_EB,0.094_EB/)
ABSF(1,10,11:20,1) = (/0.096_EB,0.098_EB,0.097_EB,0.096_EB,0.095_EB,0.093_EB,0.093_EB,0.093_EB,0.092_EB,0.092_EB/)
ABSF(1,10,21:30,1) = (/0.098_EB,0.104_EB,0.110_EB,0.117_EB,0.123_EB,0.120_EB,0.117_EB,0.114_EB,0.111_EB,0.109_EB/)
ABSF(1,10,0:10,2) = (/0.000_EB,0.069_EB,0.079_EB,0.081_EB,0.083_EB,0.085_EB,0.087_EB,0.089_EB,0.091_EB,0.093_EB,0.095_EB/)
ABSF(1,10,11:20,2) = (/0.097_EB,0.099_EB,0.097_EB,0.096_EB,0.095_EB,0.093_EB,0.094_EB,0.094_EB,0.095_EB,0.095_EB/)
ABSF(1,10,21:30,2) = (/0.101_EB,0.107_EB,0.112_EB,0.118_EB,0.123_EB,0.121_EB,0.118_EB,0.116_EB,0.113_EB,0.110_EB/)
ABSF(1,10,0:10,3) = (/0.000_EB,0.046_EB,0.079_EB,0.081_EB,0.082_EB,0.084_EB,0.086_EB,0.088_EB,0.090_EB,0.092_EB,0.094_EB/)
ABSF(1,10,11:20,3) = (/0.096_EB,0.098_EB,0.097_EB,0.096_EB,0.095_EB,0.094_EB,0.095_EB,0.095_EB,0.096_EB,0.097_EB/)
ABSF(1,10,21:30,3) = (/0.103_EB,0.108_EB,0.114_EB,0.119_EB,0.125_EB,0.127_EB,0.129_EB,0.131_EB,0.133_EB,0.135_EB/)
ABSF(1,10,0:10,4) = (/0.000_EB,0.043_EB,0.064_EB,0.078_EB,0.091_EB,0.095_EB,0.098_EB,0.101_EB,0.105_EB,0.117_EB,0.129_EB/)
ABSF(1,10,11:20,4) = (/0.141_EB,0.153_EB,0.149_EB,0.144_EB,0.140_EB,0.136_EB,0.134_EB,0.132_EB,0.129_EB,0.127_EB/)
ABSF(1,10,21:30,4) = (/0.130_EB,0.134_EB,0.137_EB,0.140_EB,0.144_EB,0.149_EB,0.154_EB,0.160_EB,0.165_EB,0.170_EB/)
ABSF(1,10,0:10,5) = (/0.000_EB,0.046_EB,0.070_EB,0.085_EB,0.100_EB,0.097_EB,0.094_EB,0.090_EB,0.087_EB,0.107_EB,0.128_EB/)
ABSF(1,10,11:20,5) = (/0.148_EB,0.168_EB,0.168_EB,0.168_EB,0.168_EB,0.168_EB,0.163_EB,0.158_EB,0.152_EB,0.147_EB/)
ABSF(1,10,21:30,5) = (/0.145_EB,0.142_EB,0.140_EB,0.138_EB,0.135_EB,0.150_EB,0.165_EB,0.180_EB,0.194_EB,0.209_EB/)
ABSF(1,10,0:10,6) = (/0.000_EB,0.044_EB,0.073_EB,0.077_EB,0.081_EB,0.092_EB,0.103_EB,0.114_EB,0.125_EB,0.136_EB,0.148_EB/)
ABSF(1,10,11:20,6) = (/0.159_EB,0.171_EB,0.182_EB,0.192_EB,0.203_EB,0.214_EB,0.208_EB,0.203_EB,0.197_EB,0.191_EB/)
ABSF(1,10,21:30,6) = (/0.176_EB,0.160_EB,0.145_EB,0.129_EB,0.114_EB,0.140_EB,0.167_EB,0.194_EB,0.220_EB,0.247_EB/)
ABSF(1,11,0:10,1) = (/0.000_EB,0.073_EB,0.087_EB,0.087_EB,0.087_EB,0.089_EB,0.091_EB,0.093_EB,0.095_EB,0.097_EB,0.099_EB/)
ABSF(1,11,11:20,1) = (/0.100_EB,0.102_EB,0.101_EB,0.100_EB,0.099_EB,0.098_EB,0.097_EB,0.097_EB,0.097_EB,0.096_EB/)
ABSF(1,11,21:30,1) = (/0.102_EB,0.109_EB,0.115_EB,0.121_EB,0.127_EB,0.124_EB,0.121_EB,0.118_EB,0.115_EB,0.112_EB/)
ABSF(1,11,0:10,2) = (/0.000_EB,0.073_EB,0.084_EB,0.086_EB,0.089_EB,0.090_EB,0.092_EB,0.094_EB,0.095_EB,0.097_EB,0.099_EB/)
ABSF(1,11,11:20,2) = (/0.101_EB,0.103_EB,0.102_EB,0.100_EB,0.099_EB,0.097_EB,0.098_EB,0.098_EB,0.099_EB,0.099_EB/)
ABSF(1,11,21:30,2) = (/0.105_EB,0.111_EB,0.116_EB,0.122_EB,0.127_EB,0.125_EB,0.122_EB,0.119_EB,0.117_EB,0.114_EB/)
ABSF(1,11,0:10,3) = (/0.000_EB,0.052_EB,0.085_EB,0.086_EB,0.087_EB,0.089_EB,0.091_EB,0.093_EB,0.095_EB,0.097_EB,0.099_EB/)
ABSF(1,11,11:20,3) = (/0.101_EB,0.103_EB,0.102_EB,0.100_EB,0.099_EB,0.098_EB,0.099_EB,0.099_EB,0.100_EB,0.101_EB/)
ABSF(1,11,21:30,3) = (/0.107_EB,0.112_EB,0.118_EB,0.123_EB,0.129_EB,0.130_EB,0.132_EB,0.133_EB,0.134_EB,0.136_EB/)
ABSF(1,11,0:10,4) = (/0.000_EB,0.046_EB,0.069_EB,0.082_EB,0.096_EB,0.099_EB,0.102_EB,0.105_EB,0.108_EB,0.120_EB,0.131_EB/)
ABSF(1,11,11:20,4) = (/0.142_EB,0.154_EB,0.150_EB,0.147_EB,0.144_EB,0.140_EB,0.138_EB,0.135_EB,0.133_EB,0.130_EB/)
ABSF(1,11,21:30,4) = (/0.134_EB,0.137_EB,0.140_EB,0.143_EB,0.147_EB,0.152_EB,0.157_EB,0.162_EB,0.167_EB,0.172_EB/)
ABSF(1,11,0:10,5) = (/0.000_EB,0.050_EB,0.074_EB,0.088_EB,0.103_EB,0.100_EB,0.098_EB,0.095_EB,0.092_EB,0.112_EB,0.132_EB/)
ABSF(1,11,11:20,5) = (/0.152_EB,0.171_EB,0.172_EB,0.172_EB,0.172_EB,0.173_EB,0.167_EB,0.161_EB,0.155_EB,0.150_EB/)
ABSF(1,11,21:30,5) = (/0.147_EB,0.145_EB,0.143_EB,0.140_EB,0.138_EB,0.152_EB,0.166_EB,0.181_EB,0.195_EB,0.209_EB/)
ABSF(1,11,0:10,6) = (/0.000_EB,0.047_EB,0.077_EB,0.081_EB,0.085_EB,0.096_EB,0.107_EB,0.118_EB,0.129_EB,0.140_EB,0.152_EB/)
ABSF(1,11,11:20,6) = (/0.163_EB,0.174_EB,0.184_EB,0.195_EB,0.205_EB,0.215_EB,0.211_EB,0.207_EB,0.203_EB,0.199_EB/)
ABSF(1,11,21:30,6) = (/0.183_EB,0.167_EB,0.151_EB,0.135_EB,0.119_EB,0.145_EB,0.172_EB,0.198_EB,0.225_EB,0.251_EB/)
ABSF(1,12,0:10,1) = (/0.000_EB,0.077_EB,0.091_EB,0.092_EB,0.092_EB,0.094_EB,0.096_EB,0.098_EB,0.100_EB,0.101_EB,0.103_EB/)
ABSF(1,12,11:20,1) = (/0.105_EB,0.107_EB,0.105_EB,0.104_EB,0.103_EB,0.102_EB,0.101_EB,0.101_EB,0.101_EB,0.100_EB/)
ABSF(1,12,21:30,1) = (/0.107_EB,0.113_EB,0.119_EB,0.125_EB,0.131_EB,0.128_EB,0.125_EB,0.122_EB,0.119_EB,0.116_EB/)
ABSF(1,12,0:10,2) = (/0.000_EB,0.077_EB,0.088_EB,0.091_EB,0.094_EB,0.096_EB,0.097_EB,0.099_EB,0.100_EB,0.102_EB,0.104_EB/)
ABSF(1,12,11:20,2) = (/0.106_EB,0.108_EB,0.106_EB,0.105_EB,0.103_EB,0.102_EB,0.102_EB,0.102_EB,0.103_EB,0.103_EB/)
ABSF(1,12,21:30,2) = (/0.109_EB,0.114_EB,0.120_EB,0.126_EB,0.131_EB,0.129_EB,0.126_EB,0.123_EB,0.120_EB,0.117_EB/)
ABSF(1,12,0:10,3) = (/0.000_EB,0.058_EB,0.090_EB,0.091_EB,0.092_EB,0.094_EB,0.096_EB,0.098_EB,0.100_EB,0.102_EB,0.103_EB/)
ABSF(1,12,11:20,3) = (/0.105_EB,0.107_EB,0.106_EB,0.105_EB,0.103_EB,0.102_EB,0.103_EB,0.103_EB,0.104_EB,0.105_EB/)
ABSF(1,12,21:30,3) = (/0.110_EB,0.116_EB,0.122_EB,0.127_EB,0.133_EB,0.134_EB,0.134_EB,0.135_EB,0.136_EB,0.137_EB/)
ABSF(1,12,0:10,4) = (/0.000_EB,0.049_EB,0.074_EB,0.087_EB,0.100_EB,0.103_EB,0.106_EB,0.109_EB,0.111_EB,0.122_EB,0.133_EB/)
ABSF(1,12,11:20,4) = (/0.144_EB,0.155_EB,0.152_EB,0.150_EB,0.147_EB,0.144_EB,0.142_EB,0.139_EB,0.136_EB,0.134_EB/)
ABSF(1,12,21:30,4) = (/0.137_EB,0.140_EB,0.143_EB,0.146_EB,0.150_EB,0.154_EB,0.159_EB,0.164_EB,0.168_EB,0.173_EB/)
ABSF(1,12,0:10,5) = (/0.000_EB,0.053_EB,0.077_EB,0.092_EB,0.106_EB,0.104_EB,0.102_EB,0.100_EB,0.098_EB,0.117_EB,0.136_EB/)
ABSF(1,12,11:20,5) = (/0.156_EB,0.175_EB,0.175_EB,0.176_EB,0.177_EB,0.177_EB,0.171_EB,0.165_EB,0.159_EB,0.152_EB/)
ABSF(1,12,21:30,5) = (/0.150_EB,0.148_EB,0.145_EB,0.143_EB,0.141_EB,0.154_EB,0.168_EB,0.182_EB,0.195_EB,0.209_EB/)
ABSF(1,12,0:10,6) = (/0.000_EB,0.050_EB,0.081_EB,0.085_EB,0.089_EB,0.100_EB,0.111_EB,0.122_EB,0.134_EB,0.145_EB,0.156_EB/)
ABSF(1,12,11:20,6) = (/0.167_EB,0.178_EB,0.187_EB,0.197_EB,0.206_EB,0.216_EB,0.213_EB,0.211_EB,0.209_EB,0.207_EB/)
ABSF(1,12,21:30,6) = (/0.190_EB,0.174_EB,0.157_EB,0.141_EB,0.124_EB,0.150_EB,0.177_EB,0.203_EB,0.229_EB,0.256_EB/)
ABSF(1,13,0:10,1) = (/0.000_EB,0.081_EB,0.096_EB,0.097_EB,0.098_EB,0.099_EB,0.101_EB,0.103_EB,0.104_EB,0.106_EB,0.108_EB/)
ABSF(1,13,11:20,1) = (/0.109_EB,0.111_EB,0.110_EB,0.109_EB,0.107_EB,0.106_EB,0.106_EB,0.105_EB,0.105_EB,0.105_EB/)
ABSF(1,13,21:30,1) = (/0.111_EB,0.117_EB,0.123_EB,0.129_EB,0.135_EB,0.132_EB,0.129_EB,0.126_EB,0.122_EB,0.119_EB/)
ABSF(1,13,0:10,2) = (/0.000_EB,0.080_EB,0.093_EB,0.096_EB,0.100_EB,0.101_EB,0.102_EB,0.104_EB,0.105_EB,0.107_EB,0.109_EB/)
ABSF(1,13,11:20,2) = (/0.111_EB,0.112_EB,0.111_EB,0.109_EB,0.108_EB,0.106_EB,0.106_EB,0.106_EB,0.107_EB,0.107_EB/)
ABSF(1,13,21:30,2) = (/0.113_EB,0.118_EB,0.124_EB,0.130_EB,0.135_EB,0.133_EB,0.130_EB,0.127_EB,0.124_EB,0.121_EB/)
ABSF(1,13,0:10,3) = (/0.000_EB,0.064_EB,0.096_EB,0.096_EB,0.097_EB,0.099_EB,0.101_EB,0.103_EB,0.105_EB,0.106_EB,0.108_EB/)
ABSF(1,13,11:20,3) = (/0.110_EB,0.112_EB,0.110_EB,0.109_EB,0.107_EB,0.106_EB,0.107_EB,0.107_EB,0.108_EB,0.108_EB/)
ABSF(1,13,21:30,3) = (/0.114_EB,0.120_EB,0.125_EB,0.131_EB,0.137_EB,0.137_EB,0.137_EB,0.137_EB,0.138_EB,0.138_EB/)
ABSF(1,13,0:10,4) = (/0.000_EB,0.053_EB,0.080_EB,0.092_EB,0.104_EB,0.107_EB,0.109_EB,0.112_EB,0.115_EB,0.125_EB,0.135_EB/)
ABSF(1,13,11:20,4) = (/0.146_EB,0.156_EB,0.154_EB,0.152_EB,0.150_EB,0.148_EB,0.146_EB,0.143_EB,0.140_EB,0.137_EB/)
ABSF(1,13,21:30,4) = (/0.140_EB,0.143_EB,0.146_EB,0.149_EB,0.152_EB,0.157_EB,0.161_EB,0.166_EB,0.170_EB,0.175_EB/)
ABSF(1,13,0:10,5) = (/0.000_EB,0.056_EB,0.081_EB,0.095_EB,0.110_EB,0.108_EB,0.106_EB,0.105_EB,0.103_EB,0.122_EB,0.141_EB/)
ABSF(1,13,11:20,5) = (/0.159_EB,0.178_EB,0.179_EB,0.180_EB,0.181_EB,0.182_EB,0.175_EB,0.169_EB,0.162_EB,0.155_EB/)
ABSF(1,13,21:30,5) = (/0.153_EB,0.151_EB,0.148_EB,0.146_EB,0.143_EB,0.156_EB,0.169_EB,0.183_EB,0.196_EB,0.209_EB/)
ABSF(1,13,0:10,6) = (/0.000_EB,0.052_EB,0.085_EB,0.089_EB,0.093_EB,0.104_EB,0.115_EB,0.127_EB,0.138_EB,0.149_EB,0.160_EB/)
ABSF(1,13,11:20,6) = (/0.171_EB,0.181_EB,0.190_EB,0.199_EB,0.208_EB,0.216_EB,0.216_EB,0.215_EB,0.215_EB,0.214_EB/)
ABSF(1,13,21:30,6) = (/0.197_EB,0.180_EB,0.163_EB,0.146_EB,0.129_EB,0.156_EB,0.182_EB,0.208_EB,0.234_EB,0.260_EB/)
ABSF(1,14,0:10,1) = (/0.000_EB,0.085_EB,0.100_EB,0.101_EB,0.103_EB,0.105_EB,0.106_EB,0.107_EB,0.109_EB,0.110_EB,0.112_EB/)
ABSF(1,14,11:20,1) = (/0.114_EB,0.116_EB,0.114_EB,0.113_EB,0.111_EB,0.110_EB,0.110_EB,0.109_EB,0.109_EB,0.109_EB/)
ABSF(1,14,21:30,1) = (/0.115_EB,0.121_EB,0.127_EB,0.133_EB,0.140_EB,0.136_EB,0.133_EB,0.129_EB,0.126_EB,0.123_EB/)
ABSF(1,14,0:10,2) = (/0.000_EB,0.084_EB,0.098_EB,0.101_EB,0.105_EB,0.106_EB,0.108_EB,0.109_EB,0.110_EB,0.112_EB,0.113_EB/)
ABSF(1,14,11:20,2) = (/0.115_EB,0.117_EB,0.115_EB,0.114_EB,0.112_EB,0.110_EB,0.110_EB,0.111_EB,0.111_EB,0.111_EB/)
ABSF(1,14,21:30,2) = (/0.117_EB,0.122_EB,0.128_EB,0.134_EB,0.139_EB,0.137_EB,0.134_EB,0.131_EB,0.128_EB,0.125_EB/)
ABSF(1,14,0:10,3) = (/0.000_EB,0.070_EB,0.101_EB,0.102_EB,0.102_EB,0.104_EB,0.106_EB,0.108_EB,0.109_EB,0.111_EB,0.113_EB/)
ABSF(1,14,11:20,3) = (/0.114_EB,0.116_EB,0.115_EB,0.113_EB,0.112_EB,0.110_EB,0.111_EB,0.111_EB,0.112_EB,0.112_EB/)
ABSF(1,14,21:30,3) = (/0.118_EB,0.124_EB,0.129_EB,0.135_EB,0.141_EB,0.140_EB,0.140_EB,0.140_EB,0.139_EB,0.139_EB/)
ABSF(1,14,0:10,4) = (/0.000_EB,0.056_EB,0.085_EB,0.096_EB,0.108_EB,0.111_EB,0.113_EB,0.116_EB,0.118_EB,0.128_EB,0.138_EB/)
ABSF(1,14,11:20,4) = (/0.148_EB,0.157_EB,0.156_EB,0.155_EB,0.154_EB,0.153_EB,0.150_EB,0.147_EB,0.144_EB,0.141_EB/)
ABSF(1,14,21:30,4) = (/0.144_EB,0.146_EB,0.149_EB,0.152_EB,0.155_EB,0.159_EB,0.164_EB,0.168_EB,0.172_EB,0.176_EB/)
ABSF(1,14,0:10,5) = (/0.000_EB,0.059_EB,0.085_EB,0.099_EB,0.113_EB,0.112_EB,0.111_EB,0.109_EB,0.108_EB,0.126_EB,0.145_EB/)
ABSF(1,14,11:20,5) = (/0.163_EB,0.182_EB,0.183_EB,0.184_EB,0.185_EB,0.186_EB,0.179_EB,0.172_EB,0.165_EB,0.158_EB/)
ABSF(1,14,21:30,5) = (/0.156_EB,0.153_EB,0.151_EB,0.148_EB,0.146_EB,0.159_EB,0.171_EB,0.184_EB,0.196_EB,0.208_EB/)
ABSF(1,14,0:10,6) = (/0.000_EB,0.055_EB,0.088_EB,0.093_EB,0.097_EB,0.108_EB,0.120_EB,0.131_EB,0.143_EB,0.153_EB,0.164_EB/)
ABSF(1,14,11:20,6) = (/0.174_EB,0.185_EB,0.193_EB,0.201_EB,0.209_EB,0.217_EB,0.218_EB,0.220_EB,0.221_EB,0.222_EB/)
ABSF(1,14,21:30,6) = (/0.205_EB,0.187_EB,0.170_EB,0.152_EB,0.134_EB,0.161_EB,0.187_EB,0.213_EB,0.239_EB,0.265_EB/)
ABSF(1,15,0:10,1) = (/0.000_EB,0.089_EB,0.104_EB,0.106_EB,0.108_EB,0.110_EB,0.111_EB,0.112_EB,0.113_EB,0.115_EB,0.117_EB/)
ABSF(1,15,11:20,1) = (/0.118_EB,0.120_EB,0.119_EB,0.117_EB,0.116_EB,0.114_EB,0.114_EB,0.114_EB,0.113_EB,0.113_EB/)
ABSF(1,15,21:30,1) = (/0.119_EB,0.125_EB,0.132_EB,0.138_EB,0.144_EB,0.140_EB,0.137_EB,0.133_EB,0.130_EB,0.126_EB/)
ABSF(1,15,0:10,2) = (/0.000_EB,0.088_EB,0.102_EB,0.107_EB,0.111_EB,0.112_EB,0.113_EB,0.114_EB,0.115_EB,0.116_EB,0.118_EB/)
ABSF(1,15,11:20,2) = (/0.120_EB,0.122_EB,0.120_EB,0.118_EB,0.116_EB,0.114_EB,0.114_EB,0.115_EB,0.115_EB,0.115_EB/)
ABSF(1,15,21:30,2) = (/0.121_EB,0.126_EB,0.132_EB,0.138_EB,0.143_EB,0.140_EB,0.137_EB,0.134_EB,0.131_EB,0.128_EB/)
ABSF(1,15,0:10,3) = (/0.000_EB,0.076_EB,0.107_EB,0.107_EB,0.107_EB,0.109_EB,0.111_EB,0.112_EB,0.114_EB,0.116_EB,0.117_EB/)
ABSF(1,15,11:20,3) = (/0.119_EB,0.121_EB,0.119_EB,0.117_EB,0.116_EB,0.114_EB,0.115_EB,0.115_EB,0.115_EB,0.116_EB/)
ABSF(1,15,21:30,3) = (/0.122_EB,0.127_EB,0.133_EB,0.139_EB,0.145_EB,0.144_EB,0.143_EB,0.142_EB,0.141_EB,0.140_EB/)
ABSF(1,15,0:10,4) = (/0.000_EB,0.059_EB,0.090_EB,0.101_EB,0.112_EB,0.114_EB,0.117_EB,0.119_EB,0.121_EB,0.131_EB,0.140_EB/)
ABSF(1,15,11:20,4) = (/0.149_EB,0.159_EB,0.158_EB,0.158_EB,0.157_EB,0.157_EB,0.153_EB,0.150_EB,0.147_EB,0.144_EB/)
ABSF(1,15,21:30,4) = (/0.147_EB,0.150_EB,0.153_EB,0.155_EB,0.158_EB,0.162_EB,0.166_EB,0.170_EB,0.174_EB,0.178_EB/)
ABSF(1,15,0:10,5) = (/0.000_EB,0.062_EB,0.089_EB,0.103_EB,0.116_EB,0.116_EB,0.115_EB,0.114_EB,0.113_EB,0.131_EB,0.149_EB/)
ABSF(1,15,11:20,5) = (/0.167_EB,0.185_EB,0.187_EB,0.188_EB,0.190_EB,0.191_EB,0.183_EB,0.176_EB,0.168_EB,0.161_EB/)
ABSF(1,15,21:30,5) = (/0.158_EB,0.156_EB,0.154_EB,0.151_EB,0.149_EB,0.161_EB,0.173_EB,0.185_EB,0.196_EB,0.208_EB/)
ABSF(1,15,0:10,6) = (/0.000_EB,0.058_EB,0.092_EB,0.097_EB,0.101_EB,0.112_EB,0.124_EB,0.135_EB,0.147_EB,0.157_EB,0.168_EB/)
ABSF(1,15,11:20,6) = (/0.178_EB,0.189_EB,0.196_EB,0.203_EB,0.211_EB,0.218_EB,0.221_EB,0.224_EB,0.227_EB,0.230_EB/)
ABSF(1,15,21:30,6) = (/0.212_EB,0.194_EB,0.176_EB,0.158_EB,0.140_EB,0.166_EB,0.191_EB,0.217_EB,0.243_EB,0.269_EB/)
ABSF(1,16,0:10,1) = (/0.000_EB,0.093_EB,0.108_EB,0.111_EB,0.114_EB,0.115_EB,0.116_EB,0.117_EB,0.118_EB,0.120_EB,0.121_EB/)
ABSF(1,16,11:20,1) = (/0.123_EB,0.125_EB,0.123_EB,0.121_EB,0.120_EB,0.118_EB,0.118_EB,0.118_EB,0.117_EB,0.117_EB/)
ABSF(1,16,21:30,1) = (/0.123_EB,0.130_EB,0.136_EB,0.142_EB,0.148_EB,0.144_EB,0.141_EB,0.137_EB,0.133_EB,0.130_EB/)
ABSF(1,16,0:10,2) = (/0.000_EB,0.092_EB,0.107_EB,0.112_EB,0.117_EB,0.117_EB,0.118_EB,0.119_EB,0.120_EB,0.121_EB,0.123_EB/)
ABSF(1,16,11:20,2) = (/0.124_EB,0.126_EB,0.124_EB,0.122_EB,0.120_EB,0.119_EB,0.119_EB,0.119_EB,0.119_EB,0.119_EB/)
ABSF(1,16,21:30,2) = (/0.124_EB,0.130_EB,0.136_EB,0.142_EB,0.147_EB,0.144_EB,0.141_EB,0.138_EB,0.135_EB,0.132_EB/)
ABSF(1,16,0:10,3) = (/0.000_EB,0.082_EB,0.112_EB,0.112_EB,0.112_EB,0.114_EB,0.115_EB,0.117_EB,0.119_EB,0.121_EB,0.122_EB/)
ABSF(1,16,11:20,3) = (/0.124_EB,0.125_EB,0.123_EB,0.122_EB,0.120_EB,0.119_EB,0.119_EB,0.119_EB,0.119_EB,0.120_EB/)
ABSF(1,16,21:30,3) = (/0.125_EB,0.131_EB,0.137_EB,0.143_EB,0.149_EB,0.147_EB,0.145_EB,0.144_EB,0.142_EB,0.141_EB/)
ABSF(1,16,0:10,4) = (/0.000_EB,0.063_EB,0.095_EB,0.106_EB,0.116_EB,0.118_EB,0.120_EB,0.123_EB,0.125_EB,0.133_EB,0.142_EB/)
ABSF(1,16,11:20,4) = (/0.151_EB,0.160_EB,0.160_EB,0.160_EB,0.160_EB,0.161_EB,0.157_EB,0.154_EB,0.151_EB,0.147_EB/)
ABSF(1,16,21:30,4) = (/0.150_EB,0.153_EB,0.156_EB,0.158_EB,0.161_EB,0.165_EB,0.168_EB,0.172_EB,0.176_EB,0.179_EB/)
ABSF(1,16,0:10,5) = (/0.000_EB,0.065_EB,0.092_EB,0.106_EB,0.120_EB,0.119_EB,0.119_EB,0.119_EB,0.118_EB,0.136_EB,0.154_EB/)
ABSF(1,16,11:20,5) = (/0.171_EB,0.189_EB,0.191_EB,0.192_EB,0.194_EB,0.195_EB,0.188_EB,0.180_EB,0.172_EB,0.164_EB/)
ABSF(1,16,21:30,5) = (/0.161_EB,0.159_EB,0.156_EB,0.154_EB,0.151_EB,0.163_EB,0.174_EB,0.186_EB,0.197_EB,0.208_EB/)
ABSF(1,16,0:10,6) = (/0.000_EB,0.061_EB,0.096_EB,0.100_EB,0.105_EB,0.116_EB,0.128_EB,0.140_EB,0.151_EB,0.162_EB,0.172_EB/)
ABSF(1,16,11:20,6) = (/0.182_EB,0.192_EB,0.199_EB,0.205_EB,0.212_EB,0.219_EB,0.223_EB,0.228_EB,0.233_EB,0.238_EB/)
ABSF(1,16,21:30,6) = (/0.219_EB,0.201_EB,0.182_EB,0.163_EB,0.145_EB,0.171_EB,0.196_EB,0.222_EB,0.248_EB,0.274_EB/)
ABSF(1,17,0:10,1) = (/0.000_EB,0.097_EB,0.113_EB,0.116_EB,0.119_EB,0.120_EB,0.121_EB,0.122_EB,0.122_EB,0.124_EB,0.126_EB/)
ABSF(1,17,11:20,1) = (/0.127_EB,0.129_EB,0.128_EB,0.126_EB,0.124_EB,0.122_EB,0.122_EB,0.122_EB,0.122_EB,0.121_EB/)
ABSF(1,17,21:30,1) = (/0.128_EB,0.134_EB,0.140_EB,0.146_EB,0.152_EB,0.149_EB,0.145_EB,0.141_EB,0.137_EB,0.133_EB/)
ABSF(1,17,0:10,2) = (/0.000_EB,0.096_EB,0.112_EB,0.117_EB,0.122_EB,0.123_EB,0.123_EB,0.124_EB,0.124_EB,0.126_EB,0.127_EB/)
ABSF(1,17,11:20,2) = (/0.129_EB,0.131_EB,0.129_EB,0.127_EB,0.125_EB,0.123_EB,0.123_EB,0.123_EB,0.123_EB,0.123_EB/)
ABSF(1,17,21:30,2) = (/0.128_EB,0.134_EB,0.140_EB,0.146_EB,0.152_EB,0.148_EB,0.145_EB,0.142_EB,0.139_EB,0.135_EB/)
ABSF(1,17,0:10,3) = (/0.000_EB,0.088_EB,0.118_EB,0.117_EB,0.117_EB,0.118_EB,0.120_EB,0.122_EB,0.124_EB,0.125_EB,0.127_EB/)
ABSF(1,17,11:20,3) = (/0.128_EB,0.129_EB,0.128_EB,0.126_EB,0.124_EB,0.123_EB,0.123_EB,0.123_EB,0.123_EB,0.123_EB/)
ABSF(1,17,21:30,3) = (/0.129_EB,0.135_EB,0.141_EB,0.147_EB,0.153_EB,0.150_EB,0.148_EB,0.146_EB,0.144_EB,0.141_EB/)
ABSF(1,17,0:10,4) = (/0.000_EB,0.066_EB,0.100_EB,0.110_EB,0.120_EB,0.122_EB,0.124_EB,0.126_EB,0.128_EB,0.136_EB,0.144_EB/)
ABSF(1,17,11:20,4) = (/0.153_EB,0.161_EB,0.162_EB,0.163_EB,0.164_EB,0.165_EB,0.161_EB,0.158_EB,0.154_EB,0.151_EB/)
ABSF(1,17,21:30,4) = (/0.153_EB,0.156_EB,0.159_EB,0.161_EB,0.164_EB,0.167_EB,0.171_EB,0.174_EB,0.177_EB,0.181_EB/)
ABSF(1,17,0:10,5) = (/0.000_EB,0.069_EB,0.096_EB,0.110_EB,0.123_EB,0.123_EB,0.123_EB,0.123_EB,0.124_EB,0.141_EB,0.158_EB/)
ABSF(1,17,11:20,5) = (/0.175_EB,0.192_EB,0.194_EB,0.196_EB,0.198_EB,0.200_EB,0.192_EB,0.183_EB,0.175_EB,0.166_EB/)
ABSF(1,17,21:30,5) = (/0.164_EB,0.162_EB,0.159_EB,0.157_EB,0.154_EB,0.165_EB,0.176_EB,0.186_EB,0.197_EB,0.208_EB/)
ABSF(1,17,0:10,6) = (/0.000_EB,0.064_EB,0.100_EB,0.104_EB,0.109_EB,0.121_EB,0.132_EB,0.144_EB,0.156_EB,0.166_EB,0.176_EB/)
ABSF(1,17,11:20,6) = (/0.186_EB,0.196_EB,0.202_EB,0.208_EB,0.214_EB,0.219_EB,0.226_EB,0.233_EB,0.239_EB,0.246_EB/)
ABSF(1,17,21:30,6) = (/0.227_EB,0.207_EB,0.188_EB,0.169_EB,0.150_EB,0.176_EB,0.201_EB,0.227_EB,0.253_EB,0.278_EB/)
ABSF(1,18,0:10,1) = (/0.000_EB,0.100_EB,0.117_EB,0.121_EB,0.124_EB,0.125_EB,0.126_EB,0.126_EB,0.127_EB,0.129_EB,0.130_EB/)
ABSF(1,18,11:20,1) = (/0.132_EB,0.134_EB,0.132_EB,0.130_EB,0.128_EB,0.127_EB,0.126_EB,0.126_EB,0.126_EB,0.126_EB/)
ABSF(1,18,21:30,1) = (/0.132_EB,0.138_EB,0.144_EB,0.150_EB,0.157_EB,0.153_EB,0.149_EB,0.145_EB,0.141_EB,0.136_EB/)
ABSF(1,18,0:10,2) = (/0.000_EB,0.100_EB,0.116_EB,0.122_EB,0.128_EB,0.128_EB,0.128_EB,0.129_EB,0.129_EB,0.131_EB,0.132_EB/)
ABSF(1,18,11:20,2) = (/0.134_EB,0.135_EB,0.133_EB,0.131_EB,0.129_EB,0.127_EB,0.127_EB,0.127_EB,0.127_EB,0.126_EB/)
ABSF(1,18,21:30,2) = (/0.132_EB,0.138_EB,0.144_EB,0.150_EB,0.156_EB,0.152_EB,0.149_EB,0.146_EB,0.142_EB,0.139_EB/)
ABSF(1,18,0:10,3) = (/0.000_EB,0.094_EB,0.124_EB,0.123_EB,0.122_EB,0.123_EB,0.125_EB,0.127_EB,0.129_EB,0.130_EB,0.131_EB/)
ABSF(1,18,11:20,3) = (/0.133_EB,0.134_EB,0.132_EB,0.130_EB,0.129_EB,0.127_EB,0.127_EB,0.127_EB,0.127_EB,0.127_EB/)
ABSF(1,18,21:30,3) = (/0.133_EB,0.139_EB,0.145_EB,0.151_EB,0.157_EB,0.154_EB,0.151_EB,0.148_EB,0.145_EB,0.142_EB/)
ABSF(1,18,0:10,4) = (/0.000_EB,0.069_EB,0.105_EB,0.115_EB,0.125_EB,0.126_EB,0.128_EB,0.130_EB,0.131_EB,0.139_EB,0.147_EB/)
ABSF(1,18,11:20,4) = (/0.154_EB,0.162_EB,0.164_EB,0.165_EB,0.167_EB,0.169_EB,0.165_EB,0.162_EB,0.158_EB,0.154_EB/)
ABSF(1,18,21:30,4) = (/0.157_EB,0.159_EB,0.162_EB,0.164_EB,0.167_EB,0.170_EB,0.173_EB,0.176_EB,0.179_EB,0.182_EB/)
ABSF(1,18,0:10,5) = (/0.000_EB,0.072_EB,0.100_EB,0.113_EB,0.126_EB,0.127_EB,0.127_EB,0.128_EB,0.129_EB,0.145_EB,0.162_EB/)
ABSF(1,18,11:20,5) = (/0.179_EB,0.196_EB,0.198_EB,0.200_EB,0.202_EB,0.205_EB,0.196_EB,0.187_EB,0.178_EB,0.169_EB/)
ABSF(1,18,21:30,5) = (/0.167_EB,0.164_EB,0.162_EB,0.159_EB,0.157_EB,0.167_EB,0.177_EB,0.187_EB,0.198_EB,0.208_EB/)
ABSF(1,18,0:10,6) = (/0.000_EB,0.067_EB,0.104_EB,0.108_EB,0.113_EB,0.125_EB,0.137_EB,0.148_EB,0.160_EB,0.170_EB,0.180_EB/)
ABSF(1,18,11:20,6) = (/0.190_EB,0.200_EB,0.205_EB,0.210_EB,0.215_EB,0.220_EB,0.229_EB,0.237_EB,0.245_EB,0.254_EB/)
ABSF(1,18,21:30,6) = (/0.234_EB,0.214_EB,0.194_EB,0.175_EB,0.155_EB,0.181_EB,0.206_EB,0.232_EB,0.257_EB,0.283_EB/)
ABSF(1,19,0:10,1) = (/0.000_EB,0.104_EB,0.121_EB,0.125_EB,0.130_EB,0.130_EB,0.131_EB,0.131_EB,0.131_EB,0.133_EB,0.135_EB/)
ABSF(1,19,11:20,1) = (/0.137_EB,0.138_EB,0.136_EB,0.134_EB,0.133_EB,0.131_EB,0.130_EB,0.130_EB,0.130_EB,0.130_EB/)
ABSF(1,19,21:30,1) = (/0.136_EB,0.142_EB,0.148_EB,0.155_EB,0.161_EB,0.157_EB,0.153_EB,0.148_EB,0.144_EB,0.140_EB/)
ABSF(1,19,0:10,2) = (/0.000_EB,0.104_EB,0.121_EB,0.127_EB,0.133_EB,0.133_EB,0.134_EB,0.134_EB,0.134_EB,0.135_EB,0.137_EB/)
ABSF(1,19,11:20,2) = (/0.138_EB,0.140_EB,0.138_EB,0.135_EB,0.133_EB,0.131_EB,0.131_EB,0.131_EB,0.130_EB,0.130_EB/)
ABSF(1,19,21:30,2) = (/0.136_EB,0.142_EB,0.148_EB,0.154_EB,0.160_EB,0.156_EB,0.153_EB,0.149_EB,0.146_EB,0.143_EB/)
ABSF(1,19,0:10,3) = (/0.000_EB,0.100_EB,0.129_EB,0.128_EB,0.127_EB,0.128_EB,0.130_EB,0.132_EB,0.134_EB,0.135_EB,0.136_EB/)
ABSF(1,19,11:20,3) = (/0.137_EB,0.138_EB,0.137_EB,0.135_EB,0.133_EB,0.131_EB,0.131_EB,0.131_EB,0.131_EB,0.131_EB/)
ABSF(1,19,21:30,3) = (/0.137_EB,0.143_EB,0.149_EB,0.155_EB,0.161_EB,0.157_EB,0.154_EB,0.150_EB,0.147_EB,0.143_EB/)
ABSF(1,19,0:10,4) = (/0.000_EB,0.072_EB,0.110_EB,0.119_EB,0.129_EB,0.130_EB,0.132_EB,0.133_EB,0.135_EB,0.142_EB,0.149_EB/)
ABSF(1,19,11:20,4) = (/0.156_EB,0.163_EB,0.166_EB,0.168_EB,0.170_EB,0.173_EB,0.169_EB,0.165_EB,0.162_EB,0.158_EB/)
ABSF(1,19,21:30,4) = (/0.160_EB,0.163_EB,0.165_EB,0.167_EB,0.170_EB,0.173_EB,0.175_EB,0.178_EB,0.181_EB,0.184_EB/)
ABSF(1,19,0:10,5) = (/0.000_EB,0.075_EB,0.104_EB,0.117_EB,0.130_EB,0.131_EB,0.132_EB,0.133_EB,0.134_EB,0.150_EB,0.167_EB/)
ABSF(1,19,11:20,5) = (/0.183_EB,0.199_EB,0.202_EB,0.204_EB,0.207_EB,0.209_EB,0.200_EB,0.191_EB,0.181_EB,0.172_EB/)
ABSF(1,19,21:30,5) = (/0.170_EB,0.167_EB,0.164_EB,0.162_EB,0.159_EB,0.169_EB,0.179_EB,0.188_EB,0.198_EB,0.208_EB/)
ABSF(1,19,0:10,6) = (/0.000_EB,0.070_EB,0.108_EB,0.112_EB,0.117_EB,0.129_EB,0.141_EB,0.153_EB,0.165_EB,0.174_EB,0.184_EB/)
ABSF(1,19,11:20,6) = (/0.194_EB,0.203_EB,0.208_EB,0.212_EB,0.216_EB,0.221_EB,0.231_EB,0.241_EB,0.251_EB,0.262_EB/)
ABSF(1,19,21:30,6) = (/0.241_EB,0.221_EB,0.201_EB,0.180_EB,0.160_EB,0.186_EB,0.211_EB,0.237_EB,0.262_EB,0.288_EB/)
ABSF(1,20,0:10,1) = (/0.000_EB,0.108_EB,0.125_EB,0.130_EB,0.135_EB,0.135_EB,0.136_EB,0.136_EB,0.136_EB,0.138_EB,0.139_EB/)
ABSF(1,20,11:20,1) = (/0.141_EB,0.143_EB,0.141_EB,0.139_EB,0.137_EB,0.135_EB,0.135_EB,0.134_EB,0.134_EB,0.134_EB/)
ABSF(1,20,21:30,1) = (/0.140_EB,0.146_EB,0.153_EB,0.159_EB,0.165_EB,0.161_EB,0.156_EB,0.152_EB,0.148_EB,0.143_EB/)
ABSF(1,20,0:10,2) = (/0.000_EB,0.108_EB,0.126_EB,0.132_EB,0.139_EB,0.139_EB,0.139_EB,0.139_EB,0.139_EB,0.140_EB,0.141_EB/)
ABSF(1,20,11:20,2) = (/0.143_EB,0.144_EB,0.142_EB,0.140_EB,0.138_EB,0.135_EB,0.135_EB,0.135_EB,0.134_EB,0.134_EB/)
ABSF(1,20,21:30,2) = (/0.140_EB,0.146_EB,0.152_EB,0.158_EB,0.164_EB,0.160_EB,0.157_EB,0.153_EB,0.150_EB,0.146_EB/)
ABSF(1,20,0:10,3) = (/0.000_EB,0.107_EB,0.135_EB,0.133_EB,0.131_EB,0.133_EB,0.135_EB,0.137_EB,0.139_EB,0.140_EB,0.141_EB/)
ABSF(1,20,11:20,3) = (/0.142_EB,0.143_EB,0.141_EB,0.139_EB,0.137_EB,0.135_EB,0.135_EB,0.135_EB,0.135_EB,0.134_EB/)
ABSF(1,20,21:30,3) = (/0.140_EB,0.146_EB,0.153_EB,0.159_EB,0.165_EB,0.161_EB,0.157_EB,0.152_EB,0.148_EB,0.144_EB/)
ABSF(1,20,0:10,4) = (/0.000_EB,0.076_EB,0.115_EB,0.124_EB,0.133_EB,0.134_EB,0.135_EB,0.137_EB,0.138_EB,0.145_EB,0.151_EB/)
ABSF(1,20,11:20,4) = (/0.158_EB,0.164_EB,0.168_EB,0.171_EB,0.174_EB,0.177_EB,0.173_EB,0.169_EB,0.165_EB,0.161_EB/)
ABSF(1,20,21:30,4) = (/0.163_EB,0.166_EB,0.168_EB,0.170_EB,0.173_EB,0.175_EB,0.178_EB,0.180_EB,0.183_EB,0.185_EB/)
ABSF(1,20,0:10,5) = (/0.000_EB,0.078_EB,0.107_EB,0.120_EB,0.133_EB,0.134_EB,0.136_EB,0.138_EB,0.139_EB,0.155_EB,0.171_EB/)
ABSF(1,20,11:20,5) = (/0.187_EB,0.203_EB,0.206_EB,0.208_EB,0.211_EB,0.214_EB,0.204_EB,0.194_EB,0.185_EB,0.175_EB/)
ABSF(1,20,21:30,5) = (/0.172_EB,0.170_EB,0.167_EB,0.165_EB,0.162_EB,0.171_EB,0.180_EB,0.189_EB,0.199_EB,0.208_EB/)
ABSF(1,20,0:10,6) = (/0.000_EB,0.073_EB,0.111_EB,0.116_EB,0.121_EB,0.133_EB,0.145_EB,0.157_EB,0.169_EB,0.179_EB,0.188_EB/)
ABSF(1,20,11:20,6) = (/0.197_EB,0.207_EB,0.210_EB,0.214_EB,0.218_EB,0.222_EB,0.234_EB,0.245_EB,0.257_EB,0.269_EB/)
ABSF(1,20,21:30,6) = (/0.249_EB,0.228_EB,0.207_EB,0.186_EB,0.165_EB,0.191_EB,0.216_EB,0.241_EB,0.267_EB,0.292_EB/)
ABSF(2,0,0:10,1) = (/0.000_EB,0.009_EB,0.021_EB,0.023_EB,0.025_EB,0.026_EB,0.026_EB,0.027_EB,0.028_EB,0.031_EB,0.034_EB/)
ABSF(2,0,11:20,1) = (/0.037_EB,0.040_EB,0.039_EB,0.039_EB,0.038_EB,0.037_EB,0.038_EB,0.040_EB,0.041_EB,0.042_EB/)
ABSF(2,0,21:30,1) = (/0.043_EB,0.044_EB,0.045_EB,0.046_EB,0.047_EB,0.047_EB,0.048_EB,0.048_EB,0.048_EB,0.049_EB/)
ABSF(2,0,0:10,2) = (/0.000_EB,0.013_EB,0.020_EB,0.028_EB,0.036_EB,0.037_EB,0.037_EB,0.038_EB,0.038_EB,0.039_EB,0.039_EB/)
ABSF(2,0,11:20,2) = (/0.040_EB,0.040_EB,0.040_EB,0.040_EB,0.039_EB,0.039_EB,0.043_EB,0.048_EB,0.052_EB,0.056_EB/)
ABSF(2,0,21:30,2) = (/0.051_EB,0.045_EB,0.039_EB,0.034_EB,0.028_EB,0.035_EB,0.041_EB,0.047_EB,0.053_EB,0.059_EB/)
ABSF(2,0,0:10,3) = (/0.000_EB,0.016_EB,0.024_EB,0.032_EB,0.041_EB,0.041_EB,0.042_EB,0.042_EB,0.042_EB,0.043_EB,0.044_EB/)
ABSF(2,0,11:20,3) = (/0.045_EB,0.046_EB,0.044_EB,0.042_EB,0.039_EB,0.037_EB,0.047_EB,0.057_EB,0.067_EB,0.078_EB/)
ABSF(2,0,21:30,3) = (/0.077_EB,0.076_EB,0.076_EB,0.075_EB,0.075_EB,0.074_EB,0.073_EB,0.072_EB,0.071_EB,0.071_EB/)
ABSF(2,0,0:10,4) = (/0.000_EB,0.014_EB,0.025_EB,0.034_EB,0.044_EB,0.048_EB,0.051_EB,0.054_EB,0.058_EB,0.054_EB,0.051_EB/)
ABSF(2,0,11:20,4) = (/0.047_EB,0.043_EB,0.037_EB,0.031_EB,0.025_EB,0.019_EB,0.039_EB,0.059_EB,0.080_EB,0.100_EB/)
ABSF(2,0,21:30,4) = (/0.101_EB,0.103_EB,0.104_EB,0.105_EB,0.106_EB,0.107_EB,0.108_EB,0.109_EB,0.110_EB,0.110_EB/)
ABSF(2,0,0:10,5) = (/0.000_EB,0.018_EB,0.025_EB,0.029_EB,0.034_EB,0.043_EB,0.052_EB,0.061_EB,0.070_EB,0.059_EB,0.049_EB/)
ABSF(2,0,11:20,5) = (/0.038_EB,0.028_EB,0.028_EB,0.028_EB,0.028_EB,0.028_EB,0.048_EB,0.068_EB,0.089_EB,0.109_EB/)
ABSF(2,0,21:30,5) = (/0.114_EB,0.119_EB,0.124_EB,0.129_EB,0.134_EB,0.133_EB,0.131_EB,0.130_EB,0.128_EB,0.127_EB/)
ABSF(2,0,0:10,6) = (/0.000_EB,0.016_EB,0.033_EB,0.025_EB,0.018_EB,0.033_EB,0.049_EB,0.064_EB,0.080_EB,0.065_EB,0.050_EB/)
ABSF(2,0,11:20,6) = (/0.035_EB,0.020_EB,0.024_EB,0.028_EB,0.031_EB,0.035_EB,0.038_EB,0.042_EB,0.045_EB,0.048_EB/)
ABSF(2,0,21:30,6) = (/0.067_EB,0.085_EB,0.103_EB,0.121_EB,0.139_EB,0.140_EB,0.141_EB,0.142_EB,0.143_EB,0.144_EB/)
ABSF(2,1,0:10,1) = (/0.000_EB,0.020_EB,0.025_EB,0.027_EB,0.029_EB,0.028_EB,0.027_EB,0.027_EB,0.026_EB,0.030_EB,0.034_EB/)
ABSF(2,1,11:20,1) = (/0.038_EB,0.043_EB,0.042_EB,0.042_EB,0.042_EB,0.041_EB,0.042_EB,0.044_EB,0.045_EB,0.046_EB/)
ABSF(2,1,21:30,1) = (/0.046_EB,0.047_EB,0.047_EB,0.047_EB,0.048_EB,0.048_EB,0.049_EB,0.049_EB,0.050_EB,0.050_EB/)
ABSF(2,1,0:10,2) = (/0.000_EB,0.013_EB,0.023_EB,0.032_EB,0.041_EB,0.039_EB,0.038_EB,0.037_EB,0.036_EB,0.039_EB,0.041_EB/)
ABSF(2,1,11:20,2) = (/0.044_EB,0.047_EB,0.045_EB,0.044_EB,0.042_EB,0.040_EB,0.045_EB,0.049_EB,0.054_EB,0.059_EB/)
ABSF(2,1,21:30,2) = (/0.054_EB,0.050_EB,0.045_EB,0.041_EB,0.037_EB,0.041_EB,0.046_EB,0.051_EB,0.055_EB,0.060_EB/)
ABSF(2,1,0:10,3) = (/0.000_EB,0.017_EB,0.027_EB,0.036_EB,0.045_EB,0.046_EB,0.047_EB,0.049_EB,0.050_EB,0.048_EB,0.046_EB/)
ABSF(2,1,11:20,3) = (/0.044_EB,0.041_EB,0.042_EB,0.042_EB,0.043_EB,0.043_EB,0.052_EB,0.062_EB,0.071_EB,0.080_EB/)
ABSF(2,1,21:30,3) = (/0.078_EB,0.075_EB,0.072_EB,0.069_EB,0.066_EB,0.070_EB,0.074_EB,0.078_EB,0.082_EB,0.086_EB/)
ABSF(2,1,0:10,4) = (/0.000_EB,0.019_EB,0.028_EB,0.037_EB,0.046_EB,0.048_EB,0.051_EB,0.053_EB,0.056_EB,0.055_EB,0.053_EB/)
ABSF(2,1,11:20,4) = (/0.052_EB,0.050_EB,0.043_EB,0.037_EB,0.030_EB,0.023_EB,0.045_EB,0.068_EB,0.090_EB,0.112_EB/)
ABSF(2,1,21:30,4) = (/0.107_EB,0.102_EB,0.097_EB,0.092_EB,0.087_EB,0.087_EB,0.088_EB,0.089_EB,0.090_EB,0.091_EB/)
ABSF(2,1,0:10,5) = (/0.000_EB,0.019_EB,0.028_EB,0.034_EB,0.039_EB,0.046_EB,0.054_EB,0.061_EB,0.068_EB,0.061_EB,0.053_EB/)
ABSF(2,1,11:20,5) = (/0.046_EB,0.038_EB,0.036_EB,0.034_EB,0.031_EB,0.029_EB,0.049_EB,0.069_EB,0.088_EB,0.108_EB/)
ABSF(2,1,21:30,5) = (/0.113_EB,0.117_EB,0.122_EB,0.126_EB,0.131_EB,0.125_EB,0.120_EB,0.114_EB,0.109_EB,0.103_EB/)
ABSF(2,1,0:10,6) = (/0.000_EB,0.018_EB,0.043_EB,0.039_EB,0.035_EB,0.048_EB,0.060_EB,0.073_EB,0.086_EB,0.070_EB,0.055_EB/)
ABSF(2,1,11:20,6) = (/0.039_EB,0.024_EB,0.032_EB,0.041_EB,0.050_EB,0.059_EB,0.056_EB,0.053_EB,0.050_EB,0.047_EB/)
ABSF(2,1,21:30,6) = (/0.070_EB,0.094_EB,0.117_EB,0.141_EB,0.164_EB,0.169_EB,0.175_EB,0.180_EB,0.185_EB,0.190_EB/)
ABSF(2,2,0:10,1) = (/0.000_EB,0.024_EB,0.029_EB,0.031_EB,0.032_EB,0.031_EB,0.029_EB,0.028_EB,0.026_EB,0.031_EB,0.036_EB/)
ABSF(2,2,11:20,1) = (/0.041_EB,0.045_EB,0.045_EB,0.044_EB,0.043_EB,0.042_EB,0.043_EB,0.045_EB,0.046_EB,0.048_EB/)
ABSF(2,2,21:30,1) = (/0.048_EB,0.048_EB,0.049_EB,0.049_EB,0.049_EB,0.050_EB,0.050_EB,0.051_EB,0.052_EB,0.052_EB/)
ABSF(2,2,0:10,2) = (/0.000_EB,0.016_EB,0.025_EB,0.033_EB,0.041_EB,0.038_EB,0.036_EB,0.033_EB,0.030_EB,0.034_EB,0.037_EB/)
ABSF(2,2,11:20,2) = (/0.041_EB,0.044_EB,0.045_EB,0.045_EB,0.046_EB,0.047_EB,0.049_EB,0.052_EB,0.055_EB,0.058_EB/)
ABSF(2,2,21:30,2) = (/0.056_EB,0.055_EB,0.054_EB,0.053_EB,0.051_EB,0.055_EB,0.058_EB,0.061_EB,0.064_EB,0.068_EB/)
ABSF(2,2,0:10,3) = (/0.000_EB,0.018_EB,0.030_EB,0.039_EB,0.048_EB,0.048_EB,0.049_EB,0.049_EB,0.050_EB,0.050_EB,0.050_EB/)
ABSF(2,2,11:20,3) = (/0.050_EB,0.050_EB,0.048_EB,0.047_EB,0.045_EB,0.043_EB,0.051_EB,0.059_EB,0.067_EB,0.075_EB/)
ABSF(2,2,21:30,3) = (/0.078_EB,0.080_EB,0.083_EB,0.085_EB,0.088_EB,0.085_EB,0.081_EB,0.078_EB,0.075_EB,0.072_EB/)
ABSF(2,2,0:10,4) = (/0.000_EB,0.020_EB,0.031_EB,0.040_EB,0.050_EB,0.052_EB,0.054_EB,0.056_EB,0.058_EB,0.055_EB,0.052_EB/)
ABSF(2,2,11:20,4) = (/0.049_EB,0.047_EB,0.042_EB,0.037_EB,0.032_EB,0.027_EB,0.046_EB,0.065_EB,0.084_EB,0.102_EB/)
ABSF(2,2,21:30,4) = (/0.101_EB,0.100_EB,0.099_EB,0.098_EB,0.097_EB,0.100_EB,0.103_EB,0.107_EB,0.110_EB,0.113_EB/)
ABSF(2,2,0:10,5) = (/0.000_EB,0.020_EB,0.030_EB,0.037_EB,0.043_EB,0.050_EB,0.056_EB,0.063_EB,0.069_EB,0.060_EB,0.051_EB/)
ABSF(2,2,11:20,5) = (/0.041_EB,0.032_EB,0.032_EB,0.032_EB,0.032_EB,0.031_EB,0.052_EB,0.073_EB,0.094_EB,0.115_EB/)
ABSF(2,2,21:30,5) = (/0.117_EB,0.119_EB,0.121_EB,0.122_EB,0.124_EB,0.119_EB,0.113_EB,0.107_EB,0.102_EB,0.096_EB/)
ABSF(2,2,0:10,6) = (/0.000_EB,0.021_EB,0.037_EB,0.040_EB,0.043_EB,0.055_EB,0.068_EB,0.081_EB,0.094_EB,0.077_EB,0.061_EB/)
ABSF(2,2,11:20,6) = (/0.045_EB,0.028_EB,0.032_EB,0.035_EB,0.039_EB,0.043_EB,0.045_EB,0.047_EB,0.049_EB,0.051_EB/)
ABSF(2,2,21:30,6) = (/0.073_EB,0.094_EB,0.115_EB,0.137_EB,0.158_EB,0.156_EB,0.154_EB,0.152_EB,0.150_EB,0.148_EB/)
ABSF(2,3,0:10,1) = (/0.000_EB,0.028_EB,0.033_EB,0.034_EB,0.036_EB,0.034_EB,0.032_EB,0.030_EB,0.028_EB,0.033_EB,0.038_EB/)
ABSF(2,3,11:20,1) = (/0.043_EB,0.048_EB,0.047_EB,0.046_EB,0.045_EB,0.045_EB,0.046_EB,0.048_EB,0.050_EB,0.052_EB/)
ABSF(2,3,21:30,1) = (/0.052_EB,0.052_EB,0.052_EB,0.052_EB,0.052_EB,0.052_EB,0.053_EB,0.053_EB,0.054_EB,0.054_EB/)
ABSF(2,3,0:10,2) = (/0.000_EB,0.020_EB,0.030_EB,0.036_EB,0.042_EB,0.039_EB,0.036_EB,0.033_EB,0.030_EB,0.035_EB,0.039_EB/)
ABSF(2,3,11:20,2) = (/0.043_EB,0.047_EB,0.047_EB,0.047_EB,0.048_EB,0.048_EB,0.050_EB,0.053_EB,0.056_EB,0.058_EB/)
ABSF(2,3,21:30,2) = (/0.057_EB,0.056_EB,0.056_EB,0.055_EB,0.054_EB,0.056_EB,0.058_EB,0.060_EB,0.062_EB,0.064_EB/)
ABSF(2,3,0:10,3) = (/0.000_EB,0.020_EB,0.032_EB,0.042_EB,0.051_EB,0.050_EB,0.050_EB,0.050_EB,0.049_EB,0.050_EB,0.051_EB/)
ABSF(2,3,11:20,3) = (/0.052_EB,0.053_EB,0.052_EB,0.050_EB,0.049_EB,0.048_EB,0.056_EB,0.064_EB,0.072_EB,0.080_EB/)
ABSF(2,3,21:30,3) = (/0.081_EB,0.082_EB,0.084_EB,0.085_EB,0.086_EB,0.083_EB,0.081_EB,0.079_EB,0.077_EB,0.075_EB/)
ABSF(2,3,0:10,4) = (/0.000_EB,0.022_EB,0.034_EB,0.043_EB,0.053_EB,0.055_EB,0.057_EB,0.059_EB,0.062_EB,0.060_EB,0.058_EB/)
ABSF(2,3,11:20,4) = (/0.056_EB,0.054_EB,0.049_EB,0.044_EB,0.039_EB,0.034_EB,0.052_EB,0.071_EB,0.089_EB,0.107_EB/)
ABSF(2,3,21:30,4) = (/0.106_EB,0.105_EB,0.104_EB,0.103_EB,0.102_EB,0.104_EB,0.106_EB,0.108_EB,0.110_EB,0.112_EB/)
ABSF(2,3,0:10,5) = (/0.000_EB,0.022_EB,0.033_EB,0.040_EB,0.047_EB,0.054_EB,0.062_EB,0.069_EB,0.077_EB,0.066_EB,0.056_EB/)
ABSF(2,3,11:20,5) = (/0.045_EB,0.035_EB,0.035_EB,0.036_EB,0.037_EB,0.037_EB,0.057_EB,0.078_EB,0.098_EB,0.118_EB/)
ABSF(2,3,21:30,5) = (/0.120_EB,0.122_EB,0.124_EB,0.126_EB,0.128_EB,0.122_EB,0.117_EB,0.111_EB,0.105_EB,0.099_EB/)
ABSF(2,3,0:10,6) = (/0.000_EB,0.023_EB,0.040_EB,0.044_EB,0.047_EB,0.059_EB,0.070_EB,0.082_EB,0.094_EB,0.079_EB,0.064_EB/)
ABSF(2,3,11:20,6) = (/0.049_EB,0.034_EB,0.036_EB,0.038_EB,0.040_EB,0.042_EB,0.045_EB,0.048_EB,0.051_EB,0.054_EB/)
ABSF(2,3,21:30,6) = (/0.075_EB,0.097_EB,0.118_EB,0.140_EB,0.161_EB,0.157_EB,0.153_EB,0.149_EB,0.146_EB,0.142_EB/)
ABSF(2,4,0:10,1) = (/0.000_EB,0.032_EB,0.036_EB,0.038_EB,0.039_EB,0.037_EB,0.034_EB,0.031_EB,0.029_EB,0.034_EB,0.040_EB/)
ABSF(2,4,11:20,1) = (/0.045_EB,0.050_EB,0.050_EB,0.049_EB,0.048_EB,0.047_EB,0.049_EB,0.051_EB,0.053_EB,0.055_EB/)
ABSF(2,4,21:30,1) = (/0.055_EB,0.055_EB,0.055_EB,0.055_EB,0.055_EB,0.055_EB,0.055_EB,0.055_EB,0.056_EB,0.056_EB/)
ABSF(2,4,0:10,2) = (/0.000_EB,0.025_EB,0.034_EB,0.038_EB,0.042_EB,0.039_EB,0.037_EB,0.034_EB,0.031_EB,0.035_EB,0.040_EB/)
ABSF(2,4,11:20,2) = (/0.045_EB,0.050_EB,0.049_EB,0.049_EB,0.049_EB,0.049_EB,0.051_EB,0.054_EB,0.056_EB,0.059_EB/)
ABSF(2,4,21:30,2) = (/0.058_EB,0.058_EB,0.057_EB,0.057_EB,0.056_EB,0.057_EB,0.058_EB,0.059_EB,0.059_EB,0.060_EB/)
ABSF(2,4,0:10,3) = (/0.000_EB,0.023_EB,0.035_EB,0.044_EB,0.054_EB,0.053_EB,0.051_EB,0.050_EB,0.049_EB,0.051_EB,0.052_EB/)
ABSF(2,4,11:20,3) = (/0.054_EB,0.056_EB,0.055_EB,0.054_EB,0.053_EB,0.052_EB,0.060_EB,0.069_EB,0.077_EB,0.086_EB/)
ABSF(2,4,21:30,3) = (/0.085_EB,0.085_EB,0.084_EB,0.084_EB,0.083_EB,0.082_EB,0.081_EB,0.081_EB,0.080_EB,0.079_EB/)
ABSF(2,4,0:10,4) = (/0.000_EB,0.025_EB,0.037_EB,0.046_EB,0.056_EB,0.058_EB,0.060_EB,0.063_EB,0.065_EB,0.064_EB,0.064_EB/)
ABSF(2,4,11:20,4) = (/0.063_EB,0.062_EB,0.057_EB,0.051_EB,0.046_EB,0.040_EB,0.058_EB,0.076_EB,0.095_EB,0.113_EB/)
ABSF(2,4,21:30,4) = (/0.112_EB,0.111_EB,0.110_EB,0.109_EB,0.108_EB,0.109_EB,0.109_EB,0.110_EB,0.111_EB,0.111_EB/)
ABSF(2,4,0:10,5) = (/0.000_EB,0.024_EB,0.035_EB,0.043_EB,0.051_EB,0.059_EB,0.068_EB,0.076_EB,0.084_EB,0.073_EB,0.061_EB/)
ABSF(2,4,11:20,5) = (/0.049_EB,0.038_EB,0.039_EB,0.040_EB,0.042_EB,0.043_EB,0.062_EB,0.082_EB,0.102_EB,0.121_EB/)
ABSF(2,4,21:30,5) = (/0.123_EB,0.125_EB,0.127_EB,0.130_EB,0.132_EB,0.126_EB,0.120_EB,0.114_EB,0.109_EB,0.103_EB/)
ABSF(2,4,0:10,6) = (/0.000_EB,0.025_EB,0.044_EB,0.047_EB,0.051_EB,0.062_EB,0.072_EB,0.083_EB,0.094_EB,0.080_EB,0.067_EB/)
ABSF(2,4,11:20,6) = (/0.053_EB,0.040_EB,0.040_EB,0.041_EB,0.041_EB,0.042_EB,0.045_EB,0.049_EB,0.053_EB,0.056_EB/)
ABSF(2,4,21:30,6) = (/0.078_EB,0.099_EB,0.121_EB,0.143_EB,0.164_EB,0.158_EB,0.153_EB,0.147_EB,0.141_EB,0.135_EB/)
ABSF(2,5,0:10,1) = (/0.000_EB,0.036_EB,0.040_EB,0.041_EB,0.043_EB,0.040_EB,0.036_EB,0.033_EB,0.030_EB,0.036_EB,0.041_EB/)
ABSF(2,5,11:20,1) = (/0.047_EB,0.053_EB,0.052_EB,0.051_EB,0.050_EB,0.050_EB,0.052_EB,0.054_EB,0.057_EB,0.059_EB/)
ABSF(2,5,21:30,1) = (/0.059_EB,0.058_EB,0.058_EB,0.058_EB,0.058_EB,0.058_EB,0.058_EB,0.057_EB,0.057_EB,0.057_EB/)
ABSF(2,5,0:10,2) = (/0.000_EB,0.029_EB,0.039_EB,0.041_EB,0.043_EB,0.040_EB,0.037_EB,0.034_EB,0.031_EB,0.036_EB,0.042_EB/)
ABSF(2,5,11:20,2) = (/0.047_EB,0.052_EB,0.052_EB,0.051_EB,0.051_EB,0.050_EB,0.052_EB,0.055_EB,0.057_EB,0.059_EB/)
ABSF(2,5,21:30,2) = (/0.059_EB,0.059_EB,0.059_EB,0.059_EB,0.059_EB,0.058_EB,0.058_EB,0.057_EB,0.057_EB,0.056_EB/)
ABSF(2,5,0:10,3) = (/0.000_EB,0.026_EB,0.038_EB,0.047_EB,0.057_EB,0.055_EB,0.053_EB,0.050_EB,0.048_EB,0.051_EB,0.054_EB/)
ABSF(2,5,11:20,3) = (/0.056_EB,0.059_EB,0.058_EB,0.057_EB,0.056_EB,0.056_EB,0.064_EB,0.073_EB,0.082_EB,0.091_EB/)
ABSF(2,5,21:30,3) = (/0.089_EB,0.087_EB,0.085_EB,0.083_EB,0.081_EB,0.081_EB,0.082_EB,0.082_EB,0.082_EB,0.082_EB/)
ABSF(2,5,0:10,4) = (/0.000_EB,0.027_EB,0.040_EB,0.049_EB,0.059_EB,0.061_EB,0.064_EB,0.066_EB,0.069_EB,0.069_EB,0.069_EB/)
ABSF(2,5,11:20,4) = (/0.070_EB,0.070_EB,0.064_EB,0.058_EB,0.053_EB,0.047_EB,0.064_EB,0.082_EB,0.100_EB,0.118_EB/)
ABSF(2,5,21:30,4) = (/0.117_EB,0.116_EB,0.115_EB,0.114_EB,0.114_EB,0.113_EB,0.112_EB,0.112_EB,0.111_EB,0.110_EB/)
ABSF(2,5,0:10,5) = (/0.000_EB,0.026_EB,0.038_EB,0.046_EB,0.055_EB,0.064_EB,0.073_EB,0.083_EB,0.092_EB,0.079_EB,0.066_EB/)
ABSF(2,5,11:20,5) = (/0.053_EB,0.041_EB,0.043_EB,0.045_EB,0.047_EB,0.049_EB,0.067_EB,0.086_EB,0.105_EB,0.124_EB/)
ABSF(2,5,21:30,5) = (/0.126_EB,0.129_EB,0.131_EB,0.133_EB,0.136_EB,0.130_EB,0.124_EB,0.118_EB,0.112_EB,0.106_EB/)
ABSF(2,5,0:10,6) = (/0.000_EB,0.027_EB,0.047_EB,0.051_EB,0.055_EB,0.065_EB,0.075_EB,0.084_EB,0.094_EB,0.082_EB,0.070_EB/)
ABSF(2,5,11:20,6) = (/0.058_EB,0.046_EB,0.045_EB,0.044_EB,0.042_EB,0.041_EB,0.046_EB,0.050_EB,0.054_EB,0.059_EB/)
ABSF(2,5,21:30,6) = (/0.080_EB,0.102_EB,0.124_EB,0.146_EB,0.167_EB,0.160_EB,0.152_EB,0.144_EB,0.136_EB,0.129_EB/)
ABSF(2,6,0:10,1) = (/0.000_EB,0.039_EB,0.043_EB,0.045_EB,0.046_EB,0.043_EB,0.040_EB,0.036_EB,0.033_EB,0.039_EB,0.044_EB/)
ABSF(2,6,11:20,1) = (/0.050_EB,0.056_EB,0.055_EB,0.055_EB,0.054_EB,0.053_EB,0.055_EB,0.058_EB,0.060_EB,0.063_EB/)
ABSF(2,6,21:30,1) = (/0.062_EB,0.062_EB,0.062_EB,0.062_EB,0.061_EB,0.061_EB,0.061_EB,0.061_EB,0.061_EB,0.061_EB/)
ABSF(2,6,0:10,2) = (/0.000_EB,0.033_EB,0.043_EB,0.045_EB,0.046_EB,0.043_EB,0.040_EB,0.037_EB,0.034_EB,0.039_EB,0.045_EB/)
ABSF(2,6,11:20,2) = (/0.050_EB,0.055_EB,0.055_EB,0.054_EB,0.054_EB,0.053_EB,0.056_EB,0.058_EB,0.060_EB,0.063_EB/)
ABSF(2,6,21:30,2) = (/0.063_EB,0.063_EB,0.063_EB,0.062_EB,0.062_EB,0.062_EB,0.061_EB,0.061_EB,0.060_EB,0.060_EB/)
ABSF(2,6,0:10,3) = (/0.000_EB,0.028_EB,0.040_EB,0.049_EB,0.058_EB,0.056_EB,0.053_EB,0.051_EB,0.048_EB,0.051_EB,0.054_EB/)
ABSF(2,6,11:20,3) = (/0.058_EB,0.061_EB,0.060_EB,0.059_EB,0.059_EB,0.058_EB,0.066_EB,0.073_EB,0.081_EB,0.089_EB/)
ABSF(2,6,21:30,3) = (/0.087_EB,0.086_EB,0.084_EB,0.083_EB,0.082_EB,0.082_EB,0.082_EB,0.082_EB,0.082_EB,0.082_EB/)
ABSF(2,6,0:10,4) = (/0.000_EB,0.029_EB,0.042_EB,0.052_EB,0.061_EB,0.064_EB,0.067_EB,0.070_EB,0.072_EB,0.073_EB,0.073_EB/)
ABSF(2,6,11:20,4) = (/0.073_EB,0.073_EB,0.068_EB,0.063_EB,0.058_EB,0.054_EB,0.070_EB,0.086_EB,0.103_EB,0.119_EB/)
ABSF(2,6,21:30,4) = (/0.118_EB,0.118_EB,0.117_EB,0.117_EB,0.116_EB,0.115_EB,0.115_EB,0.114_EB,0.114_EB,0.113_EB/)
ABSF(2,6,0:10,5) = (/0.000_EB,0.027_EB,0.041_EB,0.049_EB,0.057_EB,0.066_EB,0.075_EB,0.083_EB,0.092_EB,0.080_EB,0.068_EB/)
ABSF(2,6,11:20,5) = (/0.055_EB,0.043_EB,0.045_EB,0.047_EB,0.048_EB,0.050_EB,0.070_EB,0.089_EB,0.109_EB,0.128_EB/)
ABSF(2,6,21:30,5) = (/0.130_EB,0.133_EB,0.135_EB,0.137_EB,0.140_EB,0.134_EB,0.128_EB,0.123_EB,0.117_EB,0.111_EB/)
ABSF(2,6,0:10,6) = (/0.000_EB,0.029_EB,0.051_EB,0.055_EB,0.060_EB,0.070_EB,0.080_EB,0.090_EB,0.100_EB,0.086_EB,0.072_EB/)
ABSF(2,6,11:20,6) = (/0.059_EB,0.045_EB,0.044_EB,0.044_EB,0.043_EB,0.043_EB,0.049_EB,0.056_EB,0.062_EB,0.069_EB/)
ABSF(2,6,21:30,6) = (/0.088_EB,0.108_EB,0.127_EB,0.146_EB,0.166_EB,0.160_EB,0.154_EB,0.149_EB,0.143_EB,0.137_EB/)
ABSF(2,7,0:10,1) = (/0.000_EB,0.042_EB,0.047_EB,0.048_EB,0.049_EB,0.046_EB,0.043_EB,0.039_EB,0.036_EB,0.042_EB,0.048_EB/)
ABSF(2,7,11:20,1) = (/0.053_EB,0.059_EB,0.059_EB,0.058_EB,0.057_EB,0.056_EB,0.059_EB,0.061_EB,0.064_EB,0.067_EB/)
ABSF(2,7,21:30,1) = (/0.066_EB,0.066_EB,0.066_EB,0.065_EB,0.065_EB,0.065_EB,0.065_EB,0.064_EB,0.064_EB,0.064_EB/)
ABSF(2,7,0:10,2) = (/0.000_EB,0.037_EB,0.046_EB,0.048_EB,0.050_EB,0.046_EB,0.043_EB,0.040_EB,0.037_EB,0.042_EB,0.048_EB/)
ABSF(2,7,11:20,2) = (/0.053_EB,0.059_EB,0.058_EB,0.058_EB,0.057_EB,0.056_EB,0.059_EB,0.062_EB,0.064_EB,0.067_EB/)
ABSF(2,7,21:30,2) = (/0.067_EB,0.066_EB,0.066_EB,0.066_EB,0.066_EB,0.065_EB,0.065_EB,0.064_EB,0.064_EB,0.063_EB/)
ABSF(2,7,0:10,3) = (/0.000_EB,0.029_EB,0.043_EB,0.051_EB,0.060_EB,0.057_EB,0.054_EB,0.051_EB,0.048_EB,0.052_EB,0.055_EB/)
ABSF(2,7,11:20,3) = (/0.059_EB,0.062_EB,0.062_EB,0.061_EB,0.061_EB,0.060_EB,0.067_EB,0.074_EB,0.080_EB,0.087_EB/)
ABSF(2,7,21:30,3) = (/0.086_EB,0.085_EB,0.084_EB,0.083_EB,0.082_EB,0.082_EB,0.082_EB,0.082_EB,0.082_EB,0.082_EB/)
ABSF(2,7,0:10,4) = (/0.000_EB,0.031_EB,0.044_EB,0.054_EB,0.064_EB,0.067_EB,0.070_EB,0.073_EB,0.076_EB,0.076_EB,0.076_EB/)
ABSF(2,7,11:20,4) = (/0.076_EB,0.076_EB,0.072_EB,0.068_EB,0.064_EB,0.060_EB,0.075_EB,0.090_EB,0.105_EB,0.120_EB/)
ABSF(2,7,21:30,4) = (/0.120_EB,0.119_EB,0.119_EB,0.119_EB,0.118_EB,0.118_EB,0.117_EB,0.117_EB,0.117_EB,0.116_EB/)
ABSF(2,7,0:10,5) = (/0.000_EB,0.028_EB,0.044_EB,0.052_EB,0.060_EB,0.068_EB,0.076_EB,0.084_EB,0.092_EB,0.080_EB,0.069_EB/)
ABSF(2,7,11:20,5) = (/0.057_EB,0.046_EB,0.047_EB,0.049_EB,0.050_EB,0.051_EB,0.072_EB,0.092_EB,0.112_EB,0.132_EB/)
ABSF(2,7,21:30,5) = (/0.135_EB,0.137_EB,0.139_EB,0.142_EB,0.144_EB,0.138_EB,0.133_EB,0.127_EB,0.122_EB,0.116_EB/)
ABSF(2,7,0:10,6) = (/0.000_EB,0.030_EB,0.055_EB,0.060_EB,0.064_EB,0.075_EB,0.085_EB,0.095_EB,0.106_EB,0.090_EB,0.075_EB/)
ABSF(2,7,11:20,6) = (/0.059_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.053_EB,0.062_EB,0.071_EB,0.079_EB/)
ABSF(2,7,21:30,6) = (/0.096_EB,0.113_EB,0.130_EB,0.147_EB,0.164_EB,0.160_EB,0.157_EB,0.153_EB,0.149_EB,0.146_EB/)
ABSF(2,8,0:10,1) = (/0.000_EB,0.045_EB,0.050_EB,0.051_EB,0.053_EB,0.049_EB,0.046_EB,0.042_EB,0.039_EB,0.045_EB,0.051_EB/)
ABSF(2,8,11:20,1) = (/0.057_EB,0.063_EB,0.062_EB,0.061_EB,0.060_EB,0.060_EB,0.062_EB,0.065_EB,0.068_EB,0.070_EB/)
ABSF(2,8,21:30,1) = (/0.070_EB,0.070_EB,0.069_EB,0.069_EB,0.069_EB,0.068_EB,0.068_EB,0.068_EB,0.068_EB,0.067_EB/)
ABSF(2,8,0:10,2) = (/0.000_EB,0.042_EB,0.050_EB,0.051_EB,0.053_EB,0.050_EB,0.046_EB,0.043_EB,0.040_EB,0.045_EB,0.051_EB/)
ABSF(2,8,11:20,2) = (/0.056_EB,0.062_EB,0.061_EB,0.061_EB,0.060_EB,0.060_EB,0.062_EB,0.065_EB,0.068_EB,0.071_EB/)
ABSF(2,8,21:30,2) = (/0.070_EB,0.070_EB,0.070_EB,0.070_EB,0.070_EB,0.069_EB,0.068_EB,0.068_EB,0.067_EB,0.067_EB/)
ABSF(2,8,0:10,3) = (/0.000_EB,0.031_EB,0.045_EB,0.053_EB,0.061_EB,0.058_EB,0.055_EB,0.051_EB,0.048_EB,0.052_EB,0.056_EB/)
ABSF(2,8,11:20,3) = (/0.060_EB,0.064_EB,0.064_EB,0.063_EB,0.063_EB,0.062_EB,0.068_EB,0.074_EB,0.080_EB,0.085_EB/)
ABSF(2,8,21:30,3) = (/0.085_EB,0.084_EB,0.084_EB,0.083_EB,0.082_EB,0.082_EB,0.082_EB,0.082_EB,0.082_EB,0.082_EB/)
ABSF(2,8,0:10,4) = (/0.000_EB,0.032_EB,0.047_EB,0.057_EB,0.067_EB,0.070_EB,0.073_EB,0.077_EB,0.080_EB,0.080_EB,0.079_EB/)
ABSF(2,8,11:20,4) = (/0.079_EB,0.079_EB,0.076_EB,0.073_EB,0.070_EB,0.067_EB,0.081_EB,0.094_EB,0.108_EB,0.121_EB/)
ABSF(2,8,21:30,4) = (/0.121_EB,0.121_EB,0.121_EB,0.121_EB,0.121_EB,0.120_EB,0.120_EB,0.120_EB,0.120_EB,0.119_EB/)
ABSF(2,8,0:10,5) = (/0.000_EB,0.030_EB,0.047_EB,0.055_EB,0.063_EB,0.070_EB,0.077_EB,0.085_EB,0.092_EB,0.081_EB,0.070_EB/)
ABSF(2,8,11:20,5) = (/0.059_EB,0.049_EB,0.050_EB,0.051_EB,0.052_EB,0.053_EB,0.074_EB,0.095_EB,0.116_EB,0.137_EB/)
ABSF(2,8,21:30,5) = (/0.139_EB,0.141_EB,0.143_EB,0.146_EB,0.148_EB,0.143_EB,0.137_EB,0.132_EB,0.127_EB,0.121_EB/)
ABSF(2,8,0:10,6) = (/0.000_EB,0.032_EB,0.060_EB,0.064_EB,0.069_EB,0.080_EB,0.090_EB,0.101_EB,0.112_EB,0.095_EB,0.077_EB/)
ABSF(2,8,11:20,6) = (/0.060_EB,0.043_EB,0.044_EB,0.044_EB,0.045_EB,0.046_EB,0.057_EB,0.068_EB,0.079_EB,0.090_EB/)
ABSF(2,8,21:30,6) = (/0.104_EB,0.119_EB,0.133_EB,0.148_EB,0.162_EB,0.161_EB,0.159_EB,0.157_EB,0.156_EB,0.154_EB/)
ABSF(2,9,0:10,1) = (/0.000_EB,0.048_EB,0.053_EB,0.055_EB,0.056_EB,0.052_EB,0.049_EB,0.046_EB,0.042_EB,0.048_EB,0.054_EB/)
ABSF(2,9,11:20,1) = (/0.060_EB,0.066_EB,0.065_EB,0.064_EB,0.064_EB,0.063_EB,0.066_EB,0.069_EB,0.071_EB,0.074_EB/)
ABSF(2,9,21:30,1) = (/0.074_EB,0.074_EB,0.073_EB,0.073_EB,0.072_EB,0.072_EB,0.072_EB,0.071_EB,0.071_EB,0.071_EB/)
ABSF(2,9,0:10,2) = (/0.000_EB,0.046_EB,0.054_EB,0.055_EB,0.056_EB,0.053_EB,0.049_EB,0.046_EB,0.043_EB,0.048_EB,0.054_EB/)
ABSF(2,9,11:20,2) = (/0.059_EB,0.065_EB,0.064_EB,0.064_EB,0.063_EB,0.063_EB,0.066_EB,0.069_EB,0.072_EB,0.075_EB/)
ABSF(2,9,21:30,2) = (/0.074_EB,0.074_EB,0.074_EB,0.074_EB,0.073_EB,0.073_EB,0.072_EB,0.071_EB,0.071_EB,0.070_EB/)
ABSF(2,9,0:10,3) = (/0.000_EB,0.032_EB,0.048_EB,0.055_EB,0.063_EB,0.059_EB,0.055_EB,0.052_EB,0.048_EB,0.052_EB,0.057_EB/)
ABSF(2,9,11:20,3) = (/0.061_EB,0.066_EB,0.066_EB,0.065_EB,0.065_EB,0.065_EB,0.069_EB,0.074_EB,0.079_EB,0.084_EB/)
ABSF(2,9,21:30,3) = (/0.084_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.082_EB,0.082_EB/)
ABSF(2,9,0:10,4) = (/0.000_EB,0.034_EB,0.049_EB,0.059_EB,0.069_EB,0.073_EB,0.077_EB,0.080_EB,0.084_EB,0.083_EB,0.083_EB/)
ABSF(2,9,11:20,4) = (/0.082_EB,0.082_EB,0.080_EB,0.078_EB,0.076_EB,0.074_EB,0.086_EB,0.098_EB,0.110_EB,0.122_EB/)
ABSF(2,9,21:30,4) = (/0.123_EB,0.123_EB,0.123_EB,0.123_EB,0.123_EB,0.123_EB,0.123_EB,0.123_EB,0.122_EB,0.122_EB/)
ABSF(2,9,0:10,5) = (/0.000_EB,0.031_EB,0.049_EB,0.057_EB,0.065_EB,0.072_EB,0.079_EB,0.085_EB,0.092_EB,0.082_EB,0.072_EB/)
ABSF(2,9,11:20,5) = (/0.061_EB,0.051_EB,0.052_EB,0.053_EB,0.054_EB,0.054_EB,0.076_EB,0.098_EB,0.119_EB,0.141_EB/)
ABSF(2,9,21:30,5) = (/0.143_EB,0.145_EB,0.148_EB,0.150_EB,0.152_EB,0.147_EB,0.142_EB,0.137_EB,0.131_EB,0.126_EB/)
ABSF(2,9,0:10,6) = (/0.000_EB,0.033_EB,0.064_EB,0.069_EB,0.073_EB,0.084_EB,0.096_EB,0.107_EB,0.118_EB,0.099_EB,0.080_EB/)
ABSF(2,9,11:20,6) = (/0.061_EB,0.042_EB,0.043_EB,0.045_EB,0.046_EB,0.047_EB,0.060_EB,0.074_EB,0.087_EB,0.100_EB/)
ABSF(2,9,21:30,6) = (/0.112_EB,0.124_EB,0.136_EB,0.148_EB,0.160_EB,0.161_EB,0.161_EB,0.162_EB,0.162_EB,0.163_EB/)
ABSF(2,10,0:10,1) = (/0.000_EB,0.051_EB,0.057_EB,0.058_EB,0.059_EB,0.056_EB,0.052_EB,0.049_EB,0.045_EB,0.051_EB,0.057_EB/)
ABSF(2,10,11:20,1) = (/0.063_EB,0.069_EB,0.068_EB,0.068_EB,0.067_EB,0.066_EB,0.069_EB,0.072_EB,0.075_EB,0.078_EB/)
ABSF(2,10,21:30,1) = (/0.078_EB,0.077_EB,0.077_EB,0.077_EB,0.076_EB,0.076_EB,0.075_EB,0.075_EB,0.074_EB,0.074_EB/)
ABSF(2,10,0:10,2) = (/0.000_EB,0.050_EB,0.057_EB,0.058_EB,0.059_EB,0.056_EB,0.053_EB,0.049_EB,0.046_EB,0.051_EB,0.057_EB/)
ABSF(2,10,11:20,2) = (/0.062_EB,0.068_EB,0.068_EB,0.067_EB,0.067_EB,0.066_EB,0.069_EB,0.072_EB,0.075_EB,0.078_EB/)
ABSF(2,10,21:30,2) = (/0.078_EB,0.078_EB,0.078_EB,0.077_EB,0.077_EB,0.076_EB,0.076_EB,0.075_EB,0.074_EB,0.073_EB/)
ABSF(2,10,0:10,3) = (/0.000_EB,0.034_EB,0.051_EB,0.057_EB,0.064_EB,0.060_EB,0.056_EB,0.052_EB,0.048_EB,0.053_EB,0.058_EB/)
ABSF(2,10,11:20,3) = (/0.063_EB,0.068_EB,0.067_EB,0.067_EB,0.067_EB,0.067_EB,0.071_EB,0.074_EB,0.078_EB,0.082_EB/)
ABSF(2,10,21:30,3) = (/0.082_EB,0.082_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.082_EB/)
ABSF(2,10,0:10,4) = (/0.000_EB,0.036_EB,0.052_EB,0.062_EB,0.072_EB,0.076_EB,0.080_EB,0.084_EB,0.087_EB,0.087_EB,0.086_EB/)
ABSF(2,10,11:20,4) = (/0.086_EB,0.085_EB,0.084_EB,0.083_EB,0.082_EB,0.081_EB,0.092_EB,0.102_EB,0.113_EB,0.124_EB/)
ABSF(2,10,21:30,4) = (/0.124_EB,0.124_EB,0.125_EB,0.125_EB,0.125_EB,0.125_EB,0.125_EB,0.125_EB,0.125_EB,0.125_EB/)
ABSF(2,10,0:10,5) = (/0.000_EB,0.033_EB,0.052_EB,0.060_EB,0.068_EB,0.074_EB,0.080_EB,0.086_EB,0.092_EB,0.082_EB,0.073_EB/)
ABSF(2,10,11:20,5) = (/0.063_EB,0.054_EB,0.055_EB,0.055_EB,0.055_EB,0.056_EB,0.078_EB,0.100_EB,0.123_EB,0.145_EB/)
ABSF(2,10,21:30,5) = (/0.147_EB,0.150_EB,0.152_EB,0.154_EB,0.156_EB,0.151_EB,0.146_EB,0.141_EB,0.136_EB,0.131_EB/)
ABSF(2,10,0:10,6) = (/0.000_EB,0.034_EB,0.068_EB,0.073_EB,0.078_EB,0.089_EB,0.101_EB,0.112_EB,0.124_EB,0.103_EB,0.082_EB/)
ABSF(2,10,11:20,6) = (/0.062_EB,0.041_EB,0.043_EB,0.045_EB,0.047_EB,0.049_EB,0.064_EB,0.079_EB,0.095_EB,0.110_EB/)
ABSF(2,10,21:30,6) = (/0.120_EB,0.130_EB,0.139_EB,0.149_EB,0.159_EB,0.161_EB,0.164_EB,0.166_EB,0.169_EB,0.171_EB/)
ABSF(2,11,0:10,1) = (/0.000_EB,0.053_EB,0.060_EB,0.061_EB,0.063_EB,0.059_EB,0.055_EB,0.052_EB,0.048_EB,0.054_EB,0.060_EB/)
ABSF(2,11,11:20,1) = (/0.066_EB,0.072_EB,0.071_EB,0.071_EB,0.070_EB,0.070_EB,0.073_EB,0.075_EB,0.078_EB,0.081_EB/)
ABSF(2,11,21:30,1) = (/0.081_EB,0.080_EB,0.080_EB,0.080_EB,0.079_EB,0.079_EB,0.078_EB,0.078_EB,0.077_EB,0.077_EB/)
ABSF(2,11,0:10,2) = (/0.000_EB,0.053_EB,0.061_EB,0.062_EB,0.063_EB,0.059_EB,0.056_EB,0.052_EB,0.049_EB,0.054_EB,0.060_EB/)
ABSF(2,11,11:20,2) = (/0.066_EB,0.071_EB,0.071_EB,0.070_EB,0.070_EB,0.069_EB,0.072_EB,0.075_EB,0.078_EB,0.081_EB/)
ABSF(2,11,21:30,2) = (/0.081_EB,0.081_EB,0.081_EB,0.080_EB,0.080_EB,0.079_EB,0.079_EB,0.078_EB,0.077_EB,0.076_EB/)
ABSF(2,11,0:10,3) = (/0.000_EB,0.037_EB,0.055_EB,0.061_EB,0.067_EB,0.063_EB,0.059_EB,0.055_EB,0.050_EB,0.055_EB,0.061_EB/)
ABSF(2,11,11:20,3) = (/0.066_EB,0.071_EB,0.071_EB,0.070_EB,0.070_EB,0.070_EB,0.074_EB,0.077_EB,0.081_EB,0.085_EB/)
ABSF(2,11,21:30,3) = (/0.085_EB,0.085_EB,0.085_EB,0.086_EB,0.086_EB,0.085_EB,0.085_EB,0.085_EB,0.085_EB,0.084_EB/)
ABSF(2,11,0:10,4) = (/0.000_EB,0.037_EB,0.054_EB,0.065_EB,0.075_EB,0.078_EB,0.081_EB,0.084_EB,0.087_EB,0.087_EB,0.087_EB/)
ABSF(2,11,11:20,4) = (/0.087_EB,0.086_EB,0.086_EB,0.085_EB,0.084_EB,0.083_EB,0.094_EB,0.104_EB,0.115_EB,0.125_EB/)
ABSF(2,11,21:30,4) = (/0.126_EB,0.126_EB,0.126_EB,0.126_EB,0.127_EB,0.126_EB,0.126_EB,0.126_EB,0.126_EB,0.125_EB/)
ABSF(2,11,0:10,5) = (/0.000_EB,0.034_EB,0.055_EB,0.063_EB,0.071_EB,0.077_EB,0.083_EB,0.090_EB,0.096_EB,0.086_EB,0.076_EB/)
ABSF(2,11,11:20,5) = (/0.067_EB,0.057_EB,0.058_EB,0.059_EB,0.059_EB,0.060_EB,0.082_EB,0.104_EB,0.125_EB,0.147_EB/)
ABSF(2,11,21:30,5) = (/0.149_EB,0.151_EB,0.153_EB,0.155_EB,0.157_EB,0.153_EB,0.148_EB,0.144_EB,0.139_EB,0.135_EB/)
ABSF(2,11,0:10,6) = (/0.000_EB,0.036_EB,0.070_EB,0.076_EB,0.081_EB,0.092_EB,0.104_EB,0.115_EB,0.126_EB,0.106_EB,0.086_EB/)
ABSF(2,11,11:20,6) = (/0.066_EB,0.046_EB,0.047_EB,0.049_EB,0.050_EB,0.052_EB,0.068_EB,0.084_EB,0.100_EB,0.116_EB/)
ABSF(2,11,21:30,6) = (/0.125_EB,0.134_EB,0.143_EB,0.153_EB,0.162_EB,0.164_EB,0.166_EB,0.168_EB,0.170_EB,0.173_EB/)
ABSF(2,12,0:10,1) = (/0.000_EB,0.056_EB,0.064_EB,0.065_EB,0.066_EB,0.062_EB,0.058_EB,0.055_EB,0.051_EB,0.057_EB,0.063_EB/)
ABSF(2,12,11:20,1) = (/0.069_EB,0.075_EB,0.074_EB,0.074_EB,0.074_EB,0.073_EB,0.076_EB,0.079_EB,0.081_EB,0.084_EB/)
ABSF(2,12,21:30,1) = (/0.084_EB,0.084_EB,0.083_EB,0.083_EB,0.083_EB,0.082_EB,0.082_EB,0.081_EB,0.080_EB,0.080_EB/)
ABSF(2,12,0:10,2) = (/0.000_EB,0.055_EB,0.064_EB,0.065_EB,0.066_EB,0.062_EB,0.059_EB,0.055_EB,0.051_EB,0.057_EB,0.063_EB/)
ABSF(2,12,11:20,2) = (/0.069_EB,0.074_EB,0.074_EB,0.074_EB,0.073_EB,0.073_EB,0.076_EB,0.079_EB,0.082_EB,0.085_EB/)
ABSF(2,12,21:30,2) = (/0.084_EB,0.084_EB,0.084_EB,0.083_EB,0.083_EB,0.082_EB,0.082_EB,0.081_EB,0.080_EB,0.079_EB/)
ABSF(2,12,0:10,3) = (/0.000_EB,0.041_EB,0.059_EB,0.064_EB,0.070_EB,0.066_EB,0.062_EB,0.057_EB,0.053_EB,0.058_EB,0.063_EB/)
ABSF(2,12,11:20,3) = (/0.069_EB,0.074_EB,0.074_EB,0.074_EB,0.073_EB,0.073_EB,0.077_EB,0.080_EB,0.084_EB,0.087_EB/)
ABSF(2,12,21:30,3) = (/0.087_EB,0.088_EB,0.088_EB,0.088_EB,0.088_EB,0.088_EB,0.087_EB,0.087_EB,0.087_EB,0.086_EB/)
ABSF(2,12,0:10,4) = (/0.000_EB,0.039_EB,0.057_EB,0.068_EB,0.078_EB,0.080_EB,0.083_EB,0.085_EB,0.087_EB,0.087_EB,0.087_EB/)
ABSF(2,12,11:20,4) = (/0.087_EB,0.088_EB,0.087_EB,0.087_EB,0.086_EB,0.086_EB,0.096_EB,0.107_EB,0.117_EB,0.127_EB/)
ABSF(2,12,21:30,4) = (/0.127_EB,0.128_EB,0.128_EB,0.128_EB,0.128_EB,0.127_EB,0.127_EB,0.126_EB,0.126_EB,0.125_EB/)
ABSF(2,12,0:10,5) = (/0.000_EB,0.036_EB,0.058_EB,0.066_EB,0.074_EB,0.080_EB,0.087_EB,0.093_EB,0.100_EB,0.090_EB,0.080_EB/)
ABSF(2,12,11:20,5) = (/0.070_EB,0.060_EB,0.061_EB,0.062_EB,0.064_EB,0.065_EB,0.086_EB,0.107_EB,0.128_EB,0.149_EB/)
ABSF(2,12,21:30,5) = (/0.151_EB,0.153_EB,0.155_EB,0.157_EB,0.159_EB,0.154_EB,0.150_EB,0.146_EB,0.142_EB,0.138_EB/)
ABSF(2,12,0:10,6) = (/0.000_EB,0.038_EB,0.072_EB,0.078_EB,0.084_EB,0.095_EB,0.107_EB,0.118_EB,0.129_EB,0.109_EB,0.090_EB/)
ABSF(2,12,11:20,6) = (/0.070_EB,0.050_EB,0.052_EB,0.053_EB,0.054_EB,0.056_EB,0.072_EB,0.089_EB,0.105_EB,0.122_EB/)
ABSF(2,12,21:30,6) = (/0.130_EB,0.139_EB,0.148_EB,0.156_EB,0.165_EB,0.167_EB,0.168_EB,0.170_EB,0.172_EB,0.174_EB/)
ABSF(2,13,0:10,1) = (/0.000_EB,0.058_EB,0.067_EB,0.068_EB,0.069_EB,0.065_EB,0.062_EB,0.058_EB,0.054_EB,0.060_EB,0.066_EB/)
ABSF(2,13,11:20,1) = (/0.072_EB,0.078_EB,0.077_EB,0.077_EB,0.077_EB,0.076_EB,0.079_EB,0.082_EB,0.085_EB,0.087_EB/)
ABSF(2,13,21:30,1) = (/0.087_EB,0.087_EB,0.086_EB,0.086_EB,0.086_EB,0.085_EB,0.085_EB,0.084_EB,0.084_EB,0.083_EB/)
ABSF(2,13,0:10,2) = (/0.000_EB,0.058_EB,0.067_EB,0.069_EB,0.070_EB,0.066_EB,0.062_EB,0.058_EB,0.054_EB,0.060_EB,0.066_EB/)
ABSF(2,13,11:20,2) = (/0.072_EB,0.078_EB,0.077_EB,0.077_EB,0.077_EB,0.076_EB,0.079_EB,0.082_EB,0.085_EB,0.088_EB/)
ABSF(2,13,21:30,2) = (/0.087_EB,0.087_EB,0.087_EB,0.086_EB,0.086_EB,0.085_EB,0.084_EB,0.084_EB,0.083_EB,0.082_EB/)
ABSF(2,13,0:10,3) = (/0.000_EB,0.044_EB,0.063_EB,0.068_EB,0.073_EB,0.069_EB,0.065_EB,0.060_EB,0.056_EB,0.061_EB,0.066_EB/)
ABSF(2,13,11:20,3) = (/0.072_EB,0.077_EB,0.077_EB,0.077_EB,0.077_EB,0.077_EB,0.080_EB,0.083_EB,0.087_EB,0.090_EB/)
ABSF(2,13,21:30,3) = (/0.090_EB,0.090_EB,0.090_EB,0.091_EB,0.091_EB,0.090_EB,0.090_EB,0.089_EB,0.089_EB,0.088_EB/)
ABSF(2,13,0:10,4) = (/0.000_EB,0.041_EB,0.060_EB,0.070_EB,0.081_EB,0.082_EB,0.084_EB,0.085_EB,0.087_EB,0.087_EB,0.088_EB/)
ABSF(2,13,11:20,4) = (/0.088_EB,0.089_EB,0.089_EB,0.089_EB,0.088_EB,0.088_EB,0.099_EB,0.109_EB,0.119_EB,0.129_EB/)
ABSF(2,13,21:30,4) = (/0.129_EB,0.129_EB,0.129_EB,0.129_EB,0.129_EB,0.128_EB,0.128_EB,0.127_EB,0.126_EB,0.126_EB/)
ABSF(2,13,0:10,5) = (/0.000_EB,0.038_EB,0.062_EB,0.069_EB,0.077_EB,0.084_EB,0.090_EB,0.097_EB,0.104_EB,0.094_EB,0.084_EB/)
ABSF(2,13,11:20,5) = (/0.073_EB,0.063_EB,0.065_EB,0.066_EB,0.068_EB,0.069_EB,0.089_EB,0.110_EB,0.130_EB,0.150_EB/)
ABSF(2,13,21:30,5) = (/0.152_EB,0.154_EB,0.156_EB,0.158_EB,0.160_EB,0.156_EB,0.152_EB,0.149_EB,0.145_EB,0.141_EB/)
ABSF(2,13,0:10,6) = (/0.000_EB,0.039_EB,0.074_EB,0.081_EB,0.088_EB,0.099_EB,0.110_EB,0.121_EB,0.132_EB,0.112_EB,0.093_EB/)
ABSF(2,13,11:20,6) = (/0.074_EB,0.055_EB,0.056_EB,0.057_EB,0.058_EB,0.059_EB,0.076_EB,0.093_EB,0.110_EB,0.128_EB/)
ABSF(2,13,21:30,6) = (/0.136_EB,0.144_EB,0.152_EB,0.160_EB,0.168_EB,0.169_EB,0.171_EB,0.172_EB,0.174_EB,0.176_EB/)
ABSF(2,14,0:10,1) = (/0.000_EB,0.060_EB,0.071_EB,0.072_EB,0.073_EB,0.069_EB,0.065_EB,0.061_EB,0.057_EB,0.063_EB,0.069_EB/)
ABSF(2,14,11:20,1) = (/0.075_EB,0.081_EB,0.080_EB,0.080_EB,0.080_EB,0.080_EB,0.082_EB,0.085_EB,0.088_EB,0.090_EB/)
ABSF(2,14,21:30,1) = (/0.090_EB,0.090_EB,0.090_EB,0.089_EB,0.089_EB,0.088_EB,0.088_EB,0.087_EB,0.087_EB,0.086_EB/)
ABSF(2,14,0:10,2) = (/0.000_EB,0.060_EB,0.071_EB,0.072_EB,0.073_EB,0.069_EB,0.065_EB,0.061_EB,0.057_EB,0.063_EB,0.069_EB/)
ABSF(2,14,11:20,2) = (/0.075_EB,0.081_EB,0.080_EB,0.080_EB,0.080_EB,0.080_EB,0.082_EB,0.085_EB,0.088_EB,0.091_EB/)
ABSF(2,14,21:30,2) = (/0.090_EB,0.090_EB,0.090_EB,0.089_EB,0.089_EB,0.088_EB,0.087_EB,0.087_EB,0.086_EB,0.085_EB/)
ABSF(2,14,0:10,3) = (/0.000_EB,0.048_EB,0.067_EB,0.071_EB,0.076_EB,0.072_EB,0.067_EB,0.063_EB,0.058_EB,0.064_EB,0.069_EB/)
ABSF(2,14,11:20,3) = (/0.075_EB,0.080_EB,0.080_EB,0.080_EB,0.080_EB,0.080_EB,0.083_EB,0.086_EB,0.089_EB,0.092_EB/)
ABSF(2,14,21:30,3) = (/0.093_EB,0.093_EB,0.093_EB,0.093_EB,0.093_EB,0.093_EB,0.092_EB,0.091_EB,0.091_EB,0.090_EB/)
ABSF(2,14,0:10,4) = (/0.000_EB,0.042_EB,0.062_EB,0.073_EB,0.084_EB,0.085_EB,0.085_EB,0.086_EB,0.087_EB,0.087_EB,0.088_EB/)
ABSF(2,14,11:20,4) = (/0.089_EB,0.090_EB,0.090_EB,0.091_EB,0.091_EB,0.091_EB,0.101_EB,0.111_EB,0.121_EB,0.131_EB/)
ABSF(2,14,21:30,4) = (/0.131_EB,0.131_EB,0.130_EB,0.130_EB,0.130_EB,0.129_EB,0.128_EB,0.127_EB,0.127_EB,0.126_EB/)
ABSF(2,14,0:10,5) = (/0.000_EB,0.040_EB,0.065_EB,0.072_EB,0.080_EB,0.087_EB,0.094_EB,0.101_EB,0.108_EB,0.098_EB,0.087_EB/)
ABSF(2,14,11:20,5) = (/0.077_EB,0.066_EB,0.068_EB,0.070_EB,0.072_EB,0.074_EB,0.093_EB,0.113_EB,0.133_EB,0.152_EB/)
ABSF(2,14,21:30,5) = (/0.154_EB,0.156_EB,0.157_EB,0.159_EB,0.161_EB,0.158_EB,0.154_EB,0.151_EB,0.148_EB,0.144_EB/)
ABSF(2,14,0:10,6) = (/0.000_EB,0.041_EB,0.076_EB,0.083_EB,0.091_EB,0.102_EB,0.113_EB,0.123_EB,0.134_EB,0.115_EB,0.097_EB/)
ABSF(2,14,11:20,6) = (/0.078_EB,0.059_EB,0.060_EB,0.061_EB,0.062_EB,0.063_EB,0.080_EB,0.098_EB,0.116_EB,0.133_EB/)
ABSF(2,14,21:30,6) = (/0.141_EB,0.148_EB,0.156_EB,0.163_EB,0.171_EB,0.172_EB,0.173_EB,0.175_EB,0.176_EB,0.177_EB/)
ABSF(2,15,0:10,1) = (/0.000_EB,0.063_EB,0.074_EB,0.075_EB,0.076_EB,0.072_EB,0.068_EB,0.064_EB,0.060_EB,0.066_EB,0.072_EB/)
ABSF(2,15,11:20,1) = (/0.078_EB,0.084_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.086_EB,0.088_EB,0.091_EB,0.093_EB/)
ABSF(2,15,21:30,1) = (/0.093_EB,0.093_EB,0.093_EB,0.092_EB,0.092_EB,0.092_EB,0.091_EB,0.090_EB,0.090_EB,0.089_EB/)
ABSF(2,15,0:10,2) = (/0.000_EB,0.063_EB,0.074_EB,0.075_EB,0.076_EB,0.072_EB,0.068_EB,0.064_EB,0.060_EB,0.066_EB,0.072_EB/)
ABSF(2,15,11:20,2) = (/0.078_EB,0.084_EB,0.084_EB,0.083_EB,0.083_EB,0.083_EB,0.086_EB,0.088_EB,0.091_EB,0.094_EB/)
ABSF(2,15,21:30,2) = (/0.093_EB,0.093_EB,0.093_EB,0.092_EB,0.092_EB,0.091_EB,0.090_EB,0.090_EB,0.089_EB,0.088_EB/)
ABSF(2,15,0:10,3) = (/0.000_EB,0.051_EB,0.071_EB,0.075_EB,0.079_EB,0.075_EB,0.070_EB,0.066_EB,0.061_EB,0.067_EB,0.072_EB/)
ABSF(2,15,11:20,3) = (/0.077_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.086_EB,0.089_EB,0.092_EB,0.095_EB/)
ABSF(2,15,21:30,3) = (/0.095_EB,0.095_EB,0.096_EB,0.096_EB,0.096_EB,0.095_EB,0.094_EB,0.094_EB,0.093_EB,0.092_EB/)
ABSF(2,15,0:10,4) = (/0.000_EB,0.044_EB,0.065_EB,0.076_EB,0.087_EB,0.087_EB,0.087_EB,0.087_EB,0.086_EB,0.088_EB,0.089_EB/)
ABSF(2,15,11:20,4) = (/0.090_EB,0.092_EB,0.092_EB,0.092_EB,0.093_EB,0.093_EB,0.103_EB,0.113_EB,0.123_EB,0.133_EB/)
ABSF(2,15,21:30,4) = (/0.133_EB,0.132_EB,0.132_EB,0.132_EB,0.131_EB,0.130_EB,0.129_EB,0.128_EB,0.127_EB,0.126_EB/)
ABSF(2,15,0:10,5) = (/0.000_EB,0.042_EB,0.068_EB,0.075_EB,0.083_EB,0.090_EB,0.097_EB,0.105_EB,0.112_EB,0.101_EB,0.091_EB/)
ABSF(2,15,11:20,5) = (/0.080_EB,0.069_EB,0.072_EB,0.074_EB,0.076_EB,0.078_EB,0.097_EB,0.116_EB,0.135_EB,0.154_EB/)
ABSF(2,15,21:30,5) = (/0.156_EB,0.157_EB,0.159_EB,0.160_EB,0.162_EB,0.159_EB,0.156_EB,0.153_EB,0.151_EB,0.148_EB/)
ABSF(2,15,0:10,6) = (/0.000_EB,0.043_EB,0.078_EB,0.086_EB,0.094_EB,0.105_EB,0.115_EB,0.126_EB,0.137_EB,0.119_EB,0.100_EB/)
ABSF(2,15,11:20,6) = (/0.082_EB,0.064_EB,0.064_EB,0.065_EB,0.066_EB,0.066_EB,0.084_EB,0.103_EB,0.121_EB,0.139_EB/)
ABSF(2,15,21:30,6) = (/0.146_EB,0.153_EB,0.160_EB,0.167_EB,0.174_EB,0.175_EB,0.176_EB,0.177_EB,0.178_EB,0.179_EB/)
ABSF(2,16,0:10,1) = (/0.000_EB,0.065_EB,0.078_EB,0.079_EB,0.080_EB,0.075_EB,0.071_EB,0.067_EB,0.063_EB,0.069_EB,0.075_EB/)
ABSF(2,16,11:20,1) = (/0.081_EB,0.086_EB,0.086_EB,0.086_EB,0.086_EB,0.086_EB,0.089_EB,0.091_EB,0.094_EB,0.096_EB/)
ABSF(2,16,21:30,1) = (/0.096_EB,0.096_EB,0.096_EB,0.096_EB,0.096_EB,0.095_EB,0.094_EB,0.093_EB,0.093_EB,0.092_EB/)
ABSF(2,16,0:10,2) = (/0.000_EB,0.065_EB,0.077_EB,0.079_EB,0.080_EB,0.076_EB,0.071_EB,0.067_EB,0.063_EB,0.069_EB,0.075_EB/)
ABSF(2,16,11:20,2) = (/0.081_EB,0.087_EB,0.087_EB,0.087_EB,0.087_EB,0.086_EB,0.089_EB,0.092_EB,0.094_EB,0.097_EB/)
ABSF(2,16,21:30,2) = (/0.097_EB,0.096_EB,0.096_EB,0.095_EB,0.095_EB,0.094_EB,0.093_EB,0.093_EB,0.092_EB,0.091_EB/)
ABSF(2,16,0:10,3) = (/0.000_EB,0.054_EB,0.075_EB,0.079_EB,0.082_EB,0.078_EB,0.073_EB,0.068_EB,0.064_EB,0.069_EB,0.075_EB/)
ABSF(2,16,11:20,3) = (/0.080_EB,0.086_EB,0.086_EB,0.086_EB,0.086_EB,0.086_EB,0.089_EB,0.092_EB,0.095_EB,0.098_EB/)
ABSF(2,16,21:30,3) = (/0.098_EB,0.098_EB,0.098_EB,0.098_EB,0.098_EB,0.098_EB,0.097_EB,0.096_EB,0.095_EB,0.094_EB/)
ABSF(2,16,0:10,4) = (/0.000_EB,0.045_EB,0.068_EB,0.079_EB,0.090_EB,0.089_EB,0.088_EB,0.087_EB,0.086_EB,0.088_EB,0.090_EB/)
ABSF(2,16,11:20,4) = (/0.091_EB,0.093_EB,0.094_EB,0.094_EB,0.095_EB,0.096_EB,0.105_EB,0.115_EB,0.125_EB,0.135_EB/)
ABSF(2,16,21:30,4) = (/0.134_EB,0.134_EB,0.133_EB,0.133_EB,0.132_EB,0.131_EB,0.130_EB,0.129_EB,0.127_EB,0.126_EB/)
ABSF(2,16,0:10,5) = (/0.000_EB,0.044_EB,0.071_EB,0.078_EB,0.085_EB,0.093_EB,0.101_EB,0.108_EB,0.116_EB,0.105_EB,0.094_EB/)
ABSF(2,16,11:20,5) = (/0.083_EB,0.072_EB,0.075_EB,0.077_EB,0.080_EB,0.082_EB,0.101_EB,0.119_EB,0.137_EB,0.156_EB/)
ABSF(2,16,21:30,5) = (/0.157_EB,0.159_EB,0.160_EB,0.162_EB,0.163_EB,0.161_EB,0.158_EB,0.156_EB,0.153_EB,0.151_EB/)
ABSF(2,16,0:10,6) = (/0.000_EB,0.045_EB,0.080_EB,0.089_EB,0.097_EB,0.108_EB,0.118_EB,0.129_EB,0.140_EB,0.122_EB,0.104_EB/)
ABSF(2,16,11:20,6) = (/0.086_EB,0.068_EB,0.069_EB,0.069_EB,0.069_EB,0.070_EB,0.089_EB,0.107_EB,0.126_EB,0.145_EB/)
ABSF(2,16,21:30,6) = (/0.151_EB,0.158_EB,0.164_EB,0.170_EB,0.177_EB,0.177_EB,0.178_EB,0.179_EB,0.180_EB,0.180_EB/)
ABSF(2,17,0:10,1) = (/0.000_EB,0.068_EB,0.081_EB,0.082_EB,0.083_EB,0.079_EB,0.074_EB,0.070_EB,0.066_EB,0.072_EB,0.077_EB/)
ABSF(2,17,11:20,1) = (/0.083_EB,0.089_EB,0.089_EB,0.090_EB,0.090_EB,0.090_EB,0.092_EB,0.095_EB,0.097_EB,0.099_EB/)
ABSF(2,17,21:30,1) = (/0.099_EB,0.099_EB,0.099_EB,0.099_EB,0.099_EB,0.098_EB,0.097_EB,0.097_EB,0.096_EB,0.095_EB/)
ABSF(2,17,0:10,2) = (/0.000_EB,0.067_EB,0.081_EB,0.082_EB,0.083_EB,0.079_EB,0.074_EB,0.070_EB,0.065_EB,0.072_EB,0.078_EB/)
ABSF(2,17,11:20,2) = (/0.084_EB,0.090_EB,0.090_EB,0.090_EB,0.090_EB,0.090_EB,0.092_EB,0.095_EB,0.097_EB,0.100_EB/)
ABSF(2,17,21:30,2) = (/0.100_EB,0.099_EB,0.099_EB,0.098_EB,0.098_EB,0.097_EB,0.096_EB,0.096_EB,0.095_EB,0.094_EB/)
ABSF(2,17,0:10,3) = (/0.000_EB,0.058_EB,0.079_EB,0.082_EB,0.085_EB,0.081_EB,0.076_EB,0.071_EB,0.066_EB,0.072_EB,0.078_EB/)
ABSF(2,17,11:20,3) = (/0.083_EB,0.089_EB,0.089_EB,0.089_EB,0.089_EB,0.090_EB,0.092_EB,0.095_EB,0.098_EB,0.100_EB/)
ABSF(2,17,21:30,3) = (/0.100_EB,0.100_EB,0.101_EB,0.101_EB,0.101_EB,0.100_EB,0.099_EB,0.098_EB,0.097_EB,0.096_EB/)
ABSF(2,17,0:10,4) = (/0.000_EB,0.047_EB,0.070_EB,0.082_EB,0.093_EB,0.091_EB,0.089_EB,0.088_EB,0.086_EB,0.088_EB,0.090_EB/)
ABSF(2,17,11:20,4) = (/0.092_EB,0.094_EB,0.095_EB,0.096_EB,0.097_EB,0.098_EB,0.108_EB,0.117_EB,0.127_EB,0.137_EB/)
ABSF(2,17,21:30,4) = (/0.136_EB,0.136_EB,0.135_EB,0.134_EB,0.134_EB,0.132_EB,0.131_EB,0.129_EB,0.128_EB,0.126_EB/)
ABSF(2,17,0:10,5) = (/0.000_EB,0.046_EB,0.074_EB,0.081_EB,0.088_EB,0.096_EB,0.104_EB,0.112_EB,0.120_EB,0.109_EB,0.098_EB/)
ABSF(2,17,11:20,5) = (/0.087_EB,0.075_EB,0.078_EB,0.081_EB,0.084_EB,0.087_EB,0.105_EB,0.122_EB,0.140_EB,0.158_EB/)
ABSF(2,17,21:30,5) = (/0.159_EB,0.160_EB,0.162_EB,0.163_EB,0.164_EB,0.162_EB,0.160_EB,0.158_EB,0.156_EB,0.154_EB/)
ABSF(2,17,0:10,6) = (/0.000_EB,0.046_EB,0.082_EB,0.091_EB,0.101_EB,0.111_EB,0.121_EB,0.132_EB,0.142_EB,0.125_EB,0.107_EB/)
ABSF(2,17,11:20,6) = (/0.090_EB,0.073_EB,0.073_EB,0.073_EB,0.073_EB,0.073_EB,0.093_EB,0.112_EB,0.131_EB,0.151_EB/)
ABSF(2,17,21:30,6) = (/0.156_EB,0.162_EB,0.168_EB,0.174_EB,0.180_EB,0.180_EB,0.181_EB,0.181_EB,0.181_EB,0.182_EB/)
ABSF(2,18,0:10,1) = (/0.000_EB,0.070_EB,0.085_EB,0.086_EB,0.086_EB,0.082_EB,0.078_EB,0.073_EB,0.069_EB,0.075_EB,0.080_EB/)
ABSF(2,18,11:20,1) = (/0.086_EB,0.092_EB,0.092_EB,0.093_EB,0.093_EB,0.093_EB,0.095_EB,0.098_EB,0.100_EB,0.102_EB/)
ABSF(2,18,21:30,1) = (/0.102_EB,0.102_EB,0.102_EB,0.102_EB,0.102_EB,0.101_EB,0.100_EB,0.100_EB,0.099_EB,0.098_EB/)
ABSF(2,18,0:10,2) = (/0.000_EB,0.070_EB,0.084_EB,0.085_EB,0.087_EB,0.082_EB,0.077_EB,0.073_EB,0.068_EB,0.074_EB,0.081_EB/)
ABSF(2,18,11:20,2) = (/0.087_EB,0.093_EB,0.093_EB,0.093_EB,0.093_EB,0.093_EB,0.096_EB,0.098_EB,0.101_EB,0.103_EB/)
ABSF(2,18,21:30,2) = (/0.103_EB,0.102_EB,0.102_EB,0.102_EB,0.101_EB,0.100_EB,0.099_EB,0.099_EB,0.098_EB,0.097_EB/)
ABSF(2,18,0:10,3) = (/0.000_EB,0.061_EB,0.083_EB,0.086_EB,0.088_EB,0.084_EB,0.079_EB,0.074_EB,0.069_EB,0.075_EB,0.081_EB/)
ABSF(2,18,11:20,3) = (/0.086_EB,0.092_EB,0.092_EB,0.092_EB,0.093_EB,0.093_EB,0.095_EB,0.098_EB,0.100_EB,0.103_EB/)
ABSF(2,18,21:30,3) = (/0.103_EB,0.103_EB,0.103_EB,0.103_EB,0.103_EB,0.102_EB,0.101_EB,0.100_EB,0.099_EB,0.098_EB/)
ABSF(2,18,0:10,4) = (/0.000_EB,0.048_EB,0.073_EB,0.084_EB,0.096_EB,0.093_EB,0.091_EB,0.088_EB,0.086_EB,0.088_EB,0.091_EB/)
ABSF(2,18,11:20,4) = (/0.093_EB,0.095_EB,0.097_EB,0.098_EB,0.099_EB,0.101_EB,0.110_EB,0.120_EB,0.129_EB,0.139_EB/)
ABSF(2,18,21:30,4) = (/0.138_EB,0.137_EB,0.136_EB,0.135_EB,0.135_EB,0.133_EB,0.131_EB,0.130_EB,0.128_EB,0.126_EB/)
ABSF(2,18,0:10,5) = (/0.000_EB,0.048_EB,0.077_EB,0.084_EB,0.091_EB,0.100_EB,0.108_EB,0.116_EB,0.124_EB,0.113_EB,0.101_EB/)
ABSF(2,18,11:20,5) = (/0.090_EB,0.079_EB,0.082_EB,0.085_EB,0.088_EB,0.091_EB,0.108_EB,0.125_EB,0.142_EB,0.159_EB/)
ABSF(2,18,21:30,5) = (/0.161_EB,0.162_EB,0.163_EB,0.164_EB,0.166_EB,0.164_EB,0.162_EB,0.161_EB,0.159_EB,0.158_EB/)
ABSF(2,18,0:10,6) = (/0.000_EB,0.048_EB,0.084_EB,0.094_EB,0.104_EB,0.114_EB,0.124_EB,0.135_EB,0.145_EB,0.128_EB,0.111_EB/)
ABSF(2,18,11:20,6) = (/0.094_EB,0.077_EB,0.077_EB,0.077_EB,0.077_EB,0.077_EB,0.097_EB,0.117_EB,0.136_EB,0.156_EB/)
ABSF(2,18,21:30,6) = (/0.162_EB,0.167_EB,0.172_EB,0.178_EB,0.183_EB,0.183_EB,0.183_EB,0.183_EB,0.183_EB,0.183_EB/)
ABSF(2,19,0:10,1) = (/0.000_EB,0.073_EB,0.088_EB,0.089_EB,0.090_EB,0.085_EB,0.081_EB,0.076_EB,0.072_EB,0.077_EB,0.083_EB/)
ABSF(2,19,11:20,1) = (/0.089_EB,0.095_EB,0.095_EB,0.096_EB,0.096_EB,0.096_EB,0.099_EB,0.101_EB,0.103_EB,0.105_EB/)
ABSF(2,19,21:30,1) = (/0.105_EB,0.105_EB,0.105_EB,0.105_EB,0.105_EB,0.104_EB,0.104_EB,0.103_EB,0.102_EB,0.101_EB/)
ABSF(2,19,0:10,2) = (/0.000_EB,0.072_EB,0.087_EB,0.089_EB,0.090_EB,0.085_EB,0.081_EB,0.076_EB,0.071_EB,0.077_EB,0.084_EB/)
ABSF(2,19,11:20,2) = (/0.090_EB,0.097_EB,0.097_EB,0.097_EB,0.096_EB,0.096_EB,0.099_EB,0.101_EB,0.104_EB,0.106_EB/)
ABSF(2,19,21:30,2) = (/0.106_EB,0.105_EB,0.105_EB,0.105_EB,0.104_EB,0.103_EB,0.102_EB,0.101_EB,0.101_EB,0.100_EB/)
ABSF(2,19,0:10,3) = (/0.000_EB,0.064_EB,0.087_EB,0.089_EB,0.091_EB,0.087_EB,0.082_EB,0.077_EB,0.072_EB,0.078_EB,0.084_EB/)
ABSF(2,19,11:20,3) = (/0.089_EB,0.095_EB,0.095_EB,0.096_EB,0.096_EB,0.096_EB,0.098_EB,0.101_EB,0.103_EB,0.105_EB/)
ABSF(2,19,21:30,3) = (/0.106_EB,0.106_EB,0.106_EB,0.106_EB,0.106_EB,0.105_EB,0.104_EB,0.102_EB,0.101_EB,0.100_EB/)
ABSF(2,19,0:10,4) = (/0.000_EB,0.050_EB,0.075_EB,0.087_EB,0.099_EB,0.096_EB,0.092_EB,0.089_EB,0.086_EB,0.088_EB,0.091_EB/)
ABSF(2,19,11:20,4) = (/0.094_EB,0.097_EB,0.098_EB,0.100_EB,0.101_EB,0.103_EB,0.112_EB,0.122_EB,0.131_EB,0.141_EB/)
ABSF(2,19,21:30,4) = (/0.140_EB,0.139_EB,0.138_EB,0.137_EB,0.136_EB,0.134_EB,0.132_EB,0.130_EB,0.128_EB,0.127_EB/)
ABSF(2,19,0:10,5) = (/0.000_EB,0.050_EB,0.080_EB,0.087_EB,0.094_EB,0.103_EB,0.111_EB,0.120_EB,0.128_EB,0.117_EB,0.105_EB/)
ABSF(2,19,11:20,5) = (/0.093_EB,0.082_EB,0.085_EB,0.089_EB,0.092_EB,0.096_EB,0.112_EB,0.129_EB,0.145_EB,0.161_EB/)
ABSF(2,19,21:30,5) = (/0.162_EB,0.163_EB,0.165_EB,0.166_EB,0.167_EB,0.166_EB,0.164_EB,0.163_EB,0.162_EB,0.161_EB/)
ABSF(2,19,0:10,6) = (/0.000_EB,0.050_EB,0.086_EB,0.096_EB,0.107_EB,0.117_EB,0.127_EB,0.137_EB,0.147_EB,0.131_EB,0.115_EB/)
ABSF(2,19,11:20,6) = (/0.098_EB,0.082_EB,0.081_EB,0.081_EB,0.081_EB,0.080_EB,0.101_EB,0.121_EB,0.142_EB,0.162_EB/)
ABSF(2,19,21:30,6) = (/0.167_EB,0.172_EB,0.176_EB,0.181_EB,0.186_EB,0.186_EB,0.186_EB,0.185_EB,0.185_EB,0.185_EB/)
ABSF(2,20,0:10,1) = (/0.000_EB,0.075_EB,0.092_EB,0.092_EB,0.093_EB,0.089_EB,0.084_EB,0.079_EB,0.074_EB,0.080_EB,0.086_EB/)
ABSF(2,20,11:20,1) = (/0.092_EB,0.098_EB,0.098_EB,0.099_EB,0.099_EB,0.100_EB,0.102_EB,0.104_EB,0.106_EB,0.108_EB/)
ABSF(2,20,21:30,1) = (/0.108_EB,0.108_EB,0.108_EB,0.108_EB,0.108_EB,0.108_EB,0.107_EB,0.106_EB,0.105_EB,0.104_EB/)
ABSF(2,20,0:10,2) = (/0.000_EB,0.075_EB,0.091_EB,0.092_EB,0.094_EB,0.089_EB,0.084_EB,0.079_EB,0.074_EB,0.080_EB,0.087_EB/)
ABSF(2,20,11:20,2) = (/0.093_EB,0.100_EB,0.100_EB,0.100_EB,0.100_EB,0.100_EB,0.102_EB,0.105_EB,0.107_EB,0.109_EB/)
ABSF(2,20,21:30,2) = (/0.109_EB,0.108_EB,0.108_EB,0.108_EB,0.107_EB,0.106_EB,0.105_EB,0.104_EB,0.104_EB,0.103_EB/)
ABSF(2,20,0:10,3) = (/0.000_EB,0.068_EB,0.091_EB,0.093_EB,0.094_EB,0.089_EB,0.085_EB,0.080_EB,0.075_EB,0.080_EB,0.086_EB/)
ABSF(2,20,11:20,3) = (/0.092_EB,0.098_EB,0.099_EB,0.099_EB,0.099_EB,0.099_EB,0.102_EB,0.104_EB,0.106_EB,0.108_EB/)
ABSF(2,20,21:30,3) = (/0.108_EB,0.108_EB,0.108_EB,0.108_EB,0.108_EB,0.107_EB,0.106_EB,0.105_EB,0.103_EB,0.102_EB/)
ABSF(2,20,0:10,4) = (/0.000_EB,0.051_EB,0.078_EB,0.090_EB,0.102_EB,0.098_EB,0.094_EB,0.090_EB,0.085_EB,0.089_EB,0.092_EB/)
ABSF(2,20,11:20,4) = (/0.095_EB,0.098_EB,0.100_EB,0.102_EB,0.104_EB,0.105_EB,0.115_EB,0.124_EB,0.133_EB,0.143_EB/)
ABSF(2,20,21:30,4) = (/0.141_EB,0.140_EB,0.139_EB,0.138_EB,0.137_EB,0.135_EB,0.133_EB,0.131_EB,0.129_EB,0.127_EB/)
ABSF(2,20,0:10,5) = (/0.000_EB,0.052_EB,0.083_EB,0.090_EB,0.097_EB,0.106_EB,0.115_EB,0.124_EB,0.133_EB,0.121_EB,0.109_EB/)
ABSF(2,20,11:20,5) = (/0.097_EB,0.085_EB,0.089_EB,0.092_EB,0.096_EB,0.100_EB,0.116_EB,0.132_EB,0.147_EB,0.163_EB/)
ABSF(2,20,21:30,5) = (/0.164_EB,0.165_EB,0.166_EB,0.167_EB,0.168_EB,0.167_EB,0.166_EB,0.166_EB,0.165_EB,0.164_EB/)
ABSF(2,20,0:10,6) = (/0.000_EB,0.051_EB,0.087_EB,0.099_EB,0.110_EB,0.120_EB,0.130_EB,0.140_EB,0.150_EB,0.134_EB,0.118_EB/)
ABSF(2,20,11:20,6) = (/0.102_EB,0.086_EB,0.086_EB,0.085_EB,0.085_EB,0.084_EB,0.105_EB,0.126_EB,0.147_EB,0.168_EB/)
ABSF(2,20,21:30,6) = (/0.172_EB,0.176_EB,0.180_EB,0.185_EB,0.189_EB,0.188_EB,0.188_EB,0.187_EB,0.187_EB,0.187_EB/)
ABSF(3,0,0:10,1) = (/0.000_EB,0.011_EB,0.020_EB,0.022_EB,0.023_EB,0.025_EB,0.026_EB,0.027_EB,0.029_EB,0.029_EB,0.029_EB/)
ABSF(3,0,11:20,1) = (/0.029_EB,0.029_EB,0.031_EB,0.033_EB,0.034_EB,0.036_EB,0.036_EB,0.036_EB,0.035_EB,0.035_EB/)
ABSF(3,0,21:30,1) = (/0.036_EB,0.037_EB,0.037_EB,0.038_EB,0.039_EB,0.039_EB,0.040_EB,0.041_EB,0.041_EB,0.042_EB/)
ABSF(3,0,0:10,2) = (/0.000_EB,0.012_EB,0.021_EB,0.026_EB,0.032_EB,0.035_EB,0.039_EB,0.043_EB,0.046_EB,0.045_EB,0.043_EB/)
ABSF(3,0,11:20,2) = (/0.042_EB,0.040_EB,0.042_EB,0.043_EB,0.045_EB,0.046_EB,0.044_EB,0.041_EB,0.039_EB,0.036_EB/)
ABSF(3,0,21:30,2) = (/0.039_EB,0.041_EB,0.043_EB,0.046_EB,0.048_EB,0.049_EB,0.049_EB,0.050_EB,0.051_EB,0.051_EB/)
ABSF(3,0,0:10,3) = (/0.000_EB,0.012_EB,0.023_EB,0.029_EB,0.034_EB,0.041_EB,0.048_EB,0.054_EB,0.061_EB,0.054_EB,0.047_EB/)
ABSF(3,0,11:20,3) = (/0.040_EB,0.034_EB,0.037_EB,0.040_EB,0.043_EB,0.047_EB,0.041_EB,0.035_EB,0.030_EB,0.024_EB/)
ABSF(3,0,21:30,3) = (/0.034_EB,0.045_EB,0.055_EB,0.065_EB,0.075_EB,0.073_EB,0.072_EB,0.070_EB,0.069_EB,0.067_EB/)
ABSF(3,0,0:10,4) = (/0.000_EB,0.016_EB,0.025_EB,0.031_EB,0.038_EB,0.045_EB,0.051_EB,0.058_EB,0.064_EB,0.058_EB,0.052_EB/)
ABSF(3,0,11:20,4) = (/0.046_EB,0.040_EB,0.039_EB,0.039_EB,0.039_EB,0.038_EB,0.035_EB,0.031_EB,0.028_EB,0.024_EB/)
ABSF(3,0,21:30,4) = (/0.037_EB,0.050_EB,0.063_EB,0.076_EB,0.089_EB,0.091_EB,0.093_EB,0.095_EB,0.097_EB,0.100_EB/)
ABSF(3,0,0:10,5) = (/0.000_EB,0.014_EB,0.029_EB,0.036_EB,0.043_EB,0.056_EB,0.069_EB,0.081_EB,0.094_EB,0.083_EB,0.071_EB/)
ABSF(3,0,11:20,5) = (/0.059_EB,0.048_EB,0.040_EB,0.033_EB,0.025_EB,0.018_EB,0.020_EB,0.023_EB,0.025_EB,0.028_EB/)
ABSF(3,0,21:30,5) = (/0.033_EB,0.037_EB,0.042_EB,0.046_EB,0.050_EB,0.065_EB,0.079_EB,0.093_EB,0.107_EB,0.121_EB/)
ABSF(3,0,0:10,6) = (/0.000_EB,0.017_EB,0.040_EB,0.043_EB,0.046_EB,0.050_EB,0.055_EB,0.060_EB,0.065_EB,0.058_EB,0.051_EB/)
ABSF(3,0,11:20,6) = (/0.043_EB,0.036_EB,0.033_EB,0.029_EB,0.025_EB,0.021_EB,0.023_EB,0.025_EB,0.026_EB,0.028_EB/)
ABSF(3,0,21:30,6) = (/0.026_EB,0.023_EB,0.021_EB,0.019_EB,0.016_EB,0.040_EB,0.063_EB,0.087_EB,0.110_EB,0.134_EB/)
ABSF(3,1,0:10,1) = (/0.000_EB,0.020_EB,0.023_EB,0.025_EB,0.028_EB,0.029_EB,0.030_EB,0.030_EB,0.031_EB,0.031_EB,0.030_EB/)
ABSF(3,1,11:20,1) = (/0.029_EB,0.028_EB,0.031_EB,0.033_EB,0.035_EB,0.037_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB/)
ABSF(3,1,21:30,1) = (/0.039_EB,0.039_EB,0.039_EB,0.040_EB,0.040_EB,0.041_EB,0.042_EB,0.042_EB,0.043_EB,0.044_EB/)
ABSF(3,1,0:10,2) = (/0.000_EB,0.016_EB,0.023_EB,0.028_EB,0.034_EB,0.038_EB,0.042_EB,0.046_EB,0.050_EB,0.047_EB,0.044_EB/)
ABSF(3,1,11:20,2) = (/0.041_EB,0.038_EB,0.040_EB,0.042_EB,0.044_EB,0.046_EB,0.042_EB,0.038_EB,0.034_EB,0.030_EB/)
ABSF(3,1,21:30,2) = (/0.033_EB,0.036_EB,0.039_EB,0.043_EB,0.046_EB,0.047_EB,0.048_EB,0.050_EB,0.051_EB,0.052_EB/)
ABSF(3,1,0:10,3) = (/0.000_EB,0.014_EB,0.026_EB,0.032_EB,0.039_EB,0.046_EB,0.053_EB,0.060_EB,0.067_EB,0.062_EB,0.056_EB/)
ABSF(3,1,11:20,3) = (/0.050_EB,0.045_EB,0.046_EB,0.046_EB,0.047_EB,0.048_EB,0.043_EB,0.039_EB,0.034_EB,0.029_EB/)
ABSF(3,1,21:30,3) = (/0.033_EB,0.037_EB,0.041_EB,0.045_EB,0.049_EB,0.052_EB,0.056_EB,0.060_EB,0.064_EB,0.067_EB/)
ABSF(3,1,0:10,4) = (/0.000_EB,0.017_EB,0.027_EB,0.034_EB,0.041_EB,0.050_EB,0.058_EB,0.067_EB,0.075_EB,0.067_EB,0.059_EB/)
ABSF(3,1,11:20,4) = (/0.051_EB,0.043_EB,0.042_EB,0.041_EB,0.040_EB,0.039_EB,0.035_EB,0.030_EB,0.026_EB,0.022_EB/)
ABSF(3,1,21:30,4) = (/0.038_EB,0.055_EB,0.071_EB,0.088_EB,0.104_EB,0.104_EB,0.103_EB,0.102_EB,0.101_EB,0.101_EB/)
ABSF(3,1,0:10,5) = (/0.000_EB,0.018_EB,0.032_EB,0.039_EB,0.046_EB,0.056_EB,0.067_EB,0.077_EB,0.088_EB,0.079_EB,0.071_EB/)
ABSF(3,1,11:20,5) = (/0.062_EB,0.054_EB,0.044_EB,0.034_EB,0.024_EB,0.014_EB,0.019_EB,0.025_EB,0.030_EB,0.035_EB/)
ABSF(3,1,21:30,5) = (/0.053_EB,0.070_EB,0.087_EB,0.104_EB,0.121_EB,0.124_EB,0.126_EB,0.128_EB,0.131_EB,0.133_EB/)
ABSF(3,1,0:10,6) = (/0.000_EB,0.020_EB,0.041_EB,0.043_EB,0.045_EB,0.052_EB,0.059_EB,0.067_EB,0.074_EB,0.065_EB,0.057_EB/)
ABSF(3,1,11:20,6) = (/0.049_EB,0.040_EB,0.036_EB,0.032_EB,0.028_EB,0.024_EB,0.023_EB,0.022_EB,0.021_EB,0.020_EB/)
ABSF(3,1,21:30,6) = (/0.013_EB,0.007_EB,0.000_EB,-0.006_EB,-0.012_EB,0.017_EB,0.047_EB,0.076_EB,0.106_EB,0.135_EB/)
ABSF(3,2,0:10,1) = (/0.000_EB,0.025_EB,0.025_EB,0.027_EB,0.030_EB,0.030_EB,0.031_EB,0.032_EB,0.033_EB,0.032_EB,0.031_EB/)
ABSF(3,2,11:20,1) = (/0.030_EB,0.029_EB,0.031_EB,0.034_EB,0.036_EB,0.039_EB,0.039_EB,0.040_EB,0.040_EB,0.040_EB/)
ABSF(3,2,21:30,1) = (/0.041_EB,0.042_EB,0.044_EB,0.045_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB/)
ABSF(3,2,0:10,2) = (/0.000_EB,0.016_EB,0.025_EB,0.031_EB,0.037_EB,0.038_EB,0.039_EB,0.040_EB,0.041_EB,0.038_EB,0.035_EB/)
ABSF(3,2,11:20,2) = (/0.033_EB,0.030_EB,0.033_EB,0.036_EB,0.039_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB/)
ABSF(3,2,21:30,2) = (/0.045_EB,0.048_EB,0.051_EB,0.054_EB,0.056_EB,0.055_EB,0.054_EB,0.053_EB,0.052_EB,0.051_EB/)
ABSF(3,2,0:10,3) = (/0.000_EB,0.016_EB,0.028_EB,0.034_EB,0.040_EB,0.047_EB,0.054_EB,0.062_EB,0.069_EB,0.062_EB,0.055_EB/)
ABSF(3,2,11:20,3) = (/0.048_EB,0.041_EB,0.042_EB,0.043_EB,0.045_EB,0.046_EB,0.045_EB,0.044_EB,0.043_EB,0.042_EB/)
ABSF(3,2,21:30,3) = (/0.049_EB,0.056_EB,0.063_EB,0.070_EB,0.077_EB,0.076_EB,0.075_EB,0.073_EB,0.072_EB,0.071_EB/)
ABSF(3,2,0:10,4) = (/0.000_EB,0.017_EB,0.032_EB,0.034_EB,0.037_EB,0.047_EB,0.057_EB,0.066_EB,0.076_EB,0.071_EB,0.066_EB/)
ABSF(3,2,11:20,4) = (/0.061_EB,0.055_EB,0.051_EB,0.046_EB,0.042_EB,0.037_EB,0.036_EB,0.034_EB,0.033_EB,0.032_EB/)
ABSF(3,2,21:30,4) = (/0.043_EB,0.054_EB,0.065_EB,0.076_EB,0.087_EB,0.091_EB,0.095_EB,0.098_EB,0.102_EB,0.106_EB/)
ABSF(3,2,0:10,5) = (/0.000_EB,0.018_EB,0.034_EB,0.043_EB,0.052_EB,0.058_EB,0.065_EB,0.072_EB,0.079_EB,0.072_EB,0.065_EB/)
ABSF(3,2,11:20,5) = (/0.058_EB,0.051_EB,0.043_EB,0.035_EB,0.027_EB,0.019_EB,0.025_EB,0.030_EB,0.036_EB,0.042_EB/)
ABSF(3,2,21:30,5) = (/0.040_EB,0.038_EB,0.037_EB,0.035_EB,0.034_EB,0.051_EB,0.069_EB,0.087_EB,0.104_EB,0.122_EB/)
ABSF(3,2,0:10,6) = (/0.000_EB,0.018_EB,0.041_EB,0.046_EB,0.050_EB,0.052_EB,0.054_EB,0.056_EB,0.059_EB,0.054_EB,0.049_EB/)
ABSF(3,2,11:20,6) = (/0.044_EB,0.039_EB,0.033_EB,0.028_EB,0.022_EB,0.017_EB,0.021_EB,0.026_EB,0.031_EB,0.035_EB/)
ABSF(3,2,21:30,6) = (/0.030_EB,0.024_EB,0.018_EB,0.012_EB,0.007_EB,0.032_EB,0.057_EB,0.082_EB,0.107_EB,0.132_EB/)
ABSF(3,3,0:10,1) = (/0.000_EB,0.027_EB,0.028_EB,0.031_EB,0.033_EB,0.034_EB,0.035_EB,0.035_EB,0.036_EB,0.034_EB,0.033_EB/)
ABSF(3,3,11:20,1) = (/0.031_EB,0.029_EB,0.032_EB,0.035_EB,0.038_EB,0.041_EB,0.041_EB,0.042_EB,0.042_EB,0.043_EB/)
ABSF(3,3,21:30,1) = (/0.043_EB,0.044_EB,0.045_EB,0.046_EB,0.047_EB,0.047_EB,0.047_EB,0.047_EB,0.047_EB,0.048_EB/)
ABSF(3,3,0:10,2) = (/0.000_EB,0.017_EB,0.027_EB,0.033_EB,0.039_EB,0.039_EB,0.040_EB,0.041_EB,0.041_EB,0.039_EB,0.036_EB/)
ABSF(3,3,11:20,2) = (/0.033_EB,0.031_EB,0.034_EB,0.037_EB,0.040_EB,0.043_EB,0.043_EB,0.043_EB,0.044_EB,0.044_EB/)
ABSF(3,3,21:30,2) = (/0.046_EB,0.048_EB,0.050_EB,0.052_EB,0.054_EB,0.054_EB,0.053_EB,0.053_EB,0.052_EB,0.052_EB/)
ABSF(3,3,0:10,3) = (/0.000_EB,0.018_EB,0.030_EB,0.036_EB,0.043_EB,0.050_EB,0.057_EB,0.065_EB,0.072_EB,0.065_EB,0.057_EB/)
ABSF(3,3,11:20,3) = (/0.050_EB,0.043_EB,0.045_EB,0.048_EB,0.050_EB,0.053_EB,0.051_EB,0.049_EB,0.046_EB,0.044_EB/)
ABSF(3,3,21:30,3) = (/0.051_EB,0.058_EB,0.065_EB,0.072_EB,0.079_EB,0.078_EB,0.078_EB,0.077_EB,0.077_EB,0.076_EB/)
ABSF(3,3,0:10,4) = (/0.000_EB,0.019_EB,0.034_EB,0.038_EB,0.041_EB,0.050_EB,0.059_EB,0.067_EB,0.076_EB,0.071_EB,0.065_EB/)
ABSF(3,3,11:20,4) = (/0.060_EB,0.054_EB,0.052_EB,0.049_EB,0.047_EB,0.045_EB,0.042_EB,0.040_EB,0.038_EB,0.035_EB/)
ABSF(3,3,21:30,4) = (/0.047_EB,0.059_EB,0.071_EB,0.082_EB,0.094_EB,0.096_EB,0.098_EB,0.100_EB,0.102_EB,0.104_EB/)
ABSF(3,3,0:10,5) = (/0.000_EB,0.020_EB,0.035_EB,0.045_EB,0.056_EB,0.063_EB,0.071_EB,0.078_EB,0.086_EB,0.078_EB,0.071_EB/)
ABSF(3,3,11:20,5) = (/0.063_EB,0.056_EB,0.046_EB,0.037_EB,0.028_EB,0.019_EB,0.024_EB,0.029_EB,0.035_EB,0.040_EB/)
ABSF(3,3,21:30,5) = (/0.045_EB,0.050_EB,0.055_EB,0.060_EB,0.065_EB,0.078_EB,0.090_EB,0.103_EB,0.115_EB,0.128_EB/)
ABSF(3,3,0:10,6) = (/0.000_EB,0.020_EB,0.045_EB,0.049_EB,0.053_EB,0.056_EB,0.059_EB,0.062_EB,0.066_EB,0.061_EB,0.057_EB/)
ABSF(3,3,11:20,6) = (/0.053_EB,0.048_EB,0.041_EB,0.033_EB,0.025_EB,0.018_EB,0.022_EB,0.027_EB,0.032_EB,0.037_EB/)
ABSF(3,3,21:30,6) = (/0.032_EB,0.028_EB,0.024_EB,0.020_EB,0.015_EB,0.038_EB,0.061_EB,0.084_EB,0.107_EB,0.130_EB/)
ABSF(3,4,0:10,1) = (/0.000_EB,0.029_EB,0.031_EB,0.034_EB,0.037_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB,0.036_EB,0.034_EB/)
ABSF(3,4,11:20,1) = (/0.032_EB,0.030_EB,0.033_EB,0.036_EB,0.039_EB,0.043_EB,0.043_EB,0.044_EB,0.044_EB,0.045_EB/)
ABSF(3,4,21:30,1) = (/0.045_EB,0.046_EB,0.047_EB,0.047_EB,0.048_EB,0.048_EB,0.048_EB,0.049_EB,0.049_EB,0.049_EB/)
ABSF(3,4,0:10,2) = (/0.000_EB,0.018_EB,0.030_EB,0.035_EB,0.040_EB,0.041_EB,0.041_EB,0.041_EB,0.042_EB,0.039_EB,0.036_EB/)
ABSF(3,4,11:20,2) = (/0.034_EB,0.031_EB,0.034_EB,0.037_EB,0.040_EB,0.043_EB,0.044_EB,0.045_EB,0.045_EB,0.046_EB/)
ABSF(3,4,21:30,2) = (/0.047_EB,0.048_EB,0.050_EB,0.051_EB,0.052_EB,0.052_EB,0.052_EB,0.053_EB,0.053_EB,0.053_EB/)
ABSF(3,4,0:10,3) = (/0.000_EB,0.020_EB,0.032_EB,0.039_EB,0.045_EB,0.053_EB,0.060_EB,0.067_EB,0.075_EB,0.067_EB,0.060_EB/)
ABSF(3,4,11:20,3) = (/0.052_EB,0.045_EB,0.048_EB,0.052_EB,0.056_EB,0.059_EB,0.056_EB,0.053_EB,0.050_EB,0.047_EB/)
ABSF(3,4,21:30,3) = (/0.053_EB,0.060_EB,0.067_EB,0.074_EB,0.081_EB,0.081_EB,0.081_EB,0.081_EB,0.081_EB,0.082_EB/)
ABSF(3,4,0:10,4) = (/0.000_EB,0.021_EB,0.036_EB,0.041_EB,0.045_EB,0.053_EB,0.060_EB,0.068_EB,0.076_EB,0.070_EB,0.064_EB/)
ABSF(3,4,11:20,4) = (/0.058_EB,0.053_EB,0.053_EB,0.053_EB,0.053_EB,0.053_EB,0.049_EB,0.046_EB,0.042_EB,0.039_EB/)
ABSF(3,4,21:30,4) = (/0.051_EB,0.064_EB,0.076_EB,0.089_EB,0.101_EB,0.101_EB,0.102_EB,0.102_EB,0.102_EB,0.103_EB/)
ABSF(3,4,0:10,5) = (/0.000_EB,0.022_EB,0.036_EB,0.048_EB,0.060_EB,0.068_EB,0.076_EB,0.084_EB,0.093_EB,0.085_EB,0.077_EB/)
ABSF(3,4,11:20,5) = (/0.069_EB,0.060_EB,0.050_EB,0.040_EB,0.029_EB,0.019_EB,0.024_EB,0.029_EB,0.033_EB,0.038_EB/)
ABSF(3,4,21:30,5) = (/0.050_EB,0.062_EB,0.073_EB,0.085_EB,0.097_EB,0.104_EB,0.111_EB,0.119_EB,0.126_EB,0.133_EB/)
ABSF(3,4,0:10,6) = (/0.000_EB,0.021_EB,0.048_EB,0.052_EB,0.055_EB,0.060_EB,0.064_EB,0.068_EB,0.073_EB,0.069_EB,0.066_EB/)
ABSF(3,4,11:20,6) = (/0.062_EB,0.058_EB,0.048_EB,0.038_EB,0.028_EB,0.019_EB,0.023_EB,0.028_EB,0.033_EB,0.038_EB/)
ABSF(3,4,21:30,6) = (/0.035_EB,0.033_EB,0.030_EB,0.027_EB,0.024_EB,0.045_EB,0.065_EB,0.086_EB,0.107_EB,0.128_EB/)
ABSF(3,5,0:10,1) = (/0.000_EB,0.030_EB,0.034_EB,0.037_EB,0.041_EB,0.041_EB,0.041_EB,0.041_EB,0.041_EB,0.038_EB,0.035_EB/)
ABSF(3,5,11:20,1) = (/0.033_EB,0.030_EB,0.034_EB,0.037_EB,0.041_EB,0.044_EB,0.045_EB,0.046_EB,0.046_EB,0.047_EB/)
ABSF(3,5,21:30,1) = (/0.047_EB,0.048_EB,0.048_EB,0.049_EB,0.049_EB,0.049_EB,0.050_EB,0.050_EB,0.050_EB,0.051_EB/)
ABSF(3,5,0:10,2) = (/0.000_EB,0.019_EB,0.033_EB,0.037_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB,0.039_EB,0.037_EB/)
ABSF(3,5,11:20,2) = (/0.034_EB,0.032_EB,0.035_EB,0.038_EB,0.041_EB,0.044_EB,0.045_EB,0.046_EB,0.047_EB,0.048_EB/)
ABSF(3,5,21:30,2) = (/0.048_EB,0.049_EB,0.049_EB,0.050_EB,0.050_EB,0.051_EB,0.052_EB,0.052_EB,0.053_EB,0.054_EB/)
ABSF(3,5,0:10,3) = (/0.000_EB,0.022_EB,0.034_EB,0.041_EB,0.048_EB,0.056_EB,0.063_EB,0.070_EB,0.077_EB,0.070_EB,0.062_EB/)
ABSF(3,5,11:20,3) = (/0.054_EB,0.047_EB,0.051_EB,0.056_EB,0.061_EB,0.066_EB,0.062_EB,0.057_EB,0.053_EB,0.049_EB/)
ABSF(3,5,21:30,3) = (/0.056_EB,0.062_EB,0.069_EB,0.076_EB,0.083_EB,0.084_EB,0.084_EB,0.085_EB,0.086_EB,0.087_EB/)
ABSF(3,5,0:10,4) = (/0.000_EB,0.023_EB,0.038_EB,0.044_EB,0.049_EB,0.056_EB,0.062_EB,0.069_EB,0.076_EB,0.070_EB,0.064_EB/)
ABSF(3,5,11:20,4) = (/0.057_EB,0.051_EB,0.054_EB,0.056_EB,0.058_EB,0.061_EB,0.056_EB,0.051_EB,0.047_EB,0.042_EB/)
ABSF(3,5,21:30,4) = (/0.055_EB,0.068_EB,0.082_EB,0.095_EB,0.108_EB,0.106_EB,0.105_EB,0.104_EB,0.103_EB,0.101_EB/)
ABSF(3,5,0:10,5) = (/0.000_EB,0.023_EB,0.037_EB,0.050_EB,0.064_EB,0.073_EB,0.082_EB,0.091_EB,0.099_EB,0.091_EB,0.082_EB/)
ABSF(3,5,11:20,5) = (/0.074_EB,0.065_EB,0.054_EB,0.042_EB,0.031_EB,0.019_EB,0.023_EB,0.028_EB,0.032_EB,0.037_EB/)
ABSF(3,5,21:30,5) = (/0.055_EB,0.073_EB,0.092_EB,0.110_EB,0.129_EB,0.131_EB,0.133_EB,0.135_EB,0.137_EB,0.139_EB/)
ABSF(3,5,0:10,6) = (/0.000_EB,0.023_EB,0.051_EB,0.055_EB,0.058_EB,0.064_EB,0.069_EB,0.075_EB,0.080_EB,0.077_EB,0.074_EB/)
ABSF(3,5,11:20,6) = (/0.071_EB,0.068_EB,0.056_EB,0.044_EB,0.032_EB,0.019_EB,0.025_EB,0.030_EB,0.035_EB,0.040_EB/)
ABSF(3,5,21:30,6) = (/0.038_EB,0.037_EB,0.035_EB,0.034_EB,0.033_EB,0.051_EB,0.070_EB,0.088_EB,0.107_EB,0.125_EB/)
ABSF(3,6,0:10,1) = (/0.000_EB,0.032_EB,0.036_EB,0.040_EB,0.044_EB,0.044_EB,0.044_EB,0.043_EB,0.043_EB,0.040_EB,0.037_EB/)
ABSF(3,6,11:20,1) = (/0.035_EB,0.032_EB,0.036_EB,0.039_EB,0.043_EB,0.047_EB,0.047_EB,0.048_EB,0.048_EB,0.049_EB/)
ABSF(3,6,21:30,1) = (/0.050_EB,0.050_EB,0.051_EB,0.051_EB,0.052_EB,0.052_EB,0.053_EB,0.053_EB,0.054_EB,0.054_EB/)
ABSF(3,6,0:10,2) = (/0.000_EB,0.023_EB,0.036_EB,0.040_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.041_EB,0.038_EB/)
ABSF(3,6,11:20,2) = (/0.036_EB,0.033_EB,0.036_EB,0.039_EB,0.043_EB,0.046_EB,0.047_EB,0.048_EB,0.049_EB,0.050_EB/)
ABSF(3,6,21:30,2) = (/0.050_EB,0.051_EB,0.052_EB,0.052_EB,0.053_EB,0.054_EB,0.054_EB,0.055_EB,0.056_EB,0.057_EB/)
ABSF(3,6,0:10,3) = (/0.000_EB,0.023_EB,0.035_EB,0.043_EB,0.051_EB,0.057_EB,0.063_EB,0.068_EB,0.074_EB,0.067_EB,0.060_EB/)
ABSF(3,6,11:20,3) = (/0.052_EB,0.045_EB,0.050_EB,0.055_EB,0.060_EB,0.065_EB,0.061_EB,0.058_EB,0.055_EB,0.051_EB/)
ABSF(3,6,21:30,3) = (/0.057_EB,0.062_EB,0.067_EB,0.073_EB,0.078_EB,0.079_EB,0.080_EB,0.081_EB,0.082_EB,0.083_EB/)
ABSF(3,6,0:10,4) = (/0.000_EB,0.024_EB,0.041_EB,0.047_EB,0.053_EB,0.060_EB,0.066_EB,0.073_EB,0.079_EB,0.074_EB,0.068_EB/)
ABSF(3,6,11:20,4) = (/0.062_EB,0.057_EB,0.058_EB,0.059_EB,0.061_EB,0.062_EB,0.057_EB,0.053_EB,0.048_EB,0.044_EB/)
ABSF(3,6,21:30,4) = (/0.057_EB,0.070_EB,0.084_EB,0.097_EB,0.110_EB,0.109_EB,0.108_EB,0.107_EB,0.105_EB,0.104_EB/)
ABSF(3,6,0:10,5) = (/0.000_EB,0.024_EB,0.041_EB,0.054_EB,0.067_EB,0.075_EB,0.084_EB,0.093_EB,0.101_EB,0.093_EB,0.084_EB/)
ABSF(3,6,11:20,5) = (/0.076_EB,0.068_EB,0.058_EB,0.048_EB,0.038_EB,0.029_EB,0.032_EB,0.036_EB,0.039_EB,0.043_EB/)
ABSF(3,6,21:30,5) = (/0.060_EB,0.077_EB,0.095_EB,0.112_EB,0.129_EB,0.132_EB,0.134_EB,0.136_EB,0.139_EB,0.141_EB/)
ABSF(3,6,0:10,6) = (/0.000_EB,0.024_EB,0.053_EB,0.057_EB,0.062_EB,0.067_EB,0.072_EB,0.077_EB,0.082_EB,0.079_EB,0.075_EB/)
ABSF(3,6,11:20,6) = (/0.072_EB,0.069_EB,0.056_EB,0.044_EB,0.032_EB,0.020_EB,0.026_EB,0.032_EB,0.038_EB,0.044_EB/)
ABSF(3,6,21:30,6) = (/0.041_EB,0.039_EB,0.036_EB,0.034_EB,0.031_EB,0.050_EB,0.069_EB,0.088_EB,0.107_EB,0.126_EB/)
ABSF(3,7,0:10,1) = (/0.000_EB,0.034_EB,0.039_EB,0.043_EB,0.047_EB,0.047_EB,0.046_EB,0.046_EB,0.045_EB,0.042_EB,0.040_EB/)
ABSF(3,7,11:20,1) = (/0.037_EB,0.034_EB,0.038_EB,0.041_EB,0.045_EB,0.049_EB,0.049_EB,0.050_EB,0.050_EB,0.051_EB/)
ABSF(3,7,21:30,1) = (/0.052_EB,0.052_EB,0.053_EB,0.054_EB,0.055_EB,0.055_EB,0.056_EB,0.056_EB,0.057_EB,0.057_EB/)
ABSF(3,7,0:10,2) = (/0.000_EB,0.027_EB,0.039_EB,0.042_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB,0.043_EB,0.040_EB/)
ABSF(3,7,11:20,2) = (/0.037_EB,0.034_EB,0.037_EB,0.041_EB,0.044_EB,0.048_EB,0.049_EB,0.050_EB,0.051_EB,0.052_EB/)
ABSF(3,7,21:30,2) = (/0.053_EB,0.053_EB,0.054_EB,0.055_EB,0.056_EB,0.056_EB,0.057_EB,0.058_EB,0.059_EB,0.059_EB/)
ABSF(3,7,0:10,3) = (/0.000_EB,0.025_EB,0.037_EB,0.046_EB,0.055_EB,0.059_EB,0.063_EB,0.067_EB,0.071_EB,0.064_EB,0.057_EB/)
ABSF(3,7,11:20,3) = (/0.050_EB,0.044_EB,0.049_EB,0.054_EB,0.058_EB,0.063_EB,0.061_EB,0.058_EB,0.056_EB,0.054_EB/)
ABSF(3,7,21:30,3) = (/0.058_EB,0.062_EB,0.066_EB,0.070_EB,0.074_EB,0.075_EB,0.076_EB,0.077_EB,0.078_EB,0.079_EB/)
ABSF(3,7,0:10,4) = (/0.000_EB,0.025_EB,0.043_EB,0.051_EB,0.058_EB,0.064_EB,0.070_EB,0.077_EB,0.083_EB,0.078_EB,0.073_EB/)
ABSF(3,7,11:20,4) = (/0.067_EB,0.062_EB,0.063_EB,0.063_EB,0.063_EB,0.064_EB,0.059_EB,0.054_EB,0.050_EB,0.045_EB/)
ABSF(3,7,21:30,4) = (/0.059_EB,0.072_EB,0.086_EB,0.099_EB,0.113_EB,0.112_EB,0.111_EB,0.109_EB,0.108_EB,0.107_EB/)
ABSF(3,7,0:10,5) = (/0.000_EB,0.025_EB,0.045_EB,0.057_EB,0.069_EB,0.078_EB,0.086_EB,0.095_EB,0.103_EB,0.095_EB,0.086_EB/)
ABSF(3,7,11:20,5) = (/0.078_EB,0.070_EB,0.062_EB,0.054_EB,0.046_EB,0.038_EB,0.041_EB,0.044_EB,0.046_EB,0.049_EB/)
ABSF(3,7,21:30,5) = (/0.065_EB,0.081_EB,0.097_EB,0.114_EB,0.130_EB,0.133_EB,0.135_EB,0.138_EB,0.141_EB,0.144_EB/)
ABSF(3,7,0:10,6) = (/0.000_EB,0.026_EB,0.055_EB,0.060_EB,0.065_EB,0.070_EB,0.074_EB,0.079_EB,0.084_EB,0.080_EB,0.077_EB/)
ABSF(3,7,11:20,6) = (/0.073_EB,0.069_EB,0.057_EB,0.045_EB,0.032_EB,0.020_EB,0.027_EB,0.034_EB,0.041_EB,0.048_EB/)
ABSF(3,7,21:30,6) = (/0.044_EB,0.040_EB,0.037_EB,0.033_EB,0.029_EB,0.049_EB,0.068_EB,0.088_EB,0.108_EB,0.127_EB/)
ABSF(3,8,0:10,1) = (/0.000_EB,0.036_EB,0.042_EB,0.046_EB,0.050_EB,0.050_EB,0.049_EB,0.048_EB,0.047_EB,0.044_EB,0.042_EB/)
ABSF(3,8,11:20,1) = (/0.039_EB,0.036_EB,0.040_EB,0.043_EB,0.047_EB,0.051_EB,0.051_EB,0.052_EB,0.052_EB,0.053_EB/)
ABSF(3,8,21:30,1) = (/0.054_EB,0.055_EB,0.056_EB,0.057_EB,0.058_EB,0.058_EB,0.059_EB,0.060_EB,0.060_EB,0.061_EB/)
ABSF(3,8,0:10,2) = (/0.000_EB,0.031_EB,0.042_EB,0.045_EB,0.048_EB,0.048_EB,0.048_EB,0.048_EB,0.049_EB,0.045_EB,0.042_EB/)
ABSF(3,8,11:20,2) = (/0.039_EB,0.035_EB,0.039_EB,0.043_EB,0.046_EB,0.050_EB,0.051_EB,0.052_EB,0.053_EB,0.054_EB/)
ABSF(3,8,21:30,2) = (/0.055_EB,0.056_EB,0.057_EB,0.058_EB,0.059_EB,0.059_EB,0.060_EB,0.061_EB,0.061_EB,0.062_EB/)
ABSF(3,8,0:10,3) = (/0.000_EB,0.026_EB,0.038_EB,0.048_EB,0.058_EB,0.060_EB,0.063_EB,0.065_EB,0.067_EB,0.061_EB,0.055_EB/)
ABSF(3,8,11:20,3) = (/0.049_EB,0.042_EB,0.047_EB,0.052_EB,0.057_EB,0.062_EB,0.060_EB,0.059_EB,0.057_EB,0.056_EB/)
ABSF(3,8,21:30,3) = (/0.059_EB,0.061_EB,0.064_EB,0.067_EB,0.069_EB,0.070_EB,0.071_EB,0.072_EB,0.073_EB,0.074_EB/)
ABSF(3,8,0:10,4) = (/0.000_EB,0.026_EB,0.046_EB,0.054_EB,0.062_EB,0.068_EB,0.074_EB,0.080_EB,0.086_EB,0.082_EB,0.077_EB/)
ABSF(3,8,11:20,4) = (/0.072_EB,0.068_EB,0.067_EB,0.066_EB,0.066_EB,0.065_EB,0.060_EB,0.056_EB,0.051_EB,0.047_EB/)
ABSF(3,8,21:30,4) = (/0.060_EB,0.074_EB,0.088_EB,0.102_EB,0.116_EB,0.114_EB,0.113_EB,0.112_EB,0.111_EB,0.110_EB/)
ABSF(3,8,0:10,5) = (/0.000_EB,0.027_EB,0.049_EB,0.060_EB,0.072_EB,0.080_EB,0.088_EB,0.097_EB,0.105_EB,0.097_EB,0.088_EB/)
ABSF(3,8,11:20,5) = (/0.080_EB,0.072_EB,0.066_EB,0.060_EB,0.054_EB,0.048_EB,0.050_EB,0.052_EB,0.053_EB,0.055_EB/)
ABSF(3,8,21:30,5) = (/0.070_EB,0.085_EB,0.100_EB,0.115_EB,0.130_EB,0.134_EB,0.137_EB,0.140_EB,0.143_EB,0.146_EB/)
ABSF(3,8,0:10,6) = (/0.000_EB,0.027_EB,0.057_EB,0.063_EB,0.068_EB,0.073_EB,0.077_EB,0.082_EB,0.086_EB,0.082_EB,0.078_EB/)
ABSF(3,8,11:20,6) = (/0.074_EB,0.070_EB,0.057_EB,0.045_EB,0.033_EB,0.020_EB,0.028_EB,0.036_EB,0.044_EB,0.052_EB/)
ABSF(3,8,21:30,6) = (/0.047_EB,0.042_EB,0.037_EB,0.032_EB,0.028_EB,0.048_EB,0.068_EB,0.088_EB,0.108_EB,0.128_EB/)
ABSF(3,9,0:10,1) = (/0.000_EB,0.038_EB,0.045_EB,0.049_EB,0.054_EB,0.052_EB,0.051_EB,0.050_EB,0.049_EB,0.046_EB,0.044_EB/)
ABSF(3,9,11:20,1) = (/0.041_EB,0.038_EB,0.042_EB,0.045_EB,0.049_EB,0.053_EB,0.053_EB,0.054_EB,0.054_EB,0.054_EB/)
ABSF(3,9,21:30,1) = (/0.056_EB,0.057_EB,0.058_EB,0.059_EB,0.061_EB,0.061_EB,0.062_EB,0.063_EB,0.064_EB,0.064_EB/)
ABSF(3,9,0:10,2) = (/0.000_EB,0.036_EB,0.045_EB,0.047_EB,0.050_EB,0.050_EB,0.050_EB,0.051_EB,0.051_EB,0.047_EB,0.044_EB/)
ABSF(3,9,11:20,2) = (/0.040_EB,0.036_EB,0.040_EB,0.044_EB,0.048_EB,0.052_EB,0.053_EB,0.054_EB,0.055_EB,0.056_EB/)
ABSF(3,9,21:30,2) = (/0.057_EB,0.058_EB,0.059_EB,0.060_EB,0.061_EB,0.062_EB,0.063_EB,0.063_EB,0.064_EB,0.064_EB/)
ABSF(3,9,0:10,3) = (/0.000_EB,0.028_EB,0.040_EB,0.050_EB,0.061_EB,0.062_EB,0.062_EB,0.063_EB,0.064_EB,0.058_EB,0.052_EB/)
ABSF(3,9,11:20,3) = (/0.047_EB,0.041_EB,0.046_EB,0.051_EB,0.056_EB,0.061_EB,0.060_EB,0.059_EB,0.059_EB,0.058_EB/)
ABSF(3,9,21:30,3) = (/0.060_EB,0.061_EB,0.062_EB,0.064_EB,0.065_EB,0.066_EB,0.067_EB,0.068_EB,0.069_EB,0.070_EB/)
ABSF(3,9,0:10,4) = (/0.000_EB,0.027_EB,0.048_EB,0.058_EB,0.067_EB,0.073_EB,0.078_EB,0.084_EB,0.090_EB,0.086_EB,0.082_EB/)
ABSF(3,9,11:20,4) = (/0.077_EB,0.073_EB,0.072_EB,0.070_EB,0.068_EB,0.066_EB,0.062_EB,0.057_EB,0.053_EB,0.048_EB/)
ABSF(3,9,21:30,4) = (/0.062_EB,0.076_EB,0.090_EB,0.104_EB,0.118_EB,0.117_EB,0.116_EB,0.115_EB,0.114_EB,0.113_EB/)
ABSF(3,9,0:10,5) = (/0.000_EB,0.028_EB,0.053_EB,0.063_EB,0.074_EB,0.082_EB,0.091_EB,0.099_EB,0.107_EB,0.099_EB,0.091_EB/)
ABSF(3,9,11:20,5) = (/0.082_EB,0.074_EB,0.070_EB,0.066_EB,0.062_EB,0.058_EB,0.059_EB,0.059_EB,0.060_EB,0.061_EB/)
ABSF(3,9,21:30,5) = (/0.075_EB,0.089_EB,0.103_EB,0.117_EB,0.131_EB,0.135_EB,0.138_EB,0.142_EB,0.145_EB,0.149_EB/)
ABSF(3,9,0:10,6) = (/0.000_EB,0.029_EB,0.060_EB,0.065_EB,0.071_EB,0.076_EB,0.080_EB,0.084_EB,0.088_EB,0.084_EB,0.079_EB/)
ABSF(3,9,11:20,6) = (/0.075_EB,0.070_EB,0.058_EB,0.045_EB,0.033_EB,0.020_EB,0.029_EB,0.038_EB,0.047_EB,0.056_EB/)
ABSF(3,9,21:30,6) = (/0.050_EB,0.044_EB,0.038_EB,0.032_EB,0.026_EB,0.046_EB,0.067_EB,0.088_EB,0.108_EB,0.129_EB/)
ABSF(3,10,0:10,1) = (/0.000_EB,0.040_EB,0.047_EB,0.052_EB,0.057_EB,0.055_EB,0.054_EB,0.053_EB,0.051_EB,0.049_EB,0.046_EB/)
ABSF(3,10,11:20,1) = (/0.043_EB,0.040_EB,0.044_EB,0.047_EB,0.051_EB,0.055_EB,0.055_EB,0.056_EB,0.056_EB,0.056_EB/)
ABSF(3,10,21:30,1) = (/0.058_EB,0.059_EB,0.061_EB,0.062_EB,0.063_EB,0.064_EB,0.065_EB,0.066_EB,0.067_EB,0.068_EB/)
ABSF(3,10,0:10,2) = (/0.000_EB,0.040_EB,0.047_EB,0.050_EB,0.052_EB,0.052_EB,0.052_EB,0.053_EB,0.053_EB,0.049_EB,0.045_EB/)
ABSF(3,10,11:20,2) = (/0.041_EB,0.038_EB,0.042_EB,0.046_EB,0.050_EB,0.054_EB,0.055_EB,0.056_EB,0.057_EB,0.058_EB/)
ABSF(3,10,21:30,2) = (/0.059_EB,0.060_EB,0.062_EB,0.063_EB,0.064_EB,0.065_EB,0.065_EB,0.066_EB,0.067_EB,0.067_EB/)
ABSF(3,10,0:10,3) = (/0.000_EB,0.030_EB,0.041_EB,0.053_EB,0.064_EB,0.063_EB,0.062_EB,0.062_EB,0.061_EB,0.055_EB,0.050_EB/)
ABSF(3,10,11:20,3) = (/0.045_EB,0.039_EB,0.044_EB,0.049_EB,0.054_EB,0.059_EB,0.060_EB,0.060_EB,0.060_EB,0.060_EB/)
ABSF(3,10,21:30,3) = (/0.060_EB,0.060_EB,0.060_EB,0.060_EB,0.060_EB,0.062_EB,0.063_EB,0.064_EB,0.065_EB,0.066_EB/)
ABSF(3,10,0:10,4) = (/0.000_EB,0.028_EB,0.051_EB,0.061_EB,0.071_EB,0.077_EB,0.082_EB,0.088_EB,0.093_EB,0.090_EB,0.086_EB/)
ABSF(3,10,11:20,4) = (/0.082_EB,0.079_EB,0.076_EB,0.073_EB,0.071_EB,0.068_EB,0.063_EB,0.059_EB,0.054_EB,0.050_EB/)
ABSF(3,10,21:30,4) = (/0.064_EB,0.078_EB,0.092_EB,0.107_EB,0.121_EB,0.120_EB,0.119_EB,0.118_EB,0.116_EB,0.115_EB/)
ABSF(3,10,0:10,5) = (/0.000_EB,0.029_EB,0.056_EB,0.067_EB,0.077_EB,0.085_EB,0.093_EB,0.101_EB,0.108_EB,0.101_EB,0.093_EB/)
ABSF(3,10,11:20,5) = (/0.085_EB,0.077_EB,0.074_EB,0.072_EB,0.070_EB,0.068_EB,0.067_EB,0.067_EB,0.067_EB,0.067_EB/)
ABSF(3,10,21:30,5) = (/0.080_EB,0.093_EB,0.106_EB,0.119_EB,0.132_EB,0.136_EB,0.140_EB,0.143_EB,0.147_EB,0.151_EB/)
ABSF(3,10,0:10,6) = (/0.000_EB,0.030_EB,0.062_EB,0.068_EB,0.074_EB,0.078_EB,0.082_EB,0.086_EB,0.091_EB,0.086_EB,0.081_EB/)
ABSF(3,10,11:20,6) = (/0.076_EB,0.071_EB,0.058_EB,0.046_EB,0.033_EB,0.021_EB,0.031_EB,0.040_EB,0.050_EB,0.060_EB/)
ABSF(3,10,21:30,6) = (/0.053_EB,0.046_EB,0.038_EB,0.031_EB,0.024_EB,0.045_EB,0.066_EB,0.087_EB,0.109_EB,0.130_EB/)
ABSF(3,11,0:10,1) = (/0.000_EB,0.041_EB,0.050_EB,0.055_EB,0.060_EB,0.059_EB,0.057_EB,0.056_EB,0.054_EB,0.051_EB,0.048_EB/)
ABSF(3,11,11:20,1) = (/0.045_EB,0.042_EB,0.046_EB,0.050_EB,0.054_EB,0.058_EB,0.058_EB,0.058_EB,0.059_EB,0.059_EB/)
ABSF(3,11,21:30,1) = (/0.060_EB,0.062_EB,0.063_EB,0.064_EB,0.066_EB,0.067_EB,0.067_EB,0.068_EB,0.069_EB,0.070_EB/)
ABSF(3,11,0:10,2) = (/0.000_EB,0.042_EB,0.050_EB,0.053_EB,0.055_EB,0.055_EB,0.055_EB,0.055_EB,0.055_EB,0.051_EB,0.047_EB/)
ABSF(3,11,11:20,2) = (/0.043_EB,0.039_EB,0.044_EB,0.048_EB,0.052_EB,0.056_EB,0.057_EB,0.058_EB,0.059_EB,0.060_EB/)
ABSF(3,11,21:30,2) = (/0.062_EB,0.063_EB,0.064_EB,0.066_EB,0.067_EB,0.067_EB,0.068_EB,0.068_EB,0.069_EB,0.070_EB/)
ABSF(3,11,0:10,3) = (/0.000_EB,0.031_EB,0.044_EB,0.055_EB,0.067_EB,0.066_EB,0.065_EB,0.063_EB,0.062_EB,0.057_EB,0.052_EB/)
ABSF(3,11,11:20,3) = (/0.047_EB,0.041_EB,0.046_EB,0.051_EB,0.056_EB,0.061_EB,0.062_EB,0.062_EB,0.062_EB,0.063_EB/)
ABSF(3,11,21:30,3) = (/0.063_EB,0.063_EB,0.063_EB,0.063_EB,0.063_EB,0.064_EB,0.065_EB,0.067_EB,0.068_EB,0.069_EB/)
ABSF(3,11,0:10,4) = (/0.000_EB,0.030_EB,0.052_EB,0.062_EB,0.073_EB,0.078_EB,0.083_EB,0.089_EB,0.094_EB,0.090_EB,0.087_EB/)
ABSF(3,11,11:20,4) = (/0.083_EB,0.079_EB,0.077_EB,0.075_EB,0.072_EB,0.070_EB,0.066_EB,0.062_EB,0.058_EB,0.053_EB/)
ABSF(3,11,21:30,4) = (/0.067_EB,0.081_EB,0.095_EB,0.109_EB,0.123_EB,0.121_EB,0.120_EB,0.119_EB,0.118_EB,0.116_EB/)
ABSF(3,11,0:10,5) = (/0.000_EB,0.030_EB,0.058_EB,0.069_EB,0.080_EB,0.087_EB,0.095_EB,0.103_EB,0.110_EB,0.103_EB,0.095_EB/)
ABSF(3,11,11:20,5) = (/0.087_EB,0.080_EB,0.078_EB,0.075_EB,0.073_EB,0.071_EB,0.070_EB,0.070_EB,0.070_EB,0.069_EB/)
ABSF(3,11,21:30,5) = (/0.082_EB,0.095_EB,0.108_EB,0.121_EB,0.133_EB,0.137_EB,0.141_EB,0.144_EB,0.148_EB,0.152_EB/)
ABSF(3,11,0:10,6) = (/0.000_EB,0.031_EB,0.064_EB,0.071_EB,0.078_EB,0.082_EB,0.086_EB,0.090_EB,0.094_EB,0.089_EB,0.084_EB/)
ABSF(3,11,11:20,6) = (/0.079_EB,0.074_EB,0.061_EB,0.048_EB,0.035_EB,0.022_EB,0.032_EB,0.042_EB,0.052_EB,0.061_EB/)
ABSF(3,11,21:30,6) = (/0.057_EB,0.052_EB,0.047_EB,0.043_EB,0.038_EB,0.057_EB,0.076_EB,0.095_EB,0.114_EB,0.133_EB/)
ABSF(3,12,0:10,1) = (/0.000_EB,0.043_EB,0.052_EB,0.058_EB,0.064_EB,0.062_EB,0.060_EB,0.058_EB,0.057_EB,0.053_EB,0.050_EB/)
ABSF(3,12,11:20,1) = (/0.047_EB,0.043_EB,0.048_EB,0.052_EB,0.056_EB,0.060_EB,0.060_EB,0.061_EB,0.061_EB,0.061_EB/)
ABSF(3,12,21:30,1) = (/0.063_EB,0.064_EB,0.065_EB,0.067_EB,0.068_EB,0.069_EB,0.070_EB,0.071_EB,0.071_EB,0.072_EB/)
ABSF(3,12,0:10,2) = (/0.000_EB,0.044_EB,0.053_EB,0.056_EB,0.059_EB,0.059_EB,0.058_EB,0.058_EB,0.058_EB,0.054_EB,0.049_EB/)
ABSF(3,12,11:20,2) = (/0.045_EB,0.041_EB,0.046_EB,0.050_EB,0.054_EB,0.058_EB,0.059_EB,0.061_EB,0.062_EB,0.063_EB/)
ABSF(3,12,21:30,2) = (/0.064_EB,0.066_EB,0.067_EB,0.068_EB,0.070_EB,0.070_EB,0.071_EB,0.071_EB,0.072_EB,0.072_EB/)
ABSF(3,12,0:10,3) = (/0.000_EB,0.032_EB,0.047_EB,0.058_EB,0.069_EB,0.068_EB,0.067_EB,0.065_EB,0.064_EB,0.059_EB,0.054_EB/)
ABSF(3,12,11:20,3) = (/0.049_EB,0.044_EB,0.048_EB,0.053_EB,0.058_EB,0.063_EB,0.063_EB,0.064_EB,0.065_EB,0.065_EB/)
ABSF(3,12,21:30,3) = (/0.066_EB,0.066_EB,0.066_EB,0.066_EB,0.066_EB,0.067_EB,0.068_EB,0.069_EB,0.071_EB,0.072_EB/)
ABSF(3,12,0:10,4) = (/0.000_EB,0.031_EB,0.053_EB,0.064_EB,0.074_EB,0.079_EB,0.084_EB,0.090_EB,0.095_EB,0.091_EB,0.087_EB/)
ABSF(3,12,11:20,4) = (/0.084_EB,0.080_EB,0.078_EB,0.076_EB,0.074_EB,0.072_EB,0.068_EB,0.065_EB,0.061_EB,0.057_EB/)
ABSF(3,12,21:30,4) = (/0.071_EB,0.084_EB,0.098_EB,0.111_EB,0.125_EB,0.123_EB,0.122_EB,0.120_EB,0.119_EB,0.118_EB/)
ABSF(3,12,0:10,5) = (/0.000_EB,0.031_EB,0.060_EB,0.072_EB,0.083_EB,0.090_EB,0.097_EB,0.104_EB,0.112_EB,0.105_EB,0.097_EB/)
ABSF(3,12,11:20,5) = (/0.090_EB,0.083_EB,0.081_EB,0.079_EB,0.076_EB,0.074_EB,0.073_EB,0.073_EB,0.072_EB,0.071_EB/)
ABSF(3,12,21:30,5) = (/0.084_EB,0.097_EB,0.110_EB,0.122_EB,0.135_EB,0.139_EB,0.142_EB,0.145_EB,0.149_EB,0.152_EB/)
ABSF(3,12,0:10,6) = (/0.000_EB,0.032_EB,0.066_EB,0.074_EB,0.082_EB,0.086_EB,0.090_EB,0.094_EB,0.098_EB,0.093_EB,0.088_EB/)
ABSF(3,12,11:20,6) = (/0.082_EB,0.077_EB,0.064_EB,0.051_EB,0.037_EB,0.024_EB,0.034_EB,0.043_EB,0.053_EB,0.063_EB/)
ABSF(3,12,21:30,6) = (/0.061_EB,0.058_EB,0.056_EB,0.054_EB,0.052_EB,0.069_EB,0.086_EB,0.103_EB,0.120_EB,0.137_EB/)
ABSF(3,13,0:10,1) = (/0.000_EB,0.045_EB,0.055_EB,0.061_EB,0.067_EB,0.065_EB,0.063_EB,0.061_EB,0.059_EB,0.056_EB,0.052_EB/)
ABSF(3,13,11:20,1) = (/0.049_EB,0.045_EB,0.050_EB,0.054_EB,0.058_EB,0.062_EB,0.063_EB,0.063_EB,0.064_EB,0.064_EB/)
ABSF(3,13,21:30,1) = (/0.065_EB,0.066_EB,0.068_EB,0.069_EB,0.070_EB,0.071_EB,0.072_EB,0.073_EB,0.074_EB,0.075_EB/)
ABSF(3,13,0:10,2) = (/0.000_EB,0.046_EB,0.055_EB,0.059_EB,0.063_EB,0.062_EB,0.061_EB,0.060_EB,0.060_EB,0.056_EB,0.052_EB/)
ABSF(3,13,11:20,2) = (/0.047_EB,0.043_EB,0.048_EB,0.052_EB,0.056_EB,0.060_EB,0.062_EB,0.063_EB,0.065_EB,0.066_EB/)
ABSF(3,13,21:30,2) = (/0.067_EB,0.068_EB,0.070_EB,0.071_EB,0.072_EB,0.073_EB,0.073_EB,0.074_EB,0.074_EB,0.074_EB/)
ABSF(3,13,0:10,3) = (/0.000_EB,0.034_EB,0.050_EB,0.061_EB,0.072_EB,0.070_EB,0.069_EB,0.067_EB,0.066_EB,0.061_EB,0.056_EB/)
ABSF(3,13,11:20,3) = (/0.051_EB,0.046_EB,0.050_EB,0.055_EB,0.060_EB,0.064_EB,0.065_EB,0.066_EB,0.067_EB,0.068_EB/)
ABSF(3,13,21:30,3) = (/0.068_EB,0.068_EB,0.069_EB,0.069_EB,0.069_EB,0.070_EB,0.071_EB,0.072_EB,0.073_EB,0.074_EB/)
ABSF(3,13,0:10,4) = (/0.000_EB,0.032_EB,0.055_EB,0.065_EB,0.075_EB,0.080_EB,0.085_EB,0.090_EB,0.096_EB,0.092_EB,0.088_EB/)
ABSF(3,13,11:20,4) = (/0.084_EB,0.081_EB,0.079_EB,0.077_EB,0.076_EB,0.074_EB,0.071_EB,0.068_EB,0.064_EB,0.061_EB/)
ABSF(3,13,21:30,4) = (/0.074_EB,0.087_EB,0.100_EB,0.113_EB,0.127_EB,0.125_EB,0.123_EB,0.122_EB,0.120_EB,0.119_EB/)
ABSF(3,13,0:10,5) = (/0.000_EB,0.032_EB,0.062_EB,0.074_EB,0.085_EB,0.092_EB,0.099_EB,0.106_EB,0.113_EB,0.107_EB,0.100_EB/)
ABSF(3,13,11:20,5) = (/0.093_EB,0.086_EB,0.084_EB,0.082_EB,0.080_EB,0.078_EB,0.076_EB,0.075_EB,0.074_EB,0.073_EB/)
ABSF(3,13,21:30,5) = (/0.086_EB,0.099_EB,0.111_EB,0.124_EB,0.137_EB,0.140_EB,0.143_EB,0.146_EB,0.150_EB,0.153_EB/)
ABSF(3,13,0:10,6) = (/0.000_EB,0.034_EB,0.068_EB,0.077_EB,0.086_EB,0.090_EB,0.094_EB,0.098_EB,0.102_EB,0.097_EB,0.091_EB/)
ABSF(3,13,11:20,6) = (/0.086_EB,0.080_EB,0.067_EB,0.053_EB,0.039_EB,0.026_EB,0.035_EB,0.045_EB,0.055_EB,0.064_EB/)
ABSF(3,13,21:30,6) = (/0.065_EB,0.065_EB,0.065_EB,0.066_EB,0.066_EB,0.081_EB,0.096_EB,0.111_EB,0.126_EB,0.140_EB/)
ABSF(3,14,0:10,1) = (/0.000_EB,0.047_EB,0.058_EB,0.064_EB,0.071_EB,0.068_EB,0.066_EB,0.064_EB,0.062_EB,0.058_EB,0.055_EB/)
ABSF(3,14,11:20,1) = (/0.051_EB,0.047_EB,0.052_EB,0.056_EB,0.060_EB,0.065_EB,0.065_EB,0.066_EB,0.066_EB,0.066_EB/)
ABSF(3,14,21:30,1) = (/0.068_EB,0.069_EB,0.070_EB,0.071_EB,0.072_EB,0.073_EB,0.074_EB,0.075_EB,0.076_EB,0.077_EB/)
ABSF(3,14,0:10,2) = (/0.000_EB,0.048_EB,0.058_EB,0.062_EB,0.066_EB,0.065_EB,0.064_EB,0.063_EB,0.062_EB,0.058_EB,0.054_EB/)
ABSF(3,14,11:20,2) = (/0.049_EB,0.045_EB,0.050_EB,0.054_EB,0.058_EB,0.063_EB,0.064_EB,0.066_EB,0.067_EB,0.069_EB/)
ABSF(3,14,21:30,2) = (/0.070_EB,0.071_EB,0.072_EB,0.074_EB,0.075_EB,0.075_EB,0.076_EB,0.076_EB,0.076_EB,0.077_EB/)
ABSF(3,14,0:10,3) = (/0.000_EB,0.035_EB,0.053_EB,0.064_EB,0.075_EB,0.073_EB,0.071_EB,0.069_EB,0.067_EB,0.062_EB,0.057_EB/)
ABSF(3,14,11:20,3) = (/0.053_EB,0.048_EB,0.052_EB,0.057_EB,0.062_EB,0.066_EB,0.067_EB,0.068_EB,0.069_EB,0.070_EB/)
ABSF(3,14,21:30,3) = (/0.071_EB,0.071_EB,0.071_EB,0.072_EB,0.072_EB,0.073_EB,0.074_EB,0.075_EB,0.076_EB,0.077_EB/)
ABSF(3,14,0:10,4) = (/0.000_EB,0.034_EB,0.056_EB,0.066_EB,0.077_EB,0.081_EB,0.086_EB,0.091_EB,0.096_EB,0.093_EB,0.089_EB/)
ABSF(3,14,11:20,4) = (/0.085_EB,0.081_EB,0.080_EB,0.079_EB,0.078_EB,0.077_EB,0.074_EB,0.071_EB,0.068_EB,0.065_EB/)
ABSF(3,14,21:30,4) = (/0.078_EB,0.090_EB,0.103_EB,0.116_EB,0.128_EB,0.127_EB,0.125_EB,0.123_EB,0.121_EB,0.120_EB/)
ABSF(3,14,0:10,5) = (/0.000_EB,0.033_EB,0.064_EB,0.076_EB,0.088_EB,0.095_EB,0.102_EB,0.108_EB,0.115_EB,0.109_EB,0.102_EB/)
ABSF(3,14,11:20,5) = (/0.096_EB,0.090_EB,0.088_EB,0.085_EB,0.083_EB,0.081_EB,0.079_EB,0.078_EB,0.076_EB,0.075_EB/)
ABSF(3,14,21:30,5) = (/0.088_EB,0.101_EB,0.113_EB,0.126_EB,0.139_EB,0.142_EB,0.145_EB,0.147_EB,0.150_EB,0.153_EB/)
ABSF(3,14,0:10,6) = (/0.000_EB,0.035_EB,0.070_EB,0.080_EB,0.090_EB,0.094_EB,0.098_EB,0.102_EB,0.106_EB,0.100_EB,0.095_EB/)
ABSF(3,14,11:20,6) = (/0.089_EB,0.084_EB,0.070_EB,0.056_EB,0.041_EB,0.027_EB,0.037_EB,0.047_EB,0.056_EB,0.066_EB/)
ABSF(3,14,21:30,6) = (/0.069_EB,0.071_EB,0.074_EB,0.077_EB,0.080_EB,0.093_EB,0.105_EB,0.118_EB,0.131_EB,0.144_EB/)
ABSF(3,15,0:10,1) = (/0.000_EB,0.049_EB,0.060_EB,0.067_EB,0.074_EB,0.072_EB,0.069_EB,0.067_EB,0.064_EB,0.061_EB,0.057_EB/)
ABSF(3,15,11:20,1) = (/0.053_EB,0.049_EB,0.054_EB,0.058_EB,0.063_EB,0.067_EB,0.068_EB,0.068_EB,0.069_EB,0.069_EB/)
ABSF(3,15,21:30,1) = (/0.070_EB,0.071_EB,0.072_EB,0.074_EB,0.075_EB,0.076_EB,0.077_EB,0.077_EB,0.078_EB,0.079_EB/)
ABSF(3,15,0:10,2) = (/0.000_EB,0.050_EB,0.061_EB,0.065_EB,0.070_EB,0.068_EB,0.067_EB,0.066_EB,0.064_EB,0.060_EB,0.056_EB/)
ABSF(3,15,11:20,2) = (/0.052_EB,0.047_EB,0.052_EB,0.056_EB,0.060_EB,0.065_EB,0.066_EB,0.068_EB,0.070_EB,0.071_EB/)
ABSF(3,15,21:30,2) = (/0.073_EB,0.074_EB,0.075_EB,0.076_EB,0.077_EB,0.078_EB,0.078_EB,0.079_EB,0.079_EB,0.079_EB/)
ABSF(3,15,0:10,3) = (/0.000_EB,0.037_EB,0.056_EB,0.067_EB,0.077_EB,0.075_EB,0.073_EB,0.071_EB,0.069_EB,0.064_EB,0.059_EB/)
ABSF(3,15,11:20,3) = (/0.055_EB,0.050_EB,0.054_EB,0.059_EB,0.063_EB,0.068_EB,0.069_EB,0.070_EB,0.072_EB,0.073_EB/)
ABSF(3,15,21:30,3) = (/0.073_EB,0.074_EB,0.074_EB,0.074_EB,0.075_EB,0.076_EB,0.077_EB,0.078_EB,0.079_EB,0.080_EB/)
ABSF(3,15,0:10,4) = (/0.000_EB,0.035_EB,0.057_EB,0.068_EB,0.078_EB,0.083_EB,0.087_EB,0.092_EB,0.097_EB,0.093_EB,0.089_EB/)
ABSF(3,15,11:20,4) = (/0.086_EB,0.082_EB,0.081_EB,0.080_EB,0.080_EB,0.079_EB,0.076_EB,0.074_EB,0.071_EB,0.069_EB/)
ABSF(3,15,21:30,4) = (/0.081_EB,0.093_EB,0.106_EB,0.118_EB,0.130_EB,0.128_EB,0.127_EB,0.125_EB,0.123_EB,0.121_EB/)
ABSF(3,15,0:10,5) = (/0.000_EB,0.035_EB,0.066_EB,0.079_EB,0.091_EB,0.097_EB,0.104_EB,0.110_EB,0.117_EB,0.111_EB,0.105_EB/)
ABSF(3,15,11:20,5) = (/0.099_EB,0.093_EB,0.091_EB,0.089_EB,0.086_EB,0.084_EB,0.082_EB,0.081_EB,0.079_EB,0.077_EB/)
ABSF(3,15,21:30,5) = (/0.090_EB,0.102_EB,0.115_EB,0.128_EB,0.141_EB,0.143_EB,0.146_EB,0.148_EB,0.151_EB,0.154_EB/)
ABSF(3,15,0:10,6) = (/0.000_EB,0.036_EB,0.072_EB,0.083_EB,0.093_EB,0.098_EB,0.102_EB,0.106_EB,0.110_EB,0.104_EB,0.098_EB/)
ABSF(3,15,11:20,6) = (/0.093_EB,0.087_EB,0.072_EB,0.058_EB,0.044_EB,0.029_EB,0.039_EB,0.048_EB,0.058_EB,0.067_EB/)
ABSF(3,15,21:30,6) = (/0.073_EB,0.078_EB,0.083_EB,0.088_EB,0.094_EB,0.104_EB,0.115_EB,0.126_EB,0.137_EB,0.148_EB/)
ABSF(3,16,0:10,1) = (/0.000_EB,0.051_EB,0.063_EB,0.070_EB,0.078_EB,0.075_EB,0.072_EB,0.070_EB,0.067_EB,0.063_EB,0.059_EB/)
ABSF(3,16,11:20,1) = (/0.055_EB,0.051_EB,0.056_EB,0.060_EB,0.065_EB,0.070_EB,0.070_EB,0.071_EB,0.071_EB,0.071_EB/)
ABSF(3,16,21:30,1) = (/0.073_EB,0.074_EB,0.075_EB,0.076_EB,0.077_EB,0.078_EB,0.079_EB,0.080_EB,0.081_EB,0.081_EB/)
ABSF(3,16,0:10,2) = (/0.000_EB,0.051_EB,0.063_EB,0.068_EB,0.073_EB,0.072_EB,0.070_EB,0.068_EB,0.067_EB,0.062_EB,0.058_EB/)
ABSF(3,16,11:20,2) = (/0.054_EB,0.049_EB,0.054_EB,0.058_EB,0.062_EB,0.067_EB,0.069_EB,0.071_EB,0.072_EB,0.074_EB/)
ABSF(3,16,21:30,2) = (/0.075_EB,0.077_EB,0.078_EB,0.079_EB,0.080_EB,0.080_EB,0.081_EB,0.081_EB,0.081_EB,0.082_EB/)
ABSF(3,16,0:10,3) = (/0.000_EB,0.038_EB,0.059_EB,0.069_EB,0.080_EB,0.078_EB,0.075_EB,0.073_EB,0.070_EB,0.066_EB,0.061_EB/)
ABSF(3,16,11:20,3) = (/0.057_EB,0.052_EB,0.056_EB,0.061_EB,0.065_EB,0.070_EB,0.071_EB,0.072_EB,0.074_EB,0.075_EB/)
ABSF(3,16,21:30,3) = (/0.076_EB,0.076_EB,0.077_EB,0.077_EB,0.077_EB,0.079_EB,0.080_EB,0.081_EB,0.082_EB,0.083_EB/)
ABSF(3,16,0:10,4) = (/0.000_EB,0.036_EB,0.059_EB,0.069_EB,0.079_EB,0.084_EB,0.088_EB,0.093_EB,0.098_EB,0.094_EB,0.090_EB/)
ABSF(3,16,11:20,4) = (/0.086_EB,0.082_EB,0.082_EB,0.082_EB,0.081_EB,0.081_EB,0.079_EB,0.077_EB,0.075_EB,0.073_EB/)
ABSF(3,16,21:30,4) = (/0.085_EB,0.096_EB,0.108_EB,0.120_EB,0.132_EB,0.130_EB,0.128_EB,0.126_EB,0.124_EB,0.122_EB/)
ABSF(3,16,0:10,5) = (/0.000_EB,0.036_EB,0.068_EB,0.081_EB,0.094_EB,0.100_EB,0.106_EB,0.112_EB,0.118_EB,0.113_EB,0.107_EB/)
ABSF(3,16,11:20,5) = (/0.102_EB,0.096_EB,0.094_EB,0.092_EB,0.090_EB,0.088_EB,0.085_EB,0.083_EB,0.081_EB,0.079_EB/)
ABSF(3,16,21:30,5) = (/0.092_EB,0.104_EB,0.117_EB,0.130_EB,0.143_EB,0.145_EB,0.147_EB,0.149_EB,0.152_EB,0.154_EB/)
ABSF(3,16,0:10,6) = (/0.000_EB,0.037_EB,0.074_EB,0.085_EB,0.097_EB,0.101_EB,0.106_EB,0.110_EB,0.114_EB,0.108_EB,0.102_EB/)
ABSF(3,16,11:20,6) = (/0.096_EB,0.090_EB,0.075_EB,0.060_EB,0.046_EB,0.031_EB,0.040_EB,0.050_EB,0.059_EB,0.069_EB/)
ABSF(3,16,21:30,6) = (/0.076_EB,0.084_EB,0.092_EB,0.100_EB,0.107_EB,0.116_EB,0.125_EB,0.134_EB,0.143_EB,0.151_EB/)
ABSF(3,17,0:10,1) = (/0.000_EB,0.052_EB,0.065_EB,0.073_EB,0.081_EB,0.078_EB,0.075_EB,0.073_EB,0.070_EB,0.065_EB,0.061_EB/)
ABSF(3,17,11:20,1) = (/0.057_EB,0.053_EB,0.058_EB,0.062_EB,0.067_EB,0.072_EB,0.073_EB,0.073_EB,0.074_EB,0.074_EB/)
ABSF(3,17,21:30,1) = (/0.075_EB,0.076_EB,0.077_EB,0.078_EB,0.079_EB,0.080_EB,0.081_EB,0.082_EB,0.083_EB,0.084_EB/)
ABSF(3,17,0:10,2) = (/0.000_EB,0.053_EB,0.066_EB,0.072_EB,0.077_EB,0.075_EB,0.073_EB,0.071_EB,0.069_EB,0.064_EB,0.060_EB/)
ABSF(3,17,11:20,2) = (/0.056_EB,0.051_EB,0.056_EB,0.060_EB,0.065_EB,0.069_EB,0.071_EB,0.073_EB,0.075_EB,0.077_EB/)
ABSF(3,17,21:30,2) = (/0.078_EB,0.079_EB,0.080_EB,0.082_EB,0.083_EB,0.083_EB,0.083_EB,0.084_EB,0.084_EB,0.084_EB/)
ABSF(3,17,0:10,3) = (/0.000_EB,0.040_EB,0.062_EB,0.072_EB,0.083_EB,0.080_EB,0.077_EB,0.075_EB,0.072_EB,0.067_EB,0.063_EB/)
ABSF(3,17,11:20,3) = (/0.058_EB,0.054_EB,0.058_EB,0.063_EB,0.067_EB,0.071_EB,0.073_EB,0.074_EB,0.076_EB,0.078_EB/)
ABSF(3,17,21:30,3) = (/0.078_EB,0.079_EB,0.079_EB,0.080_EB,0.080_EB,0.081_EB,0.082_EB,0.084_EB,0.085_EB,0.086_EB/)
ABSF(3,17,0:10,4) = (/0.000_EB,0.038_EB,0.060_EB,0.070_EB,0.080_EB,0.085_EB,0.089_EB,0.094_EB,0.099_EB,0.095_EB,0.091_EB/)
ABSF(3,17,11:20,4) = (/0.087_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.083_EB,0.082_EB,0.080_EB,0.078_EB,0.076_EB/)
ABSF(3,17,21:30,4) = (/0.088_EB,0.100_EB,0.111_EB,0.123_EB,0.134_EB,0.132_EB,0.130_EB,0.127_EB,0.125_EB,0.123_EB/)
ABSF(3,17,0:10,5) = (/0.000_EB,0.037_EB,0.070_EB,0.083_EB,0.097_EB,0.102_EB,0.108_EB,0.114_EB,0.120_EB,0.115_EB,0.110_EB/)
ABSF(3,17,11:20,5) = (/0.105_EB,0.099_EB,0.097_EB,0.095_EB,0.093_EB,0.091_EB,0.088_EB,0.086_EB,0.083_EB,0.081_EB/)
ABSF(3,17,21:30,5) = (/0.094_EB,0.106_EB,0.119_EB,0.132_EB,0.144_EB,0.146_EB,0.148_EB,0.150_EB,0.152_EB,0.154_EB/)
ABSF(3,17,0:10,6) = (/0.000_EB,0.039_EB,0.076_EB,0.088_EB,0.101_EB,0.105_EB,0.109_EB,0.114_EB,0.118_EB,0.112_EB,0.106_EB/)
ABSF(3,17,11:20,6) = (/0.099_EB,0.093_EB,0.078_EB,0.063_EB,0.048_EB,0.033_EB,0.042_EB,0.051_EB,0.061_EB,0.070_EB/)
ABSF(3,17,21:30,6) = (/0.080_EB,0.091_EB,0.101_EB,0.111_EB,0.121_EB,0.128_EB,0.135_EB,0.141_EB,0.148_EB,0.155_EB/)
ABSF(3,18,0:10,1) = (/0.000_EB,0.054_EB,0.068_EB,0.076_EB,0.085_EB,0.081_EB,0.078_EB,0.075_EB,0.072_EB,0.068_EB,0.063_EB/)
ABSF(3,18,11:20,1) = (/0.059_EB,0.055_EB,0.060_EB,0.065_EB,0.070_EB,0.075_EB,0.075_EB,0.076_EB,0.076_EB,0.077_EB/)
ABSF(3,18,21:30,1) = (/0.078_EB,0.079_EB,0.080_EB,0.081_EB,0.082_EB,0.082_EB,0.083_EB,0.084_EB,0.085_EB,0.086_EB/)
ABSF(3,18,0:10,2) = (/0.000_EB,0.055_EB,0.069_EB,0.075_EB,0.081_EB,0.078_EB,0.076_EB,0.073_EB,0.071_EB,0.067_EB,0.062_EB/)
ABSF(3,18,11:20,2) = (/0.058_EB,0.053_EB,0.058_EB,0.062_EB,0.067_EB,0.071_EB,0.073_EB,0.076_EB,0.078_EB,0.080_EB/)
ABSF(3,18,21:30,2) = (/0.081_EB,0.082_EB,0.083_EB,0.084_EB,0.085_EB,0.086_EB,0.086_EB,0.086_EB,0.086_EB,0.087_EB/)
ABSF(3,18,0:10,3) = (/0.000_EB,0.041_EB,0.064_EB,0.075_EB,0.085_EB,0.082_EB,0.080_EB,0.077_EB,0.074_EB,0.069_EB,0.065_EB/)
ABSF(3,18,11:20,3) = (/0.060_EB,0.056_EB,0.060_EB,0.064_EB,0.069_EB,0.073_EB,0.075_EB,0.077_EB,0.078_EB,0.080_EB/)
ABSF(3,18,21:30,3) = (/0.081_EB,0.081_EB,0.082_EB,0.083_EB,0.083_EB,0.084_EB,0.085_EB,0.086_EB,0.087_EB,0.089_EB/)
ABSF(3,18,0:10,4) = (/0.000_EB,0.039_EB,0.061_EB,0.072_EB,0.082_EB,0.086_EB,0.091_EB,0.095_EB,0.099_EB,0.095_EB,0.091_EB/)
ABSF(3,18,11:20,4) = (/0.087_EB,0.083_EB,0.084_EB,0.084_EB,0.085_EB,0.085_EB,0.084_EB,0.083_EB,0.082_EB,0.080_EB/)
ABSF(3,18,21:30,4) = (/0.091_EB,0.103_EB,0.114_EB,0.125_EB,0.136_EB,0.134_EB,0.131_EB,0.129_EB,0.126_EB,0.124_EB/)
ABSF(3,18,0:10,5) = (/0.000_EB,0.038_EB,0.072_EB,0.086_EB,0.100_EB,0.105_EB,0.111_EB,0.116_EB,0.121_EB,0.117_EB,0.112_EB/)
ABSF(3,18,11:20,5) = (/0.107_EB,0.103_EB,0.101_EB,0.099_EB,0.096_EB,0.094_EB,0.091_EB,0.089_EB,0.086_EB,0.083_EB/)
ABSF(3,18,21:30,5) = (/0.095_EB,0.108_EB,0.121_EB,0.134_EB,0.146_EB,0.148_EB,0.150_EB,0.151_EB,0.153_EB,0.155_EB/)
ABSF(3,18,0:10,6) = (/0.000_EB,0.040_EB,0.078_EB,0.091_EB,0.105_EB,0.109_EB,0.113_EB,0.117_EB,0.122_EB,0.115_EB,0.109_EB/)
ABSF(3,18,11:20,6) = (/0.103_EB,0.097_EB,0.081_EB,0.065_EB,0.050_EB,0.034_EB,0.044_EB,0.053_EB,0.062_EB,0.072_EB/)
ABSF(3,18,21:30,6) = (/0.084_EB,0.097_EB,0.110_EB,0.123_EB,0.135_EB,0.140_EB,0.145_EB,0.149_EB,0.154_EB,0.158_EB/)
ABSF(3,19,0:10,1) = (/0.000_EB,0.056_EB,0.070_EB,0.079_EB,0.088_EB,0.085_EB,0.081_EB,0.078_EB,0.075_EB,0.070_EB,0.066_EB/)
ABSF(3,19,11:20,1) = (/0.061_EB,0.056_EB,0.062_EB,0.067_EB,0.072_EB,0.077_EB,0.078_EB,0.078_EB,0.079_EB,0.079_EB/)
ABSF(3,19,21:30,1) = (/0.080_EB,0.081_EB,0.082_EB,0.083_EB,0.084_EB,0.085_EB,0.086_EB,0.086_EB,0.087_EB,0.088_EB/)
ABSF(3,19,0:10,2) = (/0.000_EB,0.057_EB,0.071_EB,0.078_EB,0.084_EB,0.081_EB,0.079_EB,0.076_EB,0.073_EB,0.069_EB,0.064_EB/)
ABSF(3,19,11:20,2) = (/0.060_EB,0.055_EB,0.060_EB,0.064_EB,0.069_EB,0.073_EB,0.076_EB,0.078_EB,0.080_EB,0.083_EB/)
ABSF(3,19,21:30,2) = (/0.084_EB,0.085_EB,0.086_EB,0.087_EB,0.088_EB,0.088_EB,0.088_EB,0.089_EB,0.089_EB,0.089_EB/)
ABSF(3,19,0:10,3) = (/0.000_EB,0.043_EB,0.067_EB,0.078_EB,0.088_EB,0.085_EB,0.082_EB,0.078_EB,0.075_EB,0.071_EB,0.067_EB/)
ABSF(3,19,11:20,3) = (/0.062_EB,0.058_EB,0.062_EB,0.066_EB,0.070_EB,0.075_EB,0.077_EB,0.079_EB,0.081_EB,0.083_EB/)
ABSF(3,19,21:30,3) = (/0.083_EB,0.084_EB,0.085_EB,0.085_EB,0.086_EB,0.087_EB,0.088_EB,0.089_EB,0.090_EB,0.091_EB/)
ABSF(3,19,0:10,4) = (/0.000_EB,0.040_EB,0.063_EB,0.073_EB,0.083_EB,0.087_EB,0.092_EB,0.096_EB,0.100_EB,0.096_EB,0.092_EB/)
ABSF(3,19,11:20,4) = (/0.088_EB,0.084_EB,0.085_EB,0.086_EB,0.087_EB,0.088_EB,0.087_EB,0.086_EB,0.085_EB,0.084_EB/)
ABSF(3,19,21:30,4) = (/0.095_EB,0.106_EB,0.116_EB,0.127_EB,0.138_EB,0.135_EB,0.133_EB,0.130_EB,0.128_EB,0.125_EB/)
ABSF(3,19,0:10,5) = (/0.000_EB,0.039_EB,0.074_EB,0.088_EB,0.102_EB,0.108_EB,0.113_EB,0.118_EB,0.123_EB,0.119_EB,0.115_EB/)
ABSF(3,19,11:20,5) = (/0.110_EB,0.106_EB,0.104_EB,0.102_EB,0.100_EB,0.098_EB,0.094_EB,0.091_EB,0.088_EB,0.085_EB/)
ABSF(3,19,21:30,5) = (/0.097_EB,0.110_EB,0.123_EB,0.135_EB,0.148_EB,0.150_EB,0.151_EB,0.152_EB,0.154_EB,0.155_EB/)
ABSF(3,19,0:10,6) = (/0.000_EB,0.041_EB,0.080_EB,0.094_EB,0.109_EB,0.113_EB,0.117_EB,0.121_EB,0.126_EB,0.119_EB,0.113_EB/)
ABSF(3,19,11:20,6) = (/0.106_EB,0.100_EB,0.084_EB,0.068_EB,0.052_EB,0.036_EB,0.045_EB,0.055_EB,0.064_EB,0.073_EB/)
ABSF(3,19,21:30,6) = (/0.088_EB,0.104_EB,0.119_EB,0.134_EB,0.149_EB,0.152_EB,0.154_EB,0.157_EB,0.159_EB,0.162_EB/)
ABSF(3,20,0:10,1) = (/0.000_EB,0.058_EB,0.073_EB,0.082_EB,0.092_EB,0.088_EB,0.085_EB,0.081_EB,0.078_EB,0.073_EB,0.068_EB/)
ABSF(3,20,11:20,1) = (/0.063_EB,0.058_EB,0.064_EB,0.069_EB,0.074_EB,0.079_EB,0.080_EB,0.081_EB,0.081_EB,0.082_EB/)
ABSF(3,20,21:30,1) = (/0.083_EB,0.083_EB,0.084_EB,0.085_EB,0.086_EB,0.087_EB,0.088_EB,0.089_EB,0.090_EB,0.091_EB/)
ABSF(3,20,0:10,2) = (/0.000_EB,0.059_EB,0.074_EB,0.081_EB,0.088_EB,0.085_EB,0.082_EB,0.079_EB,0.075_EB,0.071_EB,0.066_EB/)
ABSF(3,20,11:20,2) = (/0.062_EB,0.057_EB,0.062_EB,0.066_EB,0.071_EB,0.076_EB,0.078_EB,0.080_EB,0.083_EB,0.085_EB/)
ABSF(3,20,21:30,2) = (/0.086_EB,0.087_EB,0.088_EB,0.089_EB,0.090_EB,0.091_EB,0.091_EB,0.091_EB,0.091_EB,0.092_EB/)
ABSF(3,20,0:10,3) = (/0.000_EB,0.044_EB,0.070_EB,0.081_EB,0.091_EB,0.087_EB,0.084_EB,0.080_EB,0.077_EB,0.073_EB,0.069_EB/)
ABSF(3,20,11:20,3) = (/0.064_EB,0.060_EB,0.064_EB,0.068_EB,0.072_EB,0.076_EB,0.079_EB,0.081_EB,0.083_EB,0.085_EB/)
ABSF(3,20,21:30,3) = (/0.086_EB,0.087_EB,0.087_EB,0.088_EB,0.089_EB,0.090_EB,0.091_EB,0.092_EB,0.093_EB,0.094_EB/)
ABSF(3,20,0:10,4) = (/0.000_EB,0.042_EB,0.064_EB,0.074_EB,0.084_EB,0.088_EB,0.093_EB,0.097_EB,0.101_EB,0.097_EB,0.093_EB/)
ABSF(3,20,11:20,4) = (/0.089_EB,0.085_EB,0.086_EB,0.087_EB,0.089_EB,0.090_EB,0.089_EB,0.089_EB,0.088_EB,0.088_EB/)
ABSF(3,20,21:30,4) = (/0.098_EB,0.109_EB,0.119_EB,0.130_EB,0.140_EB,0.137_EB,0.134_EB,0.132_EB,0.129_EB,0.126_EB/)
ABSF(3,20,0:10,5) = (/0.000_EB,0.041_EB,0.076_EB,0.091_EB,0.105_EB,0.110_EB,0.115_EB,0.120_EB,0.125_EB,0.121_EB,0.117_EB/)
ABSF(3,20,11:20,5) = (/0.113_EB,0.109_EB,0.107_EB,0.105_EB,0.103_EB,0.101_EB,0.097_EB,0.094_EB,0.090_EB,0.087_EB/)
ABSF(3,20,21:30,5) = (/0.099_EB,0.112_EB,0.125_EB,0.137_EB,0.150_EB,0.151_EB,0.152_EB,0.153_EB,0.155_EB,0.156_EB/)
ABSF(3,20,0:10,6) = (/0.000_EB,0.042_EB,0.082_EB,0.097_EB,0.112_EB,0.117_EB,0.121_EB,0.125_EB,0.129_EB,0.123_EB,0.116_EB/)
ABSF(3,20,11:20,6) = (/0.110_EB,0.103_EB,0.087_EB,0.070_EB,0.054_EB,0.038_EB,0.047_EB,0.056_EB,0.065_EB,0.075_EB/)
ABSF(3,20,21:30,6) = (/0.092_EB,0.110_EB,0.128_EB,0.145_EB,0.163_EB,0.164_EB,0.164_EB,0.165_EB,0.165_EB,0.166_EB/)
ABSF(4,0,0:10,1) = (/0.000_EB,0.010_EB,0.018_EB,0.020_EB,0.021_EB,0.023_EB,0.024_EB,0.026_EB,0.027_EB,0.027_EB,0.027_EB/)
ABSF(4,0,11:20,1) = (/0.027_EB,0.027_EB,0.028_EB,0.029_EB,0.031_EB,0.032_EB,0.033_EB,0.034_EB,0.034_EB,0.035_EB/)
ABSF(4,0,21:30,1) = (/0.035_EB,0.035_EB,0.035_EB,0.035_EB,0.034_EB,0.035_EB,0.036_EB,0.037_EB,0.038_EB,0.039_EB/)
ABSF(4,0,0:10,2) = (/0.000_EB,0.013_EB,0.019_EB,0.022_EB,0.025_EB,0.028_EB,0.032_EB,0.036_EB,0.040_EB,0.040_EB,0.040_EB/)
ABSF(4,0,11:20,2) = (/0.040_EB,0.040_EB,0.040_EB,0.040_EB,0.040_EB,0.040_EB,0.038_EB,0.037_EB,0.036_EB,0.035_EB/)
ABSF(4,0,21:30,2) = (/0.032_EB,0.029_EB,0.027_EB,0.024_EB,0.021_EB,0.024_EB,0.026_EB,0.028_EB,0.031_EB,0.033_EB/)
ABSF(4,0,0:10,3) = (/0.000_EB,0.014_EB,0.022_EB,0.029_EB,0.037_EB,0.026_EB,0.015_EB,0.004_EB,-0.007_EB,0.006_EB,0.020_EB/)
ABSF(4,0,11:20,3) = (/0.034_EB,0.048_EB,0.054_EB,0.061_EB,0.068_EB,0.074_EB,0.064_EB,0.055_EB,0.045_EB,0.035_EB/)
ABSF(4,0,21:30,3) = (/0.031_EB,0.027_EB,0.023_EB,0.019_EB,0.015_EB,0.021_EB,0.028_EB,0.035_EB,0.042_EB,0.049_EB/)
ABSF(4,0,0:10,4) = (/0.000_EB,0.014_EB,0.028_EB,0.034_EB,0.041_EB,0.048_EB,0.054_EB,0.061_EB,0.067_EB,0.068_EB,0.068_EB/)
ABSF(4,0,11:20,4) = (/0.069_EB,0.069_EB,0.067_EB,0.065_EB,0.063_EB,0.061_EB,0.068_EB,0.074_EB,0.080_EB,0.086_EB/)
ABSF(4,0,21:30,4) = (/0.068_EB,0.051_EB,0.033_EB,0.016_EB,-0.002_EB,0.010_EB,0.021_EB,0.032_EB,0.044_EB,0.055_EB/)
ABSF(4,0,0:10,5) = (/0.000_EB,0.014_EB,0.033_EB,0.041_EB,0.049_EB,0.056_EB,0.063_EB,0.070_EB,0.077_EB,0.076_EB,0.076_EB/)
ABSF(4,0,11:20,5) = (/0.076_EB,0.076_EB,0.056_EB,0.037_EB,0.018_EB,-0.002_EB,0.020_EB,0.042_EB,0.064_EB,0.086_EB/)
ABSF(4,0,21:30,5) = (/0.085_EB,0.083_EB,0.082_EB,0.080_EB,0.078_EB,0.087_EB,0.095_EB,0.104_EB,0.112_EB,0.121_EB/)
ABSF(4,0,0:10,6) = (/0.000_EB,0.017_EB,0.038_EB,0.042_EB,0.046_EB,0.054_EB,0.062_EB,0.070_EB,0.079_EB,0.083_EB,0.088_EB/)
ABSF(4,0,11:20,6) = (/0.092_EB,0.097_EB,0.082_EB,0.068_EB,0.053_EB,0.038_EB,0.052_EB,0.066_EB,0.080_EB,0.094_EB/)
ABSF(4,0,21:30,6) = (/0.096_EB,0.097_EB,0.098_EB,0.099_EB,0.101_EB,0.106_EB,0.111_EB,0.116_EB,0.121_EB,0.126_EB/)
ABSF(4,1,0:10,1) = (/0.000_EB,0.015_EB,0.019_EB,0.021_EB,0.024_EB,0.025_EB,0.026_EB,0.027_EB,0.028_EB,0.027_EB,0.027_EB/)
ABSF(4,1,11:20,1) = (/0.027_EB,0.026_EB,0.028_EB,0.029_EB,0.031_EB,0.033_EB,0.033_EB,0.034_EB,0.035_EB,0.036_EB/)
ABSF(4,1,21:30,1) = (/0.036_EB,0.037_EB,0.037_EB,0.037_EB,0.037_EB,0.037_EB,0.037_EB,0.037_EB,0.037_EB,0.037_EB/)
ABSF(4,1,0:10,2) = (/0.000_EB,0.012_EB,0.021_EB,0.025_EB,0.029_EB,0.032_EB,0.035_EB,0.038_EB,0.041_EB,0.042_EB,0.043_EB/)
ABSF(4,1,11:20,2) = (/0.044_EB,0.045_EB,0.044_EB,0.044_EB,0.043_EB,0.043_EB,0.041_EB,0.039_EB,0.038_EB,0.036_EB/)
ABSF(4,1,21:30,2) = (/0.033_EB,0.031_EB,0.028_EB,0.026_EB,0.023_EB,0.030_EB,0.036_EB,0.043_EB,0.049_EB,0.055_EB/)
ABSF(4,1,0:10,3) = (/0.000_EB,0.014_EB,0.023_EB,0.030_EB,0.037_EB,0.042_EB,0.046_EB,0.050_EB,0.054_EB,0.053_EB,0.051_EB/)
ABSF(4,1,11:20,3) = (/0.050_EB,0.049_EB,0.055_EB,0.060_EB,0.065_EB,0.071_EB,0.063_EB,0.055_EB,0.047_EB,0.040_EB/)
ABSF(4,1,21:30,3) = (/0.035_EB,0.029_EB,0.024_EB,0.019_EB,0.014_EB,0.026_EB,0.038_EB,0.049_EB,0.061_EB,0.073_EB/)
ABSF(4,1,0:10,4) = (/0.000_EB,0.016_EB,0.030_EB,0.038_EB,0.047_EB,0.052_EB,0.058_EB,0.064_EB,0.069_EB,0.070_EB,0.071_EB/)
ABSF(4,1,11:20,4) = (/0.071_EB,0.072_EB,0.072_EB,0.073_EB,0.073_EB,0.073_EB,0.062_EB,0.051_EB,0.039_EB,0.028_EB/)
ABSF(4,1,21:30,4) = (/0.022_EB,0.016_EB,0.009_EB,0.003_EB,-0.003_EB,0.017_EB,0.036_EB,0.055_EB,0.075_EB,0.094_EB/)
ABSF(4,1,0:10,5) = (/0.000_EB,0.015_EB,0.034_EB,0.041_EB,0.047_EB,0.056_EB,0.065_EB,0.073_EB,0.082_EB,0.074_EB,0.066_EB/)
ABSF(4,1,11:20,5) = (/0.059_EB,0.051_EB,0.053_EB,0.054_EB,0.056_EB,0.058_EB,0.069_EB,0.080_EB,0.091_EB,0.102_EB/)
ABSF(4,1,21:30,5) = (/0.096_EB,0.090_EB,0.085_EB,0.079_EB,0.074_EB,0.083_EB,0.093_EB,0.103_EB,0.112_EB,0.122_EB/)
ABSF(4,1,0:10,6) = (/0.000_EB,0.019_EB,0.036_EB,0.041_EB,0.047_EB,0.055_EB,0.062_EB,0.070_EB,0.078_EB,0.080_EB,0.082_EB/)
ABSF(4,1,11:20,6) = (/0.085_EB,0.087_EB,0.086_EB,0.084_EB,0.083_EB,0.082_EB,0.089_EB,0.095_EB,0.102_EB,0.108_EB/)
ABSF(4,1,21:30,6) = (/0.103_EB,0.097_EB,0.091_EB,0.086_EB,0.080_EB,0.083_EB,0.085_EB,0.088_EB,0.091_EB,0.093_EB/)
ABSF(4,2,0:10,1) = (/0.000_EB,0.022_EB,0.024_EB,0.026_EB,0.027_EB,0.028_EB,0.028_EB,0.029_EB,0.029_EB,0.029_EB,0.029_EB/)
ABSF(4,2,11:20,1) = (/0.028_EB,0.028_EB,0.029_EB,0.031_EB,0.032_EB,0.034_EB,0.035_EB,0.036_EB,0.036_EB,0.037_EB/)
ABSF(4,2,21:30,1) = (/0.037_EB,0.037_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB,0.039_EB,0.040_EB,0.040_EB,0.041_EB/)
ABSF(4,2,0:10,2) = (/0.000_EB,0.013_EB,0.023_EB,0.027_EB,0.032_EB,0.033_EB,0.034_EB,0.035_EB,0.037_EB,0.036_EB,0.035_EB/)
ABSF(4,2,11:20,2) = (/0.035_EB,0.034_EB,0.036_EB,0.037_EB,0.039_EB,0.041_EB,0.041_EB,0.041_EB,0.041_EB,0.041_EB/)
ABSF(4,2,21:30,2) = (/0.039_EB,0.037_EB,0.035_EB,0.033_EB,0.031_EB,0.036_EB,0.041_EB,0.047_EB,0.052_EB,0.057_EB/)
ABSF(4,2,0:10,3) = (/0.000_EB,0.017_EB,0.025_EB,0.031_EB,0.037_EB,0.027_EB,0.016_EB,0.006_EB,-0.004_EB,0.009_EB,0.022_EB/)
ABSF(4,2,11:20,3) = (/0.035_EB,0.048_EB,0.056_EB,0.064_EB,0.072_EB,0.080_EB,0.072_EB,0.064_EB,0.056_EB,0.048_EB/)
ABSF(4,2,21:30,3) = (/0.039_EB,0.030_EB,0.021_EB,0.012_EB,0.002_EB,0.017_EB,0.032_EB,0.047_EB,0.061_EB,0.076_EB/)
ABSF(4,2,0:10,4) = (/0.000_EB,0.016_EB,0.032_EB,0.039_EB,0.046_EB,0.052_EB,0.058_EB,0.064_EB,0.069_EB,0.071_EB,0.072_EB/)
ABSF(4,2,11:20,4) = (/0.073_EB,0.074_EB,0.078_EB,0.081_EB,0.085_EB,0.088_EB,0.086_EB,0.083_EB,0.081_EB,0.078_EB/)
ABSF(4,2,21:30,4) = (/0.075_EB,0.072_EB,0.070_EB,0.067_EB,0.064_EB,0.069_EB,0.074_EB,0.078_EB,0.083_EB,0.088_EB/)
ABSF(4,2,0:10,5) = (/0.000_EB,0.015_EB,0.038_EB,0.043_EB,0.049_EB,0.056_EB,0.063_EB,0.070_EB,0.077_EB,0.066_EB,0.054_EB/)
ABSF(4,2,11:20,5) = (/0.042_EB,0.031_EB,0.034_EB,0.036_EB,0.039_EB,0.042_EB,0.048_EB,0.054_EB,0.060_EB,0.066_EB/)
ABSF(4,2,21:30,5) = (/0.076_EB,0.085_EB,0.095_EB,0.105_EB,0.115_EB,0.112_EB,0.108_EB,0.105_EB,0.101_EB,0.098_EB/)
ABSF(4,2,0:10,6) = (/0.000_EB,0.018_EB,0.039_EB,0.045_EB,0.051_EB,0.059_EB,0.067_EB,0.075_EB,0.083_EB,0.078_EB,0.073_EB/)
ABSF(4,2,11:20,6) = (/0.068_EB,0.064_EB,0.052_EB,0.040_EB,0.029_EB,0.017_EB,0.038_EB,0.059_EB,0.080_EB,0.101_EB/)
ABSF(4,2,21:30,6) = (/0.091_EB,0.082_EB,0.073_EB,0.064_EB,0.055_EB,0.069_EB,0.082_EB,0.096_EB,0.110_EB,0.124_EB/)
ABSF(4,3,0:10,1) = (/0.000_EB,0.024_EB,0.026_EB,0.028_EB,0.030_EB,0.030_EB,0.031_EB,0.031_EB,0.031_EB,0.030_EB,0.030_EB/)
ABSF(4,3,11:20,1) = (/0.029_EB,0.028_EB,0.030_EB,0.032_EB,0.034_EB,0.035_EB,0.036_EB,0.037_EB,0.038_EB,0.039_EB/)
ABSF(4,3,21:30,1) = (/0.040_EB,0.040_EB,0.040_EB,0.040_EB,0.040_EB,0.040_EB,0.041_EB,0.041_EB,0.041_EB,0.042_EB/)
ABSF(4,3,0:10,2) = (/0.000_EB,0.015_EB,0.024_EB,0.028_EB,0.032_EB,0.033_EB,0.034_EB,0.035_EB,0.036_EB,0.035_EB,0.034_EB/)
ABSF(4,3,11:20,2) = (/0.033_EB,0.032_EB,0.034_EB,0.036_EB,0.038_EB,0.040_EB,0.040_EB,0.041_EB,0.041_EB,0.041_EB/)
ABSF(4,3,21:30,2) = (/0.040_EB,0.039_EB,0.038_EB,0.037_EB,0.035_EB,0.039_EB,0.042_EB,0.046_EB,0.049_EB,0.053_EB/)
ABSF(4,3,0:10,3) = (/0.000_EB,0.018_EB,0.027_EB,0.034_EB,0.040_EB,0.034_EB,0.029_EB,0.024_EB,0.018_EB,0.026_EB,0.034_EB/)
ABSF(4,3,11:20,3) = (/0.042_EB,0.050_EB,0.057_EB,0.064_EB,0.071_EB,0.078_EB,0.071_EB,0.064_EB,0.056_EB,0.049_EB/)
ABSF(4,3,21:30,3) = (/0.041_EB,0.033_EB,0.025_EB,0.018_EB,0.010_EB,0.023_EB,0.037_EB,0.050_EB,0.064_EB,0.077_EB/)
ABSF(4,3,0:10,4) = (/0.000_EB,0.018_EB,0.033_EB,0.041_EB,0.049_EB,0.054_EB,0.060_EB,0.065_EB,0.071_EB,0.071_EB,0.072_EB/)
ABSF(4,3,11:20,4) = (/0.072_EB,0.072_EB,0.075_EB,0.078_EB,0.082_EB,0.085_EB,0.084_EB,0.083_EB,0.081_EB,0.080_EB/)
ABSF(4,3,21:30,4) = (/0.076_EB,0.073_EB,0.069_EB,0.065_EB,0.061_EB,0.067_EB,0.074_EB,0.081_EB,0.087_EB,0.094_EB/)
ABSF(4,3,0:10,5) = (/0.000_EB,0.017_EB,0.040_EB,0.045_EB,0.050_EB,0.058_EB,0.065_EB,0.073_EB,0.080_EB,0.073_EB,0.065_EB/)
ABSF(4,3,11:20,5) = (/0.058_EB,0.050_EB,0.048_EB,0.046_EB,0.044_EB,0.042_EB,0.050_EB,0.057_EB,0.065_EB,0.072_EB/)
ABSF(4,3,21:30,5) = (/0.080_EB,0.087_EB,0.095_EB,0.103_EB,0.110_EB,0.109_EB,0.108_EB,0.107_EB,0.106_EB,0.105_EB/)
ABSF(4,3,0:10,6) = (/0.000_EB,0.019_EB,0.043_EB,0.047_EB,0.052_EB,0.060_EB,0.069_EB,0.077_EB,0.086_EB,0.082_EB,0.078_EB/)
ABSF(4,3,11:20,6) = (/0.074_EB,0.070_EB,0.061_EB,0.052_EB,0.042_EB,0.033_EB,0.045_EB,0.057_EB,0.069_EB,0.080_EB/)
ABSF(4,3,21:30,6) = (/0.079_EB,0.077_EB,0.076_EB,0.075_EB,0.073_EB,0.084_EB,0.094_EB,0.104_EB,0.115_EB,0.125_EB/)
ABSF(4,4,0:10,1) = (/0.000_EB,0.025_EB,0.027_EB,0.030_EB,0.033_EB,0.033_EB,0.033_EB,0.033_EB,0.033_EB,0.032_EB,0.031_EB/)
ABSF(4,4,11:20,1) = (/0.029_EB,0.028_EB,0.030_EB,0.032_EB,0.035_EB,0.037_EB,0.038_EB,0.039_EB,0.040_EB,0.042_EB/)
ABSF(4,4,21:30,1) = (/0.042_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB/)
ABSF(4,4,0:10,2) = (/0.000_EB,0.017_EB,0.025_EB,0.029_EB,0.033_EB,0.033_EB,0.034_EB,0.035_EB,0.036_EB,0.035_EB,0.033_EB/)
ABSF(4,4,11:20,2) = (/0.032_EB,0.031_EB,0.033_EB,0.035_EB,0.037_EB,0.039_EB,0.040_EB,0.040_EB,0.041_EB,0.042_EB/)
ABSF(4,4,21:30,2) = (/0.041_EB,0.041_EB,0.041_EB,0.040_EB,0.040_EB,0.042_EB,0.043_EB,0.045_EB,0.047_EB,0.049_EB/)
ABSF(4,4,0:10,3) = (/0.000_EB,0.019_EB,0.029_EB,0.036_EB,0.043_EB,0.042_EB,0.042_EB,0.041_EB,0.040_EB,0.043_EB,0.047_EB/)
ABSF(4,4,11:20,3) = (/0.050_EB,0.053_EB,0.059_EB,0.065_EB,0.071_EB,0.077_EB,0.070_EB,0.063_EB,0.056_EB,0.049_EB/)
ABSF(4,4,21:30,3) = (/0.043_EB,0.036_EB,0.030_EB,0.024_EB,0.017_EB,0.029_EB,0.042_EB,0.054_EB,0.066_EB,0.079_EB/)
ABSF(4,4,0:10,4) = (/0.000_EB,0.019_EB,0.035_EB,0.043_EB,0.052_EB,0.057_EB,0.062_EB,0.067_EB,0.073_EB,0.072_EB,0.071_EB/)
ABSF(4,4,11:20,4) = (/0.070_EB,0.069_EB,0.072_EB,0.075_EB,0.078_EB,0.081_EB,0.081_EB,0.082_EB,0.082_EB,0.083_EB/)
ABSF(4,4,21:30,4) = (/0.078_EB,0.073_EB,0.068_EB,0.063_EB,0.058_EB,0.066_EB,0.074_EB,0.083_EB,0.091_EB,0.099_EB/)
ABSF(4,4,0:10,5) = (/0.000_EB,0.018_EB,0.042_EB,0.047_EB,0.052_EB,0.060_EB,0.067_EB,0.075_EB,0.083_EB,0.080_EB,0.077_EB/)
ABSF(4,4,11:20,5) = (/0.073_EB,0.070_EB,0.063_EB,0.056_EB,0.049_EB,0.042_EB,0.051_EB,0.060_EB,0.070_EB,0.079_EB/)
ABSF(4,4,21:30,5) = (/0.084_EB,0.089_EB,0.095_EB,0.100_EB,0.105_EB,0.107_EB,0.108_EB,0.109_EB,0.110_EB,0.112_EB/)
ABSF(4,4,0:10,6) = (/0.000_EB,0.020_EB,0.047_EB,0.050_EB,0.053_EB,0.062_EB,0.070_EB,0.079_EB,0.088_EB,0.085_EB,0.083_EB/)
ABSF(4,4,11:20,6) = (/0.080_EB,0.077_EB,0.070_EB,0.063_EB,0.056_EB,0.049_EB,0.052_EB,0.055_EB,0.057_EB,0.060_EB/)
ABSF(4,4,21:30,6) = (/0.066_EB,0.073_EB,0.079_EB,0.085_EB,0.091_EB,0.098_EB,0.105_EB,0.113_EB,0.120_EB,0.127_EB/)
ABSF(4,5,0:10,1) = (/0.000_EB,0.027_EB,0.029_EB,0.032_EB,0.036_EB,0.036_EB,0.036_EB,0.036_EB,0.035_EB,0.033_EB,0.032_EB/)
ABSF(4,5,11:20,1) = (/0.030_EB,0.028_EB,0.030_EB,0.033_EB,0.036_EB,0.038_EB,0.040_EB,0.041_EB,0.042_EB,0.044_EB/)
ABSF(4,5,21:30,1) = (/0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.043_EB,0.043_EB,0.043_EB/)
ABSF(4,5,0:10,2) = (/0.000_EB,0.018_EB,0.026_EB,0.030_EB,0.033_EB,0.034_EB,0.034_EB,0.035_EB,0.036_EB,0.034_EB,0.032_EB/)
ABSF(4,5,11:20,2) = (/0.031_EB,0.029_EB,0.031_EB,0.034_EB,0.036_EB,0.038_EB,0.039_EB,0.040_EB,0.041_EB,0.042_EB/)
ABSF(4,5,21:30,2) = (/0.043_EB,0.043_EB,0.043_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.045_EB,0.045_EB,0.045_EB/)
ABSF(4,5,0:10,3) = (/0.000_EB,0.019_EB,0.031_EB,0.039_EB,0.046_EB,0.050_EB,0.054_EB,0.058_EB,0.062_EB,0.061_EB,0.059_EB/)
ABSF(4,5,11:20,3) = (/0.058_EB,0.056_EB,0.061_EB,0.065_EB,0.070_EB,0.075_EB,0.069_EB,0.062_EB,0.056_EB,0.050_EB/)
ABSF(4,5,21:30,3) = (/0.045_EB,0.040_EB,0.035_EB,0.030_EB,0.025_EB,0.036_EB,0.047_EB,0.058_EB,0.069_EB,0.080_EB/)
ABSF(4,5,0:10,4) = (/0.000_EB,0.020_EB,0.036_EB,0.046_EB,0.055_EB,0.060_EB,0.064_EB,0.069_EB,0.074_EB,0.072_EB,0.071_EB/)
ABSF(4,5,11:20,4) = (/0.069_EB,0.067_EB,0.070_EB,0.072_EB,0.075_EB,0.077_EB,0.079_EB,0.081_EB,0.083_EB,0.085_EB/)
ABSF(4,5,21:30,4) = (/0.079_EB,0.073_EB,0.067_EB,0.061_EB,0.055_EB,0.065_EB,0.075_EB,0.085_EB,0.095_EB,0.105_EB/)
ABSF(4,5,0:10,5) = (/0.000_EB,0.020_EB,0.044_EB,0.049_EB,0.054_EB,0.062_EB,0.070_EB,0.078_EB,0.085_EB,0.087_EB,0.088_EB/)
ABSF(4,5,11:20,5) = (/0.089_EB,0.090_EB,0.078_EB,0.066_EB,0.054_EB,0.042_EB,0.053_EB,0.064_EB,0.075_EB,0.085_EB/)
ABSF(4,5,21:30,5) = (/0.088_EB,0.091_EB,0.094_EB,0.097_EB,0.100_EB,0.104_EB,0.108_EB,0.111_EB,0.115_EB,0.118_EB/)
ABSF(4,5,0:10,6) = (/0.000_EB,0.022_EB,0.051_EB,0.052_EB,0.053_EB,0.063_EB,0.072_EB,0.081_EB,0.091_EB,0.089_EB,0.087_EB/)
ABSF(4,5,11:20,6) = (/0.086_EB,0.084_EB,0.079_EB,0.074_EB,0.069_EB,0.065_EB,0.059_EB,0.052_EB,0.046_EB,0.040_EB/)
ABSF(4,5,21:30,6) = (/0.054_EB,0.068_EB,0.082_EB,0.096_EB,0.110_EB,0.113_EB,0.117_EB,0.121_EB,0.124_EB,0.128_EB/)
ABSF(4,6,0:10,1) = (/0.000_EB,0.029_EB,0.031_EB,0.035_EB,0.039_EB,0.039_EB,0.038_EB,0.038_EB,0.037_EB,0.035_EB,0.033_EB/)
ABSF(4,6,11:20,1) = (/0.031_EB,0.029_EB,0.032_EB,0.034_EB,0.037_EB,0.039_EB,0.041_EB,0.043_EB,0.044_EB,0.046_EB/)
ABSF(4,6,21:30,1) = (/0.046_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB,0.045_EB,0.045_EB,0.045_EB,0.045_EB/)
ABSF(4,6,0:10,2) = (/0.000_EB,0.021_EB,0.029_EB,0.033_EB,0.037_EB,0.037_EB,0.037_EB,0.037_EB,0.038_EB,0.036_EB,0.034_EB/)
ABSF(4,6,11:20,2) = (/0.032_EB,0.030_EB,0.032_EB,0.035_EB,0.037_EB,0.039_EB,0.040_EB,0.041_EB,0.042_EB,0.043_EB/)
ABSF(4,6,21:30,2) = (/0.044_EB,0.044_EB,0.045_EB,0.045_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB/)
ABSF(4,6,0:10,3) = (/0.000_EB,0.020_EB,0.033_EB,0.040_EB,0.047_EB,0.051_EB,0.054_EB,0.057_EB,0.060_EB,0.059_EB,0.057_EB/)
ABSF(4,6,11:20,3) = (/0.056_EB,0.054_EB,0.059_EB,0.063_EB,0.068_EB,0.072_EB,0.067_EB,0.062_EB,0.057_EB,0.052_EB/)
ABSF(4,6,21:30,3) = (/0.048_EB,0.043_EB,0.038_EB,0.034_EB,0.029_EB,0.039_EB,0.050_EB,0.060_EB,0.071_EB,0.081_EB/)
ABSF(4,6,0:10,4) = (/0.000_EB,0.021_EB,0.038_EB,0.047_EB,0.057_EB,0.062_EB,0.066_EB,0.071_EB,0.076_EB,0.074_EB,0.072_EB/)
ABSF(4,6,11:20,4) = (/0.071_EB,0.069_EB,0.071_EB,0.074_EB,0.076_EB,0.079_EB,0.080_EB,0.081_EB,0.083_EB,0.084_EB/)
ABSF(4,6,21:30,4) = (/0.077_EB,0.070_EB,0.063_EB,0.055_EB,0.048_EB,0.059_EB,0.071_EB,0.082_EB,0.093_EB,0.104_EB/)
ABSF(4,6,0:10,5) = (/0.000_EB,0.021_EB,0.046_EB,0.051_EB,0.055_EB,0.063_EB,0.071_EB,0.079_EB,0.086_EB,0.087_EB,0.088_EB/)
ABSF(4,6,11:20,5) = (/0.089_EB,0.090_EB,0.079_EB,0.067_EB,0.055_EB,0.043_EB,0.052_EB,0.061_EB,0.070_EB,0.079_EB/)
ABSF(4,6,21:30,5) = (/0.083_EB,0.086_EB,0.090_EB,0.094_EB,0.097_EB,0.101_EB,0.105_EB,0.109_EB,0.112_EB,0.116_EB/)
ABSF(4,6,0:10,6) = (/0.000_EB,0.023_EB,0.053_EB,0.055_EB,0.056_EB,0.065_EB,0.074_EB,0.083_EB,0.092_EB,0.090_EB,0.087_EB/)
ABSF(4,6,11:20,6) = (/0.085_EB,0.082_EB,0.079_EB,0.076_EB,0.073_EB,0.070_EB,0.065_EB,0.061_EB,0.056_EB,0.052_EB/)
ABSF(4,6,21:30,6) = (/0.063_EB,0.074_EB,0.086_EB,0.097_EB,0.108_EB,0.112_EB,0.116_EB,0.119_EB,0.123_EB,0.127_EB/)
ABSF(4,7,0:10,1) = (/0.000_EB,0.030_EB,0.033_EB,0.038_EB,0.042_EB,0.041_EB,0.041_EB,0.040_EB,0.039_EB,0.037_EB,0.034_EB/)
ABSF(4,7,11:20,1) = (/0.032_EB,0.030_EB,0.033_EB,0.035_EB,0.038_EB,0.041_EB,0.042_EB,0.044_EB,0.046_EB,0.047_EB/)
ABSF(4,7,21:30,1) = (/0.048_EB,0.048_EB,0.048_EB,0.048_EB,0.048_EB,0.048_EB,0.047_EB,0.047_EB,0.047_EB,0.046_EB/)
ABSF(4,7,0:10,2) = (/0.000_EB,0.023_EB,0.032_EB,0.036_EB,0.040_EB,0.040_EB,0.040_EB,0.040_EB,0.039_EB,0.037_EB,0.035_EB/)
ABSF(4,7,11:20,2) = (/0.033_EB,0.031_EB,0.033_EB,0.036_EB,0.038_EB,0.040_EB,0.041_EB,0.043_EB,0.044_EB,0.045_EB/)
ABSF(4,7,21:30,2) = (/0.045_EB,0.046_EB,0.046_EB,0.047_EB,0.047_EB,0.047_EB,0.047_EB,0.047_EB,0.047_EB,0.047_EB/)
ABSF(4,7,0:10,3) = (/0.000_EB,0.021_EB,0.035_EB,0.041_EB,0.048_EB,0.051_EB,0.053_EB,0.056_EB,0.059_EB,0.057_EB,0.056_EB/)
ABSF(4,7,11:20,3) = (/0.054_EB,0.053_EB,0.057_EB,0.061_EB,0.065_EB,0.069_EB,0.066_EB,0.062_EB,0.058_EB,0.055_EB/)
ABSF(4,7,21:30,3) = (/0.051_EB,0.046_EB,0.042_EB,0.038_EB,0.034_EB,0.043_EB,0.053_EB,0.063_EB,0.072_EB,0.082_EB/)
ABSF(4,7,0:10,4) = (/0.000_EB,0.022_EB,0.039_EB,0.049_EB,0.059_EB,0.064_EB,0.068_EB,0.073_EB,0.078_EB,0.076_EB,0.074_EB/)
ABSF(4,7,11:20,4) = (/0.072_EB,0.071_EB,0.073_EB,0.075_EB,0.077_EB,0.080_EB,0.081_EB,0.082_EB,0.082_EB,0.083_EB/)
ABSF(4,7,21:30,4) = (/0.075_EB,0.067_EB,0.058_EB,0.050_EB,0.042_EB,0.054_EB,0.066_EB,0.079_EB,0.091_EB,0.104_EB/)
ABSF(4,7,0:10,5) = (/0.000_EB,0.022_EB,0.048_EB,0.053_EB,0.057_EB,0.065_EB,0.072_EB,0.080_EB,0.087_EB,0.088_EB,0.089_EB/)
ABSF(4,7,11:20,5) = (/0.090_EB,0.090_EB,0.079_EB,0.068_EB,0.056_EB,0.045_EB,0.052_EB,0.059_EB,0.066_EB,0.073_EB/)
ABSF(4,7,21:30,5) = (/0.077_EB,0.081_EB,0.086_EB,0.090_EB,0.095_EB,0.098_EB,0.102_EB,0.106_EB,0.110_EB,0.114_EB/)
ABSF(4,7,0:10,6) = (/0.000_EB,0.024_EB,0.055_EB,0.057_EB,0.059_EB,0.068_EB,0.076_EB,0.085_EB,0.094_EB,0.091_EB,0.087_EB/)
ABSF(4,7,11:20,6) = (/0.084_EB,0.080_EB,0.079_EB,0.078_EB,0.076_EB,0.075_EB,0.072_EB,0.069_EB,0.066_EB,0.064_EB/)
ABSF(4,7,21:30,6) = (/0.072_EB,0.081_EB,0.089_EB,0.098_EB,0.106_EB,0.110_EB,0.114_EB,0.118_EB,0.122_EB,0.126_EB/)
ABSF(4,8,0:10,1) = (/0.000_EB,0.032_EB,0.035_EB,0.040_EB,0.046_EB,0.044_EB,0.043_EB,0.042_EB,0.041_EB,0.038_EB,0.036_EB/)
ABSF(4,8,11:20,1) = (/0.034_EB,0.031_EB,0.034_EB,0.037_EB,0.039_EB,0.042_EB,0.044_EB,0.046_EB,0.048_EB,0.049_EB/)
ABSF(4,8,21:30,1) = (/0.049_EB,0.050_EB,0.050_EB,0.050_EB,0.050_EB,0.050_EB,0.049_EB,0.049_EB,0.048_EB,0.048_EB/)
ABSF(4,8,0:10,2) = (/0.000_EB,0.026_EB,0.034_EB,0.039_EB,0.044_EB,0.043_EB,0.043_EB,0.042_EB,0.041_EB,0.039_EB,0.037_EB/)
ABSF(4,8,11:20,2) = (/0.034_EB,0.032_EB,0.034_EB,0.037_EB,0.039_EB,0.041_EB,0.042_EB,0.044_EB,0.045_EB,0.046_EB/)
ABSF(4,8,21:30,2) = (/0.047_EB,0.047_EB,0.048_EB,0.048_EB,0.048_EB,0.049_EB,0.049_EB,0.049_EB,0.049_EB,0.049_EB/)
ABSF(4,8,0:10,3) = (/0.000_EB,0.022_EB,0.036_EB,0.043_EB,0.049_EB,0.051_EB,0.053_EB,0.055_EB,0.057_EB,0.055_EB,0.054_EB/)
ABSF(4,8,11:20,3) = (/0.053_EB,0.051_EB,0.055_EB,0.059_EB,0.063_EB,0.067_EB,0.064_EB,0.062_EB,0.060_EB,0.057_EB/)
ABSF(4,8,21:30,3) = (/0.053_EB,0.050_EB,0.046_EB,0.042_EB,0.038_EB,0.047_EB,0.056_EB,0.065_EB,0.074_EB,0.083_EB/)
ABSF(4,8,0:10,4) = (/0.000_EB,0.024_EB,0.041_EB,0.051_EB,0.061_EB,0.066_EB,0.070_EB,0.075_EB,0.080_EB,0.078_EB,0.076_EB/)
ABSF(4,8,11:20,4) = (/0.074_EB,0.072_EB,0.074_EB,0.077_EB,0.079_EB,0.081_EB,0.081_EB,0.082_EB,0.082_EB,0.083_EB/)
ABSF(4,8,21:30,4) = (/0.073_EB,0.064_EB,0.054_EB,0.044_EB,0.035_EB,0.049_EB,0.062_EB,0.076_EB,0.090_EB,0.103_EB/)
ABSF(4,8,0:10,5) = (/0.000_EB,0.024_EB,0.050_EB,0.054_EB,0.058_EB,0.066_EB,0.073_EB,0.081_EB,0.088_EB,0.089_EB,0.089_EB/)
ABSF(4,8,11:20,5) = (/0.090_EB,0.091_EB,0.080_EB,0.068_EB,0.057_EB,0.046_EB,0.051_EB,0.056_EB,0.061_EB,0.066_EB/)
ABSF(4,8,21:30,5) = (/0.071_EB,0.076_EB,0.081_EB,0.087_EB,0.092_EB,0.096_EB,0.100_EB,0.104_EB,0.108_EB,0.112_EB/)
ABSF(4,8,0:10,6) = (/0.000_EB,0.025_EB,0.058_EB,0.059_EB,0.061_EB,0.070_EB,0.079_EB,0.087_EB,0.096_EB,0.092_EB,0.087_EB/)
ABSF(4,8,11:20,6) = (/0.083_EB,0.079_EB,0.079_EB,0.079_EB,0.079_EB,0.080_EB,0.079_EB,0.078_EB,0.077_EB,0.076_EB/)
ABSF(4,8,21:30,6) = (/0.081_EB,0.087_EB,0.093_EB,0.099_EB,0.105_EB,0.109_EB,0.113_EB,0.117_EB,0.121_EB,0.125_EB/)
ABSF(4,9,0:10,1) = (/0.000_EB,0.034_EB,0.037_EB,0.043_EB,0.049_EB,0.047_EB,0.046_EB,0.044_EB,0.042_EB,0.040_EB,0.037_EB/)
ABSF(4,9,11:20,1) = (/0.035_EB,0.032_EB,0.035_EB,0.038_EB,0.041_EB,0.043_EB,0.045_EB,0.047_EB,0.049_EB,0.051_EB/)
ABSF(4,9,21:30,1) = (/0.051_EB,0.052_EB,0.052_EB,0.052_EB,0.052_EB,0.052_EB,0.051_EB,0.051_EB,0.050_EB,0.050_EB/)
ABSF(4,9,0:10,2) = (/0.000_EB,0.028_EB,0.037_EB,0.042_EB,0.048_EB,0.047_EB,0.045_EB,0.044_EB,0.043_EB,0.041_EB,0.038_EB/)
ABSF(4,9,11:20,2) = (/0.036_EB,0.033_EB,0.035_EB,0.038_EB,0.040_EB,0.042_EB,0.044_EB,0.045_EB,0.046_EB,0.048_EB/)
ABSF(4,9,21:30,2) = (/0.048_EB,0.048_EB,0.049_EB,0.049_EB,0.050_EB,0.050_EB,0.050_EB,0.050_EB,0.050_EB,0.050_EB/)
ABSF(4,9,0:10,3) = (/0.000_EB,0.023_EB,0.038_EB,0.044_EB,0.050_EB,0.051_EB,0.053_EB,0.054_EB,0.055_EB,0.054_EB,0.052_EB/)
ABSF(4,9,11:20,3) = (/0.051_EB,0.050_EB,0.053_EB,0.057_EB,0.061_EB,0.064_EB,0.063_EB,0.062_EB,0.061_EB,0.059_EB/)
ABSF(4,9,21:30,3) = (/0.056_EB,0.053_EB,0.049_EB,0.046_EB,0.043_EB,0.051_EB,0.059_EB,0.067_EB,0.076_EB,0.084_EB/)
ABSF(4,9,0:10,4) = (/0.000_EB,0.025_EB,0.043_EB,0.053_EB,0.063_EB,0.068_EB,0.072_EB,0.077_EB,0.082_EB,0.080_EB,0.078_EB/)
ABSF(4,9,11:20,4) = (/0.076_EB,0.074_EB,0.076_EB,0.078_EB,0.080_EB,0.082_EB,0.082_EB,0.082_EB,0.082_EB,0.082_EB/)
ABSF(4,9,21:30,4) = (/0.071_EB,0.061_EB,0.050_EB,0.039_EB,0.028_EB,0.043_EB,0.058_EB,0.073_EB,0.088_EB,0.103_EB/)
ABSF(4,9,0:10,5) = (/0.000_EB,0.025_EB,0.053_EB,0.056_EB,0.060_EB,0.067_EB,0.075_EB,0.082_EB,0.089_EB,0.090_EB,0.090_EB/)
ABSF(4,9,11:20,5) = (/0.090_EB,0.091_EB,0.080_EB,0.069_EB,0.058_EB,0.048_EB,0.051_EB,0.054_EB,0.057_EB,0.060_EB/)
ABSF(4,9,21:30,5) = (/0.066_EB,0.071_EB,0.077_EB,0.083_EB,0.089_EB,0.093_EB,0.097_EB,0.101_EB,0.106_EB,0.110_EB/)
ABSF(4,9,0:10,6) = (/0.000_EB,0.026_EB,0.060_EB,0.062_EB,0.064_EB,0.072_EB,0.081_EB,0.089_EB,0.098_EB,0.092_EB,0.087_EB/)
ABSF(4,9,11:20,6) = (/0.082_EB,0.077_EB,0.079_EB,0.081_EB,0.083_EB,0.085_EB,0.085_EB,0.086_EB,0.087_EB,0.088_EB/)
ABSF(4,9,21:30,6) = (/0.091_EB,0.094_EB,0.097_EB,0.100_EB,0.103_EB,0.107_EB,0.111_EB,0.116_EB,0.120_EB,0.124_EB/)
ABSF(4,10,0:10,1) = (/0.000_EB,0.036_EB,0.039_EB,0.046_EB,0.052_EB,0.050_EB,0.048_EB,0.046_EB,0.044_EB,0.041_EB,0.039_EB/)
ABSF(4,10,11:20,1) = (/0.036_EB,0.034_EB,0.036_EB,0.039_EB,0.042_EB,0.045_EB,0.047_EB,0.049_EB,0.051_EB,0.053_EB/)
ABSF(4,10,21:30,1) = (/0.053_EB,0.053_EB,0.054_EB,0.054_EB,0.054_EB,0.054_EB,0.053_EB,0.052_EB,0.052_EB,0.051_EB/)
ABSF(4,10,0:10,2) = (/0.000_EB,0.031_EB,0.040_EB,0.046_EB,0.051_EB,0.050_EB,0.048_EB,0.046_EB,0.045_EB,0.042_EB,0.039_EB/)
ABSF(4,10,11:20,2) = (/0.037_EB,0.034_EB,0.036_EB,0.039_EB,0.041_EB,0.043_EB,0.045_EB,0.046_EB,0.048_EB,0.049_EB/)
ABSF(4,10,21:30,2) = (/0.049_EB,0.050_EB,0.050_EB,0.051_EB,0.051_EB,0.051_EB,0.051_EB,0.051_EB,0.051_EB,0.051_EB/)
ABSF(4,10,0:10,3) = (/0.000_EB,0.024_EB,0.040_EB,0.045_EB,0.051_EB,0.051_EB,0.052_EB,0.053_EB,0.053_EB,0.052_EB,0.051_EB/)
ABSF(4,10,11:20,3) = (/0.049_EB,0.048_EB,0.052_EB,0.055_EB,0.058_EB,0.062_EB,0.062_EB,0.062_EB,0.062_EB,0.062_EB/)
ABSF(4,10,21:30,3) = (/0.059_EB,0.056_EB,0.053_EB,0.050_EB,0.047_EB,0.055_EB,0.062_EB,0.070_EB,0.077_EB,0.085_EB/)
ABSF(4,10,0:10,4) = (/0.000_EB,0.026_EB,0.044_EB,0.055_EB,0.065_EB,0.070_EB,0.074_EB,0.079_EB,0.083_EB,0.082_EB,0.080_EB/)
ABSF(4,10,11:20,4) = (/0.078_EB,0.076_EB,0.078_EB,0.079_EB,0.081_EB,0.083_EB,0.082_EB,0.082_EB,0.082_EB,0.081_EB/)
ABSF(4,10,21:30,4) = (/0.069_EB,0.057_EB,0.046_EB,0.034_EB,0.022_EB,0.038_EB,0.054_EB,0.070_EB,0.086_EB,0.102_EB/)
ABSF(4,10,0:10,5) = (/0.000_EB,0.026_EB,0.055_EB,0.058_EB,0.062_EB,0.069_EB,0.076_EB,0.083_EB,0.090_EB,0.090_EB,0.091_EB/)
ABSF(4,10,11:20,5) = (/0.091_EB,0.091_EB,0.080_EB,0.070_EB,0.060_EB,0.049_EB,0.050_EB,0.051_EB,0.052_EB,0.053_EB/)
ABSF(4,10,21:30,5) = (/0.060_EB,0.066_EB,0.073_EB,0.079_EB,0.086_EB,0.090_EB,0.095_EB,0.099_EB,0.103_EB,0.108_EB/)
ABSF(4,10,0:10,6) = (/0.000_EB,0.027_EB,0.062_EB,0.064_EB,0.066_EB,0.075_EB,0.083_EB,0.091_EB,0.099_EB,0.093_EB,0.087_EB/)
ABSF(4,10,11:20,6) = (/0.081_EB,0.075_EB,0.079_EB,0.082_EB,0.086_EB,0.090_EB,0.092_EB,0.095_EB,0.097_EB,0.099_EB/)
ABSF(4,10,21:30,6) = (/0.100_EB,0.100_EB,0.101_EB,0.101_EB,0.101_EB,0.106_EB,0.110_EB,0.114_EB,0.119_EB,0.123_EB/)
ABSF(4,11,0:10,1) = (/0.000_EB,0.037_EB,0.042_EB,0.048_EB,0.055_EB,0.052_EB,0.050_EB,0.048_EB,0.046_EB,0.043_EB,0.041_EB/)
ABSF(4,11,11:20,1) = (/0.038_EB,0.035_EB,0.038_EB,0.041_EB,0.043_EB,0.046_EB,0.048_EB,0.050_EB,0.053_EB,0.055_EB/)
ABSF(4,11,21:30,1) = (/0.055_EB,0.055_EB,0.055_EB,0.056_EB,0.056_EB,0.055_EB,0.055_EB,0.055_EB,0.054_EB,0.054_EB/)
ABSF(4,11,0:10,2) = (/0.000_EB,0.032_EB,0.042_EB,0.047_EB,0.053_EB,0.052_EB,0.050_EB,0.048_EB,0.047_EB,0.044_EB,0.041_EB/)
ABSF(4,11,11:20,2) = (/0.038_EB,0.036_EB,0.038_EB,0.040_EB,0.043_EB,0.045_EB,0.047_EB,0.048_EB,0.050_EB,0.051_EB/)
ABSF(4,11,21:30,2) = (/0.052_EB,0.052_EB,0.052_EB,0.053_EB,0.053_EB,0.054_EB,0.054_EB,0.054_EB,0.054_EB,0.054_EB/)
ABSF(4,11,0:10,3) = (/0.000_EB,0.025_EB,0.041_EB,0.048_EB,0.054_EB,0.054_EB,0.054_EB,0.054_EB,0.054_EB,0.053_EB,0.052_EB/)
ABSF(4,11,11:20,3) = (/0.050_EB,0.049_EB,0.052_EB,0.055_EB,0.058_EB,0.062_EB,0.062_EB,0.062_EB,0.062_EB,0.062_EB/)
ABSF(4,11,21:30,3) = (/0.060_EB,0.057_EB,0.055_EB,0.053_EB,0.050_EB,0.057_EB,0.064_EB,0.070_EB,0.077_EB,0.084_EB/)
ABSF(4,11,0:10,4) = (/0.000_EB,0.027_EB,0.046_EB,0.056_EB,0.067_EB,0.071_EB,0.076_EB,0.080_EB,0.085_EB,0.083_EB,0.081_EB/)
ABSF(4,11,11:20,4) = (/0.079_EB,0.077_EB,0.079_EB,0.080_EB,0.082_EB,0.084_EB,0.083_EB,0.083_EB,0.082_EB,0.081_EB/)
ABSF(4,11,21:30,4) = (/0.070_EB,0.059_EB,0.048_EB,0.037_EB,0.026_EB,0.040_EB,0.054_EB,0.069_EB,0.083_EB,0.098_EB/)
ABSF(4,11,0:10,5) = (/0.000_EB,0.027_EB,0.056_EB,0.060_EB,0.064_EB,0.071_EB,0.078_EB,0.085_EB,0.092_EB,0.092_EB,0.091_EB/)
ABSF(4,11,11:20,5) = (/0.091_EB,0.091_EB,0.080_EB,0.068_EB,0.057_EB,0.046_EB,0.048_EB,0.050_EB,0.052_EB,0.053_EB/)
ABSF(4,11,21:30,5) = (/0.059_EB,0.065_EB,0.070_EB,0.076_EB,0.081_EB,0.087_EB,0.093_EB,0.099_EB,0.105_EB,0.110_EB/)
ABSF(4,11,0:10,6) = (/0.000_EB,0.028_EB,0.063_EB,0.065_EB,0.068_EB,0.076_EB,0.084_EB,0.092_EB,0.099_EB,0.094_EB,0.089_EB/)
ABSF(4,11,11:20,6) = (/0.084_EB,0.078_EB,0.080_EB,0.082_EB,0.085_EB,0.087_EB,0.088_EB,0.090_EB,0.091_EB,0.092_EB/)
ABSF(4,11,21:30,6) = (/0.094_EB,0.095_EB,0.097_EB,0.098_EB,0.099_EB,0.102_EB,0.105_EB,0.108_EB,0.111_EB,0.114_EB/)
ABSF(4,12,0:10,1) = (/0.000_EB,0.039_EB,0.044_EB,0.050_EB,0.057_EB,0.055_EB,0.052_EB,0.050_EB,0.048_EB,0.045_EB,0.042_EB/)
ABSF(4,12,11:20,1) = (/0.040_EB,0.037_EB,0.040_EB,0.042_EB,0.045_EB,0.048_EB,0.050_EB,0.052_EB,0.054_EB,0.056_EB/)
ABSF(4,12,21:30,1) = (/0.056_EB,0.057_EB,0.057_EB,0.057_EB,0.058_EB,0.057_EB,0.057_EB,0.057_EB,0.056_EB,0.056_EB/)
ABSF(4,12,0:10,2) = (/0.000_EB,0.034_EB,0.044_EB,0.049_EB,0.055_EB,0.053_EB,0.052_EB,0.050_EB,0.048_EB,0.046_EB,0.043_EB/)
ABSF(4,12,11:20,2) = (/0.040_EB,0.037_EB,0.040_EB,0.042_EB,0.045_EB,0.047_EB,0.049_EB,0.050_EB,0.052_EB,0.053_EB/)
ABSF(4,12,21:30,2) = (/0.054_EB,0.054_EB,0.055_EB,0.055_EB,0.056_EB,0.056_EB,0.056_EB,0.056_EB,0.057_EB,0.057_EB/)
ABSF(4,12,0:10,3) = (/0.000_EB,0.026_EB,0.043_EB,0.050_EB,0.056_EB,0.056_EB,0.056_EB,0.056_EB,0.055_EB,0.054_EB,0.052_EB/)
ABSF(4,12,11:20,3) = (/0.051_EB,0.049_EB,0.052_EB,0.055_EB,0.058_EB,0.062_EB,0.062_EB,0.062_EB,0.062_EB,0.063_EB/)
ABSF(4,12,21:30,3) = (/0.061_EB,0.059_EB,0.057_EB,0.055_EB,0.053_EB,0.059_EB,0.065_EB,0.071_EB,0.077_EB,0.083_EB/)
ABSF(4,12,0:10,4) = (/0.000_EB,0.028_EB,0.047_EB,0.058_EB,0.068_EB,0.073_EB,0.077_EB,0.082_EB,0.086_EB,0.084_EB,0.082_EB/)
ABSF(4,12,11:20,4) = (/0.079_EB,0.077_EB,0.079_EB,0.081_EB,0.083_EB,0.085_EB,0.084_EB,0.083_EB,0.082_EB,0.081_EB/)
ABSF(4,12,21:30,4) = (/0.071_EB,0.061_EB,0.050_EB,0.040_EB,0.029_EB,0.042_EB,0.055_EB,0.068_EB,0.081_EB,0.093_EB/)
ABSF(4,12,0:10,5) = (/0.000_EB,0.028_EB,0.058_EB,0.062_EB,0.066_EB,0.073_EB,0.080_EB,0.087_EB,0.094_EB,0.093_EB,0.092_EB/)
ABSF(4,12,11:20,5) = (/0.092_EB,0.091_EB,0.079_EB,0.067_EB,0.055_EB,0.043_EB,0.045_EB,0.048_EB,0.051_EB,0.053_EB/)
ABSF(4,12,21:30,5) = (/0.058_EB,0.063_EB,0.067_EB,0.072_EB,0.076_EB,0.084_EB,0.091_EB,0.098_EB,0.106_EB,0.113_EB/)
ABSF(4,12,0:10,6) = (/0.000_EB,0.029_EB,0.063_EB,0.067_EB,0.070_EB,0.077_EB,0.085_EB,0.092_EB,0.099_EB,0.095_EB,0.091_EB/)
ABSF(4,12,11:20,6) = (/0.086_EB,0.082_EB,0.082_EB,0.082_EB,0.083_EB,0.083_EB,0.084_EB,0.084_EB,0.085_EB,0.085_EB/)
ABSF(4,12,21:30,6) = (/0.088_EB,0.090_EB,0.093_EB,0.095_EB,0.097_EB,0.099_EB,0.100_EB,0.102_EB,0.103_EB,0.105_EB/)
ABSF(4,13,0:10,1) = (/0.000_EB,0.040_EB,0.046_EB,0.053_EB,0.060_EB,0.057_EB,0.055_EB,0.052_EB,0.050_EB,0.047_EB,0.044_EB/)
ABSF(4,13,11:20,1) = (/0.041_EB,0.039_EB,0.041_EB,0.044_EB,0.047_EB,0.050_EB,0.052_EB,0.054_EB,0.056_EB,0.058_EB/)
ABSF(4,13,21:30,1) = (/0.058_EB,0.058_EB,0.059_EB,0.059_EB,0.059_EB,0.059_EB,0.059_EB,0.059_EB,0.059_EB,0.058_EB/)
ABSF(4,13,0:10,2) = (/0.000_EB,0.036_EB,0.046_EB,0.051_EB,0.057_EB,0.055_EB,0.053_EB,0.052_EB,0.050_EB,0.047_EB,0.045_EB/)
ABSF(4,13,11:20,2) = (/0.042_EB,0.039_EB,0.041_EB,0.044_EB,0.047_EB,0.049_EB,0.051_EB,0.052_EB,0.054_EB,0.055_EB/)
ABSF(4,13,21:30,2) = (/0.056_EB,0.056_EB,0.057_EB,0.057_EB,0.058_EB,0.058_EB,0.059_EB,0.059_EB,0.059_EB,0.060_EB/)
ABSF(4,13,0:10,3) = (/0.000_EB,0.027_EB,0.044_EB,0.052_EB,0.059_EB,0.058_EB,0.058_EB,0.057_EB,0.057_EB,0.055_EB,0.053_EB/)
ABSF(4,13,11:20,3) = (/0.052_EB,0.050_EB,0.053_EB,0.056_EB,0.059_EB,0.061_EB,0.062_EB,0.062_EB,0.063_EB,0.063_EB/)
ABSF(4,13,21:30,3) = (/0.062_EB,0.060_EB,0.059_EB,0.057_EB,0.056_EB,0.061_EB,0.067_EB,0.072_EB,0.077_EB,0.082_EB/)
ABSF(4,13,0:10,4) = (/0.000_EB,0.029_EB,0.049_EB,0.059_EB,0.070_EB,0.074_EB,0.079_EB,0.083_EB,0.087_EB,0.085_EB,0.082_EB/)
ABSF(4,13,11:20,4) = (/0.080_EB,0.078_EB,0.080_EB,0.082_EB,0.084_EB,0.086_EB,0.085_EB,0.084_EB,0.083_EB,0.081_EB/)
ABSF(4,13,21:30,4) = (/0.072_EB,0.062_EB,0.053_EB,0.043_EB,0.033_EB,0.045_EB,0.056_EB,0.067_EB,0.078_EB,0.089_EB/)
ABSF(4,13,0:10,5) = (/0.000_EB,0.029_EB,0.059_EB,0.064_EB,0.068_EB,0.075_EB,0.082_EB,0.088_EB,0.095_EB,0.094_EB,0.093_EB/)
ABSF(4,13,11:20,5) = (/0.092_EB,0.091_EB,0.078_EB,0.065_EB,0.052_EB,0.039_EB,0.043_EB,0.046_EB,0.050_EB,0.054_EB/)
ABSF(4,13,21:30,5) = (/0.057_EB,0.061_EB,0.064_EB,0.068_EB,0.072_EB,0.080_EB,0.089_EB,0.098_EB,0.107_EB,0.116_EB/)
ABSF(4,13,0:10,6) = (/0.000_EB,0.030_EB,0.064_EB,0.068_EB,0.072_EB,0.079_EB,0.086_EB,0.093_EB,0.099_EB,0.096_EB,0.092_EB/)
ABSF(4,13,11:20,6) = (/0.088_EB,0.085_EB,0.084_EB,0.083_EB,0.081_EB,0.080_EB,0.080_EB,0.079_EB,0.079_EB,0.078_EB/)
ABSF(4,13,21:30,6) = (/0.082_EB,0.085_EB,0.089_EB,0.092_EB,0.096_EB,0.096_EB,0.096_EB,0.096_EB,0.096_EB,0.096_EB/)
ABSF(4,14,0:10,1) = (/0.000_EB,0.042_EB,0.048_EB,0.055_EB,0.062_EB,0.060_EB,0.057_EB,0.054_EB,0.052_EB,0.049_EB,0.046_EB/)
ABSF(4,14,11:20,1) = (/0.043_EB,0.040_EB,0.043_EB,0.046_EB,0.049_EB,0.051_EB,0.053_EB,0.055_EB,0.057_EB,0.059_EB/)
ABSF(4,14,21:30,1) = (/0.060_EB,0.060_EB,0.060_EB,0.061_EB,0.061_EB,0.061_EB,0.061_EB,0.061_EB,0.061_EB,0.061_EB/)
ABSF(4,14,0:10,2) = (/0.000_EB,0.038_EB,0.048_EB,0.053_EB,0.058_EB,0.057_EB,0.055_EB,0.054_EB,0.052_EB,0.049_EB,0.046_EB/)
ABSF(4,14,11:20,2) = (/0.043_EB,0.040_EB,0.043_EB,0.046_EB,0.048_EB,0.051_EB,0.053_EB,0.054_EB,0.056_EB,0.058_EB/)
ABSF(4,14,21:30,2) = (/0.058_EB,0.058_EB,0.059_EB,0.059_EB,0.060_EB,0.060_EB,0.061_EB,0.062_EB,0.062_EB,0.063_EB/)
ABSF(4,14,0:10,3) = (/0.000_EB,0.028_EB,0.046_EB,0.054_EB,0.062_EB,0.061_EB,0.060_EB,0.059_EB,0.058_EB,0.056_EB,0.054_EB/)
ABSF(4,14,11:20,3) = (/0.052_EB,0.051_EB,0.053_EB,0.056_EB,0.059_EB,0.061_EB,0.062_EB,0.063_EB,0.063_EB,0.064_EB/)
ABSF(4,14,21:30,3) = (/0.063_EB,0.062_EB,0.061_EB,0.060_EB,0.059_EB,0.063_EB,0.068_EB,0.073_EB,0.077_EB,0.082_EB/)
ABSF(4,14,0:10,4) = (/0.000_EB,0.030_EB,0.050_EB,0.061_EB,0.071_EB,0.076_EB,0.080_EB,0.084_EB,0.088_EB,0.086_EB,0.083_EB/)
ABSF(4,14,11:20,4) = (/0.081_EB,0.078_EB,0.081_EB,0.083_EB,0.085_EB,0.087_EB,0.086_EB,0.084_EB,0.083_EB,0.081_EB/)
ABSF(4,14,21:30,4) = (/0.073_EB,0.064_EB,0.055_EB,0.046_EB,0.037_EB,0.047_EB,0.056_EB,0.066_EB,0.075_EB,0.085_EB/)
ABSF(4,14,0:10,5) = (/0.000_EB,0.030_EB,0.061_EB,0.065_EB,0.070_EB,0.077_EB,0.083_EB,0.090_EB,0.097_EB,0.096_EB,0.094_EB/)
ABSF(4,14,11:20,5) = (/0.093_EB,0.091_EB,0.078_EB,0.064_EB,0.050_EB,0.036_EB,0.040_EB,0.045_EB,0.049_EB,0.054_EB/)
ABSF(4,14,21:30,5) = (/0.056_EB,0.059_EB,0.061_EB,0.064_EB,0.067_EB,0.077_EB,0.088_EB,0.098_EB,0.108_EB,0.119_EB/)
ABSF(4,14,0:10,6) = (/0.000_EB,0.031_EB,0.065_EB,0.069_EB,0.074_EB,0.080_EB,0.087_EB,0.093_EB,0.099_EB,0.097_EB,0.094_EB/)
ABSF(4,14,11:20,6) = (/0.091_EB,0.088_EB,0.085_EB,0.083_EB,0.080_EB,0.077_EB,0.076_EB,0.074_EB,0.073_EB,0.071_EB/)
ABSF(4,14,21:30,6) = (/0.076_EB,0.080_EB,0.085_EB,0.089_EB,0.094_EB,0.092_EB,0.091_EB,0.090_EB,0.088_EB,0.087_EB/)
ABSF(4,15,0:10,1) = (/0.000_EB,0.043_EB,0.050_EB,0.058_EB,0.065_EB,0.062_EB,0.059_EB,0.056_EB,0.053_EB,0.051_EB,0.048_EB/)
ABSF(4,15,11:20,1) = (/0.045_EB,0.042_EB,0.045_EB,0.048_EB,0.050_EB,0.053_EB,0.055_EB,0.057_EB,0.059_EB,0.061_EB/)
ABSF(4,15,21:30,1) = (/0.061_EB,0.062_EB,0.062_EB,0.063_EB,0.063_EB,0.063_EB,0.063_EB,0.063_EB,0.063_EB,0.063_EB/)
ABSF(4,15,0:10,2) = (/0.000_EB,0.040_EB,0.050_EB,0.055_EB,0.060_EB,0.058_EB,0.057_EB,0.055_EB,0.054_EB,0.051_EB,0.048_EB/)
ABSF(4,15,11:20,2) = (/0.045_EB,0.042_EB,0.045_EB,0.047_EB,0.050_EB,0.053_EB,0.055_EB,0.056_EB,0.058_EB,0.060_EB/)
ABSF(4,15,21:30,2) = (/0.060_EB,0.061_EB,0.061_EB,0.062_EB,0.062_EB,0.063_EB,0.063_EB,0.064_EB,0.065_EB,0.065_EB/)
ABSF(4,15,0:10,3) = (/0.000_EB,0.030_EB,0.047_EB,0.056_EB,0.065_EB,0.063_EB,0.062_EB,0.060_EB,0.059_EB,0.057_EB,0.055_EB/)
ABSF(4,15,11:20,3) = (/0.053_EB,0.051_EB,0.054_EB,0.056_EB,0.059_EB,0.061_EB,0.062_EB,0.063_EB,0.063_EB,0.064_EB/)
ABSF(4,15,21:30,3) = (/0.064_EB,0.063_EB,0.063_EB,0.062_EB,0.062_EB,0.066_EB,0.069_EB,0.073_EB,0.077_EB,0.081_EB/)
ABSF(4,15,0:10,4) = (/0.000_EB,0.031_EB,0.052_EB,0.062_EB,0.073_EB,0.077_EB,0.081_EB,0.085_EB,0.090_EB,0.087_EB,0.084_EB/)
ABSF(4,15,11:20,4) = (/0.082_EB,0.079_EB,0.081_EB,0.084_EB,0.086_EB,0.089_EB,0.087_EB,0.085_EB,0.083_EB,0.081_EB/)
ABSF(4,15,21:30,4) = (/0.073_EB,0.065_EB,0.057_EB,0.049_EB,0.041_EB,0.049_EB,0.057_EB,0.065_EB,0.073_EB,0.081_EB/)
ABSF(4,15,0:10,5) = (/0.000_EB,0.031_EB,0.062_EB,0.067_EB,0.072_EB,0.079_EB,0.085_EB,0.092_EB,0.099_EB,0.097_EB,0.095_EB/)
ABSF(4,15,11:20,5) = (/0.093_EB,0.092_EB,0.077_EB,0.062_EB,0.047_EB,0.033_EB,0.038_EB,0.043_EB,0.048_EB,0.054_EB/)
ABSF(4,15,21:30,5) = (/0.055_EB,0.057_EB,0.059_EB,0.060_EB,0.062_EB,0.074_EB,0.086_EB,0.098_EB,0.110_EB,0.122_EB/)
ABSF(4,15,0:10,6) = (/0.000_EB,0.032_EB,0.066_EB,0.071_EB,0.075_EB,0.081_EB,0.087_EB,0.093_EB,0.099_EB,0.097_EB,0.095_EB/)
ABSF(4,15,11:20,6) = (/0.093_EB,0.091_EB,0.087_EB,0.083_EB,0.078_EB,0.074_EB,0.072_EB,0.069_EB,0.067_EB,0.064_EB/)
ABSF(4,15,21:30,6) = (/0.070_EB,0.075_EB,0.081_EB,0.086_EB,0.092_EB,0.089_EB,0.086_EB,0.084_EB,0.081_EB,0.078_EB/)
ABSF(4,16,0:10,1) = (/0.000_EB,0.045_EB,0.052_EB,0.060_EB,0.068_EB,0.065_EB,0.061_EB,0.058_EB,0.055_EB,0.052_EB,0.050_EB/)
ABSF(4,16,11:20,1) = (/0.047_EB,0.044_EB,0.046_EB,0.049_EB,0.052_EB,0.055_EB,0.057_EB,0.059_EB,0.061_EB,0.062_EB/)
ABSF(4,16,21:30,1) = (/0.063_EB,0.063_EB,0.064_EB,0.064_EB,0.065_EB,0.065_EB,0.065_EB,0.065_EB,0.065_EB,0.066_EB/)
ABSF(4,16,0:10,2) = (/0.000_EB,0.042_EB,0.053_EB,0.057_EB,0.062_EB,0.060_EB,0.059_EB,0.057_EB,0.056_EB,0.053_EB,0.050_EB/)
ABSF(4,16,11:20,2) = (/0.047_EB,0.044_EB,0.046_EB,0.049_EB,0.052_EB,0.055_EB,0.057_EB,0.058_EB,0.060_EB,0.062_EB/)
ABSF(4,16,21:30,2) = (/0.062_EB,0.063_EB,0.063_EB,0.064_EB,0.064_EB,0.065_EB,0.066_EB,0.067_EB,0.067_EB,0.068_EB/)
ABSF(4,16,0:10,3) = (/0.000_EB,0.031_EB,0.049_EB,0.058_EB,0.067_EB,0.065_EB,0.064_EB,0.062_EB,0.060_EB,0.058_EB,0.056_EB/)
ABSF(4,16,11:20,3) = (/0.054_EB,0.052_EB,0.054_EB,0.057_EB,0.059_EB,0.061_EB,0.062_EB,0.063_EB,0.064_EB,0.065_EB/)
ABSF(4,16,21:30,3) = (/0.065_EB,0.065_EB,0.065_EB,0.065_EB,0.065_EB,0.068_EB,0.071_EB,0.074_EB,0.077_EB,0.080_EB/)
ABSF(4,16,0:10,4) = (/0.000_EB,0.032_EB,0.053_EB,0.064_EB,0.074_EB,0.079_EB,0.083_EB,0.087_EB,0.091_EB,0.088_EB,0.085_EB/)
ABSF(4,16,11:20,4) = (/0.082_EB,0.080_EB,0.082_EB,0.085_EB,0.087_EB,0.090_EB,0.088_EB,0.086_EB,0.084_EB,0.082_EB/)
ABSF(4,16,21:30,4) = (/0.074_EB,0.067_EB,0.060_EB,0.052_EB,0.045_EB,0.051_EB,0.058_EB,0.064_EB,0.070_EB,0.076_EB/)
ABSF(4,16,0:10,5) = (/0.000_EB,0.032_EB,0.064_EB,0.069_EB,0.074_EB,0.081_EB,0.087_EB,0.094_EB,0.100_EB,0.098_EB,0.096_EB/)
ABSF(4,16,11:20,5) = (/0.094_EB,0.092_EB,0.076_EB,0.061_EB,0.045_EB,0.029_EB,0.035_EB,0.042_EB,0.048_EB,0.054_EB/)
ABSF(4,16,21:30,5) = (/0.054_EB,0.055_EB,0.056_EB,0.057_EB,0.057_EB,0.071_EB,0.084_EB,0.098_EB,0.111_EB,0.124_EB/)
ABSF(4,16,0:10,6) = (/0.000_EB,0.033_EB,0.066_EB,0.072_EB,0.077_EB,0.083_EB,0.088_EB,0.094_EB,0.099_EB,0.098_EB,0.097_EB/)
ABSF(4,16,11:20,6) = (/0.096_EB,0.095_EB,0.089_EB,0.083_EB,0.077_EB,0.071_EB,0.067_EB,0.064_EB,0.061_EB,0.057_EB/)
ABSF(4,16,21:30,6) = (/0.064_EB,0.070_EB,0.077_EB,0.083_EB,0.090_EB,0.086_EB,0.081_EB,0.077_EB,0.073_EB,0.069_EB/)
ABSF(4,17,0:10,1) = (/0.000_EB,0.046_EB,0.055_EB,0.062_EB,0.070_EB,0.067_EB,0.064_EB,0.060_EB,0.057_EB,0.054_EB,0.051_EB/)
ABSF(4,17,11:20,1) = (/0.048_EB,0.045_EB,0.048_EB,0.051_EB,0.054_EB,0.056_EB,0.058_EB,0.060_EB,0.062_EB,0.064_EB/)
ABSF(4,17,21:30,1) = (/0.065_EB,0.065_EB,0.066_EB,0.066_EB,0.067_EB,0.067_EB,0.067_EB,0.067_EB,0.068_EB,0.068_EB/)
ABSF(4,17,0:10,2) = (/0.000_EB,0.044_EB,0.055_EB,0.059_EB,0.063_EB,0.062_EB,0.060_EB,0.059_EB,0.058_EB,0.054_EB,0.051_EB/)
ABSF(4,17,11:20,2) = (/0.048_EB,0.045_EB,0.048_EB,0.051_EB,0.054_EB,0.057_EB,0.059_EB,0.060_EB,0.062_EB,0.064_EB/)
ABSF(4,17,21:30,2) = (/0.064_EB,0.065_EB,0.065_EB,0.066_EB,0.066_EB,0.067_EB,0.068_EB,0.069_EB,0.070_EB,0.071_EB/)
ABSF(4,17,0:10,3) = (/0.000_EB,0.032_EB,0.050_EB,0.060_EB,0.070_EB,0.068_EB,0.066_EB,0.063_EB,0.061_EB,0.059_EB,0.057_EB/)
ABSF(4,17,11:20,3) = (/0.054_EB,0.052_EB,0.055_EB,0.057_EB,0.059_EB,0.061_EB,0.062_EB,0.063_EB,0.064_EB,0.065_EB/)
ABSF(4,17,21:30,3) = (/0.065_EB,0.066_EB,0.066_EB,0.067_EB,0.068_EB,0.070_EB,0.072_EB,0.075_EB,0.077_EB,0.079_EB/)
ABSF(4,17,0:10,4) = (/0.000_EB,0.033_EB,0.055_EB,0.065_EB,0.076_EB,0.080_EB,0.084_EB,0.088_EB,0.092_EB,0.089_EB,0.086_EB/)
ABSF(4,17,11:20,4) = (/0.083_EB,0.080_EB,0.083_EB,0.086_EB,0.088_EB,0.091_EB,0.089_EB,0.086_EB,0.084_EB,0.082_EB/)
ABSF(4,17,21:30,4) = (/0.075_EB,0.069_EB,0.062_EB,0.055_EB,0.049_EB,0.054_EB,0.058_EB,0.063_EB,0.067_EB,0.072_EB/)
ABSF(4,17,0:10,5) = (/0.000_EB,0.033_EB,0.065_EB,0.071_EB,0.076_EB,0.083_EB,0.089_EB,0.096_EB,0.102_EB,0.099_EB,0.097_EB/)
ABSF(4,17,11:20,5) = (/0.094_EB,0.092_EB,0.075_EB,0.059_EB,0.043_EB,0.026_EB,0.033_EB,0.040_EB,0.047_EB,0.054_EB/)
ABSF(4,17,21:30,5) = (/0.053_EB,0.053_EB,0.053_EB,0.053_EB,0.052_EB,0.067_EB,0.082_EB,0.097_EB,0.112_EB,0.127_EB/)
ABSF(4,17,0:10,6) = (/0.000_EB,0.034_EB,0.067_EB,0.073_EB,0.079_EB,0.084_EB,0.089_EB,0.094_EB,0.099_EB,0.099_EB,0.099_EB/)
ABSF(4,17,11:20,6) = (/0.098_EB,0.098_EB,0.090_EB,0.083_EB,0.075_EB,0.068_EB,0.063_EB,0.059_EB,0.055_EB,0.050_EB/)
ABSF(4,17,21:30,6) = (/0.058_EB,0.065_EB,0.073_EB,0.080_EB,0.088_EB,0.082_EB,0.077_EB,0.071_EB,0.066_EB,0.060_EB/)
ABSF(4,18,0:10,1) = (/0.000_EB,0.048_EB,0.057_EB,0.065_EB,0.073_EB,0.069_EB,0.066_EB,0.063_EB,0.059_EB,0.056_EB,0.053_EB/)
ABSF(4,18,11:20,1) = (/0.050_EB,0.047_EB,0.050_EB,0.053_EB,0.055_EB,0.058_EB,0.060_EB,0.062_EB,0.064_EB,0.066_EB/)
ABSF(4,18,21:30,1) = (/0.066_EB,0.067_EB,0.067_EB,0.068_EB,0.068_EB,0.069_EB,0.069_EB,0.070_EB,0.070_EB,0.070_EB/)
ABSF(4,18,0:10,2) = (/0.000_EB,0.046_EB,0.057_EB,0.061_EB,0.065_EB,0.064_EB,0.062_EB,0.061_EB,0.059_EB,0.056_EB,0.053_EB/)
ABSF(4,18,11:20,2) = (/0.050_EB,0.047_EB,0.050_EB,0.053_EB,0.056_EB,0.059_EB,0.061_EB,0.062_EB,0.064_EB,0.066_EB/)
ABSF(4,18,21:30,2) = (/0.067_EB,0.067_EB,0.068_EB,0.068_EB,0.069_EB,0.070_EB,0.071_EB,0.072_EB,0.073_EB,0.074_EB/)
ABSF(4,18,0:10,3) = (/0.000_EB,0.033_EB,0.052_EB,0.062_EB,0.073_EB,0.070_EB,0.067_EB,0.065_EB,0.062_EB,0.060_EB,0.057_EB/)
ABSF(4,18,11:20,3) = (/0.055_EB,0.053_EB,0.055_EB,0.057_EB,0.059_EB,0.061_EB,0.062_EB,0.063_EB,0.064_EB,0.065_EB/)
ABSF(4,18,21:30,3) = (/0.066_EB,0.067_EB,0.068_EB,0.069_EB,0.070_EB,0.072_EB,0.074_EB,0.075_EB,0.077_EB,0.079_EB/)
ABSF(4,18,0:10,4) = (/0.000_EB,0.033_EB,0.056_EB,0.067_EB,0.077_EB,0.081_EB,0.085_EB,0.089_EB,0.093_EB,0.090_EB,0.087_EB/)
ABSF(4,18,11:20,4) = (/0.084_EB,0.081_EB,0.084_EB,0.086_EB,0.089_EB,0.092_EB,0.089_EB,0.087_EB,0.084_EB,0.082_EB/)
ABSF(4,18,21:30,4) = (/0.076_EB,0.070_EB,0.064_EB,0.059_EB,0.053_EB,0.056_EB,0.059_EB,0.062_EB,0.065_EB,0.068_EB/)
ABSF(4,18,0:10,5) = (/0.000_EB,0.034_EB,0.067_EB,0.073_EB,0.078_EB,0.085_EB,0.091_EB,0.097_EB,0.104_EB,0.101_EB,0.098_EB/)
ABSF(4,18,11:20,5) = (/0.095_EB,0.092_EB,0.075_EB,0.057_EB,0.040_EB,0.023_EB,0.031_EB,0.038_EB,0.046_EB,0.054_EB/)
ABSF(4,18,21:30,5) = (/0.052_EB,0.051_EB,0.050_EB,0.049_EB,0.048_EB,0.064_EB,0.081_EB,0.097_EB,0.114_EB,0.130_EB/)
ABSF(4,18,0:10,6) = (/0.000_EB,0.035_EB,0.068_EB,0.074_EB,0.081_EB,0.086_EB,0.090_EB,0.095_EB,0.099_EB,0.100_EB,0.100_EB/)
ABSF(4,18,11:20,6) = (/0.101_EB,0.101_EB,0.092_EB,0.083_EB,0.074_EB,0.064_EB,0.059_EB,0.054_EB,0.049_EB,0.043_EB/)
ABSF(4,18,21:30,6) = (/0.052_EB,0.060_EB,0.069_EB,0.077_EB,0.086_EB,0.079_EB,0.072_EB,0.065_EB,0.058_EB,0.051_EB/)
ABSF(4,19,0:10,1) = (/0.000_EB,0.049_EB,0.059_EB,0.067_EB,0.076_EB,0.072_EB,0.068_EB,0.065_EB,0.061_EB,0.058_EB,0.055_EB/)
ABSF(4,19,11:20,1) = (/0.052_EB,0.049_EB,0.052_EB,0.054_EB,0.057_EB,0.060_EB,0.062_EB,0.063_EB,0.065_EB,0.067_EB/)
ABSF(4,19,21:30,1) = (/0.068_EB,0.068_EB,0.069_EB,0.069_EB,0.070_EB,0.071_EB,0.071_EB,0.072_EB,0.072_EB,0.073_EB/)
ABSF(4,19,0:10,2) = (/0.000_EB,0.048_EB,0.059_EB,0.063_EB,0.067_EB,0.065_EB,0.064_EB,0.063_EB,0.061_EB,0.058_EB,0.055_EB/)
ABSF(4,19,11:20,2) = (/0.052_EB,0.048_EB,0.051_EB,0.055_EB,0.058_EB,0.061_EB,0.062_EB,0.064_EB,0.066_EB,0.068_EB/)
ABSF(4,19,21:30,2) = (/0.069_EB,0.069_EB,0.070_EB,0.070_EB,0.071_EB,0.072_EB,0.073_EB,0.074_EB,0.076_EB,0.077_EB/)
ABSF(4,19,0:10,3) = (/0.000_EB,0.034_EB,0.053_EB,0.064_EB,0.076_EB,0.072_EB,0.069_EB,0.066_EB,0.063_EB,0.061_EB,0.058_EB/)
ABSF(4,19,11:20,3) = (/0.056_EB,0.053_EB,0.055_EB,0.057_EB,0.059_EB,0.061_EB,0.062_EB,0.064_EB,0.065_EB,0.066_EB/)
ABSF(4,19,21:30,3) = (/0.067_EB,0.069_EB,0.070_EB,0.072_EB,0.073_EB,0.074_EB,0.075_EB,0.076_EB,0.077_EB,0.078_EB/)
ABSF(4,19,0:10,4) = (/0.000_EB,0.034_EB,0.058_EB,0.068_EB,0.079_EB,0.083_EB,0.087_EB,0.091_EB,0.095_EB,0.091_EB,0.088_EB/)
ABSF(4,19,11:20,4) = (/0.085_EB,0.081_EB,0.084_EB,0.087_EB,0.090_EB,0.093_EB,0.090_EB,0.087_EB,0.085_EB,0.082_EB/)
ABSF(4,19,21:30,4) = (/0.077_EB,0.072_EB,0.067_EB,0.062_EB,0.057_EB,0.058_EB,0.059_EB,0.061_EB,0.062_EB,0.063_EB/)
ABSF(4,19,0:10,5) = (/0.000_EB,0.036_EB,0.068_EB,0.074_EB,0.080_EB,0.087_EB,0.093_EB,0.099_EB,0.105_EB,0.102_EB,0.099_EB/)
ABSF(4,19,11:20,5) = (/0.095_EB,0.092_EB,0.074_EB,0.056_EB,0.038_EB,0.020_EB,0.028_EB,0.037_EB,0.045_EB,0.054_EB/)
ABSF(4,19,21:30,5) = (/0.052_EB,0.049_EB,0.047_EB,0.045_EB,0.043_EB,0.061_EB,0.079_EB,0.097_EB,0.115_EB,0.133_EB/)
ABSF(4,19,0:10,6) = (/0.000_EB,0.036_EB,0.068_EB,0.076_EB,0.083_EB,0.087_EB,0.091_EB,0.095_EB,0.099_EB,0.101_EB,0.102_EB/)
ABSF(4,19,11:20,6) = (/0.103_EB,0.105_EB,0.094_EB,0.083_EB,0.072_EB,0.061_EB,0.055_EB,0.049_EB,0.043_EB,0.036_EB/)
ABSF(4,19,21:30,6) = (/0.046_EB,0.055_EB,0.065_EB,0.074_EB,0.084_EB,0.075_EB,0.067_EB,0.059_EB,0.051_EB,0.042_EB/)
ABSF(4,20,0:10,1) = (/0.000_EB,0.051_EB,0.061_EB,0.070_EB,0.078_EB,0.074_EB,0.070_EB,0.067_EB,0.063_EB,0.060_EB,0.057_EB/)
ABSF(4,20,11:20,1) = (/0.054_EB,0.050_EB,0.053_EB,0.056_EB,0.059_EB,0.061_EB,0.063_EB,0.065_EB,0.067_EB,0.069_EB/)
ABSF(4,20,21:30,1) = (/0.069_EB,0.070_EB,0.071_EB,0.071_EB,0.072_EB,0.072_EB,0.073_EB,0.074_EB,0.074_EB,0.075_EB/)
ABSF(4,20,0:10,2) = (/0.000_EB,0.050_EB,0.061_EB,0.065_EB,0.068_EB,0.067_EB,0.066_EB,0.064_EB,0.063_EB,0.060_EB,0.056_EB/)
ABSF(4,20,11:20,2) = (/0.053_EB,0.050_EB,0.053_EB,0.056_EB,0.059_EB,0.063_EB,0.064_EB,0.066_EB,0.068_EB,0.070_EB/)
ABSF(4,20,21:30,2) = (/0.071_EB,0.071_EB,0.072_EB,0.072_EB,0.073_EB,0.074_EB,0.076_EB,0.077_EB,0.078_EB,0.080_EB/)
ABSF(4,20,0:10,3) = (/0.000_EB,0.036_EB,0.055_EB,0.067_EB,0.078_EB,0.075_EB,0.071_EB,0.068_EB,0.064_EB,0.062_EB,0.059_EB/)
ABSF(4,20,11:20,3) = (/0.057_EB,0.054_EB,0.056_EB,0.058_EB,0.059_EB,0.061_EB,0.062_EB,0.064_EB,0.065_EB,0.066_EB/)
ABSF(4,20,21:30,3) = (/0.068_EB,0.070_EB,0.072_EB,0.074_EB,0.076_EB,0.076_EB,0.077_EB,0.077_EB,0.077_EB,0.077_EB/)
ABSF(4,20,0:10,4) = (/0.000_EB,0.035_EB,0.059_EB,0.070_EB,0.080_EB,0.084_EB,0.088_EB,0.092_EB,0.096_EB,0.092_EB,0.089_EB/)
ABSF(4,20,11:20,4) = (/0.085_EB,0.082_EB,0.085_EB,0.088_EB,0.091_EB,0.094_EB,0.091_EB,0.088_EB,0.085_EB,0.082_EB/)
ABSF(4,20,21:30,4) = (/0.077_EB,0.073_EB,0.069_EB,0.065_EB,0.061_EB,0.060_EB,0.060_EB,0.060_EB,0.059_EB,0.059_EB/)
ABSF(4,20,0:10,5) = (/0.000_EB,0.037_EB,0.070_EB,0.076_EB,0.082_EB,0.089_EB,0.095_EB,0.101_EB,0.107_EB,0.103_EB,0.100_EB/)
ABSF(4,20,11:20,5) = (/0.096_EB,0.092_EB,0.073_EB,0.054_EB,0.035_EB,0.016_EB,0.026_EB,0.035_EB,0.044_EB,0.054_EB/)
ABSF(4,20,21:30,5) = (/0.051_EB,0.047_EB,0.044_EB,0.041_EB,0.038_EB,0.058_EB,0.077_EB,0.097_EB,0.116_EB,0.136_EB/)
ABSF(4,20,0:10,6) = (/0.000_EB,0.037_EB,0.069_EB,0.077_EB,0.085_EB,0.088_EB,0.092_EB,0.096_EB,0.099_EB,0.102_EB,0.104_EB/)
ABSF(4,20,11:20,6) = (/0.106_EB,0.108_EB,0.095_EB,0.083_EB,0.070_EB,0.058_EB,0.051_EB,0.044_EB,0.036_EB,0.029_EB/)
ABSF(4,20,21:30,6) = (/0.040_EB,0.050_EB,0.061_EB,0.071_EB,0.082_EB,0.072_EB,0.062_EB,0.053_EB,0.043_EB,0.034_EB/)
ABSF(5,0,0:10,1) = (/0.000_EB,0.010_EB,0.018_EB,0.020_EB,0.023_EB,0.025_EB,0.027_EB,0.029_EB,0.031_EB,0.030_EB,0.029_EB/)
ABSF(5,0,11:20,1) = (/0.028_EB,0.027_EB,0.028_EB,0.030_EB,0.031_EB,0.032_EB,0.033_EB,0.033_EB,0.034_EB,0.035_EB/)
ABSF(5,0,21:30,1) = (/0.035_EB,0.035_EB,0.036_EB,0.036_EB,0.037_EB,0.036_EB,0.036_EB,0.036_EB,0.036_EB,0.036_EB/)
ABSF(5,0,0:10,2) = (/0.000_EB,0.013_EB,0.020_EB,0.025_EB,0.031_EB,0.033_EB,0.036_EB,0.038_EB,0.041_EB,0.040_EB,0.040_EB/)
ABSF(5,0,11:20,2) = (/0.039_EB,0.039_EB,0.039_EB,0.040_EB,0.040_EB,0.041_EB,0.041_EB,0.041_EB,0.041_EB,0.041_EB/)
ABSF(5,0,21:30,2) = (/0.039_EB,0.038_EB,0.036_EB,0.034_EB,0.032_EB,0.030_EB,0.027_EB,0.024_EB,0.021_EB,0.018_EB/)
ABSF(5,0,0:10,3) = (/0.000_EB,0.016_EB,0.023_EB,0.028_EB,0.032_EB,0.035_EB,0.039_EB,0.042_EB,0.046_EB,0.047_EB,0.049_EB/)
ABSF(5,0,11:20,3) = (/0.050_EB,0.052_EB,0.046_EB,0.041_EB,0.035_EB,0.030_EB,0.030_EB,0.030_EB,0.031_EB,0.031_EB/)
ABSF(5,0,21:30,3) = (/0.032_EB,0.032_EB,0.032_EB,0.033_EB,0.033_EB,0.034_EB,0.034_EB,0.035_EB,0.035_EB,0.036_EB/)
ABSF(5,0,0:10,4) = (/0.000_EB,0.015_EB,0.030_EB,0.033_EB,0.036_EB,0.039_EB,0.043_EB,0.046_EB,0.050_EB,0.054_EB,0.058_EB/)
ABSF(5,0,11:20,4) = (/0.062_EB,0.067_EB,0.057_EB,0.048_EB,0.039_EB,0.030_EB,0.026_EB,0.022_EB,0.019_EB,0.015_EB/)
ABSF(5,0,21:30,4) = (/0.014_EB,0.013_EB,0.013_EB,0.012_EB,0.011_EB,0.009_EB,0.006_EB,0.004_EB,0.002_EB,-0.000_EB/)
ABSF(5,0,0:10,5) = (/0.000_EB,0.016_EB,0.035_EB,0.042_EB,0.049_EB,0.051_EB,0.052_EB,0.054_EB,0.056_EB,0.059_EB,0.063_EB/)
ABSF(5,0,11:20,5) = (/0.067_EB,0.070_EB,0.058_EB,0.046_EB,0.034_EB,0.021_EB,0.016_EB,0.010_EB,0.004_EB,-0.002_EB/)
ABSF(5,0,21:30,5) = (/0.000_EB,0.003_EB,0.005_EB,0.007_EB,0.010_EB,0.005_EB,0.001_EB,-0.004_EB,-0.009_EB,-0.013_EB/)
ABSF(5,0,0:10,6) = (/0.000_EB,0.021_EB,0.038_EB,0.043_EB,0.047_EB,0.045_EB,0.043_EB,0.041_EB,0.039_EB,0.050_EB,0.061_EB/)
ABSF(5,0,11:20,6) = (/0.072_EB,0.083_EB,0.077_EB,0.072_EB,0.066_EB,0.061_EB,0.042_EB,0.023_EB,0.004_EB,-0.016_EB/)
ABSF(5,0,21:30,6) = (/-0.013_EB,-0.010_EB,-0.007_EB,-0.004_EB,-0.001_EB,0.000_EB,0.002_EB,0.004_EB,0.005_EB,0.007_EB/)
ABSF(5,1,0:10,1) = (/0.000_EB,0.015_EB,0.022_EB,0.023_EB,0.025_EB,0.026_EB,0.027_EB,0.028_EB,0.029_EB,0.028_EB,0.027_EB/)
ABSF(5,1,11:20,1) = (/0.026_EB,0.026_EB,0.027_EB,0.029_EB,0.030_EB,0.032_EB,0.033_EB,0.034_EB,0.035_EB,0.036_EB/)
ABSF(5,1,21:30,1) = (/0.036_EB,0.037_EB,0.037_EB,0.037_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB/)
ABSF(5,1,0:10,2) = (/0.000_EB,0.014_EB,0.021_EB,0.026_EB,0.031_EB,0.034_EB,0.036_EB,0.038_EB,0.040_EB,0.039_EB,0.039_EB/)
ABSF(5,1,11:20,2) = (/0.039_EB,0.039_EB,0.039_EB,0.039_EB,0.039_EB,0.039_EB,0.038_EB,0.038_EB,0.037_EB,0.037_EB/)
ABSF(5,1,21:30,2) = (/0.036_EB,0.036_EB,0.036_EB,0.036_EB,0.036_EB,0.034_EB,0.032_EB,0.030_EB,0.028_EB,0.026_EB/)
ABSF(5,1,0:10,3) = (/0.000_EB,0.017_EB,0.026_EB,0.028_EB,0.031_EB,0.035_EB,0.039_EB,0.044_EB,0.048_EB,0.048_EB,0.048_EB/)
ABSF(5,1,11:20,3) = (/0.048_EB,0.048_EB,0.044_EB,0.039_EB,0.035_EB,0.030_EB,0.031_EB,0.033_EB,0.034_EB,0.035_EB/)
ABSF(5,1,21:30,3) = (/0.037_EB,0.038_EB,0.039_EB,0.041_EB,0.042_EB,0.035_EB,0.029_EB,0.023_EB,0.016_EB,0.010_EB/)
ABSF(5,1,0:10,4) = (/0.000_EB,0.017_EB,0.030_EB,0.033_EB,0.035_EB,0.038_EB,0.041_EB,0.044_EB,0.047_EB,0.051_EB,0.054_EB/)
ABSF(5,1,11:20,4) = (/0.058_EB,0.062_EB,0.053_EB,0.044_EB,0.035_EB,0.027_EB,0.026_EB,0.025_EB,0.025_EB,0.024_EB/)
ABSF(5,1,21:30,4) = (/0.025_EB,0.027_EB,0.028_EB,0.029_EB,0.030_EB,0.024_EB,0.018_EB,0.012_EB,0.006_EB,-0.000_EB/)
ABSF(5,1,0:10,5) = (/0.000_EB,0.017_EB,0.035_EB,0.043_EB,0.051_EB,0.053_EB,0.054_EB,0.056_EB,0.058_EB,0.061_EB,0.065_EB/)
ABSF(5,1,11:20,5) = (/0.069_EB,0.072_EB,0.062_EB,0.051_EB,0.040_EB,0.029_EB,0.019_EB,0.008_EB,-0.002_EB,-0.012_EB/)
ABSF(5,1,21:30,5) = (/-0.010_EB,-0.008_EB,-0.005_EB,-0.003_EB,-0.001_EB,0.003_EB,0.008_EB,0.012_EB,0.016_EB,0.020_EB/)
ABSF(5,1,0:10,6) = (/0.000_EB,0.022_EB,0.041_EB,0.045_EB,0.049_EB,0.046_EB,0.044_EB,0.041_EB,0.039_EB,0.050_EB,0.061_EB/)
ABSF(5,1,11:20,6) = (/0.072_EB,0.084_EB,0.075_EB,0.067_EB,0.059_EB,0.051_EB,0.035_EB,0.019_EB,0.003_EB,-0.013_EB/)
ABSF(5,1,21:30,6) = (/-0.013_EB,-0.013_EB,-0.012_EB,-0.012_EB,-0.012_EB,-0.013_EB,-0.015_EB,-0.016_EB,-0.017_EB,-0.019_EB/)
ABSF(5,2,0:10,1) = (/0.000_EB,0.021_EB,0.025_EB,0.026_EB,0.028_EB,0.029_EB,0.030_EB,0.031_EB,0.032_EB,0.030_EB,0.029_EB/)
ABSF(5,2,11:20,1) = (/0.028_EB,0.026_EB,0.028_EB,0.029_EB,0.030_EB,0.032_EB,0.033_EB,0.034_EB,0.035_EB,0.036_EB/)
ABSF(5,2,21:30,1) = (/0.036_EB,0.037_EB,0.037_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB/)
ABSF(5,2,0:10,2) = (/0.000_EB,0.016_EB,0.023_EB,0.029_EB,0.034_EB,0.037_EB,0.039_EB,0.041_EB,0.044_EB,0.041_EB,0.039_EB/)
ABSF(5,2,11:20,2) = (/0.037_EB,0.035_EB,0.035_EB,0.036_EB,0.036_EB,0.036_EB,0.036_EB,0.036_EB,0.037_EB,0.037_EB/)
ABSF(5,2,21:30,2) = (/0.037_EB,0.037_EB,0.037_EB,0.036_EB,0.036_EB,0.035_EB,0.033_EB,0.031_EB,0.030_EB,0.028_EB/)
ABSF(5,2,0:10,3) = (/0.000_EB,0.015_EB,0.028_EB,0.031_EB,0.034_EB,0.038_EB,0.042_EB,0.045_EB,0.049_EB,0.050_EB,0.052_EB/)
ABSF(5,2,11:20,3) = (/0.053_EB,0.055_EB,0.048_EB,0.041_EB,0.034_EB,0.027_EB,0.028_EB,0.029_EB,0.030_EB,0.031_EB/)
ABSF(5,2,21:30,3) = (/0.032_EB,0.033_EB,0.034_EB,0.035_EB,0.036_EB,0.035_EB,0.034_EB,0.033_EB,0.032_EB,0.031_EB/)
ABSF(5,2,0:10,4) = (/0.000_EB,0.018_EB,0.033_EB,0.035_EB,0.037_EB,0.041_EB,0.045_EB,0.048_EB,0.052_EB,0.055_EB,0.058_EB/)
ABSF(5,2,11:20,4) = (/0.062_EB,0.065_EB,0.057_EB,0.050_EB,0.042_EB,0.035_EB,0.028_EB,0.022_EB,0.016_EB,0.010_EB/)
ABSF(5,2,21:30,4) = (/0.013_EB,0.015_EB,0.018_EB,0.021_EB,0.024_EB,0.016_EB,0.008_EB,0.000_EB,-0.008_EB,-0.016_EB/)
ABSF(5,2,0:10,5) = (/0.000_EB,0.017_EB,0.036_EB,0.044_EB,0.053_EB,0.054_EB,0.055_EB,0.056_EB,0.057_EB,0.061_EB,0.066_EB/)
ABSF(5,2,11:20,5) = (/0.071_EB,0.076_EB,0.060_EB,0.045_EB,0.030_EB,0.015_EB,0.007_EB,-0.002_EB,-0.010_EB,-0.018_EB/)
ABSF(5,2,21:30,5) = (/-0.011_EB,-0.004_EB,0.004_EB,0.011_EB,0.018_EB,0.016_EB,0.015_EB,0.013_EB,0.012_EB,0.010_EB/)
ABSF(5,2,0:10,6) = (/0.000_EB,0.023_EB,0.041_EB,0.045_EB,0.049_EB,0.047_EB,0.045_EB,0.042_EB,0.040_EB,0.051_EB,0.062_EB/)
ABSF(5,2,11:20,6) = (/0.073_EB,0.084_EB,0.078_EB,0.072_EB,0.065_EB,0.059_EB,0.040_EB,0.021_EB,0.002_EB,-0.017_EB/)
ABSF(5,2,21:30,6) = (/-0.016_EB,-0.014_EB,-0.012_EB,-0.010_EB,-0.009_EB,-0.009_EB,-0.010_EB,-0.010_EB,-0.010_EB,-0.011_EB/)
ABSF(5,3,0:10,1) = (/0.000_EB,0.023_EB,0.026_EB,0.028_EB,0.029_EB,0.030_EB,0.031_EB,0.032_EB,0.033_EB,0.031_EB,0.030_EB/)
ABSF(5,3,11:20,1) = (/0.028_EB,0.027_EB,0.028_EB,0.029_EB,0.031_EB,0.032_EB,0.033_EB,0.035_EB,0.036_EB,0.037_EB/)
ABSF(5,3,21:30,1) = (/0.038_EB,0.038_EB,0.039_EB,0.039_EB,0.040_EB,0.040_EB,0.040_EB,0.040_EB,0.040_EB,0.040_EB/)
ABSF(5,3,0:10,2) = (/0.000_EB,0.017_EB,0.024_EB,0.029_EB,0.033_EB,0.035_EB,0.038_EB,0.040_EB,0.042_EB,0.040_EB,0.038_EB/)
ABSF(5,3,11:20,2) = (/0.036_EB,0.033_EB,0.034_EB,0.034_EB,0.034_EB,0.035_EB,0.036_EB,0.036_EB,0.037_EB,0.038_EB/)
ABSF(5,3,21:30,2) = (/0.038_EB,0.038_EB,0.038_EB,0.038_EB,0.039_EB,0.037_EB,0.036_EB,0.035_EB,0.034_EB,0.032_EB/)
ABSF(5,3,0:10,3) = (/0.000_EB,0.017_EB,0.029_EB,0.034_EB,0.039_EB,0.042_EB,0.045_EB,0.047_EB,0.050_EB,0.051_EB,0.053_EB/)
ABSF(5,3,11:20,3) = (/0.054_EB,0.055_EB,0.048_EB,0.041_EB,0.035_EB,0.028_EB,0.029_EB,0.030_EB,0.031_EB,0.033_EB/)
ABSF(5,3,21:30,3) = (/0.034_EB,0.035_EB,0.037_EB,0.038_EB,0.039_EB,0.038_EB,0.036_EB,0.034_EB,0.032_EB,0.030_EB/)
ABSF(5,3,0:10,4) = (/0.000_EB,0.019_EB,0.034_EB,0.036_EB,0.039_EB,0.042_EB,0.046_EB,0.050_EB,0.054_EB,0.056_EB,0.059_EB/)
ABSF(5,3,11:20,4) = (/0.062_EB,0.065_EB,0.058_EB,0.051_EB,0.044_EB,0.037_EB,0.033_EB,0.030_EB,0.026_EB,0.022_EB/)
ABSF(5,3,21:30,4) = (/0.023_EB,0.024_EB,0.024_EB,0.025_EB,0.026_EB,0.019_EB,0.012_EB,0.005_EB,-0.002_EB,-0.009_EB/)
ABSF(5,3,0:10,5) = (/0.000_EB,0.019_EB,0.038_EB,0.046_EB,0.054_EB,0.055_EB,0.057_EB,0.058_EB,0.060_EB,0.064_EB,0.068_EB/)
ABSF(5,3,11:20,5) = (/0.072_EB,0.076_EB,0.064_EB,0.052_EB,0.040_EB,0.029_EB,0.018_EB,0.008_EB,-0.002_EB,-0.013_EB/)
ABSF(5,3,21:30,5) = (/-0.007_EB,-0.001_EB,0.005_EB,0.011_EB,0.017_EB,0.016_EB,0.014_EB,0.012_EB,0.011_EB,0.009_EB/)
ABSF(5,3,0:10,6) = (/0.000_EB,0.024_EB,0.043_EB,0.047_EB,0.051_EB,0.048_EB,0.046_EB,0.043_EB,0.041_EB,0.052_EB,0.063_EB/)
ABSF(5,3,11:20,6) = (/0.074_EB,0.085_EB,0.079_EB,0.072_EB,0.066_EB,0.059_EB,0.041_EB,0.023_EB,0.005_EB,-0.013_EB/)
ABSF(5,3,21:30,6) = (/-0.011_EB,-0.009_EB,-0.007_EB,-0.006_EB,-0.004_EB,-0.004_EB,-0.005_EB,-0.006_EB,-0.007_EB,-0.007_EB/)
ABSF(5,4,0:10,1) = (/0.000_EB,0.025_EB,0.028_EB,0.029_EB,0.030_EB,0.031_EB,0.032_EB,0.033_EB,0.033_EB,0.032_EB,0.030_EB/)
ABSF(5,4,11:20,1) = (/0.029_EB,0.027_EB,0.028_EB,0.029_EB,0.031_EB,0.032_EB,0.033_EB,0.035_EB,0.037_EB,0.039_EB/)
ABSF(5,4,21:30,1) = (/0.039_EB,0.040_EB,0.040_EB,0.041_EB,0.041_EB,0.041_EB,0.041_EB,0.042_EB,0.042_EB,0.042_EB/)
ABSF(5,4,0:10,2) = (/0.000_EB,0.018_EB,0.025_EB,0.028_EB,0.032_EB,0.034_EB,0.037_EB,0.039_EB,0.041_EB,0.039_EB,0.036_EB/)
ABSF(5,4,11:20,2) = (/0.034_EB,0.032_EB,0.032_EB,0.032_EB,0.033_EB,0.033_EB,0.035_EB,0.036_EB,0.038_EB,0.039_EB/)
ABSF(5,4,21:30,2) = (/0.039_EB,0.040_EB,0.040_EB,0.040_EB,0.041_EB,0.040_EB,0.039_EB,0.038_EB,0.037_EB,0.037_EB/)
ABSF(5,4,0:10,3) = (/0.000_EB,0.018_EB,0.029_EB,0.036_EB,0.043_EB,0.045_EB,0.048_EB,0.050_EB,0.052_EB,0.053_EB,0.053_EB/)
ABSF(5,4,11:20,3) = (/0.054_EB,0.055_EB,0.048_EB,0.042_EB,0.035_EB,0.028_EB,0.030_EB,0.031_EB,0.033_EB,0.034_EB/)
ABSF(5,4,21:30,3) = (/0.036_EB,0.038_EB,0.039_EB,0.041_EB,0.043_EB,0.040_EB,0.038_EB,0.035_EB,0.033_EB,0.030_EB/)
ABSF(5,4,0:10,4) = (/0.000_EB,0.019_EB,0.035_EB,0.038_EB,0.040_EB,0.044_EB,0.048_EB,0.051_EB,0.055_EB,0.058_EB,0.060_EB/)
ABSF(5,4,11:20,4) = (/0.063_EB,0.066_EB,0.059_EB,0.052_EB,0.046_EB,0.039_EB,0.038_EB,0.037_EB,0.036_EB,0.035_EB/)
ABSF(5,4,21:30,4) = (/0.033_EB,0.032_EB,0.030_EB,0.029_EB,0.027_EB,0.021_EB,0.016_EB,0.010_EB,0.004_EB,-0.002_EB/)
ABSF(5,4,0:10,5) = (/0.000_EB,0.020_EB,0.040_EB,0.047_EB,0.054_EB,0.056_EB,0.058_EB,0.060_EB,0.062_EB,0.066_EB,0.069_EB/)
ABSF(5,4,11:20,5) = (/0.073_EB,0.076_EB,0.068_EB,0.059_EB,0.051_EB,0.042_EB,0.030_EB,0.017_EB,0.005_EB,-0.008_EB/)
ABSF(5,4,21:30,5) = (/-0.003_EB,0.002_EB,0.007_EB,0.012_EB,0.016_EB,0.015_EB,0.013_EB,0.012_EB,0.010_EB,0.008_EB/)
ABSF(5,4,0:10,6) = (/0.000_EB,0.024_EB,0.044_EB,0.048_EB,0.053_EB,0.050_EB,0.047_EB,0.045_EB,0.042_EB,0.053_EB,0.064_EB/)
ABSF(5,4,11:20,6) = (/0.075_EB,0.086_EB,0.079_EB,0.073_EB,0.066_EB,0.059_EB,0.042_EB,0.026_EB,0.009_EB,-0.008_EB/)
ABSF(5,4,21:30,6) = (/-0.006_EB,-0.004_EB,-0.002_EB,-0.001_EB,0.001_EB,0.000_EB,-0.001_EB,-0.002_EB,-0.003_EB,-0.004_EB/)
ABSF(5,5,0:10,1) = (/0.000_EB,0.026_EB,0.029_EB,0.030_EB,0.031_EB,0.032_EB,0.033_EB,0.033_EB,0.034_EB,0.033_EB,0.031_EB/)
ABSF(5,5,11:20,1) = (/0.029_EB,0.028_EB,0.029_EB,0.030_EB,0.031_EB,0.032_EB,0.034_EB,0.036_EB,0.038_EB,0.040_EB/)
ABSF(5,5,21:30,1) = (/0.041_EB,0.041_EB,0.042_EB,0.042_EB,0.043_EB,0.043_EB,0.043_EB,0.043_EB,0.043_EB,0.043_EB/)
ABSF(5,5,0:10,2) = (/0.000_EB,0.019_EB,0.025_EB,0.028_EB,0.031_EB,0.033_EB,0.035_EB,0.038_EB,0.040_EB,0.037_EB,0.035_EB/)
ABSF(5,5,11:20,2) = (/0.032_EB,0.030_EB,0.030_EB,0.031_EB,0.031_EB,0.032_EB,0.034_EB,0.036_EB,0.038_EB,0.040_EB/)
ABSF(5,5,21:30,2) = (/0.041_EB,0.041_EB,0.042_EB,0.042_EB,0.043_EB,0.042_EB,0.042_EB,0.042_EB,0.041_EB,0.041_EB/)
ABSF(5,5,0:10,3) = (/0.000_EB,0.020_EB,0.030_EB,0.039_EB,0.048_EB,0.049_EB,0.050_EB,0.052_EB,0.053_EB,0.054_EB,0.054_EB/)
ABSF(5,5,11:20,3) = (/0.055_EB,0.055_EB,0.049_EB,0.042_EB,0.035_EB,0.028_EB,0.030_EB,0.032_EB,0.034_EB,0.035_EB/)
ABSF(5,5,21:30,3) = (/0.038_EB,0.040_EB,0.042_EB,0.044_EB,0.046_EB,0.043_EB,0.040_EB,0.036_EB,0.033_EB,0.030_EB/)
ABSF(5,5,0:10,4) = (/0.000_EB,0.020_EB,0.037_EB,0.039_EB,0.042_EB,0.046_EB,0.049_EB,0.053_EB,0.056_EB,0.059_EB,0.061_EB/)
ABSF(5,5,11:20,4) = (/0.064_EB,0.067_EB,0.060_EB,0.054_EB,0.047_EB,0.041_EB,0.043_EB,0.044_EB,0.046_EB,0.048_EB/)
ABSF(5,5,21:30,4) = (/0.044_EB,0.040_EB,0.036_EB,0.033_EB,0.029_EB,0.024_EB,0.019_EB,0.015_EB,0.010_EB,0.005_EB/)
ABSF(5,5,0:10,5) = (/0.000_EB,0.022_EB,0.042_EB,0.049_EB,0.055_EB,0.058_EB,0.060_EB,0.063_EB,0.065_EB,0.068_EB,0.071_EB/)
ABSF(5,5,11:20,5) = (/0.074_EB,0.077_EB,0.072_EB,0.066_EB,0.061_EB,0.056_EB,0.042_EB,0.027_EB,0.012_EB,-0.003_EB/)
ABSF(5,5,21:30,5) = (/0.001_EB,0.005_EB,0.008_EB,0.012_EB,0.016_EB,0.014_EB,0.012_EB,0.011_EB,0.009_EB,0.007_EB/)
ABSF(5,5,0:10,6) = (/0.000_EB,0.025_EB,0.045_EB,0.050_EB,0.055_EB,0.052_EB,0.049_EB,0.046_EB,0.043_EB,0.054_EB,0.065_EB/)
ABSF(5,5,11:20,6) = (/0.076_EB,0.087_EB,0.080_EB,0.073_EB,0.066_EB,0.060_EB,0.044_EB,0.028_EB,0.012_EB,-0.004_EB/)
ABSF(5,5,21:30,6) = (/-0.002_EB,0.000_EB,0.002_EB,0.004_EB,0.006_EB,0.005_EB,0.004_EB,0.002_EB,0.001_EB,-0.000_EB/)
ABSF(5,6,0:10,1) = (/0.000_EB,0.027_EB,0.031_EB,0.032_EB,0.033_EB,0.033_EB,0.034_EB,0.035_EB,0.035_EB,0.034_EB,0.032_EB/)
ABSF(5,6,11:20,1) = (/0.031_EB,0.029_EB,0.030_EB,0.031_EB,0.031_EB,0.032_EB,0.034_EB,0.037_EB,0.039_EB,0.041_EB/)
ABSF(5,6,21:30,1) = (/0.042_EB,0.042_EB,0.043_EB,0.043_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.045_EB,0.045_EB/)
ABSF(5,6,0:10,2) = (/0.000_EB,0.020_EB,0.028_EB,0.030_EB,0.032_EB,0.034_EB,0.036_EB,0.039_EB,0.041_EB,0.038_EB,0.036_EB/)
ABSF(5,6,11:20,2) = (/0.033_EB,0.031_EB,0.031_EB,0.031_EB,0.032_EB,0.032_EB,0.035_EB,0.037_EB,0.039_EB,0.042_EB/)
ABSF(5,6,21:30,2) = (/0.042_EB,0.043_EB,0.043_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.043_EB/)
ABSF(5,6,0:10,3) = (/0.000_EB,0.021_EB,0.031_EB,0.040_EB,0.048_EB,0.050_EB,0.051_EB,0.053_EB,0.054_EB,0.054_EB,0.054_EB/)
ABSF(5,6,11:20,3) = (/0.054_EB,0.054_EB,0.049_EB,0.043_EB,0.037_EB,0.032_EB,0.033_EB,0.035_EB,0.036_EB,0.038_EB/)
ABSF(5,6,21:30,3) = (/0.039_EB,0.041_EB,0.043_EB,0.045_EB,0.046_EB,0.043_EB,0.040_EB,0.037_EB,0.034_EB,0.031_EB/)
ABSF(5,6,0:10,4) = (/0.000_EB,0.021_EB,0.038_EB,0.041_EB,0.044_EB,0.047_EB,0.051_EB,0.054_EB,0.058_EB,0.061_EB,0.064_EB/)
ABSF(5,6,11:20,4) = (/0.067_EB,0.070_EB,0.063_EB,0.056_EB,0.050_EB,0.043_EB,0.043_EB,0.043_EB,0.043_EB,0.043_EB/)
ABSF(5,6,21:30,4) = (/0.040_EB,0.038_EB,0.035_EB,0.032_EB,0.029_EB,0.026_EB,0.022_EB,0.018_EB,0.014_EB,0.010_EB/)
ABSF(5,6,0:10,5) = (/0.000_EB,0.023_EB,0.044_EB,0.051_EB,0.057_EB,0.059_EB,0.062_EB,0.064_EB,0.067_EB,0.069_EB,0.072_EB/)
ABSF(5,6,11:20,5) = (/0.075_EB,0.078_EB,0.073_EB,0.068_EB,0.062_EB,0.057_EB,0.042_EB,0.027_EB,0.012_EB,-0.003_EB/)
ABSF(5,6,21:30,5) = (/0.001_EB,0.005_EB,0.009_EB,0.013_EB,0.017_EB,0.016_EB,0.015_EB,0.013_EB,0.012_EB,0.011_EB/)
ABSF(5,6,0:10,6) = (/0.000_EB,0.026_EB,0.047_EB,0.052_EB,0.057_EB,0.054_EB,0.050_EB,0.047_EB,0.044_EB,0.056_EB,0.067_EB/)
ABSF(5,6,11:20,6) = (/0.078_EB,0.089_EB,0.082_EB,0.075_EB,0.068_EB,0.061_EB,0.046_EB,0.031_EB,0.015_EB,0.000_EB/)
ABSF(5,6,21:30,6) = (/0.002_EB,0.004_EB,0.006_EB,0.007_EB,0.009_EB,0.007_EB,0.005_EB,0.004_EB,0.002_EB,-0.000_EB/)
ABSF(5,7,0:10,1) = (/0.000_EB,0.029_EB,0.032_EB,0.033_EB,0.034_EB,0.035_EB,0.035_EB,0.036_EB,0.037_EB,0.035_EB,0.034_EB/)
ABSF(5,7,11:20,1) = (/0.032_EB,0.031_EB,0.031_EB,0.032_EB,0.032_EB,0.033_EB,0.035_EB,0.037_EB,0.040_EB,0.042_EB/)
ABSF(5,7,21:30,1) = (/0.043_EB,0.043_EB,0.044_EB,0.044_EB,0.045_EB,0.045_EB,0.045_EB,0.046_EB,0.046_EB,0.047_EB/)
ABSF(5,7,0:10,2) = (/0.000_EB,0.021_EB,0.030_EB,0.032_EB,0.033_EB,0.035_EB,0.037_EB,0.039_EB,0.041_EB,0.039_EB,0.036_EB/)
ABSF(5,7,11:20,2) = (/0.034_EB,0.031_EB,0.032_EB,0.032_EB,0.032_EB,0.033_EB,0.035_EB,0.038_EB,0.040_EB,0.043_EB/)
ABSF(5,7,21:30,2) = (/0.043_EB,0.044_EB,0.044_EB,0.045_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB,0.046_EB/)
ABSF(5,7,0:10,3) = (/0.000_EB,0.021_EB,0.033_EB,0.041_EB,0.049_EB,0.050_EB,0.052_EB,0.053_EB,0.055_EB,0.054_EB,0.054_EB/)
ABSF(5,7,11:20,3) = (/0.053_EB,0.053_EB,0.048_EB,0.044_EB,0.039_EB,0.035_EB,0.036_EB,0.037_EB,0.039_EB,0.040_EB/)
ABSF(5,7,21:30,3) = (/0.041_EB,0.042_EB,0.044_EB,0.045_EB,0.047_EB,0.044_EB,0.041_EB,0.039_EB,0.036_EB,0.033_EB/)
ABSF(5,7,0:10,4) = (/0.000_EB,0.022_EB,0.039_EB,0.043_EB,0.046_EB,0.049_EB,0.053_EB,0.056_EB,0.059_EB,0.063_EB,0.066_EB/)
ABSF(5,7,11:20,4) = (/0.070_EB,0.074_EB,0.066_EB,0.059_EB,0.052_EB,0.045_EB,0.043_EB,0.042_EB,0.040_EB,0.039_EB/)
ABSF(5,7,21:30,4) = (/0.037_EB,0.035_EB,0.033_EB,0.032_EB,0.030_EB,0.027_EB,0.024_EB,0.021_EB,0.018_EB,0.014_EB/)
ABSF(5,7,0:10,5) = (/0.000_EB,0.024_EB,0.047_EB,0.052_EB,0.058_EB,0.061_EB,0.063_EB,0.066_EB,0.068_EB,0.071_EB,0.073_EB/)
ABSF(5,7,11:20,5) = (/0.076_EB,0.079_EB,0.074_EB,0.069_EB,0.063_EB,0.058_EB,0.043_EB,0.028_EB,0.012_EB,-0.003_EB/)
ABSF(5,7,21:30,5) = (/0.002_EB,0.006_EB,0.010_EB,0.015_EB,0.019_EB,0.018_EB,0.017_EB,0.016_EB,0.015_EB,0.014_EB/)
ABSF(5,7,0:10,6) = (/0.000_EB,0.027_EB,0.048_EB,0.053_EB,0.059_EB,0.055_EB,0.052_EB,0.049_EB,0.046_EB,0.057_EB,0.068_EB/)
ABSF(5,7,11:20,6) = (/0.080_EB,0.091_EB,0.084_EB,0.077_EB,0.070_EB,0.063_EB,0.048_EB,0.033_EB,0.019_EB,0.004_EB/)
ABSF(5,7,21:30,6) = (/0.006_EB,0.007_EB,0.009_EB,0.010_EB,0.012_EB,0.009_EB,0.007_EB,0.005_EB,0.002_EB,-0.000_EB/)
ABSF(5,8,0:10,1) = (/0.000_EB,0.030_EB,0.034_EB,0.035_EB,0.035_EB,0.036_EB,0.037_EB,0.037_EB,0.038_EB,0.037_EB,0.035_EB/)
ABSF(5,8,11:20,1) = (/0.034_EB,0.032_EB,0.032_EB,0.033_EB,0.033_EB,0.033_EB,0.035_EB,0.038_EB,0.040_EB,0.043_EB/)
ABSF(5,8,21:30,1) = (/0.043_EB,0.044_EB,0.044_EB,0.045_EB,0.045_EB,0.046_EB,0.047_EB,0.047_EB,0.048_EB,0.048_EB/)
ABSF(5,8,0:10,2) = (/0.000_EB,0.022_EB,0.032_EB,0.034_EB,0.035_EB,0.037_EB,0.038_EB,0.040_EB,0.042_EB,0.040_EB,0.037_EB/)
ABSF(5,8,11:20,2) = (/0.035_EB,0.032_EB,0.032_EB,0.033_EB,0.033_EB,0.033_EB,0.036_EB,0.038_EB,0.041_EB,0.044_EB/)
ABSF(5,8,21:30,2) = (/0.044_EB,0.045_EB,0.046_EB,0.046_EB,0.047_EB,0.047_EB,0.048_EB,0.048_EB,0.048_EB,0.049_EB/)
ABSF(5,8,0:10,3) = (/0.000_EB,0.022_EB,0.035_EB,0.042_EB,0.049_EB,0.051_EB,0.052_EB,0.054_EB,0.056_EB,0.055_EB,0.054_EB/)
ABSF(5,8,11:20,3) = (/0.053_EB,0.052_EB,0.048_EB,0.045_EB,0.042_EB,0.038_EB,0.039_EB,0.040_EB,0.041_EB,0.042_EB/)
ABSF(5,8,21:30,3) = (/0.043_EB,0.044_EB,0.045_EB,0.046_EB,0.047_EB,0.044_EB,0.042_EB,0.040_EB,0.037_EB,0.035_EB/)
ABSF(5,8,0:10,4) = (/0.000_EB,0.023_EB,0.041_EB,0.044_EB,0.048_EB,0.051_EB,0.054_EB,0.058_EB,0.061_EB,0.065_EB,0.069_EB/)
ABSF(5,8,11:20,4) = (/0.073_EB,0.077_EB,0.069_EB,0.062_EB,0.054_EB,0.046_EB,0.043_EB,0.040_EB,0.037_EB,0.034_EB/)
ABSF(5,8,21:30,4) = (/0.033_EB,0.033_EB,0.032_EB,0.031_EB,0.031_EB,0.028_EB,0.026_EB,0.024_EB,0.021_EB,0.019_EB/)
ABSF(5,8,0:10,5) = (/0.000_EB,0.025_EB,0.049_EB,0.054_EB,0.060_EB,0.062_EB,0.065_EB,0.067_EB,0.069_EB,0.072_EB,0.075_EB/)
ABSF(5,8,11:20,5) = (/0.077_EB,0.080_EB,0.075_EB,0.070_EB,0.064_EB,0.059_EB,0.044_EB,0.028_EB,0.013_EB,-0.003_EB/)
ABSF(5,8,21:30,5) = (/0.002_EB,0.007_EB,0.011_EB,0.016_EB,0.021_EB,0.020_EB,0.019_EB,0.019_EB,0.018_EB,0.018_EB/)
ABSF(5,8,0:10,6) = (/0.000_EB,0.028_EB,0.050_EB,0.055_EB,0.061_EB,0.057_EB,0.054_EB,0.051_EB,0.048_EB,0.059_EB,0.070_EB/)
ABSF(5,8,11:20,6) = (/0.082_EB,0.093_EB,0.086_EB,0.079_EB,0.072_EB,0.064_EB,0.050_EB,0.036_EB,0.022_EB,0.008_EB/)
ABSF(5,8,21:30,6) = (/0.009_EB,0.010_EB,0.012_EB,0.013_EB,0.014_EB,0.012_EB,0.009_EB,0.006_EB,0.003_EB,0.000_EB/)
ABSF(5,9,0:10,1) = (/0.000_EB,0.031_EB,0.036_EB,0.036_EB,0.037_EB,0.037_EB,0.038_EB,0.039_EB,0.039_EB,0.038_EB,0.036_EB/)
ABSF(5,9,11:20,1) = (/0.035_EB,0.034_EB,0.034_EB,0.034_EB,0.034_EB,0.033_EB,0.036_EB,0.039_EB,0.041_EB,0.044_EB/)
ABSF(5,9,21:30,1) = (/0.044_EB,0.045_EB,0.045_EB,0.046_EB,0.046_EB,0.047_EB,0.048_EB,0.048_EB,0.049_EB,0.050_EB/)
ABSF(5,9,0:10,2) = (/0.000_EB,0.023_EB,0.035_EB,0.035_EB,0.036_EB,0.038_EB,0.039_EB,0.041_EB,0.043_EB,0.040_EB,0.038_EB/)
ABSF(5,9,11:20,2) = (/0.035_EB,0.033_EB,0.033_EB,0.033_EB,0.033_EB,0.033_EB,0.036_EB,0.039_EB,0.042_EB,0.045_EB/)
ABSF(5,9,21:30,2) = (/0.046_EB,0.046_EB,0.047_EB,0.048_EB,0.049_EB,0.049_EB,0.050_EB,0.050_EB,0.051_EB,0.051_EB/)
ABSF(5,9,0:10,3) = (/0.000_EB,0.022_EB,0.037_EB,0.043_EB,0.050_EB,0.051_EB,0.053_EB,0.055_EB,0.056_EB,0.055_EB,0.053_EB/)
ABSF(5,9,11:20,3) = (/0.052_EB,0.050_EB,0.048_EB,0.046_EB,0.044_EB,0.042_EB,0.042_EB,0.043_EB,0.043_EB,0.044_EB/)
ABSF(5,9,21:30,3) = (/0.045_EB,0.045_EB,0.046_EB,0.046_EB,0.047_EB,0.045_EB,0.043_EB,0.041_EB,0.039_EB,0.037_EB/)
ABSF(5,9,0:10,4) = (/0.000_EB,0.024_EB,0.042_EB,0.046_EB,0.050_EB,0.053_EB,0.056_EB,0.059_EB,0.062_EB,0.067_EB,0.071_EB/)
ABSF(5,9,11:20,4) = (/0.076_EB,0.080_EB,0.072_EB,0.064_EB,0.056_EB,0.048_EB,0.044_EB,0.039_EB,0.034_EB,0.029_EB/)
ABSF(5,9,21:30,4) = (/0.030_EB,0.030_EB,0.031_EB,0.031_EB,0.031_EB,0.030_EB,0.028_EB,0.027_EB,0.025_EB,0.024_EB/)
ABSF(5,9,0:10,5) = (/0.000_EB,0.025_EB,0.051_EB,0.056_EB,0.061_EB,0.064_EB,0.066_EB,0.068_EB,0.071_EB,0.073_EB,0.076_EB/)
ABSF(5,9,11:20,5) = (/0.079_EB,0.081_EB,0.076_EB,0.071_EB,0.066_EB,0.060_EB,0.044_EB,0.029_EB,0.013_EB,-0.003_EB/)
ABSF(5,9,21:30,5) = (/0.002_EB,0.007_EB,0.012_EB,0.017_EB,0.022_EB,0.022_EB,0.022_EB,0.022_EB,0.021_EB,0.021_EB/)
ABSF(5,9,0:10,6) = (/0.000_EB,0.029_EB,0.051_EB,0.057_EB,0.062_EB,0.059_EB,0.056_EB,0.053_EB,0.049_EB,0.061_EB,0.072_EB/)
ABSF(5,9,11:20,6) = (/0.084_EB,0.095_EB,0.088_EB,0.081_EB,0.073_EB,0.066_EB,0.052_EB,0.039_EB,0.025_EB,0.012_EB/)
ABSF(5,9,21:30,6) = (/0.013_EB,0.014_EB,0.015_EB,0.016_EB,0.017_EB,0.014_EB,0.010_EB,0.007_EB,0.004_EB,0.000_EB/)
ABSF(5,10,0:10,1) = (/0.000_EB,0.032_EB,0.037_EB,0.038_EB,0.038_EB,0.039_EB,0.039_EB,0.040_EB,0.041_EB,0.039_EB,0.038_EB/)
ABSF(5,10,11:20,1) = (/0.037_EB,0.035_EB,0.035_EB,0.035_EB,0.034_EB,0.034_EB,0.037_EB,0.039_EB,0.042_EB,0.045_EB/)
ABSF(5,10,21:30,1) = (/0.045_EB,0.046_EB,0.046_EB,0.047_EB,0.047_EB,0.048_EB,0.049_EB,0.050_EB,0.051_EB,0.051_EB/)
ABSF(5,10,0:10,2) = (/0.000_EB,0.024_EB,0.037_EB,0.037_EB,0.038_EB,0.039_EB,0.040_EB,0.042_EB,0.043_EB,0.041_EB,0.039_EB/)
ABSF(5,10,11:20,2) = (/0.036_EB,0.034_EB,0.034_EB,0.034_EB,0.034_EB,0.033_EB,0.037_EB,0.040_EB,0.043_EB,0.046_EB/)
ABSF(5,10,21:30,2) = (/0.047_EB,0.048_EB,0.048_EB,0.049_EB,0.050_EB,0.051_EB,0.052_EB,0.052_EB,0.053_EB,0.054_EB/)
ABSF(5,10,0:10,3) = (/0.000_EB,0.023_EB,0.038_EB,0.044_EB,0.050_EB,0.052_EB,0.054_EB,0.055_EB,0.057_EB,0.055_EB,0.053_EB/)
ABSF(5,10,11:20,3) = (/0.051_EB,0.049_EB,0.048_EB,0.047_EB,0.046_EB,0.045_EB,0.045_EB,0.046_EB,0.046_EB,0.046_EB/)
ABSF(5,10,21:30,3) = (/0.046_EB,0.047_EB,0.047_EB,0.047_EB,0.047_EB,0.045_EB,0.044_EB,0.042_EB,0.040_EB,0.039_EB/)
ABSF(5,10,0:10,4) = (/0.000_EB,0.025_EB,0.043_EB,0.047_EB,0.052_EB,0.055_EB,0.058_EB,0.061_EB,0.064_EB,0.069_EB,0.074_EB/)
ABSF(5,10,11:20,4) = (/0.079_EB,0.084_EB,0.075_EB,0.067_EB,0.058_EB,0.050_EB,0.044_EB,0.037_EB,0.031_EB,0.025_EB/)
ABSF(5,10,21:30,4) = (/0.026_EB,0.028_EB,0.029_EB,0.030_EB,0.032_EB,0.031_EB,0.031_EB,0.030_EB,0.029_EB,0.029_EB/)
ABSF(5,10,0:10,5) = (/0.000_EB,0.026_EB,0.053_EB,0.058_EB,0.063_EB,0.065_EB,0.067_EB,0.070_EB,0.072_EB,0.075_EB,0.077_EB/)
ABSF(5,10,11:20,5) = (/0.080_EB,0.082_EB,0.077_EB,0.072_EB,0.067_EB,0.061_EB,0.045_EB,0.029_EB,0.013_EB,-0.003_EB/)
ABSF(5,10,21:30,5) = (/0.002_EB,0.008_EB,0.013_EB,0.019_EB,0.024_EB,0.024_EB,0.024_EB,0.024_EB,0.025_EB,0.025_EB/)
ABSF(5,10,0:10,6) = (/0.000_EB,0.030_EB,0.053_EB,0.059_EB,0.064_EB,0.061_EB,0.058_EB,0.054_EB,0.051_EB,0.062_EB,0.074_EB/)
ABSF(5,10,11:20,6) = (/0.086_EB,0.097_EB,0.090_EB,0.082_EB,0.075_EB,0.068_EB,0.055_EB,0.042_EB,0.029_EB,0.016_EB/)
ABSF(5,10,21:30,6) = (/0.016_EB,0.017_EB,0.018_EB,0.019_EB,0.020_EB,0.016_EB,0.012_EB,0.008_EB,0.004_EB,0.000_EB/)
ABSF(5,11,0:10,1) = (/0.000_EB,0.033_EB,0.039_EB,0.039_EB,0.040_EB,0.040_EB,0.041_EB,0.042_EB,0.042_EB,0.041_EB,0.039_EB/)
ABSF(5,11,11:20,1) = (/0.038_EB,0.036_EB,0.036_EB,0.036_EB,0.035_EB,0.035_EB,0.038_EB,0.041_EB,0.043_EB,0.046_EB/)
ABSF(5,11,21:30,1) = (/0.046_EB,0.047_EB,0.048_EB,0.048_EB,0.049_EB,0.050_EB,0.051_EB,0.052_EB,0.053_EB,0.053_EB/)
ABSF(5,11,0:10,2) = (/0.000_EB,0.026_EB,0.039_EB,0.039_EB,0.039_EB,0.041_EB,0.042_EB,0.043_EB,0.045_EB,0.042_EB,0.040_EB/)
ABSF(5,11,11:20,2) = (/0.037_EB,0.035_EB,0.035_EB,0.035_EB,0.035_EB,0.035_EB,0.038_EB,0.041_EB,0.044_EB,0.047_EB/)
ABSF(5,11,21:30,2) = (/0.048_EB,0.049_EB,0.050_EB,0.051_EB,0.051_EB,0.052_EB,0.053_EB,0.054_EB,0.055_EB,0.056_EB/)
ABSF(5,11,0:10,3) = (/0.000_EB,0.024_EB,0.040_EB,0.045_EB,0.051_EB,0.052_EB,0.054_EB,0.056_EB,0.058_EB,0.056_EB,0.054_EB/)
ABSF(5,11,11:20,3) = (/0.052_EB,0.050_EB,0.048_EB,0.047_EB,0.046_EB,0.045_EB,0.046_EB,0.046_EB,0.047_EB,0.048_EB/)
ABSF(5,11,21:30,3) = (/0.048_EB,0.048_EB,0.048_EB,0.048_EB,0.049_EB,0.047_EB,0.046_EB,0.045_EB,0.044_EB,0.042_EB/)
ABSF(5,11,0:10,4) = (/0.000_EB,0.026_EB,0.044_EB,0.049_EB,0.053_EB,0.056_EB,0.059_EB,0.063_EB,0.066_EB,0.071_EB,0.075_EB/)
ABSF(5,11,11:20,4) = (/0.080_EB,0.085_EB,0.076_EB,0.068_EB,0.059_EB,0.051_EB,0.045_EB,0.039_EB,0.033_EB,0.027_EB/)
ABSF(5,11,21:30,4) = (/0.029_EB,0.030_EB,0.032_EB,0.033_EB,0.035_EB,0.034_EB,0.033_EB,0.033_EB,0.032_EB,0.032_EB/)
ABSF(5,11,0:10,5) = (/0.000_EB,0.027_EB,0.055_EB,0.059_EB,0.064_EB,0.066_EB,0.069_EB,0.071_EB,0.073_EB,0.076_EB,0.079_EB/)
ABSF(5,11,11:20,5) = (/0.081_EB,0.084_EB,0.078_EB,0.073_EB,0.067_EB,0.062_EB,0.046_EB,0.031_EB,0.015_EB,-0.001_EB/)
ABSF(5,11,21:30,5) = (/0.005_EB,0.010_EB,0.016_EB,0.021_EB,0.027_EB,0.026_EB,0.026_EB,0.025_EB,0.025_EB,0.025_EB/)
ABSF(5,11,0:10,6) = (/0.000_EB,0.031_EB,0.054_EB,0.060_EB,0.066_EB,0.062_EB,0.059_EB,0.056_EB,0.052_EB,0.064_EB,0.075_EB/)
ABSF(5,11,11:20,6) = (/0.086_EB,0.098_EB,0.090_EB,0.082_EB,0.073_EB,0.065_EB,0.053_EB,0.041_EB,0.029_EB,0.016_EB/)
ABSF(5,11,21:30,6) = (/0.017_EB,0.018_EB,0.019_EB,0.020_EB,0.021_EB,0.018_EB,0.014_EB,0.011_EB,0.007_EB,0.004_EB/)
ABSF(5,12,0:10,1) = (/0.000_EB,0.034_EB,0.041_EB,0.041_EB,0.041_EB,0.042_EB,0.043_EB,0.043_EB,0.044_EB,0.042_EB,0.041_EB/)
ABSF(5,12,11:20,1) = (/0.039_EB,0.037_EB,0.037_EB,0.037_EB,0.037_EB,0.036_EB,0.039_EB,0.042_EB,0.045_EB,0.047_EB/)
ABSF(5,12,21:30,1) = (/0.048_EB,0.048_EB,0.049_EB,0.050_EB,0.050_EB,0.051_EB,0.052_EB,0.053_EB,0.054_EB,0.055_EB/)
ABSF(5,12,0:10,2) = (/0.000_EB,0.028_EB,0.040_EB,0.041_EB,0.041_EB,0.042_EB,0.043_EB,0.045_EB,0.046_EB,0.043_EB,0.041_EB/)
ABSF(5,12,11:20,2) = (/0.038_EB,0.036_EB,0.036_EB,0.036_EB,0.036_EB,0.036_EB,0.039_EB,0.042_EB,0.045_EB,0.048_EB/)
ABSF(5,12,21:30,2) = (/0.049_EB,0.050_EB,0.051_EB,0.052_EB,0.053_EB,0.054_EB,0.055_EB,0.056_EB,0.056_EB,0.057_EB/)
ABSF(5,12,0:10,3) = (/0.000_EB,0.025_EB,0.042_EB,0.046_EB,0.051_EB,0.053_EB,0.055_EB,0.056_EB,0.058_EB,0.056_EB,0.054_EB/)
ABSF(5,12,11:20,3) = (/0.052_EB,0.050_EB,0.049_EB,0.048_EB,0.046_EB,0.045_EB,0.046_EB,0.047_EB,0.048_EB,0.049_EB/)
ABSF(5,12,21:30,3) = (/0.049_EB,0.050_EB,0.050_EB,0.050_EB,0.050_EB,0.050_EB,0.049_EB,0.048_EB,0.047_EB,0.046_EB/)
ABSF(5,12,0:10,4) = (/0.000_EB,0.026_EB,0.045_EB,0.050_EB,0.054_EB,0.057_EB,0.061_EB,0.064_EB,0.068_EB,0.072_EB,0.077_EB/)
ABSF(5,12,11:20,4) = (/0.081_EB,0.085_EB,0.077_EB,0.068_EB,0.060_EB,0.051_EB,0.046_EB,0.040_EB,0.035_EB,0.029_EB/)
ABSF(5,12,21:30,4) = (/0.031_EB,0.032_EB,0.034_EB,0.036_EB,0.037_EB,0.037_EB,0.036_EB,0.036_EB,0.035_EB,0.035_EB/)
ABSF(5,12,0:10,5) = (/0.000_EB,0.027_EB,0.056_EB,0.061_EB,0.066_EB,0.068_EB,0.070_EB,0.072_EB,0.075_EB,0.077_EB,0.080_EB/)
ABSF(5,12,11:20,5) = (/0.083_EB,0.086_EB,0.080_EB,0.074_EB,0.068_EB,0.063_EB,0.047_EB,0.032_EB,0.017_EB,0.002_EB/)
ABSF(5,12,21:30,5) = (/0.007_EB,0.013_EB,0.018_EB,0.024_EB,0.029_EB,0.028_EB,0.027_EB,0.026_EB,0.025_EB,0.024_EB/)
ABSF(5,12,0:10,6) = (/0.000_EB,0.032_EB,0.055_EB,0.061_EB,0.067_EB,0.064_EB,0.060_EB,0.057_EB,0.054_EB,0.065_EB,0.076_EB/)
ABSF(5,12,11:20,6) = (/0.087_EB,0.099_EB,0.090_EB,0.081_EB,0.072_EB,0.063_EB,0.052_EB,0.040_EB,0.029_EB,0.017_EB/)
ABSF(5,12,21:30,6) = (/0.019_EB,0.020_EB,0.021_EB,0.022_EB,0.023_EB,0.020_EB,0.017_EB,0.014_EB,0.011_EB,0.007_EB/)
ABSF(5,13,0:10,1) = (/0.000_EB,0.036_EB,0.043_EB,0.043_EB,0.043_EB,0.044_EB,0.044_EB,0.045_EB,0.046_EB,0.044_EB,0.042_EB/)
ABSF(5,13,11:20,1) = (/0.040_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB,0.040_EB,0.043_EB,0.046_EB,0.048_EB/)
ABSF(5,13,21:30,1) = (/0.049_EB,0.050_EB,0.051_EB,0.051_EB,0.052_EB,0.053_EB,0.054_EB,0.055_EB,0.056_EB,0.057_EB/)
ABSF(5,13,0:10,2) = (/0.000_EB,0.030_EB,0.042_EB,0.043_EB,0.043_EB,0.044_EB,0.045_EB,0.046_EB,0.047_EB,0.044_EB,0.042_EB/)
ABSF(5,13,11:20,2) = (/0.040_EB,0.037_EB,0.037_EB,0.037_EB,0.037_EB,0.037_EB,0.040_EB,0.043_EB,0.046_EB,0.049_EB/)
ABSF(5,13,21:30,2) = (/0.050_EB,0.051_EB,0.052_EB,0.053_EB,0.054_EB,0.055_EB,0.056_EB,0.057_EB,0.058_EB,0.059_EB/)
ABSF(5,13,0:10,3) = (/0.000_EB,0.025_EB,0.043_EB,0.047_EB,0.052_EB,0.053_EB,0.055_EB,0.057_EB,0.059_EB,0.057_EB,0.054_EB/)
ABSF(5,13,11:20,3) = (/0.052_EB,0.050_EB,0.049_EB,0.048_EB,0.047_EB,0.045_EB,0.047_EB,0.048_EB,0.049_EB,0.051_EB/)
ABSF(5,13,21:30,3) = (/0.051_EB,0.051_EB,0.051_EB,0.052_EB,0.052_EB,0.052_EB,0.051_EB,0.051_EB,0.050_EB,0.050_EB/)
ABSF(5,13,0:10,4) = (/0.000_EB,0.027_EB,0.047_EB,0.051_EB,0.055_EB,0.059_EB,0.062_EB,0.066_EB,0.069_EB,0.074_EB,0.078_EB/)
ABSF(5,13,11:20,4) = (/0.082_EB,0.086_EB,0.078_EB,0.069_EB,0.060_EB,0.052_EB,0.047_EB,0.042_EB,0.037_EB,0.032_EB/)
ABSF(5,13,21:30,4) = (/0.033_EB,0.035_EB,0.037_EB,0.038_EB,0.040_EB,0.039_EB,0.039_EB,0.039_EB,0.038_EB,0.038_EB/)
ABSF(5,13,0:10,5) = (/0.000_EB,0.028_EB,0.057_EB,0.062_EB,0.067_EB,0.069_EB,0.072_EB,0.074_EB,0.076_EB,0.079_EB,0.082_EB/)
ABSF(5,13,11:20,5) = (/0.085_EB,0.087_EB,0.081_EB,0.075_EB,0.069_EB,0.063_EB,0.048_EB,0.034_EB,0.019_EB,0.004_EB/)
ABSF(5,13,21:30,5) = (/0.010_EB,0.015_EB,0.021_EB,0.026_EB,0.032_EB,0.030_EB,0.029_EB,0.027_EB,0.026_EB,0.024_EB/)
ABSF(5,13,0:10,6) = (/0.000_EB,0.033_EB,0.057_EB,0.062_EB,0.068_EB,0.065_EB,0.062_EB,0.058_EB,0.055_EB,0.066_EB,0.077_EB/)
ABSF(5,13,11:20,6) = (/0.088_EB,0.099_EB,0.090_EB,0.080_EB,0.070_EB,0.061_EB,0.050_EB,0.040_EB,0.029_EB,0.018_EB/)
ABSF(5,13,21:30,6) = (/0.020_EB,0.021_EB,0.022_EB,0.023_EB,0.024_EB,0.022_EB,0.019_EB,0.016_EB,0.014_EB,0.011_EB/)
ABSF(5,14,0:10,1) = (/0.000_EB,0.037_EB,0.045_EB,0.045_EB,0.044_EB,0.045_EB,0.046_EB,0.047_EB,0.048_EB,0.045_EB,0.043_EB/)
ABSF(5,14,11:20,1) = (/0.041_EB,0.039_EB,0.039_EB,0.039_EB,0.039_EB,0.039_EB,0.042_EB,0.044_EB,0.047_EB,0.050_EB/)
ABSF(5,14,21:30,1) = (/0.051_EB,0.051_EB,0.052_EB,0.053_EB,0.053_EB,0.055_EB,0.056_EB,0.057_EB,0.058_EB,0.059_EB/)
ABSF(5,14,0:10,2) = (/0.000_EB,0.032_EB,0.044_EB,0.044_EB,0.045_EB,0.045_EB,0.046_EB,0.047_EB,0.048_EB,0.046_EB,0.043_EB/)
ABSF(5,14,11:20,2) = (/0.041_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB,0.038_EB,0.041_EB,0.044_EB,0.047_EB,0.050_EB/)
ABSF(5,14,21:30,2) = (/0.051_EB,0.052_EB,0.053_EB,0.054_EB,0.055_EB,0.057_EB,0.058_EB,0.059_EB,0.060_EB,0.061_EB/)
ABSF(5,14,0:10,3) = (/0.000_EB,0.026_EB,0.045_EB,0.048_EB,0.052_EB,0.054_EB,0.056_EB,0.057_EB,0.059_EB,0.057_EB,0.055_EB/)
ABSF(5,14,11:20,3) = (/0.053_EB,0.051_EB,0.049_EB,0.048_EB,0.047_EB,0.046_EB,0.047_EB,0.049_EB,0.050_EB,0.052_EB/)
ABSF(5,14,21:30,3) = (/0.052_EB,0.053_EB,0.053_EB,0.053_EB,0.054_EB,0.054_EB,0.054_EB,0.053_EB,0.053_EB,0.053_EB/)
ABSF(5,14,0:10,4) = (/0.000_EB,0.028_EB,0.048_EB,0.052_EB,0.056_EB,0.060_EB,0.064_EB,0.068_EB,0.071_EB,0.075_EB,0.079_EB/)
ABSF(5,14,11:20,4) = (/0.083_EB,0.087_EB,0.078_EB,0.070_EB,0.061_EB,0.052_EB,0.048_EB,0.043_EB,0.038_EB,0.034_EB/)
ABSF(5,14,21:30,4) = (/0.035_EB,0.037_EB,0.039_EB,0.041_EB,0.043_EB,0.042_EB,0.042_EB,0.042_EB,0.041_EB,0.041_EB/)
ABSF(5,14,0:10,5) = (/0.000_EB,0.029_EB,0.058_EB,0.063_EB,0.068_EB,0.071_EB,0.073_EB,0.075_EB,0.077_EB,0.080_EB,0.083_EB/)
ABSF(5,14,11:20,5) = (/0.086_EB,0.089_EB,0.083_EB,0.076_EB,0.070_EB,0.064_EB,0.050_EB,0.035_EB,0.021_EB,0.007_EB/)
ABSF(5,14,21:30,5) = (/0.012_EB,0.018_EB,0.023_EB,0.029_EB,0.034_EB,0.032_EB,0.030_EB,0.028_EB,0.026_EB,0.024_EB/)
ABSF(5,14,0:10,6) = (/0.000_EB,0.034_EB,0.058_EB,0.064_EB,0.070_EB,0.066_EB,0.063_EB,0.059_EB,0.056_EB,0.067_EB,0.078_EB/)
ABSF(5,14,11:20,6) = (/0.089_EB,0.100_EB,0.090_EB,0.079_EB,0.069_EB,0.058_EB,0.049_EB,0.039_EB,0.029_EB,0.019_EB/)
ABSF(5,14,21:30,6) = (/0.021_EB,0.022_EB,0.023_EB,0.025_EB,0.026_EB,0.024_EB,0.021_EB,0.019_EB,0.017_EB,0.015_EB/)
ABSF(5,15,0:10,1) = (/0.000_EB,0.038_EB,0.047_EB,0.046_EB,0.046_EB,0.047_EB,0.048_EB,0.049_EB,0.049_EB,0.047_EB,0.045_EB/)
ABSF(5,15,11:20,1) = (/0.042_EB,0.040_EB,0.040_EB,0.040_EB,0.040_EB,0.040_EB,0.043_EB,0.046_EB,0.048_EB,0.051_EB/)
ABSF(5,15,21:30,1) = (/0.052_EB,0.053_EB,0.053_EB,0.054_EB,0.055_EB,0.056_EB,0.058_EB,0.059_EB,0.060_EB,0.061_EB/)
ABSF(5,15,0:10,2) = (/0.000_EB,0.034_EB,0.046_EB,0.046_EB,0.046_EB,0.047_EB,0.048_EB,0.048_EB,0.049_EB,0.047_EB,0.044_EB/)
ABSF(5,15,11:20,2) = (/0.042_EB,0.039_EB,0.039_EB,0.039_EB,0.039_EB,0.039_EB,0.042_EB,0.045_EB,0.048_EB,0.051_EB/)
ABSF(5,15,21:30,2) = (/0.052_EB,0.053_EB,0.054_EB,0.056_EB,0.057_EB,0.058_EB,0.059_EB,0.060_EB,0.062_EB,0.063_EB/)
ABSF(5,15,0:10,3) = (/0.000_EB,0.027_EB,0.047_EB,0.050_EB,0.052_EB,0.054_EB,0.056_EB,0.058_EB,0.060_EB,0.058_EB,0.055_EB/)
ABSF(5,15,11:20,3) = (/0.053_EB,0.051_EB,0.050_EB,0.048_EB,0.047_EB,0.046_EB,0.048_EB,0.050_EB,0.051_EB,0.053_EB/)
ABSF(5,15,21:30,3) = (/0.054_EB,0.054_EB,0.055_EB,0.055_EB,0.056_EB,0.056_EB,0.056_EB,0.056_EB,0.057_EB,0.057_EB/)
ABSF(5,15,0:10,4) = (/0.000_EB,0.029_EB,0.049_EB,0.053_EB,0.058_EB,0.062_EB,0.065_EB,0.069_EB,0.073_EB,0.077_EB,0.080_EB/)
ABSF(5,15,11:20,4) = (/0.084_EB,0.088_EB,0.079_EB,0.070_EB,0.062_EB,0.053_EB,0.049_EB,0.044_EB,0.040_EB,0.036_EB/)
ABSF(5,15,21:30,4) = (/0.038_EB,0.040_EB,0.041_EB,0.043_EB,0.045_EB,0.045_EB,0.045_EB,0.045_EB,0.045_EB,0.044_EB/)
ABSF(5,15,0:10,5) = (/0.000_EB,0.030_EB,0.060_EB,0.065_EB,0.070_EB,0.072_EB,0.074_EB,0.076_EB,0.079_EB,0.082_EB,0.085_EB/)
ABSF(5,15,11:20,5) = (/0.088_EB,0.091_EB,0.084_EB,0.078_EB,0.071_EB,0.064_EB,0.051_EB,0.037_EB,0.023_EB,0.009_EB/)
ABSF(5,15,21:30,5) = (/0.015_EB,0.020_EB,0.026_EB,0.031_EB,0.037_EB,0.034_EB,0.032_EB,0.029_EB,0.027_EB,0.024_EB/)
ABSF(5,15,0:10,6) = (/0.000_EB,0.035_EB,0.059_EB,0.065_EB,0.071_EB,0.068_EB,0.064_EB,0.061_EB,0.057_EB,0.068_EB,0.079_EB/)
ABSF(5,15,11:20,6) = (/0.090_EB,0.101_EB,0.090_EB,0.078_EB,0.067_EB,0.056_EB,0.047_EB,0.038_EB,0.029_EB,0.020_EB/)
ABSF(5,15,21:30,6) = (/0.022_EB,0.023_EB,0.025_EB,0.026_EB,0.027_EB,0.026_EB,0.024_EB,0.022_EB,0.020_EB,0.018_EB/)
ABSF(5,16,0:10,1) = (/0.000_EB,0.039_EB,0.049_EB,0.048_EB,0.048_EB,0.049_EB,0.049_EB,0.050_EB,0.051_EB,0.049_EB,0.046_EB/)
ABSF(5,16,11:20,1) = (/0.043_EB,0.041_EB,0.041_EB,0.041_EB,0.041_EB,0.042_EB,0.044_EB,0.047_EB,0.050_EB,0.052_EB/)
ABSF(5,16,21:30,1) = (/0.053_EB,0.054_EB,0.055_EB,0.056_EB,0.057_EB,0.058_EB,0.059_EB,0.061_EB,0.062_EB,0.063_EB/)
ABSF(5,16,0:10,2) = (/0.000_EB,0.036_EB,0.048_EB,0.048_EB,0.048_EB,0.049_EB,0.049_EB,0.050_EB,0.050_EB,0.048_EB,0.045_EB/)
ABSF(5,16,11:20,2) = (/0.043_EB,0.041_EB,0.041_EB,0.041_EB,0.041_EB,0.041_EB,0.043_EB,0.046_EB,0.049_EB,0.052_EB/)
ABSF(5,16,21:30,2) = (/0.053_EB,0.054_EB,0.056_EB,0.057_EB,0.058_EB,0.059_EB,0.061_EB,0.062_EB,0.063_EB,0.065_EB/)
ABSF(5,16,0:10,3) = (/0.000_EB,0.028_EB,0.048_EB,0.051_EB,0.053_EB,0.055_EB,0.057_EB,0.058_EB,0.060_EB,0.058_EB,0.056_EB/)
ABSF(5,16,11:20,3) = (/0.053_EB,0.051_EB,0.050_EB,0.049_EB,0.047_EB,0.046_EB,0.048_EB,0.050_EB,0.053_EB,0.055_EB/)
ABSF(5,16,21:30,3) = (/0.055_EB,0.056_EB,0.056_EB,0.057_EB,0.057_EB,0.058_EB,0.059_EB,0.059_EB,0.060_EB,0.060_EB/)
ABSF(5,16,0:10,4) = (/0.000_EB,0.029_EB,0.050_EB,0.054_EB,0.059_EB,0.063_EB,0.067_EB,0.071_EB,0.075_EB,0.078_EB,0.082_EB/)
ABSF(5,16,11:20,4) = (/0.085_EB,0.089_EB,0.080_EB,0.071_EB,0.062_EB,0.053_EB,0.050_EB,0.046_EB,0.042_EB,0.038_EB/)
ABSF(5,16,21:30,4) = (/0.040_EB,0.042_EB,0.044_EB,0.046_EB,0.048_EB,0.048_EB,0.048_EB,0.048_EB,0.048_EB,0.048_EB/)
ABSF(5,16,0:10,5) = (/0.000_EB,0.030_EB,0.061_EB,0.066_EB,0.071_EB,0.073_EB,0.076_EB,0.078_EB,0.080_EB,0.083_EB,0.086_EB/)
ABSF(5,16,11:20,5) = (/0.089_EB,0.092_EB,0.086_EB,0.079_EB,0.072_EB,0.065_EB,0.052_EB,0.038_EB,0.025_EB,0.012_EB/)
ABSF(5,16,21:30,5) = (/0.017_EB,0.023_EB,0.028_EB,0.034_EB,0.039_EB,0.036_EB,0.033_EB,0.030_EB,0.027_EB,0.024_EB/)
ABSF(5,16,0:10,6) = (/0.000_EB,0.036_EB,0.060_EB,0.066_EB,0.072_EB,0.069_EB,0.065_EB,0.062_EB,0.059_EB,0.069_EB,0.080_EB/)
ABSF(5,16,11:20,6) = (/0.091_EB,0.102_EB,0.090_EB,0.078_EB,0.066_EB,0.053_EB,0.045_EB,0.037_EB,0.029_EB,0.021_EB/)
ABSF(5,16,21:30,6) = (/0.023_EB,0.024_EB,0.026_EB,0.027_EB,0.029_EB,0.028_EB,0.026_EB,0.025_EB,0.023_EB,0.022_EB/)
ABSF(5,17,0:10,1) = (/0.000_EB,0.040_EB,0.051_EB,0.050_EB,0.049_EB,0.050_EB,0.051_EB,0.052_EB,0.053_EB,0.050_EB,0.047_EB/)
ABSF(5,17,11:20,1) = (/0.044_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB,0.043_EB,0.045_EB,0.048_EB,0.051_EB,0.054_EB/)
ABSF(5,17,21:30,1) = (/0.055_EB,0.055_EB,0.056_EB,0.057_EB,0.058_EB,0.060_EB,0.061_EB,0.062_EB,0.064_EB,0.065_EB/)
ABSF(5,17,0:10,2) = (/0.000_EB,0.038_EB,0.049_EB,0.050_EB,0.050_EB,0.050_EB,0.051_EB,0.051_EB,0.051_EB,0.049_EB,0.047_EB/)
ABSF(5,17,11:20,2) = (/0.044_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB,0.042_EB,0.045_EB,0.047_EB,0.050_EB,0.053_EB/)
ABSF(5,17,21:30,2) = (/0.054_EB,0.056_EB,0.057_EB,0.058_EB,0.060_EB,0.061_EB,0.062_EB,0.064_EB,0.065_EB,0.066_EB/)
ABSF(5,17,0:10,3) = (/0.000_EB,0.028_EB,0.050_EB,0.052_EB,0.053_EB,0.055_EB,0.057_EB,0.059_EB,0.061_EB,0.058_EB,0.056_EB/)
ABSF(5,17,11:20,3) = (/0.054_EB,0.051_EB,0.050_EB,0.049_EB,0.048_EB,0.046_EB,0.049_EB,0.051_EB,0.054_EB,0.056_EB/)
ABSF(5,17,21:30,3) = (/0.057_EB,0.057_EB,0.058_EB,0.059_EB,0.059_EB,0.060_EB,0.061_EB,0.062_EB,0.063_EB,0.064_EB/)
ABSF(5,17,0:10,4) = (/0.000_EB,0.030_EB,0.051_EB,0.056_EB,0.060_EB,0.064_EB,0.068_EB,0.073_EB,0.077_EB,0.080_EB,0.083_EB/)
ABSF(5,17,11:20,4) = (/0.086_EB,0.089_EB,0.080_EB,0.072_EB,0.063_EB,0.054_EB,0.050_EB,0.047_EB,0.044_EB,0.040_EB/)
ABSF(5,17,21:30,4) = (/0.042_EB,0.044_EB,0.046_EB,0.048_EB,0.050_EB,0.051_EB,0.051_EB,0.051_EB,0.051_EB,0.051_EB/)
ABSF(5,17,0:10,5) = (/0.000_EB,0.031_EB,0.062_EB,0.067_EB,0.073_EB,0.075_EB,0.077_EB,0.079_EB,0.081_EB,0.085_EB,0.088_EB/)
ABSF(5,17,11:20,5) = (/0.091_EB,0.094_EB,0.087_EB,0.080_EB,0.073_EB,0.066_EB,0.053_EB,0.040_EB,0.027_EB,0.014_EB/)
ABSF(5,17,21:30,5) = (/0.020_EB,0.025_EB,0.031_EB,0.036_EB,0.042_EB,0.038_EB,0.035_EB,0.031_EB,0.027_EB,0.024_EB/)
ABSF(5,17,0:10,6) = (/0.000_EB,0.037_EB,0.062_EB,0.068_EB,0.074_EB,0.070_EB,0.067_EB,0.063_EB,0.060_EB,0.071_EB,0.081_EB/)
ABSF(5,17,11:20,6) = (/0.092_EB,0.103_EB,0.090_EB,0.077_EB,0.064_EB,0.051_EB,0.044_EB,0.037_EB,0.030_EB,0.022_EB/)
ABSF(5,17,21:30,6) = (/0.024_EB,0.026_EB,0.027_EB,0.029_EB,0.031_EB,0.030_EB,0.029_EB,0.028_EB,0.027_EB,0.026_EB/)
ABSF(5,18,0:10,1) = (/0.000_EB,0.041_EB,0.052_EB,0.052_EB,0.051_EB,0.052_EB,0.053_EB,0.054_EB,0.055_EB,0.052_EB,0.049_EB/)
ABSF(5,18,11:20,1) = (/0.046_EB,0.043_EB,0.043_EB,0.043_EB,0.044_EB,0.044_EB,0.047_EB,0.049_EB,0.052_EB,0.055_EB/)
ABSF(5,18,21:30,1) = (/0.056_EB,0.057_EB,0.058_EB,0.059_EB,0.060_EB,0.061_EB,0.063_EB,0.064_EB,0.066_EB,0.067_EB/)
ABSF(5,18,0:10,2) = (/0.000_EB,0.040_EB,0.051_EB,0.051_EB,0.052_EB,0.052_EB,0.052_EB,0.052_EB,0.053_EB,0.050_EB,0.048_EB/)
ABSF(5,18,11:20,2) = (/0.045_EB,0.043_EB,0.043_EB,0.043_EB,0.043_EB,0.043_EB,0.046_EB,0.048_EB,0.051_EB,0.054_EB/)
ABSF(5,18,21:30,2) = (/0.055_EB,0.057_EB,0.058_EB,0.059_EB,0.061_EB,0.062_EB,0.064_EB,0.065_EB,0.067_EB,0.068_EB/)
ABSF(5,18,0:10,3) = (/0.000_EB,0.029_EB,0.051_EB,0.053_EB,0.054_EB,0.056_EB,0.058_EB,0.059_EB,0.061_EB,0.059_EB,0.057_EB/)
ABSF(5,18,11:20,3) = (/0.054_EB,0.052_EB,0.050_EB,0.049_EB,0.048_EB,0.047_EB,0.049_EB,0.052_EB,0.055_EB,0.057_EB/)
ABSF(5,18,21:30,3) = (/0.058_EB,0.059_EB,0.060_EB,0.060_EB,0.061_EB,0.062_EB,0.064_EB,0.065_EB,0.066_EB,0.067_EB/)
ABSF(5,18,0:10,4) = (/0.000_EB,0.031_EB,0.052_EB,0.057_EB,0.061_EB,0.066_EB,0.070_EB,0.074_EB,0.078_EB,0.081_EB,0.084_EB/)
ABSF(5,18,11:20,4) = (/0.087_EB,0.090_EB,0.081_EB,0.072_EB,0.063_EB,0.054_EB,0.051_EB,0.049_EB,0.046_EB,0.043_EB/)
ABSF(5,18,21:30,4) = (/0.045_EB,0.047_EB,0.049_EB,0.051_EB,0.053_EB,0.053_EB,0.053_EB,0.054_EB,0.054_EB,0.054_EB/)
ABSF(5,18,0:10,5) = (/0.000_EB,0.032_EB,0.064_EB,0.069_EB,0.074_EB,0.076_EB,0.078_EB,0.080_EB,0.083_EB,0.086_EB,0.089_EB/)
ABSF(5,18,11:20,5) = (/0.093_EB,0.096_EB,0.088_EB,0.081_EB,0.074_EB,0.066_EB,0.054_EB,0.041_EB,0.029_EB,0.016_EB/)
ABSF(5,18,21:30,5) = (/0.022_EB,0.028_EB,0.033_EB,0.039_EB,0.044_EB,0.040_EB,0.036_EB,0.032_EB,0.028_EB,0.024_EB/)
ABSF(5,18,0:10,6) = (/0.000_EB,0.038_EB,0.063_EB,0.069_EB,0.075_EB,0.071_EB,0.068_EB,0.065_EB,0.061_EB,0.072_EB,0.082_EB/)
ABSF(5,18,11:20,6) = (/0.093_EB,0.104_EB,0.090_EB,0.076_EB,0.062_EB,0.049_EB,0.042_EB,0.036_EB,0.030_EB,0.023_EB/)
ABSF(5,18,21:30,6) = (/0.025_EB,0.027_EB,0.029_EB,0.030_EB,0.032_EB,0.032_EB,0.031_EB,0.030_EB,0.030_EB,0.029_EB/)
ABSF(5,19,0:10,1) = (/0.000_EB,0.042_EB,0.054_EB,0.053_EB,0.052_EB,0.053_EB,0.054_EB,0.055_EB,0.056_EB,0.053_EB,0.050_EB/)
ABSF(5,19,11:20,1) = (/0.047_EB,0.043_EB,0.044_EB,0.044_EB,0.045_EB,0.045_EB,0.048_EB,0.051_EB,0.053_EB,0.056_EB/)
ABSF(5,19,21:30,1) = (/0.057_EB,0.058_EB,0.059_EB,0.060_EB,0.061_EB,0.063_EB,0.064_EB,0.066_EB,0.068_EB,0.069_EB/)
ABSF(5,19,0:10,2) = (/0.000_EB,0.042_EB,0.053_EB,0.053_EB,0.053_EB,0.054_EB,0.054_EB,0.054_EB,0.054_EB,0.051_EB,0.049_EB/)
ABSF(5,19,11:20,2) = (/0.046_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.044_EB,0.047_EB,0.049_EB,0.052_EB,0.055_EB/)
ABSF(5,19,21:30,2) = (/0.056_EB,0.058_EB,0.059_EB,0.061_EB,0.062_EB,0.064_EB,0.065_EB,0.067_EB,0.068_EB,0.070_EB/)
ABSF(5,19,0:10,3) = (/0.000_EB,0.030_EB,0.053_EB,0.054_EB,0.054_EB,0.056_EB,0.058_EB,0.060_EB,0.062_EB,0.059_EB,0.057_EB/)
ABSF(5,19,11:20,3) = (/0.055_EB,0.052_EB,0.051_EB,0.049_EB,0.048_EB,0.047_EB,0.050_EB,0.053_EB,0.056_EB,0.059_EB/)
ABSF(5,19,21:30,3) = (/0.060_EB,0.060_EB,0.061_EB,0.062_EB,0.063_EB,0.064_EB,0.066_EB,0.068_EB,0.069_EB,0.071_EB/)
ABSF(5,19,0:10,4) = (/0.000_EB,0.031_EB,0.053_EB,0.058_EB,0.062_EB,0.067_EB,0.071_EB,0.076_EB,0.080_EB,0.083_EB,0.086_EB/)
ABSF(5,19,11:20,4) = (/0.088_EB,0.091_EB,0.082_EB,0.073_EB,0.064_EB,0.055_EB,0.052_EB,0.050_EB,0.047_EB,0.045_EB/)
ABSF(5,19,21:30,4) = (/0.047_EB,0.049_EB,0.051_EB,0.054_EB,0.056_EB,0.056_EB,0.056_EB,0.057_EB,0.057_EB,0.057_EB/)
ABSF(5,19,0:10,5) = (/0.000_EB,0.032_EB,0.065_EB,0.070_EB,0.075_EB,0.077_EB,0.080_EB,0.082_EB,0.084_EB,0.087_EB,0.091_EB/)
ABSF(5,19,11:20,5) = (/0.094_EB,0.098_EB,0.090_EB,0.082_EB,0.075_EB,0.067_EB,0.055_EB,0.043_EB,0.031_EB,0.019_EB/)
ABSF(5,19,21:30,5) = (/0.024_EB,0.030_EB,0.036_EB,0.041_EB,0.047_EB,0.042_EB,0.038_EB,0.033_EB,0.028_EB,0.023_EB/)
ABSF(5,19,0:10,6) = (/0.000_EB,0.039_EB,0.064_EB,0.070_EB,0.076_EB,0.073_EB,0.069_EB,0.066_EB,0.062_EB,0.073_EB,0.083_EB/)
ABSF(5,19,11:20,6) = (/0.094_EB,0.104_EB,0.090_EB,0.075_EB,0.061_EB,0.046_EB,0.041_EB,0.035_EB,0.030_EB,0.024_EB/)
ABSF(5,19,21:30,6) = (/0.026_EB,0.028_EB,0.030_EB,0.032_EB,0.034_EB,0.034_EB,0.033_EB,0.033_EB,0.033_EB,0.033_EB/)
ABSF(5,20,0:10,1) = (/0.000_EB,0.044_EB,0.056_EB,0.055_EB,0.054_EB,0.055_EB,0.056_EB,0.057_EB,0.058_EB,0.055_EB,0.051_EB/)
ABSF(5,20,11:20,1) = (/0.048_EB,0.044_EB,0.045_EB,0.045_EB,0.046_EB,0.047_EB,0.049_EB,0.052_EB,0.055_EB,0.058_EB/)
ABSF(5,20,21:30,1) = (/0.059_EB,0.060_EB,0.061_EB,0.062_EB,0.063_EB,0.064_EB,0.066_EB,0.068_EB,0.069_EB,0.071_EB/)
ABSF(5,20,0:10,2) = (/0.000_EB,0.044_EB,0.055_EB,0.055_EB,0.055_EB,0.055_EB,0.055_EB,0.055_EB,0.055_EB,0.052_EB,0.050_EB/)
ABSF(5,20,11:20,2) = (/0.047_EB,0.045_EB,0.045_EB,0.045_EB,0.045_EB,0.045_EB,0.048_EB,0.051_EB,0.053_EB,0.056_EB/)
ABSF(5,20,21:30,2) = (/0.057_EB,0.059_EB,0.060_EB,0.062_EB,0.064_EB,0.065_EB,0.067_EB,0.068_EB,0.070_EB,0.072_EB/)
ABSF(5,20,0:10,3) = (/0.000_EB,0.031_EB,0.055_EB,0.055_EB,0.055_EB,0.057_EB,0.058_EB,0.060_EB,0.062_EB,0.060_EB,0.057_EB/)
ABSF(5,20,11:20,3) = (/0.055_EB,0.052_EB,0.051_EB,0.050_EB,0.048_EB,0.047_EB,0.050_EB,0.054_EB,0.057_EB,0.060_EB/)
ABSF(5,20,21:30,3) = (/0.061_EB,0.062_EB,0.063_EB,0.064_EB,0.064_EB,0.066_EB,0.068_EB,0.071_EB,0.073_EB,0.075_EB/)
ABSF(5,20,0:10,4) = (/0.000_EB,0.032_EB,0.054_EB,0.059_EB,0.064_EB,0.068_EB,0.073_EB,0.077_EB,0.082_EB,0.084_EB,0.087_EB/)
ABSF(5,20,11:20,4) = (/0.089_EB,0.092_EB,0.083_EB,0.074_EB,0.065_EB,0.055_EB,0.053_EB,0.051_EB,0.049_EB,0.047_EB/)
ABSF(5,20,21:30,4) = (/0.049_EB,0.052_EB,0.054_EB,0.056_EB,0.058_EB,0.059_EB,0.059_EB,0.060_EB,0.060_EB,0.060_EB/)
ABSF(5,20,0:10,5) = (/0.000_EB,0.033_EB,0.066_EB,0.071_EB,0.077_EB,0.079_EB,0.081_EB,0.083_EB,0.085_EB,0.089_EB,0.092_EB/)
ABSF(5,20,11:20,5) = (/0.096_EB,0.099_EB,0.091_EB,0.083_EB,0.075_EB,0.067_EB,0.056_EB,0.044_EB,0.033_EB,0.021_EB/)
ABSF(5,20,21:30,5) = (/0.027_EB,0.033_EB,0.038_EB,0.044_EB,0.049_EB,0.044_EB,0.039_EB,0.034_EB,0.029_EB,0.023_EB/)
ABSF(5,20,0:10,6) = (/0.000_EB,0.040_EB,0.065_EB,0.071_EB,0.077_EB,0.074_EB,0.071_EB,0.067_EB,0.064_EB,0.074_EB,0.084_EB/)
ABSF(5,20,11:20,6) = (/0.095_EB,0.105_EB,0.090_EB,0.075_EB,0.059_EB,0.044_EB,0.039_EB,0.035_EB,0.030_EB,0.025_EB/)
ABSF(5,20,21:30,6) = (/0.027_EB,0.029_EB,0.031_EB,0.033_EB,0.035_EB,0.035_EB,0.036_EB,0.036_EB,0.036_EB,0.037_EB/)

END SUBROUTINE DEFINE_QREF_ARRAYS


END MODULE PHYSICAL_FUNCTIONS


!> \brief A collection of routines that manage memory

MODULE MEMORY_FUNCTIONS

USE PRECISION_PARAMETERS
USE MESH_VARIABLES
USE COMP_FUNCTIONS, ONLY : SHUTDOWN
IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS

!> \brief Constructs an error message if there was an error during allocatiion of an array
!> \param CodeSect Character string containing the subroutine or function the allocation statement is in
!> \param VarName Character string containing the name of the array being allocated
!> \param IZERO Error value returned by the ALLOCATE statment where 0 is no error

SUBROUTINE ChkMemErr(CodeSect,VarName,IZERO)

! Memory checking routine

CHARACTER(*), INTENT(IN) :: CodeSect, VarName
INTEGER IZERO
CHARACTER(MESSAGE_LENGTH) MESSAGE

IF (IZERO==0) RETURN

WRITE(MESSAGE,'(4A)') 'ERROR: Memory allocation failed for ', TRIM(VarName),' in the routine ',TRIM(CodeSect)
CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)

END SUBROUTINE ChkMemErr


!> \brief Changes the allocation of an array with DIMENSION 1
!> \param P Original array
!> \param N1 Lower bound of new allocation
!> \param N2 Upper bound of new allocation

FUNCTION REALLOCATE(P,N1,N2)

! Resize the array P

REAL(EB), POINTER, DIMENSION(:) :: P, REALLOCATE
INTEGER, INTENT(IN) :: N1,N2
INTEGER :: NOLD, IERR
CHARACTER(MESSAGE_LENGTH) :: MESSAGE
ALLOCATE(REALLOCATE(N1:N2), STAT=IERR)
IF (IERR /= 0) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: Memory allocation failed in REALLOCATE'
   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
ENDIF
IF (.NOT. ASSOCIATED(P)) RETURN
NOLD = MIN(SIZE(P), N2-N1+1)
REALLOCATE(N1:NOLD+N1-1) = P(N1:NOLD+N1-1)  ! Restore the contents of the reallocated array
DEALLOCATE(P)

END FUNCTION REALLOCATE


!> \brief Change the allocation of a real array with DIMENSION 1
!> \param P Original array
!> \param N1 Lower bound of old/new allocation
!> \param N2 Upper bound of old allocation
!> \param N3 Upper bound of new allocation

SUBROUTINE REALLOCATE_REAL_ARRAY(P,N1,N2,N3)

REAL(EB), ALLOCATABLE, DIMENSION(:) :: P,DUMMY
INTEGER, INTENT(IN) :: N1,N2,N3

ALLOCATE(DUMMY(N1:N3))
DUMMY(N1:N2) = P(N1:N2)
CALL MOVE_ALLOC(DUMMY,P)

END SUBROUTINE REALLOCATE_REAL_ARRAY


!> \brief Change the allocation of an integer array with DIMENSION 1
!> \param P Original array
!> \param N1 Lower bound of old/new allocation
!> \param N2 Upper bound of old allocation
!> \param N3 Upper bound of new allocation

SUBROUTINE REALLOCATE_INTEGER_ARRAY(P,N1,N2,N3)

INTEGER, ALLOCATABLE, DIMENSION(:) :: P,DUMMY
INTEGER, INTENT(IN) :: N1,N2,N3

ALLOCATE(DUMMY(N1:N3))
IF (ALLOCATED(P)) DUMMY(N1:N2) = P(N1:N2)
CALL MOVE_ALLOC(DUMMY,P)

END SUBROUTINE REALLOCATE_INTEGER_ARRAY


!> \brief Changes the allocation of a string array
!> \param P Original array
!> \param CLEN Length of string
!> \param N1 Lower bound of new allocation
!> \param N2 Upper bound of new allocation

FUNCTION REALLOCATE_CHARACTER_ARRAY(P,CLEN,N1,N2)

! Resize the character array P

INTEGER, INTENT(IN) :: N1,N2,CLEN
CHARACTER(CLEN), POINTER, DIMENSION(:) :: P, REALLOCATE_CHARACTER_ARRAY
INTEGER :: NOLD, IERR
CHARACTER(MESSAGE_LENGTH) :: MESSAGE
ALLOCATE(REALLOCATE_CHARACTER_ARRAY(N1:N2), STAT=IERR)
IF (IERR /= 0) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: Memory allocation failed in REALLOCATE_CHARACTER_ARRAY'
   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
ENDIF
IF (.NOT. ASSOCIATED(P)) RETURN
NOLD = MIN(SIZE(P), N2-N1+1)
REALLOCATE_CHARACTER_ARRAY(N1:NOLD+N1-1) = P(N1:NOLD+N1-1)  ! Restore the contents of the reallocated array
DEALLOCATE(P)

END FUNCTION REALLOCATE_CHARACTER_ARRAY


!> \brief Changes the allocation of an array with DIMENSION 2
!> \param P Original array
!> \param M1 Lower bound of first dimension of new allocation
!> \param M2 Upper bound of first dimension of new allocation
!> \param N1 Lower bound of second dimension of new allocation
!> \param N2 Upper bound of second dimension of new allocation

FUNCTION REALLOCATE2D(P,M1,M2,N1,N2)

REAL(EB), POINTER, DIMENSION(:,:) :: P, REALLOCATE2D
INTEGER, INTENT(IN) :: M1,M2,N1,N2
INTEGER :: MOLD,NOLD,IERR
CHARACTER(MESSAGE_LENGTH) :: MESSAGE
ALLOCATE(REALLOCATE2D(M1:M2,N1:N2), STAT=IERR)
IF (IERR /= 0) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: Memory allocation failed in REALLOCATE2D'
   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
ENDIF
IF (.NOT. ASSOCIATED(P)) RETURN
MOLD = MIN(SIZE(P,DIM=1), M2-M1+1)
NOLD = MIN(SIZE(P,DIM=2), N2-N1+1)
REALLOCATE2D(M1:MOLD+M1-1,N1:NOLD+N1-1) = P(M1:MOLD+M1-1,N1:NOLD+N1-1)  ! Restore the contents of the reallocated array
DEALLOCATE(P)

END FUNCTION REALLOCATE2D


!> \brief Changes the allocation of the array STRING by adding 100 new entries
!> \param NM Mesh number

SUBROUTINE RE_ALLOCATE_STRINGS(NM)

CHARACTER(MESH_STRING_LENGTH), ALLOCATABLE, DIMENSION(:) :: DUMMY
INTEGER, INTENT(IN) :: NM
TYPE (MESH_TYPE), POINTER :: M=>NULL()

M=>MESHES(NM)
ALLOCATE(DUMMY(1:M%N_STRINGS))
DUMMY = M%STRING
DEALLOCATE(M%STRING)
ALLOCATE(M%STRING(M%N_STRINGS_MAX+100))
M%STRING(1:M%N_STRINGS) = DUMMY(1:M%N_STRINGS)
M%N_STRINGS_MAX = M%N_STRINGS_MAX+100
DEALLOCATE(DUMMY)

END SUBROUTINE RE_ALLOCATE_STRINGS


!> \brief Re-allocate the derived type array P_ZONE

SUBROUTINE REALLOCATE_P_ZONE(N1,N2)

INTEGER, INTENT(IN) :: N1,N2
TYPE (P_ZONE_TYPE), DIMENSION(:), ALLOCATABLE :: P_ZONE_DUMMY

ALLOCATE(P_ZONE_DUMMY(1:N2))
IF (ALLOCATED(P_ZONE)) THEN
   P_ZONE_DUMMY(1:N1) = P_ZONE(1:N1)
   DEALLOCATE(P_ZONE)
ENDIF
ALLOCATE(P_ZONE(1:N2))
P_ZONE(1:N2) = P_ZONE_DUMMY(1:N2)
DEALLOCATE(P_ZONE_DUMMY)

END SUBROUTINE REALLOCATE_P_ZONE


!> \brief Allocate or reallocate the derived type array CELL
!> \param NM Mesh index
!> \param N1 Current upper bound of the array CELL
!> \param N2 New upper bound of the array CELL

SUBROUTINE REALLOCATE_CELL(NM,N1,N2)

USE GLOBAL_CONSTANTS, ONLY: CELL_COUNT,CELL_COUNT_INTEGERS,CELL_COUNT_LOGICALS
TYPE(CELL_TYPE), ALLOCATABLE, DIMENSION(:) :: CELL_DUMMY
INTEGER, ALLOCATABLE, DIMENSION(:) :: INTEGER_DUMMY
LOGICAL, ALLOCATABLE, DIMENSION(:) :: LOGICAL_DUMMY
TYPE(MESH_TYPE), POINTER :: M
INTEGER, INTENT(IN) :: NM,N1,N2
INTEGER :: IC,LC,CELL_COUNT_INTEGERS_OLD,CELL_COUNT_LOGICALS_OLD

M => MESHES(NM)

IF (.NOT.ALLOCATED(M%CELL)) THEN
   ALLOCATE(M%CELL(0:N1))
   CELL_COUNT(NM) = N1
   CALL PACK_CELL(NM,.TRUE.,N_INTEGERS=IC,N_LOGICALS=LC)
   CELL_COUNT_INTEGERS(NM) = IC
   CELL_COUNT_LOGICALS(NM) = LC
   ALLOCATE(M%CELL_INTEGERS(1:IC))
   ALLOCATE(M%CELL_LOGICALS(1:LC))
ENDIF

IF (.NOT.ALLOCATED(M%CELL_ILW)) ALLOCATE(M%CELL_ILW(1:N1,1:3))

IF (N1==N2) RETURN

ALLOCATE(CELL_DUMMY(0:N2))
CELL_DUMMY(0:N1) = M%CELL(0:N1)
CALL MOVE_ALLOC(CELL_DUMMY,M%CELL)
CELL_COUNT(NM) = N2

CELL_COUNT_INTEGERS_OLD = SIZE(M%CELL_INTEGERS)
CELL_COUNT_LOGICALS_OLD = SIZE(M%CELL_LOGICALS)

CALL PACK_CELL(NM,.TRUE.,N_INTEGERS=IC,N_LOGICALS=LC)

ALLOCATE(INTEGER_DUMMY(1:IC))
INTEGER_DUMMY(1:CELL_COUNT_INTEGERS_OLD) = M%CELL_INTEGERS(1:CELL_COUNT_INTEGERS_OLD)
CALL MOVE_ALLOC(INTEGER_DUMMY,M%CELL_INTEGERS)
CELL_COUNT_INTEGERS(NM) = IC

ALLOCATE(LOGICAL_DUMMY(1:LC))
LOGICAL_DUMMY(1:CELL_COUNT_LOGICALS_OLD) = M%CELL_LOGICALS(1:CELL_COUNT_LOGICALS_OLD)
CALL MOVE_ALLOC(LOGICAL_DUMMY,M%CELL_LOGICALS)
CELL_COUNT_LOGICALS(NM) = LC

DEALLOCATE(M%CELL_ILW)
ALLOCATE(M%CELL_ILW(1:N2,1:3))

END SUBROUTINE REALLOCATE_CELL


!> \brief Allocate or reallocate the derived type array EDGE
!> \param NM Mesh index
!> \param N1 Current upper bound of the array EDGE
!> \param N2 New upper bound of the array EDGE

SUBROUTINE REALLOCATE_EDGE(NM,N1,N2)

TYPE(EDGE_TYPE), ALLOCATABLE, DIMENSION(:) :: EDGE_DUMMY
TYPE(MESH_TYPE), POINTER :: M
INTEGER, INTENT(IN) :: NM,N1,N2

M => MESHES(NM)

IF (.NOT.ALLOCATED(M%EDGE)) ALLOCATE(M%EDGE(0:N1))

IF (N1==N2) RETURN

ALLOCATE(EDGE_DUMMY(0:N2))
EDGE_DUMMY(0:N1) = M%EDGE(0:N1)
CALL MOVE_ALLOC(EDGE_DUMMY,M%EDGE)

END SUBROUTINE REALLOCATE_EDGE


!> \brief Create the appropriate storage space for derived types WALL, CFACE, or LAGRANGIAN_PARTICLE
!> \param NM Index to the current mesh
!> \param LP_INDEX Index within the LAGRANGIAN_PARTICLE array
!> \param LPC_INDEX Index within the LAGRANGIAN_PARTICLE_CLASS array
!> \param SURF_INDEX Index within the SURFACE array
!> \param WALL_INDEX Index within the WALL array
!> \param CFACE_INDEX Index within the CFACE array
!> \param THIN_WALL_INDEX Index within the THIN_WALL array

SUBROUTINE ALLOCATE_STORAGE(NM,LP_INDEX,LPC_INDEX,SURF_INDEX,WALL_INDEX,CFACE_INDEX,THIN_WALL_INDEX)

USE GLOBAL_CONSTANTS, ONLY: PYROLYSIS_PREDICTED
INTEGER, INTENT(IN):: NM,SURF_INDEX
INTEGER, INTENT(IN), OPTIONAL :: LP_INDEX,LPC_INDEX,WALL_INDEX,CFACE_INDEX,THIN_WALL_INDEX
INTEGER :: N_NEW_STORAGE_SLOTS,OD_INDEX,TD_INDEX,B1_INDEX,B2_INDEX,BR_INDEX,BC_INDEX,I,NAS
INTEGER, ALLOCATABLE, DIMENSION(:) :: DUMMY
LOGICAL :: ALREADY_ALLOCATED
TYPE (LAGRANGIAN_PARTICLE_TYPE), ALLOCATABLE, DIMENSION(:) :: LP_DUMMY
TYPE (WALL_TYPE), ALLOCATABLE, DIMENSION(:) :: WALL_DUMMY
TYPE (THIN_WALL_TYPE), ALLOCATABLE, DIMENSION(:) :: THIN_WALL_DUMMY
TYPE (CFACE_TYPE), ALLOCATABLE, DIMENSION(:) :: CFACE_DUMMY
TYPE (BOUNDARY_COORD_TYPE), ALLOCATABLE, DIMENSION(:) :: BC_DUMMY
TYPE (BOUNDARY_PROP1_TYPE), ALLOCATABLE, DIMENSION(:) :: B1_DUMMY
TYPE (BOUNDARY_PROP2_TYPE), ALLOCATABLE, DIMENSION(:) :: B2_DUMMY
TYPE (BOUNDARY_RADIA_TYPE), ALLOCATABLE, DIMENSION(:) :: BR_DUMMY
TYPE (BOUNDARY_ONE_D_TYPE), ALLOCATABLE, DIMENSION(:) :: OD_DUMMY
TYPE (BOUNDARY_THR_D_TYPE), ALLOCATABLE, DIMENSION(:) :: TD_DUMMY
TYPE (MESH_TYPE), POINTER :: M
TYPE (SURFACE_TYPE), POINTER :: SF
TYPE (WALL_TYPE), POINTER :: WC
TYPE (THIN_WALL_TYPE), POINTER :: TW
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE (CFACE_TYPE), POINTER :: CFA
TYPE (STORAGE_TYPE), POINTER :: OS_DUMMY

M => MESHES(NM)

! Check to see if there is enough space within the LAGRANGIAN_PARTICLE, WALL, or CFACE arrays

IF (PRESENT(LP_INDEX)) THEN

   LPC => LAGRANGIAN_PARTICLE_CLASS(LPC_INDEX)
   SF => SURFACE(LPC%SURF_INDEX)
   N_NEW_STORAGE_SLOTS = LPC%NEW_PARTICLE_INCREMENT
   IF (LP_INDEX>M%NLPDIM) THEN
      ALLOCATE(LP_DUMMY(M%NLPDIM))
      LP_DUMMY(1:M%NLPDIM) = M%LAGRANGIAN_PARTICLE(1:M%NLPDIM)
      DEALLOCATE(M%LAGRANGIAN_PARTICLE)
      ALLOCATE(M%LAGRANGIAN_PARTICLE(M%NLPDIM+N_NEW_STORAGE_SLOTS))
      M%LAGRANGIAN_PARTICLE(1:M%NLPDIM) = LP_DUMMY
      DEALLOCATE(LP_DUMMY)
      M%NLPDIM = M%NLPDIM + N_NEW_STORAGE_SLOTS
   ENDIF
   LP => M%LAGRANGIAN_PARTICLE(LP_INDEX)
   BC_INDEX = LP%BC_INDEX
   OD_INDEX = LP%OD_INDEX
   B1_INDEX = LP%B1_INDEX
   B2_INDEX = LP%B2_INDEX
   BR_INDEX = LP%BR_INDEX
   IF (LPC%INCLUDE_BOUNDARY_COORD_TYPE) CALL ALLOCATE_BOUNDARY_COORD_ARRAYS
   IF (LPC%INCLUDE_BOUNDARY_ONE_D_TYPE) CALL ALLOCATE_BOUNDARY_ONE_D_ARRAYS
   IF (LPC%INCLUDE_BOUNDARY_PROP1_TYPE) CALL ALLOCATE_BOUNDARY_PROP1_ARRAYS
   IF (LPC%INCLUDE_BOUNDARY_PROP2_TYPE) CALL ALLOCATE_BOUNDARY_PROP2_ARRAYS
   IF (LPC%INCLUDE_BOUNDARY_RADIA_TYPE) CALL ALLOCATE_BOUNDARY_RADIA_ARRAYS
   LP%BC_INDEX = BC_INDEX
   LP%OD_INDEX = OD_INDEX
   LP%B1_INDEX = B1_INDEX
   LP%B2_INDEX = B2_INDEX
   LP%BR_INDEX = BR_INDEX

ELSEIF (PRESENT(WALL_INDEX)) THEN

   SF => SURFACE(SURF_INDEX)
   N_NEW_STORAGE_SLOTS = 1000
   IF (WALL_INDEX>M%N_WALL_CELLS_DIM) THEN
      ALLOCATE(WALL_DUMMY(0:M%N_WALL_CELLS_DIM))
      WALL_DUMMY(0:M%N_WALL_CELLS_DIM) = M%WALL(0:M%N_WALL_CELLS_DIM)
      DEALLOCATE(M%WALL)
      ALLOCATE(M%WALL(0:M%N_WALL_CELLS_DIM+N_NEW_STORAGE_SLOTS))
      M%WALL(0:M%N_WALL_CELLS_DIM) = WALL_DUMMY(0:M%N_WALL_CELLS_DIM)
      DEALLOCATE(WALL_DUMMY)
      M%N_WALL_CELLS_DIM = M%N_WALL_CELLS_DIM + N_NEW_STORAGE_SLOTS
   ENDIF
   WC => M%WALL(WALL_INDEX)
   BC_INDEX = WC%BC_INDEX
   OD_INDEX = WC%OD_INDEX
   TD_INDEX = WC%TD_INDEX
   B1_INDEX = WC%B1_INDEX
   B2_INDEX = WC%B2_INDEX
   BR_INDEX = WC%BR_INDEX
   IF (SF%INCLUDE_BOUNDARY_COORD_TYPE) CALL ALLOCATE_BOUNDARY_COORD_ARRAYS
   IF (SF%INCLUDE_BOUNDARY_ONE_D_TYPE) CALL ALLOCATE_BOUNDARY_ONE_D_ARRAYS
   IF (SF%INCLUDE_BOUNDARY_THR_D_TYPE) CALL ALLOCATE_BOUNDARY_THR_D_ARRAYS
   IF (SF%INCLUDE_BOUNDARY_PROP1_TYPE) CALL ALLOCATE_BOUNDARY_PROP1_ARRAYS
   IF (SF%INCLUDE_BOUNDARY_PROP2_TYPE) CALL ALLOCATE_BOUNDARY_PROP2_ARRAYS
   IF (SF%INCLUDE_BOUNDARY_RADIA_TYPE) CALL ALLOCATE_BOUNDARY_RADIA_ARRAYS
   WC%BC_INDEX = BC_INDEX
   WC%OD_INDEX = OD_INDEX
   WC%TD_INDEX = TD_INDEX
   WC%B1_INDEX = B1_INDEX
   WC%B2_INDEX = B2_INDEX
   WC%BR_INDEX = BR_INDEX

   WC%N_REALS = 0
   WC%N_INTEGERS = 0
   WC%N_LOGICALS = 0
   CALL PACK_WALL(NM,OS_DUMMY,WC,SURF_INDEX,WC%N_REALS,WC%N_INTEGERS,WC%N_LOGICALS,UNPACK_IT=.FALSE.,COUNT_ONLY=.TRUE.)

ELSEIF (PRESENT(THIN_WALL_INDEX)) THEN

   SF => SURFACE(SURF_INDEX)
   N_NEW_STORAGE_SLOTS = 1000
   IF (THIN_WALL_INDEX>M%N_THIN_WALL_CELLS_DIM) THEN
      ALLOCATE(THIN_WALL_DUMMY(1:M%N_THIN_WALL_CELLS_DIM))
      THIN_WALL_DUMMY(1:M%N_THIN_WALL_CELLS_DIM) = M%THIN_WALL(1:M%N_THIN_WALL_CELLS_DIM)
      DEALLOCATE(M%THIN_WALL)
      ALLOCATE(M%THIN_WALL(1:M%N_THIN_WALL_CELLS_DIM+N_NEW_STORAGE_SLOTS))
      M%THIN_WALL(1:M%N_THIN_WALL_CELLS_DIM) = THIN_WALL_DUMMY(1:M%N_THIN_WALL_CELLS_DIM)
      DEALLOCATE(THIN_WALL_DUMMY)
      M%N_THIN_WALL_CELLS_DIM = M%N_THIN_WALL_CELLS_DIM + N_NEW_STORAGE_SLOTS
   ENDIF
   TW => M%THIN_WALL(THIN_WALL_INDEX)
   BC_INDEX = TW%BC_INDEX
   OD_INDEX = TW%OD_INDEX
   TD_INDEX = TW%TD_INDEX
   B1_INDEX = TW%B1_INDEX
   IF (SF%INCLUDE_BOUNDARY_COORD_TYPE) CALL ALLOCATE_BOUNDARY_COORD_ARRAYS
   IF (SF%INCLUDE_BOUNDARY_PROP1_TYPE) CALL ALLOCATE_BOUNDARY_PROP1_ARRAYS
   IF (SF%INCLUDE_BOUNDARY_ONE_D_TYPE) CALL ALLOCATE_BOUNDARY_ONE_D_ARRAYS
   IF (SF%INCLUDE_BOUNDARY_THR_D_TYPE) CALL ALLOCATE_BOUNDARY_THR_D_ARRAYS
   TW%BC_INDEX = BC_INDEX
   TW%OD_INDEX = OD_INDEX
   TW%TD_INDEX = TD_INDEX
   TW%B1_INDEX = B1_INDEX

   CALL PACK_THIN_WALL(NM,OS_DUMMY,TW,SURF_INDEX,TW%N_REALS,TW%N_INTEGERS,TW%N_LOGICALS,UNPACK_IT=.FALSE.,COUNT_ONLY=.TRUE.)

ELSEIF (PRESENT(CFACE_INDEX)) THEN

   SF => SURFACE(SURF_INDEX)
   N_NEW_STORAGE_SLOTS = 1000
   IF (CFACE_INDEX>M%N_CFACE_CELLS_DIM) THEN
      ALLOCATE(CFACE_DUMMY(0:M%N_CFACE_CELLS_DIM))
      CFACE_DUMMY(0:M%N_CFACE_CELLS_DIM) = M%CFACE(0:M%N_CFACE_CELLS_DIM)
      DEALLOCATE(M%CFACE)
      ALLOCATE(M%CFACE(0:M%N_CFACE_CELLS_DIM+N_NEW_STORAGE_SLOTS))
      M%CFACE(0:M%N_CFACE_CELLS_DIM) = CFACE_DUMMY(0:M%N_CFACE_CELLS_DIM)
      DEALLOCATE(CFACE_DUMMY)
      M%N_CFACE_CELLS_DIM = M%N_CFACE_CELLS_DIM + N_NEW_STORAGE_SLOTS
   ENDIF
   CFA => M%CFACE(CFACE_INDEX)
   BC_INDEX = CFA%BC_INDEX
   OD_INDEX = CFA%OD_INDEX
   B1_INDEX = CFA%B1_INDEX
   B2_INDEX = CFA%B2_INDEX
   BR_INDEX = CFA%BR_INDEX
   IF (SF%INCLUDE_BOUNDARY_COORD_TYPE) CALL ALLOCATE_BOUNDARY_COORD_ARRAYS
   IF (SF%INCLUDE_BOUNDARY_ONE_D_TYPE) CALL ALLOCATE_BOUNDARY_ONE_D_ARRAYS
   IF (SF%INCLUDE_BOUNDARY_PROP1_TYPE) CALL ALLOCATE_BOUNDARY_PROP1_ARRAYS
   IF (SF%INCLUDE_BOUNDARY_PROP2_TYPE) CALL ALLOCATE_BOUNDARY_PROP2_ARRAYS
   IF (SF%INCLUDE_BOUNDARY_RADIA_TYPE) CALL ALLOCATE_BOUNDARY_RADIA_ARRAYS
   CFA%CFACE_INDEX = CFACE_INDEX
   CFA%BC_INDEX = BC_INDEX
   CFA%OD_INDEX = OD_INDEX
   CFA%B1_INDEX = B1_INDEX
   CFA%B2_INDEX = B2_INDEX
   CFA%BR_INDEX = BR_INDEX

   CALL PACK_CFACE(NM,OS_DUMMY,CFA,SURF_INDEX,CFA%N_REALS,CFA%N_INTEGERS,CFA%N_LOGICALS,UNPACK_IT=.FALSE.,COUNT_ONLY=.TRUE.)

ENDIF

CONTAINS

SUBROUTINE ALLOCATE_BOUNDARY_COORD_ARRAYS

IF (BC_INDEX==0 .AND. M%N_BOUNDARY_COORD_DIM>0) THEN
   NAS = M%NEXT_AVAILABLE_BOUNDARY_COORD_SLOT
   DO I=NAS,M%N_BOUNDARY_COORD_DIM
      IF (M%BOUNDARY_COORD_OCCUPANCY(I)==0) EXIT
      M%NEXT_AVAILABLE_BOUNDARY_COORD_SLOT = M%NEXT_AVAILABLE_BOUNDARY_COORD_SLOT + 1
   ENDDO
   BC_INDEX = M%NEXT_AVAILABLE_BOUNDARY_COORD_SLOT
ENDIF

IF (BC_INDEX==0 .OR. BC_INDEX>M%N_BOUNDARY_COORD_DIM) THEN  ! There are no open slots for boundary coordinates
   ALLOCATE(BC_DUMMY(1:M%N_BOUNDARY_COORD_DIM+N_NEW_STORAGE_SLOTS))
   IF (M%N_BOUNDARY_COORD_DIM>0) BC_DUMMY(1:M%N_BOUNDARY_COORD_DIM) = M%BOUNDARY_COORD(1:M%N_BOUNDARY_COORD_DIM)
   CALL MOVE_ALLOC(BC_DUMMY,M%BOUNDARY_COORD)
   ALLOCATE(DUMMY(1:M%N_BOUNDARY_COORD_DIM+N_NEW_STORAGE_SLOTS)) ; DUMMY = 0
   IF (M%N_BOUNDARY_COORD_DIM>0) DUMMY(1:M%N_BOUNDARY_COORD_DIM) = M%BOUNDARY_COORD_OCCUPANCY(1:M%N_BOUNDARY_COORD_DIM)
   CALL MOVE_ALLOC(DUMMY,M%BOUNDARY_COORD_OCCUPANCY)
   M%N_BOUNDARY_COORD_DIM = M%N_BOUNDARY_COORD_DIM + N_NEW_STORAGE_SLOTS
   IF (BC_INDEX==0) BC_INDEX = 1
ENDIF

M%BOUNDARY_COORD_OCCUPANCY(BC_INDEX) = 1

END SUBROUTINE ALLOCATE_BOUNDARY_COORD_ARRAYS


SUBROUTINE ALLOCATE_BOUNDARY_ONE_D_ARRAYS

TYPE(BOUNDARY_ONE_D_TYPE), POINTER :: ONE_D

IF (OD_INDEX==0 .AND. M%N_BOUNDARY_ONE_D_DIM>0) THEN
   NAS = M%NEXT_AVAILABLE_BOUNDARY_ONE_D_SLOT
   DO I=NAS,M%N_BOUNDARY_ONE_D_DIM
      IF (M%BOUNDARY_ONE_D_OCCUPANCY(I)==0) EXIT
      M%NEXT_AVAILABLE_BOUNDARY_ONE_D_SLOT = M%NEXT_AVAILABLE_BOUNDARY_ONE_D_SLOT + 1
   ENDDO
   OD_INDEX = M%NEXT_AVAILABLE_BOUNDARY_ONE_D_SLOT
ENDIF

IF (OD_INDEX==0 .OR. OD_INDEX>M%N_BOUNDARY_ONE_D_DIM) THEN  ! There are no open slots for boundary coordinates
   ALLOCATE(OD_DUMMY(1:M%N_BOUNDARY_ONE_D_DIM+N_NEW_STORAGE_SLOTS))
   IF (M%N_BOUNDARY_ONE_D_DIM>0) OD_DUMMY(1:M%N_BOUNDARY_ONE_D_DIM) = M%BOUNDARY_ONE_D(1:M%N_BOUNDARY_ONE_D_DIM)
   CALL MOVE_ALLOC(OD_DUMMY,M%BOUNDARY_ONE_D)
   ALLOCATE(DUMMY(1:M%N_BOUNDARY_ONE_D_DIM+N_NEW_STORAGE_SLOTS)) ; DUMMY = 0
   IF (M%N_BOUNDARY_ONE_D_DIM>0) DUMMY(1:M%N_BOUNDARY_ONE_D_DIM) = M%BOUNDARY_ONE_D_OCCUPANCY(1:M%N_BOUNDARY_ONE_D_DIM)
   CALL MOVE_ALLOC(DUMMY,M%BOUNDARY_ONE_D_OCCUPANCY)
   M%N_BOUNDARY_ONE_D_DIM = M%N_BOUNDARY_ONE_D_DIM + N_NEW_STORAGE_SLOTS
   IF (OD_INDEX==0) OD_INDEX = 1
ENDIF

M%BOUNDARY_ONE_D_OCCUPANCY(OD_INDEX) = 1

ONE_D => M%BOUNDARY_ONE_D(OD_INDEX)

IF (ONE_D%SURF_INDEX==SURF_INDEX) THEN
   IF (.NOT. PRESENT(LPC_INDEX)) RETURN
   IF (LAGRANGIAN_PARTICLE_CLASS(LPC_INDEX)%MONODISPERSE) RETURN
ENDIF

! Assume that most ONE_D arrays will have the same dimensions of the corresponding SURF_INDEX. If not, that will be adjusted later

ONE_D%SURF_INDEX  = SURF_INDEX
ONE_D%N_CELLS_MAX = SF%N_CELLS_MAX
ONE_D%N_CELLS_INI = SF%N_CELLS_INI
IF (SURFACE(SURF_INDEX)%PYROLYSIS_MODEL==PYROLYSIS_PREDICTED) ONE_D%N_CELLS_OLD = SF%N_CELLS_MAX
ONE_D%N_LAYERS    = SF%N_LAYERS
ONE_D%N_MATL      = SF%N_MATL
ONE_D%N_LPC       = SF%N_LPC

! If not already done, allocate the ONE_D arrays using the bounds that have just been set.

CALL REALLOCATE_BOUNDARY_ONE_D(ONE_D)

! Set initial values for ONE_D variables

CALL INITIALIZE_BOUNDARY_ONE_D(NM,OD_INDEX,SURF_INDEX)

END SUBROUTINE ALLOCATE_BOUNDARY_ONE_D_ARRAYS


SUBROUTINE ALLOCATE_BOUNDARY_THR_D_ARRAYS

TYPE(BOUNDARY_THR_D_TYPE), POINTER :: THR_D

IF (TD_INDEX==0 .AND. M%N_BOUNDARY_THR_D_DIM>0) THEN
   NAS = M%NEXT_AVAILABLE_BOUNDARY_THR_D_SLOT
   DO I=NAS,M%N_BOUNDARY_THR_D_DIM
      IF (M%BOUNDARY_THR_D_OCCUPANCY(I)==0) EXIT
      M%NEXT_AVAILABLE_BOUNDARY_THR_D_SLOT = M%NEXT_AVAILABLE_BOUNDARY_THR_D_SLOT + 1
   ENDDO
   TD_INDEX = M%NEXT_AVAILABLE_BOUNDARY_THR_D_SLOT
ENDIF

IF (TD_INDEX==0 .OR. TD_INDEX>M%N_BOUNDARY_THR_D_DIM) THEN  ! There are no open slots for boundary coordinates
   ALLOCATE(TD_DUMMY(1:M%N_BOUNDARY_THR_D_DIM+N_NEW_STORAGE_SLOTS))
   IF (M%N_BOUNDARY_THR_D_DIM>0) TD_DUMMY(1:M%N_BOUNDARY_THR_D_DIM) = M%BOUNDARY_THR_D(1:M%N_BOUNDARY_THR_D_DIM)
   CALL MOVE_ALLOC(TD_DUMMY,M%BOUNDARY_THR_D)
   ALLOCATE(DUMMY(1:M%N_BOUNDARY_THR_D_DIM+N_NEW_STORAGE_SLOTS)) ; DUMMY = 0
   IF (M%N_BOUNDARY_THR_D_DIM>0) DUMMY(1:M%N_BOUNDARY_THR_D_DIM) = M%BOUNDARY_THR_D_OCCUPANCY(1:M%N_BOUNDARY_THR_D_DIM)
   CALL MOVE_ALLOC(DUMMY,M%BOUNDARY_THR_D_OCCUPANCY)
   M%N_BOUNDARY_THR_D_DIM = M%N_BOUNDARY_THR_D_DIM + N_NEW_STORAGE_SLOTS
   IF (TD_INDEX==0) TD_INDEX = 1
ENDIF

M%BOUNDARY_THR_D_OCCUPANCY(TD_INDEX) = 1

THR_D => M%BOUNDARY_THR_D(TD_INDEX)

END SUBROUTINE ALLOCATE_BOUNDARY_THR_D_ARRAYS


SUBROUTINE ALLOCATE_BOUNDARY_PROP1_ARRAYS

USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,N_SURFACE_DENSITY_SPECIES
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1

IF (B1_INDEX==0 .AND. M%N_BOUNDARY_PROP1_DIM>0) THEN
   NAS = M%NEXT_AVAILABLE_BOUNDARY_PROP1_SLOT
   DO I=NAS,M%N_BOUNDARY_PROP1_DIM
      IF (M%BOUNDARY_PROP1_OCCUPANCY(I)==0) EXIT
      M%NEXT_AVAILABLE_BOUNDARY_PROP1_SLOT = M%NEXT_AVAILABLE_BOUNDARY_PROP1_SLOT + 1
   ENDDO
   B1_INDEX = M%NEXT_AVAILABLE_BOUNDARY_PROP1_SLOT
ENDIF

IF (B1_INDEX==0 .OR. B1_INDEX>M%N_BOUNDARY_PROP1_DIM) THEN  ! There are no open slots for boundary coordinates
   ALLOCATE(B1_DUMMY(1:M%N_BOUNDARY_PROP1_DIM+N_NEW_STORAGE_SLOTS))
   IF (M%N_BOUNDARY_PROP1_DIM>0) B1_DUMMY(1:M%N_BOUNDARY_PROP1_DIM) = M%BOUNDARY_PROP1(1:M%N_BOUNDARY_PROP1_DIM)
   CALL MOVE_ALLOC(B1_DUMMY,M%BOUNDARY_PROP1)
   ALLOCATE(DUMMY(1:M%N_BOUNDARY_PROP1_DIM+N_NEW_STORAGE_SLOTS)) ; DUMMY = 0
   IF (M%N_BOUNDARY_PROP1_DIM>0) DUMMY(1:M%N_BOUNDARY_PROP1_DIM) = M%BOUNDARY_PROP1_OCCUPANCY(1:M%N_BOUNDARY_PROP1_DIM)
   CALL MOVE_ALLOC(DUMMY,M%BOUNDARY_PROP1_OCCUPANCY)
   M%N_BOUNDARY_PROP1_DIM = M%N_BOUNDARY_PROP1_DIM + N_NEW_STORAGE_SLOTS
   IF (B1_INDEX==0) B1_INDEX = 1
ENDIF

M%BOUNDARY_PROP1_OCCUPANCY(B1_INDEX) = 1

B1 => M%BOUNDARY_PROP1(B1_INDEX)

IF (B1%SURF_INDEX==SURF_INDEX) THEN
   ALREADY_ALLOCATED = .TRUE.
ELSE
   ALREADY_ALLOCATED = .FALSE.
ENDIF

B1%SURF_INDEX = SURF_INDEX

IF (.NOT.ALREADY_ALLOCATED) THEN
   IF (ALLOCATED(B1%M_DOT_G_PP_ACTUAL)) DEALLOCATE(B1%M_DOT_G_PP_ACTUAL)
   ALLOCATE(B1%M_DOT_G_PP_ACTUAL(N_TRACKED_SPECIES))
   IF (ALLOCATED(B1%M_DOT_G_PP_ADJUST)) DEALLOCATE(B1%M_DOT_G_PP_ADJUST)
   ALLOCATE(B1%M_DOT_G_PP_ADJUST(N_TRACKED_SPECIES))
   IF (ALLOCATED(B1%ZZ_F)) DEALLOCATE(B1%ZZ_F)
   ALLOCATE(B1%ZZ_F(N_TRACKED_SPECIES))
   IF (ALLOCATED(B1%ZZ_G)) DEALLOCATE(B1%ZZ_G)
   ALLOCATE(B1%ZZ_G(N_TRACKED_SPECIES))
   IF (ALLOCATED(B1%RHO_D_F)) DEALLOCATE(B1%RHO_D_F)
   ALLOCATE(B1%RHO_D_F(N_TRACKED_SPECIES))
   IF (ALLOCATED(B1%RHO_D_DZDN_F)) DEALLOCATE(B1%RHO_D_DZDN_F)
   ALLOCATE(B1%RHO_D_DZDN_F(N_TRACKED_SPECIES))
   IF (ALLOCATED(B1%AWM_AEROSOL)) DEALLOCATE(B1%AWM_AEROSOL)
   ALLOCATE(B1%AWM_AEROSOL(N_SURFACE_DENSITY_SPECIES))
   IF (ALLOCATED(B1%QDOTPP_INT)) DEALLOCATE(B1%QDOTPP_INT)
   ALLOCATE(B1%QDOTPP_INT(SURFACE(SURF_INDEX)%N_QDOTPP_REF))
   IF (ALLOCATED(B1%Q_IN_SMOOTH_INT)) DEALLOCATE(B1%Q_IN_SMOOTH_INT)
   ALLOCATE(B1%Q_IN_SMOOTH_INT(SURFACE(SURF_INDEX)%N_THICK_REF))

   CALL INITIALIZE_BOUNDARY_PROP1(NM,B1_INDEX,OD_INDEX,SURF_INDEX)
ENDIF

END SUBROUTINE ALLOCATE_BOUNDARY_PROP1_ARRAYS


SUBROUTINE ALLOCATE_BOUNDARY_PROP2_ARRAYS

USE GLOBAL_CONSTANTS, ONLY: TMPA,N_LP_ARRAY_INDICES
TYPE(BOUNDARY_PROP2_TYPE), POINTER :: B2

IF (B2_INDEX==0 .AND. M%N_BOUNDARY_PROP2_DIM>0) THEN
   NAS = M%NEXT_AVAILABLE_BOUNDARY_PROP2_SLOT
   DO I=NAS,M%N_BOUNDARY_PROP2_DIM
      IF (M%BOUNDARY_PROP2_OCCUPANCY(I)==0) EXIT
      M%NEXT_AVAILABLE_BOUNDARY_PROP2_SLOT = M%NEXT_AVAILABLE_BOUNDARY_PROP2_SLOT + 1
   ENDDO
   B2_INDEX = M%NEXT_AVAILABLE_BOUNDARY_PROP2_SLOT
ENDIF

IF (B2_INDeX==0 .OR. B2_INDEX>M%N_BOUNDARY_PROP2_DIM) THEN  ! There are no open slots for boundary coordinates
   ALLOCATE(B2_DUMMY(1:M%N_BOUNDARY_PROP2_DIM+N_NEW_STORAGE_SLOTS))
   IF (M%N_BOUNDARY_PROP2_DIM>0) B2_DUMMY(1:M%N_BOUNDARY_PROP2_DIM) = M%BOUNDARY_PROP2(1:M%N_BOUNDARY_PROP2_DIM)
   CALL MOVE_ALLOC(B2_DUMMY,M%BOUNDARY_PROP2)
   ALLOCATE(DUMMY(1:M%N_BOUNDARY_PROP2_DIM+N_NEW_STORAGE_SLOTS)) ; DUMMY = 0
   IF (M%N_BOUNDARY_PROP2_DIM>0) DUMMY(1:M%N_BOUNDARY_PROP2_DIM) = M%BOUNDARY_PROP2_OCCUPANCY(1:M%N_BOUNDARY_PROP2_DIM)
   CALL MOVE_ALLOC(DUMMY,M%BOUNDARY_PROP2_OCCUPANCY)
   M%N_BOUNDARY_PROP2_DIM = M%N_BOUNDARY_PROP2_DIM + N_NEW_STORAGE_SLOTS
   IF (B2_INDEX==0) B2_INDEX = 1
ENDIF

M%BOUNDARY_PROP2_OCCUPANCY(B2_INDEX) = 1

B2 => M%BOUNDARY_PROP2(B2_INDEX)

IF (B2%SURF_INDEX==SURF_INDEX) THEN
   ALREADY_ALLOCATED = .TRUE.
ELSE
   ALREADY_ALLOCATED = .FALSE.
ENDIF

IF (.NOT.ALREADY_ALLOCATED) THEN
   IF (ALLOCATED(B2%A_LP_MPUA)) DEALLOCATE(B2%A_LP_MPUA)
   ALLOCATE(B2%A_LP_MPUA(N_LP_ARRAY_INDICES))
   IF (ALLOCATED(B2%LP_CPUA)) DEALLOCATE(B2%LP_CPUA)
   ALLOCATE(B2%LP_CPUA(N_LP_ARRAY_INDICES))
   IF (ALLOCATED(B2%LP_EMPUA)) DEALLOCATE(B2%LP_EMPUA)
   ALLOCATE(B2%LP_EMPUA(N_LP_ARRAY_INDICES))
   IF (ALLOCATED(B2%LP_MPUA)) DEALLOCATE(B2%LP_MPUA)
   ALLOCATE(B2%LP_MPUA(N_LP_ARRAY_INDICES))
   IF (ALLOCATED(B2%LP_TEMP)) DEALLOCATE(B2%LP_TEMP)
   ALLOCATE(B2%LP_TEMP(N_LP_ARRAY_INDICES))
ENDIF

B2%A_LP_MPUA = 0._EB
B2%LP_CPUA = 0._EB
B2%LP_EMPUA = 0._EB
B2%LP_MPUA = 0._EB
B2%LP_TEMP = TMPA

B2%SURF_INDEX = SURF_INDEX

END SUBROUTINE ALLOCATE_BOUNDARY_PROP2_ARRAYS


SUBROUTINE ALLOCATE_BOUNDARY_RADIA_ARRAYS

USE GLOBAL_CONSTANTS, ONLY: TMPA4,SIGMA,PI,NUMBER_SPECTRAL_BANDS,NUMBER_RADIATION_ANGLES
TYPE(BOUNDARY_RADIA_TYPE), POINTER :: BR
INTEGER :: NN

IF (BR_INDEX==0 .AND. M%N_BOUNDARY_RADIA_DIM>0) THEN
   NAS = M%NEXT_AVAILABLE_BOUNDARY_RADIA_SLOT
   DO I=NAS,M%N_BOUNDARY_RADIA_DIM
      IF (M%BOUNDARY_RADIA_OCCUPANCY(I)==0) EXIT
      M%NEXT_AVAILABLE_BOUNDARY_RADIA_SLOT = M%NEXT_AVAILABLE_BOUNDARY_RADIA_SLOT + 1
   ENDDO
   BR_INDEX = M%NEXT_AVAILABLE_BOUNDARY_RADIA_SLOT
ENDIF

IF (BR_INDEX==0 .OR. BR_INDEX>M%N_BOUNDARY_RADIA_DIM) THEN  ! There are no open slots for boundary coordinates
   ALLOCATE(BR_DUMMY(1:M%N_BOUNDARY_RADIA_DIM+N_NEW_STORAGE_SLOTS))
   IF (M%N_BOUNDARY_RADIA_DIM>0) BR_DUMMY(1:M%N_BOUNDARY_RADIA_DIM) = M%BOUNDARY_RADIA(1:M%N_BOUNDARY_RADIA_DIM)
   CALL MOVE_ALLOC(BR_DUMMY,M%BOUNDARY_RADIA)
   ALLOCATE(DUMMY(1:M%N_BOUNDARY_RADIA_DIM+N_NEW_STORAGE_SLOTS)) ; DUMMY = 0
   IF (M%N_BOUNDARY_RADIA_DIM>0) DUMMY(1:M%N_BOUNDARY_RADIA_DIM) = M%BOUNDARY_RADIA_OCCUPANCY(1:M%N_BOUNDARY_RADIA_DIM)
   CALL MOVE_ALLOC(DUMMY,M%BOUNDARY_RADIA_OCCUPANCY)
   M%N_BOUNDARY_RADIA_DIM = M%N_BOUNDARY_RADIA_DIM + N_NEW_STORAGE_SLOTS
   IF (BR_INDEX==0) BR_INDEX = 1
ENDIF

M%BOUNDARY_RADIA_OCCUPANCY(BR_INDEX) = 1

BR => M%BOUNDARY_RADIA(BR_INDEX)

IF (.NOT.ALLOCATED(BR%BAND)) ALLOCATE(BR%BAND(NUMBER_SPECTRAL_BANDS))
DO NN=1,NUMBER_SPECTRAL_BANDS
   IF (.NOT.ALLOCATED(BR%BAND(NN)%ILW)) ALLOCATE(BR%BAND(NN)%ILW(NUMBER_RADIATION_ANGLES))
   BR%BAND(NN)%ILW = 0._EB
ENDDO

IF (.NOT.ALLOCATED(BR%IL)) ALLOCATE(BR%IL(NUMBER_SPECTRAL_BANDS))
BR%IL = SIGMA*TMPA4/PI

END SUBROUTINE ALLOCATE_BOUNDARY_RADIA_ARRAYS

END SUBROUTINE ALLOCATE_STORAGE


!> \brief Pack a single LAGRANGIAN_PARTICLE into REAL, INTEGER, and LOGICAL 1-D arrays
!> \param NM Mesh index
!> \param OS Pointer to storage array
!> \param LP Pointer to particle
!> \param LPC_INDEX Class index of the particle
!> \param RC Counter of real variables
!> \param IC Counter of integer variables
!> \param LC Counter of logical variables
!> \param UNPACK_IT Logical indicating whether the data is to be packed or unpacked
!> \param COUNT_ONLY Logical indicating whether the aim is simply to count variables

SUBROUTINE PACK_PARTICLE(NM,OS,LP,LPC_INDEX,RC,IC,LC,UNPACK_IT,COUNT_ONLY)

USE COMP_OPERATORS, ONLY: EQUATE
INTEGER, INTENT(IN) :: NM,LPC_INDEX
LOGICAL, INTENT(IN) :: UNPACK_IT,COUNT_ONLY
INTEGER, INTENT(INOUT) :: IC,RC,LC
INTEGER :: BC_INDEX,OD_INDEX,B1_INDEX,B2_INDEX,BR_INDEX
TYPE(LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE(STORAGE_TYPE), POINTER :: OS
TYPE(SURFACE_TYPE), POINTER :: SF

LPC => LAGRANGIAN_PARTICLE_CLASS(LPC_INDEX)
SF => SURFACE(LPC%SURF_INDEX)

BC_INDEX = LP%BC_INDEX
OD_INDEX = LP%OD_INDEX
B1_INDEX = LP%B1_INDEX
B2_INDEX = LP%B2_INDEX
BR_INDEX = LP%BR_INDEX

! Assign integer pointers and values

IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),LP%TAG,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),LP%CLASS_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),LP%ORIENTATION_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),LP%WALL_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),LP%DUCT_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),LP%DUCT_CELL_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),LP%INIT_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),LP%CFACE_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),LP%PROP_INDEX,UNPACK_IT)

! Assign logical pointers and values

LC=LC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%LOGICALS(LC),LP%SHOW,UNPACK_IT)
LC=LC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%LOGICALS(LC),LP%SPLAT,UNPACK_IT)
LC=LC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%LOGICALS(LC),LP%EMBER,UNPACK_IT)
LC=LC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%LOGICALS(LC),LP%PATH_PARTICLE,UNPACK_IT)

! Assign real pointers and values

RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%U,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%V,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%W,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%PWT,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%RE,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%MASS,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%T_INSERT,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%DX,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%DY,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%DZ,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%LENGTH,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%C_DRAG,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%RADIUS,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%ACCEL_X,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%ACCEL_Y,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%ACCEL_Z,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),LP%RVC,UNPACK_IT)

! Boundary Coordinates

IF (LPC%INCLUDE_BOUNDARY_COORD_TYPE) CALL PACK_BOUNDARY_COORD(NM,IC,RC,OS,BC_INDEX,UNPACK_IT,COUNT_ONLY)
IF (LPC%INCLUDE_BOUNDARY_ONE_D_TYPE) CALL PACK_BOUNDARY_ONE_D(NM,IC,RC,LC,OS,OD_INDEX,UNPACK_IT,COUNT_ONLY)
IF (LPC%INCLUDE_BOUNDARY_PROP1_TYPE) CALL PACK_BOUNDARY_PROP1(NM,IC,RC,LC,OS,B1_INDEX,UNPACK_IT,COUNT_ONLY,LPC%SURF_INDEX)
IF (LPC%INCLUDE_BOUNDARY_PROP2_TYPE) CALL PACK_BOUNDARY_PROP2(NM,IC,RC,OS,B2_INDEX,UNPACK_IT,COUNT_ONLY)
IF (LPC%INCLUDE_BOUNDARY_RADIA_TYPE) CALL PACK_BOUNDARY_RADIA(NM,RC,OS,BR_INDEX,UNPACK_IT,COUNT_ONLY)

END SUBROUTINE PACK_PARTICLE


!> \brief Pack WALL components into REAL, INTEGER, and LOGICAL 1-D arrays
!> \param NM Mesh index
!> \param OS Pointer to storage array
!> \param WC Pointer to wall cell
!> \param SURF_INDEX Index of the surface type
!> \param RC Counter of real variables
!> \param IC Counter of integer variables
!> \param LC Counter of logical variables
!> \param UNPACK_IT Logical indicating whether the data is to be packed or unpacked
!> \param COUNT_ONLY Logical indicating whether the aim is simply to count variables

SUBROUTINE PACK_WALL(NM,OS,WC,SURF_INDEX,RC,IC,LC,UNPACK_IT,COUNT_ONLY)

USE COMP_OPERATORS, ONLY: EQUATE
INTEGER, INTENT(IN) :: NM,SURF_INDEX
LOGICAL, INTENT(IN) :: UNPACK_IT
INTEGER, INTENT(INOUT) :: IC,RC,LC
INTEGER :: BC_INDEX,OD_INDEX,B1_INDEX,B2_INDEX,BR_INDEX
LOGICAL, INTENT(IN) :: COUNT_ONLY
TYPE(WALL_TYPE), POINTER :: WC
TYPE(STORAGE_TYPE), POINTER :: OS
TYPE(SURFACE_TYPE), POINTER :: SF

! Point to the storage arrays that are sized for this particular surface type

SF => SURFACE(SURF_INDEX)

! Point to the wall cell that was just set up and extract the various boundary indices

BC_INDEX = WC%BC_INDEX
OD_INDEX = WC%OD_INDEX
B1_INDEX = WC%B1_INDEX
B2_INDEX = WC%B2_INDEX
BR_INDEX = WC%BR_INDEX

! Pack or unpack the integer values

IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),WC%SURF_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),WC%BOUNDARY_TYPE,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),WC%OBST_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),WC%VENT_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),WC%JD11_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),WC%JD12_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),WC%JD21_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),WC%JD22_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),WC%CUT_FACE_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),WC%N_REALS,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),WC%N_INTEGERS,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),WC%N_LOGICALS,UNPACK_IT)

! Pack or unpack the appropriate derived type variables tied to this wall cell

IF (SF%INCLUDE_BOUNDARY_COORD_TYPE) CALL PACK_BOUNDARY_COORD(NM,IC,RC,OS,BC_INDEX,UNPACK_IT,COUNT_ONLY)
IF (SF%INCLUDE_BOUNDARY_ONE_D_TYPE) CALL PACK_BOUNDARY_ONE_D(NM,IC,RC,LC,OS,OD_INDEX,UNPACK_IT,COUNT_ONLY)
IF (SF%INCLUDE_BOUNDARY_PROP1_TYPE) CALL PACK_BOUNDARY_PROP1(NM,IC,RC,LC,OS,B1_INDEX,UNPACK_IT,COUNT_ONLY,SURF_INDEX)
IF (SF%INCLUDE_BOUNDARY_PROP2_TYPE) CALL PACK_BOUNDARY_PROP2(NM,IC,RC,OS,B2_INDEX,UNPACK_IT,COUNT_ONLY)
IF (SF%INCLUDE_BOUNDARY_RADIA_TYPE) CALL PACK_BOUNDARY_RADIA(NM,RC,OS,BR_INDEX,UNPACK_IT,COUNT_ONLY)

END SUBROUTINE PACK_WALL


!> \brief Pack THIN_WALL components into REAL, INTEGER, and LOGICAL 1-D arrays
!> \param NM Mesh index
!> \param OS Pointer to storage array
!> \param TW Pointer to thin wall cell
!> \param SURF_INDEX Index of the surface type
!> \param RC Counter of real variables
!> \param IC Counter of integer variables
!> \param LC Counter of logical variables
!> \param UNPACK_IT Logical indicating whether the data is to be packed or unpacked
!> \param COUNT_ONLY Logical indicating whether the aim is simply to count variables

SUBROUTINE PACK_THIN_WALL(NM,OS,TW,SURF_INDEX,RC,IC,LC,UNPACK_IT,COUNT_ONLY)

USE COMP_OPERATORS, ONLY: EQUATE
INTEGER, INTENT(IN) :: NM,SURF_INDEX
LOGICAL, INTENT(IN) :: UNPACK_IT
INTEGER :: BC_INDEX,OD_INDEX
INTEGER, INTENT(INOUT) :: RC,IC,LC
LOGICAL, INTENT(IN) :: COUNT_ONLY
TYPE(THIN_WALL_TYPE), POINTER :: TW
TYPE(STORAGE_TYPE), POINTER :: OS
TYPE(SURFACE_TYPE), POINTER :: SF

! Point to the storage arrays that are sized for this particular surface type

SF => SURFACE(SURF_INDEX)

! Point to the wall cell that was just set up and extract the various boundary indices

BC_INDEX = TW%BC_INDEX
OD_INDEX = TW%OD_INDEX

! Pack or unpack the integer values

IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),TW%SURF_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),TW%BOUNDARY_TYPE,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),TW%OBST_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),TW%IEC,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),TW%WALL_INDEX_M,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),TW%WALL_INDEX_P,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),TW%N_REALS,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),TW%N_INTEGERS,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),TW%N_LOGICALS,UNPACK_IT)

! Pack or unpack the appropriate derived type variables tied to this wall cell

IF (SF%INCLUDE_BOUNDARY_COORD_TYPE) CALL PACK_BOUNDARY_COORD(NM,IC,RC,OS,BC_INDEX,UNPACK_IT,COUNT_ONLY)
IF (SF%INCLUDE_BOUNDARY_ONE_D_TYPE) CALL PACK_BOUNDARY_ONE_D(NM,IC,RC,LC,OS,OD_INDEX,UNPACK_IT,COUNT_ONLY)

END SUBROUTINE PACK_THIN_WALL


!> \brief Pack CFACE components into REAL, INTEGER, and LOGICAL 1-D arrays
!> \param NM Mesh index
!> \param OS Pointer to storage array
!> \param CFA Pointer to CFACE
!> \param SURF_INDEX Index of the surface type
!> \param RC Counter of real variables
!> \param IC Counter of integer variables
!> \param LC Counter of logical variables
!> \param UNPACK_IT Logical indicating whether the data is to be packed or unpacked
!> \param COUNT_ONLY Logical indicating whether the aim is simply to count variables

SUBROUTINE PACK_CFACE(NM,OS,CFA,SURF_INDEX,RC,IC,LC,UNPACK_IT,COUNT_ONLY)

USE COMP_OPERATORS, ONLY: EQUATE
INTEGER, INTENT(IN) :: NM,SURF_INDEX
LOGICAL, INTENT(IN) :: UNPACK_IT,COUNT_ONLY
INTEGER, INTENT(INOUT) :: RC,IC,LC
INTEGER :: BC_INDEX,OD_INDEX,B1_INDEX,B2_INDEX,BR_INDEX
TYPE(CFACE_TYPE), POINTER :: CFA
TYPE(STORAGE_TYPE), POINTER :: OS
TYPE(SURFACE_TYPE), POINTER :: SF

SF => SURFACE(SURF_INDEX)

BC_INDEX = CFA%BC_INDEX
OD_INDEX = CFA%OD_INDEX
B1_INDEX = CFA%B1_INDEX
B2_INDEX = CFA%B2_INDEX
BR_INDEX = CFA%BR_INDEX

! Assign integer pointers and values

IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),CFA%SURF_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),CFA%VENT_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),CFA%BOUNDARY_TYPE,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),CFA%CUT_FACE_IND1,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),CFA%CUT_FACE_IND2,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),CFA%N_REALS,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),CFA%N_INTEGERS,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),CFA%N_LOGICALS,UNPACK_IT)

! Assign and initialize reals

RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),CFA%AREA,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),CFA%DUNDT,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),CFA%RSUM_G,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),CFA%MU_G,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),CFA%PRES_BXN,UNPACK_IT)

! Assign and initialize logicals

IF (SF%INCLUDE_BOUNDARY_COORD_TYPE) CALL PACK_BOUNDARY_COORD(NM,IC,RC,OS,BC_INDEX,UNPACK_IT,COUNT_ONLY)
IF (SF%INCLUDE_BOUNDARY_ONE_D_TYPE) CALL PACK_BOUNDARY_ONE_D(NM,IC,RC,LC,OS,OD_INDEX,UNPACK_IT,COUNT_ONLY)
IF (SF%INCLUDE_BOUNDARY_PROP1_TYPE) CALL PACK_BOUNDARY_PROP1(NM,IC,RC,LC,OS,B1_INDEX,UNPACK_IT,COUNT_ONLY,SURF_INDEX)
IF (SF%INCLUDE_BOUNDARY_PROP2_TYPE) CALL PACK_BOUNDARY_PROP2(NM,IC,RC,OS,B2_INDEX,UNPACK_IT,COUNT_ONLY)
IF (SF%INCLUDE_BOUNDARY_RADIA_TYPE) CALL PACK_BOUNDARY_RADIA(NM,RC,OS,BR_INDEX,UNPACK_IT,COUNT_ONLY)

END SUBROUTINE PACK_CFACE


!> \brief Pack BOUNDARY_COORD into REAL, INTEGER, and LOGICAL 1-D arrays for a WALL, CFACE, or LAGRANGIAN_PARTICLE
!> \param NM Mesh index
!> \param IC Integer counter
!> \param RC Real counter
!> \param OS STORAGE_ARRAY name
!> \param BC_INDEX Index of the BOUNDARY_COORD array
!> \param UNPACK_IT Flag indicating whether the data is to be packed into the 1-D array or unpacked from it
!> \param COUNT_ONLY Flag signifying that only a variable count is to be done; no packing

SUBROUTINE PACK_BOUNDARY_COORD(NM,IC,RC,OS,BC_INDEX,UNPACK_IT,COUNT_ONLY)

USE COMP_OPERATORS, ONLY: EQUATE
INTEGER, INTENT(IN) :: NM,BC_INDEX
INTEGER, INTENT(INOUT) :: IC,RC
LOGICAL, INTENT(IN) :: UNPACK_IT,COUNT_ONLY
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(STORAGE_TYPE), POINTER :: OS

IF (.NOT.COUNT_ONLY) BC => MESHES(NM)%BOUNDARY_COORD(BC_INDEX)

IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),BC%II,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),BC%JJ,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),BC%KK,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),BC%IIG,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),BC%JJG,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),BC%KKG,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),BC%IOR,UNPACK_IT)

RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),BC%X,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),BC%Y,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),BC%Z,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),BC%X1,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),BC%X2,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),BC%Y1,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),BC%Y2,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),BC%Z1,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),BC%Z2,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),BC%NVEC(1),UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),BC%NVEC(2),UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),BC%NVEC(3),UNPACK_IT)

END SUBROUTINE PACK_BOUNDARY_COORD


!> \brief Pack BOUNDARY_ONE_D into REAL, INTEGER, and LOGICAL 1-D arrays for a WALL, CFACE, or LAGRANGIAN_PARTICLE
!> \param NM Mesh index
!> \param IC Integer counter
!> \param RC Real counter
!> \param LC Logical counter
!> \param OS STORAGE_ARRAY name
!> \param OD_INDEX Index of the BOUNDARY_ONE_D array
!> \param UNPACK_IT Flag indicating whether the data is to be packed into the 1-D array or unpacked from it
!> \param COUNT_ONLY Flag signifying that only a variable count is to be done; no packing

SUBROUTINE PACK_BOUNDARY_ONE_D(NM,IC,RC,LC,OS,OD_INDEX,UNPACK_IT,COUNT_ONLY)

USE GLOBAL_CONSTANTS, ONLY: CHECK_BOUNDARY_ONE_D_ARRAYS
USE COMP_OPERATORS, ONLY: EQUATE
INTEGER, INTENT(IN) :: NM,OD_INDEX
INTEGER, INTENT(INOUT) :: IC,RC,LC
LOGICAL, INTENT(IN) :: UNPACK_IT,COUNT_ONLY
INTEGER :: I1,NN,NL
TYPE(BOUNDARY_ONE_D_TYPE), POINTER :: ONE_D
TYPE(STORAGE_TYPE), POINTER :: OS

ONE_D => MESHES(NM)%BOUNDARY_ONE_D(OD_INDEX)

! Integer scalars

IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),ONE_D%SURF_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),ONE_D%N_CELLS_MAX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),ONE_D%N_CELLS_INI,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),ONE_D%N_CELLS_OLD,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),ONE_D%N_LAYERS,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),ONE_D%N_MATL,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),ONE_D%N_LPC,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),ONE_D%BACK_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),ONE_D%BACK_MESH,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),ONE_D%BACK_SURF,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),ONE_D%PYROLYSIS_MODEL,UNPACK_IT)

! Check if the array bounds are appropriate

IF (UNPACK_IT .AND. CHECK_BOUNDARY_ONE_D_ARRAYS) CALL REALLOCATE_BOUNDARY_ONE_D(ONE_D)

DO NL=1,ONE_D%N_LAYERS
   IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC) , ONE_D%N_LAYER_CELLS(NL) , UNPACK_IT)
   IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC) , ONE_D%REMESH_NWP(NL) , UNPACK_IT)
   IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC) , ONE_D%N_LAYER_CELLS_MAX(NL) , UNPACK_IT)
   IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC) , ONE_D%RAMP_IHS_INDEX(NL) , UNPACK_IT)
ENDDO

DO NN=1,ONE_D%N_MATL
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC) , ONE_D%MATL_INDEX(NN) , UNPACK_IT)
ENDDO

RC=RC+1                             ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC)   ,ONE_D%LAYER_DIVIDE_DEPTH        , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_MATL-1      ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%M_DOT_S_PP(1:RC-I1+1)  , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_CELLS_MAX   ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%X(0:RC-I1)             , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_CELLS_OLD-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%DX_OLD(1:RC-I1+1)      , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_CELLS_MAX+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%TMP(0:RC-I1)           , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_CELLS_MAX+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%DELTA_TMP(0:RC-I1)     , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_CELLS_MAX-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%RHO_C_S(1:RC-I1+1)     , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_CELLS_MAX+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%K_S(0:RC-I1)           , UNPACK_IT)

I1=RC+1 ; RC=I1+ONE_D%N_LAYERS-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%LAYER_THICKNESS(1:RC-I1+1)    , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_LAYERS-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%LAYER_THICKNESS_OLD(1:RC-I1+1), UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_LAYERS-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%MIN_LAYER_THICKNESS(1:RC-I1+1), UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_LAYERS-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%MIN_DIFFUSIVITY(1:RC-I1+1)    , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_LAYERS-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%DDSUM(1:RC-I1+1)              , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_LAYERS-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%SMALLEST_CELL_SIZE(1:RC-I1+1) , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_LAYERS-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%STRETCH_FACTOR(1:RC-I1+1)     , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_LAYERS-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%HEAT_SOURCE(1:RC-I1+1)        , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_LAYERS-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%CELL_SIZE_FACTOR(1:RC-I1+1)   , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_LAYERS-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%CELL_SIZE(1:RC-I1+1)          , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_LPC-1    ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%PART_MASS(1:RC-I1+1)          , UNPACK_IT)
I1=RC+1 ; RC=I1+ONE_D%N_LPC-1    ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC),ONE_D%PART_ENTHALPY(1:RC-I1+1)      , UNPACK_IT)

LC=LC+1                          ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%LOGICALS(LC)    , ONE_D%INTERNAL_RADIATION    , UNPACK_IT)
I1=LC+1 ; LC=I1+ONE_D%N_LAYERS-1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%LOGICALS(I1:LC) , ONE_D%HT3D_LAYER(1:LC-I1+1) , UNPACK_IT)

DO NN=1,ONE_D%N_MATL
   I1=RC+1 ; RC=I1+ONE_D%N_LAYERS-1
   IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC) , ONE_D%MATL_COMP(NN)%MASS_FRACTION(1:RC-I1+1) , UNPACK_IT)
   I1=RC+1 ; RC=I1+ONE_D%N_CELLS_MAX+1
   IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:RC) , ONE_D%MATL_COMP(NN)%RHO(0:RC-I1) , UNPACK_IT)
ENDDO

END SUBROUTINE PACK_BOUNDARY_ONE_D


SUBROUTINE REALLOCATE_BOUNDARY_ONE_D(ONE_D)

TYPE(BOUNDARY_ONE_D_TYPE), POINTER :: ONE_D
INTEGER :: NN

IF (ALLOCATED(ONE_D%N_LAYER_CELLS)) DEALLOCATE(ONE_D%N_LAYER_CELLS)           ; ALLOCATE(ONE_D%N_LAYER_CELLS(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%REMESH_NWP)) DEALLOCATE(ONE_D%REMESH_NWP)                 ; ALLOCATE(ONE_D%REMESH_NWP(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%N_LAYER_CELLS_MAX)) DEALLOCATE(ONE_D%N_LAYER_CELLS_MAX)   ; ALLOCATE(ONE_D%N_LAYER_CELLS_MAX(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%RAMP_IHS_INDEX)) DEALLOCATE(ONE_D%RAMP_IHS_INDEX)         ; ALLOCATE(ONE_D%RAMP_IHS_INDEX(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%MATL_INDEX)) DEALLOCATE(ONE_D%MATL_INDEX)                 ; ALLOCATE(ONE_D%MATL_INDEX(ONE_D%N_MATL))
IF (ALLOCATED(ONE_D%M_DOT_S_PP)) DEALLOCATE(ONE_D%M_DOT_S_PP) ; ALLOCATE(ONE_D%M_DOT_S_PP(ONE_D%N_MATL)) ; ONE_D%M_DOT_S_PP=0._EB
IF (ALLOCATED(ONE_D%X)) DEALLOCATE(ONE_D%X)                                   ; ALLOCATE(ONE_D%X(0:ONE_D%N_CELLS_MAX))
IF (ALLOCATED(ONE_D%DX_OLD)) DEALLOCATE(ONE_D%DX_OLD)                         ; ALLOCATE(ONE_D%DX_OLD(ONE_D%N_CELLS_OLD))
IF (ALLOCATED(ONE_D%TMP)) DEALLOCATE(ONE_D%TMP)                               ; ALLOCATE(ONE_D%TMP(0:ONE_D%N_CELLS_MAX+1))
IF (ALLOCATED(ONE_D%DELTA_TMP)) DEALLOCATE(ONE_D%DELTA_TMP)                   ; ALLOCATE(ONE_D%DELTA_TMP(0:ONE_D%N_CELLS_MAX+1))
IF (ALLOCATED(ONE_D%LAYER_THICKNESS)) DEALLOCATE(ONE_D%LAYER_THICKNESS)       ; ALLOCATE(ONE_D%LAYER_THICKNESS(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%LAYER_THICKNESS_OLD)) DEALLOCATE(ONE_D%LAYER_THICKNESS_OLD) 
   ALLOCATE(ONE_D%LAYER_THICKNESS_OLD(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%MIN_LAYER_THICKNESS)) DEALLOCATE(ONE_D%MIN_LAYER_THICKNESS)       
   ALLOCATE(ONE_D%MIN_LAYER_THICKNESS(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%HT3D_LAYER)) DEALLOCATE(ONE_D%HT3D_LAYER)                 ; ALLOCATE(ONE_D%HT3D_LAYER(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%MIN_DIFFUSIVITY)) DEALLOCATE(ONE_D%MIN_DIFFUSIVITY)       ; ALLOCATE(ONE_D%MIN_DIFFUSIVITY(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%RHO_C_S)) DEALLOCATE(ONE_D%RHO_C_S)                       ; ALLOCATE(ONE_D%RHO_C_S(ONE_D%N_CELLS_MAX))
IF (ALLOCATED(ONE_D%K_S)) DEALLOCATE(ONE_D%K_S)                               ; ALLOCATE(ONE_D%K_S(0:ONE_D%N_CELLS_MAX+1))
IF (ALLOCATED(ONE_D%DDSUM)) DEALLOCATE(ONE_D%DDSUM)                           ; ALLOCATE(ONE_D%DDSUM(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%SMALLEST_CELL_SIZE)) DEALLOCATE(ONE_D%SMALLEST_CELL_SIZE) ; ALLOCATE(ONE_D%SMALLEST_CELL_SIZE(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%PART_MASS)) DEALLOCATE(ONE_D%PART_MASS)                   ; ALLOCATE(ONE_D%PART_MASS(ONE_D%N_LPC))
IF (ALLOCATED(ONE_D%PART_ENTHALPY)) DEALLOCATE(ONE_D%PART_ENTHALPY)           ; ALLOCATE(ONE_D%PART_ENTHALPY(ONE_D%N_LPC))
IF (ALLOCATED(ONE_D%MATL_COMP)) DEALLOCATE(ONE_D%MATL_COMP)                   ; ALLOCATE(ONE_D%MATL_COMP(ONE_D%N_MATL))
IF (ALLOCATED(ONE_D%STRETCH_FACTOR)) DEALLOCATE(ONE_D%STRETCH_FACTOR)         ; ALLOCATE(ONE_D%STRETCH_FACTOR(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%HEAT_SOURCE)) DEALLOCATE(ONE_D%HEAT_SOURCE)               ; ALLOCATE(ONE_D%HEAT_SOURCE(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%CELL_SIZE_FACTOR)) DEALLOCATE(ONE_D%CELL_SIZE_FACTOR)     ; ALLOCATE(ONE_D%CELL_SIZE_FACTOR(ONE_D%N_LAYERS))
IF (ALLOCATED(ONE_D%CELL_SIZE)) DEALLOCATE(ONE_D%CELL_SIZE)                   ; ALLOCATE(ONE_D%CELL_SIZE(ONE_D%N_LAYERS))

DO NN=1,ONE_D%N_MATL
   IF (ALLOCATED(ONE_D%MATL_COMP(NN)%MASS_FRACTION)) DEALLOCATE(ONE_D%MATL_COMP(NN)%MASS_FRACTION)
   ALLOCATE(ONE_D%MATL_COMP(NN)%MASS_FRACTION(ONE_D%N_LAYERS))
   IF (ALLOCATED(ONE_D%MATL_COMP(NN)%RHO)) DEALLOCATE(ONE_D%MATL_COMP(NN)%RHO)
   ALLOCATE(ONE_D%MATL_COMP(NN)%RHO(0:ONE_D%N_CELLS_MAX+1))
ENDDO

END SUBROUTINE REALLOCATE_BOUNDARY_ONE_D


SUBROUTINE INITIALIZE_BOUNDARY_ONE_D(NM,OD_INDEX,SURF_INDEX)

USE GLOBAL_CONSTANTS, ONLY: RADIATION
INTEGER, INTENT(IN) :: NM,OD_INDEX,SURF_INDEX
TYPE(BOUNDARY_ONE_D_TYPE), POINTER :: ONE_D
TYPE(SURFACE_TYPE), POINTER :: SF
INTEGER :: NN,I,II
REAL(EB) :: RAMP_POSITION,TF,TB

ONE_D => MESHES(NM)%BOUNDARY_ONE_D(OD_INDEX)
SF => SURFACE(SURF_INDEX)

ONE_D%N_LAYER_CELLS(1:ONE_D%N_LAYERS) = SF%N_LAYER_CELLS(1:SF%N_LAYERS)
IF (ALLOCATED(ONE_D%REMESH_NWP)) ONE_D%REMESH_NWP(1:ONE_D%N_LAYERS) = SF%N_LAYER_CELLS(1:SF%N_LAYERS)
ONE_D%N_LAYER_CELLS_MAX(1:ONE_D%N_LAYERS) = SF%N_LAYER_CELLS_MAX(1:SF%N_LAYERS)
ONE_D%RAMP_IHS_INDEX(1:ONE_D%N_LAYERS) = SF%RAMP_IHS_INDEX(1:SF%N_LAYERS)
IF (ALLOCATED(ONE_D%MATL_INDEX) .AND. ALLOCATED(SF%MATL_INDEX)) ONE_D%MATL_INDEX(1:ONE_D%N_MATL) = SF%MATL_INDEX(1:SF%N_MATL)
ONE_D%M_DOT_S_PP = 0._EB
ONE_D%X=0._EB ; ONE_D%X(0:ONE_D%N_CELLS_INI) = SF%X_S(0:ONE_D%N_CELLS_INI)
ONE_D%DX_OLD=0._EB
DO I=1,MIN(ONE_D%N_CELLS_OLD,ONE_D%N_CELLS_INI)
   ONE_D%DX_OLD(I) = ONE_D%X(I)-ONE_D%X(I-1)
ENDDO
IF (ALLOCATED(ONE_D%LAYER_THICKNESS_OLD)) ONE_D%LAYER_THICKNESS_OLD(1:ONE_D%N_LAYERS) = SF%LAYER_THICKNESS(1:SF%N_LAYERS)
IF (SF%RAMP_T_I_INDEX > 0) THEN
   !NOTE: Replicating EVALUATE_RAMP since MODULE MATH_FUNCTIONS uses the MODULE which contains this routine
   DO I=1,ONE_D%N_CELLS_INI
      RAMP_POSITION = &
           MAX(0._EB,MIN(RAMPS(SF%RAMP_T_I_INDEX)%SPAN,ONE_D%X(I) - RAMPS(SF%RAMP_T_I_INDEX)%T_MIN))
      II = MIN(UBOUND(RAMPS(SF%RAMP_T_I_INDEX)%INTERPOLATED_DATA,1),&
           MAX(LBOUND(RAMPS(SF%RAMP_T_I_INDEX)%INTERPOLATED_DATA,1),NINT(RAMP_POSITION*RAMPS(SF%RAMP_T_I_INDEX)%RDT)))
      ONE_D%TMP(I) = RAMPS(SF%RAMP_T_I_INDEX)%INTERPOLATED_DATA(II)
   ENDDO
   RAMP_POSITION = &
         MAX(0._EB,MIN(RAMPS(SF%RAMP_T_I_INDEX)%SPAN,0._EB - RAMPS(SF%RAMP_T_I_INDEX)%T_MIN))
   II = MIN(UBOUND(RAMPS(SF%RAMP_T_I_INDEX)%INTERPOLATED_DATA,1),&
        MAX(LBOUND(RAMPS(SF%RAMP_T_I_INDEX)%INTERPOLATED_DATA,1),NINT(RAMP_POSITION*RAMPS(SF%RAMP_T_I_INDEX)%RDT)))
   TF = RAMPS(SF%RAMP_T_I_INDEX)%INTERPOLATED_DATA(II)
   ONE_D%TMP(0) = TF * 2._EB - ONE_D%TMP(1)
   RAMP_POSITION = &
         MAX(0._EB,MIN(RAMPS(SF%RAMP_T_I_INDEX)%SPAN,SUM(SF%LAYER_THICKNESS) - RAMPS(SF%RAMP_T_I_INDEX)%T_MIN))
   II = MIN(UBOUND(RAMPS(SF%RAMP_T_I_INDEX)%INTERPOLATED_DATA,1),&
        MAX(LBOUND(RAMPS(SF%RAMP_T_I_INDEX)%INTERPOLATED_DATA,1),NINT(RAMP_POSITION*RAMPS(SF%RAMP_T_I_INDEX)%RDT)))
   TB = RAMPS(SF%RAMP_T_I_INDEX)%INTERPOLATED_DATA(II)
   ONE_D%TMP(ONE_D%N_CELLS_INI+1) = TB * 2._EB - ONE_D%TMP(ONE_D%N_CELLS_INI)
ELSE
   ONE_D%TMP = SF%TMP_INNER
ENDIF
ONE_D%DELTA_TMP = 0._EB
ONE_D%LAYER_THICKNESS(1:ONE_D%N_LAYERS) = SF%LAYER_THICKNESS(1:SF%N_LAYERS)
ONE_D%MIN_LAYER_THICKNESS(1:ONE_D%N_LAYERS) = SF%MIN_LAYER_THICKNESS(1:SF%N_LAYERS)
ONE_D%HT3D_LAYER(1:ONE_D%N_LAYERS) = SF%HT3D_LAYER(1:SF%N_LAYERS)
ONE_D%MIN_DIFFUSIVITY(1:ONE_D%N_LAYERS) = SF%MIN_DIFFUSIVITY(1:SF%N_LAYERS)
ONE_D%STRETCH_FACTOR(1:ONE_D%N_LAYERS)  = SF%STRETCH_FACTOR(1:SF%N_LAYERS)
ONE_D%HEAT_SOURCE(1:ONE_D%N_LAYERS)  = SF%HEAT_SOURCE(1:SF%N_LAYERS)
ONE_D%CELL_SIZE(1:ONE_D%N_LAYERS)  = SF%CELL_SIZE(1:SF%N_LAYERS)
ONE_D%CELL_SIZE_FACTOR(1:ONE_D%N_LAYERS)  = SF%CELL_SIZE_FACTOR(1:SF%N_LAYERS)
ONE_D%RHO_C_S = 1.E6_EB
ONE_D%K_S = 0._EB
ONE_D%DDSUM(1:SF%N_LAYERS) = SF%DDSUM(1:SF%N_LAYERS)
ONE_D%SMALLEST_CELL_SIZE(1:SF%N_LAYERS) = SF%SMALLEST_CELL_SIZE(1:SF%N_LAYERS)
ONE_D%PART_MASS = 0._EB
ONE_D%PART_ENTHALPY = 0._EB
DO NN=1,ONE_D%N_MATL
   ONE_D%MATL_COMP(NN)%RHO(0:ONE_D%N_CELLS_INI+1) = SF%RHO_0(0:ONE_D%N_CELLS_INI+1,NN)
   IF (RADIATION .AND. MATERIAL(ONE_D%MATL_INDEX(NN))%KAPPA_S<4.9E4_EB) ONE_D%INTERNAL_RADIATION = .TRUE.
ENDDO
IF (RADIATION .AND. ANY(SF%KAPPA_S(1:SF%N_LAYERS)>0._EB)) ONE_D%INTERNAL_RADIATION = .TRUE.
ONE_D%BACK_MESH = NM
ONE_D%BACK_SURF = SURF_INDEX
ONE_D%PYROLYSIS_MODEL = SF%PYROLYSIS_MODEL

END SUBROUTINE INITIALIZE_BOUNDARY_ONE_D


SUBROUTINE INITIALIZE_BOUNDARY_PROP1(NM,B1_INDEX,OD_INDEX,SURF_INDEX)

USE PHYSICAL_FUNCTIONS, ONLY: GET_EMISSIVITY
USE GLOBAL_CONSTANTS, ONLY: SIGMA,TMPA,RHOA,TMPA4,N_TRACKED_SPECIES, RADIATION
INTEGER, INTENT(IN) :: NM,B1_INDEX,SURF_INDEX,OD_INDEX
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(BOUNDARY_ONE_D_TYPE), POINTER :: ONE_D
TYPE(SURFACE_TYPE), POINTER :: SF

B1 => MESHES(NM)%BOUNDARY_PROP1(B1_INDEX)
SF => SURFACE(SURF_INDEX)

IF (OD_INDEX==0 .OR. SF%EMISSIVITY_SPECIFIED) THEN
   B1%EMISSIVITY      = SF%EMISSIVITY
ELSE
   ONE_D => MESHES(NM)%BOUNDARY_ONE_D(OD_INDEX)
   CALL GET_EMISSIVITY(ONE_D,1,B1%EMISSIVITY)
ENDIF
IF (RADIATION) THEN
   B1%Q_RAD_IN        = B1%EMISSIVITY*SIGMA*TMPA4
   B1%Q_RAD_OUT       = B1%EMISSIVITY*SIGMA*TMPA4
ELSE
   B1%Q_RAD_IN        = 0._EB
   B1%Q_RAD_OUT       = 0._EB
ENDIF
B1%T_IGN           = SF%T_IGN
B1%TMP_F           = SF%TMP_FRONT
B1%TMP_G           = TMPA
B1%TMP_B           = SF%TMP_BACK
B1%RHO_F           = RHOA
B1%RHO_G           = RHOA
B1%TMP_F_OLD       = SF%TMP_FRONT
B1%BURN_DURATION   = SF%BURN_DURATION

B1%M_DOT_G_PP_ACTUAL = 0._EB
B1%M_DOT_G_PP_ADJUST = 0._EB
B1%ZZ_F(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
B1%ZZ_G(1:N_TRACKED_SPECIES) = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
B1%RHO_D_F = 0._EB
B1%RHO_D_DZDN_F = 0._EB
B1%AWM_AEROSOL = 0._EB
B1%Q_IN_SMOOTH = 0._EB
B1%QDOTPP_INT= 0._EB
B1%Q_IN_SMOOTH_INT = 0._EB

END SUBROUTINE INITIALIZE_BOUNDARY_PROP1


!> \brief Pack BOUNDARY_PROP1 into REAL, INTEGER, and LOGICAL 1-D arrays for a WALL, CFACE, or LAGRANGIAN_PARTICLE
!> \param NM Mesh index
!> \param IC Integer counter
!> \param RC Real counter
!> \param LC Logical counter
!> \param OS STORAGE_ARRAY name
!> \param B1_INDEX Index of the BOUNDARY_PROP1 array
!> \param UNPACK_IT Flag indicating whether the data is to be packed into the 1-D array or unpacked from it
!> \param COUNT_ONLY Flag signifying that only a variable count is to be done; no packing
!> \param SURF_INDEX Index of the SURFACE type

SUBROUTINE PACK_BOUNDARY_PROP1(NM,IC,RC,LC,OS,B1_INDEX,UNPACK_IT,COUNT_ONLY,SURF_INDEX)

USE COMP_OPERATORS, ONLY: EQUATE
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,N_SURFACE_DENSITY_SPECIES
INTEGER, INTENT(IN) :: NM,B1_INDEX,SURF_INDEX
INTEGER, INTENT(INOUT) :: IC,RC,LC
LOGICAL, INTENT(IN) :: UNPACK_IT,COUNT_ONLY
INTEGER :: I1,I2
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(STORAGE_TYPE), POINTER :: OS

IF (.NOT.COUNT_ONLY) B1 => MESHES(NM)%BOUNDARY_PROP1(B1_INDEX)

! Integers

IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),B1%SURF_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),B1%PRESSURE_ZONE,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),B1%NODE_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),B1%N_SUBSTEPS,UNPACK_IT)

! Logicals

LC=LC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%LOGICALS(LC),B1%BURNAWAY,UNPACK_IT)
LC=LC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%LOGICALS(LC),B1%LAYER_REMOVED,UNPACK_IT)

! Reals

RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%AREA,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%HEAT_TRANS_COEF,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%Q_CON_F,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%Q_RAD_IN,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%Q_RAD_OUT,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%EMISSIVITY,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%AREA_ADJUST,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%T_IGN,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%TMP_F,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%TMP_G,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%TMP_B,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%U_NORMAL,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%U_NORMAL_S,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%U_NORMAL_0,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%U_TANG,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%U_IMPACT,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%RHO_F,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%RHO_G,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%RDN,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%K_G,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%Q_DOT_G_PP,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%Q_DOT_O2_PP,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%Q_CONDENSE,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%TMP_F_OLD,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%BURN_DURATION,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%Q_IN_SMOOTH,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%T_MATL_PART,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%B_NUMBER,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%M_DOT_PART_ACTUAL,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%Q_LEAK,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B1%VEL_ERR_NEW,UNPACK_IT)

I2 = RC  ! I1 and I2 continue counting reals

I1 = I2+1 ; I2 = I1 + N_TRACKED_SPECIES - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2) , B1%M_DOT_G_PP_ACTUAL(1:I2-I1+1) , UNPACK_IT)

I1 = I2+1 ; I2 = I1 + N_TRACKED_SPECIES - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2) , B1%M_DOT_G_PP_ADJUST(1:I2-I1+1) , UNPACK_IT)

I1 = I2+1 ; I2 = I1 + N_TRACKED_SPECIES - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2) , B1%ZZ_F(1:I2-I1+1) , UNPACK_IT)

I1 = I2+1 ; I2 = I1 + N_TRACKED_SPECIES - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2) , B1%ZZ_G(1:I2-I1+1) , UNPACK_IT)

I1 = I2+1 ; I2 = I1 + N_TRACKED_SPECIES - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2) , B1%RHO_D_F(1:I2-I1+1) , UNPACK_IT)

I1 = I2+1 ; I2 = I1 + N_TRACKED_SPECIES - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2) , B1%RHO_D_DZDN_F(1:I2-I1+1) , UNPACK_IT)

I1 = I2+1 ; I2 = I1 + N_SURFACE_DENSITY_SPECIES - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2) , B1%AWM_AEROSOL(1:I2-I1+1) , UNPACK_IT)

I1 = I2+1 ; I2 = I1 + SURFACE(SURF_INDEX)%N_QDOTPP_REF - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2) , B1%QDOTPP_INT(1:I2-I1+1) , UNPACK_IT)

I1 = I2+1 ; I2 = I1 + SURFACE(SURF_INDEX)%N_THICK_REF - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2) , B1%Q_IN_SMOOTH_INT(1:I2-I1+1) , UNPACK_IT)


RC = I2

END SUBROUTINE PACK_BOUNDARY_PROP1


!> \brief Pack or unpack BOUNDARY_PROP2 components into 1-D arrays
!> \param NM Mesh number
!> \param IC Integer Counter
!> \param RC Real Counter
!> \param OS Storage array
!> \param B2_INDEX Index within BOUNDARY_PROP2 array
!> \param UNPACK_IT Flag indicating whether to pack or unpack
!> \param COUNT_ONLY Flag signifying that only a variable count is to be done; no packing

SUBROUTINE PACK_BOUNDARY_PROP2(NM,IC,RC,OS,B2_INDEX,UNPACK_IT,COUNT_ONLY)

USE COMP_OPERATORS, ONLY: EQUATE
USE GLOBAL_CONSTANTS, ONLY: N_LP_ARRAY_INDICES
INTEGER, INTENT(IN) :: NM,B2_INDEX
INTEGER, INTENT(INOUT) :: RC,IC
LOGICAL, INTENT(IN) :: UNPACK_IT,COUNT_ONLY
INTEGER :: I1,I2
TYPE(BOUNDARY_PROP2_TYPE), POINTER :: B2
TYPE(STORAGE_TYPE), POINTER :: OS

IF (.NOT.COUNT_ONLY) B2 => MESHES(NM)%BOUNDARY_PROP2(B2_INDEX)

IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),B2%SURF_INDEX,UNPACK_IT)
IC=IC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%INTEGERS(IC),B2%HEAT_TRANSFER_REGIME,UNPACK_IT)

RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B2%U_TAU,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B2%Y_PLUS,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B2%Z_STAR,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B2%PHI_LS,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B2%WORK1,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B2%WORK2,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B2%WORK3,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B2%K_SUPPRESSION,UNPACK_IT)
RC=RC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(RC),B2%V_DEP,UNPACK_IT)

I2 = RC

I1 = I2+1 ; I2 = I1 + N_LP_ARRAY_INDICES - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2),B2%A_LP_MPUA(1:I2-I1+1),UNPACK_IT)

I1 = I2+1 ; I2 = I1 + N_LP_ARRAY_INDICES - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2),B2%LP_CPUA(1:I2-I1+1),UNPACK_IT)

I1 = I2+1 ; I2 = I1 + N_LP_ARRAY_INDICES - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2),B2%LP_EMPUA(1:I2-I1+1),UNPACK_IT)

I1 = I2+1 ; I2 = I1 + N_LP_ARRAY_INDICES - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2),B2%LP_MPUA(1:I2-I1+1),UNPACK_IT)

I1 = I2+1 ; I2 = I1 + N_LP_ARRAY_INDICES - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2),B2%LP_TEMP(1:I2-I1+1),UNPACK_IT)

RC = I2

END SUBROUTINE PACK_BOUNDARY_PROP2


!> \brief Pack or unpack BOUNDARY_RADIA components into 1-D arrays
!> \param NM Mesh number
!> \param RC Real Counter
!> \param OS Storage array
!> \param BR_INDEX Index within BOUNDARY_RADIA array
!> \param UNPACK_IT Flag indicating whether to pack or unpack
!> \param COUNT_ONLY Flag signifying that only a variable count is to be done; no packing

SUBROUTINE PACK_BOUNDARY_RADIA(NM,RC,OS,BR_INDEX,UNPACK_IT,COUNT_ONLY)

USE COMP_OPERATORS, ONLY: EQUATE
USE GLOBAL_CONSTANTS, ONLY: NUMBER_SPECTRAL_BANDS,NUMBER_RADIATION_ANGLES
INTEGER, INTENT(IN) :: NM,BR_INDEX
INTEGER, INTENT(INOUT) :: RC
LOGICAL, INTENT(IN) :: UNPACK_IT,COUNT_ONLY
INTEGER :: I1,I2,NN
TYPE(BOUNDARY_RADIA_TYPE), POINTER :: BR
TYPE(STORAGE_TYPE), POINTER :: OS

IF (.NOT.COUNT_ONLY) BR => MESHES(NM)%BOUNDARY_RADIA(BR_INDEX)

I2 = RC

DO NN=1,NUMBER_SPECTRAL_BANDS
   I1 = I2+1 ; I2 = I1 + NUMBER_RADIATION_ANGLES - 1
   IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2),BR%BAND(NN)%ILW(1:I2-I1+1),UNPACK_IT)
ENDDO

I1 = I2+1 ; I2 = I1 + NUMBER_SPECTRAL_BANDS - 1
IF (.NOT.COUNT_ONLY) CALL EQUATE(OS%REALS(I1:I2),BR%IL(1:I2-I1+1),UNPACK_IT)

RC = I2

END SUBROUTINE PACK_BOUNDARY_RADIA


!> \brief Pack or unpack CELL
!> \param NM Mesh number
!> \param UNPACK_IT Logical indicating whether to pack or unpack
!> \param N_INTEGERS Number of integers packed
!> \param N_LOGICALS Number of logicals packed

SUBROUTINE PACK_CELL(NM,UNPACK_IT,N_INTEGERS,N_LOGICALS)

USE COMP_OPERATORS, ONLY: EQUATE
USE GLOBAL_CONSTANTS, ONLY: CELL_COUNT
INTEGER, INTENT(IN) :: NM
LOGICAL, INTENT(IN) :: UNPACK_IT
INTEGER, INTENT(OUT), OPTIONAL :: N_INTEGERS,N_LOGICALS
INTEGER :: IC,LC,ICC,N
LOGICAL :: COUNT_ONLY
TYPE(MESH_TYPE), POINTER :: M

M => MESHES(NM)

IF (PRESENT(N_INTEGERS) .OR. PRESENT(N_LOGICALS)) THEN
   COUNT_ONLY = .TRUE.
ELSE
   COUNT_ONLY = .FALSE.
ENDIF

IC=0 ; LC=0

DO ICC=1,CELL_COUNT(NM)

   LC=LC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(M%CELL_LOGICALS(LC),M%CELL(ICC)%SOLID,UNPACK_IT)
   LC=LC+1 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(M%CELL_LOGICALS(LC),M%CELL(ICC)%EXTERIOR,UNPACK_IT)

   IC=IC+1  ; IF (.NOT.COUNT_ONLY) CALL EQUATE(M%CELL_INTEGERS(IC)      ,M%CELL(ICC)%OBST_INDEX,UNPACK_IT)
   IC=IC+1  ; IF (.NOT.COUNT_ONLY) CALL EQUATE(M%CELL_INTEGERS(IC)      ,M%CELL(ICC)%I,UNPACK_IT)
   IC=IC+1  ; IF (.NOT.COUNT_ONLY) CALL EQUATE(M%CELL_INTEGERS(IC)      ,M%CELL(ICC)%J,UNPACK_IT)
   IC=IC+1  ; IF (.NOT.COUNT_ONLY) CALL EQUATE(M%CELL_INTEGERS(IC)      ,M%CELL(ICC)%K,UNPACK_IT)
   IC=IC+12 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(M%CELL_INTEGERS(IC-11:IC),M%CELL(ICC)%EDGE_INDEX(1:12),UNPACK_IT)
   IC=IC+7  ; IF (.NOT.COUNT_ONLY) CALL EQUATE(M%CELL_INTEGERS(IC- 6:IC),M%CELL(ICC)%WALL_INDEX(-3:3),UNPACK_IT)
   IC=IC+7  ; IF (.NOT.COUNT_ONLY) CALL EQUATE(M%CELL_INTEGERS(IC- 6:IC),M%CELL(ICC)%SURF_INDEX(-3:3),UNPACK_IT)
   DO N=1,3
      IC=IC+7 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(M%CELL_INTEGERS(IC-6:IC),M%CELL(ICC)%THIN_WALL_INDEX(-3:3,N),UNPACK_IT)
      IC=IC+7 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(M%CELL_INTEGERS(IC-6:IC),M%CELL(ICC)%THIN_SURF_INDEX(-3:3,N),UNPACK_IT)
      IC=IC+7 ; IF (.NOT.COUNT_ONLY) CALL EQUATE(M%CELL_INTEGERS(IC-6:IC),M%CELL(ICC)%THIN_OBST_INDEX(-3:3,N),UNPACK_IT)
   ENDDO

ENDDO

! If this is a COUNT_ONLY, return the total number of variables stored

IF (PRESENT(N_INTEGERS)) N_INTEGERS = IC
IF (PRESENT(N_LOGICALS)) N_LOGICALS = LC

END SUBROUTINE PACK_CELL


!> \brief Make a particle suitable for reuse
!> \param NM Mesh number
!> \param LP_INDEX Particle index

SUBROUTINE NULLIFY_PARTICLE(NM,LP_INDEX)

INTEGER, INTENT(IN) :: NM,LP_INDEX
TYPE(MESH_TYPE), POINTER :: M
TYPE(LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP

M => MESHES(NM)
LP => M%LAGRANGIAN_PARTICLE(LP_INDEX)

IF (LP%BC_INDEX>0) THEN
   M%BOUNDARY_COORD_OCCUPANCY(LP%BC_INDEX) = 0
   M%NEXT_AVAILABLE_BOUNDARY_COORD_SLOT = MIN(LP%BC_INDEX,M%NEXT_AVAILABLE_BOUNDARY_COORD_SLOT)
   LP%BC_INDEX = 0
ENDIF
IF (LP%OD_INDEX>0) THEN
   M%BOUNDARY_ONE_D_OCCUPANCY(LP%OD_INDEX) = 0
   M%NEXT_AVAILABLE_BOUNDARY_ONE_D_SLOT = MIN(LP%OD_INDEX,M%NEXT_AVAILABLE_BOUNDARY_ONE_D_SLOT)
   LP%OD_INDEX = 0
ENDIF
IF (LP%B1_INDEX>0) THEN
   M%BOUNDARY_PROP1_OCCUPANCY(LP%B1_INDEX) = 0
   M%NEXT_AVAILABLE_BOUNDARY_PROP1_SLOT = MIN(LP%B1_INDEX,M%NEXT_AVAILABLE_BOUNDARY_PROP1_SLOT)
   LP%B1_INDEX = 0
ENDIF
IF (LP%B2_INDEX>0) THEN
   M%BOUNDARY_PROP2_OCCUPANCY(LP%B2_INDEX) = 0
   M%NEXT_AVAILABLE_BOUNDARY_PROP2_SLOT = MIN(LP%B2_INDEX,M%NEXT_AVAILABLE_BOUNDARY_PROP2_SLOT)
   LP%B2_INDEX = 0
ENDIF
IF (LP%BR_INDEX>0) THEN
   M%BOUNDARY_RADIA_OCCUPANCY(LP%BR_INDEX) = 0
   M%NEXT_AVAILABLE_BOUNDARY_RADIA_SLOT = MIN(LP%BR_INDEX,M%NEXT_AVAILABLE_BOUNDARY_RADIA_SLOT)
   LP%BR_INDEX = 0
ENDIF

END SUBROUTINE NULLIFY_PARTICLE


!> \brief Increase the size of the storage arrays used to hold WALL, CFACE, LAGRANGIAN_PARTICLE data
!> \param OS Storage variable
!> \param N_REALS Current dimension of OS%REALS
!> \param N_INTEGERS Current dimension of OS%INTEGERS
!> \param N_LOGICALS Current dimension of OS%LOGICALS
!> \param N_REALS_NEW New dimension of OS%REALS
!> \param N_INTEGERS_NEW New dimension of OS%INTEGERS
!> \param N_LOGICALS_NEW New dimension of OS%LOGICALS

SUBROUTINE REALLOCATE_STORAGE_ARRAYS(OS,N_REALS,N_INTEGERS,N_LOGICALS,N_REALS_NEW,N_INTEGERS_NEW,N_LOGICALS_NEW)

INTEGER, OPTIONAL, INTENT(IN) :: N_REALS,N_INTEGERS,N_LOGICALS,N_REALS_NEW,N_INTEGERS_NEW,N_LOGICALS_NEW
TYPE (STORAGE_TYPE), POINTER :: OS
REAL(EB), ALLOCATABLE, DIMENSION(:) :: DUMMY_REALS
INTEGER, ALLOCATABLE, DIMENSION(:) :: DUMMY_INTEGERS
LOGICAL, ALLOCATABLE, DIMENSION(:) :: DUMMY_LOGICALS

IF (PRESENT(N_REALS) .AND. ALLOCATED(OS%REALS)) THEN
   ALLOCATE(DUMMY_REALS(N_REALS_NEW))
   DUMMY_REALS(1:N_REALS) = OS%REALS(1:N_REALS)
   CALL MOVE_ALLOC(DUMMY_REALS,OS%REALS)
ELSEIF (PRESENT(N_REALS)) THEN
   ALLOCATE(OS%REALS(N_REALS_NEW))
ENDIF

IF (PRESENT(N_INTEGERS) .AND. ALLOCATED(OS%INTEGERS)) THEN
   ALLOCATE(DUMMY_INTEGERS(N_INTEGERS_NEW))
   DUMMY_INTEGERS(1:N_INTEGERS) = OS%INTEGERS(1:N_INTEGERS)
   CALL MOVE_ALLOC(DUMMY_INTEGERS,OS%INTEGERS)
ELSEIF (PRESENT(N_INTEGERS)) THEN
   ALLOCATE(OS%INTEGERS(N_INTEGERS_NEW))
ENDIF

IF (PRESENT(N_LOGICALS) .AND. ALLOCATED(OS%LOGICALS)) THEN
   ALLOCATE(DUMMY_LOGICALS(N_LOGICALS_NEW))
   DUMMY_LOGICALS(1:N_LOGICALS) = OS%LOGICALS(1:N_LOGICALS)
   CALL MOVE_ALLOC(DUMMY_LOGICALS,OS%LOGICALS)
ELSEIF (PRESENT(N_LOGICALS)) THEN
   ALLOCATE(OS%LOGICALS(N_LOGICALS_NEW))
ENDIF

END SUBROUTINE REALLOCATE_STORAGE_ARRAYS


!> \brief Find the index of the LAGRANGIAN_PARTICLE with the given TAG. If not found, return 0.
!> \param NM Mesh number
!> \param LP_TAG A unique integer assigned to each LAGRANGIAN_PARTICLE, regardless of mesh
!> \param LP_INDEX The index of the LAGRANGIAN_PARTICLE of the input class and tag

SUBROUTINE GET_LAGRANGIAN_PARTICLE_INDEX(NM,LP_TAG,LP_INDEX)

INTEGER, INTENT(IN) :: NM,LP_TAG
INTEGER, INTENT(OUT) :: LP_INDEX
INTEGER :: IP

LP_INDEX = 0
DO IP=1,MESHES(NM)%NLP
   IF (MESHES(NM)%LAGRANGIAN_PARTICLE(IP)%TAG==LP_TAG) THEN
      LP_INDEX = IP
      EXIT
   ENDIF
ENDDO

END SUBROUTINE GET_LAGRANGIAN_PARTICLE_INDEX


!> \brief Broadcast basic geometry info from each mesh, NM, to its neighbors in MPI_COMM_NEIGHBORS(NM)

SUBROUTINE EXCHANGE_GEOMETRY_INFO

USE GLOBAL_CONSTANTS, ONLY: N_MPI_PROCESSES,MY_RANK,LOWER_MESH_INDEX,UPPER_MESH_INDEX,PROCESS,NMESHES,MPI_COMM_NEIGHBORS,&
                            MPI_COMM_NEIGHBORS_ROOT,CELL_COUNT,CELL_COUNT_INTEGERS,CELL_COUNT_LOGICALS,&
                            COUNTS,DISPLS
#ifdef WITHOUT_MPIF08
USE MPI
#else
USE MPI_F08
#endif
INTEGER :: NM,NOM,IERR,SNODE,N
TYPE(MESH_TYPE), POINTER ::  M,M4

IF (N_MPI_PROCESSES==1) RETURN

! Ensure that all MPI processes have consistent bounds for the CELL and related arrays

CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,CELL_COUNT(1:NMESHES),COUNTS,DISPLS,MPI_INTEGER,MPI_COMM_WORLD,IERR)
CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,CELL_COUNT_INTEGERS(1:NMESHES),COUNTS,DISPLS,MPI_INTEGER,MPI_COMM_WORLD,IERR)
CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,CELL_COUNT_LOGICALS(1:NMESHES),COUNTS,DISPLS,MPI_INTEGER,MPI_COMM_WORLD,IERR)

! Copy the real, integer, and logical components of CELL into 1-D arrays for the purpose of MPI broadcast

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL PACK_CELL(NM,.FALSE.)
ENDDO

! Allocate CELL and related arrays in neighboring meshes

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   M => MESHES(NM)
   DO N=1,M%N_NEIGHBORING_MESHES
      NOM = M%NEIGHBORING_MESH(N)
      SNODE = PROCESS(NOM)
      IF (SNODE==MY_RANK) CYCLE
      M4 => MESHES(NOM)
      IF (.NOT.ALLOCATED(M4%CELL_INDEX))    ALLOCATE(M4%CELL_INDEX(0:M4%IBP1,0:M4%JBP1,0:M4%KBP1))
      CALL REALLOCATE_CELL(NOM,UBOUND(M4%CELL,DIM=1),CELL_COUNT(NOM))
   ENDDO
ENDDO

! Broadcast CELL to other MPI processes

DO NM=1,NMESHES
   IF (MPI_COMM_NEIGHBORS(NM)==MPI_COMM_NULL) CYCLE
   M => MESHES(NM)
   CALL MPI_BCAST(M%CELL_INDEX(0,0,0),SIZE(M%CELL_INDEX),MPI_INTEGER,MPI_COMM_NEIGHBORS_ROOT(NM),MPI_COMM_NEIGHBORS(NM),IERR)
   CALL MPI_BCAST(M%CELL_INTEGERS(1),SIZE(M%CELL_INTEGERS),MPI_INTEGER,MPI_COMM_NEIGHBORS_ROOT(NM),MPI_COMM_NEIGHBORS(NM),IERR)
   CALL MPI_BCAST(M%CELL_LOGICALS(1),SIZE(M%CELL_LOGICALS),MPI_LOGICAL,MPI_COMM_NEIGHBORS_ROOT(NM),MPI_COMM_NEIGHBORS(NM),IERR)
ENDDO

! Copy CELL components out of the 1-D arrays

DO NM=1,NMESHES
   IF (PROCESS(NM)==MY_RANK) CYCLE
   IF (MPI_COMM_NEIGHBORS(NM)==MPI_COMM_NULL) CYCLE
   CALL PACK_CELL(NM,.TRUE.)
ENDDO

END SUBROUTINE EXCHANGE_GEOMETRY_INFO

END MODULE MEMORY_FUNCTIONS



!> \brief Functions and subroutines for manipulating geometry

MODULE GEOMETRY_FUNCTIONS

USE PRECISION_PARAMETERS
USE MESH_VARIABLES
USE GLOBAL_CONSTANTS
IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS


!> \brief Determine the mesh index and cell indices of an input point.
!> \param XX x-coordinate of the point (m)
!> \param YY y-coordinate of the point (m)
!> \param ZZ z-coordinate of the point (m)
!> \param NOM Number of the other mesh where the point is located
!> \param IIO I-index of the cell in the other mesh where the point is located
!> \param JJO J-index of the cell in the other mesh where the point is located
!> \param KKO K-index of the cell in the other mesh where the point is located
!> \param XXI (Optional) continuous form of I index
!> \param YYJ (Optional) continuous form of J index
!> \param ZZK (Optional) continuous form of K index

SUBROUTINE SEARCH_OTHER_MESHES(XX,YY,ZZ,NOM,IIO,JJO,KKO,XXI,YYJ,ZZK)

REAL(EB), INTENT(IN) :: XX,YY,ZZ
REAL(EB), OPTIONAL :: XXI,YYJ,ZZK
REAL(EB) :: XI,YJ,ZK
INTEGER, INTENT(OUT) :: NOM,IIO,JJO,KKO
TYPE (MESH_TYPE), POINTER :: M2=>NULL()

IF (PRESENT(XXI)) XXI = 0._EB
IF (PRESENT(YYJ)) YYJ = 0._EB
IF (PRESENT(ZZK)) ZZK = 0._EB

OTHER_MESH_LOOP: DO NOM=1,NMESHES
   M2=>MESHES(NOM)
   IF (XX>=M2%XS .AND. XX<=M2%XF .AND.  YY>=M2%YS .AND. YY<=M2%YF .AND. ZZ>=M2%ZS .AND. ZZ<=M2%ZF) THEN
      IF (ALLOCATED(M2%CELLSI)) THEN
         XI  = MAX( 1._EB , MIN( REAL(M2%IBP1,EB)*ONE_M_EPS , M2%CELLSI(FLOOR((XX-M2%XS)*M2%RDXINT)) + 1._EB ) )
         YJ  = MAX( 1._EB , MIN( REAL(M2%JBP1,EB)*ONE_M_EPS , M2%CELLSJ(FLOOR((YY-M2%YS)*M2%RDYINT)) + 1._EB ) )
         ZK  = MAX( 1._EB , MIN( REAL(M2%KBP1,EB)*ONE_M_EPS , M2%CELLSK(FLOOR((ZZ-M2%ZS)*M2%RDZINT)) + 1._EB ) )
         IIO = FLOOR(XI)
         JJO = FLOOR(YJ)
         KKO = FLOOR(ZK)
         IF (PRESENT(XXI)) XXI = XI - 1._EB
         IF (PRESENT(YYJ)) YYJ = YJ - 1._EB
         IF (PRESENT(ZZK)) ZZK = ZK - 1._EB
      ELSE
         IIO = 0  ! The mesh if found, but no detailed information is available to the current process
         JJO = 0
         KKO = 0
      ENDIF
      RETURN
   ENDIF
ENDDO OTHER_MESH_LOOP

NOM = 0

END SUBROUTINE SEARCH_OTHER_MESHES


!> \brief Block or unblock a cell of a given mesh
!> \param NM Mesh number
!> \param I1 Lower bound of I cell indices of the region to be blocked or unblocked
!> \param I2 Upper bound of I cell indices of the region to be blocked or unblocked
!> \param J1 Lower bound of J cell indices of the region to be blocked or unblocked
!> \param J2 Upper bound of J cell indices of the region to be blocked or unblocked
!> \param K1 Lower bound of K cell indices of the region to be blocked or unblocked
!> \param K2 Upper bound of K cell indices of the region to be blocked or unblocked
!> \param IVAL Indicator if the cell is to blocked (1) or unblocked (0)
!> \param OBST_INDEX Index of the OBSTruction blocking the cell

SUBROUTINE BLOCK_CELL(NM,I1,I2,J1,J2,K1,K2,IVAL,OBST_INDEX)

INTEGER :: NM,I1,I2,J1,J2,K1,K2,IVAL,I,J,K,OBST_INDEX,IC
TYPE (MESH_TYPE), POINTER :: M=>NULL()

M => MESHES(NM)
DO K=K1,K2
   DO J=J1,J2
      DO I=I1,I2
         IC = M%CELL_INDEX(I,J,K)
         SELECT CASE(IVAL)
            CASE(0)
               M%CELL(IC)%SOLID   = .FALSE.
               M%CELL(IC)%OBST_INDEX = 0
            CASE(1)
               M%CELL(IC)%SOLID   = .TRUE.
               M%CELL(IC)%OBST_INDEX = OBST_INDEX
         END SELECT
         IF (OBST_INDEX==0) M%CELL(IC)%EXTERIOR = .TRUE.
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE BLOCK_CELL


!> \brief Determine if a point is in the interior of at least one mesh
!> \param XX x-coordinate of the point (m)
!> \param YY y-coordinate of the point (m)
!> \param ZZ z-coordinate of the point (m)

LOGICAL FUNCTION INTERIOR(XX,YY,ZZ)

INTEGER NM
REAL(EB), INTENT(IN) :: XX,YY,ZZ

INTERIOR = .FALSE.

DO NM=1,NMESHES
   IF (XX>MESHES(NM)%XS .AND. XX<MESHES(NM)%XF .AND. &
       YY>MESHES(NM)%YS .AND. YY<MESHES(NM)%YF .AND. &
       ZZ>MESHES(NM)%ZS .AND. ZZ<MESHES(NM)%ZF) INTERIOR = .TRUE.
ENDDO

END FUNCTION INTERIOR


!> \brief Assign the pressure zone for all cells connected to a given input point
!> \param NM Mesh number
!> \param II i-coordinate of the point (m)
!> \param JJ j-coordinate of the point (m)
!> \param KK k-coordinate of the point (m)
!> \param I_ZONE Index of the pressure zone to assign to the connected cells
!> \param I_ZONE_OVERLAP Index of a pressure zone that overlaps I_ZONE

SUBROUTINE ASSIGN_PRESSURE_ZONE(NM,II,JJ,KK,I_ZONE,I_ZONE_OVERLAP)

USE COMP_FUNCTIONS, ONLY : SHUTDOWN
INTEGER, INTENT(IN) :: NM,I_ZONE,II,JJ,KK
INTEGER, INTENT(OUT) :: I_ZONE_OVERLAP
INTEGER :: NN,IOR,IC,IC_OLD,IIN,JJN,KKN,III,JJJ,KKK,Q_N
INTEGER, ALLOCATABLE, DIMENSION(:) :: Q_I,Q_J,Q_K
CHARACTER(MESSAGE_LENGTH) :: MESSAGE
TYPE (MESH_TYPE), POINTER :: M
TYPE (OBSTRUCTION_TYPE), POINTER :: OB

M => MESHES(NM)

I_ZONE_OVERLAP = -1
IC = M%CELL_INDEX(II,JJ,KK)

IF (M%CELL(IC)%SOLID) THEN
   WRITE(MESSAGE,'(A,I3,A)')  'ERROR: XYZ point for ZONE ',I_ZONE,' is inside a solid obstruction. Choose another XYZ.'
   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
   RETURN
ENDIF

IF (.NOT.M%CELL(IC)%SOLID .AND. M%PRESSURE_ZONE(II,JJ,KK)>=0 .AND.  M%PRESSURE_ZONE(II,JJ,KK)/=I_ZONE) THEN
   I_ZONE_OVERLAP = M%PRESSURE_ZONE(II,JJ,KK)
   RETURN
ELSEIF (M%PRESSURE_ZONE(II,JJ,KK)>=0) THEN
   RETURN
ELSE
   M%PRESSURE_ZONE(II,JJ,KK) = I_ZONE
ENDIF

! Add the first entry to "queue" of cells that need a pressure zone number

ALLOCATE(Q_I(M%IBAR*M%JBAR*M%KBAR))
ALLOCATE(Q_J(M%IBAR*M%JBAR*M%KBAR))
ALLOCATE(Q_K(M%IBAR*M%JBAR*M%KBAR))

Q_I(1) = II
Q_J(1) = JJ
Q_K(1) = KK
Q_N    = 1

! Look to all cells adjacent to the starting cell and determine if they are in the ZONE as well.
! Repeat process until all cells are found.

SORT_QUEUE: DO

   IF (Q_N<1) EXIT SORT_QUEUE

   III = Q_I(Q_N)
   JJJ = Q_J(Q_N)
   KKK = Q_K(Q_N)
   Q_N = Q_N - 1

   SEARCH_LOOP: DO IOR=-3,3

      IF (IOR==0) CYCLE SEARCH_LOOP

      IC  = M%CELL_INDEX(III,JJJ,KKK)

      IF (M%CELL(IC)%SOLID) THEN
         IF (.NOT.M%OBSTRUCTION(M%CELL(IC)%OBST_INDEX)%REMOVABLE) CYCLE SEARCH_LOOP  ! Do not search within a non-removable solid
      ENDIF
      IF (CC_IBM) THEN
         IF(M%CCVAR(III,JJJ,KKK,1)==1) CYCLE SEARCH_LOOP  ! Cycle if cell is of type CC_SOLID
      ENDIF

      SELECT CASE(IOR)
         CASE(-1)
            IF (III==1) CYCLE SEARCH_LOOP
            IIN = III-1
            JJN = JJJ
            KKN = KKK
         CASE( 1)
            IF (III==M%IBAR) CYCLE SEARCH_LOOP
            IIN = III+1
            JJN = JJJ
            KKN = KKK
         CASE(-2)
            IF (JJJ==1) CYCLE SEARCH_LOOP
            IIN = III
            JJN = JJJ-1
            KKN = KKK
         CASE( 2)
            IF (JJJ==M%JBAR) CYCLE SEARCH_LOOP
            IIN = III
            JJN = JJJ+1
            KKN = KKK
         CASE(-3)
            IF (KKK==1) CYCLE SEARCH_LOOP
            IIN = III
            JJN = JJJ
            KKN = KKK-1
         CASE( 3)
            IF (KKK==M%KBAR) CYCLE SEARCH_LOOP
            IIN = III
            JJN = JJJ
            KKN = KKK+1
      END SELECT

      ! Look for thin obstructions bordering the current cell

      DO NN=1,M%N_OBST
         OB=>M%OBSTRUCTION(NN)
         SELECT CASE(IOR)
            CASE(-1)
               IF (IIN==  OB%I1.AND.IIN==  OB%I2.AND.JJN>OB%J1.AND.JJN<=OB%J2.AND.KKN>OB%K1.AND.KKN<=OB%K2) CYCLE SEARCH_LOOP
            CASE( 1)
               IF (IIN-1==OB%I1.AND.IIN-1==OB%I2.AND.JJN>OB%J1.AND.JJN<=OB%J2.AND.KKN>OB%K1.AND.KKN<=OB%K2) CYCLE SEARCH_LOOP
            CASE(-2)
               IF (JJN==  OB%J1.AND.JJN==  OB%J2.AND.IIN>OB%I1.AND.IIN<=OB%I2.AND.KKN>OB%K1.AND.KKN<=OB%K2) CYCLE SEARCH_LOOP
            CASE( 2)
               IF (JJN-1==OB%J1.AND.JJN-1==OB%J2.AND.IIN>OB%I1.AND.IIN<=OB%I2.AND.KKN>OB%K1.AND.KKN<=OB%K2) CYCLE SEARCH_LOOP
            CASE(-3)
               IF (KKN==  OB%K1.AND.KKN==  OB%K2.AND.IIN>OB%I1.AND.IIN<=OB%I2.AND.JJN>OB%J1.AND.JJN<=OB%J2) CYCLE SEARCH_LOOP
            CASE( 3)
               IF (KKN-1==OB%K1.AND.KKN-1==OB%K2.AND.IIN>OB%I1.AND.IIN<=OB%I2.AND.JJN>OB%J1.AND.JJN<=OB%J2) CYCLE SEARCH_LOOP
         END SELECT
      ENDDO

      ! If the last cell was solid and the current cell is not solid, stop the current directional march.

      IC_OLD = IC
      IC = M%CELL_INDEX(IIN,JJN,KKN)

      IF (M%CELL(IC_OLD)%SOLID .AND. .NOT.M%CELL(IC)%SOLID) CYCLE SEARCH_LOOP

      ! If the current cell is not solid, but it is assigned another ZONE index, mark it as an overlap error and return

      IF (CC_IBM) THEN

         ! Here CC_CGSC=1, CC_SOLID=1:
         IF (M%CCVAR(III,JJJ,KKK,1)==1 .AND. M%CCVAR(IIN,JJN,KKN,1)/=1) CYCLE SEARCH_LOOP
         IF (M%CCVAR(III,JJJ,KKK,1)==0 .AND. M%CCVAR(IIN,JJN,KKN,1)==0) THEN
            IF(IOR>0 .AND. M%FCVAR(III,JJJ,KKK,1,ABS(IOR))==1) CYCLE SEARCH_LOOP
            IF(IOR<0 .AND. M%FCVAR( IIN, JJN, KKN,1,ABS(IOR))==1) CYCLE SEARCH_LOOP
         ENDIF

         ! Cell not SOLID for OBSTS, or GEOM cell not CC_SOLID:
         IF (.NOT.M%CELL(IC)%SOLID .AND. M%CCVAR(IIN,JJN,KKN,1)/=1 .AND. &
            M%PRESSURE_ZONE(IIN,JJN,KKN)>=0 .AND.  M%PRESSURE_ZONE(IIN,JJN,KKN)/=I_ZONE) THEN
            I_ZONE_OVERLAP = M%PRESSURE_ZONE(IIN,JJN,KKN)
            EXIT SORT_QUEUE
         ENDIF

      ELSEIF (.NOT.M%CELL(IC)%SOLID .AND. M%PRESSURE_ZONE(IIN,JJN,KKN)>=0 .AND.  M%PRESSURE_ZONE(IIN,JJN,KKN)/=I_ZONE) THEN

         I_ZONE_OVERLAP = M%PRESSURE_ZONE(IIN,JJN,KKN)
         EXIT SORT_QUEUE

      ENDIF

      ! If the current cell is unassigned, assign the cell the ZONE index, I_ZONE, and then add this cell to the
      ! queue so that further searches might originate from it.

      IF (M%PRESSURE_ZONE(IIN,JJN,KKN)<0) THEN
         M%PRESSURE_ZONE(IIN,JJN,KKN) = I_ZONE
         Q_N      = Q_N+1
         Q_I(Q_N) = IIN
         Q_J(Q_N) = JJN
         Q_K(Q_N) = KKN
      ENDIF

   ENDDO SEARCH_LOOP

ENDDO SORT_QUEUE

DEALLOCATE(Q_I)
DEALLOCATE(Q_J)
DEALLOCATE(Q_K)

END SUBROUTINE ASSIGN_PRESSURE_ZONE


!> \brief Determine the number of 1-D cells in a layer of solid material
!> \param DIFFUSIVITY \f$ k/(\rho c) \; (\hbox{m}^2/\hbox{s}) \f$, used to determine cell size
!> \param LAYER_THICKNESS Thickness of the material layer (m)
!> \param STRETCH_FACTOR A real array indicating the amount to stretch the second and subsequent cells in the interior
!> \param CELL_SIZE_FACTOR A real number indicating how much to shrink the first cell at the surface
!> \param CELL_SIZE Width of a uniform cell. This makes STRETCH_FACTOR, CELL_SIZE_FACTOR moot.
!> \param N_LAYER_CELLS_MAX Maximum number of cells to assign to the layer
!> \param N_CELLS Number of cells in the layer
!> \param DX_MIN Minimum cell size
!> \param DDSUM Divisor of LAYER_THICKNESS used to determine N_LAYER_CELLS

SUBROUTINE GET_N_LAYER_CELLS(DIFFUSIVITY,LAYER_THICKNESS,STRETCH_FACTOR,CELL_SIZE_FACTOR,CELL_SIZE,&
                             N_LAYER_CELLS_MAX,N_CELLS,DX_MIN,DDSUM)

INTEGER, INTENT(OUT)  :: N_CELLS
INTEGER, INTENT(IN) :: N_LAYER_CELLS_MAX
REAL(EB), INTENT(OUT) :: DX_MIN, DDSUM
REAL(EB), INTENT(IN) :: DIFFUSIVITY,LAYER_THICKNESS,STRETCH_FACTOR,CELL_SIZE_FACTOR,CELL_SIZE
INTEGER  :: N, I

DX_MIN = 0._EB
IF (LAYER_THICKNESS<=TWO_EPSILON_EB) THEN
   N_CELLS = 0
   DX_MIN = 0._EB
   RETURN
ENDIF

IF (CELL_SIZE>0._EB) THEN
   N_CELLS = MAX(1,NINT(LAYER_THICKNESS/CELL_SIZE))
   DDSUM   = REAL(N_CELLS)
   DX_MIN  = LAYER_THICKNESS/DDSUM
   RETURN
ENDIF

IF (N_LAYER_CELLS_MAX==1) THEN
   N_CELLS = 1
   DX_MIN = LAYER_THICKNESS
   RETURN
ENDIF

SHRINK_LOOP: DO N=1,N_LAYER_CELLS_MAX-1
   DDSUM = 0._EB
   SUM_LOOP: DO I=1,N
      DDSUM = DDSUM + STRETCH_FACTOR**(MIN(I-1,N-I))
   ENDDO SUM_LOOP
   IF ((LAYER_THICKNESS/DDSUM < CELL_SIZE_FACTOR*SQRT(DIFFUSIVITY)) .OR. (N==N_LAYER_CELLS_MAX-1)) THEN
      N_CELLS = N
      DX_MIN = LAYER_THICKNESS/DDSUM
      EXIT SHRINK_LOOP
   ENDIF
ENDDO SHRINK_LOOP

END SUBROUTINE GET_N_LAYER_CELLS


!> \brief Determine the coordinates of the 1-D interior cells
!> \param N_CELLS Number of cells in the entire collection of layers
!> \param N_CELLS_OLD Former value of N_CELLS
!> \param N_LAYERS Number of layers
!> \param N_LAYER_CELLS Array holding the number of cells in each layer
!> \param N_LAYER_CELLS_OLD Former value of N_LAYER_CELLS
!> \param SMALLEST_CELL_SIZE Array holding the sizes of the bounding cells in each layer (m)
!> \param STRETCH_FACTOR Node stretching factor of each layer
!> \param REMESH_LAYER Logical array indicating if current layer should be remeshed
!> \param X_S Array containing boundaries of the interior cells
!> \param X_S_OLD Array containing old boundaries of the interior cells used when skipping remesh for a layer
!> \param LAYER_THICKNESS Array containing layer thicknesses

SUBROUTINE GET_WALL_NODE_COORDINATES(N_CELLS,N_CELLS_OLD,N_LAYERS,N_LAYER_CELLS,N_LAYER_CELLS_OLD,SMALLEST_CELL_SIZE,&
                                     STRETCH_FACTOR,REMESH_LAYER,X_S,X_S_OLD,LAYER_THICKNESS)

INTEGER, INTENT(IN) :: N_CELLS,N_CELLS_OLD,N_LAYERS, N_LAYER_CELLS(N_LAYERS),N_LAYER_CELLS_OLD(N_LAYERS)
REAL(EB), INTENT(IN) :: SMALLEST_CELL_SIZE(N_LAYERS),STRETCH_FACTOR(N_LAYERS),X_S_OLD(0:N_CELLS_OLD),LAYER_THICKNESS(N_LAYERS)
REAL(EB), INTENT(OUT) :: X_S(0:N_CELLS)
LOGICAL, INTENT(IN) :: REMESH_LAYER(N_LAYERS)

INTEGER I, II, NL, IL, I_START
REAL(EB) DX_S,DX_SUM

II = 0
IL = 0
X_S(0) = 0._EB
I_START = 0

DO NL=1,N_LAYERS
   IF (NL > 1) I_START = I_START+N_LAYER_CELLS_OLD(NL-1)
   DX_SUM = 0._EB
   DO I=1,N_LAYER_CELLS(NL)-1
      II = II + 1
      IF (REMESH_LAYER(NL)) THEN
         DX_S = SMALLEST_CELL_SIZE(NL)*STRETCH_FACTOR(NL)**(MIN(I-1,N_LAYER_CELLS(NL)-I))
         DX_SUM = DX_SUM + DX_S
      ELSE
         DX_S = X_S_OLD(I_START+I) - X_S_OLD(I_START+I-1)
         DX_SUM = DX_SUM + DX_S
      ENDIF
      X_S(II) = X_S(II-1) + DX_S
   ENDDO
   IF (N_LAYER_CELLS(NL) > 0) THEN
      II = II + 1
      DX_S = LAYER_THICKNESS(NL) - DX_SUM
      X_S(II) = X_S(II-1) + DX_S
   ENDIF
   IL = IL + N_LAYER_CELLS_OLD(NL)
ENDDO

END SUBROUTINE GET_WALL_NODE_COORDINATES


!> \brief Determine the internal coordinates of the 1-D grid inside a solid.
!> \param N_CELLS Number of cells in the entire 1-D array
!> \param N_LAYERS Number of layers
!> \param N_LAYER_CELLS Array holding the number of cells in each layer
!> \param LAYER_THICKNESS Array of layer thicknesses (m)
!> \param GEOMETRY Indicates whether the surface is Cartesian, cylindrical, or spherical
!> \param X_S Array holding node coordinates (m)
!> \param LAYER_DIVIDE Number of layers that off gas to the front surface
!> \param DX Array of cell sizes (m)
!> \param RDX Reciprocal of DX (1/m)
!> \param RDXN Reciprocal of the distance from cell center to cell center (1/m)
!> \param DX_WGT Ratio of cell sizes: DX_WGT(I)=DX(I)*RDXN(I)*2
!> \param DXF Size of cell nearest the front surface (m)
!> \param DXB Size of cell nearest the back surface (m)
!> \param LAYER_INDEX Array of indices indicating the layer to which each interior cell belongs
!> \param MF_FRAC Array containing the fraction of each cells mass that is assigned to the front surface
!> \param INNER_RADIUS Inner radius of hollow cylinder or sphere (m)
!> \param X_DIVIDE Depth at which pyrolyzates move to back side (m)

SUBROUTINE GET_WALL_NODE_WEIGHTS(N_CELLS,N_LAYERS,N_LAYER_CELLS, &
         LAYER_THICKNESS,GEOMETRY,X_S,LAYER_DIVIDE,DX,RDX,RDXN,DX_WGT,DXF,DXB,LAYER_INDEX,MF_FRAC,INNER_RADIUS,X_DIVIDE)

! Get the wall internal coordinates

INTEGER, INTENT(IN)     :: N_CELLS, N_LAYERS, N_LAYER_CELLS(N_LAYERS),GEOMETRY
REAL(EB), INTENT(IN)    :: X_S(0:N_CELLS),LAYER_THICKNESS(1:N_LAYERS),LAYER_DIVIDE,INNER_RADIUS
INTEGER, INTENT(OUT)    :: LAYER_INDEX(0:N_CELLS+1)
REAL(EB), INTENT(OUT)   :: DX(1:N_CELLS),RDX(0:N_CELLS+1),RDXN(0:N_CELLS),DX_WGT(0:N_CELLS),DXF,DXB, &
                           MF_FRAC(1:N_CELLS),X_DIVIDE

INTEGER :: I, II, NL, I_GRAD
REAL(EB) :: R, THICKNESS

   THICKNESS = SUM(LAYER_THICKNESS)

   SELECT CASE(GEOMETRY)
      CASE(SURF_CARTESIAN)
         I_GRAD = 0
      CASE(SURF_CYLINDRICAL,SURF_INNER_CYLINDRICAL)
         I_GRAD = 1
      CASE(SURF_SPHERICAL)
         I_GRAD = 2
   END SELECT

   II = 0
   DO NL=1,N_LAYERS
      LAYER_INDEX(II+1:II+N_LAYER_CELLS(NL)) = NL
      II = II + N_LAYER_CELLS(NL)
   ENDDO
   LAYER_INDEX(0) = LAYER_INDEX(1)
   LAYER_INDEX(N_CELLS+1) = LAYER_INDEX(N_CELLS)
   DXF = X_S(1)       - X_S(0)
   DXB = X_S(N_CELLS) - X_S(N_CELLS-1)

   ! Compute dx_weight for each node (dx_weight is the ratio of cell size to the combined size of the current and next cell)

   DO I=1,N_CELLS-1
      DX_WGT(I) = (X_S(I)-X_S(I-1))/(X_S(I+1)-X_S(I-1))
   ENDDO
   DX_WGT(0)       = 0.5_EB
   DX_WGT(N_CELLS) = 0.5_EB

   ! Compute dx and 1/dx for each node (dx is the distance from cell boundary to cell boundary)

   DO I=1,N_CELLS
      DX(I)  = X_S(I)-X_S(I-1)
      RDX(I) = 1._EB/DX(I)
   ENDDO

   ! Adjust 1/dx_n to 1/(r dr) for cylindrical case and 1/(r^2 dr) for spherical

   SELECT CASE(GEOMETRY)
      CASE(SURF_CYLINDRICAL)
         DO I=1,N_CELLS
            R = 0.5_EB*RDX(I)*((INNER_RADIUS+THICKNESS-X_S(I-1))**2-(INNER_RADIUS+THICKNESS-X_S(I))**2)
            RDX(I) = RDX(I)/R**I_GRAD
         ENDDO
      CASE(SURF_INNER_CYLINDRICAL)
         DO I=1,N_CELLS
            R = 0.5_EB*RDX(I)*((INNER_RADIUS+X_S(I))**2-(INNER_RADIUS+X_S(I-1))**2)
            RDX(I) = RDX(I)/R**I_GRAD
         ENDDO
      CASE(SURF_SPHERICAL)
         DO I=1,N_CELLS
            R = SQRT(ONTH*RDX(I)*((INNER_RADIUS+THICKNESS-X_S(I-1))**3-(INNER_RADIUS+THICKNESS-X_S(I))**3))
            RDX(I) = RDX(I)/R**I_GRAD
         ENDDO
   END SELECT
   RDX(0)         = RDX(1)
   RDX(N_CELLS+1) = RDX(N_CELLS)

   ! Compute 1/dx_n for each node (dx_n is the distance from cell center to cell center)

   DO I=1,N_CELLS-1
      RDXN(I) = 2._EB/(X_S(I+1)-X_S(I-1))
   ENDDO
   RDXN(0)       = 1._EB/(X_S(1)-X_S(0))
   RDXN(N_CELLS) = 1._EB/(X_S(N_CELLS)-X_S(N_CELLS-1))

   ! Adjust 1/dx_n to r/dr for cylindrical case and r^2/dr for spaherical

   IF (GEOMETRY == SURF_INNER_CYLINDRICAL) THEN
      DO I=0,N_CELLS
         R = INNER_RADIUS+X_S(I)
         RDXN(I) = RDXN(I)*R**I_GRAD
      ENDDO
   ELSEIF (GEOMETRY /= SURF_CARTESIAN) THEN
      DO I=0,N_CELLS
         R = INNER_RADIUS+THICKNESS-X_S(I)
         RDXN(I) = RDXN(I)*R**I_GRAD
      ENDDO
   ENDIF

   ! Compute mass flux fraction array (array numbers indicate the fraction of mass flux that is added to the front

   IF (LAYER_DIVIDE >= REAL(N_LAYERS,EB)) THEN

      MF_FRAC = 1._EB
      X_DIVIDE = THICKNESS

   ELSE

      MF_FRAC = 0._EB

      IF (LAYER_DIVIDE>=0._EB) THEN
         X_DIVIDE = 0._EB
         DO NL=1,N_LAYERS
            IF (LAYER_DIVIDE>=REAL(NL,EB)) THEN
               X_DIVIDE  = X_DIVIDE + LAYER_THICKNESS(NL)
            ELSE
               X_DIVIDE  = X_DIVIDE + MOD(LAYER_DIVIDE,1.0_EB)*LAYER_THICKNESS(NL)
               EXIT
            ENDIF
         ENDDO
      ELSE
         X_DIVIDE = 0.5_EB*THICKNESS
      ENDIF

      II = 0
      DIVILOOP: DO NL=1,N_LAYERS
         DO I=1,N_LAYER_CELLS(NL)
            II = II + 1
            IF (X_S(II) < X_DIVIDE) THEN
               MF_FRAC(II) = 1._EB
            ELSEIF (X_S(II-1) < X_DIVIDE) THEN
               MF_FRAC(II) = (X_DIVIDE-X_S(II-1))/DX(II)
               EXIT DIVILOOP
            ENDIF
         ENDDO
      ENDDO DIVILOOP

   ENDIF

END SUBROUTINE GET_WALL_NODE_WEIGHTS


!> \brief Determine weighting factors for 1-D solid cells
!> \param GEOMETRY Indicator of surface geometry: Cartesian, cylindrical, or spherical
!> \param NWP Number of interior cells
!> \param NWP_NEW Number of interior cells after shrinkage or swelling
!> \param INNER_RADIUS Inner radius of hollow cylinder or sphere (m)
!> \param X_S Array of interior cell edge positions (m)
!> \param X_S_NEW Array of interior cell edge positions after shrinkage or swelling (m)
!> \param INT_WGT Array of weighting factors for new arrangement of interior cells

SUBROUTINE GET_INTERPOLATION_WEIGHTS(GEOMETRY,NWP,NWP_NEW,INNER_RADIUS,X_S,X_S_NEW,INT_WGT)

INTEGER, INTENT(IN)  :: GEOMETRY,NWP,NWP_NEW
REAL(EB), INTENT(IN) :: X_S(0:NWP), X_S_NEW(0:NWP_NEW), INNER_RADIUS
REAL(EB), INTENT(OUT) :: INT_WGT(NWP_NEW,NWP)

REAL(EB) XUP,XLOW,XUP_NEW,XLOW_NEW,VOL_NEW,VOL,THICKNESS
INTEGER I_NEW, I_OLD, I_GRAD


SELECT CASE(GEOMETRY)
   CASE(SURF_CARTESIAN)
      I_GRAD = 1
   CASE(SURF_CYLINDRICAL,SURF_INNER_CYLINDRICAL)
      I_GRAD = 2
   CASE(SURF_SPHERICAL)
      I_GRAD = 3
END SELECT

I_OLD = 1
INT_WGT = 0._EB
THICKNESS = X_S(NWP)
POINT_LOOP: DO I_NEW=1,NWP_NEW
   XLOW_NEW = X_S_NEW(I_NEW-1)
   XUP_NEW = X_S_NEW(I_NEW)
   OLD_POINT_LOOP: DO
      XLOW = X_S(I_OLD-1)
      XUP  = X_S(I_OLD)
      IF (GEOMETRY==SURF_INNER_CYLINDRICAL) THEN
         VOL = (INNER_RADIUS+XUP)**I_GRAD-(INNER_RADIUS+XLOW)**I_GRAD
         VOL_NEW = (INNER_RADIUS+MAX(XUP_NEW,XUP))**I_GRAD-(INNER_RADIUS+MIN(XLOW_NEW,XLOW))**I_GRAD
      ELSE
         VOL = (THICKNESS+INNER_RADIUS-XLOW)**I_GRAD-(THICKNESS+INNER_RADIUS-XUP)**I_GRAD
         VOL_NEW = (THICKNESS+INNER_RADIUS-MAX(XLOW_NEW,XLOW))**I_GRAD-(THICKNESS+INNER_RADIUS-MIN(XUP_NEW,XUP))**I_GRAD
      ENDIF
      IF (VOL > 0._EB) INT_WGT(I_NEW,I_OLD) = MAX(0._EB,MIN(1._EB,VOL_NEW/VOL))
      IF (XUP >= XUP_NEW .OR. I_OLD==NWP) EXIT OLD_POINT_LOOP
      I_OLD = I_OLD+1
   ENDDO OLD_POINT_LOOP

ENDDO POINT_LOOP

END SUBROUTINE GET_INTERPOLATION_WEIGHTS


!> \brief Interpolates old array of wall properties to reflect new wall noding
!>
!> \param N_CELLS Maximum of the previous or current number of wall nodes
!> \param NWP Previous number of wall nodes
!> \param NWP_NEW Current number of wall nodes
!> \param INT_WGT array of weighting factors mapping old wall nodes to new wall nodes
!> \param ARR On input this is the old array and on output this is the new array

SUBROUTINE INTERPOLATE_WALL_ARRAY(N_CELLS,NWP,NWP_NEW,INT_WGT,ARR)

INTEGER, INTENT(IN)  :: N_CELLS,NWP,NWP_NEW
REAL(EB), INTENT(IN) :: INT_WGT(NWP_NEW,NWP)
REAL(EB), INTENT(INOUT ):: ARR(N_CELLS)
REAL(EB) :: TMP(N_CELLS)

INTEGER I,J

TMP = ARR
ARR = 0._EB
DO I = 1,NWP_NEW
  DO J = 1,NWP
        ARR(I) = ARR(I) + INT_WGT(I,J) * TMP(J)
  ENDDO
ENDDO

END SUBROUTINE INTERPOLATE_WALL_ARRAY


!> \brief Choose a random point from within a rectangular volume
!> \param XX x-coordinate of the point (m)
!> \param YY y-coordinate of the point (m)
!> \param ZZ z-coordinate of the point (m)
!> \param X1 Lower x-coordinate of the volume (m)
!> \param X2 Upper x-coordinate of the volume (m)
!> \param Y1 Lower y-coordinate of the volume (m)
!> \param Y2 Upper y-coordinate of the volume (m)
!> \param Z1 Lower z-coordinate of the volume (m)
!> \param Z2 Upper z-coordinate of the volume (m)

SUBROUTINE RANDOM_RECTANGLE(XX,YY,ZZ,X1,X2,Y1,Y2,Z1,Z2)

REAL(EB), INTENT(IN) :: X1,X2,Y1,Y2,Z1,Z2
REAL(EB), INTENT(OUT) :: XX,YY,ZZ
REAL(EB) :: RN

CALL RANDOM_NUMBER(RN)
XX = X1 + RN*(X2-X1)
CALL RANDOM_NUMBER(RN)
YY = Y1 + RN*(Y2-Y1)
CALL RANDOM_NUMBER(RN)
ZZ = Z1 + RN*(Z2-Z1)

END SUBROUTINE RANDOM_RECTANGLE


!> \brief Choose a random point from within a vertically oriented cone
!> \param NM Mesh number
!> \param XX x-coordinate of the point (m)
!> \param YY y-coordinate of the point (m)
!> \param ZZ z-coordinate of the point (m)
!> \param X0 x-coordinate of the center of the code base (m)
!> \param Y0 y-coordinate of the center of the code base (m)
!> \param Z0 z-coordinate of the center of the code base (m)
!> \param RR0 Radius of the base (m)
!> \param RRI Inner radius of the base of the cone (if hollow) (m)
!> \param HH0 Height of the cone (m)
!> \param G_FACTOR Geometry indicator: 1 for a cone, 0 for a cylinder

SUBROUTINE RANDOM_CONE(NM,XX,YY,ZZ,X0,Y0,Z0,RR0,RRI,HH0,G_FACTOR)

INTEGER, INTENT(IN) :: NM,G_FACTOR
REAL(EB), INTENT(IN) :: X0,Y0,Z0,RR0,RRI,HH0
REAL(EB), INTENT(OUT) :: XX,YY,ZZ
REAL(EB) :: THETA,RR,RN,RR1,RR2
TYPE (MESH_TYPE), POINTER :: M

M => MESHES(NM)
SEARCH_LOOP: DO
   CALL RANDOM_NUMBER(RN)
   ZZ = Z0 + HH0*RN
   CALL RANDOM_NUMBER(RN)
   THETA = TWOPI*RN
   CALL RANDOM_NUMBER(RN)
   RR = SQRT(RN)*RR0
   XX = X0 + RR*COS(THETA)
   YY = Y0 + RR*SIN(THETA)
   IF (G_FACTOR==1) THEN
      RR1=RRI*(1._EB-(ZZ-Z0)/HH0)
      RR2=RR0*(1._EB-(ZZ-Z0)/HH0)
   ELSE
      RR1=RRI
      RR2=RR0
   ENDIF
   IF ((RR**2>RR2**2) .OR. (RR**2<RR1**2)) CYCLE SEARCH_LOOP
   IF (XX>=M%XS .AND. XX<=M%XF .AND. YY>=M%YS .AND. YY<=M%YF .AND. ZZ>=M%ZS .AND. ZZ<=M%ZF) EXIT SEARCH_LOOP
ENDDO SEARCH_LOOP

END SUBROUTINE RANDOM_CONE


!> \brief Choose a random point on a horizontally-oriented ring
!> \param NM Mesh number
!> \param XX x-coordinate of the point (m)
!> \param YY y-coordinate of the point (m)
!> \param X0 x-coordinate of the center of the ring (m)
!> \param Y0 y-coordinate of the center of the ring (m)
!> \param RR0 Radius of the ring (m)

SUBROUTINE RANDOM_RING(NM,XX,YY,X0,Y0,RR0)

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: X0,Y0,RR0
REAL(EB), INTENT(OUT) :: XX,YY
REAL(EB) :: THETA,RN
TYPE (MESH_TYPE), POINTER :: M

M => MESHES(NM)
SEARCH_LOOP: DO
   CALL RANDOM_NUMBER(RN)
   THETA = TWOPI*RN
   XX = X0 + RR0*COS(THETA)
   YY = Y0 + RR0*SIN(THETA)
   IF (XX>=M%XS .AND. XX<=M%XF .AND. YY>=M%YS .AND. YY<=M%YF) EXIT SEARCH_LOOP
ENDDO SEARCH_LOOP

END SUBROUTINE RANDOM_RING


!> \brief Assign an ordered point on a horizontally-oriented ring
!> \param XX x-coordinate of the point (m)
!> \param YY y-coordinate of the point (m)
!> \param X0 x-coordinate of the center of the ring (m)
!> \param Y0 y-coordinate of the center of the ring (m)
!> \param RR0 Radius of the ring (m)
!> \param NP Number of the current point
!> \param NP_TOTAL Total number of points in the ring

SUBROUTINE UNIFORM_RING(XX,YY,X0,Y0,RR0,NP,NP_TOTAL)

INTEGER, INTENT(IN) :: NP,NP_TOTAL
REAL(EB), INTENT(IN) :: X0,Y0,RR0
REAL(EB), INTENT(OUT) :: XX,YY
REAL(EB) :: THETA

THETA = TWOPI*REAL(NP,EB)/REAL(NP_TOTAL,EB)
XX = X0 + RR0*COS(THETA)
YY = Y0 + RR0*SIN(THETA)

END SUBROUTINE UNIFORM_RING


!> \brief Compute the volume of the intersection of a given rectangular volume and a given mesh
!> \param NM Mesh number
!> \param X1 Lower x-coordinate of the rectangular volume (m)
!> \param X2 Upper x-coordinate of the rectangular volume (m)
!> \param Y1 Lower y-coordinate of the rectangular volume (m)
!> \param Y2 Upper y-coordinate of the rectangular volume (m)
!> \param Z1 Lower z-coordinate of the rectangular volume (m)
!> \param Z2 Upper z-coordinate of the rectangular volume (m)
!> \param VOLUME_ADJUST Multiplicative factor to account for mismatch between block volume and FDS mesh volume

SUBROUTINE BLOCK_MESH_INTERSECTION_VOLUME(NM,X1,X2,Y1,Y2,Z1,Z2,VOLUME_ADJUST)

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: X1,X2,Y1,Y2,Z1,Z2
REAL(EB), INTENT(OUT) :: VOLUME_ADJUST
REAL(EB) :: X1_C,X2_C,Y1_C,Y2_C,Z1_C,Z2_C,ACTUAL_VOLUME,FDS_VOLUME
INTEGER :: I,J,K
TYPE(MESH_TYPE), POINTER :: M

VOLUME_ADJUST = 1._EB

IF (MY_RANK/=PROCESS(NM)) RETURN

M => MESHES(NM)

IF (X1>=M%XF .OR. X2<=M%XS .OR. Y1>=M%YF .OR. Y2<=M%YS .OR. Z1>=M%ZF .OR. Z2<=M%ZS) RETURN

X1_C = MAX(X1,M%XS)
X2_C = MIN(X2,M%XF)
Y1_C = MAX(Y1,M%YS)
Y2_C = MIN(Y2,M%YF)
Z1_C = MAX(Z1,M%ZS)
Z2_C = MIN(Z2,M%ZF)

ACTUAL_VOLUME = (X2_C-X1_C)*(Y2_C-Y1_C)*(Z2_C-Z1_C)
FDS_VOLUME = 0._EB

DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         IF (M%XC(I)<X1_C .OR. M%XC(I)>X2_C .OR. M%YC(J)<Y1_C .OR. M%YC(J)>Y2_C .OR. M%ZC(K)<Z1_C .OR. M%ZC(K)>Z2_C) CYCLE
         FDS_VOLUME = FDS_VOLUME + M%DX(I)*M%DY(J)*M%DZ(K)
      ENDDO
   ENDDO
ENDDO

IF (FDS_VOLUME>TWO_EPSILON_EB) VOLUME_ADJUST = ACTUAL_VOLUME/FDS_VOLUME

END SUBROUTINE BLOCK_MESH_INTERSECTION_VOLUME


!> \brief Estimate area of intersection of circle and grid cell
!> \param X0 x-coordinate of the center of the circle (m)
!> \param Y0 y-coordinate of the center of the circle (m)
!> \param RAD Radius of the circle (m)
!> \param X1 Lower x-coordinate of the grid cell (m)
!> \param X2 Upper x-coordinate of the grid cell (m)
!> \param Y1 Lower y-coordinate of the grid cell (m)
!> \param Y2 Upper y-coordinate of the grid cell (m)

REAL(EB) FUNCTION CIRCLE_CELL_INTERSECTION_AREA(X0,Y0,RAD,X1,X2,Y1,Y2)

REAL(EB), INTENT(IN) :: X0,Y0,RAD,X1,X2,Y1,Y2
INTEGER :: FC
REAL(EB) :: XC,YC,XP1,XP2,YP1,YP2
REAL(EB) :: RC,R2,THETA

R2 = RAD**2
FC = 0

! No overlap

IF ((X2 < X0-RAD) .OR. (X1 > X0+RAD) .OR. (Y2 < Y0-RAD) .OR. (Y1 > Y0+RAD)) THEN
   CIRCLE_CELL_INTERSECTION_AREA = 0._EB
   RETURN
ENDIF

! Count corners inside circle
RC = (X1-X0)**2+(Y1-Y0)**2
IF (RC <= R2) FC=IBSET(FC,0)
RC = (X1-X0)**2+(Y2-Y0)**2
IF (RC <= R2) FC=IBSET(FC,1)
RC = (X2-X0)**2+(Y2-Y0)**2
IF (RC <= R2) FC=IBSET(FC,2)
RC = (X2-X0)**2+(Y1-Y0)**2
IF (RC <= R2) FC=IBSET(FC,3)

SELECT CASE(FC)
   ! Grid cell surrounds the circle: area of circle
   CASE (0)
      CIRCLE_CELL_INTERSECTION_AREA = PI*R2
   ! One corner in the circle: chord + area of triangle
   CASE(1,2,4,8)
      SELECT CASE(FC)
         CASE(1)
            XC=X1
            YC=Y1
            YP1=SQRT(R2-(X1-X0)**2)+Y0
            XP1=SQRT(R2-(Y1-Y0)**2)+X0
         CASE(2)
            XC=X1
            YC=Y2
            YP1=Y0-SQRT(R2-(X1-X0)**2)
            XP1=SQRT(R2-(Y2-Y0)**2)+X0
         CASE(4)
            XC=X2
            YC=Y2
            YP1=Y0-SQRT(R2-(X2-X0)**2)
            XP1=X0-SQRT(R2-(Y2-Y0)**2)
         CASE(8)
            XC=X2
            YC=Y1
            YP1=SQRT(R2-(X2-X0)**2)+Y0
            XP1=X0-SQRT(R2-(Y1-Y0)**2)
      END SELECT
      THETA = 2._EB*ASIN(0.5_EB*SQRT((XC-XP1)**2+(YC-YP1)**2)/RAD)
      CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA)) + 0.5_EB*ABS(XC-XP1)*ABS(YC-YP1)
   ! Two corners in the circle: chord + trapezoid
   CASE(3)
      XP1=SQRT(R2-(Y1-Y0)**2)+X0
      XP2=SQRT(R2-(Y2-Y0)**2)+X0
      THETA = 2._EB*ASIN(0.5_EB*SQRT((XP1-XP2)**2+(Y2-Y1)**2)/RAD)
      CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA)) + 0.5_EB*(Y2-Y1)*(XP1+XP2-2._EB*X1)
   CASE(6)
      YP1=Y0-SQRT(R2-(X1-X0)**2)
      YP2=Y0-SQRT(R2-(X2-X0)**2)
      THETA = 2._EB*ASIN(0.5_EB*SQRT((X2-X1)**2+(YP2-YP1)**2)/RAD)
      CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA)) + 0.5_EB*(X2-X1)*(2._EB*Y2-YP1-YP2)
   CASE(9)
      YP1=SQRT(R2-(X1-X0)**2)+Y0
      YP2=SQRT(R2-(X2-X0)**2)+Y0
      THETA = 2._EB*ASIN(0.5_EB*SQRT((X2-X1)**2+(YP2-YP1)**2)/RAD)
      CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA)) + 0.5_EB*(X2-X1)*(YP1+YP2-2._EB*Y1)
   CASE(12)
      XP1=X0-SQRT(R2-(Y1-Y0)**2)
      XP2=X0-SQRT(R2-(Y2-Y0)**2)
      THETA = 2._EB*ASIN(0.5_EB*SQRT((XP1-XP2)**2+(Y2-Y1)**2)/RAD)
      CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA)) + 0.5_EB*(Y2-Y1)*(2._EB*X2-XP1-XP2)
   ! Three corners in the circle: chord + irregular pentagon (two rectangles and a triangle)
   CASE(7)
      YP1=Y0-SQRT(R2-(X2-X0)**2)
      XP1=SQRT(R2-(Y1-Y0)**2)+X0
      THETA = 2._EB*ASIN(0.5_EB*SQRT((XP1-X2)**2+(YP1-Y1)**2)/RAD)
      CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA)) + (Y2-Y1)*(XP1-X1) + (X2-XP1)*(Y2-YP1) + &
         0.5_EB*(X2-XP1)*(YP1-Y1)
   CASE(11)
      YP1=SQRT(R2-(X2-X0)**2)+Y0
      XP1=SQRT(R2-(Y2-Y0)**2)+X0
      THETA = 2._EB*ASIN(0.5_EB*SQRT((XP1-X2)**2+(YP1-Y2)**2)/RAD)
      CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA)) + (Y2-Y1)*(XP1-X1) + (X2-XP1)*(YP1-Y1) + &
         0.5_EB*(X2-XP1)*(Y2-YP1)
   CASE(13)
      YP1=SQRT(R2-(X1-X0)**2)+Y0
      XP1=X0-SQRT(R2-(Y2-Y0)**2)
      THETA = 2._EB*ASIN(0.5_EB*SQRT((XP1-X1)**2+(YP1-Y2)**2)/RAD)
      CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA)) + (Y2-Y1)*(X2-XP1) + (XP1-X1)*(YP1-Y1) + &
         0.5_EB*(XP1-X1)*(Y2-YP1)
   CASE(14)
      YP1=Y0-SQRT(R2-(X1-X0)**2)
      XP1=X0-SQRT(R2-(Y1-Y0)**2)
      THETA = 2._EB*ASIN(0.5_EB*SQRT((XP1-X1)**2+(YP1-Y1)**2)/RAD)
      CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA)) + (Y2-Y1)*(X2-XP1) + (XP1-X1)*(Y2-YP1) + &
         0.5_EB*(XP1-X1)*(YP1-Y1)
   ! Entire grid cell is in circle: area of grid cell
   CASE(15)
      CIRCLE_CELL_INTERSECTION_AREA = (X2-X1)*(Y2-Y1)
END SELECT

END FUNCTION CIRCLE_CELL_INTERSECTION_AREA


!> \brief Calculate volume of the interaction of MESH NM with a cone
!> \param NM Mesh number
!> \param X0 x-coordinate of the center of the base of the cone (m)
!> \param Y0 y-coordinate of the center of the base of the cone (m)
!> \param Z0 z-coordinate of the center of the base of the cone (m)
!> \param RR0 Radius of the base of the cone (m)
!> \param RRI Inner radius of the base of the cone (if hollow) (m)
!> \param HH0 Height of the cone (m)
!> \param G_FACTOR Geometry indicator: 1 for a cone, 0 for a cylinder

REAL(EB) FUNCTION CONE_MESH_INTERSECTION_VOLUME(NM,X0,Y0,Z0,RR0,RRI,HH0,G_FACTOR)

INTEGER, INTENT(IN) :: NM,G_FACTOR
REAL(EB), INTENT(IN) :: X0,Y0,Z0,RR0,RRI,HH0
INTEGER :: I,J,K,NX,NY,NZ
REAL(EB) :: XX,YY,ZZ,R2,X_MIN,X_MAX,Y_MIN,Y_MAX,Z_MIN,Z_MAX,DX,DY,DZ,SHORT_DIMENSION
TYPE (MESH_TYPE), POINTER :: M

CONE_MESH_INTERSECTION_VOLUME = 0._EB
M => MESHES(NM)

X_MIN = MAX(M%XS,X0-RR0)
X_MAX = MIN(M%XF,X0+RR0)
Y_MIN = MAX(M%YS,Y0-RR0)
Y_MAX = MIN(M%YF,Y0+RR0)
Z_MIN = MAX(M%ZS,Z0)
Z_MAX = MIN(M%ZF,Z0+HH0)
IF (X_MAX<=X_MIN .OR. Y_MAX<=Y_MIN .OR. Z_MAX<=Z_MIN) RETURN
SHORT_DIMENSION = MIN(X_MAX-X_MIN,Y_MAX-Y_MIN,Z_MAX-Z_MIN)
NX = CEILING(50*(X_MAX-X_MIN)/SHORT_DIMENSION)
NY = CEILING(50*(Y_MAX-Y_MIN)/SHORT_DIMENSION)
NZ = CEILING(50*(Z_MAX-Z_MIN)/SHORT_DIMENSION)
DX = (X_MAX-X_MIN)/REAL(NX,EB)
DY = (Y_MAX-Y_MIN)/REAL(NY,EB)
DZ = (Z_MAX-Z_MIN)/REAL(NZ,EB)

DO K=1,NZ
   ZZ = Z_MIN + (K-0.5_EB)*DZ
   IF (ZZ<Z0 .OR. ZZ>Z0+HH0) CYCLE
   DO J=1,NY
      YY = Y_MIN + (J-0.5_EB)*DY
      DO I=1,NX
         XX = X_MIN + (I-0.5_EB)*DX
         R2 = (XX-X0)**2+(YY-Y0)**2
         IF ((R2<(RR0*(1._EB-G_FACTOR*(ZZ-Z0)/HH0))**2) .AND. (R2>(RRI*(1._EB-G_FACTOR*(ZZ-Z0)/HH0))**2)) &
            CONE_MESH_INTERSECTION_VOLUME = CONE_MESH_INTERSECTION_VOLUME + DX*DY*DZ
      ENDDO
   ENDDO
ENDDO

END FUNCTION CONE_MESH_INTERSECTION_VOLUME


!> \brief Calculate arc length of the interaction of MESH NM with a ring
!> \param NM Mesh number
!> \param X0 x-coordinate of the center of the ring (m)
!> \param Y0 y-coordinate of the center of the ring (m)
!> \param RR0 Radius of the ring (m)

REAL(EB) FUNCTION RING_MESH_INTERSECTION_ARC(NM,X0,Y0,RR0)

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: X0,Y0,RR0
REAL(EB) :: XX,YY,THETA
INTEGER :: IR,N_ITER
TYPE (MESH_TYPE), POINTER :: M

RING_MESH_INTERSECTION_ARC = 0._EB
N_ITER = 1000
M => MESHES(NM)

DO IR=1,N_ITER
   THETA = TWOPI*REAL(IR,EB)/REAL(N_ITER,EB)
   XX = X0 + RR0*COS(THETA)
   YY = Y0 + RR0*SIN(THETA)
   IF (XX>M%XS .AND. XX<M%XF .AND. YY>M%YS .AND. YY<M%YF) &
      RING_MESH_INTERSECTION_ARC = RING_MESH_INTERSECTION_ARC + 1
ENDDO
RING_MESH_INTERSECTION_ARC = RING_MESH_INTERSECTION_ARC/REAL(N_ITER,EB)*TWOPI*RR0

END FUNCTION RING_MESH_INTERSECTION_ARC

!> \brief Rotate and translate a point using the parameters of MOVE input line
!> \param X x-coordinate of the point (m)
!> \param Y y-coordinate of the point (m)
!> \param Z z-coordinate of the point (m)
!> \param MOVE_INDEX Index of the set of parameters on a MOVE input line
!> \param MODE Mode of transformation: 1 for a point, 2 for a vector

SUBROUTINE TRANSFORM_COORDINATES(X,Y,Z,MOVE_INDEX,MODE)

INTEGER, INTENT(IN) :: MOVE_INDEX,MODE
REAL(EB), INTENT(INOUT) :: X,Y,Z
REAL(EB) :: M(3,3),UP(3,1),S(3,3),UUT(3,3),IDENTITY(3,3),X_VECTOR(3,1),X_VECTOR_0(3,1),X_VECTOR_1(4,1),SCL(3,3)
TYPE(MOVEMENT_TYPE), POINTER :: MV

MV => MOVEMENT(MOVE_INDEX)

! If T34 present, use Matrix T3x4 for transformation: Note it takes precedence over other transformation parameters.

IF (SUM(ABS(MV%T34(1:3,1:4))) > TWO_EPSILON_EB) THEN
   X_VECTOR_1 = RESHAPE( (/X,Y,Z,1._EB/), (/4, 1/) )
   X_VECTOR   = MATMUL(MV%T34, X_VECTOR_1)
   X = X_VECTOR(1,1)
   Y = X_VECTOR(2,1)
   Z = X_VECTOR(3,1)
   RETURN
ENDIF

X_VECTOR = RESHAPE( (/X,Y,Z/),(/3,1/) )
X_VECTOR_0 = RESHAPE( (/MV%X0,MV%Y0,MV%Z0/),(/3,1/) )
UP = RESHAPE(MV%AXIS,(/3,1/))
S =  RESHAPE( (/  0.0_EB, -UP(3,1),  UP(2,1),&
                 UP(3,1),   0.0_EB, -UP(1,1),&
                -UP(2,1),  UP(1,1),  0.0_EB  /),(/3,3/))
UUT = MATMUL(UP,TRANSPOSE(UP))
IDENTITY = RESHAPE ((/ 1.0_EB,0.0_EB,0.0_EB,&
                       0.0_EB,1.0_EB,0.0_EB,&
                       0.0_EB,0.0_EB,1.0_EB /),(/3,3/))
M = UUT + COS(MV%ROTATION_ANGLE*DEG2RAD)*(IDENTITY - UUT) + SIN(MV%ROTATION_ANGLE*DEG2RAD)*S

SCL= RESHAPE((/ MV%SCALE*MV%SCALEX,0.0_EB,0.0_EB, &
                0.0_EB,MV%SCALE*MV%SCALEY,0.0_EB, &
                0.0_EB,0.0_EB,MV%SCALE*MV%SCALEZ /),(/3,3/))

! Scaling takes precedence over rotation transformation (applied first on the local axes):

SELECT CASE(MODE)
   CASE(1)
      X_VECTOR = X_VECTOR_0 + MATMUL(MATMUL(M,SCL),X_VECTOR-X_VECTOR_0)
      X = X_VECTOR(1,1) + MV%DX
      Y = X_VECTOR(2,1) + MV%DY
      Z = X_VECTOR(3,1) + MV%DZ
   CASE(2)
      X_VECTOR = MATMUL(MATMUL(M,SCL),X_VECTOR)
      X = X_VECTOR(1,1)
      Y = X_VECTOR(2,1)
      Z = X_VECTOR(3,1)
END SELECT

END SUBROUTINE TRANSFORM_COORDINATES

END MODULE GEOMETRY_FUNCTIONS


!> \brief Coordinate transformation functions

MODULE TRAN

USE PRECISION_PARAMETERS
IMPLICIT NONE (TYPE,EXTERNAL)
TYPE TRAN_TYPE
   REAL(EB), POINTER, DIMENSION(:,:) :: C1=>NULL(),C2=>NULL(),C3=>NULL(),CCSTORE=>NULL(),PCSTORE=>NULL()
   INTEGER, POINTER, DIMENSION(:,:) :: IDERIVSTORE=>NULL()
   INTEGER NOC(3),ITRAN(3),NOCMAX
END TYPE TRAN_TYPE
TYPE (TRAN_TYPE), ALLOCATABLE, TARGET, DIMENSION(:) :: TRANS


CONTAINS


!> \brief Transformation function that maps uniform numerical coordinate to non-uniform physical coordinate
!> \param X Coordinate on uniform numerical grid (m)
!> \param IC Coordinate index, 1, 2, or 3
!> \param NM Mesh number

REAL(EB) FUNCTION G(X,IC,NM)

REAL(EB), INTENT(IN) :: X
INTEGER, INTENT(IN)  :: IC,NM
INTEGER :: I,II,N
TYPE (TRAN_TYPE), POINTER :: T=>NULL()

T => TRANS(NM)

N = T%NOC(IC)
IF (N==0) THEN
   G = X
   RETURN
ENDIF

SELECT CASE(T%ITRAN(IC))
   CASE(1)
      G = 0._EB
      DO I=1,N+1
         G = G + T%C1(I,IC)*X**I
      ENDDO
   CASE(2)
      ILOOP: DO I=1,N+1
         II = I
         IF (X<=T%C1(I,IC)) EXIT ILOOP
      ENDDO ILOOP
      G = T%C2(II-1,IC) + T%C3(II,IC)*(X-T%C1(II-1,IC))
END SELECT

END FUNCTION G


!> \brief First derivative of the transformation function that maps uniform numerical coordinate to non-uniform physical coordinate
!> \param X Coordinate on uniform numerical grid (m)
!> \param IC Coordinate index, 1, 2, or 3
!> \param NM Mesh number

REAL(EB) FUNCTION GP(X,IC,NM)

REAL(EB), INTENT(IN) :: X
INTEGER, INTENT(IN)  :: IC,NM
INTEGER :: I,II,N
TYPE (TRAN_TYPE), POINTER :: T=>NULL()

T => TRANS(NM)
N =  T%NOC(IC)
IF (N==0) THEN
   GP = 1._EB
   RETURN
ENDIF

SELECT CASE(T%ITRAN(IC))
   CASE(1)
      GP = 0._EB
      DO I=1,N+1
         GP = GP + I*T%C1(I,IC)*X**(I-1)
      ENDDO
   CASE(2)
      ILOOP: DO I=1,N+1
         II = I
         IF (X<=T%C1(I,IC)) EXIT ILOOP
      ENDDO ILOOP
      GP = T%C3(II,IC)
   CASE DEFAULT
      GP = 0._EB
END SELECT

END FUNCTION GP


!> \brief Inverse of the transformation function that maps uniform numerical coordinate to non-uniform physical coordinate
!> \param Z Coordinate on non-uniform physical grid (m)
!> \param IC Coordinate index, 1, 2, or 3
!> \param NM Mesh number

REAL(EB) FUNCTION GINV(Z,IC,NM)

REAL(EB) :: GF
INTEGER :: N,IT,II,I
REAL(EB), INTENT(IN) :: Z
INTEGER, INTENT(IN)  :: IC,NM
TYPE (TRAN_TYPE), POINTER :: T=>NULL()

T => TRANS(NM)
GINV = Z
N = T%NOC(IC)
IF (N==0) RETURN

SELECT CASE(T%ITRAN(IC))
   CASE(1)
      LOOP1: DO IT=1,10
         GF = G(GINV,IC,NM)-Z
         IF (ABS(GF)<0.00001_EB) EXIT LOOP1
         GINV = GINV - GF/GP(GINV,IC,NM)
      ENDDO LOOP1
   CASE(2)
      ILOOP: DO I=1,N+1
         II = I
         IF (Z<=T%C2(I,IC)) EXIT ILOOP
      ENDDO ILOOP
      GINV = T%C1(II-1,IC) + (Z-T%C2(II-1,IC))/T%C3(II,IC)
END SELECT

END FUNCTION GINV


!> \brief Determine the continuous and discrete form of the cell indices of a given point
!> \param X x-coordinate of the point (m)
!> \param Y y-coordinate of the point (m)
!> \param Z z-coordinate of the point (m)
!> \param NM Mesh number
!> \param XI Real number indicating the cell index of the point in the x direction
!> \param YJ Real number indicating the cell index of the point in the y direction
!> \param ZK Real number indicating the cell index of the point in the z direction
!> \param I Cell index of the point in the x direction
!> \param J Cell index of the point in the y direction
!> \param K Cell index of the point in the z direction
!> \details For example, if the mesh is the unit cube with 10 cm cells, for the point (X,Y,Z)=(0.23,0.46,0.66),
!> GET_IJK would return (XI,JK,ZK)=(2.3,4.6,6.6) and (I,J,K)=(3,5,7).

SUBROUTINE GET_IJK(X,Y,Z,NM,XI,YJ,ZK,I,J,K)

USE MESH_VARIABLES
REAL(EB),INTENT(IN) :: X,Y,Z
INTEGER,INTENT(IN) :: NM
REAL(EB), INTENT(OUT) :: XI,YJ,ZK
INTEGER, INTENT(OUT) :: I,J,K

XI = MESHES(NM)%CELLSI(  MIN( MESHES(NM)%CELLSI_HI, MAX(MESHES(NM)%CELLSI_LO,FLOOR((X-MESHES(NM)%XS)*MESHES(NM)%RDXINT)) )  )
YJ = MESHES(NM)%CELLSJ(  MIN( MESHES(NM)%CELLSJ_HI, MAX(MESHES(NM)%CELLSJ_LO,FLOOR((Y-MESHES(NM)%YS)*MESHES(NM)%RDYINT)) )  )
ZK = MESHES(NM)%CELLSK(  MIN( MESHES(NM)%CELLSK_HI, MAX(MESHES(NM)%CELLSK_LO,FLOOR((Z-MESHES(NM)%ZS)*MESHES(NM)%RDZINT)) )  )
I = FLOOR(XI+1._EB)
J = FLOOR(YJ+1._EB)
K = FLOOR(ZK+1._EB)

END SUBROUTINE GET_IJK


END MODULE TRAN



!> \brief Module for various OpenMP functions

MODULE OPENMP_FDS

USE GLOBAL_CONSTANTS, ONLY : OPENMP_AVAILABLE_THREADS, USE_OPENMP
!$ USE OMP_LIB
IMPLICIT NONE (TYPE,EXTERNAL)
PUBLIC OPENMP_INIT, OPENMP_PRINT_STATUS

CONTAINS

!> \brief Set the control flag USE_OPENMP if OpenMP is used.

SUBROUTINE OPENMP_INIT

!$OMP PARALLEL
!$OMP MASTER
!$ IF (OMP_GET_NUM_THREADS() /= 0) THEN
!$    USE_OPENMP = .TRUE.
!$    OPENMP_AVAILABLE_THREADS = OMP_GET_NUM_THREADS()
!$ ENDIF
!$OMP END MASTER
!$OMP BARRIER
!$OMP END PARALLEL

END SUBROUTINE OPENMP_INIT


!> \brief Write OpenMP status to standard error.

SUBROUTINE OPENMP_PRINT_STATUS
!   USE GLOBAL_CONSTANTS, ONLY : LU_ERR, MY_RANK, N_MPI_PROCESSES, VERBOSE
!   INTEGER :: THREAD_ID

!   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(THREAD_ID)

!   THREAD_ID = 0
!   !$ THREAD_ID = OMP_GET_THREAD_NUM()

!   !$OMP CRITICAL
!   IF (USE_OPENMP .AND. OPENMP_AVAILABLE_THREADS>1 .AND. VERBOSE) WRITE(LU_ERR,91) " OpenMP thread ",THREAD_ID," of ",&
!      OPENMP_AVAILABLE_THREADS-1," assigned to MPI process ",MY_RANK," of ",N_MPI_PROCESSES-1
!   IF (.NOT.USE_OPENMP .AND. VERBOSE) WRITE(LU_ERR,92) " MPI process ",MY_RANK," of ",N_MPI_PROCESSES-1
!   !$OMP END CRITICAL

!   !$OMP END PARALLEL

!   91 FORMAT(A,I3,A,I3,A,I6,A,I6)
!   92 FORMAT(A,I6,A,I6)

END SUBROUTINE OPENMP_PRINT_STATUS

END MODULE OPENMP_FDS



MODULE MISC_FUNCTIONS

USE PRECISION_PARAMETERS
IMPLICIT NONE (TYPE,EXTERNAL)

CONTAINS


!> \brief Determine the index of a RAMP with a given name (ID)
!> \param ID The name of the ramp function
!> \param TYPE The kind of ramp
!> \param RAMP_INDEX The index of the ramp
!> \param DUPLICATE_RAMP Optional logical parameter indicating if the ramp has already been processed

SUBROUTINE GET_RAMP_INDEX(ID,TYPE,RAMP_INDEX,DUPLICATE_RAMP)

USE GLOBAL_CONSTANTS, ONLY: N_RAMP,RAMP_ID,RAMP_TYPE
USE MEMORY_FUNCTIONS, ONLY: REALLOCATE_CHARACTER_ARRAY
CHARACTER(*), INTENT(IN) :: ID,TYPE
INTEGER, INTENT(OUT) :: RAMP_INDEX
INTEGER :: NR
LOGICAL, INTENT(IN), OPTIONAL :: DUPLICATE_RAMP

IF (ID=='null') THEN
   RAMP_INDEX = 0
   RETURN
ENDIF

IF (.NOT.PRESENT(DUPLICATE_RAMP)) THEN
   SEARCH: DO NR=1,N_RAMP
      IF (ID==RAMP_ID(NR)) THEN
         RAMP_INDEX = NR
         RETURN
      ENDIF
   ENDDO SEARCH
ENDIF

IF (N_RAMP>=SIZE(RAMP_ID)) THEN
   RAMP_ID   => REALLOCATE_CHARACTER_ARRAY(RAMP_ID,  LABEL_LENGTH,1,SIZE(RAMP_ID)+100)
   RAMP_TYPE => REALLOCATE_CHARACTER_ARRAY(RAMP_TYPE,LABEL_LENGTH,1,SIZE(RAMP_ID)+100)
ENDIF

N_RAMP                = N_RAMP + 1
RAMP_INDEX            = N_RAMP
RAMP_ID(RAMP_INDEX)   = ID
RAMP_TYPE(RAMP_INDEX) = TYPE

END SUBROUTINE GET_RAMP_INDEX


SUBROUTINE WRITE_SUMMARY_INFO(LU,INPUT_FILE_INCLUDED)

USE GLOBAL_CONSTANTS
USE COMP_FUNCTIONS, ONLY : GET_DATE
#ifdef WITHOUT_MPIF08
USE MPI
#else
USE MPI_F08
#endif
USE ISO_FORTRAN_ENV
CHARACTER(LABEL_LENGTH) :: DATE
INTEGER, INTENT(IN) :: LU
INTEGER :: MPIVERSION,MPISUBVERSION,MPILIBLENGTH,IERR
CHARACTER(LEN=MPI_MAX_LIBRARY_VERSION_STRING) :: MPILIBVERSION
LOGICAL, INTENT(IN) :: INPUT_FILE_INCLUDED

CALL GET_DATE(DATE)

WRITE(LU,'(/A/)')      ' Fire Dynamics Simulator'
WRITE(LU,'(A,A)')      ' Current Date     : ',TRIM(DATE)
WRITE(LU,'(A,A)')      ' Revision         : ',TRIM(GITHASH_PP)
WRITE(LU,'(A,A)')      ' Revision Date    : ',TRIM(GITDATE_PP)
WRITE(LU,'(A,A)')      ' Compiler         : ',TRIM(COMPILER_VERSION())
WRITE(LU,'(A,A/)')     ' Compilation Date : ',TRIM(BUILDDATE_PP)

IF (INPUT_FILE_INCLUDED) THEN
                   WRITE(LU,'(A,I0)')  ' Number of MPI Processes:  ',N_MPI_PROCESSES
   IF (USE_OPENMP) WRITE(LU,'(A,I0)')  ' Number of OpenMP Threads: ',OPENMP_AVAILABLE_THREADS
ENDIF

CALL MPI_GET_LIBRARY_VERSION(MPILIBVERSION,MPILIBLENGTH,IERR)
CALL MPI_GET_VERSION(MPIVERSION,MPISUBVERSION,IERR)
WRITE(LU,'(/A,I1,A,I1)') ' MPI version: ',MPIVERSION,'.',MPISUBVERSION
WRITE(LU,'(A,A)') ' MPI library version: ',TRIM(MPILIBVERSION)
#ifdef HYPRE_PP
WRITE(LU,'(A,A)') ' Hypre library version: ',TRIM(HYPRE_PP)
#else
WRITE(LU,'(A)') ' Hypre library: not used '
#endif
#ifdef SUNDIALS_PP
WRITE(LU,'(A,A)') ' Sundials library version: ',TRIM(SUNDIALS_PP)
#else
WRITE(LU,'(A)') ' Sundials library: not used'
#endif

END SUBROUTINE WRITE_SUMMARY_INFO


!> \brief Write VERBOSE diagnostic output

SUBROUTINE VERBOSE_PRINTOUT(DIAGNOSTIC_MESSAGE)

USE COMP_FUNCTIONS, ONLY : GET_DATE_ISO_8601
USE GLOBAL_CONSTANTS, ONLY: CPU_TIME_START,LU_ERR,MY_RANK
REAL :: CPUTIME
CHARACTER(*), INTENT(IN) :: DIAGNOSTIC_MESSAGE
CHARACTER(LABEL_LENGTH) :: DATE

CALL GET_DATE_ISO_8601(DATE)
CALL CPU_TIME(CPUTIME)
WRITE(LU_ERR,'(1X,A,I0,1X,A,A,F12.3,3X,A)') 'RANK=',MY_RANK,[CHARACTER(LEN=50)::DIAGNOSTIC_MESSAGE],&
                                            ' CPU Time:',CPUTIME-CPU_TIME_START,TRIM(DATE)

END SUBROUTINE VERBOSE_PRINTOUT


!> \brief Finds the device or control function assoicated with an input
!> \param NAME Namelist ID
!> \param CTRL_ID String containing name of control function.
!> \param DEVC_ID String containing name of device.
!> \param DEVICE_INDEX Integer index locating DEVC_ID in the array of devices.
!> \param CONTROL_INDEX Integer index locating CTRL_ID in the array of control functions.
!> \param INPUT_INDEX The current count of inputs of type NAME.

SUBROUTINE SEARCH_CONTROLLER(NAME,CTRL_ID,DEVC_ID,DEVICE_INDEX,CONTROL_INDEX,INPUT_INDEX)

USE DEVICE_VARIABLES, ONLY: DEVICE,N_DEVC
USE CONTROL_VARIABLES, ONLY: CONTROL,N_CTRL
USE COMP_FUNCTIONS, ONLY: SHUTDOWN
CHARACTER(MESSAGE_LENGTH) :: MESSAGE
CHARACTER(*), INTENT (IN) :: NAME,CTRL_ID,DEVC_ID
INTEGER :: I, DEVICE_INDEX,CONTROL_INDEX
INTEGER , INTENT(IN) :: INPUT_INDEX

! There cannot be both a device and controller for any given entity

IF (DEVC_ID /= 'null' .AND. CTRL_ID /='null') THEN
   WRITE(MESSAGE,'(A,A,1X,I3,A)')  'ERROR: ',TRIM(NAME),INPUT_INDEX,' has both a device (DEVC) and a control (CTRL) specified'
   CALL SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)
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


!> \brief Return the index of a material
!> \param ID Material identifier

INTEGER FUNCTION GET_MATL_INDEX(ID)

USE GLOBAL_CONSTANTS, ONLY: N_MATL
USE TYPES, ONLY: MATERIAL
CHARACTER(LABEL_LENGTH), INTENT(IN) :: ID
INTEGER :: N

DO N=1,N_MATL
   IF (MATERIAL(N)%ID/=ID) CYCLE
   GET_MATL_INDEX = N
   RETURN
ENDDO
GET_MATL_INDEX = 0

END FUNCTION GET_MATL_INDEX


!> \brief Return the index of a surface
!> \param ID Surface identifier

INTEGER FUNCTION GET_SURF_INDEX(ID)

USE GLOBAL_CONSTANTS, ONLY: N_SURF
USE TYPES, ONLY: SURFACE
CHARACTER(LABEL_LENGTH), INTENT(IN) :: ID
INTEGER :: N

DO N=1,N_SURF
   IF (SURFACE(N)%ID/=ID) CYCLE
   GET_SURF_INDEX = N
   RETURN
ENDDO
GET_SURF_INDEX = 0

END FUNCTION GET_SURF_INDEX


!> \brief Returns true if current MPI process MY_RANK controls mesh NM or its neighbors.
!> \param NM Mesh number.

LOGICAL FUNCTION PROCESS_MESH_NEIGHBORHOOD(NM)

USE MESH_VARIABLES
USE GLOBAL_CONSTANTS, ONLY: MY_RANK,PROCESS
INTEGER, INTENT(IN) :: NM
INTEGER :: N

PROCESS_MESH_NEIGHBORHOOD = .FALSE.
DO N=1,MESHES(NM)%N_NEIGHBORING_MESHES
   IF (MY_RANK==PROCESS(MESHES(NM)%NEIGHBORING_MESH(N))) PROCESS_MESH_NEIGHBORHOOD = .TRUE.
ENDDO
IF (MY_RANK==PROCESS(NM)) PROCESS_MESH_NEIGHBORHOOD = .TRUE.

END FUNCTION PROCESS_MESH_NEIGHBORHOOD


!> \brief Find the appropriate SPEC or SMIX index for the given SPEC_ID

SUBROUTINE GET_SPEC_OR_SMIX_INDEX(SPEC_ID,Y_INDX,Z_INDX)
USE TYPES, ONLY: SPECIES, SPECIES_MIXTURE
USE GLOBAL_CONSTANTS, ONLY: N_SPECIES, N_TRACKED_SPECIES
CHARACTER(*), INTENT(IN) :: SPEC_ID
INTEGER, INTENT(OUT) :: Y_INDX,Z_INDX
INTEGER :: NS

Y_INDX = -999
Z_INDX = -999

DO NS=1,N_SPECIES
   IF (TRIM(SPEC_ID)==TRIM(SPECIES(NS)%ID) .OR. TRIM(SPEC_ID)==TRIM(SPECIES(NS)%ALT_ID)) THEN
      Y_INDX = NS
      EXIT
    ENDIF
ENDDO

DO NS=1,N_TRACKED_SPECIES
   IF (TRIM(SPEC_ID)==TRIM(SPECIES_MIXTURE(NS)%ID) .OR. TRIM(SPEC_ID)==TRIM(SPECIES_MIXTURE(NS)%ALT_ID)) THEN
      Z_INDX = NS
      RETURN
   ENDIF
ENDDO

END SUBROUTINE GET_SPEC_OR_SMIX_INDEX

!> \brief Accumulate string MYSTR(LEN=STRING_SIZE) into allocatable string ACCSTR(LEN=ACCSTR_T_LEN).
!> \brief Note: An extra End of Line character is added at the end of MYSTR.
!> \param STRING_SIZE Length of MYSTR
!> \param MYSTR String to be added
!> \param ACCSTR Accumulated string
!> \param ACCSTR_T_LEN Total length of ACCSTR
!> \param ACCSTR_USE_LEN Length of ACCSTR used so far

SUBROUTINE ACCUMULATE_STRING(STRING_SIZE,MYSTR,ACCSTR,ACCSTR_T_LEN,ACCSTR_USE_LEN)
   INTEGER, INTENT(IN) :: STRING_SIZE
   CHARACTER(LEN=STRING_SIZE), INTENT(IN) :: MYSTR
   CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT) :: ACCSTR
   INTEGER, INTENT(INOUT) :: ACCSTR_T_LEN, ACCSTR_USE_LEN

   ! Local variables:
   INTEGER :: MYSTR_LEN
   CHARACTER(LEN=:), ALLOCATABLE :: AUXSTR

   ! Length of string to be added:
   MYSTR_LEN=LEN_TRIM(MYSTR)+1;
   ! Check ACCSTR_T_LEN (size of ACCSTR) and allocate more memory if needed: NOTE ACCSTR_T_LEN can be changed.
   IF(ACCSTR_T_LEN<ACCSTR_USE_LEN+MYSTR_LEN) THEN
      ACCSTR_T_LEN=ACCSTR_T_LEN+MYSTR_LEN+STRING_SIZE; ALLOCATE(CHARACTER(LEN=ACCSTR_T_LEN)::AUXSTR)
      IF(ALLOCATED(ACCSTR)) AUXSTR(1:ACCSTR_USE_LEN)=ACCSTR(1:ACCSTR_USE_LEN);
      CALL MOVE_ALLOC(AUXSTR,ACCSTR)
   ENDIF
   ! Write into ACCSTR:
   WRITE(ACCSTR(ACCSTR_USE_LEN+1:ACCSTR_USE_LEN+MYSTR_LEN),'(A,A)') TRIM(MYSTR),NEW_LINE(' ')
   ACCSTR_USE_LEN=ACCSTR_USE_LEN+MYSTR_LEN
END SUBROUTINE ACCUMULATE_STRING

END MODULE MISC_FUNCTIONS
   
