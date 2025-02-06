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
USE MPI_F08
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
USE QREF_DATA
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

I_HRRPUA = MIN(14,INT(HRRPUA*0.000005_EB))
F_HRRPUA = MIN(15._EB,HRRPUA*0.000005_EB) - REAL(I_HRRPUA,EB)
IF (F_HRRPUA < TWO_EPSILON_EB) THEN
   Q_REF_FIT = Q_INC
   RETURN
ENDIF
I_HOC = MIN(4,MAX(1,INT(HOC*1.E-7_EB)))
F_HOC = MIN(5._EB,MAX(1._EB,HOC*1.E-7_EB)) - REAL(I_HOC,EB)
I_Y_S = MIN(9,INT(Y_S*50._EB))
F_Y_S = MIN(10._EB,Y_S*50._EB) - REAL(I_Y_S,EB)
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
USE MPI_F08
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

IF ((X2 <= X0-RAD) .OR. (X1 >= X0+RAD) .OR. (Y2 <= Y0-RAD) .OR. (Y1 >= Y0+RAD) .OR. RAD<TWO_EPSILON_EB) THEN
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
   ! No corners either full area or circle minus one or more chords
   CASE (0)
      ! First four cases are more than half the cirlce outside the rectangle.
      ! Can only have one edge where this is the case without having a corner inside. This edge makes a chord.
      ! Intersection area is the area of the chord.
      IF (X2 <= X0 .AND. X2 > X0-RAD) THEN
         THETA = 2._EB*ACOS((X0-X2)/RAD)         
         CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA))
         RETURN
      ENDIF
      IF (Y2 <= Y0 .AND. Y2 > Y0-RAD) THEN
         THETA = 2._EB*ACOS((Y0-Y2)/RAD)         
         CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA))
         RETURN
      ENDIF
      IF (X1 >= X0 .AND. X1 < X0+RAD) THEN
         THETA = 2._EB*ACOS((X1-X0)/RAD)         
         CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA))
         RETURN
      ENDIF
      IF (Y1 >= Y0 .AND. Y1 < Y0+RAD) THEN
         THETA = 2._EB*ACOS((Y1-Y0)/RAD)         
         CIRCLE_CELL_INTERSECTION_AREA = 0.5_EB*R2*(THETA-SIN(THETA))
         RETURN
      ENDIF
      ! Remaining cases is where one or more rectangle edges are inside but the corners are outside
      ! With no corner inside each internal edge is a chord
      ! Intersection area is the circe area minus the area of each chord
      CIRCLE_CELL_INTERSECTION_AREA = PI*R2
      IF (X1 > X0 - RAD) THEN
         THETA = 2._EB*ACOS((X0-X1)/RAD)
         CIRCLE_CELL_INTERSECTION_AREA = CIRCLE_CELL_INTERSECTION_AREA - 0.5_EB*R2*(THETA-SIN(THETA))
      ENDIF
      IF (Y1 > Y0 - RAD) THEN
         THETA = 2._EB*ACOS((Y0-Y1)/RAD)
         CIRCLE_CELL_INTERSECTION_AREA = CIRCLE_CELL_INTERSECTION_AREA - 0.5_EB*R2*(THETA-SIN(THETA))
      ENDIF
      IF (X2 < X0 + RAD) THEN
         THETA = 2._EB*ACOS((X2-X0)/RAD)
         CIRCLE_CELL_INTERSECTION_AREA = CIRCLE_CELL_INTERSECTION_AREA - 0.5_EB*R2*(THETA-SIN(THETA))
      ENDIF
      IF (Y2 < Y0 + RAD) THEN
         THETA = 2._EB*ACOS((Y2-Y0)/RAD)
         CIRCLE_CELL_INTERSECTION_AREA = CIRCLE_CELL_INTERSECTION_AREA - 0.5_EB*R2*(THETA-SIN(THETA))
      ENDIF
      RETURN
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
INTEGER :: K,NZ
REAL(EB) :: X_MIN,X_MAX,Y_MIN,Y_MAX,Z_MIN,Z_MAX,DZ,Z1,Z2,&
            OUTER_BASE_AREA,INNER_BASE_AREA,OUTER_TOP_AREA,INNER_TOP_AREA,HH_IN,&
            RR0_BASE,RRI_BASE,RR0_TOP,RRI_TOP
TYPE (MESH_TYPE), POINTER :: M

CONE_MESH_INTERSECTION_VOLUME = 0._EB
M => MESHES(NM)

X_MIN = MAX(M%XS,X0-RR0)
X_MAX = MIN(M%XF,X0+RR0)
Y_MIN = MAX(M%YS,Y0-RR0)
Y_MAX = MIN(M%YF,Y0+RR0)
Z_MIN = MAX(M%ZS,Z0)
Z_MAX = MIN(M%ZF,Z0+HH0)
HH_IN = Z_MAX-Z_MIN

! Case 1: volume fully outside the mesh
IF (X_MAX<=X_MIN .OR. Y_MAX<=Y_MIN .OR. Z_MAX<=Z_MIN) RETURN

! Case 2: base area fully inside the mesh
IF ((X_MIN-M%XS)>TWO_EPSILON_EB .AND. (M%XF-X_MAX)>TWO_EPSILON_EB &
   .AND. (Y_MIN-M%YS)>TWO_EPSILON_EB .AND. (M%YF-Y_MAX)>TWO_EPSILON_EB) THEN
   IF (G_FACTOR==0) THEN
      ! Subtract inner from outer cylinder volume
      CONE_MESH_INTERSECTION_VOLUME = PI*(RR0**2._EB-RRI**2._EB)*HH_IN
   ELSEIF (G_FACTOR==1) THEN
      RR0_BASE = RR0*(1-(Z_MIN-Z0)/HH0)
      RRI_BASE = RRI*(1-(Z_MIN-Z0)/HH0)
      RR0_TOP = RR0*(1-(Z_MAX-Z0)/HH0)
      RRI_TOP = RRI*(1-(Z_MAX-Z0)/HH0)
      ! Subtract inner from outer truncated cone volume
      CONE_MESH_INTERSECTION_VOLUME = PI*HH_IN/3._EB*&
         (RR0_BASE**2._EB+RR0_BASE*RR0_TOP+RR0_TOP**2._EB-RRI_BASE**2._EB-RRI_BASE*RRI_TOP-RRI_TOP**2._EB)
   ENDIF
   RETURN
ENDIF

! Case 3: partial hollow cylinder intersection volume
IF (G_FACTOR==0) THEN
   OUTER_BASE_AREA = CIRCLE_CELL_INTERSECTION_AREA(X0,Y0,RR0,M%XS,M%XF,M%YS,M%YF)
   INNER_BASE_AREA = CIRCLE_CELL_INTERSECTION_AREA(X0,Y0,RRI,M%XS,M%XF,M%YS,M%YF)
   CONE_MESH_INTERSECTION_VOLUME = (OUTER_BASE_AREA-INNER_BASE_AREA)*HH_IN
   RETURN
ENDIF

! Case 4: partial hollow cone intersection volume (trapezoidal integration)
NZ = 100
DZ = HH_IN/REAL(NZ,EB)

DO K=1,NZ
   Z1 = Z_MIN + (K-1)*DZ
   Z2 = Z1 + DZ
   RR0_BASE = RR0*(1-(Z1-Z0)/HH0)
   RRI_BASE = RRI*(1-(Z1-Z0)/HH0)
   RR0_TOP = RR0*(1-(Z2-Z0)/HH0)
   RRI_TOP = RRI*(1-(Z2-Z0)/HH0)
   OUTER_BASE_AREA = CIRCLE_CELL_INTERSECTION_AREA(X0,Y0,RR0_BASE,M%XS,M%XF,M%YS,M%YF)
   INNER_BASE_AREA = CIRCLE_CELL_INTERSECTION_AREA(X0,Y0,RRI_BASE,M%XS,M%XF,M%YS,M%YF)
   OUTER_TOP_AREA = CIRCLE_CELL_INTERSECTION_AREA(X0,Y0,RR0_TOP,M%XS,M%XF,M%YS,M%YF)
   INNER_TOP_AREA = CIRCLE_CELL_INTERSECTION_AREA(X0,Y0,RRI_TOP,M%XS,M%XF,M%YS,M%YF)

   CONE_MESH_INTERSECTION_VOLUME = CONE_MESH_INTERSECTION_VOLUME + &
      DZ*0.5_EB*(OUTER_BASE_AREA+OUTER_TOP_AREA-INNER_BASE_AREA-INNER_TOP_AREA)
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
  USE GLOBAL_CONSTANTS, ONLY : LU_ERR, MY_RANK, N_MPI_PROCESSES, VERBOSE
  INTEGER :: THREAD_ID

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(THREAD_ID)

  THREAD_ID = 0
  !$ THREAD_ID = OMP_GET_THREAD_NUM()

  !$OMP CRITICAL
  IF (USE_OPENMP .AND. OPENMP_AVAILABLE_THREADS>1 .AND. VERBOSE) WRITE(LU_ERR,91) " OpenMP thread ",THREAD_ID," of ",&
     OPENMP_AVAILABLE_THREADS-1," assigned to MPI process ",MY_RANK," of ",N_MPI_PROCESSES-1
  IF (.NOT.USE_OPENMP .AND. VERBOSE) WRITE(LU_ERR,92) " MPI process ",MY_RANK," of ",N_MPI_PROCESSES-1
  !$OMP END CRITICAL

  !$OMP END PARALLEL

  91 FORMAT(A,I3,A,I3,A,I6,A,I6)
  92 FORMAT(A,I6,A,I6)

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
USE MPI_F08
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
   
