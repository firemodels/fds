PROGRAM TWOWAYCOUPLING

! Program to convert exposed surface elementos into ANSYS SURF152
USE IFPORT
IMPLICIT NONE

CHARACTER(256) FDSRUN,ANSYSRUN
CHARACTER(256) CHID,INPUT,ANSYS_STEP1,ANSYS_STEP2
CHARACTER(8) ARG
INTEGER :: OWNER_INDEX,I,J,K,STEPS
REAL TINT,TBEG,TEND,HC,TFINAL,DT
LOGICAL(4) COMMAND

!********TRY TO GET INFORMATION THROUGH ARGUMENTS*********
CALL GET_COMMAND_ARGUMENT(1,INPUT)
CALL GET_COMMAND_ARGUMENT(2,FDSRUN)
IF (LEN_TRIM(FDSRUN)/=0) THEN
   CHID=INPUT
   CALL GET_COMMAND_ARGUMENT(3,ANSYSRUN)
   CALL GET_COMMAND_ARGUMENT(4,ARG)
   READ(ARG,*) TFINAL
   CALL GET_COMMAND_ARGUMENT(5,ARG)
   READ(ARG,*) TINT
   CALL GET_COMMAND_ARGUMENT(6,ARG)
   READ(ARG,*) DT
   CALL GET_COMMAND_ARGUMENT(7,ARG)
   READ(ARG,*) HC
   CALL GET_COMMAND_ARGUMENT(8,ARG)
   READ(ARG,*) STEPS   
   CALL GET_COMMAND_ARGUMENT(9,ANSYS_STEP1)
   IF (STEPS==2) CALL GET_COMMAND_ARGUMENT(10,ANSYS_STEP2)   
!********TRY TO GET INFORMATION THROUGH INPUTFILE*********
ELSEIF (LEN_TRIM(INPUT)/=0) THEN
   OPEN(3,FILE=TRIM(INPUT),STATUS='OLD',FORM='FORMATTED')
   READ(3,'(a)') CHID
   READ(3,'(a)') FDSRUN
   READ(3,'(a)') ANSYSRUN 
   READ(3,*) TFINAL, TINT
   READ(3,*) DT
   READ(3,*) HC
   READ(3,*) STEPS
   READ(3,*) ANSYS_STEP1
   IF (STEPS==2) READ(3,*) ANSYS_STEP2
   CLOSE(3)
ENDIF
!***********GET INFORMATION THROUGH USER INPUT************
IF (LEN_TRIM(CHID)==0) THEN
   WRITE(6,*) ' Enter Job ID string (CHID):'
   READ(*,'(a)') CHID
   WRITE(6,*) ' Enter FDS executable path'
   READ(*,'(a)') FDSRUN
   WRITE(6,*) ' Enter ANSYS executable path'
   READ(*,'(a)') ANSYSRUN
   WRITE(6,*) ' Enter simulation ending time (s)'
   READ(*,*) TFINAL   
   WRITE(6,*) ' Enter boundary file time interval (s)'
   READ(*,*) TINT
   WRITE(6,*) ' Enter FDS-FEM communication time interval (s)'
   READ(*,*) DT
   WRITE(6,*) ' Specify the convective heat transfer coefficient (W/(m².K))'
   READ(*,*) HC
   WRITE(6,*) ' Specify the number of steps of the ANSYS simulation'
   READ(*,*) STEPS
   WRITE(6,*) ' Enter the name of the ANSYS model for step 1'
   READ(*,'(a)') ANSYS_STEP1
   IF (STEPS==2)THEN
      WRITE(6,*) ' Enter the name of the ANSYS model for step 2'
      READ(*,'(a)') ANSYS_STEP2
   ENDIF
ENDIF
!*******************************************************
OPEN(10,FILE=TRIM(CHID)//'.gc',ACTION='WRITE',FORM='UNFORMATTED')
OWNER_INDEX=0
WRITE(10) OWNER_INDEX
CLOSE(10)  

! start fds in backgound
IF (SYS==1) COMMAND=SYSTEMQQ (TRIM(FDSRUN)//' '//TRIM(CHID)//'.fds &')
IF (SYS==2) COMMAND=SYSTEMQQ ('start /b '//TRIM(FDSRUN)//' '//TRIM(CHID)//'.fds')

! main loop
MAIN_LOOP:DO J=0,INT(TFINAL/DT)-1
             TBEG=J*DT
             TEND=TBEG+DT
! wait for FDS flag
BOUNDARY_CHECK_LOOP: DO I=1,240 ! 300 loops x 15 sec = 60min
   OPEN(10,FILE=TRIM(CHID)//'.gc',ACTION='READ',FORM='UNFORMATTED')
   READ(10) OWNER_INDEX
   IF (OWNER_INDEX/=2) THEN
      CLOSE(10)
      IF (I==1) THEN
         WRITE (6,'(A)') ' COUPLING: waiting for FDS results ... '
         K=1
      ELSE 
         IF (K==4) THEN
            K=0
            WRITE(6,'(A,I2,A)') ' COUPLING: waiting ... ', (I-1)/4,' min'
         ELSE
            K=K+1
         ENDIF
      ENDIF
      CALL SLEEP(15)
      IF (I==240) THEN
         WRITE (6,*) ' ERROR: .be FILE WAS NOT UPDATED BY FDS'
         STOP
      ENDIF
   ELSEIF (OWNER_INDEX==2) THEN
      CLOSE(10)
      EXIT BOUNDARY_CHECK_LOOP
   ENDIF
ENDDO BOUNDARY_CHECK_LOOP

! extract data from .be file
CALL INTERFACE (CHID,TBEG,TEND,TINT,HC,DT,STEPS,ANSYS_STEP1,ANSYS_STEP2)

! run ansys from TBEG to TEND
COMMAND=SYSTEMQQ (TRIM(ANSYSRUN)//' -j '//TRIM(CHID)//'_ansys -i run_case.dat')

! write .gc file to update geometry at FDS
CALL WRITE_GC (CHID,TBEG,TEND)

END DO MAIN_LOOP
END PROGRAM