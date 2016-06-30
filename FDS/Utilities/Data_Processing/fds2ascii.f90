#ifndef GITHASH_PP
#define GITHASH_PP "unknown"
#endif
#ifndef GITDATE_PP
#define GITDATE_PP "unknown"
#endif
#ifndef BUILDDATE_PP
#define BUILDDATE_PP "unknown"
#endif
PROGRAM FDS2ASCII

! Program to convert various FDS output files to ASCII

IMPLICIT NONE

INTERFACE
SUBROUTINE PARSE(BUFFER,SB_TOKS,SE_TOKS,N_TOKS)
IMPLICIT NONE
CHARACTER(*), INTENT(INOUT) :: BUFFER
INTEGER, DIMENSION(*), INTENT(OUT) :: SB_TOKS, SE_TOKS
INTEGER, INTENT(OUT) :: N_TOKS
END SUBROUTINE PARSE
END INTERFACE

CHARACTER(255), PARAMETER :: f2aversion='2.1.0'
INTEGER, PARAMETER :: FB = SELECTED_REAL_KIND(6)
INTEGER, PARAMETER :: FILE_DIM = 500
INTEGER :: IERR, NMESHES, NM, NOC, I, J, K, L
INTEGER :: IDUM, IFILE, NSAM, NV, MV
INTEGER :: I1, I2, J1, J2, K1, K2, I3, J3, K3
INTEGER :: I10, I20, J10, J20, K10, K20
INTEGER :: NCOUNT, IOR_INPUT, NPATCH, IJBAR, JKBAR
INTEGER :: II, NXP, NYP, NZP, N
REAL(FB) :: XS, XF, YS, YF, ZS, ZF, TIME
REAL(FB) :: D1, D2, D3, D4, TBEG, TEND
CHARACTER(256) :: AUTO_SLICE_LABEL
INTEGER :: AUTO_SLICE_FLAG, N_AUTO_SLICES, IAUTO, IZTEMP
INTEGER, DIMENSION(1) :: IZMIN1
INTEGER :: IZMIN
INTEGER, DIMENSION(FILE_DIM) :: AUTO_SLICE_LISTS
REAL(FB), DIMENSION(FILE_DIM) :: AUTO_SLICE_Z
REAL(FB) :: AX1, AY1, AZ1, AZ2
REAL(FB) :: EPS, ZMIN
 
TYPE MESH_TYPE
   REAL(FB), POINTER, DIMENSION(:) :: X,Y,Z
   REAL(FB) :: D1,D2,D3,D4
   INTEGER :: IBAR,JBAR,KBAR,IERR,NXP,NYP,NZP
END TYPE MESH_TYPE
 
TYPE (MESH_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: MESH
TYPE(MESH_TYPE), POINTER :: M
 
REAL(FB), DIMENSION(60) :: SUM
INTEGER, ALLOCATABLE, DIMENSION(:) :: IOR,I1B,I2B,J1B,J2B,K1B,K2B
INTEGER, ALLOCATABLE, DIMENSION(:) :: AUTO_SLICES
REAL(FB), ALLOCATABLE, DIMENSION(:,:,:,:) :: Q
REAL(FB), ALLOCATABLE, DIMENSION(:,:,:) :: F
LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: ALREADY_USED
LOGICAL :: NEW_PLOT3D=.TRUE.
INTEGER IOR_SLCF
CHARACTER(2) :: ANS
CHARACTER(4) :: CHOICE
CHARACTER(30) :: UNITJUNK
CHARACTER(256) GRIDFILE,QNAME,CHID,QFILE,JUNK,FRMT,OUTFILE,SLCF_LABEL_DUMMY
CHARACTER(256), DIMENSION(FILE_DIM) :: PL3D_FILE,SLCF_FILE,BNDF_FILE,SLCF_TEXT,BNDF_TEXT,SLCF_UNIT,BNDF_UNIT,SLCF_LABEL
CHARACTER(256), DIMENSION(FILE_DIM) :: SLICE_LABELS
CHARACTER(20), DIMENSION(FILE_DIM) :: BNDF_TYPE
CHARACTER(20) :: BNDF_TYPE_CHOSEN
INTEGER :: NSLICE_LABELS, SLICE_EXIST
REAL(FB), DIMENSION(FILE_DIM) :: X1, X2, Y1, Y2, Z1, Z2
INTEGER,  DIMENSION(60) :: IB,IS
INTEGER,  DIMENSION(FILE_DIM) :: PL3D_MESH,SLCF_MESH,BNDF_MESH
REAL(FB), DIMENSION(FILE_DIM) :: PL3D_TIME
LOGICAL, DIMENSION(FILE_DIM) :: FILE_EXISTS
INTEGER :: RCODE
INTEGER :: NFILES_EXIST
REAL(FB) :: ZOFFSET
INTEGER :: ZOFFSET_FLAG
LOGICAL :: EXISTS
INTEGER :: LU_IN, NARGS, IARG, LENSTRING, BATCHMODE
CHARACTER(256) :: BUFFER, FILEIN
INTEGER, DIMENSION(256) :: SB_TOKS, SE_TOKS
INTEGER :: N_TOKS
INTEGER :: ERROR_STATUS
CHARACTER(256) :: ARG

! Optional arguments in the call sequence

NARGS = IARGC()
IF(NARGS.GT.0)THEN
  DO I = 1, NARGS
    CALL GETARG(I,ARG)
    IF(ARG.EQ."-h".OR.ARG.EQ."-H")THEN
      CALL USAGE2(f2aversion)
      STOP
    ENDIF
    IF(ARG.EQ."-v".OR.ARG.EQ."-V")THEN
      CALL VERSION2(f2aversion)
      STOP
    ENDIF
  END DO
ENDIF

! Set a few default values

ZOFFSET=0.0
ZOFFSET_FLAG=0
EPS=0.001
BATCHMODE=0
FILEIN='STDIN'

! Parse command line arguments
      
IF (FILEIN.EQ.'STDIN') THEN
   LU_IN=5
ELSE
   LU_IN=7
   INQUIRE(FILE=FILEIN,EXIST=EXISTS)
   IF (.NOT.EXISTS) THEN
      WRITE(6,*)"*** Fatal error: The file: ",TRIM(FILEIN), " does not exist"
      STOP
   ENDIF
   OPEN(LU_IN,FILE=FILEIN,STATUS='OLD',FORM='FORMATTED')
ENDIF
      
IF (BATCHMODE==0) WRITE(6,*) ' Enter Job ID string (CHID):'

READ(LU_IN,'(a)') CHID
 
! Check to see if the .smv file exists
 
GRIDFILE = TRIM(CHID)//'.smv'
INQUIRE(FILE=GRIDFILE,EXIST=EXISTS)
IF (.NOT.EXISTS) THEN
   WRITE(6,*)"*** Fatal error: The file: ",TRIM(GRIDFILE)," does not exist"
   STOP
ENDIF

! Open the .smv file

OPEN(11,FILE=GRIDFILE,STATUS='OLD',FORM='FORMATTED')
 
 
! Determine the number of meshes

REWIND(11)
 
CALL SEARCH('NMESHES',7,11,IERR)
IF (IERR.EQ.1) THEN
   WRITE(6,*) ' WARNING: Assuming 1 mesh'
   NMESHES = 1
ELSE
   READ(11,*) NMESHES
ENDIF

ALLOCATE(MESH(NMESHES))
 
! Get the coordinates of the meshes

REWIND(11)
 
READ_SMV: DO NM=1,NMESHES
   M=>MESH(NM)
   CALL SEARCH('GRID',4,11,IERR)
   READ(11,*) M%IBAR,M%JBAR,M%KBAR
   ALLOCATE(M%X(0:M%IBAR))
   ALLOCATE(M%Y(0:M%JBAR))
   ALLOCATE(M%Z(0:M%KBAR))
 
   CALL SEARCH('TRNX',4,11,IERR)
   READ(11,*) NOC
   DO I=1,NOC
      READ(11,*)
      ENDDO
      DO I=0,M%IBAR
      READ(11,*) IDUM,M%X(I)
   ENDDO
 
   CALL SEARCH('TRNY',4,11,IERR)
   READ(11,*) NOC
   DO I=1,NOC
      READ(11,*)
      ENDDO
      DO J=0,M%JBAR
      READ(11,*) IDUM,M%Y(J)
   ENDDO
 
   CALL SEARCH('TRNZ',4,11,IERR)
   READ(11,*) NOC
   DO I=1,NOC
      READ(11,*)
      ENDDO
      DO K=0,M%KBAR
      READ(11,*) IDUM,M%Z(K)
   ENDDO
 
ENDDO READ_SMV
 
! Determine what kind of file to work with

IF (BATCHMODE.EQ.0) THEN
   write(6,*) ' What type of file to parse?'
   write(6,*) ' PL3D file? Enter 1'
   write(6,*) ' SLCF file? Enter 2'
   write(6,*) ' BNDF file? Enter 3'
ENDIF

READ(LU_IN,*) IFILE
 
! Get sampling factor for data

IF (BATCHMODE.EQ.0) THEN
   write(6,*) ' Enter Sampling Factor for Data?'
   write(6,*) ' (1 for all data, 2 for every other point, etc.)'
ENDIF

READ(LU_IN,*) NSAM
 
! Determine whether to limit domain size

! ANS(1:1) may be 'y', 'n' or 'z' (or upper case equivalents)
! ANS(2:2) may be 'a' or ' '
 
! y - domain size is limited
! n - domain size is not limited
! z - domain size is not limited z levels are offset by zoffset
! a - slice files are selected based on slice file type and location

IF (BATCHMODE.EQ.0)THEN
   write(6,*) ' Domain selection:'
   write(6,*) '   y - domain size is limited'
   write(6,*) '   n - domain size is not limited'
   write(6,*) '   z - domain size is not limited and z levels are offset'
   write(6,*) '   ya, na or za - slice files are selected based on type and location.'
   write(6,*) '       The y, n, z prefix are defined as before.  '
ENDIF
READ(LU_IN,'(A)') ANS
CALL TOUPPER(ANS,ANS)
IF (ANS(1:1).EQ.'Y') THEN
   IF (BATCHMODE.EQ.0)write(6,*) ' Enter min/max x, y and z'
   READ(LU_IN,*) XS,XF,YS,YF,ZS,ZF
ELSE
   XS = -100000.
   XF =  100000.
   YS = -100000.
   YF =  100000.
   ZS = -100000.
   ZF =  100000.
ENDIF
       
! If ANS is z or Z then subtract zoffset from z data (for multi-level cases)

IF (ANS(1:1).EQ.'Z') THEN
   ZOFFSET_FLAG=1
ELSE
   ZOFFSET_FLAG=0
ENDIF

IF(LEN(TRIM(ANS)).GT.1)THEN
   IF (ANS(2:2).EQ.'A')THEN
      AUTO_SLICE_FLAG=1
   ELSE
      AUTO_SLICE_FLAG=0
   ENDIF
ELSE
   AUTO_SLICE_FLAG=0
ENDIF
 
! Select which kind of file to work on

EXTRA_PLOT3D_FILES: DO
 
   FILETYPE: SELECT CASE(IFILE)
 
      CASE(1) FILETYPE  ! Read a PLOT3D file and print out the values in ASCII text

      PL3D_MESH = 1
      REWIND(11)
      NFILES_EXIST=0
      SEARCH_PL3D: DO I=1,FILE_DIM
      CALL SEARCH('PL3D',4,11,IERR)
      IF (IERR.EQ.1) EXIT SEARCH_PL3D
      BACKSPACE(11) 
      READ(11,*) JUNK,PL3D_TIME(I),PL3D_MESH(I)
      READ(11,'(A)') PL3D_FILE(I)
      DO II=1,5
      IS(II) = II
      READ(11,'(A)') SLCF_TEXT(II)
      READ(11,*)
      READ(11,'(A)') SLCF_UNIT(II)
      ENDDO
      OPEN(12,FILE=TRIM(PL3D_FILE(I)),FORM='UNFORMATTED',STATUS='OLD',IOSTAT=RCODE)
      CLOSE(12)
      FILE_EXISTS(I)=.TRUE.
      IF(RCODE.NE.0)THEN
        FILE_EXISTS(I)=.FALSE.
        CYCLE
      ENDIF
      NFILES_EXIST = NFILES_EXIST + 1
      IF(BATCHMODE.EQ.0)THEN
        write(6,'(I3,3X,A,A,I2,A,F5.0)') I,TRIM(PL3D_FILE(I)),', MESH ',PL3D_MESH(I),', Time: ',PL3D_TIME(I)
      ENDIF
      ENDDO SEARCH_PL3D
 
      IF(NFILES_EXIST.EQ.0)THEN
        write(6,*)"There are no PLOT3D files to convert"
        STOP
      ENDIF
      
      IF(BATCHMODE.EQ.0)THEN
       write(6,'(/A)') ' Choose PLOT3D .q file by its index: (0 to end)'
      ENDIF
      READ(LU_IN,*) I
      IF(I.NE.0)THEN
        IF(.NOT.FILE_EXISTS(I))THEN
          write(6,*)"The file, ",TRIM(PL3D_FILE(I)),", does not exist"
          CYCLE
        ENDIF
      ENDIF
      IF (I.EQ.0) STOP
 
      QFILE = PL3D_FILE(I)
      NM = PL3D_MESH(I)
      M=>MESH(NM)
 
      I1 = 0
      I2 = M%IBAR
      J1 = 0
      J2 = M%JBAR
      K1 = 0
      K2 = M%KBAR
      NV = 5
      IF (NEW_PLOT3D) ALLOCATE(Q(0:M%IBAR,0:M%JBAR,0:M%KBAR,5))
 
      OPEN(12,FILE=QFILE,FORM='UNFORMATTED',STATUS='OLD')
      IF(BATCHMODE.EQ.0)write(6,*) ' Reading PLOT3D data file...'
      READ(12) NXP,NYP,NZP
      READ(12) D1,D2,D3,D4
      READ(12) ((((Q(I,J,K,N),I=0,M%IBAR),J=0,M%JBAR),K=0,M%KBAR),N=1,5)
      CLOSE(12)
 

      CASE(2) FILETYPE  ! Slice files
 

         NSLICE_LABELS=0
         SLCF_LABEL = 'null'

         IF (ZOFFSET_FLAG.EQ.0) THEN
            IF (BATCHMODE.EQ.0) write(6,*) ' Enter starting and ending time for averaging (s)'
            READ(LU_IN,*) TBEG,TEND
         ELSE
            IF (BATCHMODE.EQ.0) write(6,*) ' Enter starting and ending time for averaging (s) and zoffset (m)'
            READ(LU_IN,*) TBEG,TEND,ZOFFSET
         ENDIF

         SLCF_MESH = 1

         REWIND(11)
         NFILES_EXIST=0

         SEARCH_SLCF: DO I=1,FILE_DIM
            CALL SEARCH2('SLCF',4,'SLCC',4,11,IERR,CHOICE)
            IF (IERR.EQ.1) EXIT SEARCH_SLCF
            BACKSPACE(11)
            READ(11,*) JUNK,SLCF_MESH(I)
            READ(11,'(A)') SLCF_FILE(I)
            READ(11,'(A)') SLCF_TEXT(I)
   
            ! Create unique list of slice types      
   
            SLICE_EXIST=0
            DO II=1, NSLICE_LABELS
              IF(TRIM(SLCF_TEXT(I)).EQ.SLICE_LABELS(II))THEN
                SLICE_EXIST=1
                EXIT
              ENDIF
            ENDDO
            IF (SLICE_EXIST.EQ.0) THEN
               NSLICE_LABELS=NSLICE_LABELS+1
               SLICE_LABELS(NSLICE_LABELS)=TRIM(SLCF_TEXT(I))
            ENDIF
   
            READ(11,*) 
            READ(11,'(A)') SLCF_UNIT(I)
            OPEN(12,FILE=SLCF_FILE(I),FORM='UNFORMATTED',STATUS='OLD', IOSTAT=RCODE)
            IF (RCODE.NE.0) THEN
               CLOSE(12)
               CYCLE
            ENDIF
   
            NFILES_EXIST=NFILES_EXIST+1
   
            READ(12) UNITJUNK
            READ(12) UNITJUNK
            READ(12) UNITJUNK
            READ(12) I1,I2,J1,J2,K1,K2
            CLOSE(12)
   
            NM=SLCF_MESH(I)
            M=>MESH(NM)
            X1(I)=M%X(I1)
            X2(I)=M%X(I2)
            Y1(I)=M%Y(J1)
            Y2(I)=M%Y(J2)
            Z1(I)=M%Z(K1)
            Z2(I)=M%Z(K2)
   
            IF(AUTO_SLICE_FLAG.EQ.0)THEN
               IF (BATCHMODE.EQ.0)THEN
                  write(6,'(I3,1x,A,1x,A)')I,TRIM(SLCF_TEXT(I)),TRIM(SLCF_FILE(I))
                  write(6,'(3x,A,6(1x,f8.2))')'slice bounds:',M%X(I1),M%X(I2),M%Y(J1),M%Y(J2),M%Z(K1),M%Z(K2)
               ENDIF
            ENDIF
   
         ENDDO SEARCH_SLCF
    
         IF (NFILES_EXIST.EQ.0)THEN
            WRITE(6,*)"There are no slice files to convert"
            STOP
         ENDIF

         AUTO_SLICE_IF: IF (AUTO_SLICE_FLAG.EQ.0) THEN

            IF(BATCHMODE.EQ.0) write(6,*)'How many variables to read:'
            READ(LU_IN,*) NV
            N_AUTO_SLICES=1

         ELSE AUTO_SLICE_IF

            N_AUTO_SLICES=0
            AUTO_SLICE_FLAG=1
            IF (BATCHMODE.EQ.0)THEN
               write(6,*)' Enter slice file type index'
               DO II=1, NSLICE_LABELS
                  write(6,'(2x,I3,1x,A)')II,TRIM(SLICE_LABELS(II))
               ENDDO
            ENDIF
            READ(LU_IN,'(I3)')II

            AUTO_SLICE_LABEL=TRIM(SLICE_LABELS(II))
            boundloop: DO II=1,NFILES_EXIST
               IF(TRIM(AUTO_SLICE_LABEL).NE.TRIM(SLCF_TEXT(II))) CYCLE boundloop
               IF(X2(II)-X1(II).GT.EPS)CYCLE boundloop
               IF(Y2(II)-Y1(II).GT.EPS)CYCLE boundloop
               IF(BATCHMODE.EQ.0) write(6,'(4(a,f8.3))')"x=",X1(II)," y=",Y1(II)," zmin=",Z1(II)," zmax=",Z2(II)
            ENDDO boundloop
        
            IF (BATCHMODE.EQ.0) write(6,*)"Enter 1D SLICE location (x y zmin zmax)"
            READ(LU_IN,*)AX1, AY1, AZ1, AZ2
            EPS=0.001

            ! Create list of slices that match requested type and location also record the elevation (z) of each 
            ! slice so that they may be later output in increasing order

            AUTO_LOOP1: DO I=1,NFILES_EXIST
               IF (TRIM(SLCF_TEXT(I)).EQ.TRIM(AUTO_SLICE_LABEL))THEN
                  IF(X2(I)-X1(I).GT.EPS)CYCLE
                  IF(Y2(I)-Y1(I).GT.EPS)CYCLE
                  IF(AX1+EPS.LT.X1(I))CYCLE
                  IF(AX1-EPS.GT.X2(I))CYCLE
                  IF(AY1+EPS.LT.Y1(I))CYCLE
                  IF(AY1-EPS.GT.Y2(I))CYCLE
                  N_AUTO_SLICES = N_AUTO_SLICES + 1
                  AUTO_SLICE_LISTS(N_AUTO_SLICES)=I
                  AUTO_SLICE_Z(N_AUTO_SLICES) = Z1(I)
               ENDIF
            ENDDO AUTO_LOOP1

         ENDIF AUTO_SLICE_IF
 
         SUM = 0.
 
         IF (AUTO_SLICE_FLAG.EQ.1) THEN
            ALLOCATE(Q(0:1,0:1,0:1,1))
            ALLOCATE(F(0:1,0:1,0:1))
            NV = 1

            ! Sort slice list        

            DO I = 1, N_AUTO_SLICES
               IZMIN1= MINLOC(AUTO_SLICE_Z(I:N_AUTO_SLICES))
               IZMIN = IZMIN1(1) + I - 1
               IF (IZMIN.NE.I)THEN
                  ZMIN=AUTO_SLICE_Z(IZMIN)
                  AUTO_SLICE_Z(IZMIN)=AUTO_SLICE_Z(I)
                  AUTO_SLICE_Z(I) = ZMIN
                  IZTEMP = AUTO_SLICE_LISTS(IZMIN)
                  AUTO_SLICE_LISTS(IZMIN)=AUTO_SLICE_LISTS(I)
                  AUTO_SLICE_LISTS(I)=IZTEMP
               ENDIF
            ENDDO
         ELSE
            N_AUTO_SLICES=1
         ENDIF

      AUTO_LIST: DO IAUTO=1, N_AUTO_SLICES
      VARLOOP: DO MV=1,NV
 
      IF (AUTO_SLICE_FLAG.EQ.0)THEN
         SLCF_LABEL_DUMMY=' '
 ! note:  The following two lines of FORTRAN was my attempt at getting EOR to work.
 !        I got compile errors. (My "research" said that you can only use EOR with non-advancing IO)
 !  see http://www.pcc.qub.ac.uk/tec/courses/f77tof90/stu-notes/f90studentMIF_7.html for an example
         !IF (BATCHMODE.EQ.0) write(6,'(A,I2)',ADVANCE='NO') ' Enter index for variable',MV
         !READ(LU_IN,*,EOR=200,ADVANCE='NO') I,SLCF_LABEL_DUMMY
         IF (BATCHMODE.EQ.0) write(6,'(A,I2)') ' Enter index for variable',MV
         READ(LU_IN,'(A)',IOSTAT=ERROR_STATUS) BUFFER
         IF(ERROR_STATUS.NE.0)THEN
             WRITE(6,*)"*** fatal error: read of variable index failed"
             STOP
         ENDIF
         CALL PARSE(BUFFER,SB_TOKS,SE_TOKS,N_TOKS)
         IF(N_TOKS.GE.1)THEN
            READ(BUFFER(SB_TOKS(1):SE_TOKS(1)),*)I
            IF(N_TOKS.GT.1)SLCF_LABEL_DUMMY=BUFFER(SB_TOKS(2):SE_TOKS(N_TOKS))
         ELSE
            WRITE(6,*)"*** fatal error: index for variable ",MV," not entered"
            STOP
         ENDIF
         200 CONTINUE
         SLCF_LABEL(I) = SLCF_LABEL_DUMMY
         IF (SLCF_LABEL(I)=='null' .OR. SLCF_LABEL(I)==' ') SLCF_LABEL(I) = SLCF_TEXT(I)
      ELSE
         I=AUTO_SLICE_LISTS(IAUTO)
      ENDIF

      IS(MV) = I
      QFILE = SLCF_FILE(I)
 
      IF (MV.EQ.1) THEN
      NM = SLCF_MESH(I)
      M=>MESH(NM)
      IF(AUTO_SLICE_FLAG.EQ.1)THEN
        DEALLOCATE(Q)
        DEALLOCATE(F)
      ENDIF
      ALLOCATE(Q(0:M%IBAR,0:M%JBAR,0:M%KBAR,NV))
      ALLOCATE(F(0:M%IBAR,0:M%JBAR,0:M%KBAR))
      F = 0.
      Q = 0.
      ELSE
      IF (SLCF_MESH(I).NE.NM) THEN
        WRITE(6,*) ' ERROR: All slices must have the same mesh'
        STOP
        ENDIF
      ENDIF
 
      OPEN(12,FILE=QFILE,FORM='UNFORMATTED',STATUS='OLD')
 
      READ(12)
      READ(12)
      READ(12)
      READ(12) I1,I2,J1,J2,K1,K2                    
 
      IF (MV.EQ.1) THEN
         I10=I1 ; I20=I2 ; J10=J1 ; J20=J2 ; K10=K1 ; K20=K2
         IF (I1.EQ.I2) IOR_SLCF = 1
         IF (J1.EQ.J2) IOR_SLCF = 2
         IF (K1.EQ.K2) IOR_SLCF = 3
      ELSE
         IF (I1.EQ.I2 .AND. I10.EQ.I20) THEN
            I1=I10
            I2=I20
            ENDIF
         IF (J1.EQ.J2 .AND. J10.EQ.J20) THEN
            J1=J10
            J2=J20
            ENDIF
         IF (K1.EQ.K2 .AND. K10.EQ.K20) THEN
            K1=K10
            K2=K20
            ENDIF
         IF (((I1.NE.I10.OR.I2.NE.I20).AND.(I10.NE.I20)) .OR. &
             ((J1.NE.J10.OR.J2.NE.J20).AND.(J10.NE.J20)) .OR. &
             ((K1.NE.K10.OR.K2.NE.K20).AND.(K10.NE.K20))) THEN
            WRITE(6,*) ' ERROR: Slice files are incompatible'
            STOP
         ENDIF
      ENDIF
 
      NCOUNT = 0
      READ_LOOP: DO
      READ(12,END=99) TIME
      READ(12,END=99) (((F(I,J,K),I=I1,I2),J=J1,J2),K=K1,K2)
      IF (TIME.LT.TBEG) CYCLE READ_LOOP
      IF (TIME.GT.TEND) EXIT READ_LOOP
      NCOUNT = NCOUNT + 1
 
      Q(I1:I2,J1:J2,K1:K2,MV) = Q(I1:I2,J1:J2,K1:K2,MV)+F(I1:I2,J1:J2,K1:K2)
      ENDDO READ_LOOP

99    CLOSE(12)
 
      IF (NCOUNT.EQ.0) NCOUNT=1
 
      DO K=K1,K2
      DO J=J1,J2
      DO I=I1,I2
      Q(I,J,K,MV) = Q(I,J,K,MV)/REAL(NCOUNT)
      ENDDO
      ENDDO
      ENDDO
 
      SELECT CASE(IOR_SLCF)
      CASE(1)
      DO K=K1+1,K2
      DO J=J1+1,J2
      SUM(MV) = SUM(MV) + 0.25*(Q(I1,J,K,MV)+Q(I1,J-1,K,MV)+Q(I1,J,K-1,MV)+Q(I1,J-1,K-1,MV))*(M%Y(J)-M%Y(J-1))*(M%Z(K)-M%Z(K-1))
      ENDDO
      ENDDO
      CASE(2)
      DO K=K1+1,K2
      DO I=I1+1,I2
      SUM(MV) = SUM(MV) + 0.25*(Q(I,J1,K,MV)+Q(I-1,J1,K,MV)+Q(I,J1,K-1,MV)+Q(I-1,J1,K-1,MV))*(M%X(I)-M%X(I-1))*(M%Z(K)-M%Z(K-1))
      ENDDO
      ENDDO
      CASE(3)
      DO J=J1+1,J2
      DO I=I1+1,I2
      SUM(MV) = SUM(MV) + 0.25*(Q(I,J,K1,MV)+Q(I-1,J,K1,MV)+Q(I,J-1,K1,MV)+Q(I-1,J-1,K1,MV))*(M%X(I)-M%X(I-1))*(M%Y(J)-M%Y(J-1))
      ENDDO
      ENDDO
      END SELECT
 
      IF(BATCHMODE.EQ.0)write(6,'(A,A,A,ES12.4)') ' Integral of ',TRIM(SLCF_TEXT(IS(MV))),' = ',SUM(MV)
 
      ENDDO VARLOOP

      ! WRITE OUT SLICE DATA IF AUTO_SLICE_FLAG IS SET
      !(WRITE IT OUT HERE RATHER THAN later SO WE DON'T HAVE TO SAVE EXTRA DATA)
         
      IF (AUTO_SLICE_FLAG.EQ.1) THEN
        IF(IAUTO.EQ.1)THEN
          IF(BATCHMODE.EQ.0)write(6,*) 'Enter output file name:'
          READ(LU_IN,'(A)') OUTFILE
          OPEN(44,FILE=OUTFILE,FORM='FORMATTED',STATUS='UNKNOWN')
 
          write(6,*) ' Writing to file...      ',TRIM(OUTFILE)
          ZMIN=-1000000000000.0
        ENDIF
 
        I3 = I2 - I1 + 1
        J3 = J2 - J1 + 1
        K3 = K2 - K1 + 1
 
! One-dimensional section file
!    NOTE: IN AUTO_SLICE MODE ONLY VERTICAL 1D SLICES ARE SUPPORTED
 
        IF (I1.EQ.I2 .AND. J1.EQ.J2 .AND. K1.NE.K2) then
          IF(IAUTO.EQ.1)THEN
            WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV,"(A,','),A)"
            WRITE(44,FRMT) 'Z',(TRIM(SLCF_LABEL(IS(L))),L=1,NV)
            WRITE(44,FRMT) 'm',(TRIM(SLCF_UNIT(IS(L))),L=1,NV)
            WRITE(FRMT,'(A,I2.2,A)') "(",NV,"(E12.5,','),E12.5)"
          ENDIF
          LOOP1A: DO K=K1,K2,NSAM
            IF (M%Z(K).GT.ZF .OR. M%Z(K).LT.ZS) CYCLE LOOP1A
            IF (ZOFFSET_FLAG.EQ.1.AND.M%Z(K)-ZOFFSET.LT.0.0)CYCLE LOOP1A
            IF(M%Z(K).LT.AZ1)CYCLE LOOP1A
            IF(M%Z(K).GT.AZ2)EXIT LOOP1A
            IF(M%Z(K)-ZOFFSET.GT.ZMIN)THEN
              WRITE(44,FRMT) M%Z(K)-ZOFFSET,(Q(I2,J2,K,L),L=1,NV)
            ENDIF
            ZMIN=M%Z(k)-ZOFFSET
          enddo LOOP1A
        endif
        ENDIF
      ENDDO AUTO_LIST
  

      CASE(3) FILETYPE

         BNDF_MESH = 1
         REWIND(11)
         NFILES_EXIST=0

         SEARCH_BNDF: DO I=1,FILE_DIM
            CALL SEARCH2('BNDF',4,'BNDC',4,11,IERR,CHOICE)
            IF (IERR.EQ.1) EXIT SEARCH_BNDF
            BACKSPACE(11)
            READ(11,*) JUNK,BNDF_MESH(I)
            READ(11,'(A)') BNDF_FILE(I)
            READ(11,'(A)') BNDF_TEXT(I)
            READ(11,*)
            READ(11,'(A)') BNDF_UNIT(I)
            OPEN(12,FILE=BNDF_FILE(I),FORM='UNFORMATTED',STATUS='OLD', IOSTAT=RCODE)
            CLOSE(12)
            IF (RCODE.NE.0) CYCLE
            NFILES_EXIST=NFILES_EXIST+1
            IF (CHOICE=='BNDC') BNDF_TYPE(I) = 'CENTERED'
            IF (CHOICE=='BNDF') BNDF_TYPE(I) = 'STAGGERED'
            IF (BATCHMODE.EQ.0) write(6,'(I3,3X,A,I2,A,A)') I,'MESH ',BNDF_MESH(I), ', ',TRIM(BNDF_TEXT(I))
         ENDDO SEARCH_BNDF
    
         IF (NFILES_EXIST.EQ.0)THEN
            IF(BATCHMODE.EQ.0) write(6,*)"There are no boundary files to convert"
            STOP
         ENDIF
         
         IF (BATCHMODE.EQ.0) write(6,*) ' Enter starting and ending time for averaging (s)'
         READ(LU_IN,*) TBEG,TEND
         IF (BATCHMODE.EQ.0) write(6,*) ' Enter orientation: (plus or minus 1, 2 or 3)'
         READ(LU_IN,*) IOR_INPUT
    
         IF(BATCHMODE.EQ.0)write(6,*) ' Enter number of variables'
         READ(LU_IN,*) NV
    
         BVARLOOP: DO MV=1,NV
   
            IF (BATCHMODE.EQ.0) write(6,'(A,I2)') ' Enter boundary file index for variable',MV
            READ(LU_IN,*) I
            IB(MV) = I
            QFILE = BNDF_FILE(I)
            OPEN(12,FILE=QFILE,FORM='UNFORMATTED',STATUS='OLD', IOSTAT=RCODE)
            IF (RCODE.NE.0) THEN
              CLOSE(12)
              CYCLE
            ENDIF
    
            IF (MV.EQ.1) THEN
               NM = BNDF_MESH(I)
               M=>MESH(NM)
               ALLOCATE(Q(0:M%IBAR,0:M%JBAR,0:M%KBAR,NV))
               Q = 0.
               ALLOCATE(ALREADY_USED(0:M%IBAR,0:M%JBAR,0:M%KBAR))
               BNDF_TYPE_CHOSEN = BNDF_TYPE(I)
            ELSE
               IF (BNDF_MESH(I).NE.NM) THEN
                  WRITE(6,*) ' ERROR: All boundary files must have the same mesh'
                  STOP
               ENDIF
               IF (BNDF_TYPE(I).NE.BNDF_TYPE_CHOSEN) THEN
                  IF (BNDF_TYPE_CHOSEN=='CENTERED') WRITE(6,*) ' ERROR: All boundary files must be CELL_CENTERED'
                  IF (BNDF_TYPE_CHOSEN/='CENTERED') WRITE(6,*) ' ERROR: No boundary files can be CELL_CENTERED'
                  STOP
               ENDIF
            ENDIF
    
            READ(12)
            READ(12)
            READ(12)
            READ(12) NPATCH
    
            IF (MV.EQ.1) THEN
               ALLOCATE(IOR(1:NPATCH))
               ALLOCATE(I1B(1:NPATCH))
               ALLOCATE(I2B(1:NPATCH))
               ALLOCATE(J1B(1:NPATCH))
               ALLOCATE(J2B(1:NPATCH))
               ALLOCATE(K1B(1:NPATCH))
               ALLOCATE(K2B(1:NPATCH))
            ENDIF
    
            DO I=1,NPATCH
               READ(12) I1B(I),I2B(I),J1B(I),J2B(I),K1B(I),K2B(I),IOR(I)
            ENDDO
    
            IF (MV.EQ.1) THEN
               IJBAR = MAX(M%IBAR,M%JBAR)
               JKBAR = MAX(M%JBAR,M%KBAR)
               ALLOCATE(F(0:IJBAR,0:JKBAR,NPATCH))
            ENDIF
            F = 0.
    
            NCOUNT = 0
            READ_BLOOP: DO
               READ(12,END=199) TIME
               DO II=1,NPATCH
                  SELECT CASE(ABS(IOR(II)))
                     CASE(1)
                        IF (BNDF_TYPE_CHOSEN=='STAGGERED')READ(12,END=199) ((F(J,K,II),J=J1B(II),J2B(II)),K=K1B(II),K2B(II))
                        IF (BNDF_TYPE_CHOSEN=='CENTERED') READ(12,END=199) ((F(J,K,II),J=J1B(II)+1,J2B(II)+1),K=K1B(II)+1,K2B(II)+1)
                     CASE(2)
                        IF (BNDF_TYPE_CHOSEN=='STAGGERED')READ(12,END=199) ((F(I,K,II),I=I1B(II),I2B(II)),K=K1B(II),K2B(II))
                        IF (BNDF_TYPE_CHOSEN=='CENTERED') READ(12,END=199) ((F(I,K,II),I=I1B(II)+1,I2B(II)+1),K=K1B(II)+1,K2B(II)+1)
                     CASE(3)
                        IF (BNDF_TYPE_CHOSEN=='STAGGERED')READ(12,END=199) ((F(I,J,II),I=I1B(II),I2B(II)),J=J1B(II),J2B(II))
                        IF (BNDF_TYPE_CHOSEN=='CENTERED') READ(12,END=199) ((F(I,J,II),I=I1B(II)+1,I2B(II)+1),J=J1B(II)+1,J2B(II)+1)
                  END SELECT
               ENDDO
               IF (TIME.LT.TBEG) CYCLE READ_BLOOP
               IF (TIME.GT.TEND) EXIT READ_BLOOP
       
               NCOUNT = NCOUNT + 1
               ALREADY_USED = .FALSE.
    
               REC_PATCH: DO II=1,NPATCH
                  IF (IOR(II).NE.IOR_INPUT) CYCLE REC_PATCH
                  SELECT CASE(ABS(IOR(II)))
                     CASE(1)
                        DO K=K1B(II),K2B(II)
                           DO J=J1B(II),J2B(II)
                              IF (.NOT.ALREADY_USED(I1B(II),J,K)) THEN
                                 Q(I1B(II),J,K,MV) = Q(I1B(II),J,K,MV) + F(J,K,II)
                                 ALREADY_USED(I1B(II),J,K) = .TRUE.
                              ENDIF
                           ENDDO
                        ENDDO
                     CASE(2)
                        DO K=K1B(II),K2B(II)
                           DO I=I1B(II),I2B(II)
                              IF (.NOT.ALREADY_USED(I,J1B(II),K)) THEN
                                 Q(I,J1B(II),K,MV) = Q(I,J1B(II),K,MV) + F(I,K,II)
                                 ALREADY_USED(I,J1B(II),K) = .TRUE.
                              ENDIF
                           ENDDO
                        ENDDO
                     CASE(3)
                        DO J=J1B(II),J2B(II)
                           DO I=I1B(II),I2B(II)
                              IF (.NOT.ALREADY_USED(I,J,K1B(II))) THEN
                                 Q(I,J,K1B(II),MV) = Q(I,J,K1B(II),MV) + F(I,J,II)
                                 ALREADY_USED(I,J,K1B(II)) = .TRUE.
                              ENDIF
                           ENDDO
                        ENDDO
                  END SELECT
               ENDDO REC_PATCH
       
            ENDDO READ_BLOOP
      
               199 CLOSE(12)
    
               DO K=0,M%KBAR
                  DO J=0,M%JBAR
                     DO I=0,M%IBAR
                        Q(I,J,K,MV) = Q(I,J,K,MV)/REAL(NCOUNT)
                     ENDDO
                  ENDDO
               ENDDO
    
         ENDDO BVARLOOP
    
   END SELECT FILETYPE
   
   ! AUTO_SLICE DATA WAS ALREADY GENERATED

   IF(AUTO_SLICE_FLAG.EQ.1) STOP

   ! Write out the data to an ASCII file

   IF(BATCHMODE.EQ.0)write(6,*) 'Enter output file name:'
   READ(LU_IN,'(A)') OUTFILE
   OPEN(44,FILE=OUTFILE,FORM='FORMATTED',STATUS='UNKNOWN')
 
   write(6,*) ' Writing to file...      ',TRIM(OUTFILE)
 
   SELECT CASE(IFILE)
 
      CASE(1:2)
 
      I3 = I2 - I1 + 1
      J3 = J2 - J1 + 1
      K3 = K2 - K1 + 1

      ! One-dimensional section file

      IF (I1.EQ.I2 .AND. J1.EQ.J2 .AND. K1.NE.K2) then
        WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV,"(A,','),A)"
        WRITE(44,FRMT) 'Z',(TRIM(SLCF_LABEL(IS(L))),L=1,NV)
        WRITE(44,FRMT) 'm',(TRIM(SLCF_UNIT(IS(L))),L=1,NV)
        WRITE(FRMT,'(A,I2.2,A)') "(",NV,"(E12.5,','),E12.5)"
        LOOP1: DO K=K1,K2,NSAM
        IF (M%Z(K).GT.ZF .OR. M%Z(K).LT.ZS) CYCLE LOOP1
        IF (ZOFFSET_FLAG.EQ.1.AND.M%Z(K)-ZOFFSET.LT.0.0) CYCLE LOOP1
        write(44,FRMT) M%Z(K)-ZOFFSET,(Q(I2,J2,K,L),L=1,NV)
        enddo LOOP1
        endif

        if(i1.eq.i2.and.j1.ne.j2.and.k1.eq.k2) then
        WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV,"(A,','),A)"
        WRITE(44,FRMT) 'Y',(TRIM(SLCF_LABEL(IS(L))),L=1,NV)
        WRITE(44,FRMT) 'm',(TRIM(SLCF_UNIT(IS(L))),L=1,NV)
        WRITE(FRMT,'(A,I2.2,A)') "(",NV,"(E12.5,','),E12.5)"
        LOOP2: DO J=J1,J2,NSAM
        IF (M%Y(J).GT.YF .OR. M%Y(J).LT.YS) CYCLE LOOP2
        write(44,FRMT) M%y(j),(q(i2,j,k2,l),l=1,nv)
        enddo LOOP2
        endif

        if(i1.ne.i2.and.j1.eq.j2.and.k1.eq.k2) then
        WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV,"(A,','),A)"
        WRITE(44,FRMT) 'X',(TRIM(SLCF_LABEL(IS(L))),L=1,NV)
        WRITE(44,FRMT) 'm',(TRIM(SLCF_UNIT(IS(L))),L=1,NV)
        WRITE(FRMT,'(A,I2.2,A)') "(",NV,"(E12.5,','),E12.5)"
        LOOP3: DO I=I1,I2,NSAM
        IF (M%X(I).GT.XF .OR. M%X(I).LT.XS) CYCLE LOOP3
        write(44,FRMT) M%x(i),(q(i,j2,k2,l),l=1,nv)
        enddo LOOP3
        endif

        ! Two-dimensional section file

        if(i1.eq.i2.and.j1.ne.j2.and.k1.ne.k2) then
        WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV+1,"(A,','),A)"
        WRITE(44,FRMT) 'Y','Z',(TRIM(SLCF_LABEL(IS(L))),L=1,NV)
        WRITE(44,FRMT) 'm','m',(TRIM(SLCF_UNIT(IS(L))),L=1,NV)
        WRITE(FRMT,'(A,I2.2,A)') "(",NV+1,"(E12.5,','),E12.5)"
        DO K=K1,K2,NSAM
        LOOP4: DO J=J1,J2,NSAM
        IF (M%Y(J).GT.YF .OR. M%Y(J).LT.YS) CYCLE LOOP4
        IF (M%Z(K).GT.ZF .OR. M%Z(K).LT.ZS) CYCLE LOOP4
        write(44,FRMT) M%y(j),M%z(k)-zoffset,(q(i2,j,k,l),l=1,nv)
        enddo LOOP4
        enddo
        endif

        if (j1.eq.j2.and.i1.ne.i2.and.k1.ne.k2) then
        WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV+1,"(A,','),A)"
        WRITE(44,FRMT) 'X','Z',(TRIM(SLCF_LABEL(IS(L))),L=1,NV)
        WRITE(44,FRMT) 'm','m',(TRIM(SLCF_UNIT(IS(L))),L=1,NV)
        WRITE(FRMT,'(A,I2.2,A)') "(",NV+1,"(E12.5,','),E12.5)"
        DO K=K1,K2,NSAM
        LOOP5: DO I=I1,I2,NSAM
        IF (M%X(I).GT.XF .OR. M%X(I).LT.XS) CYCLE LOOP5
        IF (M%Z(K).GT.ZF .OR. M%Z(K).LT.ZS) CYCLE LOOP5
        write(44,FRMT) M%x(i),M%z(k)-zoffset,(q(i,j2,k,l),l=1,nv)
        enddo LOOP5
        enddo
        endif

        if(k1.eq.k2.and.i1.ne.i2.and.j1.ne.j2) then
        WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV+1,"(A,','),A)"
        WRITE(44,FRMT) 'X','Y',(TRIM(SLCF_LABEL(IS(L))),L=1,NV)
        WRITE(44,FRMT) 'm','m',(TRIM(SLCF_UNIT(IS(L))),L=1,NV)
        WRITE(FRMT,'(A,I2.2,A)') "(",NV+1,"(E12.5,','),E12.5)"
        DO J=J1,J2,NSAM
        LOOP6: DO I=I1,I2,NSAM
        IF (M%X(I).GT.XF .OR. M%X(I).LT.XS) CYCLE LOOP6
        IF (M%Y(J).GT.YF .OR. M%Y(J).LT.YS) CYCLE LOOP6
        write(44,FRMT) M%x(i),M%y(j),(q(i,j,k2,l),l=1,nv)
        enddo LOOP6
        enddo
        endif

        ! Three-dimensional section file

        if (i1.ne.i2.and.j1.ne.j2.and.k1.ne.k2) then
        WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV+2,"(A,','),A)"
        WRITE(44,FRMT) 'X','Y','Z',(TRIM(SLCF_LABEL(IS(L))),L=1,NV)
        WRITE(44,FRMT) 'm','m','m',(TRIM(SLCF_UNIT(IS(L))),L=1,NV)
        WRITE(FRMT,'(A,I2.2,A)') "(",NV+2,"(E12.5,','),E12.5)"
        DO K=K1,K2,NSAM
        DO J=J1,J2,NSAM
        LOOP7: DO I=I1,I2,NSAM
        IF (M%X(I).GT.XF .OR. M%X(I).LT.XS) CYCLE LOOP7
        IF (M%Y(J).GT.YF .OR. M%Y(J).LT.YS) CYCLE LOOP7
        IF (M%Z(K).GT.ZF .OR. M%Z(K).LT.ZS) CYCLE LOOP7
        write(44,FRMT) M%x(i),M%y(j),M%z(k),(q(i,j,k,l),l=1,nv)
        enddo LOOP7
        enddo
        enddo
        endif
 
      CASE(3)  ! Write out boundary file output
 
         PATCHES: DO II=1,NPATCH
            IF (IOR(II).NE.IOR_INPUT) CYCLE PATCHES
            IF (M%X(I1B(II)).GT.XF .OR. M%X(I2B(II)).LT.XS) CYCLE PATCHES
            IF (M%Y(J1B(II)).GT.YF .OR. M%Y(J2B(II)).LT.YS) CYCLE PATCHES
            IF (M%Z(K1B(II)).GT.ZF .OR. M%Z(K2B(II)).LT.ZS) CYCLE PATCHES
            WRITE(44,'(A,I4,5(F7.2,A),F8.2)') 'Patch',II, &
               M%X(I1B(II)),'<x<',M%X(I2B(II)),',  ', &
               M%Y(J1B(II)),'<y<',M%Y(J2B(II)),',  ', &
               M%Z(K1B(II)),'<z<',M%Z(K2B(II)) 
            WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV+2,"(A,','),A)"
            WRITE(44,FRMT) 'X','Y','Z',(TRIM(BNDF_TEXT(IB(L))),L=1,NV)
            WRITE(44,FRMT) 'm','m','m',(TRIM(BNDF_UNIT(IB(L))),L=1,NV)
            WRITE(FRMT,'(A,I2.2,A)') "(",NV+2,"(E12.5,','),E12.5)"
            IF (BNDF_TYPE_CHOSEN=='STAGGERED') THEN
               DO K=K1B(II),K2B(II),NSAM
                  DO J=J1B(II),J2B(II),NSAM
                     DO I=I1B(II),I2B(II),NSAM
                        IF (M%X(I).GT.XF .OR. M%X(I).LT.XS) CYCLE 
                        IF (M%Y(J).GT.YF .OR. M%Y(J).LT.YS) CYCLE
                        IF (M%Z(K).GT.ZF .OR. M%Z(K).LT.ZS) CYCLE
                        WRITE(44,FRMT) M%X(I),M%Y(J),M%Z(K),(Q(I,J,K,L),L=1,NV)
                     ENDDO 
                  ENDDO
               ENDDO
            ELSE
               IF (ABS(IOR_INPUT)==1) I1B(II) = I1B(II)-1
               IF (ABS(IOR_INPUT)==2) J1B(II) = J1B(II)-1
               IF (ABS(IOR_INPUT)==3) K1B(II) = K1B(II)-1
               DO K=K1B(II)+1,K2B(II),NSAM
                  DO J=J1B(II)+1,J2B(II),NSAM
                     DO I=I1B(II)+1,I2B(II),NSAM
                        IF (M%X(I).GT.XF .OR. M%X(I).LT.XS) CYCLE
                        IF (M%Y(J).GT.YF .OR. M%Y(J).LT.YS) CYCLE
                        IF (M%Z(K).GT.ZF .OR. M%Z(K).LT.ZS) CYCLE
                        IF (ABS(IOR_INPUT)==1) WRITE(44,FRMT) M%X(I),0.5*(M%Y(J-1)+M%Y(J)),0.5*(M%Z(K-1)+M%Z(K)),(Q(I,J,K,L),L=1,NV)
                        IF (ABS(IOR_INPUT)==2) WRITE(44,FRMT) 0.5*(M%X(I-1)+M%X(I)),M%Y(J),0.5*(M%Z(K-1)+M%Z(K)),(Q(I,J,K,L),L=1,NV)
                        IF (ABS(IOR_INPUT)==3) WRITE(44,FRMT) 0.5*(M%X(I-1)+M%X(I)),0.5*(M%Y(J-1)+M%Y(J)),M%Z(K),(Q(I,J,K,L),L=1,NV)
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF
         ENDDO PATCHES
 
   END SELECT
 
   CLOSE(44)
 
   IF (IFILE.NE.1) EXIT EXTRA_PLOT3D_FILES

   NEW_PLOT3D = .FALSE.

ENDDO EXTRA_PLOT3D_FILES

STOP
END PROGRAM FDS2ASCII

! *********************** SEARCH *******************************

SUBROUTINE SEARCH(STRING,LENGTH,LU,IERR)

IMPLICIT NONE
CHARACTER(*), INTENT(IN) :: STRING
INTEGER, INTENT(OUT) :: IERR
INTEGER, INTENT(IN) :: LU, LENGTH
CHARACTER(20) :: JUNK

SEARCH_LOOP: DO 
   READ(LU,'(A)',END=10) JUNK
   IF (JUNK(1:LENGTH).EQ.STRING(1:LENGTH)) EXIT SEARCH_LOOP
ENDDO SEARCH_LOOP

IERR = 0
RETURN

10 IERR = 1
RETURN

END SUBROUTINE SEARCH

SUBROUTINE SEARCH2(STRING,LENGTH,STRING2,LENGTH2,LU,IERR,CHOICE)

IMPLICIT NONE
CHARACTER(*), INTENT(IN) :: STRING,STRING2
CHARACTER(*), INTENT(OUT) :: CHOICE
INTEGER, INTENT(OUT) :: IERR
INTEGER, INTENT(IN) :: LU, LENGTH, LENGTH2
CHARACTER(20) :: JUNK

SEARCH_LOOP: DO 
   READ(LU,'(A)',END=10) JUNK
   IF (JUNK(1:LENGTH).EQ.STRING(1:LENGTH).OR.JUNK(1:LENGTH2).EQ.STRING2(1:LENGTH2)) THEN
      IF (JUNK(1:LENGTH) .EQ.STRING(1:LENGTH))   CHOICE = JUNK(1:LENGTH)
      IF (JUNK(1:LENGTH2).EQ.STRING2(1:LENGTH2)) CHOICE = JUNK(1:LENGTH2)
      EXIT SEARCH_LOOP
   ENDIF
ENDDO SEARCH_LOOP

IERR = 0
RETURN

10 IERR = 1
RETURN

END SUBROUTINE SEARCH2

! *********************** USAGE2 *******************************

SUBROUTINE USAGE2(f2aversion)
IMPLICIT NONE

CHARACTER(255), intent(in) :: f2aversion
INTEGER :: lastchar
        
        
WRITE(6,*)"fds2ascii ",trim(f2aversion)," ",TRIM(GITHASH_PP)
WRITE(6,*)""
WRITE(6,*)"  Convert boundary, slice or plot3d data generated"
WRITE(6,*)"  by FDS to an ascii spreadsheet file."
WRITE(6,*)""
WRITE(6,*)"Usage:"
WRITE(6,*)"  fds2ascii [-h] [-v] [input]"
WRITE(6,*)""
WRITE(6,*)"   -h    - print out this message"
WRITE(6,*)"   -v    - print out version info"
WRITE(6,*)"   input - read input from a file named input or from"
WRITE(6,*)"           the console if no file is specified"
END SUBROUTINE USAGE2

! *********************** VERSION2 *******************************

SUBROUTINE VERSION2(f2aversion)
IMPLICIT NONE

CHARACTER(255), intent(in) :: f2aversion

CHARACTER(60) :: DATE

WRITE(6,'(/A/)')      ' fds2ascii'
WRITE(6,'(A,A)')      ' Version          : ',TRIM(f2aversion)
WRITE(6,'(A,A)')      ' Revision         : ',TRIM(GITHASH_PP)
WRITE(6,'(A,A)')      ' Revision Date    : ',TRIM(GITDATE_PP)
WRITE(6,'(A,A/)')     ' Compilation Date : ',TRIM(BUILDDATE_PP)


END SUBROUTINE VERSION2

! *********************** PARSE *******************************

SUBROUTINE PARSE(BUFFER,SB_TOKS,SE_TOKS,N_TOKS)

! parse buffer into a series of tokens each separated by blanks

IMPLICIT NONE
CHARACTER(*), INTENT(INOUT) :: BUFFER
INTEGER, DIMENSION(*), INTENT(OUT) :: SB_TOKS, SE_TOKS
INTEGER, INTENT(OUT) :: N_TOKS
INTEGER :: I, INTOK, INQUOTE, LENBUF
CHARACTER(LEN=1) :: C

N_TOKS=0
LENBUF = LEN(TRIM(BUFFER))
IF(LENBUF.EQ.0)RETURN ! buffer is blank so there are no tokens
INTOK=0
INQUOTE=0
DO I = 1, LENBUF
  IF(INTOK.EQ.0)THEN  
     IF(BUFFER(I:I).NE.' ')THEN  ! beginning of a new token since previous char
       INTOK=1                   ! was not in a token and this one is
       N_TOKS=N_TOKS + 1
       SB_TOKS(N_TOKS)=I
       IF(BUFFER(I:I).EQ."'")INQUOTE=1
     ENDIF
  ELSE
     IF(INQUOTE.EQ.1)THEN
        IF(BUFFER(I:I).EQ."'")THEN
           SE_TOKS(N_TOKS)=I
           INTOK=0
           INQUOTE=0
        ENDIF
     ENDIF
     IF(BUFFER(I:I).EQ.' ')THEN
       SE_TOKS(N_TOKS)=I-1       ! previous char was in a token, this one is not
       INTOK=0                   ! so previous char is end of token
     ENDIF
  ENDIF
END DO
IF(BUFFER(LENBUF:LENBUF).NE.' ')SE_TOKS(N_TOKS)=LENBUF ! last char in buffer is not blank
                                                       ! so it is end of last token

! strip out single or double quotes if present
DO I = 1, N_TOKS
   C = BUFFER(SB_TOKS(I):SB_TOKS(I))
   IF(C.EQ."'")SB_TOKS(I)=SB_TOKS(I)+1
   C = BUFFER(SE_TOKS(I):SE_TOKS(I))
   IF(C.EQ."'")SE_TOKS(I)=SE_TOKS(I)-1
   IF(SE_TOKS(I).LT.SB_TOKS(I))THEN
     SE_TOKS(I)=SB_TOKS(I)
     BUFFER(SB_TOKS(I):SE_TOKS(I))=' '
   ENDIF
END DO
END SUBROUTINE PARSE

! *********************** TOUPPER *******************************

SUBROUTINE TOUPPER(BUFFERIN, BUFFEROUT)
CHARACTER(LEN=*), INTENT(IN) :: BUFFERIN
CHARACTER(LEN=*), INTENT(OUT) :: BUFFEROUT
CHARACTER(LEN=1) :: C

INTEGER :: LENBUF, I

LENBUF=MIN(LEN(TRIM(BUFFERIN)),LEN(BUFFEROUT))
DO I = 1, LENBUF
   C = BUFFERIN(I:I)
   IF(C.GE.'a'.AND.C.LE.'z')C=CHAR(ICHAR(C)+ICHAR('A')-ICHAR('a'))
    BUFFEROUT(I:I)=C
END DO

END SUBROUTINE TOUPPER
