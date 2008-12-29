      PROGRAM FDS2ASCII
      IMPLICIT NONE
C $Id$  
C $Revision$
C $Date$
 
C
C Program to convert various FDS output files to ASCII
C
      INTEGER, PARAMETER :: FB = SELECTED_REAL_KIND(6)
      INTEGER :: IERR, VERSION, NMESHES, NM, NOC, I, J, K, L
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
      INTEGER, DIMENSION(500) :: AUTO_SLICE_LISTS
      REAL(FB), DIMENSION(500) :: AUTO_SLICE_Z
      REAL(FB) :: AX1, AY1, AZ1, AZ2
      REAL(FB) :: EPS, ZMIN
C
      TYPE MESH_TYPE
      REAL(FB), POINTER, DIMENSION(:) :: X,Y,Z
      REAL(FB) :: D1,D2,D3,D4
      INTEGER :: IBAR,JBAR,KBAR,IERR,NXP,NYP,NZP
      END TYPE MESH_TYPE
C
      TYPE (MESH_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: MESH
      TYPE(MESH_TYPE), POINTER :: M
C
      REAL(FB), DIMENSION(60) :: SUM
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IOR,I1B,I2B,J1B,J2B,K1B,K2B
      INTEGER, ALLOCATABLE, DIMENSION(:) :: AUTO_SLICES
      REAL(FB), ALLOCATABLE, DIMENSION(:,:,:,:) :: Q
      REAL(FB), ALLOCATABLE, DIMENSION(:,:,:) :: F
      LOGICAL :: NEW_PLOT3D=.TRUE.
      INTEGER IOR_SLCF
      CHARACTER(2) ANS
      CHARACTER(30) :: UNITJUNK
      CHARACTER(256) GRIDFILE,QNAME,CHID,QFILE,JUNK,FRMT,OUTFILE
      CHARACTER(256), DIMENSION(500) :: PL3D_FILE,SLCF_FILE,BNDF_FILE,
     .                                 SLCF_TEXT,BNDF_TEXT,
     .                                 SLCF_UNIT,BNDF_UNIT
      CHARACTER(256), DIMENSION(500) :: SLICE_LABELS
      INTEGER :: NSLICE_LABELS, SLICE_EXIST
      REAL(FB), DIMENSION(500) :: X1, X2, Y1, Y2, Z1, Z2
      INTEGER,  DIMENSION(60) :: IB,IS
      INTEGER,  DIMENSION(500) :: PL3D_MESH,SLCF_MESH,BNDF_MESH
      REAL(FB), DIMENSION(500) :: PL3D_TIME
      LOGICAL, DIMENSION(500) :: FILE_EXISTS
      INTEGER :: RCODE
      INTEGER :: NFILES_EXIST
      REAL(FB) :: ZOFFSET
      INTEGER :: ZOFFSET_FLAG
      LOGICAL :: EXISTS

      ZOFFSET=0.0
      ZOFFSET_FLAG=0
      EPS=0.001
C
      WRITE(6,*) ' Enter Job ID string (CHID):'
      READ(5,'(a)') CHID
      GRIDFILE = TRIM(CHID)//'.smv'
C
C Open grid file
C
      INQUIRE(FILE=GRIDFILE,EXIST=EXISTS)
      IF(.NOT.EXISTS)THEN
        WRITE(6,*)"*** Fatal error: The file: ",TRIM(GRIDFILE),
     .        " does not exist"
        STOP
      ENDIF

      OPEN(11,FILE=GRIDFILE,STATUS='OLD',FORM='FORMATTED')
C
      CALL SEARCH('VERSION',7,11,IERR)
      IF (IERR.EQ.1) THEN
         WRITE(6,*) ' WARNING: Assuming FDS version 2 or less'
         VERSION = 2.
      ELSE
         READ(11,*) VERSION
      ENDIF
C
      REWIND(11)
C
      CALL SEARCH('NMESHES',7,11,IERR)
      IF (IERR.EQ.1) THEN
         WRITE(6,*) ' WARNING: Assuming 1 mesh'
         NMESHES = 1
      ELSE
         READ(11,*) NMESHES
      ENDIF
      ALLOCATE(MESH(NMESHES))
C
      REWIND(11)
C
      READ_SMV: DO NM=1,NMESHES
      M=>MESH(NM)
      CALL SEARCH('GRID',4,11,IERR)
      READ(11,*) M%IBAR,M%JBAR,M%KBAR
      ALLOCATE(M%X(0:M%IBAR))
      ALLOCATE(M%Y(0:M%JBAR))
      ALLOCATE(M%Z(0:M%KBAR))
C
      CALL SEARCH('TRNX',4,11,IERR)
      READ(11,*) NOC
      DO I=1,NOC
      READ(11,*)
      ENDDO
      DO I=0,M%IBAR
      READ(11,*) IDUM,M%X(I)
      ENDDO
C
      CALL SEARCH('TRNY',4,11,IERR)
      READ(11,*) NOC
      DO I=1,NOC
      READ(11,*)
      ENDDO
      DO J=0,M%JBAR
      READ(11,*) IDUM,M%Y(J)
      ENDDO
C
      CALL SEARCH('TRNZ',4,11,IERR)
      READ(11,*) NOC
      DO I=1,NOC
      READ(11,*)
      ENDDO
      DO K=0,M%KBAR
      READ(11,*) IDUM,M%Z(K)
      ENDDO
C
      ENDDO READ_SMV
C
      WRITE(6,*) ' What type of file to parse?'
      WRITE(6,*) ' PL3D file? Enter 1'
      WRITE(6,*) ' SLCF file? Enter 2'
      WRITE(6,*) ' BNDF file? Enter 3'
      READ(5,*) IFILE
C
      WRITE(6,*) ' Enter Sampling Factor for Data?'
      WRITE(6,*) ' (1 for all data, 2 for every other point, etc.)'
      READ(5,*) NSAM
C
C*** ANS(1:1) may be 'y', 'n' or 'z' (or upper case equivalents)
C    ANS(2:2) may be 'a' or ' '
C
C     y - domain size is limited
C     n - domain size is not limited
C     z - domain size is not limted z levels are offset by zoffset
C     a - slice files are selected based on slice file type and location

      WRITE(6,*) ' Limit the domain size? (y or n)'
      READ(5,'(A)') ANS
      IF (ANS(1:1).EQ.'y' .OR. ANS(1:1).EQ.'Y') THEN
         WRITE(6,*) ' Enter min/max x, y and z'
         READ(5,*) XS,XF,YS,YF,ZS,ZF
      ELSE
         XS = -100000.
         XF =  100000.
         YS = -100000.
         YF =  100000.
         ZS = -100000.
         ZF =  100000.
      ENDIF
C      
C*** if ANS is z or Z then subtract zoffset from z data (for multi-level cases)
      IF (ANS(1:1).EQ.'z' .OR. ANS(1:1).EQ.'Z') THEN
         ZOFFSET_FLAG=1
        ELSE
         ZOFFSET_FLAG=0
      ENDIF
      IF (ANS(2:2).eq.'a'.or.ANS(2:2).eq.'A')THEN
         AUTO_SLICE_FLAG=1
        ELSE
         AUTO_SLICE_FLAG=0
      ENDIF
C
      EXTRA_PLOT3D_FILES: DO
C
      FILETYPE: SELECT CASE(IFILE)
C
      CASE(1) FILETYPE
C
C Read a PLOT3D file and print out the values in ASCII text
C
      PL3D_MESH = 1
      REWIND(11)
      NFILES_EXIST=0
      SEARCH_PL3D: DO I=1,500
      CALL SEARCH('PL3D',4,11,IERR)
      IF (IERR.EQ.1) EXIT SEARCH_PL3D
      BACKSPACE(11) 
      IF (VERSION.LE.2.) READ(11,*) JUNK,PL3D_TIME(I)
      IF (VERSION.GT.2.) READ(11,*) JUNK,PL3D_TIME(I),PL3D_MESH(I)
      READ(11,'(A)') PL3D_FILE(I)
      DO II=1,5
      IS(II) = II
      READ(11,'(A)') SLCF_TEXT(II)
      READ(11,*)
      READ(11,'(A)') SLCF_UNIT(II)
      ENDDO
      OPEN(12,FILE=TRIM(PL3D_FILE(I)),FORM='UNFORMATTED',STATUS='OLD',
     .                          IOSTAT=RCODE)
      CLOSE(12)
      FILE_EXISTS(I)=.TRUE.
      IF(RCODE.NE.0)THEN
        FILE_EXISTS(I)=.FALSE.
        CYCLE
      ENDIF
      NFILES_EXIST = NFILES_EXIST + 1
      WRITE(6,'(I3,3X,A,A,I2,A,F5.0)') I,TRIM(PL3D_FILE(I)),
     .    ', MESH ',PL3D_MESH(I),', Time: ',PL3D_TIME(I)
      ENDDO SEARCH_PL3D
C
      IF(NFILES_EXIST.EQ.0)THEN
        WRITE(6,*)"There are no PLOT3D files to convert"
        GOTO 500
      ENDIF
      
      WRITE(6,'(/A)') ' Choose PLOT3D .q file by its index: (0 to end)'
      READ(5,*) I
      IF(I.NE.0)THEN
        IF(.NOT.FILE_EXISTS(I))THEN
          WRITE(6,*)"The file, ",TRIM(PL3D_FILE(I)),", does not exist"
          CYCLE
        ENDIF
      ENDIF
      IF (I.EQ.0) GOTO 500
C
      QFILE = PL3D_FILE(I)
      NM = PL3D_MESH(I)
      M=>MESH(NM)
C
      I1 = 0
      I2 = M%IBAR
      J1 = 0
      J2 = M%JBAR
      K1 = 0
      K2 = M%KBAR
      NV = 5
      IF (NEW_PLOT3D) ALLOCATE(Q(0:M%IBAR,0:M%JBAR,0:M%KBAR,5))
C
      OPEN(12,FILE=QFILE,FORM='UNFORMATTED',STATUS='OLD')
      WRITE(6,*) ' Reading PLOT3D data file...'
      READ(12) NXP,NYP,NZP
      READ(12) D1,D2,D3,D4
      READ(12) ((((Q(I,J,K,N),I=0,M%IBAR),J=0,M%JBAR),K=0,M%KBAR),N=1,5)
      CLOSE(12)
C
      CASE(2) FILETYPE
C
C The following program averages up to six variables per section dumped 
C to a slice file.  The data is then dumped to ASCII file(s) in the form
C of the coordinate(s) and averaged variables.
C Example:  for two-dimensional data, the form is  x,y,t,u,v, where x,y
C are the coordinates, and t,u,v are averaged scalar variables.
C
      WRITE(6,*) ' Enter starting and ending time for averaging (s)'
      NSLICE_LABELS=0
      IF(ZOFFSET_FLAG.EQ.0)THEN
        READ(5,*) TBEG,TEND
      ELSE
        READ(5,*) TBEG,TEND,ZOFFSET
      ENDIF
      SLCF_MESH = 1
      REWIND(11)
      NFILES_EXIST=0
      SEARCH_SLCF: DO I=1,500
      CALL SEARCH('SLCF',4,11,IERR)
      IF (IERR.EQ.1) EXIT SEARCH_SLCF
      BACKSPACE(11)
      IF (VERSION.LE.2.) READ(11,*) JUNK
      IF (VERSION.GT.2.) READ(11,*) JUNK,SLCF_MESH(I)
      READ(11,'(A)') SLCF_FILE(I)
      READ(11,'(A)') SLCF_TEXT(I)

C*** create unique list of slice types      
      SLICE_EXIST=0
      DO II=1, NSLICE_LABELS
        IF(TRIM(SLCF_TEXT(I)).EQ.SLICE_LABELS(II))THEN
          SLICE_EXIST=1
          EXIT
        ENDIF
      ENDDO
      IF(SLICE_EXIST.EQ.0)THEN
        NSLICE_LABELS=NSLICE_LABELS+1
        SLICE_LABELS(NSLICE_LABELS)=TRIM(SLCF_TEXT(I))
      ENDIF
      READ(11,*) 
      READ(11,'(A)') SLCF_UNIT(I)
      OPEN(12,FILE=SLCF_FILE(I),FORM='UNFORMATTED',STATUS='OLD',
     .                          IOSTAT=RCODE)
      IF(RCODE.NE.0)THEN
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
        write(6,'(I3,1x,A,1x,A)')I,TRIM(SLCF_TEXT(I)),TRIM(SLCF_FILE(I))
        write(6,'(3x,A,6(1x,f8.2))')'slice bounds:',
     .  M%X(I1),M%X(I2),M%Y(J1),M%Y(J2),M%Z(K1),M%Z(K2)
      ENDIF
      ENDDO SEARCH_SLCF
C
      IF(NFILES_EXIST.EQ.0)THEN
        WRITE(6,*)"There are no slice files to convert"
        GOTO 500
      ENDIF

      IF(AUTO_SLICE_FLAG.EQ.0)THEN
        WRITE(6,*)'How many variables to read: (6 max)'
        READ(5,*) NV
        N_AUTO_SLICES=1
       ELSE
        N_AUTO_SLICES=0
        AUTO_SLICE_FLAG=1
        WRITE(6,*)' Enter slice file type index'
        DO II=1, NSLICE_LABELS
          WRITE(6,'(I3,1x,A)')II,TRIM(SLICE_LABELS(II))
        ENDDO
        READ(5,'(I3)')II
        AUTO_SLICE_LABEL=TRIM(SLICE_LABELS(II))
        boundloop: DO II=1,NFILES_EXIST
          IF(TRIM(AUTO_SLICE_LABEL).NE.TRIM(SLCF_TEXT(II)))THEN
            CYCLE boundloop
          ENDIF
          IF(X2(II)-X1(II).GT.EPS)CYCLE boundloop
          IF(Y2(II)-Y1(II).GT.EPS)CYCLE boundloop
          WRITE(6,'(4(a,f8.3))')"x=",X1(II)," y=",Y1(II),
     .                          " zmin=",Z1(II)," zmax=",Z2(II)
        ENDDO boundloop
        
        WRITE(6,*)"Enter 1D SLICE bounds (x y zmin zmax)"
        READ(5,*)AX1, AY1, AZ1, AZ2
        EPS=0.001
C*** create list of slices that match requested type and location
C         also record the elevation (z) of each slice so that they
C         may be later output in increasing order
        AUTO_LOOP1: DO I=1,NFILES_EXIST
          IF(TRIM(SLCF_TEXT(I)).EQ.TRIM(AUTO_SLICE_LABEL))THEN
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
      ENDIF
C
      SUM = 0.
C
      IF(AUTO_SLICE_FLAG.EQ.1)THEN
        ALLOCATE(Q(0:1,0:1,0:1,1))
        ALLOCATE(F(0:1,0:1,0:1))
        NV = 1
C*** sort slice list        
        DO I = 1, N_AUTO_SLICES
          IZMIN1= MINLOC(AUTO_SLICE_Z(I:N_AUTO_SLICES))
          IZMIN = IZMIN1(1) + I - 1
          IF(IZMIN.NE.I)THEN
            ZMIN=AUTO_SLICE_Z(IZMIN)
            AUTO_SLICE_Z(IZMIN)=AUTO_SLICE_Z(I)
            AUTO_SLICE_Z(I) = ZMIN
            IZTEMP = AUTO_SLICE_LISTS(IZMIN)
            AUTO_SLICE_LISTS(IZMIN)=AUTO_SLICE_LISTS(I)
            AUTO_SLICE_LISTS(I)=IZTEMP
          ENDIF
        END DO
       ELSE
        N_AUTO_SLICES=1
      ENDIF
      AUTO_LIST: DO IAUTO=1, N_AUTO_SLICES
      VARLOOP: DO MV=1,NV
C
      IF(AUTO_SLICE_FLAG.EQ.0)THEN
        WRITE(6,'(A,I2)') ' Enter index for variable',MV
        READ(5,*) I
      ELSE
        I=AUTO_SLICE_LISTS(IAUTO)
      ENDIF
      IS(MV) = I
      QFILE = SLCF_FILE(I)
C
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
C
      OPEN(12,FILE=QFILE,FORM='UNFORMATTED',STATUS='OLD')
C
      READ(12)
      READ(12)
      READ(12)
      READ(12) I1,I2,J1,J2,K1,K2                    
C
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
      IF (((I1.NE.I10.OR.I2.NE.I20).AND.(I10.NE.I20)) .OR.
     .    ((J1.NE.J10.OR.J2.NE.J20).AND.(J10.NE.J20)) .OR.
     .    ((K1.NE.K10.OR.K2.NE.K20).AND.(K10.NE.K20))) THEN
         WRITE(6,*) ' ERROR: Slice files are incompatible'
         STOP
         ENDIF
      ENDIF
C
      NCOUNT = 0
      READ_LOOP: DO
      READ(12,END=99) TIME
      READ(12,END=99) (((F(I,J,K),I=I1,I2),J=J1,J2),K=K1,K2)
      IF (TIME.LT.TBEG) CYCLE READ_LOOP
      IF (TIME.GT.TEND) EXIT READ_LOOP
      NCOUNT = NCOUNT + 1
C
      Q(I1:I2,J1:J2,K1:K2,MV) = Q(I1:I2,J1:J2,K1:K2,MV)+
     .                          F(I1:I2,J1:J2,K1:K2)
      ENDDO READ_LOOP

99    CLOSE(12)
C
      IF (NCOUNT.EQ.0) NCOUNT=1
C
      DO K=K1,K2
      DO J=J1,J2
      DO I=I1,I2
      Q(I,J,K,MV) = Q(I,J,K,MV)/REAL(NCOUNT)
      ENDDO
      ENDDO
      ENDDO
C
      SELECT CASE(IOR_SLCF)
      CASE(1)
      DO K=K1+1,K2
      DO J=J1+1,J2
      SUM(MV) = SUM(MV) + 0.25*(Q(I1,J,K,MV)+Q(I1,J-1,K,MV)+
     .                          Q(I1,J,K-1,MV)+Q(I1,J-1,K-1,MV))*
     .                  (M%Y(J)-M%Y(J-1))*(M%Z(K)-M%Z(K-1))
      ENDDO
      ENDDO
      CASE(2)
      DO K=K1+1,K2
      DO I=I1+1,I2
      SUM(MV) = SUM(MV) + 0.25*(Q(I,J1,K,MV)+Q(I-1,J1,K,MV)+
     .                          Q(I,J1,K-1,MV)+Q(I-1,J1,K-1,MV))*
     .                  (M%X(I)-M%X(I-1))*(M%Z(K)-M%Z(K-1))
      ENDDO
      ENDDO
      CASE(3)
      DO J=J1+1,J2
      DO I=I1+1,I2
      SUM(MV) = SUM(MV) + 0.25*(Q(I,J,K1,MV)+Q(I-1,J,K1,MV)+
     .                          Q(I,J-1,K1,MV)+Q(I-1,J-1,K1,MV))*
     .                  (M%X(I)-M%X(I-1))*(M%Y(J)-M%Y(J-1))
      ENDDO
      ENDDO
      END SELECT
C
      WRITE(6,'(A,A,A,ES12.4)') 
     .    ' Integral of ',TRIM(SLCF_TEXT(IS(MV))),' = ',SUM(MV)
C
      ENDDO VARLOOP
C*** WRITE OUT SLICE DATA IF AUTO_SLICE_FLAG IS SET
C       (WRITE IT OUT HERE RATHER THAN later SO WE DON'T HAVE TO SAVE EXTRA DATA)
c        
      IF(AUTO_SLICE_FLAG.EQ.1)THEN
        IF(IAUTO.EQ.1)THEN
          WRITE(6,*) 'Enter output file name:'
          READ(5,'(A)') OUTFILE
          OPEN(44,FILE=OUTFILE,FORM='FORMATTED',STATUS='UNKNOWN')
C
          WRITE(6,*) ' Writing to file...      ',TRIM(OUTFILE)
          ZMIN=-1000000000000.0
        ENDIF
C
        I3 = I2 - I1 + 1
        J3 = J2 - J1 + 1
        K3 = K2 - K1 + 1
C
C One-dimensional section file
C    NOTE: IN AUTO_SLICE MODE ONLY VERTICAL 1D SLICES ARE SUPPORTED
C
        IF (I1.EQ.I2 .AND. J1.EQ.J2 .AND. K1.NE.K2) then
          IF(IAUTO.EQ.1)THEN
            WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV,"(A,','),A)"
            WRITE(44,FRMT) 'Z',(TRIM(SLCF_TEXT(IS(L))),L=1,NV)
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
C
C
      CASE(3) FILETYPE
C
      BNDF_MESH = 1
      REWIND(11)
      NFILES_EXIST=0
      SEARCH_BNDF: DO I=1,500
      CALL SEARCH('BNDF',4,11,IERR)
      IF (IERR.EQ.1) EXIT SEARCH_BNDF
      BACKSPACE(11)
      IF (VERSION.LE.2.) READ(11,*) JUNK
      IF (VERSION.GT.2.) READ(11,*) JUNK,BNDF_MESH(I)
      READ(11,'(A)') BNDF_FILE(I)
      READ(11,'(A)') BNDF_TEXT(I)
      READ(11,*)
      READ(11,'(A)') BNDF_UNIT(I)
      OPEN(12,FILE=BNDF_FILE(I),FORM='UNFORMATTED',STATUS='OLD',
     .                   IOSTAT=RCODE)
      CLOSE(12)
      IF(RCODE.NE.0)THEN
        CYCLE
      ENDIF
      NFILES_EXIST=NFILES_EXIST+1
      WRITE(6,'(I3,3X,A,I2,A,A)') I,'MESH ',BNDF_MESH(I),
     .     ', ',TRIM(BNDF_TEXT(I))
      ENDDO SEARCH_BNDF
C
      IF(NFILES_EXIST.EQ.0)THEN
        WRITE(6,*)"There are no boundary files to convert"
        GOTO 500
      ENDIF
      
      WRITE(6,*) ' Enter starting and ending time for averaging (s)'
      READ(5,*) TBEG,TEND
      WRITE(6,*) ' Enter orientation: (plus or minus 1, 2 or 3)'
      READ(5,*) IOR_INPUT
      IF (VERSION.LE.2.0) IOR_INPUT = ABS(IOR_INPUT)
C
      WRITE(6,*) ' Enter number of variables'
      READ(5,*) NV
C
      BVARLOOP: DO MV=1,NV
C
      WRITE(6,'(A,I2)') ' Enter boundary file index for variable',MV
      READ(5,*) I
      IB(MV) = I
      QFILE = BNDF_FILE(I)
      OPEN(12,FILE=QFILE,FORM='UNFORMATTED',STATUS='OLD',
     .                   IOSTAT=RCODE)
      IF(RCODE.NE.0)THEN
        CLOSE(12)
        CYCLE
      ENDIF
C
      IF (MV.EQ.1) THEN
      NM = BNDF_MESH(I)
      M=>MESH(NM)
      ALLOCATE(Q(0:M%IBAR,0:M%JBAR,0:M%KBAR,NV))
      Q = 0.
      ELSE
      IF (BNDF_MESH(I).NE.NM) THEN
        WRITE(6,*) ' ERROR: All boundary files must have the same mesh'
        STOP
        ENDIF
      ENDIF
C
      READ(12)
      READ(12)
      READ(12)
      READ(12) NPATCH
C
      IF (MV.EQ.1) THEN
      ALLOCATE(IOR(1:NPATCH))
      ALLOCATE(I1B(1:NPATCH))
      ALLOCATE(I2B(1:NPATCH))
      ALLOCATE(J1B(1:NPATCH))
      ALLOCATE(J2B(1:NPATCH))
      ALLOCATE(K1B(1:NPATCH))
      ALLOCATE(K2B(1:NPATCH))
      ENDIF
C
      DO I=1,NPATCH
      IF (VERSION.LE.2.0) THEN
      READ(12) I1B(I),I2B(I),J1B(I),J2B(I),K1B(I),K2B(I)
      IF (I1B(I).EQ.I2B(I)) IOR(I) = 1
      IF (J1B(I).EQ.J2B(I)) IOR(I) = 2
      IF (K1B(I).EQ.K2B(I)) IOR(I) = 3
      ELSE
      READ(12) I1B(I),I2B(I),J1B(I),J2B(I),K1B(I),K2B(I),IOR(I)
      ENDIF
      ENDDO
C
      IF (MV.EQ.1) THEN
      IJBAR = MAX(M%IBAR,M%JBAR)
      JKBAR = MAX(M%JBAR,M%KBAR)
      ALLOCATE(F(0:IJBAR,0:JKBAR,NPATCH))
      ENDIF
      F = 0.
C
      NCOUNT = 0
      READ_BLOOP: DO
      READ(12,END=199) TIME
      DO II=1,NPATCH
      SELECT CASE(ABS(IOR(II)))
      CASE(1)
      READ(12,END=199) ((F(J,K,II),J=J1B(II),J2B(II)),K=K1B(II),K2B(II))
      CASE(2)
      READ(12,END=199) ((F(I,K,II),I=I1B(II),I2B(II)),K=K1B(II),K2B(II))
      CASE(3)
      READ(12,END=199) ((F(I,J,II),I=I1B(II),I2B(II)),J=J1B(II),J2B(II))
      END SELECT
      ENDDO
      IF (TIME.LT.TBEG) CYCLE READ_BLOOP
      IF (TIME.GT.TEND) EXIT READ_BLOOP
C
      NCOUNT = NCOUNT + 1
C
      REC_PATCH: DO II=1,NPATCH
      IF (IOR(II).NE.IOR_INPUT) CYCLE REC_PATCH
      SELECT CASE(ABS(IOR(II)))
      CASE(1)
      DO K=K1B(II),K2B(II)
      DO J=J1B(II),J2B(II)
      Q(I1B(II),J,K,MV) = Q(I1B(II),J,K,MV) + F(J,K,II)
      ENDDO
      ENDDO
      CASE(2)
      DO K=K1B(II),K2B(II)
      DO I=I1B(II),I2B(II)
      Q(I,J1B(II),K,MV) = Q(I,J1B(II),K,MV) + F(I,K,II)
      ENDDO
      ENDDO
      CASE(3)
      DO J=J1B(II),J2B(II)
      DO I=I1B(II),I2B(II)
      Q(I,J,K1B(II),MV) = Q(I,J,K1B(II),MV) + F(I,J,II)
      ENDDO
      ENDDO
      END SELECT
      ENDDO REC_PATCH
C
      ENDDO READ_BLOOP

199   CLOSE(12)
C
      DO K=0,M%KBAR
      DO J=0,M%JBAR
      DO I=0,M%IBAR
      Q(I,J,K,MV) = Q(I,J,K,MV)/REAL(NCOUNT)
      ENDDO
      ENDDO
      ENDDO
C
      ENDDO BVARLOOP
C
      END SELECT FILETYPE

C*** AUTO_SLICE DATA WAS ALREADY GENERATED
      IF(AUTO_SLICE_FLAG.EQ.1)GOTO 500
C
C Write out the data to an ASCII file
C
      WRITE(6,*) 'Enter output file name:'
      READ(5,'(A)') OUTFILE
C     WRITE(QNAME,'(A,A)') TRIM(CHID),'_fds2ascii.csv'
      OPEN(44,FILE=OUTFILE,FORM='FORMATTED',STATUS='UNKNOWN')
C
      WRITE(6,*) ' Writing to file...      ',TRIM(OUTFILE)
C
      SELECT CASE(IFILE)
C
      CASE(1:2)
C
      I3 = I2 - I1 + 1
      J3 = J2 - J1 + 1
      K3 = K2 - K1 + 1
C
C One-dimensional section file
C
      IF (I1.EQ.I2 .AND. J1.EQ.J2 .AND. K1.NE.K2) then
        WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV,"(A,','),A)"
        WRITE(44,FRMT) 'Z',(TRIM(SLCF_TEXT(IS(L))),L=1,NV)
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
        WRITE(44,FRMT) 'Y',(TRIM(SLCF_TEXT(IS(L))),L=1,NV)
        WRITE(44,FRMT) 'm',(TRIM(SLCF_UNIT(IS(L))),L=1,NV)
        WRITE(FRMT,'(A,I2.2,A)') "(",NV,"(E12.5,','),E12.5)"
        LOOP2: DO J=J1,J2,NSAM
        IF (M%Y(J).GT.YF .OR. M%Y(J).LT.YS) CYCLE LOOP2
        write(44,FRMT) M%y(j),(q(i2,j,k2,l),l=1,nv)
        enddo LOOP2
        endif

        if(i1.ne.i2.and.j1.eq.j2.and.k1.eq.k2) then
        WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV,"(A,','),A)"
        WRITE(44,FRMT) 'X',(TRIM(SLCF_TEXT(IS(L))),L=1,NV)
        WRITE(44,FRMT) 'm',(TRIM(SLCF_UNIT(IS(L))),L=1,NV)
        WRITE(FRMT,'(A,I2.2,A)') "(",NV,"(E12.5,','),E12.5)"
        LOOP3: DO I=I1,I2,NSAM
        IF (M%X(I).GT.XF .OR. M%X(I).LT.XS) CYCLE LOOP3
        write(44,FRMT) M%x(i),(q(i,j2,k2,l),l=1,nv)
        enddo LOOP3
        endif

c ... Two-dimensional section file
        if(i1.eq.i2.and.j1.ne.j2.and.k1.ne.k2) then
        WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV+1,"(A,','),A)"
        WRITE(44,FRMT) 'Y','Z',(TRIM(SLCF_TEXT(IS(L))),L=1,NV)
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
        WRITE(44,FRMT) 'X','Z',(TRIM(SLCF_TEXT(IS(L))),L=1,NV)
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
        WRITE(44,FRMT) 'X','Y',(TRIM(SLCF_TEXT(IS(L))),L=1,NV)
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

c ... Three-dimensional section file
        if (i1.ne.i2.and.j1.ne.j2.and.k1.ne.k2) then
        WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV+2,"(A,','),A)"
        WRITE(44,FRMT) 'X','Y','Z',(TRIM(SLCF_TEXT(IS(L))),L=1,NV)
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
C
      CASE(3)
C
      PATCHES: DO II=1,NPATCH
      IF (IOR(II).NE.IOR_INPUT) CYCLE PATCHES
      IF (M%X(I1B(II)).GT.XF .OR. M%X(I2B(II)).LT.XS) CYCLE PATCHES
      IF (M%Y(J1B(II)).GT.YF .OR. M%Y(J2B(II)).LT.YS) CYCLE PATCHES
      IF (M%Z(K1B(II)).GT.ZF .OR. M%Z(K2B(II)).LT.ZS) CYCLE PATCHES
      WRITE(44,'(A,I4,5(F7.2,A),F8.2)') 'Patch',II,
     .   M%X(I1B(II)),'<x<',M%X(I2B(II)),',  ',
     .   M%Y(J1B(II)),'<y<',M%Y(J2B(II)),',  ',
     .   M%Z(K1B(II)),'<z<',M%Z(K2B(II))
      WRITE(FRMT,'(A,I2.2,A)') "(1X,",NV+2,"(A,','),A)"
      WRITE(44,FRMT) 'X','Y','Z',(TRIM(BNDF_TEXT(IB(L))),L=1,NV)
      WRITE(44,FRMT) 'm','m','m',(TRIM(BNDF_UNIT(IB(L))),L=1,NV)
      WRITE(FRMT,'(A,I2.2,A)') "(",NV+2,"(E12.5,','),E12.5)"
      DO K=K1B(II),K2B(II),NSAM
      DO J=J1B(II),J2B(II),NSAM
      ILOOP: DO I=I1B(II),I2B(II),NSAM
      IF (M%X(I).GT.XF .OR. M%X(I).LT.XS) CYCLE ILOOP
      IF (M%Y(J).GT.YF .OR. M%Y(J).LT.YS) CYCLE ILOOP
      IF (M%Z(K).GT.ZF .OR. M%Z(K).LT.ZS) CYCLE ILOOP
      WRITE(44,FRMT) M%X(I),M%Y(J),M%Z(K),(Q(I,J,K,L),L=1,NV)
      ENDDO ILOOP
      ENDDO
      ENDDO
      ENDDO PATCHES
C
      END SELECT
C
      CLOSE(44)
C
      IF (IFILE.NE.1) EXIT EXTRA_PLOT3D_FILES
C
      NEW_PLOT3D = .FALSE.
      ENDDO EXTRA_PLOT3D_FILES
C
  500 CONTINUE
C
      STOP
      END


      SUBROUTINE SEARCH(STRING,LENGTH,LU,IERR)
      IMPLICIT NONE
C
      CHARACTER(*), INTENT(IN) :: STRING
      INTEGER, INTENT(OUT) :: IERR
      INTEGER, INTENT(IN) :: LU, LENGTH
      CHARACTER(20) :: JUNK
C
      SEARCH_LOOP: DO 
      READ(LU,'(A)',END=10) JUNK
      IF (JUNK(1:LENGTH).EQ.STRING(1:LENGTH)) EXIT SEARCH_LOOP
      ENDDO SEARCH_LOOP
C
      IERR = 0
      RETURN
C
   10 IERR = 1
      RETURN
C
      END SUBROUTINE SEARCH
