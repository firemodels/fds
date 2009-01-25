      PROGRAM FDS2ASCII
C
C Program to convert various FDS output files to ASCII
C
      INTEGER, PARAMETER :: FB = SELECTED_REAL_KIND(6)
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
      REAL(FB), ALLOCATABLE, DIMENSION(:,:,:,:) :: Q
      REAL(FB), ALLOCATABLE, DIMENSION(:,:,:) :: F
      LOGICAL :: NEW_PLOT3D=.TRUE.
      INTEGER IOR_SLCF
      CHARACTER(1) ANS
      CHARACTER(30) :: UNITJUNK
      CHARACTER(40) GRIDFILE,QNAME,CHID,QFILE,JUNK,FRMT,OUTFILE
      CHARACTER(40), DIMENSION(500) :: PL3D_FILE,SLCF_FILE,BNDF_FILE,
     .                                 SLCF_TEXT,BNDF_TEXT,
     .                                 SLCF_UNIT,BNDF_UNIT,
     .                                 SLCF_LABEL
      INTEGER,  DIMENSION(60) :: IB,IS
      INTEGER,  DIMENSION(500) :: PL3D_MESH,SLCF_MESH,BNDF_MESH
      REAL(FB), DIMENSION(500) :: PL3D_TIME
C
      WRITE(6,*) ' Enter Job ID string (CHID):'
      READ(5,'(a)') CHID
      GRIDFILE = TRIM(CHID)//'.smv'
C
C Open grid file
C
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
         XS = -100000.
         XF =  100000.
         YS = -100000.
         YF =  100000.
         ZS = -100000.
         ZF =  100000.
C
      BNDF_MESH = 1
      REWIND(11)
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
      WRITE(6,'(I3,3X,A,I2,A,A)') I,'MESH ',BNDF_MESH(I),
     .     ', ',TRIM(BNDF_TEXT(I))
      ENDDO SEARCH_BNDF
C
      WRITE(6,*) ' Enter starting and ending time for averaging (s)'
      READ(5,*) TBEG,TEND
C
      NV=1
C
      BVARLOOP: DO MV=1,NV
C
      WRITE(6,'(A,I2)') ' Enter boundary file index for variable',MV
      READ(5,*) I
      IB(MV) = I
      QFILE = BNDF_FILE(I)
      OPEN(12,FILE=QFILE,FORM='UNFORMATTED',STATUS='OLD')
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
      READ(12) I1B(I),I2B(I),J1B(I),J2B(I),K1B(I),K2B(I),IOR(I)
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
C Write out the data to an ASCII file
C
      WRITE(6,*) 'Enter output file name:'
      READ(5,'(A)') OUTFILE
C     WRITE(QNAME,'(A,A)') TRIM(CHID),'_fds2ascii.csv'
      OPEN(44,FILE=OUTFILE,FORM='FORMATTED',STATUS='UNKNOWN')
C
      WRITE(6,*) ' Writing to file...      ',OUTFILE
C
      WRITE(44,'(A)') 'Height,Left 29,Left 16,Left 3,
     .Right 3,Right 16,Right 29'
      DO K=1,M%KBAR
      WRITE(44,"(6(E12.5,','),E12.5)") M%Z(K),Q(15,0,K,1),Q(8,0,K,1),
     . Q(2,0,K,1),Q(0,2,K,1),Q(0,8,K,1),Q(0,15,K,1)
      ENDDO
C
      CLOSE(44)
C
  500 CONTINUE
C
      STOP
      END


      SUBROUTINE SEARCH(STRING,LENGTH,LU,IERR)
C
      CHARACTER(*), INTENT(IN) :: STRING
      INTEGER, INTENT(OUT) :: IERR
      CHARACTER(20) :: JUNK
      INTEGER LU,LENGTH
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
