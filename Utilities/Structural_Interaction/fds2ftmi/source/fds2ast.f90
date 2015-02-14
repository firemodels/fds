SUBROUTINE FDS2AST (CHID,NOMASTER,C_S,TBEG,TEND,TINT,COUNTER,VARIABLE,N)

! Program to read the Adiabatic Surface Temperature from bndf files for a specific 
! location, based on fds2ascii (part of FDS package)
! This program is part of a Coupling Procedure between fire and structural models
!                                                                       - SILVA, JC

! $Id$  
! $Revision$
! $Date$
 
IMPLICIT NONE

INTEGER, PARAMETER :: FB = SELECTED_REAL_KIND(6)
INTEGER :: IERR, NMESHES, NM, NOC, I, J, K, IOR_LOOP
INTEGER :: IDUM, IFILE, NSAM, NV, MV
!INTEGER :: I1, I2, J1, J2, K1, K2, I3, J3, K3
INTEGER :: IOR_INPUT, NPATCH, IJBAR, JKBAR,II
REAL(FB) :: XS, XF, YS, YF, ZS, ZF, TIME

TYPE MESH_TYPE
   REAL(FB), POINTER, DIMENSION(:) :: X,Y,Z
   REAL(FB) :: D1,D2,D3,D4
   INTEGER :: IBAR,JBAR,KBAR,IERR
END TYPE MESH_TYPE
 
TYPE (MESH_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: MESH
TYPE (MESH_TYPE), POINTER :: M
 
INTEGER, ALLOCATABLE, DIMENSION(:) :: IOR,I1B,I2B,J1B,J2B,K1B,K2B
REAL(FB), ALLOCATABLE, DIMENSION(:,:,:,:) :: Q
REAL(FB), ALLOCATABLE, DIMENSION(:,:,:) :: F
LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: ALREADY_USED
CHARACTER(2) :: ANS
CHARACTER(4) :: CHOICE
CHARACTER(256) GRIDFILE,CHID,QFILE,JUNK
CHARACTER(256), DIMENSION(500) :: BNDF_FILE,BNDF_TEXT,BNDF_UNIT
CHARACTER(20), DIMENSION(500) :: BNDF_TYPE
CHARACTER(20) :: BNDF_TYPE_CHOSEN
!REAL(FB), DIMENSION(500) :: X1, X2, Y1, Y2, Z1, Z2
INTEGER,  DIMENSION(60) :: IB
INTEGER,  DIMENSION(500) :: BNDF_MESH
INTEGER :: RCODE
INTEGER :: NFILES_EXIST
LOGICAL :: EXISTS
INTEGER :: BATCHMODE
!CHARACTER(256) :: BUFFER
!INTEGER :: ERROR_STATUS
!CHARACTER(256) :: ARG

INTEGER COUNTER
INTEGER V,SIZE,VCOUNT(7),P,POND,VARIABLE,WARNING,IN,VAR
REAL Xa,Ya,Za,C_SIZE,TINT,A,B,C,C_S,NX,NY,NZ
REAL TBEG,TEND
REAL NOMASTER(COUNTER,4),N(COUNTER,4)
REAL, ALLOCATABLE, DIMENSION(:,:) :: M_AST
REAL, ALLOCATABLE, DIMENSION(:) :: V_TIME,TAST
CHARACTER(8) OUTFILE
CHARACTER(30) VARIABLE_KIND
!REAL MED

!REAL(FB) :: tb_init,te_init,te_nmeshes,tb_nmeshes,te_fds2ast,tb_fds2ast

!OPEN(50,FILE='cputime.dat',FORM='FORMATTED',STATUS='UNKNOWN')
!CALL CPU_TIME (tb_init)

! Set a few default values
BATCHMODE=0

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
!   DO I=0,M%IBAR
!   M%X(I)=0
!   ENDDO
!   DO I=0,M%JBAR
!   M%Y(I)=0
!   ENDDO
!   DO I=0,M%KBAR
!   M%Z(I)=0
!   ENDDO
    
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
 
IFILE=3 
NSAM=1 
ANS='y'
CALL TOUPPER(ANS,ANS)
! **************************************************
! Beggining of the inputs for the new code - SILVA, JC
!   WRITE(6,*) ' Enter XYZ position'
!   READ(LU_IN,*) X,Y,Z
!   WRITE(6,*) ' Enter Cell Size'
!   READ(LU_IN,*) C_SIZE
            
         BNDF_MESH = 1
         REWIND(11)
         NFILES_EXIST=0

         SEARCH_BNDF: DO I=1,500
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
            !IF (BATCHMODE.EQ.0) write(6,'(I3,3X,A,I2,A,A)') I,'MESH ',BNDF_MESH(I), ', ',TRIM(BNDF_TEXT(I))
         ENDDO SEARCH_BNDF
    
         IF (NFILES_EXIST.EQ.0)THEN
            IF(BATCHMODE.EQ.0) write(6,*)"There are no boundary files to convert"
            STOP
         ENDIF
         
! Rest of the input parameters - SILVA, JC         
!         IF (BATCHMODE.EQ.0) write(6,*) ' Enter starting and ending time (s)'
!         READ(LU_IN,*) TBEG,TEND
!         IF (BATCHMODE.EQ.0) write(6,*) ' Enter time interval (s)'
!         READ(LU_IN,*) TINT
!         IF (BATCHMODE.EQ.0) write(6,*) ' Enter output file name:'
!         READ(LU_IN,'(A)') OUTFILE

!CALL CPU_TIME (te_init)
!WRITE (50,*) 'Time of initialization was ', te_init-tb_init, ' seconds'

LOOP_FDS2AST:DO P=1,COUNTER
!CALL CPU_TIME (tb_fds2ast)
    Xa=NOMASTER(P,2)
    Ya=NOMASTER(P,3)
    Za=NOMASTER(P,4)
    WRITE (OUTFILE,'(g8.0)') INT(NOMASTER(P,1))
    PRINT *, OUTFILE
    C_SIZE=C_S
    NX=N(P,2)
    NY=N(P,3)
    NZ=N(P,4)

         
!$omp parallel 
!$omp do
VARIABLE_NUMBER: DO VAR=1,VARIABLE
! 500 - Point of the loop in case of the point does not reach a meaning value - SILVA, JC         
500 CONTINUE
   XS = Xa-(C_SIZE/2)
   XF = Xa+(C_SIZE/2)
   YS = Ya-(C_SIZE/2)
   YF = Ya+(C_SIZE/2)
   ZS = Za-(C_SIZE/2)
   ZF = Za+(C_SIZE/2)

! SIZE is the lenght of the vectors - SILVA, JC
             SIZE=((TEND-TBEG)/TINT)+1
             ALLOCATE (V_TIME(SIZE))
             ALLOCATE (TAST(SIZE))
             ALLOCATE (M_AST(SIZE,7))
!             DO I=1,SIZE
!                V_TIME(I)=0.0
!                TAST(I)=0.0
!                DO J=1,7
!                M_AST(I,J)=0.0
!                ENDDO
!             ENDDO
             IF (VAR.EQ.1) VCOUNT=0.0
             V_TIME=0.0
             TAST=0.0
             M_AST=0.0

             
! IN_LOOP is a loop at index file - SILVA, JC
IN_LOOP: DO IN=VAR,VARIABLE*NMESHES,VARIABLE

!CALL CPU_TIME (tb_nmeshes)
!*********
         !IF (BATCHMODE.EQ.0) write(6,*) ' Enter orientation: (plus or minus 1, 2 or 3)'
         !READ(LU_IN,*) IOR_INPUT
         IOR_LOOPING: DO IOR_LOOP=1,7
            IOR_INPUT=IOR_LOOP-4
            IF (IOR_INPUT.EQ.0) CYCLE
            IF (IOR_INPUT.EQ.-3 .AND. NZ.GE.0) CYCLE
            IF (IOR_INPUT.EQ.-2 .AND. NY.GE.0) CYCLE
            IF (IOR_INPUT.EQ.-1 .AND. NX.GE.0) CYCLE
            IF (IOR_INPUT.EQ.1 .AND. NX.LE.0) CYCLE
            IF (IOR_INPUT.EQ.2 .AND. NY.LE.0) CYCLE
            IF (IOR_INPUT.EQ.3 .AND. NZ.LE.0) CYCLE

!***********************
            NV=1
            MV=1
            IB(MV) = IN
            QFILE = BNDF_FILE(IN)
            OPEN(12+VAR,FILE=QFILE,FORM='UNFORMATTED',STATUS='OLD', IOSTAT=RCODE)
            IF (RCODE.NE.0) THEN
              CLOSE(12+VAR)
              CYCLE
            ENDIF

            NM = BNDF_MESH(IN)
            M=>MESH(NM)
            ALLOCATE(Q(0:M%IBAR,0:M%JBAR,0:M%KBAR,NV))
            Q = 0.0
!            DO I=0,M%IBAR
!                DO J=0,M%JBAR
!                    DO K=0,M%KBAR
!                    Q(I,J,K,NV)=0.0
!                    ENDDO
!                ENDDO
!             ENDDO
            ALLOCATE(ALREADY_USED(0:M%IBAR,0:M%JBAR,0:M%KBAR))
            BNDF_TYPE_CHOSEN = BNDF_TYPE(IN)

            READ(12+VAR)
            READ(12+VAR)
            READ(12+VAR)
            READ(12+VAR) NPATCH
    
            ALLOCATE(IOR(1:NPATCH))
            ALLOCATE(I1B(1:NPATCH))
            ALLOCATE(I2B(1:NPATCH))
            ALLOCATE(J1B(1:NPATCH))
            ALLOCATE(J2B(1:NPATCH))
            ALLOCATE(K1B(1:NPATCH))
            ALLOCATE(K2B(1:NPATCH))
    
            DO I=1,NPATCH
               READ(12+VAR) I1B(I),I2B(I),J1B(I),J2B(I),K1B(I),K2B(I),IOR(I)
            ENDDO
    
            IJBAR = MAX(M%IBAR,M%JBAR)
            JKBAR = MAX(M%JBAR,M%KBAR)
            ALLOCATE(F(0:IJBAR,0:JKBAR,NPATCH))
!            F = 0.
!            DO I=0,IJBAR
!                DO J=0,JKBAR
!                    DO K=1,NPATCH
!                    F(I,J,K)=0.0
!                    ENDDO
!                ENDDO
!             ENDDO

            READ_BLOOP: DO
               READ(12+VAR,END=199) TIME
               DO II=1,NPATCH
                  SELECT CASE(ABS(IOR(II)))
                   CASE(1)
                   IF (BNDF_TYPE_CHOSEN=='STAGGERED')READ(12+VAR,END=199) ((F(J,K,II),J=J1B(II),J2B(II)),K=K1B(II),K2B(II))
                   IF (BNDF_TYPE_CHOSEN=='CENTERED') READ(12+VAR,END=199) ((F(J,K,II),J=J1B(II)+1,J2B(II)+1),K=K1B(II)+1,K2B(II)+1)
                   CASE(2)
                   IF (BNDF_TYPE_CHOSEN=='STAGGERED')READ(12+VAR,END=199) ((F(I,K,II),I=I1B(II),I2B(II)),K=K1B(II),K2B(II))
                   IF (BNDF_TYPE_CHOSEN=='CENTERED') READ(12+VAR,END=199) ((F(I,K,II),I=I1B(II)+1,I2B(II)+1),K=K1B(II)+1,K2B(II)+1)
                   CASE(3)
                   IF (BNDF_TYPE_CHOSEN=='STAGGERED')READ(12+VAR,END=199) ((F(I,J,II),I=I1B(II),I2B(II)),J=J1B(II),J2B(II))
                   IF (BNDF_TYPE_CHOSEN=='CENTERED') READ(12+VAR,END=199) ((F(I,J,II),I=I1B(II)+1,I2B(II)+1),J=J1B(II)+1,J2B(II)+1)
                  END SELECT
               ENDDO
               IF (TIME.LT.TBEG) CYCLE READ_BLOOP
               IF (INT(TIME).GT.TEND) EXIT READ_BLOOP
                      
               ALREADY_USED = .FALSE.
!                DO I=0,M%IBAR
!                    DO J=0,M%JBAR
!                        DO K=0,M%KBAR
!                        ALREADY_USED(I,J,K) = .FALSE.
!                        ENDDO
!                    ENDDO
!                ENDDO

               REC_PATCH: DO II=1,NPATCH
                  IF (IOR(II).NE.IOR_INPUT) CYCLE REC_PATCH
                  IF (M%X(I1B(II)).GT.XF .OR. M%X(I2B(II)).LT.XS) CYCLE REC_PATCH
                  IF (M%Y(J1B(II)).GT.YF .OR. M%Y(J2B(II)).LT.YS) CYCLE REC_PATCH
                  IF (M%Z(K1B(II)).GT.ZF .OR. M%Z(K2B(II)).LT.ZS) CYCLE REC_PATCH
                  SELECT CASE(ABS(IOR(II)))
                     CASE(1)
                        DO K=K1B(II),K2B(II)
                           DO J=J1B(II),J2B(II)
                              IF (.NOT. ALREADY_USED(I1B(II),J,K)) THEN
!                                 Q(I1B(II),J,K,MV) = Q(I1B(II),J,K,MV) + F(J,K,II)
                                 Q(I1B(II),J,K,MV) = F(J,K,II)
                                 ALREADY_USED(I1B(II),J,K) = .TRUE.
                              ENDIF
                           ENDDO
                        ENDDO
                     CASE(2)
                        DO K=K1B(II),K2B(II)
                           DO I=I1B(II),I2B(II)
                              IF (.NOT. ALREADY_USED(I,J1B(II),K)) THEN
!                                 Q(I,J1B(II),K,MV) = Q(I,J1B(II),K,MV) + F(I,K,II)
                                 Q(I,J1B(II),K,MV) = F(I,K,II)
                                 ALREADY_USED(I,J1B(II),K) = .TRUE.
                              ENDIF
                           ENDDO
                        ENDDO
                     CASE(3)
                        DO J=J1B(II),J2B(II)
                           DO I=I1B(II),I2B(II)
                              IF (.NOT. ALREADY_USED(I,J,K1B(II))) THEN
!                                 Q(I,J,K1B(II),MV) = Q(I,J,K1B(II),MV) + F(I,J,II)
                                 Q(I,J,K1B(II),MV) =  F(I,J,II)
                                 ALREADY_USED(I,J,K1B(II)) = .TRUE.
                              ENDIF
                           ENDDO
                        ENDDO
                  END SELECT
               ENDDO REC_PATCH

            V=(TIME/TINT)+1
            V_TIME(V)=TIME
            PATCHES: DO II=1,NPATCH
               DO K=K1B(II),K2B(II),NSAM
                  DO J=J1B(II),J2B(II),NSAM
                     DO I=I1B(II),I2B(II),NSAM
                        IF (M%X(I).GT.XF .OR. M%X(I).LT.XS) CYCLE 
                        IF (M%Y(J).GT.YF .OR. M%Y(J).LT.YS) CYCLE
                        IF (M%Z(K).GT.ZF .OR. M%Z(K).LT.ZS) CYCLE
                        M_AST(V,IOR_LOOP)=M_AST(V,IOR_LOOP)+Q(I,J,K,NV)
                        IF (VAR.EQ.1 .AND. TIME.EQ.TBEG .AND. Q(I,J,K,NV).GT.0) VCOUNT(IOR_LOOP)=VCOUNT(IOR_LOOP)+1
                     ENDDO 
                  ENDDO
               ENDDO
             ENDDO PATCHES
     
            ENDDO READ_BLOOP
      
               199 CLOSE(12+VAR)
    
            DEALLOCATE(Q)
            DEALLOCATE(ALREADY_USED)
            DEALLOCATE(IOR)
            DEALLOCATE(I1B)
            DEALLOCATE(I2B)
            DEALLOCATE(J1B)
            DEALLOCATE(J2B)
            DEALLOCATE(K1B)
            DEALLOCATE(K2B)
            DEALLOCATE(F)

         ENDDO IOR_LOOPING

    ENDDO IN_LOOP 

!CALL CPU_TIME (te_nmeshes)
!WRITE (50,*) 'Time of internal loop was ', te_nmeshes-tb_nmeshes, ' seconds'
!************************************
! Writing the vector time (V_TIME) and tast matrix - SILVA, JC
IF (VAR.EQ.2) THEN
    OPEN(100+VAR,FILE=TRIM(OUTFILE)//'H.dat',FORM='FORMATTED',STATUS='UNKNOWN')
ELSE
    OPEN(100+VAR,FILE=TRIM(OUTFILE)//'.dat',FORM='FORMATTED',STATUS='UNKNOWN')
ENDIF
!
!WRITE (6,'(F8.3,F8.3,F8.3,F8.3,F8.3,F8.3,F8.3)') C_SIZE,M_AST(1,1),M_AST(1,2),M_AST(1,3),M_AST(1,5), &
!M_AST(1,6),M_AST(1,7)
!WRITE (6,'(I8,I8,I8,I8,I8,I8)') VCOUNT(1),VCOUNT(2),VCOUNT(3),VCOUNT(5),VCOUNT(6),VCOUNT(7)
!
IF (VCOUNT(1).GT.1) THEN
    DO I=1,V
        M_AST(I,1)=M_AST(I,1)/VCOUNT(1)
    END DO
END IF
IF (VCOUNT(2).GT.1) THEN
    DO I=1,V
        M_AST(I,2)=M_AST(I,2)/VCOUNT(2)
    END DO
END IF
IF (VCOUNT(3).GT.1) THEN
    DO I=1,V
        M_AST(I,3)=M_AST(I,3)/VCOUNT(3)
    END DO
END IF
IF (VCOUNT(5).GT.1) THEN
    DO I=1,V
        M_AST(I,5)=M_AST(I,5)/VCOUNT(5)
    END DO
END IF
IF (VCOUNT(6).GT.1) THEN
    DO I=1,V
        M_AST(I,6)=M_AST(I,6)/VCOUNT(6)
    END DO
END IF
IF (VCOUNT(7).GT.1) THEN
    DO I=1,V
        M_AST(I,7)=M_AST(I,7)/VCOUNT(7)
    END DO
END IF
!
VARIABLE_KIND=TRIM(BNDF_TEXT(VAR))
IF (VARIABLE_KIND==' ADIABATIC SURFACE TEMPERATURE') THEN
!
        IF (NZ.GT.0) THEN
            DO I=1,V
            M_AST(I,1)=0.D0
            END DO
        END IF
        IF (NZ.LT.0) THEN 
            DO I=1,V
            M_AST(I,7)=0.D0
            END DO
        END IF
        IF (NZ.EQ.0) THEN 
            DO I=1,V
            M_AST(I,1)=0.D0
            M_AST(I,7)=0.D0
            END DO
        END IF
!        
        IF (NY.GT.0) THEN
            DO I=1,V
            M_AST(I,2)=0.D0
            END DO
        END IF
        IF (NY.LT.0) THEN 
            DO I=1,V
            M_AST(I,6)=0.D0
            END DO
        END IF
        IF (NY.EQ.0) THEN 
            DO I=1,V
            M_AST(I,2)=0.D0
            M_AST(I,6)=0.D0
            END DO
        END IF
!
        IF (NX.GT.0) THEN
            DO I=1,V
            M_AST(I,3)=0.D0
            END DO
        END IF
        IF (NX.LT.0) THEN 
            DO I=1,V
            M_AST(I,5)=0.D0
            END DO
        END IF
        IF (NX.EQ.0) THEN 
            DO I=1,V
            M_AST(I,3)=0.D0
            M_AST(I,5)=0.D0
            END DO
        END IF
!
!
        IF (M_AST(2,7).LT.19 .AND. M_AST(2,7).GT.0) THEN
            DO I=1,V
            M_AST(I,7)=20.0
            END DO
        END IF 
        IF (M_AST(2,1).LT.19 .AND. M_AST(2,1).GT.0) THEN
            DO I=1,V
            M_AST(I,1)=20.0
            END DO
        END IF 
        IF (M_AST(2,2).LT.19 .AND. M_AST(2,2).GT.0) THEN
            DO I=1,V
            M_AST(I,2)=20.0
            END DO
        END IF 
        IF (M_AST(2,6).LT.19 .AND. M_AST(2,6).GT.0) THEN
            DO I=1,V
            M_AST(I,6)=20.0
            END DO
        END IF 
        IF (M_AST(2,5).LT.19 .AND. M_AST(2,5).GT.0) THEN
            DO I=1,V
            M_AST(I,5)=20.0
            END DO
        END IF 
        IF (M_AST(2,3).LT.19 .AND. M_AST(2,3).GT.0) THEN
            DO I=1,V
            M_AST(I,3)=20.0
            END DO
        END IF 
        IF (M_AST(1,7).LT.20 .AND. M_AST(1,7).GT.0) THEN
            M_AST(1,7)=20.0
        END IF 
        IF (M_AST(1,1).LT.20 .AND. M_AST(1,1).GT.0) THEN
            M_AST(1,1)=20.0
        END IF 
        IF (M_AST(1,2).LT.20 .AND. M_AST(1,2).GT.0) THEN
            M_AST(1,2)=20.0
        END IF 
        IF (M_AST(1,6).LT.20 .AND. M_AST(1,6).GT.0) THEN
            M_AST(1,6)=20.0
        END IF 
        IF (M_AST(1,5).LT.20 .AND. M_AST(1,5).GT.0) THEN
            M_AST(1,5)=20.0
        END IF 
        IF (M_AST(1,3).LT.20 .AND. M_AST(1,3).GT.0) THEN
            M_AST(1,3)=20.0
        END IF 
!         
C=SQRT((NX**2)+(NY**2)+(NZ**2))
!
        IF (M_AST(1,1).EQ.20 .AND. M_AST(1,2).EQ.20 .AND. M_AST(1,3).EQ.20 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.20 .AND. M_AST(1,2).EQ.20 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.20 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.20 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.20 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.20 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.20 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.20 .AND. M_AST(1,6).EQ.20 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.20 .AND. M_AST(1,2).EQ.20 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.20 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.20 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF  
        IF (M_AST(1,1).EQ.20 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.20 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF  
        IF (M_AST(1,1).EQ.20 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.20 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF  
        IF (M_AST(1,1).EQ.20 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ))
            B=SQRT((M_AST(I,1)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.20 .AND. M_AST(1,3).EQ.20 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.20) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.20 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.20 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.20) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.20 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.20 .AND. M_AST(1,7).EQ.20) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.20 .AND. M_AST(1,6).EQ.20 .AND. M_AST(1,7).EQ.20) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.20 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.20) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.20 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.20) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF  
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.20 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.20) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF  
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.20 .AND. M_AST(1,7).EQ.20) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF  
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.20) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ))
            B=SQRT((M_AST(I,7)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.20 .AND. M_AST(1,3).EQ.20 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.20 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.20 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.20 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.20 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.20 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.20 .AND. M_AST(1,6).EQ.20 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.20 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((M_AST(I,6)*NY))
            B=SQRT((M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.20 .AND. &
        M_AST(1,5).EQ.0 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(1,1).EQ.0 .AND. M_AST(1,2).EQ.0 .AND. M_AST(1,3).EQ.0 .AND. &
        M_AST(1,5).EQ.20 .AND. M_AST(1,6).EQ.0 .AND. M_AST(1,7).EQ.0) THEN
            DO I=1,V
            A=((M_AST(I,5)*NX))
            B=SQRT((M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
!********************************
        IF (TAST(1).GT.0) THEN
!            TAST(1)=20.0
            ELSE 
            IF (C_SIZE.GE.(2.5*C_S)) THEN
                WRITE (6,'(F8.3,F8.3,F8.3,F8.3,F8.3,F8.3,F8.3)') C_SIZE,M_AST(1,1),M_AST(1,2),M_AST(1,3),M_AST(1,5), &
                M_AST(1,6),M_AST(1,7)
                WRITE (6,'(F8.5,F8.5,F8.5,F8.5,F8.5,F8.5)') NX,NY,NZ
121             WRITE(6,*) ' WARNING: C_SIZE has reached the maximum for',OUTFILE,'!',C_SIZE
                WRITE(6,*) ' (1) if you want to ignore this node'
                WRITE(6,*) ' (2) if you want to stop the code'
                WRITE(6,*) ' (3) if you want to continue increasing C_SIZE'
                READ(*,*) WARNING
                IF (WARNING.EQ.1) THEN
                    DEALLOCATE (V_TIME)
                    DEALLOCATE (TAST)
                    DEALLOCATE (M_AST)
                    GO TO 100
                END IF
                IF (WARNING.EQ.2) STOP
                IF (WARNING.NE.1 .AND. WARNING.NE.2 .AND. WARNING.NE.3) THEN
                WRITE(6,*) ' WRONG ANSWER! '
                GO TO 121
                END IF
            END IF
            DEALLOCATE (V_TIME)
            DEALLOCATE (TAST)
            DEALLOCATE (M_AST)
            C_SIZE=C_SIZE*1.1
            GOTO 500
        END IF            
END IF
!***************************************
IF (VARIABLE_KIND==' WALL TEMPERATURE             ') THEN
        POND=0
        IF (M_AST(2,1).GT.20) POND=POND+1
        IF (M_AST(2,2).GT.20) POND=POND+1
        IF (M_AST(2,3).GT.20) POND=POND+1
        IF (M_AST(2,5).GT.20) POND=POND+1
        IF (M_AST(2,6).GT.20) POND=POND+1
        IF (M_AST(2,7).GT.20) POND=POND+1
            IF (POND.NE.0) THEN
                DO I=1,V
                TAST(I)=(M_AST(I,1)+M_AST(I,2)+M_AST(I,3)+M_AST(I,4)+M_AST(I,5)+M_AST(I,6))/POND
                ENDDO
            END IF
!*************
        IF (TAST(1).GE.20) THEN
            TAST(1)=20.0
        ELSE 
            IF (C_SIZE.GE.(2.5*C_S)) THEN
122             WRITE(6,*) ' WARNING: C_SIZE has reached the maximum for',OUTFILE,'!',C_SIZE
                WRITE(6,*) ' (1) if you want to ignore this node'
                WRITE(6,*) ' (2) if you want to stop the code'
                WRITE(6,*) ' (3) if you want to continue increasing C_SIZE'
                READ(*,*) WARNING
                IF (WARNING.EQ.1) THEN
                    DEALLOCATE (V_TIME)
                    DEALLOCATE (TAST)
                    DEALLOCATE (M_AST)
                    GO TO 100
                END IF
                IF (WARNING.EQ.2) STOP
                IF (WARNING.NE.1 .AND. WARNING.NE.2 .AND. WARNING.NE.3) THEN
                WRITE(6,*) ' WRONG ANSWER! '
                GO TO 122
                END IF
            END IF
        DEALLOCATE (V_TIME)
        DEALLOCATE (TAST)
        DEALLOCATE (M_AST)
        C_SIZE=C_SIZE*1.1
        GOTO 500
        END IF 
END IF

IF (VARIABLE_KIND==' HEAT TRANSFER COEFFICIENT    ') THEN
        IF (NZ.GT.0) THEN
            DO I=1,V
            M_AST(I,1)=0.D0
            END DO
        END IF
        IF (NZ.LT.0) THEN 
            DO I=1,V
            M_AST(I,7)=0.D0
            END DO
        END IF
        IF (NZ.EQ.0) THEN 
            DO I=1,V
            M_AST(I,1)=0.D0
            M_AST(I,7)=0.D0
            END DO
        END IF
!        
        IF (NY.GT.0) THEN
            DO I=1,V
            M_AST(I,2)=0.D0
            END DO
        END IF
        IF (NY.LT.0) THEN 
            DO I=1,V
            M_AST(I,6)=0.D0
            END DO
        END IF
        IF (NY.EQ.0) THEN 
            DO I=1,V
            M_AST(I,2)=0.D0
            M_AST(I,6)=0.D0
            END DO
        END IF
!
        IF (NX.GT.0) THEN
            DO I=1,V
            M_AST(I,3)=0.D0
            END DO
        END IF
        IF (NX.LT.0) THEN 
            DO I=1,V
            M_AST(I,5)=0.D0
            END DO
        END IF
        IF (NX.EQ.0) THEN 
            DO I=1,V
            M_AST(I,3)=0.D0
            M_AST(I,5)=0.D0
            END DO
        END IF
!
C=SQRT((NX**2)+(NY**2)+(NZ**2))
!
        IF (M_AST(10,1).GT.0 .AND. M_AST(10,2).GT.0 .AND. M_AST(10,3).GT.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).GT.0 .AND. M_AST(10,2).GT.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).GT.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).GT.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).GT.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).GT.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).GT.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).GT.0 .AND. M_AST(10,6).GT.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).GT.0 .AND. M_AST(10,2).GT.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).GT.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).GT.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF  
        IF (M_AST(10,1).GT.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).GT.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF  
        IF (M_AST(10,1).GT.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).GT.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF  
        IF (M_AST(10,1).GT.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,1)*NZ))
            B=SQRT((M_AST(I,1)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).GT.0 .AND. M_AST(10,3).GT.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).GT.0) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).GT.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).GT.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).GT.0) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).GT.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).GT.0 .AND. M_AST(10,7).GT.0) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).GT.0 .AND. M_AST(10,6).GT.0 .AND. M_AST(10,7).GT.0) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).GT.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).GT.0) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).GT.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).GT.0) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF  
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).GT.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).GT.0) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF  
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).GT.0 .AND. M_AST(10,7).GT.0) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF  
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).GT.0) THEN
            DO I=1,V
            A=((M_AST(I,7)*NZ))
            B=SQRT((M_AST(I,7)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).GT.0 .AND. M_AST(10,3).GT.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).GT.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).GT.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).GT.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).GT.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).GT.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).GT.0 .AND. M_AST(10,6).GT.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).GT.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((M_AST(I,6)*NY))
            B=SQRT((M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).GT.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).GT.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).EQ.0) THEN
            DO I=1,V
            A=((M_AST(I,5)*NX))
            B=SQRT((M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
            END DO
        END IF
        IF (M_AST(10,1).EQ.0 .AND. M_AST(10,2).EQ.0 .AND. M_AST(10,3).EQ.0 .AND. &
        M_AST(10,5).EQ.0 .AND. M_AST(10,6).EQ.0 .AND. M_AST(10,7).EQ.0) THEN
            TAST=0.D0
        END IF
!********************************
        IF (TAST(2).GE.0 .AND. TAST(2).LT.1000) THEN
            TAST(1)=0.0
        ELSE 
            IF (C_SIZE.GE.(2.5*C_S)) THEN
123             WRITE(6,*) ' WARNING: C_SIZE has reached the maximum for',OUTFILE,'!',C_SIZE
                WRITE(6,*) ' (1) if you want to ignore this node'
                WRITE(6,*) ' (2) if you want to stop the code'
                WRITE(6,*) ' (3) if you want to continue increasing C_SIZE'
                READ(*,*) WARNING
                IF (WARNING.EQ.1) THEN
                    DEALLOCATE (V_TIME)
                    DEALLOCATE (TAST)
                    DEALLOCATE (M_AST)
                    GO TO 100
                END IF
                IF (WARNING.EQ.2) STOP
                IF (WARNING.NE.1 .AND. WARNING.NE.2 .AND. WARNING.NE.3) THEN
                WRITE(6,*) ' WRONG ANSWER! '
                GO TO 123
                END IF
            END IF
        DEALLOCATE (V_TIME)
        DEALLOCATE (TAST)
        DEALLOCATE (M_AST)
        C_SIZE=C_SIZE*1.1
        GOTO 500
        END IF
END IF
!****************
!CALL CPU_TIME (te_fds2ast)
!WRITE (50,*) 'Time of FDS2AST was ', te_fds2ast-tb_fds2ast, ' seconds'
    DO I=1,V
    WRITE (100+VAR,'(E12.5)', ADVANCE='NO') V_TIME(I)
    WRITE (100+VAR,'(E12.5)', ADVANCE='YES') TAST(I)
    ENDDO
!**********************
!*** Set an average value for a steady simulation (last 50 results)
!    MED=0.0
!    DO I=0,49
!      MED=MED+TAST(V-I)
!    ENDDO
!    MED=MED/(49+1)
!    WRITE (100+VAR,'(E12.5)', ADVANCE='NO') V_TIME(V)+TINT
!    WRITE (100+VAR,'(E12.5)', ADVANCE='YES') MED          
!    WRITE (100+VAR,'(E12.5)', ADVANCE='NO') 18000.0
!    WRITE (100+VAR,'(E12.5)', ADVANCE='YES') MED      
!**********************
    IF (VARIABLE_KIND==' ADIABATIC SURFACE TEMPERATURE') THEN
      WRITE (6,'(A, F8.3)') 'Cell_Size  ',C_SIZE
      WRITE (6,'(A, F8.3, F8.3, F8.3, F8.3, F8.3, F8.3)') 'AST_Matrix ',M_AST(1,1),M_AST(1,2),M_AST(1,3),M_AST(1,5), &
      M_AST(1,6),M_AST(1,7)
      WRITE (6,'(A, F8.5, F8.5, F8.5)') 'Normal_Vector',NX,NY,NZ
    END IF
    DEALLOCATE (V_TIME)
    DEALLOCATE (TAST)
    DEALLOCATE (M_AST)
    CLOSE (100+VAR)
100 CONTINUE
    ENDDO VARIABLE_NUMBER
!$omp end do
!$omp end parallel
    ENDDO LOOP_FDS2AST
RETURN
END SUBROUTINE FDS2AST

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

! *********************** SEARCH2 *******************************

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
