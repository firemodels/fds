SUBROUTINE FDS2AST (CHID,NOMASTER,C_S,TBEG,TEND,TINT,COUNTER,VARIABLE,N,N_AVERAGE,T_AVERAGE)

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
INTEGER,  DIMENSION(60) :: IB
INTEGER,  DIMENSION(500) :: BNDF_MESH
INTEGER :: RCODE
INTEGER :: NFILES_EXIST
LOGICAL :: EXISTS
INTEGER :: BATCHMODE

INTEGER COUNTER
INTEGER V,SIZE,VCOUNT(7),P,POND,VARIABLE,VARIAB,WARNING,IN,VAR,BEG_VAR
REAL Xa,Ya,Za,C_SIZE,TINT,A,B,C,C_S,NX,NY,NZ
REAL TBEG,TEND,SUM
REAL NOMASTER(COUNTER,4),N(COUNTER,4)
REAL, ALLOCATABLE, DIMENSION(:,:) :: M_AST
REAL, ALLOCATABLE, DIMENSION(:) :: V_TIME,TAST,MED_TAST
CHARACTER(8) OUTFILE
CHARACTER(30) VARIABLE_KIND
REAL MED
INTEGER N_AVERAGE,T_AVERAGE,I_AVERAGE,COMEBACK,NUM_INT_CSIZE(7),H_NULL(3),N_DIR

! Set a few default values
BATCHMODE=0
VARIAB=VARIABLE

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
IF (IERR==1) THEN
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
!   
   CALL SEARCH('TRNX',4,11,IERR)
   READ(11,*) NOC
   DO I=1,NOC
      READ(11,*)
      ENDDO
      DO I=0,M%IBAR
      READ(11,*) IDUM,M%X(I)
   ENDDO
! 
   CALL SEARCH('TRNY',4,11,IERR)
   READ(11,*) NOC
   DO I=1,NOC
      READ(11,*)
      ENDDO
      DO J=0,M%JBAR
      READ(11,*) IDUM,M%Y(J)
   ENDDO
! 
   CALL SEARCH('TRNZ',4,11,IERR)
   READ(11,*) NOC
   DO I=1,NOC
      READ(11,*)
      ENDDO
      DO K=0,M%KBAR
      READ(11,*) IDUM,M%Z(K)
   ENDDO
! 
ENDDO READ_SMV
! 
IFILE=3 
NSAM=1 
ANS='y'
CALL TOUPPER(ANS,ANS)
! **************************************************
BNDF_MESH = 1
REWIND(11)
NFILES_EXIST=0
!
SEARCH_BNDF: DO I=1,500
   CALL SEARCH2('BNDF',4,'BNDC',4,11,IERR,CHOICE)
   IF (IERR==1) EXIT SEARCH_BNDF
      BACKSPACE(11)
      READ(11,*) JUNK,BNDF_MESH(I)
      READ(11,'(A)') BNDF_FILE(I)
      READ(11,'(A)') BNDF_TEXT(I)
      READ(11,*)
      READ(11,'(A)') BNDF_UNIT(I)
      OPEN(12,FILE=BNDF_FILE(I),FORM='UNFORMATTED',STATUS='OLD', IOSTAT=RCODE)
      CLOSE(12)
   IF (RCODE/=0) CYCLE
      NFILES_EXIST=NFILES_EXIST+1
      IF (CHOICE=='BNDC') BNDF_TYPE(I) = 'CENTERED'
      IF (CHOICE=='BNDF') BNDF_TYPE(I) = 'STAGGERED'
     !IF (BATCHMODE==0) write(6,'(I3,3X,A,I2,A,A)') I,'MESH ',BNDF_MESH(I), ', ',TRIM(BNDF_TEXT(I))
ENDDO SEARCH_BNDF
!
IF (NFILES_EXIST==0)THEN
   IF(BATCHMODE==0) write(6,*)"There are no boundary files to convert"
   STOP
ENDIF
!
BEG_VAR=1
IF (VARIAB==3) THEN
   DO I=1,500
      IF (TRIM(BNDF_TEXT(I))==' NET HEAT FLUX                ') THEN
         BEG_VAR=I
         VARIAB=I
         EXIT
      ENDIF
   ENDDO
ENDIF
!
LOOP_FDS2AST:DO P=1,COUNTER
   Xa=NOMASTER(P,2)
   Ya=NOMASTER(P,3)
   Za=NOMASTER(P,4)
   WRITE (OUTFILE,'(g8.0)') INT(NOMASTER(P,1))
   PRINT *, OUTFILE
   C_SIZE=C_S
   NX=N(P,2)
   NY=N(P,3)
   NZ=N(P,4)
   IF (ABS(NX)<(MAX(ABS(NY),ABS(NZ))/50)) NX=0
   IF (ABS(NY)<(MAX(ABS(NX),ABS(NZ))/50)) NY=0
   IF (ABS(NZ)<(MAX(ABS(NY),ABS(NX))/50)) NZ=0
   N_DIR=0.D0
   IF (NX/=0) N_DIR=N_DIR+1
   IF (NY/=0) N_DIR=N_DIR+1
   IF (NZ/=0) N_DIR=N_DIR+1
   !
   XS = Xa-(C_SIZE/2)
   XF = Xa+(C_SIZE/2)
   YS = Ya-(C_SIZE/2)
   YF = Ya+(C_SIZE/2)
   ZS = Za-(C_SIZE/2)
   ZF = Za+(C_SIZE/2)
   NUM_INT_CSIZE=1
   VARIABLE_NUMBER: DO VAR=BEG_VAR,VARIAB
!     500 - Point of the loop in case of the point does not reach a meaning value - SILVA, JC         
   500 CONTINUE   
   COMEBACK=0
!     SIZE is the lenght of the vectors - SILVA, JC
   SIZE=((TEND-TBEG)/TINT)+1
   ALLOCATE (V_TIME(SIZE))
   ALLOCATE (TAST(SIZE))
   ALLOCATE (M_AST(SIZE,7))
   IF (TRIM(BNDF_TEXT(VAR))==' ADIABATIC SURFACE TEMPERATURE') VCOUNT=0.0
   IF (TRIM(BNDF_TEXT(VAR))==' NET HEAT FLUX                ') VCOUNT=0.0
   !VCOUNT=0.0
   V_TIME=0.0
   TAST=0.0
   M_AST=0.0
!     IN_LOOP is a loop at index file - SILVA, JC
   IN_LOOP: DO IN=VAR,VARIAB*NMESHES,VARIAB
!*********
      IOR_LOOPING: DO IOR_LOOP=1,7
         IOR_INPUT=IOR_LOOP-4
         IF (IOR_INPUT==0) CYCLE
         IF (IOR_INPUT==-3 .AND. NZ>=0) CYCLE
         IF (IOR_INPUT==-2 .AND. NY>=0) CYCLE
         IF (IOR_INPUT==-1 .AND. NX>=0) CYCLE
         IF (IOR_INPUT==1 .AND. NX<=0) CYCLE
         IF (IOR_INPUT==2 .AND. NY<=0) CYCLE
         IF (IOR_INPUT==3 .AND. NZ<=0) CYCLE
         NV=1
         MV=1
         IB(MV) = IN
         QFILE = BNDF_FILE(IN)
         OPEN(12+VAR,FILE=QFILE,FORM='UNFORMATTED',STATUS='OLD', IOSTAT=RCODE)
         IF (RCODE/=0) THEN
            CLOSE(12+VAR)
            CYCLE
         ENDIF
!
      NM = BNDF_MESH(IN)
      M=>MESH(NM)
      ALLOCATE(Q(0:M%IBAR,0:M%JBAR,0:M%KBAR,NV))
      Q = 0.0
      ALLOCATE(ALREADY_USED(0:M%IBAR,0:M%JBAR,0:M%KBAR))
      BNDF_TYPE_CHOSEN = BNDF_TYPE(IN)
!
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
!    
      DO I=1,NPATCH
         READ(12+VAR) I1B(I),I2B(I),J1B(I),J2B(I),K1B(I),K2B(I),IOR(I)
      ENDDO
!    
      IJBAR = MAX(M%IBAR,M%JBAR)
      JKBAR = MAX(M%JBAR,M%KBAR)
      ALLOCATE(F(0:IJBAR,0:JKBAR,NPATCH))
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
         IF (TIME<TBEG) CYCLE READ_BLOOP
         IF (INT(TIME)>TEND) EXIT READ_BLOOP
                    
         ALREADY_USED = .FALSE.
            REC_PATCH: DO II=1,NPATCH
               IF (IOR(II)/=IOR_INPUT) CYCLE REC_PATCH
               IF (M%X(I1B(II))>XF .OR. M%X(I2B(II))<XS) CYCLE REC_PATCH
               IF (M%Y(J1B(II))>YF .OR. M%Y(J2B(II))<YS) CYCLE REC_PATCH
               IF (M%Z(K1B(II))>ZF .OR. M%Z(K2B(II))<ZS) CYCLE REC_PATCH
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
!                               Q(I,J1B(II),K,MV) = Q(I,J1B(II),K,MV) + F(I,K,II)
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
!
            V=(TIME/TINT)+1
            V_TIME(V)=TIME
               PATCHES: DO II=1,NPATCH
                  DO K=K1B(II),K2B(II),NSAM
                     DO J=J1B(II),J2B(II),NSAM
                        DO I=I1B(II),I2B(II),NSAM
                           IF (M%X(I)>XF .OR. M%X(I)<XS) CYCLE 
                           IF (M%Y(J)>YF .OR. M%Y(J)<YS) CYCLE
                           IF (M%Z(K)>ZF .OR. M%Z(K)<ZS) CYCLE
                           M_AST(V,IOR_LOOP)=M_AST(V,IOR_LOOP)+Q(I,J,K,NV)
                           !IF (TRIM(BNDF_TEXT(VAR))==' HEAT TRANSFER COEFFICIENT    ') THEN
                           !   IF (TIME==TBEG .AND. Q(I,J,K,NV)==0) VCOUNT(IOR_LOOP)=VCOUNT(IOR_LOOP)+1
                           !ENDIF
                           IF (TRIM(BNDF_TEXT(VAR))==' ADIABATIC SURFACE TEMPERATURE') THEN
                              IF (TIME==TBEG .AND. Q(I,J,K,NV)>0) VCOUNT(IOR_LOOP)=VCOUNT(IOR_LOOP)+1
                           ENDIF
                           IF (TRIM(BNDF_TEXT(VAR))==' NET HEAT FLUX                ') THEN   
                              IF (TIME==TBEG .AND. Q(I,J,K,NV)/=0) VCOUNT(IOR_LOOP)=VCOUNT(IOR_LOOP)+1
                           ENDIF
                        ENDDO 
                     ENDDO
                  ENDDO
               ENDDO PATCHES
!   
      ENDDO READ_BLOOP
!      
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
!
   IF (VCOUNT(1)>1) THEN
      DO I=1,V
         M_AST(I,1)=M_AST(I,1)/VCOUNT(1)
      END DO
   END IF
   IF (VCOUNT(2)>1) THEN
      DO I=1,V
         M_AST(I,2)=M_AST(I,2)/VCOUNT(2)
      END DO
   END IF
   IF (VCOUNT(3)>1) THEN
      DO I=1,V
      M_AST(I,3)=M_AST(I,3)/VCOUNT(3)
      END DO
   END IF
   IF (VCOUNT(5)>1) THEN
      DO I=1,V
         M_AST(I,5)=M_AST(I,5)/VCOUNT(5)
      END DO
   END IF
   IF (VCOUNT(6)>1) THEN
      DO I=1,V
         M_AST(I,6)=M_AST(I,6)/VCOUNT(6)
      END DO
   END IF
   IF (VCOUNT(7)>1) THEN
      DO I=1,V
         M_AST(I,7)=M_AST(I,7)/VCOUNT(7)
      END DO
   END IF
!
   VARIABLE_KIND=TRIM(BNDF_TEXT(VAR))
   IF (VARIABLE_KIND==' ADIABATIC SURFACE TEMPERATURE') THEN
!
      IF (NZ>0) THEN
         DO I=1,V
            M_AST(I,1)=0.D0
         END DO
      END IF
      IF (NZ<0) THEN 
         DO I=1,V
            M_AST(I,7)=0.D0
         END DO
      END IF
      IF (NZ==0) THEN 
          DO I=1,V
             M_AST(I,1)=0.D0
             M_AST(I,7)=0.D0
          END DO
      END IF
!        
      IF (NY>0) THEN
         DO I=1,V
            M_AST(I,2)=0.D0
         END DO
      END IF
      IF (NY<0) THEN 
         DO I=1,V
            M_AST(I,6)=0.D0
         END DO
      END IF
      IF (NY==0) THEN 
         DO I=1,V
            M_AST(I,2)=0.D0
            M_AST(I,6)=0.D0
         END DO
      END IF
!
      IF (NX>0) THEN
         DO I=1,V
            M_AST(I,3)=0.D0
         END DO
      END IF
      IF (NX<0) THEN 
         DO I=1,V
            M_AST(I,5)=0.D0
         END DO
      END IF
      IF (NX==0) THEN 
         DO I=1,V
            M_AST(I,3)=0.D0
            M_AST(I,5)=0.D0
         END DO
      END IF
!
      IF (M_AST(2,7)<19 .AND. M_AST(2,7)>0) THEN
         DO I=1,V
            M_AST(I,7)=20.0
         END DO
      END IF 
      IF (M_AST(2,1)<19 .AND. M_AST(2,1)>0) THEN
         DO I=1,V
            M_AST(I,1)=20.0
         END DO
        END IF 
      IF (M_AST(2,2)<19 .AND. M_AST(2,2)>0) THEN
         DO I=1,V
            M_AST(I,2)=20.0
        END DO
      END IF 
      IF (M_AST(2,6)<19 .AND. M_AST(2,6)>0) THEN
         DO I=1,V
            M_AST(I,6)=20.0
         END DO
      END IF 
      IF (M_AST(2,5)<19 .AND. M_AST(2,5)>0) THEN
         DO I=1,V
            M_AST(I,5)=20.0
         END DO
      END IF 
      IF (M_AST(2,3)<19 .AND. M_AST(2,3)>0) THEN
         DO I=1,V
            M_AST(I,3)=20.0
         END DO
      END IF 
      IF (M_AST(1,7)<20 .AND. M_AST(1,7)>0) THEN
         M_AST(1,7)=20.0
      END IF 
      IF (M_AST(1,1)<20 .AND. M_AST(1,1)>0) THEN
         M_AST(1,1)=20.0
      END IF 
      IF (M_AST(1,2)<20 .AND. M_AST(1,2)>0) THEN
         M_AST(1,2)=20.0
      END IF 
      IF (M_AST(1,6)<20 .AND. M_AST(1,6)>0) THEN
         M_AST(1,6)=20.0
      END IF 
      IF (M_AST(1,5)<20 .AND. M_AST(1,5)>0) THEN
         M_AST(1,5)=20.0
      END IF 
      IF (M_AST(1,3)<20 .AND. M_AST(1,3)>0) THEN
         M_AST(1,3)=20.0
      END IF 
!   
      IF (NZ<0) THEN
         IF (M_AST(1,1)==0) THEN
            ZS = ZS-(0.05*(ZF-ZS))
            ZF = ZF+(0.05*(ZF-ZS))
            COMEBACK=1
            NUM_INT_CSIZE(1)=NUM_INT_CSIZE(1)+1
         ENDIF
      END IF
!
      IF (NY<0) THEN
         IF (M_AST(1,2)==0) THEN
            YS = YS-(0.05*(YF-YS))
            YF = YF+(0.05*(YF-YS))
            COMEBACK=1
            NUM_INT_CSIZE(2)=NUM_INT_CSIZE(2)+1
         ENDIF
      END IF
! 
      IF (NX<0) THEN
         IF (M_AST(1,3)==0) THEN
            XS = XS-(0.05*(XF-XS))
            XF = XF+(0.05*(XF-XS))      
            COMEBACK=1    
            NUM_INT_CSIZE(3)=NUM_INT_CSIZE(3)+1  
         ENDIF
      END IF 
!
      IF (NZ>0) THEN
         IF (M_AST(1,7)==0) THEN
            ZS = ZS-(0.05*(ZF-ZS))
            ZF = ZF+(0.05*(ZF-ZS))    
            COMEBACK=1       
            NUM_INT_CSIZE(7)=NUM_INT_CSIZE(7)+1 
         ENDIF
      END IF
!
      IF (NY>0) THEN
         IF (M_AST(1,6)==0) THEN
            YS = YS-(0.05*(YF-YS))
            YF = YF+(0.05*(YF-YS))   
            COMEBACK=1         
            NUM_INT_CSIZE(6)=NUM_INT_CSIZE(6)+1
         ENDIF
      END IF
! 
      IF (NX>0) THEN
         IF (M_AST(1,5)==0) THEN
            XS = XS-(0.05*(XF-XS))
            XF = XF+(0.05*(XF-XS))  
            COMEBACK=1          
            NUM_INT_CSIZE(5)=NUM_INT_CSIZE(5)+1
         ENDIF
      END IF
!
      IF (MAX(NUM_INT_CSIZE(1),NUM_INT_CSIZE(2),NUM_INT_CSIZE(3),NUM_INT_CSIZE(5),NUM_INT_CSIZE(6),NUM_INT_CSIZE(7))>=5) THEN
         XS = XS-(0.05*(XF-XS))
         XF = XF+(0.05*(XF-XS)) 
         YS = YS-(0.05*(YF-YS))
         YF = YF+(0.05*(YF-YS))
         ZS = ZS-(0.05*(ZF-ZS))
         ZF = ZF+(0.05*(ZF-ZS))  
      ENDIF                 
!                    
      IF (COMEBACK==1) THEN
         DEALLOCATE (V_TIME)
         DEALLOCATE (TAST)
         DEALLOCATE (M_AST)
         GOTO 500      
      ENDIF     
! 
      C=SQRT((NX**2)+(NY**2)+(NZ**2))
!
      IF (NZ<0 .AND. NY<0 .AND. NX<0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY<0 .AND. NX>0) THEN              
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY>0 .AND. NX<0) THEN              
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY>0 .AND. NX>0) THEN              
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY<0 .AND. NX==0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY==0 .AND. NX<0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ<0 .AND. NY==0 .AND. NX>0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ<0 .AND. NY>0 .AND. NX==0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ<0 .AND. NY==0 .AND. NX==0) THEN 
         DO I=1,V
            A=((-M_AST(I,1)*NZ))
            B=SQRT((M_AST(I,1)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY<0 .AND. NX<0) THEN         
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY<0 .AND. NX>0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY>0 .AND. NX<0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY>0 .AND. NX>0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY<0 .AND. NX==0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY==0 .AND. NX<0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ>0 .AND. NY==0 .AND. NX>0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ>0 .AND. NY>0 .AND. NX==0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ>0 .AND. NY==0 .AND. NX==0) THEN 
         DO I=1,V
            A=((M_AST(I,7)*NZ))
            B=SQRT((M_AST(I,7)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY<0 .AND. NX<0) THEN 
         DO I=1,V
            A=((-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY<0 .AND. NX>0) THEN
         DO I=1,V
            A=((-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY<0 .AND. NX==0) THEN
         DO I=1,V
            A=((-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY>0 .AND. NX<0) THEN
         DO I=1,V
            A=((M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY>0 .AND. NX>0) THEN
         DO I=1,V
            A=((M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY>0 .AND. NX==0) THEN
         DO I=1,V
            A=((M_AST(I,6)*NY))
            B=SQRT((M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY==0 .AND. NX<0) THEN
         DO I=1,V
            A=((-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY==0 .AND. NX>0) THEN
         DO I=1,V
            A=((M_AST(I,5)*NX))
            B=SQRT((M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
!********************************
      IF (TAST(1)>0) THEN
         TAST(1)=20
      ELSE
         XS = XS-(0.05*(XF-XS))
         XF = XF+(0.05*(XF-XS)) 
         YS = YS-(0.05*(YF-YS))
         YF = YF+(0.05*(YF-YS)) 
         ZS = ZS-(0.05*(ZF-ZS))
         ZF = ZF+(0.05*(ZF-ZS))    
         DEALLOCATE (V_TIME)
         DEALLOCATE (TAST)
         DEALLOCATE (M_AST)
         GOTO 500
      END IF  
   !
   END IF
!***************************************
   IF (VARIABLE_KIND==' WALL TEMPERATURE             ') THEN
      POND=0
      IF (M_AST(2,1)>20) POND=POND+1
      IF (M_AST(2,2)>20) POND=POND+1
      IF (M_AST(2,3)>20) POND=POND+1
      IF (M_AST(2,5)>20) POND=POND+1
      IF (M_AST(2,6)>20) POND=POND+1
      IF (M_AST(2,7)>20) POND=POND+1
      IF (POND/=0) THEN
         DO I=1,V
            TAST(I)=(M_AST(I,1)+M_AST(I,2)+M_AST(I,3)+M_AST(I,4)+M_AST(I,5)+M_AST(I,6))/POND
         ENDDO
      END IF
   END IF
!************************************************************
   IF (VARIABLE_KIND==' HEAT TRANSFER COEFFICIENT    ') THEN     
      IF (NZ>0) THEN
         DO I=1,V
            M_AST(I,1)=0.D0
         END DO
      END IF
      IF (NZ<0) THEN 
         DO I=1,V
            M_AST(I,7)=0.D0
         END DO
      END IF
      IF (NZ==0) THEN 
         DO I=1,V
            M_AST(I,1)=0.D0
            M_AST(I,7)=0.D0
         END DO
      END IF
!      
      IF (NY>0) THEN
         DO I=1,V
            M_AST(I,2)=0.D0
         END DO
      END IF
      IF (NY<0) THEN 
         DO I=1,V
            M_AST(I,6)=0.D0
         END DO
      END IF
      IF (NY==0) THEN 
         DO I=1,V
            M_AST(I,2)=0.D0
            M_AST(I,6)=0.D0
         END DO
      END IF
!
      IF (NX>0) THEN
         DO I=1,V
            M_AST(I,3)=0.D0
         END DO
      END IF
      IF (NX<0) THEN 
         DO I=1,V
            M_AST(I,5)=0.D0
         END DO
      END IF
      IF (NX==0) THEN 
         DO I=1,V
            M_AST(I,3)=0.D0
            M_AST(I,5)=0.D0
         END DO
      END IF
!
      H_NULL=0.d0
      IF (NZ<0) THEN
         SUM=0.d0
         DO I=1,V
            SUM=SUM+ABS(M_AST(I,1))
         ENDDO
         IF (SUM==0) H_NULL(3)=1
      END IF
!
      IF (NY<0) THEN
         SUM=0.d0
         DO I=1,V
            SUM=SUM+ABS(M_AST(I,2))
         ENDDO
         IF (SUM==0) H_NULL(2)=1
      END IF
! 
      IF (NX<0) THEN
         SUM=0.d0
         DO I=1,V
            SUM=SUM+ABS(M_AST(I,3))
         ENDDO
         IF (SUM==0) H_NULL(1)=1
      END IF 
!
      IF (NZ>0) THEN
         SUM=0.d0
         DO I=1,V
            SUM=SUM+ABS(M_AST(I,7))
         ENDDO
         IF (SUM==0) H_NULL(3)=1
      END IF
!
      IF (NY>0) THEN
         SUM=0.d0
         DO I=1,V
            SUM=SUM+ABS(M_AST(I,6))
         ENDDO
         IF (SUM==0) H_NULL(2)=1
      END IF
! 
      IF (NX>0) THEN
         SUM=0.d0
         DO I=1,V
            SUM=SUM+ABS(M_AST(I,5))
         ENDDO
         IF (SUM==0) H_NULL(1)=1
      END IF
!
      IF (NZ/=0 .AND. H_NULL(3)==1) COMEBACK=COMEBACK+1
      IF (NY/=0 .AND. H_NULL(2)==1) COMEBACK=COMEBACK+1
      IF (NX/=0 .AND. H_NULL(1)==1) COMEBACK=COMEBACK+1
      IF (COMEBACK==N_DIR) THEN
         DO I=1,V
            TAST(I)=0.d0
         END DO
         GOTO 600      
      ENDIF          
!
      C=SQRT((NX**2)+(NY**2)+(NZ**2))
!
      IF (NZ<0 .AND. NY<0 .AND. NX<0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY<0 .AND. NX>0) THEN              
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY>0 .AND. NX<0) THEN              
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY>0 .AND. NX>0) THEN              
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY<0 .AND. NX==0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY==0 .AND. NX<0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ<0 .AND. NY==0 .AND. NX>0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ<0 .AND. NY>0 .AND. NX==0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ<0 .AND. NY==0 .AND. NX==0) THEN 
         DO I=1,V
            A=((-M_AST(I,1)*NZ))
            B=SQRT((M_AST(I,1)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY<0 .AND. NX<0) THEN         
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY<0 .AND. NX>0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY>0 .AND. NX<0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY>0 .AND. NX>0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY<0 .AND. NX==0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY==0 .AND. NX<0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ>0 .AND. NY==0 .AND. NX>0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ>0 .AND. NY>0 .AND. NX==0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ>0 .AND. NY==0 .AND. NX==0) THEN 
         DO I=1,V
            A=((M_AST(I,7)*NZ))
            B=SQRT((M_AST(I,7)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY<0 .AND. NX<0) THEN 
         DO I=1,V
            A=((-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY<0 .AND. NX>0) THEN
         DO I=1,V
            A=((-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY<0 .AND. NX==0) THEN
         DO I=1,V
            A=((-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY>0 .AND. NX<0) THEN
         DO I=1,V
            A=((M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY>0 .AND. NX>0) THEN
         DO I=1,V
            A=((M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY>0 .AND. NX==0) THEN
         DO I=1,V
            A=((M_AST(I,6)*NY))
            B=SQRT((M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY==0 .AND. NX<0) THEN
         DO I=1,V
            A=((-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY==0 .AND. NX>0) THEN
         DO I=1,V
            A=((M_AST(I,5)*NX))
            B=SQRT((M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
!********************************
      IF (TAST(2)>=0 .AND. TAST(2)<1000) THEN
         TAST(1)=0.0
!      ELSE 
!         XS = XS-(0.05*(XF-XS))
!         XF = XF+(0.05*(XF-XS)) 
!         YS = YS-(0.05*(YF-YS))
!         YF = YF+(0.05*(YF-YS)) 
!         ZS = ZS-(0.05*(ZF-ZS))
!         ZF = ZF+(0.05*(ZF-ZS))     
!         DEALLOCATE (V_TIME)
!         DEALLOCATE (TAST)
!         DEALLOCATE (M_AST)
!         GOTO 500  
      END IF
!   
   END IF
!***************************************************************
   IF (VARIABLE_KIND==' NET HEAT FLUX                ') THEN
      IF (NZ>0) THEN
         DO I=1,V
            M_AST(I,1)=0.D0
         END DO
      END IF
      IF (NZ<0) THEN 
         DO I=1,V
            M_AST(I,7)=0.D0
         END DO
      END IF
      IF (NZ==0) THEN 
         DO I=1,V
            M_AST(I,1)=0.D0
            M_AST(I,7)=0.D0
         END DO
      END IF
!        
      IF (NY>0) THEN
         DO I=1,V
            M_AST(I,2)=0.D0
         END DO
      END IF
      IF (NY<0) THEN 
         DO I=1,V
            M_AST(I,6)=0.D0
         END DO
      END IF
      IF (NY==0) THEN 
         DO I=1,V
            M_AST(I,2)=0.D0
            M_AST(I,6)=0.D0
         END DO
      END IF
!
      IF (NX>0) THEN
         DO I=1,V
            M_AST(I,3)=0.D0
         END DO
      END IF
      IF (NX<0) THEN 
         DO I=1,V
            M_AST(I,5)=0.D0
         END DO
      END IF
      IF (NX==0) THEN 
         DO I=1,V
            M_AST(I,3)=0.D0
            M_AST(I,5)=0.D0
         END DO
      END IF
!
      IF (NZ<0) THEN
         SUM=0.d0
         DO I=1,V
            SUM=SUM+ABS(M_AST(I,1))
         ENDDO
         IF (SUM==0) THEN
            ZS = ZS-(0.05*(ZF-ZS))
            ZF = ZF+(0.05*(ZF-ZS))            
            COMEBACK=1
         ENDIF
      END IF
!
      IF (NY<0) THEN
         SUM=0.d0
         DO I=1,V
            SUM=SUM+ABS(M_AST(I,2))
         ENDDO
         IF (SUM==0) THEN
            YS = YS-(0.05*(YF-YS))
            YF = YF+(0.05*(YF-YS))            
            COMEBACK=1
         ENDIF
      END IF
! 
      IF (NX<0) THEN
         SUM=0.d0
         DO I=1,V
            SUM=SUM+ABS(M_AST(I,3))
         ENDDO
         IF (SUM==0) THEN
            XS = XS-(0.05*(XF-XS))
            XF = XF+(0.05*(XF-XS))            
            COMEBACK=1
         ENDIF
      END IF 
!
      IF (NZ>0) THEN
         SUM=0.d0
         DO I=1,V
            SUM=SUM+ABS(M_AST(I,7))
         ENDDO
         IF (SUM==0) THEN
            ZS = ZS-(0.05*(ZF-ZS))
            ZF = ZF+(0.05*(ZF-ZS))            
            COMEBACK=1
         ENDIF
      END IF
!
      IF (NY>0) THEN
         SUM=0.d0
         DO I=1,V
            SUM=SUM+ABS(M_AST(I,6))
         ENDDO
         IF (SUM==0) THEN
            YS = YS-(0.05*(YF-YS))
            YF = YF+(0.05*(YF-YS))            
            COMEBACK=1
         ENDIF
      END IF
! 
      IF (NX>0) THEN
         SUM=0.d0
         DO I=1,V
            SUM=SUM+ABS(M_AST(I,5))
         ENDDO
         IF (SUM==0) THEN
            XS = XS-(0.05*(XF-XS))
            XF = XF+(0.05*(XF-XS))            
            COMEBACK=1
         ENDIF
      END IF
      IF (COMEBACK==1) THEN
         DEALLOCATE (V_TIME)
         DEALLOCATE (TAST)
         DEALLOCATE (M_AST)
         GOTO 500      
      ENDIF
!        
      C=SQRT((NX**2)+(NY**2)+(NZ**2))
!
      IF (NZ<0 .AND. NY<0 .AND. NX<0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY<0 .AND. NX>0) THEN              
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY>0 .AND. NX<0) THEN              
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY>0 .AND. NX>0) THEN              
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY<0 .AND. NX==0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ<0 .AND. NY==0 .AND. NX<0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ<0 .AND. NY==0 .AND. NX>0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ<0 .AND. NY>0 .AND. NX==0) THEN      
         DO I=1,V
            A=((-M_AST(I,1)*NZ)+(M_AST(I,6)*NY))
            B=SQRT((M_AST(I,1)**2)+(M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ<0 .AND. NY==0 .AND. NX==0) THEN 
         DO I=1,V
            A=((-M_AST(I,1)*NZ))
            B=SQRT((M_AST(I,1)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY<0 .AND. NX<0) THEN         
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY<0 .AND. NX>0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY>0 .AND. NX<0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY>0 .AND. NX>0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY<0 .AND. NX==0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ>0 .AND. NY==0 .AND. NX<0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ>0 .AND. NY==0 .AND. NX>0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ>0 .AND. NY>0 .AND. NX==0) THEN  
         DO I=1,V
            A=((M_AST(I,7)*NZ)+(M_AST(I,6)*NY))
            B=SQRT((M_AST(I,7)**2)+(M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF  
      IF (NZ>0 .AND. NY==0 .AND. NX==0) THEN 
         DO I=1,V
            A=((M_AST(I,7)*NZ))
            B=SQRT((M_AST(I,7)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY<0 .AND. NX<0) THEN 
         DO I=1,V
            A=((-M_AST(I,2)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,2)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY<0 .AND. NX>0) THEN
         DO I=1,V
            A=((-M_AST(I,2)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,2)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY<0 .AND. NX==0) THEN
         DO I=1,V
            A=((-M_AST(I,2)*NY))
            B=SQRT((M_AST(I,2)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY>0 .AND. NX<0) THEN
         DO I=1,V
            A=((M_AST(I,6)*NY)+(-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,6)**2)+(M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY>0 .AND. NX>0) THEN
         DO I=1,V
            A=((M_AST(I,6)*NY)+(M_AST(I,5)*NX))
            B=SQRT((M_AST(I,6)**2)+(M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY>0 .AND. NX==0) THEN
         DO I=1,V
            A=((M_AST(I,6)*NY))
            B=SQRT((M_AST(I,6)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY==0 .AND. NX<0) THEN
         DO I=1,V
            A=((-M_AST(I,3)*NX))
            B=SQRT((M_AST(I,3)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
      IF (NZ==0 .AND. NY==0 .AND. NX>0) THEN
         DO I=1,V
            A=((M_AST(I,5)*NX))
            B=SQRT((M_AST(I,5)**2))
            TAST(I)=(A/(B*C))*B
         END DO
      END IF
!********************************
      SUM=0.d0
      DO I=1,V
         SUM=SUM+ABS(TAST(I))
      ENDDO
      IF (SUM>0) THEN
         TAST(1)=0.0
      ELSE
         XS = XS-(0.05*(XF-XS))
         XF = XF+(0.05*(XF-XS)) 
         YS = YS-(0.05*(YF-YS))
         YF = YF+(0.05*(YF-YS)) 
         ZS = ZS-(0.05*(ZF-ZS))
         ZF = ZF+(0.05*(ZF-ZS))   
         DEALLOCATE (V_TIME)
         DEALLOCATE (TAST)
         DEALLOCATE (M_AST)   
         GOTO 500                 
     END IF 
!
   END IF
!
!*** Set an time averaged function to represent the evaluation of the variables (T_AVERAGE)
   600 CONTINUE
   IF (T_AVERAGE==0) THEN
      DO I=1,V
         IF (VARIABLE_KIND==' ADIABATIC SURFACE TEMPERATURE') THEN
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", I, ",0),", V_TIME(I)
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", I, ",1),", TAST(I)
         ENDIF
         IF (VARIABLE_KIND==' HEAT TRANSFER COEFFICIENT    ') THEN
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", I, ",0),", V_TIME(I)
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", I, ",1),", TAST(I)
         ENDIF
         IF (VARIABLE_KIND==' NET HEAT FLUX                ') THEN
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", I, ",0),", V_TIME(I)
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", I, ",1),", TAST(I)*1000
         ENDIF
      ENDDO
   ELSE
      I_AVERAGE=T_AVERAGE/TINT
      ALLOCATE (MED_TAST((V-1)/I_AVERAGE))
      DO I=1,(V-1)/I_AVERAGE
         MED=0.0
         DO J=I_AVERAGE*(I-1),(I_AVERAGE*I)-1
            MED=MED+TAST(J+2)
         ENDDO
         MED=MED/I_AVERAGE
         MED_TAST(I)=MED
      ENDDO
      IF (VARIABLE_KIND==' ADIABATIC SURFACE TEMPERATURE') THEN
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", 1, ",0),", V_TIME(1)
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", 1, ",1),", 20.0
         DO I=2,((V-1)/I_AVERAGE)
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", I, ",0),", V_TIME((I-1)*I_AVERAGE+1)
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", I, ",1),", (MED_TAST(I-1)+MED_TAST(I))/2
         ENDDO
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+1, ",0),", V_TIME(V)
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+1, ",1),", MED_TAST((V-1)/I_AVERAGE)+ &
                                                                    (MED_TAST((V-1)/I_AVERAGE)-MED_TAST((V-1)/I_AVERAGE-1))/2
      ENDIF
      IF (VARIABLE_KIND==' HEAT TRANSFER COEFFICIENT    ') THEN
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", 1, ",0),", V_TIME(1)
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", 1, ",1),", 0.0
         DO I=2,((V-1)/I_AVERAGE)
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", I, ",0),", V_TIME((I-1)*I_AVERAGE+1)
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", I, ",1),", (MED_TAST(I-1)+MED_TAST(I))/2
         ENDDO
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+1, ",0),", V_TIME(V)
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+1, ",1),", MED_TAST((V-1)/I_AVERAGE)+ & 
                                                                    (MED_TAST((V-1)/I_AVERAGE)-MED_TAST((V-1)/I_AVERAGE-1))/2
      ENDIF
      IF (VARIABLE_KIND==' NET HEAT FLUX                ') THEN
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", 1, ",0),", V_TIME(1)
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", 1, ",1),", 0.0
         DO I=2,((V-1)/I_AVERAGE)
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", I, ",0),", V_TIME((I-1)*I_AVERAGE+1)
            WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", I, ",1),", ((MED_TAST(I-1)+MED_TAST(I))/2)*1000
         ENDDO
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+1, ",0),", V_TIME(V)
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+1, ",1),", (MED_TAST((V-1)/I_AVERAGE)+ &
                                                                 (MED_TAST((V-1)/I_AVERAGE)-MED_TAST((V-1)/I_AVERAGE-1))/2)*1000
      ENDIF 
   ENDIF 
!**********************
!*** Set an average value for a steady simulation (last N_AVERAGE results)
   IF (N_AVERAGE==0) THEN
      CONTINUE
   ELSE IF (T_AVERAGE==0) THEN
      MED=0.0
      DO I=0,N_AVERAGE-1
         MED=MED+TAST(V-I)
      ENDDO
      MED=MED/(N_AVERAGE)
      IF (VARIABLE_KIND==' ADIABATIC SURFACE TEMPERATURE') THEN    
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", V, ",0),", V_TIME(V)+TINT
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", V, ",1),", MED
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", V+1, ",0),", 18000.0
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", V+1, ",1),", MED
      ENDIF
      IF (VARIABLE_KIND==' HEAT TRANSFER COEFFICIENT    ') THEN
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", V, ",0),", V_TIME(V)+TINT
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", V, ",1),", MED
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", V+1, ",0),", 18000.0
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", V+1, ",1),", MED
      ENDIF
      IF (VARIABLE_KIND==' NET HEAT FLUX                ') THEN    
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", V, ",0),", V_TIME(V)+TINT
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", V, ",1),", MED*1000
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", V+1, ",0),", 18000.0
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", V+1, ",1),", MED*1000
      ENDIF
   ELSE IF (T_AVERAGE/=0) THEN
      MED=0.0
      DO I=0,N_AVERAGE-1
         MED=MED+MED_TAST(((V-1)/I_AVERAGE)-I)
      ENDDO
      MED=MED/(N_AVERAGE)
      IF (VARIABLE_KIND==' ADIABATIC SURFACE TEMPERATURE') THEN    
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+2, ",0),", V_TIME(V)+T_AVERAGE
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+2, ",1),", MED
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+2+1, ",0),", 18000.0
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+2+1, ",1),", MED
      ENDIF
      IF (VARIABLE_KIND==' HEAT TRANSFER COEFFICIENT    ') THEN
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+2, ",0),", V_TIME(V)+T_AVERAGE
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+2, ",1),", MED
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+2+1, ",0),", 18000.0
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,H", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+2+1, ",1),", MED
      ENDIF
      IF (VARIABLE_KIND==' NET HEAT FLUX                ') THEN    
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+2, ",0),", V_TIME(V)+T_AVERAGE
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+2, ",1),", MED*1000
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+2+1, ",0),", 18000.0
         WRITE (70,'(A,A,A,I8,A,E12.5)') "*set,A", TRIM(OUTFILE), "(", (V-1)/I_AVERAGE+2+1, ",1),", MED*1000
      ENDIF    
   ENDIF    
!**********************
   IF (VARIABLE_KIND==' ADIABATIC SURFACE TEMPERATURE') THEN
      WRITE (6,'(A, F8.3, A, F8.3, A, F8.3)') 'Interface point    ', Xa, '        ', Ya, &
      '        ', Za
      WRITE (6,'(A, F8.3, A, F8.3, A, F8.3)') 'Cell dimension     ', (XF-XS), '        ', (YF-YS), &
      '        ', (ZF-ZS)
      WRITE (6,'(A, F8.3, F8.3, F8.3, F8.3, F8.3, F8.3)') 'AST Results    ', M_AST(5,3), M_AST(5,5), M_AST(5,2), &
      M_AST(5,6), M_AST(5,1), M_AST(5,7)
      IF (VARIABLE==1) WRITE (6,'(A, F8.5, A, F8.5, A, F8.5)') 'Normal_Vector       ',NX, '        ', NY, '        ', NZ
   END IF
   IF (VARIABLE_KIND==' HEAT TRANSFER COEFFICIENT    ') THEN
      WRITE (6,'(A, F8.3, A, F8.3, A, F8.3)') 'Cell dimension     ', (XF-XS), '        ', (YF-YS), &
      '        ', (ZF-ZS)
      WRITE (6,'(A, F8.3, F8.3, F8.3, F8.3, F8.3, F8.3)') ' h  Results    ', M_AST(5,3), M_AST(5,5), M_AST(5,2), &
      M_AST(5,6), M_AST(5,1), M_AST(5,7)
      IF (VARIABLE==2) WRITE (6,'(A, F8.5, A, F8.5, A, F8.5)') 'Normal_Vector       ',NX, '        ', NY, '        ', NZ
   END IF    
   IF (VARIABLE_KIND==' NET HEAT FLUX                ') THEN
      WRITE (6,'(A, F8.3, A, F8.3, A, F8.3)') 'Interface point    ', Xa, '        ', Ya, &
      '        ', Za
      WRITE (6,'(A, F8.3, A, F8.3, A, F8.3)') 'Cell dimension     ', (XF-XS), '        ', (YF-YS), &
      '        ', (ZF-ZS)
      WRITE (6,'(A, F8.3, F8.3, F8.3, F8.3, F8.3, F8.3)') 'qnet Results   ', M_AST(5,3), M_AST(5,5), M_AST(5,2), &
      M_AST(5,6), M_AST(5,1), M_AST(5,7)
      WRITE (6,'(A, F8.5, A, F8.5, A, F8.5)') 'Normal_Vector       ',NX, '        ', NY, '        ', NZ
   END IF
   DEALLOCATE (V_TIME)
   DEALLOCATE (TAST)
   DEALLOCATE (M_AST)
   IF (T_AVERAGE/=0) DEALLOCATE (MED_TAST)  
    !CLOSE (100+VAR)
   100 CONTINUE
   ENDDO VARIABLE_NUMBER
!
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
   IF (JUNK(1:LENGTH)==STRING(1:LENGTH)) EXIT SEARCH_LOOP
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
   IF (JUNK(1:LENGTH)==STRING(1:LENGTH).OR.JUNK(1:LENGTH2)==STRING2(1:LENGTH2)) THEN
      IF (JUNK(1:LENGTH) ==STRING(1:LENGTH))   CHOICE = JUNK(1:LENGTH)
      IF (JUNK(1:LENGTH2)==STRING2(1:LENGTH2)) CHOICE = JUNK(1:LENGTH2)
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
   IF(C>='a'.AND.C<='z')C=CHAR(ICHAR(C)+ICHAR('A')-ICHAR('a'))
    BUFFEROUT(I:I)=C
END DO

END SUBROUTINE TOUPPER
