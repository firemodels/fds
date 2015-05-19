PROGRAM interface

! Program to convert exposed surface elementos into ANSYS SURF152

IMPLICIT NONE

CHARACTER(256) NODES,ELMS,INTFILE
INTEGER NUMEL,NNODE,I,J,NO1,NO2,NO3,NO4,LINHAS,VARIABLE
INTEGER SHELL,EL,LAYER,HIGHNODE,N_OR,ROUND,CODE,N_AVERAGE
REAL SOMAX,SOMAY,SOMAZ,A(3),B(3),C(3)
REAL, ALLOCATABLE, DIMENSION(:,:) :: NOS,NOMASTER,N
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ELEMENTOS
REAL C_SIZE,TINT,TBEG,TEND,VEC_MULT,HC
CHARACTER(256) CHID,INPUT
CHARACTER(8) INTFILE2,INTFILE3,FILE_END,ARG

C_SIZE=0.0
TINT=0.0
TBEG=0.0
TEND=0.0
SHELL=0
ROUND=1


! Open node file        
!WRITE(6,*) ' Enter node file ID :'
!READ(*,'(a)') NODES
NODES='nodes'
OPEN(1,FILE=TRIM(NODES)//'.dat',STATUS='OLD',FORM='FORMATTED')
!WRITE(6,*) ' Enter number of nodes :'
READ(1,*) NNODE,HIGHNODE
ALLOCATE (NOS(NNODE,4))

! Open element file
!WRITE(6,*) ' Enter element file ID :'
!READ(*,'(a)') ELMS
ELMS='elements'
OPEN(2,FILE=TRIM(ELMS)//'.dat',STATUS='OLD',FORM='FORMATTED')
!WRITE(6,*) ' Enter number of elements :'
!READ(2,'(I8)',ADVANCE='NO') NUMEL
READ(2,*) NUMEL,SHELL,LAYER,EL

! Allocate the pointers
ALLOCATE (NOMASTER(NUMEL,4))
NOMASTER=0.D0
ALLOCATE (ELEMENTOS(NUMEL,9))
ELEMENTOS=0.D0
ALLOCATE (N(NUMEL,4))
N=0.D0

!Shell elements?
!WRITE(6,*) ' Shell elements? (1)YES (2) NO :'
!READ(2,'(I8)',ADVANCE='NO') SHELL

!IF (SHELL.EQ.1) THEN
!    WRITE(6,*) ' Attach AST to which layer? (1)Top (2)Bottom :'
!     READ(2,'(I8)',ADVANCE='NO') LAYER
!END IF
!
!WRITE(6,*) ' Write number of SURF152 element at ansys :'
!READ(2,*) EL
!
!Normal Orientation (1) for outside the solid and (-1) for inside the solid
    N_OR=1
!

!********TRY TO GET INFORMATION THROUGH ARGUMENTS*********
CALL GET_COMMAND_ARGUMENT(1,INPUT)
CALL GET_COMMAND_ARGUMENT(2,ARG)
IF (LEN_TRIM(ARG)/=0) THEN
CHID=INPUT
READ(ARG,*) C_SIZE
CALL GET_COMMAND_ARGUMENT(3,ARG)
READ(ARG,*) VARIABLE
IF (VARIABLE==1) THEN
    CALL GET_COMMAND_ARGUMENT(4,ARG)
    READ(ARG,*) HC
    CALL GET_COMMAND_ARGUMENT(5,ARG)
    READ(ARG,*) TBEG
    CALL GET_COMMAND_ARGUMENT(6,ARG)
    READ(ARG,*) TEND
    CALL GET_COMMAND_ARGUMENT(7,ARG)
    READ(ARG,*) TINT
    CALL GET_COMMAND_ARGUMENT(8,ARG)
    READ(ARG,*) N_AVERAGE
    CALL GET_COMMAND_ARGUMENT(9,INTFILE)
ELSE
    CALL GET_COMMAND_ARGUMENT(4,ARG)
    READ(ARG,*) TBEG
    CALL GET_COMMAND_ARGUMENT(5,ARG)
    READ(ARG,*) TEND
    CALL GET_COMMAND_ARGUMENT(6,ARG)
    READ(ARG,*) TINT
    CALL GET_COMMAND_ARGUMENT(7,ARG)
    READ(ARG,*) N_AVERAGE
    CALL GET_COMMAND_ARGUMENT(8,INTFILE)
ENDIF
!********TRY TO GET INFORMATION THROUGH INPUTFILE*********
ELSEIF (LEN_TRIM(INPUT)/=0) THEN
OPEN(3,FILE=TRIM(INPUT),STATUS='OLD',FORM='FORMATTED')
READ(3,'(a)') CHID
READ(3,*) C_SIZE
READ(3,*) VARIABLE,HC
READ(3,*) TBEG,TEND,TINT
READ(3,*) N_AVERAGE
READ(3,'(a)') INTFILE
CLOSE(3)
ENDIF
!***********GET INFORMATION THROUGH USER INPUT************
IF (LEN_TRIM(CHID)==0) THEN
WRITE(6,*) ' Enter Job ID string (CHID):'
READ(*,'(a)') CHID
! Beggining of the inputs for the new code 
   WRITE(6,*) ' Enter Cell Size (m)'
   READ(*,*) C_SIZE
   WRITE(6,*) ' Enter number of variables '
   READ(*,*) VARIABLE
   IF (VARIABLE==1) THEN
        WRITE(6,*) ' Specify the convective heat transfer coefficient (W/(m².K))'
        READ(*,*) HC
   ENDIF
! Rest of the input parameters          
         WRITE(6,*) ' Enter starting and ending time (s)'
         READ(*,*) TBEG,TEND
         WRITE(6,*) ' Enter time interval (s)'
         READ(*,*) TINT
         WRITE(6,*) ' Enter number of average points to consider constant steady-state exposure conditions (0 if NO average)'
         READ(*,*) N_AVERAGE
END IF         
!*********************************************************
! Read nodes and elements
READ_NODES: DO I=1,NNODE
            READ (1,*) NOS(I,1),NOS(I,2),NOS(I,3),NOS(I,4)
ENDDO READ_NODES

10  CONTINUE

READ_ELMS: DO I=1,NUMEL
            READ (2,*) ELEMENTOS(I,1),ELEMENTOS(I,2),ELEMENTOS(I,3),ELEMENTOS(I,4), &
            ELEMENTOS(I,5),ELEMENTOS(I,6),ELEMENTOS(I,7),ELEMENTOS(I,8),ELEMENTOS(I,9)
ENDDO READ_ELMS

CREATE_NOMASTER: DO I=1,NUMEL
HIGHNODE=HIGHNODE+1
            NO1=ELEMENTOS(I,2)
            NO2=ELEMENTOS(I,3)
            NO3=ELEMENTOS(I,4)
            NO4=ELEMENTOS(I,5)
            SOMAX=0.D0
            SOMAY=0.D0
            SOMAZ=0.D0
                    DO J=1,NNODE
                    IF (NO1.EQ.NOS(J,1)) THEN
                    SOMAX=SOMAX+NOS(J,2)
                    SOMAY=SOMAY+NOS(J,3)
                    SOMAZ=SOMAZ+NOS(J,4)
                    A(1)=NOS(J,2)
                    A(2)=NOS(J,3)
                    A(3)=NOS(J,4)
                    ELSE IF (NO2.EQ.NOS(J,1)) THEN
                    SOMAX=SOMAX+NOS(J,2)
                    SOMAY=SOMAY+NOS(J,3)
                    SOMAZ=SOMAZ+NOS(J,4)
                    B(1)=NOS(J,2)
                    B(2)=NOS(J,3)
                    B(3)=NOS(J,4)
                    ELSE IF (NO3.EQ.NOS(J,1)) THEN
                    SOMAX=SOMAX+NOS(J,2)
                    SOMAY=SOMAY+NOS(J,3)
                    SOMAZ=SOMAZ+NOS(J,4)
                    C(1)=NOS(J,2)
                    C(2)=NOS(J,3)
                    C(3)=NOS(J,4)
                    ELSE IF (NO4.NE.NO3 .AND. NO4.EQ.NOS(J,1)) THEN
                    SOMAX=SOMAX+NOS(J,2)
                    SOMAY=SOMAY+NOS(J,3)
                    SOMAZ=SOMAZ+NOS(J,4)
                    END IF
                    END DO
            IF (NO4.EQ.NO3) THEN
            NOMASTER(I,1)=HIGHNODE
            NOMASTER(I,2)=SOMAX/3
            NOMASTER(I,3)=SOMAY/3
            NOMASTER(I,4)=SOMAZ/3
            ELSE 
            NOMASTER(I,1)=HIGHNODE
            NOMASTER(I,2)=SOMAX/4
            NOMASTER(I,3)=SOMAY/4
            NOMASTER(I,4)=SOMAZ/4
            END IF

        N(I,1)=NOMASTER(I,1)
        N(I,2)=N_OR*((B(2)-A(2))*(C(3)-A(3))-(B(3)-A(3))*(C(2)-A(2)))
        N(I,3)=N_OR*((B(3)-A(3))*(C(1)-A(1))-(B(1)-A(1))*(C(3)-A(3)))
        N(I,4)=N_OR*((B(1)-A(1))*(C(2)-A(2))-(B(2)-A(2))*(C(1)-A(1)))

        IF (SHELL.EQ.1) THEN
            VEC_MULT=(C_SIZE/2)/(SQRT((N(I,2)**2)+(N(I,3)**2)+(N(I,4)**2)))
            NOMASTER(I,2)=NOMASTER(I,2)+(N_OR*(N(I,2)*VEC_MULT))
            NOMASTER(I,3)=NOMASTER(I,3)+(N_OR*(N(I,3)*VEC_MULT))
            NOMASTER(I,4)=NOMASTER(I,4)+(N_OR*(N(I,4)*VEC_MULT))
        END IF

ENDDO CREATE_NOMASTER

IF (ROUND.EQ.1) THEN
IF (LEN_TRIM(INTFILE)==0) THEN
WRITE(6,*) ' Enter output file name :'
READ(*,'(a)') INTFILE
END IF
OPEN(70,FILE=TRIM(INTFILE)//'.dat',FORM='FORMATTED',STATUS='UNKNOWN')
END IF
!***** SELECT CASE - WHICH FEM CODE ARE YOU USING! *****
CODE =1
! 1- ANSYS, 2- .....
SELECT CASE (CODE)
CASE (1) ! ANSYS
!************** CREATE EXTRA NODES *********************
   WRITE(70,'(A)') '/PREP7' 
DO I=1,NUMEL
   WRITE(70,'(A, I8, A, F7.3, A, F7.3, A, F7.3, A)') "N,", INT(NOMASTER(I,1)), ",", NOMASTER(I,2), ",", &
   NOMASTER(I,3), ",", NOMASTER(I,4), ",,,,"
!  N,"NODE NUMBER","COORD X","COORD Y","COORD Z",,,,
ENDDO
!***********************************
!!WORKING!!
!**************** CREATE ELEMENT TYPE SURF152 ********************
!IF (SHELL.NE.1) EL=5
   WRITE(70,'(A)') "!*"
   WRITE(70,'(A, I8, A)') "ET,", EL, ",SURF152"
   WRITE(70,'(A)') "!*"
   WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",1,0"
   WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",2,0"
   WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",3,0"
   WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",4,0"
   WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",5,1"
   WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",6,0"
   WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",7,0"
   WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",8,4"
   WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",9,1"
! IF IT IS SHELL ELEMENTS,
! 1 FOR TOP LAYER AND 2 FOR BOTTOM LAYER !
IF (SHELL.EQ.1) THEN
   IF (LAYER.EQ.1) WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",11,1"
   IF (LAYER.EQ.2)  WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",11,2"
END IF
!
   WRITE(70,'(A)') "!*"
   WRITE(70,'(A)') "!*"
   WRITE(70,'(A, I8, A)') "R,", EL, ",1,5.67E-08, , , ,"
   WRITE(70,'(A)') "RMORE, , , ,"
   WRITE(70,'(A)') "RMORE, , ,"
   WRITE(70,'(A)') "!*"
!****************************
!!WORKING!!
!************* CREATE SURF152 ELEMENTS ********
PRINT *, 'LOOP_SURF152'
    WRITE(70,'(A)') "!*"
    WRITE(70,'(A, I8)') "TYPE,", EL
    WRITE(70,'(A, I8)') "MAT,", 1
    WRITE(70,'(A, I8)') "REAL,", EL
    WRITE(70,'(A, I8)') "ESYS,", 0
    WRITE(70,'(A)') "SECNUM,"
    WRITE(70,'(A)') "TSHAP,LINE"
    WRITE(70,'(A)') "!*"
LOOP_SURF152:DO I=1,NUMEL
    IF (ELEMENTOS(I,6).EQ.0) THEN
    WRITE(70,'(A, I8)') "nsel,S,node,,", ELEMENTOS(I,2)
    WRITE(70,'(A, I8)') "nsel,A,node,,", ELEMENTOS(I,3)
    WRITE(70,'(A, I8)') "nsel,A,node,,", ELEMENTOS(I,4)
    WRITE(70,'(A, I8)') "nsel,A,node,,", ELEMENTOS(I,5)
    WRITE(70,'(A, I8)') "ESURF,", INT(NOMASTER(I,1))
    WRITE(70,'(A)') "!*"
    ELSE
    WRITE(70,'(A, I8)') "nsel,S,node,,", ELEMENTOS(I,2)
    WRITE(70,'(A, I8)') "nsel,A,node,,", ELEMENTOS(I,3)
    WRITE(70,'(A, I8)') "nsel,A,node,,", ELEMENTOS(I,4)
    WRITE(70,'(A, I8)') "nsel,A,node,,", ELEMENTOS(I,5)
    WRITE(70,'(A, I8)') "nsel,A,node,,", ELEMENTOS(I,6)
    WRITE(70,'(A, I8)') "nsel,A,node,,", ELEMENTOS(I,7)
    WRITE(70,'(A, I8)') "nsel,A,node,,", ELEMENTOS(I,8)
    WRITE(70,'(A, I8)') "nsel,A,node,,", ELEMENTOS(I,9)
    WRITE(70,'(A, I8)') "ESURF,", INT(NOMASTER(I,1))
    WRITE(70,'(A)') "!*"
    END IF
    WRITE(70,'(A)') "ALLSELL,ALL"
    WRITE(70,'(A)') "!*"
ENDDO LOOP_SURF152
!*******************************
!!WORKING!!
!*********** CREATING TABLES (.DAT) **********
PRINT *, 'LOOP_TABELAS'
LOOP_TABELAS:DO I=1,NUMEL
    WRITE (INTFILE2,'(g8.0)') INT(NOMASTER(I,1))
    PRINT *, INTFILE2
    IF (N_AVERAGE==0) THEN
        LINHAS=INT((TEND-TBEG)/TINT)
    ELSE
        LINHAS=INT((TEND-TBEG)/TINT)+2
    ENDIF
    WRITE(70,'(A, A, A, I8, A)') "*DIM,A", INTFILE2, ",TABLE,", LINHAS, ",1,1,TIME,TEMP," 
    WRITE(70,'(A)') "!*"
    IF (VARIABLE.EQ.2) WRITE(70,'(A, A, A, I8, A)') "*DIM,H", INTFILE2, ",TABLE,", LINHAS, ",1,1,TIME,," 
    IF (VARIABLE.EQ.2) WRITE(70,'(A)') "!*"
    !WRITE(70,'(A, A, A, A, A)') "*TREAD,A", INTFILE2, ",'", INTFILE2, "','dat',' ', ,"
    !WRITE(70,'(A)') "!*"
    !IF (VARIABLE.EQ.2) WRITE(70,'(A, A, A, A, A)') "*TREAD,H", INTFILE2, ",'", INTFILE2, "H','dat',' ', ,"
    !IF (VARIABLE.EQ.2) WRITE(70,'(A)') "!*"
ENDDO LOOP_TABELAS    
!****************************
!!WORKING!!
!200 CONTINUE
!****** LOOP AT FDS2AST TO EXTRACT RESULTS FROM BNDF FILES *****
PRINT *, 'LOOP_FDS2AST'
CALL FDS2AST (CHID,NOMASTER,C_SIZE,TBEG,TEND,TINT,NUMEL,VARIABLE,N,N_AVERAGE)
!*******************
!!WORKING!!
!*********** APPLYING TABLES AS NODAL LOADS **********
    WRITE(70,'(A)') "/PREP7"
    PRINT *, 'LOOP_CARGAS'
LOOP_CARGAS:DO I=1,NUMEL
    WRITE (INTFILE2,'(g8.0)') INT(NOMASTER(I,1))
    PRINT *, INTFILE2
    WRITE(70,'(A)') "!*"
    WRITE(70,'(A, A, A, A, A)') "D,",INTFILE2,", , %A", INTFILE2, "% , , , ,TEMP, , , , ,"
ENDDO LOOP_CARGAS
!*******************************
!!WORKING!!
!*********** APPLYING TABLES AS CONVECTIVE HEAT TRANSFER COEFFICIENT **********
IF (VARIABLE.EQ.2) THEN
    WRITE(70,'(A)') "/PREP7"
    PRINT *, 'LOOP_HEAT_TRANSFER'
LOOP_HEAT2:DO I=1,NUMEL
    WRITE (INTFILE2,'(g8.0)') INT(NOMASTER(I,1))
    PRINT *, INTFILE2
    WRITE (INTFILE3,'(g8.0)') INT(ELEMENTOS(I,1))
    PRINT *, INTFILE3
    WRITE(70,'(A,A,A,A,A)') "SFE,",INTFILE3,",1,CONV,0,%H",INTFILE2,"%"
ENDDO LOOP_HEAT2
ENDIF
IF (VARIABLE.EQ.1) THEN
    WRITE(70,'(A)') "/PREP7"
    PRINT *, 'LOOP_HEAT_TRANSFER'
LOOP_HEAT1:DO I=1,NUMEL
    WRITE (INTFILE2,'(g8.0)') INT(NOMASTER(I,1))
    PRINT *, INTFILE2
    WRITE (INTFILE3,'(g8.0)') INT(ELEMENTOS(I,1))
    PRINT *, INTFILE3
    WRITE(70,'(A,A,A,F7.3)') "SFE,",INTFILE3,",1,CONV,0,",HC
ENDDO LOOP_HEAT1
ENDIF
!*******************************
READ(2,'(a)') FILE_END
IF (FILE_END=='END') THEN
WRITE(6,*) FILE_END
ELSE
    DEALLOCATE (NOMASTER)
    DEALLOCATE (N)
    DEALLOCATE (ELEMENTOS)
    BACKSPACE 2
    READ(2,*) NUMEL,SHELL,LAYER,EL
    ALLOCATE (NOMASTER(NUMEL,4))
    NOMASTER=0.D0
    ALLOCATE (ELEMENTOS(NUMEL,9))
    ELEMENTOS=0.D0
    ALLOCATE (N(NUMEL,4))
    N=0.D0
    ROUND=2
    GO TO 10 
END IF
!*******************************
! SPACE LEFT TO OTHER CODES IMPLEMENTATION
! JUST CHANGE THE OUPUT FORMAT AND KEEP THE 
! CODE STRUCTURE
!*******************************
END SELECT
END PROGRAM interface 