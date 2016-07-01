SUBROUTINE INTERFACE(CHID,TBEG,TEND,TINT,HC,DT,STEPS,ANSYS_STEP1,ANSYS_STEP2)

! Program to convert exposed surface elementos into ANSYS SURF152

IMPLICIT NONE

CHARACTER(256) NODES,ELMS!,INTFILE
INTEGER NUMEL,NNODE,I,J,NO1,NO2,NO3,NO4,LINHAS,VARIABLE,SHELL,EL,LAYER,HIGHNODE,N_OR,ROUND,L_LAYER
INTEGER ST,ENDI,STEPS
REAL SOMAX,SOMAY,SOMAZ,A(3),B(3),C(3)
REAL, ALLOCATABLE, DIMENSION(:,:) :: NOS,NOMASTER
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ELEMENTOS
REAL TINT,TBEG,TEND,HC,DT
CHARACTER(256) CHID,ANSYS_STEP1,ANSYS_STEP2
CHARACTER(8) INTFILE2,INTFILE3

SHELL=0
ROUND=1
LAYER=0
L_LAYER=0
ST=1

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
ENDI=NUMEL

! Allocate the pointers
ALLOCATE (NOMASTER(NUMEL,4))
NOMASTER=0.D0
ALLOCATE (ELEMENTOS(NUMEL,9))
ELEMENTOS=0.D0

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
!************* DATA FOR SUBROUTINE FDS2AST ***************
!WRITE(6,*) ' Enter Job ID string (CHID):'
!READ(*,'(a)') CHID
! Beggining of the inputs for the new code 
!   WRITE(6,*) ' Enter number of variables '
!   READ(*,*) VARIABLE
VARIABLE=1
!         WRITE(6,*) ' Enter starting and ending time (s)'
!         READ(*,*) TBEG,TEND
!         WRITE(6,*) ' Enter time interval (s)'
!         READ(*,*) TINT
!         WRITE(6,*) ' Specify the convective heat transfer coefficient (W/(m².K))'
!         READ(*,*) HC
!*********************************************************

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

!        N(I,1)=NOMASTER(I,1)
!        N(I,2)=N_OR*((B(2)-A(2))*(C(3)-A(3))-(B(3)-A(3))*(C(2)-A(2)))
!        N(I,3)=N_OR*((B(3)-A(3))*(C(1)-A(1))-(B(1)-A(1))*(C(3)-A(3)))
!        N(I,4)=N_OR*((B(1)-A(1))*(C(2)-A(2))-(B(2)-A(2))*(C(1)-A(1)))

!        IF (SHELL.EQ.1) THEN
!            VEC_MULT=(C_SIZE/2)/(SQRT((N(I,2)**2)+(N(I,3)**2)+(N(I,4)**2)))
!            NOMASTER(I,2)=NOMASTER(I,2)+(N_OR*(N(I,2)*VEC_MULT))
!            NOMASTER(I,3)=NOMASTER(I,3)+(N_OR*(N(I,3)*VEC_MULT))
!            NOMASTER(I,4)=NOMASTER(I,4)+(N_OR*(N(I,4)*VEC_MULT))
!        END IF

ENDDO CREATE_NOMASTER

OPEN(70,FILE=TRIM(CHID)//'_thermal_boundary.dat',FORM='FORMATTED',STATUS='UNKNOWN')
OPEN(3,FILE='run_case.dat',FORM='FORMATTED',STATUS='UNKNOWN')
WRITE(3,'(A)') "/batch"
IF (TBEG==0) THEN
  WRITE(3,'(A,A,A)') "resume,", TRIM(ANSYS_STEP1), ",db"
ELSE
  WRITE(3,'(A)') "resume"
ENDIF
WRITE(3,'(A)') "!"
WRITE(3,'(A)') "/output,ansys_output,dat"
WRITE(3,'(A,A,A)') "/input,", TRIM(CHID),'_thermal_boundary,dat'
WRITE(3,'(A)') "ncnv,2"
WRITE(3,'(A)') "/solu"
WRITE(3,'(A)') "solve"
WRITE(3,'(A)') "finish"

IF (STEPS.EQ.2) THEN
   OPEN(71,FILE=TRIM(CHID)//'_mechanical_run.dat',FORM='FORMATTED',STATUS='UNKNOWN')
   WRITE(3,'(A)') "save"
   WRITE(3,'(A,A,A)') "/input,", TRIM(CHID),'_mechanical_run,dat'
   WRITE(3,'(A)') "/solu"
   WRITE(3,'(A)') "solve"
   WRITE(3,'(A)') "finish"
END IF
WRITE(3,'(A,A,A)') "*set,DISPNAME,'", TRIM(CHID), "_ansys_displacements'"
WRITE(3,'(A)') "/input,plot_displacements,rsts"
WRITE(3,'(A)') "save"
!********** Writing results into displacements file ***********
OPEN(5,FILE='plot_displacements.rsts',FORM='FORMATTED',STATUS='UNKNOWN')
WRITE(5,'(A)') '/post1'
WRITE(5,'(A)') 'esel,S,type,,1'                               ! Select just the solid elements (exclude surface effect elements)
WRITE(5,'(A)') 'nsle,r,corner'                                ! Select just the corner nodes (exclude midside nodes)
WRITE(5,'(A)') '*get,NCOUNT,node,,count'                      ! Get total number of selected nodes
WRITE(5,'(A)') '*get,ECOUNT,elem,,count'                      ! Get total number of selected elements
WRITE(5,'(A)') '*dim,NARRAY,array,NCOUNT,3'                   ! Create NCOUNT x 3 array
WRITE(5,'(A)') '*dim,EARRAY,array,ECOUNT,8'                   ! Create ECOUNT x 8 array
WRITE(5,'(A)') '*dim,N2V,array,NCOUNT,1'                      ! Create NCOUNT x 1 array
WRITE(5,'(A)') '*dim,ELM,array,ECOUNT,1'                      ! Create ECOUNT x 1 array
WRITE(5,'(A)') '*get,NUMBER,active,,set,nset'                 ! Get number of substeps
WRITE(5,'(A)') '*cfopen,DISPNAME,dat'                         ! Create file called “displacements.dat”
WRITE(5,'(A)') 'set,,, ,,, ,NUMBER'
WRITE(5,'(A)') '*get,CONV,active,,solu,cnvg'                  ! Get the solution status - 0 unconverged, 1 converged
WRITE(5,'(A)') '*if,CONV,EQ,0,THEN'
WRITE(5,'(A)') '   SET,LAST'
WRITE(5,'(A)') '   SET,PREVIOUS'
WRITE(5,'(A)') '*endif'
WRITE(5,'(A)') '*get,TIME,active,,set,TIME'                   ! Get the time of the current substep
WRITE(5,'(A)') '*vget,N2V(1),node,,nlist'                     ! Fill N2V with node ID, VERTs ID is the vector line
WRITE(5,'(A)') '*do,I,1,NCOUNT'
WRITE(5,'(A)') '   *get,NARRAY(I,1),node,N2V(I),U,x'             ! Fill first column with x-coord.
WRITE(5,'(A)') '   *get,NARRAY(I,2),node,N2V(I),U,y'             ! Fill second column with y-coord.
WRITE(5,'(A)') '   *get,NARRAY(I,3),node,N2V(I),U,z'             ! Fill third column with z-coord.
WRITE(5,'(A)') '*enddo'
WRITE(5,'(A)') '*do,I,1,ncount'
WRITE(5,'(A)') '   NARRAY(I,1)=NARRAY(I,1)'       
WRITE(5,'(A)') '   NARRAY(I,2)=NARRAY(I,2)'
WRITE(5,'(A)') '   NARRAY(I,3)=NARRAY(I,3)'
WRITE(5,'(A)') '*enddo'
WRITE(5,'(A)') '*vwrite,TIME'                                 ! Writing Nodes'
WRITE(5,'(A)') '%8.3F'
WRITE(5,'(A)') '*vwrite,N2V(1),NARRAY(1,1),NARRAY(1,2),NARRAY(1,3)'            ! Write NARRAY to file
WRITE(5,'(A)') '%8I%8.5F%8.5F%8.5F' 
WRITE(5,'(A)') '*cfclose'                                     ! Close file called “displacements.dat”
WRITE(5,'(A)') '*del,NCOUNT'                                  ! Delete created variables
WRITE(5,'(A)') '*del,ECOUNT'
WRITE(5,'(A)') '*del,NARRAY'
WRITE(5,'(A)') '*del,EARRAY'
WRITE(5,'(A)') '*del,N2V'
WRITE(5,'(A)') '*del,ELM'
WRITE(5,'(A)') '*del,NUMBER'
WRITE(5,'(A)') '*del,TIME'
CLOSE(5)

IF (TBEG.NE.0) GO TO 15
!************** CREATE EXTRA NODES *********************
   WRITE(70,'(A)') '/PREP7' 
DO I=1,NUMEL
   WRITE(70,'(A, I8, A, F7.3, A, F7.3, A, F7.3, A)') "N,", INT(NOMASTER(I,1)), ",", NOMASTER(I,2), ",", &
   NOMASTER(I,3), ",", NOMASTER(I,4), ",,,,"
!  N,"NODE NUMBER","COORD X","COORD Y","COORD Z",,,,
ENDDO
15  CONTINUE
!***********************************
!!WORKING!!
!*******************
IF (LAYER.EQ.3) THEN 
L_LAYER=1
END IF
!**********
IF (TBEG.NE.0) THEN
   WRITE(70,'(A)') "/prep7"
   WRITE(70,'(A)') "antype,,rest,last,last,0"
   WRITE(70,'(A)') "!*"
ELSE
   WRITE(70,'(A)') "/SOLU"
   WRITE(70,'(A)') "antype,4"
   WRITE(70,'(A)') "!*"
   WRITE(70,'(A)') "TRNOPT,FULL"
   WRITE(70,'(A)') "TLUMPM,0"
   WRITE(70,'(A,I8,A,I8,A,I8)') "DELTIM,", INT(DT/2),",", INT(DT/10),",", INT(DT)
   WRITE(70,'(A)') "OUTRES,ERASE"
   WRITE(70,'(A,I8)') "OUTRES,BASI,", 10   
ENDIF
   WRITE(70,'(A, I8)')"time,", INT(TEND) 
   WRITE(70,'(A)') "!*"
   WRITE(70,'(A)') "/PREP7"
!
IF (STEPS.EQ.2) THEN
   WRITE(71,'(A, A, A)') "/filname,", TRIM(CHID), "_ansys_mechanical ,0" 
   IF (TBEG==0) THEN
      WRITE(71,'(A,A,A)') "resume,", TRIM(ANSYS_STEP2), ",db"
      WRITE(71,'(A)') "/SOLU"
      WRITE(71,'(A)') "ncnv,2"
      WRITE(71,'(A)') "antype,0,new"
      WRITE(71,'(A)') "!*"   
   ELSE
      WRITE(71,'(A)') "resume"
      WRITE(71,'(A)') "/SOLU"
      WRITE(71,'(A)') "ncnv,2"
      WRITE(71,'(A)') "antype,,rest,last,last,0"
      WRITE(71,'(A)') "!*"
   END IF
   WRITE(71,'(A)') "/solu"
   WRITE(71,'(A,I8,A,I8,A,I8)') "DELTIM,",INT(DT/2),",",INT(DT/10),",",INT(DT)  
   WRITE(71,'(A)') "OUTRES,ERASE"
   WRITE(71,'(A,I8)') "OUTRES,BASI,", 10  
   WRITE(71,'(A, I8)')"time,", INT(TEND) 
   WRITE(71,'(A)') "!*"
   WRITE(71,'(A, I8, A, A, A)') "ldread,temp,,,",INT(TEND),",,", TRIM(CHID),"_ansys,rth "
   CLOSE(71)
ENDIF
20 CONTINUE
IF (TBEG.NE.0) GOTO 25
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
IF (SHELL.EQ.1 .AND. LAYER.EQ.3) THEN
   IF (L_LAYER.EQ.1) WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",11,1"
   IF (L_LAYER.EQ.2)  WRITE(70,'(A, I8, A)') "KEYOPT,", EL, ",11,2"
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
IF (LAYER.EQ.3 .AND. L_LAYER.EQ.1) THEN 
ST=1
ENDI=NUMEL/2
END IF
IF (LAYER.EQ.3 .AND. L_LAYER.EQ.2) THEN 
ST=(NUMEL/2)+1
ENDI=NUMEL
END IF

LOOP_SURF152:DO I=ST,ENDI
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
25 CONTINUE
PRINT *, 'LOOP_TABLES'
IF (TBEG.EQ.0) THEN
LOOP_TABLES_BEG:DO I=ST,ENDI
    WRITE (INTFILE2,'(g8.0)') INT(NOMASTER(I,1))
!    PRINT *, INTFILE2
    LINHAS=INT((TEND-TBEG)/TINT)+1
    WRITE(70,'(A, A, A, I8, A)') "*DIM,A", INTFILE2, ",TABLE,", LINHAS, ",1,1,TIME,TEMP," 
    WRITE(70,'(A)') "!*"
    IF (VARIABLE.EQ.2) WRITE(70,'(A, A, A, I8, A)') "*DIM,H", INTFILE2, ",TABLE,", LINHAS, ",1,1,TIME,," 
    IF (VARIABLE.EQ.2) WRITE(70,'(A)') "!*"
    !WRITE(70,'(A, A, A, A, A)') "*TREAD,A", INTFILE2, ",'", INTFILE2, "','dat',' ', ,"
    !WRITE(70,'(A)') "!*"
    !IF (VARIABLE.EQ.2) WRITE(70,'(A, A, A, A, A)') "*TREAD,H", INTFILE2, ",'", INTFILE2, "H','dat',' ', ,"
    !IF (VARIABLE.EQ.2) WRITE(70,'(A)') "!*"
ENDDO LOOP_TABLES_BEG  
ENDIF
!  
IF (TBEG.NE.0) THEN
LOOP_TABLES:DO I=ST,ENDI
    WRITE (INTFILE2,'(g8.0)') INT(NOMASTER(I,1))
!    PRINT *, INTFILE2
    !WRITE(70,'(A, A, A, A, A)') "*TREAD,A", INTFILE2, ",'", INTFILE2, "','dat',' ', ,"
    !WRITE(70,'(A)') "!*"
    !IF (VARIABLE.EQ.2) WRITE(70,'(A, A, A, A, A)') "*TREAD,H", INTFILE2, ",'", INTFILE2, "H','dat',' ', ,"
    !IF (VARIABLE.EQ.2) WRITE(70,'(A)') "!*"
ENDDO LOOP_TABLES 
ENDIF     
!****************************
!!WORKING!!
!*********** APPLYING TABLES AS NODAL LOADS **********
    WRITE(70,'(A)') "/PREP7"
    PRINT *, 'LOOP_LOADS'
LOOP_LOADS:DO I=ST,ENDI
    WRITE (INTFILE2,'(g8.0)') INT(NOMASTER(I,1))
!    PRINT *, INTFILE2
    WRITE(70,'(A)') "!*"
    WRITE(70,'(A, A, A, A, A)') "D,",INTFILE2,", , %A", INTFILE2, "% , , , ,TEMP, , , , ,"
ENDDO LOOP_LOADS
!*******************************
!!WORKING!!
!*********** APPLYING CONVECTIVE HEAT TRANSFER COEFFICIENT **********
    WRITE(70,'(A)') "/PREP7"
    PRINT *, 'LOOP_HEAT_TRANSFER'
LOOP_HEAT:DO I=ST,ENDI
    WRITE (INTFILE3,'(g8.0)') INT(ELEMENTOS(I,1))
!    PRINT *, INTFILE3
    WRITE(70,'(A,A,A,F7.3)') "SFE,",INTFILE3,",1,CONV,0,",HC
ENDDO LOOP_HEAT
!*******************************
IF (LAYER.EQ.3 .AND. L_LAYER.EQ.1) THEN
    L_LAYER=2
    EL=EL+1
    ROUND=ROUND+1
    GO TO 20
END IF    


!****** LOOP NO FDS2AST PARA RETIRAR OS RESULTADOS DO FDS *****
PRINT *, 'LOOP_FDS2AST'
    CALL BE2FTMI (CHID,NOMASTER,TBEG,TEND,TINT,NUMEL,SHELL)
!*******************
!!WORKING!!


!READ(2,'(a)') FILE_END
!IF (FILE_END=='END') THEN
!WRITE(6,*) FILE_END
!ELSE
!    DEALLOCATE (NOMASTER)
!    DEALLOCATE (ELEMENTOS)
!    BACKSPACE 2
!    READ(2,*) NUMEL,SHELL,LAYER,EL
!    ALLOCATE (NOMASTER(NUMEL,4))
!    NOMASTER=0.D0
!    ALLOCATE (ELEMENTOS(NUMEL,9))
!    ELEMENTOS=0.D0
!    ROUND=2
!    GO TO 10 
!END IF
!**********************************
CLOSE (1)
CLOSE (2)
CLOSE (3)
CLOSE (70)
!**********************************
END SUBROUTINE INTERFACE