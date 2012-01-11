! Routines related to unstructured geometry and immersed boundary methods

MODULE COMPLEX_GEOMETRY

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: CHECKREAD,SHUTDOWN
USE MEMORY_FUNCTIONS, ONLY: ChkMemErr

IMPLICIT NONE

CHARACTER(255), PARAMETER :: geomid='$Id$'
CHARACTER(255), PARAMETER :: geomrev='$Revision$'
CHARACTER(255), PARAMETER :: geomdate='$Date$'

PRIVATE
PUBLIC :: INIT_IBM,TRILINEAR,GETX,GETU,GETGRAD,INIT_FACE,GET_REV_geom, &
          READ_GEOM,READ_VERT,READ_FACE,READ_VOLU,GET_VELO_IBM,GET_CUTCELL_AREA
 
CONTAINS

SUBROUTINE READ_GEOM

CHARACTER(30) :: ID,SHAPE,TFILE
REAL(EB) :: XYZ(3),ORIENTATION(3),XB(6),RADIUS,VELOCITY(3),OMEGA,RGB(3),PIXELS,ROUGHNESS
LOGICAL :: HOLE,TWO_SIDED
INTEGER :: IOS,IZERO,NG
TYPE(GEOMETRY_TYPE), POINTER :: G=>NULL()
NAMELIST /GEOM/ ID,SHAPE,XB,XYZ,ORIENTATION,RADIUS,VELOCITY,OMEGA,HOLE,RGB,PIXELS,TWO_SIDED,TFILE,ROUGHNESS

N_GEOM=0
REWIND(LU_INPUT)
COUNT_GEOM_LOOP: DO
   CALL CHECKREAD('GEOM',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_GEOM_LOOP
   READ(LU_INPUT,NML=GEOM,END=11,ERR=12,IOSTAT=IOS)
   N_GEOM=N_GEOM+1
   12 IF (IOS>0) CALL SHUTDOWN('ERROR: problem with GEOM line')
ENDDO COUNT_GEOM_LOOP
11 REWIND(LU_INPUT)

IF (N_GEOM==0) RETURN

! Allocate GEOMETRY array

ALLOCATE(GEOMETRY(N_GEOM),STAT=IZERO)
CALL ChkMemErr('READ','GEOMETRY',IZERO)

READ_GEOM_LOOP: DO NG=1,N_GEOM
   G=>GEOMETRY(NG)
   
   CALL CHECKREAD('GEOM',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_GEOM_LOOP
   
   ! Set defaults
   
   ID = 'geom'
   TFILE = 'null'
   SHAPE = 'SPHERE'
   HOLE = .FALSE.
   XYZ = 0._EB
   ORIENTATION = (/0._EB,0._EB,1._EB/)
   XB = 0._EB
   RADIUS = 1._EB
   VELOCITY = 0._EB
   OMEGA = 0._EB
   RGB = (/192,192,192/)
   PIXELS = 1.0
   TWO_SIDED = .FALSE.
   ROUGHNESS = 0._EB
   
   ! Read the GEOM line
   
   READ(LU_INPUT,GEOM,END=35)
   
   G%X1 = XB(1)
   G%X2 = XB(2)
   G%Y1 = XB(3)
   G%Y2 = XB(4)
   G%Z1 = XB(5)
   G%Z2 = XB(6)
   G%X0 = XYZ(1)
   G%Y0 = XYZ(2)
   G%Z0 = XYZ(3)
   G%X  = XYZ(1)
   G%Y  = XYZ(2)
   G%Z  = XYZ(3)
   G%XOR = ORIENTATION(1)
   G%YOR = ORIENTATION(2)
   G%ZOR = ORIENTATION(3)
   G%U0 = VELOCITY(1)
   G%V0 = VELOCITY(2)
   G%W0 = VELOCITY(3)
   G%OMEGA = OMEGA
   G%OMEGA_X = G%OMEGA*G%XOR
   G%OMEGA_Y = G%OMEGA*G%YOR
   G%OMEGA_Z = G%OMEGA*G%ZOR
   G%RADIUS = RADIUS
   G%SHAPE = SHAPE
   G%HOLE = HOLE
   G%RGB = RGB
   G%PIXELS = PIXELS
   G%TWO_SIDED = TWO_SIDED
   G%ID = ID
   G%TFILE = TFILE
   G%ROUGHNESS = ROUGHNESS
   
   IF (ABS(G%U0)>0._EB .OR. ABS(G%V0)>0._EB .OR. ABS(G%W0)>0._EB) G%TRANSLATE = .TRUE.
   IF (ABS(G%OMEGA)>0._EB) G%ROTATE = .TRUE.
   
   SELECT CASE(G%SHAPE)
      CASE('BOX')
         G%ISHAPE = IBOX
         G%X0 = 0.5_EB*(G%X1+G%X2)
         G%Y0 = 0.5_EB*(G%Y1+G%Y2)
         G%Z0 = 0.5_EB*(G%Z1+G%Z2)
         G%SMVOBJECT = 'movingbox'
      CASE('SPHERE')
         G%ISHAPE = ISPHERE
         G%SMVOBJECT = 'movingsphere'
      CASE('CYLINDER')
         G%ISHAPE = ICYLINDER
         G%SMVOBJECT = 'tube'
      CASE('PLANE')
         G%ISHAPE = IPLANE
         G%SMVOBJECT = 'plane'
      CASE DEFAULT
         CALL SHUTDOWN('ERROR: unrecognized SHAPE on GEOM line')
   END SELECT
   
   ! Allocate bounding box arrays
   
   ALLOCATE(G%MIN_I(NMESHES),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','MIN_I',IZERO)
   ALLOCATE(G%MAX_I(NMESHES),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','MAX_I',IZERO)
   ALLOCATE(G%MIN_J(NMESHES),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','MIN_J',IZERO)
   ALLOCATE(G%MAX_J(NMESHES),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','MAX_J',IZERO)
   ALLOCATE(G%MIN_K(NMESHES),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','MIN_K',IZERO)
   ALLOCATE(G%MAX_K(NMESHES),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','MAX_K',IZERO)
   
ENDDO READ_GEOM_LOOP
35 REWIND(LU_INPUT)

END SUBROUTINE READ_GEOM


SUBROUTINE READ_VERT

REAL(EB) :: X(3)=0._EB
INTEGER :: I,IOS,IZERO
NAMELIST /VERT/ X

N_VERT=0
REWIND(LU_INPUT)
COUNT_VERT_LOOP: DO
   CALL CHECKREAD('VERT',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_VERT_LOOP
   READ(LU_INPUT,NML=VERT,END=14,ERR=15,IOSTAT=IOS)
   N_VERT=N_VERT+1
   15 IF (IOS>0) CALL SHUTDOWN('ERROR: problem with VERT line')
ENDDO COUNT_VERT_LOOP
14 REWIND(LU_INPUT)

IF (N_VERT==0) RETURN

! Allocate VERTEX array

ALLOCATE(VERTEX(N_VERT),STAT=IZERO)
CALL ChkMemErr('READ','VERT',IZERO)

READ_VERT_LOOP: DO I=1,N_VERT
   
   CALL CHECKREAD('VERT',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_VERT_LOOP
   
   ! Read the VERT line
   
   READ(LU_INPUT,VERT,END=36)
   
   VERTEX(I)%X = X(1)
   VERTEX(I)%Y = X(2)
   VERTEX(I)%Z = X(3)

ENDDO READ_VERT_LOOP
36 REWIND(LU_INPUT)

END SUBROUTINE READ_VERT


SUBROUTINE READ_FACE

INTEGER :: N(3),I,IOS,IZERO,NNN,NS
LOGICAL :: EX
CHARACTER(30) :: SURF_ID='INERT'
CHARACTER(100) :: MESSAGE
NAMELIST /FACE/ N,SURF_ID

N_FACE=0
REWIND(LU_INPUT)
COUNT_FACE_LOOP: DO
   CALL CHECKREAD('FACE',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_FACE_LOOP
   READ(LU_INPUT,NML=FACE,END=16,ERR=17,IOSTAT=IOS)
   N_FACE=N_FACE+1
   16 IF (IOS>0) CALL SHUTDOWN('ERROR: problem with FACE line')
ENDDO COUNT_FACE_LOOP
17 REWIND(LU_INPUT)

IF (N_FACE==0) RETURN

! Allocate FACET array

ALLOCATE(FACET(N_FACE),STAT=IZERO)
CALL ChkMemErr('READ','FACE',IZERO)

READ_FACE_LOOP: DO I=1,N_FACE
   
   CALL CHECKREAD('FACE',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_FACE_LOOP
   
   ! Read the FACE line
   
   READ(LU_INPUT,FACE,END=37)
   
   FACET(I)%VERTEX(1) = N(1)
   FACET(I)%VERTEX(2) = N(2)
   FACET(I)%VERTEX(3) = N(3)

   ! Check the SURF_ID against the list of SURF's

   EX = .FALSE.
   DO NS=0,N_SURF
      IF (TRIM(SURF_ID)==SURFACE(NS)%ID) EX = .TRUE.
   ENDDO
   IF (.NOT.EX) THEN
      WRITE(MESSAGE,'(A,A,A)') 'ERROR: SURF_ID ',TRIM(SURF_ID),' not found'
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   ! Assign SURF_INDEX, Index of the Boundary Condition

   FACET(I)%SURF_ID = TRIM(SURF_ID)
   FACET(I)%SURF_INDEX = DEFAULT_SURF_INDEX
   DO NNN=0,N_SURF
      IF (SURF_ID==SURFACE(NNN)%ID) FACET(I)%SURF_INDEX = NNN
   ENDDO

ENDDO READ_FACE_LOOP
37 REWIND(LU_INPUT)

CALL INIT_FACE

END SUBROUTINE READ_FACE


SUBROUTINE READ_VOLU

INTEGER :: N(4),I,IOS,IZERO
CHARACTER(30) :: MATL_ID
NAMELIST /VOLU/ N,MATL_ID

N_VOLU=0
REWIND(LU_INPUT)
COUNT_VOLU_LOOP: DO
   CALL CHECKREAD('VOLU',LU_INPUT,IOS)
   IF (IOS==1) EXIT COUNT_VOLU_LOOP
   READ(LU_INPUT,NML=VOLU,END=18,ERR=19,IOSTAT=IOS)
   N_VOLU=N_VOLU+1
   18 IF (IOS>0) CALL SHUTDOWN('ERROR: problem with VOLU line')
ENDDO COUNT_VOLU_LOOP
19 REWIND(LU_INPUT)

IF (N_VOLU==0) RETURN

! Allocate VOLUME array

ALLOCATE(VOLUME(N_VOLU),STAT=IZERO)
CALL ChkMemErr('READ','VOLU',IZERO)

READ_VOLU_LOOP: DO I=1,N_VOLU
   
   CALL CHECKREAD('VOLU',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_VOLU_LOOP
   
   ! Read the VOLU line
   
   READ(LU_INPUT,VOLU,END=38)
   
   VOLUME(I)%VERTEX(1) = N(1)
   VOLUME(I)%VERTEX(2) = N(2)
   VOLUME(I)%VERTEX(3) = N(3)
   VOLUME(I)%VERTEX(4) = N(4)
   VOLUME(I)%MATL_ID = TRIM(MATL_ID)

ENDDO READ_VOLU_LOOP
38 REWIND(LU_INPUT)

END SUBROUTINE READ_VOLU


SUBROUTINE INIT_IBM(T,NM)
USE COMP_FUNCTIONS, ONLY: GET_FILE_NUMBER
IMPLICIT NONE
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
INTEGER :: I,J,K,N,IERR,I_MIN,I_MAX,J_MIN,J_MAX,K_MIN,K_MAX,IC,DUMMY_INTEGER,IZERO,LU
TYPE (MESH_TYPE), POINTER :: M
REAL(EB) :: BB(6),V1(3),V2(3),V3(3)
LOGICAL :: EX,OP
CHARACTER(60) :: FN
REAL(FB) :: DUMMY_FB_REAL

M => MESHES(NM)

CUTCELL_TEST: IF (ABS(T-T_BEGIN)<ZERO_P) THEN

FACE_LOOP: DO N=1,N_FACE

   V1 = (/VERTEX(FACET(N)%VERTEX(1))%X,VERTEX(FACET(N)%VERTEX(1))%Y,VERTEX(FACET(N)%VERTEX(1))%Z/)
   V2 = (/VERTEX(FACET(N)%VERTEX(2))%X,VERTEX(FACET(N)%VERTEX(2))%Y,VERTEX(FACET(N)%VERTEX(2))%Z/)
   V3 = (/VERTEX(FACET(N)%VERTEX(3))%X,VERTEX(FACET(N)%VERTEX(3))%Y,VERTEX(FACET(N)%VERTEX(3))%Z/)

   BB(1) = MIN(V1(1),V2(1),V3(1))
   BB(2) = MAX(V1(1),V2(1),V3(1))
   BB(3) = MIN(V1(2),V2(2),V3(2))
   BB(4) = MAX(V1(2),V2(2),V3(2))
   BB(5) = MIN(V1(3),V2(3),V3(3))
   BB(6) = MAX(V1(3),V2(3),V3(3))

   I_MIN = MAX(1,FLOOR((BB(1)-M%XS)/M%DX(1))-1)
   J_MIN = MAX(1,FLOOR((BB(3)-M%YS)/M%DY(1))-1)
   K_MIN = MAX(1,FLOOR((BB(5)-M%ZS)/M%DZ(1))-1)

   I_MAX = MIN(M%IBAR,CEILING((BB(2)-M%XS)/M%DX(1))+1)
   J_MAX = MIN(M%JBAR,CEILING((BB(4)-M%YS)/M%DY(1))+1)
   K_MAX = MIN(M%KBAR,CEILING((BB(6)-M%ZS)/M%DZ(1))+1)

   DO K=K_MIN,K_MAX
      DO J=J_MIN,J_MAX
         DO I=I_MIN,I_MAX

            IC = (K-1)*M%IBAR*M%JBAR + (J-1)*M%IBAR + I

            BB(1) = M%X(I-1)
            BB(2) = M%X(I)
            BB(3) = M%Y(J-1)
            BB(4) = M%Y(J)
            BB(5) = M%Z(K-1)
            BB(6) = M%Z(K)
            CALL TRIANGLE_BOX_INTERSECT(IERR,V1,V2,V3,BB)
            IF (IERR==1) CALL CUTCELL_INSERT(IC,FACET(N)%CUTCELL_LIST)

         ENDDO
      ENDDO
   ENDDO

ENDDO FACE_LOOP

ENDIF CUTCELL_TEST

!print *
!LL=>FACET(1)%CUTCELL_LIST
!IF ( ASSOCIATED(LL) ) THEN
!    END_OF_LIST=.FALSE.
!    DO WHILE (.NOT.END_OF_LIST)
!       print *, LL%INDEX, LL%AREA
!       LL=>LL%NEXT
!       IF ( .NOT.ASSOCIATED(LL) ) THEN
!          print *,'done printing linked list!'
!          END_OF_LIST=.TRUE.
!       ENDIF
!    ENDDO
!ENDIF

! Read boundary condition from file

SURF_LOOP: DO N=1,N_SURF

   IF (TRIM(SURFACE(N)%BC_FILENAME)=='null') CYCLE SURF_LOOP

   IF (ABS(T-T_BEGIN)<ZERO_P .AND. LU_SFBC<0) THEN

      FN = TRIM(SURFACE(N)%BC_FILENAME)
      INQUIRE(FILE=FN,EXIST=EX,OPENED=OP,NUMBER=LU)
      IF (.NOT.EX) CALL SHUTDOWN('Error: boundary condition file does not exist.')
      IF (OP) CLOSE(LU)
      LU_SFBC = GET_FILE_NUMBER()
      OPEN(LU_SFBC,FILE=FN,ACTION='READ',FORM='UNFORMATTED')

      READ(LU_SFBC) DUMMY_INTEGER ! one
      READ(LU_SFBC) DUMMY_INTEGER ! version

      IF (ALLOCATED(FB_REAL_FACE_VALS_ARRAY)) DEALLOCATE(FB_REAL_FACE_VALS_ARRAY)
      ALLOCATE(FB_REAL_FACE_VALS_ARRAY(N_FACE),STAT=IZERO)
      CALL ChkMemErr('INIT_IBM','FB_REAL_FACE_VALS_ARRAY',IZERO)

   ENDIF

   READ(LU_SFBC) DUMMY_FB_REAL ! stime

   IF (T>=REAL(DUMMY_FB_REAL,EB)) THEN
      ! n_vert_s_vals,n_vert_d_vals,n_face_s_vals,n_face_d_vals
      READ(LU_SFBC) DUMMY_INTEGER,DUMMY_INTEGER,DUMMY_INTEGER,DUMMY_INTEGER
      READ(LU_SFBC) (FB_REAL_FACE_VALS_ARRAY(I),I=1,N_FACE)
      FACET%TMP_F = REAL(FB_REAL_FACE_VALS_ARRAY,EB) + TMPM
   ELSE
      BACKSPACE LU_SFBC
   ENDIF

ENDDO SURF_LOOP

END SUBROUTINE INIT_IBM


SUBROUTINE INIT_FACE
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
IMPLICIT NONE

INTEGER :: I
REAL(EB) :: N_VEC(3),N_LENGTH,U_VEC(3),V_VEC(3),V1(3),V2(3),V3(3)
TYPE(FACET_TYPE), POINTER :: FACE
REAL(EB), PARAMETER :: TOL=1.E-10_EB

DO I=1,N_FACE

   FACE=>FACET(I)

   V1 = (/VERTEX(FACE%VERTEX(1))%X,VERTEX(FACE%VERTEX(1))%Y,VERTEX(FACE%VERTEX(1))%Z/)
   V2 = (/VERTEX(FACE%VERTEX(2))%X,VERTEX(FACE%VERTEX(2))%Y,VERTEX(FACE%VERTEX(2))%Z/)
   V3 = (/VERTEX(FACE%VERTEX(3))%X,VERTEX(FACE%VERTEX(3))%Y,VERTEX(FACE%VERTEX(3))%Z/)

   U_VEC = V2-V1
   V_VEC = V3-V1

   CALL CROSS_PRODUCT(N_VEC,U_VEC,V_VEC)
   N_LENGTH = SQRT(DOT_PRODUCT(N_VEC,N_VEC))

   IF (N_LENGTH>TOL) THEN
      FACE%NVEC = N_VEC/N_LENGTH
   ELSE
      FACE%NVEC = 0._EB
   ENDIF

   FACE%AREA  = TRIANGLE_AREA(V1,V2,V3)
   FACE%TMP_F = SURFACE(FACE%SURF_INDEX)%TMP_FRONT
   FACE%TMP_G = TMPA

ENDDO

END SUBROUTINE INIT_FACE


! http://www.sdsc.edu/~tkaiser/f90.html#Linked lists
RECURSIVE SUBROUTINE INSERT(ITEM,ROOT) 
   IMPLICIT NONE 
   TYPE(LINKED_LIST_TYPE), POINTER :: ROOT 
   INTEGER :: ITEM,IZERO
   IF (.NOT.ASSOCIATED(ROOT)) THEN 
      ALLOCATE(ROOT,STAT=IZERO)
      CALL ChkMemErr('GEOM','ROOT',IZERO)
      NULLIFY(ROOT%NEXT) 
      ROOT%INDEX = ITEM 
   ELSE 
      CALL INSERT(ITEM,ROOT%NEXT) 
   ENDIF 
END SUBROUTINE INSERT


RECURSIVE SUBROUTINE CUTCELL_INSERT(ITEM,ROOT) 
   IMPLICIT NONE 
   TYPE(CUTCELL_LINKED_LIST_TYPE), POINTER :: ROOT 
   INTEGER :: ITEM,IZERO
   IF (.NOT.ASSOCIATED(ROOT)) THEN 
      ALLOCATE(ROOT,STAT=IZERO)
      CALL ChkMemErr('GEOM','ROOT',IZERO)
      NULLIFY(ROOT%NEXT) 
      ROOT%INDEX = ITEM
      ROOT%AREA = 1._EB ! Charles Luo call here, cutcell AREA of FACET(N) with cell index ITEM
   ELSE 
      CALL CUTCELL_INSERT(ITEM,ROOT%NEXT) 
   ENDIF 
END SUBROUTINE CUTCELL_INSERT


SUBROUTINE TRIANGLE_BOX_INTERSECT(IERR,V1,V2,V3,BB)
IMPLICIT NONE

INTEGER, INTENT(OUT) :: IERR
REAL(EB), INTENT(IN) :: V1(3),V2(3),V3(3),BB(6)
REAL(EB) :: PLANE(4),P0(3),P1(3)

IERR=0

!! Filter small triangles
!
!A_TRI = TRIANGLE_AREA(V1,V2,V3)
!A_BB  = MIN( (BB(2)-BB(1))*(BB(4)-BB(3)), (BB(2)-BB(1))*(BB(6)-BB(5)), (BB(4)-BB(3))*(BB(6)-BB(5)) )
!IF (A_TRI < 0.01*A_BB) RETURN

! Are vertices outside of bounding planes?

IF (MAX(V1(1),V2(1),V3(1))<BB(1)) RETURN
IF (MIN(V1(1),V2(1),V3(1))>BB(2)) RETURN
IF (MAX(V1(2),V2(2),V3(2))<BB(3)) RETURN
IF (MIN(V1(2),V2(2),V3(2))>BB(4)) RETURN
IF (MAX(V1(3),V2(3),V3(3))<BB(5)) RETURN
IF (MIN(V1(3),V2(3),V3(3))>BB(6)) RETURN

! Any vertices inside bounding box?

IF ( V1(1)>=BB(1).AND.V1(1)<=BB(2) .AND. &
     V1(2)>=BB(3).AND.V1(2)<=BB(4) .AND. &
     V1(3)>=BB(5).AND.V1(3)<=BB(6) ) THEN
   IERR=1
   RETURN
ENDIF
IF ( V2(1)>=BB(1).AND.V2(1)<=BB(2) .AND. &
     V2(2)>=BB(3).AND.V2(2)<=BB(4) .AND. &
     V2(3)>=BB(5).AND.V2(3)<=BB(6) ) THEN
   IERR=1
   RETURN
ENDIF
IF ( V3(1)>=BB(1).AND.V3(1)<=BB(2) .AND. &
     V3(2)>=BB(3).AND.V3(2)<=BB(4) .AND. &
     V3(3)>=BB(5).AND.V3(3)<=BB(6) ) THEN
   IERR=1
   RETURN
ENDIF

! There are a couple other trivial rejection tests we could employ.
! But for now we jump straight to line segment--plane intersection.

! Test edge V1,V2 for intersection with each face of box
PLANE = (/-1._EB,0._EB,0._EB, BB(1)/); CALL LINE_PLANE_INTERSECT(IERR,V1,V2,PLANE,BB,-1); IF (IERR==1) RETURN
PLANE = (/ 1._EB,0._EB,0._EB,-BB(2)/); CALL LINE_PLANE_INTERSECT(IERR,V1,V2,PLANE,BB, 1); IF (IERR==1) RETURN
PLANE = (/0._EB,-1._EB,0._EB, BB(3)/); CALL LINE_PLANE_INTERSECT(IERR,V1,V2,PLANE,BB,-2); IF (IERR==1) RETURN
PLANE = (/0._EB, 1._EB,0._EB,-BB(4)/); CALL LINE_PLANE_INTERSECT(IERR,V1,V2,PLANE,BB, 2); IF (IERR==1) RETURN
PLANE = (/0._EB,0._EB,-1._EB, BB(5)/); CALL LINE_PLANE_INTERSECT(IERR,V1,V2,PLANE,BB,-3); IF (IERR==1) RETURN
PLANE = (/0._EB,0._EB, 1._EB,-BB(6)/); CALL LINE_PLANE_INTERSECT(IERR,V1,V2,PLANE,BB, 3); IF (IERR==1) RETURN

! Test edge V2,V3 for intersection with each face of box
PLANE = (/-1._EB,0._EB,0._EB, BB(1)/); CALL LINE_PLANE_INTERSECT(IERR,V2,V3,PLANE,BB,-1); IF (IERR==1) RETURN
PLANE = (/ 1._EB,0._EB,0._EB,-BB(2)/); CALL LINE_PLANE_INTERSECT(IERR,V2,V3,PLANE,BB, 1); IF (IERR==1) RETURN
PLANE = (/0._EB,-1._EB,0._EB, BB(3)/); CALL LINE_PLANE_INTERSECT(IERR,V2,V3,PLANE,BB,-2); IF (IERR==1) RETURN
PLANE = (/0._EB, 1._EB,0._EB,-BB(4)/); CALL LINE_PLANE_INTERSECT(IERR,V2,V3,PLANE,BB, 2); IF (IERR==1) RETURN
PLANE = (/0._EB,0._EB,-1._EB, BB(5)/); CALL LINE_PLANE_INTERSECT(IERR,V2,V3,PLANE,BB,-3); IF (IERR==1) RETURN
PLANE = (/0._EB,0._EB, 1._EB,-BB(6)/); CALL LINE_PLANE_INTERSECT(IERR,V2,V3,PLANE,BB, 3); IF (IERR==1) RETURN

! Test edge V3,V1 for intersection with each face of box
PLANE = (/-1._EB,0._EB,0._EB, BB(1)/); CALL LINE_PLANE_INTERSECT(IERR,V3,V1,PLANE,BB,-1); IF (IERR==1) RETURN
PLANE = (/ 1._EB,0._EB,0._EB,-BB(2)/); CALL LINE_PLANE_INTERSECT(IERR,V3,V1,PLANE,BB, 1); IF (IERR==1) RETURN
PLANE = (/0._EB,-1._EB,0._EB, BB(3)/); CALL LINE_PLANE_INTERSECT(IERR,V3,V1,PLANE,BB,-2); IF (IERR==1) RETURN
PLANE = (/0._EB, 1._EB,0._EB,-BB(4)/); CALL LINE_PLANE_INTERSECT(IERR,V3,V1,PLANE,BB, 2); IF (IERR==1) RETURN
PLANE = (/0._EB,0._EB,-1._EB, BB(5)/); CALL LINE_PLANE_INTERSECT(IERR,V3,V1,PLANE,BB,-3); IF (IERR==1) RETURN
PLANE = (/0._EB,0._EB, 1._EB,-BB(6)/); CALL LINE_PLANE_INTERSECT(IERR,V3,V1,PLANE,BB, 3); IF (IERR==1) RETURN

! The remaining possibility for tri-box intersection is that the corner of the box pokes through
! the triangle such that neither the vertices nor the edges of tri intersect any of the box faces.
! In this case the diagonal of the box corner intersects the plane formed by the tri.  The diagonal
! is defined as the line segment from point P0 to P1, formed from the corners of the bounding box.

! Test the four box diagonals:

P0 = (/BB(1),BB(3),BB(5)/)
P1 = (/BB(2),BB(4),BB(6)/)
CALL LINE_SEGMENT_TRIANGLE_INTERSECT(IERR,V1,V2,V3,P0,P1); IF (IERR==1) RETURN

P0 = (/BB(2),BB(3),BB(5)/)
P1 = (/BB(1),BB(4),BB(6)/)
CALL LINE_SEGMENT_TRIANGLE_INTERSECT(IERR,V1,V2,V3,P0,P1); IF (IERR==1) RETURN

P0 = (/BB(1),BB(3),BB(6)/)
P1 = (/BB(2),BB(4),BB(5)/)
CALL LINE_SEGMENT_TRIANGLE_INTERSECT(IERR,V1,V2,V3,P0,P1); IF (IERR==1) RETURN

P0 = (/BB(1),BB(4),BB(5)/)
P1 = (/BB(2),BB(3),BB(6)/)
CALL LINE_SEGMENT_TRIANGLE_INTERSECT(IERR,V1,V2,V3,P0,P1); IF (IERR==1) RETURN

! test commit from Charles Luo

END SUBROUTINE TRIANGLE_BOX_INTERSECT


REAL(EB) FUNCTION TRIANGLE_AREA(V1,V2,V3)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT,NORM2
IMPLICIT NONE

REAL(EB), INTENT(IN) :: V1(3),V2(3),V3(3)
REAL(EB) :: N(3),R1(3),R2(3)

R1 = V2-V1
R2 = V3-V1
CALL CROSS_PRODUCT(N,R1,R2)

TRIANGLE_AREA = 0.5_EB*NORM2(N)

END FUNCTION TRIANGLE_AREA


SUBROUTINE LINE_SEGMENT_TRIANGLE_INTERSECT(IERR,V1,V2,V3,P0,P1)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
IMPLICIT NONE

INTEGER, INTENT(OUT) :: IERR
REAL(EB), INTENT(IN) :: V1(3),V2(3),V3(3),P0(3),P1(3)
REAL(EB) :: E1(3),E2(3),S(3),Q(3),U,V,TMP,T,D(3),P(3)
REAL(EB), PARAMETER :: EPS=1.E-10_EB

IERR=0

! Schneider and Eberly, Section 11.1

D = P1-P0

E1 = V2-V1
E2 = V3-V1

CALL CROSS_PRODUCT(P,D,E2)

TMP = DOT_PRODUCT(P,E1)

IF ( ABS(TMP)<EPS ) RETURN

TMP = 1._EB/TMP
S = P0-V1

U = TMP*DOT_PRODUCT(S,P)
IF (U<0._EB .OR. U>1._EB) RETURN

CALL CROSS_PRODUCT(Q,S,E1)
V = TMP*DOT_PRODUCT(D,Q)
IF (V<0._EB .OR. (U+V)>1._EB) RETURN

T = TMP*DOT_PRODUCT(E2,Q)
!XI = P0 + T*D ! the intersection point

IF (T>=0._EB .AND. T<=1._EB) IERR=1

END SUBROUTINE LINE_SEGMENT_TRIANGLE_INTERSECT


SUBROUTINE LINE_PLANE_INTERSECT(IERR,P0,P1,PP,BB,IOR)
USE MATH_FUNCTIONS, ONLY: NORM2
IMPLICIT NONE

INTEGER, INTENT(OUT) :: IERR
REAL(EB), INTENT(IN) :: P0(3),P1(3),PP(4),BB(6)
INTEGER, INTENT(IN) :: IOR
REAL(EB) :: D(3),T,DENOM, Q0(3)
REAL(EB), PARAMETER :: EPS=1.E-10_EB

IERR=0
Q0=-999._EB
T=0._EB

D = P1-P0
DENOM = DOT_PRODUCT(PP(1:3),D)

IF (ABS(DENOM)>EPS) THEN
   T = -( DOT_PRODUCT(PP(1:3),P0)+PP(4) )/DENOM
   IF (T>=0._EB .AND. T<=1._EB) THEN
      Q0 = P0 + T*D ! instersection point
      IF (POINT_IN_BOX_2D(Q0,BB,IOR)) IERR=1
   ENDIF
ENDIF

END SUBROUTINE LINE_PLANE_INTERSECT


LOGICAL FUNCTION POINT_IN_BOX_2D(P,BB,IOR)
IMPLICIT NONE

REAL(EB), INTENT(IN) :: P(3),BB(6)
INTEGER, INTENT(IN) :: IOR

POINT_IN_BOX_2D=.FALSE.

SELECT CASE(ABS(IOR))
   CASE(1) ! YZ plane
      IF ( P(2)>=BB(3).AND.P(2)<=BB(4) .AND. &
           P(3)>=BB(5).AND.P(3)<=BB(6) ) POINT_IN_BOX_2D=.TRUE.
   CASE(2) ! XZ plane
      IF ( P(1)>=BB(1).AND.P(1)<=BB(2) .AND. &
           P(3)>=BB(5).AND.P(3)<=BB(6) ) POINT_IN_BOX_2D=.TRUE.
   CASE(3) ! XY plane
      IF ( P(1)>=BB(1).AND.P(1)<=BB(2) .AND. &
           P(2)>=BB(3).AND.P(2)<=BB(4) ) POINT_IN_BOX_2D=.TRUE.
END SELECT

END FUNCTION POINT_IN_BOX_2D


LOGICAL FUNCTION POINT_IN_TETRAHEDRON(XP,V1,V2,V3,V4,BB)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
IMPLICIT NONE

REAL(EB), INTENT(IN) :: XP(3),V1(3),V2(3),V3(3),V4(3),BB(6)
REAL(EB) :: U_VEC(3),V_VEC(3),N_VEC(3),Q_VEC(3),R_VEC(3)
INTEGER :: I

! In this routine, we test all four faces of the tet volume defined by the points X(i),Y(i),Z(i); i=1:4.
! If the point is on the negative side of all the faces, it is inside the volume.

POINT_IN_TETRAHEDRON=.FALSE.

! first test bounding box

IF (XP(1)<BB(1)) RETURN
IF (XP(1)>BB(2)) RETURN
IF (XP(2)<BB(3)) RETURN
IF (XP(2)>BB(4)) RETURN
IF (XP(3)<BB(5)) RETURN
IF (XP(3)>BB(6)) RETURN

POINT_IN_TETRAHEDRON=.TRUE.

FACE_LOOP: DO I=1,4

   SELECT CASE(I)
      CASE(1)
         ! vertex ordering = 1,2,3,4
         Q_VEC = XP-(/V1(1),V1(2),V1(3)/) ! form a vector from a point on the triangular surface to the point XP
         R_VEC = (/V4(1),V4(2),V4(3)/)-(/V1(1),V1(2),V1(3)/) ! vector from the tri to other point of volume defining inside
         U_VEC = (/V2(1)-V1(1),V2(2)-V1(2),V2(3)-V1(3)/) ! vectors forming the sides of the triangle
         V_VEC = (/V3(1)-V1(1),V3(2)-V1(2),V3(3)-V1(3)/)
      CASE(2)
         ! vertex ordering = 1,3,4,2
         Q_VEC = XP-(/V1(1),V1(2),V1(3)/)
         R_VEC = (/V2(1),V2(2),V2(3)/)-(/V1(1),V1(2),V1(3)/)
         U_VEC = (/V3(1)-V1(1),V3(2)-V1(2),V3(3)-V1(3)/)
         V_VEC = (/V4(1)-V1(1),V4(2)-V1(2),V4(3)-V1(3)/)
      CASE(3)
         ! vertex ordering = 1,4,2,3
         Q_VEC = XP-(/V1(1),V1(2),V1(3)/)
         R_VEC = (/V2(1),V2(2),V2(3)/)-(/V1(1),V1(2),V1(3)/)
         U_VEC = (/V4(1)-V1(1),V4(2)-V1(2),V4(3)-V1(3)/)
         V_VEC = (/V2(1)-V1(1),V2(2)-V1(2),V2(3)-V1(3)/)
      CASE(4)
         ! vertex ordering = 2,4,3,1
         Q_VEC = XP-(/V2(1),V2(2),V2(3)/)
         R_VEC = (/V1(1),V1(2),V1(3)/)-(/V2(1),V2(2),V2(3)/)
         U_VEC = (/V4(1)-V2(1),V4(2)-V2(2),V4(3)-V2(3)/)
         V_VEC = (/V3(1)-V2(1),V3(2)-V2(2),V3(3)-V2(3)/)
   END SELECT

   ! if the sign of the dot products are equal, the point is inside, else it is outside and we return

   IF ( ABS( SIGN(1._EB,DOT_PRODUCT(Q_VEC,N_VEC))-SIGN(1._EB,DOT_PRODUCT(R_VEC,N_VEC)) )>ZERO_P ) THEN
      POINT_IN_TETRAHEDRON=.FALSE.
      RETURN
   ENDIF

ENDDO FACE_LOOP

END FUNCTION POINT_IN_TETRAHEDRON


LOGICAL FUNCTION POINT_IN_POLYHEDRON(XP,BB)
IMPLICIT NONE

REAL(EB) :: XP(3),BB(6),XX(3),YY(3),ZZ(3),RAY_DIRECTION(3)
INTEGER :: I,J,N_INTERSECTIONS,IRAY
REAL(EB), PARAMETER :: EPS=1.E-6_EB

! Schneider and Eberly, Geometric Tools for Computer Graphics, Morgan Kaufmann, 2003. Section 13.4

POINT_IN_POLYHEDRON=.FALSE.

! test global bounding box

IF ( XP(1)<BB(1) .OR. XP(1)>BB(2) ) RETURN
IF ( XP(2)<BB(3) .OR. XP(2)>BB(4) ) RETURN
IF ( XP(3)<BB(5) .OR. XP(3)>BB(6) ) RETURN

N_INTERSECTIONS=0

RAY_DIRECTION = (/0._EB,0._EB,1._EB/)

FACE_LOOP: DO I=1,N_FACE

   ! test bounding box
   XX(1) = VERTEX(FACET(I)%VERTEX(1))%X
   XX(2) = VERTEX(FACET(I)%VERTEX(2))%X
   XX(3) = VERTEX(FACET(I)%VERTEX(3))%X

   IF (XP(1)<MINVAL(XX)) CYCLE FACE_LOOP
   IF (XP(1)>MAXVAL(XX)) CYCLE FACE_LOOP

   YY(1) = VERTEX(FACET(I)%VERTEX(1))%Y
   YY(2) = VERTEX(FACET(I)%VERTEX(2))%Y
   YY(3) = VERTEX(FACET(I)%VERTEX(3))%Y

   IF (XP(2)<MINVAL(YY)) CYCLE FACE_LOOP
   IF (XP(2)>MAXVAL(YY)) CYCLE FACE_LOOP

   ZZ(1) = VERTEX(FACET(I)%VERTEX(1))%Z
   ZZ(2) = VERTEX(FACET(I)%VERTEX(2))%Z
   ZZ(3) = VERTEX(FACET(I)%VERTEX(3))%Z

   IF (XP(3)>MAXVAL(ZZ)) CYCLE FACE_LOOP

   RAY_TEST_LOOP: DO J=1,3
      IRAY = RAY_TRIANGLE_INTERSECT(I,XP,RAY_DIRECTION)
      SELECT CASE(IRAY)
         CASE(0)
            ! does not intersect
            EXIT RAY_TEST_LOOP
         CASE(1)
            ! ray intersects triangle
            N_INTERSECTIONS=N_INTERSECTIONS+1
            EXIT RAY_TEST_LOOP
         CASE(2)
            ! ray intersects edge, try new ray (shift origin)
            IF (J==1) XP=XP+(/EPS,0._EB,0._EB/) ! shift in x direction
            IF (J==2) XP=XP+(/0._EB,EPS,0._EB/) ! shift in y direction
            IF (J==3) WRITE(LU_ERR,*) 'WARNING: ray test failed'
      END SELECT
   ENDDO RAY_TEST_LOOP

ENDDO FACE_LOOP

IF ( MOD(N_INTERSECTIONS,2)/=0 ) POINT_IN_POLYHEDRON=.TRUE.

END FUNCTION POINT_IN_POLYHEDRON


LOGICAL FUNCTION POINT_IN_TRIANGLE(P,V1,V2,V3)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
IMPLICIT NONE

REAL(EB), INTENT(IN) :: P(3),V1(3),V2(3),V3(3)
REAL(EB) :: E(3),E1(3),E2(3),N(3),R(3),Q(3)
INTEGER :: I
REAL(EB), PARAMETER :: EPS=1.E-10_EB

! This routine tests whether the projection of P, in the plane normal
! direction, onto to the plane defined by the triangle (V1,V2,V3) is
! inside the triangle.

POINT_IN_TRIANGLE=.TRUE. ! start by assuming the point is inside

! compute face normal
E1 = V2-V1
E2 = V3-V1
CALL CROSS_PRODUCT(N,E1,E2)

EDGE_LOOP: DO I=1,3
   SELECT CASE(I)
      CASE(1)
         E = V2-V1
         R = P-V1
      CASE(2)
         E = V3-V2
         R = P-V2
      CASE(3)
         E = V1-V3
         R = P-V3
   END SELECT
   CALL CROSS_PRODUCT(Q,E,R)
   IF ( DOT_PRODUCT(Q,N) < -EPS ) THEN
      POINT_IN_TRIANGLE=.FALSE.
      RETURN
   ENDIF
ENDDO EDGE_LOOP

END FUNCTION POINT_IN_TRIANGLE


INTEGER FUNCTION RAY_TRIANGLE_INTERSECT(TRI,XP,D)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
IMPLICIT NONE

INTEGER, INTENT(IN) :: TRI
REAL(EB), INTENT(IN) :: XP(3),D(3)
REAL(EB) :: E1(3),E2(3),P(3),S(3),Q(3),U,V,TMP,V1(3),V2(3),V3(3),T !,XI(3)
REAL(EB), PARAMETER :: EPS=1.E-10_EB

! Schneider and Eberly, Section 11.1

V1(1) = VERTEX(FACET(TRI)%VERTEX(1))%X
V1(2) = VERTEX(FACET(TRI)%VERTEX(1))%Y
V1(3) = VERTEX(FACET(TRI)%VERTEX(1))%Z

V2(1) = VERTEX(FACET(TRI)%VERTEX(2))%X
V2(2) = VERTEX(FACET(TRI)%VERTEX(2))%Y
V2(3) = VERTEX(FACET(TRI)%VERTEX(2))%Z

V3(1) = VERTEX(FACET(TRI)%VERTEX(3))%X
V3(2) = VERTEX(FACET(TRI)%VERTEX(3))%Y
V3(3) = VERTEX(FACET(TRI)%VERTEX(3))%Z

E1 = V2-V1
E2 = V3-V1

CALL CROSS_PRODUCT(P,D,E2)

TMP = DOT_PRODUCT(P,E1)

IF ( ABS(TMP)<EPS ) THEN
   RAY_TRIANGLE_INTERSECT=0
   RETURN
ENDIF

TMP = 1._EB/TMP
S = XP-V1

U = TMP*DOT_PRODUCT(S,P)
IF (U<-EPS .OR. U>(1._EB+EPS)) THEN
   ! ray does not intersect triangle
   RAY_TRIANGLE_INTERSECT=0
   RETURN
ENDIF

IF (U<EPS .OR. U>(1._EB-EPS)) THEN
   ! ray intersects edge
   RAY_TRIANGLE_INTERSECT=2
   RETURN
ENDIF

CALL CROSS_PRODUCT(Q,S,E1)
V = TMP*DOT_PRODUCT(D,Q)
IF (V<-EPS .OR. (U+V)>(1._EB+EPS)) THEN
   ! ray does not intersect triangle
   RAY_TRIANGLE_INTERSECT=0
   RETURN
ENDIF

IF (V<EPS .OR. (U+V)>(1._EB-EPS)) THEN
   ! ray intersects edge
   RAY_TRIANGLE_INTERSECT=2
   RETURN
ENDIF

T = TMP*DOT_PRODUCT(E2,Q)
!XI = XP + T*D ! the intersection point

IF (T>0._EB) THEN
   RAY_TRIANGLE_INTERSECT=1
ELSE
   RAY_TRIANGLE_INTERSECT=0
ENDIF
RETURN

END FUNCTION RAY_TRIANGLE_INTERSECT


REAL(EB) FUNCTION TRILINEAR(UU,DXI,LL)
IMPLICIT NONE

REAL(EB), INTENT(IN) :: UU(0:1,0:1,0:1),DXI(3),LL(3)

! Comments:
!
! see http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/index.html
! with appropriate scaling. LL is length of side.
!
!                       UU(1,1,1)
!        z /----------/
!        ^/          / |
!        ------------  |    Particle position
!        |          |  |
!  LL(3) |   o<-----|------- DXI = [DXI(1),DXI(2),DXI(3)]
!        |          | /        
!        |          |/      Particle property at XX = TRILINEAR
!        ------------> x
!        ^
!        |
!   X0 = [0,0,0]
!
!    UU(0,0,0)
!
!===========================================================

TRILINEAR = UU(0,0,0)*(LL(1)-DXI(1))*(LL(2)-DXI(2))*(LL(3)-DXI(3)) +    &
            UU(1,0,0)*DXI(1)*(LL(2)-DXI(2))*(LL(3)-DXI(3)) +            &
            UU(0,1,0)*(LL(1)-DXI(1))*DXI(2)*(LL(3)-DXI(3)) +            &
            UU(0,0,1)*(LL(1)-DXI(1))*(LL(2)-DXI(2))*DXI(3) +            &
            UU(1,0,1)*DXI(1)*(LL(2)-DXI(2))*DXI(3) +                    &
            UU(0,1,1)*(LL(1)-DXI(1))*DXI(2)*DXI(3) +                    &
            UU(1,1,0)*DXI(1)*DXI(2)*(LL(3)-DXI(3)) +                    &
            UU(1,1,1)*DXI(1)*DXI(2)*DXI(3)

TRILINEAR = TRILINEAR/(LL(1)*LL(2)*LL(3))

END FUNCTION TRILINEAR


SUBROUTINE GETX(XI,XU,NG)
IMPLICIT NONE

REAL(EB), INTENT(OUT) :: XI(3)
REAL(EB), INTENT(IN) :: XU(3)
INTEGER, INTENT(IN) :: NG
TYPE(GEOMETRY_TYPE), POINTER :: G
REAL(EB) :: PP(3),RU(3),RUMAG,DR,DP

G => GEOMETRY(NG)
SELECT CASE(G%ISHAPE)
   CASE(IBOX)
      XI = XU
      IF (XU(1)<G%X1) XI(1) = XU(1) + (XU(1)-G%X1)
      IF (XU(1)>G%X2) XI(1) = XU(1) + (XU(1)-G%X2)
      IF (XU(2)<G%Y1) XI(2) = XU(2) + (XU(2)-G%Y1)
      IF (XU(2)>G%Y2) XI(2) = XU(2) + (XU(2)-G%Y2)
      IF (XU(3)<G%Z1) XI(3) = XU(3) + (XU(3)-G%Z1)
      IF (XU(3)>G%Z2) XI(3) = XU(3) + (XU(3)-G%Z2)
   CASE(ISPHERE)
      RU     = XU-(/G%X,G%Y,G%Z/)
      RUMAG  = SQRT(DOT_PRODUCT(RU,RU))
      DR     = RUMAG-G%RADIUS
      XI     = XU + DR*RU/RUMAG
   CASE(ICYLINDER)
      ! at the moment, cylinder must be aligned with an axis
      IF (ABS(G%XOR-1._EB)<EPSILON_EB) RU = (/0._EB,XU(2),XU(3)/)-(/0._EB,G%Y,G%Z/)
      IF (ABS(G%YOR-1._EB)<EPSILON_EB) RU = (/XU(1),0._EB,XU(3)/)-(/G%X,0._EB,G%Z/)
      IF (ABS(G%ZOR-1._EB)<EPSILON_EB) RU = (/XU(1),XU(2),0._EB/)-(/G%X,G%Y,0._EB/)
      RUMAG  = SQRT(DOT_PRODUCT(RU,RU))
      DR     = RUMAG-G%RADIUS
      XI     = XU + DR*RU/RUMAG
   CASE(IPLANE)
      PP = (/G%X,G%Y,G%Z/)           ! point in the plane
      DP = DOT_PRODUCT(G%NN,XU-PP)   ! signed distance to plane
      XI = XU + DP*G%NN
END SELECT

END SUBROUTINE GETX


SUBROUTINE GETU(U_DATA,DXI,XI_IN,I_VEL,NM)
IMPLICIT NONE

REAL(EB), INTENT(OUT) :: U_DATA(0:1,0:1,0:1),DXI(3)
REAL(EB), INTENT(IN) :: XI_IN(3)
INTEGER, INTENT(IN) :: I_VEL,NM
TYPE(MESH_TYPE), POINTER :: M
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
INTEGER :: II,JJ,KK
!CHARACTER(100) :: MESSAGE
REAL(EB) :: XI(3)

M=>MESHES(NM)
IF (PREDICTOR) THEN
   UU => M%U
   VV => M%V
   WW => M%W
ELSE
   UU => M%US
   VV => M%VS
   WW => M%WS
ENDIF

!II = INDU(1)
!JJ = INDU(2)
!KK = INDU(3)
!
!IF (XI(1)<XU(1)) THEN
!   N=CEILING((XU(1)-XI(1))/M%DX(II))
!   II=MAX(0,II-N)
!   DXI(1)=XI(1)-(XU(1)-REAL(N,EB)*M%DX(II))
!ELSE
!   N=FLOOR((XI(1)-XU(1))/M%DX(II))
!   II=MIN(IBP1,II+N)
!   DXI(1)=XI(1)-(XU(1)+REAL(N,EB)*M%DX(II))
!ENDIF
!
!IF (XI(2)<XU(2)) THEN
!   N=CEILING((XU(2)-XI(2))/M%DY(JJ))
!   JJ=MAX(0,JJ-N)
!   DXI(2)=XI(2)-(XU(2)-REAL(N,EB)*M%DY(JJ))
!ELSE
!   N=FLOOR((XI(2)-XU(2))/M%DY(JJ))
!   JJ=MIN(JBP1,JJ+N)
!   DXI(2)=XI(2)-(XU(2)+REAL(N,EB)*M%DY(JJ))
!ENDIF
!
!IF (XI(3)<XU(3)) THEN
!   N=CEILING((XU(3)-XI(3))/M%DZ(KK))
!   KK=MAX(0,KK-N)
!   DXI(3)=XI(3)-(XU(3)-REAL(N,EB)*M%DZ(KK))
!ELSE
!   N=FLOOR((XI(3)-XU(3))/M%DZ(KK))
!   KK=MIN(KBP1,KK+N)
!   DXI(3)=XI(3)-(XU(3)+REAL(N,EB)*M%DZ(KK))
!ENDIF

XI(1) = MAX(M%XS,MIN(M%XF,XI_IN(1)))
XI(2) = MAX(M%YS,MIN(M%YF,XI_IN(2)))
XI(3) = MAX(M%ZS,MIN(M%ZF,XI_IN(3)))

SELECT CASE(I_VEL)
   CASE(1)
      II = FLOOR((XI(1)-M%XS)/M%DX(1))
      JJ = FLOOR((XI(2)-M%YS)/M%DY(1)+0.5_EB)
      KK = FLOOR((XI(3)-M%ZS)/M%DZ(1)+0.5_EB)
      DXI(1) = XI(1) - M%X(II)
      DXI(2) = XI(2) - M%YC(JJ)
      DXI(3) = XI(3) - M%ZC(KK)
   CASE(2)
      II = FLOOR((XI(1)-M%XS)/M%DX(1)+0.5_EB)
      JJ = FLOOR((XI(2)-M%YS)/M%DY(1))
      KK = FLOOR((XI(3)-M%ZS)/M%DZ(1)+0.5_EB)
      DXI(1) = XI(1) - M%XC(II)
      DXI(2) = XI(2) - M%Y(JJ)
      DXI(3) = XI(3) - M%ZC(KK)
   CASE(3)
      II = FLOOR((XI(1)-M%XS)/M%DX(1)+0.5_EB)
      JJ = FLOOR((XI(2)-M%YS)/M%DY(1)+0.5_EB)
      KK = FLOOR((XI(3)-M%ZS)/M%DZ(1))
      DXI(1) = XI(1) - M%XC(II)
      DXI(2) = XI(2) - M%YC(JJ)
      DXI(3) = XI(3) - M%Z(KK)
   CASE(4)
      II = FLOOR((XI(1)-M%XS)/M%DX(1)+0.5_EB)
      JJ = FLOOR((XI(2)-M%YS)/M%DY(1)+0.5_EB)
      KK = FLOOR((XI(3)-M%ZS)/M%DZ(1)+0.5_EB)
      DXI(1) = XI(1) - M%XC(II)
      DXI(2) = XI(2) - M%YC(JJ)
      DXI(3) = XI(3) - M%ZC(KK)
END SELECT

DXI = MAX(0._EB,DXI)

!IF (ANY(DXI<0._EB)) THEN
!   WRITE(MESSAGE,'(A)') 'ERROR: DXI<0 in GETU'
!   CALL SHUTDOWN(MESSAGE)
!ENDIF
!IF (DXI(1)>M%DX(II) .OR. DXI(2)>M%DY(JJ) .OR. DXI(3)>M%DZ(KK)) THEN
!   WRITE(MESSAGE,'(A)') 'ERROR: DXI>DX in GETU'
!   CALL SHUTDOWN(MESSAGE)
!ENDIF

SELECT CASE(I_VEL)
   CASE(1)
      U_DATA(0,0,0) = UU(II,JJ,KK)
      U_DATA(1,0,0) = UU(II+1,JJ,KK)
      U_DATA(0,1,0) = UU(II,JJ+1,KK)
      U_DATA(0,0,1) = UU(II,JJ,KK+1)
      U_DATA(1,0,1) = UU(II+1,JJ,KK+1)
      U_DATA(0,1,1) = UU(II,JJ+1,KK+1)
      U_DATA(1,1,0) = UU(II+1,JJ+1,KK)
      U_DATA(1,1,1) = UU(II+1,JJ+1,KK+1)
   CASE(2)
      U_DATA(0,0,0) = VV(II,JJ,KK)
      U_DATA(1,0,0) = VV(II+1,JJ,KK)
      U_DATA(0,1,0) = VV(II,JJ+1,KK)
      U_DATA(0,0,1) = VV(II,JJ,KK+1)
      U_DATA(1,0,1) = VV(II+1,JJ,KK+1)
      U_DATA(0,1,1) = VV(II,JJ+1,KK+1)
      U_DATA(1,1,0) = VV(II+1,JJ+1,KK)
      U_DATA(1,1,1) = VV(II+1,JJ+1,KK+1)
   CASE(3)
      U_DATA(0,0,0) = WW(II,JJ,KK)
      U_DATA(1,0,0) = WW(II+1,JJ,KK)
      U_DATA(0,1,0) = WW(II,JJ+1,KK)
      U_DATA(0,0,1) = WW(II,JJ,KK+1)
      U_DATA(1,0,1) = WW(II+1,JJ,KK+1)
      U_DATA(0,1,1) = WW(II,JJ+1,KK+1)
      U_DATA(1,1,0) = WW(II+1,JJ+1,KK)
      U_DATA(1,1,1) = WW(II+1,JJ+1,KK+1)
   CASE(4) ! viscosity
      U_DATA(0,0,0) = M%MU(II,JJ,KK)
      U_DATA(1,0,0) = M%MU(II+1,JJ,KK)
      U_DATA(0,1,0) = M%MU(II,JJ+1,KK)
      U_DATA(0,0,1) = M%MU(II,JJ,KK+1)
      U_DATA(1,0,1) = M%MU(II+1,JJ,KK+1)
      U_DATA(0,1,1) = M%MU(II,JJ+1,KK+1)
      U_DATA(1,1,0) = M%MU(II+1,JJ+1,KK)
      U_DATA(1,1,1) = M%MU(II+1,JJ+1,KK+1)
END SELECT

END SUBROUTINE GETU


SUBROUTINE GETGRAD(G_DATA,DXI,XI,XU,INDU,COMP_I,COMP_J,NM)
IMPLICIT NONE

REAL(EB), INTENT(OUT) :: G_DATA(0:1,0:1,0:1),DXI(3)
REAL(EB), INTENT(IN) :: XI(3),XU(3)
INTEGER, INTENT(IN) :: INDU(3),COMP_I,COMP_J,NM
TYPE(MESH_TYPE), POINTER :: M
REAL(EB), POINTER, DIMENSION(:,:,:) :: DUDX
INTEGER :: II,JJ,KK,N
CHARACTER(100) :: MESSAGE

M=>MESHES(NM)

IF (COMP_I==1 .AND. COMP_J==1) DUDX => M%WORK5
IF (COMP_I==1 .AND. COMP_J==2) DUDX => M%IBM_SAVE1
IF (COMP_I==1 .AND. COMP_J==3) DUDX => M%IBM_SAVE2
IF (COMP_I==2 .AND. COMP_J==1) DUDX => M%IBM_SAVE3
IF (COMP_I==2 .AND. COMP_J==2) DUDX => M%WORK6
IF (COMP_I==2 .AND. COMP_J==3) DUDX => M%IBM_SAVE4
IF (COMP_I==3 .AND. COMP_J==1) DUDX => M%IBM_SAVE5
IF (COMP_I==3 .AND. COMP_J==2) DUDX => M%IBM_SAVE6
IF (COMP_I==3 .AND. COMP_J==3) DUDX => M%WORK7

II = INDU(1)
JJ = INDU(2)
KK = INDU(3)

IF (XI(1)<XU(1)) THEN
   N=CEILING((XU(1)-XI(1))/M%DX(II))
   II=MAX(0,II-N)
   DXI(1)=XI(1)-(XU(1)-REAL(N,EB)*M%DX(II))
ELSE
   N=FLOOR((XI(1)-XU(1))/M%DX(II))
   II=MIN(IBP1,II+N)
   DXI(1)=XI(1)-(XU(1)+REAL(N,EB)*M%DX(II))
ENDIF

IF (XI(2)<XU(2)) THEN
   N=CEILING((XU(2)-XI(2))/M%DY(JJ))
   JJ=MAX(0,JJ-N)
   DXI(2)=XI(2)-(XU(2)-REAL(N,EB)*M%DY(JJ))
ELSE
   N=FLOOR((XI(2)-XU(2))/M%DY(JJ))
   JJ=MIN(JBP1,JJ+N)
   DXI(2)=XI(2)-(XU(2)+REAL(N,EB)*M%DY(JJ))
ENDIF

IF (XI(3)<XU(3)) THEN
   N=CEILING((XU(3)-XI(3))/M%DZ(KK))
   KK=MAX(0,KK-N)
   DXI(3)=XI(3)-(XU(3)-REAL(N,EB)*M%DZ(KK))
ELSE
   N=FLOOR((XI(3)-XU(3))/M%DZ(KK))
   KK=MIN(KBP1,KK+N)
   DXI(3)=XI(3)-(XU(3)+REAL(N,EB)*M%DZ(KK))
ENDIF

IF (ANY(DXI<0._EB)) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: DXI<0 in GETGRAD'
   CALL SHUTDOWN(MESSAGE)
ENDIF
IF (DXI(1)>M%DX(II) .OR. DXI(2)>M%DY(JJ) .OR. DXI(3)>M%DZ(KK)) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: DXI>DX in GETGRAD'
   CALL SHUTDOWN(MESSAGE)
ENDIF

G_DATA(0,0,0) = DUDX(II,JJ,KK)
G_DATA(1,0,0) = DUDX(II+1,JJ,KK)
G_DATA(0,1,0) = DUDX(II,JJ+1,KK)
G_DATA(0,0,1) = DUDX(II,JJ,KK+1)
G_DATA(1,0,1) = DUDX(II+1,JJ,KK+1)
G_DATA(0,1,1) = DUDX(II,JJ+1,KK+1)
G_DATA(1,1,0) = DUDX(II+1,JJ+1,KK)
G_DATA(1,1,1) = DUDX(II+1,JJ+1,KK+1)

END SUBROUTINE GETGRAD


SUBROUTINE GET_VELO_IBM(VELO_IBM,IERR,VELO_INDEX,XVELO,TRI_INDEX,IBM_INDEX,DXC,NM)
USE MATH_FUNCTIONS, ONLY: NORM2
IMPLICIT NONE

REAL(EB), INTENT(OUT) :: VELO_IBM
INTEGER, INTENT(OUT) :: IERR
REAL(EB), INTENT(IN) :: XVELO(3),DXC(3)
INTEGER, INTENT(IN) :: VELO_INDEX,TRI_INDEX,IBM_INDEX,NM
REAL(EB) :: N(3),R(3),V1(3),V2(3),V3(3),T,U_DATA(0:1,0:1,0:1),XI(3),DXI(3)
REAL(EB), PARAMETER :: EPS=1.E-10_EB

IERR=0
VELO_IBM=0._EB

V1 = (/VERTEX(FACET(TRI_INDEX)%VERTEX(1))%X,VERTEX(FACET(TRI_INDEX)%VERTEX(1))%Y,VERTEX(FACET(TRI_INDEX)%VERTEX(1))%Z/)
V2 = (/VERTEX(FACET(TRI_INDEX)%VERTEX(2))%X,VERTEX(FACET(TRI_INDEX)%VERTEX(2))%Y,VERTEX(FACET(TRI_INDEX)%VERTEX(2))%Z/)
V3 = (/VERTEX(FACET(TRI_INDEX)%VERTEX(3))%X,VERTEX(FACET(TRI_INDEX)%VERTEX(3))%Y,VERTEX(FACET(TRI_INDEX)%VERTEX(3))%Z/)
N = FACET(TRI_INDEX)%NVEC

R = XVELO-V1
IF ( NORM2(R)<EPS ) R = XVELO-V2 ! select a different vertex

T = DOT_PRODUCT(R,N)

IF (IBM_INDEX==0 .AND. T<EPS) RETURN ! the velocity point is on or interior to the surface

IF (IBM_INDEX==1) THEN
   XI = XVELO + T*N
   CALL GETU(U_DATA,DXI,XI,VELO_INDEX,NM)
   VELO_IBM = 0.5_EB*TRILINEAR(U_DATA,DXI,DXC)
   RETURN
ENDIF

IERR=1

END SUBROUTINE GET_VELO_IBM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Cut-cell subroutines by Charles Luo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TRI_PLANE_BOX_INTERSECT(NP,PC,V1,V2,V3,BB)
USE MATH_FUNCTIONS
IMPLICIT NONE
! get the intersection points (cooridnates) of the BB's 12 edges and the plane of the trianlge
! regular intersection polygons with 0, 3, 4, 5, or 6 corners
! irregular intersection case (corner, edge, or face intersection) should also be ok.

INTEGER, INTENT(OUT) :: NP
REAL(EB), INTENT(OUT) :: PC(18) ! max 6 points but maybe repeated at the vertices
REAL(EB), INTENT(IN) :: V1(3),V2(3),V3(3),BB(6)
REAL(EB) :: P0(3),P1(3),Q(3),PC_TMP(24)
INTEGER :: I,J,IERR,IERR2

NP = 0
EDGE_LOOP: DO I=1,12
    SELECT CASE(I) 
        CASE(1)
            P0(1)=BB(1)
            P0(2)=BB(3)
            P0(3)=BB(5)
            P1(1)=BB(2)
            P1(2)=BB(3)
            P1(3)=BB(5)
        CASE(2)
            P0(1)=BB(2)
            P0(2)=BB(3)
            P0(3)=BB(5)
            P1(1)=BB(2)
            P1(2)=BB(4)
            P1(3)=BB(5)
        CASE(3)
            P0(1)=BB(2)
            P0(2)=BB(4)
            P0(3)=BB(5)
            P1(1)=BB(1)
            P1(2)=BB(4)
            P1(3)=BB(5)
        CASE(4)
            P0(1)=BB(1)
            P0(2)=BB(4)
            P0(3)=BB(5)
            P1(1)=BB(1)
            P1(2)=BB(3)
            P1(3)=BB(5)
        CASE(5)
            P0(1)=BB(1)
            P0(2)=BB(3)
            P0(3)=BB(6)
            P1(1)=BB(2)
            P1(2)=BB(3)
            P1(3)=BB(6)
        CASE(6)
            P0(1)=BB(2)
            P0(2)=BB(3)
            P0(3)=BB(6)
            P1(1)=BB(2)
            P1(2)=BB(4)
            P1(3)=BB(6)
        CASE(7)
            P0(1)=BB(2)
            P0(2)=BB(4)
            P0(3)=BB(6)
            P1(1)=BB(1)
            P1(2)=BB(4)
            P1(3)=BB(6)
        CASE(8)
            P0(1)=BB(1)
            P0(2)=BB(4)
            P0(3)=BB(6)
            P1(1)=BB(1)
            P1(2)=BB(3)
            P1(3)=BB(6)
        CASE(9)
            P0(1)=BB(1)
            P0(2)=BB(3)
            P0(3)=BB(5)
            P1(1)=BB(1)
            P1(2)=BB(3)
            P1(3)=BB(6)
        CASE(10)
            P0(1)=BB(2)
            P0(2)=BB(3)
            P0(3)=BB(5)
            P1(1)=BB(2)
            P1(2)=BB(3)
            P1(3)=BB(6)
        CASE(11)
            P0(1)=BB(2)
            P0(2)=BB(4)
            P0(3)=BB(5)
            P1(1)=BB(2)
            P1(2)=BB(4)
            P1(3)=BB(6)
        CASE(12)
            P0(1)=BB(1)
            P0(2)=BB(4)
            P0(3)=BB(5)
            P1(1)=BB(1)
            P1(2)=BB(4)
            P1(3)=BB(6)
    END SELECT 
    CALL LINE_SEG_TRI_PLANE_INTERSECT(IERR,IERR2,Q,V1,V2,V3,P0,P1)
    
    IF (IERR==1) THEN
        NP=NP+1
        DO J=1,3
            PC_TMP((NP-1)*3+J)=Q(J)
        ENDDO
    ENDIF
ENDDO EDGE_LOOP

! For more than 3 intersection points
! they have to be sorted in order to create a convex polygon
CALL ELIMATE_REPEATED_POINTS(NP,PC_TMP)
DO I=1,NP*3
   PC(I) = PC_TMP(I)
ENDDO
IF (NP > 3) THEN 
    CALL SORT_POLYGON_CORNERS(NP,V1,V2,V3,PC)
ENDIF

RETURN
END SUBROUTINE TRI_PLANE_BOX_INTERSECT


SUBROUTINE SORT_POLYGON_CORNERS(NP,V1,V2,V3,PC)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
IMPLICIT NONE
! Sort all the corners of a polygon
! Ref: Gernot Hoffmann, Cube Plane Intersection.

INTEGER, INTENT(IN) :: NP
REAL(EB), INTENT(INOUT) :: PC(3*NP)
REAL(EB), INTENT(IN) :: V1(3),V2(3),V3(3)
REAL(EB) :: MEAN_VALUE(3),POLY_NORM(3),R1,R2,TMP(3),U(3),W(3)
INTEGER :: I,J,K,IOR,NA,NB

IF (NP <=3 ) RETURN

U = V2-V1
W = V3-V1
CALL CROSS_PRODUCT(POLY_NORM,U,W)

DO I=1,3
    MEAN_VALUE(I) = 0._EB
    DO J=1,NP
        MEAN_VALUE(I) = MEAN_VALUE(I) + PC((J-1)*3+I)/REAL(NP)
    ENDDO
ENDDO

!get normal of ploygan 
IF (ABS(POLY_NORM(1)) >= ABS(POLY_NORM(2)) .AND. ABS(POLY_NORM(1)) >= ABS(POLY_NORM(2)) ) THEN
    IOR = 1
    NA = 2
    NB = 3
ELSE IF (ABS(POLY_NORM(2)) >= ABS(POLY_NORM(3)) ) THEN
    IOR = 2
    NA = 1
    NB = 3
ELSE
    IOR = 3
    NA = 1
    NB = 2
ENDIF

DO I=1,NP-1
    R1 = ATAN2(PC((I-1)*3+NB)-MEAN_VALUE(NB), PC((I-1)*3+NA)-MEAN_VALUE(NA))
    DO J=I+1, NP
        R2 = ATAN2(PC((J-1)*3+NB)-MEAN_VALUE(NB), PC((J-1)*3+NA)-MEAN_VALUE(NA))
        IF (R2 < R1) THEN
            DO K=1,3
                TMP(K) = PC((J-1)*3+K)
                PC((J-1)*3+K) = PC((I-1)*3+K)
                PC((I-1)*3+K) = TMP(K)
                R1 = R2
            ENDDO
        ENDIF
    ENDDO
ENDDO
    
RETURN
END SUBROUTINE SORT_POLYGON_CORNERS


SUBROUTINE TRIANGLE_POLYGON_POINTS(IERR,NXP,XPC,V1,V2,V3,NP,PC,BB)
IMPLICIT NONE
! Calculate the intersection points of a triangle and a polygon, if intersected.
! http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm

INTEGER, INTENT(IN) :: NP
INTEGER, INTENT(OUT) :: NXP,IERR
REAL(EB), INTENT(OUT) :: XPC(27)
REAL(EB), INTENT(IN) :: V1(3),V2(3),V3(3),PC(18),BB(6)
INTEGER :: I,J,K
REAL(EB) :: U(3),V(3),W(3),S1P0(3),XC(3)
REAL(EB) :: A,B,C,D,E,DD,SC,TC
REAL(EB), PARAMETER :: EPS=1.E-10_EB
!LOGICAL :: POINT_IN_BB, POINT_IN_TRIANGLE

IERR = 0
SC = 0._EB
TC = 0._EB
NXP = 0
TRIANGLE_LOOP: DO I=1,3
    SELECT CASE(I)
    CASE(1)
        U = V2-V1
        S1P0 = V1
    CASE(2)
        U = V3-V2
        S1P0 = V2
    CASE(3)
        U = V1-V3
        S1P0 = V3
    END SELECT
    
    POLYGON_LOOP: DO J=1,NP
        IF (J < NP) THEN
            DO K=1,3
                V(K) = PC(J*3+K)-PC((J-1)*3+K)
            ENDDO
        ELSE IF (J == NP) THEN        
            DO K=1,3
                V(K) = PC(K)-PC((J-1)*3+K)
            ENDDO
        ENDIF
        
        DO K=1,3
            W(K) = S1P0(K)-PC((J-1)*3+K)
        ENDDO
        
        A = DOT_PRODUCT(U,U)
        B = DOT_PRODUCT(U,V)
        C = DOT_PRODUCT(V,V)
        D = DOT_PRODUCT(U,W)
        E = DOT_PRODUCT(V,W)
        DD = A*C-B*B
        
        IF (DD < EPS) THEN ! almost parallel
            IERR = 0
            EXIT
        ELSE 
            SC = (B*E-C*D)/DD
            TC = (A*E-B*D)/DD
            IF (SC>0._EB .AND. SC<1._EB .AND. TC>0._EB .AND. TC<1._EB ) THEN
                NXP = NXP+1
                XC = S1P0+SC*U
                DO K=1,3
                    XPC((NXP-1)*3+K) = XC(K)
                ENDDO
            ENDIF
        ENDIF
                
    ENDDO POLYGON_LOOP
ENDDO TRIANGLE_LOOP

!WRITE(LU_ERR,*) 'A', NXP
! add triangle vertices in polygon
DO I=1,3
    SELECT CASE(I)
    CASE(1)
        V = V1
    CASE(2)
        V = V2
    CASE(3)
        V = V3
    END SELECT
    
    IF (POINT_IN_BB(V,BB)) THEN
        NXP = NXP+1
        DO K=1,3
            XPC((NXP-1)*3+K) = V(K)
        ENDDO
    ENDIF
ENDDO

!WRITE(LU_ERR,*) 'B', NXP
! add polygon vertices in triangle
DO I=1,NP
    DO J=1,3
        V(J) = PC((I-1)*3+J)
    ENDDO
    IF (POINT_IN_TRIANGLE(V,V1,V2,V3)) THEN
        NXP = NXP+1
        DO J=1,3
            XPC((NXP-1)*3+J) = V(J)
        ENDDO
    ENDIF
ENDDO

!WRITE(LU_ERR,*) 'C', NXP

CALL ELIMATE_REPEATED_POINTS(NXP,XPC)

!WRITE(LU_ERR,*) 'D', NXP

IF (NXP > 3) THEN 
    CALL SORT_POLYGON_CORNERS(NXP,V1,V2,V3,XPC)
ENDIF

!WRITE(LU_ERR,*) 'E', NXP

IF (NXP >= 1) THEN
    IERR = 1 ! index for intersecting
ELSE
    IERR = 0
ENDIF

RETURN
END SUBROUTINE TRIANGLE_POLYGON_POINTS


SUBROUTINE ELIMATE_REPEATED_POINTS(NP,PC)
USE MATH_FUNCTIONS, ONLY:NORM2
IMPLICIT NONE

INTEGER, INTENT(INOUT):: NP
REAL(EB), INTENT(INOUT) :: PC(3*NP)
INTEGER :: NP2,I,J,K
REAL(EB), PARAMETER :: EPS=1.E-6_EB
REAL(EB) :: U(3),V(3),W(3)

I = 1
DO WHILE (I <= NP-1)
    DO K=1,3
        U(K) = PC(3*(I-1)+K)
    ENDDO
    
    J = I+1
    NP2 = NP
    DO WHILE (J <= NP2)
        DO K=1,3
            V(K) = PC(3*(J-1)+K)
        ENDDO
        W = U-V
        IF (NORM2(W) <= EPS) THEN
            DO K=3*J+1,3*NP
                PC(K-3) = PC(K)
            ENDDO
            NP = NP-1
            J = J-1
        ENDIF
        J = J+1
        IF (J > NP) EXIT
    ENDDO
    I = I+1
ENDDO

RETURN
END SUBROUTINE ELIMATE_REPEATED_POINTS



LOGICAL FUNCTION POINT_IN_BB(V1,BB)
IMPLICIT NONE

REAL(EB), INTENT(IN) :: V1(3),BB(6)

POINT_IN_BB=.FALSE.
IF ( V1(1)>=BB(1).AND.V1(1)<=BB(2) .AND. &
     V1(2)>=BB(3).AND.V1(2)<=BB(4) .AND. &
     V1(3)>=BB(5).AND.V1(3)<=BB(6) ) THEN
   POINT_IN_BB=.TRUE.
   RETURN
ENDIF

RETURN
END FUNCTION POINT_IN_BB


SUBROUTINE LINE_SEG_TRI_PLANE_INTERSECT(IERR,IERR2,Q,V1,V2,V3,P0,P1)
USE MATH_FUNCTIONS, ONLY:CROSS_PRODUCT
IMPLICIT NONE

INTEGER, INTENT(OUT) :: IERR
REAL(EB), INTENT(OUT) :: Q(3)
REAL(EB), INTENT(IN) :: V1(3),V2(3),V3(3),P0(3),P1(3)
REAL(EB) :: E1(3),E2(3),S(3),U,V,TMP,T,D(3),P(3)
REAL(EB), PARAMETER :: EPS=1.E-10_EB
INTEGER :: IERR2

IERR  = 0
IERR2 = 1
! IERR=1:  line segment intersect with the plane
! IERR2=1: the intersection point is in the triangle

! Schneider and Eberly, Section 11.1

D = P1-P0

E1 = V2-V1
E2 = V3-V1

CALL CROSS_PRODUCT(P,D,E2)

TMP = DOT_PRODUCT(P,E1)

IF ( ABS(TMP)<EPS ) RETURN

TMP = 1._EB/TMP
S = P0-V1

U = TMP*DOT_PRODUCT(S,P)
IF (U<0._EB .OR. U>1._EB) IERR2=0

CALL CROSS_PRODUCT(Q,S,E1)
V = TMP*DOT_PRODUCT(D,Q)
IF (V<0._EB .OR. (U+V)>1._EB) IERR2=0

T = TMP*DOT_PRODUCT(E2,Q)
Q = P0 + T*D ! the intersection point

IF (T>=0._EB .AND. T<=1._EB) IERR=1

END SUBROUTINE LINE_SEG_TRI_PLANE_INTERSECT


REAL(EB) FUNCTION POLYGON_AREA(NP,PC)
IMPLICIT NONE
! Calculate the area of a polygon

INTEGER, INTENT(IN) :: NP
REAL(EB), INTENT(IN) :: PC(3*NP)
INTEGER :: I,K
REAL(EB) :: V1(3),V2(3),V3(3)
    
POLYGON_AREA = 0._EB
V3 = 0._EB ! mean of the polygon

DO I=1,NP
    DO K=1,3
        V3(K) = V3(K)+PC((I-1)*3+K)/NP
    ENDDO
ENDDO
DO I=1,NP
    IF (I < NP) THEN
        DO K=1,3
            V1(K) = PC((I-1)*3+K)
            V2(K) = PC(I*3+K)
        ENDDO
    ELSE IF (I == NP) THEN        
        DO K=1,3
            V1(K) = PC((I-1)*3+K)
            V2(K) = PC(K)
        ENDDO
    ENDIF
    POLYGON_AREA = POLYGON_AREA+TRIANGLE_AREA(V1,V2,V3)
ENDDO

RETURN
END FUNCTION POLYGON_AREA


SUBROUTINE GET_CUTCELL_AREA()
IMPLICIT NONE

REAL(EB) :: V1(3),V2(3),V3(3),BB(6)
REAL(EB) :: PC(18), AREA,AREA0, XPCTMP(27)
INTEGER :: IERR,NP,NXP
REAL(EB), ALLOCATABLE, DIMENSION(:) :: XPC
TYPE (MESH_TYPE), POINTER :: M

V1(1) = VERTEX(FACET(1)%VERTEX(1))%X
V1(2) = VERTEX(FACET(1)%VERTEX(1))%Y
V1(3) = VERTEX(FACET(1)%VERTEX(1))%Z

V2(1) = VERTEX(FACET(1)%VERTEX(2))%X
V2(2) = VERTEX(FACET(1)%VERTEX(2))%Y
V2(3) = VERTEX(FACET(1)%VERTEX(2))%Z

V3(1) = VERTEX(FACET(1)%VERTEX(3))%X
V3(2) = VERTEX(FACET(1)%VERTEX(3))%Y
V3(3) = VERTEX(FACET(1)%VERTEX(3))%Z

M => MESHES(1)
BB(1) = M%X(1)
BB(2) = M%X(2)
BB(3) = M%Y(1)
BB(4) = M%Y(2)
BB(5) = M%Z(1)
BB(6) = M%Z(2)

AREA0 = TRIANGLE_AREA(V1,V2,V3)
    
CALL TRIANGLE_BOX_INTERSECT(IERR,V1,V2,V3,BB)

IF (IERR == 0) THEN
    WRITE(LU_ERR,*) 'The triangle is not intersecting with the BBox'
    
ELSE IF (IERR == 1) THEN
    ! if the triangle intersects with the BB
    ! next to determine the in intersection points
    
    CALL TRI_PLANE_BOX_INTERSECT(NP,PC,V1,V2,V3,BB)
    WRITE(LU_ERR,*) 'The intermeidate polygon has NP vertices, NP=', NP
!    DO I=1,NP
!        WRITE(LU_ERR,*) (PC((I-1)*3+J),J=1,3)
!    ENDDO
    
    ! get the intersection points, then calculate the area
    CALL TRIANGLE_POLYGON_POINTS(IERR,NXP,XPCTMP,V1,V2,V3,NP,PC,BB)
    WRITE(LU_ERR,*) 'The final polygon has NXP vertices,    NXP = ', NXP
    ALLOCATE(XPC(3*NXP))
    XPC = XPCTMP
!    INTP_LOOP: DO I=1,NXP
!        WRITE(LU_ERR,*) (XPC((I-1)*3+J),J=1,3)
!    ENDDO INTP_LOOP
    
    IF (IERR == 1)  AREA = POLYGON_AREA(NXP,XPC)
    WRITE(LU_ERR,*) 'The intersection area is = ', AREA

    WRITE(LU_ERR,*) 'The triangle     area is = ', AREA0
    
!    CALL WRITE_MATLAB_VISU(V1,V2,V3,NXP,XPC,NP,PC)
ENDIF

RETURN
END SUBROUTINE GET_CUTCELL_AREA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End Cut-cell subroutines by Charles Luo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE GET_REV_geom(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') geomrev(INDEX(geomrev,':')+1:LEN_TRIM(geomrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') geomdate

END SUBROUTINE GET_REV_geom


END MODULE COMPLEX_GEOMETRY
