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
PUBLIC :: INIT_IBM,TRILINEAR,GETU,GETGRAD,INIT_FACE,GET_REV_geom, &
          READ_GEOM,READ_VERT,READ_FACE,READ_VOLU,LINKED_LIST_INSERT
 
CONTAINS

SUBROUTINE READ_GEOM

CHARACTER(LABEL_LENGTH) :: ID,SHAPE,TFILE,SURF_ID
CHARACTER(60) :: BNDC_FILENAME,GEOC_FILENAME
REAL(EB) :: XYZ(3),ORIENTATION(3),XB(6),RADIUS,VELOCITY(3),OMEGA,RGB(3),PIXELS,ROUGHNESS
LOGICAL :: HOLE,TWO_SIDED
INTEGER :: IOS,IZERO,N
TYPE(GEOMETRY_TYPE), POINTER :: G=>NULL()
NAMELIST /GEOM/ ID,SHAPE,XB,XYZ,ORIENTATION,RADIUS,VELOCITY,OMEGA,HOLE,RGB,PIXELS,TWO_SIDED,TFILE,ROUGHNESS,&
                BNDC_FILENAME,DT_GEOC,DT_BNDC,GEOC_FILENAME,SURF_ID

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

READ_GEOM_LOOP: DO N=1,N_GEOM
   G=>GEOMETRY(N)
   
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
   DT_GEOC =  1._EB
   DT_BNDC =  1._EB
   SURF_ID = 'INERT'
   BNDC_FILENAME = 'null'
   GEOC_FILENAME = 'null'

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
   G%BNDC_FILENAME = BNDC_FILENAME
   G%GEOC_FILENAME = GEOC_FILENAME
   
   IF (ABS(G%U0)>0._EB .OR. ABS(G%V0)>0._EB .OR. ABS(G%W0)>0._EB) G%TRANSLATE = .TRUE.
   IF (ABS(G%OMEGA)>0._EB) G%ROTATE = .TRUE.
   
   SELECT CASE(G%SHAPE)
      CASE('COMPLEX')
         G%ISHAPE = ICOMPLEX
      CASE DEFAULT
         CALL SHUTDOWN('ERROR: unrecognized SHAPE on GEOM line')
   END SELECT
   
   ! Allocate bounding box arrays
   
   ALLOCATE(G%MIN_I(NMESHES),STAT=IZERO)
   CALL ChkMemErr('READ_GEOM','MIN_I',IZERO)
   ALLOCATE(G%MAX_I(NMESHES),STAT=IZERO)
   CALL ChkMemErr('READ_GEOM','MAX_I',IZERO)
   ALLOCATE(G%MIN_J(NMESHES),STAT=IZERO)
   CALL ChkMemErr('READ_GEOM','MIN_J',IZERO)
   ALLOCATE(G%MAX_J(NMESHES),STAT=IZERO)
   CALL ChkMemErr('READ_GEOM','MAX_J',IZERO)
   ALLOCATE(G%MIN_K(NMESHES),STAT=IZERO)
   CALL ChkMemErr('READ_GEOM','MIN_K',IZERO)
   ALLOCATE(G%MAX_K(NMESHES),STAT=IZERO)
   CALL ChkMemErr('READ_GEOM','MAX_K',IZERO)
   
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
CHARACTER(LABEL_LENGTH) :: SURF_ID='INERT'
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
CALL ChkMemErr('READ','FACET',IZERO)

READ_FACE_LOOP: DO I=1,N_FACE
   
   CALL CHECKREAD('FACE',LU_INPUT,IOS)
   IF (IOS==1) EXIT READ_FACE_LOOP
   
   ! Read the FACE line
   
   READ(LU_INPUT,FACE,END=37)

   IF (ANY(N>N_VERT)) THEN
      WRITE(MESSAGE,'(A,A,A)') 'ERROR: problem with FACE line ',TRIM(CHAR(I)),', N>N_VERT' 
      CALL SHUTDOWN(MESSAGE)
   ENDIF
   
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

   ! Allocate 1D arrays

   IF (N_TRACKED_SPECIES>0) THEN
      ALLOCATE(FACET(I)%RHODW(N_TRACKED_SPECIES),STAT=IZERO)
      CALL ChkMemErr('READ','FACET%RHODW',IZERO)
      ALLOCATE(FACET(I)%ZZ_F(N_TRACKED_SPECIES),STAT=IZERO)
      CALL ChkMemErr('READ','FACET%ZZ_F',IZERO)
   ENDIF

ENDDO READ_FACE_LOOP
37 REWIND(LU_INPUT)

CALL INIT_FACE

END SUBROUTINE READ_FACE


SUBROUTINE INIT_FACE
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
IMPLICIT NONE

INTEGER :: I,IZERO
REAL(EB) :: N_VEC(3),N_LENGTH,U_VEC(3),V_VEC(3),V1(3),V2(3),V3(3)
TYPE(FACET_TYPE), POINTER :: FC
TYPE(SURFACE_TYPE), POINTER :: SF
REAL(EB), PARAMETER :: TOL=1.E-10_EB

DO I=1,N_FACE

   FC=>FACET(I)

   V1 = (/VERTEX(FC%VERTEX(1))%X,VERTEX(FC%VERTEX(1))%Y,VERTEX(FC%VERTEX(1))%Z/)
   V2 = (/VERTEX(FC%VERTEX(2))%X,VERTEX(FC%VERTEX(2))%Y,VERTEX(FC%VERTEX(2))%Z/)
   V3 = (/VERTEX(FC%VERTEX(3))%X,VERTEX(FC%VERTEX(3))%Y,VERTEX(FC%VERTEX(3))%Z/)

   U_VEC = V2-V1
   V_VEC = V3-V1

   CALL CROSS_PRODUCT(N_VEC,U_VEC,V_VEC)
   N_LENGTH = SQRT(DOT_PRODUCT(N_VEC,N_VEC))

   IF (N_LENGTH>TOL) THEN
      FC%NVEC = N_VEC/N_LENGTH
   ELSE
      FC%NVEC = 0._EB
   ENDIF

   FC%AW = TRIANGLE_AREA(V1,V2,V3)
   IF (SURFACE(FC%SURF_INDEX)%TMP_FRONT>0._EB) THEN
      FC%TMP_F = SURFACE(FC%SURF_INDEX)%TMP_FRONT
   ELSE
      FC%TMP_F = TMPA
   ENDIF
   FC%TMP_G = TMPA
   FC%BOUNDARY_TYPE = SOLID_BOUNDARY

   IF (RADIATION) THEN
      SF => SURFACE(FC%SURF_INDEX)
      ALLOCATE(FC%ILW(1:SF%NRA,1:SF%NSB),STAT=IZERO)
      CALL ChkMemErr('INIT_FACE','FACET%ILW',IZERO)
   ENDIF

ENDDO

! Surface work arrays

IF (RADIATION) THEN
   ALLOCATE(FACE_WORK1(N_FACE),STAT=IZERO)
   CALL ChkMemErr('INIT_FACE','FACE_WORK1',IZERO)
   ALLOCATE(FACE_WORK2(N_FACE),STAT=IZERO)
   CALL ChkMemErr('INIT_FACE','FACE_WORK2',IZERO)
ENDIF

END SUBROUTINE INIT_FACE


SUBROUTINE READ_VOLU

INTEGER :: N(4),I,IOS,IZERO
CHARACTER(LABEL_LENGTH) :: MATL_ID
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
USE PHYSICAL_FUNCTIONS, ONLY: LES_FILTER_WIDTH_FUNCTION
IMPLICIT NONE
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
INTEGER :: I,J,K,N,IERR,IERR1,IERR2,I_MIN,I_MAX,J_MIN,J_MAX,K_MIN,K_MAX,IC,IOR,IIG,JJG,KKG
INTEGER :: NP,NXP,IZERO,LU,CUTCELL_COUNT,GEOM_TYPE,OWNER_INDEX,VERSION,NV,NF,DUMMY_INTEGER
TYPE (MESH_TYPE), POINTER :: M
REAL(EB) :: BB(6),V1(3),V2(3),V3(3),AREA,PC(18),XPC(60),V_POLYGON_CENTROID(3),VC,AREA_CHECK,&
            X_MIN,X_MAX,Y_MIN,Y_MAX,Z_MIN,Z_MAX
LOGICAL :: EX,OP
CHARACTER(60) :: FN
CHARACTER(100) :: MESSAGE
REAL(FB) :: DUMMY_FB_REAL,TIME_STRU
REAL(EB), PARAMETER :: CUTCELL_TOLERANCE=1.E-10_EB,MIN_AREA=1.E-16_EB,TOL=1.E-9_EB
!LOGICAL :: END_OF_LIST
TYPE(GEOMETRY_TYPE), POINTER :: G
TYPE(FACET_TYPE), POINTER :: FC=>NULL()
TYPE(CUTCELL_LINKED_LIST_TYPE), POINTER :: CL=>NULL()

IF (EVACUATION_ONLY(NM)) RETURN
M => MESHES(NM)

! primitive geometry loop

IF (N_GEOM>0) THEN
   IF (ABS(T-T_BEGIN)<TWO_EPSILON_EB .OR. ANY(GEOMETRY%TRANSLATE) .OR. ANY(GEOMETRY%ROTATE)) THEN
      M%U_MASK=1
      M%V_MASK=1
      M%W_MASK=1
      M%P_MASK=1
   ENDIF
ENDIF

GEOM_LOOP: DO N=1,N_GEOM

   G => GEOMETRY(N)

   IF ( ABS(T-T_BEGIN)>TWO_EPSILON_EB .AND. (.NOT. G%TRANSLATE) .AND. (.NOT. G%ROTATE) ) CYCLE GEOM_LOOP
   
   ! acceleration (not implemented yet)
   
   G%U = G%U0
   G%V = G%V0
   G%W = G%W0
   
   ! translation (linear for now)
   
   G%X = G%X0 + G%U*T
   G%Y = G%Y0 + G%V*T
   G%Z = G%Z0 + G%W*T
        
   DELTA_IBM = LES_FILTER_WIDTH_FUNCTION(M%DX(1),M%DY(1),M%DZ(1))

   ! find bounding box

   G%MIN_I(NM) = M%IBAR
   G%MIN_J(NM) = M%JBAR
   G%MIN_K(NM) = M%KBAR

   IF (X_MIN>=M%XS .AND. X_MIN<=M%XF) G%MIN_I(NM) = MAX(0,FLOOR((X_MIN-DELTA_IBM-M%XS)/M%DX(1))-1)
   IF (Y_MIN>=M%YS .AND. Y_MIN<=M%YF) G%MIN_J(NM) = MAX(0,FLOOR((Y_MIN-DELTA_IBM-M%YS)/M%DY(1))-1)
   IF (Z_MIN>=M%ZS .AND. Z_MIN<=M%ZF) G%MIN_K(NM) = MAX(0,FLOOR((Z_MIN-DELTA_IBM-M%ZS)/M%DZ(1))-1)
   
   G%MAX_I(NM) = 0
   G%MAX_J(NM) = 0
   G%MAX_K(NM) = 0

   IF (X_MAX>=M%XS .AND. X_MAX<=M%XF) G%MAX_I(NM) = MIN(M%IBAR,CEILING((X_MAX+DELTA_IBM-M%XS)/M%DX(1))+1)
   IF (Y_MAX>=M%YS .AND. Y_MAX<=M%YF) G%MAX_J(NM) = MIN(M%JBAR,CEILING((Y_MAX+DELTA_IBM-M%YS)/M%DY(1))+1)
   IF (Z_MAX>=M%ZS .AND. Z_MAX<=M%ZF) G%MAX_K(NM) = MIN(M%KBAR,CEILING((Z_MAX+DELTA_IBM-M%ZS)/M%DZ(1))+1)
  
   
   IF ( G%MAX_I(NM)<G%MIN_I(NM) .OR. &
        G%MAX_J(NM)<G%MIN_J(NM) .OR. &
        G%MAX_K(NM)<G%MIN_K(NM) ) CYCLE GEOM_LOOP
   
ENDDO GEOM_LOOP

! reinitialize complex geometry from geometry coordinate (.gc) file based on DT_GEOC frequency

IF (T<GEOC_CLOCK) RETURN
GEOC_CLOCK = GEOC_CLOCK + DT_GEOC
IF (ABS(T-T_BEGIN)<TWO_EPSILON_EB) THEN
   IZERO = 0
   IF (.NOT.ALLOCATED(M%CUTCELL_INDEX)) ALLOCATE(M%CUTCELL_INDEX(0:IBP1,0:JBP1,0:KBP1),STAT=IZERO) 
   CALL ChkMemErr('INIT_GEOM','CUTCELL_INDEX',IZERO) 
ENDIF

M%CUTCELL_INDEX = 0
CUTCELL_COUNT = 0

GEOC_LOOP: DO N=1,N_GEOM
   IF (TRIM(GEOMETRY(N)%GEOC_FILENAME)=='null') CYCLE
   
   FN = TRIM(GEOMETRY(N)%GEOC_FILENAME)
   INQUIRE(FILE=FN,EXIST=EX,OPENED=OP,NUMBER=LU)
   IF (.NOT.EX) CALL SHUTDOWN('ERROR: GEOMetry Coordinate file does not exist.')
   IF (OP) CLOSE(LU)
   LU_GEOC = GET_FILE_NUMBER()
      
   GEOC_CHECK_LOOP: DO I=1,30
      OPEN(LU_GEOC,FILE=FN,ACTION='READ',FORM='UNFORMATTED')
      READ(LU_GEOC) OWNER_INDEX     ! always 1, written by Abq
      READ(LU_GEOC) VERSION
      READ(LU_GEOC) TIME_STRU,GEOM_TYPE ! stime,gem_type
      IF ( ABS(T-REAL(TIME_STRU,EB)) > DT_BNDC*0.1_EB) THEN
         CLOSE(LU_GEOC)
         WRITE(LU_ERR,'(4X,A,I2,A)')  'waiting ... ', I,' s'
         CALL SLEEP(1)
      ELSE
         WRITE(LU_ERR,'(4X,A,F10.2,A,F10.2)')  'GEOM was updated at ',T,' s, GEOM Time:', TIME_STRU
         EXIT GEOC_CHECK_LOOP
      ENDIF
      IF (I==30) CALL SHUTDOWN('ERROR: BNDC FILE WAS NOT UPDATED BY STRUCTURE CODE')
   ENDDO GEOC_CHECK_LOOP
   
   IF (T>=REAL(TIME_STRU,EB)) THEN
      SELECT CASE(GEOM_TYPE)
         CASE(0)
            READ(LU_GEOC) NV,NF
            
            IF (NV>0 .AND. .NOT.ALLOCATED(FB_REAL_VERT_ARRAY)) THEN
               ALLOCATE(FB_REAL_VERT_ARRAY(NV*3),STAT=IZERO)
               CALL ChkMemErr('INIT_IBM','FB_REAL_VERT_ARRAY',IZERO)
            ENDIF
            IF (NF>0) THEN
               IF (.NOT.ALLOCATED(INT_FACE_VALS_ARRAY)) THEN
                  ALLOCATE(INT_FACE_VALS_ARRAY(NF*3),STAT=IZERO)
                  CALL ChkMemErr('INIT_IBM','FB_REAL_FACE_ARRAY',IZERO)
               ENDIF
               IF (.NOT.ALLOCATED(INT_SURF_VALS_ARRAY)) THEN
                  ALLOCATE(INT_SURF_VALS_ARRAY(NF),STAT=IZERO)
                  CALL ChkMemErr('INIT_IBM','INT_SURF_VALS_ARRAY',IZERO)
               ENDIF
            ENDIF
      
            IF (NV>0) READ(LU_GEOC) ((FB_REAL_VERT_ARRAY((I-1)*3+J),J=1,3),I=1,NV)
            DO I=1,NV
               VERTEX(I)%X = REAL(FB_REAL_VERT_ARRAY((I-1)*3+1),EB)
               VERTEX(I)%Y = REAL(FB_REAL_VERT_ARRAY((I-1)*3+2),EB)
               VERTEX(I)%Z = REAL(FB_REAL_VERT_ARRAY((I-1)*3+3),EB)
!                IF (T<10) VERTEX(I)%Z = VERTEX(I)%Z+0.03
!                IF (T>10 .AND. T<20) VERTEX(I)%Y = VERTEX(I)%Y+0.1
!                IF (T>20 .AND. T<30) VERTEX(I)%Z = VERTEX(I)%Z-0.02
!                IF (T>30 .AND. T<60) VERTEX(I)%Y = VERTEX(I)%Y-0.05
            ENDDO
            IF (NF>0) READ(LU_GEOC) ((INT_FACE_VALS_ARRAY((I-1)*3+J),J=1,3),I=1,NF)
            IF (NF>0) READ(LU_GEOC) (INT_SURF_VALS_ARRAY(I),I=1,NF)
            ! assuming the surface triangles have the same vertices
         CASE(1)
            READ(LU_GEOC) (DUMMY_FB_REAL,I=1,8)
      END SELECT
   ELSE
      CALL SHUTDOWN('ERROR: GEOMetry Coordinate file was not updated')
   ENDIF
ENDDO GEOC_LOOP

FACE_LOOP: DO N=1,N_FACE

   ! re-initialize the cutcell linked list
   IF (ASSOCIATED(FACET(N)%CUTCELL_LIST)) CALL CUTCELL_DESTROY(FACET(N)%CUTCELL_LIST)

   V1 = (/VERTEX(FACET(N)%VERTEX(1))%X,VERTEX(FACET(N)%VERTEX(1))%Y,VERTEX(FACET(N)%VERTEX(1))%Z/)
   V2 = (/VERTEX(FACET(N)%VERTEX(2))%X,VERTEX(FACET(N)%VERTEX(2))%Y,VERTEX(FACET(N)%VERTEX(2))%Z/)
   V3 = (/VERTEX(FACET(N)%VERTEX(3))%X,VERTEX(FACET(N)%VERTEX(3))%Y,VERTEX(FACET(N)%VERTEX(3))%Z/)
   FACET(N)%AW = TRIANGLE_AREA(V1,V2,V3)
   
   BB(1) = MIN(V1(1),V2(1),V3(1))
   BB(2) = MAX(V1(1),V2(1),V3(1))
   BB(3) = MIN(V1(2),V2(2),V3(2))
   BB(4) = MAX(V1(2),V2(2),V3(2))
   BB(5) = MIN(V1(3),V2(3),V3(3))
   BB(6) = MAX(V1(3),V2(3),V3(3))
   
   I_MIN = MAX(1,FLOOR((BB(1)-M%XS)/M%DX(1))-1) ! assumes uniform grid for now
   J_MIN = MAX(1,FLOOR((BB(3)-M%YS)/M%DY(1))-1)
   K_MIN = MAX(1,FLOOR((BB(5)-M%ZS)/M%DZ(1))-1)

   I_MAX = MIN(M%IBAR,CEILING((BB(2)-M%XS)/M%DX(1))+1)
   J_MAX = MIN(M%JBAR,CEILING((BB(4)-M%YS)/M%DY(1))+1)
   K_MAX = MIN(M%KBAR,CEILING((BB(6)-M%ZS)/M%DZ(1))+1)

   DO K=K_MIN,K_MAX
      DO J=J_MIN,J_MAX
         DO I=I_MIN,I_MAX

            BB(1) = M%X(I-1)
            BB(2) = M%X(I)
            BB(3) = M%Y(J-1)
            BB(4) = M%Y(J)
            BB(5) = M%Z(K-1)
            BB(6) = M%Z(K)
            CALL TRIANGLE_BOX_INTERSECT(IERR,V1,V2,V3,BB)
            
            IF (IERR==1) THEN
               CALL TRIANGLE_ON_CELL_SURF(IERR1,FACET(N)%NVEC,V1,M%XC(I),M%YC(J),M%ZC(K),M%DX(I),M%DY(J),M%DZ(K))  
               IF (IERR1==-1) CYCLE ! remove the possibility of double counting
                              
               CALL TRI_PLANE_BOX_INTERSECT(NP,PC,V1,V2,V3,BB)
               CALL TRIANGLE_POLYGON_POINTS(IERR2,NXP,XPC,V1,V2,V3,NP,PC,BB)
               IF (IERR2 == 1)  THEN                  
                  AREA = POLYGON_AREA(NXP,XPC)
                  IF (AREA > MIN_AREA) THEN
                     
                     ! check if the cutcell area needs to be assigned to a neighbor cell
                     V_POLYGON_CENTROID = POLYGON_CENTROID(NXP,XPC)
                     CALL POLYGON_CLOSE_TO_EDGE(IOR,FACET(N)%NVEC,V_POLYGON_CENTROID,&
                                                M%XC(I),M%YC(J),M%ZC(K),M%DX(I),M%DY(J),M%DZ(K))
                     IF (IOR/=0) THEN ! assign the cutcell area to a neighbor cell
                        SELECT CASE(IOR)
                           CASE(1)
                              IF (I==M%IBAR) THEN 
                                 IOR=0
                                 EXIT
                              ENDIF
                              IF (M%CUTCELL_INDEX(I+1,J,K)==0) THEN
                                 CUTCELL_COUNT = CUTCELL_COUNT+1
                                 IC = CUTCELL_COUNT
                                 M%CUTCELL_INDEX(I+1,J,K) = IC
                               ELSE
                                 IC = M%CUTCELL_INDEX(I+1,J,K)
                               ENDIF
                           CASE(-1)
                              IF (I==1) THEN 
                                 IOR=0
                                 EXIT
                              ENDIF
                              IF (M%CUTCELL_INDEX(I-1,J,K)==0) THEN
                                 CUTCELL_COUNT = CUTCELL_COUNT+1
                                 IC = CUTCELL_COUNT
                                 M%CUTCELL_INDEX(I-1,J,K) = IC
                               ELSE
                                 IC = M%CUTCELL_INDEX(I-1,J,K)
                               ENDIF
                           CASE(2)
                              IF (J==M%JBAR) THEN 
                                 IOR=0
                                 EXIT
                              ENDIF
                              IF (M%CUTCELL_INDEX(I,J+1,K)==0) THEN
                                 CUTCELL_COUNT = CUTCELL_COUNT+1
                                 IC = CUTCELL_COUNT
                                 M%CUTCELL_INDEX(I,J+1,K) = IC
                               ELSE
                                 IC = M%CUTCELL_INDEX(I,J+1,K)
                               ENDIF
                           CASE(-2)
                              IF (J==1) THEN 
                                 IOR=0
                                 EXIT
                              ENDIF
                              IF (M%CUTCELL_INDEX(I,J-1,K)==0) THEN
                                 CUTCELL_COUNT = CUTCELL_COUNT+1
                                 IC = CUTCELL_COUNT
                                 M%CUTCELL_INDEX(I,J-1,K) = IC
                               ELSE
                                 IC = M%CUTCELL_INDEX(I,J-1,K)
                               ENDIF
                           CASE(3)
                              IF (K==M%KBAR) THEN 
                                 IOR=0
                                 EXIT
                              ENDIF
                              IF (M%CUTCELL_INDEX(I,J,K+1)==0) THEN
                                 CUTCELL_COUNT = CUTCELL_COUNT+1
                                 IC = CUTCELL_COUNT
                                 M%CUTCELL_INDEX(I,J,K+1) = IC
                               ELSE
                                 IC = M%CUTCELL_INDEX(I,J,K+1)
                               ENDIF
                           CASE(-3)
                              IF (K==1) THEN 
                                 IOR=0
                                 EXIT
                              ENDIF
                              IF (M%CUTCELL_INDEX(I,J,K-1)==0) THEN
                                 CUTCELL_COUNT = CUTCELL_COUNT+1
                                 IC = CUTCELL_COUNT
                                 M%CUTCELL_INDEX(I,J,K-1) = IC
                               ELSE
                                 IC = M%CUTCELL_INDEX(I,J,K-1)
                               ENDIF
                        END SELECT
                     ENDIF
                     
                     IF (IOR==0) THEN 
                        IF (M%CUTCELL_INDEX(I,J,K)==0) THEN
                           CUTCELL_COUNT = CUTCELL_COUNT+1
                           IC = CUTCELL_COUNT
                           M%CUTCELL_INDEX(I,J,K) = IC
                        ELSE
                           IC = M%CUTCELL_INDEX(I,J,K)
                        ENDIF
                     ENDIF
                                            
                     CALL CUTCELL_INSERT(IC,AREA,FACET(N)%CUTCELL_LIST)
                  ENDIF
               ENDIF
            ENDIF

         ENDDO
      ENDDO
   ENDDO

ENDDO FACE_LOOP

! Create arrays to hold cutcell indices

CUTCELL_INDEX_IF: IF (CUTCELL_COUNT>0) THEN

   IF (ALLOCATED(M%I_CUTCELL)) DEALLOCATE(M%I_CUTCELL)
   IF (ALLOCATED(M%J_CUTCELL)) DEALLOCATE(M%J_CUTCELL)
   IF (ALLOCATED(M%K_CUTCELL)) DEALLOCATE(M%K_CUTCELL)

   ALLOCATE(M%I_CUTCELL(CUTCELL_COUNT),STAT=IZERO) 
   CALL ChkMemErr('READ','I_CUTCELL',IZERO) 
   M%I_CUTCELL = -1
   
   ALLOCATE(M%J_CUTCELL(CUTCELL_COUNT),STAT=IZERO) 
   CALL ChkMemErr('READ','J_CUTCELL',IZERO) 
   M%J_CUTCELL = -1
   
   ALLOCATE(M%K_CUTCELL(CUTCELL_COUNT),STAT=IZERO) 
   CALL ChkMemErr('READ','K_CUTCELL',IZERO) 
   M%K_CUTCELL = -1
 
   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
            IC = M%CUTCELL_INDEX(I,J,K)
            IF (IC>0) THEN
               M%I_CUTCELL(IC) = I
               M%J_CUTCELL(IC) = J
               M%K_CUTCELL(IC) = K
            ENDIF
         ENDDO
      ENDDO
   ENDDO

ENDIF CUTCELL_INDEX_IF

!CL=>FACET(1)%CUTCELL_LIST
!IF ( ASSOCIATED(CL) ) THEN
!    END_OF_LIST=.FALSE.
!    DO WHILE (.NOT.END_OF_LIST)
!       print *, CL%INDEX, CL%AREA
!       CL=>CL%NEXT
!       IF ( .NOT.ASSOCIATED(CL) ) THEN
!          print *,'done printing linked list!'
!          END_OF_LIST=.TRUE.
!       ENDIF
!    ENDDO
!ENDIF

! Set up any face parameters related to cut cells and check area sums

DO N=1,N_FACE

   FC=>FACET(N)
   CL=>FC%CUTCELL_LIST

   FC%DN=0._EB
   FC%RDN=0._EB
   AREA_CHECK=0._EB

   CUTCELL_LOOP: DO
      IF ( .NOT. ASSOCIATED(CL) ) EXIT CUTCELL_LOOP ! if the next index does not exist, exit the loop
      IC = CL%INDEX
      IIG = M%I_CUTCELL(IC)
      JJG = M%J_CUTCELL(IC)
      KKG = M%K_CUTCELL(IC)
      VC = M%DX(IIG)*M%DY(JJG)*M%DZ(KKG)

      FC%DN = FC%DN + CL%AREA*VC

      AREA_CHECK = AREA_CHECK + CL%AREA

      CL=>CL%NEXT ! point to the next index in the linked list
   ENDDO CUTCELL_LOOP

   IF (ABS(FC%AW-AREA_CHECK)>CUTCELL_TOLERANCE) THEN
      WRITE(MESSAGE,'(A,1I4)') 'ERROR: cutcell area checksum failed for facet ',N
      CALL SHUTDOWN(MESSAGE)
   ENDIF

   IF (FC%DN>CUTCELL_TOLERANCE) THEN
      FC%DN = FC%DN**ONTH ! wall normal length scale
      FC%RDN = 1._EB/FC%DN
   ENDIF

ENDDO

! Read boundary condition from file

BNDC_LOOP: DO N=1,N_GEOM
   IF (TRIM(GEOMETRY(N)%BNDC_FILENAME)=='null') CYCLE
   
   FN = TRIM(GEOMETRY(N)%BNDC_FILENAME)
   INQUIRE(FILE=FN,EXIST=EX,OPENED=OP,NUMBER=LU)
   IF (.NOT.EX) CALL SHUTDOWN('Error: boundary condition file does not exist.')
   IF (OP) CLOSE(LU)
   LU_BNDC = GET_FILE_NUMBER()
      
   BNDC_CHECK_LOOP: DO I=1,30
      OPEN(LU_BNDC,FILE=FN,ACTION='READ',FORM='UNFORMATTED')
      READ(LU_BNDC) OWNER_INDEX     ! 1 means written by Abq, 0 by FDS
      READ(LU_BNDC) VERSION         ! version
      READ(LU_BNDC) TIME_STRU       ! stime
      IF (OWNER_INDEX /=1 .OR. T-REAL(TIME_STRU,EB)>DT_BNDC*0.1_EB) THEN
         CLOSE(LU_BNDC)
         WRITE(LU_ERR,'(4X,A,F10.2,A)')  'BNDC not updated at ',T,' s'
         CALL SLEEP(1)
      ELSE
         WRITE(LU_ERR,'(4X,A,F10.2,A)')  'BNDC was updated at ',T,' s'
         EXIT BNDC_CHECK_LOOP
      ENDIF
      IF (I==30) THEN
         IF (OWNER_INDEX == 0) THEN
            CALL SHUTDOWN("ERROR: BNDC FILE WAS NOT UPDATED BY STRUCTURE CODE")
         ELSE
            CALL SHUTDOWN("ERROR: TIME MARKS do not match")
         ENDIF
      ENDIF
   ENDDO BNDC_CHECK_LOOP
   
   IF (ALLOCATED(FB_REAL_FACE_VALS_ARRAY)) DEALLOCATE(FB_REAL_FACE_VALS_ARRAY)
   ALLOCATE(FB_REAL_FACE_VALS_ARRAY(N_FACE),STAT=IZERO)
   CALL ChkMemErr('INIT_IBM','FB_REAL_FACE_VALS_ARRAY',IZERO)

   IF (T>=REAL(TIME_STRU,EB)) THEN
      ! n_vert_s_vals,n_vert_d_vals,n_face_s_vals,n_face_d_vals
      READ(LU_BNDC) (DUMMY_INTEGER,I=1,4)
      READ(LU_BNDC) (FB_REAL_FACE_VALS_ARRAY(I),I=1,N_FACE)
      DO I=1,N_FACE
         FC=>FACET(I)
         FC%TMP_F = REAL(FB_REAL_FACE_VALS_ARRAY(I),EB) + TMPM
      ENDDO
      IBM_FEM_COUPLING=.TRUE. ! immersed boundary method / finite-element method coupling
   ELSE
      BACKSPACE LU_BNDC
   ENDIF

ENDDO BNDC_LOOP

END SUBROUTINE INIT_IBM


! http://www.sdsc.edu/~tkaiser/f90.html#Linked lists
RECURSIVE SUBROUTINE LINKED_LIST_INSERT(ITEM,ROOT) 
   IMPLICIT NONE 
   TYPE(LINKED_LIST_TYPE), POINTER :: ROOT 
   INTEGER :: ITEM,IZERO
   IF (.NOT.ASSOCIATED(ROOT)) THEN 
      ALLOCATE(ROOT,STAT=IZERO)
      CALL ChkMemErr('GEOM','ROOT',IZERO)
      NULLIFY(ROOT%NEXT) 
      ROOT%INDEX = ITEM 
   ELSE 
      CALL LINKED_LIST_INSERT(ITEM,ROOT%NEXT) 
   ENDIF 
END SUBROUTINE LINKED_LIST_INSERT


RECURSIVE SUBROUTINE CUTCELL_INSERT(ITEM,AREA,ROOT) 
   IMPLICIT NONE 
   TYPE(CUTCELL_LINKED_LIST_TYPE), POINTER :: ROOT 
   INTEGER :: ITEM,IZERO
   REAL(EB):: AREA
   IF (.NOT.ASSOCIATED(ROOT)) THEN 
      ALLOCATE(ROOT,STAT=IZERO)
      CALL ChkMemErr('GEOM','ROOT',IZERO)
      NULLIFY(ROOT%NEXT) 
      ROOT%INDEX = ITEM
      ROOT%AREA = AREA 
   ELSE 
      CALL CUTCELL_INSERT(ITEM,AREA,ROOT%NEXT) 
   ENDIF 
END SUBROUTINE CUTCELL_INSERT


SUBROUTINE CUTCELL_DESTROY(ROOT)
IMPLICIT NONE
TYPE(CUTCELL_LINKED_LIST_TYPE), POINTER :: ROOT,CURRENT
DO WHILE (ASSOCIATED(ROOT))
  CURRENT => ROOT
  ROOT => CURRENT%NEXT
  DEALLOCATE(CURRENT)
ENDDO
RETURN
END SUBROUTINE CUTCELL_DESTROY


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

   IF ( ABS( SIGN(1._EB,DOT_PRODUCT(Q_VEC,N_VEC))-SIGN(1._EB,DOT_PRODUCT(R_VEC,N_VEC)) )>TWO_EPSILON_EB ) THEN
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
REAL(EB), PARAMETER :: EPS=1.E-16_EB

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
REAL(EB) :: XX,YY,ZZ

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

XX = DXI(1)/LL(1)
YY = DXI(2)/LL(2)
ZZ = DXI(3)/LL(3)

TRILINEAR = UU(0,0,0)*(1._EB-XX)*(1._EB-YY)*(1._EB-ZZ) + &
            UU(1,0,0)*XX*(1._EB-YY)*(1._EB-ZZ) +         & 
            UU(0,1,0)*(1._EB-XX)*YY*(1._EB-ZZ) +         &
            UU(0,0,1)*(1._EB-XX)*(1._EB-YY)*ZZ +         &
            UU(1,0,1)*XX*(1._EB-YY)*ZZ +                 &
            UU(0,1,1)*(1._EB-XX)*YY*ZZ +                 &
            UU(1,1,0)*XX*YY*(1._EB-ZZ) +                 & 
            UU(1,1,1)*XX*YY*ZZ

END FUNCTION TRILINEAR


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
REAL(EB) :: P0(3),P1(3),Q(3),PC_TMP(60)
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
IF (NP > 3) THEN 
   CALL SORT_POLYGON_CORNERS(NP,V1,V2,V3,PC_TMP)
ENDIF
DO I=1,NP*3
   PC(I) = PC_TMP(I)
ENDDO

RETURN
END SUBROUTINE TRI_PLANE_BOX_INTERSECT


SUBROUTINE SORT_POLYGON_CORNERS(NP,V1,V2,V3,PC)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
IMPLICIT NONE
! Sort all the corners of a polygon
! Ref: Gernot Hoffmann, Cube Plane Intersection.

INTEGER, INTENT(IN) :: NP
REAL(EB), INTENT(INOUT) :: PC(60)
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
IF (ABS(POLY_NORM(1)) >= ABS(POLY_NORM(2)) .AND. ABS(POLY_NORM(1)) >= ABS(POLY_NORM(3)) ) THEN
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
REAL(EB), INTENT(OUT) :: XPC(60)
REAL(EB), INTENT(IN) :: V1(3),V2(3),V3(3),PC(18),BB(6)
INTEGER :: I,J,K
REAL(EB) :: U(3),V(3),W(3),S1P0(3),XC(3)
REAL(EB) :: A,B,C,D,E,DD,SC,TC
REAL(EB), PARAMETER :: EPS=1.E-20_EB,TOL=1.E-12_EB
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
      ELSE
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
         CYCLE
      ELSE 
         SC = (B*E-C*D)/DD
         TC = (A*E-B*D)/DD
         IF (SC>-TOL .AND. SC<1._EB+TOL .AND. TC>-TOL .AND. TC<1._EB+TOL ) THEN
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
REAL(EB), INTENT(INOUT) :: PC(60)
INTEGER :: NP2,I,J,K
REAL(EB), PARAMETER :: EPS=1.E-16_EB
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
REAL(EB), PARAMETER :: EPS=1.E-10_EB,TOL=1.E-15
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

IF (T>=0._EB-TOL .AND. T<=1._EB+TOL) IERR=1

RETURN
END SUBROUTINE LINE_SEG_TRI_PLANE_INTERSECT


REAL(EB) FUNCTION POLYGON_AREA(NP,PC)
IMPLICIT NONE
! Calculate the area of a polygon

INTEGER, INTENT(IN) :: NP
REAL(EB), INTENT(IN) :: PC(60)
INTEGER :: I,K
REAL(EB) :: V1(3),V2(3),V3(3)
    
POLYGON_AREA = 0._EB
V3 = POLYGON_CENTROID(NP,PC)

DO I=1,NP
   IF (I < NP) THEN
      DO K=1,3
         V1(K) = PC((I-1)*3+K)
         V2(K) = PC(I*3+K)
      ENDDO
   ELSE
      DO K=1,3
         V1(K) = PC((I-1)*3+K)
         V2(K) = PC(K)
      ENDDO
   ENDIF
   POLYGON_AREA = POLYGON_AREA+TRIANGLE_AREA(V1,V2,V3)
ENDDO

RETURN
END FUNCTION POLYGON_AREA


REAL(EB) FUNCTION POLYGON_CENTROID(NP,PC)
IMPLICIT NONE
! Calculate the centroid of polygon vertices

DIMENSION :: POLYGON_CENTROID(3)
INTEGER, INTENT(IN) :: NP
REAL(EB), INTENT(IN) :: PC(60)
INTEGER :: I,K

POLYGON_CENTROID = 0._EB
DO I=1,NP
   DO K=1,3
      POLYGON_CENTROID(K) = POLYGON_CENTROID(K)+PC((I-1)*3+K)/NP
   ENDDO
ENDDO

RETURN
END FUNCTION POLYGON_CENTROID


SUBROUTINE TRIANGLE_ON_CELL_SURF(IERR1,N_VEC,V,XC,YC,ZC,DX,DY,DZ)
USE MATH_FUNCTIONS, ONLY:NORM2
IMPLICIT NONE

INTEGER, INTENT(OUT) :: IERR1
REAL(EB), INTENT(IN) :: N_VEC(3),V(3),XC,YC,ZC,DX,DY,DZ
REAL(EB) :: DIST(3),TOL=1.E-15_EB

IERR1 = 1
DIST = 0._EB
!IF (NORM2(N_VEC)>1._EB) N_VEC = N_VEC/NORM2(N_VEC)

IF (N_VEC(1)==1._EB .OR. N_VEC(1)==-1._EB) THEN
   DIST(1) = XC-V(1)
   IF ( ABS(ABS(DIST(1))-DX*0.5_EB)<TOL .AND. DOT_PRODUCT(DIST,N_VEC)<0._EB) THEN
      IERR1 = -1
   ENDIF
   RETURN
ENDIF

IF (N_VEC(2)==1._EB .OR. N_VEC(2)==-1._EB) THEN
   DIST(2) = YC-V(2)
   IF ( ABS(ABS(DIST(2))-DY*0.5_EB)<TOL .AND. DOT_PRODUCT(DIST,N_VEC)<0._EB) THEN
      IERR1 = -1
   ENDIF
   RETURN
ENDIF

IF (N_VEC(3)==1._EB .OR. N_VEC(3)==-1._EB) THEN
   DIST(3) = ZC-V(3)
   IF ( ABS(ABS(DIST(3))-DZ*0.5_EB)<TOL .AND. DOT_PRODUCT(DIST,N_VEC)<0._EB) THEN
      IERR1 = -1
   ENDIF
   RETURN
ENDIF

RETURN
END SUBROUTINE TRIANGLE_ON_CELL_SURF


SUBROUTINE POLYGON_CLOSE_TO_EDGE(IOR,N_VEC,V,XC,YC,ZC,DX,DY,DZ)
IMPLICIT NONE
INTEGER, INTENT(OUT) :: IOR
REAL(EB), INTENT(IN) :: N_VEC(3),V(3),XC,YC,ZC,DX,DY,DZ
REAL(EB) :: DIST(3),DMAX
REAL(EB), PARAMETER :: TOLERANCE=0.01_EB

IOR = 0
DIST(1) = XC-V(1)
DIST(2) = YC-V(2)
DIST(3) = ZC-V(3)

IF (ABS(DIST(1)/DX) >= ABS(DIST(2)/DY) .AND. ABS(DIST(1)/DX) >= ABS(DIST(3)/DZ)) THEN
   DMAX = ABS(DIST(1)/DX*2._EB)
   IF (DMAX < (1._EB-TOLERANCE) .OR. DOT_PRODUCT(DIST,N_VEC) > 0._EB) RETURN
   IF (DIST(1) < 0._EB) THEN
      IOR = 1
   ELSE
      IOR = -1
   ENDIF
ELSEIF (ABS(DIST(2)/DY) >= ABS(DIST(3)/DZ)) THEN
   DMAX = ABS(DIST(2)/DY*2._EB)
   IF (DMAX < (1._EB-TOLERANCE) .OR. DOT_PRODUCT(DIST,N_VEC) > 0._EB) RETURN
   IF (DIST(2) < 0._EB) THEN
      IOR = 2
   ELSE
      IOR = -2
   ENDIF
ELSE
   DMAX = ABS(DIST(3)/DZ*2._EB)
   IF (DMAX < (1._EB-TOLERANCE) .OR. DOT_PRODUCT(DIST,N_VEC) > 0._EB) RETURN
   IF (DIST(3) < 0._EB) THEN
      IOR = 3
   ELSE
      IOR = -3
   ENDIF
ENDIF
   
RETURN
END SUBROUTINE POLYGON_CLOSE_TO_EDGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End cut cell subroutines by Charles Luo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE GET_REV_geom(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') geomrev(INDEX(geomrev,':')+2:LEN_TRIM(geomrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') geomdate

END SUBROUTINE GET_REV_geom


END MODULE COMPLEX_GEOMETRY
