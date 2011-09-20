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
          READ_GEOM,READ_VERT,READ_FACE,READ_VOLU
 
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

INTEGER :: N(3),I,IOS,IZERO
CHARACTER(30) :: SURF_ID
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
   FACET(I)%SURF_ID = TRIM(SURF_ID)

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
IMPLICIT NONE
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
INTEGER :: I,J,K,NG,NV,I_MIN,I_MAX,J_MIN,J_MAX,K_MIN,K_MAX
TYPE (MESH_TYPE), POINTER :: M
TYPE (GEOMETRY_TYPE), POINTER :: G
REAL(EB) :: DELTA,RP,XU(3),PP(3),DP,TIME,TOL=1.E-9_EB, &
            X_MIN,Y_MIN,Z_MIN,X_MAX,Y_MAX,Z_MAX,XX(4),YY(4),ZZ(4),XP(3),BB(6)

IF (ICYC>1 .AND. N_GEOM==0) RETURN

TIME = T
M => MESHES(NM)

M%U_MASK=1
M%V_MASK=1
M%W_MASK=1
M%P_MASK=1

! geometry loop

GEOM_LOOP: DO NG=1,N_GEOM

   G => GEOMETRY(NG)

   IF (ICYC>1 .AND. (.NOT. G%TRANSLATE) .AND. (.NOT. G%ROTATE)) CYCLE GEOM_LOOP
   
   ! acceleration (not implemented yet)
   
   G%U = G%U0
   G%V = G%V0
   G%W = G%W0
   
   ! translation (linear for now)
   
   G%X = G%X0 + G%U*TIME
   G%Y = G%Y0 + G%V*TIME
   G%Z = G%Z0 + G%W*TIME
        
   IF (TWO_D) THEN
      DELTA = MIN(M%DX(1),M%DZ(1))
   ELSE
      DELTA = MIN(M%DX(1),M%DY(1),M%DZ(1))
   ENDIF
   
   IBM_UVWMAX = MAXVAL((/ABS(G%U),ABS(G%V),ABS(G%W)/))/DELTA

   ! find bounding box

   SELECT_SHAPE: SELECT CASE(G%ISHAPE)
      CASE(IBOX) SELECT_SHAPE
         G%X1 = G%X1 + G%U*DT
         G%X2 = G%X2 + G%U*DT
         G%Y1 = G%Y1 + G%V*DT
         G%Y2 = G%Y2 + G%V*DT
         G%Z1 = G%Z1 + G%W*DT
         G%Z2 = G%Z2 + G%W*DT
         X_MIN = G%X1
         X_MAX = G%X2
         Y_MIN = G%Y1
         Y_MAX = G%Y2
         Z_MIN = G%Z1
         Z_MAX = G%Z2
         G%HL(1) = 0.5_EB*(X_MAX-X_MIN) + TOL
         G%HL(2) = 0.5_EB*(X_MAX-X_MIN) + TOL
         G%HL(3) = 0.5_EB*(X_MAX-X_MIN) + TOL
      CASE(ISPHERE) SELECT_SHAPE
         X_MIN = G%X-G%RADIUS
         Y_MIN = G%Y-G%RADIUS
         Z_MIN = G%Z-G%RADIUS
         X_MAX = G%X+G%RADIUS
         Y_MAX = G%Y+G%RADIUS
         Z_MAX = G%Z+G%RADIUS
         IBM_UVWMAX = IBM_UVWMAX + G%RADIUS*MAXVAL((/ABS(G%OMEGA_X),ABS(G%OMEGA_Y),ABS(G%OMEGA_Z)/))*RDX(1)
      CASE(ICYLINDER) SELECT_SHAPE
         G%HL(1) = 0.5_EB*(G%X2-G%X1)
         G%HL(2) = 0.5_EB*(G%Y2-G%Y1)
         G%HL(3) = 0.5_EB*(G%Z2-G%Z1)
         
         IF (ABS(G%XOR-1._EB)<EPSILON_EB) THEN ! cylinder aligned with x axis
            X_MIN = G%X-G%HL(1)
            Y_MIN = G%Y-G%RADIUS
            Z_MIN = G%Z-G%RADIUS
            X_MAX = G%X+G%HL(1)
            Y_MAX = G%Y+G%RADIUS
            Z_MAX = G%Z+G%RADIUS
         ENDIF
         
         IF (ABS(G%YOR-1._EB)<EPSILON_EB) THEN ! cylinder aligned with y axis
            X_MIN = G%X-G%RADIUS
            Y_MIN = G%Y-G%HL(2)
            Z_MIN = G%Z-G%RADIUS
            X_MAX = G%X+G%RADIUS
            Y_MAX = G%Y+G%HL(2)
            Z_MAX = G%Z+G%RADIUS
         ENDIF
         
         IF (ABS(G%ZOR-1._EB)<EPSILON_EB) THEN ! cylinder aligned with z axis
            X_MIN = G%X-G%RADIUS
            Y_MIN = G%Y-G%RADIUS
            Z_MIN = G%Z-G%HL(3)
            X_MAX = G%X+G%RADIUS
            Y_MAX = G%Y+G%RADIUS
            Z_MAX = G%Z+G%HL(3)
         ENDIF
         
      CASE(IPLANE)
         X_MIN = M%XS
         Y_MIN = M%YS
         Z_MIN = M%ZS
         X_MAX = M%XF
         Y_MAX = M%YF
         Z_MAX = M%ZF
         PP   = (/G%X,G%Y,G%Z/)
         G%NN = (/G%XOR,G%YOR,G%ZOR/) - PP          ! normal vector to plane
         G%NN = G%NN/SQRT(DOT_PRODUCT(G%NN,G%NN))   ! unit normal
   END SELECT SELECT_SHAPE

   G%MIN_I(NM) = M%IBAR
   G%MIN_J(NM) = M%JBAR
   G%MIN_K(NM) = M%KBAR

   IF (X_MIN>=M%XS .AND. X_MIN<=M%XF) G%MIN_I(NM) = MAX(0,FLOOR((X_MIN-M%XS)/M%DX(1))-1)
   IF (Y_MIN>=M%YS .AND. Y_MIN<=M%YF) G%MIN_J(NM) = MAX(0,FLOOR((Y_MIN-M%YS)/M%DY(1))-1)
   IF (Z_MIN>=M%ZS .AND. Z_MIN<=M%ZF) G%MIN_K(NM) = MAX(0,FLOOR((Z_MIN-M%ZS)/M%DZ(1))-1)
   
   G%MAX_I(NM) = 0
   G%MAX_J(NM) = 0
   G%MAX_K(NM) = 0

   IF (X_MAX>=M%XS .AND. X_MAX<=M%XF) G%MAX_I(NM) = MIN(M%IBAR,CEILING((X_MAX-M%XS)/M%DX(1))+1)
   IF (Y_MAX>=M%YS .AND. Y_MAX<=M%YF) G%MAX_J(NM) = MIN(M%JBAR,CEILING((Y_MAX-M%YS)/M%DY(1))+1)
   IF (Z_MAX>=M%ZS .AND. Z_MAX<=M%ZF) G%MAX_K(NM) = MIN(M%KBAR,CEILING((Z_MAX-M%ZS)/M%DZ(1))+1)
   
   IF (TWO_D) THEN
      G%MIN_J(NM)=1
      G%MAX_J(NM)=1
   ENDIF
   
   IF ( G%MAX_I(NM)<G%MIN_I(NM) .OR. &
        G%MAX_J(NM)<G%MIN_J(NM) .OR. &
        G%MAX_K(NM)<G%MIN_K(NM) ) CYCLE GEOM_LOOP
   
   ! mask cells

   DO K=G%MIN_K(NM),G%MAX_K(NM)
      DO J=G%MIN_J(NM),G%MAX_J(NM)
         DO I=G%MIN_I(NM),G%MAX_I(NM)
         
            MASK_SHAPE: SELECT CASE(G%ISHAPE)
            
               CASE(IBOX) MASK_SHAPE
               
                  !! this will not work for overlapping geometry, but use for now
                  !M%U_MASK(I,J,K)=1
                  !M%V_MASK(I,J,K)=1
                  !M%W_MASK(I,J,K)=1
                  !M%P_MASK(I,J,K)=1
                  
                  ! see if the point is inside geometry
                  IF (ABS( M%X(I)-G%X)<G%HL(1) .AND. &
                      ABS(M%YC(J)-G%Y)<G%HL(2) .AND. &
                      ABS(M%ZC(K)-G%Z)<G%HL(3)) M%U_MASK(I,J,K) = -1
                  
                  IF (ABS(M%XC(I)-G%X)<G%HL(1) .AND. &
                      ABS( M%Y(J)-G%Y)<G%HL(2) .AND. &
                      ABS(M%ZC(K)-G%Z)<G%HL(3)) M%V_MASK(I,J,K) = -1
                  
                  IF (ABS(M%XC(I)-G%X)<G%HL(1) .AND. &
                      ABS(M%YC(J)-G%Y)<G%HL(2) .AND. &
                      ABS( M%Z(K)-G%Z)<G%HL(3)) M%W_MASK(I,J,K) = -1
                  
                  ! see if the point is in surface layer
                  IF (X_MAX<M%X(I) .AND. M%X(I)<X_MAX+DELTA) M%U_MASK(I,J,K) = 0
                  IF (Y_MAX<M%Y(J) .AND. M%Y(J)<Y_MAX+DELTA) M%V_MASK(I,J,K) = 0
                  IF (Z_MAX<M%Z(K) .AND. M%Z(K)<Z_MAX+DELTA) M%W_MASK(I,J,K) = 0
                  
                  IF (X_MIN-DELTA<M%X(I) .AND. M%X(I)<X_MIN) M%U_MASK(I,J,K) = 0
                  IF (Y_MIN-DELTA<M%Y(J) .AND. M%Y(J)<Y_MIN) M%V_MASK(I,J,K) = 0
                  IF (Z_MIN-DELTA<M%Z(K) .AND. M%Z(K)<Z_MIN) M%W_MASK(I,J,K) = 0
                  
               CASE(ISPHERE) MASK_SHAPE
               
                  !M%U_MASK(I,J,K)=1
                  !M%V_MASK(I,J,K)=1
                  !M%W_MASK(I,J,K)=1
                  !M%P_MASK(I,J,K)=1
               
                  RP = SQRT( (M%X(I)-G%X)**2+(M%YC(J)-G%Y)**2+(M%ZC(K)-G%Z)**2 )
                  IF (RP-G%RADIUS < DELTA) M%U_MASK(I,J,K) = 0
                  IF (RP-G%RADIUS < 0._EB) M%U_MASK(I,J,K) = -1
                  
                  RP = SQRT( (M%XC(I)-G%X)**2+(M%Y(J)-G%Y)**2+(M%ZC(K)-G%Z)**2 )
                  IF (RP-G%RADIUS < DELTA) M%V_MASK(I,J,K) = 0
                  IF (RP-G%RADIUS < 0._EB) M%V_MASK(I,J,K) = -1
                  
                  RP = SQRT( (M%XC(I)-G%X)**2+(M%YC(J)-G%Y)**2+(M%Z(K)-G%Z)**2 )
                  IF (RP-G%RADIUS < DELTA) M%W_MASK(I,J,K) = 0
                  IF (RP-G%RADIUS < 0._EB) M%W_MASK(I,J,K) = -1
                  
                  RP = SQRT( (M%XC(I)-G%X)**2+(M%YC(J)-G%Y)**2+(M%ZC(K)-G%Z)**2 )
                  IF (RP-G%RADIUS < DELTA) M%P_MASK(I,J,K) = 0
                  IF (RP-G%RADIUS < 0._EB) M%P_MASK(I,J,K) = -1
                  
               CASE(ICYLINDER) MASK_SHAPE
               
                  !M%U_MASK(I,J,K)=1
                  !M%V_MASK(I,J,K)=1
                  !M%W_MASK(I,J,K)=1
                  !M%P_MASK(I,J,K)=1
               
                  CYLINDER_Y: IF (ABS(G%YOR-1._EB)<EPSILON_EB) THEN
                     RP = SQRT( (M%X(I)-G%X)**2+(M%ZC(K)-G%Z)**2 )
                     IF (RP-G%RADIUS < DELTA) M%U_MASK(I,J,K) = 0
                     IF (RP-G%RADIUS < 0._EB) M%U_MASK(I,J,K) = -1
                  
                     RP = SQRT( (M%XC(I)-G%X)**2+(M%ZC(K)-G%Z)**2 )
                     IF (RP-G%RADIUS < DELTA) M%V_MASK(I,J,K) = 0
                     IF (RP-G%RADIUS < 0._EB) M%V_MASK(I,J,K) = -1
                  
                     RP = SQRT( (M%XC(I)-G%X)**2+(M%Z(K)-G%Z)**2 )
                     IF (RP-G%RADIUS < DELTA) M%W_MASK(I,J,K) = 0
                     IF (RP-G%RADIUS < 0._EB) M%W_MASK(I,J,K) = -1
                  
                     RP = SQRT( (M%XC(I)-G%X)**2+(M%ZC(K)-G%Z)**2 )
                     IF (RP-G%RADIUS < DELTA) M%P_MASK(I,J,K) = 0
                     IF (RP-G%RADIUS < 0._EB) M%P_MASK(I,J,K) = -1
                  ENDIF CYLINDER_Y
                  
                  CYLINDER_Z: IF (ABS(G%ZOR-1._EB)<EPSILON_EB) THEN
                     RP = SQRT( (M%X(I)-G%X)**2+(M%YC(J)-G%Y)**2 )
                     IF (RP-G%RADIUS < DELTA) M%U_MASK(I,J,K) = 0
                     IF (RP-G%RADIUS < 0._EB) M%U_MASK(I,J,K) = -1
                  
                     RP = SQRT( (M%XC(I)-G%X)**2+(M%Y(J)-G%Y)**2 )
                     IF (RP-G%RADIUS < DELTA) M%V_MASK(I,J,K) = 0
                     IF (RP-G%RADIUS < 0._EB) M%V_MASK(I,J,K) = -1
                  
                     RP = SQRT( (M%XC(I)-G%X)**2+(M%YC(J)-G%Y)**2 )
                     IF (RP-G%RADIUS < DELTA) M%W_MASK(I,J,K) = 0
                     IF (RP-G%RADIUS < 0._EB) M%W_MASK(I,J,K) = -1
                  
                     RP = SQRT( (M%XC(I)-G%X)**2+(M%YC(J)-G%Y)**2 )
                     IF (RP-G%RADIUS < DELTA) M%P_MASK(I,J,K) = 0
                     IF (RP-G%RADIUS < 0._EB) M%P_MASK(I,J,K) = -1
                     
                  ENDIF CYLINDER_Z

                  
               CASE(IPLANE) MASK_SHAPE
               
                  !M%U_MASK(I,J,K)=1
                  !M%V_MASK(I,J,K)=1
                  !M%W_MASK(I,J,K)=1
                  !M%P_MASK(I,J,K)=1
               
                  ! see Section 10.3 Schneider and Eberly
                  
                  XU = (/M%X(I),M%YC(J),M%ZC(K)/)        ! point of interest
                  IF (G%TWO_SIDED) THEN
                     DP = ABS(DOT_PRODUCT(G%NN,XU-PP))   ! distance to plane
                  ELSE
                     DP = DOT_PRODUCT(G%NN,XU-PP)        ! signed distance to plane
                  ENDIF
                  IF (DP<DELTA) M%U_MASK(I,J,K) = 0
                  IF (DP<0._EB) M%U_MASK(I,J,K) = -1
                  
                  XU = (/M%XC(I),M%Y(J),M%ZC(K)/)
                  IF (G%TWO_SIDED) THEN
                     DP = ABS(DOT_PRODUCT(G%NN,XU-PP))
                  ELSE
                     DP = DOT_PRODUCT(G%NN,XU-PP)
                  ENDIF
                  IF (DP<DELTA) M%V_MASK(I,J,K) = 0
                  IF (DP<0._EB) M%V_MASK(I,J,K) = -1
                  
                  XU = (/M%XC(I),M%YC(J),M%Z(K)/)
                  IF (G%TWO_SIDED) THEN
                     DP = ABS(DOT_PRODUCT(G%NN,XU-PP))
                  ELSE
                     DP = DOT_PRODUCT(G%NN,XU-PP)
                  ENDIF
                  IF (DP<DELTA) M%W_MASK(I,J,K) = 0
                  IF (DP<0._EB) M%W_MASK(I,J,K) = -1
                  
            END SELECT MASK_SHAPE
      
         ENDDO
      ENDDO
   ENDDO

ENDDO GEOM_LOOP

! unstructured geometry

VOLUME_LOOP: DO NV=1,N_VOLU

   XX(1) = VERTEX(VOLUME(NV)%VERTEX(1))%X
   XX(2) = VERTEX(VOLUME(NV)%VERTEX(2))%X
   XX(3) = VERTEX(VOLUME(NV)%VERTEX(3))%X
   XX(4) = VERTEX(VOLUME(NV)%VERTEX(4))%X

   YY(1) = VERTEX(VOLUME(NV)%VERTEX(1))%Y
   YY(2) = VERTEX(VOLUME(NV)%VERTEX(2))%Y
   YY(3) = VERTEX(VOLUME(NV)%VERTEX(3))%Y
   YY(4) = VERTEX(VOLUME(NV)%VERTEX(4))%Y

   ZZ(1) = VERTEX(VOLUME(NV)%VERTEX(1))%Z
   ZZ(2) = VERTEX(VOLUME(NV)%VERTEX(2))%Z
   ZZ(3) = VERTEX(VOLUME(NV)%VERTEX(3))%Z
   ZZ(4) = VERTEX(VOLUME(NV)%VERTEX(4))%Z

   ! bounding box

   X_MIN = MINVAL(XX)
   X_MAX = MAXVAL(XX)

   Y_MIN = MINVAL(YY)
   Y_MAX = MAXVAL(YY)

   Z_MIN = MINVAL(ZZ)
   Z_MAX = MAXVAL(ZZ)

   I_MIN = M%IBAR
   J_MIN = M%JBAR
   K_MIN = M%KBAR

   IF (X_MIN>=M%XS .AND. X_MIN<=M%XF) I_MIN = MAX(0,FLOOR((X_MIN-M%XS)/M%DX(1))-1)
   IF (Y_MIN>=M%YS .AND. Y_MIN<=M%YF) J_MIN = MAX(0,FLOOR((Y_MIN-M%YS)/M%DY(1))-1)
   IF (Z_MIN>=M%ZS .AND. Z_MIN<=M%ZF) K_MIN = MAX(0,FLOOR((Z_MIN-M%ZS)/M%DZ(1))-1)
   
   I_MAX = 0
   J_MAX = 0
   K_MAX = 0

   IF (X_MAX>=M%XS .AND. X_MAX<=M%XF) I_MAX = MIN(M%IBAR,CEILING((X_MAX-M%XS)/M%DX(1))+1)
   IF (Y_MAX>=M%YS .AND. Y_MAX<=M%YF) J_MAX = MIN(M%JBAR,CEILING((Y_MAX-M%YS)/M%DY(1))+1)
   IF (Z_MAX>=M%ZS .AND. Z_MAX<=M%ZF) K_MAX = MIN(M%KBAR,CEILING((Z_MAX-M%ZS)/M%DZ(1))+1)
   
   IF (TWO_D) THEN
      J_MIN=1
      J_MAX=1
   ENDIF
   
   IF ( I_MAX<I_MIN .OR. &
        J_MAX<J_MIN .OR. &
        K_MAX<K_MIN ) CYCLE VOLUME_LOOP

   ! mask cells

   DO K=K_MIN,K_MAX
      DO J=J_MIN,J_MAX
         DO I=I_MIN,I_MAX

            ! test pcell (cell center)
            XP = (/M%XC(I),M%YC(J),M%ZC(K)/)
            IF ( POINT_IN_TETRAHEDRON(XP,XX,YY,ZZ) ) M%P_MASK(I,J,K)=-1

            ! test ucell
            XP = (/M%X(I),M%YC(J),M%ZC(K)/)
            IF ( POINT_IN_TETRAHEDRON(XP,XX,YY,ZZ) ) M%U_MASK(I,J,K)=-1

            ! test vcell
            XP = (/M%XC(I),M%Y(J),M%ZC(K)/)
            IF ( POINT_IN_TETRAHEDRON(XP,XX,YY,ZZ) ) M%V_MASK(I,J,K)=-1

            ! test wcell
            XP = (/M%XC(I),M%YC(J),M%Z(K)/)
            IF ( POINT_IN_TETRAHEDRON(XP,XX,YY,ZZ) ) M%W_MASK(I,J,K)=-1 

         ENDDO
      ENDDO
   ENDDO

ENDDO VOLUME_LOOP

! point in general polyhedron

RAY_TEST: IF (.TRUE.) THEN

! bounding box (will use better data structure later)

BB(1) = MINVAL(VERTEX%X)
BB(2) = MAXVAL(VERTEX%X)
BB(3) = MINVAL(VERTEX%Y)
BB(4) = MAXVAL(VERTEX%Y)
BB(5) = MINVAL(VERTEX%Z)
BB(6) = MAXVAL(VERTEX%Z)

DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         ! test pcell (cell center)
         XP = (/M%XC(I),M%YC(J),M%ZC(K)/)
         IF ( POINT_IN_POLYHEDRON(XP,BB) ) M%P_MASK(I,J,K)=-1
      ENDDO
   ENDDO
ENDDO

DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         ! test ucell
         XP = (/M%X(I),M%YC(J),M%ZC(K)/)
         IF ( POINT_IN_POLYHEDRON(XP,BB) ) M%U_MASK(I,J,K)=-1
      ENDDO
   ENDDO
ENDDO

DO K=1,M%KBAR
   DO J=0,M%JBAR
      DO I=1,M%IBAR
         ! test vcell
         XP = (/M%XC(I),M%Y(J),M%ZC(K)/)
         IF ( POINT_IN_POLYHEDRON(XP,BB) ) M%V_MASK(I,J,K)=-1
      ENDDO
   ENDDO
ENDDO

DO K=0,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         ! test wcell
         XP = (/M%XC(I),M%YC(J),M%Z(K)/)
         IF ( POINT_IN_POLYHEDRON(XP,BB) ) M%W_MASK(I,J,K)=-1 
      ENDDO
   ENDDO
ENDDO

ENDIF RAY_TEST

END SUBROUTINE INIT_IBM


LOGICAL FUNCTION POINT_IN_TETRAHEDRON(XP,XX,YY,ZZ)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
IMPLICIT NONE

REAL(EB), INTENT(IN) :: XP(3),XX(4),YY(4),ZZ(4)
REAL(EB) :: U_VEC(3),V_VEC(3),N_VEC(3),Q_VEC(3),R_VEC(3)
INTEGER :: I,N(4,4)

! In this routine, we test all four faces of the tet volume defined by the points X(i),Y(i),Z(i); i=1:4.
! If the point is on the negative side of all the faces, it is inside the volume.

! define the vertex ordering for each face (store later)

N(1,:) = (/1,2,3,4/)
N(2,:) = (/1,3,4,2/)
N(3,:) = (/1,4,2,3/)
N(4,:) = (/2,4,3,1/)

POINT_IN_TETRAHEDRON=.TRUE. ! start by assuming the point is inside

FACE_LOOP: DO I=1,4

   ! vectors forming the sides of the triangle

   U_VEC = (/XX(N(I,2))-XX(N(I,1)),YY(N(I,2))-YY(N(I,1)),ZZ(N(I,2))-ZZ(N(I,1))/)
   V_VEC = (/XX(N(I,3))-XX(N(I,1)),YY(N(I,3))-YY(N(I,1)),ZZ(N(I,3))-ZZ(N(I,1))/)

   CALL CROSS_PRODUCT(N_VEC,U_VEC,V_VEC)

   ! form a vector from a point on the triangular surface to the point XP

   Q_VEC = XP-(/XX(N(I,1)),YY(N(I,1)),ZZ(N(I,1))/)

   ! also form a vector from the triangular surface to the other point on the volume defining inside

   R_VEC = (/XX(N(I,4)),YY(N(I,4)),ZZ(N(I,4))/)-(/XX(N(I,1)),YY(N(I,1)),ZZ(N(I,1))/)

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


LOGICAL FUNCTION POINT_IN_TRIANGLE(XP,XX,YY)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
IMPLICIT NONE

REAL(EB), INTENT(IN) :: XP(2),XX(3),YY(3)
REAL(EB) :: V_VEC(2),N_VEC(2),Q_VEC(2),R_VEC(2)
INTEGER :: I,N(3,3)

! This routine is similar to POINT_IN_TETRAHEDRON

N(1,:) = (/1,2,3/)
N(2,:) = (/2,3,1/)
N(3,:) = (/3,1,2/)

POINT_IN_TRIANGLE=.TRUE. ! start by assuming the point is inside

EDGE_LOOP: DO I=1,3

   ! vector along the direction of edge I

   V_VEC = (/XX(N(I,2))-XX(N(I,1)),YY(N(I,2))-YY(N(I,1))/)

   ! find vector normal to edge

   IF (ABS(V_VEC(2))>ZERO_P) THEN
      N_VEC = (/1._EB, -V_VEC(1)/V_VEC(2)/)
   ELSE
      N_VEC = (/0._EB, 1._EB/)
   ENDIF

   print *,N_VEC

   ! form a vector from a point on the edge to the point XP

   Q_VEC = XP-(/XX(N(I,1)),YY(N(I,1))/)

   ! also form a vector from the edge to the other point on the triangle defining inside

   R_VEC = (/XX(N(I,3)),YY(N(I,3))/)-(/XX(N(I,1)),YY(N(I,1))/)

   ! if the sign of the dot products are equal, the point is inside, else it is outside and we return

   IF ( ABS( SIGN(1._EB,DOT_PRODUCT(Q_VEC,N_VEC))-SIGN(1._EB,DOT_PRODUCT(R_VEC,N_VEC)) )>ZERO_P ) THEN
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


SUBROUTINE INIT_FACE
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
IMPLICIT NONE

INTEGER :: I
REAL(EB) :: N_VEC(3),N_LENGTH,U_VEC(3),V_VEC(3)
TYPE(FACET_TYPE), POINTER :: F
REAL(EB), PARAMETER :: TOL=1.E-10_EB

DO I=1,N_FACE

   F=>FACET(I)

   U_VEC(1) = VERTEX(F%VERTEX(2))%X - VERTEX(F%VERTEX(1))%X
   U_VEC(2) = VERTEX(F%VERTEX(2))%Y - VERTEX(F%VERTEX(1))%Y
   U_VEC(3) = VERTEX(F%VERTEX(2))%Z - VERTEX(F%VERTEX(1))%Z

   V_VEC(1) = VERTEX(F%VERTEX(3))%X - VERTEX(F%VERTEX(1))%X
   V_VEC(2) = VERTEX(F%VERTEX(3))%Y - VERTEX(F%VERTEX(1))%Y
   V_VEC(3) = VERTEX(F%VERTEX(3))%Z - VERTEX(F%VERTEX(1))%Z

   CALL CROSS_PRODUCT(N_VEC,U_VEC,V_VEC)
   N_LENGTH = SQRT(DOT_PRODUCT(N_VEC,N_VEC))

   IF (N_LENGTH>TOL) THEN
      F%NVEC(1) = N_VEC(1)/N_LENGTH
      F%NVEC(2) = N_VEC(2)/N_LENGTH
      F%NVEC(3) = N_VEC(3)/N_LENGTH
   ELSE
      CALL SHUTDOWN('ERROR: Invalid facet')
   ENDIF

ENDDO

END SUBROUTINE INIT_FACE


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


SUBROUTINE GETU(U_DATA,DXI,XI,XU,INDU,I_VEL,NM)
IMPLICIT NONE

REAL(EB), INTENT(OUT) :: U_DATA(0:1,0:1,0:1),DXI(3)
REAL(EB), INTENT(IN) :: XI(3),XU(3)
INTEGER, INTENT(IN) :: INDU(3),I_VEL,NM
TYPE(MESH_TYPE), POINTER :: M
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
INTEGER :: II,JJ,KK,N
CHARACTER(100) :: MESSAGE

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
   WRITE(MESSAGE,'(A)') 'ERROR: DXI<0 in GETU'
   CALL SHUTDOWN(MESSAGE)
ENDIF
IF (DXI(1)>M%DX(II) .OR. DXI(2)>M%DY(JJ) .OR. DXI(3)>M%DZ(KK)) THEN
   WRITE(MESSAGE,'(A)') 'ERROR: DXI>DX in GETU'
   CALL SHUTDOWN(MESSAGE)
ENDIF

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


SUBROUTINE GET_REV_geom(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') geomrev(INDEX(geomrev,':')+1:LEN_TRIM(geomrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') geomdate

END SUBROUTINE GET_REV_geom


END MODULE COMPLEX_GEOMETRY
