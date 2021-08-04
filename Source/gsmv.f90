
!  +++++++++++++++++++++++ BOXTETRA_ROUTINES ++++++++++++++++++++++++++

MODULE BOXTETRA_ROUTINES
USE PRECISION_PARAMETERS
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT
IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE

INTEGER, DIMENSION(0:3,0:5), TARGET :: BOX_PLANE2VERT
INTEGER, DIMENSION(0:3,0:5), TARGET :: BOX_PLANE2EDGE
INTEGER, DIMENSION(0:1,0:11), TARGET :: BOX_EDGE2VERT

INTEGER, DIMENSION(0:2,0:3), TARGET :: TETRA_PLANE2VERT
INTEGER, DIMENSION(0:2,0:3), TARGET :: TETRA_PLANE2EDGE
INTEGER, DIMENSION(0:1,0:5) :: TETRA_EDGE2VERT

INTEGER, PARAMETER :: MIN_X=0, MAX_X=1, MIN_Y=2, MAX_Y=3, MIN_Z=4, MAX_Z=5

INTEGER :: I, J

!       6-----------7
!      /.          /|
!    /  .        /  |
!   4-----------5   |
!   |   .       |   |
!   |   .       |   |
!   |   2.......|...3
!   |  .        |  /
!   | .         | /
!   |.          |/
!   0-----------1

! BOX_PLANE2VERT(edge,plane)
DATA ( (BOX_PLANE2VERT(I,J), I=0,3),J=0,5) /&
  0,2,4,6,&
  1,3,5,7,&
  0,1,4,5,&
  2,3,6,7,&
  0,1,2,3,&
  4,5,6,7 &
  /

!       ------7------
!      /.           /
!     2 .         3 |
!    /  .        /  |
!   ------6------   |
!   |  10       |  11
!   |   .       |   |
!   |   .....5..|...|
!   8  .        9  /
!   | 0         | 1
!   |.          |/
!   |----4------|

DATA ( (BOX_PLANE2EDGE(I,J), I=0,3),J=0,5) /&
  0,2,8,10,&
  1,3,9,11,&
  4,6,8,9,&
  5,7,10,11,&
  0,1,4,5,&
  2,3,6,7 &
  /

!       6-----7-----7
!      /.           /
!     2 .         3 |
!    /  .        /  |
!   4-----6-----5   |
!   |  10       |  11
!   |   .       |   |
!   |   2....5..|...3
!   8  .        9  /
!   | 0         | 1
!   |.          |/
!   0----4------1
! planes: 0-left 1-right 2-front 3-back 4-bottom 5-top
! edges: 0-bottom left  1-bottom right 2-top left   3-top right
!        4-bottom front 5-bottom back  6-top front  7-top back
!        8-front left   9-front right 10-back left 11-back right
! vertices: 0-bottom left front 1-bottom right front 2-bottom left back 3-bottom right back
!           4-top left front    5-top right front    6-top left back    7-top right back

 DATA ( (BOX_EDGE2VERT(I,J), I=0,1), J=0,11) /&
  0,2,  1,3,  4,6,  5,7,&
  0,1,  2,3,  4,5,  6,7,&
  0,4,  1,5,  2,6,  3,7&
  /

!           3
!          /.\
!         / . \
!        /  5  \
!       /   .   \
!      3    2    4
!     /   .   .   \
!    /  2       1  \
!   / .           . \
!  0-------0---------1

DATA ( (TETRA_PLANE2VERT(I,J), I=0,2),J=0,3) /&
  0,3,1,&
  1,3,2,&
  0,2,3,&
  0,1,2&
  /

DATA ( (TETRA_PLANE2EDGE(I,J), I=0,2),J=0,3) /&
  0,3,4,&
  1,4,5,&
  2,5,3,& !double check (was 2 3 5)
  0,1,2&
  /

 DATA ( (TETRA_EDGE2VERT(I,J), I=0,1), J=0,5) /&
  0,1,  1,2,  2,0,&
  0,3,  1,3,  2,3&
  /

PUBLIC TETRAHEDRON_VOLUME, REMOVE_DUPLICATE_VERTS, GEOMCLIP, DECIMATE_FB

CONTAINS

!  ------------------ TETRAHEDRON_VOLUME ------------------------

!              D1
!             /|\
!            / | \
!           /  |  \
!          /   |   \
!         /    |    \
!        /     B4    \
!       /     . .     \
!      /     .    .    \
!     /    .        .   \
!    /   .            .  \
!   /  .               .  \
!  / .                    .\
! C2------------------------A3

REAL(EB) FUNCTION TETRAHEDRON_VOLUME(A,B,C,D)

! determine the volume of a tetrahedron formed from vertices A, B, C and D

REAL(EB), DIMENSION(0:2), INTENT(IN) :: A, B, C, D
REAL(EB), DIMENSION(0:2) :: AMC, BMC, DMC, ACROSSB

AMC = A - C
BMC = B - C
DMC = D - C
CALL CROSS_PRODUCT(ACROSSB,AMC, BMC)
TETRAHEDRON_VOLUME = DOT_PRODUCT(ACROSSB,DMC)/6.0_EB
END FUNCTION TETRAHEDRON_VOLUME

! ---------------------------- REMOVE_DUPLICATE_VERTS ----------------------------------------

SUBROUTINE REMOVE_DUPLICATE_VERTS(N_VERTS,N_FACES,N_VOLUS,&
                            MAX_VERTS,MAX_FACES,MAX_VOLUS,FIRST_FACE_INDEX,&
                            VERTS,FACES,VOLUS,VERT_EPS)
INTEGER, INTENT(INOUT) :: N_VERTS
INTEGER, INTENT(IN) :: N_FACES, N_VOLUS
INTEGER, INTENT(IN) :: MAX_VERTS, MAX_FACES, MAX_VOLUS
INTEGER, INTENT(IN) :: FIRST_FACE_INDEX
INTEGER, TARGET, INTENT(INOUT) :: VOLUS(4*MAX_VOLUS)
REAL(EB), TARGET, INTENT(INOUT) :: VERTS(3*MAX_VERTS)
INTEGER, TARGET, INTENT(INOUT) :: FACES(3*MAX_FACES)
REAL(EB), INTENT(IN) :: VERT_EPS

REAL(EB), POINTER, DIMENSION(:) :: VI, VJ
REAL(EB), DIMENSION(3) :: VIMVJ
REAL(EB) :: NORM_VI, NORM_VJ, NORM_VIMVJ

INTEGER :: I, J, K
INTEGER :: OFFSET

OFFSET = 1-FIRST_FACE_INDEX
I = 1
DO WHILE (I<=N_VERTS)
   VI=>VERTS(3*I-2:3*I)
   NORM_VI = MAX(ABS(VI(1)),ABS(VI(2)),ABS(VI(3)))
   J = I+1
   DO WHILE (J<=N_VERTS)
      VJ=>VERTS(3*J-2:3*J)
      VIMVJ = VI-VJ
      NORM_VJ = MAX(ABS(VJ(1)),ABS(VJ(2)),ABS(VJ(3)))
      NORM_VIMVJ = MAX(ABS(VIMVJ(1)),ABS(VIMVJ(2)),ABS(VIMVJ(3)))
      IF (NORM_VIMVJ <= MAX(1.0_EB,NORM_VI,NORM_VJ)*VERT_EPS) THEN
         ! vertex I and J are the same
         ! first copy index J -> I in face list
         ! next copy index N_VERTS -> J in face list
         ! finally reduce N_VERTS by 1

         DO K = 1, 3*N_FACES
           IF (FACES(K)==J-OFFSET)FACES(K)=I-OFFSET
         ENDDO
         DO K = 1, 4*N_VOLUS
           IF (VOLUS(K)==J-OFFSET)VOLUS(K)=I-OFFSET
         ENDDO
         VJ(1:3)=VERTS(3*N_VERTS-2:3*N_VERTS)

         DO K = 1, 3*N_FACES
           IF (FACES(K)==N_VERTS-OFFSET)FACES(K)=J-OFFSET
         ENDDO
         DO K = 1, 4*N_VOLUS
           IF (VOLUS(K)==N_VERTS-OFFSET)VOLUS(K)=J-OFFSET
         ENDDO
         N_VERTS=N_VERTS-1
         CYCLE
      ENDIF
      J=J+1
   ENDDO
   I=I+1
ENDDO
END SUBROUTINE REMOVE_DUPLICATE_VERTS

!  ------------------ IN_BOX3 ------------------------

INTEGER FUNCTION IN_BOX3(XB,V)
REAL(EB) :: XB(6), V(3)

IF (XB(1)<=V(1) .AND.V(1)<=XB(2) .AND.&
    XB(3)<=V(2) .AND.V(2)<=XB(4) .AND.&
    XB(5)<=V(3) .AND.V(3)<=XB(6)) THEN
   IN_BOX3=1
   RETURN
ENDIF
IN_BOX3=0
END FUNCTION

!  ------------------ DISTANCE3 ------------------------

REAL(EB) FUNCTION DISTANCE3(V1,V2)
REAL(EB), INTENT(IN), DIMENSION(3) :: V1, V2
REAL(EB) :: DX, DY, DZ

DX = V1(1)-V2(1)
DY = V1(2)-V2(2)
DZ = V1(3)-V2(3)
DISTANCE3 = SQRT(DX*DX+DY*DY+DZ*DZ)
END FUNCTION DISTANCE3

!  ------------------ GET_VERTTYPE ------------------------

SUBROUTINE GET_VERTTYPE(NVERTS, TRIANGLES, NTRIANGLES, VERT_TYPE)
INTEGER, INTENT(IN) :: NVERTS, NTRIANGLES
INTEGER, INTENT(OUT), DIMENSION(NVERTS) :: VERT_TYPE
INTEGER, INTENT(IN), DIMENSION(3*NTRIANGLES), TARGET :: TRIANGLES

! classify each vertex in a geometry as either interior or exterior
! a vertex VI is interior if the vertices connected to I form a cycle or loop
! ie VI is connected to v1, v2, v3 and v1 -> v2 -> v3 -> v1
! if they don't form a loop then it is an exterior vertex

! exterior vertices (connected to a blockage or mesh boundary) won't be moved or deleted

INTEGER, DIMENSION(NVERTS) :: VERT_COUNT, TRIANGLE_COUNT
INTEGER :: I, J, K, VERTJ_INDEX, TRIJ_INDEX, MAXCOUNT, VERTK_INDEX
INTEGER, POINTER, DIMENSION(:) :: TRII
INTEGER, DIMENSION(:,:), ALLOCATABLE :: VERT_TRILIST

! count number of triangles connected to each vertex

TRIANGLE_COUNT(1:NVERTS) = 0
DO I = 1, NTRIANGLES
   TRII=>TRIANGLES(3*I-2:3*I)
   DO J = 1, 3
     VERTJ_INDEX = TRII(J)
     IF(VERTJ_INDEX>=1 .AND. VERTJ_INDEX<=NVERTS)TRIANGLE_COUNT(VERTJ_INDEX) = TRIANGLE_COUNT(VERTJ_INDEX) + 1
   END DO
END DO

MAXCOUNT = TRIANGLE_COUNT(1)
DO I = 2, NVERTS
  MAXCOUNT = MAX(MAXCOUNT,TRIANGLE_COUNT(I))
END DO

! construct a list of triangles connected to each vertex
! VERT_TRILIST(I,1) contains number of triangles connected to vertex I
! VERT_TRILIST(I,2-> ) contains the triangle indices

ALLOCATE(VERT_TRILIST(NVERTS,MAXCOUNT+1))
VERT_TRILIST(1:NVERTS,1)=0
DO I = 1, NTRIANGLES
   TRII=>TRIANGLES(3*I-2:3*I)
   DO J = 1, 3
     VERTJ_INDEX = TRII(J)
     IF(VERTJ_INDEX>=1 .AND. VERTJ_INDEX<=NVERTS)THEN
        VERT_TRILIST(VERTJ_INDEX,1) = VERT_TRILIST(VERTJ_INDEX,1) + 1 ! bump triangle count by 1 for vertex j
        K = VERT_TRILIST(VERTJ_INDEX,1)+1
        VERT_TRILIST(VERTJ_INDEX,K) = I  ! put triangle index into triangle count + 1
     ENDIF
   END DO
END DO

! count vertices connected to each vertex

VERT_TYPE(1:NVERTS) = 1
OUTER: DO I = 1, NVERTS
  VERT_COUNT(1:NVERTS) = 0
  DO J = 1, VERT_TRILIST(I,1) ! loop over triangles connected to vertex I
    TRIJ_INDEX = VERT_TRILIST(I,1+J)
    DO K = 1, 3 ! loop over vertices of triangle J
      VERTK_INDEX = TRIANGLES(3*TRIJ_INDEX-3+K)
      IF(VERTK_INDEX/=I)VERT_COUNT(VERTK_INDEX) = VERT_COUNT(VERTK_INDEX) + 1
    END DO
  END DO
  DO J = 1, NVERTS
    IF(VERT_COUNT(J)==1)THEN ! consider all vertices that are connected to vertex I
      VERT_TYPE(I) = 0       ! if all of these neighbors have two neighbors among this set then I is interior
      CYCLE OUTER            ! if at least one of these neighbors has only one neigbor among this set then I is on the exterior
    ENDIF
  END DO
END DO OUTER

DEALLOCATE(VERT_TRILIST)

END SUBROUTINE GET_VERTTYPE

!  ------------------ DECIMATE_FB ------------------------

SUBROUTINE DECIMATE_FB(VERTS, NVERTS, FACES, NFACES, MESH_BOUNDS_FB, DELTA)
INTEGER, INTENT(INOUT) :: NVERTS, NFACES
REAL(FB), INTENT(INOUT), DIMENSION(3*NVERTS) :: VERTS
INTEGER, INTENT(INOUT), DIMENSION(3*NFACES) :: FACES
REAL(FB), INTENT(IN) :: DELTA
REAL(FB), DIMENSION(6), INTENT(IN) :: MESH_BOUNDS_FB

REAL(EB), DIMENSION(3*NVERTS) :: VERTS_EB
REAL(EB) :: DELTA_EB
REAL(EB), DIMENSION(6) :: MESH_BOUNDS_EB

DELTA_EB = REAL(DELTA,EB)
VERTS_EB(1:3*NVERTS) = REAL(VERTS(1:3*NVERTS),EB)
MESH_BOUNDS_EB(1:6) = REAL(MESH_BOUNDS_FB(1:6),EB)
CALL DECIMATE(VERTS_EB,NVERTS,FACES,NFACES,MESH_BOUNDS_EB,DELTA_EB)
VERTS(1:3*NVERTS) = REAL(VERTS_EB(1:3*NVERTS),FB)
RETURN
END SUBROUTINE DECIMATE_FB

!  ------------------ AVERAGE_VERTS2 ------------------------

SUBROUTINE AVERAGE_VERTS2(V1,V1TYPE,V2,V2TYPE,MESH_BOUNDS,VAVG)
REAL(EB), DIMENSION(3), INTENT(IN) :: V1, V2
INTEGER, INTENT(IN) :: V1TYPE, V2TYPE
REAL(EB), DIMENSION(6), INTENT(IN) :: MESH_BOUNDS
REAL(EB), DIMENSION(3), INTENT(OUT) :: VAVG

REAL(EB) :: BOXEPS

BOXEPS = 0.001

IF (V1TYPE == 0) THEN
  VAVG(1:3) = V1(1:3)
  RETURN
ENDIF

IF (V2TYPE == 0) THEN
  VAVG(1:3) = V2(1:3)
  RETURN
ENDIF

IF (ABS(V1(1)-MESH_BOUNDS(1))<BOXEPS .OR. ABS(V1(1)-MESH_BOUNDS(2))<BOXEPS) THEN
  VAVG(1) = V1(1)
ELSE IF (ABS(V2(1)-MESH_BOUNDS(1))<BOXEPS .OR. ABS(V2(1)-MESH_BOUNDS(2))<BOXEPS) THEN
  VAVG(1) = V2(1)
ELSE
  VAVG(1) = (V1(1)+V2(1))/2.0
ENDIF

IF (ABS(V1(2)-MESH_BOUNDS(3))<BOXEPS .OR. ABS(V1(2)-MESH_BOUNDS(4))<BOXEPS) THEN
  VAVG(2) = V1(2)
ELSE IF (ABS(V2(2)-MESH_BOUNDS(3))<BOXEPS .OR. ABS(V2(2)-MESH_BOUNDS(4))<BOXEPS) THEN
  VAVG(2) = V2(2)
ELSE
  VAVG(2) = (V1(2)+V2(2))/2.0
ENDIF

IF ( ABS( V1(3)-MESH_BOUNDS(5) )<BOXEPS .OR. ABS(V1(3)-MESH_BOUNDS(6))<BOXEPS) THEN
  VAVG(3) = V1(3)
ELSE IF (ABS(V2(3)-MESH_BOUNDS(5))<BOXEPS .OR. ABS(V2(3)-MESH_BOUNDS(6))<BOXEPS) THEN
  VAVG(3) = V2(3)
ELSE
  VAVG(3) = (V1(3)+V2(3))/2.0
ENDIF

END SUBROUTINE AVERAGE_VERTS2

!  ------------------ AVERAGE_VERTS3 ------------------------

SUBROUTINE AVERAGE_VERTS3(V1,V1TYPE,V2,V2TYPE,V3,V3TYPE,MESH_BOUNDS,VAVG)
REAL(EB), DIMENSION(3), INTENT(IN) :: V1, V2, V3
INTEGER, INTENT(IN) :: V1TYPE, V2TYPE, V3TYPE
REAL(EB), DIMENSION(6), INTENT(IN) :: MESH_BOUNDS
REAL(EB), DIMENSION(3), INTENT(OUT) :: VAVG

REAL(EB) :: BOXEPS

BOXEPS = 0.001

IF (V1TYPE==0) THEN
  VAVG(1:3) = V1(1:3)
  RETURN
ENDIF

IF (V2TYPE==0) THEN
  VAVG(1:3) = V2(1:3)
  RETURN
ENDIF

IF (V3TYPE==0) THEN
  VAVG(1:3) = V3(1:3)
  RETURN
ENDIF

IF (ABS(V1(1)-MESH_BOUNDS(1))<BOXEPS .OR. ABS(V1(1)-MESH_BOUNDS(2))<BOXEPS) THEN
  VAVG(1) = V1(1)
ELSE IF (ABS(V2(1)-MESH_BOUNDS(1))<BOXEPS .OR. ABS(V2(1)-MESH_BOUNDS(2))<BOXEPS) THEN
  VAVG(1) = V2(1)
ELSE IF (ABS(V3(1)-MESH_BOUNDS(1))<BOXEPS .OR. ABS(V3(1)-MESH_BOUNDS(2))<BOXEPS) THEN
  VAVG(1) = V3(1)
ELSE
  VAVG(1) = (V1(1)+V2(1)+V3(1))/3.0
ENDIF

IF (ABS(V1(2)-MESH_BOUNDS(3))<BOXEPS .OR. ABS(V1(2)-MESH_BOUNDS(4))<BOXEPS) THEN
  VAVG(2) = V1(2)
ELSE IF (ABS(V2(2)-MESH_BOUNDS(3))<BOXEPS .OR. ABS(V2(2)-MESH_BOUNDS(4))<BOXEPS) THEN
  VAVG(2) = V2(2)
ELSE IF (ABS(V3(2)-MESH_BOUNDS(3))<BOXEPS .OR. ABS(V3(2)-MESH_BOUNDS(4))<BOXEPS) THEN
  VAVG(2) = V3(2)
ELSE
  VAVG(2) = (V1(2)+V2(2)+V3(2))/3.0
ENDIF

IF (ABS(V1(3)-MESH_BOUNDS(5))<BOXEPS .OR. ABS(V1(3)-MESH_BOUNDS(6))<BOXEPS) THEN
  VAVG(3) = V1(3)
ELSE IF (ABS(V2(3)-MESH_BOUNDS(5))<BOXEPS .OR. ABS(V2(3)-MESH_BOUNDS(6))<BOXEPS) THEN
  VAVG(3) = V2(3)
ELSE IF (ABS(V3(3)-MESH_BOUNDS(5))<BOXEPS .OR. ABS(V3(3)-MESH_BOUNDS(6))<BOXEPS) THEN
  VAVG(3) = V3(3)
ELSE
  VAVG(3) = (V1(3)+V2(3)+V3(3))/3.0
ENDIF

END SUBROUTINE AVERAGE_VERTS3

!  ------------------ DECIMATE ------------------------

SUBROUTINE DECIMATE(VERTS, NVERTS, FACES, NFACES, MESH_BOUNDS, DELTA)
! This routine reduces the size of a geometry by
!  1) merging vertices that are "close" together
!  2) eliminating redundent vertices
!  3) eliminating "singular" triangles

INTEGER, INTENT(INOUT) :: NVERTS, NFACES
REAL(EB), INTENT(INOUT), DIMENSION(3*NVERTS), TARGET :: VERTS
INTEGER, INTENT(INOUT), DIMENSION(3*NFACES), TARGET :: FACES
REAL(EB), INTENT(IN) :: DELTA
REAL(EB), DIMENSION(6), INTENT(IN) :: MESH_BOUNDS

INTEGER, PARAMETER :: V_MERGED=-1, V_DISCARD=0, V_ORIGINAL=1

INTEGER, DIMENSION(NVERTS) :: VERT_STATE, VERT_MAP, VERT_TYPE
REAL(EB), POINTER, DIMENSION(:) :: V1, V2, V3, VERTFROM, VERTTO
INTEGER, POINTER, DIMENSION(:) :: TRI_I, TRI_FROM, TRI_TO
INTEGER, DIMENSION(3) :: TRI_NEW
REAL(EB) :: D12, D13, D23
INTEGER :: I, IFROM, ITO, ITER, MAX_ITER
LOGICAL :: HAVE_SMALL
REAL(EB), DIMENSION(3) :: VAVG

HAVE_SMALL = .TRUE.
MAX_ITER = 4
ITER = 0
DO WHILE (HAVE_SMALL .AND. ITER<MAX_ITER) ! iterate until no further changes are made (or 10 times whichever comes first)
   HAVE_SMALL = .FALSE.

   ! VERT_STATE
   !    V_MERGE =   -1  -  merged vertex
   !    V_DISCARD =  0  -  discard vertex
   !    V_ORIGINAL = 1  -  vertex kept and not changed

   VERT_STATE(1:NVERTS) = V_ORIGINAL

   DO I = 1, NVERTS
     VERT_MAP(I) = I
   END DO

   ITER = ITER + 1

   CALL GET_VERTTYPE(NVERTS, FACES, NFACES, VERT_TYPE)

! combine vertices that are close together

   TRI_LOOP: DO I = 1, NFACES
      TRI_I=>FACES(3*I-2:3*I)
      V1=>VERTS(3*TRI_I(1)-2:3*TRI_I(1))
      V2=>VERTS(3*TRI_I(2)-2:3*TRI_I(2))
      V3=>VERTS(3*TRI_I(3)-2:3*TRI_I(3))

      IF (VERT_STATE(TRI_I(1))/=V_ORIGINAL .OR. &  ! only look at triangles that have not changed
          VERT_STATE(TRI_I(2))/=V_ORIGINAL .OR. &
          VERT_STATE(TRI_I(3))/=V_ORIGINAL) CYCLE TRI_LOOP

      D12 = DISTANCE3(V1,V2)
      D13 = DISTANCE3(V1,V3)
      D23 = DISTANCE3(V2,V3)
      IF (D12>DELTA .AND.D13>DELTA .AND.D23>DELTA) CYCLE TRI_LOOP ! triangle too large, do not combine verts

      HAVE_SMALL = .TRUE.
      IF (D12<DELTA .AND.D13>DELTA .AND.D23>DELTA) THEN ! combine verts 1 and 2 leave 3 alone
         VERT_STATE(TRI_I(1)) = V_MERGED
         VERT_STATE(TRI_I(2)) = V_DISCARD
         CALL AVERAGE_VERTS2(V1,VERT_TYPE(TRI_I(1)),V2,VERT_TYPE(TRI_I(2)),MESH_BOUNDS,VAVG)
         V1(1:3) = VAVG(1:3)
         VERT_MAP(TRI_I(2)) = TRI_I(1)
         TRI_I(2) = TRI_I(1)
      ELSE IF (D13<DELTA .AND.D12>DELTA .AND.D23>DELTA) THEN ! combine verts 1 and 3
         VERT_STATE(TRI_I(1)) = V_MERGED
         VERT_STATE(TRI_I(3)) = V_DISCARD
         CALL AVERAGE_VERTS2(V1,VERT_TYPE(TRI_I(1)),V3,VERT_TYPE(TRI_I(3)),MESH_BOUNDS,VAVG)
         V1(1:3) = VAVG(1:3)
         VERT_MAP(TRI_I(3)) = TRI_I(1)
         TRI_I(3) = TRI_I(1)
      ELSE IF (D23<DELTA .AND.D12>DELTA .AND.D13>DELTA) THEN ! combine verts 2 and 3
         VERT_STATE(TRI_I(2)) = V_MERGED
         VERT_STATE(TRI_I(3)) = V_DISCARD
         CALL AVERAGE_VERTS2(V2,VERT_TYPE(TRI_I(2)),V3,VERT_TYPE(TRI_I(3)),MESH_BOUNDS,VAVG)
         V2(1:3) = VAVG(1:3)
         VERT_MAP(TRI_I(3)) = TRI_I(2)
         TRI_I(3) = TRI_I(2)
      ELSE  ! combine verts 1, 2 and 3
         VERT_STATE(TRI_I(1))= V_MERGED
         VERT_STATE(TRI_I(2))= V_DISCARD
         VERT_STATE(TRI_I(3))= V_DISCARD
         CALL AVERAGE_VERTS3(V1, VERT_TYPE(TRI_I(1)), V2, VERT_TYPE(TRI_I(2)), V3, VERT_TYPE(TRI_I(3)), MESH_BOUNDS,VAVG)
         V1(1:3) = VAVG(1:3)
         VERT_MAP(TRI_I(2)) = TRI_I(1)
         VERT_MAP(TRI_I(3)) = TRI_I(1)
         TRI_I(2) = TRI_I(1)
         TRI_I(3) = TRI_I(1)
      ENDIF
   ENDDO TRI_LOOP

   ! remap triangle vertices

   DO I = 1, NFACES
      TRI_I=>FACES(3*I-2:3*I)
      TRI_I(1) = VERT_MAP(TRI_I(1))
      TRI_I(2) = VERT_MAP(TRI_I(2))
      TRI_I(3) = VERT_MAP(TRI_I(3))
   END DO

   ! construct new vertex list skipping over vertices that have been removed

   ITO = 0
   DO IFROM = 1, NVERTS
      IF (VERT_STATE(IFROM)/=V_DISCARD) THEN
        ITO = ITO + 1
        VERTFROM=>VERTS(3*IFROM-2:3*IFROM)
        VERTTO=>VERTS(3*ITO-2:3*ITO)
        VERTTO(1:3)=VERTFROM(1:3)
      ENDIF
      VERT_MAP(IFROM) = ITO
   ENDDO
   NVERTS = ITO

   ! eliminate singular triangles (as a result of merged vertices)

   ITO = 0
   DO IFROM = 1, NFACES
      TRI_FROM=>FACES(3*IFROM-2:3*IFROM)
      IF (TRI_FROM(1)/=TRI_FROM(2) .AND. TRI_FROM(1)/=TRI_FROM(3) .AND. TRI_FROM(2)/=TRI_FROM(3)) THEN
        TRI_NEW(1)=VERT_MAP(TRI_FROM(1))
        TRI_NEW(2)=VERT_MAP(TRI_FROM(2))
        TRI_NEW(3)=VERT_MAP(TRI_FROM(3))
        IF (TRI_NEW(1)/=0 .AND. TRI_NEW(2)/=0 .AND. TRI_NEW(3)/=0) THEN
          ITO=ITO+1
          TRI_TO=>FACES(3*ITO-2:3*ITO)
          TRI_TO(1)=TRI_NEW(1)
          TRI_TO(2)=TRI_NEW(2)
          TRI_TO(3)=TRI_NEW(3)
        ENDIF
      ENDIF
   ENDDO
   NFACES = ITO
ENDDO

END SUBROUTINE DECIMATE

!  ------------------ GEOMCLIP ------------------------

SUBROUTINE GEOMCLIP(VERTS, NVERTS, FACES, NFACES, XB)
INTEGER, INTENT(INOUT) :: NVERTS, NFACES
REAL(EB), INTENT(INOUT), DIMENSION(3*NVERTS), TARGET :: VERTS
INTEGER, INTENT(INOUT), DIMENSION(3*NFACES), TARGET :: FACES
REAL(EB), DIMENSION(6), INTENT(IN) :: XB

INTEGER, POINTER, DIMENSION(:) :: TRII
REAL(EB), POINTER, DIMENSION(:) :: V1, V2, V3
INTEGER, DIMENSION(NVERTS) :: HAVE_VERT, VERT_MAP
INTEGER, DIMENSION(NFACES) :: HAVE_TRIANGLE
INTEGER :: NT, NV

! find triangles that are not within the box defined by XB

HAVE_VERT(1:NVERTS) = 0
HAVE_TRIANGLE(1:NFACES) = 0
DO I = 1, NFACES
  TRII(1:3) => FACES(3*I-2:3*I)
  V1(1:3) => VERTS(3*TRII(1)-2:3*TRII(1))
  V2(1:3) => VERTS(3*TRII(2)-2:3*TRII(2))
  V3(1:3) => VERTS(3*TRII(3)-2:3*TRII(3))
  IF (IN_BOX3(XB,V1)==1 .OR. IN_BOX3(XB,V2)==1 .OR. IN_BOX3(XB,V3)==1) THEN
    HAVE_VERT(TRII(1)) = 1
    HAVE_VERT(TRII(2)) = 1
    HAVE_VERT(TRII(3)) = 1
    HAVE_TRIANGLE(I) = 1
  ENDIF
END DO

! construct a map pointing from old vertex index to new vertex index

NV = 0
DO I = 1, NVERTS
  IF (HAVE_VERT(I) == 1) THEN
    NV = NV + 1
    IF (I /= NV) THEN
       VERTS(3*NV-2:3*NV) = VERTS(3*I-2:3*I)
    ENDIF
  ENDIF
  VERT_MAP(I) = NV
END DO

! define new vertex indices

NT=0
DO I = 1, NFACES
  IF (HAVE_TRIANGLE(I)==1) THEN
    TRII(1:3) => FACES(3*I-2:3*I)
    TRII(1) = VERT_MAP(TRII(1))
    TRII(2) = VERT_MAP(TRII(2))
    TRII(3) = VERT_MAP(TRII(3))
    NT = NT + 1
    IF (NT /= I)FACES(3*NT-2:3*NT) = TRII(1:3)
  ENDIF
END DO

NFACES = NT
NVERTS = NV

END SUBROUTINE GEOMCLIP

END MODULE BOXTETRA_ROUTINES

