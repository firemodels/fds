MODULE SCARC_SOLVER
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
 
IMPLICIT NONE
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! communication structure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHARACTER (40) :: SCARC_MSG
INTEGER        :: SCARC_METHOD=0, SCARC_PRECON=4, SCARC_NIT=1000, SCARC_DEBUG=0
INTEGER        :: SCARC_COMPARE=0, SCARC_PERIODIC=0, SCARC_CASE=0
REAL (EB)      :: SCARC_EPS=1.E-16_EB, SCARC_EPS_REL=1.E-2_EB

INTEGER :: SCARC_COUNT2 = 0
 
INTEGER :: NumMaster, NumSlaves, NMESHES_SQRT
INTEGER :: myMasterId, mySlaveId, myself
INTEGER :: MpiGroupWorld, MpiGroupSlaves
INTEGER :: myMPI_COMM_WORLD, myMpi_COMM_SLAVES
INTEGER :: NumExclude, RankExclude (4)
 
INTEGER :: NCOM_TYPE_GLOBAL= 1, &
           NCOM_TYPE_LOCAL = 2
INTEGER :: NUPDATE_FULL    = 1, &
           NUPDATE_MEDIUM  = 2, &
           NUPDATE_VEL     = 3, &
           NUPDATE_FV      = 4, &
           NUPDATE_Q       = 5, &
           NUPDATE_H       = 6

INTEGER :: NEXCHANGE_EXTENDED = 1
INTEGER :: NNLEVEL = 10, NLEVEL, ILEVEL, NLEV1 = 1

INTEGER :: SNODE, RNODE, NDIAG
INTEGER :: NFACE, NEDGE, NVRTX
 
INTEGER :: N_REQ_MV, TAG_MV
INTEGER :: N_REQ_FACE, TAG_FACE
INTEGER :: N_REQ_EDGE, TAG_EDGE
INTEGER :: N_REQ_VRTX, TAG_VRTX
 
INTEGER :: IJKW_FACE (6, 12), IJKW_EDGE (12, 12), IJKW_VRTX (8, 12)
 
 
INTEGER, ALLOCATABLE, DIMENSION (:)    :: REQ_MV, REQ_FACE, REQ_EDGE, REQ_VRTX
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: NBR_FACE, TAGS_FACE
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: NBR_EDGE, TAGS_EDGE
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: NBR_VRTX, TAGS_VRTX
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: TAGS_MV, STATUSES_MV
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: NIC_FACE, STATUSES_FACE
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: NIC_EDGE, STATUSES_EDGE
INTEGER, ALLOCATABLE, DIMENSION (:, :) :: NIC_VRTX, STATUSES_VRTX
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scarc type on 'other' mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE OSCARC_TYPE
 
   REAL (EB), POINTER, DIMENSION (:) :: SENDBUF_MV, RECVBUF_MV
   REAL (EB), POINTER, DIMENSION (:) :: SENDBUF_FACE, RECVBUF_FACE
   REAL (EB), POINTER, DIMENSION (:) :: SENDBUF_EDGE, RECVBUF_EDGE
   REAL (EB), POINTER, DIMENSION (:) :: SENDBUF_VRTX, RECVBUF_VRTX
 
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Y_MV, R_MV, G_MV
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Z_FACE, Z_EDGE, Z_VRTX
 
END TYPE OSCARC_TYPE
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scarc type on 'global' mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE SCARC_TYPE
 
! some method parameter
   LOGICAL :: BFIRST, BPRECON, BGSADI, BFIRST_CG
 
   REAL (EB), POINTER, DIMENSION (:) :: SENDBUF_MV, SENDBUF_FACE, SENDBUF_EDGE, &
    SENDBUF_VRTX
   REAL (EB), POINTER, DIMENSION (:) :: RECVBUF_MV, RECVBUF_FACE, RECVBUF_EDGE, &
    RECVBUF_VRTX
 
! number of local and global grid cells
   INTEGER :: NCELLS_LOCAL, NCELLS_GLOBAL
   INTEGER, ALLOCATABLE, DIMENSION (:) :: NCELLS_LOCAL_ALL, NCELLS_GLOBAL_ALL
 
   INTEGER :: NPRECON
 
! local and global scalar products
   REAL (EB) :: ASUB
   REAL (EB) :: SP_LOCAL, SP_GLOBAL
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: SP_LOCAL_ALL, SP_GLOBAL_ALL
 
! full matrix and auxiliary vectors for global defect correction
   REAL (EB), POINTER, DIMENSION (:, :) :: AGLOB, ALOC
   REAL (EB), POINTER, DIMENSION (:, :, :) :: Z
   REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y, F, D, G, R, TMP
   REAL (EB), POINTER, DIMENSION (:) :: DIAG, LOW, UP, PERIOD
   REAL (EB), POINTER, DIMENSION (:, :) :: BXS0, BXF0, BYS0, BYF0, BZS0, BZF0
 
   REAL (EB) :: RES, RESIN, OMEGA, CAPPA
 
   REAL (EB) :: DXMIN, DYMIN, DZMIN
   REAL (EB) :: DXMAX, DYMAX, DZMAX
 
   REAL (EB) :: ALPHA, SIGMA0, SIGMA1, GAMMA
 
   TYPE (OSCARC_TYPE), POINTER, DIMENSION (:) :: OSCARC
 
END TYPE SCARC_TYPE
 
TYPE (SCARC_TYPE), SAVE, DIMENSION (:), ALLOCATABLE, TARGET :: SCARC
 
 
TYPE (SCARC_TYPE), POINTER :: S, S4
TYPE (OSCARC_TYPE), POINTER :: S2, S3
 
CONTAINS
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCARC : Initialize Scarc structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE (MYID)
   INTEGER :: MYID
 
   IF (SCARC_DEBUG .GE. 1) THEN
      SCARC_MSG = 'msg/  _scarc'
      WRITE (SCARC_MSG(5:6), '(i2.2)') MYID + 1
      OPEN (9, FILE=SCARC_MSG)
   ENDIF
 
   IF (TWO_D) THEN
      CALL SCARC_INITIALIZE_NEIGHBORS2D
      CALL SCARC_INITIALIZE_MATRICES2D (MYID+1)
      CALL SCARC_INITIALIZE_MESH_EXCHANGE2D (MYID+1)
      CALL SCARC_INITIALIZE_CG2D (MYID+1)
   ELSE
      CALL SCARC_INITIALIZE_NEIGHBORS3D
      CALL SCARC_INITIALIZE_MATRICES3D (MYID+1)
      CALL SCARC_INITIALIZE_MESH_EXCHANGE3D (MYID+1)
      CALL SCARC_INITIALIZE_CG3D (MYID+1)
   ENDIF
!ALLOCATE(BCOPY(0:IBAR+1,0:JBAR+1,0:KBAR+1))
 
END SUBROUTINE SCARC_INITIALIZE
 
 
!------------------------------------------------------------------------------------
! determine neighborship/communication structure for data exchange for 2D-case
!------------------------------------------------------------------------------------
SUBROUTINE SCARC_INITIALIZE_NEIGHBORS2D
 
INTEGER :: NM, NOM, III
INTEGER :: INBR1, INBR2, INBR3, JNBR1, JNBR2
INTEGER :: IFACE1, IFACE2, JFACE1, JFACE2
INTEGER :: IEDGE, IEDGE1, IEDGE2, IVRTX
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER, DIMENSION (2, 12) :: FACE_MAP = reshape ( (/ 1, 2, 1, 3, 1, 4, 1, 5, 2, &
 5, 2, 3, 3, 4, 4, 5, 2, 6, 3, 6, 4, 6, 5, 6 /), (/ 2, 12 /))
INTEGER, DIMENSION (3, 8) :: VRTX_MAP = reshape ( (/ 4, 1, 1, 1, 2, 1, 2, 3, 1, &
 3, 4, 1, 12, 9, 6, 9, 10, 6, 10, 11, 6, 11, 12, 6 /), (/ 3, 8 /))
LOGICAL FOUND
TYPE (MESH_TYPE), POINTER :: M, M1
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'Starting SCARC_INITIALIZE_NEIGHBORS2D'
   WRITE (9,*) 'NMESHES=', NMESHES
   WRITE (9,*) 'NMESHES=', NMESHES
   CALL flush (9)
ENDIF
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  2D: find adjacent neighbors via common edges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NFACE = 4
ALLOCATE (NBR_FACE(NMESHES, NFACE))
NBR_FACE = 0
 
MESH_LOOP1: DO NM = 1, NMESHES
   M => MESHES (NM)
   OTHER_MESH_LOOP1: DO NOM = 1, NMESHES
 
      IF (NIC(NOM, NM) == 0) CYCLE OTHER_MESH_LOOP1
 
      IF (SCARC_DEBUG .GE. 2) THEN
         WRITE (9,*) '----->  NM=', NM, ': NOM=', NOM
         WRITE (9,*) 'NIC(', NOM, ',', NM, ')=', NIC (NOM, NM)
      ENDIF
 
      M1 => MESHES (NOM)
 
      IF (MATCH(M%ZS, M1%ZF)) THEN ! ZS: lower edge in z-direction
         IF (MATCH(M%XS, M1%XS) .AND. MATCH(M%XF, M1%XF)) NBR_FACE (NM, 1) = NOM
      ELSE IF (MATCH(M%XF, M1%XS)) THEN ! XF: right edge in x-direction
         IF (MATCH(M%ZS, M1%ZS) .AND. MATCH(M%ZF, M1%ZF)) NBR_FACE (NM, 2) = NOM
      ELSE IF (MATCH(M%ZF, M1%ZS)) THEN ! ZF: upper edge in z-direction
         IF (MATCH(M%XS, M1%XS) .AND. MATCH(M%XF, M1%XF)) NBR_FACE (NM, 3) = NOM
      ELSE IF (MATCH(M%XS, M1%XF)) THEN ! XS: left edge in x-direction
         IF (MATCH(M%ZS, M1%ZS) .AND. MATCH(M%ZF, M1%ZF)) NBR_FACE (NM, 4) = NOM
      ENDIF
 
   ENDDO OTHER_MESH_LOOP1
ENDDO MESH_LOOP1
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) '============================================================='
   WRITE (9,*) 'NBR_FACE:'
   WRITE (9, '(4i6)') ((NBR_FACE(NM, III), III=1, NFACE), NM=1, NMESHES)
ENDIF
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine diagonal neighbours in 2D: only along vertices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NVRTX = 4
ALLOCATE (NBR_VRTX(NMESHES, NVRTX))
NBR_VRTX = 0
 
 
MESH_LOOP2: DO NM = 1, NMESHES
 
   IF (SCARC_DEBUG .GE. 2) write (9,*) '======================= LOOKING AT ', NM
 
   VRTX_LOOP2D: DO IFACE1 = 1, NFACE
 
!
! get both neighbors at two joining edges of NM and determine
! their common neighbor which is a diagonal neighbor of NM
!
      IFACE2 = Mod (IFACE1, 4) + 1
 
      INBR1 = NBR_FACE (NM, IFACE1)
      INBR2 = NBR_FACE (NM, IFACE2)
 
      IF (INBR1 == 0 .OR. INBR2 == 0) CYCLE VRTX_LOOP2D
 
      IF (SCARC_DEBUG .GT. 2) THEN
         WRITE (9,*) 'IFACE1=', IFACE1
         WRITE (9,*) 'IFACE2=', IFACE2
      ENDIF
 
      NBR1_FACE_LOOP2D: DO JFACE1 = 1, NFACE
 
         JNBR1 = NBR_FACE (INBR1, JFACE1)
         IF (JNBR1 == 0 .OR. JNBR1 == NM) CYCLE NBR1_FACE_LOOP2D
 
         NBR2_FACE_LOOP2D: DO JFACE2 = 1, NFACE
            JNBR2 = NBR_FACE (INBR2, JFACE2)
            IF (JNBR2 == 0 .OR. JNBR2 == NM) CYCLE NBR2_FACE_LOOP2D
 
            IF (SCARC_DEBUG .GT. 2) THEN
               WRITE (9,*) 'INBR1=', INBR1, ':  JFACE1=', JFACE1, ':  JNBR1 =', &
                JNBR1
               WRITE (9,*) 'INBR2=', INBR2, ':  JFACE2=', JFACE2, ':  JNBR2 =', &
                JNBR2
            ENDIF
 
            IF (JNBR1 == JNBR2) THEN
               NBR_VRTX (NM, IFACE2) = JNBR1
               IF (SCARC_DEBUG .GE. 6) write (9,*) 'NBR_VRTX(', NM, ',', IFACE2, &
                ')=', NBR_VRTX (NM, IFACE2)
               EXIT NBR1_FACE_LOOP2D
            ENDIF
 
         ENDDO NBR2_FACE_LOOP2D
      ENDDO NBR1_FACE_LOOP2D
 
   ENDDO VRTX_LOOP2D
 
ENDDO MESH_LOOP2
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) '=========================================='
   WRITE (9,*) '=========================================='
   WRITE (9,*) '=========================================='
   WRITE (9,*) 'NBR_FACE:'
   WRITE (9, '(4i4)') ((NBR_FACE(NM, III), III=1, NFACE), NM=1, NMESHES)
   WRITE (9,*) 'NBR_VRTX:'
   WRITE (9, '(4i4)') ((NBR_VRTX(NM, III), III=1, NVRTX), NM=1, NMESHES)
ENDIF
 
END SUBROUTINE SCARC_INITIALIZE_NEIGHBORS2D
 
 
 
!------------------------------------------------------------------------------------
! determine neighborship/communication structure for data exchange for 3D case
!------------------------------------------------------------------------------------
SUBROUTINE SCARC_INITIALIZE_NEIGHBORS3D
 
INTEGER :: NM, NOM, III
INTEGER :: INBR1, INBR2, INBR3, JNBR1, JNBR2
INTEGER :: IFACE1, IFACE2, JFACE1, JFACE2
INTEGER :: IEDGE, IEDGE1, IEDGE2, IVRTX
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER, DIMENSION (2, 12) :: FACE_MAP = reshape ( (/ 1, 2, 1, 3, 1, 4, 1, 5, 2, &
 5, 2, 3, 3, 4, 4, 5, 2, 6, 3, 6, 4, 6, 5, 6 /), (/ 2, 12 /))
INTEGER, DIMENSION (3, 8) :: VRTX_MAP = reshape ( (/ 4, 1, 1, 1, 2, 1, 2, 3, 1, &
 3, 4, 1, 12, 9, 6, 9, 10, 6, 10, 11, 6, 11, 12, 6 /), (/ 3, 8 /))
LOGICAL FOUND
TYPE (MESH_TYPE), POINTER :: M, M1
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'Starting SCARC_INITIALIZE_NEIGHBORS3D'
   WRITE (9,*) 'NMESHES=', NMESHES
   WRITE (9,*) 'NMESHES=', NMESHES
ENDIF
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 3D: find adjacent neighbors via common faces and store them in NBR_FACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NFACE = 6
ALLOCATE (NBR_FACE(NMESHES, NFACE))
NBR_FACE = 0
 
MESH_LOOP1: DO NM = 1, NMESHES
   M => MESHES (NM)
   OTHER_MESH_LOOP1: DO NOM = 1, NMESHES
 
      IF (NIC(NOM, NM) == 0) CYCLE OTHER_MESH_LOOP1
 
      IF (SCARC_DEBUG .GE. 2) THEN
         WRITE (9,*) '----->  NM=', NM, ': NOM=', NOM
         WRITE (9,*) 'NIC(', NOM, ',', NM, ')=', NIC (NOM, NM)
      ENDIF
 
      M1 => MESHES (NOM)
 
      IF (SCARC_DEBUG .GE. 2) write (9, '(a,i2,a,6i6)') '0:NBR_FACE(', NM, '):', &
       (NBR_FACE(NM, III), III=1, 6)
 
      IF (MATCH(M%ZS, M1%ZF)) THEN ! ZS: lower face in z-direction
         IF (MATCH(M%XS, M1%XS) .AND. MATCH(M%XF, M1%XF) .AND. MATCH(M%YS, &
          M1%YS) .AND. MATCH(M%YF, M1%YF)) NBR_FACE (NM, 1) = NOM
         IF (SCARC_DEBUG .GE. 6) write (9,*) '1:NBR_FACE(', NM, ',', 3, ')=', &
          NBR_FACE (NM, 1)
      ELSE IF (MATCH(M%YS, M1%YF)) THEN ! YS: front face in y-direction
         IF (MATCH(M%XS, M1%XS) .AND. MATCH(M%XF, M1%XF) .AND. MATCH(M%ZS, &
          M1%ZS) .AND. MATCH(M%ZF, M1%ZF)) NBR_FACE (NM, 2) = NOM
         IF (SCARC_DEBUG .GE. 6) write (9,*) '2:NBR_FACE(', NM, ',', 3, ')=', &
          NBR_FACE (NM, 2)
      ELSE IF (MATCH(M%XF, M1%XS)) THEN ! XF: right face in x-direction
         IF (MATCH(M%YS, M1%YS) .AND. MATCH(M%YF, M1%YF) .AND. MATCH(M%ZS, &
          M1%ZS) .AND. MATCH(M%ZF, M1%ZF)) NBR_FACE (NM, 3) = NOM
         IF (SCARC_DEBUG .GE. 6) write (9,*) '3.2:NBR_FACE(', NM, ',', 3, ')=', &
          NBR_FACE (NM, 3)
      ELSE IF (MATCH(M%YF, M1%YS)) THEN ! YF: back face in y-direction
         IF (MATCH(M%XS, M1%XS) .AND. MATCH(M%XF, M1%XF) .AND. MATCH(M%ZS, &
          M1%ZS) .AND. MATCH(M%ZF, M1%ZF)) NBR_FACE (NM, 4) = NOM
         IF (SCARC_DEBUG .GE. 6) write (9,*) '4:NBR_FACE(', NM, ',', 3, ')=', &
          NBR_FACE (NM, 4)
      ELSE IF (MATCH(M%XS, M1%XF)) THEN ! XS: left face in x-direction
         IF (MATCH(M%YS, M1%YS) .AND. MATCH(M%YF, M1%YF) .AND. MATCH(M%ZS, &
          M1%ZS) .AND. MATCH(M%ZF, M1%ZF)) NBR_FACE (NM, 5) = NOM
         IF (SCARC_DEBUG .GE. 6) write (9,*) '5:NBR_FACE(', NM, ',', 3, ')=', &
          NBR_FACE (NM, 5)
      ELSE IF (MATCH(M%ZF, M1%ZS)) THEN ! ZF: upper face in z-direction
         IF (MATCH(M%XS, M1%XS) .AND. MATCH(M%XF, M1%XF) .AND. MATCH(M%YS, &
          M1%YS) .AND. MATCH(M%YF, M1%YF)) NBR_FACE (NM, 6) = NOM
         IF (SCARC_DEBUG .GE. 6) write (9,*) '6:NBR_FACE(', NM, ',', 3, ')=', &
          NBR_FACE (NM, 6)
      ENDIF
 
      IF (SCARC_DEBUG .GE. 2) write (9, '(a,i2,a,6i6)') '1:NBR_FACE(', NM, '):', &
       (NBR_FACE(NM, III), III=1, 6)
 
 
   ENDDO OTHER_MESH_LOOP1
ENDDO MESH_LOOP1
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) '=========================================='
   WRITE (9,*) '=========================================='
   WRITE (9,*) '=========================================='
   WRITE (9,*) 'NBR_FACE:'
   WRITE (9, '(6i6)') ((NBR_FACE(NM, III), III=1, NFACE), NM=1, NMESHES)
ENDIF
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine diagonal neighbours in 3D along edges and vertices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NEDGE = 12
ALLOCATE (NBR_EDGE(NMESHES, NEDGE))
NBR_EDGE = 0
 
NVRTX = 8
ALLOCATE (NBR_VRTX(NMESHES, NVRTX))
NBR_VRTX = 0
 
 
MESH_LOOP2: DO NM = 1, NMESHES
   IF (SCARC_DEBUG .GE. 2) THEN
      WRITE (9,*) 'FACE_MAP(:,1) =', FACE_MAP (1, 1), FACE_MAP (2, 1)
      WRITE (9,*) 'FACE_MAP(:,2) =', FACE_MAP (1, 2), FACE_MAP (2, 2)
      WRITE (9,*) 'FACE_MAP(:,3) =', FACE_MAP (1, 3), FACE_MAP (2, 3)
      WRITE (9,*) 'FACE_MAP(:,4) =', FACE_MAP (1, 4), FACE_MAP (2, 4)
      WRITE (9,*) 'FACE_MAP(:,5) =', FACE_MAP (1, 5), FACE_MAP (2, 5)
      WRITE (9,*) 'FACE_MAP(:,6) =', FACE_MAP (1, 6), FACE_MAP (2, 6)
      WRITE (9,*) 'FACE_MAP(:,7) =', FACE_MAP (1, 7), FACE_MAP (2, 7)
      WRITE (9,*) 'FACE_MAP(:,8) =', FACE_MAP (1, 8), FACE_MAP (2, 8)
      WRITE (9,*) 'FACE_MAP(:,9) =', FACE_MAP (1, 9), FACE_MAP (2, 9)
      WRITE (9,*) 'FACE_MAP(:,10)=', FACE_MAP (1, 10), FACE_MAP (2, 10)
      WRITE (9,*) 'FACE_MAP(:,11)=', FACE_MAP (1, 11), FACE_MAP (2, 11)
      WRITE (9,*) 'FACE_MAP(:,12)=', FACE_MAP (1, 12), FACE_MAP (2, 12)
      WRITE (9,*)
      WRITE (9,*) 'VRTX_MAP(:,1) =', VRTX_MAP (1, 1), VRTX_MAP (2, 1)
      WRITE (9,*) 'VRTX_MAP(:,2) =', VRTX_MAP (1, 2), VRTX_MAP (2, 2)
      WRITE (9,*) 'VRTX_MAP(:,3) =', VRTX_MAP (1, 3), VRTX_MAP (2, 3)
      WRITE (9,*) 'VRTX_MAP(:,4) =', VRTX_MAP (1, 4), VRTX_MAP (2, 4)
      WRITE (9,*) 'VRTX_MAP(:,5) =', VRTX_MAP (1, 5), VRTX_MAP (2, 5)
      WRITE (9,*) 'VRTX_MAP(:,6) =', VRTX_MAP (1, 6), VRTX_MAP (2, 6)
      WRITE (9,*) 'VRTX_MAP(:,7) =', VRTX_MAP (1, 7), VRTX_MAP (2, 7)
      WRITE (9,*) 'VRTX_MAP(:,8) =', VRTX_MAP (1, 8), VRTX_MAP (2, 8)
      WRITE (9,*)
      WRITE (9,*) '======================= LOOKING AT ', NM
   ENDIF
 
!
! get both neighbors at two joining edges of NM and determine
! their common neighbor which is a diagonal neighbor of NM
!
   EDGE_LOOP3D: DO IEDGE = 1, NEDGE
 
      IFACE1 = FACE_MAP (1, IEDGE)
      IFACE2 = FACE_MAP (2, IEDGE)
 
      INBR1 = NBR_FACE (NM, IFACE1)
      INBR2 = NBR_FACE (NM, IFACE2)
 
      IF (INBR1 == 0 .OR. INBR2 == 0) CYCLE EDGE_LOOP3D
 
      IF (SCARC_DEBUG .GE. 6) THEN
         WRITE (9,*) '---------IEDGE=', IEDGE
         WRITE (9,*) 'IFACE1=', IFACE1
         WRITE (9,*) 'IFACE2=', IFACE2
      ENDIF
 
      FACE_NBR1_LOOP3D: DO JFACE1 = 1, NFACE
 
         JNBR1 = NBR_FACE (INBR1, JFACE1)
         IF (JNBR1 == 0 .OR. JNBR1 == NM) CYCLE FACE_NBR1_LOOP3D
 
         FACE_NBR2_LOOP3D: DO JFACE2 = 1, NFACE
 
            JNBR2 = NBR_FACE (INBR2, JFACE2)
            IF (JNBR2 == 0 .OR. JNBR2 == NM) CYCLE FACE_NBR2_LOOP3D
 
            IF (SCARC_DEBUG .GE. 6) write (9,*) 'INBR1=', INBR1, ':  JFACE1=', &
             JFACE1, ':  JNBR1 =', JNBR1
            IF (SCARC_DEBUG .GE. 6) write (9,*) 'INBR2=', INBR2, ':  JFACE2=', &
             JFACE2, ':  JNBR2 =', JNBR2
 
            IF (JNBR1 == JNBR2) THEN
               NBR_EDGE (NM, IEDGE) = JNBR1
               IF (SCARC_DEBUG .GE. 6) write (9,*) 'NBR_EDGE(', NM, ',', IEDGE, &
                ')=', NBR_VRTX (NM, IEDGE)
               EXIT FACE_NBR1_LOOP3D
            ENDIF
 
         ENDDO FACE_NBR2_LOOP3D
      ENDDO FACE_NBR1_LOOP3D
 
   ENDDO EDGE_LOOP3D
!
! get both neighbors at two joining edges of NM and determine
! their common neighbor which is a diagonal neighbor of NM
!
   VRTX_LOOP3D: DO IVRTX = 1, NVRTX
 
      IEDGE1 = VRTX_MAP (1, IVRTX)
      IEDGE2 = VRTX_MAP (2, IVRTX)
      IFACE1 = VRTX_MAP (3, IVRTX)
 
      INBR1 = NBR_EDGE (NM, IEDGE1)
      INBR2 = NBR_EDGE (NM, IEDGE2)
      INBR3 = NBR_FACE (NM, IFACE1)
 
      IF (INBR1 == 0 .OR. INBR2 == 0) CYCLE VRTX_LOOP3D
 
      IF (SCARC_DEBUG .GE. 6) write (9,*) '=====================================&
     &==IVRTX=', IVRTX
      IF (SCARC_DEBUG .GE. 6) write (9,*) 'IEDGE1=', IEDGE1
      IF (SCARC_DEBUG .GE. 6) write (9,*) 'IEDGE2=', IEDGE2
 
      VRTX_NBR1_LOOP3D: DO JFACE1 = 1, NFACE
 
         JNBR1 = NBR_FACE (INBR1, JFACE1)
         IF (JNBR1 == 0 .OR. JNBR1 == NM .OR. JNBR1 == INBR3) CYCLE &
          VRTX_NBR1_LOOP3D
 
         VRTX_NBR2_LOOP3D: DO JFACE2 = 1, NFACE
 
            JNBR2 = NBR_FACE (INBR2, JFACE2)
            IF (JNBR2 == 0 .OR. JNBR2 == NM .OR. JNBR1 == INBR3) CYCLE &
             VRTX_NBR2_LOOP3D
 
            IF (SCARC_DEBUG .GE. 6) write (9,*) 'INBR1=', INBR1, ':  JFACE1=', &
             JFACE1, ':  JNBR1 =', JNBR1
            IF (SCARC_DEBUG .GE. 6) write (9,*) 'INBR2=', INBR2, ':  JFACE2=', &
             JFACE2, ':  JNBR2 =', JNBR2
 
            IF (JNBR1 == JNBR2) THEN
               NBR_VRTX (NM, IVRTX) = JNBR1
               IF (SCARC_DEBUG .GE. 6) write (9,*) 'NBR_VRTX(', NM, ',', IVRTX, &
                ')=', NBR_VRTX (NM, IVRTX)
               EXIT VRTX_NBR1_LOOP3D
            ENDIF
 
         ENDDO VRTX_NBR2_LOOP3D
      ENDDO VRTX_NBR1_LOOP3D
 
   ENDDO VRTX_LOOP3D
 
ENDDO MESH_LOOP2
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) '=========================================='
   WRITE (9,*) '=========================================='
   WRITE (9,*) '=========================================='
   WRITE (9,*) 'NBR_FACE:'
   WRITE (9, '(6i4)') ((NBR_FACE(NM, III), III=1, NFACE), NM=1, NMESHES)
   WRITE (9,*) 'NBR_EDGE:'
   WRITE (9, '(12i4)') ((NBR_EDGE(NM, III), III=1, NEDGE), NM=1, NMESHES)
   WRITE (9,*) 'NBR_VRTX:'
   WRITE (9, '( 8i4)') ((NBR_VRTX(NM, III), III=1, NVRTX), NM=1, NMESHES)
ENDIF
 
 
END SUBROUTINE SCARC_INITIALIZE_NEIGHBORS3D
 
 
 
 
 
 
!------------------------------------------------------------------------------------
! determine neighborship/communication structure for data exchange in 2D
!------------------------------------------------------------------------------------
SUBROUTINE SCARC_INITIALIZE_MESH_EXCHANGE2D (NM)
 
INTEGER :: IM, NM, NOM, III
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: IFACE, IEDGE, IVRTX
INTEGER :: I1, I2, J1, J2, K1, K2
INTEGER :: IERR = 0
LOGICAL FOUND
 
INCLUDE 'mpif.h'
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'Starting SCARC_INITIALIZE_MESH_EXCHANGE2D'
   WRITE (9,*) 'NM=', NM
   CALL flush (9)
ENDIF
 
ALLOCATE (STATUSES_MV(MPI_STATUS_SIZE, NMESHES*NMESHES*20))
ALLOCATE (STATUSES_FACE(MPI_STATUS_SIZE, NMESHES*NMESHES*20))
ALLOCATE (STATUSES_VRTX(MPI_STATUS_SIZE, NMESHES*NMESHES*20))
 
ALLOCATE (REQ_MV(NMESHES*NMESHES*20))
ALLOCATE (REQ_FACE(NMESHES*NMESHES*20))
ALLOCATE (REQ_VRTX(NMESHES*NMESHES*20))
REQ_MV = MPI_REQUEST_NULL
REQ_FACE = MPI_REQUEST_NULL
REQ_VRTX = MPI_REQUEST_NULL
 
ALLOCATE (TAGS_MV(NMESHES, NMESHES))
ALLOCATE (TAGS_FACE(NMESHES, NMESHES))
ALLOCATE (TAGS_VRTX(NMESHES, NMESHES))
TAGS_MV = 0
TAGS_FACE = 0
TAGS_VRTX = 0
 
ALLOCATE (NIC_FACE(NMESHES, NMESHES))
ALLOCATE (NIC_VRTX(NMESHES, NMESHES))
NIC_FACE = 0
NIC_VRTX = 0
 
IF (SCARC_DEBUG .GE. 2) write (9,*) 'After Allocate'
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Create unique tags arrays for face, edge and vertex exchanges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TAG_MV = 0
TAG_FACE = 0
TAG_VRTX = 0
 
MESH_LOOP1: DO IM = 1, NMESHES
   DO NOM = IM, NMESHES
 
      TAG_MV = TAG_MV + 1
      TAG_FACE = TAG_FACE + 1
      TAG_VRTX = TAG_VRTX + 1
 
      TAGS_MV (IM, NOM) = TAG_MV
      TAGS_MV (NOM, IM) = TAG_MV
 
      TAGS_FACE (IM, NOM) = TAG_FACE
      TAGS_FACE (NOM, IM) = TAG_FACE
 
      TAGS_VRTX (IM, NOM) = TAG_VRTX
      TAGS_VRTX (NOM, IM) = TAG_VRTX
 
   ENDDO
ENDDO MESH_LOOP1
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'TAGS_FACE:'
   WRITE (9, '(4i16)') ((TAGS_FACE(I1, J1), I1=1, NMESHES), J1=1, NMESHES)
ENDIF
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize structures for face neighbor communication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
S => SCARC (NM)
MV_NBR_LOOP: DO NOM = 1, NMESHES
 
   IF (NOM == NM) CYCLE MV_NBR_LOOP
   IF (NIC(NM, NOM) == 0 .AND. NIC(NOM, NM) == 0) CYCLE MV_NBR_LOOP
 
   IMIN = I_MIN (NOM, NM)
   IMAX = I_MAX (NOM, NM)
   JMIN = J_MIN (NOM, NM)
   JMAX = J_MAX (NOM, NM)
   KMIN = K_MIN (NOM, NM)
   KMAX = K_MAX (NOM, NM)
 
   IF (SCARC_DEBUG .GE. 2) THEN
      WRITE (9,*) '---NM=', NM
      WRITE (9,*) '---NOM=', NOM
      WRITE (9,*) 'IMIN=', IMIN
      WRITE (9,*) 'IMAX=', IMAX
      WRITE (9,*) 'JMIN=', JMIN
      WRITE (9,*) 'JMAX=', JMAX
      WRITE (9,*) 'KMIN=', KMIN
      WRITE (9,*) 'KMAX=', KMAX
      WRITE (9,*) 'allocate SCARC(', NM, ')%OSCARC(', NOM, ')'
      CALL flush (9)
   ENDIF
 
!JMIN=1
!JMAX=1
 
!ALLOCATE(S%OSCARC(NOM)%Y_MV(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE (S%OSCARC(NOM)%Y_MV(0:IBAR+1, 1, 0:KBAR+1))
   S%OSCARC(NOM)%Y_MV = 0.0_EB
 
!ALLOCATE(S%OSCARC(NOM)%R_MV(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE (S%OSCARC(NOM)%R_MV(0:IBAR+1, 1, 0:KBAR+1))
   S%OSCARC(NOM)%R_MV = 0.0_EB
 
!ALLOCATE(S%OSCARC(NOM)%G_MV(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE (S%OSCARC(NOM)%G_MV(0:IBAR+1, 1, 0:KBAR+1))
   S%OSCARC(NOM)%G_MV = 0.0_EB
 
ENDDO MV_NBR_LOOP
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize communication structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!
!!! Initialize structures for extended face neigbours communication
!!!
IF (NEXCHANGE_EXTENDED == 1) THEN ! exchange includes ghost cells
   IJKW_FACE (1, 1:12) = (/ 0, 1, 1, IBAR + 1, 1, 1, 0, 1, 0, IBAR + 1, 1, 0 /)
   IJKW_FACE (2, 1:12) = (/ IBAR, 1, 0, IBAR, 1, KBAR + 1, IBAR + 1, 1, 0, IBAR &
    + 1, 1, KBAR + 1 /)
   IJKW_FACE (3, 1:12) = (/ 0, 1, KBAR, IBAR + 1, 1, KBAR, 0, 1, KBAR + 1, IBAR &
    + 1, 1, KBAR + 1 /)
   IJKW_FACE (4, 1:12) = (/ 1, 1, 0, 1, 1, KBAR + 1, 0, 1, 0, 0, 1, KBAR + 1 /)
ELSE ! exchange excludes ghost cells
   IJKW_FACE (1, 1:12) = (/ 1, 1, 1, IBAR, 1, 1, 1, 1, 0, IBAR, 1, 0 /)
   IJKW_FACE (2, 1:12) = (/ IBAR, 1, 1, IBAR, 1, KBAR, IBAR + 1, 1, 1, IBAR + 1, &
    1, KBAR /)
   IJKW_FACE (3, 1:12) = (/ 1, 1, KBAR, IBAR, 1, KBAR, 1, 1, KBAR + 1, IBAR, 1, &
    KBAR + 1 /)
   IJKW_FACE (4, 1:12) = (/ 1, 1, 1, 1, 1, KBAR, 0, 1, 1, 0, 1, KBAR /)
ENDIF
 
FACE_INIT_LOOP2D: DO IFACE = 1, NFACE
 
   NOM = NBR_FACE (NM, IFACE)!  neighbor at edge IFACE
   IF (NOM == 0) CYCLE FACE_INIT_LOOP2D
 
   I1 = IJKW_FACE (IFACE, 1)
   I2 = IJKW_FACE (IFACE, 4)
   K1 = IJKW_FACE (IFACE, 3)
   K2 = IJKW_FACE (IFACE, 6)
 
   NIC_FACE (NM, NOM) = (I2-I1+1) * (K2-K1+1)
   NIC_FACE (NOM, NM) = NIC_FACE (NM, NOM)
 
   ALLOCATE (S%OSCARC(NOM)%Z_FACE(I1:I2, 1, K1:K2))
   S%OSCARC(NOM)%Z_FACE = 0.0_EB
 
   IF (SCARC_DEBUG .GE. 2) THEN
      WRITE (9,*) 'NIC_FACE(', NM, ',', NOM, ')=', NIC_FACE (NM, NOM)
      WRITE (9,*) 'NIC_FACE(', NOM, ',', NM, ')=', NIC_FACE (NOM, NM)
      WRITE (9,*) 'IFACE=', IEDGE, ': NOM=', NOM
      WRITE (9,*) 'IJKW_FACE:'
      WRITE (9, '(12i4)') (IJKW_FACE(IFACE, III), III=1, 12)
      WRITE (9,*) 'I1=', I1
      WRITE (9,*) 'I2=', I2
      WRITE (9,*) 'K1=', K1
      WRITE (9,*) 'K2=', K2
   ENDIF
ENDDO FACE_INIT_LOOP2D
 
 
!!!
!!! Initialize structures for diagonal-vertex neigbours communication
!!!
IJKW_VRTX (1, 1:12) = (/ 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0 /)
IJKW_VRTX (2, 1:12) = (/ IBAR, 1, 1, IBAR, 1, 1, IBAR + 1, 1, 0, IBAR + 1, 1, 0 &
 /)
IJKW_VRTX (3, 1:12) = (/ IBAR, 1, KBAR, IBAR, 1, KBAR, IBAR + 1, 1, KBAR + 1, &
 IBAR + 1, 1, KBAR + 1 /)
IJKW_VRTX (4, 1:12) = (/ 1, 1, KBAR, 1, 1, KBAR, 0, 1, KBAR + 1, 0, 1, KBAR + 1 &
 /)
 
VRTX_INIT_LOOP2D: DO IVRTX = 1, NVRTX
 
   NOM = NBR_VRTX (NM, IVRTX)!  diagonal neighbor at vertex IVRTX
   IF (NOM == 0) CYCLE VRTX_INIT_LOOP2D
 
   NIC_VRTX (NM, NOM) = 1
   NIC_VRTX (NOM, NM) = 1
   IF (SCARC_DEBUG .GE. 2) write (9,*) 'IVRTX=', IVRTX, ': NIC_VRTX(', NOM, ',', &
    NM, ')=', NIC_VRTX (NM, NOM)
 
   ALLOCATE (S%OSCARC(NOM)%Z_VRTX(1, 1, 1))
 
ENDDO VRTX_INIT_LOOP2D
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'IJKW_FACE:'
   DO IFACE = 1, NFACE
      WRITE (9, '(12i4)') (IJKW_FACE(IFACE, III), III=1, 12)
   ENDDO
   WRITE (9,*) 'IJKW_VRTX:'
   DO IVRTX = 1, NVRTX
      WRITE (9, '(12i4)') (IJKW_VRTX(IVRTX, III), III=1, 12)
   ENDDO
   CALL flush (9)
ENDIF
 
 
CALL SCARC_RECEIVE (0, 0)
CALL SCARC_EXCHANGE (0, 0)
CALL flush (9)
 
IF (SCARC_DEBUG .GE. 2) write (9,*) 'Leaving SCARC_INITIALIZE_MESH_EXCHANGE2D'
 
 
END SUBROUTINE SCARC_INITIALIZE_MESH_EXCHANGE2D
 
 
 
 
!------------------------------------------------------------------------------------
! determine neighborship/communication structure for data exchange in 3D
!------------------------------------------------------------------------------------
SUBROUTINE SCARC_INITIALIZE_MESH_EXCHANGE3D (NM)
 
INTEGER :: IM, NM, NOM, III
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: IFACE, IEDGE, IVRTX
INTEGER :: I1, I2, J1, J2, K1, K2
INTEGER :: IERR = 0
!INTEGER, DIMENSION(:):: DISPLS, DISPLS2D, COUNTS2D
LOGICAL FOUND
TYPE (MESH_TYPE), POINTER :: M, M1
 
INCLUDE 'mpif.h'
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'Starting SCARC_INITIALIZE_MESH_EXCHANGE3D'
   WRITE (9,*) 'NM=', NM
ENDIF
CALL flush (9)
 
ALLOCATE (REQ_MV(NMESHES*NMESHES*20))
REQ_MV = MPI_REQUEST_NULL
 
ALLOCATE (REQ_FACE(NMESHES*NMESHES*20))
REQ_FACE = MPI_REQUEST_NULL
 
ALLOCATE (REQ_EDGE(NMESHES*NMESHES*20))
REQ_EDGE = MPI_REQUEST_NULL
 
ALLOCATE (REQ_VRTX(NMESHES*NMESHES*20))
REQ_VRTX = MPI_REQUEST_NULL
 
ALLOCATE (STATUSES_MV(MPI_STATUS_SIZE, NMESHES*NMESHES*20))
ALLOCATE (STATUSES_FACE(MPI_STATUS_SIZE, NMESHES*NMESHES*20))
ALLOCATE (STATUSES_EDGE(MPI_STATUS_SIZE, NMESHES*NMESHES*20))
ALLOCATE (STATUSES_VRTX(MPI_STATUS_SIZE, NMESHES*NMESHES*20))
 
ALLOCATE (TAGS_MV(NMESHES, NMESHES))
ALLOCATE (TAGS_FACE(NMESHES, NMESHES))
ALLOCATE (TAGS_EDGE(NMESHES, NMESHES))
ALLOCATE (TAGS_VRTX(NMESHES, NMESHES))
 
TAGS_MV = 0
TAGS_FACE = 0
TAGS_EDGE = 0
TAGS_VRTX = 0
 
ALLOCATE (NIC_FACE(NMESHES, NMESHES))
NIC_FACE = 0
 
ALLOCATE (NIC_EDGE(NMESHES, NMESHES))
NIC_EDGE = 0
 
ALLOCATE (NIC_VRTX(NMESHES, NMESHES))
NIC_VRTX = 0
 
IF (SCARC_DEBUG .GE. 2) write (9,*) 'After Allocate'
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Create unique tags arrays for face, edge and vertex exchanges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TAG_MV = 0
TAG_FACE = 0
TAG_EDGE = 0
TAG_VRTX = 0
 
MESH_LOOP1: DO IM = 1, NMESHES
   DO NOM = IM, NMESHES
 
      TAG_MV = TAG_MV + 1
      TAG_FACE = TAG_FACE + 1
      TAG_EDGE = TAG_EDGE + 1
      TAG_VRTX = TAG_VRTX + 1
 
      TAGS_MV (IM, NOM) = TAG_MV
      TAGS_MV (NOM, IM) = TAG_MV
 
      TAGS_FACE (IM, NOM) = TAG_FACE
      TAGS_FACE (NOM, IM) = TAG_FACE
 
      TAGS_EDGE (IM, NOM) = TAG_EDGE
      TAGS_EDGE (NOM, IM) = TAG_EDGE
 
      TAGS_VRTX (IM, NOM) = TAG_VRTX
      TAGS_VRTX (NOM, IM) = TAG_VRTX
 
      IF (SCARC_DEBUG .GE. 2) write (9,*) 'TAGS_FACE(', IM, ',', NOM, ')=', &
       TAGS_FACE (IM, NOM)
      CALL flush (9)
 
   ENDDO
ENDDO MESH_LOOP1
 
 
 
S => SCARC (NM)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize structures for face neighbor communication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FACE_NBR_LOOP: DO NOM = 1, NMESHES
 
   IF (NOM == NM) CYCLE FACE_NBR_LOOP
   IF (NIC(NM, NOM) == 0 .AND. NIC(NOM, NM) == 0) CYCLE FACE_NBR_LOOP
 
   IMIN = I_MIN (NOM, NM)
   IMAX = I_MAX (NOM, NM)
   JMIN = J_MIN (NOM, NM)
   JMAX = J_MAX (NOM, NM)
   KMIN = K_MIN (NOM, NM)
   KMAX = K_MAX (NOM, NM)
 
   IF (SCARC_DEBUG .GE. 2) THEN
      WRITE (9,*) '---NM=', NM
      WRITE (9,*) '---NOM=', NOM
      WRITE (9,*) 'IMIN=', IMIN
      WRITE (9,*) 'IMAX=', IMAX
      WRITE (9,*) 'JMIN=', JMIN
      WRITE (9,*) 'JMAX=', JMAX
      WRITE (9,*) 'KMIN=', KMIN
      WRITE (9,*) 'KMAX=', KMAX
      WRITE (9,*) 'allocate SCARC(', NM, ')%OSCARC(', NOM, ')'
      CALL flush (9)
   ENDIF
 
!ALLOCATE(S%OSCARC(NOM)%Y_MV(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE (S%OSCARC(NOM)%Y_MV(0:IBP1, 0:JBP1, 0:KBP1))
   S%OSCARC(NOM)%Y_MV = 0.0_EB
 
!ALLOCATE(S%OSCARC(NOM)%R_MV(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE (S%OSCARC(NOM)%R_MV(0:IBP1, 0:JBP1, 0:KBP1))
   S%OSCARC(NOM)%R_MV = 0.0_EB
 
!ALLOCATE(S%OSCARC(NOM)%G_MV(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))
   ALLOCATE (S%OSCARC(NOM)%G_MV(0:IBP1, 0:JBP1, 0:KBP1))
   S%OSCARC(NOM)%G_MV = 0.0_EB
 
ENDDO FACE_NBR_LOOP
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize communication structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!
!!! Initialize structures for extended face neigbours communication
!!!
IF (NEXCHANGE_EXTENDED == 1) THEN ! exchange includes ghost cells
   IJKW_FACE (1, 1:12) = (/ 0, 0, 1, IBAR + 1, JBAR + 1, 1, 0, 0, 0, IBAR + 1, &
    JBAR + 1, 0 /)
   IJKW_FACE (2, 1:12) = (/ 0, 1, 0, IBAR + 1, 1, KBAR + 1, 0, 0, 0, IBAR + 1, &
    0, KBAR + 1 /)
   IJKW_FACE (3, 1:12) = (/ IBAR, 0, 0, IBAR, JBAR + 1, KBAR + 1, IBAR + 1, 0, &
    0, IBAR + 1, JBAR + 1, KBAR + 1 /)
   IJKW_FACE (4, 1:12) = (/ 0, JBAR, 0, IBAR + 1, JBAR, KBAR + 1, 0, JBAR + 1, &
    0, IBAR + 1, JBAR + 1, KBAR + 1 /)
   IJKW_FACE (5, 1:12) = (/ 1, 0, 0, 1, JBAR + 1, KBAR + 1, 0, 0, 0, 0, JBAR + &
    1, KBAR + 1 /)
   IJKW_FACE (6, 1:12) = (/ 0, 0, KBAR, IBAR + 1, JBAR + 1, KBAR, 0, 0, KBAR + &
    1, IBAR + 1, JBAR + 1, KBAR + 1 /)
ELSE ! exchange excludes ghost cells
   IJKW_FACE (1, 1:12) = (/ 1, 1, 1, IBAR, JBAR, 1, 1, 1, 0, IBAR, JBAR, 0 /)
   IJKW_FACE (2, 1:12) = (/ 1, 1, 1, IBAR, 1, KBAR, 1, 0, 1, IBAR, 0, KBAR /)
   IJKW_FACE (3, 1:12) = (/ IBAR, 1, 1, IBAR, JBAR, KBAR, IBAR + 1, 1, 1, IBAR + &
    1, JBAR, KBAR /)
   IJKW_FACE (4, 1:12) = (/ 1, JBAR, 1, IBAR, JBAR, KBAR, 1, JBAR + 1, 1, IBAR, &
    JBAR + 1, KBAR /)
   IJKW_FACE (5, 1:12) = (/ 1, 1, 1, 1, JBAR, KBAR, 0, 1, 1, 0, JBAR, KBAR /)
   IJKW_FACE (6, 1:12) = (/ 1, 1, KBAR, IBAR, JBAR, KBAR, 1, 1, KBAR + 1, IBAR, &
    JBAR, KBAR + 1 /)
ENDIF
 
FACE_NBR_LOOP2: DO IFACE = 1, NFACE
 
   NOM = NBR_FACE (NM, IFACE)!  diagonal neighbor at edge IFACE
   IF (NOM == 0) CYCLE FACE_NBR_LOOP2
 
   I1 = IJKW_FACE (IFACE, 1)
   I2 = IJKW_FACE (IFACE, 4)
   J1 = IJKW_FACE (IFACE, 2)
   J2 = IJKW_FACE (IFACE, 5)
   K1 = IJKW_FACE (IFACE, 3)
   K2 = IJKW_FACE (IFACE, 6)
 
   NIC_FACE (NM, NOM) = (I2-I1+1) * (J2-J1+1) * (K2-K1+1)
   NIC_FACE (NOM, NM) = NIC_FACE (NM, NOM)
 
   ALLOCATE (S%OSCARC(NOM)%Z_FACE(I1:I2, J1:J2, K1:K2))
 
   IF (SCARC_DEBUG .GE. 2) THEN
      WRITE (9,*) 'NIC_FACE(', NOM, ',', NM, ')=', NIC_FACE (NM, NOM)
      WRITE (9,*) 'IFACE=', IEDGE, ': NOM=', NOM
      WRITE (9,*) 'IJKW_FACE:'
      WRITE (9, '(12i4)') (IJKW_FACE(IFACE, III), III=1, 12)
      WRITE (9,*) 'I1=', I1
      WRITE (9,*) 'I2=', I2
      WRITE (9,*) 'J1=', J1
      WRITE (9,*) 'J2=', J2
      WRITE (9,*) 'K1=', K1
      WRITE (9,*) 'K2=', K2
      CALL flush (9)
   ENDIF
 
ENDDO FACE_NBR_LOOP2
 
 
!!!
!!! Initialize structures for diagonal-edge neigbours communication
!!!
IF (NEXCHANGE_EXTENDED == 1) THEN ! exchange includes ghost cells
   IJKW_EDGE (1, 1:12) = (/ 0, 1, 1, IBAR + 1, 1, 1, 0, 0, 0, IBAR + 1, 0, 0 /)
   IJKW_EDGE (2, 1:12) = (/ IBAR, 0, 1, IBAR, JBAR + 1, 1, IBAR + 1, 0, 0, IBAR &
    + 1, JBAR + 1, 0 /)
   IJKW_EDGE (3, 1:12) = (/ 0, JBAR, 1, IBAR + 1, JBAR, 1, 0, JBAR + 1, 0, IBAR &
    + 1, JBAR + 1, 0 /)
   IJKW_EDGE (4, 1:12) = (/ 1, 0, 1, 1, JBAR + 1, 1, 0, 0, 0, 0, JBAR + 1, 0 /)
   IJKW_EDGE (5, 1:12) = (/ 1, 1, 0, 1, 1, KBAR + 1, 0, 0, 0, 0, 0, KBAR + 1 /)
   IJKW_EDGE (6, 1:12) = (/ IBAR, 1, 0, IBAR, 1, KBAR + 1, IBAR + 1, 0, 0, IBAR &
    + 1, 0, KBAR + 1 /)
   IJKW_EDGE (7, 1:12) = (/ IBAR, JBAR, 0, IBAR, JBAR, KBAR + 1, IBAR + 1, JBAR &
    + 1, 0, IBAR + 1, JBAR + 1, KBAR + 1 /)
   IJKW_EDGE (8, 1:12) = (/ 1, JBAR, 0, 1, JBAR, KBAR + 1, 0, JBAR + 1, 0, 0, &
    JBAR + 1, KBAR + 1 /)
   IJKW_EDGE (9, 1:12) = (/ 0, 1, KBAR, IBAR + 1, 1, KBAR, 0, 0, KBAR + 1, IBAR &
    + 1, 0, KBAR + 1 /)
   IJKW_EDGE (10, 1:12) = (/ IBAR, 0, KBAR, IBAR, JBAR + 1, KBAR, IBAR + 1, 0, &
    KBAR + 1, IBAR + 1, JBAR + 1, KBAR + 1 /)
   IJKW_EDGE (11, 1:12) = (/ 0, JBAR, KBAR, IBAR + 1, JBAR, KBAR, 0, JBAR + 1, &
    KBAR + 1, IBAR + 1, JBAR + 1, KBAR + 1 /)
   IJKW_EDGE (12, 1:12) = (/ 1, 0, KBAR, 1, JBAR + 1, KBAR, 0, 0, KBAR + 1, 0, &
    JBAR + 1, KBAR + 1 /)
ELSE ! exchange excludes ghost cells
   IJKW_EDGE (1, 1:12) = (/ 1, 1, 1, IBAR, 1, 1, 1, 0, 0, IBAR, 0, 0 /)
   IJKW_EDGE (2, 1:12) = (/ IBAR, 1, 1, IBAR, JBAR, 1, IBAR + 1, 1, 0, IBAR + 1, &
    JBAR, 0 /)
   IJKW_EDGE (3, 1:12) = (/ 1, JBAR, 1, IBAR, JBAR, 1, 1, JBAR + 1, 0, IBAR, &
    JBAR + 1, 0 /)
   IJKW_EDGE (4, 1:12) = (/ 1, 1, 1, 1, JBAR, 1, 0, 1, 0, 0, JBAR, 0 /)
   IJKW_EDGE (5, 1:12) = (/ 1, 1, 1, 1, 1, KBAR, 0, 0, 1, 0, 0, KBAR /)
   IJKW_EDGE (6, 1:12) = (/ IBAR, 1, 1, IBAR, 1, KBAR, IBAR + 1, 0, 1, IBAR + 1, &
    0, KBAR /)
   IJKW_EDGE (7, 1:12) = (/ IBAR, JBAR, 1, IBAR, JBAR, KBAR, IBAR + 1, JBAR + 1, &
    1, IBAR + 1, JBAR + 1, KBAR /)
   IJKW_EDGE (8, 1:12) = (/ 1, JBAR, 1, 1, JBAR, KBAR, 0, JBAR + 1, 1, 0, JBAR + &
    1, KBAR /)
   IJKW_EDGE (9, 1:12) = (/ 1, 1, KBAR, IBAR, 1, KBAR, 1, 0, KBAR + 1, IBAR, 0, &
    KBAR + 1 /)
   IJKW_EDGE (10, 1:12) = (/ IBAR, 1, KBAR, IBAR, JBAR, KBAR, IBAR + 1, 1, KBAR &
    + 1, IBAR + 1, JBAR, KBAR + 1 /)
   IJKW_EDGE (11, 1:12) = (/ 1, JBAR, KBAR, IBAR, JBAR, KBAR, 1, JBAR + 1, KBAR &
    + 1, IBAR, JBAR + 1, KBAR + 1 /)
   IJKW_EDGE (12, 1:12) = (/ 1, 1, KBAR, 1, JBAR, KBAR, 0, 1, KBAR + 1, 0, JBAR, &
    KBAR + 1 /)
ENDIF
 
 
EDGE_NBR_LOOP: DO IEDGE = 1, NEDGE
 
   NOM = NBR_EDGE (NM, IEDGE)!  diagonal neighbor at edge IEDGE
   IF (NOM == 0) CYCLE EDGE_NBR_LOOP
 
   I1 = IJKW_EDGE (IEDGE, 1)
   I2 = IJKW_EDGE (IEDGE, 4)
   J1 = IJKW_EDGE (IEDGE, 2)
   J2 = IJKW_EDGE (IEDGE, 5)
   K1 = IJKW_EDGE (IEDGE, 3)
   K2 = IJKW_EDGE (IEDGE, 6)
 
   NIC_EDGE (NM, NOM) = (I2-I1+1) * (J2-J1+1) * (K2-K1+1)
   NIC_EDGE (NOM, NM) = NIC_EDGE (NM, NOM)
 
   ALLOCATE (S%OSCARC(NOM)%Z_EDGE(I1:I2, J1:J2, K1:K2))
 
   IF (SCARC_DEBUG .GE. 2) THEN
      WRITE (9,*) 'NIC_EDGE(', NOM, ',', NM, ')=', NIC_EDGE (NM, NOM)
      WRITE (9,*) 'IEDGE=', IEDGE, ': NOM=', NOM
      WRITE (9,*) 'IJKW_EDGE:'
      WRITE (9, '(12i4)') (IJKW_EDGE(IEDGE, III), III=1, 12)
      WRITE (9,*) 'I1=', I1
      WRITE (9,*) 'I2=', I2
      WRITE (9,*) 'J1=', J1
      WRITE (9,*) 'J2=', J2
      WRITE (9,*) 'K1=', K1
      WRITE (9,*) 'K2=', K2
      CALL flush (9)
   ENDIF
 
ENDDO EDGE_NBR_LOOP
 
 
 
!!!
!!! Initialize structures for diagonal-vertex neigbours communication
!!!
IJKW_VRTX (1, 1:12) = (/ 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0 /)
IJKW_VRTX (2, 1:12) = (/ IBAR, 1, 1, IBAR, 1, 1, IBAR + 1, 0, 0, IBAR + 1, 0, 0 &
 /)
IJKW_VRTX (3, 1:12) = (/ IBAR, JBAR, 1, IBAR, JBAR, 1, IBAR + 1, JBAR + 1, 0, &
 IBAR + 1, JBAR + 1, 0 /)
IJKW_VRTX (4, 1:12) = (/ 1, JBAR, 1, 1, JBAR, 1, 0, JBAR + 1, 0, 0, JBAR + 1, 0 &
 /)
IJKW_VRTX (5, 1:12) = (/ 1, 1, KBAR, 1, 1, KBAR, 0, 0, KBAR + 1, 0, 0, KBAR + 1 &
 /)
IJKW_VRTX (6, 1:12) = (/ IBAR, 1, KBAR, IBAR, 1, KBAR, IBAR + 1, 0, KBAR + 1, &
 IBAR + 1, 0, KBAR + 1 /)
IJKW_VRTX (7, 1:12) = (/ IBAR, JBAR, KBAR, IBAR, JBAR, KBAR, IBAR + 1, JBAR + 1, &
 KBAR + 1, IBAR + 1, JBAR + 1, KBAR + 1 /)
IJKW_VRTX (8, 1:12) = (/ 1, JBAR, KBAR, 1, JBAR, KBAR, 0, JBAR + 1, KBAR + 1, 0, &
 JBAR + 1, KBAR + 1 /)
 
 
VRTX_NBR_LOOP: DO IVRTX = 1, NVRTX
 
   NOM = NBR_VRTX (NM, IVRTX)!  diagonal neighbor at vertex IVRTX
   IF (NOM == 0) CYCLE VRTX_NBR_LOOP
 
   NIC_VRTX (NM, NOM) = 1
   NIC_VRTX (NOM, NM) = 1
 
   IF (SCARC_DEBUG .GE. 2) write (9,*) 'IVRTX=', IVRTX, ': NIC_VRTX(', NOM, ',', &
    NM, ')=', NIC_VRTX (NM, NOM)
 
   ALLOCATE (S%OSCARC(NOM)%Z_VRTX(1, 1, 1))
 
ENDDO VRTX_NBR_LOOP
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'IJKW_FACE:'
   DO IFACE = 1, NFACE
      WRITE (9, '(12i4)') (IJKW_FACE(IFACE, III), III=1, 12)
   ENDDO
   WRITE (9,*) 'IJKW_EDGE:'
   DO IEDGE = 1, NEDGE
      WRITE (9, '(12i4)') (IJKW_EDGE(IEDGE, III), III=1, 12)
   ENDDO
   WRITE (9,*) 'IJKW_VRTX:'
   DO IVRTX = 1, NVRTX
      WRITE (9, '(12i4)') (IJKW_VRTX(IVRTX, III), III=1, 12)
   ENDDO
   CALL flush (9)
ENDIF
 
CALL SCARC_RECEIVE (0, 0)
CALL SCARC_EXCHANGE (0, 0)
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'Leaving SCARC_INITIALIZE_MESH_EXCHANGE'
   CALL flush (9)
ENDIF
 
 
END SUBROUTINE SCARC_INITIALIZE_MESH_EXCHANGE3D
 
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Allocate full matrix for the usual 5-point-stencil (2D)
!!! 'full' means that it corresponds to the 1-process matrix of the serial problem
!!!
!!! Note: The matrix is stored "band wise", namely diagonal for diagonal.
!!! Although the sub- and superdiagonals are shorter than the main diagonal,
!!! all bands are stored in the same length and are simply filled with zeros
!!! at redundant positions. The 'wasting' of this storage space is justified
!!! by the possibility to use a much more efficient matrix-vector-multiplication
!!! (which doesn't have to use an expensive referencing logic)
!!!
!!! 1. band:  ibar-th subdiagonal for the southern '1'-entry
!!! 2. band:  1-th subdiagonal for the western  '1'-entry
!!! 3. band:  main diagonal with the '-4'-entry
!!! 4. band:  1-th superdiagonal for the eastern  '1'-entry
!!! 5. band:  ibar-th superdiagonal for the northern '1'-entry
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE_MATRICES2D (NM)
 
INTEGER :: NM, NOM
INTEGER :: I, K, IC, IC0, IDIAG, IW, IM, IZERO, IERR, III
TYPE (MESH_TYPE), POINTER :: M, M2
 
INCLUDE 'mpif.h'
 
IERR = 0
 
!!!
!!! Allocation of SCARC and OSCARC structures
!!!
ALLOCATE (SCARC(NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES2D', 'SCARC', IZERO)
 
S => SCARC (NM)
 
ALLOCATE (S%OSCARC(NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES2D', 'OSCARC', IZERO)
 
!!! CG-parameters for global and local solves
S%BFIRST_CG = .TRUE.
S%BFIRST = .TRUE.
S%BPRECON = .TRUE.
S%BGSADI = .FALSE.
 
!!!
!!! Initialization of SCARC variables
!!!
S%NCELLS_LOCAL = IBAR * KBAR ! number of local and global grid cells
 
ALLOCATE (S%NCELLS_LOCAL_ALL(1:NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES2D', 'NCELLS_LOCAL_ALL', IZERO)
 
S%NCELLS_LOCAL_ALL = 0
S%NCELLS_LOCAL_ALL (NM) = S%NCELLS_LOCAL
 
ALLOCATE (S%NCELLS_GLOBAL_ALL(1:NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES2D', 'NCELLS_GLOBAL_ALL', IZERO)
 
CALL MPI_ALLREDUCE (S%NCELLS_LOCAL_ALL(1), S%NCELLS_GLOBAL_ALL(1), NMESHES, &
 MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
 
S%NCELLS_GLOBAL = 0
DO IM = 1, NMESHES
   S%NCELLS_GLOBAL = S%NCELLS_GLOBAL + S%NCELLS_GLOBAL_ALL(IM)
ENDDO
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) '===================== NM=', NM
   WRITE (9,*) 'NMESHES=', NMESHES
   WRITE (9,*) 'NCELLS_LOCAL=', S%NCELLS_LOCAL
   WRITE (9,*) 'NCELLS_LOCAL_ALL=', (S%NCELLS_LOCAL_ALL(I), I=1, NMESHES)
   WRITE (9,*) 'NCELLS_GLOBAL_ALL=', (S%NCELLS_GLOBAL_ALL(I), I=1, NMESHES)
   WRITE (9,*) 'NCELLS_GLOBAL=', S%NCELLS_GLOBAL
   CALL flush (9)
ENDIF
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only temporary: check if singles grids are equidistant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
S%DXMIN = 1000.0_EB
S%DXMAX = 0.0_EB
DO I = 0, IBAR + 1
   S%DXMIN = Min (DX(I), S%DXMIN)
   S%DXMAX = Max (DX(I), S%DXMAX)
!if (SCARC_DEBUG.ge.2) write(9,*) 'DXMIN=',S%DXMIN,': DXMAX=',S%DXMAX
ENDDO
 
S%DZMIN = 1000.0_EB
S%DZMAX = 0.0_EB
DO I = 0, KBAR + 1
   S%DZMIN = Min (DZ(I), S%DZMIN)
   S%DZMAX = Max (DZ(I), S%DZMAX)
!if (SCARC_DEBUG.ge.2) write(9,*) 'DZMIN=',S%DZMIN,': DZMAX=',S%DZMAX
ENDDO
 
 
!write(*,'(i3,a,e24.18,a,e24.18)') NM,':DXMAX=',S%DXMAX,': DXMIN=',S%DXMIN
!write(*,'(i3,a,e24.18,a,e24.18)') NM,':DZMAX=',S%DZMAX,': DZMIN=',S%DZMIN
IF (Abs(S%DXMAX-S%DXMIN) .GE. 1.0E-6_EB) THEN
   WRITE (*, '(i3,a,e24.18,a,e24.18,a,e24.18)') NM, ':S%DXMIN=', S%DXMIN, ': S%D&
  &XMAX=', S%DXMAX, ' DIFF=', Abs (S%DXMAX-S%DXMIN)
   WRITE (*,*) 'x-grid-size not equidistant!'
ENDIF
IF (Abs(S%DZMAX-S%DZMIN) .GE. 1.0E-6_EB) THEN
   WRITE (*, '(i3,a,e24.18,a,e24.18,a,e24.18)') NM, ':S%DZMIN=', S%DZMIN, ': S%D&
  &ZMAX=', S%DZMAX, ' DIFF=', Abs (S%DZMAX-S%DZMIN)
   WRITE (*,*) 'z-grid-size not equidistant!'
ENDIF
IF (Abs(S%DXMAX-S%DZMAX) .GE. 1.0E-6_EB) THEN
   WRITE (*, '(i3,a,e24.18,a,e24.18,a,e24.18)') NM, ':S%DXMAX=', S%DXMAX, ': S%D&
  &ZMAX=', S%DZMAX, ' DIFF=', Abs (S%DXMAX-S%DZMAX)
   WRITE (*,*) 'x- and z-grid-size not equidistant!'
ENDIF
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize matrix diagonals for global matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NDIAG = 5
 
!!! Allocate full matrix corresponding to the band-wise storage technique
ALLOCATE (S%AGLOB(1:S%NCELLS_LOCAL, 1:NDIAG), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES2D', 'S%AGLOB', IZERO)
S%AGLOB = 0.0_EB
 
!!! Allocate local matrix corresponding to the band-wise storage technique
ALLOCATE (S%ALOC(1:S%NCELLS_LOCAL, 1:NDIAG), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES2D', 'S%ALOC', IZERO)
 
!!! Allocate Matrix Period if periodic boundary conditions are set
ALLOCATE(S%PERIOD(1:NEWC),STAT=IZERO)
CALL ChkMemErr('SCARC_INIT2D','S%PERIOD',IZERO)
S%PERIOD=0.0_EB
 
 
 
DO K = 1, KBAR
   DO I = 1, IBAR
 
      IC = (K-1) * IBAR + I
 
      IF (K > 1) S%AGLOB(IC, 1) = 1.0_EB ! lower subdiagonal corresponding to z-direction
      IF (I > 1) S%AGLOB(IC, 2) = 1.0_EB ! lower subdiagonal corresponding to x-direction
 
      S%AGLOB (IC, 3) = - 4.0_EB
 
      IF (I < IBAR) S%AGLOB(IC, 4) = 1.0_EB ! upper subdiagonal corresponding to x-direction
      IF (K < KBAR) S%AGLOB(IC, 5) = 1.0_EB ! upper subdiagonal corresponding to z-direction
 
 
      IF (SCARC_DEBUG .GE. 4) THEN
         WRITE (9,*) 'D1: AG(', IC, ',1)=', S%AGLOB(IC, 1)
         WRITE (9,*) 'D2: AG(', IC, ',2)=', S%AGLOB(IC, 2)
         WRITE (9,*) 'D3: AG(', IC, ',3)=', S%AGLOB(IC, 3)
         WRITE (9,*) 'D4: AG(', IC, ',4)=', S%AGLOB(IC, 4)
         WRITE (9,*) 'D5: AG(', IC, ',5)=', S%AGLOB(IC, 5)
      ENDIF
 
   ENDDO
ENDDO
 
!!! local matrix AL is first set to global matrix AG
S%ALOC = S%AGLOB
 
!!! set subdiaogonal value (needed for matvec-communication)
S%ASUB = 1.0_EB
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCARC: end periodic boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
IF (SCARC_PERIODIC==1) THEN
 
if (SCARC_DEBUG.ge.2) write(9,*) 'IN SCARC_PERIODIC:'
!  IF (NMESHES==1) THEN
!    DO K=1,KBAR
!      IC=(K-1)*IBAR+1
!      S%AGLOB(IC,3)=S%AGLOB(IC,3)-1.0_EB
!      !write(9,*) 'D:AGLOB(',IC,',3)=',S%AGLOB(IC,3)
!    ENDDO
!    DO K=1,KBAR
!      IC=K*IBAR
!      S%AGLOB(IC,3)=S%AGLOB(IC,3)-1.0_EB
!      !write(9,*) 'D:AGLOB(',IC,',3)=',S%AGLOB(IC,3)
!    ENDDO
!    DO I=1,IBAR
!      IC=I
!      S%AGLOB(IC,3)=S%AGLOB(IC,3)-1.0_EB
!      !write(9,*) 'D:AGLOB(',IC,',3)=',S%AGLOB(IC,3)
!    ENDDO
!    DO I=1,IBAR
!      IC=(KBAR-1)*IBAR+I
!      S%AGLOB(IC,3)=S%AGLOB(IC,3)-1.0_EB
!      !write(9,*) 'D:AGLOB(',IC,',3)=',S%AGLOB(IC,3)
!    ENDDO
 
 
!  !------------------------------------------------------------
!  ! square MxM-topology
!  !------------------------------------------------------------
 
    IF (NMESHES==1) THEN
      NMESHES_SQRT=1
    ELSE IF (NMESHES==4) THEN
      NMESHES_SQRT=2
    ELSE IF (NMESHES==9) THEN
      NMESHES_SQRT=3
    ELSE IF (NMESHES==16) THEN
      NMESHES_SQRT=4
    ELSE IF (NMESHES==25) THEN
      NMESHES_SQRT=5
    ELSE IF (NMESHES==36) THEN
      NMESHES_SQRT=6
    ELSE IF (NMESHES==49) THEN
      NMESHES_SQRT=7
    ELSE IF (NMESHES==64) THEN
      NMESHES_SQRT=8
    ENDIF

if (SCARC_DEBUG.ge.2) then
   write(9,*) 'NMESHES_SQRT=',NMESHES_SQRT
   write(9,*) 'NM=',NM
endif

    IF (NMESHES==1.OR.NM.LE.NMESHES_SQRT) THEN              !xmin --> left edge Diric
      DO K=1,KBAR
        IC=(K-1)*IBAR+1
        S%AGLOB(IC,3)=S%AGLOB(IC,3)-1.0_EB
        IW=WALL_INDEX(CELL_INDEX(1,1,K),-1)
        S%PERIOD(IW)=100
        !write(9,*) 'D:AGLOB(',IC,',3)=',S%AGLOB(IC,3),': S%PERIOD(',IW,')=100'
      ENDDO
    ENDIF
    IF (NMESHES==1.OR.NM.GT.NMESHES-NMESHES_SQRT) THEN      !xmax --> right edge Diric
      !write(9,*) ' right edge, NM=',NM,' NMESHES -NMESHES_SQRT:',NMESHES-NMESHES_SQRT
      DO K=1,KBAR
        IC=K*IBAR
        S%AGLOB(IC,3)=S%AGLOB(IC,3)-1.0_EB
        IW=WALL_INDEX(CELL_INDEX(IBAR,1,K),1)
        S%PERIOD(IW)=100
        !write(9,*) 'D:AGLOB(',IC,',3)=',S%AGLOB(IC,3),': S%PERIOD(',IW,')=100'
      ENDDO
    ENDIF
    IF (NMESHES==1.OR.MOD(NM,NMESHES_SQRT)==1) THEN         !zmin --> lower edge Diric
      !write(9,*) ' lower edge'
      DO I=1,IBAR
        IC=I
        S%AGLOB(IC,3)=S%AGLOB(IC,3)-1.0_EB
        IW=WALL_INDEX(CELL_INDEX(I,1,1),-3)
       !write(9,*) 'CELLINDEX=',CELL_INDEX(I,1,1)
        !write(9,*) 'WALLINDEX=',WALL_INDEX(CELL_INDEX(I,1,1),-3)
        S%PERIOD(IW)=100
        !write(9,*) 'D:AGLOB(',IC,',3)=',S%AGLOB(IC,3),': S%PERIOD(',IW,')=100'
      ENDDO
    ENDIF
    IF (NMESHES==1.OR.MOD(NM,NMESHES_SQRT)==0) THEN         !zmax --> upper edge Diric
      !write(9,*) ' upper edge'
      DO I=1,IBAR
        IC=(KBAR-1)*IBAR+I
        S%AGLOB(IC,3)=S%AGLOB(IC,3)-1.0_EB
        IW=WALL_INDEX(CELL_INDEX(I,1,KBAR),3)
        S%PERIOD(IW)=100
        !write(9,*) 'D:AGLOB(',IC,',3)=',S%AGLOB(IC,3),': S%PERIOD(',IW,')=100'
      ENDDO
    ENDIF

 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! SCARC: end periodic boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
else
 
 
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! adapt main diagonal to boundary type
!!! at interior boundaries full matrix information is kept
!!! at exterior boundaries matrix entries are corrected corresponding
!!! to the number of adjacent cells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! IOR = 1
!!!
 
DO K = 1, KBAR
 
   IC = (K-1) * IBAR + 1
 
   IF (LBC == 3 .OR. LBC == 4 .OR. LBC == 6) THEN ! Neumann boundary
      S%AGLOB (IC, 3) = S%AGLOB(IC, 3) + 1.0_EB
      S%ALOC (IC, 3) = S%ALOC(IC, 3) + 1.0_EB
      IF (SCARC_DEBUG .GE. 4) write (9,*) 'N 1:AG(', IC, ',4)=', S%AGLOB(IC, 3)
   ELSE IF (LBC == 1 .OR. LBC == 2 .OR. LBC == 5) THEN ! Dirichlet boundary
      IF (NBR_FACE(NM, 4) == 0) S%AGLOB(IC, 3) = S%AGLOB(IC, 3) - 1.0_EB
      S%ALOC (IC, 3) = S%ALOC(IC, 3) - 1.0_EB
      IF (SCARC_DEBUG .GE. 4) write (9, '(a,i3,a,f12.6,a,i3,a,i3)') 'D 1:AG(', &
       IC, ',4)=', S%AGLOB(IC, 3)
   ENDIF
 
ENDDO
 
!!!
!!! IOR = -1
!!!
DO K = 1, KBAR
 
   IC = K * IBAR
 
   IF (LBC == 2 .OR. LBC == 3 .OR. LBC == 6) THEN ! Neumann boundary
      S%AGLOB (IC, 3) = S%AGLOB(IC, 3) + 1.0_EB
      S%ALOC (IC, 3) = S%ALOC(IC, 3) + 1.0_EB
      IF (SCARC_DEBUG .GE. 4) write (9,*) 'N-1:AG(', IC, ',4)=', S%AGLOB(IC, 3)
   ELSE IF (LBC == 1 .OR. LBC == 4 .OR. LBC == 5) THEN ! Dirichlet boundary
      IF (NBR_FACE(NM, 2) == 0) S%AGLOB(IC, 3) = S%AGLOB(IC, 3) - 1.0_EB
      S%ALOC (IC, 3) = S%ALOC(IC, 3) - 1.0_EB
      IF (SCARC_DEBUG .GE. 4) write (9, '(a,i3,a,f12.6,a,i3,a,i3)') 'D-1:AG(', &
       IC, ',4)=', S%AGLOB(IC, 3)
   ENDIF
 
ENDDO
 
 
!!!
!!! IOR = 3
!!!
DO I = 1, IBAR
 
   IC = I
 
   IF (NBC == 3 .OR. NBC == 4) THEN ! Neumann boundary
      S%AGLOB (IC, 3) = S%AGLOB(IC, 3) + 1.0_EB
      S%ALOC (IC, 3) = S%ALOC(IC, 3) + 1.0_EB
      IF (SCARC_DEBUG .GE. 4) write (9,*) 'N 3:AG(', IC, ',4)=', S%AGLOB(IC, 3)
   ELSE IF (NBC == 1 .OR. NBC == 2) THEN ! Dirichlet boundary
      IF (NBR_FACE(NM, 1) == 0) S%AGLOB(IC, 3) = S%AGLOB(IC, 3) - 1.0_EB
      S%ALOC (IC, 3) = S%ALOC(IC, 3) - 1.0_EB
      IF (SCARC_DEBUG .GE. 4) write (9, '(a,i3,a,f12.6,a,i3,a,i3)') 'D 3:AG(', &
       IC, ',4)=', S%AGLOB(IC, 3)
   ENDIF
 
ENDDO
 
 
!!!
!!! IOR = -3
!!!
DO I = 1, IBAR
 
   IC = (KBAR-1) * IBAR + I
 
   IF (NBC == 3 .OR. NBC == 2) THEN ! Neumann boundary
      S%AGLOB (IC, 3) = S%AGLOB(IC, 3) + 1.0_EB
      S%ALOC (IC, 3) = S%ALOC(IC, 3) + 1.0_EB
      IF (SCARC_DEBUG .GE. 4) write (9,*) 'N-3:AG(', IC, ',4)=', S%AGLOB(IC, 3)
   ELSE IF (NBC == 1 .OR. NBC == 4) THEN ! Dirichlet boundary
      IF (NBR_FACE(NM, 3) == 0) S%AGLOB(IC, 3) = S%AGLOB(IC, 3) - 1.0_EB
      S%ALOC (IC, 3) = S%ALOC(IC, 3) - 1.0_EB
      IF (SCARC_DEBUG .GE. 4) write (9, '(a,i3,a,f12.6,a,i3,a,i3)') 'D-3:AG(', &
       IC, ',4)=', S%AGLOB(IC, 3)
   ENDIF
 
ENDDO
 
ENDIF
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'LBC=', LBC
   WRITE (9,*) 'MBC=', MBC
   WRITE (9,*) 'NBC=', NBC
   WRITE (9,*)
   WRITE (9,*) 'ASUB=', S%ASUB
   WRITE (9,*)
   WRITE (9,*) 'GLOBAL MATRIX AG:'
   WRITE (9,*) '1. Diag | 2. Diag | 3. Diag | 4. Diag | 5. Diag '
   DO IC = 1, S%NCELLS_LOCAL
      WRITE (9, '(5f12.4)') (S%AGLOB(IC, I), I=1, 5)
      IF (Mod(IC, IBAR) == 0) write (9,*) '-------------------------------------&
     &--------------------'
   ENDDO
   WRITE (9,*)
   WRITE (9,*) 'LOCAL MATRIX AL:'
   WRITE (9,*) '1. Diag | 2. Diag | 3. Diag | 4. Diag | 5. Diag '
   DO IC = 1, S%NCELLS_LOCAL
      WRITE (9, '(5f12.4)') (S%ALOC(IC, I), I=1, 5)
      IF (Mod(IC, IBAR) == 0) write (9,*) '-------------------------------------&
     &--------------------'
   ENDDO
   WRITE (9,*) 'DXMAX=', S%DXMAX
   WRITE (9,*) 'DXMAX^2=', S%DXMAX ** 2
ENDIF
 
 
!!!
!!! scale matrix with 1/h**2
!!! at the moment, where all the grids are equidistant, h = DXMAX
!!!
S%AGLOB = S%AGLOB / (S%DXMAX**2)
S%ALOC = S%ALOC / (S%DXMAX**2)
 
!!! ONLY TEMPORARY - works only for equidistant meshes !!!!!!!!!!
S%ASUB = 1.0_EB / (S%DXMAX**2)
!if (S%AGLOB(1,4)/=0) S%ASUB=S%AGLOB(1,4)
 
 
IF (SCARC_DEBUG .GE. 2) write (9,*) 'S%ASUB=', S%ASUB
 
 
END SUBROUTINE SCARC_INITIALIZE_MATRICES2D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  3D-initialization of SCARC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE_MATRICES3D (NM)
 
INTEGER :: NM, NOM, idim
INTEGER :: I, J, K, IC, IDIAG, IW, IM, IZERO, IERR, III, IPTR
TYPE (MESH_TYPE), POINTER :: M
 
INCLUDE 'mpif.h'
 
IERR = 0
 
ALLOCATE (SCARC(NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES3D', 'SCARC', IZERO)
 
S => SCARC (NM)
 
IF (NMESHES == 1) THEN
 
   WRITE (*,*) '================================================================&
  &====='
   WRITE (*,*) '=== SCARC-INFO: single-mesh topology'
   WRITE (*,*) '================================================================&
  &====='
 
ELSE
 
   IF (NM == 1) THEN
      WRITE (*,*) '=============================================================&
     &========'
      WRITE (*,*) '=== SCARC-INFO: multi-mesh topology'
      WRITE (*,*) '=== SCARC-INFO: NMESHES =', NMESHES
      WRITE (*,*) '=============================================================&
     &========'
   ENDIF
 
ENDIF
 
!!!
!!! Allocation of SCARC structures
!!!
 
!!! for all neighboring 'other meshes'
ALLOCATE (S%OSCARC(NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES3D', 'OSCARC', IZERO)
 
 
!!!
!!! Initialize some SCARC variables
!!!
 
!!! CG-parameters for global and local solves
S%BFIRST_CG = .TRUE.
S%BFIRST = .TRUE.
S%BPRECON = .TRUE.
S%BGSADI = .FALSE.
 
!!! number of local and global grid cells
S%NCELLS_LOCAL = IBAR * JBAR * KBAR
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'STARTING SCARC_INITIALIZE_MATRICES3D'
   WRITE (9,*) 'NCELLS_LOCAL=', S%NCELLS_LOCAL
   WRITE (9,*) 'NMESHES=', NMESHES
   WRITE (9,*) 'IBAR=', IBAR
   WRITE (9,*) 'JBAR=', JBAR
   WRITE (9,*) 'KBAR=', KBAR
   CALL flush (9)
ENDIF
 
ALLOCATE (S%NCELLS_LOCAL_ALL(1:NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES3D', 'NCELLS_LOCAL_ALL', IZERO)
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'After allocate, Izero=', IZERO
   CALL flush (9)
ENDIF
 
S%NCELLS_LOCAL_ALL = 0
S%NCELLS_LOCAL_ALL (NM) = S%NCELLS_LOCAL
 
ALLOCATE (S%NCELLS_GLOBAL_ALL(1:NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES3D', 'NCELLS_GLOBAL_ALL', IZERO)
 
CALL MPI_ALLREDUCE (S%NCELLS_LOCAL_ALL(1), S%NCELLS_GLOBAL_ALL(1), NMESHES, &
 MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
 
S%NCELLS_GLOBAL = 0
DO IM = 1, NMESHES
   S%NCELLS_GLOBAL = S%NCELLS_GLOBAL + S%NCELLS_GLOBAL_ALL(IM)
ENDDO
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'BFIRST=', S%BFIRST
   WRITE (9,*) 'NMESHES=', NMESHES
   WRITE (9,*) 'NCELLS_LOCAL=', S%NCELLS_LOCAL
   WRITE (9,*) 'NCELLS_LOCAL_ALL=', (S%NCELLS_LOCAL_ALL(I), I=1, NMESHES)
   WRITE (9,*) 'NCELLS_GLOBAL_ALL=', (S%NCELLS_GLOBAL_ALL(I), I=1, NMESHES)
   WRITE (9,*) 'NCELLS_GLOBAL=', S%NCELLS_GLOBAL
   WRITE (9,*) 'NMESHES=', NMESHES
   WRITE (9,*) 'NM=', NM
ENDIF
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! only temporary: check if singles grids are equidistant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
S%DXMIN = 1000.0_EB
S%DXMAX = 0.0_EB
DO I = 0, IBAR + 1
   S%DXMIN = Min (DX(I), S%DXMIN)
   S%DXMAX = Max (DX(I), S%DXMAX)
   IF (SCARC_DEBUG .GE. 2) write (9,*) 'DXMIN=', S%DXMIN, ': DXMAX=', S%DXMAX
ENDDO
 
S%DYMIN = 1000.0_EB
S%DYMAX = 0.0_EB
DO I = 0, JBAR + 1
   S%DYMIN = Min (DY(I), S%DYMIN)
   S%DYMAX = Max (DY(I), S%DYMAX)
   IF (SCARC_DEBUG .GE. 2) write (9,*) 'DYMIN=', S%DYMIN, ': DYMAX=', S%DYMAX
ENDDO
 
S%DZMIN = 1000.0_EB
S%DZMAX = 0.0_EB
DO I = 0, KBAR + 1
   S%DZMIN = Min (DZ(I), S%DZMIN)
   S%DZMAX = Max (DZ(I), S%DZMAX)
   IF (SCARC_DEBUG .GE. 2) write (9,*) 'DZMIN=', S%DZMIN, ': DZMAX=', S%DZMAX
ENDDO
 
 
!write(*,'(i3,a,e24.18,a,e24.18)') NM,':DXMAX=',S%DXMAX,': DXMIN=',S%DXMIN
!write(*,'(i3,a,e24.18,a,e24.18)') NM,':DZMAX=',S%DZMAX,': DZMIN=',S%DZMIN
IF (Abs(S%DXMAX-S%DXMIN) .GE. 1.0E-6_EB) THEN
   WRITE (*, '(i3,a,e24.18,a,e24.18,a,e24.18)') NM, ':S%DXMIN=', S%DXMIN, ': S%D&
  &XMAX=', S%DXMAX, ' DIFF=', Abs (S%DXMAX-S%DXMIN)
   WRITE (*,*) 'x-grid-size not equidistant!'
ENDIF
IF (Abs(S%DYMAX-S%DYMIN) .GE. 1.0E-6_EB) THEN
   WRITE (*, '(i3,a,e24.18,a,e24.18,a,e24.18)') NM, ':S%DXMIN=', S%DYMIN, ': S%D&
  &YMAX=', S%DYMAX, ' DIFF=', Abs (S%DYMAX-S%DYMIN)
   WRITE (*,*) 'y-grid-size not equidistant!'
ENDIF
IF (Abs(S%DZMAX-S%DZMIN) .GE. 1.0E-6_EB) THEN
   WRITE (*, '(i3,a,e24.18,a,e24.18,a,e24.18)') NM, ':S%DZMIN=', S%DZMIN, ': S%D&
  &ZMAX=', S%DZMAX, ' DIFF=', Abs (S%DZMAX-S%DZMIN)
   WRITE (*,*) 'z-grid-size not equidistant!'
ENDIF
IF (Abs(S%DXMAX-S%DZMAX) .GE. 1.0E-6_EB) THEN
   WRITE (*, '(i3,a,e24.18,a,e24.18,a,e24.18)') NM, ':S%DXMAX=', S%DXMAX, ': S%D&
  &ZMAX=', S%DZMAX, ' DIFF=', Abs (S%DXMAX-S%DZMAX)
   WRITE (*,*) 'x- and z-grid-size not equidistant!'
ENDIF
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize matrix diagonals for global matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NDIAG = 7
 
!!! Allocate full matrix corresponding to the band-wise storage technique
ALLOCATE (S%AGLOB(1:S%NCELLS_LOCAL, 1:NDIAG), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES3D', 'S%AGLOB', IZERO)
S%AGLOB = 0.0_EB
 
ALLOCATE (S%ALOC(1:S%NCELLS_LOCAL, 1:NDIAG), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES3D', 'S%ALOC', IZERO)
S%ALOC = 0.0_EB
 
IF (SCARC_DEBUG .GE. 2) write (9,*) 'NEWC=', NEWC
ALLOCATE (S%PERIOD(1:NEWC), STAT=IZERO)
CALL CHKMEMERR ('SCARC_INITIALIZE_MATRICES3D', 'S%PERIOD', IZERO)
S%PERIOD = 0.0_EB
 
DO K = 1, KBAR
   DO J = 1, JBAR
      DO I = 1, IBAR
 
         IC = (K-1) * IBAR * JBAR + (J-1) * IBAR + I
 
         IF (K > 1) S%AGLOB(IC, 1) = 1.0_EB ! lower subdiagonal corresponding to kbar
         IF (J > 1) S%AGLOB(IC, 2) = 1.0_EB ! lower subdiagonal corresponding to jbar
         IF (I > 1) S%AGLOB(IC, 3) = 1.0_EB ! lower subdiagonal corresponding to ibar
 
         S%AGLOB (IC, 4) = - 6.0_EB
 
         IF (I < IBAR) S%AGLOB(IC, 5) = 1.0_EB ! upper subdiagonal corresponding to ibar
         IF (J < JBAR) S%AGLOB(IC, 6) = 1.0_EB ! upper subdiagonal corresponding to jbar
         IF (K < KBAR) S%AGLOB(IC, 7) = 1.0_EB ! upper subdiagonal corresponding to kbar
 
 
         IF (SCARC_DEBUG .GE. 2) THEN
            WRITE (9,*) 'D1: AGLOB(', IC, ',1)=', S%AGLOB(IC, 1)
            WRITE (9,*) 'D2: AGLOB(', IC, ',2)=', S%AGLOB(IC, 2)
            WRITE (9,*) 'D3: AGLOB(', IC, ',3)=', S%AGLOB(IC, 3)
            WRITE (9,*) 'D4: AGLOB(', IC, ',4)=', S%AGLOB(IC, 4)
            WRITE (9,*) 'D5: AGLOB(', IC, ',5)=', S%AGLOB(IC, 5)
            WRITE (9,*) 'D6: AGLOB(', IC, ',6)=', S%AGLOB(IC, 6)
            WRITE (9,*) 'D7: AGLOB(', IC, ',7)=', S%AGLOB(IC, 7)
         ENDIF
 
      ENDDO
   ENDDO
ENDDO
 
S%ALOC = S%AGLOB
 
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'NBR_FACE: '
   WRITE (9, '(6i4)') (NBR_FACE(NM, III), III=1, 6)
   WRITE (9,*) 'BEFORE BOUNDARY SETTING'
   WRITE (9,*) 'LBC=', LBC
   WRITE (9,*) 'MBC=', MBC
   WRITE (9,*) 'NBC=', NBC
   WRITE (9,*)
   WRITE (9,*) 'ASUB=', S%ASUB
   WRITE (9,*)
   WRITE (9,*) 'GLOBAL MATRIX AGLOB:'
   WRITE (9,*) '1. Diag | 2. Diag | 3. Diag | 4. Diag | 5. Diag | 6. Diag | 7. D&
  &iag '
   DO IC = 1, S%NCELLS_LOCAL
      WRITE (9, '(7f10.4)') (S%AGLOB(IC, I), I=1, 7)
      IF (Mod(IC, IBAR) == 0) write (9,*) '-------------------------------------&
     &--------------------'
   ENDDO
   WRITE (9,*) 'LOCAL MATRIX AGLOB:'
   WRITE (9,*) '1. Diag | 2. Diag | 3. Diag | 4. Diag | 5. Diag | 6. Diag | 7. D&
  &iag '
   DO IC = 1, S%NCELLS_LOCAL
      WRITE (9, '(7f10.4)') (S%ALOC(IC, I), I=1, 7)
      IF (Mod(IC, IBAR) == 0) write (9,*) '-------------------------------------&
     &--------------------'
   ENDDO
ENDIF
 
 
 
!!!
!!! adapt main diagonal to boundary type
!!! at interior boundaries full matrix information is kept
!!! at exterior boundaries matrix entries are corrected corresponding
!!! to the number of adjacent cells
!!!
 
!!!
!!! IOR = 1
!!!
DO K = 1, KBAR
   DO J = 1, JBAR
 
      IC = (K-1) * IBAR * JBAR + (J-1) * IBAR + 1
 
      IF (LBC == 3 .OR. LBC == 4 .OR. LBC == 6) THEN ! Neumann boundary
         S%AGLOB (IC, 4) = S%AGLOB(IC, 4) + 1.0_EB
         S%ALOC (IC, 4) = S%ALOC(IC, 4) + 1.0_EB
         IF (SCARC_DEBUG .GE. 2) write (9,*) 'N 1:AGLOB(', IC, ',4)=', &
          S%AGLOB(IC, 4)
      ELSE IF (LBC == 1 .OR. LBC == 2 .OR. LBC == 5) THEN ! Dirichlet boundary
!IW=WALL_INDEX(CELL_INDEX(1,J,K),-1)
!IF (BOUNDARY_TYPE(IW) /= INTERPOLATED_BOUNDARY) S%AGLOB(IC,4)=S%AGLOB(IC,4)-1.0_EB
         IF (NBR_FACE(NM, 5) == 0) S%AGLOB(IC, 4) = S%AGLOB(IC, 4) - 1.0_EB
         S%ALOC (IC, 4) = S%ALOC(IC, 4) - 1.0_EB
         IF (SCARC_DEBUG .GE. 2) write (9, '(a,i3,a,f12.6,a,i3,a,i3)') 'D 1:AGLO&
        &B(', IC, ',4)=', S%AGLOB(IC, 4), ' B-TYPE(', IW, ')=', BOUNDARY_TYPE &
          (IW)
      ENDIF
 
   ENDDO
   IF (SCARC_DEBUG .GE. 2) write (9,*) '------------------------------------'
ENDDO
 
 
!!!
!!! IOR = -1
!!!
DO K = 1, KBAR
   DO J = 1, JBAR
 
      IC = (K-1) * IBAR * JBAR + J * IBAR
 
      IF (LBC == 2 .OR. LBC == 3 .OR. LBC == 6) THEN ! Neumann boundary
         S%AGLOB (IC, 4) = S%AGLOB(IC, 4) + 1.0_EB
         S%ALOC (IC, 4) = S%ALOC(IC, 4) + 1.0_EB
         IF (SCARC_DEBUG .GE. 2) write (9,*) 'N-1:AGLOB(', IC, ',4)=', &
          S%AGLOB(IC, 4)
      ELSE IF (LBC == 1 .OR. LBC == 4 .OR. LBC == 5) THEN ! Dirichlet boundary
!IW=WALL_INDEX(CELL_INDEX(IBAR,J,K),1)
!IF (BOUNDARY_TYPE(IW) /= INTERPOLATED_BOUNDARY) S%AGLOB(IC,4)=S%AGLOB(IC,4)-1.0_EB
         IF (NBR_FACE(NM, 3) == 0) S%AGLOB(IC, 4) = S%AGLOB(IC, 4) - 1.0_EB
         S%ALOC (IC, 4) = S%ALOC(IC, 4) - 1.0_EB
         IF (SCARC_DEBUG .GE. 2) write (9, '(a,i3,a,f12.6,a,i3,a,i3)') 'D-1:AGLO&
        &B(', IC, ',4)=', S%AGLOB(IC, 4), ' B-TYPE(', IW, ')=', BOUNDARY_TYPE &
          (IW)
      ENDIF
 
   ENDDO
   IF (SCARC_DEBUG .GE. 2) write (9,*) '------------------------------------'
ENDDO
 
 
 
!!!
!!! IOR = 2
!!!
DO K = 1, KBAR
   DO I = 1, IBAR
 
      IC = (K-1) * IBAR * JBAR + I
 
      IF (MBC == 3 .OR. MBC == 4) THEN ! Neumann boundary
         S%AGLOB (IC, 4) = S%AGLOB(IC, 4) + 1.0_EB
         S%ALOC (IC, 4) = S%ALOC(IC, 4) + 1.0_EB
         IF (SCARC_DEBUG .GE. 2) write (9,*) 'N 2:AGLOB(', IC, ',4)=', &
          S%AGLOB(IC, 4)
      ELSE IF (MBC == 1 .OR. MBC == 2) THEN ! Dirichlet boundary
!IW=WALL_INDEX(CELL_INDEX(I,1,K),-2)
!IF (BOUNDARY_TYPE(IW) /= INTERPOLATED_BOUNDARY) S%AGLOB(IC,4)=S%AGLOB(IC,4)-1.0_EB
         IF (NBR_FACE(NM, 2) == 0) S%AGLOB(IC, 4) = S%AGLOB(IC, 4) - 1.0_EB
         S%ALOC (IC, 4) = S%ALOC(IC, 4) - 1.0_EB
         IF (SCARC_DEBUG .GE. 2) write (9, '(a,i3,a,f12.6,a,i3,a,i3)') 'D 2:AGLO&
        &B(', IC, ',4)=', S%AGLOB(IC, 4), ' B-TYPE(', IW, ')=', BOUNDARY_TYPE &
          (IW)
      ENDIF
 
   ENDDO
   IF (SCARC_DEBUG .GE. 2) write (9,*) '------------------------------------'
ENDDO
 
DO IC = 1, S%NCELLS_LOCAL
   IF (SCARC_DEBUG .GE. 2) write (9, '(7f10.4)') (S%AGLOB(IC, I), I=1, 7)
!if (mod(IC,IBAR)==0) write(9,*) '333---------------------------------------------------------'
ENDDO
 
 
!!!
!!! IOR = -2
!!!
DO K = 1, KBAR
   DO I = 1, IBAR
 
      IC = (K-1) * IBAR * JBAR + (JBAR-1) * IBAR + I
 
      IF (MBC == 3 .OR. MBC == 2) THEN ! Neumann boundary
         S%AGLOB (IC, 4) = S%AGLOB(IC, 4) + 1.0_EB
         S%ALOC (IC, 4) = S%ALOC(IC, 4) + 1.0_EB
         IF (SCARC_DEBUG .GE. 2) write (9,*) 'N-2:AGLOB(', IC, ',4)=', &
          S%AGLOB(IC, 4)
      ELSE IF (MBC == 1 .OR. MBC == 4) THEN ! Dirichlet boundary
!IW=WALL_INDEX(CELL_INDEX(I,JBAR,K),2)
!IF (BOUNDARY_TYPE(IW) /= INTERPOLATED_BOUNDARY) S%AGLOB(IC,4)=S%AGLOB(IC,4)-1.0_EB
         IF (NBR_FACE(NM, 4) == 0) S%AGLOB(IC, 4) = S%AGLOB(IC, 4) - 1.0_EB
         S%ALOC (IC, 4) = S%ALOC(IC, 4) - 1.0_EB
         IF (SCARC_DEBUG .GE. 2) write (9, '(a,i3,a,f12.6,a,i3,a,i3)') 'D-2:AGLO&
        &B(', IC, ',4)=', S%AGLOB(IC, 4), ' B-TYPE(', IW, ')=', BOUNDARY_TYPE &
          (IW)
      ENDIF
 
   ENDDO
   IF (SCARC_DEBUG .GE. 2) write (9,*) '------------------------------------'
ENDDO
 
 
 
 
!!!
!!! IOR = 3
!!!
DO J = 1, JBAR
   DO I = 1, IBAR
 
      IC = (J-1) * IBAR + I
 
      IF (NBC == 3 .OR. NBC == 4) THEN ! Neumann boundary
         S%AGLOB (IC, 4) = S%AGLOB(IC, 4) + 1.0_EB
         S%ALOC (IC, 4) = S%ALOC(IC, 4) + 1.0_EB
         IF (SCARC_DEBUG .GE. 2) write (9,*) 'N 3:AGLOB(', IC, ',4)=', &
          S%AGLOB(IC, 4)
      ELSE IF (NBC == 1 .OR. NBC == 2) THEN ! Dirichlet boundary
!IW=WALL_INDEX(CELL_INDEX(I,J,1),-3)
!IF (BOUNDARY_TYPE(IW) /= INTERPOLATED_BOUNDARY) S%AGLOB(IC,4)=S%AGLOB(IC,4)-1.0_EB
         IF (NBR_FACE(NM, 1) == 0) S%AGLOB(IC, 4) = S%AGLOB(IC, 4) - 1.0_EB
         S%ALOC (IC, 4) = S%ALOC(IC, 4) - 1.0_EB
         IF (SCARC_DEBUG .GE. 2) write (9, '(a,i3,a,f12.6,a,i3,a,i3)') 'D 3:AGLO&
        &B(', IC, ',4)=', S%AGLOB(IC, 4), ' B-TYPE(', IW, ')=', BOUNDARY_TYPE &
          (IW)
      ENDIF
 
   ENDDO
   IF (SCARC_DEBUG .GE. 2) write (9,*) '------------------------------------'
ENDDO
 
 
 
!!!
!!! IOR = -3
!!!
DO J = 1, JBAR
   DO I = 1, IBAR
 
      IC = (KBAR-1) * IBAR * JBAR + (J-1) * IBAR + I
 
      IF (NBC == 3 .OR. NBC == 2) THEN ! Neumann boundary
         S%AGLOB (IC, 4) = S%AGLOB(IC, 4) + 1.0_EB
         S%ALOC (IC, 4) = S%ALOC(IC, 4) + 1.0_EB
         IF (SCARC_DEBUG .GE. 2) write (9,*) 'N-3:AGLOB(', IC, ',4)=', &
          S%AGLOB(IC, 4)
      ELSE IF (NBC == 1 .OR. NBC == 4) THEN ! Dirichlet boundary
!IW=WALL_INDEX(CELL_INDEX(I,J,KBAR),3)
!IF (BOUNDARY_TYPE(IW) /= INTERPOLATED_BOUNDARY) S%AGLOB(IC,4)=S%AGLOB(IC,4)-1.0_EB
         IF (NBR_FACE(NM, 6) == 0) S%AGLOB(IC, 4) = S%AGLOB(IC, 4) - 1.0_EB
         S%ALOC (IC, 4) = S%ALOC(IC, 4) - 1.0_EB
         IF (SCARC_DEBUG .GE. 2) write (9, '(a,i3,a,f12.6,a,i3,a,i3)') 'D-3:AGLO&
        &B(', IC, ',4)=', S%AGLOB(IC, 4), ' B-TYPE(', IW, ')=', BOUNDARY_TYPE &
          (IW)
      ENDIF
 
   ENDDO
   IF (SCARC_DEBUG .GE. 2) write (9,*) '------------------------------------'
ENDDO
 
 
 
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'LBC=', LBC
   WRITE (9,*) 'MBC=', MBC
   WRITE (9,*) 'NBC=', NBC
   WRITE (9,*)
   WRITE (9,*) 'ASUB=', S%ASUB
   WRITE (9,*)
   WRITE (9,*) 'AFTER: GLOBAL MATRIX AGLOB:'
   WRITE (9,*) '1. Diag | 2. Diag | 3. Diag | 4. Diag | 5. Diag | 6. Diag | 7. D&
  &iag '
ENDIF
 
DO IC = 1, S%NCELLS_LOCAL
   IF (SCARC_DEBUG .GE. 2) write (9, '(7f10.4)') (S%AGLOB(IC, I), I=1, 7)
!if (mod(IC,IBAR)==0) write(9,*) '---------------------------------------------------------'
ENDDO
 
 
 
IF (SCARC_DEBUG .GE. 2) write (9,*)
IF (SCARC_DEBUG .GE. 2) write (9,*) 'AFTER: LOCAL MATRIX ALOC:'
IF (SCARC_DEBUG .GE. 2) write (9,*) '1. Diag | 2. Diag | 3. Diag | 4. Diag | 5. &
  &Diag | 6. Diag | 7. Diag '
DO IC = 1, S%NCELLS_LOCAL
   IF (SCARC_DEBUG .GE. 2) write (9, '(7f10.4)') (S%ALOC(IC, I), I=1, 7)
!if (mod(IC,IBAR)==0) write(9,*) '---------------------------------------------------------'
ENDDO
 
IF (SCARC_DEBUG .GE. 2) write (9,*) 'DXMAX=', S%DXMAX
IF (SCARC_DEBUG .GE. 2) write (9,*) 'DXMAX^2=', S%DXMAX ** 2
 
!!!
!!! scale matrix with 1/h**2
!!! at the moment, where all the grids are equidistant, h = DXMAX
!!!
S%ASUB = 1.0_EB
 
!write(*,*) 'Achtung: Matrix ist NICHT skaliert !!!'
S%AGLOB = S%AGLOB / (S%DXMAX**2)
S%ALOC = S%ALOC / (S%DXMAX**2)
S%ASUB = S%ASUB / (S%DXMAX**2)
 
!!! ONLY TEMPORARY - works only for equidistant meshes !!!!!!!!!!
 
 
IF (SCARC_DEBUG .GE. 2) write (9,*) 'S%ASUB=', S%ASUB
 
 
END SUBROUTINE SCARC_INITIALIZE_MATRICES3D
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize global cg-method based on global possion-matrix S%AGLOB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE_CG2D (NM)
 
INTEGER :: NM, NOM, IZERO, I
 
 
IF (SCARC_DEBUG .GE. 2) write (9,*) '=========== STARTING SCARC_INITIALIZE_CG2D'
CALL flush (9)
 
S => SCARC (NM)
IZERO = 0
 
!!! Allocate auxiliary vectors for defect correction method
ALLOCATE (S%X(IBAR, 1, KBAR), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'X', IZERO)
S%X = 0.0_EB
 
ALLOCATE (S%F(IBAR, 1, KBAR), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'B', IZERO)
S%F = 0.0_EB
 
ALLOCATE (S%D(IBAR, 1, KBAR), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'D', IZERO)
S%D = 0.0_EB
 
ALLOCATE (S%R(IBAR, 1, KBAR), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'R', IZERO)
S%R = 0.0_EB
 
ALLOCATE (S%G(IBAR, 1, KBAR), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'G', IZERO)
S%G = 0.0_EB
 
ALLOCATE (S%Y(IBAR, 1, KBAR), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'Y', IZERO)
S%Y = 0.0_EB
 
ALLOCATE (S%Z(0:IBP1, 0:JBP1, 0:KBP1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'Z', IZERO)
S%Z = 0.0_EB
 
 
ALLOCATE (S%BXS0(1, KBP1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BXS0', IZERO)
S%BXS0 = 0.0_EB
 
ALLOCATE (S%BXF0(1, KBP1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BXF0', IZERO)
S%BXF0 = 0.0_EB
 
ALLOCATE (S%BZS0(IBP1, 1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BZS0', IZERO)
S%BZS0 = 0.0_EB
 
ALLOCATE (S%BZF0(IBP1, 1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BZF0', IZERO)
S%BZF0 = 0.0_EB
 
 
!!! Allocate vectors for global scalar products
ALLOCATE (S%SP_LOCAL_ALL(1:NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'SP_LOCAL_ALL', IZERO)
S%SP_LOCAL_ALL = 0.0_EB
 
ALLOCATE (S%SP_GLOBAL_ALL(1:NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'SP_GLOBAL_ALL', IZERO)
S%SP_GLOBAL_ALL = 0.0_EB
 
 
 
S%BFIRST = .FALSE.
IF (SCARC_DEBUG .GE. 2) write (9,*) 'INITIALIZE_CG2D: STOP'
CALL flush (9)
 
END SUBROUTINE SCARC_INITIALIZE_CG2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize global cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INITIALIZE_CG3D (NM)
 
INTEGER :: NM, NOM, IZERO
 
IF (SCARC_DEBUG .GE. 2) write (9,*) '=========== STARTING SCARC_INITIALIZE_CG3D'
CALL flush (9)
 
S => SCARC (NM)
IZERO = 0
 
!!! Allocate auxiliary vectors for defect correction method
ALLOCATE (S%X(0:IBAR+1, 0:JBAR+1, 0:KBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'X', IZERO)
S%X = 0.0_EB
 
ALLOCATE (S%F(0:IBAR+1, 0:JBAR+1, 0:KBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'B', IZERO)
S%F = 0.0_EB
 
ALLOCATE (S%D(0:IBAR+1, 0:JBAR+1, 0:KBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'D', IZERO)
S%D = 0.0_EB
 
ALLOCATE (S%R(0:IBAR+1, 0:JBAR+1, 0:KBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'R', IZERO)
S%R = 0.0_EB
 
ALLOCATE (S%G(0:IBAR+1, 0:JBAR+1, 0:KBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'G', IZERO)
S%G = 0.0_EB
 
ALLOCATE (S%Y(0:IBAR+1, 0:JBAR+1, 0:KBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'Y', IZERO)
S%Y = 0.0_EB
 
ALLOCATE (S%Z(0:IBAR+1, 0:JBAR+1, 0:KBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'Z', IZERO)
S%Z = 0.0_EB
 
ALLOCATE (S%TMP(0:IBAR+1, 0:JBAR+1, 0:KBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'TMP', IZERO)
S%TMP = 0.0_EB
 
ALLOCATE (S%BXS0(JBAR+1, KBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BXS0', IZERO)
S%BXS0 = 0.0_EB
 
ALLOCATE (S%BXF0(JBAR+1, KBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BXF0', IZERO)
S%BXF0 = 0.0_EB
 
ALLOCATE (S%BYS0(IBAR+1, KBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BXS0', IZERO)
S%BYS0 = 0.0_EB
 
ALLOCATE (S%BYF0(IBAR+1, KBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BXF0', IZERO)
S%BYF0 = 0.0_EB
 
ALLOCATE (S%BZS0(IBAR+1, JBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BZS0', IZERO)
S%BZS0 = 0.0_EB
 
ALLOCATE (S%BZF0(IBAR+1, JBAR+1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BZF0', IZERO)
S%BZF0 = 0.0_EB
 
 
!!! Allocate vectors for global scalar products
ALLOCATE (S%SP_LOCAL_ALL(1:NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'SP_LOCAL_ALL', IZERO)
S%SP_LOCAL_ALL = 0.0_EB
 
ALLOCATE (S%SP_GLOBAL_ALL(1:NMESHES), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'SP_GLOBAL_ALL', IZERO)
S%SP_GLOBAL_ALL = 0.0_EB
 
 
S%BFIRST = .FALSE.
IF (SCARC_DEBUG .GE. 2) write (9,*) 'INITIALIZE_CG3D: STOP'
CALL flush (9)
 
END SUBROUTINE SCARC_INITIALIZE_CG3D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global cg-method based on global possion-matrix
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE SCARC_CG2D(NM, ITRN, POIS_PTB, SAVE1, WORK, HX)
SUBROUTINE SCARC_CG2D (BXS, BXF, BZS, BZF, LDIMF, PRHS, NM, POIS_PTB, SAVE1, WORK, &
 HX)
 
INTEGER :: LDIMF, NM
INTEGER :: ITE, IREL, ICONV, ITYPE, I, J, K, IZERO, NOM, IERR, IPRECON, III, JJJ
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, I1, J1, II, KK
 
REAL (EB), DIMENSION (:) :: WORK
REAL (EB), DIMENSION (0:) :: HX
REAL (EB), DIMENSION (-3:) :: SAVE1
REAL (EB) :: POIS_PTB
REAL (EB) PRHS (LDIMF,*)
REAL (EB) BXS (*), BXF (*), BZS (*), BZF (*)
!REAL(EB), DIMENSION(:,:):: BXS, BXF, BYS, BYF, BZS, BZF
!REAL(EB), DIMENSION(:)  :: SAVE1, WORK, HX
 
TYPE (SCARC_TYPE), POINTER :: S
TYPE (MESH_TYPE),  POINTER :: M
 
INCLUDE 'mpif.h'
 
S => SCARC (NM)
M  => MESHES(NM)
 
 
N_REQ_MV = 0
IERR = 0
 
 
! initialize working vectors
S%G = 0.0_EB
S%Y = 0.0_EB
S%R = 0.0_EB
S%D = 0.0_EB
S%X = 0.0_EB
 

 
!do k=1,kbar
!  do i=1,ibar
!    S%X(I,1,K)=M%H(I,1,K)
!  enddo
!enddo
 
!S%F=0.0_EB
DO K = 1, KBAR
   DO I = 1, IBAR
!S%F(I,1,K)=M%PRHS(I,1,K)
      S%F (I, 1, K) = PRHS (I, K)
   ENDDO
ENDDO
 
 
IF (SCARC_DEBUG .GE. 3) THEN
   WRITE (9,*) 'Starting SCARC_CG2D'
   WRITE (9,*) 'IBAR=:', IBAR
   WRITE (9,*) 'JBAR=:', JBAR
   WRITE (9,*) 'KBAR=:', KBAR
   WRITE (9,*) 'SCARC_NIT=:', SCARC_NIT
   WRITE (9,*) 'SCARC_EPS=:', SCARC_EPS
   WRITE (9,*) 'SCARC_EPS_REL=:', SCARC_EPS_REL
   WRITE (9, '(a,10f12.6)') '1: S%SP_LOCAL_ALL=', (S%SP_LOCAL_ALL(I), I=1, &
    NMESHES)
   CALL flush (9)
ENDIF
CALL SCARC_COMPARE_SINGLE (S%X, 'SARC', 1, 'X1    ',0)
CALL SCARC_COMPARE_SINGLE (S%F, 'SARC', 1, 'F1    ',0)
 
 
!!!
!!! initialize some method parameters
!!!
ITYPE = NCOM_TYPE_GLOBAL
IREL = 0
ICONV = 0
IPRECON = 2
 
!S%NIT    =5000
!S%EPS    =1.0E-10_EB
!S%EPS_REL=1.0E-2_EB
 
ALLOCATE (S%BXS0(1, KBAR), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BXS0', IZERO)
S%BXS0 = 0.0_EB
 
ALLOCATE (S%BXF0(1, KBAR), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BXF0', IZERO)
S%BXF0 = 0.0_EB
 
ALLOCATE (S%BZS0(IBAR, 1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BZS0', IZERO)
S%BZS0 = 0.0_EB
 
ALLOCATE (S%BZF0(IBAR, 1), STAT=IZERO)
CALL CHKMEMERR ('SCARC', 'BZF0', IZERO)
S%BZF0 = 0.0_EB
 
 
!!!
!!! initialize boundary values
!!!
IF (NMESHES == 1) THEN
   DO KK = 1, KBAR
      S%BXS0 (1, KK) = BXS (KK)
      S%BXF0 (1, KK) = BXF (KK)
   ENDDO
   DO II = 1, IBAR
      S%BZS0 (II, 1) = BZS (II)
      S%BZF0 (II, 1) = BZF (II)
   ENDDO
ELSE
   S%BXS0 = 0.0_EB
   S%BXF0 = 0.0_EB
   S%BZS0 = 0.0_EB
   S%BZF0 = 0.0_EB
   IF (NBR_FACE(NM, 4) == 0) THEN ! IOR= 1: external boundary
      DO KK = 1, KBAR
         S%BXS0 (1, KK) = BXS (KK)
      ENDDO
   ENDIF
   IF (NBR_FACE(NM, 2) == 0) THEN ! IOR=-1: external boundary
      DO KK = 1, KBAR
         S%BXF0 (1, KK) = BXF (KK)
      ENDDO
   ENDIF
   IF (NBR_FACE(NM, 1) == 0) THEN ! IOR= 3: external boundary
      DO II = 1, IBAR
         S%BZS0 (II, 1) = BZS (II)
      ENDDO
   ENDIF
   IF (NBR_FACE(NM, 3) == 0) THEN ! IOR=-3: external boundary
      DO II = 1, IBAR
         S%BZF0 (II, 1) = BZF (II)
      ENDDO
   ENDIF
ENDIF
 
CALL SCARC_SETBDRY2D (S%F, S%BXS0, S%BXF0, S%BZS0, S%BZF0, SAVE1)
 
 
if (scarc_debug.ge.2) then
CALL SCARC_COMPARE_SINGLE (S%X, 'SARC', 1, 'X0    ',0)
CALL SCARC_COMPARE_SINGLE (S%F, 'SARC', 1, 'F0    ',0)
CALL SCARC_COMPARE_SINGLE (S%R, 'SARC', 1, 'R0    ',0)
WRITE (9,100) SHAPE(S%D)
100 FORMAT (1X,'The shape of the array is:          ', 7I6)

WRITE (9,110) SIZE(S%D)
110 FORMAT (1X,'The shape of the array is:          ', I10)

WRITE (9,120) LBOUND(S%D)
120 FORMAT (1X,'The lower bounds of the array are: ', 7I6)

WRITE (9,130) UBOUND(S%D)
130 FORMAT (1X,'The upper bounds of the array are: ', 7I6)
endif

!!!
!!! calculate initial defect Y = B - A*X and get l2-norm of it
!!!
CALL SCARC_MATVEC_MUL2D (S%AGLOB, S%X, S%R, 1.0_EB, 0.0_EB, NM, ITYPE, 3)
CALL SCARC_VECADD2D (S%F, S%R,-1.0_EB, 1.0_EB)
CALL SCARC_L2NORM2D (S%R, S%RES, NM, ITYPE)
S%RESIN = S%RES
 
IF (SCARC_DEBUG .GE. 3) THEN
   WRITE (9,*) 'Initial DEFECT '
   WRITE (9,*) 'RES=', S%RES
   WRITE (9,*) 'S%X:'
   WRITE (9, '(a,10f12.6)') '2: S%SP_LOCAL_ALL=', (S%SP_LOCAL_ALL(I), I=1, &
    NMESHES)
   CALL flush (9)
CALL SCARC_COMPARE_SINGLE (S%X, 'SARC', 1, 'X3    ',0)
CALL SCARC_COMPARE_SINGLE (S%F, 'SARC', 1, 'F3    ',0)
CALL SCARC_COMPARE_SINGLE (S%R, 'SARC', 1, 'R3    ',0)
ENDIF
 
IF (S%RESIN.lt.1.0E-16) THEN
   ICONV=1
   ITE=0
   S%CAPPA=0.0_EB
   goto 12345
ENDIF

 
IF (S%BPRECON) THEN
 
   IF (SCARC_PRECON == 1) THEN
      WRITE (*,*) 'IPRECON=1 not yet implemented!'
   ELSE IF (SCARC_PRECON == 2) THEN
      CALL SCARC_VECCOPY2D (S%R, S%G)
      CALL SCARC_PRECON_JACOBI2D (S%AGLOB, S%G)
   ELSE IF (SCARC_PRECON == 3) THEN
      CALL SCARC_VECCOPY2D (S%R, S%G)
      CALL SCARC_PRECON_GS2D (S%AGLOB, S%G)
   ELSE IF (SCARC_PRECON == 4) THEN
      CALL SCARC_VECCOPY2D (S%R, S%G)
      CALL SCARC_PRECON_SSOR2D (S%AGLOB, S%G)
   ENDIF
 
   CALL SCARC_SCALPROD2D (S%R, S%G, S%SIGMA0, NM)
 
ELSE
   S%SIGMA0 = S%RESIN * S%RESIN
ENDIF
 
IF (SCARC_DEBUG .GE. 3) THEN
   WRITE (9,*) 'After initial Precon '
   WRITE (9,*) 'Sigma0=', S%SIGMA0
   WRITE (9,*) 'S%X:'
   WRITE (9, '(a,10f12.6)') '4: S%SP_LOCAL_ALL=', (S%SP_LOCAL_ALL(I), I=1, &
    NMESHES)
   CALL SCARC_COMPARE_SINGLE (S%R, 'SARC', 1, 'R4    ',0)
   CALL SCARC_COMPARE_SINGLE (S%G, 'SARC', 1, 'G4    ',0)
   CALL flush (9)
ENDIF
 
CALL SCARC_VECADD2D (S%G, S%D,-1.0_EB, 0.0_EB)
 
IF (SCARC_DEBUG .GE. 3) THEN
   WRITE (9,*) 'After initial vecadd S%D'
   WRITE (9, '(a,10f12.6)') '5: S%SP_LOCAL_ALL=', (S%SP_LOCAL_ALL(I), I=1, &
    NMESHES)
   CALL flush (9)
ENDIF
CALL SCARC_COMPARE_SINGLE (S%D, 'SARC', 1, 'D4    ',0)
 
 
!!!
!!! start defect correction loop
!!!
CG_LOOP: DO ITE = 1, SCARC_NIT
 
   IF (SCARC_DEBUG .GE. 3) THEN
      WRITE (9,*) '=========================== STARTING CG-LOOP ', ITE
   ENDIF
 
 
!!!  calculate new defect and get L2-norm of it
   IF (SCARC_DEBUG.GE.1) THEN
      CALL SCARC_COMPARE_SINGLE (S%D, 'SARC', 1, 'D10   ',0)
      CALL SCARC_COMPARE_SINGLE (S%Y, 'SARC', 1, 'Y10   ',0)
   ENDIF

   CALL SCARC_MATVEC_MUL2D (S%AGLOB, S%D, S%Y, 1.0_EB, 0.0_EB, NM, ITYPE, 5)
   CALL SCARC_SCALPROD2D (S%D, S%Y, S%ALPHA, NM)

   IF (SCARC_DEBUG.GE.3) THEN
      CALL SCARC_COMPARE_SINGLE (S%D, 'SARC', 1, 'D11   ',0)
      CALL SCARC_COMPARE_SINGLE (S%Y, 'SARC', 1, 'Y11   ',0)
   ENDIF
 
   S%ALPHA = S%SIGMA0 / S%ALPHA
 
   IF (SCARC_DEBUG.GE.3) THEN
      CALL SCARC_COMPARE_SINGLE (S%D, 'SARC', 1, 'D12   ',0)
      CALL SCARC_COMPARE_SINGLE (S%X, 'SARC', 1, 'X12   ',0)
      CALL SCARC_COMPARE_SINGLE (S%Y, 'SARC', 1, 'Y12   ',0)
   ENDIF
   CALL SCARC_VECADD2D (S%D, S%X, S%ALPHA, 1.0_EB)
   CALL SCARC_VECADD2D (S%Y, S%R, S%ALPHA, 1.0_EB)
   CALL SCARC_L2NORM2D (S%R, S%RES, NM, ITYPE)

   IF (SCARC_DEBUG.GE.3) THEN
      CALL SCARC_COMPARE_SINGLE (S%D, 'SARC', 1, 'D13   ',0)
      CALL SCARC_COMPARE_SINGLE (S%X, 'SARC', 1, 'X13   ',0)
      CALL SCARC_COMPARE_SINGLE (S%Y, 'SARC', 1, 'Y13   ',0)
   ENDIF
 
 
   IF (SCARC_DEBUG .GE. 3) THEN
      WRITE (9,*) 'New DEFECT in iteration ', ITE, ': ALPHA=', S%ALPHA
      WRITE (9,*) 'New Iterate in iteration ', ITE
      WRITE (9, '(a10,f22.16)') 'ALPHA=', S%ALPHA
      WRITE (9, '(a10,f22.16)') 'SIGMA0=', S%SIGMA0
      WRITE (9, '(a10,f22.16)') 'RES=', S%RES
      WRITE (9, '(a,10f12.6)') '6: S%SP_LOCAL_ALL=', (S%SP_LOCAL_ALL(I), I=1, &
       NMESHES)
      CALL SCARC_COMPARE_SINGLE (S%Y, 'SARC', 1, 'Y5    ',0)
      CALL SCARC_COMPARE_SINGLE (S%R, 'SARC', 1, 'R5    ',0)
      CALL flush (9)
   ENDIF
   CALL SCARC_COMPARE_SINGLE (S%X, 'SARC', 1, 'XITE   ',0)
 
 
!if (MOD(ITE,1)==0) write(*,'(a,i4,a,f25.16,a,f25.16)') '### SCARC_CG, ite=',ITE,' RES=',S%RES,'  ####'
   IF (NM == 1 .AND. SCARC_DEBUG .GE. 1) write (9, '(a,i4,a,f25.16,a,f25.16)') '&
  &### SCARC_CG, ite=', ITE, ' RES=', S%RES, '  ####'
 
!!!  stop in case of convergence or divergence
   IF (IREL .EQ. 1) THEN
      IF (S%RES .LE. S%RESIN*SCARC_EPS) ICONV = 1
   ELSE
      IF (S%RES .LE. SCARC_EPS .AND. S%RES .LE. S%RESIN*SCARC_EPS_REL) ICONV = 1
   ENDIF
   IF (S%RES .GT. 100000000000.0_EB) ICONV = - 1
   IF (Abs(ICONV) == 1) EXIT CG_LOOP
 
 
!!! preconditioning
   IF (S%BPRECON) THEN
 
      IF (SCARC_PRECON == 1) THEN
         CALL SCARC_RECEIVE (2, 0)
         CALL SCARC_EXCHANGE (2, 0)
         WRITE (*,*) 'IPRECON=1 not yet implemented!'
!     CALL SCARC_CG0(S%R,NM)
      ELSE IF (SCARC_PRECON == 2) THEN
         CALL SCARC_VECCOPY2D (S%R, S%G)
         CALL SCARC_PRECON_JACOBI2D (S%AGLOB, S%G)
      ELSE IF (SCARC_PRECON == 3) THEN
         CALL SCARC_VECCOPY2D (S%R, S%G)
         CALL SCARC_PRECON_GS2D (S%AGLOB, S%G)
      ELSE IF (SCARC_PRECON == 4) THEN
         CALL SCARC_VECCOPY2D (S%R, S%G)
         CALL SCARC_PRECON_SSOR2D (S%AGLOB, S%G)
      ENDIF
 
      CALL SCARC_SCALPROD2D (S%R, S%G, S%SIGMA1, NM)
 
   ELSE
      S%SIGMA1 = S%RES * S%RES
   ENDIF
 
   CALL SCARC_COMPARE_SINGLE (S%G, 'SARC', 1, 'G6    ',0)
 
   S%GAMMA = S%SIGMA1 / S%SIGMA0
   S%SIGMA0 = S%SIGMA1
 
   IF (SCARC_DEBUG .GE. 3) THEN
      WRITE (9,*) 'After Precon in iteration ', ITE
      WRITE (9,*) 'Sigma0=', S%SIGMA0
      WRITE (9,*) 'Gamma =', S%GAMMA
      CALL SCARC_COMPARE_SINGLE (S%R, 'SARC', 1, 'R4    ',0)
      CALL flush (9)
   ENDIF
 
   CALL SCARC_VECADD2D (S%G, S%D,-1.0_EB, S%GAMMA)
 
   IF (SCARC_DEBUG .GE. 3) CALL SCARC_COMPARE_SINGLE (S%D, 'SARC', 1, 'D7    ',0)
 
ENDDO CG_LOOP
 
IF (SCARC_DEBUG .GE. 3) write (9,*) 'AFTER CG_LOOP'
 
 
!!! 'bad end'
IF (ICONV ==-1) THEN
   ITE = - 1
   S%CAPPA = 1.0_EB
 
!!! 'good end'
ELSE
   IF (S%RESIN .GE. 1.0-70_EB) THEN
      S%RESIN = S%RES / S%RESIN
      S%CAPPA = S%RESIN ** (1.0_EB/ITE)
   ELSE IF (ITE == 0) THEN
      S%CAPPA = 0.0_EB
   ENDIF
ENDIF
 
 
12345 CONTINUE
IF (SCARC_DEBUG .GE. 1) THEN
   WRITE (9, '(a,i4,a,f12.6,a,f12.5)') 'SCARC_CG: Ite=', ITE, ': Res=', S%RES, '&
  &: Kappa=', S%CAPPA
   CALL flush (9)
ENDIF
IF (NM == 1) write (*, '(a20,i4,a,e12.6,a,e12.6)') 'SCARC_CG: Ite=', ITE, ': Res&
  &=', S%RES, ': Kappa=', S%CAPPA
CALL flush (6)
 
 
!do k=0,kbar+1
!  do i=0,ibar+1
!    M%H(I,1,K)=S%X(I,1,K)
!  enddo
!enddo
 
DO K = 1, KBAR
   DO I = 1, IBAR
      PRHS (I, K) = S%X(I, 1, K)
   ENDDO
ENDDO
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'Leaving SCARC_CG2D'
   CALL flush (9)
ENDIF
CALL SCARC_COMPARE_SINGLE (S%X, 'SARC', 1, 'XFINAL',0)


 
END SUBROUTINE SCARC_CG2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform global cg-method based on global possion-matrix in 3D
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CG3D (NM, ITRN, POIS_PTB, SAVE1, WORK, HX)
 
INTEGER :: ITRN, NM
INTEGER :: ITE, IREL, ICONV, ITYPE, I, J, K, IZERO, NOM, IERR, IPRECON, III, JJJ
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
 
REAL (EB), DIMENSION (:) :: WORK
REAL (EB), DIMENSION (0:) :: HX
REAL (EB), DIMENSION (-3:) :: SAVE1
REAL (EB) :: POIS_PTB
!REAL(EB), POINTER, DIMENSION(:,:,:) :: H, PRHS
!REAL(EB), POINTER, DIMENSION(:,:):: BXS, BXF, BYS, BYF, BZS, BZF
!REAL(EB), DIMENSION(:,:):: BXS, BXF, BYS, BYF, BZS, BZF
!REAL(EB), DIMENSION(:)  :: SAVE1, WORK, HX
 
TYPE (MESH_TYPE), POINTER :: M
 
INCLUDE 'mpif.h'
 
S => SCARC (NM)
M => MESHES (NM)
 
N_REQ_MV = 0
IERR = 0
 
 
! initialize working vectors
S%G = 0.0_EB
S%Y = 0.0_EB
S%R = 0.0_EB
S%D = 0.0_EB
 
 
DO K = 0, KBAR + 1
   DO J = 0, JBAR + 1
      DO I = 0, IBAR + 1
         S%X (I, J, K) = M%H(I, J, K)
      ENDDO
   ENDDO
ENDDO
 
S%F = 0.0_EB
DO K = 1, KBAR + 1
   DO J = 1, JBAR + 1
      DO I = 1, IBAR + 1
         S%F (I, J, K) = M%PRHS(I, J, K)
      ENDDO
   ENDDO
ENDDO
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'Starting SCARC_CG3D'
   WRITE (9,*) 'IBAR=:', IBAR
   WRITE (9,*) 'JBAR=:', JBAR
   WRITE (9,*) 'KBAR=:', KBAR
   CALL SCARC_COMPARE_SINGLE (M%H, 'SARC', 1, 'H1    ',0)
   CALL SCARC_COMPARE_SINGLE (S%X, 'SARC', 1, 'X1    ',0)
   CALL SCARC_COMPARE_SINGLE (S%F, 'SARC', 1, 'F1    ',0)
   CALL flush (9)
ENDIF
 
!!!
!!! initialize some method parameters
!!!
ITYPE = NCOM_TYPE_GLOBAL
IREL = 0
ICONV = 0
IPRECON = 2
 
!S%NIT    =5000
!S%EPS    =1.0E-16_EB
!S%EPS_REL=1.0E-2_EB
 
 
!!!
!!! initialize boundary values
!!!
IF (NMESHES == 1) THEN
 
   S%BXS0 = M%BXS
   S%BXF0 = M%BXF
   S%BYS0 = M%BYS
   S%BYF0 = M%BYF
   S%BZS0 = M%BZS
   S%BZF0 = M%BZF
 
ELSE
 
   S%BXS0 = 0.0_EB
   S%BXF0 = 0.0_EB
   S%BYS0 = 0.0_EB
   S%BYF0 = 0.0_EB
   S%BZS0 = 0.0_EB
   S%BZF0 = 0.0_EB
 
   IF (NBR_FACE(NM, 5) == 0) S%BXS0 = M%BXS ! IOR= 1: external boundary
   IF (NBR_FACE(NM, 3) == 0) S%BXF0 = M%BXF ! IOR=-1: external boundary
   IF (NBR_FACE(NM, 2) == 0) S%BYS0 = M%BYS ! IOR= 2: external boundary
   IF (NBR_FACE(NM, 4) == 0) S%BYF0 = M%BYF ! IOR=-2: external boundary
   IF (NBR_FACE(NM, 1) == 0) S%BZS0 = M%BZS ! IOR= 3: external boundary
   IF (NBR_FACE(NM, 6) == 0) S%BZF0 = M%BZF ! IOR=-3: external boundary
 
ENDIF
 
CALL SCARC_SETBDRY3D (S%F, S%BXS0, S%BXF0, S%BYS0, S%BYF0, S%BZS0, S%BZF0, &
 SAVE1)
 
 
!!!
!!! calculate initial defect Y = B - A*X and get l2-norm of it
!!!
CALL SCARC_MATVEC_MUL3D (S%AGLOB, S%X, S%R, 1.0_EB, 0.0_EB, NM, ITYPE, 3)
CALL SCARC_VECADD3D (S%F, S%R,-1.0_EB, 1.0_EB)
CALL SCARC_L2NORM3D (S%R, S%RES, NM, ITYPE)
S%RESIN = S%RES
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'Initial DEFECT '
   WRITE (9,*) 'RES=', S%RES
   WRITE (9,*) 'S%X:'
   CALL SCARC_COMPARE_SINGLE (S%X, 'SARC', 1, 'X3    ',0)
   CALL SCARC_COMPARE_SINGLE (S%F, 'SARC', 1, 'F3    ',0)
   CALL SCARC_COMPARE_SINGLE (S%R, 'SARC', 1, 'R3    ',0)
   CALL flush (9)
ENDIF
 
 
IF (S%BPRECON) THEN
 
   IF (SCARC_PRECON == 1) THEN
      WRITE (*,*) 'IPRECON=1 not yet implemented!'
   ELSE IF (SCARC_PRECON == 2) THEN
      CALL SCARC_VECCOPY3D (S%R, S%G)
      CALL SCARC_PRECON_JACOBI3D (S%AGLOB, S%G)
   ELSE IF (SCARC_PRECON == 3) THEN
      CALL SCARC_VECCOPY3D (S%R, S%G)
      CALL SCARC_PRECON_GS3D (S%AGLOB, S%G)
   ELSE IF (SCARC_PRECON == 4) THEN
      CALL SCARC_VECCOPY3D (S%R, S%G)
      CALL SCARC_PRECON_SSOR3D (S%AGLOB, S%G)
   ENDIF
 
   CALL SCARC_SCALPROD3D (S%R, S%G, S%SIGMA0, NM)
 
ELSE
   S%SIGMA0 = S%RESIN * S%RESIN
ENDIF
 
CALL SCARC_VECADD3D (S%G, S%D,-1.0_EB, 0.0_EB)
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'After initial Precon '
   WRITE (9,*) 'Sigma0=', S%SIGMA0
   WRITE (9,*) 'S%X:'
   CALL SCARC_COMPARE_SINGLE (S%R, 'SARC', 1, 'R4    ',0)
   CALL SCARC_COMPARE_SINGLE (S%G, 'SARC', 1, 'G4    ',0)
   CALL SCARC_COMPARE_SINGLE (S%D, 'SARC', 1, 'D4    ',0)
   CALL flush (9)
ENDIF
 
!!!
!!! start defect correction loop
!!!
CG_LOOP: DO ITE = 1, SCARC_NIT
 
   IF (SCARC_DEBUG .GE. 2) write (9,*) '=========================== STARTING CG-&
  &LOOP ', ITE
 
!!!  calculate new defect and get L2-norm of it
   CALL SCARC_MATVEC_MUL3D (S%AGLOB, S%D, S%Y, 1.0_EB, 0.0_EB, NM, ITYPE, 5)
 
   IF (SCARC_DEBUG .GE. 2) CALL SCARC_COMPARE_SINGLE (S%Y, 'SARC', 1, 'Y4    ',0)
 
   CALL SCARC_SCALPROD3D (S%D, S%Y, S%ALPHA, NM)
   S%ALPHA = S%SIGMA0 / S%ALPHA
 
   IF (SCARC_DEBUG .GE. 2) CALL SCARC_COMPARE_SINGLE (S%Y, 'SARC', 1, 'Y5    ',0)
 
   CALL SCARC_VECADD3D (S%D, S%X, S%ALPHA, 1.0_EB)
   CALL SCARC_VECADD3D (S%Y, S%R, S%ALPHA, 1.0_EB)
 
   CALL SCARC_L2NORM3D (S%R, S%RES, NM, ITYPE)
 
   IF (SCARC_DEBUG .GE. 2) THEN
      WRITE (9,*) 'New DEFECT in iteration ', ITE, ': ALPHA=', S%ALPHA
      WRITE (9,*) 'New Iterate in iteration ', ITE
      WRITE (9,*) 'RES=', S%RES
      CALL flush (9)
   ENDIF
   CALL SCARC_COMPARE_SINGLE (S%X, 'SARC', 1, 'X5    ',0)
   CALL SCARC_COMPARE_SINGLE (S%R, 'SARC', 1, 'R5    ',0)
 
   IF (NM == 1 .AND. SCARC_DEBUG .GE. 2) write (9, '(a,i4,a,f25.16,a,f25.16)') '&
  &### SCARC_CG, ite=', ITE, ' RES=', S%RES, '  ####'
   IF (Mod(ITE, 200) == 0 .AND. SCARC_DEBUG == 1) write (9, '(a,i4,a,f25.16,a,f2&
  &5.16)') '### SCARC_CG, ite=', ITE, ' RES=', S%RES, '  ####'
 
!!!  stop in case of convergence or divergence
   IF (IREL .EQ. 1) THEN
      IF (S%RES .LE. S%RESIN*SCARC_EPS) ICONV = 1
   ELSE
      IF (S%RES .LE. SCARC_EPS .AND. S%RES .LE. S%RESIN*SCARC_EPS_REL) ICONV = 1
   ENDIF
   IF (S%RES .GT. 100000.0_EB) ICONV = - 1
   IF (Abs(ICONV) == 1) EXIT CG_LOOP
 
!!! preconditioning
   IF (S%BPRECON) THEN
 
      IF (SCARC_PRECON == 1) THEN
         CALL SCARC_RECEIVE (2, 0)
         CALL SCARC_EXCHANGE (2, 0)
         WRITE (*,*) 'IPRECON=1 not yet implemented!'
!     CALL SCARC_CG0(S%R,NM)
      ELSE IF (SCARC_PRECON == 2) THEN
         CALL SCARC_VECCOPY3D (S%R, S%G)
         CALL SCARC_PRECON_JACOBI3D (S%AGLOB, S%G)
      ELSE IF (SCARC_PRECON == 3) THEN
         CALL SCARC_VECCOPY3D (S%R, S%G)
         CALL SCARC_PRECON_GS3D (S%AGLOB, S%G)
      ELSE IF (SCARC_PRECON == 4) THEN
         CALL SCARC_VECCOPY3D (S%R, S%G)
         CALL SCARC_PRECON_SSOR3D (S%AGLOB, S%G)
      ENDIF
 
      CALL SCARC_SCALPROD3D (S%R, S%G, S%SIGMA1, NM)
 
   ELSE
      S%SIGMA1 = S%RES * S%RES
   ENDIF
 
   IF (SCARC_DEBUG .GE. 3) CALL SCARC_COMPARE_SINGLE (S%G, 'SARC', 1, 'G6    ',0)
 
   S%GAMMA = S%SIGMA1 / S%SIGMA0
   S%SIGMA0 = S%SIGMA1
 
 
   CALL SCARC_VECADD3D (S%G, S%D,-1.0_EB, S%GAMMA)
 
   IF (SCARC_DEBUG .GE. 2) THEN
      WRITE (9,*) 'After Precon in iteration ', ITE
      WRITE (9,*) 'Sigma0=', S%SIGMA0
      WRITE (9,*) 'Gamma =', S%GAMMA
      CALL SCARC_COMPARE_SINGLE (S%R, 'SARC', 1, 'R7    ',0)
      CALL SCARC_COMPARE_SINGLE (S%D, 'SARC', 1, 'D7    ',0)
      CALL flush (9)
   ENDIF
 
ENDDO CG_LOOP
 
IF (SCARC_DEBUG .GE. 2) write (9,*) 'AFTER CG_LOOP'
 
 
!!! 'bad end'
IF (ICONV ==-1) THEN
   ITE = - 1
   S%CAPPA = 1.0_EB
 
!!! 'good end'
ELSE
   IF (S%RESIN .GE. 1.0-70_EB) THEN
      S%RESIN = S%RES / S%RESIN
      S%CAPPA = S%RESIN ** (1.0_EB/ITE)
   ELSE IF (ITE == 0) THEN
      S%CAPPA = 0.0_EB
   ENDIF
ENDIF
 
 
IF (NM == 1) write (*, '(a20,i4,a,e12.6,a,e12.6)') 'SCARC_CG: Ite=', ITE, ': Res&
  &=', S%RES, ': Kappa=', S%CAPPA
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9, '(a,i4,a,f12.6,a,f12.5)') 'SCARC_CG: Ite=', ITE, ': Res=', S%RES, '&
  &: Kappa=', S%CAPPA
   CALL flush (9)
ENDIF
CALL flush (6)
 
 
DO K = 0, KBAR + 1
   DO J = 0, JBAR + 1
      DO I = 0, IBAR + 1
         M%H (I, J, K) = S%X(I, J, K)
      ENDDO
   ENDDO
ENDDO
 
DO K = 1, KBAR + 1
   DO J = 1, JBAR + 1
      DO I = 1, IBAR + 1
         M%PRHS (I, J, K) = S%X(I, J, K)
      ENDDO
   ENDDO
ENDDO
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) '-------SOL: S%X:'
   CALL SCARC_COMPARE_SINGLE (S%X, 'SARC', 1, 'XFINAL',0)
   CALL flush (9)
ENDIF
 
END SUBROUTINE SCARC_CG3D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform matrix-vector multiplication with full matrix
!!! corresponding to band-wise storage technique
!!!
!!!  Y = A1 * S%AGLOB * X + A2 * Y
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MATVEC_MUL2D (A, X, Y, A1, A2, NM, ITYPE, IMV)
 
REAL (EB) :: A1, A2, X1, X2, X3, X4, X5
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: NM, ITYPE
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: I, J, K, LL, IC, IW, IWW, IERR, NOM, IMV, KK, jj, II
TYPE (MESH_TYPE), POINTER :: M
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print debug messages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG .GT. 2) THEN
   WRITE (9,*) 'BEFORE MATVEC_MUL: '
   WRITE (9,*) 'TWO_D=', TWO_D
   WRITE (9,*) 'KBAR=', KBAR
   WRITE (9,*) 'JBAR=', JBAR
   WRITE (9,*) 'IBAR=', IBAR
   CALL SCARC_COMPARE_SINGLE (X, 'SARC', 1, 'matx1 ',0)
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'maty1 ',0)
WRITE (9,100) SHAPE(X)
100 FORMAT (1X,'The shape of the array is:          ', 7I6)

WRITE (9,110) SIZE(X)
110 FORMAT (1X,'The shape of the array is:          ', I10)

WRITE (9,120) LBOUND(X)
120 FORMAT (1X,'The lower bounds of the array are: ', 7I6)

WRITE (9,130) UBOUND(X)
130 FORMAT (1X,'The upper bounds of the array are: ', 7I6)

write(9,*) 'X(1,1,1)=',X(1,1,1)
ENDIF
 
!X=REAL(NM,EB)
!Y=0.0_EB
!X=1.0_EB
!Y=0.0_EB
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform 2D-matrix-vector-multiplication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO K = 1, KBAR
   DO I = 1, IBAR
      IC = (K-1) * IBAR + I
      IF (K>1) THEN 
        X1=X(I,1,K-1)
      ELSE
        X1=0.0_EB
      ENDIF
      IF (I>1) THEN 
        X2=X(I-1,1,K)
      ELSE
        X2=0.0_EB
      ENDIF
      X3=X(I,1,K)
      IF (I<IBAR) THEN 
        X4=X(I+1,1,K)
      ELSE
        X4=0.0_EB
      ENDIF
      IF (K<KBAR) THEN 
        X5=X(I,1,K+1)
      ELSE
        X5=0.0_EB
      ENDIF


      Y (I, 1, K) = A1 * (  A(IC, 1) * X1  &
                          + A(IC, 2) * X2  &
                          + A(IC, 3) * X3  &
                          + A(IC, 4) * X4  &
                          + A(IC, 5) * X5 ) &
                  + A2 * Y (I, 1, K)
if (scarc_debug.ge.2) then
write(9,*) 'X1=',X1
write(9,*) 'X2=',X2
write(9,*) 'X3=',X3
write(9,*) 'X4=',X4
write(9,*) 'X5=',X5
write(9,*) 'Y(',I,',',k,')=',Y(I,1,k)
endif
      !Y (I, 1, K) = A1 * (  A(IC, 1) * X(I  , 1, K-1)    &
      !                    + A(IC, 2) * X(I-1, 1, K  )    & 
      !                    + A(IC, 3) * X(I  , 1, K  )    &
      !                    + A(IC, 4) * X(I+1, 1, K  )    &
      !                    + A(IC, 5) * X(I  , 1, K+1) )  &
      !            + A2 * Y (I, 1, K)

   ENDDO
ENDDO
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print debug messages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG .GT. 2) THEN
   WRITE (9,*) '------->> HALLO BEFORE SCARC_EXCHANGE: Y:'
   CALL SCARC_COMPARE_SINGLE(Y,'SARC',1,'maty2 ',0)
   CALL flush (9)
ENDIF
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform data exchange at interior boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF ((ITYPE == NCOM_TYPE_GLOBAL) .AND. (NMESHES > 1)) THEN
   N_REQ_MV = 0
   CALL SCARC_RECEIVE (1, IMV)
   CALL SCARC_EXCHANGE (1, IMV)
ENDIF
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print debug messages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG .GT. 2) THEN
   WRITE (9,*) '------->> HALLO AFTER SCARC_EXCHANGE: Y:'
   CALL flush (9)
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'maty1 ',0)
ENDIF
 
END SUBROUTINE SCARC_MATVEC_MUL2D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Perform matrix-vector multiplication with full matrix
!!! corresponding to band-wise storage technique
!!!
!!!  Y = A1 * S%AGLOB * X + A2 * Y
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MATVEC_MUL3D (A, X, Y, A1, A2, NM, ITYPE, IMV)
 
REAL (EB) :: A1, A2
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: NM, ITYPE
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: I, J, K, LL, IC, IW, IWW, IERR, NOM, IMV, KK, jj, II
TYPE (MESH_TYPE), POINTER :: M
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print debug messages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG .GT. 2) THEN
   WRITE (9,*) 'BEFORE MATVEC_MUL: '
!  CALL SCARC_COMPARE_SINGLE(X,'SARC',1,'matX  ',0)
!  CALL SCARC_COMPARE_SINGLE(Y,'SARC',1,'maty  ',0)
   WRITE (9,*) 'TWO_D=', TWO_D
   WRITE (9,*) 'KBAR=', KBAR
   WRITE (9,*) 'JBAR=', JBAR
   WRITE (9,*) 'IBAR=', IBAR
ENDIF
 
!X=REAL(NM,EB)
!Y=0.0_EB
!X=1.0_EB
!Y=0.0_EB
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform 3D-matrix-vector-multiplication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO K = 1, KBAR
   IF (SCARC_DEBUG .GT. 2) write (9,*) 'K=', K, ' ==============================&
  &======='
   DO J = 1, JBAR
      DO I = 1, IBAR
         IC = (K-1) * IBAR * JBAR + (J-1) * IBAR + I
!write(9,'(8f22.16)') y(i,j,k),(A(IC,kk),kk=1,7)
         Y (I, J, K) = A1 * (A(IC, 1)*X(I, J, K-1)+A(IC, 2)*X(I, J-1, K)+A(IC, &
          3)*X(I-1, J, K)+A(IC, 4)*X(I, J, K)+A(IC, 5)*X(I+1, J, K)+A(IC, &
          6)*X(I, J+1, K)+A(IC, 7)*X(I, J, K+1)) + A2 * Y (I, J, K)
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,i3,a,7f6.0,a,f22.1&
        &6)') '(', I, ',', J, ',', K, '):A(', IC, '):', (A(IC, LL), LL=1, 7), '-&
        &-> Y=', Y (I, J, K)
      ENDDO
      IF (SCARC_DEBUG .GT. 2) write (9,*) '-------------------------------------&
     &'
   ENDDO
ENDDO
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print debug messages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG .GT. 2) THEN
   WRITE (9,*) '------->> HALLO BEFORE SCARC_EXCHANGE: Y:'
!  CALL SCARC_COMPARE_SINGLE(Y,'SARC',1,'maty2 ',0)
   CALL flush (9)
ENDIF
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform data exchange at interior boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF ((ITYPE == NCOM_TYPE_GLOBAL) .AND. (NMESHES > 1)) THEN
   N_REQ_MV = 0
   CALL SCARC_RECEIVE (1, IMV)
   CALL SCARC_EXCHANGE (1, IMV)
ENDIF
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! print debug messages !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SCARC_DEBUG .GT. 2) THEN
   WRITE (9,*) '------->> HALLO AFTER SCARC_EXCHANGE: Y:'
!  CALL SCARC_COMPARE_SINGLE(Y,'SARC',1,'maty3 ',0)
   CALL flush (9)
ENDIF
 
END SUBROUTINE SCARC_MATVEC_MUL3D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! add two vectors:
!!!
!!! Y = A1 * X + A2 * Y
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_VECADD2D (X, Y, A1, A2)
 
REAL (EB) :: A1, A2
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: I, J, K
 
 
IF (SCARC_DEBUG .GE. 6) THEN
   WRITE (9,*) 'BEFORE VECADD'
   WRITE (9,*) 'A1=', A1
   WRITE (9,*) 'A2=', A2
   WRITE (9,*) ' X:'
   DO K = KBAR + 1, 1, - 1
      WRITE (9, '(5f10.5)') (X(I, 1, K), I=1, IBAR+1)
   ENDDO
   WRITE (9,*) ' Y:'
   DO K = KBAR + 1, 1, - 1
      WRITE (9, '(5f10.5)') (Y(I, 1, K), I=1, IBAR+1)
   ENDDO
ENDIF
 
 
IF (A2 == 0.0_EB) THEN
   IF (A1 /= 1.0_EB) THEN
      DO K = 1, KBAR
         DO I = 1, IBAR
            Y (I, 1, K) = A1 * X (I, 1, K)
 !WRITE(9,*) 'Y(',I,',1,',K,')=',Y(I,1,K)
         ENDDO
      ENDDO
   ELSE
      DO K = 1, KBAR
         DO I = 1, IBAR
            Y (I, 1, K) = X (I, 1, K)
 !WRITE(9,*) 'Y(',I,',1,',K,')=',Y(I,1,K)
         ENDDO
      ENDDO
   ENDIF
ELSE
   IF (A1 /= 1.0_EB) THEN
      DO K = 1, KBAR
         DO I = 1, IBAR
            Y (I, 1, K) = A1 * X (I, 1, K) + A2 * Y (I, 1, K)
 !WRITE(9,*) 'Y(',I,',1,',K,')=',Y(I,1,K)
         ENDDO
      ENDDO
   ELSE
      DO K = 1, KBAR
         DO I = 1, IBAR
            Y (I, 1, K) = X (I, 1, K) + A2 * Y (I, 1, K)
 !WRITE(9,*) 'Y(',I,',1,',K,')=',Y(I,1,K)
         ENDDO
      ENDDO
   ENDIF
ENDIF
 
 
IF (SCARC_DEBUG .GE. 6) THEN
   WRITE (9,*) 'AFTER VECADD'
   WRITE (9,*) ' Y:'
   DO K = KBAR + 1, 1, - 1
      WRITE (9, '(5f10.5)') (Y(I, 1, K), I=1, IBAR+1)
   ENDDO
ENDIF
 
END SUBROUTINE SCARC_VECADD2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! add two vectors:
!!!
!!! Y = A1 * X + A2 * Y
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_VECADD3D (X, Y, A1, A2)
 
REAL (EB) :: A1, A2
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: I, J, K
 
 
IF (SCARC_DEBUG .GE. 6) THEN
   WRITE (9,*) 'BEFORE VECADD'
   WRITE (9,*) 'A1=', A1
   WRITE (9,*) 'A2=', A2
   WRITE (9,*) ' X:'
   DO K = KBAR + 1, 1, - 1
      WRITE (9, '(5f10.5)') ((X(I, J, K), I=1, IBAR+1), J=JBAR+1, 1,-1)
   ENDDO
   WRITE (9,*) ' Y:'
   DO K = KBAR + 1, 1, - 1
      WRITE (9, '(5f10.5)') ((Y(I, J, K), I=1, IBAR+1), J=JBAR+1, 1,-1)
   ENDDO
ENDIF
 
 
IF (A2 == 0.0_EB) THEN
   IF (A1 /= 1.0_EB) THEN
      DO K = 1, KBAR
         DO J = 1, JBAR
            DO I = 1, IBAR
               Y (I, J, K) = A1 * X (I, J, K)
!WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO K = 1, KBAR
         DO J = 1, JBAR
            DO I = 1, IBAR
               Y (I, J, K) = X (I, J, K)
!WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ELSE
   IF (A1 /= 1.0_EB) THEN
      DO K = 1, KBAR
         DO J = 1, JBAR
            DO I = 1, IBAR
               Y (I, J, K) = A1 * X (I, J, K) + A2 * Y (I, J, K)
!WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO K = 1, KBAR
         DO J = 1, KBAR
            DO I = 1, IBAR
               Y (I, J, K) = X (I, J, K) + A2 * Y (I, J, K)
!WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDIF
 
 
IF (SCARC_DEBUG .GE. 6) THEN
   WRITE (9,*) 'AFTER VECADD'
   WRITE (9,*) ' Y:'
   DO K = KBAR + 1, 1, - 1
      WRITE (9, '(5f10.5)') ((Y(I, J, K), I=1, IBAR+1), J=JBAR+1, 1,-1)
   ENDDO
ENDIF
 
 
END SUBROUTINE SCARC_VECADD3D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! copy vector
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_VECCOPY2D (X, Y)
 
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: I, J, K
 
 
IF (SCARC_DEBUG .GE. 6) THEN
   WRITE (9,*) 'BEFORE VECCOPY2D'
   WRITE (9,*) ' X:'
   DO K = KBAR, 1, - 1
      WRITE (9, '(4f10.4)') (X(I, 1, K), I=1, IBAR)
   ENDDO
   WRITE (9,*) ' Y:'
   DO K = KBAR, 1, - 1
      WRITE (9, '(4f10.4)') (Y(I, 1, K), I=1, IBAR)
   ENDDO
ENDIF
 
 
DO K = 1, KBAR
   DO I = 1, IBAR
      Y (I, 1, K) = X (I, 1, K)
!    WRITE(9,*) 'Y(',I,',1,',K,')=',Y(I,1,K)
   ENDDO
ENDDO
 
 
IF (SCARC_DEBUG .GE. 6) THEN
   WRITE (9,*) 'AFTER VECCOPY'
   WRITE (9,*) ' Y:'
   DO K = KBAR, 1, - 1
      WRITE (9, '(4f10.4)') (Y(I, 1, K), I=1, IBAR)
   ENDDO
ENDIF
 
END SUBROUTINE SCARC_VECCOPY2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! copy vector
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_VECCOPY3D (X, Y)
 
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: I, J, K
 
 
IF (SCARC_DEBUG .GE. 6) THEN
   WRITE (9,*) 'BEFORE VECCOPY'
   WRITE (9,*) ' X:'
   DO K = KBAR, 1, - 1
      WRITE (9, '(4f10.4)') ((X(I, J, K), I=1, IBAR), J=JBAR, 1,-1)
   ENDDO
   WRITE (9,*) ' Y:'
   DO K = KBAR, 1, - 1
      WRITE (9, '(4f10.4)') ((Y(I, J, K), I=1, IBAR), J=JBAR, 1,-1)
   ENDDO
ENDIF
 
 
DO K = 1, KBAR
   DO J = 1, JBAR
      DO I = 1, IBAR
         Y (I, J, K) = X (I, J, K)
!      WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
      ENDDO
   ENDDO
ENDDO
 
 
IF (SCARC_DEBUG .GE. 6) THEN
   WRITE (9,*) 'AFTER VECCOPY3D'
   WRITE (9,*) ' Y:'
   DO K = KBAR, 1, - 1
      WRITE (9, '(4f10.4)') ((Y(I, J, K), I=1, IBAR), J=JBAR, 1,-1)
   ENDDO
ENDIF
 
END SUBROUTINE SCARC_VECCOPY3D
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Compute scalarproduct of vector X and vector Y
!!!
!!! SP = (X,Y)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SCALPROD2D (X, Y, SP, NM)
 
REAL (EB) :: SP
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: I, J, K, ITYPE, IM, NM, IERR
 
INCLUDE 'mpif.h'
 
S => SCARC (NM)
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'SCAL_PROD2D'
   WRITE (9,*) 'NMESHES=', NMESHES
   WRITE (9,*) 'IBAR=', IBAR
   WRITE (9,*) 'KBAR=', KBAR
   CALL SCARC_COMPARE_SINGLE (X, 'SARC', 1, 'X      ',0)
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y      ',0)
ENDIF
 
 
 
SP = 0.0_EB
DO K = 1, KBAR
   DO I = 1, IBAR
      SP = SP + X (I, 1, K) * Y (I, 1, K)
      IF (SCARC_DEBUG .GE. 4) write (9, 1001) I, K, X (I, 1, K), Y (I, 1, K), SP
   ENDDO
ENDDO
 
IF (SCARC_DEBUG .GE. 4) write (9,*) 'SP=', SP
 
IF (NMESHES > 1) THEN
 
   S%SP_LOCAL_ALL = 0.0_EB
   S%SP_LOCAL_ALL (NM) = SP
 
   S%SP_GLOBAL_ALL = 0.0_EB
   CALL MPI_ALLREDUCE (S%SP_LOCAL_ALL(1), S%SP_GLOBAL_ALL(1), NMESHES, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
 
   IF (SCARC_DEBUG .GE. 6) write (9,*) 'S%SP_LOCAL_ALL='
   IF (SCARC_DEBUG .GE. 6) write (9, '(f22.16)') (S%SP_LOCAL_ALL(I), I=1, &
    NMESHES)
   IF (SCARC_DEBUG .GE. 6) write (9,*) 'S%SP_GLOBAL_ALL='
   IF (SCARC_DEBUG .GE. 6) write (9, '(f22.16)') (S%SP_GLOBAL_ALL(I), I=1, &
    NMESHES)
 
   SP = 0.0_EB
   DO IM = 1, NMESHES
      SP = SP + S%SP_GLOBAL_ALL(IM)
   ENDDO
 
ENDIF
IF (SCARC_DEBUG .GE. 6) write (9,*) 'SCALPROD2D: SP=', SP
 
1001     FORMAT ('sp=sp+[', i4, ',1,', i4, ']:[', f22.16, '*', f22.16, ']:', f22.16)
1002     FORMAT ('sp=sp+[', i4, ',', i4, ',', i4, ']:[', f12.6, ']:', f24.12)
 
END SUBROUTINE SCARC_SCALPROD2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Compute scalarproduct of vector X and vector Y
!!!
!!! SP = (X,Y)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SCALPROD3D (X, Y, SP, NM)
 
REAL (EB) :: SP
REAL (EB), POINTER, DIMENSION (:, :, :) :: X, Y
INTEGER :: I, J, K, ITYPE, IM, NM, IERR
 
INCLUDE 'mpif.h'
 
IF (SCARC_DEBUG .GE. 6) THEN
   WRITE (9,*) 'SCAL_PROD3D'
   CALL SCARC_COMPARE_SINGLE (X, 'SARC', 1, 'X      ',0)
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y      ',0)
ENDIF
 
 
SP = 0.0_EB
DO K = 1, KBAR
   DO J = 1, JBAR
      DO I = 1, IBAR
         SP = SP + X (I, J, K) * Y (I, J, K)
         IF (SCARC_DEBUG .GE. 4) write (9, 1002) I, J, K, X (I, J, K) * Y (I, J, &
          K), SP
      ENDDO
   ENDDO
ENDDO
 
IF (SCARC_DEBUG .GT. 2) write (9,*) 'SP=', SP
 
IF (NMESHES > 1) THEN
 
   S%SP_LOCAL_ALL = 0.0_EB
   S%SP_LOCAL_ALL (NM) = SP
 
   S%SP_GLOBAL_ALL = 0.0_EB
   CALL MPI_ALLREDUCE (S%SP_LOCAL_ALL(1), S%SP_GLOBAL_ALL(1), NMESHES, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
 
   IF (SCARC_DEBUG .GE. 6) write (9,*) 'S%SP_LOCAL_ALL='
   IF (SCARC_DEBUG .GE. 6) write (9, '(f22.16)') (S%SP_LOCAL_ALL(I), I=1, &
    NMESHES)
   IF (SCARC_DEBUG .GE. 6) write (9,*) 'S%SP_GLOBAL_ALL='
   IF (SCARC_DEBUG .GE. 6) write (9, '(f22.16)') (S%SP_GLOBAL_ALL(I), I=1, &
    NMESHES)
 
   SP = 0.0_EB
   DO IM = 1, NMESHES
      SP = SP + S%SP_GLOBAL_ALL(IM)
   ENDDO
   IF (SCARC_DEBUG .GE. 6) write (9,*) 'SP=', SP
 
ENDIF
 
1001     FORMAT ('sp=sp+[', i4, ',', i4, ']:[', f12.6, '*', f12.6, ']:', f12.6)
1002     FORMAT ('sp=sp+[', i4, ',', i4, ',', i4, ']:[', f12.6, ']:', f24.12)
 
END SUBROUTINE SCARC_SCALPROD3D
 
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Compute L2-norm of vector X:   SP = ||X||
!!! depending on ITYPE:
!!! locally: only local cells are considered
!!! globally: all cells (globally) are considered
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_L2NORM2D (X, SP, NM, ITYPE)
 
REAL (EB) :: SP
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: NM, ITYPE
INTEGER :: IC, IM, IERR, I, J, K
 
INCLUDE 'mpif.h'
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) '========= L2NORM , ITYPE=', ITYPE
   WRITE (9,*) 'NMESHES=', NMESHES
   WRITE (9,*) 'NM=', NM
   CALL SCARC_COMPARE_SINGLE (X, 'SARC', 1, 'X     ',0)
ENDIF
 
S => SCARC (NM)
 
!!! build local scalar product (x,x)
S%SP_LOCAL = 0.0_EB
DO K = 1, KBAR
   DO I = 1, IBAR
      S%SP_LOCAL = S%SP_LOCAL + X (I, 1, K) * X (I, 1, K)
      IF (SCARC_DEBUG .GE. 4) write (9, 1002) I, K, X (I, 1, K), X (I, 1, K), &
       S%SP_LOCAL
   ENDDO
ENDDO
 
IF (SCARC_DEBUG .GE. 2) write (9,*) 'S%SP_LOCAL=', S%SP_LOCAL
 
!!! scale with number of local cells
IF (ITYPE == NCOM_TYPE_LOCAL) THEN
 
   SP = Sqrt (S%SP_LOCAL/REAL(S%NCELLS_LOCAL, EB))
 
   IF (SCARC_DEBUG .GE. 2) THEN
      WRITE (9,*) 'S%NCELLS_LOCAL=', S%NCELLS_LOCAL
      WRITE (9,*) 'S%SP_LOCAL=', S%SP_LOCAL
      WRITE (9,*) 'NMESHES=', NMESHES
   ENDIF
 
!!! sum up globally and scale with global number of nodes
ELSE IF (ITYPE == NCOM_TYPE_GLOBAL) THEN
 
   IF (NMESHES > 1) THEN
 
      S%SP_LOCAL_ALL = 0.0_EB
      S%SP_LOCAL_ALL (NM) = S%SP_LOCAL
 
      S%SP_GLOBAL_ALL = 0.0_EB
      CALL MPI_ALLREDUCE (S%SP_LOCAL_ALL(1), S%SP_GLOBAL_ALL(1), NMESHES, &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
 
      IF (SCARC_DEBUG .GE. 2) write (9,*) 'S%SP_LOCAL_ALL='
      IF (SCARC_DEBUG .GE. 2) write (9, '(f22.16)') (S%SP_LOCAL_ALL(I), I=1, &
       NMESHES)
      IF (SCARC_DEBUG .GE. 2) write (9,*) 'S%SP_GLOBAL_ALL='
      IF (SCARC_DEBUG .GE. 2) write (9, '(f22.16)') (S%SP_GLOBAL_ALL(I), I=1, &
       NMESHES)
 
      S%SP_GLOBAL = 0.0_EB
      DO IM = 1, NMESHES
         S%SP_GLOBAL = S%SP_GLOBAL + S%SP_GLOBAL_ALL(IM)
         IF (SCARC_DEBUG .GE. 4) write (9,*) 'S%SP_GLOBAL=', S%SP_GLOBAL
      ENDDO
   ELSE
      S%SP_GLOBAL = S%SP_LOCAL
   ENDIF
 
   SP = Sqrt (S%SP_GLOBAL/REAL(S%NCELLS_GLOBAL, EB))
 
   IF (SCARC_DEBUG .GE. 6) write (9,*) 'S%NCELLS_GLOBAL=', S%NCELLS_GLOBAL
   IF (SCARC_DEBUG .GE. 6) write (9,*) 'L2NORM2D  : SP=', SP
 
ELSE
 
   WRITE (*,*) 'Wrong type for SCARC_L2NORM2D ', ITYPE
 
ENDIF
 
1002     FORMAT ('sp=sp+[', i4, ',', i4, ']:[', f22.16, '*', f22.16, ']:', f22.16)
END SUBROUTINE SCARC_L2NORM2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Compute L2-norm of vector X:   SP = ||X||
!!! depending on ITYPE:
!!! locally: only local cells are considered
!!! globally: all cells (globally) are considered
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_L2NORM3D (X, SP, NM, ITYPE)
 
REAL (EB) :: SP
REAL (EB), POINTER, DIMENSION (:, :, :) :: X
INTEGER :: NM, ITYPE
INTEGER :: IC, IM, IERR, I, J, K
 
INCLUDE 'mpif.h'
 
IF (SCARC_DEBUG .GE. 2) write (9,*) '========= L2NORM , ITYPE=', ITYPE
IF (SCARC_DEBUG .GE. 2) write (9,*) 'NMESHES=', NMESHES
 
S => SCARC (NM)
 
!!! build local scalar product (x,x)
S%SP_LOCAL = 0.0_EB
DO K = 1, KBAR
   DO J = 1, JBAR
      DO I = 1, IBAR
         S%SP_LOCAL = S%SP_LOCAL + X (I, J, K) * X (I, J, K)
!if (SCARC_DEBUG.ge.2) write(9,1002) I,J,K,X(I,J,K)*X(I,J,K),SP
      ENDDO
   ENDDO
ENDDO
 
IF (SCARC_DEBUG .GE. 2) write (9,*) 'S%SP_LOCAL=', S%SP_LOCAL
 
!!! scale with number of local cells
IF (ITYPE == NCOM_TYPE_LOCAL) THEN
 
   SP = Sqrt (S%SP_LOCAL/REAL(S%NCELLS_LOCAL, EB))
 
   IF (SCARC_DEBUG .GE. 2) write (9,*) 'S%NCELLS_LOCAL=', S%NCELLS_LOCAL
   IF (SCARC_DEBUG .GE. 2) write (9,*) 'S%SP_LOCAL=', S%SP_LOCAL
   IF (SCARC_DEBUG .GE. 2) write (9,*) 'NMESHES=', NMESHES
 
!!! sum up globally and scale with global number of nodes
ELSE IF (ITYPE == NCOM_TYPE_GLOBAL) THEN
 
   IF (NMESHES > 1) THEN
 
      S%SP_LOCAL_ALL = 0.0_EB
      S%SP_LOCAL_ALL (NM) = S%SP_LOCAL
 
      S%SP_GLOBAL_ALL = 0.0_EB
      CALL MPI_ALLREDUCE (S%SP_LOCAL_ALL(1), S%SP_GLOBAL_ALL(1), NMESHES, &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
 
      IF (SCARC_DEBUG .GE. 2) write (9,*) 'S%SP_LOCAL_ALL='
      IF (SCARC_DEBUG .GE. 2) write (9, '(f22.16)') (S%SP_LOCAL_ALL(I), I=1, &
       NMESHES)
      IF (SCARC_DEBUG .GE. 2) write (9,*) 'S%SP_GLOBAL_ALL='
      IF (SCARC_DEBUG .GE. 2) write (9, '(f22.16)') (S%SP_GLOBAL_ALL(I), I=1, &
       NMESHES)
 
      S%SP_GLOBAL = 0.0_EB
      DO IM = 1, NMESHES
         S%SP_GLOBAL = S%SP_GLOBAL + S%SP_GLOBAL_ALL(IM)
         IF (SCARC_DEBUG .GE. 2) write (9,*) 'S%SP_GLOBAL=', S%SP_GLOBAL
      ENDDO
   ELSE
      S%SP_GLOBAL = S%SP_LOCAL
   ENDIF
 
   SP = Sqrt (S%SP_GLOBAL/REAL(S%NCELLS_GLOBAL, EB))
 
   IF (SCARC_DEBUG .GE. 2) write (9,*) 'S%NCELLS_GLOBAL=', S%NCELLS_GLOBAL
   IF (SCARC_DEBUG .GE. 2) write (9,*) 'SP=', SP
 
ELSE
 
   WRITE (*,*) 'Wrong type for SCARC_L2NORM ', ITYPE
 
ENDIF
 
END SUBROUTINE SCARC_L2NORM3D
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Jacobi preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_JACOBI2D (A, Y)
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: YOLD
INTEGER :: I, J, K, IC
 
IF (SCARC_DEBUG .GE. 4) THEN
   WRITE (9,*) 'BEFORE PRECON_JACOBI2D: Y'
   WRITE (9,*) 'Y:'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
 
DO K = 1, KBAR
   DO I = 1, IBAR
      IC = (K-1) * IBAR + I
      Y (I, 1, K) = Y (I, 1, K) / A (IC, 3)
      IF (SCARC_DEBUG .GE. 6) write (9, '(a,i3,a,f8.2,a,f18.12)') 'A(', IC, ')=',&
        A (IC, 3), ': Y=', Y (I, 1, K)
   ENDDO
   IF (SCARC_DEBUG .GE. 6) write (9,*) '========================================&
  &======'
ENDDO
 
 
IF (SCARC_DEBUG .GE. 4) THEN
   WRITE (9,*)
   WRITE (9,*) 'AFTER  PRECON_JACOBI: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
END SUBROUTINE SCARC_PRECON_JACOBI2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Jacobi preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_JACOBI3D (A, Y)
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: YOLD
INTEGER :: I, J, K, IC
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'BEFORE PRECON_JACOBI3D: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
 
DO K = 1, KBAR
   DO J = 1, JBAR
      DO I = 1, IBAR
         IC = (K-1) * IBAR * JBAR + (J-1) * IBAR + I
         YOLD = Y (I, J, K)
         IF (TWO_D) THEN
            Y (I, J, K) = Y (I, J, K) / A (IC, 3)
         ELSE
            Y (I, J, K) = Y (I, J, K) / A (IC, 4)
            IF (SCARC_DEBUG .GE. 6) write (9, '(a,i3,a,f8.2,a,f18.12)') 'A(', &
             IC, ')=', A (IC, 4), ': Y=', Y (I, J, K)
         ENDIF
      ENDDO
      IF (SCARC_DEBUG .GE. 6) write (9,*) '-------------'
   ENDDO
   IF (SCARC_DEBUG .GE. 6) write (9,*) '========================================&
  &======'
ENDDO
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'AFTER  PRECON_JACOBI3D: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
END SUBROUTINE SCARC_PRECON_JACOBI3D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! GS preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_GS2D (A, Y)
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: YOLD, AUX, OMEGA
INTEGER :: I, J, K, IC, IDIAG
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'BEFORE PRECON_GS2D: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
 
OMEGA = 1.0_EB
 
DO K = 1, KBAR
   DO I = 1, IBAR
      IC = (K-1) * IBAR + I
 
      AUX = AUX + A (IC, 1) * Y (I, 1, K-1) + A (IC, 2) * Y (I-1, 1, K)
      Y (I, 1, K) = (Y(I, 1, K)-AUX*OMEGA) / A (IC, 3)
      IF (SCARC_DEBUG .GE. 6) write (9, '(a,i3,a,f8.2,a,f18.12)') 'A(', IC, ')=',&
        A (IC, 3), ': Y=', Y (I, J, K)
   ENDDO
   IF (SCARC_DEBUG .GE. 6) write (9,*) '========================================&
  &======'
ENDDO
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'AFTER  PRECON_GS: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
END SUBROUTINE SCARC_PRECON_GS2D
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! GS preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_GS3D (A, Y)
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: YOLD, AUX, OMEGA
INTEGER :: I, J, K, IC, IDIAG
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'BEFORE PRECON_GS3D: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
OMEGA = 1.0_EB
 
DO K = 1, KBAR
   DO J = 1, JBAR
      DO I = 1, IBAR
         IC = (K-1) * IBAR * JBAR + (J-1) * IBAR + I
         AUX = 0.0_EB
         AUX = A (IC, 1) * Y (I, J, K-1) + A (IC, 2) * Y (I, J-1, K) + A (IC, 3) &
          * Y (I-1, J, K)
         Y (I, J, K) = (Y(I, J, K)-OMEGA*AUX) / A (IC, 4)
         IF (SCARC_DEBUG .GE. 6) write (9, '(a,i3,a,f8.2,a,f18.12,a,f18.12)') 'A&
        &(', IC, ')=', A (IC, 4), ': AUX=', AUX, ': Y=', Y (I, J, K)
      ENDDO
      IF (SCARC_DEBUG .GE. 6) write (9,*) '-------------'
   ENDDO
ENDDO
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'AFTER  PRECON_GS3D: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
END SUBROUTINE SCARC_PRECON_GS3D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! SSOR preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_SSOR2D (A, Y)
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: YOLD, AUX, OMEGA
INTEGER :: I, J, K, IC, IDIAG
 
IF (SCARC_DEBUG .GE. 4) THEN
   WRITE (9,*) 'BEFORE PRECON_SSOR2D: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
OMEGA = 1.5_EB
 
DO K = 1, KBAR
   DO I = 1, IBAR
      IC = (K-1) * IBAR + I
      AUX = A (IC, 1) * Y (I, 1, K-1) + A (IC, 2) * Y (I-1, 1, K)
      Y (I, 1, K) = (Y(I, 1, K)-AUX*OMEGA) / A (IC, 3)
      IF (SCARC_DEBUG .GE. 4) write (9, '(a,i3,a,f8.2,a,f18.12,a,f18.12)') 'A(', &
       IC, ')=', A (IC, 3), ': AUX=', AUX, ': Y=', Y (I, 1, K)
   ENDDO
   IF (SCARC_DEBUG .GE. 4) write (9,*) '========================================&
  &======'
ENDDO
 
DO K = KBAR, 1, - 1
   DO I = IBAR, 1, - 1
      IC = (K-1) * IBAR + I
      AUX = A (IC, 5) * Y (I, 1, K+1) + A (IC, 4) * Y (I+1, 1, K)
      Y (I, 1, K) = Y (I, 1, K) - AUX * OMEGA / A (IC, 3)
      IF (SCARC_DEBUG .GE. 4) write (9, '(a,i3,a,f8.2,a,f18.12,a,f18.12)') 'A(', &
       IC, ')=', A (IC, 3), ': AUX=', AUX, ': Y=', Y (I, 1, K)
   ENDDO
   IF (SCARC_DEBUG .GE. 4) write (9,*) '========================================&
  &======'
ENDDO
 
IF (SCARC_DEBUG .GE. 4) THEN
   WRITE (9,*) 'AFTER  PRECON_SSOR: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
END SUBROUTINE SCARC_PRECON_SSOR2D
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! SSOR preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_SSOR3D (A, Y)
REAL (EB), POINTER, DIMENSION (:, :) :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: Y
REAL (EB) :: YOLD, AUX, OMEGA
INTEGER :: I, J, K, IC, IDIAG
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'BEFORE PRECON_SSOR: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
 
OMEGA = 1.6_EB
 
DO K = 1, KBAR
   DO J = 1, JBAR
      DO I = 1, IBAR
         IC = (K-1) * IBAR * JBAR + (J-1) * IBAR + I
 
         AUX = A (IC, 1) * Y (I, J, K-1) + A (IC, 2) * Y (I, J-1, K) + A (IC, 3) &
          * Y (I-1, J, K)
         Y (I, J, K) = (Y(I, J, K)-AUX*OMEGA) / A (IC, 4)
         IF (SCARC_DEBUG .GE. 6) write (9, '(a,i3,a,f8.2,a,f18.12,a,f18.12)') 'A&
        &(', IC, ')=', A (IC, 4), ': AUX=', AUX, ': Y=', Y (I, J, K)
      ENDDO
      IF (SCARC_DEBUG .GE. 6) write (9,*) '-------------'
   ENDDO
   IF (SCARC_DEBUG .GE. 6) write (9,*) '========================================&
  &======'
ENDDO
 
DO K = KBAR, 1, - 1
   DO J = JBAR, 1, - 1
      DO I = IBAR, 1, - 1
 
         IC = (K-1) * IBAR * JBAR + (J-1) * IBAR + I
 
         AUX = A (IC, 7) * Y (I, J, K+1) + A (IC, 6) * Y (I, J+1, K) + A (IC, 5) &
          * Y (I+1, J, K)
 
         Y (I, J, K) = Y (I, J, K) - AUX * OMEGA / A (IC, 4)
         IF (SCARC_DEBUG .GE. 6) write (9, '(a,i3,a,f8.2,a,f18.12,a,f18.12)') 'A&
        &(', IC, ')=', A (IC, 4), ': AUX=', AUX, ': Y=', Y (I, J, K)
      ENDDO
      IF (SCARC_DEBUG .GE. 6) write (9,*) '-------------'
   ENDDO
   IF (SCARC_DEBUG .GE. 6) write (9,*) '========================================&
  &======'
ENDDO
 
 
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'AFTER  PRECON_SSOR3D: Y'
   CALL SCARC_COMPARE_SINGLE (Y, 'SARC', 1, 'Y     ',0)
ENDIF
 
END SUBROUTINE SCARC_PRECON_SSOR3D
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Initialize GSADI preconditioner
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECON_GSADI_INIT
INTEGER :: IC, NC
 
NC = S%NCELLS_LOCAL
 
!!! copy diagonal
DO IC = 1, NC
   S%DIAG (IC) = S%AGLOB(IC, 3)
ENDDO
 
!!! copy subdiagonal
DO IC = 1, NC
   S%LOW (IC) = S%AGLOB(IC, 2)
ENDDO
 
!!! copy superdiagonal
DO IC = 1, NC
   S%UP (IC) = S%AGLOB(IC, 4)
ENDDO
 
!!! perform LU-factorization of matrix da
DO IC = 2, NC - 1
   S%LOW (IC) = S%LOW(IC) / S%DIAG(IC-1)
   S%DIAG (IC) = S%DIAG(IC) - S%LOW(IC) * S%UP(IC-1)
ENDDO
 
S%LOW (NC) = S%LOW(NC) / S%DIAG(NC-1)
S%DIAG (NC) = S%DIAG(NC) - S%LOW(NC) * S%UP(NC-1)
 
!!! replace diagonal values diag(i) by 1/diag(i)
!!! multiply U with 1/diag(i)
 
DO IC = 1, NC
   S%DIAG (IC) = 1.0_EB / S%DIAG(IC)
   S%UP (IC) = S%DIAG(IC) * S%UP(IC)
ENDDO
 
END SUBROUTINE SCARC_PRECON_GSADI_INIT
 
 
 
 
 
 
SUBROUTINE SCARC_RECEIVE (CODE, IMV)
 
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: I, J, K, LL, IC, IW, IWW, IERR, NM, NOM, CODE, IMV, III
INTEGER :: IFACE, IEDGE, IVRTX
TYPE (MESH_TYPE), POINTER :: M, M4
TYPE (OMESH_TYPE), POINTER :: M2, M3
 
INCLUDE 'mpif.h'
 
 
IF (SCARC_DEBUG .GE. 6) THEN
   WRITE (9,*) 'Entering SCARC_RECEIVE: CODE=', CODE, ': IMV=', IMV
   WRITE (9,*) 'NIC_FACE:'
   WRITE (9, '(4i6)') ((NIC_FACE(I, J), I=1, NMESHES), J=1, NMESHES)
   WRITE (9,*) 'TAGS_FACE:'
   WRITE (9, '(4i16)') ((TAGS_FACE(I, J), I=1, NMESHES), J=1, NMESHES)
   flush (9)
ENDIF
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CODE_IF: IF (CODE == 0) THEN
 
   INIT_LOOP: DO NM = 1, NMESHES
 
      IF (PROCESS(NM) /= MYID) CYCLE INIT_LOOP
      IF (SCARC_DEBUG .GE. 6) write (9,*) 'SCARC_RECEIVE: INIT: NM=', NM
      CALL flush (9)
 
!!!
!!! initialize matrix-vector-multiplication
!!!
      INIT_MV_LOOP: DO NOM = 1, NMESHES
 
         SNODE = PROCESS (NOM)
 
         IF (PROCESS(NM) == SNODE) CYCLE INIT_MV_LOOP
         IF (NIC(NM, NOM) == 0 .AND. NIC(NOM, NM) == 0) CYCLE INIT_MV_LOOP
 
         S3 => SCARC(NM)%OSCARC(NOM)
         TAG_MV = TAGS_MV (NM, NOM)
 
         IF (NIC(NM, NOM) > 0) THEN
            ALLOCATE (S3%RECVBUF_MV(Max(NIC(NM, NOM), NIC(NOM, NM))*2+1))
            S3%RECVBUF_MV = 0.0_EB
         ENDIF
 
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) 'SCARC_RECEIVE: INIT MV: TAG_MV=', TAG_MV
            WRITE (9,*) 'SCARC_RECEIVE: INIT MV: NOM=', NOM
            WRITE (9,*) 'SCARC_RECEIVE: ALLOCATING RECVBUF_MV in Length ', SIZE &
             (S3%RECVBUF_MV)
            CALL flush (9)
         ENDIF
 
      ENDDO INIT_MV_LOOP
 
 
!!!
!!! initialize extended face-communication
!!!
      INIT_FACE_LOOP: DO IFACE = 1, NFACE
 
         NOM = NBR_FACE (NM, IFACE)
         IF (NOM == 0) CYCLE INIT_FACE_LOOP
 
         IF (SCARC_DEBUG .GE. 6) write (9,*) 'SCARC_RECEIVE: INIT FACE: NOM=', &
          NOM
         CALL flush (9)
 
         SNODE = PROCESS (NOM)
         S3 => SCARC(NM)%OSCARC(NOM)
 
         IF (NIC_FACE(NM, NOM) > 0) THEN
            ALLOCATE (S3%RECVBUF_FACE(NIC_FACE(NM, NOM)+1))
            S3%RECVBUF_FACE = 0.0_EB
         ENDIF
 
      ENDDO INIT_FACE_LOOP
 
 
!!!
!!! initialize edge-communication (only 3D-case)
!!!
      IF ( .NOT. TWO_D) THEN
 
         INIT_EDGE_LOOP: DO IEDGE = 1, NEDGE
 
            NOM = NBR_EDGE (NM, IEDGE)
            IF (NOM == 0) CYCLE INIT_EDGE_LOOP
 
            IF (SCARC_DEBUG .GE. 6) write (9,*) 'SCARC_RECEIVE: INIT EDGE: NOM=',&
              NOM
            CALL flush (9)
 
            SNODE = PROCESS (NOM)
            S3 => SCARC(NM)%OSCARC(NOM)
 
            IF (NIC_EDGE(NM, NOM) > 0) THEN
               ALLOCATE (S3%RECVBUF_EDGE(NIC_EDGE(NM, NOM)+1))
               S3%RECVBUF_EDGE = 0.0_EB
            ENDIF
 
         ENDDO INIT_EDGE_LOOP
      ENDIF
 
 
!!!
!!! initialize vertex-communication
!!!
      INIT_VRTX_LOOP: DO IVRTX = 1, NVRTX
 
         NOM = NBR_VRTX (NM, IVRTX)
         IF (NOM == 0) CYCLE INIT_VRTX_LOOP
 
         IF (SCARC_DEBUG .GE. 6) write (9,*) 'SCARC_RECEIVE: INIT VRTX: NOM=', &
          NOM
         CALL flush (9)
 
         SNODE = PROCESS (NOM)
         S3 => SCARC(NM)%OSCARC(NOM)
 
         IF (NIC_VRTX(NM, NOM) > 0) THEN
            ALLOCATE (S3%RECVBUF_VRTX(NIC_VRTX(NM, NOM)+1))
            S3%RECVBUF_VRTX = 0.0_EB
         ENDIF
 
      ENDDO INIT_VRTX_LOOP
 
   ENDDO INIT_LOOP
 
 
   IF (SCARC_DEBUG .GE. 6) write (9,*) 'SCARC_RECEIVE: INIT LEAVING'
   CALL flush (9)
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Receive communication data for matrix-vector-multiplication (pure face communication)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE IF (CODE == 1) THEN
 
   MV_LOOP: DO NM = 1, NMESHES
 
      IF (PROCESS(NM) /= MYID) CYCLE MV_LOOP
 
      MV_NBR_LOOP: DO NOM = 1, NMESHES
 
         SNODE = PROCESS (NOM)
 
         IF (PROCESS(NM) == SNODE) CYCLE MV_NBR_LOOP
         IF (NIC(NM, NOM) == 0 .AND. NIC(NOM, NM) == 0) CYCLE MV_NBR_LOOP
 
         S3 => SCARC(NM)%OSCARC(NOM)
 
         TAG_MV = TAGS_MV (NM, NOM)
         N_REQ_MV = N_REQ_MV + 1
 
         CALL MPI_IRECV (S3%RECVBUF_MV(1), SIZE(S3%RECVBUF_MV), &
          MPI_DOUBLE_PRECISION, SNODE, TAG_MV, MPI_COMM_WORLD, REQ_MV(N_REQ_MV), &
          IERR)
 
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) 'NEXT_NEIGHBOR_RECV: ------- NOM  =', NOM
            WRITE (9,*) 'N_REQ_MV=', N_REQ_MV
            WRITE (9,*) 'REQ_MV(', N_REQ_MV, ')=', REQ_MV (N_REQ_MV)
            WRITE (9,*) 'TAG_MV=', TAG_MV
            WRITE (9,*) 'SNODE=', SNODE
            WRITE (9,*) 'NIC(', NOM, ',', NM, ')=', NIC (NOM, NM)
            WRITE (9,*) 'NIC(', NM, ',', NOM, ')=', NIC (NM, NOM)
            WRITE (9,*) NM, 'Reading ', SIZE (S3%RECVBUF_MV), ' data from ', NOM
         ENDIF
         CALL flush (9)
 
      ENDDO MV_NBR_LOOP
 
   ENDDO MV_LOOP
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Receive communication data for extended boundary update (including diagonal communication)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE IF (CODE == 2) THEN
 
 
   IF (SCARC_DEBUG .GE. 6) write (9,*) 'Halloooooooooooooooo'
 
   BDRY_LOOP: DO NM = 1, NMESHES
 
      IF (PROCESS(NM) /= MYID) CYCLE BDRY_LOOP
 
 
!!!---------------------------------------------------------------------------------------
!!! receive extended face-communication data
!!!---------------------------------------------------------------------------------------
      FACE_LOOP: DO IFACE = 1, NFACE
 
         NOM = NBR_FACE (NM, IFACE)
         IF (SCARC_DEBUG .GE. 6) write (9,*) 'NM=', NM
         IF (SCARC_DEBUG .GE. 6) write (9,*) 'NOM=', NOM
 
         IF (NOM == 0) CYCLE FACE_LOOP
 
         SNODE = PROCESS (NOM)
 
         S3 => SCARC(NM)%OSCARC(NOM)
         TAG_FACE = TAGS_FACE (NM, NOM)
 
         N_REQ_FACE = N_REQ_FACE + 1
 
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) '=================>   FACE_LOOP:'
            WRITE (9,*) 'CODE =', CODE
            WRITE (9,*) 'NM   =', NM
            WRITE (9,*) 'NOM  =', NOM
            WRITE (9,*) 'TAG_FACE  =', TAG_FACE
            WRITE (9,*) 'N_REQ_FACE=', N_REQ_FACE
            WRITE (9,*) 'REQ_FACE(', N_REQ_FACE, ')=', REQ_EDGE (N_REQ_FACE)
            WRITE (9,*) 'SNODE=', SNODE
            WRITE (9,*) 'NIC_FACE(', NOM, ',', NM, ')=', NIC_FACE (NOM, NM)
            WRITE (9,*) 'NIC_FACE(', NM, ',', NOM, ')=', NIC_FACE (NM, NOM)
            WRITE (9,*) NM, 'Reading ', SIZE (S3%RECVBUF_FACE), ' data from ', &
             NOM
            WRITE (9,*) 'RECVBUF_FACE:'
            WRITE (9, '(f25.16)') (S3%RECVBUF_FACE(I), I=1, &
             SIZE(S3%RECVBUF_FACE))
            CALL flush (9)
         ENDIF
 
         CALL MPI_IRECV (S3%RECVBUF_FACE(1), SIZE(S3%RECVBUF_FACE), &
          MPI_DOUBLE_PRECISION, SNODE, TAG_FACE, MPI_COMM_WORLD, &
          REQ_FACE(N_REQ_FACE), IERR)
 
      ENDDO FACE_LOOP
 
 
!!!---------------------------------------------------------------------------------------
!!! receive edge-communication data (only 3D-case)
!!!---------------------------------------------------------------------------------------
      IF ( .NOT. TWO_D) THEN
 
         EDGE_LOOP: DO IEDGE = 1, NEDGE
 
            NOM = NBR_EDGE (NM, IEDGE)
            IF (NOM == 0) CYCLE EDGE_LOOP
 
            SNODE = PROCESS (NOM)
            S3 => SCARC(NM)%OSCARC(NOM)
            TAG_EDGE = TAGS_EDGE (NM, NOM)
 
            N_REQ_EDGE = N_REQ_EDGE + 1
            CALL MPI_IRECV (S3%RECVBUF_EDGE(1), SIZE(S3%RECVBUF_EDGE), &
             MPI_DOUBLE_PRECISION, SNODE, TAG_EDGE, MPI_COMM_WORLD, &
             REQ_EDGE(N_REQ_EDGE), IERR)
 
            IF (SCARC_DEBUG .GE. 6) THEN
               WRITE (9,*) '=================>   EDGE_LOOP:'
               WRITE (9,*) 'NM   =', NM
               WRITE (9,*) 'NOM  =', NOM
               WRITE (9,*) 'N_REQ_EDGE=', N_REQ_EDGE
               WRITE (9,*) 'REQ_EDGE(', N_REQ_EDGE, ')=', REQ_EDGE (N_REQ_EDGE)
               WRITE (9,*) 'SNODE=', SNODE
               WRITE (9,*) 'NIC_EDGE(', NOM, ',', NM, ')=', NIC_EDGE (NOM, NM)
               WRITE (9,*) 'NIC_EDGE(', NM, ',', NOM, ')=', NIC_EDGE (NM, NOM)
               WRITE (9,*) NM, 'Reading ', SIZE (S3%RECVBUF_EDGE), ' data from ',&
                 NOM
               WRITE (9,*) 'RECVBUF_EDGE:'
               WRITE (9, '(f25.16)') (S3%RECVBUF_EDGE(I), I=1, &
                SIZE(S3%RECVBUF_EDGE))
               CALL flush (9)
            ENDIF
 
         ENDDO EDGE_LOOP
 
      ENDIF
 
 
!!!---------------------------------------------------------------------------------------
!!! receive vertex-communication data
!!!---------------------------------------------------------------------------------------
      VRTX_LOOP: DO IVRTX = 1, NVRTX
 
         NOM = NBR_VRTX (NM, IVRTX)
         IF (NOM == 0) CYCLE VRTX_LOOP
 
         SNODE = PROCESS (NOM)
         S3 => SCARC(NM)%OSCARC(NOM)
         TAG_VRTX = TAGS_VRTX (NM, NOM)
 
         N_REQ_VRTX = N_REQ_VRTX + 1
 
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) '=================>   VRTX_NEIGHBOR_IF'
            WRITE (9,*) 'NM   =', NM
            WRITE (9,*) 'NOM  =', NOM
            WRITE (9,*) 'N_REQ_VRTX=', N_REQ_VRTX
            WRITE (9,*) 'REQ_VRTX(', N_REQ_VRTX, ')=', REQ_VRTX (N_REQ_VRTX)
            WRITE (9,*) 'SNODE=', SNODE
            WRITE (9,*) 'NIC_VRTX(', NOM, ',', NM, ')=', NIC_VRTX (NOM, NM)
            WRITE (9,*) 'NIC_VRTX(', NM, ',', NOM, ')=', NIC_VRTX (NM, NOM)
            WRITE (9,*) NM, 'Reading ', SIZE (S3%RECVBUF_VRTX), ' data from ', &
             NOM
            WRITE (9,*) 'RECVBUF_VRTX:'
            WRITE (9, '(f25.16)') (S3%RECVBUF_VRTX(I), I=1, &
             SIZE(S3%RECVBUF_VRTX))
            CALL flush (9)
         ENDIF
         CALL MPI_IRECV (S3%RECVBUF_VRTX(1), SIZE(S3%RECVBUF_VRTX), &
          MPI_DOUBLE_PRECISION, SNODE, TAG_VRTX, MPI_COMM_WORLD, &
          REQ_VRTX(N_REQ_VRTX), IERR)
 
 
      ENDDO VRTX_LOOP
 
 
   ENDDO BDRY_LOOP
 
ENDIF CODE_IF
END SUBROUTINE SCARC_RECEIVE
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Scarc_Exchange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_EXCHANGE (CODE, IMV)
 
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX
INTEGER :: I, J, K, i4, J4, K4, LL, IC, IW, IWW, IERR, NM, NOM, CODE, IDIAG, &
 IMV, AAA, III, II, jj, KK
INTEGER :: IFACE, IEDGE, IVRTX
INTEGER :: iii0, jjj0, I1, J1
TYPE (MESH_TYPE), POINTER :: M, M4
TYPE (OMESH_TYPE), POINTER :: M2, M3
 
INCLUDE 'mpif.h'
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CODE_IF: IF (CODE == 0) THEN
 
   INIT_LOOP: DO NM = 1, NMESHES
 
      IF (PROCESS(NM) /= MYID) CYCLE INIT_LOOP
 
!!!----------------------------------------------------------------------------------------
!!! initialize matrix-vector-multiplication
!!!----------------------------------------------------------------------------------------
      INIT_MV_LOOP: DO NOM = 1, NMESHES
 
         IF (NIC(NM, NOM) == 0 .AND. NIC(NOM, NM) == 0) CYCLE INIT_MV_LOOP
 
         SNODE = PROCESS (NOM)
         RNODE = PROCESS (NM)
 
         S3 => SCARC(NM)%OSCARC(NOM)
 
         IF (NIC(NOM, NM) > 0 .AND. RNODE /= SNODE) THEN
            ALLOCATE (S3%SENDBUF_MV(Max(NIC(NOM, NM), NIC(NM, NOM))*2+1))
            S3%SENDBUF_MV = 0.0_EB
         ENDIF
 
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) '------->> 1:INIT_MV_LOOP'
            WRITE (9,*) 'IMV=', IMV
            WRITE (9,*) 'CODE=', CODE
            WRITE (9,*) 'NM=', NM
            WRITE (9,*) 'NOM=', NOM
            WRITE (9,*) 'SNODE=', SNODE
            WRITE (9,*) 'RNODE=', RNODE
            WRITE (9,*) 'NIC(', NM, ',', NOM, ')=', NIC (NM, NOM)
            WRITE (9,*) 'SCARC_EXCHANGE: ALLOCATING SENDBUF_MV in Length ', SIZE &
             (S3%SENDBUF_MV)
            CALL flush (9)
         ENDIF
 
      ENDDO INIT_MV_LOOP
 
 
!!!----------------------------------------------------------------------------------------
!!! initialize extended face-communication
!!!----------------------------------------------------------------------------------------
      INIT_FACE_LOOP: DO IFACE = 1, NFACE
 
         NOM = NBR_FACE (NM, IFACE)
         IF (NOM == 0) CYCLE INIT_FACE_LOOP
 
         SNODE = PROCESS (NOM)
         RNODE = PROCESS (NM)
 
         S3 => SCARC(NM)%OSCARC(NOM)
 
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) '------->> 1:INIT_FACE_LOOP'
            WRITE (9,*) 'IMV=', IMV
            WRITE (9,*) 'CODE=', CODE
            WRITE (9,*) 'NM=', NM
            WRITE (9,*) 'NOM=', NOM
            WRITE (9,*) 'SNODE=', SNODE
            WRITE (9,*) 'RNODE=', RNODE
            WRITE (9,*) 'NIC_FACE(', NM, ',', NOM, ')=', NIC_FACE (NM, NOM)
            CALL flush (9)
         ENDIF
 
         IF (NIC_FACE(NOM, NM) > 0 .AND. RNODE /= SNODE) THEN
            ALLOCATE (S3%SENDBUF_FACE(NIC_FACE(NM, NOM)+1))
            S3%SENDBUF_FACE = 0.0_EB
         ENDIF
 
      ENDDO INIT_FACE_LOOP
 
 
!!!----------------------------------------------------------------------------------------
!!! initialize edge-communication (only 3D)
!!!----------------------------------------------------------------------------------------
      IF ( .NOT. TWO_D) THEN
         INIT_EDGE_LOOP: DO IEDGE = 1, NEDGE
 
            NOM = NBR_EDGE (NM, IEDGE)
            IF (NOM == 0) CYCLE INIT_EDGE_LOOP
 
            SNODE = PROCESS (NOM)
            RNODE = PROCESS (NM)
 
            S3 => SCARC(NM)%OSCARC(NOM)
 
            IF (SCARC_DEBUG .GE. 6) THEN
               WRITE (9,*) '------->> 1:INIT_EDGE_LOOP'
               WRITE (9,*) 'IMV=', IMV
               WRITE (9,*) 'CODE=', CODE
               WRITE (9,*) 'NM=', NM
               WRITE (9,*) 'NOM=', NOM
               WRITE (9,*) 'SNODE=', SNODE
               WRITE (9,*) 'RNODE=', RNODE
               WRITE (9,*) 'NIC_EDGE(', NM, ',', NOM, ')=', NIC_EDGE (NM, NOM)
               CALL flush (9)
            ENDIF
 
            IF (NIC_EDGE(NOM, NM) > 0 .AND. RNODE /= SNODE) THEN
               ALLOCATE (S3%SENDBUF_EDGE(NIC_EDGE(NM, NOM)+1))
               S3%SENDBUF_EDGE = 0.0_EB
            ENDIF
 
         ENDDO INIT_EDGE_LOOP
      ENDIF
 
 
 
!!!----------------------------------------------------------------------------------------
!!! initialize vertex-communication
!!!----------------------------------------------------------------------------------------
      INIT_VRTX_LOOP: DO IVRTX = 1, NVRTX
 
         NOM = NBR_VRTX (NM, IVRTX)
         IF (NOM == 0) CYCLE INIT_VRTX_LOOP
 
         SNODE = PROCESS (NOM)
         RNODE = PROCESS (NM)
 
         S3 => SCARC(NM)%OSCARC(NOM)
 
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) '------->> 1:INIT_VRTX_LOOP'
            WRITE (9,*) 'IMV=', IMV
            WRITE (9,*) 'CODE=', CODE
            WRITE (9,*) 'NM=', NM
            WRITE (9,*) 'NOM=', NOM
            WRITE (9,*) 'SNODE=', SNODE
            WRITE (9,*) 'RNODE=', RNODE
            WRITE (9,*) 'NIC_VRTX(', NM, ',', NOM, ')=', NIC_VRTX (NM, NOM)
            CALL flush (9)
         ENDIF
 
         IF (NIC_VRTX(NOM, NM) > 0 .AND. RNODE /= SNODE) THEN
            ALLOCATE (S3%SENDBUF_VRTX(NIC_VRTX(NM, NOM)+1))
            S3%SENDBUF_VRTX = 0.0_EB
         ENDIF
 
      ENDDO INIT_VRTX_LOOP
 
   ENDDO INIT_LOOP
 
 
   IF (SCARC_DEBUG .GE. 6) write (9,*) 'LEAVING SCARC_EXCHANGE: INIT'
   CALL flush (9)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Exchange communication data for matrix-vector-multiplication (pure face communication)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE IF (CODE == 1) THEN
 
!!!----------------------------------------------------------------------------------------
!!! Pack exchange data for matrix-vector-multiplication
!!!----------------------------------------------------------------------------------------
   MV_PACK_LOOP: DO NM = 1, NMESHES
 
      IF (PROCESS(NM) /= MYID) CYCLE MV_PACK_LOOP
 
      MV_PACK_NBR_LOOP: DO NOM = 1, NMESHES
 
         SNODE = PROCESS (NOM)
         RNODE = PROCESS (NM)
 
         IF (NIC(NM, NOM) == 0 .AND. NIC(NOM, NM) == 0) CYCLE MV_PACK_NBR_LOOP
 
         S => SCARC (NM)
         S3 => SCARC(NM)%OSCARC(NOM)
         M3 => MESHES(NM)%OMESH(NOM)
         M4 => MESHES (NOM)
 
         TAG_MV = TAGS_MV (NM, NOM)
 
         IMIN = I_MIN (NOM, NM)
         IMAX = I_MAX (NOM, NM)
         JMIN = J_MIN (NOM, NM)
         JMAX = J_MAX (NOM, NM)
         KMIN = K_MIN (NOM, NM)
         KMAX = K_MAX (NOM, NM)
 
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) '------->> 1:IN SCARC_EXCHANGE: MV_PACK_NBR_LOOP'
            WRITE (9,*) 'IMV=', IMV
            WRITE (9,*) 'CODE=', CODE
            WRITE (9,*) 'NM=', NM
            WRITE (9,*) 'NOM=', NOM
            WRITE (9,*) 'IBAR=', IBAR
            WRITE (9,*) 'JBAR=', JBAR
            WRITE (9,*) 'KBAR=', KBAR
            WRITE (9,*) 'SNODE=', SNODE
            WRITE (9,*) 'RNODE=', RNODE
            WRITE (9,*) 'IMIN=', IMIN
            WRITE (9,*) 'IMAX=', IMAX
            WRITE (9,*) 'JMIN=', JMIN
            WRITE (9,*) 'JMAX=', JMAX
            WRITE (9,*) 'KMIN=', KMIN
            WRITE (9,*) 'KMAX=', KMAX
            WRITE (9,*) 'S  zeigt auf SCARC(', NM, ')'
            WRITE (9,*) 'S3 zeigt auf SCARC(', NM, ')%OSCARC(', NOM, ')'
            WRITE (9,*) 'M3 zeigt auf MESHES(', NM, ')%OMESH(', NOM, ')'
            WRITE (9,*) 'M4 zeigt auf MESHES(', NOM, ')'
            WRITE (9,*) 'M4%NEWC=', M4%NEWC
            CALL flush (9)
         ENDIF
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) 'X1:NIC_FACE:'
            WRITE (9, '(4i16)') ((NIC_FACE(I1, J1), I1=1, NMESHES), J1=1, &
             NMESHES)
            WRITE (9,*) 'X1:TAGS_FACE:'
            WRITE (9, '(4i16)') ((TAGS_FACE(I1, J1), I1=1, NMESHES), J1=1, &
             NMESHES)
         ENDIF
 
 
         IF (RNODE /= SNODE) THEN
 
            LL = 0
            IWW = 0
 
            MV_SENDBUF_LOOP: DO IW = 1, M4%NEWC
               IF (M3%IJKW(9, IW) /= NM) CYCLE MV_SENDBUF_LOOP
               DO K = M3%IJKW(12, IW), M3%IJKW(15, IW)
                  DO J = M3%IJKW(11, IW), M3%IJKW(14, IW)
                     DO I = M3%IJKW(10, IW), M3%IJKW(13, IW)
                        IWW = IWW + 1
                        S3%SENDBUF_MV (LL+1) = REAL (IW, EB)
                        IF (IMV .GE. 5) THEN
                           S3%SENDBUF_MV (LL+2) = S%D(I, J, K)
         !write(9,*) 'Sending D:',S%D(I,J,K)
                        ELSE IF (IMV .EQ. 3) THEN
                           S3%SENDBUF_MV (LL+2) = S%X(I, J, K)
         !write(9,*) 'Sending X:',S%X(I,J,K)
                        ELSE
                           WRITE (*,*) 'Wrong type for matrix-vector-multiplicat&
                          &ion'
                           STOP
                        ENDIF
!!!!! ACHTUNG: FALLS KEIN STANDARD-5-PUNKTE-STERN, muss hier AG(.)*X(.,.,.) gesendet werden!
                        IF (SCARC_DEBUG .GE. 6) write (9, '(a,i3,a,i3,a,i3,a,i3,&
                       &a,f22.16,a,f22.16)') 'MATVEC_MUL: (', I, ',', J, ',', K, &
                         '):   SENDBUF_MV(', LL + 1, ')=', S3%SENDBUF_MV(LL+1), &
                         ',', S3%SENDBUF_MV(LL+2)
                        CALL flush (9)
                        LL = LL + 2
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO MV_SENDBUF_LOOP
 
            IF (SCARC_DEBUG .GE. 6) THEN
               WRITE (9,*) 'X2:NIC_FACE:'
               WRITE (9, '(4i16)') ((NIC_FACE(I1, J1), I1=1, NMESHES), J1=1, &
                NMESHES)
               WRITE (9,*) 'X2:TAGS_FACE:'
               WRITE (9, '(4i16)') ((TAGS_FACE(I1, J1), I1=1, NMESHES), J1=1, &
                NMESHES)
            ENDIF
 
            S3%SENDBUF_MV (IWW*2+1) = - 999.0_EB
            N_REQ_MV = N_REQ_MV + 1
 
            CALL MPI_ISEND (S3%SENDBUF_MV(1), IWW*2+1, MPI_DOUBLE_PRECISION, &
             SNODE, TAG_MV, MPI_COMM_WORLD, REQ_MV(N_REQ_MV), IERR)
 
            IF (SCARC_DEBUG .GE. 6) THEN
               WRITE (9,*) 'SENDBUF_MV(', IWW + 2 + 1, ')=', &
                S3%SENDBUF_MV(IWW*2+1)
               WRITE (9,*) 'N_REQ_MV=', N_REQ_MV
               WRITE (9,*) 'REQ_MV(', N_REQ_MV, ')=', REQ_MV (N_REQ_MV)
               WRITE (9,*) 'S3 zeigt auf SCARC(', NM, ')%OSCARC(', NOM, ')'
               CALL flush (9)
            ENDIF
 
 
         ELSE
 
            IF (SCARC_DEBUG .GE. 2) THEN
               WRITE (9,*) 'ELSEHALLO: FILLING SCARC(', NOM, ')%OSCARC(', NM, ')&
              &'
               WRITE (9,*) '------------------'
               CALL flush (9)
            ENDIF
 
            S2 => SCARC(NOM)%OSCARC(NM)
            S2%Y_MV (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX) = S%Y(IMIN:IMAX, &
             JMIN:JMAX, KMIN:KMAX)
 
            IF (SCARC_DEBUG .GE. 2) THEN
               WRITE (9,*) '------------------'
               WRITE (9,*) 'S2%Y(', IMIN, ',', JMIN, ',', KMIN, ')=', &
                S2%Y_MV(IMIN, JMIN, KMIN)
               WRITE (9,*) 'S2%Y(', IMAX, ',', JMAX, ',', KMAX, ')=', &
                S2%Y_MV(IMAX, JMAX, KMAX)
               CALL flush (9)
            ENDIF
 
         ENDIF
 
      ENDDO MV_PACK_NBR_LOOP
   ENDDO MV_PACK_LOOP
 
   IF (SCARC_DEBUG .GE. 6) THEN
      WRITE (9,*) 'X3:NIC_FACE:'
      WRITE (9, '(4i16)') ((NIC_FACE(I1, J1), I1=1, NMESHES), J1=1, NMESHES)
      WRITE (9,*) 'X3:TAGS_FACE:'
      WRITE (9, '(4i16)') ((TAGS_FACE(I1, J1), I1=1, NMESHES), J1=1, NMESHES)
   ENDIF
 
 
!!!----------------------------------------------------------------------------------------
!!! Wait for complete matrix-vector-multiplication exchange
!!!----------------------------------------------------------------------------------------
   CALL MPI_WAITALL (N_REQ_MV, REQ_MV(1:N_REQ_MV), STATUSES_MV, IERR)
 
 
!!!----------------------------------------------------------------------------------------
!!! Unpack exchange data for matrix-vector-multiplication
!!!----------------------------------------------------------------------------------------
   LL = 0
   MV_UNPACK_LOOP: DO NOM = 1, NMESHES
 
      IF (PROCESS(NOM) /= MYID) CYCLE MV_UNPACK_LOOP
 
      MV_UNPACK_NBR_LOOP: DO NM = 1, NMESHES
 
         SNODE = PROCESS (NOM)
         RNODE = PROCESS (NM)
 
         IF (NIC(NM, NOM) == 0 .AND. NIC(NOM, NM) == 0) CYCLE MV_UNPACK_NBR_LOOP
 
         S2 => SCARC(NOM)%OSCARC(NM)
         S4 => SCARC (NOM)
 
         M4 => MESHES (NOM)
 
         IF (SCARC_DEBUG .GE. 2) THEN
            WRITE (9,*) '--------- unpacking ...'
            WRITE (9,*) 'NM=', NM
            WRITE (9,*) 'NOM=', NOM
            WRITE (9,*) 'CODE=', CODE
            WRITE (9,*) 'NIC(', NOM, ',', NM, ')=', NIC (NOM, NM)
            WRITE (9,*) 'RNODE=', RNODE
            WRITE (9,*) 'SNODE=', SNODE
            WRITE (9,*) 'S2 zeigt auf SCARC(', NOM, ')%OSCARC(', NM, ')'
            WRITE (9,*) 'S4 zeigt auf SCARC(', NOM, ')'
            WRITE (9,*) 'RECVBUF_MV: AFTER WAITALLL'
            WRITE (9, '(2f22.16)') (S2%RECVBUF_MV(I), I=1, SIZE(S2%RECVBUF_MV))
            CALL flush (9)
         ENDIF
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) 'X4:NIC_FACE:'
            WRITE (9, '(4i16)') ((NIC_FACE(I1, J1), I1=1, NMESHES), J1=1, &
             NMESHES)
            WRITE (9,*) 'X4:TAGS_FACE:'
            WRITE (9, '(4i16)') ((TAGS_FACE(I1, J1), I1=1, NMESHES), J1=1, &
             NMESHES)
         ENDIF
 
 
         MV_UNPACK_IF: IF (NIC(NOM, NM) > 0 .AND. RNODE /= SNODE) THEN
            LL = 0
 
            MV_RECVBUF_LOOP: DO
               IW = Nint (S2%RECVBUF_MV(LL+1))
               IF (IW ==-999) EXIT MV_RECVBUF_LOOP
               i4 = M4%IJKW (6, IW)
               J4 = M4%IJKW (7, IW)
               K4 = M4%IJKW (8, IW)
 
 
 
               DO I = M4%IJKW(10, IW), M4%IJKW(13, IW)
                  DO J = M4%IJKW(11, IW), M4%IJKW(14, IW)
                     DO K = M4%IJKW(12, IW), M4%IJKW(15, IW)
 
                        IF (SCARC_DEBUG .GE. 6) THEN
                           WRITE (9,*) '????????????????????????????????????????&
                          &????????'
                           WRITE (9,*) 'I4=', i4, ',J4=', J4, ',K4=', K4
                           WRITE (9,*) 'I =', I, ',J =', J, ',K =', K
                           WRITE (9,*) '????????????????????????????????????????&
                          &????????'
                        ENDIF
                        IF (IMV == 1 .OR. IMV == 5) THEN
                           IF (SCARC_DEBUG .GE. 6) write (9,*) '1:IMV=', IMV
                           S2%Y_MV (I, J, K) = S2%RECVBUF_MV(LL+2)
                           IF (SCARC_DEBUG .GE. 6) THEN
                              WRITE (9,*) 'X41A:NIC_FACE:'
                              WRITE (9, '(4i16)') ((NIC_FACE(I1, J1), I1=1, &
                               NMESHES), J1=1, NMESHES)
                           ENDIF
 
                        ELSE IF (IMV == 2 .OR. IMV == 6) THEN
                           IF (SCARC_DEBUG .GE. 6) write (9,*) '2:IMV=', IMV
                           S2%G_MV (I, J, K) = S2%RECVBUF_MV(LL+2)
                           IF (SCARC_DEBUG .GE. 6) THEN
                              WRITE (9,*) 'X41B:NIC_FACE:'
                              WRITE (9, '(4i16)') ((NIC_FACE(I1, J1), I1=1, &
                               NMESHES), J1=1, NMESHES)
                           ENDIF
 
                        ELSE IF (IMV == 3) THEN
                           IF (SCARC_DEBUG .GE. 6) write (9,*) '3:IMV=', IMV, ':&
                            RECVBUF_MV(', LL + 2, ')=', S2%RECVBUF_MV(LL+2)
                           S2%R_MV (I, J, K) = S2%RECVBUF_MV(LL+2)
                           IF (SCARC_DEBUG .GE. 6) THEN
                              WRITE (9,*) 'X41C:NIC_FACE:'
                              WRITE (9, '(4i16)') ((NIC_FACE(I1, J1), I1=1, &
                               NMESHES), J1=1, NMESHES)
                           ENDIF
 
                        ENDIF
 
                        IF (SCARC_DEBUG .GE. 6) THEN
         !write(9,*) 'S4%PERIOD(',IW,')=',S4%PERIOD(IW)
                           WRITE (9,*) '========================================&
                          &'
                           IF (IMV == 3) THEN
                              WRITE (9, '(a,i3,a,i2,a,i2,a,i2,a,f22.16)') 'IW=', &
                               IW, ':S2%R_MV(', I, ',', J, ',', K, ')=', &
                               S2%R_MV(I, J, K)
                              WRITE (9, '(a,i3,a,i2,a,i2,a,i2,a,f22.16)') 'IW=', &
                               IW, ':BEFORE :S4%R(', i4, ',', J4, ',', K4, ')=', &
                               S4%R(i4, J4, K4)
                           ELSE
                              WRITE (9, '(a,i3,a,i2,a,i2,a,i2,a,f22.16)') 'IW=', &
                               IW, ':S2%Y_MV(', I, ',', J, ',', K, ')=', &
                               S2%Y_MV(I, J, K)
                              WRITE (9, '(a,i3,a,i2,a,i2,a,i2,a,f22.16)') 'IW=', &
                               IW, ':BEFORE :S4%Y(', i4, ',', J4, ',', K4, ')=', &
                               S4%Y(i4, J4, K4)
                           ENDIF
                           CALL flush (9)
                        ENDIF
 
       if (S4%PERIOD(IW)/=100) then
                        IF (IMV == 1 .OR. IMV == 5) THEN
                           IF (SCARC_DEBUG .GE. 6) THEN
                              WRITE (9,*) '1A:IMV=', IMV
                              WRITE (9,*) 'J =', J
                              WRITE (9,*) 'J4=', J4
                              WRITE (9,*) 'SIZE(S4%Y)=', SIZE (S4%Y)
                              WRITE (9,*) 'SIZE(S2%Y_MV)=', SIZE (S2%Y_MV)
                           ENDIF
 
                           S4%Y (i4, J4, K4) = S4%Y(i4, J4, K4) + S4%ASUB * &
                            S2%Y_MV(I, J, K)
 
                           IF (SCARC_DEBUG .GE. 6) THEN
                              WRITE (9,*) '1A:IMV=', IMV, ' Y_MV=', S2%Y_MV(I, &
                               J, K), ': ASUB=', S4%ASUB
                              WRITE (9,*) 'Y:', S4%Y(i4, J4, K4)
                           ENDIF
                           IF (SCARC_DEBUG .GE. 6) THEN
                              WRITE (9, '(a,i3,a,i2,a,i2,a,i2,a,f22.16)') 'IW=', &
                               IW, ':AFTER  :S4%Y(', i4, ',', J4, ',', K4, ')=', &
                               S4%Y(i4, J4, K4)
                           ENDIF
                        ELSE IF (IMV == 2 .OR. IMV == 6) THEN
                           IF (SCARC_DEBUG .GE. 6) write (9,*) '2A:IMV=', IMV
                           S4%G (i4, J4, K4) = S4%G(i4, J4, K4) + S4%ASUB * &
                            S2%G_MV(I, J, K)
                        ELSE IF (IMV == 3) THEN
                           S4%R (i4, J4, K4) = S4%R(i4, J4, K4) + S4%ASUB * &
                            S2%R_MV(I, J, K)
 
                           IF (SCARC_DEBUG .GE. 6) THEN
                              WRITE (9,*) '3A:IMV=', IMV, ' R_MV=', S2%R_MV(I, &
                               J, K), ': ASUB=', S4%ASUB
                              WRITE (9,*) 'R:', S4%R(i4, J4, K4)
                           ENDIF
                           IF (SCARC_DEBUG .GE. 6) THEN
                              WRITE (9, '(a,i3,a,i2,a,i2,a,i2,a,f22.16)') 'IW=', &
                               IW, ':AFTER  :S4%R(', i4, ',', J4, ',', K4, ')=', &
                               S4%R(i4, J4, K4)
                           ENDIF
 
                        ENDIF
 
       endif
                        IF (SCARC_DEBUG .GE. 6) THEN
                           WRITE (9,*) 'X42:NIC_FACE:'
                           WRITE (9, '(4i16)') ((NIC_FACE(I1, J1), I1=1, &
                            NMESHES), J1=1, NMESHES)
                           WRITE (9,*) 'X42:TAGS_FACE:'
                           WRITE (9, '(4i16)') ((TAGS_FACE(I1, J1), I1=1, &
                            NMESHES), J1=1, NMESHES)
                        ENDIF
 
                        LL = LL + 2
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO MV_RECVBUF_LOOP
         ENDIF MV_UNPACK_IF
!call show_vector3(S4%R,'EX:R  ')
 
      ENDDO MV_UNPACK_NBR_LOOP
   ENDDO MV_UNPACK_LOOP
 
   IF (SCARC_DEBUG .GE. 6) THEN
      WRITE (9,*) 'X5:NIC_FACE:'
      WRITE (9, '(4i16)') ((NIC_FACE(I1, J1), I1=1, NMESHES), J1=1, NMESHES)
      WRITE (9,*) 'X5:TAGS_FACE:'
      WRITE (9, '(4i16)') ((TAGS_FACE(I1, J1), I1=1, NMESHES), J1=1, NMESHES)
   ENDIF
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Exchange communication data for extended boundary update (including diagonal communication)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE IF (CODE == 2) THEN
 
 
!!!----------------------------------------------------------------------------------------
!!! pack extended face-communication data
!!!----------------------------------------------------------------------------------------
   FACE_PACK_LOOP: DO NM = 1, NMESHES
 
      IF (PROCESS(NM) /= MYID) CYCLE FACE_PACK_LOOP
 
      FACE_PACK_NBR_LOOP: DO IFACE = 1, NFACE
 
         NOM = NBR_FACE (NM, IFACE)
         IF (NOM == 0) CYCLE FACE_PACK_NBR_LOOP
 
         SNODE = PROCESS (NOM)
         RNODE = PROCESS (NM)
 
         S => SCARC (NM)
 
         S3 => SCARC(NM)%OSCARC(NOM)
         S4 => SCARC (NOM)
 
         M3 => MESHES(NM)%OMESH(NOM)
         M4 => MESHES (NOM)
 
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) '------->> 1:IN SCARC_EXCHANGE: FACE_PACK_NBR_LOOP'
            WRITE (9,*) 'IMV=', IMV
            WRITE (9,*) 'CODE=', CODE
            WRITE (9,*) 'NM=', NM
            WRITE (9,*) 'NOM=', NOM
            WRITE (9,*) 'SNODE=', SNODE
            WRITE (9,*) 'RNODE=', RNODE
            CALL flush (9)
         ENDIF
 
         TAG_FACE = TAGS_FACE (NM, NOM)
 
         IF (RNODE /= SNODE) THEN
 
            LL = 0
            IWW = 0
 
            DO K = IJKW_FACE (IFACE, 3), IJKW_FACE (IFACE, 6)
               DO J = IJKW_FACE (IFACE, 2), IJKW_FACE (IFACE, 5)
                  DO I = IJKW_FACE (IFACE, 1), IJKW_FACE (IFACE, 4)
                     IWW = IWW + 1
                     S3%SENDBUF_FACE (LL+1) = S%Z(I, J, K)
                     IF (SCARC_DEBUG .GE. 6) write (9, '(a,i3,a,i3,a,i3,a,i3,a,f&
                    &22.16)') 'BDRY_UPDATE_FACE: (', I, ',', J, ',', K, '):   SE&
                    &NDBUF_FACE(', LL + 1, ')=', S3%SENDBUF_FACE(LL+1)
                     LL = LL + 1
                  ENDDO
               ENDDO
            ENDDO
 
            S3%SENDBUF_FACE (IWW+1) = - 999.0_EB
            N_REQ_FACE = N_REQ_FACE + 1
 
            IF (SCARC_DEBUG .GE. 6) THEN
               WRITE (9,*) 'SENDBUF_FACE(', IWW + 1, ')=', &
                S3%SENDBUF_FACE(IWW+1)
               WRITE (9,*) 'N_REQ_FACE=', N_REQ_FACE
               CALL flush (9)
            ENDIF
 
            CALL MPI_ISEND (S3%SENDBUF_FACE(1), IWW+1, MPI_DOUBLE_PRECISION, &
             SNODE, TAG_FACE, MPI_COMM_WORLD, REQ_FACE(N_REQ_FACE), IERR)
 
         ENDIF
 
      ENDDO FACE_PACK_NBR_LOOP
   ENDDO FACE_PACK_LOOP
 
 
!!!----------------------------------------------------------------------------------------
!!! Waiting for complete face data exchange
!!!----------------------------------------------------------------------------------------
   CALL MPI_WAITALL (N_REQ_FACE, REQ_FACE(1:N_REQ_FACE), STATUSES_FACE, IERR)
 
 
   IF (SCARC_DEBUG .GE. 6) THEN
      WRITE (9,*) 'REQ_FACE(', N_REQ_FACE, ')=', REQ_FACE (N_REQ_FACE)
      flush (9)
   ENDIF
 
 
!!!----------------------------------------------------------------------------------------
!!! pack edge-communication data (only 3D-case)
!!!----------------------------------------------------------------------------------------
   IF ( .NOT. TWO_D) THEN
      EDGE_PACK_LOOP: DO NM = 1, NMESHES
 
         IF (PROCESS(NM) /= MYID) CYCLE EDGE_PACK_LOOP
 
         EDGE_PACK_NBR_LOOP: DO IEDGE = 1, NEDGE
 
            NOM = NBR_EDGE (NM, IEDGE)
            IF (NOM == 0) CYCLE EDGE_PACK_NBR_LOOP
 
            SNODE = PROCESS (NOM)
            RNODE = PROCESS (NM)
 
            S => SCARC (NM)
            S3 => SCARC(NM)%OSCARC(NOM)
            S4 => SCARC (NOM)
 
            IF (SCARC_DEBUG .GE. 6) THEN
               WRITE (9,*) '------->> 1:IN SCARC_EXCHANGE: EDGE_PACK_LOOP'
               WRITE (9,*) 'IMV=', IMV
               WRITE (9,*) 'CODE=', CODE
               WRITE (9,*) 'NM=', NM
               WRITE (9,*) 'NOM=', NOM
               WRITE (9,*) 'SNODE=', SNODE
               WRITE (9,*) 'RNODE=', RNODE
               CALL flush (9)
            ENDIF
 
            M3 => MESHES(NM)%OMESH(NOM)
            M4 => MESHES (NOM)
 
            TAG_EDGE = TAGS_EDGE (NM, NOM)
 
            IF (RNODE /= SNODE) THEN
 
               LL = 0
               IWW = 0
 
               DO K = IJKW_EDGE (IEDGE, 3), IJKW_EDGE (IEDGE, 6)
                  DO J = IJKW_EDGE (IEDGE, 2), IJKW_EDGE (IEDGE, 5)
                     DO I = IJKW_EDGE (IEDGE, 1), IJKW_EDGE (IEDGE, 4)
                        IWW = IWW + 1
                        S3%SENDBUF_EDGE (LL+1) = S%Z(I, J, K)
                        IF (SCARC_DEBUG .GE. 6) write (9, '(a,i3,a,i3,a,i3,a,i3,&
                       &a,f22.16)') 'BDRY_UPDATE_EDGE: (', I, ',', J, ',', K, ')&
                       &:   SENDBUF_EDGE(', LL + 1, ')=', S3%SENDBUF_EDGE(LL+1)
                        LL = LL + 1
                     ENDDO
                  ENDDO
               ENDDO
 
               S3%SENDBUF_EDGE (IWW+1) = - 999.0_EB
               N_REQ_EDGE = N_REQ_EDGE + 1
 
               IF (SCARC_DEBUG .GE. 6) THEN
                  WRITE (9,*) 'SENDBUF_EDGE(', IWW + 1, ')=', &
                   S3%SENDBUF_EDGE(IWW+1)
                  WRITE (9,*) 'N_REQ_EDGE=', N_REQ_EDGE
                  CALL flush (9)
               ENDIF
 
               CALL MPI_ISEND (S3%SENDBUF_EDGE(1), IWW+1, MPI_DOUBLE_PRECISION, &
                SNODE, TAG_EDGE, MPI_COMM_WORLD, REQ_EDGE(N_REQ_EDGE), IERR)
 
            ENDIF
 
         ENDDO EDGE_PACK_NBR_LOOP
      ENDDO EDGE_PACK_LOOP
 
 
!!! ----------------------------------------------------------------------------------------
!!! Waiting for complete edge data exchange
!!!----------------------------------------------------------------------------------------
      CALL MPI_WAITALL (N_REQ_EDGE, REQ_EDGE(1:N_REQ_EDGE), STATUSES_EDGE, IERR)
 
      IF (SCARC_DEBUG .GE. 6) THEN
         WRITE (9,*) 'REQ_EDGE(', N_REQ_EDGE, ')=', REQ_EDGE (N_REQ_EDGE)
         flush (9)
      ENDIF
 
   ENDIF
 
 
!!!----------------------------------------------------------------------------------------
!!! pack vertex-communication data
!!!----------------------------------------------------------------------------------------
   VRTX_PACK_LOOP: DO NM = 1, NMESHES
 
      IF (PROCESS(NM) /= MYID) CYCLE VRTX_PACK_LOOP
 
      VRTX_PACK_NBR_LOOP: DO IVRTX = 1, NVRTX
 
         NOM = NBR_VRTX (NM, IVRTX)
         IF (NOM == 0) CYCLE VRTX_PACK_NBR_LOOP
 
         SNODE = PROCESS (NOM)
         RNODE = PROCESS (NM)
 
         S => SCARC (NM)
 
         S3 => SCARC(NM)%OSCARC(NOM)
         S4 => SCARC (NOM)
 
         M3 => MESHES(NM)%OMESH(NOM)
         M4 => MESHES (NOM)
 
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) '------->> 1:IN SCARC_EXCHANGE: VRTX_PACK_LOOP'
            WRITE (9,*) 'IMV=', IMV
            WRITE (9,*) 'CODE=', CODE
            WRITE (9,*) 'NM=', NM
            WRITE (9,*) 'NOM=', NOM
            WRITE (9,*) 'SNODE=', SNODE
            WRITE (9,*) 'RNODE=', RNODE
         ENDIF
 
         TAG_VRTX = TAGS_VRTX (NM, NOM)
 
         IF (RNODE /= SNODE) THEN
 
            LL = 0
            IWW = 0
 
            DO K = IJKW_VRTX (IVRTX, 3), IJKW_VRTX (IVRTX, 6)
               DO J = IJKW_VRTX (IVRTX, 2), IJKW_VRTX (IVRTX, 5)
                  DO I = IJKW_VRTX (IVRTX, 1), IJKW_VRTX (IVRTX, 4)
                     IWW = IWW + 1
                     S3%SENDBUF_VRTX (LL+1) = S%Z(I, J, K)
                     IF (SCARC_DEBUG .GE. 6) write (9, '(a,i3,a,i3,a,i3,a,i3,a,f&
                    &22.16)') 'BDRY_UPDATE_VRTX: (', I, ',', J, ',', K, '):   SE&
                    &NDBUF_VRTX(', LL + 1, ')=', S3%SENDBUF_VRTX(LL+1)
                     LL = LL + 1
                  ENDDO
               ENDDO
            ENDDO
 
            S3%SENDBUF_VRTX (IWW+1) = - 999.0_EB
            N_REQ_VRTX = N_REQ_VRTX + 1
 
            IF (SCARC_DEBUG .GE. 6) THEN
               WRITE (9,*) 'SENDBUF_VRTX(', IWW + 1, ')=', &
                S3%SENDBUF_VRTX(IWW+1)
               WRITE (9,*) 'N_REQ_VRTX=', N_REQ_VRTX
            ENDIF
 
            CALL MPI_ISEND (S3%SENDBUF_VRTX(1), IWW+1, MPI_DOUBLE_PRECISION, &
             SNODE, TAG_VRTX, MPI_COMM_WORLD, REQ_VRTX(N_REQ_VRTX), IERR)
 
         ENDIF
 
      ENDDO VRTX_PACK_NBR_LOOP
   ENDDO VRTX_PACK_LOOP
 
 
!!!----------------------------------------------------------------------------------------
!!! Waiting for complete vertex data exchange
!!!----------------------------------------------------------------------------------------
   CALL MPI_WAITALL (N_REQ_VRTX, REQ_VRTX(1:N_REQ_VRTX), STATUSES_VRTX, IERR)
 
   IF (SCARC_DEBUG .GE. 6) THEN
      WRITE (9,*) 'REQ_VRTX(', N_REQ_VRTX, ')=', REQ_VRTX (N_REQ_VRTX)
      flush (9)
   ENDIF
 
   UNPACK_LOOP: DO NM = 1, NMESHES
 
 
      IF (PROCESS(NM) /= MYID) CYCLE UNPACK_LOOP
 
!!!----------------------------------------------------------------------------------------
!!! Unpack extended face-communication data
!!!----------------------------------------------------------------------------------------
      FACE_UNPACK_LOOP: DO IFACE = 1, NFACE
 
         NOM = NBR_FACE (NM, IFACE)
         IF (NOM == 0) CYCLE FACE_UNPACK_LOOP
 
         S => SCARC (NM)
         S2 => SCARC(NM)%OSCARC(NOM)
 
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) 'IFACE=', IFACE
            WRITE (9,*) 'NM   =', NM
            WRITE (9,*) 'NOM  =', NOM
            WRITE (9,*) 'S2 zeigt auf SCARC(', NM, ')%OSCARC(', NOM, ')'
            WRITE (9,*) 'S4 zeigt auf SCARC(', NM, ')'
            WRITE (9,*) 'IJKW_FACE(:'
            WRITE (9, '(12i4)') (IJKW_FACE(IFACE, III), III=1, 12)
            CALL flush (9)
         ENDIF
 
         LL = 0
         DO K = IJKW_FACE (IFACE, 9), IJKW_FACE (IFACE, 12)
            DO J = IJKW_FACE (IFACE, 8), IJKW_FACE (IFACE, 11)
               DO I = IJKW_FACE (IFACE, 7), IJKW_FACE (IFACE, 10)
 
                  IF (SCARC_DEBUG .GE. 6) write (9, '(a,i2,a,i2,a,i2,a,f22.16)') &
                   'BEFORE :S%Z(', I, ',', J, ',', K, ')=', S%Z(I, J, K)
                  S%Z (I, J, K) = S2%RECVBUF_FACE(LL+1)
 
                  IF (SCARC_DEBUG .GE. 6) write (9, '(a,i2,a,i2,a,i2,a,f22.16)') &
                   'AFTER  :S%Z(', I, ',', J, ',', K, ')=', S%Z(I, J, K)
 
                  LL = LL + 1
               ENDDO
            ENDDO
         ENDDO
 
      ENDDO FACE_UNPACK_LOOP
 
 
!!!----------------------------------------------------------------------------------------
!!! Unpack extended edge-communication data (only 3D-case)
!!!----------------------------------------------------------------------------------------
      IF ( .NOT. TWO_D) THEN
 
         EDGE_UNPACK_LOOP: DO IEDGE = 1, NEDGE
 
            NOM = NBR_EDGE (NM, IEDGE)
            IF (NOM == 0) CYCLE EDGE_UNPACK_LOOP
 
            S => SCARC (NM)
            S2 => SCARC(NM)%OSCARC(NOM)
 
            IF (SCARC_DEBUG .GE. 6) THEN
               WRITE (9,*) 'IEDGE=', IEDGE
               WRITE (9,*) 'NM   =', NM
               WRITE (9,*) 'NOM  =', NOM
               WRITE (9,*) 'S2 zeigt auf SCARC(', NM, ')%OSCARC(', NOM, ')'
               WRITE (9,*) 'S4 zeigt auf SCARC(', NM, ')'
               WRITE (9,*) 'IJKW_EDGE(:'
               WRITE (9, '(12i4)') (IJKW_EDGE(IEDGE, III), III=1, 12)
               CALL flush (9)
            ENDIF
 
            LL = 0
            DO K = IJKW_EDGE (IEDGE, 9), IJKW_EDGE (IEDGE, 12)
               DO J = IJKW_EDGE (IEDGE, 8), IJKW_EDGE (IEDGE, 11)
                  DO I = IJKW_EDGE (IEDGE, 7), IJKW_EDGE (IEDGE, 10)
 
                     IF (SCARC_DEBUG .GE. 6) write (9, '(a,i2,a,i2,a,i2,a,f22.16&
                    &)') 'BEFORE :S%Z(', I, ',', J, ',', K, ')=', S%Z(I, J, K)
                     S%Z (I, J, K) = S2%RECVBUF_EDGE(LL+1)
 
                     IF (SCARC_DEBUG .GE. 6) write (9, '(a,i2,a,i2,a,i2,a,f22.16&
                    &)') 'AFTER  :S%Z(', I, ',', J, ',', K, ')=', S%Z(I, J, K)
 
                     LL = LL + 1
                  ENDDO
               ENDDO
            ENDDO
 
         ENDDO EDGE_UNPACK_LOOP
 
      ENDIF
 
!!!----------------------------------------------------------------------------------------
!!! Unpack extended vertex-communication data
!!!----------------------------------------------------------------------------------------
      VRTX_UNPACK_LOOP: DO IVRTX = 1, NVRTX
 
         NOM = NBR_VRTX (NM, IVRTX)
 
         IF (NOM == 0) CYCLE VRTX_UNPACK_LOOP
 
         S => SCARC (NM)
         S2 => SCARC(NM)%OSCARC(NOM)
 
         IF (SCARC_DEBUG .GE. 6) THEN
            WRITE (9,*) 'IVRTX=', IVRTX
            WRITE (9,*) 'NM   =', NM
            WRITE (9,*) 'NOM  =', NOM
            WRITE (9,*) 'S2 zeigt auf SCARC(', NM, ')%OSCARC(', NOM, ')'
            WRITE (9,*) 'S4 zeigt auf SCARC(', NM, ')'
            WRITE (9,*) 'IJKW_VRTX(:'
            WRITE (9, '(12i4)') (IJKW_VRTX(IVRTX, III), III=1, NVRTX)
            CALL flush (9)
         ENDIF
 
         LL = 0
         DO K = IJKW_VRTX (IVRTX, 9), IJKW_VRTX (IVRTX, 12)
            DO J = IJKW_VRTX (IVRTX, 8), IJKW_VRTX (IVRTX, 11)
               DO I = IJKW_VRTX (IVRTX, 7), IJKW_VRTX (IVRTX, 10)
 
                  IF (SCARC_DEBUG .GE. 6) write (9, '(a,i2,a,i2,a,i2,a,f22.16)') &
                   'BEFORE :S%Z(', I, ',', J, ',', K, ')=', S%Z(I, J, K)
                  S%Z (I, J, K) = S2%RECVBUF_VRTX(LL+1)
 
                  IF (SCARC_DEBUG .GE. 6) write (9, '(a,i2,a,i2,a,i2,a,f22.16)') &
                   'AFTER  :S%Z(', I, ',', J, ',', K, ')=', S%Z(I, J, K)
 
                  LL = LL + 1
               ENDDO
            ENDDO
         ENDDO
 
      ENDDO VRTX_UNPACK_LOOP
 
   ENDDO UNPACK_LOOP
 
!if (SCARC_DEBUG.ge.6) CALL SCARC_COMPARE_SINGLE(S%Z,'SARC',1,'EXCHZ ',0)
 
ENDIF CODE_IF
 
 
 
 
 
 
END SUBROUTINE SCARC_EXCHANGE
 
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! 2D: Set correct boundary conditions for S%AGLOB at exterior boundaries
!!! This is done corresponding to pois.f90
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETBDRY2D (F, BXS, BXF, BYS, BYF, SAVE)
REAL (EB), POINTER, DIMENSION (:, :, :) :: F
REAL (EB), POINTER, DIMENSION (:, :) :: BXS, BXF, BYS, BYF
REAL (EB), DIMENSION (-3:) :: SAVE
REAL (EB) :: DDX, DDY, DLYRCP, TWDYSQ
INTEGER :: L, LP, M, MP, IA, ID, I, J, K
 
 
DDX = SAVE (2)
L = SAVE (3)
LP = SAVE (4)
 
DDY = SAVE (5)
M = SAVE (6)
MP = SAVE (7)
 
DLYRCP = 1._EB / DDY
TWDYSQ = 2._EB / DDY ** 2
 
!                               ALLOCATE SAVE ARRAY
 
IA = 9
ID = IA + 3 * L
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) '========================Starting SCARC_SETBDRY2D'
   WRITE (9,*)
   WRITE (9,*) 'F(.,.)'
   WRITE (9, '(4f10.4)') ((F(I, 1, J), I=1, IBAR), J=KBAR, 1,-1)
   WRITE (9,*) 'DDX=', DDX
   WRITE (9,*) 'L  =', L
   WRITE (9,*) 'LP =', LP
   WRITE (9,*)
   WRITE (9,*) 'DDY=', DDY
   WRITE (9,*) 'M  =', M
   WRITE (9,*) 'MP =', MP
   WRITE (9,*)
   WRITE (9,*) 'DLYRCP=', DLYRCP
   WRITE (9,*) 'TWDYSQ=', TWDYSQ
   WRITE (9,*)
   WRITE (9,*) 'IA=', IA
   WRITE (9,*) 'ID=', ID
   WRITE (9,*)
   WRITE (9,*) 'SAVE(IA)=', SAVE (IA)
   WRITE (9,*) 'SAVE(ID-1)=', SAVE (ID-1)
   WRITE (9,*)
   WRITE (9,*) 'BXS:'
   WRITE (9,*) (BXS(1, J), J=1, M)
   WRITE (9,*)
   WRITE (9,*) 'BXF:'
   WRITE (9,*) (BXF(1, J), J=1, M)
   WRITE (9,*)
   WRITE (9,*) 'BYS:'
   WRITE (9,*) (BYS(I, 1), I=1, L)
   WRITE (9,*)
   WRITE (9,*) 'BYF:'
   WRITE (9,*) (BYF(I, 1), I=1, L)
   WRITE (9,*)
   WRITE (9,*) 'BEFORE CHANGE: F(.,.)'
   WRITE (9, '(4f10.4)') ((F(I, 1, J), I=1, IBAR), J=KBAR, 1,-1)
   WRITE (9,*)
   WRITE (9,*) 'IBAR  =', IBAR
   WRITE (9,*) 'L  =', L
   WRITE (9,*) 'M  =', M
ENDIF
 
 
!                               ENTER BOUNDARY DATA FOR X-BOUNDARIES
 
IF (LP == 2 .OR. LP == 3) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO J = 1, M
      F (1, 1, J) = F (1, 1, J) - 2._EB * BXS (1, J) * SAVE (IA)
      IF (SCARC_DEBUG .GE. 2) write (9, '(a,i2,a,i2,a,f14.8,a,f14.8,a,f14.8)') '&
     &A1:F(', 1, ',', J, ') =', F (1, 1, J), ' - 2*', BXS (1, J), '*', SAVE (IA)
 
   ENDDO
ENDIF
 
IF (LP == 4 .OR. LP == 5) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO J = 1, M
      F (1, 1, J) = F (1, 1, J) + SAVE (IA) * DDX * BXS (1, J)
      IF (SCARC_DEBUG .GE. 2) write (9, '(a,i2,a,i2,a,f14.8,a,f14.8,a,f14.8,a,f1&
     &4.8)') 'A2:F(', 1, ',', J, ') =', F (1, 1, J), ' + ', SAVE (IA), '*', DDX, &
       '*', BXS (1, J)
 
   ENDDO
ENDIF
 
IF (LP == 2 .OR. LP == 5) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO J = 1, M
      F (L, 1, J) = F (L, 1, J) - 2._EB * BXF (1, J) * SAVE (ID-1)
      IF (SCARC_DEBUG .GE. 2) write (9, '(a,i2,a,i2,a,f14.8,a,f14.8,a,f14.8)') '&
     &A3:F(', L, ',', J, ') =', F (L, 1, J), ' - 2*', BXF (1, J), '*', SAVE &
       (ID-1)
 
   ENDDO
ENDIF
 
IF (LP == 3 .OR. LP == 4) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO J = 1, M
      F (L, 1, J) = F (L, 1, J) - SAVE (ID-1) * DDX * BXF (1, J)
      IF (SCARC_DEBUG .GE. 2) write (9, '(a,i2,a,i2,a,f14.8,a,f14.8,a,f14.8,a,f1&
     &4.8)') 'A4:F(', L, ',', J, ') =', F (L, 1, J), ' - ', SAVE (ID-1), '*', &
       DDX, '*', BXF (1, J)
 
   ENDDO
ENDIF
 
 
!                               ENTER BOUNDARY DATA FOR Y-BOUNDARIES
 
IF (MP == 2 .OR. MP == 3) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO I = 1, L
      F (I, 1, 1) = F (I, 1, 1) - BYS (I, 1) * TWDYSQ
      IF (SCARC_DEBUG .GE. 2) write (9, '(a,i2,a,i2,a,f14.8,a,f14.8,a,f14.8)') '&
     &B1:F(', I, ',', 1, ') =', F (I, 1, 1), ' - ', BYS (I, 1), '*', TWDYSQ
 
   ENDDO
ENDIF
 
IF (MP == 4 .OR. MP == 5) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO I = 1, L
      F (I, 1, 1) = F (I, 1, 1) + BYS (I, 1) * DLYRCP
      IF (SCARC_DEBUG .GE. 2) write (9, '(a,i2,a,i2,a,f14.8,a,f14.8,a,f14.8)') '&
     &B2:F(', I, ',', 1, ') =', F (I, 1, 1), ' + ', BYS (I, 1), '*', DLYRCP
 
   ENDDO
ENDIF
 
IF (MP == 2 .OR. MP == 5) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO I = 1, L
      F (I, 1, M) = F (I, 1, M) - BYF (I, 1) * TWDYSQ
      IF (SCARC_DEBUG .GE. 2) write (9, '(a,i2,a,i2,a,f14.8,a,f14.8,a,f14.8)') '&
     &B3:F(', I, ',', M, ') =', F (I, 1, M), ' - ', BYF (I, 1), '*', TWDYSQ
 
   ENDDO
ENDIF
 
IF (MP == 3 .OR. MP == 4) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO I = 1, L
      F (I, 1, M) = F (I, 1, M) - BYF (I, 1) * DLYRCP
      IF (SCARC_DEBUG .GE. 2) write (9, '(a,i2,a,i2,a,f14.8,a,f14.8,a,f14.8)') '&
     &B4:F(', I, ',', M, ') =', F (I, 1, M), ' - ', BYF (I, 1), '*', DLYRCP
 
   ENDDO
ENDIF
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*)
   WRITE (9,*) 'AFTER CHANGE: F(.,.)'
   WRITE (9,*) 'IBAR=', IBAR
   WRITE (9,*) 'KBAR=', KBAR
   WRITE (9, '(4f10.4)') ((F(I, 1, J), I=1, IBAR), J=KBAR, 1,-1)
ENDIF
 
END SUBROUTINE SCARC_SETBDRY2D
 
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! 3D: Set correct boundary conditions for S%AGLOB at exterior boundaries
!!! This is done corresponding to pois.f90
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETBDRY3D (F, BDXS, BDXF, BDYS, BDYF, BDZS, BDZF, SAVE)
REAL (EB), POINTER, DIMENSION (:, :, :) :: F
REAL (EB), POINTER, DIMENSION (:, :) :: BDXS, BDXF, BDYS, BDYF, BDZS, BDZF
REAL (EB), DIMENSION (-3:) :: SAVE
REAL (EB) :: DDX, DDY, DDZ
REAL (EB) :: DLYRCP, TWDYSQ
REAL (EB) :: DLZRCP, TWDZSQ
INTEGER :: L, LP, M, MP, N, NP, IA, ID, I, J, K
 
 
DDX = SAVE (2)
L = SAVE (3)
LP = SAVE (4)
 
DDY = SAVE (5)
M = SAVE (6)
MP = SAVE (7)
 
DDZ = SAVE (8)
N = SAVE (9)
NP = SAVE (10)
 
DLYRCP = 1._EB / DDY
TWDYSQ = 2._EB / DDY ** 2
 
DLZRCP = 1._EB / DDZ
TWDZSQ = 2._EB / DDZ ** 2
 
!                               ALLOCATE SAVE ARRAY
 
IA = 12
ID = IA + 3 * L
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) '========================Starting SCARC_SETBDRY3D'
   WRITE (9,*)
   WRITE (9,*) 'DDX=', DDX
   WRITE (9,*) 'L  =', L
   WRITE (9,*) 'LP =', LP
   WRITE (9,*)
   WRITE (9,*) 'DDY=', DDY
   WRITE (9,*) 'M  =', M
   WRITE (9,*) 'MP =', MP
   WRITE (9,*)
   WRITE (9,*) 'DDZ=', DDZ
   WRITE (9,*) 'N  =', N
   WRITE (9,*) 'NP =', NP
   WRITE (9,*)
   WRITE (9,*) 'DLYRCP=', DLYRCP
   WRITE (9,*) 'TWDYSQ=', TWDYSQ
   WRITE (9,*)
   WRITE (9,*) 'DLZRCP=', DLZRCP
   WRITE (9,*) 'TWDZSQ=', TWDZSQ
   WRITE (9,*)
   WRITE (9,*) 'IA=', IA
   WRITE (9,*) 'ID=', ID
   WRITE (9,*)
   WRITE (9,*) 'SAVE(IA)=', SAVE (IA)
   WRITE (9,*) 'SAVE(ID-1)=', SAVE (ID-1)
   WRITE (9,*)
   WRITE (9,*) 'BDXS:'
   WRITE (9, '(4f9.4)') ((BDXS(J, K), J=1, M), K=1, N)
   WRITE (9,*)
   WRITE (9,*) 'BDXF:'
   WRITE (9, '(4f9.4)') ((BDXF(J, K), J=1, M), K=1, N)
   WRITE (9,*)
   WRITE (9,*) 'BDYS:'
   WRITE (9, '(4f9.4)') ((BDYS(I, K), I=1, L), K=1, N)
   WRITE (9,*)
   WRITE (9,*) 'BDYF:'
   WRITE (9, '(4f9.4)') ((BDYF(I, K), I=1, L), K=1, N)
   WRITE (9,*)
   WRITE (9,*) 'BDZS:'
   WRITE (9, '(4f9.4)') ((BDZS(I, J), I=1, L), J=1, M)
   WRITE (9,*)
   WRITE (9,*) 'BDZF:'
   WRITE (9, '(4f9.4)') ((BDZF(I, J), I=1, L), J=1, M)
   WRITE (9,*)
   WRITE (9,*) 'BEFORE CHANGE: F(.,.)'
   DO K = KBAR, 1, - 1
      WRITE (9, '(4f9.4)') ((F(I, J, K), I=1, IBAR), J=1, JBAR)
      WRITE (9,*) '--------------------------------------------'
   ENDDO
   WRITE (9,*)
   WRITE (9,*) 'IBAR  =', IBAR
   WRITE (9,*) 'L  =', L
   WRITE (9,*) 'M  =', M
   WRITE (9,*) 'N  =', N
ENDIF
 
 
!                               ENTER BOUNDARY DATA FOR X-BOUNDARIES
 
IF (LP == 2 .OR. LP == 3) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO K = 1, N
      DO J = 1, M
         F (1, J, K) = F (1, J, K) - 2._EB * BDXS (J, K) * SAVE (IA)
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,f14.8,a,f14.8,a,f1&
        &4.8)') 'A1:F(', 1, ',', J, ',', K, ') =', F (1, J, K), ' - 2*', BDXS &
          (J, K), '*', SAVE (IA)
      ENDDO
   ENDDO
ENDIF
 
IF (LP == 4 .OR. LP == 5) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO K = 1, N
      DO J = 1, M
         F (1, J, K) = F (1, J, K) + SAVE (IA) * DDX * BDXS (J, K)
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,f14.8,a,f14.8,a,f1&
        &4.8,a,f14.8)') 'A2:F(', 1, ',', J, ',', K, ') =', F (1, J, K), ' + ', &
          SAVE (IA), '*', DDX, '*', BDXS (J, K)
      ENDDO
   ENDDO
ENDIF
 
IF (LP == 2 .OR. LP == 5) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO K = 1, N
      DO J = 1, M
         F (L, J, K) = F (L, J, K) - 2._EB * BDXF (J, K) * SAVE (ID-1)
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,f14.8,a,f14.8,a,f1&
        &4.8)') 'A3:F(', L, ',', J, ',', K, ') =', F (L, J, K), ' - 2*', BDXF &
          (J, K), '*', SAVE (ID-1)
      ENDDO
   ENDDO
ENDIF
 
IF (LP == 3 .OR. LP == 4) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO K = 1, N
      DO J = 1, M
         F (L, J, K) = F (L, J, K) - SAVE (ID-1) * DDX * BDXF (J, K)
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,f14.8,a,f14.8,a,f1&
        &4.8,a,f14.8)') 'A4:F(', L, ',', J, ',', K, ') =', F (L, J, K), ' - ', &
          SAVE (ID-1), '*', DDX, '*', BDXF (J, K)
      ENDDO
   ENDDO
ENDIF
 
 
!                               ENTER BOUNDARY DATA FOR Y-BOUNDARIES
 
IF (MP == 2 .OR. MP == 3) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO K = 1, N
      DO I = 1, L
         F (I, 1, K) = F (I, 1, K) - BDYS (I, K) * TWDYSQ
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,f14.8,a,f14.8,a,f1&
        &4.8)') 'B1:F(', I, ',', 1, ',', K, ') =', F (I, 1, K), ' - ', BDYS (I, &
          K), '*', TWDYSQ
      ENDDO
   ENDDO
ENDIF
 
IF (MP == 4 .OR. MP == 5) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO K = 1, N
      DO I = 1, L
         F (I, 1, K) = F (I, 1, K) + BDYS (I, K) * DLYRCP
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,f14.8,a,f14.8,a,f1&
        &4.8)') 'B2:F(', I, ',', 1, ',', K, ') =', F (I, 1, K), ' - ', BDYS (I, &
          K), '*', DLYRCP
      ENDDO
   ENDDO
ENDIF
 
IF (MP == 2 .OR. MP == 5) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO K = 1, N
      DO I = 1, L
         F (I, M, K) = F (I, M, K) - BDYF (I, K) * TWDYSQ
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,f14.8,a,f14.8,a,f1&
        &4.8)') 'B3:F(', I, ',', M, ',', K, ') =', F (I, M, K), ' - ', BDYF (I, &
          K), '*', TWDYSQ
      ENDDO
   ENDDO
ENDIF
 
IF (MP == 3 .OR. MP == 4) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO K = 1, N
      DO I = 1, L
         F (I, M, K) = F (I, M, K) - BDYF (I, K) * DLYRCP
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,f14.8,a,f14.8,a,f1&
        &4.8)') 'B4:F(', I, ',', M, ',', K, ') =', F (I, M, K), ' - ', BDYF (I, &
          K), '*', DLYRCP
      ENDDO
   ENDDO
ENDIF
 
 
!                               ENTER BOUNDARY DATA FOR Z-BOUNDARIES
 
IF (NP == 2 .OR. NP == 3) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO J = 1, M
      DO I = 1, L
         F (I, J, 1) = F (I, J, 1) - BDZS (I, J) * TWDZSQ
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,f14.8,a,f14.8,a,f1&
        &4.8)') 'C1:F(', I, ',', J, ',', 1, ') =', F (I, J, 1), ' - ', BDZS (I, &
          J), '*', TWDZSQ
      ENDDO
   ENDDO
ENDIF
 
IF (NP == 4 .OR. NP == 5) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO J = 1, M
      DO I = 1, L
         F (I, J, 1) = F (I, J, 1) + BDZS (I, J) * DLZRCP
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,f14.8,a,f14.8,a,f1&
        &4.8)') 'C2:F(', I, ',', J, ',', 1, ') =', F (I, J, 1), ' - ', BDZS (I, &
          J), '*', DLZRCP
      ENDDO
   ENDDO
ENDIF
 
IF (NP == 2 .OR. NP == 5) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO J = 1, M
      DO I = 1, L
         F (I, J, N) = F (I, J, N) - BDZF (I, J) * TWDZSQ
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,f14.8,a,f14.8,a,f1&
        &4.8)') 'C3:F(', I, ',', J, ',', N, ') =', F (I, J, N), ' - ', BDZF (I, &
          J), '*', TWDZSQ
      ENDDO
   ENDDO
ENDIF
 
IF (NP == 3 .OR. NP == 4) THEN
   IF (SCARC_DEBUG .GE. 2) write (9,*)
   DO J = 1, M
      DO I = 1, L
         F (I, J, N) = F (I, J, N) - BDZF (I, J) * DLZRCP
         IF (SCARC_DEBUG .GT. 2) write (9, '(a,i2,a,i2,a,i2,a,f14.8,a,f14.8,a,f1&
        &4.8)') 'C4:F(', I, ',', J, ',', N, ') =', F (I, J, N), ' - ', BDZF (I, &
          J), '*', DLZRCP
      ENDDO
   ENDDO
ENDIF
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*)
   WRITE (9,*) 'AFTER CHANGE: F(.,.)'
   WRITE (9,*) 'IBAR=', IBAR
   WRITE (9,*) 'JBAR=', JBAR
   WRITE (9,*) 'KBAR=', KBAR
   DO K = N, 1, - 1
      WRITE (9, '(4f18.12)') ((F(I, J, K), I=1, L), J=1, M)
      WRITE (9,*) '--------------------------------------------'
   ENDDO
   WRITE (9,*) 'STOPPING SETBDRY3D'
ENDIF
 
END SUBROUTINE SCARC_SETBDRY3D
 
 
SUBROUTINE SCARC_UPDATE (MYID, ITYPE)
INTEGER :: MYID, ITYPE, NM
 
NM = MYID + 1
 
IF (ITYPE == NUPDATE_FULL) THEN
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%U, NM, 'U     ')
   IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%V, NM, 'V     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%W, NM, 'W     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%H, NM, 'H     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%HS, NM, 'HS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%D, NM, 'D     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%RHO, NM, 'RHO   ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%MU, NM, 'MU    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%US, NM, 'US    ')
   IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%VS, NM, 'VS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%WS, NM, 'WS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DS, NM, 'DS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%RHOS, NM, 'RHOS  ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DS, NM, 'DS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DDDT, NM, 'DDDT  ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%Q, NM, 'Q     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%QR, NM, 'QR    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%TMP, NM, 'TMP   ')
ELSE IF (ITYPE == NUPDATE_MEDIUM) THEN
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%U, NM, 'U     ')
   IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%V, NM, 'V     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%W, NM, 'W     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%H, NM, 'H     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%HS, NM, 'HS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%D, NM, 'D     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%US, NM, 'US    ')
   IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%VS, NM, 'VS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%WS, NM, 'WS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DS, NM, 'DS    ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DS, NM, 'DS   ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%DDDT, NM, 'DDDT   ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%TMP, NM, 'TMP    ')
ELSE IF (ITYPE == NUPDATE_VEL) THEN
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%U, NM, 'U     ')
   IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%V, NM, 'V     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%W, NM, 'W     ')
ELSE IF (ITYPE == NUPDATE_FV) THEN
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVX, NM, 'FVX   ')
   IF ( .NOT. TWO_D) CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVY, NM, 'FVY   ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%FVZ, NM, 'FVZ   ')
ELSE IF (ITYPE == NUPDATE_Q) THEN
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%Q, NM, 'Q     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%QR, NM, 'QR    ')
ELSE IF (ITYPE == NUPDATE_H) THEN
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%H , NM, 'H     ')
   CALL SCARC_UPDATE_QUANTITY (MESHES(NM)%HS, NM, 'HS    ')
ENDIF
 
 
END SUBROUTINE SCARC_UPDATE
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  update interior boundary values along vertical directions
!!!  to have consistent values all over the domain
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE_QUANTITY (Z, NM, CNAME)
REAL (EB), POINTER, DIMENSION (:, :, :) :: Z
INTEGER :: NM, I, J, K, IERR, III, JJJ, KKK
REAL (EB) :: VERTEX (0:7), VERTEX_ALL (0:7, 0:100)
CHARACTER (4) :: CNAME
 
INCLUDE 'mpif.h'
 
 
IF (NMESHES > 1) THEN
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF (TWO_D) THEN
 
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            S%Z (I, 1, K) = Z (I, 1, K)
         ENDDO
      ENDDO
 
      IF (SCARC_DEBUG .GE. 2) THEN
         WRITE (9,*) '//////////////////////////////////////////////////////////'
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: ==================================='
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: NMESHES=', NMESHES
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: NM=', NM
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: CNAME=', CNAME
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: TWO_D=', TWO_D
!write(9,*) 'NBR_VRTX:'
!write(9,'( 4i4)') ((NBR_VRTX(JJJ,III),III=1,NVRTX),JJJ=1,NMESHES)
         WRITE (9,*) 'NBR_VRTX(', NM, '):'
         WRITE (9, '( 4i4)') (NBR_VRTX(NM, III), III=1, NVRTX)
!write(9,*) 'NBR_FACE:'
!write(9,'( 4i4)') ((NBR_FACE(JJJ,III),III=1,4),JJJ=1,NMESHES)
         WRITE (9,*) 'NBR_FACE(', NM, '):'
         WRITE (9, '( 4i4)') (NBR_FACE(NM, III), III=1, 4)
         WRITE (9,*)
         WRITE (9,*) '-------SCARC_COPYBDRY : BEFORE EXCHANGE Z'
         DO K = 5, 0, - 1
            IF (NM == 2) THEN
               WRITE (9, '(a,i3,a,4f22.16)') 'k=', K, ' : ', (Z(I, 1, K), I=14, IBAR+1)
            ELSE
               WRITE (9, '(a,i3,a,4f22.16)') 'k=', K, ' : ', (Z(I, 1, K), I=0, 3)
            ENDIF
  !write(9,*) '---------------------------------------'
         ENDDO
         WRITE (9,*) '-------SCARC_COPYBDRY : BEFORE EXCHANGE S%Z'
         DO K = 5, 0, - 1
            IF (NM == 2) THEN
               WRITE (9, '(a,i3,a,4f22.16)') 'k=', K, ' : ', (S%Z(I, 1, K), I=14, IBAR+1)
            ELSE
               WRITE (9, '(a,i3,a,4f22.16)') 'k=', K, ' : ', (S%Z(I, 1, K), I=0, 3)
            ENDIF
  !write(9,*) '---------------------------------------'
         ENDDO
         CALL flush (9)
      ENDIF
 
 
      N_REQ_FACE = 0
      N_REQ_VRTX = 0
 
!if (CNAME(1:2)=='H4') THEN
!close(9)
!stop
!ENDIF
 
      CALL SCARC_RECEIVE (2, 0)
      CALL SCARC_EXCHANGE (2, 0)
 
 
      DO K = 0, KBAR + 1
         DO I = 0, IBAR + 1
            Z (I, 1, K) = S%Z(I, 1, K)
         ENDDO
      ENDDO
 
      IF (SCARC_DEBUG .GE. 2) THEN
         WRITE (9,*) '-------SCARC_COPYBDRY : AFTER EXCHANGE S%Z'
         DO K = 5, 0, - 1
            IF (NM == 2) THEN
               WRITE (9, '(a,i3,a,4f22.16)') 'k=', K, ' : ', (S%Z(I, 1, K), I=14, IBAR+1)
            ELSE
               WRITE (9, '(a,i3,a,4f22.16)') 'k=', K, ' : ', (S%Z(I, 1, K), I=0, 3)
            ENDIF
  !write(9,*) '---------------------------------------'
         ENDDO
         CALL flush (9)
      ENDIF
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ELSE
 
      DO K = 0, KBAR + 1
         DO J = 0, JBAR + 1
            DO I = 0, IBAR + 1
               S%Z (I, J, K) = Z (I, J, K)
            ENDDO
         ENDDO
      ENDDO
 
      IF (SCARC_DEBUG .GE. 6) THEN
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: ==================================='
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: NMESHES=', NMESHES
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: NM=', NM
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: CNAME=', CNAME
         WRITE (9,*) 'SCARC_UPDATE_QUANTITY: TWO_D=', TWO_D
         WRITE (9,*) 'NBR_VRTX:'
         WRITE (9, '( 8i4)') ((NBR_VRTX(JJJ, III), III=1, NVRTX), JJJ=1, NMESHES)
         WRITE (9,*) 'NBR_VRTX(', NM, '):'
         WRITE (9, '( 8i4)') (NBR_VRTX(NM, III), III=1, NVRTX)
         WRITE (9,*) 'NBR_EDGE:'
         WRITE (9, '( 12i4)') ((NBR_EDGE(JJJ, III), III=1, 12), JJJ=1, NMESHES)
         WRITE (9,*) 'NBR_EDGE(', NM, '):'
         WRITE (9, '( 12i4)') (NBR_EDGE(NM, III), III=1, 12)
         WRITE (9,*)
         WRITE (9,*) '-------SCARC_COPYBDRY : BEFORE EXCHANGE Z'
         DO K = KBAR + 1, 0, - 1
            WRITE (9, '(6f10.5)') ((Z(I, J, K), I=0, IBAR+1), J=JBAR+1, 0,-1)
            WRITE (9,*) '---------------------------------------'
         ENDDO
         WRITE (9,*) '-------SCARC_COPYBDRY : BEFORE EXCHANGE S%Z'
         DO K = KBAR + 1, 0, - 1
            WRITE (9, '(6f10.5)') ((S%Z(I, J, K), I=0, IBAR+1), J=JBAR+1, 0,-1)
            WRITE (9,*) '---------------------------------------'
         ENDDO
         CALL flush (9)
      ENDIF
 
 
      N_REQ_FACE = 0
      N_REQ_EDGE = 0
      N_REQ_VRTX = 0
 
      CALL SCARC_RECEIVE (2, 0)
      CALL SCARC_EXCHANGE (2, 0)
 
 
      DO K = 0, KBAR + 1
         DO J = 0, JBAR + 1
            DO I = 0, IBAR + 1
               Z (I, J, K) = S%Z(I, J, K)
            ENDDO
         ENDDO
      ENDDO
 
      IF (SCARC_DEBUG .GE. 6) THEN
         WRITE (9,*) '-------SCARC_COPYBDRY : AFTER EXCHANGE S%Y'
         IF (NMESHES == 1) THEN
            DO K = KBAR + 1, 0, - 1
               WRITE (9, '(10f10.5)') ((S%Z(I, J, K), I=0, IBAR+1), J=JBAR+1, 0,-1)
               WRITE (9,*) '---------------------------------------'
            ENDDO
         ELSE
            DO K = KBAR + 1, 0, - 1
               WRITE (9, '(6f10.5)') ((S%Z(I, J, K), I=0, IBAR+1), J=JBAR+1, 0,-1)
               WRITE (9,*) '---------------------------------------'
            ENDDO
         ENDIF
         CALL flush (9)
      ENDIF
 
   ENDIF
 
ENDIF
 
END SUBROUTINE SCARC_UPDATE_QUANTITY
 
 
 
 
SUBROUTINE SCARC_VECADD_EXTENDED (X, Y, A1, A2)
 
REAL (EB) :: A1, A2
REAL (EB), DIMENSION (0:, 0:, 0:) :: X, Y
INTEGER :: I, J, K
 
 
IF (SCARC_DEBUG .GT. 3) THEN
   WRITE (9,*) 'BEFORE VECADD_EXTENDED'
   WRITE (9,*) 'A1=', A1
   WRITE (9,*) 'A2=', A2
   WRITE (9,*) 'X='
   DO K = KBAR + 1, 0, - 1
      WRITE (9, '(6f10.5)') ((X(I, J, K), I=0, IBAR+1), J=0, JBAR+1)
      WRITE (9,*) '---------------------------------------'
   ENDDO
   WRITE (9,*) 'Y='
   DO K = KBAR + 1, 0, - 1
      WRITE (9, '(6f10.5)') ((Y(I, J, K), I=0, IBAR+1), J=0, JBAR+1)
      WRITE (9,*) '---------------------------------------'
   ENDDO
ENDIF
 
 
 
IF (A2 == 0.0_EB) THEN
   IF (A1 /= 1.0_EB) THEN
      DO K = 0, KBAR + 1
         DO J = 0, JBAR + 1
            DO I = 0, IBAR + 1
               Y (I, J, K) = A1 * X (I, J, K)
!WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO K = 0, KBAR + 1
         DO J = 0, JBAR + 1
            DO I = 0, IBAR + 1
               Y (I, J, K) = X (I, J, K)
!WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ELSE
   IF (A1 /= 1.0_EB) THEN
      DO K = 0, KBAR + 1
         DO J = 0, JBAR + 1
            DO I = 0, IBAR + 1
               Y (I, J, K) = A1 * X (I, J, K) + A2 * Y (I, J, K)
!WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO K = 0, KBAR + 1
         DO J = 0, JBAR + 1
            DO I = 0, IBAR + 1
               Y (I, J, K) = X (I, J, K) + A2 * Y (I, J, K)
!WRITE(9,*) 'Y(',I,',J,',K,')=',Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDIF
 
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   WRITE (9,*) 'AFTER VECADD_EXTENDED'
   DO K = KBAR + 1, 0, - 1
      WRITE (9, '(6f10.5)') ((Y(I, J, K), I=0, IBAR+1), J=0, JBAR+1)
      WRITE (9,*) '---------------------------------------'
   ENDDO
ENDIF
 
END SUBROUTINE SCARC_VECADD_EXTENDED
 
 
 
LOGICAL FUNCTION MATCH (VAL1, VAL2)
REAL (EB) :: VAL1, VAL2, TOL
TOL = 1.0E-8_EB
MATCH = .FALSE.
IF (Abs(VAL1-VAL2) <= TOL) MATCH = .TRUE.
END FUNCTION MATCH
 
 
 
SUBROUTINE SCARC_COMPARE_ALL (MYID, CROUTINE, INUM, ILINE)
INTEGER :: MYID, INUM, ILINE
CHARACTER (4) :: CROUTINE
 
call flush(6)
call flush(9)

CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%H, CROUTINE, INUM, 'H     ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%MU, CROUTINE, INUM, 'MU    ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%DDDT, CROUTINE, INUM, 'DDDT  ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%D, CROUTINE, INUM, 'D     ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%DS, CROUTINE, INUM, 'DS    ',ILINE)

CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%U, CROUTINE, INUM, 'U      ',ILINE)
IF ( .NOT. TWO_D) CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%V, CROUTINE, INUM, 'V      ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%W, CROUTINE, INUM, 'W      ',ILINE)
 
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%US, CROUTINE, INUM, 'US     ',ILINE)
IF ( .NOT. TWO_D) CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%VS, CROUTINE, INUM, 'VS     ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%WS, CROUTINE, INUM, 'WS     ',ILINE)
 
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%FVX, CROUTINE, INUM, 'FVX    ',ILINE)
IF ( .NOT. TWO_D) CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%FVY, CROUTINE, INUM, 'FVY    ',ILINE)
CALL SCARC_COMPARE_SINGLE (MESHES(MYID+1)%FVZ, CROUTINE, INUM, 'FVZ    ',ILINE)
 
END SUBROUTINE SCARC_COMPARE_ALL
 
 
SUBROUTINE SCARC_COMPARE_SINGLE (X, CROUTINE, INUM, CNAME, ILINE)
REAL (EB), DIMENSION (:, :, :) :: X
INTEGER :: INUM, ILINE
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
IF (SCARC_COMPARE .GE. 1 .AND. CROUTINE(1:4) /= "SARC") THEN

   SCARC_MSG = "msg/  /    _      _      "
   WRITE (SCARC_MSG(5:6), '(i2.2)') MYID + 1
   WRITE (SCARC_MSG(8:11), '(i4.4)') SCARC_COUNT2
   SCARC_MSG (13:16) = CROUTINE (1:4)
   WRITE (SCARC_MSG(17:18), '(i2.2)') INUM
   SCARC_MSG (20:25) = CNAME (1:6)
   SCARC_COUNT2 = SCARC_COUNT2 + 1
   OPEN (1000, FILE=SCARC_MSG)
 
   IF (TWO_D) THEN
      WRITE (1000,*) '=========================================================='
      WRITE (1000,'(a20,a4)') 'ROUTINE          : ',CROUTINE(1:4)
      WRITE (1000,'(a20,a6)') 'VECTOR           : ',CNAME(1:6)
      WRITE (1000,'(a20,i6)') 'CODE-LINE        : ',ILINE
      WRITE (1000,'(a20,i6)') 'INTERNAL POSITION: ',INUM
      WRITE (1000,'(a20,i6)') 'INTERNAL NUMBER  : ',SCARC_COUNT2-1
      WRITE (1000,*) '=========================================================='
      IF (NMESHES == 1) THEN
         DO KK = KBAR, 1, - 1
            WRITE (1000, '(8f22.16)') (X(II, 1, KK), II=1, IBAR)
         ENDDO
      ELSE
         DO KK = KBAR, 1, - 1
            WRITE (1000, '(4f22.16)') (X(II, 1, KK), II=1, IBAR)
         ENDDO
      ENDIF
   ELSE
      WRITE (1000,*) '=========================================================='
      WRITE (1000,*) '===        SCARC: ', CNAME
      WRITE (1000,*) '=========================================================='
      IF (NMESHES == 1) THEN
         DO KK = KBAR, 1, - 1
            WRITE (1000, '(8f22.16)') ((X(II, jj, KK), II=1, IBAR), jj=JBAR,1,-1)
            WRITE (1000,*) '----------------------------------------'
         ENDDO
      ELSE
         DO KK = KBAR, 1, - 1
            WRITE (1000, '(4f22.16)') ((X(II, jj, KK), II=1, IBAR), jj=JBAR,1,-1)
            WRITE (1000,*) '----------------------------------------'
         ENDDO
      ENDIF
   ENDIF

   CALL flush (1000)
   CLOSE (1000)
ENDIF
 
IF (SCARC_DEBUG .GE. 2) THEN
   IF (TWO_D) THEN
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': SCARC_COUNT2= ',SCARC_COUNT2-1
      WRITE (9,*) '============================================================='
      write(9,*) 'IBAR=',IBAR
      write(9,*) 'KBAR=',KBAR
      IF (NMESHES == 1) THEN
         DO KK = KBAR, 1, - 1
            WRITE (9, '(a,i3,a,8e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=1,IBAR)
            !WRITE (9,*) 'k=', KK, ' : ', (X(II, 1, KK), II=1,IBAR)
         ENDDO
      ELSE
         DO KK = KBAR, 1, - 1
            WRITE (9, '(a,i3,a,8e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=1, IBAR)
            !WRITE (9,*) 'k=', KK, ' : ', (X(II, 1, KK), II=1, IBAR)
         ENDDO
      ENDIF
   ELSE
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': SCARC_COUNT2= ',SCARC_COUNT2-1
      WRITE (9,*) '============================================================='
      IF (NMESHES == 1) THEN
         DO KK = KBAR, 1, - 1
            WRITE (9, '(8e12.3)') ((X(II, jj, KK), II=1, IBAR), jj=JBAR, 1,-1)
            WRITE (9,*) '----------------------------------------'
         ENDDO
      ELSE
         DO KK = KBAR, 1, - 1
            WRITE (9, '(4e12.3)') ((X(II, jj, KK), II=1, IBAR), jj=JBAR, 1,-1)
            WRITE (9,*) '----------------------------------------'
         ENDDO
      ENDIF
   ENDIF
   CALL flush (9)
ENDIF
 
 
 
END SUBROUTINE SCARC_COMPARE_SINGLE
 
 
SUBROUTINE SCARC_COMPARE_SINGLE0(X, CROUTINE, INUM, CNAME, ILINE)
REAL (EB), DIMENSION (0:,0 :,0 :) :: X
INTEGER :: INUM, ILINE
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
 
IF (SCARC_DEBUG .GE. 2) THEN
   IF (TWO_D) THEN
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE
      WRITE (9,*) '============================================================='
      write(9,*) 'IBAR=',IBAR
      write(9,*) 'KBAR=',KBAR
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0,IBAR+1)
            !WRITE (9,*) 'k=', KK, ' : ', (X(II, 1, KK), II=1,IBAR)
         ENDDO
      ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,6e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0, IBAR+1)
            !WRITE (9,*) 'k=', KK, ' : ', (X(II, 1, KK), II=1, IBAR)
         ENDDO
      ENDIF
   ELSE
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': SCARC_COUNT2= ',SCARC_COUNT2-1
      WRITE (9,*) '============================================================='
      IF (NMESHES == 1) THEN
         DO KK = KBAR, 1, - 1
            WRITE (9, '(8e12.3)') ((X(II, jj, KK), II=1, IBAR), jj=JBAR, 1,-1)
            WRITE (9,*) '----------------------------------------'
         ENDDO
      ELSE
         DO KK = KBAR, 1, - 1
            WRITE (9, '(4e12.3)') ((X(II, jj, KK), II=1, IBAR), jj=JBAR, 1,-1)
            WRITE (9,*) '----------------------------------------'
         ENDDO
      ENDIF
   ENDIF
   CALL flush (9)
ENDIF
 
 
 
END SUBROUTINE SCARC_COMPARE_SINGLE0

SUBROUTINE SCARC_COMPARE_SINGLE1(X, CROUTINE, INUM, CNAME, ILINE)
REAL (EB), DIMENSION (0:,0 :,0 :) :: X
INTEGER :: INUM, ILINE
INTEGER :: II, jj, KK
CHARACTER (4) :: CROUTINE
CHARACTER (6) :: CNAME
 
 
IF (SCARC_DEBUG .GE. 2) THEN
      WRITE (9,*) '============================================================='
      WRITE (9,*) '===   COMPARE= ',CNAME,': ROUTINE= ',CROUTINE
      WRITE (9,*) '============================================================='
      WRITE (9,*) '-------------------- J = 0 ---------------------------------'
      write(9,*) 'IBAR=',IBAR
      write(9,*) 'KBAR=',KBAR
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 0, KK), II=0,IBAR+1)
         ENDDO
      ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,6e12.3)') 'k=', KK, ' : ', (X(II, 0, KK), II=0, IBAR+1)
         ENDDO
      ENDIF
      WRITE (9,*) '-------------------- J = 1 ---------------------------------'
      write(9,*) 'IBAR=',IBAR
      write(9,*) 'KBAR=',KBAR
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0,IBAR+1)
         ENDDO
      ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,6e12.3)') 'k=', KK, ' : ', (X(II, 1, KK), II=0, IBAR+1)
         ENDDO
      ENDIF
      WRITE (9,*) '-------------------- J = 2 ---------------------------------'
      write(9,*) 'IBAR=',IBAR
      write(9,*) 'KBAR=',KBAR
      IF (NMESHES == 1) THEN
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,10e12.3)') 'k=', KK, ' : ', (X(II, 2, KK), II=0,IBAR+1)
         ENDDO
      ELSE
         DO KK = KBAR+1, 0, - 1
            WRITE (9, '(a,i3,a,6e12.3)') 'k=', KK, ' : ', (X(II, 2, KK), II=0, IBAR+1)
         ENDDO
      ENDIF
   CALL flush (9)
ENDIF
 
 
 
END SUBROUTINE SCARC_COMPARE_SINGLE1
 
 
SUBROUTINE SCARC_ANALYTICAL_SOLUTION(NM)
IMPLICIT NONE
! Initialize flow variables with an analytical solution of the governing equations

INTEGER, INTENT(IN) :: NM
INTEGER :: I,J,K
REAL(EB) :: UU,WW, SCAL, SCAL2, SCAL4, RHO_R, RADIUS, VNU, THETA
REAL(EB) :: SHIFT_X, SHIFT_Y, SHIFT_Z, U_MAX

CALL POINT_TO_MESH(NM)

SELECT CASE(SCARC_CASE)

!-------------------------------------------------------------------------------------
! CD_NSA_2D
!-------------------------------------------------------------------------------------
   CASE(0)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               U(I,J,K) = 1._EB - 2._EB*COS(X(I))*SIN(ZC(K))
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               W(I,J,K) = 1._EB + 2._EB*SIN(XC(I))*COS(Z(K))
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
               UU = 1._EB - 2._EB*COS(XC(I))*SIN(ZC(K))
               WW = 1._EB + 2._EB*SIN(XC(I))*COS(ZC(K))
               !H(I,J,K) = -( COS(2._EB*XC(I)) + COS(2._EB*ZC(K)) ) + 0.5_EB*(UU**2+WW**2)    
            ENDDO
         ENDDO
      ENDDO

!-------------------------------------------------------------------------------------
! VD_NSA_2D
! Bitte alles noch mal nachprüfen !!!!!!!!!!!!!
!-------------------------------------------------------------------------------------
   CASE(1)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               U(I,J,K) = 1._EB - 2._EB*COS(X(I))*SIN(ZC(K))
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               W(I,J,K) = 1._EB + 2._EB*SIN(XC(I))*COS(Z(K))
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
               UU = 1._EB - 2._EB*COS(XC(I))*SIN(ZC(K))
               WW = 1._EB + 2._EB*SIN(XC(I))*COS(ZC(K))
               !H(I,J,K) = -( COS(2._EB*XC(I)) + COS(2._EB*ZC(K)) ) + 0.5_EB*(UU**2+WW**2)   
               RHO(I,J,K) = 1.0000004463756595_EB + 1._EB * COS(XC(I))**2 * COS(ZC(K))**2

            ENDDO
         ENDDO
      ENDDO


!-------------------------------------------------------------------------------------
! CD_VA_2D
!-------------------------------------------------------------------------------------
   CASE(2)
              WRITE(*,*) '------------------- alter Testfall --------------------------'
              U_MAX = 1.0_EB
              DO K=1,KBAR
                 DO J=1,JBAR
                    DO I=0,IBAR
                       SHIFT_X= X(I) - 0.5_EB
                       SHIFT_Y= 0.0_EB
                       SHIFT_Z= ZC(K) - 0.5_EB
                       RADIUS = SQRT( (SHIFT_X)**2 + (SHIFT_Z)**2 )
        
                       IF (0.0_EB.LE.RADIUS.AND.RADIUS.LT.0.2_EB) THEN
                          VNU = 5.0_EB * RADIUS
                       ELSE IF (RADIUS.LT.0.4_EB) THEN
                         VNU = 2.0_EB - 5.0_EB * RADIUS
                       ELSE IF (RADIUS.GE.0.4_EB) THEN
                          VNU = 0.0_EB
                       ENDIF
                       U(I,J,K) = U_MAX - VNU * SHIFT_Z / RADIUS
                    ENDDO
                 ENDDO
              ENDDO
              DO K=0,KBAR
                 DO J=1,JBAR
                    DO I=1,IBAR
                       SHIFT_X= XC(I) - 0.5_EB
                       SHIFT_Y= 0.0_EB
                       SHIFT_Z= Z(K) - 0.5_EB
                       RADIUS = SQRT( (SHIFT_X)**2 + (SHIFT_Z)**2 )

                       IF (0.0_EB.LE.RADIUS.AND.RADIUS.LT.0.2_EB) THEN
                          VNU=5.0_EB*RADIUS
                       ELSE IF (RADIUS.LT.0.4_EB) THEN
                          VNU=2.0_EB-5.0_EB*RADIUS
                       ELSE IF (RADIUS.GE.0.4_EB) THEN
                          VNU=0.0_EB
                       ENDIF
                       W(I,J,K) = VNU * SHIFT_X / RADIUS
                    ENDDO
                 ENDDO
              ENDDO
              DO K=0,KBP1
                DO J=0,JBP1
                   DO I=0,IBP1
                     RHO(I,J,K) = 1.0_EB
                  ENDDO
               ENDDO
             ENDDO


!-------------------------------------------------------------------------------------
! ZM_GRAV_ADVECTED_2D
! Bitte alles noch mal nachprüfen !!!!!!!!!!!!!
!-------------------------------------------------------------------------------------
   CASE(3)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               U(I,J,K) = 1.00_EB 
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               W(I,J,K) = 1.25_EB
            ENDDO
         ENDDO
      ENDDO
      DO K=0,KBP1
         DO J=0,JBP1
            DO I=0,IBP1
               RADIUS=SQRT((XC(I)-0.75_EB)**2+(ZC(K)-0.75_EB)**2)
               IF (0.0_EB.LE.RADIUS.AND.RADIUS.LT.0.5_EB) THEN
                  RHO_R=1.0_EB+0.05_EB*(SIN(RADIUS/0.5_EB*PI+0.5_EB*PI)+1.0_EB)
               ELSE IF (0.5_EB.LE.RADIUS) THEN
                  RHO_R=1.0_EB
               ENDIF
               RHO(I,J,K) = RHO_R
            ENDDO
         ENDDO
      ENDDO
   END SELECT


END SUBROUTINE SCARC_ANALYTICAL_SOLUTION


 
SUBROUTINE SCARC_DUMP_QUANTITIES(NM,T,ITYPE)
IMPLICIT NONE
! Initialize flow variables with an analytical solution of the governing equations

INTEGER, INTENT(IN) :: NM, ITYPE
REAL(EB), INTENT(IN) :: T
REAL(EB), POINTER, DIMENSION(:,:,:) :: U, V, W, RHO, H
REAL(EB), POINTER, DIMENSION(:) :: RDXN, RDYN, RDZN, XC, YC, ZC
INTEGER, POINTER, DIMENSION(:,:,:) :: PRESSURE_ZONE
REAL(EB)::  VORTY(1000)
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR
REAL(EB), POINTER, DIMENSION(:) :: X, Y, Z, P_0, RHO_AVG !! ACHTUNG; FALSCH !!!
!REAL(EB), POINTER :: RHO_AVG
REAL(EB):: UU, VV, WW, VEL2, PRES
INTEGER, POINTER:: IIBAR, JJBAR, KKBAR
INTEGER:: II, JJ, KK, III, JJJ, KKK, MOUT, ILEN, IZERO=0
CHARACTER(40):: COUT

TYPE (MESH_TYPE), POINTER :: M


M=>MESHES(NM)

U=>M%U
V=>M%V
W=>M%W

H=>M%H

PBAR=>M%PBAR
RHO=>M%RHO
RHO_AVG=0.0_EB

P_0=>M%P_0

X=>M%X
Y=>M%Y
Z=>M%Z

XC=>M%XC
YC=>M%YC
ZC=>M%ZC
   
IIBAR=>M%IBAR
JJBAR=>M%JBAR
KKBAR=>M%KBAR

RDXN=>M%RDXN
RDYN=>M%RDYN
RDZN=>M%RDZN
   

PRESSURE_ZONE=>M%PRESSURE_ZONE

MOUT=99

! ------------------------------------------------------------------------
! print out current time
! ------------------------------------------------------------------------
IF (ITYPE==1) THEN
   COUT='dump/  _time_first'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE IF (ITYPE==2) THEN
   COUT='dump/  _time_second'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE
   COUT='dump/  _time_final'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ENDIF
OPEN(MOUT,file=COUT)
WRITE(MOUT,*) T
CLOSE(MOUT)


! ------------------------------------------------------------------------
! print out pressure
! ------------------------------------------------------------------------
IF (ITYPE==1) THEN
   COUT='dump/  _pressure_first'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE IF (ITYPE==2) THEN
   COUT='dump/  _pressure_second'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE
   COUT='dump/  _pressure_final'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ENDIF
OPEN(MOUT,file=COUT)
DO KK=1,KKBAR
   DO JJ=1,JJBAR
      DO II=1,IIBAR
         UU   = 0.5_EB*(U(MIN(IIBAR,II),JJ,KK)+U(MAX(0,II-1),JJ,KK))
         VV   = 0.5_EB*(V(II,MIN(JJBAR,JJ),KK)+V(II,MAX(0,JJ-1),KK))
         WW   = 0.5_EB*(W(II,JJ,MIN(KKBAR,KK))+W(II,JJ,MAX(0,KK-1)))
         VEL2 = UU**2+VV**2+WW**2
       !  PRES = PBAR(KK,PRESSURE_ZONE(II,JJ,KK)) + RHO_AVG*(H(II,JJ,KK)-.5_EB*VEL2) - P_0(KK)
         WRITE(MOUT,*) PRES
      ENDDO
   ENDDO
ENDDO
CLOSE(MOUT)


! ------------------------------------------------------------------------
! print out h
! ------------------------------------------------------------------------
IF (ITYPE==1) THEN
   COUT='dump/  _h_first'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE IF (ITYPE==2) THEN
   COUT='dump/  _h_second'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE
   COUT='dump/  _h_final'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ENDIF
OPEN(MOUT,file=COUT)
DO KK=1,KKBAR
   DO JJ=1,JJBAR
      DO II=1,IIBAR
         WRITE(MOUT,*) H(II,JJ,KK)
      ENDDO
   ENDDO
ENDDO
CLOSE(MOUT)


! ------------------------------------------------------------------------
! print out u-velocity
! ------------------------------------------------------------------------
IF (ITYPE==1) THEN
   COUT='dump/  _uvelocity_first'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE IF (ITYPE==2) THEN
   COUT='dump/  _uvelocity_second'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE
   COUT='dump/  _uvelocity_final'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ENDIF
OPEN(MOUT,file=COUT)
DO KK=1,KKBAR
   DO JJ=1,JJBAR
      DO II=0,IIBAR
         WRITE(MOUT,*) U(II,JJ,KK)
      ENDDO
   ENDDO
ENDDO
CLOSE(MOUT)


! ------------------------------------------------------------------------
! print out v-velocity
! ------------------------------------------------------------------------
IF (ITYPE==1) THEN
   COUT='dump/  _vvelocity_first'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE IF (ITYPE==2) THEN
   COUT='dump/  _vvelocity_second'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE
   COUT='dump/  _vvelocity_final'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ENDIF
OPEN(MOUT,file=COUT)
DO KK=1,KKBAR
   DO JJ=0,JJBAR
      DO II=1,IIBAR
         WRITE(MOUT,*) V(II,JJ,KK)
      ENDDO
   ENDDO
ENDDO
CLOSE(MOUT)


! ------------------------------------------------------------------------
! print out w-velocity
! ------------------------------------------------------------------------
IF (ITYPE==1) THEN
   COUT='dump/  _wvelocity_first'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE IF (ITYPE==2) THEN
   COUT='dump/  _wvelocity_second'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE
   COUT='dump/  _wvelocity_final'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ENDIF
OPEN(MOUT,file=COUT)
DO KK=0,KKBAR
   DO JJ=1,JJBAR
      DO II=1,IIBAR
         WRITE(MOUT,*) W(II,JJ,KK)
      ENDDO
   ENDDO
ENDDO
CLOSE(MOUT)

! ------------------------------------------------------------------------
! print out y-vorticity
! ------------------------------------------------------------------------
IF (ITYPE==1) THEN
   COUT='dump/  _yvorticity_first'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE IF (ITYPE==2) THEN
   COUT='dump/  _yvorticity_second'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE
   COUT='dump/  _yvorticity_final'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ENDIF
OPEN(MOUT,file=COUT)
JJ=1
DO KK=1,KKBAR
   DO II=1,IIBAR
      JJJ = MAX(1,MIN(JJ,JBAR))
      VORTY(II) = (U(II,JJJ,KK+1)-U(II,JJJ,KK))*RDZN(KK) - (W(II+1,JJJ,KK)-W(II,JJJ,KK))*RDXN(II)
      WRITE(MOUT,*) VORTY(II)
   ENDDO
ENDDO
CLOSE(MOUT)


! ------------------------------------------------------------------------
! print out density
! ------------------------------------------------------------------------
IF (ITYPE==1) THEN
   COUT='dump/  _density_first'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE IF (ITYPE==2) THEN
   COUT='dump/  _density_second'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ELSE
   COUT='dump/  _density_final'
   WRITE(COUT(6:7)  ,'(I2.2)') NM
ENDIF
OPEN(MOUT,file=COUT)
DO KK=1,KKBAR
   DO JJ=1,JJBAR
      DO II=1,IIBAR
         WRITE(MOUT,*) RHO(II,JJ,KK)
      ENDDO
   ENDDO
ENDDO
CLOSE(MOUT)


END SUBROUTINE SCARC_DUMP_QUANTITIES


 
SUBROUTINE SCARC_DUMP_VORTICTY(T,NM,ICYC)
IMPLICIT NONE
! Initialize flow variables with an analytical solution of the governing equations

INTEGER, INTENT(IN) :: NM, ICYC
LOGICAL:: BDUMP_DIAG, BDUMP_ALL
REAL(EB), INTENT(IN) :: T
REAL(EB):: VORTX(1000), VORTY(1000), VORTZ(1000), EPS
REAL(EB), POINTER, DIMENSION(:,:,:) :: U, V, W
REAL(EB), POINTER, DIMENSION(:) :: RDXN, RDYN, RDZN, X, Y, Z, XC, YC, ZC
INTEGER, POINTER:: IIBAR, JJBAR, KKBAR
INTEGER:: II, JJ, KK, III, JJJ, KKK, MVORTX, MVORTY, MVORTZ, IZERO=0
LOGICAL:: BVORT=.TRUE.
CHARACTER(40):: CVORTX, CVORTY, CVORTZ

TYPE (MESH_TYPE), POINTER :: M

!write(*,*) 'DUMP_VORT: ICYC=',ICYC


BDUMP_DIAG=.FALSE.
IF (NMESHES==1) THEN
   BDUMP_DIAG=.TRUE.
ELSE IF (NMESHES==4) THEN
   IF (NM== 1) BDUMP_DIAG=.TRUE.
   IF (NM== 4) BDUMP_DIAG=.TRUE.
ELSE IF (NMESHES==9) THEN
   IF (NM== 1) BDUMP_DIAG=.TRUE.
   IF (NM== 5) BDUMP_DIAG=.TRUE.
   IF (NM== 9) BDUMP_DIAG=.TRUE.
ELSE IF (NMESHES==16) THEN
   IF (NM== 1) BDUMP_DIAG=.TRUE.
   IF (NM== 6) BDUMP_DIAG=.TRUE.
   IF (NM==11) BDUMP_DIAG=.TRUE.
   IF (NM==16) BDUMP_DIAG=.TRUE.
ELSE IF (NMESHES==25) THEN
   IF (NM== 1) BDUMP_DIAG=.TRUE.
   IF (NM== 7) BDUMP_DIAG=.TRUE.
   IF (NM==13) BDUMP_DIAG=.TRUE.
   IF (NM==19) BDUMP_DIAG=.TRUE.
   IF (NM==25) BDUMP_DIAG=.TRUE.
ELSE IF (NMESHES==64) THEN
   IF (NM== 1) BDUMP_DIAG=.TRUE.
   IF (NM==10) BDUMP_DIAG=.TRUE.
   IF (NM==19) BDUMP_DIAG=.TRUE.
   IF (NM==28) BDUMP_DIAG=.TRUE.
   IF (NM==37) BDUMP_DIAG=.TRUE.
   IF (NM==46) BDUMP_DIAG=.TRUE.
   IF (NM==55) BDUMP_DIAG=.TRUE.
   IF (NM==64) BDUMP_DIAG=.TRUE.
ENDIF


BDUMP_DIAG=.FALSE.
IF (BDUMP_DIAG) THEN

   M=>MESHES(NM)
   
   U=>M%U
   V=>M%V
   W=>M%W
   
   X=>M%X
   Y=>M%Y
   Z=>M%Z
   
   XC=>M%XC
   YC=>M%YC
   ZC=>M%ZC
   
   RDXN=>M%RDXN
   RDYN=>M%RDYN
   RDZN=>M%RDZN
   
   IIBAR=>M%IBAR
   JJBAR=>M%JBAR
   KKBAR=>M%KBAR
   
   
   !MVORTX=97
   !CVORTX='dump/  _vortx_cyc    '
   !WRITE(CVORTX(6:7)  ,'(I2.2)') NM
   !WRITE(CVORTX(18:21),'(I4.4)') ICYC
   !OPEN(MVORTX,file=CVORTX)
   
   MVORTY=98
   CVORTY='dump/  _vorty_diag_cyc    '
   WRITE(CVORTY(6:7)  ,'(I2.2)') NM
   WRITE(CVORTY(23:26),'(I4.4)') ICYC
   OPEN(MVORTY,file=CVORTY)
   WRITE(MVORTY,*) 'Time=',T,': Cycle=',ICYC
   
   !MVORTZ=99
   !CVORTZ='dump/  _vortz_cyc    '
   !WRITE(CVORTZ(6:7)  ,'(I2.2)') NM
   !WRITE(CVORTZ(18:21),'(I4.4)') ICYC
   !OPEN(MVORTZ,file=CVORTZ)
   
   
   EPS=1.0E-8_EB
   
   JJ=1
   
   DO KK=1,KKBAR
      DO II=1,IIBAR
   

         IF (ABS(X(II)-Z(KK)).LT.EPS) THEN
   
            ! VORTICITY X
            !III = MAX(1,MIN(II,IBAR))
            !VORTX(II) = (W(III,JJ+1,KK)-W(III,JJ,KK))*RDYN(JJ) - (V(III,JJ,KK+1)-V(III,JJ,KK))*RDZN(KK)
            !WRITE(MVORTX,*) VORTX(II)
   
            ! VORTICITY Z
            JJJ = MAX(1,MIN(JJ,JBAR))
            VORTY(II) = (U(II,JJJ,KK+1)-U(II,JJJ,KK))*RDZN(KK) - (W(II+1,JJJ,KK)-W(II,JJJ,KK))*RDXN(II)
            WRITE(MVORTY,*) xC(II),',',VORTY(II)
            WRITE(*,*) xC(II),',',VORTY(II)
   
            ! VORTICITY Z
            !KKK = MAX(1,MIN(KK,KBAR))
            !VORTZ(II) = (V(II+1,JJ,KKK)-V(II,JJ,KKK))*RDXN(II) - (U(II,JJ+1,KKK)-U(II,JJ,KKK))*RDYN(JJ)
            !WRITE(MVORTZ,*) VORTZ(II)
   
         ENDIF
   
      ENDDO
   ENDDO
   
   CLOSE(MVORTY)

ENDIF

BDUMP_ALL=.TRUE.
IF (BDUMP_ALL) THEN

   MVORTY=198
   CVORTY='dump/  _vorty_all_cyc    '
   WRITE(CVORTY(6:7)  ,'(I2.2)') NM
   WRITE(CVORTY(22:25),'(I4.4)') ICYC
   OPEN(MVORTY,file=CVORTY)
   WRITE(MVORTY,*) 'Time=',T,': Cycle=',ICYC
   
write(*,*) 'SCARC-DUMP:1', MVORTY

   !MVORTZ=99
   M=>MESHES(NM)
   
   U=>M%U
   V=>M%V
   W=>M%W
   
   X=>M%X
   Y=>M%Y
   Z=>M%Z
   
   XC=>M%XC
   YC=>M%YC
   ZC=>M%ZC
   
   RDXN=>M%RDXN
   RDYN=>M%RDYN
   RDZN=>M%RDZN
   
   IIBAR=>M%IBAR
   JJBAR=>M%JBAR
   KKBAR=>M%KBAR

   JJ=1
   DO KK=1,KKBAR
      DO II=1,IIBAR
         JJJ = MAX(1,MIN(JJ,JBAR))
         VORTY(II) = (U(II,JJJ,KK+1)-U(II,JJJ,KK))*RDZN(KK) - (W(II+1,JJJ,KK)-W(II,JJJ,KK))*RDXN(II)
         WRITE(MVORTY,*) xC(II),ZC(KK),VORTY(II)
         WRITE(*,*) xC(II),',',VORTY(II)
      ENDDO
   ENDDO
   
   CLOSE(MVORTY)

ENDIF
   
   
END SUBROUTINE SCARC_DUMP_VORTICTY






 
 
END MODULE SCARC_SOLVER
 
 
