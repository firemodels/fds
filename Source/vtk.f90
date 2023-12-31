!> \brief Routines for handling vtk output

MODULE VTK

USE MESH_VARIABLES
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE OUTPUT_DATA

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC INITIALIZE_VTK,WRITE_VTK_SLICE_GEOMETRY,WRITE_VTK_SLICE_CELLS,WRITE_VTK_SLICE_DATA,FINALIZE_VTK,&
       WRITE_VTK_SLICE_WRAPPER,BUILD_VTK_GAS_PHASE_GEOMETRY,WRITE_VTK_SM3D_WRAPPER

CONTAINS

SUBROUTINE INITIALIZE_VTK(LU,FN,MESH_TOPOLOGY)

INTEGER, INTENT(IN) :: LU
CHARACTER(LEN=*), INTENT(IN) :: MESH_TOPOLOGY,FN
OPEN(LU,FILE=FN,FORM='FORMATTED')
WRITE(LU,'(A)',ADVANCE='YES') '<?xml version="1.0"?>'
WRITE(LU,'(A,A,A)',ADVANCE='YES') '<VTKFile type="',TRIM(MESH_TOPOLOGY),'" version="1.0" byte_order="LittleEndian">'
WRITE(LU,'(A,A,A)',ADVANCE='YES') '  <',TRIM(MESH_TOPOLOGY),'>'
CLOSE(LU)
END SUBROUTINE INITIALIZE_VTK

SUBROUTINE WRITE_VTK_SLICE_GEOMETRY(LU, FN, NP, NC, X_PTS, Y_PTS, Z_PTS, FORMAT)

INTEGER, INTENT(IN) :: LU, NC, NP
REAL(FB), DIMENSION(:), INTENT(IN) :: X_PTS, Y_PTS, Z_PTS
CHARACTER(LEN=*), INTENT(IN) :: FN, FORMAT
INTEGER :: I

OPEN(LU,FILE=FN,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
WRITE(LU,'(A,I0,A,I0,A)',ADVANCE='YES') '    <Piece NumberOfPoints="',NP,'" NumberOfCells="',NC,'">'
WRITE(LU,'(A)',ADVANCE='YES') '      <Points>'
WRITE(LU,'(A,A,A)',ADVANCE='YES') &
   '        <DataArray type="Float32" NumberOfComponents="3" Name="Points" format="',FORMAT,'">'

WRITE(LU,'(A)',ADVANCE='NO') '          '
IF (TRIM(ADJUSTL(FORMAT)) .eq. 'ascii') THEN
   DO I=1,NP
      WRITE(LU,'(E15.8,A,E15.8,A,E15.8,A)',ADVANCE='NO') X_PTS(I),' ',Y_PTS(I),' ',Z_PTS(I),' '
   ENDDO
ENDIF
WRITE(LU,'(A)',ADVANCE='YES') ' '
WRITE(LU,'(A)',ADVANCE='YES') '        </DataArray>'
WRITE(LU,'(A)',ADVANCE='YES') '      </Points>'
CLOSE(LU)

END SUBROUTINE WRITE_VTK_SLICE_GEOMETRY

SUBROUTINE WRITE_VTK_SLICE_CELLS(LU, FN, CONNECT, OFFSETS, CELL_TYPE, FORMAT)

INTEGER, INTENT(IN) :: LU
INTEGER, DIMENSION(:), INTENT(IN) :: CONNECT, OFFSETS
INTEGER(IB4), DIMENSION(:), INTENT(IN) :: CELL_TYPE
CHARACTER(LEN=*), INTENT(IN) :: FN, FORMAT
INTEGER :: I

! Open cells section
OPEN(LU,FILE=FN,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
WRITE(LU,'(A)',ADVANCE='YES') '      <Cells>'

! Write connectivity
WRITE(LU,'(A,A,A)',ADVANCE='YES') &
   '        <DataArray type="Int32" NumberOfComponents="1" Name="connectivity" format="',FORMAT,'">'

WRITE(LU,'(A)',ADVANCE='NO') '          '
IF (TRIM(ADJUSTL(FORMAT)) .EQ. 'ascii') THEN
   DO I=1,SIZE(CONNECT)
      WRITE(LU,'(I0,A)',ADVANCE='NO') CONNECT(I), ' '
   ENDDO
ENDIF
WRITE(LU,'(A)',ADVANCE='YES') ' '
WRITE(LU,'(A)',ADVANCE='YES') '        </DataArray>'

! Write offsets
WRITE(LU,'(A,A,A)',ADVANCE='YES') &
   '        <DataArray type="Int32" NumberOfComponents="1" Name="offsets" format="',FORMAT,'">'

WRITE(LU,'(A)',ADVANCE='NO') '          '
IF (TRIM(ADJUSTL(FORMAT)) .EQ. 'ascii') THEN
   DO I=1,SIZE(OFFSETS)
      WRITE(LU,'(I0,A)',ADVANCE='NO') OFFSETS(I), ' '
   ENDDO
ENDIF
WRITE(LU,'(A)',ADVANCE='YES') ' '
WRITE(LU,'(A)',ADVANCE='YES') '        </DataArray>'

! Write cell types
WRITE(LU,'(A,A,A)',ADVANCE='YES') &
   '        <DataArray type="Int8" NumberOfComponents="1" Name="types" format="',FORMAT,'">'

WRITE(LU,'(A)',ADVANCE='NO') '          '
IF (TRIM(ADJUSTL(FORMAT)) .EQ. 'ascii') THEN
   DO I=1,SIZE(CELL_TYPE)
      WRITE(LU,'(I0,A)',ADVANCE='NO') CELL_TYPE(I), ' '
   ENDDO
ENDIF
WRITE(LU,'(A)',ADVANCE='YES') ' '
WRITE(LU,'(A)',ADVANCE='YES') '        </DataArray>'

! Close cell section
WRITE(LU,'(A)',ADVANCE='YES') '      </Cells>'

! Open point data section
WRITE(LU,'(A)',ADVANCE='YES') '      <PointData>'

CLOSE(LU)

END SUBROUTINE WRITE_VTK_SLICE_CELLS

SUBROUTINE WRITE_VTK_SLICE_DATA(LU, FN, DATA, DATA_NAME, FORMAT)

INTEGER, INTENT(IN) :: LU
CHARACTER(LEN=*), INTENT(IN) :: FORMAT, DATA_NAME, FN
CHARACTER(LEN=:), ALLOCATABLE :: CODE
REAL(FB), DIMENSION(:) :: DATA
INTEGER :: I, NBITS_DATA, NSTRINGS_DATA
CHARACTER(LEN=:), ALLOCATABLE :: BIT_DATA

!real(R16P),                    intent(in)  :: n       !< Number to be encoded.
integer(I1P),     allocatable              :: nI1P(:) !< One byte integer array containing n.
integer(I4P)                               :: padd    !< Number of padding characters ('=').

integer(I4P) :: BYR16P = 2_I1P




! Open data section
OPEN(LU,FILE=FN,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')

! Write data
WRITE(LU,'(A,A,A,A,A)',ADVANCE='YES') &
   '        <DataArray type="Float32" NumberOfComponents="1" Name="',&
  TRIM(DATA_NAME),'" format="',FORMAT,'">'

WRITE(LU,'(A)',ADVANCE='NO') '          '
IF (TRIM(ADJUSTL(FORMAT)) .EQ. 'ascii') THEN
   DO I=1,SIZE(DATA)
      WRITE(LU,'(E15.8,A)',ADVANCE='NO') DATA(I), ' '
   ENDDO
ELSEIF (TRIM(ADJUSTL(FORMAT)) .EQ. 'binary') THEN
   CLOSE(LU)
   OPEN(LU,FILE=FN,FORM='UNFORMATTED',STATUS='OLD',POSITION='APPEND')
   !NBITS_DATA = SIZE(DATA) * SIZEOF(FB)*2
   !NSTRINGS_DATA = INT(FLOAT(NBITS_DATA)/3+0.9)
   
   allocate(nI1P(1:((BYR16P+2)/3)*3)) ; nI1P = 0_I1P
   code = repeat(' ',((BYR16P+2)/3)*4)
   nI1P = transfer(DATA,nI1P)
   padd = mod((BYR16P),3) ; if (padd>0_I4P) padd = 3_I4P - padd


   
   CALL ENCODE_BITS(BITS=NI1P,PADD=PADD,CODE=CODE)
   
   WRITE(*,*) "SIZES", SIZE(DATA), NBITS_DATA
   WRITE(*,*) PADD
   WRITE(*,*) CODE
   !ALLOCATE(CHARACTER(LEN=NBITS_DATA) :: BIT_DATA)
   !WRITE(BIT_DATA,'(B64)') DATA
   !WRITE(*,*) BIT_DATA
   !DEALLOCATE(BIT_DATA)
   !DO I=0,NSTRINGS_DATA
   !   STRINGBITS = DATA(I*6:(I+1)*6)
   !   WRITE(*,*), "DATA ", I, " ", STRINGBITS
   !ENDDO
   !WRITE(LU_SL3D_VTK(NM)) NBITS_DATA  , (DATA(I),I=1,SIZE(DATA))
   CLOSE(LU)
   OPEN(LU,FILE=FN,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
   !CODE=ENCODE_ASCII_DATAARRAY1_RANK1_FB(DATA)
   !WRITE(LU_SL3D_VTK(NM)) CODE
ELSEIF (TRIM(ADJUSTL(FORMAT)) .EQ. 'appended') THEN
   CLOSE(LU)
   OPEN(LU,FILE=FN,FORM='UNFORMATTED',STATUS='OLD',POSITION='APPEND')
   NBITS_DATA = SIZE(DATA) * SIZEOF(FB)
   WRITE(LU) NBITS_DATA  , (DATA(I),I=1,SIZE(DATA))
   CLOSE(LU)
   OPEN(LU,FILE=FN,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
   !CODE=ENCODE_ASCII_DATAARRAY1_RANK1_FB(DATA)
   !WRITE(LU_SL3D_VTK(NM)) CODE
ENDIF
WRITE(LU,'(A)',ADVANCE='YES') ' '

! Close data section
WRITE(LU,'(A)',ADVANCE='YES') '        </DataArray>'
CLOSE(LU)

END SUBROUTINE WRITE_VTK_SLICE_DATA

SUBROUTINE FINALIZE_VTK(LU,FN,MESH_TOPOLOGY)

INTEGER, INTENT(IN) :: LU
CHARACTER(LEN=*), INTENT(IN) :: MESH_TOPOLOGY,FN

! Close VTK file
OPEN(LU,FILE=FN,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
WRITE(LU,'(A)',ADVANCE='YES') '      </PointData>'
WRITE(LU,'(A)',ADVANCE='YES') '    </Piece>'
WRITE(LU,'(A,A,A)',ADVANCE='YES') '  </',TRIM(MESH_TOPOLOGY),'>'
WRITE(LU,'(A)',ADVANCE='YES') '</VTKFile>'

CLOSE(LU)
END SUBROUTINE FINALIZE_VTK

SUBROUTINE WRITE_VTK_SLICE_WRAPPER(T,NMESHES,MESH_TOPOLOGY)
INTEGER, INTENT(IN) :: NMESHES
CHARACTER(LEN=*), INTENT(IN) :: MESH_TOPOLOGY
REAL(EB), INTENT(IN) :: T
TYPE (SLICE_TYPE), POINTER :: SL
INTEGER :: I,NQT,I1,I2,J1,J2,K1,K2,IQ,ITM,ITM1
REAL(FB) :: STIME
REAL(EB) :: TT
CHARACTER(200) :: TMP_FILE

! Generate filename
STIME = REAL(T_BEGIN + (T-T_BEGIN)*TIME_SHRINK_FACTOR,FB)
TT   = T_BEGIN + (T-T_BEGIN)*TIME_SHRINK_FACTOR
ITM  = INT(TT)
ITM1 = NINT(ABS(TT-ITM)*100)
IF (ITM1==100) THEN
   ITM = ITM+1
   ITM1 = 0
ENDIF
WRITE(FN_SL3D_VTK(NMESHES+1),'(A,A,A,I8.8,I2.2,A)') "./results/",TRIM(CHID),'_SL3D_',ITM,ITM1,'.pvtu'

! First part of header before quantities
OPEN(LU_SL3D_VTK(NMESHES+1),FILE=FN_SL3D_VTK(NMESHES+1),FORM='FORMATTED')
WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '<?xml version="1.0"?>'
WRITE(LU_SL3D_VTK(NMESHES+1),'(A,A,A)',ADVANCE='YES') &
   '<VTKFile type="',TRIM(MESH_TOPOLOGY),'" version="1.0" byte_order="LittleEndian">'
WRITE(LU_SL3D_VTK(NMESHES+1),'(A,A,A)',ADVANCE='YES') '  <',TRIM(MESH_TOPOLOGY),' GhostLevel="0">'
WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '    <PPointData>'

! Add PointData arrays
NQT = MESHES(1)%N_SLCF
QUANTITY_LOOP2: DO IQ=1,NQT
   SL => SLICE(IQ)
   I1  = SL%I1
   I2  = SL%I2
   J1  = SL%J1
   J2  = SL%J2
   K1  = SL%K1
   K2  = SL%K2
   IF (I2-I1==0 .OR. J2-J1==0 .OR. K2-K1==0) CYCLE QUANTITY_LOOP2
   WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='NO') '      <PDataArray type="Float32" NumberOfComponents="1" Name='
   WRITE(LU_SL3D_VTK(NMESHES+1),'(A,A,A)',ADVANCE='NO') '"',TRIM(SL%SMOKEVIEW_LABEL(1:30)),'" format='
   IF (VTK_BINARY) THEN
      WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '"binary"/>'
   ELSE
      WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '"ascii"/>'
   ENDIF
ENDDO QUANTITY_LOOP2

! Close PointData
WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '    </PPointData>'

! Add Point arrays
WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '    <PPoints>'
WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='NO') '      <PDataArray type="Float32" NumberOfComponents="3" Name="Points" format='
IF (VTK_BINARY) THEN
   WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '"binary"/>'
ELSE
   WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '"ascii"/>'
ENDIF
WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '    </PPoints>'

! Add CellData
WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '    <PCellData></PCellData>'

! Write pieces
DO I=1,NMESHES
   WRITE(TMP_FILE,'(A,A,I0,A,I8.8,I2.2,A)') TRIM(CHID),'_SL3D_',I,'_',ITM,ITM1,'.vtu'
   WRITE(LU_SL3D_VTK(NMESHES+1),'(A,A,A)',ADVANCE='YES') '    <Piece Source="',TRIM(TMP_FILE),'"/>'
ENDDO

! Finalize file
WRITE(LU_SL3D_VTK(NMESHES+1),'(A,A,A)',ADVANCE='YES') '  </',TRIM(MESH_TOPOLOGY),'>'
WRITE(LU_SL3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '</VTKFile>'
CLOSE(LU_SL3D_VTK(NMESHES+1))
END SUBROUTINE WRITE_VTK_SLICE_WRAPPER

















SUBROUTINE WRITE_VTK_SM3D_WRAPPER(T,NMESHES,MESH_TOPOLOGY)
INTEGER, INTENT(IN) :: NMESHES
CHARACTER(LEN=*), INTENT(IN) :: MESH_TOPOLOGY
REAL(EB), INTENT(IN) :: T
TYPE (SLICE_TYPE), POINTER :: SL
INTEGER :: I,NQT,I1,I2,J1,J2,K1,K2,IQ,ITM,ITM1,N
REAL(FB) :: STIME
REAL(EB) :: TT
CHARACTER(200) :: TMP_FILE
TYPE(SMOKE3D_TYPE), POINTER :: S3

! Generate filename
STIME = REAL(T_BEGIN + (T-T_BEGIN)*TIME_SHRINK_FACTOR,FB)
TT   = T_BEGIN + (T-T_BEGIN)*TIME_SHRINK_FACTOR
ITM  = INT(TT)
ITM1 = NINT(ABS(TT-ITM)*100)
IF (ITM1==100) THEN
   ITM = ITM+1
   ITM1 = 0
ENDIF
WRITE(FN_SMOKE3D_VTK(NMESHES+1),'(A,A,A,A,A,I8.8,I2.2,A)') "./results/",TRIM(CHID),'_','SM3D','_',ITM,ITM1,'.pvtu'

! First part of header before quantities
OPEN(LU_SMOKE3D_VTK(NMESHES+1),FILE=FN_SMOKE3D_VTK(NMESHES+1),FORM='FORMATTED')
WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '<?xml version="1.0"?>'
WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A,A,A)',ADVANCE='YES') &
   '<VTKFile type="',TRIM(MESH_TOPOLOGY),'" version="1.0" byte_order="LittleEndian">'
WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A,A,A)',ADVANCE='YES') '  <',TRIM(MESH_TOPOLOGY),' GhostLevel="0">'
WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '    <PPointData>'

! Add PointData arrays
DATA_FILE_LOOP2: DO N=1,N_SMOKE3D

   S3 => SMOKE3D_FILE(N)
   WRITE(*,*) "SMOKE3D QUANTITY", S3%SMOKEVIEW_LABEL(1:30)
   WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='NO') '      <PDataArray type="Float32" NumberOfComponents="1" Name='
   WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A,A,A)',ADVANCE='NO') '"',TRIM(S3%SMOKEVIEW_LABEL(1:30)),'" format='
   IF (VTK_BINARY) THEN
      WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '"binary"/>'
   ELSE
      WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '"ascii"/>'
   ENDIF
ENDDO DATA_FILE_LOOP2

! Close PointData
WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '    </PPointData>'

! Add Point arrays
WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '    <PPoints>'
WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='NO') '      <PDataArray type="Float32" NumberOfComponents="3" Name="Points" format='
IF (VTK_BINARY) THEN
   WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '"binary"/>'
ELSE
   WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '"ascii"/>'
ENDIF
WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '    </PPoints>'

! Add CellData
WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '    <PCellData></PCellData>'

! Write pieces
DO I=1,NMESHES
   WRITE(TMP_FILE,'(A,A,I0,A,I8.8,I2.2,A)') TRIM(CHID),'_SM3D_',I,'_',ITM,ITM1,'.vtu'
   WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A,A,A)',ADVANCE='YES') '    <Piece Source="',TRIM(TMP_FILE),'"/>'
ENDDO

! Finalize file
WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A,A,A)',ADVANCE='YES') '  </',TRIM(MESH_TOPOLOGY),'>'
WRITE(LU_SMOKE3D_VTK(NMESHES+1),'(A)',ADVANCE='YES') '</VTKFile>'
CLOSE(LU_SMOKE3D_VTK(NMESHES+1))
END SUBROUTINE WRITE_VTK_SM3D_WRAPPER



SUBROUTINE BUILD_VTK_GAS_PHASE_GEOMETRY(NM, MESH_TOPOLOGY, FN, LU)

INTEGER :: NX, NY, NZ, NC, NP, I, J, K, IFACT, JFACT, KFACT
REAL(FB), ALLOCATABLE, DIMENSION(:) :: QQ_PACK, X_PTS, Y_PTS, Z_PTS
INTEGER, ALLOCATABLE, DIMENSION(:) :: CONNECT, OFFSETS
INTEGER(IB4), ALLOCATABLE, DIMENSION(:) :: CELL_TYPE
CHARACTER(LEN=*), INTENT(IN) :: MESH_TOPOLOGY, FN
INTEGER, INTENT(IN) :: NM, LU

NX = SIZE(MESHES(NM)%X)
NY = SIZE(MESHES(NM)%Y)
NZ = SIZE(MESHES(NM)%Z)

NP = NX*NY*NZ
NC = (NX-1)*(NY-1)*(NZ-1)

! Fill point data
ALLOCATE(X_PTS(NP))
ALLOCATE(Y_PTS(NP))
ALLOCATE(Z_PTS(NP))
IFACT = 1
DO K = 0, NZ-1
   DO J = 0, NY-1
      DO I = 0, NX-1
         X_PTS(IFACT)=MESHES(NM)%X(I)
         Y_PTS(IFACT)=MESHES(NM)%Y(J)
         Z_PTS(IFACT)=MESHES(NM)%Z(K)
         IFACT = IFACT + 1
      ENDDO
   ENDDO
ENDDO

! Fill cell data
ALLOCATE(CONNECT(NC*8))
DO I = 1, NX-1
   IFACT = (I-1)
   DO J = 1, NY-1
      JFACT = (J-1)*(NX-1)
      DO K = 1, NZ-1
         KFACT = (K - 1)*(NY-1)*(NX-1)
         CONNECT((IFACT+JFACT+KFACT)*8+1) = (K-1)*(NY*NX) + (J-1)*NX + I-1
         CONNECT((IFACT+JFACT+KFACT)*8+2) = (K-1)*(NY*NX) + (J-1)*NX + I
         CONNECT((IFACT+JFACT+KFACT)*8+3) = (K-1)*(NY*NX) + (J)*NX + I-1
         CONNECT((IFACT+JFACT+KFACT)*8+4) = (K-1)*(NY*NX) + (J)*NX + I
         CONNECT((IFACT+JFACT+KFACT)*8+5) = (K)*(NY*NX) + (J-1)*NX + I-1
         CONNECT((IFACT+JFACT+KFACT)*8+6) = (K)*(NY*NX) + (J-1)*NX + I
         CONNECT((IFACT+JFACT+KFACT)*8+7) = (K)*(NY*NX) + (J)*NX + I-1
         CONNECT((IFACT+JFACT+KFACT)*8+8) = (K)*(NY*NX) + (J)*NX + I
      ENDDO
   ENDDO
ENDDO

ALLOCATE(OFFSETS(NC))
ALLOCATE(CELL_TYPE(NC))

DO I=1,NC
   OFFSETS(I) = (I)*8_IB4
   CELL_TYPE(I) = 11_IB4
ENDDO

IF (VTK_BINARY) THEN
   CALL WRITE_VTK_SLICE_GEOMETRY(LU,FN,NP,NC,X_PTS,Y_PTS,Z_PTS,'ascii')
   CALL WRITE_VTK_SLICE_CELLS(LU,FN,CONNECT,OFFSETS,CELL_TYPE,'ascii')
ELSE
   CALL WRITE_VTK_SLICE_GEOMETRY(LU,FN,NP,NC,X_PTS,Y_PTS,Z_PTS,'ascii')
   CALL WRITE_VTK_SLICE_CELLS(LU,FN,CONNECT,OFFSETS,CELL_TYPE,'ascii')
ENDIF

DEALLOCATE(OFFSETS)
DEALLOCATE(CELL_TYPE)
DEALLOCATE(CONNECT)
DEALLOCATE(X_PTS)
DEALLOCATE(Y_PTS)
DEALLOCATE(Z_PTS)
         
ENDSUBROUTINE BUILD_VTK_GAS_PHASE_GEOMETRY
         

SUBROUTINE ENCODE_BITS(BITS, PADD, CODE)
!< Encode a bits stream (must be multiple of 24 bits) into base64 charcaters code (of length multiple of 4).
!<
!< The bits stream are encoded in chunks of 24 bits as the following example (in little endian order)
!<```
!< +--first octet--+-second octet--+--third octet--+
!< |7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0|
!< +-----------+---+-------+-------+---+-----------+
!< |5 4 3 2 1 0|5 4 3 2 1 0|5 4 3 2 1 0|5 4 3 2 1 0|
!< +--1.index--+--2.index--+--3.index--+--4.index--+
!<```
!< @note The 4 indexes are stored into 4 elements 8 bits array, thus 2 bits of each array element are not used.
!<
!< @note The number of paddings must be computed outside this procedure, into the calling scope.
!<
!< @warning This procedure is the backend of encoding, thus it must be never called outside the module.
integer(I1P), intent(in)  :: bits(1:)  !< Bits to be encoded.
integer(I4P), intent(in)  :: padd      !< Number of padding characters ('=').
character(*), intent(out) :: code      !< Characters code.
integer(I1P)              :: sixb(1:4) !< 6 bits slices (stored into 8 bits integer) of 24 bits input.
integer(I8P)              :: c         !< Counter.
integer(I8P)              :: e         !< Counter.
integer(I8P)              :: Nb        !< Length of bits array.
character(64) :: base64="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/" !< Base64 alphabet.

Nb=size(bits,dim=1,kind=I8P)
c = 1_I8P
do e=1_I8P,Nb,3_I8P ! loop over array elements: 3 bytes (24 bits) scanning
  sixb = 0_I1P
     call mvbits(bits(e  ),2,6,sixb(1),0)
     call mvbits(bits(e  ),0,2,sixb(2),4)
  if (e+1<=Nb) then
     call mvbits(bits(e+1),4,4,sixb(2),0)
     call mvbits(bits(e+1),0,4,sixb(3),2)
  endif
  if (e+2<=Nb) then
     call mvbits(bits(e+2),6,2,sixb(3),0)
     call mvbits(bits(e+2),0,6,sixb(4),0)
  endif
  sixb = sixb + 1_I1P
  code(c  :c  ) = base64(sixb(1):sixb(1))
  code(c+1:c+1) = base64(sixb(2):sixb(2))
  code(c+2:c+2) = base64(sixb(3):sixb(3))
  code(c+3:c+3) = base64(sixb(4):sixb(4))
  c = c + 4_I8P
enddo
if (padd>0) code(len(code)-padd+1:)=repeat('=',padd)
endsubroutine encode_bits

END MODULE VTK

!FUNCTION ENCODE_ASCII_DATAARRAY1_RANK1_FB(X) RESULT(CODE)
!!< Encode (Base64) a dataarray with 1 components of rank 1 (FB).
!REAL(FB), INTENT(IN)          :: X(1:) !< Data variable.
!CHARACTER(len=:), ALLOCATABLE :: CODE  !< Encoded base64 dataarray.
!INTEGER(IB16)                  :: N     !< Counter.
!INTEGER(IB16)                  :: L     !< Length
!INTEGER(IB16)                  :: SP    !< String pointer
!INTEGER(IB16)                  :: SIZE_N!< Dimension size
!
!SIZE_N = SIZE(X,DIM=1)
!L = DI2P+1
!SP = 0
!CODE = REPEAT(' ',L*SIZE_N)
!DO N = 1,SIZE_N
!    CODE(SP+1:SP+L) = STR(N=X(N))
!    SP = SP + L
!ENDDO
!ENDFUNCTION ENCODE_ASCII_DATAARRAY1_RANK1_FB
!
!END MODULE VTK
