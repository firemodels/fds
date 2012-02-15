program main
USE PRECISION_PARAMETERS
USE ISOSMOKE
implicit none

integer :: LU_GEOM, LU_GEOM_DATA
integer :: ONE,VERSION
real :: STIME
integer :: N_VERT_S, N_FACE_S, N_VERT_D, N_FACE_D
integer :: N_VERT_S_VALS, N_FACE_S_VALS, N_VERT_D_VALS, N_FACE_D_VALS
real(fb), dimension(4) :: Xvert_S, Yvert_S, Zvert_S
real(fb), dimension(4) :: Xvert_D, Yvert_D, Zvert_D
real(fb), dimension(4) :: ValVertStatic, ValVertDynamic, ValFaceStatic, ValFaceDynamic
integer, dimension(4) :: FACE1_S, FACE2_S, FACE3_S
integer, dimension(4) :: FACE1_D, FACE2_D, FACE3_D
integer, dimension(4) :: SURF_S, SURF_D
integer :: GEOM_TYPE
integer :: I
real :: XTRAN, YTRAN, ZTRAN
real :: XTRAN0, YTRAN0, ZTRAN0
real :: ROT_AZ, ROT_ELEV
real :: XROT0, YROT0, ZROT0
real :: xoffset, yoffset, zoffset
integer :: itime
real :: twfin, dt, dx
integer :: nsteps
integer :: pos
real(fb), pointer, dimension(:) :: test
real(fb), dimension(0:7) :: vals,tvals
real(fb) :: level
real(fb), dimension(0:1) :: x=(/0.0,1.0/),y=(/0.0,1.0/),z=(/0.0,1.0/)
integer :: HAVE_TVALS=0
  INTEGER, DIMENSION(0:23) :: NODEINDEXES=(/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23/)
  REAL(FB), DIMENSION(0:35) :: XYZV_LOCAL,TV_LOCAL
  INTEGER, DIMENSION(0:14) :: TRIS
  INTEGER :: NXYZV
  INTEGER :: NTRIS
  INTEGER, DIMENSION(0:35) :: CLOSESTNODES
  INTEGER :: j

level=0.75

do i = 0,255
  do pos = 0, 7
    if(btest(i,pos))then
      vals(pos)=1
    else
      vals(pos)=0
    endif
  end do
!  SUBROUTINE FGETISOBOX(X,Y,Z,VALS,HAVE_TVALS,TVALS,NODEINDEXES,LEVEL,XYZV_LOCAL,TV_LOCAL,NXYZV,TRIS,NTRIS,CLOSESTNODES)
  call FGETISOBOX(x,y,z,vals,HAVE_TVALS,tvals,NODEINDEXES,level,XYZV_LOCAL,TV_LOCAL,NXYZV,TRIS,NTRIS,CLOSESTNODES)
  write(6,*)"case ",i
  write(6,*)"vals: ",(vals(j),j=0,7)
  do j = 0, nxyzv-1
    write(6,*)xyzv_local(3*j),xyzv_local(3*j+1),xyzv_local(3*j+2)
  end do
  do j = 0, ntris-1
    write(6,*)tris(3*j),tris(3*j+1),tris(3*j+2)
  end do
end do

call REALLOCATE_F(test,0,10)
do i = 1, 10
  test(i)=float(i)
end do
call REALLOCATE_F(test,10,20)
do i = 11,20
  test(i)=float(2*i)
end do
do i = 1, 20
  write(6,*)i,test(i)
end do

xoffset=0.0
yoffset=3.0
zoffset=3.0
ONE=1
VERSION=0
STIME=0.0
N_VERT_S=4
N_FACE_S=4
N_VERT_D=4
N_FACE_D=4
N_VERT_S_VALS=4
N_FACE_S_VALS=4
N_VERT_D_VALS=4
N_FACE_D_VALS=4
Xvert_S(1)=0.0
Yvert_S(1)=3.0
Zvert_S(1)=0.0
Xvert_S(2)=1.0
Yvert_S(2)=3.0
Zvert_S(2)=0.0
Xvert_S(3)=1.0
Yvert_S(3)=4.0
Zvert_S(3)=0.0
Xvert_S(4)=1.0
Yvert_S(4)=4.0
Zvert_S(4)=1.0
Xvert_D(1)=Xvert_S(1)+xoffset
Yvert_D(1)=Yvert_S(1)+yoffset
Zvert_D(1)=Zvert_S(1)+zoffset
Xvert_D(2)=Xvert_S(2)+xoffset
Yvert_D(2)=Yvert_S(2)+yoffset
Zvert_D(2)=Zvert_S(2)+zoffset
Xvert_D(3)=Xvert_S(3)+xoffset
Yvert_D(3)=Yvert_S(3)+yoffset
Zvert_D(3)=Zvert_S(3)+zoffset
Xvert_D(4)=Xvert_S(4)+xoffset
Yvert_D(4)=Yvert_S(4)+yoffset
Zvert_D(4)=Zvert_S(4)+zoffset
FACE1_S(1)=1
FACE2_S(1)=3
FACE3_S(1)=2
FACE1_S(2)=1
FACE2_S(2)=2
FACE3_S(2)=4
FACE1_S(3)=2
FACE2_S(3)=3
FACE3_S(3)=4
FACE1_S(4)=3
FACE2_S(4)=1
FACE3_S(4)=4
FACE1_D(1)=1
FACE2_D(1)=3
FACE3_D(1)=2
FACE1_D(2)=1
FACE2_D(2)=2
FACE3_D(2)=4
FACE1_D(3)=2
FACE2_D(3)=3
FACE3_D(3)=4
FACE1_D(4)=3
FACE2_D(4)=1
FACE3_D(4)=4
SURF_S(1)=1
SURF_S(2)=1
SURF_S(3)=1
SURF_S(4)=1
SURF_D(1)=2
SURF_D(2)=2
SURF_D(3)=2
SURF_D(4)=2

do i = 1, 4
  ValVertStatic(i)=0.0
  ValVertDynamic(i)=0.0
  ValFaceStatic(i)=0.0
  ValFaceDynamic(i)=0.0
end do
LU_GEOM=10
LU_GEOM_DATA=11

OPEN(UNIT=LU_GEOM,FILE='geom_test.geo',FORM='UNFORMATTED')
OPEN(UNIT=LU_GEOM_DATA,FILE='geom_test_data.geo',FORM='UNFORMATTED')

! static geometry
WRITE(LU_GEOM) ONE
WRITE(LU_GEOM) VERSION

WRITE(LU_GEOM) STIME ! first time step
WRITE(LU_GEOM) N_VERT_S, N_FACE_S, N_VERT_D, N_FACE_D
IF (N_VERT_S>0) WRITE(LU_GEOM) (Xvert_S(I),Yvert_S(I),Zvert_S(I),I=1,N_VERT_S)
IF (N_VERT_D>0) WRITE(LU_GEOM) (Xvert_D(I),Yvert_D(I),Zvert_D(I),I=1,N_VERT_D)
IF (N_FACE_S>0) WRITE(LU_GEOM) (FACE1_S(I),FACE2_S(I),FACE3_S(I),I=1,N_FACE_S)
IF (N_FACE_D>0) WRITE(LU_GEOM) (FACE1_D(I),FACE2_D(I),FACE3_D(I),I=1,N_FACE_D)
IF (N_FACE_S>0) WRITE(LU_GEOM) (SURF_S(I),I=1,N_FACE_S)
IF (N_FACE_D>0) WRITE(LU_GEOM) (SURF_D(I),I=1,N_FACE_D)

WRITE(LU_GEOM_DATA) ONE
WRITE(LU_GEOM_DATA) VERSION
WRITE(LU_GEOM_DATA) STIME
WRITE(LU_GEOM_DATA) N_VERT_S_VALS,N_VERT_D_VALS,N_FACE_S_VALS,N_FACE_D_VALS
IF (N_VERT_S_VALS>0) WRITE(LU_GEOM_DATA) (ValVertStatic(I), I=1,N_VERT_S_VALS)
IF (N_VERT_D_VALS>0) WRITE(LU_GEOM_DATA) (ValVertDynamic(I),I=1,N_VERT_D_VALS)
IF (N_FACE_S_VALS>0) WRITE(LU_GEOM_DATA) (ValFaceStatic(I), I=1,N_VERT_S_VALS)
IF (N_FACE_D_VALS>0) WRITE(LU_GEOM_DATA) (ValFaceDynamic(I),I=1,N_FACE_D_VALS)

twfin=40.0
nsteps=401
dt=twfin/(nsteps-1)
dx = 16.0/(nsteps-1)
do itime=1, nsteps
stime=(itime-1)*dt
GEOM_TYPE=0
xoffset=(itime-1)*dx
yoffset=3.0
zoffset=3.0
Xvert_D(1)=Xvert_S(1)+xoffset
Yvert_D(1)=Yvert_S(1)+yoffset
Zvert_D(1)=Zvert_S(1)+zoffset
Xvert_D(2)=Xvert_S(2)+xoffset
Yvert_D(2)=Yvert_S(2)+yoffset
Zvert_D(2)=Zvert_S(2)+zoffset
Xvert_D(3)=Xvert_S(3)+xoffset
Yvert_D(3)=Yvert_S(3)+yoffset
Zvert_D(3)=Zvert_S(3)+zoffset+(itime-200)*dx
Xvert_D(4)=Xvert_S(4)+xoffset
Yvert_D(4)=Yvert_S(4)+yoffset
Zvert_D(4)=Zvert_S(4)+zoffset
! dynamic geometry
WRITE(LU_GEOM) STIME, GEOM_TYPE ! each successive time step (if there is time dependent geometry)
IF(GEOM_TYPE.EQ.0)THEN
  WRITE(LU_GEOM) N_VERT_D, N_FACE_D
  IF (N_VERT_D>0)THEN
    WRITE(LU_GEOM) (Xvert_D(I),Yvert_D(I),Zvert_D(I),I=1,N_VERT_D)
    WRITE(LU_GEOM) (FACE1_D(I),FACE2_D(I),FACE3_D(I),I=1,N_FACE_D)
    WRITE(LU_GEOM) (SURF_D(I),I=1,N_FACE_D)
  ENDIF
ELSE IF(GEOM_TYPE.EQ.1)THEN
  WRITE(LU_GEOM) Xtran, Ytran, Ztran, Xrot0, Yrot0, Zrot0, rot_az, rot_elev
ENDIF

WRITE(LU_GEOM_DATA) STIME
WRITE(LU_GEOM_DATA) N_VERT_S_VALS,N_VERT_D_VALS,N_FACE_S_VALS,N_FACE_D_VALS
IF (N_VERT_S_VALS>0) WRITE(LU_GEOM_DATA) (ValVertStatic(I), I=1,N_VERT_S_VALS)
IF (N_VERT_D_VALS>0) WRITE(LU_GEOM_DATA) (ValVertDynamic(I),I=1,N_VERT_D_VALS)
IF (N_FACE_S_VALS>0) WRITE(LU_GEOM_DATA) (ValFaceStatic(I), I=1,N_VERT_S_VALS)
IF (N_FACE_D_VALS>0) WRITE(LU_GEOM_DATA) (ValFaceDynamic(I),I=1,N_FACE_D_VALS)

end do
write(6,*)"test geometry output complete"


end program main