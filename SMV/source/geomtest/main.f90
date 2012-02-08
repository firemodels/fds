program main
implicit none

integer :: LU_GEOM
integer :: ONE,VERSION
real :: STIME
integer :: NVERT_S, NFACE_S, NVERT_D, NFACE_D
real, dimension(4) :: Xvert_S, Yvert_S, Zvert_S
real, dimension(4) :: Xvert_D, Yvert_D, Zvert_D
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
xoffset=0.0
yoffset=3.0
zoffset=3.0
ONE=1
VERSION=0
STIME=0.0
NVERT_S=4
NFACE_S=4
NVERT_D=4
NFACE_D=4
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
SURF_S(1)=2
SURF_S(2)=2
SURF_S(3)=2
SURF_S(4)=2
SURF_D(1)=3
SURF_D(2)=3
SURF_D(3)=3
SURF_D(4)=3
LU_GEOM=10

OPEN(UNIT=LU_GEOM,FILE='geom_test.geo',FORM='UNFORMATTED')

! static geometry
WRITE(LU_GEOM) ONE
WRITE(LU_GEOM) VERSION

WRITE(LU_GEOM) STIME ! first time step
WRITE(LU_GEOM) NVERT_S, NFACE_S, NVERT_D, NFACE_D
IF (NVERT_S>0) WRITE(LU_GEOM) (Xvert_S(I),Yvert_S(I),Zvert_S(I),I=1,NVERT_S)
IF (NVERT_D>0) WRITE(LU_GEOM) (Xvert_D(I),Yvert_D(I),Zvert_D(I),I=1,NVERT_D)
IF (NFACE_S>0) WRITE(LU_GEOM) (FACE1_S(I),FACE2_S(I),FACE3_S(I),I=1,NFACE_S)
IF (NFACE_D>0) WRITE(LU_GEOM) (FACE1_D(I),FACE2_D(I),FACE3_D(I),I=1,NFACE_D)
IF (NFACE_S>0) WRITE(LU_GEOM) (SURF_S(I),I=1,NFACE_S)
IF (NFACE_D>0) WRITE(LU_GEOM) (SURF_D(I),I=1,NFACE_D)

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
Zvert_D(3)=Zvert_S(3)+zoffset
Xvert_D(4)=Xvert_S(4)+xoffset
Yvert_D(4)=Yvert_S(4)+yoffset
Zvert_D(4)=Zvert_S(4)+zoffset
! dynamic geometry
WRITE(LU_GEOM) STIME, GEOM_TYPE ! each successive time step (if there is time dependent geometry)
IF(GEOM_TYPE.EQ.0)THEN
  WRITE(LU_GEOM) NVERT_D, NFACE_D
  IF (NVERT_D>0)THEN
    WRITE(LU_GEOM) (Xvert_D(I),Yvert_D(I),Zvert_D(I),I=1,NVERT_D)
    WRITE(LU_GEOM) (FACE1_D(I),FACE2_D(I),FACE3_D(I),I=1,NFACE_D)
    WRITE(LU_GEOM) (SURF_D(I),I=1,NFACE_D)
  ENDIF
ELSE IF(GEOM_TYPE.EQ.1)THEN
  WRITE(LU_GEOM) Xtran, Ytran, Ztran, Xrot0, Yrot0, Zrot0, rot_az, rot_elev
ENDIF
end do
write(6,*)"test geometry output complete"


end program main