! $Date: 2015-04-29 15:16:22 -0400 (Wed, 29 Apr 2015) $ 
! $Revision: 22563 $
! $Author: gforney $

!  WRITE(LU_GEOM) ONE
!  WRITE(LU_GEOM) VERSION
!  WRITE(LU_GEOM) STIME  ! first time step
!  WRITE(LU_GEOM) N_VERT_S, NFACE_S, NVERT_D, N_FACE_D
!  IF (N_VERT_S>0)  WRITE(LU_GEOM) (Xvert_S(I),Yvert_S(I),Zvert_S(I),I=1,N_VERT_S)
!  IF (N_VERT_D>0)  WRITE(LU_GEOM) (Xvert_D(I),Yvert_D(I),Zvert_D(I),I=1,N_VERT_D)
!  IF (N_FACE_S>0)  WRITE(LU_GEOM) (FACE1_S(I),FACE2_S(I),FACE3_S(I),I=1,N_FACE_S)
!  IF (N_FACE_D>0)  WRITE(LU_GEOM) (FACE1_D(I),FACE2_D(I),FACE3_D(I),I=1,N_FACE_D)
!  IF (N_FACE_S>0)  WRITE(LU_GEOM) (SURF_S(I),I=1,N_FACE_S)
!  IF (N_FACE_D>0)  WRITE(LU_GEOM) (SURF_D(I),I=1,N_FACE_D)

!  ------------------ geomout ------------------------ 

subroutine geomout(verts, N_VERT_S, faces, N_FACE_S)
implicit none
integer, intent(in) :: N_VERT_S, N_FACE_S
real, intent(in), dimension(3*N_VERT_S) :: verts
integer, intent(in), dimension(3*N_FACE_S) :: faces;
integer :: LU_GEOM, ONE
real :: STIME
integer :: dummy, VERSION, N_VERT_D, N_FACE_D
integer :: I

ONE=1
LU_GEOM = 40
STIME = 0.0
dummy = 0
VERSION=0
N_VERT_D=0
N_FACE_D=0

OPEN(UNIT=LU_GEOM,FILE="terrain.geom",FORM='UNFORMATTED')
WRITE(LU_GEOM) ONE
WRITE(LU_GEOM) VERSION
WRITE(LU_GEOM) STIME  ! first time step
WRITE(LU_GEOM) N_VERT_S, N_FACE_S, N_VERT_D, N_FACE_D
IF (N_VERT_S>0)  WRITE(LU_GEOM) (verts(3*I-2),verts(3*I-1),verts(3*I),I=1,N_VERT_S)
IF (N_FACE_S>0)  WRITE(LU_GEOM) (faces(3*I-2),faces(3*I-1),faces(3*I),I=1,N_FACE_S)
!IF (N_FACE_S>0)  WRITE(LU_GEOM) (1,I=1,N_FACE_S)
close(LU_GEOM)
end subroutine geomout

!  ------------------ getembeddata ------------------------ 

subroutine getembeddata(filename,ntimes,nvals,times,nstatics,ndynamics,vals,redirect_flag,error)
implicit none
character(len=*), intent(in) :: filename
integer, intent(in) :: ntimes, nvals, redirect_flag
integer, intent(out) :: error
real, intent(out), dimension(:) :: times(ntimes), vals(nvals)
integer, intent(out), dimension(:) :: nstatics(ntimes), ndynamics(ntimes)

integer :: lu20, finish
logical :: isopen,exists
integer :: i
integer :: one, itime, nvars
integer :: nvert_s, ntri_s, nvert_d, ntri_d
real :: valmin, valmax
integer :: version

lu20=20
inquire(unit=lu20,opened=isopen)

if(isopen)close(lu20)
inquire(file=trim(filename),exist=exists)
if(exists)then
#ifdef pp_SHARED
  open(unit=lu20,file=trim(filename),form="unformatted",shared,action="read")
#else
  open(unit=lu20,file=trim(filename),form="unformatted",action="read")
#endif
 else
  write(6,*)' The boundary element file name, ',trim(filename),' does not exist'
  error=1
  return
endif

error = 0
read(lu20)one
read(lu20)version
nvars=0
valmin = 1000000000000.0;
valmax = -valmin
do itime=1, ntimes
  read(lu20,iostat=finish)times(itime)
  if(redirect_flag.eq.0)write(6,10)times(itime)
10 format(" boundary element time=",f9.2)  
  if(finish.eq.0)read(lu20,iostat=finish)nvert_s, ntri_s, nvert_d, ntri_d
  nstatics(itime)=ntri_s
  if(finish.eq.0.and.ntri_s.gt.0)read(lu20,iostat=finish)(vals(nvars+i),i=1,ntri_s)
  ndynamics(itime)=ntri_d
  if(finish.eq.0.and.ntri_d.ne.0)then
    read(lu20,iostat=finish)(vals(nvars+ntri_s+i),i=1,ntri_d)
!    do i = 1, ntri_d
!      vals(nvars+ntri_s+i)=times(itime)
!    end do
  endif
  do i = 1, ntri_s+ntri_d
    if(vals(nvars+i).lt.valmin)valmin=vals(nvars+i)
    if(vals(nvars+i).gt.valmax)valmax=vals(nvars+i)
  end do
  if(finish.ne.0)return
  nvars = nvars + ntri_s + ntri_d
end do
if(redirect_flag.eq.0)write(6,*)" nvars=",nvars,"valmin=",valmin," valmax=",valmax
close(lu20)

end subroutine getembeddata

!  ------------------ getzonedata ------------------------ 

subroutine getzonedata(zonefilename,nzonet,nrooms, nfires, zonet,zoneqfire,zonepr, zoneylay,zonetl,zonetu,error)
implicit none
character(len=*) :: zonefilename
integer, intent(in) :: nrooms, nfires
integer, intent(inout) :: nzonet
real, intent(out), dimension(nrooms*nzonet) :: zonepr, zoneylay, zonetl, zonetu
real, intent(out), dimension(nfires*nzonet) :: zoneqfire
real, intent(out), dimension(nzonet) :: zonet
integer , intent(out) :: error

integer  :: file_unit
integer :: lu26,i,j,ii,ii2,idummy,version
real :: dummy, qdot
logical :: isopen, exists

call get_file_unit(file_unit,26)
lu26 = file_unit

inquire(unit=lu26,opened=isopen)

if(isopen)close(lu26)
inquire(file=trim(zonefilename),exist=exists)
if(exists)then
#ifdef pp_SHARED
  open(unit=lu26,file=trim(zonefilename),form="unformatted",shared,action="read")
#else
  open(unit=lu26,file=trim(zonefilename),form="unformatted",action="read")
#endif
 else
  write(6,*)' The zone file name, ',trim(zonefilename),' does not exist'
  error=1
  return
endif

read(lu26)version
read(lu26)idummy
read(lu26)idummy
ii = 0
ii2 = 0
do j = 1, nzonet
  read(lu26)zonet(j)
  do i = 1, nrooms
    ii = ii + 1
    read(lu26,iostat=error)zonepr(ii),zoneylay(ii),zonetl(ii),zonetu(ii)
    if(error.ne.0)then
      error = 1
      nzonet = j - 1
      close(lu26)
      return
    endif
  end do
  do i = 1, nfires
    ii2 = ii2 + 1
    read(lu26,iostat=error)dummy,qdot
    zoneqfire(ii2) = qdot
    if(error.ne.0)then
      error=1
      nzonet=j-1
      close(lu26)
      return
    endif
  end do 
end do

close(lu26)
end subroutine getzonedata

!  ------------------ getpatchdata ------------------------ 

subroutine getpatchdata(lunit,npatch,pi1,pi2,pj1,pj2,pk1,pk2,patchtime,pqq,npqq,error)
implicit none

integer, intent(in) :: npatch,lunit
integer, intent(in), dimension(*) :: pi1, pi2, pj1, pj2, pk1, pk2
real, intent(out), dimension(*) :: pqq
integer, intent(out) :: error,npqq
real, intent(out) :: patchtime

integer :: i, i1, i2, j1, j2, k1, k2, size, ibeg, iend, lu15, ii

lu15 = lunit
read(lu15,iostat=error)patchtime
if(error.ne.0)then
  close(lu15)
  return
endif
ibeg = 1
npqq=0
do i = 1, npatch
  i1 = pi1(i)
  i2 = pi2(i)
  j1 = pj1(i)
  j2 = pj2(i)
  k1 = pk1(i)
  k2 = pk2(i)
  size = (i2+1-i1)*(j2+1-j1)*(k2+1-k1)
  npqq=npqq+size
  iend = ibeg + size - 1
  read(lu15,iostat=error)(pqq(ii),ii=ibeg,iend)
  if(error.ne.0)then
    close(lu15)
    exit
  endif
  ibeg = iend + 1
end do
return

end subroutine getpatchdata

!  ------------------ getdata1 ------------------------ 

subroutine getdata1(file_unit,ipart,error)
implicit none

integer, intent(in) :: file_unit
integer, intent(out) :: ipart, error

integer :: lu10
real :: sarx, sary, swpar
integer :: i, j, k
integer ndum2
integer :: nspr, nv
integer :: ibar, jbar, kbar
real :: dummy
integer :: nb1, idummy

lu10 = file_unit
error=0

read(lu10,iostat=error) sarx,sary,swpar,ipart,ndum2
if(error.ne.0)return

read(lu10,iostat=error) ibar,jbar,kbar
if(error.ne.0)return

read(lu10,iostat=error) (dummy,i=1,ibar+1),(dummy,j=1,jbar+1),(dummy,k=1,kbar+1)
if(error.ne.0)return

read(lu10,iostat=error) nb1
if(error.ne.0)return

do i=1,nb1
  read(lu10,iostat=error) idummy, idummy, idummy, idummy, idummy, idummy,idummy
  if(error.ne.0)return
end do

read(lu10,iostat=error) nv
if(error.ne.0)return

do i=1,nv
  read(lu10,iostat=error) idummy, idummy, idummy, idummy, idummy, idummy,idummy
  if(error.ne.0)return
end do

read(lu10,iostat=error) nspr
if(error.ne.0)return

do i=1,nspr
  read(lu10,iostat=error) dummy,dummy,dummy
  if(error.ne.0)return
end do


return
end subroutine getdata1

!  ------------------ getdata2 ------------------------ 

subroutine getdata2(file_unit,xs,ys,zs,&
                    t,&
                    sprinkflag,isprink,tspr,bframe,sframe,sprframe,stimes,nspr,nmax,mxframes,nframes,&
                    settmin_p,settmax_p,tmin_p,tmax_p,frameloadstep,partpointstep, &
                    xbox0, xbox, ybox0, ybox, zbox0, zbox, &
                    offset_x, offset_y, offset_z, redirect_flag, &
                    error)
                   
implicit none
real, dimension(*), intent(out) :: t

integer(2), dimension(*) :: xs, ys, zs
integer, dimension(*) :: bframe, sframe, sprframe
character(len=1), dimension(*) :: isprink
real, dimension(*) ::  stimes,tspr
integer, intent(in) :: file_unit,nspr,nmax, mxframes, redirect_flag
integer, intent(out) :: nframes,error
integer, intent(in) :: settmin_p, settmax_p, frameloadstep, partpointstep
integer, intent(out) :: sprinkflag
real, intent(in) :: tmin_p, tmax_p
real, intent(in) :: xbox0, xbox, ybox0, ybox, zbox0, zbox
real, intent(in) :: offset_x, offset_y, offset_z

real :: stime

integer :: lu10, i
integer :: npp1, nn, n, npp2, naspr
real, allocatable, dimension(:) :: xpp, ypp, zpp, brp
integer, allocatable, dimension(:) :: iitemp
integer, allocatable, dimension(:) :: oldispr, ispr
integer :: npp1a, npp2a, nppold
integer :: factor

logical :: load, allocated
integer :: nf, npoints

lu10 = file_unit
factor=2**15

nframes = 0
bframe(1) = 0
nf = 0
n = 0
sprinkflag=0
if(nspr.gt.0)then
  allocate(oldispr(nspr),ispr(nspr))
  do i = 1, nspr
   oldispr(i)=0
   tspr(i)=99999.
  end do
endif

nppold=100
npp1=100
allocate(iitemp(npp1),xpp(npp1),ypp(npp1),zpp(npp1),brp(npp1),stat=error)

do
  npp1a = 0
  npp2a = 0
  read(lu10,iostat=error) stime,npp1,nn,(ispr(i),i=1,nspr)
  if(error.ne.0)go to 999
  naspr = 0
  do i = 1, nspr
    if(oldispr(i).eq.0.and.ispr(i).ne.0)tspr(i)=stime
    oldispr(i) = ispr(i)
    if(ispr(i).ne.0)naspr = naspr + 1
  end do

  allocated = .false.
  if(npp1.gt.0.and.npp1.gt.nppold)then
    allocated=.true.
    nppold = npp1
   else
    allocated = .false.
  endif
  if(allocated)then
    deallocate(iitemp,xpp,ypp,zpp,brp)
    allocate(iitemp(npp1),xpp(npp1),ypp(npp1),zpp(npp1),brp(npp1),stat=error)
    if(error.ne.0)go to 999
  endif
  read(lu10,iostat=error) (xpp(i),i=1,npp1),(ypp(i),i=1,npp1),(zpp(i),i=1,npp1),(brp(i),i=1,npp1)
  if(error.ne.0)go to 999
  if((settmin_p.ne.0.and.stime.lt.tmin_p).or.mod(nf,frameloadstep).ne.0)then
    load=.false.
   else
    load=.true.
  endif
  nf = nf + 1
  if(settmax_p.ne.0.and.stime.gt.tmax_p)go to 999
  if(load.and.nframes+1.gt.mxframes)go to 999

  if(load)then
    npoints = (npp1-1)/partpointstep + 1
    if(n+npoints.gt.nmax)then
      go to 999
    endif
   iitemp(1:npoints) = factor*(offset_x+xpp(1:npp1:partpointstep)-xbox0)/(xbox-xbox0)
   where(iitemp(1:npoints)<0)iitemp(1:npoints)=0
   where(iitemp(1:npoints).gt.factor-1)iitemp(1:npoints)=factor-1
   xs(n+1:n+npoints) = iitemp(1:npoints)

   iitemp(1:npoints) = factor*(offset_y+ypp(1:npp1:partpointstep)-ybox0)/(ybox-ybox0)
   where(iitemp(1:npoints)<0)iitemp(1:npoints)=0
   where(iitemp(1:npoints).gt.factor-1)iitemp(1:npoints)=factor-1
   ys(n+1:n+npoints) = iitemp(1:npoints)

   iitemp(1:npoints) = factor*(offset_z+zpp(1:npp1:partpointstep)-zbox0)/(zbox-zbox0)
   where(iitemp(1:npoints)<0)iitemp(1:npoints)=0
   where(iitemp(1:npoints).gt.factor-1)iitemp(1:npoints)=factor-1
   zs(n+1:n+npoints) = iitemp(1:npoints)

   t(n+1:n+npoints) = brp(1:npp1:partpointstep)
   isprink(n+1:n+npoints) = char(0)
   npp1a = npoints
   n = n + npp1a
  endif
  npp2 = 0


  if(naspr.ne.0)then       ! read in sprinkler data
    read(lu10,iostat=error) npp2
    if(npp2.lt.0)sprinkflag=1
    npp2 = abs(npp2)
    if(error.ne.0)go to 999
    if(npp2.gt.0.and.npp2.gt.nppold)then
      allocated=.true.
      nppold = npp2
     else
      allocated = .false.
    endif
    if(allocated)then
      deallocate(iitemp,xpp,ypp,zpp,brp)
      allocate(iitemp(npp2),xpp(npp2),ypp(npp2),zpp(npp2),brp(npp2),stat=error)
    endif
    if(sprinkflag.eq.0)then
      read(lu10,iostat=error) (xpp(i),i=1,npp2),(ypp(i),i=1,npp2),(zpp(i),i=1,npp2)
     else
      read(lu10,iostat=error) (xpp(i),i=1,npp2),(ypp(i),i=1,npp2),(zpp(i),i=1,npp2),(brp(i),i=1,npp2)
    endif
    if(load)then
      npoints = (npp2-1)/partpointstep + 1
      if(n+npoints.gt.nmax)go to 999
      if(sprinkflag.eq.1)then
        t(n+1:n+npoints) = brp(1:npp2:partpointstep)
       else
        t(n+1:n+npoints) = -1.0
      endif
      isprink(n+1:n+npoints) = char(1)
        iitemp(1:npoints) = factor*(xpp(1:npp2:partpointstep)-xbox0)/(xbox-xbox0)
        where(iitemp(1:npoints)<0)iitemp(1:npoints)=0
        where(iitemp(1:npoints).gt.factor-1)iitemp(1:npoints)=factor-1
        xs(n+1:n+npoints) = iitemp(1:npoints)

        iitemp(1:npoints) = factor*(ypp(1:npp2:partpointstep)-ybox0)/(ybox-ybox0)
        where(iitemp(1:npoints)<0)iitemp(1:npoints)=0
        where(iitemp(1:npoints).gt.factor-1)iitemp(1:npoints)=factor-1
        ys(n+1:n+npoints) = iitemp(1:npoints)

        iitemp(1:npoints) = factor*(zpp(1:npp2:partpointstep)-zbox0)/(zbox-zbox0)
        where(iitemp(1:npoints)<0)iitemp(1:npoints)=0
        where(iitemp(1:npoints).gt.factor-1)iitemp(1:npoints)=factor-1
        zs(n+1:n+npoints) = iitemp(1:npoints)
        npp2a = npoints
        n = n + npp2a
    endif
  end if
  if(error.ne.0)goto 999
  if(load)then
    nframes = nframes + 1
    stimes(nframes) = stime
    sframe(nframes) = npp1a + npp2a
    sprframe(nframes) = npp2a
    if(nframes+1.le.mxframes)bframe(nframes+1) = bframe(nframes) + sframe(nframes)
    if(npp2.eq.0)then
      if(redirect_flag.eq.0)write(6,10)stime
10 format(" particle time=",f9.2)
     else
      if(redirect_flag.eq.0)then
        write(6,*)" particle time=",stime,"particles",npp1,"droplets",npp2
        write(6,20)stime,npp2
      endif
20 format(" particle time=",f9.2," particles",i9," droplets",i9)      
    endif
  endif

end do
999 continue
error = 0
close(lu10)
return
end subroutine getdata2

!  ------------------ getdirval ------------------------ 

subroutine getdirval(is1,is2,js1,js2,ks1,ks2,idir,joff,koff)
implicit none
integer :: nxsp, nysp, nzsp
integer, intent(in) :: is1, js1, ks1
integer, intent(inout) :: is2, js2, ks2
integer, intent(out) :: idir, koff, joff
integer :: imin

nxsp = is2 + 1 - is1
nysp = js2 + 1 - js1
nzsp = ks2 + 1 - ks1
joff=0
koff=0
if(is1.ne.is2.and.js1.ne.js2.and.ks1.ne.ks2)then
  idir=1
  is2 = is1
  return
endif
imin = min(nxsp,nysp,nzsp)
if(nxsp.eq.imin)then
  idir = 1
  is2 = is1
 elseif(nysp.eq.imin)then
  idir = 2
  js2 = js1
 else
  idir = 3
  ks2 = ks1
endif
if(is1.eq.is2.and.js1.eq.js2)then
   idir=1
   joff=1
  elseif(is1.eq.is2.and.ks1.eq.ks2)then
   idir=1
   koff=1
  elseif(js1.eq.js2.and.ks1.eq.ks2)then
   idir=2
   koff=1
endif
return
end subroutine getdirval

!  ------------------ writeslicedata ------------------------ 

subroutine writeslicedata(file_unit,slicefilename,is1,is2,js1,js2,ks1,ks2,qdata,times,ntimes,redirect_flag)
implicit none


integer, intent(in) :: file_unit
character(len=*),intent(in) :: slicefilename
integer, intent(in) :: is1, is2, js1, js2, ks1, ks2, redirect_flag
real, intent(in), dimension(*) :: qdata
real, intent(in), dimension(*) :: times
integer, intent(in) :: ntimes

logical :: connected
integer :: error
character(len=30) :: longlbl, shortlbl, unitlbl
integer :: ibeg, iend, nframe
integer :: nxsp, nysp, nzsp
integer :: i,ii
inquire(unit=file_unit,opened=connected)
if(connected)close(file_unit)

open(unit=file_unit,file=trim(slicefilename),form="unformatted",action="write")

longlbl=" "
shortlbl=" "
unitlbl=" "

write(file_unit,iostat=error)longlbl
write(file_unit,iostat=error)shortlbl
write(file_unit,iostat=error)unitlbl

write(file_unit,iostat=error)is1, is2, js1, js2, ks1, ks2
nxsp = is2 + 1 - is1
nysp = js2 + 1 - js1
nzsp = ks2 + 1 - ks1
nframe=nxsp*nysp*nzsp
do i = 1, ntimes
  write(file_unit)times(i)
  if(redirect_flag.eq.0)write(6,10)times(i)
10 format(" outputting slice time=",f10.2)
  ibeg=1+(i-1)*nframe
  iend=i*nframe
  write(file_unit)(qdata(ii),ii=ibeg,iend)
end do

close(file_unit)

return
end subroutine writeslicedata

!  ------------------ getslicedata ------------------------ 

subroutine getslicedata(file_unit,slicefilename,shortlabel,&
            is1,is2,js1,js2,ks1,ks2,idir,qmin,qmax,qdata,times,nstepsmax,sliceframestep,&
            settmin_s,settmax_s,tmin_s,tmax_s,redirect_flag)
implicit none

character(len=*) :: slicefilename, shortlabel

integer, intent(in) :: file_unit
real, intent(out) :: qmin, qmax
real, intent(out), dimension(*) :: qdata
real, intent(out), dimension(*) :: times
integer, intent(out) :: idir
integer, intent(out) :: is1, is2, js1, js2, ks1, ks2
integer, intent(in) :: redirect_flag
integer, intent(inout) :: nstepsmax
integer, intent(in) :: settmin_s, settmax_s, sliceframestep
real, intent(in) :: tmin_s, tmax_s

real, dimension(:,:,:), pointer :: qq

integer :: i,j,k
integer :: lu11, nsteps
logical :: exists
integer :: ip1, ip2, jp1, jp2, kp1, kp2
integer :: nxsp, nysp, nzsp
integer :: error, istart, irowstart
real :: time, time_max
character(len=30) :: longlbl, shortlbl, unitlbl
integer :: lenshort, lenunits
character(len=3) :: blank
logical :: connected, load
integer :: ii, kk
integer :: joff, koff
integer :: count

lu11 = file_unit
joff = 0
koff = 0
inquire(unit=lu11,opened=connected)
if(connected)close(lu11)

inquire(file=trim(slicefilename),exist=exists)
if(exists)then
#ifdef pp_SHARED
  open(unit=lu11,file=trim(slicefilename),form="unformatted",shared,action="read")
#else
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read")
#endif
 else
  write(6,*)' the slice file ',trim(slicefilename),' does not exist'
  nsteps = 0
  return
endif

nsteps = 0
blank = '   '
longlbl=" "
shortlbl=" "
unitlbl=" "

read(lu11,iostat=error)longlbl
read(lu11,iostat=error)shortlbl
read(lu11,iostat=error)unitlbl

! longlabel=trim(longlbl)//char(0)

lenshort = min(len_trim(shortlabel),6)
! shortlabel=shortlbl(1:lenshort)//char(0)

lenunits = min(len_trim(unitlbl),6)
! units=unitlbl(1:lenunits)//char(0)

read(lu11,iostat=error)ip1, ip2, jp1, jp2, kp1, kp2
is1 = ip1 
is2 = ip2 
js1 = jp1 
js2 = jp2 
ks1 = kp1 
ks2 = kp2 
if(error.ne.0)then
  close(lu11)
  return
endif

nxsp = is2 + 1 - is1
nysp = js2 + 1 - js1
nzsp = ks2 + 1 - ks1  
call getdirval(is1,is2,js1,js2,ks1,ks2,idir,joff,koff)

allocate(qq(nxsp,nysp+joff,nzsp+koff))

qmin = 1.0e30
qmax = -1.0e30
count=-1
time_max=-1000000.0
do
  read(lu11,iostat=error)time
  if(error.ne.0)exit
  if((settmin_s.ne.0.and.time<tmin_s).or.time.le.time_max)then
    load = .false.
   else
    load = .true.
    time_max = time
  endif
  if(settmax_s.ne.0.and.time>tmax_s)exit
  read(lu11,iostat=error)(((qq(i,j,k),i=1,nxsp),j=1,nysp),k=1,nzsp)
  count=count+1
  if(mod(count,sliceframestep).ne.0)load = .false.
  if(koff.eq.1)then
    qq(1:nxsp,1:nysp,2)=qq(1:nxsp,1:nysp,1)
   elseif(joff.eq.1)then
    qq(1:nxsp,2,1:nzsp)=qq(1:nxsp,1,1:nzsp)
  endif
  if(error.ne.0.or.nsteps.ge.nstepsmax)go to 999
  if(.not.load)cycle
  nsteps = nsteps + 1
  times(nsteps) = time
  if(redirect_flag.eq.0)write(6,10)time
10 format(" slice time=",f9.2)  
  if(idir.eq.3)then
    istart = (nsteps-1)*nxsp*nysp
    do i = 1, nxsp
      irowstart = (i-1)*nysp
      ii = istart+irowstart
      qdata(ii+1:ii+nysp)=qq(i,1:nysp,1)
      qmax = max(qmax,maxval(qq(i,1:nysp,1)))
      qmin = min(qmin,minval(qq(i,1:nysp,1)))
    end do
   elseif(idir.eq.2)then
    istart = (nsteps-1)*nxsp*(nzsp+koff)
    do i = 1, nxsp
      irowstart = (i-1)*(nzsp+koff)
      kk = istart + irowstart
      qdata(kk+1:kk+nzsp+koff) = qq(i,1,1:nzsp+koff)
      qmax = max(qmax,maxval(qq(i,1,1:nzsp+koff)))
      qmin = min(qmin,minval(qq(i,1,1:nzsp+koff)))
    end do
   else
    istart = (nsteps-1)*(nysp+joff)*(nzsp+koff)*nxsp
    do i = 1, nxsp
    do j = 1, nysp+joff
      irowstart = (i-1)*nysp*(nzsp+koff)+(j-1)*(nzsp+koff)
      kk = istart + irowstart
      qdata(kk+1:kk+nzsp+koff) = qq(i,j,1:nzsp+koff)
      qmax = max(qmax,maxval(qq(i,j,1:nzsp+koff)))
      qmin = min(qmin,minval(qq(i,j,1:nzsp+koff)))
    end do
    end do
  endif

end do

999 continue
ks2 = ks2 + koff
js2 = js2 + joff
nstepsmax=nsteps
deallocate(qq)
close(lu11)

return
end subroutine getslicedata

!  ------------------ getsliceframe ------------------------ 

subroutine getsliceframe(lu11,is1,is2,js1,js2,ks1,ks2,time,qframe,testslice,error)
implicit none

real, intent(out), dimension(*) :: qframe
real, intent(out) :: time
integer, intent(out) :: error
integer, intent(in) :: lu11, is1, is2, js1, js2, ks1, ks2
integer, intent(in) :: testslice

integer :: i,j,k
integer :: nxsp, nysp, nzsp
real :: val,factor
integer :: index
real :: ii, jj, kk

nxsp = is2 + 1 - is1
nysp = js2 + 1 - js1
nzsp = ks2 + 1 - ks1  

read(lu11,iostat=error)time
if(error.ne.0)return
read(lu11,iostat=error)(((qframe(1+i+j*nxsp+k*nxsp*nysp),i=0,nxsp-1),j=0,nysp-1),k=0,nzsp-1)
if(testslice.eq.1.or.testslice.eq.2)then
  factor=1.0
  if(testslice.eq.2)factor=1.1
  do k = 0, nzsp-1
    kk = 2.0*((nzsp-1)/2.0-k)/(nzsp-1.0)
    do j = 0, nysp-1
      jj = 2.0*((nysp-1)/2.0-j)/(nysp-1.0)
      do i = 0, nxsp-1
        ii = 2.0*((nxsp-1)/2.0-i)/(nxsp-1.0)
        val = factor*(time-20.0)*(ii*ii + jj*jj + kk*kk)/20.0
        index = 1+i+j*nxsp+k*nxsp*nysp
        qframe(index) = val
      end do
    end do
  end do
endif

999 continue

return
end subroutine getsliceframe

!  ------------------ endian_out ------------------------ 

subroutine endianout(endianfilename)
implicit none
character(len=*) :: endianfilename
integer :: one
integer :: file_unit

file_unit=31

call get_file_unit(file_unit,file_unit)
open(unit=file_unit,file=trim(endianfilename),form="unformatted")
one=1
write(31)one
close(file_unit)
return
end subroutine endianout

!  ------------------ outsliceheader ------------------------ 

subroutine outsliceheader(slicefilename,unit,ip1, ip2, jp1, jp2, kp1, kp2, error)
implicit none

character(len=*) :: slicefilename
integer, intent(in) :: unit
integer, intent(in) :: ip1, ip2, jp1, jp2, kp1, kp2
integer, intent(out) :: error

character(len=30) :: longlbl, shortlbl, unitlbl
integer :: lu11
logical :: connected

lu11 = unit
inquire(unit=lu11,opened=connected)
if(connected)close(lu11)

open(unit=lu11,file=trim(slicefilename),form="unformatted")

longlbl= "long                          "
shortlbl="short                         "
unitlbl= "unit                          "
write(lu11,iostat=error)longlbl
write(lu11,iostat=error)shortlbl
write(lu11,iostat=error)unitlbl

write(lu11,iostat=error)ip1, ip2, jp1, jp2, kp1, kp2

end subroutine outsliceheader

!  ------------------ outsliceframe ------------------------ 

subroutine outsliceframe(lu11,is1,is2,js1,js2,ks1,ks2,time,qframe,error)
implicit none

real, intent(in), dimension(*) :: qframe
real, intent(in) :: time
integer, intent(out) :: error
integer, intent(in) :: lu11, is1, is2, js1, js2, ks1, ks2

integer :: i,j,k
integer :: nxsp, nysp, nzsp

nxsp = is2 + 1 - is1
nysp = js2 + 1 - js1
nzsp = ks2 + 1 - ks1  

write(lu11,iostat=error)time
if(error.ne.0)return
write(lu11,iostat=error)(((qframe(1+i+j*nxsp+k*nxsp*nysp),i=0,nxsp-1),j=0,nysp-1),k=0,nzsp-1)

999 continue

return
end subroutine outsliceframe

!  ------------------ outboundaryheader ------------------------ 

subroutine outboundaryheader(boundaryfilename,boundaryunitnumber,npatches,pi1,pi2,pj1,pj2,pk1,pk2,patchdir,error)
implicit none

character(len=*), intent(in) :: boundaryfilename
integer, intent(in) :: boundaryunitnumber, npatches
integer, intent(in), dimension(npatches) :: pi1, pi2, pj1, pj2, pk1, pk2, patchdir
integer, intent(out) :: error

character(len=30) :: blank
integer :: n, lu15

error=0
lu15 = boundaryunitnumber
open(unit=lu15,file=trim(boundaryfilename),form="unformatted",action="write")

blank="                              "
write(lu15)blank
write(lu15)blank
write(lu15)blank
write(lu15)npatches

do n = 1, npatches
  write(lu15)pi1(n), pi2(n), pj1(n), pj2(n), pk1(n), pk2(n), patchdir(n)
end do

return
end subroutine outboundaryheader

!  ------------------ outpatchframe ------------------------ 

subroutine outpatchframe(lunit,npatch,pi1,pi2,pj1,pj2,pk1,pk2,patchtime,pqq,error)
implicit none

integer, intent(in) :: npatch,lunit
integer, intent(in), dimension(*) :: pi1, pi2, pj1, pj2, pk1, pk2
real, intent(in), dimension(*) :: pqq
integer, intent(out) :: error
real, intent(in) :: patchtime

integer :: i, i1, i2, j1, j2, k1, k2, size, ibeg, iend, lu15, ii

error=0
lu15 = lunit
write(lu15)patchtime
ibeg = 1
do i = 1, npatch
  i1 = pi1(i)
  i2 = pi2(i)
  j1 = pj1(i)
  j2 = pj2(i)
  k1 = pk1(i)
  k2 = pk2(i)
  size = (i2+1-i1)*(j2+1-j1)*(k2+1-k1)
  iend = ibeg + size - 1
  write(lu15)(pqq(ii),ii=ibeg,iend)
  ibeg = iend + 1
end do
return

end subroutine outpatchframe

!  ------------------ getplot3dq ------------------------ 

subroutine getplot3dq(qfilename,nx,ny,nz,qq,error,isotest)
implicit none

character(len=*) :: qfilename
integer, intent(in) :: nx, ny, nz
integer, intent(out) :: error
real, dimension(nx,ny,nz,5)  :: qq
integer, intent(in) :: isotest

real :: dum1, dum2, dum3, dum4
logical :: exists
integer :: error2
real :: dummy, qval

integer :: nxpts, nypts, nzpts
integer :: i, j, k, n
real :: r2

integer :: u_in
logical :: connected

if(isotest.eq.0)then
  call get_file_unit(u_in,70)
  inquire(unit=u_in,opened=connected)
  if(connected)close(u_in)

  error=0
  inquire(file=qfilename,exist=exists)
  if(exists)then
#ifdef pp_SHARED
    open(unit=u_in,file=qfilename,form="unformatted",shared,action="read",iostat=error2)
#else
    open(unit=u_in,file=qfilename,form="unformatted",action="read",iostat=error2)
#endif
   else
    write(6,*)' The file name, ',trim(qfilename),' does not exist'
    read(5,*)dummy
    stop
  endif

  read(u_in,iostat=error)nxpts, nypts, nzpts
  if(nx.eq.nxpts.and.ny.eq.nypts.and.nz.eq.nzpts)then
    read(u_in,iostat=error)dum1, dum2, dum3, dum4
    read(u_in,iostat=error)((((qq(i,j,k,n),i=1,nxpts),j=1,nypts),k=1,nzpts),n=1,5)
   else
    error = 1
    write(6,*)" *** Fatal error in getplot3dq ***"
    write(6,*)" Grid size found in plot3d file was:",nxpts,nypts,nzpts
    write(6,*)" Was expecting:",nx,ny,nz
    stop
  endif
  close(u_in)
 else
    do i = 1, nx
    do j = 1, ny
    do k = 1, nz
      r2 = sqrt(float((i-nx/2)**2 + (j-ny/2)**2 + (k-nz/2)**2))
      qval = r2
      if(isotest.eq.1)then
        qq(i,j,k,1) = 0.0
        qq(i,j,k,2) = 0.0
        qq(i,j,k,3:5) = qval
      endif
      if(isotest.eq.2)then
        qq(i,j,k,1)=qval
        qq(i,j,k,2)=1.1*qval
        qq(i,j,k,3:5)=1.1*qval
      endif
    end do
    end do
  end do
  error = 0
endif
close(u_in)
return
end subroutine getplot3dq

!  ------------------ plot3dout ------------------------ 

subroutine plot3dout(outfile, nx, ny, nz, qout, error3)
implicit none

character(len=*), intent(in) :: outfile
integer, intent(in) :: nx, ny, nz
real, dimension(nx,ny,nz,5)  :: qout
integer, intent(out) :: error3

integer :: u_out
logical :: connected
integer :: i, j, k, n
real :: dummy

error3 = 0

call get_file_unit(u_out,70)
inquire(unit=u_out,opened=connected)
if(connected)close(u_out)

dummy = 0.0
open(unit=u_out,file=trim(outfile),form="unformatted",action="write",iostat=error3)
if(error3.ne.0)return

write(u_out,iostat=error3)nx, ny, nz
write(u_out,iostat=error3)dummy, dummy, dummy, dummy
write(u_out,iostat=error3)((((qout(i,j,k,n),i=1,nx),j=1,ny),k=1,nz),n=1,5)
close(u_out)

return
end subroutine plot3dout

SUBROUTINE color2rgb(RGB,COLOR)
implicit none

! Translate character string of a color name to RGB value

INTEGER :: RGB(3)
CHARACTER(len=*) :: COLOR

SELECT CASE(COLOR)
CASE ('ALICE BLUE');RGB = (/240,248,255/)
CASE ('ANTIQUE WHITE');RGB = (/250,235,215/)
CASE ('ANTIQUE WHITE 1');RGB = (/255,239,219/)
CASE ('ANTIQUE WHITE 2');RGB = (/238,223,204/)
CASE ('ANTIQUE WHITE 3');RGB = (/205,192,176/)
CASE ('ANTIQUE WHITE 4');RGB = (/139,131,120/)
CASE ('AQUAMARINE');RGB = (/127,255,212/)
CASE ('AQUAMARINE 1');RGB = (/118,238,198/)
CASE ('AQUAMARINE 2');RGB = (/102,205,170/)
CASE ('AQUAMARINE 3');RGB = (/69,139,116/)
CASE ('AZURE');RGB = (/240,255,255/)
CASE ('AZURE 1');RGB = (/224,238,238/)
CASE ('AZURE 2');RGB = (/193,205,205/)
CASE ('AZURE 3');RGB = (/131,139,139/)
CASE ('BANANA');RGB = (/227,207,87/)
CASE ('BEIGE');RGB = (/245,245,220/)
CASE ('BISQUE');RGB = (/255,228,196/)
CASE ('BISQUE 1');RGB = (/238,213,183/)
CASE ('BISQUE 2');RGB = (/205,183,158/)
CASE ('BISQUE 3');RGB = (/139,125,107/)
CASE ('BLACK');RGB = (/0,0,0/)
CASE ('BLANCHED ALMOND');RGB = (/255,235,205/)
CASE ('BLUE');RGB = (/0,0,255/)
CASE ('BLUE 2');RGB = (/0,0,238/)
CASE ('BLUE 3');RGB = (/0,0,205/)
CASE ('BLUE 4');RGB = (/0,0,139/)
CASE ('BLUE VIOLET');RGB = (/138,43,226/)
CASE ('BRICK');RGB = (/156,102,31/)
CASE ('BROWN');RGB = (/165,42,42/)
CASE ('BROWN 1');RGB = (/255,64,64/)
CASE ('BROWN 2');RGB = (/238,59,59/)
CASE ('BROWN 3');RGB = (/205,51,51/)
CASE ('BROWN 4');RGB = (/139,35,35/)
CASE ('BURLY WOOD');RGB = (/222,184,135/)
CASE ('BURLY WOOD 1');RGB = (/255,211,155/)
CASE ('BURLY WOOD 2');RGB = (/238,197,145/)
CASE ('BURLY WOOD 3');RGB = (/205,170,125/)
CASE ('BURLY WOOD 4');RGB = (/139,115,85/)
CASE ('BURNT ORANGE');RGB = (/204,85,0/)
CASE ('BURNT SIENNA');RGB = (/138,54,15/)
CASE ('BURNT UMBER');RGB = (/138,51,36/)
CASE ('CADET BLUE');RGB = (/95,158,160/)
CASE ('CADET BLUE 1');RGB = (/152,245,255/)
CASE ('CADET BLUE 2');RGB = (/142,229,238/)
CASE ('CADET BLUE 3');RGB = (/122,197,205/)
CASE ('CADET BLUE 4');RGB = (/83,134,139/)
CASE ('CADMIUM ORANGE');RGB = (/255,97,3/)
CASE ('CADMIUM YELLOW');RGB = (/255,153,18/)
CASE ('CARROT');RGB = (/237,145,33/)
CASE ('CHARTREUSE');RGB = (/127,255,0/)
CASE ('CHARTREUSE 1');RGB = (/118,238,0/)
CASE ('CHARTREUSE 2');RGB = (/102,205,0/)
CASE ('CHARTREUSE 3');RGB = (/69,139,0/)
CASE ('CHOCOLATE');RGB = (/210,105,30/)
CASE ('CHOCOLATE 1');RGB = (/255,127,36/)
CASE ('CHOCOLATE 2');RGB = (/238,118,33/)
CASE ('CHOCOLATE 3');RGB = (/205,102,29/)
CASE ('CHOCOLATE 4');RGB = (/139,69,19/)
CASE ('COBALT');RGB = (/61,89,171/)
CASE ('COBALT GREEN');RGB = (/61,145,64/)
CASE ('COLD GREY');RGB = (/128,138,135/)
CASE ('CORAL');RGB = (/255,127,80/)
CASE ('CORAL 1');RGB = (/255,114,86/)
CASE ('CORAL 2');RGB = (/238,106,80/)
CASE ('CORAL 3');RGB = (/205,91,69/)
CASE ('CORAL 4');RGB = (/139,62,47/)
CASE ('CORNFLOWER BLUE');RGB = (/100,149,237/)
CASE ('CORNSILK');RGB = (/255,248,220/)
CASE ('CORNSILK 1');RGB = (/238,232,205/)
CASE ('CORNSILK 2');RGB = (/205,200,177/)
CASE ('CORNSILK 3');RGB = (/139,136,120/)
CASE ('CRIMSON');RGB = (/220,20,60/)
CASE ('CYAN');RGB = (/0,255,255/)
CASE ('CYAN 2');RGB = (/0,238,238/)
CASE ('CYAN 3');RGB = (/0,205,205/)
CASE ('CYAN 4');RGB = (/0,139,139/)
CASE ('DARK GOLDENROD');RGB = (/184,134,11/)
CASE ('DARK GOLDENROD 1');RGB = (/255,185,15/)
CASE ('DARK GOLDENROD 2');RGB = (/238,173,14/)
CASE ('DARK GOLDENROD 3');RGB = (/205,149,12/)
CASE ('DARK GOLDENROD 4');RGB = (/139,101,8/)
CASE ('DARK GRAY');RGB = (/169,169,169/)
CASE ('DARK GREEN');RGB = (/0,100,0/)
CASE ('DARK KHAKI');RGB = (/189,183,107/)
CASE ('DARK OLIVE GREEN');RGB = (/85,107,47/)
CASE ('DARK OLIVE GREEN 1');RGB = (/202,255,112/)
CASE ('DARK OLIVE GREEN 2');RGB = (/188,238,104/)
CASE ('DARK OLIVE GREEN 3');RGB = (/162,205,90/)
CASE ('DARK OLIVE GREEN 4');RGB = (/110,139,61/)
CASE ('DARK ORANGE');RGB = (/255,140,0/)
CASE ('DARK ORANGE 1');RGB = (/255,127,0/)
CASE ('DARK ORANGE 2');RGB = (/238,118,0/)
CASE ('DARK ORANGE 3');RGB = (/205,102,0/)
CASE ('DARK ORANGE 4');RGB = (/139,69,0/)
CASE ('DARK ORCHID');RGB = (/153,50,204/)
CASE ('DARK ORCHID 1');RGB = (/191,62,255/)
CASE ('DARK ORCHID 2');RGB = (/178,58,238/)
CASE ('DARK ORCHID 3');RGB = (/154,50,205/)
CASE ('DARK ORCHID 4');RGB = (/104,34,139/)
CASE ('DARK SALMON');RGB = (/233,150,122/)
CASE ('DARK SEA GREEN');RGB = (/143,188,143/)
CASE ('DARK SEA GREEN 1');RGB = (/193,255,193/)
CASE ('DARK SEA GREEN 2');RGB = (/180,238,180/)
CASE ('DARK SEA GREEN 3');RGB = (/155,205,155/)
CASE ('DARK SEA GREEN 4');RGB = (/105,139,105/)
CASE ('DARK SLATE BLUE');RGB = (/72,61,139/)
CASE ('DARK SLATE GRAY');RGB = (/47,79,79/)
CASE ('DARK SLATE GRAY 1');RGB = (/151,255,255/)
CASE ('DARK SLATE GRAY 2');RGB = (/141,238,238/)
CASE ('DARK SLATE GRAY 3');RGB = (/121,205,205/)
CASE ('DARK SLATE GRAY 4');RGB = (/82,139,139/)
CASE ('DARK TURQUOISE');RGB = (/0,206,209/)
CASE ('DARK VIOLET');RGB = (/148,0,211/)
CASE ('DEEP PINK');RGB = (/255,20,147/)
CASE ('DEEP PINK 1');RGB = (/238,18,137/)
CASE ('DEEP PINK 2');RGB = (/205,16,118/)
CASE ('DEEP PINK 3');RGB = (/139,10,80/)
CASE ('DEEP SKYBLUE');RGB = (/0,191,255/)
CASE ('DEEP SKYBLUE 1');RGB = (/0,178,238/)
CASE ('DEEP SKYBLUE 2');RGB = (/0,154,205/)
CASE ('DEEP SKYBLUE 3');RGB = (/0,104,139/)
CASE ('DIM GRAY');RGB = (/105,105,105/)
CASE ('DODGERBLUE');RGB = (/30,144,255/)
CASE ('DODGERBLUE 1');RGB = (/28,134,238/)
CASE ('DODGERBLUE 2');RGB = (/24,116,205/)
CASE ('DODGERBLUE 3');RGB = (/16,78,139/)
CASE ('EGGSHELL');RGB = (/252,230,201/)
CASE ('EMERALD GREEN');RGB = (/0,201,87/)
CASE ('FIREBRICK');RGB = (/178,34,34/)
CASE ('FIREBRICK 1');RGB = (/255,48,48/)
CASE ('FIREBRICK 2');RGB = (/238,44,44/)
CASE ('FIREBRICK 3');RGB = (/205,38,38/)
CASE ('FIREBRICK 4');RGB = (/139,26,26/)
CASE ('FLESH');RGB = (/255,125,64/)
CASE ('FLORAL WHITE');RGB = (/255,250,240/)
CASE ('FOREST GREEN');RGB = (/34,139,34/)
CASE ('GAINSBORO');RGB = (/220,220,220/)
CASE ('GHOST WHITE');RGB = (/248,248,255/)
CASE ('GOLD');RGB = (/255,215,0/)
CASE ('GOLD 1');RGB = (/238,201,0/)
CASE ('GOLD 2');RGB = (/205,173,0/)
CASE ('GOLD 3');RGB = (/139,117,0/)
CASE ('GOLDENROD');RGB = (/218,165,32/)
CASE ('GOLDENROD 1');RGB = (/255,193,37/)
CASE ('GOLDENROD 2');RGB = (/238,180,34/)
CASE ('GOLDENROD 3');RGB = (/205,155,29/)
CASE ('GOLDENROD 4');RGB = (/139,105,20/)
CASE ('GRAY');RGB = (/128,128,128/)
CASE ('GRAY 1');RGB = (/3,3,3/)
CASE ('GRAY 10');RGB = (/26,26,26/)
CASE ('GRAY 11');RGB = (/28,28,28/)
CASE ('GRAY 12');RGB = (/31,31,31/)
CASE ('GRAY 13');RGB = (/33,33,33/)
CASE ('GRAY 14');RGB = (/36,36,36/)
CASE ('GRAY 15');RGB = (/38,38,38/)
CASE ('GRAY 16');RGB = (/41,41,41/)
CASE ('GRAY 17');RGB = (/43,43,43/)
CASE ('GRAY 18');RGB = (/46,46,46/)
CASE ('GRAY 19');RGB = (/48,48,48/)
CASE ('GRAY 2');RGB = (/5,5,5/)
CASE ('GRAY 20');RGB = (/51,51,51/)
CASE ('GRAY 21');RGB = (/54,54,54/)
CASE ('GRAY 22');RGB = (/56,56,56/)
CASE ('GRAY 23');RGB = (/59,59,59/)
CASE ('GRAY 24');RGB = (/61,61,61/)
CASE ('GRAY 25');RGB = (/64,64,64/)
CASE ('GRAY 26');RGB = (/66,66,66/)
CASE ('GRAY 27');RGB = (/69,69,69/)
CASE ('GRAY 28');RGB = (/71,71,71/)
CASE ('GRAY 29');RGB = (/74,74,74/)
CASE ('GRAY 3');RGB = (/8,8,8/)
CASE ('GRAY 30');RGB = (/77,77,77/)
CASE ('GRAY 31');RGB = (/79,79,79/)
CASE ('GRAY 32');RGB = (/82,82,82/)
CASE ('GRAY 33');RGB = (/84,84,84/)
CASE ('GRAY 34');RGB = (/87,87,87/)
CASE ('GRAY 35');RGB = (/89,89,89/)
CASE ('GRAY 36');RGB = (/92,92,92/)
CASE ('GRAY 37');RGB = (/94,94,94/)
CASE ('GRAY 38');RGB = (/97,97,97/)
CASE ('GRAY 39');RGB = (/99,99,99/)
CASE ('GRAY 4');RGB = (/10,10,10/)
CASE ('GRAY 40');RGB = (/102,102,102/)
CASE ('GRAY 42');RGB = (/107,107,107/)
CASE ('GRAY 43');RGB = (/110,110,110/)
CASE ('GRAY 44');RGB = (/112,112,112/)
CASE ('GRAY 45');RGB = (/115,115,115/)
CASE ('GRAY 46');RGB = (/117,117,117/)
CASE ('GRAY 47');RGB = (/120,120,120/)
CASE ('GRAY 48');RGB = (/122,122,122/)
CASE ('GRAY 49');RGB = (/125,125,125/)
CASE ('GRAY 5');RGB = (/13,13,13/)
CASE ('GRAY 50');RGB = (/127,127,127/)
CASE ('GRAY 51');RGB = (/130,130,130/)
CASE ('GRAY 52');RGB = (/133,133,133/)
CASE ('GRAY 53');RGB = (/135,135,135/)
CASE ('GRAY 54');RGB = (/138,138,138/)
CASE ('GRAY 55');RGB = (/140,140,140/)
CASE ('GRAY 56');RGB = (/143,143,143/)
CASE ('GRAY 57');RGB = (/145,145,145/)
CASE ('GRAY 58');RGB = (/148,148,148/)
CASE ('GRAY 59');RGB = (/150,150,150/)
CASE ('GRAY 6');RGB = (/15,15,15/)
CASE ('GRAY 60');RGB = (/153,153,153/)
CASE ('GRAY 61');RGB = (/156,156,156/)
CASE ('GRAY 62');RGB = (/158,158,158/)
CASE ('GRAY 63');RGB = (/161,161,161/)
CASE ('GRAY 64');RGB = (/163,163,163/)
CASE ('GRAY 65');RGB = (/166,166,166/)
CASE ('GRAY 66');RGB = (/168,168,168/)
CASE ('GRAY 67');RGB = (/171,171,171/)
CASE ('GRAY 68');RGB = (/173,173,173/)
CASE ('GRAY 69');RGB = (/176,176,176/)
CASE ('GRAY 7');RGB = (/18,18,18/)
CASE ('GRAY 70');RGB = (/179,179,179/)
CASE ('GRAY 71');RGB = (/181,181,181/)
CASE ('GRAY 72');RGB = (/184,184,184/)
CASE ('GRAY 73');RGB = (/186,186,186/)
CASE ('GRAY 74');RGB = (/189,189,189/)
CASE ('GRAY 75');RGB = (/191,191,191/)
CASE ('GRAY 76');RGB = (/194,194,194/)
CASE ('GRAY 77');RGB = (/196,196,196/)
CASE ('GRAY 78');RGB = (/199,199,199/)
CASE ('GRAY 79');RGB = (/201,201,201/)
CASE ('GRAY 8');RGB = (/20,20,20/)
CASE ('GRAY 80');RGB = (/204,204,204/)
CASE ('GRAY 81');RGB = (/207,207,207/)
CASE ('GRAY 82');RGB = (/209,209,209/)
CASE ('GRAY 83');RGB = (/212,212,212/)
CASE ('GRAY 84');RGB = (/214,214,214/)
CASE ('GRAY 85');RGB = (/217,217,217/)
CASE ('GRAY 86');RGB = (/219,219,219/)
CASE ('GRAY 87');RGB = (/222,222,222/)
CASE ('GRAY 88');RGB = (/224,224,224/)
CASE ('GRAY 89');RGB = (/227,227,227/)
CASE ('GRAY 9');RGB = (/23,23,23/)
CASE ('GRAY 90');RGB = (/229,229,229/)
CASE ('GRAY 91');RGB = (/232,232,232/)
CASE ('GRAY 92');RGB = (/235,235,235/)
CASE ('GRAY 93');RGB = (/237,237,237/)
CASE ('GRAY 94');RGB = (/240,240,240/)
CASE ('GRAY 95');RGB = (/242,242,242/)
CASE ('GRAY 97');RGB = (/247,247,247/)
CASE ('GRAY 98');RGB = (/250,250,250/)
CASE ('GRAY 99');RGB = (/252,252,252/)
CASE ('GREEN');RGB = (/0,255,0/)
CASE ('GREEN 2');RGB = (/0,238,0/)
CASE ('GREEN 3');RGB = (/0,205,0/)
CASE ('GREEN 4');RGB = (/0,139,0/)
CASE ('GREEN YELLOW');RGB = (/173,255,47/)
CASE ('HONEYDEW');RGB = (/240,255,240/)
CASE ('HONEYDEW 1');RGB = (/224,238,224/)
CASE ('HONEYDEW 2');RGB = (/193,205,193/)
CASE ('HONEYDEW 3');RGB = (/131,139,131/)
CASE ('HOT PINK');RGB = (/255,105,180/)
CASE ('HOT PINK 1');RGB = (/255,110,180/)
CASE ('HOT PINK 2');RGB = (/238,106,167/)
CASE ('HOT PINK 3');RGB = (/205,96,144/)
CASE ('HOT PINK 4');RGB = (/139,58,98/)
CASE ('INDIAN RED');RGB = (/205,92,92/)
CASE ('INDIAN RED 1');RGB = (/255,106,106/)
CASE ('INDIAN RED 2');RGB = (/238,99,99/)
CASE ('INDIAN RED 3');RGB = (/205,85,85/)
CASE ('INDIAN RED 4');RGB = (/139,58,58/)
CASE ('INDIGO');RGB = (/75,0,130/)
CASE ('IVORY');RGB = (/255,255,240/)
CASE ('IVORY 1');RGB = (/238,238,224/)
CASE ('IVORY 2');RGB = (/205,205,193/)
CASE ('IVORY 3');RGB = (/139,139,131/)
CASE ('IVORY BLACK');RGB = (/41,36,33/)
CASE ('KELLY GREEN');RGB = (/0,128,0/)
CASE ('KHAKI');RGB = (/240,230,140/)
CASE ('KHAKI 1');RGB = (/255,246,143/)
CASE ('KHAKI 2');RGB = (/238,230,133/)
CASE ('KHAKI 3');RGB = (/205,198,115/)
CASE ('KHAKI 4');RGB = (/139,134,78/)
CASE ('LAVENDER');RGB = (/230,230,250/)
CASE ('LAVENDER BLUSH');RGB = (/255,240,245/)
CASE ('LAVENDER BLUSH 1');RGB = (/238,224,229/)
CASE ('LAVENDER BLUSH 2');RGB = (/205,193,197/)
CASE ('LAVENDER BLUSH 3');RGB = (/139,131,134/)
CASE ('LAWN GREEN');RGB = (/124,252,0/)
CASE ('LEMON CHIFFON');RGB = (/255,250,205/)
CASE ('LEMON CHIFFON 1');RGB = (/238,233,191/)
CASE ('LEMON CHIFFON 2');RGB = (/205,201,165/)
CASE ('LEMON CHIFFON 3');RGB = (/139,137,112/)
CASE ('LIGHT BLUE');RGB = (/173,216,230/)
CASE ('LIGHT BLUE 1');RGB = (/191,239,255/)
CASE ('LIGHT BLUE 2');RGB = (/178,223,238/)
CASE ('LIGHT BLUE 3');RGB = (/154,192,205/)
CASE ('LIGHT BLUE 4');RGB = (/104,131,139/)
CASE ('LIGHT CORAL');RGB = (/240,128,128/)
CASE ('LIGHT CYAN');RGB = (/224,255,255/)
CASE ('LIGHT CYAN 1');RGB = (/209,238,238/)
CASE ('LIGHT CYAN 2');RGB = (/180,205,205/)
CASE ('LIGHT CYAN 3');RGB = (/122,139,139/)
CASE ('LIGHT GOLDENROD');RGB = (/255,236,139/)
CASE ('LIGHT GOLDENROD 1');RGB = (/238,220,130/)
CASE ('LIGHT GOLDENROD 2');RGB = (/205,190,112/)
CASE ('LIGHT GOLDENROD 3');RGB = (/139,129,76/)
CASE ('LIGHT GOLDENROD YELLOW');RGB = (/250,250,210/)
CASE ('LIGHT GREY');RGB = (/211,211,211/)
CASE ('LIGHT PINK');RGB = (/255,182,193/)
CASE ('LIGHT PINK 1');RGB = (/255,174,185/)
CASE ('LIGHT PINK 2');RGB = (/238,162,173/)
CASE ('LIGHT PINK 3');RGB = (/205,140,149/)
CASE ('LIGHT PINK 4');RGB = (/139,95,101/)
CASE ('LIGHT SALMON');RGB = (/255,160,122/)
CASE ('LIGHT SALMON 1');RGB = (/238,149,114/)
CASE ('LIGHT SALMON 2');RGB = (/205,129,98/)
CASE ('LIGHT SALMON 3');RGB = (/139,87,66/)
CASE ('LIGHT SEA GREEN');RGB = (/32,178,170/)
CASE ('LIGHT SKY BLUE');RGB = (/135,206,250/)
CASE ('LIGHT SKY BLUE 1');RGB = (/176,226,255/)
CASE ('LIGHT SKY BLUE 2');RGB = (/164,211,238/)
CASE ('LIGHT SKY BLUE 3');RGB = (/141,182,205/)
CASE ('LIGHT SKY BLUE 4');RGB = (/96,123,139/)
CASE ('LIGHT SLATE BLUE');RGB = (/132,112,255/)
CASE ('LIGHT SLATE GRAY');RGB = (/119,136,153/)
CASE ('LIGHT STEEL BLUE');RGB = (/176,196,222/)
CASE ('LIGHT STEEL BLUE 1');RGB = (/202,225,255/)
CASE ('LIGHT STEEL BLUE 2');RGB = (/188,210,238/)
CASE ('LIGHT STEEL BLUE 3');RGB = (/162,181,205/)
CASE ('LIGHT STEEL BLUE 4');RGB = (/110,123,139/)
CASE ('LIGHT YELLOW 1');RGB = (/255,255,224/)
CASE ('LIGHT YELLOW 2');RGB = (/238,238,209/)
CASE ('LIGHT YELLOW 3');RGB = (/205,205,180/)
CASE ('LIGHT YELLOW 4');RGB = (/139,139,122/)
CASE ('LIME GREEN');RGB = (/50,205,50/)
CASE ('LINEN');RGB = (/250,240,230/)
CASE ('MAGENTA');RGB = (/255,0,255/)
CASE ('MAGENTA 2');RGB = (/238,0,238/)
CASE ('MAGENTA 3');RGB = (/205,0,205/)
CASE ('MAGENTA 4');RGB = (/139,0,139/)
CASE ('MANGANESE BLUE');RGB = (/3,168,158/)
CASE ('MAROON');RGB = (/128,0,0/)
CASE ('MAROON 1');RGB = (/255,52,179/)
CASE ('MAROON 2');RGB = (/238,48,167/)
CASE ('MAROON 3');RGB = (/205,41,144/)
CASE ('MAROON 4');RGB = (/139,28,98/)
CASE ('MEDIUM ORCHID');RGB = (/186,85,211/)
CASE ('MEDIUM ORCHID 1');RGB = (/224,102,255/)
CASE ('MEDIUM ORCHID 2');RGB = (/209,95,238/)
CASE ('MEDIUM ORCHID 3');RGB = (/180,82,205/)
CASE ('MEDIUM ORCHID 4');RGB = (/122,55,139/)
CASE ('MEDIUM PURPLE');RGB = (/147,112,219/)
CASE ('MEDIUM PURPLE 1');RGB = (/171,130,255/)
CASE ('MEDIUM PURPLE 2');RGB = (/159,121,238/)
CASE ('MEDIUM PURPLE 3');RGB = (/137,104,205/)
CASE ('MEDIUM PURPLE 4');RGB = (/93,71,139/)
CASE ('MEDIUM SEA GREEN');RGB = (/60,179,113/)
CASE ('MEDIUM SLATE BLUE');RGB = (/123,104,238/)
CASE ('MEDIUM SPRING GREEN');RGB = (/0,250,154/)
CASE ('MEDIUM TURQUOISE');RGB = (/72,209,204/)
CASE ('MEDIUM VIOLET RED');RGB = (/199,21,133/)
CASE ('MELON');RGB = (/227,168,105/)
CASE ('MIDNIGHT BLUE');RGB = (/25,25,112/)
CASE ('MINT');RGB = (/189,252,201/)
CASE ('MINT CREAM');RGB = (/245,255,250/)
CASE ('MISTY ROSE');RGB = (/255,228,225/)
CASE ('MISTY ROSE 1');RGB = (/238,213,210/)
CASE ('MISTY ROSE 2');RGB = (/205,183,181/)
CASE ('MISTY ROSE 3');RGB = (/139,125,123/)
CASE ('MOCCASIN');RGB = (/255,228,181/)
CASE ('NAVAJO WHITE');RGB = (/255,222,173/)
CASE ('NAVAJO WHITE 1');RGB = (/238,207,161/)
CASE ('NAVAJO WHITE 2');RGB = (/205,179,139/)
CASE ('NAVAJO WHITE 3');RGB = (/139,121,94/)
CASE ('NAVY');RGB = (/0,0,128/)
CASE ('OLD LACE');RGB = (/253,245,230/)
CASE ('OLIVE');RGB = (/128,128,0/)
CASE ('OLIVE DRAB');RGB = (/192,255,62/)
CASE ('OLIVE DRAB 1');RGB = (/179,238,58/)
CASE ('OLIVE DRAB 2');RGB = (/154,205,50/)
CASE ('OLIVE DRAB 3');RGB = (/105,139,34/)
CASE ('ORANGE');RGB = (/255,128,0/)
CASE ('ORANGE 1');RGB = (/255,165,0/)
CASE ('ORANGE 2');RGB = (/238,154,0/)
CASE ('ORANGE 3');RGB = (/205,133,0/)
CASE ('ORANGE 4');RGB = (/139,90,0/)
CASE ('ORANGE RED');RGB = (/255,69,0/)
CASE ('ORANGE RED 1');RGB = (/238,64,0/)
CASE ('ORANGE RED 2');RGB = (/205,55,0/)
CASE ('ORANGE RED 3');RGB = (/139,37,0/)
CASE ('ORCHID');RGB = (/218,112,214/)
CASE ('ORCHID 1');RGB = (/255,131,250/)
CASE ('ORCHID 2');RGB = (/238,122,233/)
CASE ('ORCHID 3');RGB = (/205,105,201/)
CASE ('ORCHID 4');RGB = (/139,71,137/)
CASE ('PALE GOLDENROD');RGB = (/238,232,170/)
CASE ('PALE GREEN');RGB = (/152,251,152/)
CASE ('PALE GREEN 1');RGB = (/154,255,154/)
CASE ('PALE GREEN 2');RGB = (/144,238,144/)
CASE ('PALE GREEN 3');RGB = (/124,205,124/)
CASE ('PALE GREEN 4');RGB = (/84,139,84/)
CASE ('PALE TURQUOISE');RGB = (/187,255,255/)
CASE ('PALE TURQUOISE 1');RGB = (/174,238,238/)
CASE ('PALE TURQUOISE 2');RGB = (/150,205,205/)
CASE ('PALE TURQUOISE 3');RGB = (/102,139,139/)
CASE ('PALE VIOLET RED');RGB = (/219,112,147/)
CASE ('PALE VIOLET RED 1');RGB = (/255,130,171/)
CASE ('PALE VIOLET RED 2');RGB = (/238,121,159/)
CASE ('PALE VIOLET RED 3');RGB = (/205,104,137/)
CASE ('PALE VIOLET RED 4');RGB = (/139,71,93/)
CASE ('PAPAYA WHIP');RGB = (/255,239,213/)
CASE ('PEACH PUFF');RGB = (/255,218,185/)
CASE ('PEACH PUFF 1');RGB = (/238,203,173/)
CASE ('PEACH PUFF 2');RGB = (/205,175,149/)
CASE ('PEACH PUFF 3');RGB = (/139,119,101/)
CASE ('PEACOCK');RGB = (/51,161,201/)
CASE ('PINK');RGB = (/255,192,203/)
CASE ('PINK 1');RGB = (/255,181,197/)
CASE ('PINK 2');RGB = (/238,169,184/)
CASE ('PINK 3');RGB = (/205,145,158/)
CASE ('PINK 4');RGB = (/139,99,108/)
CASE ('PLUM');RGB = (/221,160,221/)
CASE ('PLUM 1');RGB = (/255,187,255/)
CASE ('PLUM 2');RGB = (/238,174,238/)
CASE ('PLUM 3');RGB = (/205,150,205/)
CASE ('PLUM 4');RGB = (/139,102,139/)
CASE ('POWDER BLUE');RGB = (/176,224,230/)
CASE ('PURPLE');RGB = (/128,0,128/)
CASE ('PURPLE 1');RGB = (/155,48,255/)
CASE ('PURPLE 2');RGB = (/145,44,238/)
CASE ('PURPLE 3');RGB = (/125,38,205/)
CASE ('PURPLE 4');RGB = (/85,26,139/)
CASE ('RASPBERRY');RGB = (/135,38,87/)
CASE ('RAW SIENNA');RGB = (/199,97,20/)
CASE ('RED');RGB = (/255,0,0/)
CASE ('RED 1');RGB = (/238,0,0/)
CASE ('RED 2');RGB = (/205,0,0/)
CASE ('RED 3');RGB = (/139,0,0/)
CASE ('ROSY BROWN');RGB = (/188,143,143/)
CASE ('ROSY BROWN 1');RGB = (/255,193,193/)
CASE ('ROSY BROWN 2');RGB = (/238,180,180/)
CASE ('ROSY BROWN 3');RGB = (/205,155,155/)
CASE ('ROSY BROWN 4');RGB = (/139,105,105/)
CASE ('ROYAL BLUE');RGB = (/65,105,225/)
CASE ('ROYAL BLUE 1');RGB = (/72,118,255/)
CASE ('ROYAL BLUE 2');RGB = (/67,110,238/)
CASE ('ROYAL BLUE 3');RGB = (/58,95,205/)
CASE ('ROYAL BLUE 4');RGB = (/39,64,139/)
CASE ('SALMON');RGB = (/250,128,114/)
CASE ('SALMON 1');RGB = (/255,140,105/)
CASE ('SALMON 2');RGB = (/238,130,98/)
CASE ('SALMON 3');RGB = (/205,112,84/)
CASE ('SALMON 4');RGB = (/139,76,57/)
CASE ('SANDY BROWN');RGB = (/244,164,96/)
CASE ('SAP GREEN');RGB = (/48,128,20/)
CASE ('SEA GREEN');RGB = (/84,255,159/)
CASE ('SEA GREEN 1');RGB = (/78,238,148/)
CASE ('SEA GREEN 2');RGB = (/67,205,128/)
CASE ('SEA GREEN 3');RGB = (/46,139,87/)
CASE ('SEASHELL');RGB = (/255,245,238/)
CASE ('SEASHELL 1');RGB = (/238,229,222/)
CASE ('SEASHELL 2');RGB = (/205,197,191/)
CASE ('SEASHELL 3');RGB = (/139,134,130/)
CASE ('SEPIA');RGB = (/94,38,18/)
CASE ('SIENNA');RGB = (/160,82,45/)
CASE ('SIENNA 1');RGB = (/255,130,71/)
CASE ('SIENNA 2');RGB = (/238,121,66/)
CASE ('SIENNA 3');RGB = (/205,104,57/)
CASE ('SIENNA 4');RGB = (/139,71,38/)
CASE ('SILVER');RGB = (/192,192,192/)
CASE ('SKY BLUE');RGB = (/135,206,235/)
CASE ('SKY BLUE 1');RGB = (/135,206,255/)
CASE ('SKY BLUE 2');RGB = (/126,192,238/)
CASE ('SKY BLUE 3');RGB = (/108,166,205/)
CASE ('SKY BLUE 4');RGB = (/74,112,139/)
CASE ('SLATE BLUE');RGB = (/106,90,205/)
CASE ('SLATE BLUE 1');RGB = (/131,111,255/)
CASE ('SLATE BLUE 2');RGB = (/122,103,238/)
CASE ('SLATE BLUE 3');RGB = (/105,89,205/)
CASE ('SLATE BLUE 4');RGB = (/71,60,139/)
CASE ('SLATE GRAY');RGB = (/112,128,144/)
CASE ('SLATE GRAY 1');RGB = (/198,226,255/)
CASE ('SLATE GRAY 2');RGB = (/185,211,238/)
CASE ('SLATE GRAY 3');RGB = (/159,182,205/)
CASE ('SLATE GRAY 4');RGB = (/108,123,139/)
CASE ('SNOW');RGB = (/255,250,250/)
CASE ('SNOW 1');RGB = (/238,233,233/)
CASE ('SNOW 2');RGB = (/205,201,201/)
CASE ('SNOW 3');RGB = (/139,137,137/)
CASE ('SPRING GREEN');RGB = (/0,255,127/)
CASE ('SPRING GREEN 1');RGB = (/0,238,118/)
CASE ('SPRING GREEN 2');RGB = (/0,205,102/)
CASE ('SPRING GREEN 3');RGB = (/0,139,69/)
CASE ('STEEL BLUE');RGB = (/70,130,180/)
CASE ('STEEL BLUE 1');RGB = (/99,184,255/)
CASE ('STEEL BLUE 2');RGB = (/92,172,238/)
CASE ('STEEL BLUE 3');RGB = (/79,148,205/)
CASE ('STEEL BLUE 4');RGB = (/54,100,139/)
CASE ('TAN');RGB = (/210,180,140/)
CASE ('TAN 1');RGB = (/255,165,79/)
CASE ('TAN 2');RGB = (/238,154,73/)
CASE ('TAN 3');RGB = (/205,133,63/)
CASE ('TAN 4');RGB = (/139,90,43/)
CASE ('TEAL');RGB = (/0,128,128/)
CASE ('THISTLE');RGB = (/216,191,216/)
CASE ('THISTLE 1');RGB = (/255,225,255/)
CASE ('THISTLE 2');RGB = (/238,210,238/)
CASE ('THISTLE 3');RGB = (/205,181,205/)
CASE ('THISTLE 4');RGB = (/139,123,139/)
CASE ('TOMATO');RGB = (/255,99,71/)
CASE ('TOMATO 1');RGB = (/238,92,66/)
CASE ('TOMATO 2');RGB = (/205,79,57/)
CASE ('TOMATO 3');RGB = (/139,54,38/)
CASE ('TURQUOISE');RGB = (/64,224,208/)
CASE ('TURQUOISE 1');RGB = (/0,245,255/)
CASE ('TURQUOISE 2');RGB = (/0,229,238/)
CASE ('TURQUOISE 3');RGB = (/0,197,205/)
CASE ('TURQUOISE 4');RGB = (/0,134,139/)
CASE ('TURQUOISE BLUE');RGB = (/0,199,140/)
CASE ('VIOLET');RGB = (/238,130,238/)
CASE ('VIOLET RED');RGB = (/208,32,144/)
CASE ('VIOLET RED 1');RGB = (/255,62,150/)
CASE ('VIOLET RED 2');RGB = (/238,58,140/)
CASE ('VIOLET RED 3');RGB = (/205,50,120/)
CASE ('VIOLET RED 4');RGB = (/139,34,82/)
CASE ('WARM GREY');RGB = (/128,128,105/)
CASE ('WHEAT');RGB = (/245,222,179/)
CASE ('WHEAT 1');RGB = (/255,231,186/)
CASE ('WHEAT 2');RGB = (/238,216,174/)
CASE ('WHEAT 3');RGB = (/205,186,150/)
CASE ('WHEAT 4');RGB = (/139,126,102/)
CASE ('WHITE');RGB = (/255,255,255/)
CASE ('WHITE SMOKE');RGB = (/245,245,245/)
CASE ('YELLOW');RGB = (/255,255,0/)
CASE ('YELLOW 1');RGB = (/238,238,0/)
CASE ('YELLOW 2');RGB = (/205,205,0/)
CASE ('YELLOW 3');RGB = (/139,139,0/)

CASE DEFAULT
  RGB = (/0,0,0/)
END SELECT

END SUBROUTINE COLOR2RGB

!  ------------------ funit ------------------------ 

subroutine get_file_unit(funit,first_unit)
integer, intent(in) :: first_unit
integer, intent(out) :: funit
logical :: is_open

do funit=first_unit,32767
  inquire(UNIT=funit,OPENED=is_open)
  if(is_open)cycle;
  return
end do
funit=-1
return
end subroutine get_file_unit


!  ------------------ GET_TETRABOX_VOLUME ------------------------ 

SUBROUTINE GET_TETRABOX_VOLUME_FB(BOX_BOUNDS_FB,V0_FB,V1_FB,V2_FB,V3_FB,TETRABOX_VOLUME_FB,AREAS_FB,CENTROID_FB)
USE PRECISION_PARAMETERS
USE BOXTETRA_ROUTINES
IMPLICIT NONE

! compute the volume of the intersection of a box and a tetrahedron

REAL(FB), DIMENSION(0:5), INTENT(IN) :: BOX_BOUNDS_FB
REAL(FB), DIMENSION(0:2), INTENT(IN) :: V0_FB, V1_FB, V2_FB, V3_FB
REAL(FB), INTENT(OUT) :: TETRABOX_VOLUME_FB,AREAS_FB(6), CENTROID_FB(3)

REAL(EB), DIMENSION(0:5) :: BOX_BOUNDS_EB
REAL(EB), DIMENSION(0:2) :: V0_EB, V1_EB, V2_EB, V3_EB
REAL(EB) :: TETRABOX_VOLUME_EB,AREAS_EB(6), CENTROID_EB(3)

BOX_BOUNDS_EB(0:5) = REAL(BOX_BOUNDS_FB(0:5),EB)
V0_EB(0:2) = REAL(V0_FB(0:2),FB)
V1_EB(0:2) = REAL(V1_FB(0:2),FB)
V2_EB(0:2) = REAL(V2_FB(0:2),FB)
V3_EB(0:2) = REAL(V3_FB(0:2),FB)

CALL GET_TETRABOX_VOLUME(BOX_BOUNDS_EB,V0_EB,V1_EB,V2_EB,V3_EB,TETRABOX_VOLUME_EB,AREAS_EB,CENTROID_EB)
TETRABOX_VOLUME_FB = REAL(TETRABOX_VOLUME_EB,FB)
AREAS_FB(1:6) = REAL(AREAS_EB(1:6),FB)
CENTROID_FB(1:3) = REAL(CENTROID_EB(1:3),FB)

END SUBROUTINE GET_TETRABOX_VOLUME_FB
!  ------------------ GETVERTS2 ------------------------

SUBROUTINE GETVERTS2(BOX_BOUNDS,V0,V1,V2,V3,VERTS,NVERTS,FACES,FACE_ID,WHICH_POLY,NFACES2,NPOLYS,BOX_STATE)
USE PRECISION_PARAMETERS
USE BOXTETRA_ROUTINES
IMPLICIT NONE

REAL(FB), DIMENSION(0:5), INTENT(IN) :: BOX_BOUNDS
REAL(FB), DIMENSION(0:2), INTENT(IN) :: V0, V1, V2, V3
REAL(FB), DIMENSION(0:599), INTENT(OUT) :: VERTS
INTEGER, DIMENSION(0:199), INTENT(OUT) :: FACE_ID
INTEGER, DIMENSION(0:599), INTENT(OUT) :: FACES
INTEGER, DIMENSION(0:199), INTENT(OUT) :: WHICH_POLY
INTEGER, INTENT(OUT) :: NVERTS,NFACES2,NPOLYS
INTEGER, INTENT(OUT) :: BOX_STATE(-1:5)

REAL(EB) :: VOL_BOX,  VOL_TETRA

REAL(EB), DIMENSION(0:5) :: BOX_BOUNDS_EB, BOX_BOUNDS_IJK_EB
REAL(EB), DIMENSION(0:2) :: V0_EB, V1_EB, V2_EB, V3_EB
REAL(EB),DIMENSION(0:599), TARGET :: VERTS_EB
REAL(EB) :: FACTORI, FACTORJ, FACTORK
REAL(EB) :: FACTORI2, FACTORJ2, FACTORK2

REAL(EB) :: XMIN, YMIN, ZMIN
REAL(EB) :: XMAX, YMAX, ZMAX
INTEGER :: VERT_STATE

REAL(EB) :: VOL

INTEGER :: I, J, K

REAL(EB) :: REL_ERROR
INTEGER :: IBOX_MAX, JBOX_MAX, KBOX_MAX

REAL(EB) :: AREAS_EB(6), CENTROID_EB(3)

V0_EB=REAL(V0,EB)
V1_EB=REAL(V1,EB)
V2_EB=REAL(V2,EB)
V3_EB=REAL(V3,EB)
IF (TETRAHEDRON_VOLUME(V0_EB,V1_EB,V2_EB,V3_EB).LT.0.0_EB) THEN
! swap vertices if orientation is reversed
   V0_EB=REAL(V0,EB)
   V1_EB=REAL(V2,EB)
   V2_EB=REAL(V1,EB)
   V3_EB=REAL(V3,EB)
ENDIF

BOX_BOUNDS_EB=REAL(BOX_BOUNDS,EB)

AREAS_EB(1:6)=0.0_EB
CENTROID_EB(1:3)=0.0_EB

CALL GET_VERTS(BOX_BOUNDS_EB,V0_EB,V1_EB,V2_EB,V3_EB,VERTS_EB,NVERTS,FACES,FACE_ID,WHICH_POLY,NFACES2,NPOLYS,&
   BOX_STATE,VERT_STATE,AREAS_EB)

VERTS(0:3*NVERTS-1) = REAL(VERTS_EB(0:3*NVERTS-1),FB)

XMIN = MIN(V0_EB(0), V1_EB(0), V2_EB(0), V3_EB(0))-0.10_EB
XMAX = MAX(V0_EB(0), V1_EB(0), V2_EB(0), V3_EB(0))+0.10_EB

YMIN = MIN(V0_EB(1), V1_EB(1), V2_EB(1), V3_EB(1))-0.10_EB
YMAX = MAX(V0_EB(1), V1_EB(1), V2_EB(1), V3_EB(1))+0.10_EB

ZMIN = MIN(V0_EB(2), V1_EB(2), V2_EB(2), V3_EB(2))-0.10_EB
ZMAX = MAX(V0_EB(2), V1_EB(2), V2_EB(2), V3_EB(2))+0.10_EB

IBOX_MAX = 5
JBOX_MAX = 5
KBOX_MAX = 5

VOL_BOX=0.0_EB
DO K = 0, KBOX_MAX
   FACTORK =  REAL(  K,EB)/REAL(KBOX_MAX+1,EB)
   FACTORK2 = REAL(K+1,EB)/REAL(KBOX_MAX+1,EB)
   BOX_BOUNDS_IJK_EB(4) = ZMIN*(1.0_EB-FACTORK)  + ZMAX*FACTORK
   BOX_BOUNDS_IJK_EB(5) = ZMIN*(1.0_EB-FACTORK2) + ZMAX*FACTORK2
   DO J = 0, JBOX_MAX
      FACTORJ =  REAL(  J,EB)/REAL(JBOX_MAX+1,EB)
      FACTORJ2 = REAL(J+1,EB)/REAL(JBOX_MAX+1,EB)
      BOX_BOUNDS_IJK_EB(2) = YMIN*(1.0_EB-FACTORJ)  + YMAX*FACTORJ
      BOX_BOUNDS_IJK_EB(3) = YMIN*(1.0_EB-FACTORJ2) + YMAX*FACTORJ2
      DO I = 0, IBOX_MAX
         FACTORI =  REAL(  I,EB)/REAL(IBOX_MAX+1,EB)
         FACTORI2 = REAL(I+1,EB)/REAL(IBOX_MAX+1,EB)
         BOX_BOUNDS_IJK_EB(0) = XMIN*(1.0_EB-FACTORI)  + XMAX*FACTORI
         BOX_BOUNDS_IJK_EB(1) = XMIN*(1.0_EB-FACTORI2) + XMAX*FACTORI2
         CALL GET_TETRABOX_VOLUME(BOX_BOUNDS_IJK_EB,V0_EB,V1_EB,V2_EB,V3_EB,VOL,AREAS_EB,CENTROID_EB)
         VOL_BOX = VOL_BOX + VOL
      END DO
   END DO
END DO
VOL_TETRA=ABS(TETRAHEDRON_VOLUME(V0_EB,V1_EB,V2_EB,V3_EB))
REL_ERROR = (VOL_BOX-VOL_TETRA)/VOL_TETRA
IF (ABS(REL_ERROR)<1.0E-14_EB) THEN
   WRITE(LU_OUTPUT,'(a,E15.8,a,e15.8)')"passed tetra=",VOL_TETRA," rel error=",REL_ERROR
ELSE
   WRITE(LU_OUTPUT,'(a,E15.8,a,e15.8)')"*** failure *** tetra=",VOL_TETRA," rel error=",REL_ERROR
ENDIF

RETURN
END SUBROUTINE GETVERTS2

!  ------------------ TEST_IN_TETRA ------------------------

SUBROUTINE TEST_IN_TETRA(XYZ, IN_TETRA, TETRA_STATE)
USE PRECISION_PARAMETERS
USE BOXTETRA_ROUTINES
IMPLICIT NONE

REAL(EB), INTENT(IN) :: XYZ(3)
INTEGER, INTENT(OUT) :: IN_TETRA, TETRA_STATE(4)

CALL TEST_IN_TETRA0(XYZ, IN_TETRA, TETRA_STATE)

END SUBROUTINE TEST_IN_TETRA


