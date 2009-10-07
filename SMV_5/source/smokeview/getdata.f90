!  ------------------ getzonedata ------------------------ 

subroutine getzonedata(nzonet,nrooms, nfires, zonet,zoneqfire,zonepr, zoneylay,zonetl,zonetu,error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getzonedata@40' :: getzonedata
#endif
implicit none
integer, intent(in) :: nrooms, nfires
integer, intent(inout) :: nzonet
real, intent(out), dimension(nrooms*nzonet) :: zonepr, zoneylay, zonetl, zonetu
real, intent(out), dimension(nfires*nzonet) :: zoneqfire
real, intent(out), dimension(nzonet) :: zonet
integer , intent(out) :: error

integer :: lu26,i,j,ii,ii2,idummy,version
real :: dummy, qdot

lu26 = 26

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

!  ------------------ getxyzdata ------------------------ 

subroutine getxyzdata(iblank,nx,ny,nz,error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getxyzdata@20' :: getxyzdata
#endif
implicit none

integer :: lu24
real :: dummyx, dummyy, dummyz
integer, intent(in) :: nx, ny, nz
integer, dimension(nx,ny,nz) :: iblank
integer, intent(out) :: error
integer :: i,j,k

error=0
lu24=24
read(LU24,iostat=error)(((dummyx,I=1,nx),J=1,ny),K=1,nz),&
(((dummyy,I=1,nx),J=1,ny),K=1,nz),          &
(((dummyz,I=1,nx),J=1,ny),K=1,nz),          &
(((iblank(i,j,k),i=1,nx),j=1,ny),k=1,nz)
close(lu24)
return
end subroutine getxyzdata

!  ------------------ getpatchdata ------------------------ 

subroutine getpatchdata(npatch,pi1,pi2,pj1,pj2,pk1,pk2,patchtime,pqq,error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getpatchdata@40' :: getpatchdata
#endif
implicit none


integer, intent(in) :: npatch
integer, intent(in), dimension(*) :: pi1, pi2, pj1, pj2, pk1, pk2
real, intent(out), dimension(*) :: pqq
integer, intent(out) :: error
real, intent(out) :: patchtime

integer :: i, i1, i2, j1, j2, k1, k2, size, ibeg, iend, lu15, ii

lu15 = 15
read(lu15,iostat=error)patchtime
if(error.ne.0)then
  close(lu15)
  return
endif
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

subroutine getdata1(ipart,error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getdata1@8' :: getdata1
#endif

implicit none

integer, intent(out) :: ipart, error

integer :: lu10
real :: sarx, sary, swpar
integer :: i, j, k
integer ndum2
integer :: nspr, nv
integer :: ibar, jbar, kbar
real :: dummy
integer :: nb1, idummy

lu10 = 10
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

!  ------------------ getdata2a ------------------------ 

subroutine getdata2a(nmax,nspr,x,y,z,t,stime,np,ns,error)
                   
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getdata2a@40' :: getdata2a
#endif

  implicit none

  integer, intent(in) :: nspr,nmax
  real, dimension(*), intent(out) :: x, y, z, t
  real, intent(out) ::  stime
  integer, intent(out) :: np, ns, error
  
  integer :: lu10, nn, i, naspr
  integer, dimension(:), allocatable :: ispr

  lu10 = 10

  if(nspr.gt.0)allocate(ispr(nspr))

  read(lu10,iostat=error) stime,np,nn,(ispr(i),i=1,nspr)
  if(error.ne.0)go to 999
  naspr = 0
  do i = 1, nspr
    if(ispr(i).ne.0)naspr = naspr + 1
  end do

  read(lu10,iostat=error) (x(i),i=1,np),(y(i),i=1,np),(z(i),i=1,np),(t(i),i=1,np)
  if(error.ne.0)go to 999

  ns = 0
  if(naspr.ne.0)then       ! read in sprinkler data
    read(lu10,iostat=error) ns
    if(error.ne.0)go to 999
	if(np+ns.gt.nmax)then
	  error=1
	  go to 999
	endif
    read(lu10,iostat=error) (x(i),i=np+1,np+ns),(y(i),i=np+1,np+ns),(z(i),i=np+1,np+ns)
    if(error.ne.0)go to 999
  end if

999 continue
if(nspr.gt.0)deallocate(ispr)
close(lu10)
return
end subroutine getdata2a

!  ------------------ getdata2 ------------------------ 

!STDCALL FORTgetdata2b(short *xs, short *yparts, short *zparts, float *tpart,
!                       int *bframe,int *sframe,float *stimes,int *npartpoints,int *mxframes, int *nframes,
!                      float *xbar0, float *xbar, float *ybar0, float *ybar, float *zbar0, float *zbar, int *error);
subroutine getdata2b(partfilename, x, y, z, tpart, bframe, sframe, stimes, &
                     npartpoints, npartframes, mxpoints, mxframes, error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getdata2b@56' :: getdata2b
#endif
implicit none
character(len=*), intent(in) :: partfilename
real, dimension(*), intent(out) :: x, y, z
integer, dimension(*), intent(out) :: bframe, sframe
real, intent(out), dimension(*) :: tpart, stimes
integer, intent(in) :: mxframes, mxpoints
integer, intent(out) :: npartpoints
integer, intent(out) :: npartframes,error

integer :: lu10, i, nparts, nlines, nxyzvals
logical :: connected, opened, exists
character(len=255) :: line
real :: t_or_x
real, allocatable, dimension(:) :: xyzvals
integer :: jbeg, jend, j
integer :: ipoint

lu10 = 10
error=0
npartpoints=0
npartframes=0
inquire(unit=lu10,opened=connected)
if(connected)close(lu10)

inquire(file=trim(partfilename),exist=exists)
if(exists)then
#ifdef pp_cvf
  open(unit=lu10,file=trim(partfilename),form="formatted",action="read",shared)
#else
  open(unit=lu10,file=trim(partfilename),form="formatted",action="read")
#endif
 else
  write(6,*)'The particle file name, ',trim(partfilename),' does not exist'
  error=-1
endif
do i = 1, 5
  read(lu10,'(a)')line
end do
nxyzvals=10
allocate(xyzvals(nxyzvals))
ipoint=0
do 
  read(lu10,*,end=999)t_or_x,nparts
  write(6,*)"frame: ",ipoint+1,"t =",t_or_x," npoints=",nparts
  nlines = (nparts-1)*3/6 + 1
  if(6*nlines>nxyzvals)then
    nxyzvals=6*nlines
    deallocate(xyzvals)
    allocate(xyzvals(nxyzvals))
  endif
  do i = 1, nlines
    jbeg=6*(i-1) + 1
    jend = jbeg + 5
    if(jend.gt.nparts*3)jend=nparts*3
    read(lu10,*,end=999)(xyzvals(j),j=jbeg,jend)
  end do
  do i = 1, nparts
    ipoint = ipoint + 1
    x(ipoint) = t_or_x
    y(ipoint) = xyzvals(i)
    z(ipoint) = xyzvals(nparts+i)
    tpart(ipoint) = xyzvals(2*nparts+i)
  end do
  npartpoints = npartpoints + nparts
  npartframes = npartframes + 1
  sframe(npartframes) = nparts
  stimes(npartframes) = t_or_x
end do
999 continue
  bframe(1)=0
  do i=2,npartframes
    bframe(i) = bframe(i-1)+sframe(i-1)
  end do

end subroutine getdata2b

!  ------------------ getdata2 ------------------------ 

subroutine getdata2(xs,ys,zs,&
                    t,&
                    sprinkflag,isprink,tspr,bframe,sframe,sprframe,stimes,nspr,nmax,mxframes,nframes,&
                    settmin_p,settmax_p,tmin_p,tmax_p,frameloadstep,partpointstep, &
              			xbox0, xbox, ybox0, ybox, zbox0, zbox, &
                    offset_x, offset_y, offset_z, &
      	    				error)
                   
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getdata2@128' :: getdata2
#endif

implicit none
real, dimension(*), intent(out) :: t

integer(2), dimension(*) :: xs, ys, zs
integer, dimension(*) :: bframe, sframe, sprframe
character(len=1), dimension(*) :: isprink
real, dimension(*) ::  stimes,tspr
integer, intent(in) :: nspr,nmax, mxframes
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

lu10 = 10
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
  	  write(6,10)stime
10    format("particle time=",f9.2)
  	 else
      write(6,*)"particle time=",stime,"particles",npp1,"droplets",npp2
      write(6,20)stime,npp2
20    format("particle time=",f9.2," particles",i9," droplets",i9)      
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

!  ------------------ getslicedata ------------------------ 

subroutine getslicedata(slicefilename,longlabel,shortlabel,units,&
            is1,is2,js1,js2,ks1,ks2,idir,qmin,qmax,qdata,times,nstepsmax,sliceframestep,&
			endian,settmin_s,settmax_s,tmin_s,tmax_s)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getslicedata@104' :: getslicedata
#endif

implicit none

character(len=*) :: slicefilename, longlabel, shortlabel,units

real, intent(out) :: qmin, qmax
real, intent(out), dimension(*) :: qdata
real, intent(out), dimension(*) :: times
integer, intent(out) :: idir
integer, intent(out) :: is1, is2, js1, js2, ks1, ks2
integer, intent(in) :: endian
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

lu11 = 11
joff = 0
koff = 0
inquire(unit=lu11,opened=connected)
if(connected)close(lu11)

inquire(file=trim(slicefilename),exist=exists)
if(exists)then
#ifdef pp_cvf
if(endian.eq.1)then
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read",shared,convert="BIG_ENDIAN")
 else
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read",shared)
endif
#else
#ifdef pp_LAHEY
if(endian.eq.1)then
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read",convert="BIG_ENDIAN")
 else
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read")
endif
#else
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read")
#endif
#endif
 else
  write(6,*)'the slice file ',trim(slicefilename),' does not exist'
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
  write(6,10)time
10 format("slice time=",f9.2)  
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
!!    istart = (nsteps-1)*(nysp+joff)*(nzsp+koff)
    istart = (nsteps-1)*(nysp+joff)*(nzsp+koff)*nxsp
    do i = 1, nxsp
    do j = 1, nysp+joff
!!      irowstart = (j-1)*(nzsp+koff)
      irowstart = (i-1)*nysp*(nzsp+koff)+(j-1)*(nzsp+koff)
      kk = istart + irowstart
!!      qdata(kk+1:kk+nzsp+koff) = qq(1,j,1:nzsp+koff)
!!      qmax = max(qmax,maxval(qq(1,j,1:nzsp+koff)))
!!      qmin = min(qmin,minval(qq(1,j,1:nzsp+koff)))
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

subroutine getsliceframe(lu11,is1,is2,js1,js2,ks1,ks2,time,qframe,error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getslicedata@40' :: getslicedata
#endif

implicit none

real, intent(out), dimension(*) :: qframe
real, intent(out) :: time
integer, intent(out) :: error
integer, intent(in) :: lu11, is1, is2, js1, js2, ks1, ks2

integer :: i,j,k
integer :: nxsp, nysp, nzsp

nxsp = is2 + 1 - is1
nysp = js2 + 1 - js1
nzsp = ks2 + 1 - ks1  

read(lu11,iostat=error)time
if(error.ne.0)return
read(lu11,iostat=error)(((qframe(1+i+j*nxsp+k*nxsp*nysp),i=0,nxsp-1),j=0,nysp-1),k=0,nzsp-1)

999 continue

return
end subroutine getsliceframe

!  ------------------ outsliceheader ------------------------ 

subroutine outsliceheader(slicefilename,unit,ip1, ip2, jp1, jp2, kp1, kp2, error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_outsliceheader@40' :: outsliceheader
#endif

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

longlbl= "                              "
shortlbl="                              "
unitlbl= "                              "
write(lu11,iostat=error)longlbl
write(lu11,iostat=error)shortlbl
write(lu11,iostat=error)unitlbl

write(lu11,iostat=error)ip1, ip2, jp1, jp2, kp1, kp2

end subroutine outsliceheader

!  ------------------ getsliceframe ------------------------ 

subroutine outsliceframe(lu11,is1,is2,js1,js2,ks1,ks2,time,qframe,error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_outsliceframe@40' :: outsliceframe
#endif

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

!  ------------------ getplot3dqa ------------------------ 

subroutine getplot3dqa(qfilename,nx,ny,nz,qq,error)

#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getplot3dqa@28' :: getplot3dqa
#endif
implicit none

character(len=*) :: qfilename
integer, intent(in) :: nx, ny, nz
integer, intent(out) :: error
real, dimension(nx,ny,nz,5)  :: qq
real, dimension(ny*nz) :: qbuffer

integer :: u_in, error2, i, j, k, kbeg, kend, nlines
logical :: connected, exists
character(len=255) :: line

  u_in=10
  inquire(unit=u_in,opened=connected)
  if(connected)close(u_in)
  

  error=0
  inquire(file=trim(qfilename),exist=exists)
  if(exists)then
#ifdef pp_cvf
    open(unit=u_in,file=trim(qfilename),form="formatted",action="read",shared,iostat=error2)
#else
    open(unit=u_in,file=trim(qfilename),form="formatted",action="read",iostat=error2)
#endif
   else
    write(6,*)'The file name, ',trim(qfilename),' does not exist'
    return
  endif

do i = 1, 4
  read(u_in,'(a)')line
end do
nlines = (ny*nz-1)/6 + 1
do i = 1, nx
  read(u_in,'(a)')line
  do j = 1, nlines
    kbeg = 6*(j-1) + 1
    kend = kbeg + 5
    if(kend.gt.ny*nz)kend=ny*nz
    read(u_in,*)(qbuffer(k),k=kbeg,kend)
  end do
  do j = 1, ny
    do k = 1, nz
      qq(i,j,k,1:5) = qbuffer(j+(k-1)*ny)
    end do
  end do
end do
return 

end subroutine getplot3dqa

!  ------------------ getplot3dq ------------------------ 

subroutine getplot3dq(qfilename,nx,ny,nz,qq,error,endian,isotest)

#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getplot3dq@36' :: getplot3dq
#endif
implicit none

character(len=*) :: qfilename
integer, intent(in) :: nx, ny, nz
integer, intent(out) :: error
real, dimension(nx,ny,nz,5)  :: qq
integer, intent(in) :: endian, isotest

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
  u_in=10
  inquire(unit=u_in,opened=connected)
  if(connected)close(u_in)

  error=0
  inquire(file=trim(qfilename),exist=exists)
  if(exists)then
#ifdef pp_cvf
  if(endian.eq.1)then
    open(unit=u_in,file=trim(qfilename),form="unformatted",action="read",shared,iostat=error2,convert="BIG_ENDIAN")
   else
    open(unit=u_in,file=trim(qfilename),form="unformatted",action="read",shared,iostat=error2)
  endif
#else
#ifdef pp_LAHEY
  if(endian.eq.1)then
    open(unit=u_in,file=trim(qfilename),form="unformatted",action="read",iostat=error2,convert="BIG_ENDIAN")
   else
    open(unit=u_in,file=trim(qfilename),form="unformatted",action="read",iostat=error2)
  endif
#else
    open(unit=u_in,file=trim(qfilename),form="unformatted",action="read",iostat=error2)
#endif
#endif
   else
    write(6,*)'The file name, ',trim(qfilename),' does not exist'
    read(5,*)dummy
    stop
  endif

  read(u_in,iostat=error)nxpts, nypts, nzpts
  if(nx.eq.nxpts.and.ny.eq.nypts.and.nz.eq.nzpts)then
    read(u_in,iostat=error)dum1, dum2, dum3, dum4
    read(u_in,iostat=error)((((qq(i,j,k,n),i=1,nxpts),j=1,nypts),k=1,nzpts),n=1,5)
   else
    error = 1
	  write(6,*)"*** Fatal error in getplot3dq ***"
  	write(6,*)"Grid size found in plot3d file was:",nxpts,nypts,nzpts
  	write(6,*)"Was expecting:",nx,ny,nz
  	stop
  endif
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
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_plot3dout@28' :: plot3dout
#endif
implicit none

character(len=*), intent(in) :: outfile
real, dimension(nx,ny,nz,5)  :: qout
integer, intent(in) :: nx, ny, nz
integer, intent(out) :: error3

integer :: u_out
logical :: connected
integer :: i, j, k, n
real :: dummy

error3 = 0

u_out=13
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

