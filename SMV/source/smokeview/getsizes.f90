! $Date: 2015-04-29 15:16:22 -0400 (Wed, 29 Apr 2015) $ 
! $Revision: 22563 $
! $Author: gforney $

!  ------------------ getembeddatasize ------------------------ 

subroutine getembeddatasize(filename,ntimes,nvars,error)
implicit none
character(len=*), intent(in) :: filename
integer, intent(out) :: ntimes, nvars, error

integer :: lu20, finish
logical :: isopen,exists
real :: time, dummy
integer :: i, one, version
integer :: nvert_s, nvert_d, nface_s, nface_d

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
ntimes=0
nvars=0
do 
  read(lu20,iostat=finish)time
  if(finish.eq.0)read(lu20,iostat=finish)nvert_s,nface_s,nvert_d,nface_d
  if(finish.eq.0.and.nvert_s>0)read(lu20,iostat=finish)(dummy,i=1,nvert_s)
  if(finish.eq.0.and.nvert_d>0)read(lu20,iostat=finish)(dummy,i=1,nvert_d)
  if(finish.eq.0.and.nface_s>0)read(lu20,iostat=finish)(dummy,i=1,nface_s)
  if(finish.eq.0.and.nface_d>0)read(lu20,iostat=finish)(dummy,i=1,nface_d)
  if(finish.ne.0)return
  nvars = nvars + nvert_s + nvert_d + nface_s + nface_d
  ntimes = ntimes + 1
end do
close(lu20)

end subroutine getembeddatasize

!  ------------------ endian_open ------------------------ 

integer function endian_open(file,lunit)
character(len=*), intent(in) :: file
integer, intent(in) :: lunit

logical isopen,exists
integer :: error
integer :: one

error=0
inquire(unit=lunit,opened=isopen)
if(isopen)close(lunit)

inquire(file=trim(file),exist=exists)
if(.not.exists)then
  endian_open=1
  return
endif
#ifdef pp_SHARED
open(unit=lunit,file=trim(file),form="unformatted",shared,action="read")
#else
open(unit=lunit,file=trim(file),form="unformatted",action="read")
#endif
read(lunit)one
if(one.eq.1)then
  endian_open=0
  rewind(lunit)
 else
  endian_open=1
endif
return
end function endian_open

!  ------------------ fcreate_part5sizefile ------------------------ 

subroutine fcreate_part5sizefile(part5file, part5sizefile, angle_flag, redirect_flag, error)
implicit none
character(len=*), intent(in) :: part5file, part5sizefile
integer, intent(in) :: angle_flag, redirect_flag
integer, intent(out) :: error

integer :: lu20, lu21, version, nclasses
logical :: isopen
integer, allocatable, dimension(:) :: numtypes, numpoints
character(len=30) :: dummy
integer :: i, j,one,idummy
real :: rdummy,time
integer :: endian_open

error=1
lu20=20

error=endian_open(trim(part5file),lu20)
if(error.ne.0)return

lu21=21
inquire(unit=lu21,opened=isopen)
if(isopen)close(lu21)
open(unit=lu21,file=trim(part5sizefile),form="formatted",action="write")

error=0
read(lu20,iostat=error)one
if(error.ne.0)go to 998
read(lu20,iostat=error)version
if(error.ne.0)go to 998
read(lu20,iostat=error)nclasses
if(error.ne.0)go to 998
allocate(numtypes(2*nclasses))
allocate(numpoints(nclasses))
do i = 1, nclasses
  read(lu20,iostat=error)numtypes(2*i-1),numtypes(2*i)
  if(error.ne.0)go to 999
  do j = 1, (numtypes(2*i-1)+numtypes(2*i))
    read(lu20,iostat=error)dummy
    if(error.ne.0)go to 999
    read(lu20,iostat=error)dummy
    if(error.ne.0)go to 999    
  end do
end do
do
  read(lu20,iostat=error)time
  if(redirect_flag.eq.0)write(6,10)time
10 format(" sizing particle time=",f9.2)
  if(error.ne.0)go to 999
  do i = 1, nclasses
    read(lu20,iostat=error)numpoints(i)
    if(error.ne.0)go to 999    
    if(angle_flag.eq.1)then
      read(lu20,iostat=error)(rdummy,j=1,7*numpoints(i))
     else
      read(lu20,iostat=error)(rdummy,j=1,3*numpoints(i))
    endif
    if(error.ne.0)go to 999    
    read(lu20,iostat=error)(rdummy,j=1,numpoints(i))
    if(error.ne.0)go to 999    
    if(numtypes(2*i-1)>0)then
      read(lu20,iostat=error)(rdummy,j=1,numpoints(i)*numtypes(2*i-1))
      if(error.ne.0)go to 999    
    endif
    if(numtypes(2*i)>0)then
      read(lu20,iostat=error)(idummy,j=1,numpoints(i)*numtypes(2*i))
      if(error.ne.0)go to 999    
    endif    
  end do
  write(lu21,"(e15.8)")time
  do i = 1, nclasses
    write(lu21,"(1x,i9)")numpoints(i)
  end do
end do
!     for(i=0;i<nclasses;i++){
!      FORTPART5READ(&nparts,1);
!      if(returncode==0)goto wrapup;
!      numpoints[i]=nparts;
!      skip = 4 + 4*nparts*3 + 4;
!      skip += 4 + 4*nparts + 4;
!      if(numtypes[2*i]>0)skip += 4 + 4*nparts*numtypes[2*i] + 4;
!      if(numtypes[2*i+1]>0)skip += 4 + 4*nparts*numtypes[2*i+1] + 4;
!      
!      returncode=fseek(PART5FILE,skip,SEEK_CUR);
!      if(returncode!=0)goto wrapup;
!    }


999 continue
deallocate(numpoints,numtypes)
998 continue
close(lu20)
close(lu21)

return
end subroutine fcreate_part5sizefile

!  ------------------ getzonesize ------------------------ 

subroutine getzonesize(zonefilename,nzonet,nrooms,nfires,error)
implicit none
character(len=*) :: zonefilename
integer, intent(out) :: nzonet,nrooms,nfires,error

logical :: isopen, exists
integer :: lu26, version
integer :: i
real :: dummy, dummy2
integer :: exit_all
integer :: file_unit

call get_file_unit(file_unit,26)
lu26 = file_unit
error = 0

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

nzonet = 0
read(lu26,iostat=error)version
if(error.eq.0)read(lu26,iostat=error)nrooms
if(error.eq.0)read(lu26,iostat=error)nfires
if(error.ne.0)then
  error=0
  rewind(lu26)
  return
endif
do
  exit_all=0
  read(lu26,iostat=error)dummy
  if(error.ne.0)then
    error = 0
    exit
  endif
  do i = 1, nrooms
    read(lu26,iostat=error)dummy,dummy,dummy,dummy
    if(error.eq.0)cycle
    error = 0
    exit_all=1
    exit
  end do 
  if(exit_all.eq.1)exit
  do i = 1, nfires
    read(lu26,iostat=error)dummy,dummy2
    if(error.eq.0)cycle
    error = 0
    exit_all=1
    exit
  end do
  if(exit_all.eq.1)exit
  nzonet = nzonet + 1
end do
close(lu26)
end subroutine getzonesize

!  ------------------ getpatchsizes1 ------------------------ 

subroutine getpatchsizes1(file_unit,patchfilename,patchlonglabel,patchshortlabel,patchunit, &
       npatch,headersize,error)
implicit none

character(len=*) :: patchfilename, patchlonglabel, patchshortlabel, patchunit
integer, intent(in) :: file_unit
integer, intent(out) :: npatch
integer, intent(out) :: headersize

integer, intent(out) :: error
integer :: lu15, lenshort, lenunits
logical :: exists
logical :: isopen

error=0
lu15 = file_unit
inquire(unit=lu15,opened=isopen)

if(isopen)close(lu15)
inquire(file=trim(patchfilename),exist=exists)
if(exists)then
#ifdef pp_SHARED
  open(unit=lu15,file=trim(patchfilename),form="unformatted",shared,action="read")
#else
  open(unit=lu15,file=trim(patchfilename),form="unformatted",action="read")
#endif
 else
  write(6,*)' The boundary file name, ',trim(patchfilename),' does not exist'
  error=1
  return
endif

if(error.eq.0)read(lu15,iostat=error)patchlonglabel
if(error.eq.0)read(lu15,iostat=error)patchshortlabel
if(error.eq.0)read(lu15,iostat=error)patchunit
if(error.eq.0)read(lu15,iostat=error)npatch
headersize = 3*(30+8) + 4 + 8 

patchlonglabel=trim(patchlonglabel)//char(0)
lenshort = min(len_trim(patchshortlabel),6)
patchshortlabel=patchshortlabel(1:lenshort)//char(0)
lenunits = min(len_trim(patchunit),6)
patchunit=patchunit(1:lenunits)//char(0)

return
end subroutine getpatchsizes1

!  ------------------ getpatchsizes2 ------------------------ 

subroutine getpatchsizes2(file_unit,version,npatch,npatchsize,pi1,pi2,pj1,pj2,pk1,pk2,patchdir,headersize,framesize)
implicit none
integer, intent(in) :: version, npatch,file_unit
integer, intent(out) :: npatchsize
integer, intent(out), dimension(npatch) :: pi1, pi2, pj1, pj2, pk1, pk2, patchdir
integer, intent(inout) :: headersize
integer, intent(out) :: framesize

integer :: n, lu15
integer :: i1, i2, j1, j2, k1, k2

lu15 = file_unit
npatchsize = 0
do n = 1, npatch
  if(version.eq.0)then
    read(lu15)i1, i2, j1, j2, k1, k2
   else
    read(lu15)i1, i2, j1, j2, k1, k2, patchdir(n)
  endif
  pi1(n)=i1
  pi2(n)=i2
  pj1(n)=j1
  pj2(n)=j2
  pk1(n)=k1
  pk2(n)=k2
  npatchsize = npatchsize + (i2+1-i1)*(j2+1-j1)*(k2+1-k1)
end do
headersize = headersize + npatch*(8+6*4)
if(version.eq.1)headersize = headersize + npatch*4
framesize = 8+4+8*npatch+npatchsize*4

return
end subroutine getpatchsizes2

!  ------------------ getsizesa ------------------------ 
!    FORTgetsizesa(file,&npartpoint,&npartframes,lenfile);

subroutine getsizesa(partfilename,npartpoints,npartframes)
implicit none

integer, intent(out) :: npartpoints, npartframes
character(len=*) :: partfilename

integer :: lu10
integer :: error
logical :: connected,exists
character(len=255) :: line
integer :: nparts
integer :: nlines
integer :: i
real :: dummy

lu10 = 10
error=0
npartpoints=0
npartframes=0
inquire(unit=lu10,opened=connected)
if(connected)close(lu10)

inquire(file=trim(partfilename),exist=exists)
if(exists)then
#ifdef pp_SHARED
  open(unit=lu10,file=trim(partfilename),form="formatted",shared,action="read")
#else
  open(unit=lu10,file=trim(partfilename),form="formatted",action="read")
#endif
 else
  write(6,*)' The particle file name, ',trim(partfilename),' does not exist'
  error=-1
endif
do i = 1, 5
  read(lu10,'(a)')line
end do
do 
  read(lu10,*,end=999)dummy,nparts
  nlines = (nparts-1)*3/6 + 1
  do i = 1, nlines
    read(lu10,'(a)',end=999)line
  end do
  npartpoints = npartpoints + nparts
  npartframes = npartframes + 1
end do
999 continue
close(lu10)
return
end subroutine getsizesa

!  ------------------ getsizes ------------------------ 

subroutine getsizes(file_unit,partfilename,nb,nv,nspr,mxframepoints,showstaticsmoke, error)
implicit none

integer, intent(out) :: nb, nv, nspr, mxframepoints, showstaticsmoke, error
integer, intent(in) :: file_unit

integer :: lu10, ii, i, j, k
real :: xx, yy
character(len=*) :: partfilename
logical :: exists, connected

integer :: ibar1, jbar1, kbar1

lu10 = file_unit
error=0
inquire(unit=lu10,opened=connected)
if(connected)close(lu10)

inquire(file=trim(partfilename),exist=exists)
if(exists)then
#ifdef pp_SHARED
  open(unit=lu10,file=trim(partfilename),form="unformatted",shared,action="read")
#else
  open(unit=lu10,file=trim(partfilename),form="unformatted",action="read")
#endif
 else
  write(6,*)' The particle file name, ',trim(partfilename),' does not exist'
  error=-1
endif
read(lu10,iostat=error) xx,xx,yy,ii,mxframepoints
showstaticsmoke=0
if(yy.lt.0)showstaticsmoke=1
if(error.ne.0)return

read(lu10,iostat=error) ibar1,jbar1,kbar1
if(error.ne.0)return

read(lu10,iostat=error) (xx,i=0,ibar1),(xx,j=0,jbar1),(xx,k=0,kbar1)
if(error.ne.0)return

read(lu10,iostat=error) nb
if(error.ne.0)return

do i=1,nb
  read(lu10,iostat=error) ii,ii,ii,ii,ii,ii,ii
  if(error.ne.0)return
end do

read(lu10,iostat=error) nv
if(error.ne.0)return

do i=1,nv
  read(lu10,iostat=error) ii,ii,ii,ii,ii,ii,ii
  if(error.ne.0)return
end do

read(lu10,iostat=error) nspr
if(error.ne.0)return
rewind(lu10)

return
end subroutine getsizes

!  ------------------ getsizes2 ------------------------ 

subroutine getsizes2(file_unit,settmin_p,tmin_p,settmax_p,tmax_p,nspr,frameloadstep,partpointstep,npartpoints,npartframes,error)
implicit none

integer, intent(in) :: file_unit,settmin_p, settmax_p, nspr
real, intent(in) :: tmin_p, tmax_p
integer, intent(in) :: frameloadstep, partpointstep
integer, intent(out) :: npartpoints, npartframes, error

integer :: lu10
integer :: idummy
real :: dummy
integer, allocatable, dimension(:) :: ispr
integer :: naspr
integer :: i
integer :: npoints, npoints2
integer :: npp1, npp2
real :: stime
integer :: nf
logical :: load

lu10 = file_unit

npartframes = 0
npartpoints = 0
if(nspr.gt.0)allocate(ispr(nspr))

nf = 0
do
  read(lu10,iostat=error) stime,npp1,idummy,(ispr(i),i=1,nspr)
  if(error.ne.0)go to 999
  naspr = 0
  do i = 1, nspr
    if(ispr(i).ne.0)naspr = naspr + 1
  end do

  read(lu10,iostat=error) (dummy,i=1,npp1),(dummy,i=1,npp1),(dummy,i=1,npp1),(dummy,i=1,npp1)
  if(error.ne.0)go to 999
  if((settmin_p.ne.0.and.stime.lt.tmin_p).or.mod(nf,frameloadstep).ne.0)then
    load=.false.
   else
    load=.true.
  endif
  nf = nf + 1
  if(settmax_p.ne.0.and.stime.gt.tmax_p)go to 999
  npoints = 0
  if(load)npoints = (npp1-1)/partpointstep + 1

  npp2 = 0

  npoints2 = 0
  if(naspr.ne.0)then       ! read in sprinkler data
    read(lu10,iostat=error) npp2
    if(error.ne.0)go to 999
    if(npp2.ge.0)then
      read(lu10,iostat=error) (dummy,i=1,npp2),(dummy,i=1,npp2),(dummy,i=1,npp2)
     else
      read(lu10,iostat=error) (dummy,i=1,npp2),(dummy,i=1,npp2),(dummy,i=1,npp2),(dummy,i=1,npp2)
    endif
    if(load)npoints2 = (abs(npp2)-1)/partpointstep + 1
  end if
  if(error.ne.0)goto 999
  if(load)then
    npartframes = npartframes + 1
    npartpoints = npartpoints + npoints + npoints2
  endif

end do
999 continue
error = 0
close(lu10)
return
end subroutine getsizes2

!  ------------------ getsliceparms ------------------------ 

subroutine getsliceparms(slicefilename, ip1, ip2, jp1, jp2, kp1, kp2, ni, nj, nk, slice3d, error)
implicit none

character(len=*) :: slicefilename
logical :: exists

integer, intent(inout) :: ip1, ip2, jp1, jp2, kp1, kp2
integer, intent(out) :: ni, nj, nk, slice3d, error

integer :: idir, joff, koff
logical :: connected
character(len=30) :: longlbl, shortlbl, unitlbl

integer :: lu11

if(ip1.eq.-1.or.ip2.eq.-1.or.jp1.eq.-1.or.jp2.eq.-1.or.kp1.eq.-1.or.kp2.eq.-1)then
  ip1 = 0
  ip2 = 0
  jp1 = 0
  jp2 = 0
  kp1 = 0
  kp2 = 0
  error=0
  lu11 = 11
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
    error=1
    return
  endif
  read(lu11,iostat=error)longlbl
  read(lu11,iostat=error)shortlbl
  read(lu11,iostat=error)unitlbl

  read(lu11,iostat=error)ip1, ip2, jp1, jp2, kp1, kp2
  close(lu11)
endif

ni = ip2 + 1 - ip1
nj = jp2 + 1 - jp1
nk = kp2 + 1 - kp1
if(ip1.eq.ip2.or.jp1.eq.jp2.or.kp1.eq.kp2)then
  slice3d=0
 else
  slice3d=1
endif

call getdirval(ip1,ip2,jp1,jp2,kp1,kp2,idir,joff,koff)

return
end subroutine getsliceparms

!  ------------------ getslicesizes ------------------------ 

subroutine getslicesizes(slicefilename, nslicei, nslicej, nslicek, nsteps, sliceframestep,&
   error, settmin_s, settmax_s, tmin_s, tmax_s, headersize, framesize)
implicit none

character(len=*) :: slicefilename
logical :: exists

integer, intent(out) :: nslicei, nslicej, nslicek, nsteps, error
integer, intent(in) :: settmin_s, settmax_s, sliceframestep
integer, intent(out) :: headersize, framesize
real, intent(in) :: tmin_s, tmax_s

integer :: ip1, ip2, jp1, jp2, kp1, kp2
integer :: nxsp, nysp, nzsp
integer :: i, j, k

integer :: lu11
real :: time, time_max
real, dimension(:,:,:), pointer :: qq
character(len=30) :: longlbl, shortlbl, unitlbl
logical :: connected, load
integer :: idir, joff, koff
integer :: count

error=0
lu11 = 11
nsteps = 0 
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
  error=1
  return
endif

headersize = 0
read(lu11,iostat=error)longlbl
read(lu11,iostat=error)shortlbl
read(lu11,iostat=error)unitlbl
headersize = headersize + 3*38

read(lu11,iostat=error)ip1, ip2, jp1, jp2, kp1, kp2
headersize = headersize + 6*4 + 8
if(error.ne.0)return
  
nxsp = ip2 + 1 - ip1
nysp = jp2 + 1 - jp1
nzsp = kp2 + 1 - kp1

allocate(qq(nxsp,nysp,nzsp))

call getdirval(ip1,ip2,jp1,jp2,kp1,kp2,idir,joff,koff)
nslicei = nxsp
nslicej = nysp + joff
nslicek = nzsp + koff

framesize = 4*(1+nxsp*nysp*nzsp)+16

count=-1
time_max=-1000000.0
do
  read(lu11,iostat=error)time
  if(error.ne.0)exit
  if((settmin_s.ne.0.and.time.lt.tmin_s).or.time.le.time_max)then
    load=.false.
   else
    load = .true.
    time_max=time
  endif
  if(settmax_s.ne.0.and.time.gt.tmax_s)then
    close(lu11)
    return
  endif
  read(lu11,iostat=error)(((qq(i,j,k),i=1,nxsp),j=1,nysp),k=1,nzsp)
  count = count + 1
  if(mod(count,sliceframestep).ne.0)load = .false.
  if(error.ne.0)exit
  if(load)nsteps = nsteps + 1
end do

error = 0
close(lu11)

return

end subroutine getslicesizes

!  ------------------ openpart ------------------------ 

subroutine openpart(partfilename, unit, error)
implicit none

character(len=*) :: partfilename
logical :: exists

integer, intent(in) :: unit
integer, intent(out) :: error
logical :: connected

integer :: lu11

error=0
lu11 = unit
inquire(unit=lu11,opened=connected)
if(connected)close(lu11)

inquire(file=partfilename,exist=exists)
if(exists)then
#ifdef pp_SHARED
  open(unit=lu11,file=partfilename,form="unformatted",shared,action="read")
#else
  open(unit=lu11,file=partfilename,form="unformatted",action="read")
#endif
 else
  error=1
  return
endif

return
end subroutine openpart

!  ------------------ openslice ------------------------ 

subroutine openslice(slicefilename, unitnum, is1, is2, js1, js2, ks1, ks2, error)
implicit none

character(len=*) :: slicefilename
logical :: exists

integer, intent(inout) :: unitnum
integer, intent(out) :: is1, is2, js1, js2, ks1, ks2
integer, intent(out) :: error
character(len=30) :: longlbl, shortlbl, unitlbl

integer :: lu11

error=0
lu11 = unitnum;
exists=.true.

inquire(file=slicefilename,exist=exists)
if(exists)then
#ifdef pp_SHARED
  open(unit=lu11,file=slicefilename,form="unformatted",shared,action="read")
#else
  open(unit=lu11,file=slicefilename,form="unformatted",action="read")
#endif
 else
  error=1
  return
endif
read(lu11,iostat=error)longlbl
read(lu11,iostat=error)shortlbl
read(lu11,iostat=error)unitlbl

read(lu11,iostat=error)is1, is2, js1, js2, ks1, ks2

return
end subroutine openslice

!  ------------------ closefortranfile ------------------------ 

subroutine closefortranfile(unit)
implicit none

integer, intent(in) :: unit

close(unit)

return
end subroutine closefortranfile

!  ------------------ getboundaryheader1 ------------------------ 

subroutine getboundaryheader1(boundaryfilename,boundaryunitnumber,npatch,error)
implicit none

character(len=*), intent(in) :: boundaryfilename
integer, intent(inout) :: boundaryunitnumber
integer, intent(out) :: npatch, error

character(len=30) :: patchlonglabel, patchshortlabel, patchunit

integer :: lu15
logical :: exists
logical :: isopen

error=0
call get_file_unit(boundaryunitnumber,boundaryunitnumber)
lu15 = boundaryunitnumber
inquire(unit=lu15,opened=isopen)

if(isopen)close(lu15)
inquire(file=trim(boundaryfilename),exist=exists)
if(exists)then
#ifdef pp_SHARED
  open(unit=lu15,file=trim(boundaryfilename),form="unformatted",shared,action="read")
#else
  open(unit=lu15,file=trim(boundaryfilename),form="unformatted",action="read")
#endif
 else
  write(6,*)' The boundary file name, ',trim(boundaryfilename),' does not exist'
  error=1
  return
endif

if(error.eq.0)read(lu15,iostat=error)patchlonglabel
if(error.eq.0)read(lu15,iostat=error)patchshortlabel
if(error.eq.0)read(lu15,iostat=error)patchunit
if(error.eq.0)read(lu15,iostat=error)npatch
if(error.ne.0)close(lu15)

return
end subroutine getboundaryheader1

!  ------------------ getboundaryheader2 ------------------------ 

subroutine getboundaryheader2(boundaryunitnumber,version,npatch,pi1,pi2,pj1,pj2,pk1,pk2,patchdir)
implicit none
integer, intent(in) :: boundaryunitnumber, version, npatch
integer, intent(out), dimension(npatch) :: pi1, pi2, pj1, pj2, pk1, pk2, patchdir

integer :: n, lu15
integer :: i1, i2, j1, j2, k1, k2

lu15 = boundaryunitnumber
do n = 1, npatch
  if(version.eq.0)then
    read(lu15)i1, i2, j1, j2, k1, k2
   else
    read(lu15)i1, i2, j1, j2, k1, k2, patchdir(n)
  endif
  pi1(n)=i1
  pi2(n)=i2
  pj1(n)=j1
  pj2(n)=j2
  pk1(n)=k1
  pk2(n)=k2
end do
close(lu15)

return
end subroutine getboundaryheader2

!  ------------------ openboundary ------------------------ 

subroutine openboundary(boundaryfilename,boundaryunitnumber,version,error)
implicit none

character(len=*), intent(in) :: boundaryfilename
integer, intent(in) :: boundaryunitnumber,version
integer, intent(out) :: error

character(len=30) :: patchlonglabel, patchshortlabel, patchunit

integer :: lu15
logical :: exists
logical :: isopen
integer :: npatch,n
integer :: i1, i2, j1, j2, k1, k2, patchdir

error=0
lu15 = boundaryunitnumber
inquire(unit=lu15,opened=isopen)

if(isopen)close(lu15)
inquire(file=boundaryfilename,exist=exists)
if(exists)then
#ifdef pp_SHARED
  open(unit=lu15,file=boundaryfilename,form="unformatted",shared,action="read")
#else
  open(unit=lu15,file=boundaryfilename,form="unformatted",action="read")
#endif
 else
  write(6,*)' The boundary file name, ',boundaryfilename,' does not exist'
  error=1
  return
endif

if(error.eq.0)read(lu15,iostat=error)patchlonglabel
if(error.eq.0)read(lu15,iostat=error)patchshortlabel
if(error.eq.0)read(lu15,iostat=error)patchunit
if(error.eq.0)read(lu15,iostat=error)npatch

do n = 1, npatch
  if(version.eq.0)then
    if(error.eq.0)read(lu15,iostat=error)i1, i2, j1, j2, k1, k2
   else
    if(error.eq.0)read(lu15,iostat=error)i1, i2, j1, j2, k1, k2, patchdir
  endif
end do

if(error.ne.0)close(lu15)

return
end subroutine openboundary

!  ------------------ getpartheader1 ------------------------ 

subroutine getpartheader1(unit,nclasses,fdsversion,size)
implicit none

integer, intent(in) :: unit
integer, intent(out) :: nclasses,fdsversion,size

integer :: one

read(unit)one
read(unit)fdsversion

read(unit)nclasses
size=12

return

end subroutine getpartheader1

!  ------------------ getpartheader2 ------------------------ 

subroutine getpartheader2(unit,nclasses,nquantities,size)
implicit none

integer, intent(in) :: unit,nclasses
integer, intent(out), dimension(nclasses) :: nquantities
integer, intent(out) :: size

character(len=30) :: clabel
integer :: i, j, dummy

size=0

do i = 1, nclasses
  read(unit)nquantities(i),dummy
  size=size+4+2*nquantities(i)*(4+30+4)
  do j=1, nquantities(i)
    read(unit)clabel
    read(unit)clabel
  end do
end do

return

end subroutine getpartheader2

!  ------------------ getpartdataframe ------------------------ 

subroutine getpartdataframe(unit,nclasses,nquantities,npoints,time,tagdata,pdata,size,error)
implicit none

integer, intent(in) :: unit,nclasses
integer, intent(in), dimension(nclasses) :: nquantities
integer, intent(out), dimension(nclasses) :: npoints
real, intent(out), dimension(*) :: pdata
integer, intent(out), dimension(*) :: tagdata
real, intent(out) :: time
integer, intent(out) :: size,error

integer :: pstart, pend
integer :: tagstart, tagend
integer :: i, j, nparticles

size=0
pend=0
tagend=0
error=0
read(unit,iostat=error)time
size=4
if(error.ne.0)return
do i = 1, nclasses
  read(unit,iostat=error)nparticles
  if(error.ne.0)return
  npoints(i)=nparticles
  
  pstart=pend+1
  pend=pstart+3*nparticles-1
  read(unit,iostat=error)(pdata(j),j=pstart,pend)
  if(error.ne.0)return

  tagstart = tagend + 1
  tagend = tagstart + nparticles - 1
  read(unit,iostat=error)(tagdata(j),j=tagstart,tagend)
  if(error.ne.0)return

  if(nquantities(i).gt.0)then
    pstart = pend + 1
    pend = pstart + nparticles*nquantities(i) - 1
    read(unit,iostat=error)(pdata(j),j=pstart,pend)
    if(error.ne.0)return
  endif
  size=size+4+(4*3*nparticles)+4*nparticles+4*nparticles*nquantities(i)
end do
error=0

end subroutine getpartdataframe
