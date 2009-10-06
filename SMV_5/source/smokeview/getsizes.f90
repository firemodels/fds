
!  ------------------ endian_open ------------------------ 

integer function endian_open(file,lunit,endian)
character(len=*), intent(in) :: file
integer, intent(in) :: lunit, endian

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
#ifdef pp_cvf
if(endian.eq.0)then
  open(unit=lunit,file=trim(file),form="unformatted",action="read",shared,convert="BIG_ENDIAN")
 else
  open(unit=lunit,file=trim(file),form="unformatted",action="read",shared)
endif
#elif pp_LAHEY
if(endian.eq.0)then
  open(unit=lunit,file=trim(file),form="unformatted",action="read",convert="BIG_ENDIAN")
 else
  open(unit=lunit,file=trim(file),form="unformatted",action="read")
endif
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

subroutine fcreate_part5sizefile(part5file, part5sizefile, angle_flag, error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_fcreate_part5sizefile@24' :: fcreate_part5sizefile
#endif
implicit none
character(len=*), intent(in) :: part5file, part5sizefile
integer, intent(in) :: angle_flag
integer, intent(out) :: error

integer :: lu20, lu21, version, nclasses
logical :: isopen, exists
integer, allocatable, dimension(:) :: numtypes, numpoints
character(len=30) :: dummy
integer :: i, j,one,idummy
real :: rdummy,time
integer :: endian
integer :: endian_open

error=1
lu20=20

error=endian_open(trim(part5file),lu20,1)
if(error.ne.0)then
  error=endian_open(trim(part5file),lu20,0)
endif
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
  write(lu21,"(e11.4)")time
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

!  ------------------ getisosize ------------------------ 

subroutine getisosize(isofilename,endian,nisopoints,nisotriangles,nisolevels,nisosteps,error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getisosize@32' :: getisosize
#endif
implicit none
character(len=*), intent(in) :: isofilename
integer, intent(in) :: endian
integer, intent(out) :: nisopoints, nisotriangles, nisolevels, nisosteps, error

integer :: endian2, lu20
logical :: isopen,exists
character(len=30) :: label
integer :: version,i,j, nisopoints_i, nisotriangles_i
real :: dummy
integer :: idummy,finish
integer(2) :: idummy2
character(len=1) :: idummy1


lu20=20
nisopoints=0
nisotriangles=0
nisosteps=0
inquire(unit=lu20,opened=isopen)

if(isopen)close(lu20)
inquire(file=trim(isofilename),exist=exists)
if(exists)then
#ifdef pp_cvf
endian2=0
endian2=endian
if(endian2.eq.1)then
  open(unit=lu20,file=trim(isofilename),form="unformatted",action="read",shared,convert="BIG_ENDIAN")
 else
  open(unit=lu20,file=trim(isofilename),form="unformatted",action="read",shared)
endif
#elif pp_LAHEY
endian2=0
endian2=endian
if(endian2.eq.1)then
  open(unit=lu20,file=trim(isofilename),form="unformatted",action="read",convert="BIG_ENDIAN")
 else
  open(unit=lu20,file=trim(isofilename),form="unformatted",action="read")
endif
#else	   
  open(unit=lu20,file=trim(isofilename),form="unformatted",action="read")
#endif
 else
  write(6,*)'The iso file name, ',trim(isofilename),' does not exist'
  error=1
  return
endif

error = 0
read(lu20)version
read(lu20)label
read(lu20)label
read(lu20)label
read(lu20)nisolevels
read(lu20)(dummy,i=1,nisolevels)

do 
  read(lu20,iostat=finish)dummy
  if(finish.ne.0)return
  do i = 1, nisolevels
    read(lu20,iostat=finish)nisopoints_i,nisotriangles_i
	if(finish.ne.0)return
	if(nisopoints_i.lt.256)then
      read(lu20,iostat=finish)(idummy2,j=1,3*nisopoints_i),(idummy1,j=1,nisotriangles_i)
	 elseif(nisopoints_i.ge.256.and.nisopoints_i.lt.32768)then
      read(lu20,iostat=finish)(idummy2,j=1,3*nisopoints_i),(idummy2,j=1,nisotriangles_i)
	 else
      read(lu20,iostat=finish)(idummy2,j=1,3*nisopoints_i),(idummy,j=1,nisotriangles_i)
	endif
	if(finish.ne.0)return
	nisopoints = nisopoints + nisopoints_i
	nisotriangles = nisotriangles + nisotriangles_i
  end do
  nisosteps=nisosteps+1
end do

end subroutine getisosize

!  ------------------ getzonesize ------------------------ 

subroutine getzonesize(zonefilename,nzonet,nrooms,nfires,endian,error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getzonesize@28' :: getzonesize
#endif
implicit none
character(len=*) :: zonefilename
integer, intent(out) :: nzonet,nrooms,nfires,error
integer, intent(in) :: endian

logical :: isopen, exists
integer :: lu26, version
integer :: i, endian2
real :: dummy, dummy2

lu26 = 26
error = 0

inquire(unit=lu26,opened=isopen)

if(isopen)close(lu26)
inquire(file=trim(zonefilename),exist=exists)
if(exists)then
#ifdef pp_cvf
endian2=0
endian2=endian
if(endian2.eq.1)then
  open(unit=lu26,file=trim(zonefilename),form="unformatted",action="read",shared,convert="BIG_ENDIAN")
 else
  open(unit=lu26,file=trim(zonefilename),form="unformatted",action="read",shared)
endif
#elif pp_LAHEY
endian2=0
endian2=endian
if(endian2.eq.1)then
  open(unit=lu26,file=trim(zonefilename),form="unformatted",action="read",convert="BIG_ENDIAN")
 else
  open(unit=lu26,file=trim(zonefilename),form="unformatted",action="read")
endif
#else
  open(unit=lu26,file=trim(zonefilename),form="unformatted",action="read")
#endif
 else
  write(6,*)'The zone file name, ',trim(zonefilename),' does not exist'
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
  read(lu26,iostat=error)dummy
  if(error.ne.0)then
    error = 0
    rewind(lu26)
    return
  endif
  do i = 1, nrooms
    read(lu26,iostat=error)dummy,dummy,dummy,dummy
    if(error.ne.0)then
      error = 0
      rewind(lu26)
      return
    endif
  end do 
  do i = 1, nfires
    read(lu26,iostat=error)dummy,dummy2
    if(error.ne.0)then
      error=0
      rewind(lu26)
      return
    endif
  end do
  nzonet = nzonet + 1
end do

end subroutine getzonesize


!  ------------------ getxyzsize ------------------------ 

subroutine getxyzsize(xyzfilename,ibp1,jbp1,kbp1,endian,error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getxyzsize@28' :: getxyzsize
#endif
implicit none

character(len=*) :: xyzfilename
integer, intent(in) :: endian
integer, intent(out) :: ibp1, jbp1, kbp1, error

logical isopen, exists
integer :: lu24

error=0
lu24 = 24

inquire(unit=lu24,opened=isopen)

if(isopen)close(lu24)
inquire(file=trim(xyzfilename),exist=exists)
if(exists)then
#ifdef pp_cvf
if(endian.eq.1)then
  open(unit=lu24,file=trim(xyzfilename),form="unformatted",action="read",shared,convert="BIG_ENDIAN")
 else
  open(unit=lu24,file=trim(xyzfilename),form="unformatted",action="read",shared)
endif
#elif pp_LAHEY
if(endian.eq.1)then
  open(unit=lu24,file=trim(xyzfilename),form="unformatted",action="read",convert="BIG_ENDIAN")
 else
  open(unit=lu24,file=trim(xyzfilename),form="unformatted",action="read")
endif
#else
  open(unit=lu24,file=trim(xyzfilename),form="unformatted",action="read")
#endif
 else
  write(6,*)'The xyz file name, ',trim(xyzfilename),' does not exist'
  error=1
  return
endif

read(LU24,iostat=error) IBP1,JBP1,KBP1
end subroutine getxyzsize


!  ------------------ closezone ------------------------ 

subroutine closezone()
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_closezone@0' :: closezone
#endif
implicit none

integer :: lu26
logical :: isopen

lu26 = 26
inquire(unit=lu26,opened=isopen)

if(isopen)close(lu26)

return
end subroutine closezone

!  ------------------ closeiso ------------------------ 

subroutine closeiso()
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_closeiso@0' :: closeiso
#endif
implicit none

integer :: lu20
logical :: isopen

lu20 = 20
inquire(unit=lu20,opened=isopen)

if(isopen)close(lu20)

return
end subroutine closeiso

!  ------------------ closepatch ------------------------ 

subroutine closepatch()
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_closepatch@0' :: closepatch
#endif
implicit none

integer :: lu15
logical :: isopen

lu15 = 15
inquire(unit=lu15,opened=isopen)

if(isopen)close(lu15)

return
end subroutine closepatch

!  ------------------ closepart ------------------------ 

subroutine closepart()
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_closepart@0' :: closepart
#endif
implicit none

integer :: lu10
logical :: isopen

lu10 = 10
inquire(unit=lu10,opened=isopen)

if(isopen)close(lu10)

return
end subroutine closepart

!  ------------------ getpatchsizes1 ------------------------ 

subroutine getpatchsizes1(patchfilename,patchlonglabel,patchshortlabel,patchunit, &
       endian,npatch,headersize,error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getpatchsizes1@48' :: getpatchsizes1
#endif
implicit none

character(len=*) :: patchfilename, patchlonglabel, patchshortlabel, patchunit
integer, intent(in) :: endian
integer, intent(out) :: npatch
integer, intent(out) :: headersize

integer, intent(out) :: error
integer :: lu15, lenshort, lenunits
logical :: exists
logical :: isopen

error=0
lu15 = 15
inquire(unit=lu15,opened=isopen)

if(isopen)close(lu15)
inquire(file=trim(patchfilename),exist=exists)
if(exists)then
#ifdef pp_cvf
if(endian.eq.1)then
  open(unit=lu15,file=trim(patchfilename),form="unformatted",action="read",shared,convert="BIG_ENDIAN")
 else
  open(unit=lu15,file=trim(patchfilename),form="unformatted",action="read",shared)
endif
#elif pp_LAHEY
if(endian.eq.1)then
  open(unit=lu15,file=trim(patchfilename),form="unformatted",action="read",convert="BIG_ENDIAN")
 else
  open(unit=lu15,file=trim(patchfilename),form="unformatted",action="read")
endif
#else
  open(unit=lu15,file=trim(patchfilename),form="unformatted",action="read")
#endif
 else
  write(6,*)'The boundary file name, ',trim(patchfilename),' does not exist'
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

subroutine getpatchsizes2(version,npatch,npatchsize,pi1,pi2,pj1,pj2,pk1,pk2,patchdir,headersize,framesize)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getpatchsizes2@48' :: getpatchsizes2
#endif
implicit none
integer, intent(in) :: version, npatch
integer, intent(out) :: npatchsize
integer, intent(out), dimension(npatch) :: pi1, pi2, pj1, pj2, pk1, pk2, patchdir
integer, intent(inout) :: headersize
integer, intent(out) :: framesize

integer :: n, lu15
integer :: i1, i2, j1, j2, k1, k2

lu15 = 15
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
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getsizesa@16' :: getsizesa
#endif
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

subroutine getsizes(partfilename,ibar,jbar,kbar,nb,nv,nspr,mxframepoints,endian, showstaticsmoke, error)
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getsizes@48' :: getsizes
#endif
implicit none

integer, intent(out) :: nb, nv, nspr, mxframepoints, showstaticsmoke, error
integer, intent(in) :: ibar, jbar, kbar

integer :: lu10, ii, i, j, k
integer, intent(in) :: endian
real :: xx, yy
character(len=*) :: partfilename
logical :: exists, connected

integer :: ibar1, jbar1, kbar1


lu10 = 10
error=0
inquire(unit=lu10,opened=connected)
if(connected)close(lu10)

inquire(file=trim(partfilename),exist=exists)
if(exists)then
#ifdef pp_cvf
if(endian.eq.1)then
  open(unit=lu10,file=trim(partfilename),form="unformatted",action="read",convert="BIG_ENDIAN",shared)
 else
  open(unit=lu10,file=trim(partfilename),form="unformatted",action="read",shared)
endif
#elif pp_LAHEY
if(endian.eq.1)then
  open(unit=lu10,file=trim(partfilename),form="unformatted",action="read",convert="BIG_ENDIAN")
 else
  open(unit=lu10,file=trim(partfilename),form="unformatted",action="read")
endif
#else
  open(unit=lu10,file=trim(partfilename),form="unformatted",action="read")
#endif
 else
  write(6,*)'The particle file name, ',trim(partfilename),' does not exist'
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

subroutine getsizes2(settmin_p,tmin_p,settmax_p,tmax_p,nspr,frameloadstep,partpointstep,npartpoints,npartframes,error)
                   
#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getsizes2@40' :: getsizes2
#endif

implicit none

integer, intent(in) :: settmin_p, settmax_p, nspr
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

lu10 = 10

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

subroutine getsliceparms(slicefilename, endian, ip1, ip2, jp1, jp2, kp1, kp2, slice3d, error)

#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getsliceparms@44' :: getsliceparms
#endif

implicit none

character(len=*) :: slicefilename
logical :: exists

integer, intent(out) :: ip1, ip2, jp1, jp2, kp1, kp2, slice3d, error
integer, intent(in) :: endian
integer :: idir, joff, koff
logical :: connected
character(len=30) :: longlbl, shortlbl, unitlbl


integer :: lu11

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
#ifdef pp_cvf
if(endian.eq.1)then
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read",shared,convert="BIG_ENDIAN")
 else
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read",shared)
endif
#elif pp_LAHEY
if(endian.eq.1)then
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read",convert="BIG_ENDIAN")
 else
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read")
endif
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

if(ip1.eq.ip2.or.jp1.eq.jp2.or.kp1.eq.kp2)then
  slice3d=0
 else
  slice3d=1
endif

call getdirval(ip1,ip2,jp1,jp2,kp1,kp2,idir,joff,koff)
close(lu11)

return
end subroutine getsliceparms

!  ------------------ getslicesizes ------------------------ 

subroutine getslicesizes(slicefilename, nslicei, nslicej, nslicek, nsteps, sliceframestep,&
   endian, error, settmin_s, settmax_s, tmin_s, tmax_s, &
   headersize, framesize, statfile)

#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_getslicesizes@64' :: getslicesizes
#endif

implicit none

character(len=*) :: slicefilename
logical :: exists

integer, intent(out) :: nslicei, nslicej, nslicek, nsteps, error
integer, intent(in) :: endian
integer, intent(in) :: settmin_s, settmax_s, statfile , sliceframestep
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
#ifdef pp_cvf
if(endian.eq.1)then
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read",shared,convert="BIG_ENDIAN")
 else
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read",shared)
endif
#elif pp_LAHEY
if(endian.eq.1)then
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read",convert="BIG_ENDIAN")
 else
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read")
endif
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
if(settmin_s.ne.0.or.settmax_s.ne.0.or.statfile.ne.0)then
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
endif

error = 0
close(lu11)


return

end subroutine getslicesizes

!  ------------------ openslice ------------------------ 

subroutine openslice(slicefilename, unit, endian, is1, is2, js1, js2, ks1, ks2, error)

#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_openslice@44' :: openslice
#endif

implicit none

character(len=*) :: slicefilename
logical :: exists

integer, intent(in) :: endian, unit
integer, intent(out) :: is1, is2, js1, js2, ks1, ks2
integer, intent(out) :: error
logical :: connected
character(len=30) :: longlbl, shortlbl, unitlbl


integer :: lu11

error=0
lu11 = unit
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
#elif pp_LAHEY
if(endian.eq.1)then
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read",convert="BIG_ENDIAN")
 else
  open(unit=lu11,file=trim(slicefilename),form="unformatted",action="read")
endif
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

read(lu11,iostat=error)is1, is2, js1, js2, ks1, ks2

return
end subroutine openslice

!  ------------------ closeslice ------------------------ 

subroutine closefortranfile(unit)

#ifdef pp_cvf
!DEC$ ATTRIBUTES ALIAS:'_closefortranfile@4' :: closefortranfile
#endif

implicit none

integer, intent(in) :: unit

close(unit)

return
end subroutine closefortranfile