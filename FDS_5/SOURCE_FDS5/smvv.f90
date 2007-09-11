module isodefs
implicit none

CHARACTER(255), PARAMETER :: smvvid='$Id$'
CHARACTER(255), PARAMETER :: smvvrev='$Revision$'
CHARACTER(255), PARAMETER :: smvvdate='$Date$'

interface

subroutine isoheader(isofile,isolonglabel,isoshortlabel,isounits,levels,nlevels,error)

!DEC$ ATTRIBUTES C :: ISOHEADER
!DEC$ ATTRIBUTES REFERENCE :: ISOFILE,ISOLONGLABEL,ISOSHORTLABEL
!DEC$ ATTRIBUTES REFERENCE :: ISOUNITS,LEVELS,NLEVELS,ERROR

character(len=*) :: isofile
character(len=30), intent(in) :: isolonglabel, isoshortlabel, isounits
integer, intent(in) :: nlevels
integer, intent(out) :: error
real, dimension(nlevels), intent(in) :: levels

end subroutine isoheader

subroutine tisoheader(isofile,isolonglabel,isoshortlabel,isounits,levels,nlevels,error)

!DEC$ ATTRIBUTES C :: TISOHEADER
!DEC$ ATTRIBUTES REFERENCE :: ISOFILE,ISOLONGLABEL,ISOSHORTLABEL
!DEC$ ATTRIBUTES REFERENCE :: ISOUNITS,LEVELS,NLEVELS,ERROR

character(len=*) :: isofile
character(len=30), intent(in) :: isolonglabel, isoshortlabel, isounits
integer, intent(in) :: nlevels
integer, intent(out) :: error
real, dimension(nlevels), intent(in) :: levels

end subroutine tisoheader

subroutine iso2file(isofile,t,data,iblank,level,nlevels, xplt, nx, yplt, ny, zplt, nz, isooffset, reduce_triangles, error)

!DEC$ ATTRIBUTES C :: ISO2FILE
!DEC$ ATTRIBUTES REFERENCE :: ISOFILE, T, DATA, IBLANK, LEVEL
!DEC$ ATTRIBUTES REFERENCE :: NLEVELS, XPLT, NX, YPLT, NY, ZPLT, NZ 
!DEC$ ATTRIBUTES REFERENCE :: ISOOFFSET, REDUCE_TRIANGLES, ERROR

character(len=*), intent(in) :: isofile
integer, intent(in) :: nlevels, isooffset, nx, ny, nz
integer, intent(out) :: error
real, intent(in) :: t
real, dimension(nlevels), intent(in)  :: level
real, intent(in), dimension(nx*ny*nz) :: data
integer, intent(in), dimension(nx*ny*nz) :: iblank
real, intent(in), dimension(nx) :: xplt
real, intent(in), dimension(ny) :: yplt
real, intent(in), dimension(nz) :: zplt
integer, intent(in) :: reduce_triangles

end subroutine iso2file

subroutine isot2file(isofile,t,data,data2flag,data2, iblank,level,nlevels, &
                      xplt, nx, yplt, ny, zplt, nz, isooffset, reduce_triangles, error)

!DEC$ ATTRIBUTES C :: ISOT2FILE
!DEC$ ATTRIBUTES REFERENCE :: ISOFILE, T, DATA, DATA2FLAG, DATA2
!DEC$ ATTRIBUTES REFERENCE :: IBLANK, LEVEL, NLEVELS, XPLT, NX
!DEC$ ATTRIBUTES REFERENCE :: YPLT, NY, ZPLT, NZ, ISOOFFSET
!DEC$ ATTRIBUTES REFERENCE :: REDUCE_TRIANGLES, ERROR

character(len=*), intent(in) :: isofile
integer, intent(in) :: nlevels, isooffset, nx, ny, nz, data2flag
integer, intent(out) :: error
real, intent(in) :: t
real, dimension(nlevels), intent(in)  :: level
real, intent(in), dimension(nx*ny*nz) :: data, data2
integer, intent(in), dimension(nx*ny*nz) :: iblank
real, intent(in), dimension(nx) :: xplt
real, intent(in), dimension(ny) :: yplt
real, intent(in), dimension(nz) :: zplt
integer, intent(in) :: reduce_triangles

end subroutine isot2file

end interface

CONTAINS

SUBROUTINE GET_REV_smvv(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,*) smvvrev(INDEX(smvvrev,':')+1:LEN_TRIM(smvvrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,*) smvvdate

END SUBROUTINE GET_REV_smvv

end module isodefs

module compressdefs

interface

subroutine smoke3dheader(file,is1,is2,js1,js2,ks1,ks2)

!DEC$ ATTRIBUTES C :: SMOKE3DHEADER
!DEC$ ATTRIBUTES REFERENCE :: FILE,IS1,IS2,JS1,JS2,KS1,KS2
character(len=*), intent(in) :: file
integer, intent(in) ::is1,is2,js1,js2,ks1,ks2

end subroutine smoke3dheader

subroutine smoke3dtofile(file,time,dx,extcoef,type,xyz,nx,ny,nz)

!DEC$ ATTRIBUTES C :: SMOKE3DTOFILE
!DEC$ ATTRIBUTES REFERENCE :: FILE,TIME,DX,EXTCOEF,TYPE,XYZ,NX,NY,NZ
character(len=*), intent(in) :: file
real, intent(in) :: time, dx, extcoef
integer, intent(in) :: nx,ny,nz
real, intent(in), dimension(nx*ny*nz) :: xyz
integer, intent(in) :: type
end subroutine smoke3dtofile

end interface
end module compressdefs
