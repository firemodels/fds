program triangle_test

! fire dynamics simulator, main program, multiple cpu version.

use precision_parameters
use complex_geometry, only: triangulate, valid_triangle

implicit none

integer :: i

! miscellaneous declarations
!subroutine triangulate(verts,nverts,vert_offset,faces)
!  integer, intent(in) :: nverts, vert_offset
!  real(fb), intent(in) :: verts(3*nverts)
!  integer, intent(out) :: faces(3*(nverts-2))

integer :: nverts
real(fb) :: verts(300)
integer :: vert_offset
integer :: faces(300)
logical :: tritest
integer :: iv1, iv2, iv3, dir

vert_offset=0
nverts=6

verts(1:3*nverts) = (/2.0, 0.0, 0.0, 0.0,0.0,0.0, 0.0,2.0,0.0, 1.0,2.0,0.0,  1.0,1.0,0.0, 2.0, 1.0, 0.0    /)
dir=3

do i = 1, nverts 
  iv1 = i-1
  iv2 = i
  iv3 = i+1
  if(iv1<1)iv1 = iv1 + nverts
  if(iv2>nverts)iv2 = iv2 - nverts
  if(iv3>nverts)iv3 = iv3 - nverts
  tritest = valid_triangle(dir,verts, nverts, iv1, iv2, iv3)
end do

call triangulate(dir,verts,nverts,vert_offset,faces)
do i = 1, nverts - 2
   write(0,*)"face",i,faces(3*i-2),faces(3*i-1),faces(3*i)
end do
write(0,*)"complete"

end program triangle_test
