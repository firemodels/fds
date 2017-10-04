program triangle_test

! fire dynamics simulator, main program, multiple cpu version.

use precision_parameters
use complex_geometry, only: triangulate, valid_triangle

implicit none

integer :: i, j, nverts,input_unit
integer :: vert_offset, faces(300),loctype(300)
integer :: iv1, iv2, iv3, dir
logical :: tritest
real(fb) :: verts(300)

vert_offset=0

!nverts=6
!dir=3
!verts(1:3*nverts) = (/2.0, 0.0, 0.0, 0.0,0.0,0.0, 0.0,2.0,0.0, 1.0,2.0,0.0,  1.0,1.0,0.0, 2.0, 1.0, 0.0    /)

input_unit=7
open(input_unit,file="input.txt")
read(input_unit,*)nverts,dir
do i = 0, nverts-1
   read(input_unit,*)(verts(j),j=3*i+1,3*i+3)
end do

write(0,*)"valid_triangle"
do i = 1, nverts 
  iv1 = i-1
  iv2 = i
  iv3 = i+1
  if(iv1<1)iv1 = iv1 + nverts
  if(iv2>nverts)iv2 = iv2 - nverts
  if(iv3>nverts)iv3 = iv3 - nverts
!  tritest = valid_triangle(dir,verts, nverts, iv1, iv2, iv3)
!  write(0,*)iv1,iv2,iv3,tritest
end do

write(0,*)""
write(0,*)"triangulation"
call triangulate(dir,verts,nverts,vert_offset,faces,loctype)
do i = 1, nverts - 2
   write(0,*)"face:",i,"vertices:",faces(3*i-2),faces(3*i-1),faces(3*i)," loc:",loctype(i)
end do
write(0,*)"complete"

end program triangle_test
