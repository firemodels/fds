program ns2d

character(30) :: infile(4,2),outfile(4,2)
real :: t(10000),u(10000,4,2),dx(4),rms(4,2),cnt(4,2),nu(2),pi
integer :: i,j,k

infile(1,1) = 'ns2d_8_devc.csv'
infile(2,1) = 'ns2d_16_devc.csv'
infile(3,1) = 'ns2d_32_devc.csv'
infile(4,1) = 'ns2d_64_devc.csv'
infile(1,2) = 'ns2d_8_nupt1_devc.csv'
infile(2,2) = 'ns2d_16_nupt1_devc.csv'
infile(3,2) = 'ns2d_32_nupt1_devc.csv'
infile(4,2) = 'ns2d_64_nupt1_devc.csv'

outfile(1,1) = 'ns2d_8_exact.csv'
outfile(2,1) = 'ns2d_16_exact.csv'
outfile(3,1) = 'ns2d_32_exact.csv'
outfile(4,1) = 'ns2d_64_exact.csv'
outfile(1,2) = 'ns2d_8_nupt1_exact.csv'
outfile(2,2) = 'ns2d_16_nupt1_exact.csv'
outfile(3,2) = 'ns2d_32_nupt1_exact.csv'
outfile(4,2) = 'ns2d_64_nupt1_exact.csv'

nu(1) = 0.0
nu(2) = 0.1

pi = 4.*atan(1.)
x = pi
y = pi

dx(1) = 2.*pi/8.
dx(2) = 2.*pi/16.
dx(3) = 2.*pi/32.
dx(4) = 2.*pi/64.

case_loop: do k=1,2
   resolution_loop: do j=1,4
      open(10,file=infile(j,k),form='formatted',status='old')
      open(11,file=outfile(j,k),form='formatted',status='replace')
      write(11,'(a)') 'Time,u-vel'
      read(10,*)
      read(10,*)
      rms(j,k) = 0.
      cnt(j,k) = 0.
      time_loop: do i=1,10000
         read(10,*,end=20) t(i),u(i,j,k)
         u_exact = 1. - 2.*cos(x-t(i))*sin(y-0.5*dx(j)-t(i))*exp(-2.*nu(k)*t(i))
         write(11,"(f8.3,',',f9.5)") t(i),u_exact
         rms(j,k) = rms(j,k) + (u(i,j,k)-u_exact)**2
         cnt(j,k) = cnt(j,k) + 1.
      enddo time_loop
   20 continue
      close(10)
      close(11)
      rms(j,k) = sqrt(rms(j,k)/cnt(j,k))
   enddo resolution_loop
enddo case_loop

! Write the error files

open(11,file='ns2d_error.csv',form='formatted',status='replace')
open(12,file='ns2d_nupt1_error.csv',form='formatted',status='replace')

do k=1,2
   write(10+k,'(a)') 'dx,dx^2,rms error'
   do j=1,4
      write(10+k,"(2(f8.3,','),f8.4)") dx(j),dx(j)**2,rms(j,k)
   enddo
enddo

end program
