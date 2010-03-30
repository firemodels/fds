program compute_flame_height

character(30) :: infile(16,3)
real :: z(0:200),hrrpul(200),height(16,3),qstar(16),diameter,sumold
integer :: i,n,npts

diameter = 1.13  ! Equivalent diameter of 1 m2 square

qstar(1)   = 0.1
qstar(2)   = 0.2
qstar(3)   = 0.5
qstar(4)   = 1.0
qstar(5)   = 2.0
qstar(6)   = 5.0
qstar(7)   = 10.
qstar(8)   = 20.
qstar(9)   = 50.
qstar(10)  = 100.
qstar(11)  = 200.
qstar(12)  = 500.
qstar(13)  = 1000.
qstar(14)  = 2000.
qstar(15)  = 5000.
qstar(16)  = 10000.

infile(1,1)  = 'Qs=p1_RI=05_fds2ascii.csv'
infile(2,1)  = 'Qs=p2_RI=05_fds2ascii.csv'
infile(3,1)  = 'Qs=p5_RI=05_fds2ascii.csv'
infile(4,1)  = 'Qs=1_RI=05_fds2ascii.csv'
infile(5,1)  = 'Qs=2_RI=05_fds2ascii.csv'
infile(6,1)  = 'Qs=5_RI=05_fds2ascii.csv'
infile(7,1)  = 'Qs=10_RI=05_fds2ascii.csv'
infile(8,1)  = 'Qs=20_RI=05_fds2ascii.csv'
infile(9,1)  = 'Qs=50_RI=05_fds2ascii.csv'
infile(10,1) = 'Qs=100_RI=05_fds2ascii.csv'
infile(11,1) = 'Qs=200_RI=05_fds2ascii.csv'
infile(12,1) = 'Qs=500_RI=05_fds2ascii.csv'
infile(13,1) = 'Qs=1000_RI=05_fds2ascii.csv'
infile(14,1) = 'Qs=2000_RI=05_fds2ascii.csv'
infile(15,1) = 'Qs=5000_RI=05_fds2ascii.csv'
infile(16,1) = 'Qs=10000_RI=05_fds2ascii.csv'

infile(1,2)  = 'Qs=p1_fds2ascii.csv'
infile(2,2)  = 'Qs=p2_fds2ascii.csv'
infile(3,2)  = 'Qs=p5_fds2ascii.csv'
infile(4,2)  = 'Qs=1_fds2ascii.csv'
infile(5,2)  = 'Qs=2_fds2ascii.csv'
infile(6,2)  = 'Qs=5_fds2ascii.csv'
infile(7,2)  = 'Qs=10_fds2ascii.csv'
infile(8,2)  = 'Qs=20_fds2ascii.csv'
infile(9,2)  = 'Qs=50_fds2ascii.csv'
infile(10,2) = 'Qs=100_fds2ascii.csv'
infile(11,2) = 'Qs=200_fds2ascii.csv'
infile(12,2) = 'Qs=500_fds2ascii.csv'
infile(13,2) = 'Qs=1000_fds2ascii.csv'
infile(14,2) = 'Qs=2000_fds2ascii.csv'
infile(15,2) = 'Qs=5000_fds2ascii.csv'
infile(16,2) = 'Qs=10000_fds2ascii.csv'

infile(1,3)  = 'Qs=p1_RI=20_fds2ascii.csv'
infile(2,3)  = 'Qs=p2_RI=20_fds2ascii.csv'
infile(3,3)  = 'Qs=p5_RI=20_fds2ascii.csv'
infile(4,3)  = 'Qs=1_RI=20_fds2ascii.csv'
infile(5,3)  = 'Qs=2_RI=20_fds2ascii.csv'
infile(6,3)  = 'Qs=5_RI=20_fds2ascii.csv'
infile(7,3)  = 'Qs=10_RI=20_fds2ascii.csv'
infile(8,3)  = 'Qs=20_RI=20_fds2ascii.csv'
infile(9,3)  = 'Qs=50_RI=20_fds2ascii.csv'
infile(10,3) = 'Qs=100_RI=20_fds2ascii.csv'
infile(11,3) = 'Qs=200_RI=20_fds2ascii.csv'
infile(12,3) = 'Qs=500_RI=20_fds2ascii.csv'
infile(13,3) = 'Qs=1000_RI=20_fds2ascii.csv'
infile(14,3) = 'Qs=2000_RI=20_fds2ascii.csv'
infile(15,3) = 'Qs=5000_RI=20_fds2ascii.csv'
infile(16,3) = 'Qs=10000_RI=20_fds2ascii.csv'

write(6,"(a)") "Q*,L/D (RI=5),L/D (RI=10),L/D (RI=20)"

file_loop: do n=1,16

   resolution_loop: do i=1,3

      if (i==1) npts=39
      if (i==2) npts=76
      if (i==3) npts=151

      open(10,file=infile(n,i),form='formatted',status='old')
      read(10,*)
      read(10,*)
      z    = 0.
      z(0) = 0.
      do k=1,npts
         read(10,*,end=20) z(k),hrrpul(k)
      enddo
   20 continue
      sum = 0.
      do k=1,npts
    !!   sum = sum + hrrpul(k)*(z(k)-z(k-1))
         sum = sum + hrrpul(k)
      enddo
      sum1 = 0.
      do k=1,npts
    !!   sum1 = sum1 + hrrpul(k)*(z(k)-z(k-1))
         sumold = sum1
         sum1 = sum1 + hrrpul(k)
         if (sum1/sum>0.99) then
            height(n,i) = z(k-1) + (z(k)-z(k-1))*(0.99*sum-sumold)/(sum1-sumold)
            exit
         endif
      enddo
  
   enddo resolution_loop

   write(6,"(f8.1,',',f8.2,',',f8.2,',',f8.2)") qstar(n),height(n,1)/diameter,height(n,2)/diameter,height(n,3)/diameter
   close(10)

enddo file_loop

end program
