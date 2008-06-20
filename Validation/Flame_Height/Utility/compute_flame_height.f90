program compute_flame_height

character(30) :: infile(16)
real :: z(0:100),hrrpul(100),height(16)

infile(1)  = 'Qs=p1_fds2ascii.csv'
infile(2)  = 'Qs=p2_fds2ascii.csv'
infile(3)  = 'Qs=p5_fds2ascii.csv'
infile(4)  = 'Qs=1_fds2ascii.csv'
infile(5)  = 'Qs=2_fds2ascii.csv'
infile(6)  = 'Qs=5_fds2ascii.csv'
infile(7)  = 'Qs=10_fds2ascii.csv'
infile(8)  = 'Qs=20_fds2ascii.csv'
infile(9)  = 'Qs=50_fds2ascii.csv'
infile(10) = 'Qs=100_fds2ascii.csv'
infile(11) = 'Qs=200_fds2ascii.csv'
infile(12) = 'Qs=500_fds2ascii.csv'
infile(13) = 'Qs=1000_fds2ascii.csv'
infile(14) = 'Qs=2000_fds2ascii.csv'
infile(15) = 'Qs=5000_fds2ascii.csv'
infile(16) = 'Qs=10000_fds2ascii.csv'

file_loop: do n=1,16

   open(10,file=infile(n),form='formatted',status='old')
   read(10,*)
   read(10,*)
   z(0) = 0.
   do k=1,76
      read(10,*,end=20) z(k),hrrpul(k)
   enddo
20 continue
   sum = 0.
   do k=1,76
      sum = sum + hrrpul(k)*(z(k)-z(k-1))
   enddo
   sum1 = 0.
   do k=1,76
      sum1 = sum1 + hrrpul(k)*(z(k)-z(k-1))
      if (sum1/sum>0.95) then
         height(n) = z(k)
         exit
      endif
   enddo

   write(6,"(i3,',',f8.2)") n,height(n)
   close(10)

enddo file_loop

end program
