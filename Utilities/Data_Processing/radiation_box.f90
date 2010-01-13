program radiation_box

character(50) :: infile(10),label(10)
real :: t,flux(20,10),x(20)
integer :: i,j,k

x( 1) = 2.50E-02
x( 2) = 7.50E-02
x( 3) = 1.25E-01
x( 4) = 1.75E-01
x( 5) = 2.25E-01
x( 6) = 2.75E-01
x( 7) = 3.25E-01
x( 8) = 3.75E-01
x( 9) = 4.25E-01
x(10) = 4.75E-01
x(11) = 5.25E-01
x(12) = 5.75E-01
x(13) = 6.25E-01
x(14) = 6.75E-01
x(15) = 7.25E-01
x(16) = 7.75E-01
x(17) = 8.25E-01
x(18) = 8.75E-01
x(19) = 9.25E-01
x(20) = 9.75E-01

infile(1)  = 'radiation_box__20___50_devc.csv'
infile(2)  = 'radiation_box__20__100_devc.csv'
infile(3)  = 'radiation_box__20__300_devc.csv'
infile(4)  = 'radiation_box__20_1000_devc.csv'
infile(5)  = 'radiation_box__20_2000_devc.csv'
infile(6)  = 'radiation_box_100___50_devc.csv'
infile(7)  = 'radiation_box_100__100_devc.csv'
infile(8)  = 'radiation_box_100__300_devc.csv'
infile(9)  = 'radiation_box_100_1000_devc.csv'
infile(10) = 'radiation_box_100_2000_devc.csv'

label( 1) = 'Flux_20_50'
label( 2) = 'Flux_20_100'
label( 3) = 'Flux_20_300'
label( 4) = 'Flux_20_1000'
label( 5) = 'Flux_20_2000'
label( 6) = 'Flux_100_50'
label( 7) = 'Flux_100_100'
label( 8) = 'Flux_100_300'
label( 9) = 'Flux_100_1000'
label(10) = 'Flux_100_2000'

file_loop: do n=1,10
   open(10,file=infile(n),form='formatted',status='old')
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*,end=20) t,(flux(k,n),k=1,20)
20 continue
   close(10)
enddo file_loop

open(11,file='radiation_box_devc.csv',form='formatted',status='replace')
write(11,"(10(a,','),a)") "Position",(trim(label(n)),n=1,10)
do i=1,20
   write(11,"(10(f8.3,','),f8.3)") x(i),(flux(i,j),j=1,10)
enddo

end program
