program NIST_RSE_1994

character(50) :: infile(10),label(13)
real :: t,value(13),hrr(10),value_avg(13)
integer :: i,j,k

hrr(1) =  50.
hrr(2) =  75.
hrr(3) = 100.
hrr(4) = 150.
hrr(5) = 200.
hrr(6) = 300.
hrr(7) = 400.
hrr(8) = 500.
hrr(9) = 600.

infile(1) = 'NIST_RSE_1994_50_devc.csv'
infile(2) = 'NIST_RSE_1994_75_devc.csv'
infile(3) = 'NIST_RSE_1994_100_devc.csv'
infile(4) = 'NIST_RSE_1994_150_devc.csv'
infile(5) = 'NIST_RSE_1994_200_devc.csv'
infile(6) = 'NIST_RSE_1994_300_devc.csv'
infile(7) = 'NIST_RSE_1994_400_devc.csv'
infile(8) = 'NIST_RSE_1994_500_devc.csv'
infile(9) = 'NIST_RSE_1994_600_devc.csv'

label(1)  = "O2Rear_FDS"
label(2)  = "CO2Rear_FDS"
label(3)  = "CORear_FDS"
label(4)  = "UHRear_FDS"
label(5)  = "O2Front_FDS"
label(6)  = "CO2Front_FDS"
label(7)  = "COFront_FDS"
label(8)  = "UHFront_FDS"
label(9)  = "TRSampA_FDS"
label(10) = "TRSampBB_FDS"
label(11) = "TFSampA_FDS"
label(12) = "TFSampBB_FDS"

open(11,file='NIST_RSE_1994_FDS.csv',form='formatted',status='replace')
write(11,"(12(a,','),a)") "HRR",(trim(label(i)),i=1,12)

file_loop: do n=1,9
   value_avg = 0.
   open(10,file=infile(n),form='formatted',status='old')
   read(10,*)
   read(10,*)
   do
      read(10,*,end=20) t,(value(k),k=1,12)
      value_avg(1:12) = 0.9*value_avg(1:12) + 0.1*value(1:12)
      cycle
   20 continue
      write(11,"(f5.1,',',8(f8.4,','),3(f6.1,','),f6.1)") hrr(n),(value_avg(k),k=1,12)
      exit
   enddo

   close(10)
enddo file_loop

end program
