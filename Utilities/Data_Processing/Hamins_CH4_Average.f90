!  Hamins_CH4.f90 
!
! Purpose
! =======
! This program computes the time averaged radiative profiles
! Utility program for the FDS validation case Hamins_CH4
!
! This program reads in the devc files for cases Hamins_CH4 1,5,23,21,7,19
! and computes time averaged radial and vertical profiles of radiative flux.
! The averaged proviles are written in files Hamins_CH4_#_devc_aver.csv
!
! Usage
! =====
! Compile for executable program.
! Run ./hamins_CH4 in the same directory witht the full devc files.
!

program NIST_Hamins
implicit none

! Variables
character(60) infile, outfile
integer ncd, ncol, nR, nV
logical EX
real mintime
real, dimension(:), allocatable :: M,R,Z,Ra1,Ra2,RaR,RaV

ncd = 100

allocate(M(ncd))
allocate(Ra1(ncd))
allocate(Ra2(ncd))
allocate(RaR(ncd))
allocate(RaV(ncd))
allocate(R(ncd))
allocate(Z(ncd))


!
! 10 cm burners
mintime = 3.0
ncol = 37
nR = 6
nV = 12
R = (/7.8, 10.0, 12.8, 15.0, 17.5, 20.6/)
Z = (/0.2, 1.1, 2.2, 3.3, 4.4, 6.6, 8.8, 12.2, 15.2, 20.3, 30.0, 39.3/)

! Hamins_CH4_1
infile = 'Hamins_CH4_01_devc.csv'
outfile = 'Hamins_CH4_01_devc_avg.csv'
call ReadCH4(infile,outfile,mintime,ncol,nR,nV,ncd,M,Ra1,Ra2,RaR,RaV,R,Z)

! Hamins_CH4_5
infile = 'Hamins_CH4_05_devc.csv'
outfile = 'Hamins_CH4_05_devc_avg.csv'
call ReadCH4(infile,outfile,mintime,ncol,nR,nV,ncd,M,Ra1,Ra2,RaR,RaV,R,Z)

!
! 38 cm burners
mintime = 2.5
ncol = 41
nR = 7
nV = 13
R = (/24.0, 31.2, 40.0, 51.2, 60.0, 70.0, 82.4/)
Z = (/0.008, 0.044, 0.088, 0.132, 0.176, 0.264, 0.352, 0.488, 0.608, 0.808, 1.200, 1.572, 1.890/)
Z = Z*100.0

! Hamins_CH4_21
infile = 'Hamins_CH4_21_devc.csv'
outfile = 'Hamins_CH4_21_devc_avg.csv'
call ReadCH4(infile,outfile,mintime,ncol,nR,nV,ncd,M,Ra1,Ra2,RaR,RaV,R,Z)

! Hamins_CH4_23
infile = 'Hamins_CH4_23_devc.csv'
outfile = 'Hamins_CH4_23_devc_avg.csv'
call ReadCH4(infile,outfile,mintime,ncol,nR,nV,ncd,M,Ra1,Ra2,RaR,RaV,R,Z)

!
! 100 cm burners
mintime = 3.0
ncol = 41
nR = 8
nV = 12
R = (/0.55, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.45/)
Z = (/0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5 /)
R = R*100.0
Z = Z*100.0

! Hamins_CH4_7
infile = 'Hamins_CH4_07_devc.csv'
outfile = 'Hamins_CH4_07_devc_avg.csv'
call ReadCH4(infile,outfile,mintime,ncol,nR,nV,ncd,M,Ra1,Ra2,RaR,RaV,R,Z)

! Hamins_CH4_19
infile = 'Hamins_CH4_19_devc.csv'
outfile = 'Hamins_CH4_19_devc_avg.csv'
call ReadCH4(infile,outfile,mintime,ncol,nR,nV,ncd,M,Ra1,Ra2,RaR,RaV,R,Z)

end program NIST_Hamins

subroutine ReadCH4(infile,outfile,mintime,ncol,nR,nV,ncd,M,Ra1,Ra2,RaR,RaV,R,Z)
character(60) infile, outfile
character(60) form
real mintime
integer FID, FOD,ncd, ncol, i, nt, nR, nV, c1, c2, c3, c4, c5
logical EX
real M(ncd),Ra1(ncd),Ra2(ncd),RaR(ncd),RaV(ncd),R(ncd),Z(ncd)
FID = 10
FOD = 11

c1 = 2
c2 = c1 + nR
c3 = c2 + nR
c4 = c3 + nV
c5 = c4 + nV

RaR = 0.0
RaV = 0.0

INQUIRE(file=infile,exist=EX)
If (EX) Then
   open(FID,FILE=infile,FORM='FORMATTED',action='read')
   write(*,*) 'Reading ',infile
   ! read blanks
   read(FID,*)
   read(FID,*)
   nt = 0
   do while (.not.eof(FID))
      read(FID,*) (M(i), i=1,ncol)
      If (M(1) >= mintime) Then
         nt = nt + 1
         Ra1 = 0.0
         Ra2 = 0.0
!        Angular average of radial fluxes
         Ra1(1:nR) = M(c1:c2-1)
         Ra2(1:nR) = M(c2:c3-1)
         RaR = RaR + 0.5*(Ra1+Ra2)
         Ra1 = 0.0
         Ra2 = 0.0
!        Angular average of vertical fluxes
         Ra1(1:nV) = M(c3:c4-1)
         Ra2(1:nV) = M(c4:c5-1)
         RaV = RaV + 0.5*(Ra1+Ra2)
      Endif
   enddo
   close(FID)
   ! time averages
   RaR = RaR / real(nt)
   RaV = RaV / real(nt)

Else
   write(*,*) infile, ' is missing!'
   stop
Endif

open(FID,file=trim(outfile),form='formatted',status='replace')
write(FID,"(a)") '"Radius","Radial Flux","Height","Vertical Flux"'

!write radial data
do i=1,nV
   if (i<=nR) then
      write(FID,"(3(e10.4,','),e10.4)") R(i),RaR(i),Z(i),RaV(i)
   else
      write(FID,"(A,1(e10.4,','),e10.4)") ',,',Z(i),RaV(i)
   endif
enddo
close(FID)

End subroutine ReadCH4
