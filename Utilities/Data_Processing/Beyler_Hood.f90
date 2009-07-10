!This program processes the fds devc_csv files for the Beyler_Hood validation cases
!$Id: init.f90 4308 2009-07-01 18:44:04Z mcgratta $
PROGRAM BEYLER_HOOD
IMPLICIT NONE
INTEGER :: i,k
CHARACTER(40) :: infile(12)
REAL(8) :: time(1000),fuel(4,1000),o2(4,1000),co2(4,1000),co(4,1000),n2(4,1000),h2o(4,1000),c(4,1000),&
           t(4,1000),wv(4,1000),nv(4,1000)
REAL(8) fave(12),o2ave(12),co2ave(12),coave(12)

infile(1)  = 'Beyler_Hood_propane-792-19-5_devc.csv'
infile(2)  = 'Beyler_Hood_propane-1353-19-5_devc.csv'
infile(3)  = 'Beyler_Hood_propane-1825-19-5_devc.csv'
infile(4)  = 'Beyler_Hood_propane-2430-19-5_devc.csv'
infile(5)  = 'Beyler_Hood_propane-821-19--10_devc.csv'
infile(6)  = 'Beyler_Hood_propane-1353-19--10_devc.csv'
infile(7)  = 'Beyler_Hood_propane-1825-19--10_devc.csv'
infile(8)  = 'Beyler_Hood_propane-2430-19--10_devc.csv'
infile(9)  = 'Beyler_Hood_propane-821-19-0_devc.csv'
infile(10)  = 'Beyler_Hood_propane-1353-19-0_devc.csv'
infile(11)  = 'Beyler_Hood_propane-1825-19-0_devc.csv'
infile(12)  = 'Beyler_Hood_propane-3152-19-0_devc.csv'

fave=0._8
o2ave=0._8
co2ave=0._8
coave=0._8

DO i=1,12
   OPEN(11,FILE=infile(i))
   READ(11,*)
   READ(11,*)
   READ(11,*)
   DO k=1,1000
      READ(11,*) time(k),fuel(1:4,k),o2(1:4,k),co2(1:4,k),co(1:4,k),n2(1:4,k),h2o(1:4,k),c(1:4,k),&
                 t(1:4,k),wv(1:4,k),nv(1:4,k)
   ENDDO
   DO k=901,1000
      fave(i) =fave(i)   +0.25_8*SUM(fuel(1:4,k))
      o2ave(i) =o2ave(i) +0.25_8*SUM(o2(1:4,k))
      co2ave(i)=co2ave(i)+0.25_8*SUM(co2(1:4,k))
      coave(i) =coave(i) +0.25_8*SUM(co(1:4,k))
   ENDDO
   CLOSE(11)
ENDDO
fave = fave / 100._8
o2ave = o2ave / 100._8
co2ave = co2ave / 100._8
coave = coave / 100._8

OPEN(11,FILE='Beyler_Hood_FDS.csv')
WRITE(11,*) "Q-5,CO2-5,CO-5,O2-5,C3H8-5,Q--10,CO2--10,CO--10,O2--10,C3H8--10,Q-0,CO2-0,CO-0,O2-0,C3H8-0"
WRITE(11,*) 7.92,",",co2ave(1),",",coave(1),",",o2ave(1),",",fave(1),",", &
            8.21,",",co2ave(5),",",coave(5),",",o2ave(5),",",fave(5),",", &
            8.21,",",co2ave(9),",",coave(9),",",o2ave(9),",",fave(9)
WRITE(11,*) 13.53,",",co2ave(2),",",coave(2),",",o2ave(2),",",fave(2),",", &
            13.53,",",co2ave(6),",",coave(6),",",o2ave(6),",",fave(6),",", &
            13.53,",",co2ave(10),",",coave(10),",",o2ave(10),",",fave(10)
WRITE(11,*) 18.25,",",co2ave(3),",",coave(3),",",o2ave(3),",",fave(3),",", &
            18.25,",",co2ave(7),",",coave(7),",",o2ave(7),",",fave(7),",", &
            18.25,",",co2ave(11),",",coave(11),",",o2ave(11),",",fave(11)
WRITE(11,*) 24.3,",",co2ave(4),",",coave(4),",",o2ave(4),",",fave(4),",", &
            24.3,",",co2ave(8),",",coave(8),",",o2ave(8),",",fave(8),",", &
            31.52,",",co2ave(12),",",coave(12),",",o2ave(12),",",fave(12)                        
CLOSE (11)

END PROGRAM
