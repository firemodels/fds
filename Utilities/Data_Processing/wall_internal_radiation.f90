program wall_internal_radiation

character(50) :: infile
real :: t,flux(5)
integer :: i

infile  = 'wall_internal_radiation_devc.csv'

open(10,file=infile,form='formatted',status='old')
read(10,*)
read(10,*)
read(10,*)
read(10,*)
read(10,*) t,(flux(i),i=1,5)
close(10)

open(11,file='wall_internal_radiation.tex',form='formatted',status='replace')

write(11,'(a)') "\begin{center}"
write(11,'(a)') "\begin{tabular}{|c|c|c|} \hline"
write(11,'(a)') "$\tau$      & $S(\tau)$   & FDS \\"
write(11,'(a)') "(m$^{-1}$)  & (kW/m$^2$)  & (kW/m$^2$) \\ \hline\hline"
write(11,'(a,f5.3,a)') "0.01        & 2.897       &",-flux(1)," \\"
write(11,'(a,f5.2,a)') "0.1         & 24.94       &",-flux(2)," \\"
write(11,'(a,f5.2,a)') "0.5         & 82.95       &",-flux(3)," \\"
write(11,'(a,f5.1,a)') "1.0         & 116.3       &",-flux(4)," \\"
write(11,'(a,f5.1,a)') "10.         & 149.0       &",-flux(5)," \\ \hline"
write(11,'(a)') "\end{tabular}"
write(11,'(a)') "\end{center}"

end program
