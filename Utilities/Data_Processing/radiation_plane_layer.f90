program radiation_plane_layer

character(50) :: infile(6,5),kappa(6)
real :: t,flux(6,5),exact(6)
integer :: i,j

infile(1,1)  = 'radiation_plane_layer_1_1_devc.csv'
infile(2,1)  = 'radiation_plane_layer_2_1_devc.csv'
infile(3,1)  = 'radiation_plane_layer_3_1_devc.csv'
infile(4,1)  = 'radiation_plane_layer_4_1_devc.csv'
infile(5,1)  = 'radiation_plane_layer_5_1_devc.csv'
infile(6,1)  = 'radiation_plane_layer_6_1_devc.csv'
infile(1,2)  = 'radiation_plane_layer_1_2_devc.csv'
infile(2,2)  = 'radiation_plane_layer_2_2_devc.csv'
infile(3,2)  = 'radiation_plane_layer_3_2_devc.csv'
infile(4,2)  = 'radiation_plane_layer_4_2_devc.csv'
infile(5,2)  = 'radiation_plane_layer_5_2_devc.csv'
infile(6,2)  = 'radiation_plane_layer_6_2_devc.csv'
infile(1,3)  = 'radiation_plane_layer_1_3_devc.csv'
infile(2,3)  = 'radiation_plane_layer_2_3_devc.csv'
infile(3,3)  = 'radiation_plane_layer_3_3_devc.csv'
infile(4,3)  = 'radiation_plane_layer_4_3_devc.csv'
infile(5,3)  = 'radiation_plane_layer_5_3_devc.csv'
infile(6,3)  = 'radiation_plane_layer_6_3_devc.csv'
infile(1,4)  = 'radiation_plane_layer_1_4_devc.csv'
infile(2,4)  = 'radiation_plane_layer_2_4_devc.csv'
infile(3,4)  = 'radiation_plane_layer_3_4_devc.csv'
infile(4,4)  = 'radiation_plane_layer_4_4_devc.csv'
infile(5,4)  = 'radiation_plane_layer_5_4_devc.csv'
infile(6,4)  = 'radiation_plane_layer_6_4_devc.csv'
infile(1,5)  = 'radiation_plane_layer_1_5_devc.csv'
infile(2,5)  = 'radiation_plane_layer_2_5_devc.csv'
infile(3,5)  = 'radiation_plane_layer_3_5_devc.csv'
infile(4,5)  = 'radiation_plane_layer_4_5_devc.csv'
infile(5,5)  = 'radiation_plane_layer_5_5_devc.csv'
infile(6,5)  = 'radiation_plane_layer_6_5_devc.csv'

kappa(1) = "0.01"
kappa(2) = "0.1"
kappa(3) = "0.5"
kappa(4) = "1.0"
kappa(5) = "10"
kappa(6) = "$\infty$"

exact(1) = 2.8970
exact(2) = 24.9403
exact(3) = 82.9457
exact(4) = 116.2891
exact(5) = 148.9698
exact(6) = 148.9709

do j=1,5
   do i=1,6
      open(10,file=infile(i,j),form='formatted',status='old')
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*) t,flux(i,j)
      close(10)
   enddo
enddo 

open(11,file='../../Manuals/FDS_5_Verification_Guide/SCRIPT_FIGURES/radiation_plane_layer.tex',form='formatted',status='replace')

write(11,'(a)') "\begin{center}"
write(11,'(a)') "\begin{tabular}{|c|c|c|c|c|c|c|} \hline"
write(11,'(a)') "$\tau$ & $S(\tau)$ & \multicolumn{2}{|c|}{FDS (I=20,J=20)} &"
write(11,'(a)') "\multicolumn{2}{|c|}{FDS (I=20,J=1)} & FDS (I=150) \\ \cline{3-7}"
write(11,'(a)') "(m$^{-1}$) & (kW/m$^2$)& 1 band    & 6 bands   & 1 band    & 6 bands  & 1 band \\ \hline\hline"
do i=1,6
   write(11,"(a9,6(' & ',f9.4),' \\')") trim(kappa(i)),exact(i),(flux(i,j),j=1,5)
enddo
write(11,'(a)') "\hline"
write(11,'(a)') "\end{tabular}"
write(11,'(a)') "\end{center}"

end program
