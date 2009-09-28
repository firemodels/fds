set intelbin=c:\bin

call %intelbin%\ifortvars ia32

ifort -o fds2ascii fds2ascii.f90
ifort -o flame_height flame_height.f90
ifort -o Hamins_CH4_Average Hamins_CH4_Average.f90
ifort -o layer_height layer_height.f
ifort -o NIST_RSE_1994 NIST_RSE_1994.f90
ifort -o ns2d ns2d.f90
ifort -o radiation_box radiation_box.f90
ifort -o radiation_plane_layer radiation_plane_layer.f90
ifort -o Ulster_SBI_fds2ascii Ulster_SBI_fds2ascii.f
ifort -o wall_internal_radiation wall_internal_radiation.f90
pause