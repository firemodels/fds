#!/bin/csh -f

#  This script compiles a few utility programs.
 
set platform=intel64
setenv SVNROOT ~/FDS-SMV
source $SVNROOT/FDS_Compilation/Scripts/set_fort.csh $platform
if ($status == 1) exit

ifort -o Beyler_Hood Beyler_Hood.f90
ifort -o fds2ascii fds2ascii.f90
ifort -o flame_height flame_height.f90
ifort -o Hamins_CH4_Average Hamins_CH4_Average.f90
ifort -o layer_height layer_height.f
ifort -o NIST_RSE_1994 NIST_RSE_1994.f90
ifort -o ns2d ns2d.f90
ifort -o radiation_box radiation_box.f90
ifort -o radiation_plane_layer radiation_plane_layer.f90
ifort -o wall_internal_radiation wall_internal_radiation.f90

