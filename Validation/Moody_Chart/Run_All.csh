#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results moody_dpdx=-0p01_N16  fire61 &
$RUNFDS Current_Results moody_dpdx=-0p01_N32  fire61 &
$RUNFDS Current_Results moody_dpdx=-0p01_N8   fire62 &
$RUNFDS Current_Results moody_dpdx=-100_N16   fire62 &
$RUNFDS Current_Results moody_dpdx=-100_N32   fire63 &
$RUNFDS Current_Results moody_dpdx=-100_N8    fire63 &
$RUNFDS Current_Results moody_dpdx=-1_N16     fire64 &
$RUNFDS Current_Results moody_dpdx=-1_N32     fire64 &
$RUNFDS Current_Results moody_dpdx=-1_N8      fire65 &
$RUNFDS Current_Results poiseuille_N16_mu025  fire65 &
$RUNFDS Current_Results poiseuille_N32_mu025  fire66 &
$RUNFDS Current_Results poiseuille_N64_mu0125 fire66 &
$RUNFDS Current_Results poiseuille_N64_mu025  fire67 &
$RUNFDS Current_Results poiseuille_N8_mu025   fire67 &

echo FDS cases submitted
