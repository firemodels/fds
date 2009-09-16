#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results moody_dpdx=-0p01_N16.fds  fire51 &
$RUNFDS Current_Results moody_dpdx=-0p01_N32.fds  fire51 &
$RUNFDS Current_Results moody_dpdx=-0p01_N8.fds   fire52 &
$RUNFDS Current_Results moody_dpdx=-100_N16.fds   fire52 &
$RUNFDS Current_Results moody_dpdx=-100_N32.fds   fire53 &
$RUNFDS Current_Results moody_dpdx=-100_N8.fds    fire53 &
$RUNFDS Current_Results moody_dpdx=-1_N16.fds     fire54 &
$RUNFDS Current_Results moody_dpdx=-1_N32.fds     fire54 &
$RUNFDS Current_Results moody_dpdx=-1_N8.fds      fire55 &
$RUNFDS Current_Results poiseuille_N16_mu025.fds  fire55 &
$RUNFDS Current_Results poiseuille_N32_mu025.fds  fire56 &
$RUNFDS Current_Results poiseuille_N64_mu0125.fds fire56 &
$RUNFDS Current_Results poiseuille_N64_mu025.fds  fire57 &
$RUNFDS Current_Results poiseuille_N8_mu025.fds   fire57 &

echo FDS cases submitted
