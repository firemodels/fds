#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results moody_dpdx=-0p01_N16  fire41 &
$RUNFDS Current_Results moody_dpdx=-0p01_N32  fire41 &
$RUNFDS Current_Results moody_dpdx=-0p01_N8   fire41 &
$RUNFDS Current_Results moody_dpdx=-100_N16   fire43 &
$RUNFDS Current_Results moody_dpdx=-100_N32   fire43 &
$RUNFDS Current_Results moody_dpdx=-100_N8    fire43 &
$RUNFDS Current_Results moody_dpdx=-1_N16     fire44 &
$RUNFDS Current_Results moody_dpdx=-1_N32     fire44 &
$RUNFDS Current_Results moody_dpdx=-1_N8      fire45 &
$RUNFDS Current_Results poiseuille_N16_mu025  fire45 &
$RUNFDS Current_Results poiseuille_N32_mu025  fire46 &
$RUNFDS Current_Results poiseuille_N64_mu0125 fire46 &
$RUNFDS Current_Results poiseuille_N64_mu025  fire47 &
$RUNFDS Current_Results poiseuille_N8_mu025   fire47 &

$RUNFDS Current_Results z0=p0001_dpdx=-1_N8     fire41 &
$RUNFDS Current_Results z0=p001_dpdx=-1_N8      fire41 &
$RUNFDS Current_Results z0=p01_dpdx=-1_N8       fire41 &
$RUNFDS Current_Results z0=p1_dpdx=-1_N8        fire43 &
$RUNFDS Current_Results z0=p0001_dpdx=-1_N16    fire43 &
$RUNFDS Current_Results z0=p001_dpdx=-1_N16     fire43 &
$RUNFDS Current_Results z0=p01_dpdx=-1_N16      fire44 &
$RUNFDS Current_Results z0=p0001_dpdx=-100_N8   fire44 &
$RUNFDS Current_Results z0=p001_dpdx=-100_N8    fire44 &
$RUNFDS Current_Results z0=p01_dpdx=-100_N8     fire45 &
$RUNFDS Current_Results z0=p1_dpdx=-100_N8      fire45 &
$RUNFDS Current_Results z0=p0001_dpdx=-100_N16  fire45 &
$RUNFDS Current_Results z0=p001_dpdx=-100_N16   fire46 &
$RUNFDS Current_Results z0=p01_dpdx=-100_N16    fire46 &
$RUNFDS Current_Results z0=p0001_dpdx=-p0001_N8 fire46 &
$RUNFDS Current_Results z0=p0001_dpdx=-p01_N8   fire47 &
$RUNFDS Current_Results z0=p001_dpdx=-p0001_N8  fire47 &
$RUNFDS Current_Results z0=p001_dpdx=-p01_N8    fire47 &
$RUNFDS Current_Results z0=p01_dpdx=-p0001_N8   fire47 &
$RUNFDS Current_Results z0=p01_dpdx=-p01_N8     fire46 &
$RUNFDS Current_Results z0=p1_dpdx=-p0001_N8    fire45 &
$RUNFDS Current_Results z0=p1_dpdx=-p01_N8      fire44 &
$RUNFDS Current_Results z0=p02_dpdx=-1_N16      fire43 &
$RUNFDS Current_Results z0=p02_dpdx=-100_N16    fire41 &

echo FDS cases submitted
