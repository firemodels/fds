#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Annulus                 fire51 &
$RUNFDS Current_Results Main_Control_Room       fire52 &
$RUNFDS Current_Results Pump_Room               fire52 &
$RUNFDS Current_Results Switchgear_Room_Cabinet fire53 &

setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Training/NRC_EPRI
cd $WDIR/Current_Results
mpirun n11 n11           $FDS Cable_Spreading_Room.fds >& Cable_Spreading_Room.err &
mpirun n12 n12 n13 n13   $FDS Corridor.fds             >& Corridor.err             &
mpirun n14 n14           $FDS Switchgear_Room_MCC.fds  >& Switchgear_Room_MCC.err  &
mpirun n15 n15           $FDS Turbine_Building.fds     >& Turbine_Building.err     &
cd $WDIR

