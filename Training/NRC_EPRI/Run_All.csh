#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results_2 Annulus                 fire70 &
$RUNFDS Current_Results_2 Main_Control_Room_Purge fire70 &
$RUNFDS Current_Results_2 Main_Control_Room_No_Purge fire70 &
$RUNFDS Current_Results_2 Pump_Room               fire70 &

setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Training/NRC_EPRI
cd $WDIR/Current_Results_2
mpirun n1 n1                $FDS Cable_Spreading_Room.fds     >& Cable_Spreading_Room.err     &
mpirun n2 n2 n3 n3          $FDS Corridor.fds                 >& Corridor.err                 &
mpirun n1 n1 n1 n0 n0       $FDS Switchgear_Room_Cabinet.fds  >& Switchgear_Room_Cabinet.err  &
mpirun n2 n2                $FDS Switchgear_Room_MCC.fds      >& Switchgear_Room_MCC.err      &
mpirun n3 n3                $FDS Turbine_Building.fds         >& Turbine_Building.err         &
cd $WDIR

