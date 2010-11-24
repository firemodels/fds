#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS ~/FDS/FDS5/bin/fds5_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Annulus                    fire45 &
$RUNFDS Current_Results Main_Control_Room_Purge    fire45 &
$RUNFDS Current_Results Main_Control_Room_No_Purge fire46 &
$RUNFDS Current_Results Pump_Room                  fire46 &

setenv FDS ~/FDS/FDS5/bin/fds5_mpi_linux_32
set WDIR=$SVNROOT/Training/NRC_EPRI
cd $WDIR/Current_Results
mpirun n17 n17              $FDS Cable_Spreading_Room.fds     >& Cable_Spreading_Room.err     &
mpirun n11 n11 n11 n11      $FDS Corridor.fds                 >& Corridor.err                 &
mpirun n18 n19 n19 n20 n20  $FDS Switchgear_Room_Cabinet.fds  >& Switchgear_Room_Cabinet.err  &
mpirun n30 n30              $FDS Switchgear_Room_MCC.fds      >& Switchgear_Room_MCC.err      &
mpirun n26 n26              $FDS Turbine_Building.fds         >& Turbine_Building.err         &
cd $WDIR

