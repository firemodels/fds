#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG  $QUEUE -d $INDIR Cup_acetone_ar.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_acetone_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_acetone_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_acetylene_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_acetylene_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_benzene_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_benzene_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_butane_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_butane_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_dodecane_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_ethanol_ar.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_ethanol_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_ethanol_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_ethylene_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_ethylene_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_heptane_ar.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_heptane_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_heptane_he.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_heptane_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_heptane_sf6.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_hexane_ar.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_hexane_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_hexane_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_hydrogen_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_hydrogen_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_isopropanol_ar.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_isopropanol_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_isopropanol_he.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_isopropanol_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_isopropanol_sf6.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_methane_ar.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_methane_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_methane_he.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_methane_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_methanol_ar.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_methanol_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_methanol_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_octane_ar.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_octane_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_octane_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_propane_ar.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_propane_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_propane_n2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_toluene_ar.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_toluene_co2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR Cup_toluene_n2.fds

echo FDS cases submitted
