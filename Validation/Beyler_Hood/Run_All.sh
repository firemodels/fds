#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_acetone_117_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_acetone_119_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_acetone_122_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_acetone_142_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_acetone_145_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_ethanol_106_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_ethanol_107_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_ethanol_108_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_ethanol_110_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_ethanol_115_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_isopropanol_130_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_isopropanol_132_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_isopropanol_133_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_isopropanol_136_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_isopropanol_141_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_methanol_942_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_methanol_943_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_methanol_945_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_methanol_947_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_methanol_951_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_232_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_257_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_287_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_303_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_307_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_318_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_322_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_334_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_355_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_359_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_371_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_389_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_389_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_429_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_433_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propane_445_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propylene_780_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propylene_805_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propylene_859_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propylene_870_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propylene_882_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propylene_886_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_propylene_910_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_toluene_160_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_toluene_162_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_toluene_165_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_toluene_166_lr.fds
$QFDS $DEBUG -d $INDIR $QUEUE Beyler_Hood_toluene_170_lr.fds

echo FDS cases submitted
