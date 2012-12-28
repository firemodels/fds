#!/bin/bash -f

export SVNROOT=`pwd`/../..
export BASEDIR=`pwd`
export INDIR=Current_Results
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Beyler_Hood_acetone_117_lr
$RUNFDS $INDIR Beyler_Hood_acetone_119_lr
$RUNFDS $INDIR Beyler_Hood_acetone_122_lr
$RUNFDS $INDIR Beyler_Hood_acetone_142_lr
$RUNFDS $INDIR Beyler_Hood_acetone_145_lr
$RUNFDS $INDIR Beyler_Hood_ethanol_106_lr
$RUNFDS $INDIR Beyler_Hood_ethanol_107_lr
$RUNFDS $INDIR Beyler_Hood_ethanol_108_lr
$RUNFDS $INDIR Beyler_Hood_ethanol_110_lr
$RUNFDS $INDIR Beyler_Hood_ethanol_115_lr
$RUNFDS $INDIR Beyler_Hood_isopropanol_130_lr
$RUNFDS $INDIR Beyler_Hood_isopropanol_132_lr
$RUNFDS $INDIR Beyler_Hood_isopropanol_133_lr
$RUNFDS $INDIR Beyler_Hood_isopropanol_136_lr
$RUNFDS $INDIR Beyler_Hood_isopropanol_141_lr
$RUNFDS $INDIR Beyler_Hood_methanol_942_lr
$RUNFDS $INDIR Beyler_Hood_methanol_943_lr
$RUNFDS $INDIR Beyler_Hood_methanol_945_lr
$RUNFDS $INDIR Beyler_Hood_methanol_947_lr
$RUNFDS $INDIR Beyler_Hood_methanol_951_lr
$RUNFDS $INDIR Beyler_Hood_propane_232_lr
$RUNFDS $INDIR Beyler_Hood_propane_257_lr
$RUNFDS $INDIR Beyler_Hood_propane_287_lr
$RUNFDS $INDIR Beyler_Hood_propane_303_lr
$RUNFDS $INDIR Beyler_Hood_propane_307_lr
$RUNFDS $INDIR Beyler_Hood_propane_318_lr
$RUNFDS $INDIR Beyler_Hood_propane_318_mr
$RUNFDS $INDIR Beyler_Hood_propane_318_hr
$RUNFDS $INDIR Beyler_Hood_propane_322_lr
$RUNFDS $INDIR Beyler_Hood_propane_334_lr
$RUNFDS $INDIR Beyler_Hood_propane_355_lr
$RUNFDS $INDIR Beyler_Hood_propane_359_lr
$RUNFDS $INDIR Beyler_Hood_propane_371_lr
$RUNFDS $INDIR Beyler_Hood_propane_389_lr
$RUNFDS $INDIR Beyler_Hood_propane_389_lr
$RUNFDS $INDIR Beyler_Hood_propane_389_mr
$RUNFDS $INDIR Beyler_Hood_propane_389_hr
$RUNFDS $INDIR Beyler_Hood_propane_429_lr
$RUNFDS $INDIR Beyler_Hood_propane_433_lr
$RUNFDS $INDIR Beyler_Hood_propane_445_lr
$RUNFDS $INDIR Beyler_Hood_propylene_780_lr
$RUNFDS $INDIR Beyler_Hood_propylene_805_lr
$RUNFDS $INDIR Beyler_Hood_propylene_859_lr
$RUNFDS $INDIR Beyler_Hood_propylene_870_lr
$RUNFDS $INDIR Beyler_Hood_propylene_882_lr
$RUNFDS $INDIR Beyler_Hood_propylene_886_lr
$RUNFDS $INDIR Beyler_Hood_propylene_910_lr
$RUNFDS $INDIR Beyler_Hood_toluene_160_lr
$RUNFDS $INDIR Beyler_Hood_toluene_162_lr
$RUNFDS $INDIR Beyler_Hood_toluene_165_lr
$RUNFDS $INDIR Beyler_Hood_toluene_166_lr
$RUNFDS $INDIR Beyler_Hood_toluene_170_lr


