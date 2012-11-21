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
$RUNFDS $INDIR Beyler_Hood_acetone_120_lr
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
$RUNFDS $INDIR Beyler_Hood_isopropanol_140_lr
$RUNFDS $INDIR Beyler_Hood_isopropanol_141_lr
$RUNFDS $INDIR Beyler_Hood_methanol_942_lr
$RUNFDS $INDIR Beyler_Hood_methanol_943_lr
$RUNFDS $INDIR Beyler_Hood_methanol_945_lr
$RUNFDS $INDIR Beyler_Hood_methanol_947_lr
$RUNFDS $INDIR Beyler_Hood_methanol_951_lr
$RUNFDS $INDIR Beyler_Hood_propane_232
$RUNFDS $INDIR Beyler_Hood_propane_257
$RUNFDS $INDIR Beyler_Hood_propane_287
$RUNFDS $INDIR Beyler_Hood_propane_303
$RUNFDS $INDIR Beyler_Hood_propane_307
$RUNFDS $INDIR Beyler_Hood_propane_318_lr
$RUNFDS $INDIR Beyler_Hood_propane_318_mr
$RUNFDS $INDIR Beyler_Hood_propane_318
$RUNFDS $INDIR Beyler_Hood_propane_322
$RUNFDS $INDIR Beyler_Hood_propane_334
$RUNFDS $INDIR Beyler_Hood_propane_355
$RUNFDS $INDIR Beyler_Hood_propane_359
$RUNFDS $INDIR Beyler_Hood_propane_371
$RUNFDS $INDIR Beyler_Hood_propane_389_lr
$RUNFDS $INDIR Beyler_Hood_propane_389_mr
$RUNFDS $INDIR Beyler_Hood_propane_389
$RUNFDS $INDIR Beyler_Hood_propane_429
$RUNFDS $INDIR Beyler_Hood_propane_433
$RUNFDS $INDIR Beyler_Hood_propane_445
$RUNFDS $INDIR Beyler_Hood_propylene_780_lr
$RUNFDS $INDIR Beyler_Hood_propylene_805_lr
$RUNFDS $INDIR Beyler_Hood_propylene_859_lr
$RUNFDS $INDIR Beyler_Hood_propylene_870_lr
$RUNFDS $INDIR Beyler_Hood_propylene_882_lr
$RUNFDS $INDIR Beyler_Hood_propylene_910_lr
$RUNFDS $INDIR Beyler_Hood_toluene_160_lr
$RUNFDS $INDIR Beyler_Hood_toluene_162_lr
$RUNFDS $INDIR Beyler_Hood_toluene_165_lr
$RUNFDS $INDIR Beyler_Hood_toluene_166_lr
$RUNFDS $INDIR Beyler_Hood_toluene_170_lr


