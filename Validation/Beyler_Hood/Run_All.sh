#!/bin/bash -f

export SVNROOT=`pwd`/../..
export BASEDIR=`pwd`
export INDIR=Current_Results
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Beyler_Hood_acetone_117
$RUNFDS $INDIR Beyler_Hood_acetone_118
$RUNFDS $INDIR Beyler_Hood_acetone_119
$RUNFDS $INDIR Beyler_Hood_acetone_120
$RUNFDS $INDIR Beyler_Hood_acetone_121
$RUNFDS $INDIR Beyler_Hood_acetone_122
$RUNFDS $INDIR Beyler_Hood_acetone_124
$RUNFDS $INDIR Beyler_Hood_acetone_126
$RUNFDS $INDIR Beyler_Hood_acetone_128
$RUNFDS $INDIR Beyler_Hood_acetone_129
$RUNFDS $INDIR Beyler_Hood_acetone_142
$RUNFDS $INDIR Beyler_Hood_acetone_143
$RUNFDS $INDIR Beyler_Hood_acetone_145
$RUNFDS $INDIR Beyler_Hood_ethanol_106
$RUNFDS $INDIR Beyler_Hood_ethanol_107
$RUNFDS $INDIR Beyler_Hood_ethanol_108
$RUNFDS $INDIR Beyler_Hood_ethanol_109
$RUNFDS $INDIR Beyler_Hood_ethanol_110
$RUNFDS $INDIR Beyler_Hood_ethanol_111
$RUNFDS $INDIR Beyler_Hood_ethanol_112
$RUNFDS $INDIR Beyler_Hood_ethanol_113
$RUNFDS $INDIR Beyler_Hood_ethanol_114
$RUNFDS $INDIR Beyler_Hood_ethanol_115
$RUNFDS $INDIR Beyler_Hood_ethanol_116
$RUNFDS $INDIR Beyler_Hood_isopropanol_129
$RUNFDS $INDIR Beyler_Hood_isopropanol_130
$RUNFDS $INDIR Beyler_Hood_isopropanol_131
$RUNFDS $INDIR Beyler_Hood_isopropanol_132
$RUNFDS $INDIR Beyler_Hood_isopropanol_133
$RUNFDS $INDIR Beyler_Hood_isopropanol_135
$RUNFDS $INDIR Beyler_Hood_isopropanol_136
$RUNFDS $INDIR Beyler_Hood_isopropanol_137
$RUNFDS $INDIR Beyler_Hood_isopropanol_138
$RUNFDS $INDIR Beyler_Hood_isopropanol_140
$RUNFDS $INDIR Beyler_Hood_isopropanol_141
$RUNFDS $INDIR Beyler_Hood_methanol_938
$RUNFDS $INDIR Beyler_Hood_methanol_939
$RUNFDS $INDIR Beyler_Hood_methanol_940
$RUNFDS $INDIR Beyler_Hood_methanol_941
$RUNFDS $INDIR Beyler_Hood_methanol_942
$RUNFDS $INDIR Beyler_Hood_methanol_943
$RUNFDS $INDIR Beyler_Hood_methanol_945
$RUNFDS $INDIR Beyler_Hood_methanol_946
$RUNFDS $INDIR Beyler_Hood_methanol_947
$RUNFDS $INDIR Beyler_Hood_methanol_948
$RUNFDS $INDIR Beyler_Hood_methanol_949
$RUNFDS $INDIR Beyler_Hood_methanol_950
$RUNFDS $INDIR Beyler_Hood_methanol_951
$RUNFDS $INDIR Beyler_Hood_methanol_953
$RUNFDS $INDIR Beyler_Hood_methanol_954
$RUNFDS $INDIR Beyler_Hood_methanol_955
$RUNFDS $INDIR Beyler_Hood_methanol_956
$RUNFDS $INDIR Beyler_Hood_methanol_957
$RUNFDS $INDIR Beyler_Hood_propane_220
$RUNFDS $INDIR Beyler_Hood_propane_228
$RUNFDS $INDIR Beyler_Hood_propane_232
$RUNFDS $INDIR Beyler_Hood_propane_239
$RUNFDS $INDIR Beyler_Hood_propane_241
$RUNFDS $INDIR Beyler_Hood_propane_245
$RUNFDS $INDIR Beyler_Hood_propane_253
$RUNFDS $INDIR Beyler_Hood_propane_257
$RUNFDS $INDIR Beyler_Hood_propane_269
$RUNFDS $INDIR Beyler_Hood_propane_278
$RUNFDS $INDIR Beyler_Hood_propane_282
$RUNFDS $INDIR Beyler_Hood_propane_286
$RUNFDS $INDIR Beyler_Hood_propane_287
$RUNFDS $INDIR Beyler_Hood_propane_291
$RUNFDS $INDIR Beyler_Hood_propane_295
$RUNFDS $INDIR Beyler_Hood_propane_299
$RUNFDS $INDIR Beyler_Hood_propane_303
$RUNFDS $INDIR Beyler_Hood_propane_307
$RUNFDS $INDIR Beyler_Hood_propane_311
$RUNFDS $INDIR Beyler_Hood_propane_314
$RUNFDS $INDIR Beyler_Hood_propane_318
$RUNFDS $INDIR Beyler_Hood_propane_322
$RUNFDS $INDIR Beyler_Hood_propane_326
$RUNFDS $INDIR Beyler_Hood_propane_330
$RUNFDS $INDIR Beyler_Hood_propane_334
$RUNFDS $INDIR Beyler_Hood_propane_338
$RUNFDS $INDIR Beyler_Hood_propane_351
$RUNFDS $INDIR Beyler_Hood_propane_355
$RUNFDS $INDIR Beyler_Hood_propane_359
$RUNFDS $INDIR Beyler_Hood_propane_363
$RUNFDS $INDIR Beyler_Hood_propane_367
$RUNFDS $INDIR Beyler_Hood_propane_371
$RUNFDS $INDIR Beyler_Hood_propane_385
$RUNFDS $INDIR Beyler_Hood_propane_389
$RUNFDS $INDIR Beyler_Hood_propane_393
$RUNFDS $INDIR Beyler_Hood_propane_399
$RUNFDS $INDIR Beyler_Hood_propane_403
$RUNFDS $INDIR Beyler_Hood_propane_407
$RUNFDS $INDIR Beyler_Hood_propane_412
$RUNFDS $INDIR Beyler_Hood_propane_417
$RUNFDS $INDIR Beyler_Hood_propane_421
$RUNFDS $INDIR Beyler_Hood_propane_425
$RUNFDS $INDIR Beyler_Hood_propane_429
$RUNFDS $INDIR Beyler_Hood_propane_433
$RUNFDS $INDIR Beyler_Hood_propane_437
$RUNFDS $INDIR Beyler_Hood_propane_441
$RUNFDS $INDIR Beyler_Hood_propane_445
$RUNFDS $INDIR Beyler_Hood_propylene_776
$RUNFDS $INDIR Beyler_Hood_propylene_780
$RUNFDS $INDIR Beyler_Hood_propylene_784
$RUNFDS $INDIR Beyler_Hood_propylene_792
$RUNFDS $INDIR Beyler_Hood_propylene_801
$RUNFDS $INDIR Beyler_Hood_propylene_805
$RUNFDS $INDIR Beyler_Hood_propylene_809
$RUNFDS $INDIR Beyler_Hood_propylene_813
$RUNFDS $INDIR Beyler_Hood_propylene_838
$RUNFDS $INDIR Beyler_Hood_propylene_842
$RUNFDS $INDIR Beyler_Hood_propylene_846
$RUNFDS $INDIR Beyler_Hood_propylene_850
$RUNFDS $INDIR Beyler_Hood_propylene_854
$RUNFDS $INDIR Beyler_Hood_propylene_859
$RUNFDS $INDIR Beyler_Hood_propylene_863
$RUNFDS $INDIR Beyler_Hood_propylene_867
$RUNFDS $INDIR Beyler_Hood_propylene_870
$RUNFDS $INDIR Beyler_Hood_propylene_874
$RUNFDS $INDIR Beyler_Hood_propylene_878
$RUNFDS $INDIR Beyler_Hood_propylene_882
$RUNFDS $INDIR Beyler_Hood_propylene_886
$RUNFDS $INDIR Beyler_Hood_propylene_890
$RUNFDS $INDIR Beyler_Hood_propylene_899
$RUNFDS $INDIR Beyler_Hood_propylene_903
$RUNFDS $INDIR Beyler_Hood_propylene_910
$RUNFDS $INDIR Beyler_Hood_propylene_914
$RUNFDS $INDIR Beyler_Hood_toluene_160
$RUNFDS $INDIR Beyler_Hood_toluene_161
$RUNFDS $INDIR Beyler_Hood_toluene_162
$RUNFDS $INDIR Beyler_Hood_toluene_163
$RUNFDS $INDIR Beyler_Hood_toluene_164
$RUNFDS $INDIR Beyler_Hood_toluene_165
$RUNFDS $INDIR Beyler_Hood_toluene_166
$RUNFDS $INDIR Beyler_Hood_toluene_167
$RUNFDS $INDIR Beyler_Hood_toluene_168
$RUNFDS $INDIR Beyler_Hood_toluene_170


