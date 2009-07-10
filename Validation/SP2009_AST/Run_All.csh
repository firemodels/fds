#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results SP2009_AST_Test_1 fire57 &
$RUNFDS Current_Results SP2009_AST_Test_2 fire58 &
$RUNFDS Current_Results SP2009_AST_Test_3 fire59 &

echo FDS cases submitted
