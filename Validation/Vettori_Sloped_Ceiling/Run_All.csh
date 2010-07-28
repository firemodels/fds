#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`
cp $BASEDIR/FDS_Input_Files/sge-fds.sh $BASEDIR/Current_Results

$RUNFDS Current_Results Vettori_Sloped_FSSW fire51 &
$RUNFDS Current_Results Vettori_Sloped_FSSD fire51 &
$RUNFDS Current_Results Vettori_Sloped_FSSC fire61 &
$RUNFDS Current_Results Vettori_Sloped_FSFW fire61 &
$RUNFDS Current_Results Vettori_Sloped_FSFD fire62 &
$RUNFDS Current_Results Vettori_Sloped_FSFC fire62 &
$RUNFDS Current_Results Vettori_Sloped_FOFW fire63 &
$RUNFDS Current_Results Vettori_Sloped_FOFD fire63 &
$RUNFDS Current_Results Vettori_Sloped_FOFC fire64 &
$RUNFDS Current_Results Vettori_Sloped_FOSW fire64 &
$RUNFDS Current_Results Vettori_Sloped_FOSD fire65 &
$RUNFDS Current_Results Vettori_Sloped_FOSC fire65 &
$RUNFDS Current_Results Vettori_Sloped_13SSW fire66 &
$RUNFDS Current_Results Vettori_Sloped_13SSD fire66 &
$RUNFDS Current_Results Vettori_Sloped_13SSC fire67 &
$RUNFDS Current_Results Vettori_Sloped_13SFW fire67 &
$RUNFDS Current_Results Vettori_Sloped_13SFD fire68 &
$RUNFDS Current_Results Vettori_Sloped_13SFC fire68 &
$RUNFDS Current_Results Vettori_Sloped_13OFW fire52 &
$RUNFDS Current_Results Vettori_Sloped_13OFD fire52 &
$RUNFDS Current_Results Vettori_Sloped_13OFC fire53 &
$RUNFDS Current_Results Vettori_Sloped_13OSW fire53 &
$RUNFDS Current_Results Vettori_Sloped_13OSD fire54 &
$RUNFDS Current_Results Vettori_Sloped_13OSC fire54 &
$RUNFDS Current_Results Vettori_Sloped_24SSW fire55 &
$RUNFDS Current_Results Vettori_Sloped_24SSD fire55 &
$RUNFDS Current_Results Vettori_Sloped_24SSC fire74 &
$RUNFDS Current_Results Vettori_Sloped_24SFW fire74 &
$RUNFDS Current_Results Vettori_Sloped_24SFD fire75 &
$RUNFDS Current_Results Vettori_Sloped_24SFC fire75 &
$RUNFDS Current_Results Vettori_Sloped_24OFW fire76 &
$RUNFDS Current_Results Vettori_Sloped_24OFD fire76 &
$RUNFDS Current_Results Vettori_Sloped_24OFC fire77 &
$RUNFDS Current_Results Vettori_Sloped_24OSW fire77 &
$RUNFDS Current_Results Vettori_Sloped_24OSD fire79 &
$RUNFDS Current_Results Vettori_Sloped_24OSC fire79 &

echo FDS cases submitted
