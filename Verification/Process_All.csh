#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
cd Flowfields
$SVNROOT/Utilities/Data_Processing/moody
cd ..
cd NS_Analytical_Solution
$SVNROOT/Utilities/Data_Processing/ns2d
cd ..
cd Radiation
$SVNROOT/Utilities/Data_Processing/radiation_box
$SVNROOT/Utilities/Data_Processing/radiation_plane_layer
cd ..


