#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
cd NS_Analytical_Solution
$SVNROOT/Utilities/Data_Processing/ns2d
cd ..
cd Radiation
$SVNROOT/Utilities/Data_Processing/radiation_box
$SVNROOT/Utilities/Data_Processing/radiation_plane_layer
$SVNROOT/Utilities/Data_Processing/wall_internal_radiation
cd ..


