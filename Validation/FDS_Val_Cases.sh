#!/bin/bash

# To include validation data sets to be run automatically by Validationbot,
# be sure to include the data set in ~/FDS-SMV/Validation/Process_All_Output.sh

$QFDS -m 1 -p 204 -n 8 -d  UMD_Line_Burner/FDS_Input_Files methane_dx_p3125cm.fds
$QFDS -m 1 -p 204 -n 8 -d  UMD_Line_Burner/FDS_Input_Files propane_dx_p3125cm.fds
