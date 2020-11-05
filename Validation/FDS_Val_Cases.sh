#!/bin/bash

# These cases are run by firebot 1 time step when the -V test option is used
# to test the health of a Linux cluster

$QFDS -m 1 -p 204 -n 8 -d  UMD_Line_Burner/FDS_Input_Files methane_dx_p3125cm.fds
$QFDS -m 1 -p 204 -n 8 -d  UMD_Line_Burner/FDS_Input_Files propane_dx_p3125cm.fds
