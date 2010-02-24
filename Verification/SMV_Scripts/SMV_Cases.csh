#!/bin/csh -f

# To generate a batch file used to run the below cases on Windows 
# do the following

# 1. make a copy of this file and call it FDS_Cases.bat
# 2. open SMV_Cases.bat in a text editor and change every
#    occurrence of $RUNFDS to %RUNFDS%
# 3. Then run the batch file Run_All.bat


$RUNFDS Visualization colorbar                    fire47
$RUNFDS Visualization colorconv                   fire53
$RUNFDS Visualization objects_elem                fire53
$RUNFDS Visualization objects_static              fire53
$RUNFDS Visualization objects_dynamic             fire53
$RUNFDS Visualization objects_dynamic3            fire53
$RUNFDS Visualization plume5a                     fire53
$RUNFDS Visualization plume5b                     fire53
$RUNFDS Visualization plume5c                     fire53
$RUNFDS Visualization plume5c_bounddef            fire53
$RUNFDS Visualization sillytexture                fire45
$RUNFDS Visualization script_test                 fire45
$RUNFDS Visualization smoke_sensor                fire45
$RUNFDS Visualization smoke_test                  fire46
$RUNFDS Visualization smoke_test2                 fire45
$RUNFDS Visualization sprinkler_many              fire47
$RUNFDS Visualization thouse5                     fire47
$RUNFDS Visualization transparency                fire47
$RUNFDS Wui fire_line                             fire51
