#!/bin/csh -f

# To generate a batch file used to run the below cases on Windows 
# do the following

# 1. make a copy of this file and call it FDS_Cases.bat
# 2. open FDS_Cases.bat in a text editor and change every
#    occurrence of $RUNFDS to %RUNFDS%
# 3. Then run the batch file Run_All.bat

$RUNFDS Surf_mass surf_mass_part_char_cart_fuel     fire61
$RUNFDS Surf_mass surf_mass_part_char_cart_gas      fire61
$RUNFDS Surf_mass surf_mass_part_char_cyl_fuel      fire61
$RUNFDS Surf_mass surf_mass_part_char_cyl_gas       fire61
$RUNFDS Surf_mass surf_mass_part_char_spher_fuel    fire62
$RUNFDS Surf_mass surf_mass_part_char_spher_gas     fire62
$RUNFDS Surf_mass surf_mass_part_nonchar_cart_fuel  fire62
$RUNFDS Surf_mass surf_mass_part_nonchar_cart_gas   fire62
$RUNFDS Surf_mass surf_mass_part_nonchar_cyl_fuel   fire63
$RUNFDS Surf_mass surf_mass_part_nonchar_cyl_gas    fire63
$RUNFDS Surf_mass surf_mass_part_nonchar_spher_fuel fire63
$RUNFDS Surf_mass surf_mass_part_nonchar_spher_gas  fire63
$RUNFDS Surf_mass surf_mass_vent_char_cart_fuel     fire64
$RUNFDS Surf_mass surf_mass_vent_char_cart_gas      fire64
$RUNFDS Surf_mass surf_mass_vent_char_cyl_fuel      fire64
$RUNFDS Surf_mass surf_mass_vent_char_cyl_gas       fire64
$RUNFDS Surf_mass surf_mass_vent_char_spher_fuel    fire65
$RUNFDS Surf_mass surf_mass_vent_char_spher_gas     fire65
$RUNFDS Surf_mass surf_mass_vent_nonchar_cart_fuel  fire65
$RUNFDS Surf_mass surf_mass_vent_nonchar_cart_gas   fire65
$RUNFDS Surf_mass surf_mass_vent_nonchar_cyl_fuel   fire66
$RUNFDS Surf_mass surf_mass_vent_nonchar_cyl_gas    fire66
$RUNFDS Surf_mass surf_mass_vent_nonchar_spher_fuel fire66
$RUNFDS Surf_mass surf_mass_vent_nonchar_spher_gas  fire66

