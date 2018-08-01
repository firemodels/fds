@echo off

:: Produce fds input files with terrains in both formats, ZVALS (&GEOM) and &OBSTs.
:: This routine is a specialized version of smv/Build/dem2fds/data/test3.sh, for Gatlinburg test case.
:: Should run from both Blaze and Burn.

:: Relative location of compiled dem2fds. Assumes both fds and smv repos are in the same directory.
set dem2fds=..\..\..\..\smv\Build\dem2fds\intel_win_64\dem2fds_win_64

:: terrain directory containing DEM (Digital Elevation Model) and image data. 
set terraindir=%userprofile%\terrain


%dem2fds% -fds   Gatlinburg_1000m.fds       -dir %terraindir%\gatlinburg Gatlinburg_1000m.in
%dem2fds% -fds Gatlinburg_1000m_g.fds -geom -dir %terraindir%\gatlinburg Gatlinburg_1000m.in
