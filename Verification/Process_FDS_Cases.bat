@echo off
echo.
echo Processing...
set DPDIR=..\..\Utilities\Data_Processing
cd NS_Analytical_Solution
"%DPDIR%\ns2d"
cd ..\Radiation
"%DPDIR%\radiation_box"
"%DPDIR%\radiation_plane_layer"
"%DPDIR%\wall_internal_radiation"
cd ..
pause

