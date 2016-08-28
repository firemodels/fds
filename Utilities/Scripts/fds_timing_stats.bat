@echo off

::This script checks the .out files for the FDS Verification Suite and generates
::a .csv file of the CPU time (and other metrics) called fds_timing_stats.csv

set CURDIR=%CD%
cd ..\..
set SVNROOT=%CD%

cd %SVNROOT%\FDS\Verification
::   OUT_FILES=*/*.out

:: Write header information to fds_timing_stats.csv file
echo 'FDS Case,Wall Clock Time (s),CPU Time (s),Number of Cells,Number of Time Steps,Performance Metric (1e-6)' > %SVNROOT%\Utilities\Scripts\fds_timing_stats.csv

