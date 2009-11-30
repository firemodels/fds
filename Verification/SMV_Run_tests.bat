@echo off

Rem Batch script to run FDS cases used to test Smokeview

echo About to run FDS cases used to test Smokeview
echo   (will take several hours depending on computer speed)

set wui="%CD%\Wui"
set vis="%CD%\Visualization"
set smvug="%CD%\..\Manuals\SMV_5_User_Guide\"
set logfile="%CD%\fds_run_tests.log"

echo Press CTRL c to abort
echo Press any other key to start FDS cases
pause>NUL

echo Starting Visualization verification cases > %logfile%
date /T >> %logfile%
time /T >> %logfile%
echo. >> %logfile%

echo | fds5 2> "%smvug%\figures\fds5.version"

cd %wui%
fds5 fire_line.fds

cd %vis%
fds5 colorbar.fds
fds5 colorconv.fds
fds5 devices_elem.fds
fds5 devices_vistest.fds
fds5 devcices_vistest2.fds
fds5 plume5a.fds
fds5 plume5b.fds
fds5 plume5c.fds
fds5 plume5c_bounddef.fds
fds5 sillytexture.fds
fds5 script_test.fds
fds5 smoke_sensor.fds
fds5 smoke_test.fds
fds5 smoke_test2.fds
fds5 sprinkler_many.fds
fds5 thouse5.fds
fds5 transparency.fds

echo Cases complete >> %logfile%
date /T >> %logfile%
time /T >> %logfile%



