@echo off

Rem Batch script to run FDS cases used to test Smokeview

echo About to run FDS cases used to test Smokeview
echo   (will take several hours depending on computer speed)

set wui="%CD%\Wui"
set vis="%CD%\Visualization"
set smvug="%CD%\..\Manuals\SMV_5_User_Guide\"

echo Press CTRL c to abort
echo Press any other key to start FDS cases
pause>NUL

echo | fds5 2> %smvug%\scriptfigures\fds5.version

cd %wui%
fds5 fire_line.fdv

cd %vis%
fds5 colorconv.fdv
fds5 plume5a.fdv
fds5 plume5b.fdv
fds5 plume5c.fdv
fds5 plume5c_bounddef.fdv
fds5 sillytexture.fdv
fds5 script_test.fdv
fds5 smoke_sensor.fdv
fds5 smoke_test.fdv
fds5 smoke_test2.fdv
fds5 thouse5.fdv

