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
fds5 fire_line.fds

cd %vis%
fds5 colorconv.fds
fds5 plume5a.fds
fds5 plume5b.fds
fds5 plume5c.fds
fds5 plume5c_bounddef.fds
fds5 sillytexture.fds
fds5 script_test.fds
fds5 smoke_sensor.fds
fds5 smoke_test.fds
fds5 smoke_test2.fds
fds5 thouse5.fds

