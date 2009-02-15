@echo off
echo "Press <CTRL> c to abort"
echo "Press any other key to start FDS cases"
pause
fds5 colorconv
fds5 plume5a.fds
fds5 plume5b.fds
fds5 plume5c.fds
fds5 sillytexture.fds
fds5 smoke_sensor.fds
fds5 smoke_test.fds
fds5 smoke_test2.fds
fds5 thouse5.fds