@echo off

set svnroot=d:\fds-smv
set casedir=%svnroot%\Verification\SMV_Script_Example\case

echo Press any case to begin FDS simulation
pause>NUL

d:
cd %casedir%

fds5 plume5c.fds > plume5c.err 