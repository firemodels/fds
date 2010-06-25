@echo off

set svnroot=d:\fds-smv
set casedir=%svnroot%\Verification\SMV_Script_Example\case

echo Press any key to generate images 
pause>NUL

d:
cd %casedir%

smokeview -runscript plume5c 