@ECHO OFF
set dir=%1
set infile=%2

set fulldir=%curdir%\..\%dir%

cd %fulldir%
IF EXIST %infile%.smv (
smokeview -runscript %infile%
) ELSE (
echo %infile%.smv does not exist
)
cd %curdir%

