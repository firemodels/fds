@echo off

::  batch file to copy the smokeview source git log to the the smokeview readme file

:: setup environment variables (defining where repository resides etc) 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/Scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%

%svn_drive%
set curdir=%CD%
cd %svn_root%\SMV\Source
echo git log --pretty=oneline --abbrev-commit --since=%smvlogdate% 
pause
git log --pretty=oneline --abbrev-commit --since=%smvlogdate% .>> %userprofile%\FDS-SMVwebpages\smv_readme.html
cd %curdir%
pause
