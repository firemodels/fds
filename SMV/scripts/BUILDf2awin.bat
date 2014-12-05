@echo off

::  Windows batch file to build 64 bit windows version of fds2ascii

:: setup environment variables (defining where repository resides etc) 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%

echo.
echo ---------------------------------------------
echo Building 64 bit Windows versions of fds2ascii
echo ---------------------------------------------

cd %svn_root%\Utilities\fds2ascii\intel_win_64
call %svn_root%\Utilities\Scripts\setup_intel_compilers.bat

Title Building fds2ascii for 64 bit Windows

erase *.obj *.exe
ifort -o fds2ascii_win_64.exe -O2 ..\..\Data_processing\fds2ascii.f90
pause
