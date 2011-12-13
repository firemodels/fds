@echo off
Title Cleaning FDS for 64 bit Windows 

Rem Batch file used to clean Windows 64 bit FDS build directories

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

%svn_drive%
echo.
cd %svn_root%\FDS_Compilation\Intel_Win_64
echo Cleaning intel_win_64

erase *.obj 
erase *.mod

echo.
cd %svn_root%\FDS_Compilation\mpi_intel_win_64
echo Cleaning mpi_intel_win_64

erase *.obj 
erase *.mod

pause