@echo off
Title Cleaning FDS for 32 bit Windows 

Rem Batch file used to clean 32 and 64 bit FDS build directories

set envfile=%homedrive%\%homepath%\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist
pause
call %envfile%

%svn_drive%
echo.
echo cleaning intel_win_32
cd %svn_root%\FDS_Compilation\intel_win_32
set out=intel_win_32.out
date /t | tee %out%
time /t | tee -a %out%
echo Cleaning intel_win_32 | tee -a %out%

make -f ..\makefile clean | tee -a %out%

echo.
echo cleaning mpi_intel_win_32
cd %svn_root%\FDS_Compilation\mpi_intel_win_32
set out=%svn_root%\FDS_Compilation\mpi_intel_win_32\mpi_intel_win_32.out
date /t | tee %out%
time /t | tee -a %out%
echo Cleaning mpi_intel_win_32 | tee -a %out%

make -f ..\makefile clean | tee -a %out%

pause