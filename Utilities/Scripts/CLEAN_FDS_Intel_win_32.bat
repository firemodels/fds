@echo off
Title Cleaning FDS for 32 bit Windows 

Rem Batch file used to clean 32 and 64 bit FDS build directories

set envfile=%userprofile%\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

Rem location of batch files used to set up Intel compilation environment
set intelbin=c:\bin

call %envfile%

%svn_drive%
echo.
echo cleaning Intel_Win_32
cd %svn_root%\Utilities\Makefile\Intel_Win_32
set out=intel_win_32.out
date /t | tee %out%
time /t | tee -a %out%
echo Cleaning intel_win_32 | tee -a %out%

make -f ..\makefile clean | tee -a %out%

echo.
echo cleaning Mpi_Intel_Win_32
cd %svn_root%\Utilities\Makefile\Mpi_Intel_Win_32
set out=%svn_root%\Utilities\Makefile\Mpi_Intel_Win_32\mpi_intel_win_32.out
date /t | tee %out%
time /t | tee -a %out%
echo Cleaning mpi_intel_win_32 | tee -a %out%

make -f ..\makefile clean | tee -a %out%

cd ..\..\Scripts
pause