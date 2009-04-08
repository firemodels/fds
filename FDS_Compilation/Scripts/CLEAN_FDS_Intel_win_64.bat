@echo off
Title Cleaning FDS for 64 bit Windows 

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

Rem location of batch files used to set up Intel compilation environment
set intelbin=c:\bin

call %envfile%

%svn_drive%
echo.
echo cleaning Intel_Win_64
cd %svn_root%\FDS_Compilation\Intel_Win_64
set out=intel_win_64.out
date /t | tee %out%
time /t | tee -a %out%
echo Cleaning intel_win_64 | tee -a %out%
make -f ..\makefile clean | tee -a %out%
pscp %out% %svn_logon%:%linux_svn_root%/FDS_Compilation/intel_win_64/.

echo.
echo cleaning mpi_intel_win_64
cd %svn_root%\FDS_Compilation\mpi_intel_win_64
set out=mpi_intel_win_64.out
date /t | tee %out%
time /t | tee -a %out%
echo Cleaning mpi_intel_win_64 | tee -a %out%
make -f ..\makefile clean | tee -a %out%
pscp %out% %svn_logon%:%linux_svn_root%/FDS_Compilation/mpi_intel_win_64/.

pause