@echo off
Title Building Parallel FDS for 64 bit Windows 
	       
Rem Batch file used to build a 32 bit version of FDS

set envfile="%homedrive%\%homepath%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%

%svn_drive%

call %svn_root%\FDS_Compilation\Scripts\SET_INTEL_64.bat

cd %svn_root%\FDS_Compilation\mpi_intel_win_64

set out=mpi_intel_win_64.out
echo. | tee -a %out%
date /t | tee -a  %out%
time /t | tee -a  %out%
make VPATH="../../FDS_Source" -f ..\makefile mpi_intel_win_64 | tee -a %out%
pscp %out% %svn_logon%:%linux_svn_root%/FDS_Compilation/mpi_intel_win_64/mpi_intel_win_64.out

pause
