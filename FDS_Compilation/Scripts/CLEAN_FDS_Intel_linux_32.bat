@echo off
Title Cleaning FDS for 32 bit Linux

Rem Batch file used to build a 32 bit version of FDS

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

Rem location of batch files used to set up Intel compilation environment

call %envfile%

set scriptdir=%linux_svn_root%/FDS_Compilation/Scripts
set fdsdir=%linux_svn_root%/FDS_Compilation/

plink %svn_logon% %scriptdir%/MAKE_fds_onhost.csh intel_linux_32  %fdsdir%/intel_linux_32 %linux_hostname% clean

plink %svn_logon% %scriptdir%/MAKE_fds_onhost.csh mpi_intel_linux_32  %fdsdir%/mpi_intel_linux_32 %linux_hostname% clean

pause