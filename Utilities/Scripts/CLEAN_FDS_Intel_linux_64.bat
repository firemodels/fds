@echo off
Title Cleaning FDS for 64 bit Linux

Rem Batch file used to build a 64 bit version of FDS

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

call %envfile%

set scriptdir=%linux_svn_root%/Utilities/Scripts

set target=intel_linux_64
set fdsdir=%linux_svn_root%/Utilities/Makefile/Intel_Linux_64
plink %svn_logon% %scriptdir%/MAKE_fds_onhost.csh %target% %fdsdir% %LINUXCOMPILE% clean

set target=mpi_intel_linux_64
set fdsdir=%linux_svn_root%/Utilities/Makefile/Mpi_Intel_Linux_64
plink %svn_logon% %scriptdir%/MAKE_fds_onhost.csh %target% %fdsdir% %LINUXCOMPILE% clean

pause
