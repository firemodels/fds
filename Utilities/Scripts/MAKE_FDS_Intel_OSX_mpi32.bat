@echo off
Title Building Parallel FDS for 32 bit OSX

Rem Batch file used to build a 32 bit version of FDS

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

call %envfile%

set target=mpi_intel_osx_32
set fdsdir=%linux_svn_root%/Utilities/Makefile/Mpi_Intel_OSX_32
set scriptdir=%linux_svn_root%/Utilities/Scripts

plink %svn_logon% %scriptdir%/MAKE_fds_onhost.csh %target% %fdsdir% devi1.nist.gov

pause