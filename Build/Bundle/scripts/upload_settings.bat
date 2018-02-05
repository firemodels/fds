@echo off
setlocal EnableDelayedExpansion
set platform=%1

:: batch file to install the FDS-SMV bundle on Windows, Linux or OSX systems

:: setup environment variables (defining where repository resides etc) 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use smv/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%

%svn_drive%
cd %svn_root%\fds\Build\Bundle\scripts

set bashscript=%svn_root%\fds\Build\Bundle\scripts\env_params.sh

echo #!/bin/bash > %bashscript%

:: -------FDS and smokeview version---------------

echo.  >> %bashscript%
echo # ---- FDS and smokeview version ---- >> %bashscript%
echo.  >> %bashscript%

echo export fds_version=%fds_version% >> %bashscript%
echo export smv_version=%smv_version% >> %bashscript%

:: ----------MPI version-----------------------

echo.  >> %bashscript%
echo #  ---- MPI version ---- >> %bashscript%
echo.  >> %bashscript%

echo export linux_mpi_version=%linux_mpi_version% >> %bashscript%
echo export osx_mpi_version=%osx_mpi_version% >> %bashscript%

:: ----------repo locations---------------

echo.  >> %bashscript%
echo #  ---- Linux/OSX repo locations ---- >> %bashscript%
echo.  >> %bashscript%

echo export linux_svn_root=%linux_svn_root% >> %bashscript%
echo export compiler_dir=%compiler_dir% >> %bashscript%
echo export misc_dir=%misc_dir% >> %bashscript%

echo.  >> %bashscript%
echo #  -------------------------------- >> %bashscript%
echo #  shouldn't have to change anything below >> %bashscript%

:: ---------Guide location------------------------

echo.  >> %bashscript%
echo #  ---- Guide location ---- >> %bashscript%
echo.  >> %bashscript%

echo export GUIDE_DIR=$HOME/%GUIDE_DIR% >> %bashscript%

:: ---------openmpi location----------------------

echo.  >> %bashscript%
echo #  ---- openmpi location ---- >> %bashscript%
echo.  >> %bashscript%

echo export OPENMPI_DIR=$HOME/%OPENMPI_DIR% >> %bashscript%

:: ---------bundle location----------------------

echo.  >> %bashscript%
echo #  ---- bundle location ---- >> %bashscript%
echo.  >> %bashscript%

echo export BUNDLE_DIR=$HOME/%BUNDLE_DIR% >> %bashscript%

:: ----------bot locations---------------------

echo.  >> %bashscript%
echo #  ---- bot locations ---- >> %bashscript%
echo.  >> %bashscript%

echo export firebotrepo=%firebotrepo% >> %bashscript%
echo export firebothome=%firebothome% >> %bashscript%
echo.
echo export smokebotrepo=%smokebotrepo% >> %bashscript%
echo export smokebothome=%smokebothome% >> %bashscript%

:: ----------linux hostnames---------------------

echo.  >> %bashscript%
echo #  ---- Linux login info ---- >> %bashscript%
echo.  >> %bashscript%

echo export linux_hostname=%linux_hostname% >> %bashscript%
echo export linux_username=%linux_username% >> %bashscript%
echo export linux_logon=%linux_logon% >> %bashscript%

echo.  >> %bashscript%
echo #  ---- OSX login info ---- >> %bashscript%
echo.  >> %bashscript%
echo export osx_hostname=%osx_hostname% >> %bashscript%
echo export osx_username=%osx_username% >> %bashscript%
echo export osx_logon=%osx_logon% >> %bashscript%

if "%platform%" == "linux" (
  pscp %bashscript% %linux_hostname%:FDS_SMV_ENVpc.sh
  plink %linux_logon% %linux_svn_root%/fds/Build/Bundle/scripts/dos2unix.sh FDS_SMV_ENVpc.sh
)

if "%platform%" == "osx" (
  pscp %bashscript% %osx_hostname%:FDS_SMV_ENVpc.sh
  plink %osx_logon% %linux_svn_root%/fds/Build/Bundle/scripts/dos2unix.sh FDS_SMV_ENVpc.sh
)

erase %bashscript%

pause
