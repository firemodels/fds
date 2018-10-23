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

call :title2bash # ---- FDS and smokeview version ----
call :dos2bash fds_version %fds_version%
call :dos2bash smv_version %smv_version%
call :dos2bash fds_debug   %fds_debug%
call :dos2bash FDSEDITION  %fds_edition%

call :title2bash #  ---- MPI version ----
call :dos2bash linux_mpi_version %linux_mpi_version%
call :dos2bash osx_mpi_version   %osx_mpi_version%

call :title2bash #  ---- INTEL compiler version ----
call :dos2bash INTELVERSION %INTELVERSION%

call :title2bash #  ---- Linux/OSX repo locations ----
call :dos2bash linux_svn_root %linux_svn_root%
call :dos2bash INTEL_LIB_DIR  %INTEL_LIB_DIR%
call :dos2bash INTEL_BIN_DIR  %INTEL_BIN_DIR%
call :dos2bash OS_LIB_DIR     %OS_LIB_DIR%

call :title2bash #  shouldn't have to change anything below

call :title2bash #  ---- Guide location ----
call :dos2bash GUIDE_DIR $HOME/%GUIDE_DIR%

call :title2bash #  ---- openmpi location ----
call :dos2bash OPENMPI_DIR $HOME/%OPENMPI_DIR%

call :title2bash #  ---- bundle location ----
call :dos2bash BUNDLE_DIR $HOME/%BUNDLE_DIR%

call :title2bash #  ---- bot locations ----
call :dos2bash firebotrepo %firebotrepo%
call :dos2bash firebothome %firebothome%
echo.  >> %bashscript%
call :dos2bash smokebotrepo %smokebotrepo%
call :dos2bash smokebothome %smokebothome%

call :title2bash #  ---- Linux login info ----
call :dos2bash linux_hostname `hostname`
call :dos2bash linux_username `whoami`
call :dos2bash linux_logon    $linux_username@$linux_hostname

call :title2bash #  ---- OSX login info ----
call :dos2bash osx_hostname `hostname`
call :dos2bash osx_username `whoami`
call :dos2bash osx_logon    $osx_username@$osx_hostname

if "%platform%" == "linux" (
  pscp %bashscript% %linux_hostname%:.bundle/FDS_SMV_ENVpc.sh
  plink %linux_logon% %linux_svn_root%/fds/Build/Bundle/scripts/dos2unix.sh .bundle FDS_SMV_ENVpc.sh
)

if "%platform%" == "osx" (
  pscp %bashscript% %osx_hostname%:.bundle/FDS_SMV_ENVpc.sh
  plink %osx_logon% %linux_svn_root%/fds/Build/Bundle/scripts/dos2unix.sh .bundle FDS_SMV_ENVpc.sh
)

:: -------------------------------------------------------------
 :title2bash
:: -------------------------------------------------------------

set val=%*
echo.  >> %bashscript%
echo %val% >> %bashscript%
echo.  >> %bashscript%
exit /b /0


erase %bashscript%
goto :eof

:: -------------------------------------------------------------
 :dos2bash
:: -------------------------------------------------------------

set var=%1
set val=%2
echo export %var%=%val% >> %bashscript%
exit /b /0


erase %bashscript%
goto :eof

:eof

pause
