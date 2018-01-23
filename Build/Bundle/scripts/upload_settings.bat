@echo off
setlocal EnableDelayedExpansion
set platform=%1

set bashscript=env_params.sh

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

echo #!/bin/bash > %bashscript%

echo.  >> %bashscript%
echo # ---- FDS and smokeview version ---- >> %bashscript%
echo.  >> %bashscript%

echo export fds_version=%fds_version% >> %bashscript%
echo export smv_version=%smv_version% >> %bashscript%

echo.  >> %bashscript%
echo #  ---- repo locations ---- >> %bashscript%
echo.  >> %bashscript%

echo.  >> %bashscript%
echo #  *** Linux/OSX >> %bashscript%
echo.  >> %bashscript%

echo export linux_svn_root=%linux_svn_root% >> %bashscript%
echo export compiler_dir=%compiler_dir% >> %bashscript%
echo export misc_dir=%misc_dir% >> %bashscript%

echo.  >> %bashscript%
echo #  ---- MPI library locations ---- >> %bashscript%
echo.  >> %bashscript%

echo export linux_mpi_version=%linux_mpi_version% >> %bashscript%
echo export osx_mpi_version=%osx_mpi_version% >> %bashscript%

echo.  >> %bashscript%
echo #  ---- bot locations ---- >> %bashscript%
echo.  >> %bashscript%

echo export firebotrepo=%firebotrepo% >> %bashscript%
echo export smokebotrepo=%smokebotrepo% >> %bashscript%

echo.  >> %bashscript%
echo #  ---- hostnames ---- >> %bashscript%
echo.  >> %bashscript%

echo.  >> %bashscript%
echo #  *** linux >> %bashscript%
echo.  >> %bashscript%

echo export linux_hostname=%linux_hostname% >> %bashscript%
echo export linux_username=%linux_username% >> %bashscript%
echo export linux_logon=%linux_logon% >> %bashscript%

echo.  >> %bashscript%
echo #  *** osx >> %bashscript%
echo.  >> %bashscript%
echo export osx_hostname=%osx_hostname% >> %bashscript%
echo export osx_username=%osx_username% >> %bashscript%
echo export osx_logon=%osx_logon% >> %bashscript%

sed "s/^M$//" %bashscript% > FDS_SMV_ENVpc.sh

if "%platform%" == "linux" (
pscp FDS_SMV_ENVpc.sh %linux_hostname%:.
)

if "%platform%" == "osx" (
pscp FDS_SMV_ENVpc.sh %osx_hostname%:.
)

erase %bashscript%

pause
