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

echo #!/bin/bash > FDS_SMV_ENV.sh

echo.  >> FDS_SMV_ENV.sh
echo # ---- FDS and smokeview version ---- >> FDS_SMV_ENV.sh
echo.  >> FDS_SMV_ENV.sh

echo export fds_version=%fds_version% >> FDS_SMV_ENV.sh
echo export smv_version=%smv_version% >> FDS_SMV_ENV.sh

echo.  >> FDS_SMV_ENV.sh
echo #  ---- repo locations ---- >> FDS_SMV_ENV.sh
echo.  >> FDS_SMV_ENV.sh

echo.  >> FDS_SMV_ENV.sh
echo #  *** Linux/OSX >> FDS_SMV_ENV.sh
echo.  >> FDS_SMV_ENV.sh

echo export linux_svn_root=%linux_svn_root% >> FDS_SMV_ENV.sh
echo export compiler_dir=%compiler_dir% >> FDS_SMV_ENV.sh
echo export misc_dir=%misc_dir% >> FDS_SMV_ENV.sh

echo.  >> FDS_SMV_ENV.sh
echo #  ---- MPI library locations ---- >> FDS_SMV_ENV.sh
echo.  >> FDS_SMV_ENV.sh

echo export linux_mpi_version=%linux_mpi_version% >> FDS_SMV_ENV.sh
echo export osx_mpi_version=%osx_mpi_version% >> FDS_SMV_ENV.sh

echo.  >> FDS_SMV_ENV.sh
echo #  ---- bot locations ---- >> FDS_SMV_ENV.sh
echo.  >> FDS_SMV_ENV.sh

echo export firebotrepo=%firebotrepo% >> FDS_SMV_ENV.sh
echo export smokebotrepo=%smokebotrepo% >> FDS_SMV_ENV.sh

echo.  >> FDS_SMV_ENV.sh
echo #  ---- hostnames ---- >> FDS_SMV_ENV.sh
echo.  >> FDS_SMV_ENV.sh

echo.  >> FDS_SMV_ENV.sh
echo #  *** linux >> FDS_SMV_ENV.sh
echo.  >> FDS_SMV_ENV.sh

echo export linux_hostname=%linux_hostname% >> FDS_SMV_ENV.sh
echo export linux_username=%linux_username% >> FDS_SMV_ENV.sh
echo export linux_logon=%linux_logon% >> FDS_SMV_ENV.sh

echo.  >> FDS_SMV_ENV.sh
echo #  *** osx >> FDS_SMV_ENV.sh
echo.  >> FDS_SMV_ENV.sh
echo export osx_hostname=%osx_hostname% >> FDS_SMV_ENV.sh
echo export osx_username=%osx_username% >> FDS_SMV_ENV.sh
echo export osx_logon=%osx_logon% >> FDS_SMV_ENV.sh

