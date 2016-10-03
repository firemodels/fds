@echo off

:: build FDS on Windows, Linux and OSX platforms

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
echo.

%svn_drive%

echo ----------------------------------------------------------------------------
echo building windows fds
echo. 
cd %svn_root%\fds\Build\mpi_intel_win_64
erase *.obj *.mod *.exe
call make_fds

echo ----------------------------------------------------------------------------
echo building Linux fds
echo. 
plink %linux_logon% %linux_svn_root%/smv/scripts/run_command.sh fds/Build/mpi_intel_linux_64 clean_fds.sh
plink %linux_logon% %linux_svn_root%/smv/scripts/run_command.sh fds/Build/mpi_intel_linux_64 make_fds.sh

echo ----------------------------------------------------------------------------
echo building Linux OSX
echo. 
plink %osx_logon% %linux_svn_root%/smv/scripts/run_command.sh fds/Build/mpi_intel_osx_64 clean_fds.sh
plink %osx_logon% %linux_svn_root%/smv/scripts/run_command.sh fds/Build/mpi_intel_osx_64 make_fds.sh

:eof
