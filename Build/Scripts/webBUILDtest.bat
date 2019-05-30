@echo off
set platform=%1

:: batch file to build FDS on Windows, Linux or OSX platforms

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
echo  Building test_mpi for %platform%
Title Building test_mpi for %platform%

%svn_drive%

if "%platform%" == "windows" (
  cd %svn_root%\fds\Utilities\test_mpi\impi_intel_win
  erase *.obj *.mod *.exe
  call make_test_mpi
  goto eof
)
if "%platform%" == "linux" (
  plink %plink_options% %linux_logon% %linux_svn_root%/smv/scripts/clean.sh fds/Utilities/test_mpi/impi_intel_linux
  plink %plink_options% %linux_logon% %linux_svn_root%/smv/scripts/run_command.sh fds/Utilities/test_mpi/impi_intel_linux make_test_mpi.sh
  pause
  goto eof
)
if "%platform%" == "osx" (
  plink %plink_options% %linux_logon% %linux_svn_root%/smv/scripts/clean.sh fds/Utilities/test_mpi/mpi_intel_osx
  plink %plink_options% %osx_logon% %linux_svn_root%/smv/scripts/run_command.sh fds/Utilities/test_mpi/mpi_intel_osx make_test_mpi.sh
  pause
  goto eof
)

:eof
