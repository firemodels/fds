@echo off
setlocal EnableDelayedExpansion
set platform=%1

Rem  Windows batch file to build a test Smokeview for Windows 64

Rem setup environment variables (defining where repository resides etc) 

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

call %envfile%
%svn_drive%
echo.
echo   Installing FDS-Smokeview bundle on %platform%
Title  Installing FDS-Smokeview bundle on %platform%

if "%platform%" == "windows" (
  cd %svn_root%\Utilities\uploads
  call FDS_%fds_version%-SMV_%smv_version%_win64.exe
  goto eof
)
if "%platform%" == "linux" (
  plink %linux_logon% %linux_svn_root%/SMV/scripts/run_command.sh Utilities/uploads FDS_%fds_version%-SMV_%smv_version%_linux64.sh y
  goto eof
)
if "%platform%" == "osx" (
  plink %osx_logon% %linux_svn_root%/SMV/scripts/run_command.sh Utilities/uploads FDS_%fds_version%-SMV_%smv_version%_osx64.sh y
  goto eof
)

:eof
echo.
echo installation complete
pause
