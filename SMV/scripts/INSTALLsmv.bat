@echo off
set platform=%1
set buildtype=%2

::  batch file to build test or release Smokeview on a Windows, OSX or Linux system

:: setup environment variables (defining where repository resides etc) 

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
if "%buildtype%" == "test" (
  echo.
  echo  Installing test %platform% Smokeview
  Title  Installing test %platform% Smokeview
)
if "%buildtype%" == "release" (
  echo.
  echo  Installing %platform% Smokeview
  Title  Installing %platform% Smokeview
)

%svn_drive%

if "%platform%" == "windows" (
  cd %svn_root%\SMV\uploads
  if "%buildtype%" == "test" (
    echo Running Smokeview installer:  smv_test_%smv_revision%_win64.exe
    pause
    call smv_test_%smv_revision%_win64.exe
  )
  if "%buildtype%" == "release" (
    echo Running Smokeview installer: smv_%smv_version%_win64.exe
    pause
    call smv_%smv_version%_win64.exe
  )
  pause
  goto eof
)
if "%platform%" == "linux" (
  if "%buildtype%" == "test" (
    plink %linux_logon% %linux_svn_root%/SMV/scripts/run_command.sh SMV/uploads smv_test_%smv_revision%_linux64.sh y
  )
  if "%buildtype%" == "release" (
    plink %linux_logon% %linux_svn_root%/SMV/scripts/run_command.sh SMV/uploads smv_%smv_version%_linux64.sh y
  )
  goto eof
)
if "%platform%" == "osx" (
  if "%buildtype%" == "test" (
    plink %osx_logon% %linux_svn_root%/SMV/scripts/run_command.sh SMV/uploads smv_test_%smv_revision%_osx64.sh y
  )
  if "%buildtype%" == "release" (
    plink %osx_logon% %linux_svn_root%/SMV/scripts/run_command.sh SMV/uploads smv_%smv_version%_osx64.sh y
  )
  goto eof
)

:eof
echo.
echo compilation complete
pause
