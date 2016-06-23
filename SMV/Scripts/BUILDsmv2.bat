@echo off
set whichsmv=%1

::  batch file to build a debug windows or incremental Windows/Linux/OSX smokeview

:: setup environment variables (defining where repository resides etc) 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/Scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%

%svn_drive%


if "%whichsmv%" == "win_smvd" (
  echo.
  echo  Building windows debug smokeview
  Title Building windows debug smokeview
  
  cd %svn_root%\SMV\Build\smokeview\intel_win_64
  call make_smv_db -r
  goto eof
)
if "%whichsmv%" == "win_smvtd" (
  echo.
  echo  Building windows test debug smokeview
  Title Building windows test debug smokeview

  cd %svn_root%\SMV\Build\smokeview\intel_win_64
   call make_smv_db -t
  goto eof
)
if "%whichsmv%" == "linux_smvti" (
  echo.
  echo  Building linux smokeview incrementally
  Title Building linux smokeview incrementally
  
  plink %linux_logon% %linux_svn_root%/SMV/Scripts/run_command.sh SMV/Build/smokeview/intel_linux_64  make_smv_inc.sh
  goto eof
)
if "%whichsmv%" == "osx_smvti" (
  echo.
  echo  Building OSX smokeview incrementally
  Title Building OSX smokeview incrementally

  plink %osx_logon% %linux_svn_root%/SMV/Scripts/run_command.sh SMV/Build/smokeview/intel_osx_64  make_smv_inc.sh
  goto eof
)
if "%whichsmv%" == "win_smvti" (
  echo.
  echo  Building windows smokeview incrementally
  Title Building windows smokeview incrementally

  cd %svn_root%\SMV\Build\smokeview\intel_win_64
  call make_smv_inc -t
  goto eof
)

:eof
echo.
echo compilation complete
pause
