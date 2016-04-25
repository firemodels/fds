@echo off
set platform=%1
set buildtype=%2

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
echo.
echo  Building %type% Smokeview for 64 bit %platform%
Title Building %type% Smokeview for 64 bit %platform%

%svn_drive%

set type=
if "%buildtype%" == "test" (
   set type=-t
)
if "%buildtype%" == "release" (
   set type=-r
)

if "%platform%" == "windows" (
  cd %svn_root%\SMV\Build\smokeview\intel_win_64
  call make_smv %type%
  goto eof
)
if "%platform%" == "linux" (
  plink %linux_logon% %linux_svn_root%/SMV/scripts/run_command.sh SMV/Build/smokeview/intel_linux_64 make_smv.sh %type%
  goto eof
)
if "%platform%" == "osx" (
  plink %osx_logon% %linux_svn_root%/SMV/scripts/run_command.sh SMV/Build/smokeview/intel_osx_64 make_smv.sh %type%
  goto eof
)

:eof
echo.
echo compilation complete
pause
