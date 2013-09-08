@echo off

Rem Windows batch file to upload Smokeview test files to
Rem the google download site.  This script assume that the Windows
Rem batch file, MAKEtest.bat, has already been run.

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
cd %svn_root%\SMV\for_bundle\uploads

Rem ----------------------------------------------------------
Rem should not need to edit any lines below

set level=Release-3_Maintenance

set version=%smv_version%

Rem --------------- 64 bit Windows ----------------

set glabels=Type-Installer,OpSys-Windows_64,%level%
set dplatform=64 bit Windows
set platform=win64
set summary=Smokeview %smv_version% for %dplatform% (SVN r%smv_revision%)
set exe=smv_%version%_%platform%.exe
echo.
  echo Uploading %exe% 
     %upload% --passwd-file-dir %google_password_dir% --config-dir none -s "%summary%" -p fds-smv -u %google_username% -l %glabels% %exe%

echo.
echo Upload complete
pause
