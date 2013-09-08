@echo off

Rem Windows batch file to upload 64 bit Smokeview osx files to
Rem the google download site.

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

Rem --------------- 64 bit OSX ----------------

  set glabels=Type-Archive,OpSys-OSX_64,%level%
  set dplatform=64 bit OSX
  set platform=osx64
  set summary=Smokeview %smv_version% for %dplatform% (SVN r%smv_revision%)
  set exe=smv_%smv_version%_%platform%.tar.gz
  echo.
  echo Uploading %exe% 
       %upload% --passwd-file-dir %google_password_dir% --config-dir none  -s "%summary%" -p fds-smv -u %google_username% -l %glabels% %exe%

echo.
echo Upload complete
pause
