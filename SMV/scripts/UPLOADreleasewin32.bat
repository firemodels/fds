@echo off

Rem Windows batch file to upload 32 bit Smokeview windows archive

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
cd %svn_root%\smv\for_bundle\uploads

set level=Release-3_Maintenance

Rem --------------- 32 bit Windows ----------------

set glabels=Type-Installer,OpSys-Windows,%level%
set dplatform=32 bit Windows
set platform=win32
set summary=Smokeview %smv_version% for %dplatform% (SVN r%smv_revision%)
set exe=smv_%smv_version%_%platform%.exe
echo.
  echo Uploading %exe% 
     %upload% --passwd-file-dir %google_password_dir% --config-dir none -s "%summary%" -p fds-smv -u %google_username% -l %glabels% %exe%

echo.
echo Upload complete
pause
