@echo off

Rem Windows batch file to upload FDS release to the google download site.

Rem uncomment 2 of the following 3 lines depending on what type of release is
Rem being uploaded

  Rem set level=Release-1_Major
  Rem set level=Release-2_Minor
  set level=Release-3_Maintenance

Rem directory containing googlecode password file (gc.passwd)
  set pwdir=c:\bin\

set envfile=c:\bin\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%

%svn_drive%
cd %svn_root%\Utilities\Scripts\to_google


  set glabels=Type-Installer,Opsys-Windows,%level%
  set dplatform=32 bit Windows
  set platform=win32
  set summary=FDS %fds_version% for %dplatform% (SVN r%fds_revision%)
  set exe=fds_%fds_version%_%fds_revision%_%platform%.exe
echo Uploading %exe%
echo FDS %level% version=%fds_version% revision=%fds_revision%
echo.
echo press any key to proceed with upload, CTRL c to abort
pause>NUL

  if not exist %exe% goto abort_upload
  echo.
  echo Uploading %summary% - %exe%
  echo googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
       googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%

echo.
echo Uploads complete
pause
goto:eof

:abort_upload
echo %exe% does not exist - upload to google failed
pause
