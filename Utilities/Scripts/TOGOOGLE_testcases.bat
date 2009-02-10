@echo off

Rem Windows batch file to upload FDS test cases

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


Rem set level=Release-1_Major
Rem set level=Release-2_Minor
set level=Release-3_Maintenance

Rem directory containing googlecode password file (gc.passwd)
set pwdir=c:\bin\

Rem -------- should not need to edit any lines below ---------

set winfile=fds_test_cases_%test_cases_revision%.exe
set unixfile=fds_test_cases_%test_cases_revision%.tar.gz

echo Uploading %winfile%
echo Uploading %unixfile%
pause

  set glabels=Type-Installer,Opsys-Windows,%level%
  set dplatform=Windows
  set summary=Test cases (SVN r%test_cases_revision% - Windows line endings)
  echo.
  echo Uploading %summary% - %wincase%
       googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %winfile%

  set glabels=Type-Archive,Opsys-Unix,%level%
  set dplatform=Linux/OSX
  set summary=Test cases (SVN r%test_cases_revision% - Unix line endings)
  echo.
  echo Uploading %summary% - %wincase%
       googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %unixfile%

echo.
echo Uploads complete
pause
