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

Rem un-comment 1 of the following 3 lines

Rem set level=Release-1_Major
Rem set level=Release-2_Minor
set level=Release-3_Maintenance

Rem directory containing googlecode password file (gc.passwd)
set pwdir=c:\bin\

set windocs=fds_smv_docs_%docs_revision%.exe
echo Uploading %windocs%
pause

  set glabels=Type-Installer,Opsys-Windows,%level%
  set dplatform=Windows
  set summary=FDS and Smokeview documentation for %dplatform% (SVN r%docs_revision%)
  set file=%windocs%
  echo.
  echo Uploading %summary% 
       googlecode_upload.py --passwd-file-dir %pwdir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %file%

echo -----
echo Uploads complete
pause
