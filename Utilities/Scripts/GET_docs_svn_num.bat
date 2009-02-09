@echo off

Rem Batch file used to svn info for FDS, Smokeview and the SVN test cases

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

echo.

cd %svn_root%\Manuals\All_PDF_Files
svn update
echo ---------------------------
echo All_PDF_Files svn info
echo ---------------------------
svn info


pause