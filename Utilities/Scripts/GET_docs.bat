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

set mandir=%svn_root%\Manuals\All_PDF_Files
set zipbase=fds_smv_docs_%docs_revision%
set togoogle=%svn_root%\Utilities\Scripts\to_google
set docbasedir=%togoogle%\%zipbase%\nist\
set fdsuploaddir=%docbasedir%\FDS\Documentation
set smvuploaddir=%docbasedir%\Smokeview\Documentation

mkdir %fdsuploaddir%
mkdir %smvuploaddir%

copy %mandir%\FDS_5_Technical_Reference_Guide.pdf %fdsuploaddir%\.
copy %mandir%\FDS_5_User_Guide.pdf %fdsuploaddir%\.
copy %mandir%\FDS_5_Validation_Guide.pdf %fdsuploaddir%\.
copy %mandir%\FDS_5_Verification_Guide.pdf %fdsuploaddir%\.

copy %mandir%\SMV_5_User_Guide.pdf %smvuploaddir%\.

cd %docbasedir%
if exist %togoogle%\%zipbase%.zip erase %togoogle%\%zipbase%.zip
wzzip -a -r -P %togoogle%\%zipbase%.zip *

echo.
echo creating self-extracting archive

cd %togoogle%
if exist %zipbase%.exe erase %zipbase%.exe
d:\bin\winzip\wzipse32 %zipbase%.zip -d "C:\Program Files\nist\"

echo Documentation self-extractor, %zipbase%.exe, 
echo located in %togoogle%

cd %svn_root%\Utilities\Scripts
pause