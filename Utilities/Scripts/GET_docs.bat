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
set docbasedir=%togoogle%\%zipbase%\

mkdir %docbasedir%

copy %mandir%\FDS_5_Technical_Reference_Guide.pdf %docbasedir%\.
copy %mandir%\FDS_5_User_Guide.pdf %docbasedir%\.
copy %mandir%\FDS_5_Validation_Guide.pdf %docbasedir%\.
copy %mandir%\FDS_5_Verification_Guide.pdf %docbasedir%\.

copy %mandir%\SMV_5_User_Guide.pdf %docbasedir%\.

cd %docbasedir%
if exist %togoogle%\%zipbase%.zip erase %togoogle%\%zipbase%.zip
wzzip -a -r -P %togoogle%\%zipbase%.zip *

echo.
echo creating self-extracting archive

cd %togoogle%
if exist %zipbase%.exe erase %zipbase%.exe
c:\bin\winzip\wzipse32 %zipbase%.zip -d "C:\Program Files\nist\Documentation"

echo Documentation self-extractor, %zipbase%.exe, 
echo located in %togoogle%

set LREPOS=FDS-SMV

Rem get Documentation for Linux/OSX systems

plink %svn_logon% %LREPOS%/Utilities/Scripts/get_docs.csh %docs_revision%
pscp %svn_logon%:%LREPOS%/Utilities/Scripts/to_google/%zipbase%.tar.gz %togoogle%\.

cd %svn_root%\Utilities\Scripts
pause