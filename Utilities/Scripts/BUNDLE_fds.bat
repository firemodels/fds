@echo off

Rem Script to bundle fds and smokeview into an installation file

set envfile=%homedrive%\%homepath%\fds_smv_env.bat
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

set in_pdf=%svn_root%\Manuals\All_PDF_Files
set in_fds=%svn_root%\Utilities\to_google\fds_%fds_version%_%fds_revision%_win32
set in_smv=%svn_root%\SMV_5\for_bundle\to_google\smv_%smv_version%_%smv_revision%_win32

set to_google=%svn_root%\Utilities\to_google
set basename=FDS_%fds_version%-SMV_%smv_version%_Windows
set out_bundle=%to_google%\%basename%
set out_bin=%out_bundle%\FDS5\bin
set out_doc=%out_bundle%\FDS5\Documentation
set out_examples=%out_bundle%\FDS5\Examples

set bundleinfo=%svn_root%\Utilities\Scripts\bundle_setup

if exist %out_bundle% rmdir /s /q %out_bundle%
mkdir %out_bundle%
mkdir %out_bin%
mkdir %out_doc%

echo.
echo Copying files to bin directory
copy %in_fds%\fds5.exe %out_bin%\.
copy %in_smv%\smokeview.exe %out_bin%\.
copy %in_smv%\objects.svo %out_bin%\.
copy %in_smv%\pthreadVC.dll %out_bin%\.
copy %in_smv%\smokezip.exe %out_bin%\.
copy %in_smv%\glew32.dll %out_bin%\.
copy %in_smv%\smokeview.ini %out_bin%\.

echo.
echo Copying files to Documentation directory

copy %in_pdf%\FDS_5_Technical_Reference_Guide.pdf %out_doc%\.
copy %in_pdf%\FDS_5_Validation_Guide.pdf          %out_doc%\.
copy %in_pdf%\SMV_5_User_Guide.pdf                %out_doc%\.
copy %in_pdf%\FDS_5_User_Guide.pdf                %out_doc%\.
copy %in_pdf%\FDS_5_Verification_Guide.pdf        %out_doc%\.

echo.
echo Creating Examples directory

svn -q export https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Verification %out_examples%
set turbfile=%out_examples%\Decaying_Isotropic_Turbulence
if exist %turbfile% rmdir /S /Q %turbfile%

echo. 
echo Compressing %out_bundle%

copy %bundleinfo%\setup.bat %out_bundle%\FDS5\setup.bat

cd %to_google%
if exist %basename%.zip erase %basename%.zip
cd %basename%
wzzip -a -r -P ..\%basename%.zip *

echo.
echo Creating installer
cd ..
echo Setup is about to install FDS %fds_version% and Smokeview %smv_version% > %bundleinfo%\message.txt
echo Press Setup to begin installation. > %bundleinfo%\main.txt
if exist %basename%.exe erase %basename%.exe
Rem c:\bin\winzip\wzipse32 %basename%.zip -d "c:\program files\nist\"
c:\bin\winzip\wzipse32 -setup -a %bundleinfo%\about.txt -st"FDS-Smokeview Setup" -t %bundleinfo%\main.txt -m %bundleinfo%\message.txt %basename%.zip -c FDS5\setup.bat
pause