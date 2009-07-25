@echo off

Rem Script to bundle fds and smokeview into an installation file

set envfile="%homedrive%\%homepath%"\fds_smv_env.bat
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

Rem define variables used by script

set in_pdf=%svn_root%\Manuals\All_PDF_Files
set in_fds=%svn_root%\FDS_Compilation\intel_win_32
set in_fds_mpi=%svn_root%\FDS_Compilation\mpi_intel_win_32
set in_smv=%svn_root%\SMV_5\for_bundle\

set to_google=%svn_root%\Utilities\to_google
set basename=FDS_%fds_version%-SMV_%smv_version%_Windows
set out_bundle=%to_google%\%basename%
set out_bin=%out_bundle%\FDS5\bin
set out_doc=%out_bundle%\FDS5\Documentation
set out_examples=%out_bundle%\FDS5\Examples

set bundleinfo=%svn_root%\Utilities\Scripts\bundle_setup

Rem erase the temporary bundle directory if it already exists

if exist %out_bundle% rmdir /s /q %out_bundle%
mkdir %out_bundle%
mkdir %out_bin%
mkdir %out_doc%
mkdir %out_examples%

Rem Copy FDS, Smokeview and other needed files to the bin  directory

echo.
echo Copying files to bin directory
copy %in_fds%\fds5_win_32.exe %out_bin%\fds5.exe
copy %in_fds_mpi%\fds5_win_mpi_32 %out_bin%\fds5_mpi.exe
copy %in_smv%\smokeview_release.exe %out_bin%\smokeview.exe
copy %in_smv%\devices.svo %out_bin%\.
copy %in_smv%\pthreadVC.dll %out_bin%\.
copy %in_smv%\smokezip_release.exe %out_bin%\smokezip.exe
copy %in_smv%\glew32.dll %out_bin%\.
copy %in_smv%\smokeview.ini %out_bin%\.

Rem Include documentation in the bundle only if the variable, docs_include_in_bundles,
Rem is not set to 0.  This variable is defined in the fds_smv_env.bat setup  file

if %docs_include_in_bundle% EQU 0 goto end_docs
echo.
echo Copying documentation and readme_docs.html to the Documentation directory
copy %in_pdf%\FDS_5_Technical_Reference_Guide.pdf %out_doc%\.
copy %in_pdf%\FDS_5_Validation_Guide.pdf          %out_doc%\.
copy %in_pdf%\SMV_5_Verification_Guide.pdf        %out_doc%\.
copy %in_pdf%\SMV_5_Technical_Reference_Guide.pdf %out_doc%\.
copy %in_pdf%\SMV_5_User_Guide.pdf                %out_doc%\.
copy %in_pdf%\FDS_5_User_Guide.pdf                %out_doc%\.
copy %in_pdf%\FDS_5_Verification_Guide.pdf        %out_doc%\.
:end_docs

copy %bundleinfo%\readme_docs.html %out_doc%\readme_docs.html

Rem Copy readme_examples file to Examples directory to let user download all examples

echo.
echo Copying readme_examples.html to the Examples directory
copy %bundleinfo%\readme_examples.html %out_examples%\readme_examples.html

copy %bundleinfo%\wrapup.bat %out_bundle%\FDS5\wrapup.bat
copy %bundleinfo%\set_fds5_path.exe %out_bundle%\FDS5\set_fds5_path.exe

Rem compress bundle directory

cd %to_google%
if exist %basename%.zip erase %basename%.zip
cd %basename%\fds5\
wzzip -a -r -P ..\..\%basename%.zip *

Rem create an installation file from the zipped bundle directory

echo.
echo Creating installer
cd %to_google%
echo Setup is about to install FDS %fds_version% and Smokeview %smv_version% > %bundleinfo%\message.txt
echo Press Setup to begin installation. > %bundleinfo%\main.txt
if exist %basename%.exe erase %basename%.exe
c:\bin\winzip\wzipse32 %basename%.zip -d "c:\program files\fds5" -c wrapup.bat
Rem c:\bin\winzip\wzipse32 -setup -a %bundleinfo%\about.txt -st"FDS-Smokeview Setup" -t %bundleinfo%\main.txt -mokcancel %bundleinfo%\message.txt %basename%.zip -c wrapup.bat
pause