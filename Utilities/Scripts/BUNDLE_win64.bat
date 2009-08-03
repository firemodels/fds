@echo off

Rem bundle 64 bit versions of fds and smokeview into an installation file

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
set in_fds=%svn_root%\FDS_Compilation\intel_win_64
set in_fds_mpi=%svn_root%\FDS_Compilation\mpi_intel_win_64
set in_fds2ascii=%svn_root%\Utilities\fds2ascii
set in_smv=%svn_root%\SMV_5\for_bundle\

set to_google=%svn_root%\Utilities\to_google
set basename=FDS_%fds_version%_SMV_%smv_version%_win64
set out_bundle=%to_google%\%basename%\FDS
set out_bin=%out_bundle%\FDS5\bin
set out_doc=%out_bundle%\FDS5\Documentation
set out_examples=%out_bundle%\FDS5\Examples

set fds5=fds5.exe
set fds5mpi=fds5_mpi.exe
set smokeview=smokeview.exe
set smokezip=smokezip.exe

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
copy %in_fds%\fds5_win_64.exe %out_bin%\%fds5%
copy %in_fds_mpi%\fds5_win_mpi_64.exe %out_bin%\%fds5mpi%
copy %in_fds2ascii%\fds2ascii.exe %out_bin%\fds2ascii.exe
copy %in_smv%\smokeview64_release.exe %out_bin%\%smokeview%
copy %in_smv%\devices.svo %out_bin%\.
copy %in_smv%\pthreadVC.dll %out_bin%\.
copy %in_smv%\smokezip_release.exe %out_bin%\%smokezip%
copy %in_smv%\glew32.dll %out_bin%\.
copy %in_smv%\smokeview.ini %out_bin%\.

Rem Include documentation in the bundle only if the variable, docs_include_in_bundles,
Rem is not set to 0.  This variable is defined in the fds_smv_env.bat setup  file

echo.
echo Copying FDS and Smokeview users guide and other documentation 
echo to the Documentation directory

copy %in_pdf%\FDS_5_User_Guide.pdf        %out_doc%\.
copy %in_pdf%\SMV_5_User_Guide.pdf        %out_doc%\.
copy %bundleinfo%\readme_docs.html        %out_doc%\readme_docs.html
copy "%bundleinfo%\FDS Web Site.url"      %out_doc%\.
copy "%bundleinfo%\FDS Development Web Site.url"       %out_doc%\.
copy %in_smv%\readme.html %out_doc%\readme_smokeview.html
copy %bundleinfo%\readme_fds.url %out_doc%\readme_fds.url

if %docs_include_in_bundle% EQU 0 goto end_docs
echo.
echo Copying all other documentation to the Documentation directory
copy %in_pdf%\FDS_5_Technical_Reference_Guide.pdf %out_doc%\.
copy %in_pdf%\FDS_5_Validation_Guide.pdf          %out_doc%\.
copy %in_pdf%\FDS_5_Verification_Guide.pdf        %out_doc%\.
copy %in_pdf%\SMV_5_Verification_Guide.pdf        %out_doc%\.
copy %in_pdf%\SMV_5_Technical_Reference_Guide.pdf %out_doc%\.
:end_docs

Rem Copy readme_examples file to Examples directory to let user download all examples

echo.
echo Copying readme_examples.html and the examples to the Examples directory
copy %bundleinfo%\readme_examples.html %out_examples%\readme_examples.html
svn export --force https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Verification %out_examples%

copy %bundleinfo%\wrapup_fds_install.bat %out_bundle%\FDS5\wrapup_fds_install.bat
copy %bundleinfo%\shortcut.exe %out_bundle%\FDS5\shortcut.exe
copy %bundleinfo%\set_path.exe %out_bundle%\FDS5\set_path.exe

Rem compress bundle directory

cd %to_google%
if exist %basename%.zip erase %basename%.zip
cd %basename%\fds\fds5\
wzzip -a -r -xExamples\*.csv -P ..\..\..\%basename%.zip *

Rem create an installation file from the zipped bundle directory

echo.
echo Creating installer
cd %to_google%
echo Setup is about to install FDS %fds_version% and Smokeview %smv_version% > %bundleinfo%\message.txt
echo Press Setup to begin installation. > %bundleinfo%\main.txt
if exist %basename%.exe erase %basename%.exe
c:\bin\winzip\wzipse32 %basename%.zip -a %bundleinfo%\about.txt -st"FDS-Smokeview Setup" -d "c:\Program Files\FDS\FDS5" -c wrapup_fds_install.bat
Rem c:\bin\winzip\wzipse32 -setup -a %bundleinfo%\about.txt -st"FDS-Smokeview Setup" -t %bundleinfo%\main.txt -mokcancel %bundleinfo%\message.txt %basename%.zip -c wrapup_fds_install.bat
pause