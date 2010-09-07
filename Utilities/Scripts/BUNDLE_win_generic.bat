@echo off

set fdsdir=%svn_root%\FDS_Compilation\intel_win_%platform%
set fdsmpidir=%svn_root%\FDS_Compilation\mpi_intel_win_%platform%
set basename=FDS_%fds_version%_SMV_%smv_version%_win%platform%

set in_pdf=%svn_root%\Manuals\All_PDF_Files
set in_fds2ascii=%svn_root%\Utilities\fds2ascii
set in_smokediff=%svn_root%\Utilities\smokediff
Rem set in_smokezip=%svn_root%\Utilities\smokezip
set in_smokezip=%svn_root%\SMV_5\for_bundle
set in_background=%svn_root%\Utilities\background
set in_smv=%svn_root%\SMV_5\for_bundle\

set to_google=%svn_root%\Utilities\to_google
set out_bundle=%to_google%\%basename%\FDS
set out_bin=%out_bundle%\FDS5\bin
set out_textures=%out_bin%\textures
set out_uninstall=%out_bundle%\FDS5\Uninstall
set out_doc=%out_bundle%\FDS5\Documentation
set out_guides="%out_doc%\Guides_and_Release_Notes"
set out_web="%out_doc%\FDS_on_the_Web"
set out_examples=%out_bundle%\FDS5\Examples

set manifest=%out_bin%\RELEASE_INFO.txt
set bundleinfo=%svn_root%\Utilities\Scripts\bundle_setup
set wikify=%svn_root%\Utilities\Scripts\wikify.py

Rem erase the temporary bundle directory if it already exists

if exist %out_bundle% rmdir /s /q %out_bundle%
echo making directories
mkdir %out_bundle%
mkdir %out_bin%
mkdir %out_textures%
mkdir %out_doc%
mkdir %out_guides%
mkdir %out_web%
mkdir %out_examples%
mkdir %out_uninstall%

Rem Copy FDS, Smokeview and other needed files to the bin  directory

echo.
echo Copying FDS to the bin directory
if "%platform%"=="32" copy %fdsdir%\fds5_win_%platform%.exe         %out_bin%\fds5.exe
if "%platform%"=="32" copy %fdsmpidir%\fds5_mpi_win_%platform%.exe  %out_bin%\fds5_mpi.exe
if "%platform%"=="64" copy %fdsdir%\fds5_win_%platform%.exe         %out_bin%\.
if "%platform%"=="64" copy %fdsmpidir%\fds5_mpi_win_%platform%.exe  %out_bin%\.

echo Copying Smokeview to the bin directory
copy %in_smv%\smokeview%platform%_release.exe   %out_bin%\smokeview.exe

echo Copying smokediff to the bin directory
if "%platform%"=="32" copy %in_smokediff%\intel_win_32\smokediff.exe     %out_bin%\smokediff.exe
if "%platform%"=="64" copy %in_smokediff%\intel_win_64\smokediff64.exe   %out_bin%\smokediff_win_64.exe

echo Copying smokezip to the bin directory
Rem if "%platform%"=="32" copy %in_smokezip%\intel_win_32\smokezip.exe     %out_bin%\smokezip.exe
if "%platform%"=="32" copy %in_smokezip%\smokezip32_release.exe     %out_bin%\smokezip.exe
if "%platform%"=="64" copy %in_smokezip%\intel_win_64\smokezip64.exe   %out_bin%\smokezip_win_64.exe

echo Copying fds2ascii to the bin directory
if "%platform%"=="32" copy %in_fds2ascii%\intel_win_32\fds2ascii_win_32.exe     %out_bin%\fds2ascii.exe
if "%platform%"=="64" copy %in_fds2ascii%\intel_win_64\fds2ascii_win_64.exe     %out_bin%\fds2ascii_win_64.exe

echo Copying background to the bin directory
copy %in_background%\intel_win_32\background.exe %out_bin%\background.exe

echo.
echo Creating Manifiest
echo FDS-Smokeview bundle created > %manifest%
date /t >> %manifest%
time /t >> %manifest%
echo. >> %manifest%
echo Versions:>> %manifest%
echo. >> %manifest%

if "%platform%"=="64" GOTO endif_32
echo -------------------------- >> %manifest%
echo | %out_bin%\fds5.exe 2>> %manifest%
echo. >> %manifest%
echo -------------------------- >> %manifest%
echo | %out_bin%\fds5_mpi.exe 2>> %manifest%

echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\fds2ascii.exe -v >>%manifest%

echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\smokeview.exe -v >> %manifest%
echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\smokediff.exe -v >> %manifest%
echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\smokezip.exe -v >> %manifest%

echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\background.exe -v >> %manifest%
:endif_32

if "%platform%"=="32" GOTO endif_64 
echo -------------------------- >> %manifest%
echo | %out_bin%\fds5_win_%platform%.exe 2>> %manifest%
echo. >> %manifest%
echo -------------------------- >> %manifest%
echo | %out_bin%\fds5_mpi_win_%platform%.exe 2>> %manifest%

echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\fds2ascii_win_64.exe -v >> %manifest%

echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\smokeview.exe -v >> %manifest%
echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\smokediff_win_64.exe >> %manifest%
echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\smokezip_win_64.exe >> %manifest%

echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\background.exe -v >> %manifest%
:endif_64


echo.
echo Copying auxillary files to the bin directory
copy %in_smv%\objects.svo             %out_bin%\.
copy %in_smv%\pthreadVC.dll           %out_bin%\.
copy %in_smv%\glew32.dll              %out_bin%\.
copy %in_smv%\smokeview.ini           %out_bin%\.
echo.
echo Copying textures to the bin\textures directory
copy %in_smv%\textures\*.jpg          %out_textures%\.
copy %in_smv%\textures\*.png          %out_textures%\.

echo.
echo Copying Uninstaller to Uninstall directory
copy "%bundleinfo%\uninstall_fds5.bat"             "%out_uninstall%\Uninstall.bat"
copy "%bundleinfo%\set_path%platform%.exe"         "%out_uninstall%\set_path.exe"

Rem Include documentation in the bundle only if the variable, docs_include_in_bundles,
Rem is not set to 0.  This variable is defined in the fds_smv_env.bat setup  file

echo.
echo Getting the FDS release notes from the repository
svn export --quiet --force http://fds-smv.googlecode.com/svn/wiki/FDS_Release_Notes.wiki "%bundleinfo%\FDS_Release_Notes.wiki"

echo.
echo Converting the FDS release notes from wiki to html format
"%wikify%" -r "%bundleinfo%\FDS_Release_Notes.wiki" > "%out_guides%\FDS_Release_Notes.htm"

echo.
echo Copying Documentation to the Documentation directory

copy %in_pdf%\FDS_5_User_Guide.pdf               %out_guides%\.
copy %in_pdf%\SMV_5_User_Guide.pdf               %out_guides%\.
copy %in_pdf%\SMV_5_Technical_Reference_Guide.pdf %out_guides%\.
copy %in_pdf%\FDS_5_Technical_Reference_Guide.pdf %out_guides%\.
copy "%in_smv%\readme.html"                      "%out_guides%\Smokeview_release_notes.html"

echo.
echo Copying web page shortcuts
copy "%bundleinfo%\Overview.html"             "%out_doc%\Overview.html"
copy "%bundleinfo%\FDS_Web_Site.url"          "%out_web%\Official_Web_Site.url"
copy "%bundleinfo%\Updates.url"               "%out_web%\Software_Updates.url"
copy "%bundleinfo%\Docs.url"               "%out_web%\Documentation_Updates.url"
copy "%bundleinfo%\FDS_Development_Web_Site.url" "%out_web%\Developer_Web_Site.url"
copy "%bundleinfo%\discussion_group.url"          "%out_web%\Discussion_Group.url"
copy "%bundleinfo%\issue_tracker.url"          "%out_web%\Issue_Tracker.url"

Rem Copy readme_examples file to Examples directory to let user download all examples

echo.
echo Copying readme_examples.html to the Examples directory
copy %bundleinfo%\readme_examples.html "%out_examples%\Examples notes.html"
echo.
echo Getting the Verification cases from the repository
svn export --quiet --force https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Verification %out_examples%

Rem echo.
Rem echo Make a generic copy of the installation directory
Rem if exist "%out_bundle%\..\..\FDS_SMOKEVIEW\" rmdir /s /q "%out_bundle%\..\..\FDS_SMOKEVIEW\"
Rem mkdir "%out_bundle%\..\..\FDS_SMOKEVIEW\"
Rem xcopy /E "%out_bundle%\..\FDS\*" "%out_bundle%\..\..\FDS_SMOKEVIEW\"
Rem if exist "%out_bundle%\..\..\FDS_SMOKEVIEW\FDS5\Uninstall" rmdir /s /q "%out_bundle%\..\..\FDS_SMOKEVIEW\FDS5\Uninstall"

echo.
echo Copying wrapup scripts for use in final installation
copy "%bundleinfo%\wrapup_fds_install.bat" "%out_bundle%\FDS5\wrapup_fds_install.bat"

copy "%bundleinfo%\shortcut.exe" "%out_bundle%\FDS5\shortcut.exe"
copy "%bundleinfo%\set_path.exe" "%out_bundle%\FDS5\set_path.exe"

echo.
echo Compressing FDS/Smokeview distribution

cd %to_google%
if exist %basename%.zip erase %basename%.zip
cd %basename%\fds\fds5\
wzzip -a -r -el -xExamples\*.csv -P ..\..\..\%basename%.zip *

Rem create an installation file from the zipped bundle directory

echo.
echo Creating installer
cd %to_google%
echo Setup is about to install FDS %fds_version% and Smokeview %smv_version% > %bundleinfo%\message.txt
echo Press Setup to begin installation. > %bundleinfo%\main.txt
if exist %basename%.exe erase %basename%.exe
wzipse32 %basename%.zip -runasadmin -a %bundleinfo%\about.txt -st"FDS-Smokeview Setup" -d "c:\Program Files\FDS\FDS5" -c wrapup_fds_install.bat
Rem wzipse32 -setup -a %bundleinfo%\about.txt -st"FDS-Smokeview Setup" -t %bundleinfo%\main.txt -mokcancel %bundleinfo%\message.txt %basename%.zip -c wrapup_fds_install.bat
