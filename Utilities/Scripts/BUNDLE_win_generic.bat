@echo off

set fdsversion=%fds_edition%

set fdsdir=%svn_root%\FDS_Compilation\intel_win_%platform%
set fdsmpidir=%svn_root%\FDS_Compilation\mpi_intel_win_%platform%
set basename=FDS_%fds_version%_SMV_%smv_version%_win_%platform%

set in_pdf=%svn_root%\Manuals\All_PDF_Files
set in_fds2ascii=%svn_root%\Utilities\fds2ascii
set in_smokediff=%svn_root%\Utilities\smokediff
set in_smokezip=%svn_root%\Utilities\smokezip
set in_background=%svn_root%\Utilities\background
set in_smv=%svn_root%\SMV\Build\intel_win_%platform%
set in_for_bundle=%svn_root%\SMV\for_bundle

set to_google=%svn_root%\Utilities\to_google
set basedir=%to_google%\%basename%
set out_bundle=%basedir%\FDS
set out_bin=%out_bundle%\%fdsversion%\bin
set out_textures=%out_bin%\textures
set out_uninstall=%out_bundle%\%fdsversion%\Uninstall
set out_doc=%out_bundle%\%fdsversion%\Documentation
set out_guides="%out_doc%\Guides_and_Release_Notes"
set out_web="%out_doc%\FDS_on_the_Web"
set out_examples=%out_bundle%\%fdsversion%\Examples

set manifest=%out_bin%\manifest.html
set bundleinfo=%svn_root%\Utilities\Scripts\bundle_setup

Rem erase the temporary bundle directory if it already exists

if exist %basedir% rmdir /s /q %basedir%
echo making directories
mkdir %basedir%
mkdir %out_bundle%
mkdir %out_bundle%\%fdsversion%
mkdir %out_bin%
mkdir %out_textures%
mkdir %out_doc%
mkdir %out_guides%
mkdir %out_web%
mkdir %out_examples%
mkdir %out_uninstall%

set release_version=%fdssmv_major_version%_win_%platform%
set release_version=

echo.
echo copying fds_win_%platform%.exe to fds.exe
copy %fdsdir%\fds_win_%platform%.exe         %out_bin%\fds%release_version%.exe

echo copying fds_mpi_win_%platform%.exe to fds_mpi.exe
copy %fdsmpidir%\fds_mpi_win_%platform%.exe  %out_bin%\fds_mpi.exe

echo copying smokeview_win_%platform%.exe to smokeview%release_version%.exe
copy %in_smv%\smokeview_win_%platform%.exe   %out_bin%\smokeview%release_version%.exe

echo copying smokediff_win_%platform%.exe to smokediff%release_version%.exe
copy %in_smokediff%\intel_win_%platform%\smokediff_win_%platform%.exe     %out_bin%\smokediff%release_version%.exe

echo copying smokezip_win_%platform%.exe to smokezip%release_version%.exe
copy %in_smokezip%\intel_win_%platform%\smokezip_win_%platform%.exe       %out_bin%\smokezip%release_version%.exe 

echo copying fds2ascii_win_%platform%.exe to fds2ascii%release_version%.exe
copy %in_fds2ascii%\intel_win_%platform%\fds2ascii_win_%platform%.exe     %out_bin%\fds2ascii%release_version%.exe

echo copying background.exe
copy %in_background%\intel_win_32\background.exe %out_bin%\background.exe

echo.
echo Creating Manifiest

echo ^<html^> > %manifest%
echo ^<head^> >> %manifest%
echo ^<TITLE^>Build FDS^</TITLE^> >> %manifest%
echo ^</HEAD^> >> %manifest%
echo ^<BODY BGCOLOR="#FFFFFF" ^> >> %manifest%
echo ^<pre^> >> %manifest%
echo FDS-Smokeview bundle created >> %manifest%
date /t >> %manifest%
time /t >> %manifest%
echo. >> %manifest%
echo Versions:>> %manifest%
echo. >> %manifest%

echo -------------------------- >> %manifest%
echo | %out_bin%\fds%release_version%.exe 2>> %manifest%
echo. >> %manifest%
echo -------------------------- >> %manifest%
echo | %out_bin%\fds%release_version%.exe 2>> %manifest%

echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\fds2ascii%release_version%.exe -v >>%manifest%

echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\smokeview%release_version%.exe -v >> %manifest%
echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\smokediff%release_version%.exe -v >> %manifest%
echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\smokezip%release_version%.exe -v >> %manifest%

echo. >> %manifest%
echo -------------------------- >> %manifest%
%out_bin%\background.exe -v >> %manifest%

echo ^</body^> >> %manifest%
echo ^</html^> >> %manifest%

echo.
echo Copying auxillary files to the bin directory

echo copying objects.svo
copy %in_for_bundle%\objects.svo             %out_bin%\.

if "%platform%"=="32" echo copying pthreadVC.dll
if "%platform%"=="32" copy %in_for_bundle%\pthreadVC.dll           %out_bin%\.

if "%platform%"=="64" echo copying pthreadVC2_x64.dll
if "%platform%"=="64" copy %in_for_bundle%\pthreadVC2_x64.dll         %out_bin%\.

if "%platform%"=="32" echo copying glew32.dll
if "%platform%"=="32" copy %in_for_bundle%\glew32.dll              %out_bin%\.

if "%platform%"=="64" echo copying glew32_x64.dll
if "%platform%"=="64" copy %in_for_bundle%\glew32_x64.dll              %out_bin%\.

echo copying smokeview.ini
copy %in_for_bundle%\smokeview.ini           %out_bin%\.

echo.
echo Copying textures to the bin\textures directory
copy %in_for_bundle%\textures\*.jpg          %out_textures%\.
copy %in_for_bundle%\textures\*.png          %out_textures%\.

echo.
echo Copying Uninstaller to Uninstall directory
echo copying uninstall_fds.bat
copy "%bundleinfo%\uninstall_fds.bat" "%out_uninstall%\uninstall.bat"

echo copying set_path.exe
copy "%bundleinfo%\set_path.exe"         "%out_uninstall%\set_path.exe"

echo.
echo copying FDS_Release_Notes.htm
copy "%bundleinfo%\FDS_Release_Notes.htm" "%out_guides%\FDS_Release_Notes.htm"

echo.
echo copying Documentation to the Documentation directory

echo copying FDS_User_Guide.pdf
copy %in_pdf%\FDS_User_Guide.pdf %out_guides%\.

echo copying SMV_User_Guide.pdf
copy %in_pdf%\SMV_User_Guide.pdf %out_guides%\.

echo copying SMV_Technical_Reference_Guide.pdf
copy %in_pdf%\SMV_Technical_Reference_Guide.pdf %out_guides%\.

echo copying FDS_Technical_Reference_Guide.pdf
copy %in_pdf%\FDS_Technical_Reference_Guide.pdf %out_guides%\.

echo copying readme.html
copy "%in_for_bundle%\readme.html" "%out_guides%\Smokeview_release_notes.html"

echo.
echo copying web page shortcuts

echo copying Overview.html
copy "%bundleinfo%\Overview.html"             "%out_doc%\Overview.html"

echo copying FDS_Web_Site.url
copy "%bundleinfo%\FDS_Web_Site.url"          "%out_web%\Official_Web_Site.url"

echo copying Updates.url
copy "%bundleinfo%\Updates.url"               "%out_web%\Software_Updates.url"

echo copying Docs.url
copy "%bundleinfo%\Docs.url"               "%out_web%\Documentation_Updates.url"

echo copying FDS_Development_Web_Site.url
copy "%bundleinfo%\FDS_Development_Web_Site.url" "%out_web%\Developer_Web_Site.url"

echo copying discussion_group.url
copy "%bundleinfo%\discussion_group.url"          "%out_web%\Discussion_Group.url"

echo copying issue_tracker.url
copy "%bundleinfo%\issue_tracker.url"          "%out_web%\Issue_Tracker.url"

Rem Copy readme_examples file to Examples directory to let user download all examples

echo.
echo Copying readme_examples.html to the Examples directory
copy %bundleinfo%\readme_examples.html "%out_examples%\Examples notes.html"

echo.
echo Getting the Verification cases from the repository
svn export --quiet --force https://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Verification %out_examples%

echo.
echo Copying wrapup scripts for use in final installation

echo copying wrapup_fds_install_%platform%.bat
copy "%bundleinfo%\wrapup_fds_install_%platform%.bat" "%out_bundle%\%fdsversion%\wrapup_fds_install.bat
copy "%bundleinfo%\wrapup_fds_install_gen.bat" "%out_bundle%\%fdsversion%\wrapup_fds_install_gen.bat

echo copying shortcut.exe
copy "%bundleinfo%\shortcut.exe" "%out_bundle%\%fdsversion%\shortcut.exe"

echo copying set_path.exe
copy "%bundleinfo%\set_path.exe" "%out_bundle%\%fdsversion%\set_path.exe"

echo.
echo Compressing FDS/Smokeview distribution

cd %to_google%
if exist %basename%.zip erase %basename%.zip
cd %out_bundle%\%fdsversion%\
wzzip -a -r -xExamples\*.csv -P ..\..\..\%basename%.zip * > winzip.out

Rem create an installation file from the zipped bundle directory

echo.
echo Creating installer
cd %to_google%
echo Setup is about to install FDS %fds_version% and Smokeview %smv_version% > %bundleinfo%\message.txt
echo Press Setup to begin installation. > %bundleinfo%\main.txt
if exist %basename%.exe erase %basename%.exe
wzipse32 %basename%.zip -runasadmin -a %bundleinfo%\about.txt -st"FDS %fds_version% Smokeview %smv_version% Setup" -d "c:\Program Files\FDS\%fdsversion%" -c wrapup_fds_install.bat

start explorer %manifest%
type %manifest%
