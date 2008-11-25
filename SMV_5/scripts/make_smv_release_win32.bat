@echo off
set version=%1
REM
REM This batch file creates a self unarchiving file containing
REM a release version of smokeview and associated files
REM
REM usage: 
REM  make_smv_release_win32 version
REM    where version is of the from X.Y_svn#

cd ..\for_bundle
set zipbase=smv_%version%_win32
set smvdir=to_google\%zipbase%

echo
echo filling distribution directory
mkdir %smvdir%
mkdir %smvdir%\Documentation
copy smokeview_release.exe %smvdir%\smokeview.exe
copy smokezip_release.exe %smvdir%\smokezip.exe
copy devices.svo %smvdir%\.
copy readme.html %smvdir%\Documentation\.
copy ..\..\Manuals\All_PDF_Files\SMV_5_User_Guide.pdf %smvdir%\Documentation\.

echo
echo winzipping distribution directory
cd %smvdir%
wzzip -a -r -P %zipbase%.zip *

echo
echo creating self-extracting archive
d:\bin\winzip\wzipse32 %zipbase%.zip -d "c:\program files\nist\smokeview"
copy %zipbase%.exe ..\.

cd ..\..\..\scripts