@echo off
set version=%1
REM
REM This batch file creates a self unarchiving file containing
REM a test version of smokeview and associated files
REM
REM usage: 
REM  make_smv_release_win64 version
REM    where version is of the from test_svn#

cd ..\for_bundle
set zipbase=smv_%version%_win64
set smvdir=to_google\%zipbase%

echo
echo filling distribution directory
IF EXIST %smvdir% rmdir /S /Q %smvdir%
mkdir %smvdir%
copy smokeview64_test.exe %smvdir%\smokeview.exe
copy smokezip_release.exe %smvdir%\smokezip.exe
copy devices.svo %smvdir%\.
copy glew32.dll %smvdir%\.
copy pthreadVC.dll %smvdir%\.
copy note.txt %smvdir%\.

echo
echo winzipping distribution directory
cd %smvdir%
wzzip -a %zipbase%.zip *

echo
echo creating self-extracting archive
c:\bin\winzip\wzipse32 %zipbase%.zip -d "c:\program files\nist\smokeview"
copy %zipbase%.exe ..\.

cd ..\..\..\scripts