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
IF EXIST %smvdir% rmdir /S /Q %smvdir%
mkdir %smvdir%
copy smokeview.ini %smvdir%\smokeview.ini
copy smokeview32_release.exe %smvdir%\smokeview.exe
copy smokediff32_release.exe %smvdir%\smokediff.exe
copy smokezip32_release.exe %smvdir%\smokezip.exe
copy background.exe %smvdir%\background.exe
copy objects.svo %smvdir%\.
copy glew32.dll %smvdir%\.
copy pthreadVC.dll %smvdir%\.
copy readme.html %smvdir%\release_notes.html

echo
echo winzipping distribution directory
cd %smvdir%
wzzip -a -r -P %zipbase%.zip *

echo
echo creating self-extracting archive
wzipse32 %zipbase%.zip -d "c:\program files\fds\fds5\bin"
copy %zipbase%.exe ..\.

cd ..\..\..\scripts
