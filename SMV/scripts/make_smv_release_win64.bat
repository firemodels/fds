@echo off
set version=%1
REM
REM This batch file creates a self unarchiving file containing
REM a release version of smokeview and associated files
REM
REM usage: 
REM  make_smv_release_win64 version
REM    where version is of the from X.Y_svn#

cd ..\for_bundle
set zipbase=smv_%version%_win64
set smvdir=uploads\%zipbase%

echo
echo filling distribution directory
IF EXIST %smvdir% rmdir /S /Q %smvdir%
mkdir %smvdir%
copy smokeview.ini %smvdir%\smokeview.ini
copy smokeview64_release.exe %smvdir%\smokeview.exe
copy smokezip64_release.exe %smvdir%\smokezip_win_64.exe
copy background.exe %smvdir%\background.exe
copy smokediff64_release.exe %smvdir%\smokediff_win_64.exe
copy objects.svo %smvdir%\.
copy glew32_x64.dll %smvdir%\.
copy pthreadVC2_x64.dll %smvdir%\.
copy readme.html %smvdir%\release_notes.html

echo
echo winzipping distribution directory
cd %smvdir%
wzzip -a -r -P %zipbase%.zip *

echo
echo creating self-extracting archive
wzipse32 %zipbase%.zip -d "c:\Program Files\FDS\%fds_edition%\bin"
copy %zipbase%.exe ..\.

cd ..\..\..\scripts
