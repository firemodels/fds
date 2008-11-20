@echo off
set revision=2700
REM
REM create a directory containing files needed for a
REM test version of smokeview
REM
cd ..\for_bundle
set smvdir=smv_5_test_%revision%_win32
mkdir %smvdir%
copy smokeview_test.exe %smvdir%\smokeview.exe
copy smokezip_release.exe %smvdir%\smokezip.exe
copy devices.svo %smvdir%\.
copy note.txt %smvdir%\.
wzzip -a %smvdir%.zip %smvdir%\*
d:\bin\winzip\wzipse32 %smvdir% -d "c:\program files\nist\smokeview"