@echo off

Rem Example illustrating how a script to bundle fds and smokeview for windows might be written

set REPOS=D:\fds-smv
set bundledir=D:\fds_bundle\NIST\Windows

Rem should not need to edit lines below

set fdsdir=%REPOS%\Utilities\Makefile
set scriptdir=%REPOS%\Utilities\Scripts
set pdfdir=%REPOS%\Manuals\All_PDF_Files
set smvbundle=%REPOS%\SMV_5\for_bundle
set exampledir=%bundledir%\Examples

Rem build FDS
echo Building FDS
cd %fdsdir%
echo %scriptdir%
call %scriptdir%\iclvars.bat ia32
call %scriptdir%\ifortvars.bat ia32
Rem erase *.obj
make intel_win_32

Rem Copy FDS to bundle directory

copy fds5_win_32.exe %bundledir%\FDS\fds5.exe

Rem Copy Documentation to bundle directory
echo Copying Documentation
copy %pdfdir%\FDS_5_User_Guide.pdf               %bundledir%\Documentation\.
copy %pdfdir%\FDS_5_Technical_Reference_Guide.pdf %bundledir%\Documentation\.
copy %pdfdir%\SMV_5_User_Guide.pdf               %bundledir%\Documentation\.

Rem Copy Smokeview files to bundle directory
echo Copying Smokeview files
copy %smvbundle%\smokeview_release.exe %bundledir%\Smokeview\smokeview.exe
copy %smvbundle%\smokezip_release.exe %bundledir%\Smokeview\smokezip.exe
copy %smvbundle%\devices.svo %bundledir%\Smokeview\.
copy %smvbundle%\glew32.dll %bundledir%\Smokeview\.
copy %smvbundle%\pthreadVC.dll %bundledir%\Smokeview\.
copy %smvbundle%\readme.html %bundledir%\Smokeview\.
copy %smvbundle%\smokeview.ini %bundledir%\Smokeview\.

Rem Copy examples to bundle directory

echo Copying Examples
cd %REPOS%\Test_cases

copy Atmospheric_Effects\*.fds %exampledir%\.
copy Controls\*.fds %exampledir%\.
copy Detectors\*.fds %exampledir%\.
copy Fires\*.fds %exampledir%\.
copy Flowfields\*.fds %exampledir%\.
copy Heat_Transfer\*.fds %exampledir%\.
copy Miscellaneous\*.fds %exampledir%\.
copy Pressure_Effects\*.fds %exampledir%\.
copy Pyrolysis\*.fds %exampledir%\.
copy Sprinklers_and_Sprays\*.fds %exampledir%\.
copy Timing_Benchmarks\*.fds %exampledir%\.
copy Visualization\*.fds %exampledir%\.



pause