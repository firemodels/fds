@echo off

Rem  Windows batch file to build directory tree of FDS figures.

set version_tag=2010_0618_r6356

cd to_google

set indir=..\..\..\Manuals
set outdir=fdsfigures_%version_tag%
set FDSUG=FDS_User_Guide\SCRIPT_FIGURES
set FDSVG=FDS_Verification_Guide\SCRIPT_FIGURES
set FDSTG=FDS_Technical_Reference_Guide\SCRIPT_FIGURES

echo.
echo filling FDS figures directory
IF EXIST %outdir% rmdir /S /Q %outdir%
mkdir %outdir%

mkdir %outdir%\%FDSTG%
mkdir %outdir%\%FDSUG%
mkdir %outdir%\%FDSVG%

copy %indir%\%FDSUG%\*.pdf  %outdir%\%FDSUG%
copy %indir%\%FDSUG%\*.png  %outdir%\%FDSUG%

copy %indir%\%FDSVG%\*.pdf  %outdir%\%FDSVG%
copy %indir%\%FDSVG%\*.png  %outdir%\%FDSVG%
copy %indir%\%FDSVG%\*.tex  %outdir%\%FDSVG%

echo figure bundle build finished

echo
echo winzipping figures directory
cd %outdir%
wzzip -a -r -P ..\%outdir%.zip *
cd ..
echo
echo creating self-extracting archive
wzipse32 %outdir%.zip -d "d:\fds-smv\manuals"

cd ..

pause
