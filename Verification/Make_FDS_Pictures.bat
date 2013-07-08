@echo off

echo.
echo erasing User Guide scripted figures generated previously
erase ..\Manuals\FDS_User_Guide\SCRIPT_FIGURES\*.png

echo.
echo erasing Verification scripted figures generated previously
erase ..\Manuals\FDS_Verification_Guide\SCRIPT_FIGURES\*.png

set BASEDIR=%CD%\
set RUNSMV=call "%BASEDIR%\scripts\runsmv.bat"
set SMOKEVIEW=smokeview

..\Utilities\Data_Processing\sh2bat FDS_Pictures.sh FDS_Pictures.bat

call FDS_Pictures.bat

cd %BASEDIR%\Verification
erase FDS_Pictures.bat