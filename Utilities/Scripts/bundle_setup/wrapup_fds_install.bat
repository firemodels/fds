@echo off

IF not EXIST placeholder.txt goto dircheck
echo ***error: This script is running in wrong directory.
pause
exit
:dircheck

echo.
echo FDS and Smokeview installation
echo.
echo ***Note*** This installation removes entries in the system path for
echo pre 5.4 versions of FDS and Smokeview.  Only path entries will be
echo changed, program and data files WILL NOT be removed.
echo.
echo Press any key to proceed or CTRL C to abort
pause>NUL

echo.  
echo Proceeding...

echo.
echo Associating the smv file extension with smokeview.exe

ftype smvDoc="%CD%\bin\smokeview.exe" "%%1" >Nul
assoc .smv=smvDoc>Nul

echo. 
echo Adding FDS and Smokeview shortcuts to the Start menu.
if exist "%USERPROFILE%\Start Menu\Programs\FDS5" rmdir /q /s "%USERPROFILE%\Start Menu\Programs\FDS5"

mkdir "%USERPROFILE%\Start Menu\Programs\FDS5"

mkdir "%USERPROFILE%\Start Menu\Programs\FDS5\FDS_on_the_Web"
copy "%CD%\Documentation\FDS_on_the_Web\D*"          "%USERPROFILE%\Start Menu\Programs\FDS5\FDS_on_the_Web"
copy "%CD%\Documentation\FDS_on_the_Web\O*"          "%USERPROFILE%\Start Menu\Programs\FDS5\FDS_on_the_Web"
copy "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"        %USERPROFILE%\Start Menu\Programs\FDS5\FDS_on_the_Web\Discussion_Group.url"
copy "%CD%\Documentation\FDS_on_the_Web\Issue_Tracker.url"          "%USERPROFILE%\Start Menu\Programs\FDS5\FDS_on_the_Web\Issue_Tracker.url"
copy "%CD%\Documentation\FDS_on_the_Web\Updates.url" "%USERPROFILE%\Start Menu\Programs\FDS5\FDS_on_the_Web\Updates.url"  

mkdir "%USERPROFILE%\Start Menu\Programs\FDS5\User_Guides_and_Release_Notes"
"%CD%\shortcut.exe" /F:"%USERPROFILE%\Start Menu\Programs\FDS5\User_Guides_and_Release_Notes\FDS_5_User_Guide.lnk"  /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_5_User_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%USERPROFILE%\Start Menu\Programs\FDS5\User_Guides_and_Release_Notes\FDS_Release_Notes.lnk"  /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Release_Notes.htm" /A:C >NUL
"%CD%\shortcut.exe" /F:"%USERPROFILE%\Start Menu\Programs\FDS5\User_Guides_and_Release_Notes\SMV_5_User_Guide.lnk"  /T:"%CD%\Documentation\Guides_and_Release_Notes\SMV_5_User_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%USERPROFILE%\Start Menu\Programs\FDS5\User_Guides_and_Release_Notes\Smokeview_release_notes.lnk"  /T:"%CD%\Documentation\Guides_and_Release_Notes\Smokeview_release_notes.html" /A:C >NUL

"%CD%\shortcut.exe" /F:"%USERPROFILE%\Start Menu\Programs\FDS5\Overview.lnk"  /T:"%CD%\Documentation\Overview.html" /A:C >NUL

echo.
echo Adding the directory %CD% to the user path variable for: %USERNAME%

call "%CD%"\set_path.exe -a "%CD%\bin"

echo.
echo Looking for pre 5.4 FDS/Smokeview path entries in the system path

call "%CD%"\set_path.exe -r

erase "%CD%"\set_path.exe
erase "%CD%"\shortcut.exe

echo.
echo Press any key to complete the installation.
pause>NUL
erase "%CD%"\wrapup_fds_install.bat

