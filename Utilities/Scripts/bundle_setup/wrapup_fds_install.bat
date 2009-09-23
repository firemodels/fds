@echo off

IF not EXIST placeholder.txt goto dircheck
echo ***error: This script is running in the wrong directory.
pause
exit
:dircheck

echo.
echo FDS and Smokeview installation.
echo.
echo Removing pre 5.4 FDS/Smokeview entries (if present) from the System Path.
echo Press any key to proceed or CTRL C to abort
pause>NUL

call "%CD%\set_path.exe" -s -m -b -r "nist\fds"

echo.
echo Associating the smv file extension with smokeview.exe

ftype smvDoc="%CD%\bin\smokeview.exe" "%%1" >Nul
assoc .smv=smvDoc>Nul

set FDS5START=%USERPROFILE%\Start Menu\Programs\FDS5

echo. 
echo Adding FDS and Smokeview shortcuts to the Start menu.
if exist "%FDS5START%" rmdir /q /s "%FDS5START%"

mkdir "%FDS5START%"

mkdir "%FDS5START%\FDS on the Web"
copy "%CD%\Documentation\FDS_on_the_Web\Developer_Web_Site.url" "%FDS5START%\FDS on the Web\Developer Web Site.url"
copy "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"   "%FDS5START%\FDS on the Web\Discussion Group.url"
copy "%CD%\Documentation\FDS_on_the_Web\Official_Web_Site.url"  "%FDS5START%\FDS on the Web\Official Web Site.url"
copy "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"   "%FDS5START%\FDS on the Web\Discussion Group.url"
copy "%CD%\Documentation\FDS_on_the_Web\Issue_Tracker.url"      "%FDS5START%\FDS on the Web\Issue Tracker.url"
copy "%CD%\Documentation\FDS_on_the_Web\Updates.url"            "%FDS5START%\FDS on the Web\Updates.url"  

mkdir "%FDS5START%\Guides and Release Notes"
"%CD%\shortcut.exe" /F:"%FDS5START%\Guides and Release Notes\FDS 5 User Guide.lnk"        /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_5_User_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDS5START%\Guides and Release Notes\FDS Release Notes.lnk"       /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Release_Notes.htm" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDS5START%\Guides and Release Notes\SMV 5 User Guide.lnk"        /T:"%CD%\Documentation\Guides_and_Release_Notes\SMV_5_User_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDS5START%\Guides and Release Notes\Smokeview release notes.lnk" /T:"%CD%\Documentation\Guides_and_Release_Notes\Smokeview_release_notes.html" /A:C >NUL
copy "%CD%\Documentation\Guides_and_Release_Notes\Latest_Documentation.url"            "%FDS5START%\Guides and Release Notes\Latest Documentation.url"  

"%CD%\shortcut.exe" /F:"%FDS5START%\Overview.lnk"  /T:"%CD%\Documentation\Overview.html" /A:C >NUL

echo.
echo Adding %CD%\bin to the User Path (if absent) for: %USERNAME%

call "%CD%\set_path.exe" -m -a "%CD%\bin"

erase "%CD%"\set_path.exe
erase "%CD%"\shortcut.exe

echo.
echo Installation complete.  Press any key to continue.
pause>NUL
erase "%CD%"\wrapup_fds_install.bat

