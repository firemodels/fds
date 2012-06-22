@echo off

IF not EXIST placeholder.txt goto dircheck
echo ***error: This script is running in the wrong directory.
pause
exit
:dircheck

echo.
echo Wrapping up FDS and Smokeview installation.
echo.
echo.
echo echo Removing pre 5.4 FDS/Smokeview entries (if present) from the system and user path.

call "%CD%\set_path.exe" -s -m -b -r "nist\fds"

call "%CD%\custom_env.bat"

echo.
echo Associating the smv file extension with smokeview%fdssmv_major_version%_win_64.exe

ftype smvDoc="%CD%\bin\smokeview%fdssmv_major_version%_win_64.exe" "%%1" >Nul
assoc .smv=smvDoc>Nul

set FDSSTART=%ALLUSERSPROFILE%\Start Menu\Programs\%fds_edition%
echo FDSSTART=%FDSSTART%

echo. 
echo Adding FDS and Smokeview shortcuts to the Start menu.
if exist "%FDSSTART%" rmdir /q /s "%FDSSTART%"

mkdir "%FDSSTART%"

mkdir "%FDSSTART%\FDS on the Web"
echo copying "%CD%\Documentation\FDS_on_the_Web\Software_Updates.url"
copy "%CD%\Documentation\FDS_on_the_Web\Software_Updates.url"            "%FDSSTART%\FDS on the Web\Software Updates.url"

echo copying "%CD%\Documentation\FDS_on_the_Web\Documentation_Updates.url"
copy "%CD%\Documentation\FDS_on_the_Web\Documentation_Updates.url"       "%FDSSTART%\FDS on the Web\Documentation Updates.url"

echo copying "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"
copy "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"   "%FDSSTART%\FDS on the Web\Discussion Group.url"

echo copying "%CD%\Documentation\FDS_on_the_Web\Official_Web_Site.url"
copy "%CD%\Documentation\FDS_on_the_Web\Official_Web_Site.url"  "%FDSSTART%\FDS on the Web\Official Web Site.url"

echo copying "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"
copy "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"   "%FDSSTART%\FDS on the Web\Discussion Group.url"

echo copy "%CD%\Documentation\FDS_on_the_Web\Issue_Tracker.url"
copy "%CD%\Documentation\FDS_on_the_Web\Issue_Tracker.url"      "%FDSSTART%\FDS on the Web\Issue Tracker.url"

mkdir "%FDSSTART%\Guides and Release Notes"
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS User Guide.lnk"        /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_User_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Technical Reference Guide.lnk"    /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Technical_Reference_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Release Notes.lnk"       /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Release_Notes.htm" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview User Guide.lnk"        /T:"%CD%\Documentation\Guides_and_Release_Notes\SMV_User_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview Technical Reference Guide.lnk"    /T:"%CD%\Documentation\Guides_and_Release_Notes\SMV_Technical_Reference_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview release notes.lnk" /T:"%CD%\Documentation\Guides_and_Release_Notes\Smokeview_release_notes.html" /A:C >NUL

"%CD%\shortcut.exe" /F:"%FDSSTART%\Overview.lnk"  /T:"%CD%\Documentation\Overview.html" /A:C >NUL

echo.
echo Adding %CD%\bin to the system path 

call "%CD%\set_path.exe" -s -m -a "%CD%\bin"

erase "%CD%"\set_path.exe
erase "%CD%"\shortcut.exe

echo.
echo Press any key to complete Installation.
pause>NUL
Rem erase "%CD%"\wrapup_fds_install.bat
Rem erase "%CD%"\custom_env.bat

