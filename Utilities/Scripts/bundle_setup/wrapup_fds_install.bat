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
echo Removing previous FDS/Smokeview entries from the system and user path.
call "%CD%\set_path.exe" -s -m -b -r "nist\fds"
call "%CD%\set_path.exe" -u -m -b -r "FDS\FDS5"

Rem create c:\bin directory

if exist c:\bin goto existbin
echo.
echo Creating the directory c:\bin
mkdir c:\bin
:existbin

Rem create 32 bit fds5 alias

set fds5exe="c:\Program Files\FDS\FDS5\bin\fds5.exe"
set fds5bat="c:\bin\fds5.bat"

if exist %fds5exe% (
  echo Adding fds5.bat to c:\bin
  echo @echo off > %fds5bat%
  echo %fds5exe% %%* >> %fds5bat%
)

Rem create 64 bit fds5 alias

set fds5exe="c:\Program Files\FDS\FDS5\bin\fds5_win_64.exe"

if exist %fds5exe% (
  echo Adding fds5.bat to c:\bin
  echo @echo off > %fds5bat%
  echo %fds5exe% %%* >> %fds5bat%
)


Rem create smokeview5 alias

set smv5exe="c:\Program Files\FDS\FDS5\bin\smokeview.exe"
set smv5bat="c:\bin\smokeview5.bat"

if exist %smv5exe% (
  echo Adding smokeview5.bat to c:\bin
  echo @echo off > %smv5bat%
  echo %smv5exe% %%* >> %smv5bat%
)

Rem create fds6 and smokeview6 aliases

echo.
echo Adding fds6.bat and smokeview6.bat to c:\bin
set fds6=c:\bin\fds6.bat
echo @echo off > %fds6%
echo "%CD%\bin\fds" %%* >> %fds6%

set smv6=c:\bin\smokeview6.bat
echo @echo off > %smv6%
echo "%CD%\bin\smokeview" %%* >> %smv6%

echo.
echo Adding c:\bin to the system path 
call "%CD%\set_path.exe" -s -m -a "c:\bin"

echo.
echo Adding %CD%\bin to the system path 
call "%CD%\set_path.exe" -s -m -a "%CD%\bin"

echo.
echo Associating the .smv file extension with smokeview.exe

ftype smvDoc="%CD%\bin\smokeview.exe" "%%1" >Nul
assoc .smv=smvDoc>Nul

set FDSSTART=%ALLUSERSPROFILE%\Start Menu\Programs\FDS6

echo. 
echo Adding shortcuts to the Start menu.
if exist "%FDSSTART%" rmdir /q /s "%FDSSTART%"

mkdir "%FDSSTART%"

mkdir "%FDSSTART%\FDS on the Web"
copy "%CD%\Documentation\FDS_on_the_Web\Software_Updates.url"            "%FDSSTART%\FDS on the Web\Software Updates.url"
copy "%CD%\Documentation\FDS_on_the_Web\Documentation_Updates.url"       "%FDSSTART%\FDS on the Web\Documentation Updates.url"
copy "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"   "%FDSSTART%\FDS on the Web\Discussion Group.url"
copy "%CD%\Documentation\FDS_on_the_Web\Official_Web_Site.url"  "%FDSSTART%\FDS on the Web\Official Web Site.url"
copy "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"   "%FDSSTART%\FDS on the Web\Discussion Group.url"
copy "%CD%\Documentation\FDS_on_the_Web\Issue_Tracker.url"      "%FDSSTART%\FDS on the Web\Issue Tracker.url"

mkdir "%FDSSTART%\Guides and Release Notes"
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS User Guide.lnk"        /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_User_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Technical Reference Guide.lnk"    /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Technical_Reference_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Release Notes.lnk"       /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Release_Notes.htm" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview User Guide.lnk"        /T:"%CD%\Documentation\Guides_and_Release_Notes\SMV_User_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview Technical Reference Guide.lnk"    /T:"%CD%\Documentation\Guides_and_Release_Notes\SMV_Technical_Reference_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview release notes.lnk" /T:"%CD%\Documentation\Guides_and_Release_Notes\Smokeview_release_notes.html" /A:C >NUL

"%CD%\shortcut.exe" /F:"%FDSSTART%\Overview.lnk"  /T:"%CD%\Documentation\Overview.html" /A:C >NUL
Rem "%CD%\shortcut.exe" /F:"%FDSSTART%\Uninstall.lnk"  /T:"%CD%\Uninstall\uninstall.bat" /A:C >NUL

erase "%CD%"\set_path.exe
erase "%CD%"\shortcut.exe

echo echo. >> Uninstall\Uninstall.bat
echo echo Removing directory, %CD%\bin, from the System Path >> Uninstall\Uninstall.bat
echo call "%CD%\Uninstall\set_path.exe" -s -b -r "%CD%\bin" >> Uninstall\Uninstall.bat

echo echo. >> Uninstall\Uninstall.bat
echo echo Delete the directory %CD% by hand (as administrator) to complete the removal of FDS and Smokeview >> Uninstall\Uninstall.bat
echo pause >> Uninstall\Uninstall.bat

echo.
echo Press any key to complete Installation.
pause>NUL

erase "%CD%"\wrapup_fds_install.bat

