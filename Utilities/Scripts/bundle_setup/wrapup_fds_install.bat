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
rmdir /q /s "%USERPROFILE%\Start Menu\Programs\FDS5"
mkdir "%USERPROFILE%\Start Menu\Programs\FDS5"
"%CD%\shortcut.exe" /F:"%USERPROFILE%\Start Menu\Programs\FDS5\Documentation and Updates.lnk"  /T:"%CD%"\Documentation /A:C >NUL
"%CD%\shortcut.exe" /F:"%USERPROFILE%\Start Menu\Programs\FDS5\FDS (version only).lnk"     /T:"%CD%"\bin\fds5.exe /A:C >NUL
"%CD%\shortcut.exe" /F:"%USERPROFILE%\Start Menu\Programs\FDS5\Smokeview.lnk"      /T:"%CD%"\bin\smokeview.exe /A:C >NUL

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

