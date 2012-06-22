@echo off

IF not EXIST placeholder.txt goto dircheck
echo ***error: This script is running in wrong directory.
pause
exit
:dircheck

echo.
echo Uninstall FDS and Smokeview
echo.
echo Press any key to proceed or CTRL C to abort
pause>NUL

echo.
echo Removing the association between .smv and Smokeview

assoc .smv=
echo ftype smvDoc=
Rem ftype smvDoc=

echo. 
echo Removing the FDS entry from the Start menu.
echo rmdir /q /s "%ALLUSERSPROFILE%\Start Menu\Programs\FDS6"
Rem rmdir /q /s "%ALLUSERSPROFILE%\Start Menu\Programs\FDS6"


cd ..
echo.
echo Removing %CD%\bin from the System User Path
echo call Uninstall\set_path.exe -s -m -r "%CD%\bin"
Rem call Uninstall\set_path.exe -s -m -r "%CD%\bin"

cd ..
echo.
echo Removing directory "%CD%\FDS6"
Rem rmdir /q /s "%CD%\FDS6"
echo rmdir /q /s "%CD%\FDS6"

echo.
echo Uninstall Complete.  Press any key to continue.
pause>NUL

