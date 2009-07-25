@echo off
set bindir="%CD%"\bin

echo .
echo current directory is %CD%

Rem assoc .smv=smokeview.Document>Nul
Rem ftype smokeview.Document=%bindir%\smokeview.exe %1 >Nul
echo.
echo Setting Path for FDS
"%CD%"\set_fds5_path.exe
echo.
erase "%CD%"\set_fds5_path.exe
erase "%CD%"\wrapup.bat
echo Installation Complete.  
echo Press any key to continue . . . 
pause>NUL
