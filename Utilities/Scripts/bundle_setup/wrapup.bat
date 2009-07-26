@echo off
set bindir=%CD%\bin

echo .
echo current directory is %CD%

echo .
echo associating the smt extension with smokeviewt.exe

ftype smtDoc="%bindir%\smokeviewt.exe" "%%1" >Nul
assoc .smt=smtDoc>Nul

echo.
echo Setting Path for FDS

"%CD%"\set_fds5_path.exe
echo.
erase "%CD%"\set_fds5_path.exe
erase "%CD%"\wrapup.bat
echo Installation Complete.  
echo Press any key to continue . . . 
pause>NUL
