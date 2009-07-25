@echo off
set fromdir=%CD%
set todir=%ProgramFiles%\FDS5
set bindir=%todir%\bin

echo .
echo copying files to %todir% from %fromdir%
xcopy  "%fromdir%" "%todir%" /E /C /I /Y /Q
erase "%todir%"\setup.bat
erase "%todir%"\set_fds5_path.bat

Rem assoc .smv=smokeview.Document>Nul
Rem ftype smokeview.Document=%bindir%\smokeview.exe %1 >Nul
echo.
echo Setting Path for FDS
%fromdir%\set_fds5_path.exe
echo.
echo Installation Complete.  
echo Press any key to continue . . . 
pause>NUL
