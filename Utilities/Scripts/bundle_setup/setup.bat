@echo off
set fromdir=%CD%\fds5
set todir=%ProgramFiles%\FDS5
set bindir=%todir%\bin

echo .
echo copying files to %todir%
xcopy  "%fromdir%" "%todir%" /E /C /I /Y /Q
erase "%todir%"\setup.bat


Rem assoc .smv=smokeview.Document>Nul
Rem ftype smokeview.Document=%bindir%\smokeview.exe %1 >Nul
echo.
echo Setting Path for FDS
%fromdir%\fds5_setup.exe
echo.
echo Installation Complete.  
echo Press any key to continue . . . 
pause>NUL
