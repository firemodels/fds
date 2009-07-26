@echo off

IF not EXIST placeholder.txt goto dircheck
echo ***error: wrapup script running in wrong directory.
pause
exit
:dircheck


echo .
echo associating the smt extension with smokeviewt.exe

ftype smtDoc="%CD%\bin\smokeviewt.exe" "%%1" >Nul
assoc .smt=smtDoc>Nul

echo. 
echo Setting shortcuts to FDS/Smokeview documentation and web sites.
rmdir /q /s "%USERPROFILE%\Start Menu\Programs\FDS5"
mkdir "%USERPROFILE%\Start Menu\Programs\FDS5"
"%CD%\shortcut.exe" /F:"%USERPROFILE%\Start Menu\Programs\FDS5\FDS Web Site.lnk"  /T:"%CD%\Documentation\FDS Web Site.url" /A:C
"%CD%\shortcut.exe" /F:"%USERPROFILE%\Start Menu\Programs\FDS5\FDS Development Web Site.lnk"  /T:"%CD%\Documentation\FDS Development Web Site.url" /A:C
"%CD%\shortcut.exe" /F:"%USERPROFILE%\Start Menu\Programs\FDS5\Documentation.lnk"  /T:"%CD%"\Documentation /A:C

echo.
echo Setting Path for FDS

call "%CD%"\set_fds5_path.exe

erase "%CD%"\set_fds5_path.exe
erase "%CD%"\shortcut.exe
erase "%CD%"\wrapup.bat

echo Installation Complete.  
echo Press any key to continue . . . 

pause>NUL

