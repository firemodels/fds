@echo off

echo.
echo Uninstall FDS and Smokeview
echo.
echo Press any key to proceed or CTRL C to abort
pause>NUL

echo.
echo Removing the association between .smv and Smokeview

assoc .smv=
ftype smvDoc=

echo. 
echo Removing FDS from the Start menu.
rmdir /q /s "%ALLUSERSPROFILE%\Start Menu\Programs\FDS6"


