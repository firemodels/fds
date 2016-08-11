@echo off

:: setup environment variables (defining where repository resides etc) 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%

echo.
echo *** FDS/Smokeview version and revision numbers

echo.
echo smv_version=%smv_version%
echo smv_revision=%smv_revision%
echo fds_version=%fds_version%
echo fds_revision=%fds_revision%

echo.
echo Press any key to continue
pause>NUL

echo.
echo *** FDS-Smokeview repository settings

echo.
echo svn_root=%svn_root%
echo svn_drive=%svn_drive%
echo linux_svn_root=%linux_svn_root%

echo.
echo Press any key to continue
pause>NUL

echo.
echo *** Linux User/Host names

echo linux_hostname=%linux_hostname%
echo linux_username=%linux_username%

echo.
echo *** OSX User/Host names

echo osx_hostname=%osx_hostname%
echo osx_username=%osx_username%


echo.
echo Press any key to continue
pause>NUL

echo.
echo variable display complete.
echo Press any key to continue.
pause>NUL

