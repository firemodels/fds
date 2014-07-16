@echo off

IF not EXIST placeholder.txt goto dircheck
echo ***error: This script is running in the wrong directory.
pause
exit
:dircheck

echo.
echo *** Wrapping up the FDS and Smokeview installation.
echo.
echo *** Removing previous FDS/Smokeview entries from the system and user path.
call "%CD%\set_path.exe" -s -m -b -r "nist\fds"
call "%CD%\set_path.exe" -u -m -b -r "FDS\FDS5"
call "%CD%\set_path.exe" -s -m -b -r "FDS\FDS5"
call "%CD%\set_path.exe" -u -m -b -r "FDS\FDS6"
call "%CD%\set_path.exe" -s -m -b -r "FDS\FDS6"

set SAVECD="%CD%"

cd "%CD%\.."
set SHORTCUTSDIR=%CD%\shortcuts

cd %SAVECD%

:: create shortcuts directory

echo.
if exist "%SHORTCUTSDIR%" goto existbin
echo.
echo Creating the directory %SHORTCUTSDIR%
mkdir "%SHORTCUTSDIR%"
:existbin

echo *** Adding program shortcuts to %SHORTCUTSDIR%
:: ------------ create aliases ----------------

set tempfile="%TEMP%\fds_tempfile"
set numcoresfile="%TEMP%\numcoresfile"

:: *** fds5 (32 bit)

set fds5exe="c:\Program Files\FDS\FDS5\bin\fds5.exe"
set fds5bat="%SHORTCUTSDIR%\fds5.bat"

if exist %fds5exe% (
  echo @echo off > %tempfile%
  echo %fds5exe% %%* >> %tempfile%
  copy %tempfile% "%SHORTCUTSDIR%\fds5.bat" > Nul
)

:: *** fds5 (64 bit)

set fds5exe="c:\Program Files\FDS\FDS5\bin\fds5_win_64.exe"
if exist %fds5exe% (
  echo @echo off > %tempfile%
  echo %fds5exe% %%* >> %tempfile%
  copy %tempfile% "%SHORTCUTSDIR%\fds5.bat" > Nul
)

:: *** smokeview5

set smv5exe="c:\Program Files\FDS\FDS5\bin\smokeview.exe"
set smv5bat="%SHORTCUTSDIR%\smokeview5.bat"
if exist %smv5exe% (
  echo @echo off > %tempfile%
  echo %smv5exe% %%* >> %tempfile%
  copy %tempfile% "%SHORTCUTSDIR%\smokeview5.bat" > Nul
)

:: *** fds6

echo @echo off > %tempfile%
echo "%CD%\bin\fds" %%* >> %tempfile%
copy %tempfile% "%SHORTCUTSDIR%\fds6.bat" > Nul

:: *** smokeview6

echo @echo off > %tempfile%
echo "%CD%\bin\smokeview" %%* >> %tempfile%
copy %tempfile% "%SHORTCUTSDIR%\smokeview6.bat" > Nul

:: *** smokediff6

echo @echo off > %tempfile%
echo "%CD%\bin\smokediff" %%* >> %tempfile%
copy %tempfile% %smd6% "%SHORTCUTSDIR%\smokediff6.bat" > Nul

:: *** smokezip6

echo @echo off > %tempfile%
echo "%CD%\bin\smokezip" %%* >> %tempfile%
copy %tempfile% %smd6% "%SHORTCUTSDIR%\smokezip6.bat" > Nul

:: ------------ setting up path ------------

echo.
echo *** Setting up the PATH variable

:: *** c:\...\FDS\FDS6\bin
call "%CD%\set_path.exe" -s -m -f "%CD%\bin"

:: *** c:\...\FDS\shortcuts
call "%CD%\set_path.exe" -s -m -f "%SHORTCUTSDIR%"

:: ------------- file association -------------
echo.
echo *** Associating the .smv file extension with smokeview.exe

ftype smvDoc="%CD%\bin\smokeview.exe" "%%1" >Nul
assoc .smv=smvDoc>Nul

set FDSSTART=%ALLUSERSPROFILE%\Start Menu\Programs\FDS6

:: ------------- start menu shortcuts ---------------
echo. 
echo *** Adding document shortcuts to the Start menu.
if exist "%FDSSTART%" rmdir /q /s "%FDSSTART%"

mkdir "%FDSSTART%"

mkdir "%FDSSTART%\FDS on the Web"
copy "%CD%\Documentation\FDS_on_the_Web\Software_Updates.url"            "%FDSSTART%\FDS on the Web\Software Updates.url" > Nul
copy "%CD%\Documentation\FDS_on_the_Web\Documentation_Updates.url"       "%FDSSTART%\FDS on the Web\Documentation Updates.url" > Nul
copy "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"   "%FDSSTART%\FDS on the Web\Discussion Group.url" > Nul
copy "%CD%\Documentation\FDS_on_the_Web\Official_Web_Site.url"  "%FDSSTART%\FDS on the Web\Official Web Site.url" > Nul
copy "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"   "%FDSSTART%\FDS on the Web\Discussion Group.url" > Nul
copy "%CD%\Documentation\FDS_on_the_Web\Issue_Tracker.url"      "%FDSSTART%\FDS on the Web\Issue Tracker.url" > Nul

mkdir "%FDSSTART%\Guides and Release Notes"
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Configuration Management Plan.lnk"   /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Configuration_Management_Plan.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS User Guide.lnk"                      /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_User_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Technical Reference Guide.lnk"       /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Technical_Reference_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Validation Guide.lnk"                /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Valication_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Verification Guide.lnk"              /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Verification_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Release Notes.lnk"                   /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Release_Notes.htm" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview User Guide.lnk"                /T:"%CD%\Documentation\Guides_and_Release_Notes\SMV_User_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview Technical Reference Guide.lnk" /T:"%CD%\Documentation\Guides_and_Release_Notes\SMV_Technical_Reference_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview Verification Guide.lnk"        /T:"%CD%\Documentation\Guides_and_Release_Notes\SMV_Verification_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview release notes.lnk"             /T:"%CD%\Documentation\Guides_and_Release_Notes\Smokeview_release_notes.html" /A:C >NUL

"%CD%\shortcut.exe" /F:"%FDSSTART%\Overview.lnk"  /T:"%CD%\Documentation\Overview.html" /A:C >NUL
:: "%CD%\shortcut.exe" /F:"%FDSSTART%\Uninstall.lnk"  /T:"%CD%\Uninstall\uninstall.bat" /A:C >NUL

erase "%CD%"\set_path.exe
erase "%CD%"\shortcut.exe

:: ----------- setting up openmp threads environment variable

WMIC CPU Get NumberofLogicalProcessors | more +1 > %numcoresfile%
set /p ncores=<%numcoresfile%

if %ncores% GEQ 8 (
  set nthreads=4
) else (
  if %ncores% GEQ 4 (
    set nthreads=2
  ) else (
    set nthreads=1 
  )
)
setx OMP_NUM_THREADS %nthreads%

:: ----------- setting up firewall for mpi version of FDS

::new
:: set firewall_setup="%CD%\setup_fds_firewall.bat"
:: if exist "%firewall_setup%" (
::    echo setting up firewall exceptions
::    call "%firewall_setup%"
:: )

:: ----------- setting up uninstall file

echo echo. >> Uninstall\Uninstall.bat
echo echo Removing directories, %CD%\bin and %SHORTCUTSDIR%, from the System Path >> Uninstall\Uninstall.bat
echo call "%CD%\Uninstall\set_path.exe" -s -b -r "%CD%\bin" >> Uninstall\Uninstall.bat
echo call "%CD%\Uninstall\set_path.exe" -s -b -r "%SHORTCUTSDIR%" >> Uninstall\Uninstall.bat

echo echo. >> Uninstall\Uninstall.bat
echo echo Delete the directory %CD% by hand (as administrator) to complete the removal of FDS and Smokeview >> Uninstall\Uninstall.bat
echo pause >> Uninstall\Uninstall.bat

echo.
echo *** Press any key to complete the installation.
pause>NUL

:: erase "%CD%"\wrapup_fds_install.bat

