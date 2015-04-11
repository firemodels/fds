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
call "%CD%\set_path.exe" -s -m -b -r "nist\fds" >Nul
call "%CD%\set_path.exe" -u -m -b -r "FDS\FDS5" >Nul
call "%CD%\set_path.exe" -s -m -b -r "FDS\FDS5" >Nul
call "%CD%\set_path.exe" -u -m -b -r "FDS\FDS6" >Nul
call "%CD%\set_path.exe" -s -m -b -r "FDS\FDS6" >Nul

set SAVECD="%CD%"

cd "%CD%\.."
set SHORTCUTSDIR=%CD%\shortcuts

cd %SAVECD%

:: create shortcuts directory

echo.
if exist "%SHORTCUTSDIR%" goto existbin
echo.
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

:: *** fds_mpi

echo @echo off > %tempfile%
echo "%CD%\bin\fds" %%* >> %tempfile%
copy %tempfile% "%SHORTCUTSDIR%\fds_mpi.bat" > Nul

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
echo *** Setting up the PATH variable.

:: *** c:\...\FDS\FDS6\bin
call "%CD%\set_path.exe" -s -m -f "%CD%\bin" >Nul

:: *** c:\...\FDS\shortcuts
call "%CD%\set_path.exe" -s -m -f "%SHORTCUTSDIR%" > Nul

:: ------------- file association -------------
echo.
echo *** Associating the .smv file extension with smokeview.exe

ftype smvDoc="%CD%\bin\smokeview.exe" "%%1" >Nul
assoc .smv=smvDoc>Nul

:: ------------- remove old executables -------------
echo.
echo *** Removing old executables if they exist
if exist fds6.exe erase fds6.exe
if exist fds_mpi.exe erase fds_mpi.exe

set FDSSTART=%ALLUSERSPROFILE%\Start Menu\Programs\FDS6

:: ------------- start menu shortcuts ---------------
echo. 
echo *** Adding document shortcuts to the Start menu.
if exist "%FDSSTART%" rmdir /q /s "%FDSSTART%"

mkdir "%FDSSTART%"

mkdir "%FDSSTART%\FDS on the Web"
copy "%CD%\Documentation\FDS_on_the_Web\Software_Updates.url"            "%FDSSTART%\FDS on the Web\Software Updates.url" > Nul
copy "%CD%\Documentation\FDS_on_the_Web\Documentation_Updates.url"       "%FDSSTART%\FDS on the Web\Documentation Updates.url" > Nul
copy "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"            "%FDSSTART%\FDS on the Web\Discussion Group.url" > Nul
copy "%CD%\Documentation\FDS_on_the_Web\Official_Web_Site.url"           "%FDSSTART%\FDS on the Web\Official Web Site.url" > Nul
copy "%CD%\Documentation\FDS_on_the_Web\Discussion_Group.url"            "%FDSSTART%\FDS on the Web\Discussion Group.url" > Nul
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
"%CD%\shortcut.exe" /F:"%FDSSTART%\Uninstall.lnk"  /T:"%CD%\Uninstall\uninstall.bat" /A:C >NUL

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
setx -m OMP_NUM_THREADS %nthreads% > Nul

:: ----------- setting up firewall for mpi version of FDS

:: remove smpd and hydra

smpd -remove 1>> Nul 2>&1
hydra_service -remove 1>> Nul 2>&1

copy "%CD%\bin"\hydra_service2.exe "%CD%\bin"\hydra_service.exe>Nul
erase "%CD%\bin"\hydra_service2.exe >Nul

set firewall_setup="%CD%\setup_fds_firewall.bat"
echo.
echo *** Setting up firewall exceptions.
call %firewall_setup% "%CD%\bin"

:: ----------- copy backup files to Uninstall directory

copy *.txt Uninstall > Nul
erase /q *.txt

:: ----------- setting up uninstall file

echo.
echo *** Setting up Uninstall script.
echo echo. >> Uninstall\uninstall_base.bat
echo echo Removing directories, %CD%\bin and %SHORTCUTSDIR%, from the System Path >> Uninstall\uninstall_base.bat
echo call "%CD%\Uninstall\set_path.exe" -s -b -r "%CD%\bin" >> Uninstall\uninstall_base.bat
echo call "%CD%\Uninstall\set_path.exe" -s -b -r "%SHORTCUTSDIR%" >> Uninstall\uninstall_base.bat

echo echo. >> Uninstall\uninstall_base.bat
echo echo Removing %CD% >> Uninstall\uninstall_base.bat
echo rmdir /s /q "%SHORTCUTSDIR%" >> Uninstall\uninstall_base.bat
echo rmdir /s /q "%CD%" >> Uninstall\Uninstall_base.bat
echo echo *** Uninstall complete >> Uninstall\uninstall_base.bat
echo pause>Nul >> Uninstall\uninstall_base.bat

echo "%CD%\Uninstall\uninstall.vbs" >> Uninstall\uninstall.bat
echo echo Uninstall complete >> Uninstall\uninstall.bat
echo pause >> Uninstall\uninstall.bat

set ELEVATE_APP=%CD%\Uninstall\Uninstall_base.bat
set ELEVATE_PARMS=
echo Set objShell = CreateObject("Shell.Application") >>Uninstall\uninstall.vbs
echo Set objWshShell = WScript.CreateObject("WScript.Shell") >>Uninstall\uninstall.vbs
echo Set objWshProcessEnv = objWshShell.Environment("PROCESS") >>Uninstall\uninstall.vbs
echo objShell.ShellExecute "%ELEVATE_APP%", "%ELEVATE_PARMS%", "", "runas" >>Uninstall\uninstall.vbs

echo.
echo *** Press any key to complete the installation.
pause>NUL

erase "%CD%"\setup_fds_firewall.bat >Nul
erase "%CD%"\wrapup_fds_install.bat >Nul

