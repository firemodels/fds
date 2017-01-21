@echo off

echo.
echo *** Wrapping up the FDS and Smokeview installation.
echo.
echo *** Testing whether Windows is 64 bit
echo.

if defined PROGRAMFILES(X86) (
  echo *** 64 bit Windows detected
  echo.
  echo *** Wrapup proceeding
  echo.
) else (
  echo *** Fatal error: 32 bit Windows detected
  echo     FDS and Smokeview only run on 64 bit systems
  echo     Installation wrap up aborted
  echo *** Press any key to continue.    ***
pause>NUL
  goto abort
)

echo *** Removing previous FDS/Smokeview entries from the system and user path.
call "%CD%\set_path.exe" -s -m -b -r "nist\fds" >Nul
call "%CD%\set_path.exe" -u -m -b -r "FDS\FDS5" >Nul
call "%CD%\set_path.exe" -s -m -b -r "FDS\FDS5" >Nul
call "%CD%\set_path.exe" -u -m -b -r "FDS\FDS6" >Nul
call "%CD%\set_path.exe" -s -m -b -r "FDS\FDS6" >Nul
call "%CD%\set_path.exe" -s -m -b -r "firemodels\FDS6" >Nul
call "%CD%\set_path.exe" -s -m -b -r "firemodels\SMV6" >Nul

set SAVECD="%CD%"

cd "%CD%\.."
set SMV6=%CD%\SMV6

cd %SAVECD%

:: ------------ create aliases ----------------

set numcoresfile="%TEMP%\numcoresfile"

:: ------------ setting up path ------------

echo.
echo *** Setting up the PATH variable.

:: *** c:\...\FDS\FDS6\bin
call "%CD%\set_path.exe" -s -m -f "%CD%\bin" >Nul

:: *** c:\...\FDS\SMV6
call "%CD%\set_path.exe" -s -m -f "%SMV6%" > Nul

:: ------------- file association -------------
echo.
echo *** Associating the .smv file extension with smokeview.exe

ftype smvDoc="%SMV6%\smokeview.exe" "%%1" >Nul
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

copy "%CD%\Documentation\FDS_on_the_Web\Official_Web_Site.url"     "%FDSSTART%\FDS Home Page.url"         > Nul

mkdir "%FDSSTART%\Guides and Release Notes"
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Config Management Plan.lnk"          /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Config_Management_Plan.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS User Guide.lnk"                      /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_User_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Technical Reference Guide.lnk"       /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Technical_Reference_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Validation Guide.lnk"                /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Validation_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Verification Guide.lnk"              /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Verification_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\FDS Release Notes.lnk"                   /T:"%CD%\Documentation\Guides_and_Release_Notes\FDS_Release_Notes.htm" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview User Guide.lnk"                /T:"%CD%\Documentation\Guides_and_Release_Notes\SMV_User_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview Technical Reference Guide.lnk" /T:"%CD%\Documentation\Guides_and_Release_Notes\SMV_Technical_Reference_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview Verification Guide.lnk"        /T:"%CD%\Documentation\Guides_and_Release_Notes\SMV_Verification_Guide.pdf" /A:C >NUL
"%CD%\shortcut.exe" /F:"%FDSSTART%\Guides and Release Notes\Smokeview release notes.lnk"             /T:"%CD%\Documentation\Guides_and_Release_Notes\Smokeview_release_notes.html" /A:C >NUL

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
echo echo.                                                    >> Uninstall\uninstall_base.bat

:: remove smokeview path and directory
echo if %%cfastinstalled%% == 1 goto skip2                    >> Uninstall\uninstall_base.bat
echo echo Removing directory, %SMV6%, from the System Path    >> Uninstall\uninstall_base.bat
echo call "%CD%\Uninstall\set_path.exe" -s -b -r "%SMV6%"     >> Uninstall\uninstall_base.bat
echo rmdir /s /q "%CD%\..\SMV6"                               >> Uninstall\Uninstall_base.bat
echo :skip2 >> Uninstall\uninstall_base.bat

:: remove FDS path and directory
echo echo Removing directory, %CD%\bin , from the System Path >> Uninstall\uninstall_base.bat
echo call "%CD%\Uninstall\set_path.exe" -s -b -r "%CD%\bin"   >> Uninstall\uninstall_base.bat
echo echo.                                                    >> Uninstall\uninstall_base.bat
echo echo Removing %CD%                                       >> Uninstall\uninstall_base.bat
echo rmdir /s /q "%CD%"                                       >> Uninstall\Uninstall_base.bat
echo pause                                                    >> Uninstall\Uninstall_base.bat

echo echo *** Uninstall complete                              >> Uninstall\uninstall_base.bat
echo pause>Nul                                                >> Uninstall\uninstall_base.bat

type Uninstall\uninstall_base2.bat                            >> Uninstall\uninstall_base.bat
erase Uninstall\uninstall_base2.bat

echo "%CD%\Uninstall\uninstall.vbs"                           >> Uninstall\uninstall.bat
echo echo Uninstall complete                                  >> Uninstall\uninstall.bat
echo pause                                                    >> Uninstall\uninstall.bat

set ELEVATE_APP=%CD%\Uninstall\Uninstall_base.bat
set ELEVATE_PARMS=
echo Set objShell = CreateObject("Shell.Application")                       > Uninstall\uninstall.vbs
echo Set objWshShell = WScript.CreateObject("WScript.Shell")               >> Uninstall\uninstall.vbs
echo Set objWshProcessEnv = objWshShell.Environment("PROCESS")             >> Uninstall\uninstall.vbs
echo objShell.ShellExecute "%ELEVATE_APP%", "%ELEVATE_PARMS%", "", "runas" >> Uninstall\uninstall.vbs

echo.
echo *** Press any key.  Then reboot to complete the installation.  ***
pause>NUL

erase "%CD%"\setup_fds_firewall.bat >Nul
erase wrapup_fds_install.bat >Nul
goto eof

:abort
erase "%CD%"\setup_fds_firewall.bat >Nul
erase wrapup_fds_install.bat >Nul

:eof

