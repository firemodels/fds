@echo off
set LIB_TAG=v2.32.0

::*** library and tag name are the same

set LIB_DIR=%LIB_TAG%


::*** placehoder for parsing options

set clean_hypre=

call :getopts %*
if %stopscript% == 1 exit /b

::*** make sure cmake and gcc are installed

set abort=0
set buildstatus=
call :is_file_installed cmake || set abort=1
call :is_file_installed gcc   || set abort=1
if %abort% == 1 exit /b

set CURDIR=%CD%

::*** define root directory where fds repo and libs directories are located

set FIREMODELS=..\..\..\..
cd %FIREMODELS%
set FIREMODELS=%CD%
cd %CURDIR%

set INSTALLDIR=%FIREMODELS%\libs\hypre\%LIB_DIR%

::*** erase install directory if clean option was specified

if "x%clean_hypre%" == "x" goto endif1
  if exist %INSTALLDIR% rmdir /s /q %INSTALLDIR%
:endif1

::*** if hypre library directory exists exit and use it

if not exist %INSTALLDIR% goto endif2
  set HYPRE_HOME=%INSTALLDIR%
  set buildstatus=prebuilt
  goto eof
:endif2

::*** if directory pointed to by HYPRE_HOME exists exit and use it

if "x%HYPRE_HOME%" == "x" goto endif3
if not exist %HYPRE_HOME%  goto endif3
if "x%clean_hypre%" == "x" goto endif3
  rmdir /s /q %HYPRE_HOME%
:endif3

if "x%HYPRE_HOME%" == "x" goto endif4
  if not exist %HYPRE_HOME%  goto endif4
    set buildstatus=prebuilt
    goto eof
:endif4

::*** if hypre repo does not exist exit and build fds without it

set LIB_REPO=%FIREMODELS%\hypre
if exist %LIB_REPO% goto endif5
  set HYPRE_HOME=
  set buildstatus=norepo
  goto eof
:endif5

::*** if we've gotten this far the prebuilt libraries do not exist, the repo does exist so build the hypre library

cd %CURDIR%

echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo building Hypre library version %LIB_TAG%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.

set buildstatus=build
echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo setting up Intel compilers
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
call %FIREMODELS%\fds\Build\Scripts\setup_intel_compilers.bat

cd %CURDIR%

echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo checking out tag %LIB_TAG%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
cd %LIB_REPO%
git checkout %LIB_REPO%\src\config\HYPRE_config.h.cmake.in
git checkout %LIB_TAG%

echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo changing HYPRE_FMANGLE 0 to HYPRE_FMANGLE 4
echo in the file HYPRE_config.h.cmake.in
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
powershell -Command "(Get-Content %LIB_REPO%\src\config\HYPRE_config.h.cmake.in) -replace 'HYPRE_FMANGLE 0', 'HYPRE_FMANGLE 4' | Set-Content %LIB_REPO%\src\config\HYPRE_config.h.cmake.in"

cd %CURDIR%

echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo cleaning hypre repo
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.

cd %LIB_REPO%
git clean -dxf

echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo configuring hypre version %LIB_TAG%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.

set BUILDDIR=%LIB_REPO%\src\cmbuild
cd %BUILDDIR%
cmake ..\  ^
-G "MinGW Makefiles" ^
-DCMAKE_INSTALL_PREFIX="%INSTALLDIR%" ^
-DCMAKE_C_COMPILER=icx ^
-DCMAKE_C_FLAGS="/DWIN32 -O3 /fp:precise" ^
-DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded" ^
-DCMAKE_INSTALL_LIBDIR="lib"

echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo building and installing hypre version %LIB_TAG%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
call make install

echo.
set HYPRE_HOME=%INSTALLDIR%
echo Hypre library version %LIB_TAG% built in %INSTALLDIR%
echo.

cd %CURDIR%

goto eof

:: -------------------------------------------------------------
:getopts
:: -------------------------------------------------------------
 set stopscript=0
 if (%1)==() exit /b
 set valid=0
 set arg=%1
 if /I "%1" EQU "--clean-hypre" (
    set clean_hypre=1
    set valid=1
 )
 if /I "%1" EQU "-help" (
   call :usage
   set stopscript=1
   exit /b
 )
 shift
 if %valid% == 0 (
   echo.
   echo ***Error: the input argument %arg% is invalid
   echo.
   echo Usage:
   call :usage
   set stopscript=1
   exit /b
 )
if not (%1)==() goto getopts
exit /b

:: -------------------------------------------------------------
:is_file_installed
:: -------------------------------------------------------------

  set program=%1
  where %program% 1> installed_error.txt 2>&1
  type installed_error.txt | find /i /c "Could not find" > installed_error_count.txt
  set /p nothave=<installed_error_count.txt
  erase installed_error_count.txt installed_error.txt
  if %nothave% == 1 (
    echo "***Fatal error: %program% not present"
    exit /b 1
  )
  exit /b 0

:: -------------------------------------------------------------
:have_program
:: -------------------------------------------------------------
:: same as is_file_installed except does not abort script if program is not insstalled

  set program=%1
  where %program% 1> installed_error.txt 2>&1
  type installed_error.txt | find /i /c "Could not find" > installed_error_count.txt
  set /p nothave=<installed_error_count.txt
  erase installed_error_count.txt installed_error.txt
  if %nothave% == 1 (
    exit /b 1
  )
  exit /b 0

:: -------------------------------------------------------------
:usage  
:: -------------------------------------------------------------
echo build hypre
echo. 
echo -clean-hypre    - force build of hypre library
echo -help           - display this message
exit /b

:eof
echo.
if "%buildstatus%" == "norepo"   echo HYPRE library not built, The hypre git repo does not exist
if "%buildstatus%" == "prebuilt" echo HYPRE library not built. It exists in %HYPRE_HOME%
