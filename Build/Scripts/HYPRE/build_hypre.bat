@echo off
set INSTALLDIR=hypre-2.31.0
set HYPREVERSION=v2.31.0

call :getopts %*
if %stopscript% == 1 exit /b

set have_setx=1
call :have_program setx       || set have_setx=0

set abort=0
call :is_file_installed cmake || set abort=1
call :is_file_installed gcc   || set abort=1
if %abort% == 1 exit /b

set CURDIR=%CD%

set FIREMODELS=..\..\..\..
cd %FIREMODELS%
set FIREMODELS=%CD%
cd %CURDIR%

set LIBS=%FIREMODELS%\libs
if not exist %LIBS% mkdir %LIBS%
if not exist %LIBS% echo failed to create LIBS directory
if not exist %LIBS% exit
cd %LIBS%
set LIBS=%CD%
set INSTALLDIR=%LIBS%\%INSTALLDIR%
cd %CURDIR%

echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo setting up Intel compilers
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
call %FIREMODELS%\fds\Build\Scripts\setup_intel_compilers.bat

cd %CURDIR%

set HYPRE=%FIREMODELS%\hypre

:: clone hypre repo (at same level as fds, smv etc repos) if it doesn't exist
if exist %HYPRE% goto endif1
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo cloning hypre from https://github.com/LLNL/hypre.git
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.

  cd %FIREMODELS%
  git clone https://github.com/hypre-space/hypre.git
  cd hypre
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo checking out version %HYPREVERSION%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
  git checkout %HYPREVERSION%

echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo modify HYPRE_config.h.cmake.in file
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
echo change HYPRE_FMANGLE line to #define HYPRE_FMANGLE 4
echo after saving file, press enter
  notepad %HYPRE%\src\config\HYPRE_config.h.cmake.in

  pause
  cd %CURDIR%
:endif1

echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo cleaning hypre repo
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.

cd %HYPRE%
set HYPRE=%CD%
git clean -dxf

:: configure hypre
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo configuring hypre version %HYPREVERSION%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.

set BUILDDIR=%HYPRE%\src\cmbuild
cd %BUILDDIR%
cmake ..\  ^
-G "MinGW Makefiles" ^
-DCMAKE_INSTALL_PREFIX="%INSTALLDIR%" ^
-DCMAKE_C_COMPILER=icx ^
-DCMAKE_C_FLAGS="/DWIN32 -O3 /fp:precise" ^
-DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded"

:: build and install hypre
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo building and installing hypre version %HYPREVERSION%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
call make install

if %have_setx% == 0 goto else_setx
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo setting HYPRE_HOME environment variable to %INSTALLDIR%
setx HYPRE_HOME %INSTALLDIR%
echo note: the environment variable HYPRE_HOME takes effect after opening a new command shell
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
goto endif_setx
:else_setx
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo set environment variable HYPRE_HOME to %INSTALLDIR%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
:endif_setx

echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo hypre version %HYPREVERSION% installed in %INSTALLDIR%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
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
echo -help           - display this message
exit /b

:eof