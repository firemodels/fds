@echo off
set LIB_TAG=hdf5_1.14.5

::*** library and tag name are the same

set LIB_DIR=%LIB_TAG%

::*** parse options

set clean_hdf5=
call :getopts %*
if %stopscript% == 1 exit /b

::*** make sure cmake and make are installed

set abort=0
set buildstatus=
call :is_file_installed cmake || set abort=1
call :is_file_installed make  || set abort=1
if %abort% == 1 exit /b

set CURDIR=%CD%

::*** define root directory where fds repo and libs directories are located

set FIREMODELS=..\..\..\..
cd %FIREMODELS%
set FIREMODELS=%CD%
cd %CURDIR%

set INSTALLDIR=%FIREMODELS%\libs\hdf5\%LIB_DIR%

::*** erase directory if it exists and clean option was specified

if "x%clean_hdf5%" == "x" goto endif1
  if exist %INSTALLDIR% rmdir /s /q %INSTALLDIR%
:endif1

::*** if hdf5 library directory exists exit and use it

if not exist %INSTALLDIR% goto endif2
  set HDF5_HOME=%INSTALLDIR%
  set buildstatus=prebuilt
  goto eof
:endif2

::*** hdf5 library doesn't exist, if hdf5 repo exists build hdf5 library

set LIB_REPO=%FIREMODELS%\hdf5
if exist %LIB_REPO% goto buildlib

::*** if directory pointed to by SUNDIALS_HOME exists exit and use it

if "x%HDF5_HOME%" == "x" goto else4
if not exist %HDF5_HOME%  goto else4
    set buildstatus=prebuilt
    goto endif4
:else4
  set buildstatus=norepo
:endif4
goto eof

::*** if we've gotten this far the prebuilt libraries do not exist, the repo does exist so build the hdf5 library

:buildlib
cd %LIB_REPO%

set buildstatus=build
echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo building hdf5 library version %LIB_TAG%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.

git checkout %LIB_TAG%

echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo cleaning hdf5 repo
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.

cd %LIB_REPO%
set BUILDDIR=%LIB_REPO%\BUILDDIR
git clean -dxf

mkdir %BUILDDIR%
cd %BUILDDIR%

::*** configure hdf5

echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo configuring hdf5 version %LIB_TAG%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.

cmake ..\  ^
-G "MinGW Makefiles" ^
-DCMAKE_INSTALL_PREFIX="%INSTALLDIR%" ^
-DCMAKE_C_COMPILER=%COMP_CC% ^
-DCMAKE_Fortran_COMPILER=%COMP_FC% ^
-DCMAKE_BUILD_TYPE="Release" ^
-DBUILD_SHARED_LIBS="OFF" ^
-DBUILD_TESTING="OFF" ^
-DHDF5_BUILD_TOOLS="OFF" ^
-DHDF5_BUILD_FORTRAN="ON" ^
-DHDF5_ENABLE_PARALLEL="ON" ^
-DCMAKE_MAKE_PROGRAM="%CMAKE_MAKE_PROGRAM%" ^
-DCMAKE_C_FLAGS_RELEASE="${CMAKE_C_FLAGS_RELEASE} /MT" ^
-DCMAKE_C_FLAGS_DEBUG="${CMAKE_C_FLAGS_DEBUG} /MTd"

::*** build and install hdf5

echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo building hdf5 version %LIB_TAG%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
call make 

echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo installing hdf5 version %LIB_TAG% in %INSTALLDIR%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
call make install

echo.
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo setting HDF5_HOME environment variable to %INSTALLDIR%
set HDF5_HOME=%INSTALLDIR%
echo.

echo The hdf5 library version %LIB_TAG% was built and installed in %INSTALLDIR%
echo.

echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo removing .obj and .mod files from Windows fds build directories
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
for /D %%f in (%FIREMODELS%\fds\Build\*win*) do (
  cd %%f
  erase *.obj *.mod > Nul 2> Nul
)

cd %CURDIR%

goto eof

:: -------------------------------------------------------------
:getopts
:: -------------------------------------------------------------
 set stopscript=0
 if (%1)==() exit /b
 set valid=0
 set arg=%1
if /I "%1" EQU "--clean-hdf5" (
    set clean_hdf5=1
    set valid=1
 )
if /I "%1" EQU "--help" (
   call :usage
   set stopscript=1
   exit /b
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
:usage  
:: -------------------------------------------------------------
echo build hdf5
echo. 
echo --clean-hdf5 - force build of hdf5 library
echo --help           - display this message
exit /b

:eof
echo.
if "%buildstatus%" == "norepo"   echo The hdf5 git repo does not exist, The hdf5 library was not built.  FDS will be built without it.
if "%buildstatus%" == "prebuilt" echo The hdf5 library exists. Skipping hdf5 build. FDS will be built using the
if "%buildstatus%" == "prebuilt" echo hdf5 library in %HDF5_HOME%
echo.
