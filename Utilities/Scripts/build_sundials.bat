@echo off

set abort=0

call :is_file_installed cmake || set abort=1
if %abort% == 1 exit /b

call :getopts %*
if %stopscript% == 1 exit /b

set INSTALLDIR=C:\sundials_test
set SUNDIALSVERSION=v6.7.0

set CURDIR=%CD%

echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo setting up Intel compilers
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
call ..\..\Build\Scripts\setup_intel_compilers.bat

cd %CURDIR%

set SUNDIALS=..\..\..\sundials

:: clone sundials repo (at same level as fds, smv etc repos) if it doesn't exist
if exist %SUNDIALS% goto endif1
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo cloning sundials from https://github.com/LLNL/sundials.git
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.

  cd ..\..\..
  git clone https://github.com/LLNL/sundials.git
  cd sundials
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo checking out version %SUNDIALSVERSION%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
  git checkout %SUNDIALSVERSION%
  cd %CURDIR%
:endif1

echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo cleaning sundials repo
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.

cd %SUNDIALS%
set SUNDIALS=%CD%
set BUILDDIR=%SUNDIALS%\BUILDDIR
git clean -dxf

mkdir %BUILDDIR%
cd %BUILDDIR%

:: configure sundials
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo configuring sundials version %SUNDIALSVERSION%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.

cmake ..\  ^
-G "MinGW Makefiles" ^
-DCMAKE_INSTALL_PREFIX="%INSTALLDIR%" ^
-DEXAMPLES_INSTALL_PATH="%INSTALLDIR%\examples" ^
-DBUILD_FORTRAN_MODULE_INTERFACE=ON ^
-DCMAKE_C_COMPILER=icx ^
-DCMAKE_Fortran_COMPILER=ifort ^
-DEXAMPLES_ENABLE_C=OFF ^
-DEXAMPLES_ENABLE_CXX=OFF ^
-DEXAMPLES_ENABLE_F2003=OFF ^
-DENABLE_OPENMP=ON ^
-DBUILD_SHARED_LIBS=OFF ^
-DCMAKE_C_FLAGS_RELEASE="${CMAKE_C_FLAGS_RELEASE} /MT" ^
-DCMAKE_C_FLAGS_DEBUG="${CMAKE_C_FLAGS_DEBUG} /MTd"

:: build and install sundials
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo building sundials version %SUNDIALSVERSION%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
make 

echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo installing sundials version %SUNDIALSVERSION% in %INSTALLDIR%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
make install

echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo setting SUNDIALSHOME environment variable to %INSTALLDIR%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo.
setx SUNDIALSHOME %INSTALLDIR%

echo sundials version %SUNDIALSVERSION% installed in %INSTALLDIR%

cd %CURDIR%

goto eof

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
:usage  
:: -------------------------------------------------------------
echo build sundials
echo. 
echo -help           - display this message
exit /b

:eof