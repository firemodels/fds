@echo off

set INSTALLDIR=C:\sundials_test
set SUNDIALSVERSION=v6.7.0

set CURDIR=%CD%

echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo setting up Intel compilers
echo ----------------------------------------------------------
echo ----------------------------------------------------------
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

  cd ..\..\..
  git clone https://github.com/LLNL/sundials.git
  cd sundials
echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo checking out version %SUNDIALSVERSION%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
  git checkout %SUNDIALSVERSION%
  cd %CURDIR%
:endif1

echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo cleaning sundials repo
echo ----------------------------------------------------------
echo ----------------------------------------------------------

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
make 

echo ----------------------------------------------------------
echo ----------------------------------------------------------
echo installing sundials version %SUNDIALSVERSION% in %INSTALLDIR%
echo ----------------------------------------------------------
echo ----------------------------------------------------------
make install

echo sundials version %SUNDIALSVERSION% installed in %INSTALLDIR%

cd %CURDIR%