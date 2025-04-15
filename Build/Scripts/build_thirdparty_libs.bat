@echo off

:: Select the compilers
call ..\Scripts\set_compilers.bat

set clean_hypre=
set clean_sundials=
set clean_hdf5=
set clean_fds=
set no_libs=
set stopscript=0
set FDS_BUILDDIR=%CD%

:: Define an array of allowed GPU architectures
::set allowed_gpus=("cuda")

call :getopts %*
if %stopscript% == 1 exit /b

goto eof

:: -------------------------------------------------------------
:: Function to check if a value exists in the allowed list
::is_valid_gpu
:: -------------------------------------------------------------
:: set gpu="$1"
:: for /l %%x in (1, 1, 1) do (
::   if allowed_gpus[i] EQU gpu (
::     return 0
::   )
:: )
::EXIT /b 1

:: -------------------------------------------------------------
:getopts
:: -------------------------------------------------------------
 set stopscript=0
 if (%1)==() exit /b
 set valid=0
 set arg=%1
 if /I "%1" EQU "--clean-fds" (
   set clean_fds=--clean-fds
   set valid=1
 )
 if /I "%1" EQU "--clean-hypre" (
   set clean_hypre=--clean-hypre
   set valid=1
 )
 if /I "%1" EQU "bot" (
   set valid=1
 )
 if /I "%1" EQU "--clean-sundials" (
   set clean_sundials=--clean-sundials
   set valid=1
 )
 if /I "%1" EQU "--clean-hdf5" (
   set clean_hdf5=--clean-hdf5
   set valid=1
 )
 if /I "%1" EQU "--clean-all" (
   set clean_fds=--clean-fds
   set clean_hypre=--clean-hypre
   set clean_sundials=--clean-sundials
   set clean_hdf5=--clean-hdf5
   set valid=1
 )
 if /I "%1" EQU "--with-gpu" (
   set BUILD_WITH_GPU=ON
   set GPU_HOME=%CUDA_PATH%
   set GPU_ARCH=cuda
   set with_gpu=true
   ::set gpu_arch="${1#*=}"
   :: Check if the value is in the allowed list
   ::  if (call :is_valid_gpu "$gpu_arch") (
   ::    echo "Error: Invalid option for --with-gpu. Allowed values: ${allowed_gpus[*]}." >&2
   ::    exit 1
   ::  )
   ::set BUILD_WITH_GPU=ON
   ::set GPU_ARCH="$gpu_arch"
   ::if $(GPU_ARCH) EQU cuda (
   ::   set GPU_HOME=%CUDA_PATH%
   ::)
   set valid=1
)
 if /I "%1" EQU "--no-libs" (
   set no_libs=--no-libs
   set clean_fds=--clean-fds
   set valid=1
 )
 if /I "%1" EQU "-help" (
   call :usage
   set stopscript=1
   exit /b
 )
 if /I "%1" EQU "--help" (
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
echo build 3rd party libraries
echo.
echo --clean-all      - rebuild all libraries, remove .obj and .mod files from the fds build directory
echo --clean-fds      - remove .obj and .mod files from the fds build directory
echo --clean-hypre    - rebuild hypre library
echo --clean-sundials - rebuild sundials library
echo --clean-hdf5     - rebuild hdf5 library
echo --no-libs        - build without thirdparty libs
echo --help           - display this message
echo --with-gpu arch  - build with gpu with cuda
exit /b

:eof

set CURDIR_3RDPARTY=%CD%

if "x%clean_fds%" == "x" goto endif1
  echo. 
  echo removing .obj and .mod files from %CD%
  erase *.obj *.mod > Nul 2> Nul
:endif1

set SCRIPTDIR=%~dp0
cd %SCRIPTDIR%
SET SCRIPTDIR=%CD%

if "x%no_libs%" == "x" (
    :: Call HYPRE and SUNDIALS build script
    cd %SCRIPTDIR%\HYPRE
    call build_hypre %clean_hypre%

    cd %SCRIPTDIR%\SUNDIALS
    call build_sundials %clean_sundials%

    cd %SCRIPTDIR%\HDF5
    call build_hdf5 %clean_hdf5%
) else (
    set SUNDIALS_HOME=
    set HYPRE_HOME=
    set HDF5_HOME=
    echo Building FDS without third-party libraries.
)

cd %CURDIR_3RDPARTY%
