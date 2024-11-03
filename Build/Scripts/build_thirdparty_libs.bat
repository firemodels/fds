@echo off

set clean_hypre=
set clean_sundials=
set FDS_BUILDDIR=%CD%

call :getopts %*
if %stopscript% == 1 exit /b


goto eof

:: -------------------------------------------------------------
:getopts
:: -------------------------------------------------------------
 set stopscript=0
 if (%1)==() exit /b
 set valid=0
 set arg=%1
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
 if /I "%1" EQU "--clean-all" (
   set clean_hypre=--clean-hypre
   set clean_sundials=--clean-sundials
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
:usage  
:: -------------------------------------------------------------
echo build 3rd party libraries
echo.
echo --clean-all      - rebuild all libraries
echo --clean-hypre    - rebuild hypre library
echo --clean-sundials - rebuild sundials library
echo -help            - display this message
exit /b

:eof

set CURDIR_3RDPARTY=%CD%
set SCRIPTDIR=%~dp0
cd %SCRIPTDIR%
SET SCRIPTDIR=%CD%

cd %SCRIPTDIR%\HYPRE
call build_hypre %clean_hypre%

cd %SCRIPTDIR%\SUNDIALS
call build_sundials %clean_sundials%

cd %CURDIR_3RDPARTY%