@echo off

:: Initialize compiler flags
set set_COMP_CC=0
set set_COMP_FC=0

:: Check and set C compiler
if defined FIREMODELS_CC (
    set COMP_CC=%FIREMODELS_CC%
    where /q %COMP_CC% && (
        set set_COMP_CC=1
    ) || (
        echo Warning: %FIREMODELS_CC% is not available. Searching for an alternative.
    )
)
if %set_COMP_CC%==0 (
    for %%C in (icx icc) do (
        where /q %%C && set COMP_CC=%%C && set set_COMP_CC=1 && goto :found_cc
    )
    echo Error: Neither icx nor icc is available. & exit /b 1
)
:found_cc

:: Check and set Fortran compiler
if defined FIREMODELS_FC (
    set COMP_FC=%FIREMODELS_FC%
    where /q %COMP_FC% && (
        set set_COMP_FC=1
    ) || (
        echo Warning: %FIREMODELS_FC% is not available. Searching for an alternative.
    )
)
if %set_COMP_FC%==0 (
    for %%F in (ifx ifort) do (
        where /q %%F && set COMP_FC=%%F && set set_COMP_FC=1 && goto :found_fc
    )
    echo Error: Neither ifx nor ifort is available. & exit /b 1
)
:found_fc

:: Display selected compilers
echo.
echo Third-party libs C Compiler: %COMP_CC%
echo Firemodels and Third-party libs Fortran Compiler: %COMP_FC%
echo.


::Check if make.bat or make.exe exists, and set CMAKE_MAKE_PROGRAM accordingly
::------------------------------------------------------------------------------
set CMAKE_MAKE_PROGRAM=
for /f "delims=" %%i in ('where make.bat 2^>nul') do set CMAKE_MAKE_PROGRAM=%%i
if not defined CMAKE_MAKE_PROGRAM (
    for /f "delims=" %%i in ('where make.exe 2^>nul') do set CMAKE_MAKE_PROGRAM=%%i
)
if not defined CMAKE_MAKE_PROGRAM (
    echo Error: Neither make.bat nor make.exe found in PATH.
    exit /b 1
)

echo.
echo make proram is %CMAKE_MAKE_PROGRAM%
echo.


