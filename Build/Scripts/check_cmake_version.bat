@echo off
setlocal enabledelayedexpansion

set REQUIRED_MAJOR=3
set REQUIRED_MINOR=21
set VERSION_STRING=%REQUIRED_MAJOR%.%REQUIRED_MINOR%.0

:: Check if cmake is installed
where cmake >nul 2>nul
if errorlevel 1 (
    echo Error: CMake is not installed. Please install CMake version %VERSION_STRING% or newer.
    exit /b 1
)

:: Extract the version number
for /f "tokens=3" %%v in ('cmake --version ^| findstr /i "version"') do (
    set VERSION=%%v
)

:: Parse version into major, minor, patch
for /f "tokens=1,2,3 delims=." %%a in ("!VERSION!") do (
    set MAJOR=%%a
    set MINOR=%%b
    set PATCH=%%c
)

:: Compare major version
if !MAJOR! LSS %REQUIRED_MAJOR% (
    echo Error: Installed CMake version is !VERSION!. Version %VERSION_STRING% or newer is required.
    exit /b 1
)

:: Compare minor version if major matches
if !MAJOR! EQU %REQUIRED_MAJOR% if !MINOR! LSS %REQUIRED_MINOR% (
    echo Error: Installed CMake version is !VERSION!. Version %VERSION_STRING% or newer is required.
    exit /b 1
)

echo CMake version !VERSION! is sufficient.
exit /b 0
