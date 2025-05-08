@echo off
setlocal enabledelayedexpansion

set REQUIRED_MAJOR=3
set REQUIRED_MINOR=21

:: Check if cmake is installed
where cmake >nul 2>nul
if errorlevel 1 (
    echo Error: CMake is not installed. Please install CMake version 3.21.0 or newer.
    exit /b 1
)

:: Get installed cmake version
for /f "tokens=3" %%v in ('cmake --version') do (
    set VERSION=%%v
    goto :parse_version
)

:parse_version
for /f "tokens=1,2,3 delims=." %%a in ("!VERSION!") do (
    set MAJOR=%%a
    set MINOR=%%b
    set PATCH=%%c
)

:: Compare version
if !MAJOR! LSS %REQUIRED_MAJOR% (
    echo Error: Installed CMake version is !VERSION!. Version 3.21.0 or newer is required.
    exit /b 1
)

if !MAJOR! EQU %REQUIRED_MAJOR% if !MINOR! LSS %REQUIRED_MINOR% (
    echo Error: Installed CMake version is !VERSION!. Version 3.21.0 or newer is required.
    exit /b 1
)

echo CMake version !VERSION! is sufficient.
exit /b 0
