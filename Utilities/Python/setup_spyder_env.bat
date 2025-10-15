@echo off
REM Cross-platform setup script for Python 3 environment (Windows version)
REM Usage: call setup_python3.bat
setlocal enabledelayedexpansion

REM === Function-like error handler ===
:ERROR_EXIT
echo ***Error: %~1
exit /b 1

REM === Check for python3 ===
where python >nul 2>nul
if errorlevel 1 (
    call :ERROR_EXIT "python is not installed or not in PATH"
)

REM === Save current directory ===
set CURDIR=%cd%

REM === Navigate to repo root ===
cd /d "%~dp0..\..\.."
set REPOROOT=%cd%

REM === Move to .github folder ===
cd /d "%REPOROOT%\fds\.github"
if errorlevel 1 (
    call :ERROR_EXIT "Directory not found: %REPOROOT%\fds\.github"
)

REM === Create venv if not exists ===
set VENV_DIR=fds_python_env
if not exist "%VENV_DIR%" (
    python -m venv "%VENV_DIR%"
    if errorlevel 1 (
        call :ERROR_EXIT "Failed to create virtual environment"
    )
)

REM === Activate venv ===
call "%VENV_DIR%\Scripts\activate"
if errorlevel 1 (
    call :ERROR_EXIT "Failed to activate virtual environment"
)

REM === Install dependencies ===
python -m pip install --upgrade pip
if exist requirements.txt (
    python -m pip install -r requirements.txt
    if errorlevel 1 (
        call :ERROR_EXIT "Failed to install requirements"
    )
)
python -m pip install spyder

REM === Set PYTHONPATH ===
set PYTHONPATH=%REPOROOT%\fds\Utilities\Python;%PYTHONPATH%

REM === Run test script ===
cd /d "%REPOROOT%\fds\Utilities\Python"
if not exist hello_world.py (
    call :ERROR_EXIT "hello_world.py not found"
)
python hello_world.py
if errorlevel 1 (
    call :ERROR_EXIT "hello_world.py failed"
)

REM === Return to original directory ===
cd /d "%CURDIR%"

echo.
echo Python environment setup complete.
endlocal
