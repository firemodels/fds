@echo off
REM Cross-platform setup script for Python 3 environment (Windows version)
REM Usage: call setup_python3.bat

setlocal enabledelayedexpansion


REM === Check for python ===
where python >nul 2>nul
if errorlevel 1 (
    set ERROR_MSG=python is not installed or not in PATH
    call :ERROR_EXIT
)

REM === Ensure it's Python 3 ===
for /f "delims=" %%v in ('python -c "import sys; print(sys.version_info[0])"') do (
    if not %%v==3 (
        set ERROR_MSG=Python 3 is required but not found
        call :ERROR_EXIT
    )
)

REM === Save current directory ===
set CURDIR=%cd%

REM === Navigate to repo root ===
REM NOTE: Assumes script is 3 levels deep. Change if structure changes.
cd /d "%~dp0..\..\.."
if errorlevel 1 (
    set ERROR_MSG=Failed to navigate to repo root
    call :ERROR_EXIT
)
set REPOROOT=%cd%

REM === Move to .github folder ===
cd /d "%REPOROOT%\fds\.github"
if errorlevel 1 (
    set ERROR_MSG=Directory not found: %REPOROOT%\fds\.github
    call :ERROR_EXIT
)

REM === Create venv if not exists ===
set VENV_DIR=fds_python_env
if not exist "%VENV_DIR%" (
    python -m venv "%VENV_DIR%"
    if errorlevel 1 (
        set ERROR_MSG=Failed to create virtual environment
        call :ERROR_EXIT
    )
)

REM === Check if activate.bat exists ===
if not exist "%VENV_DIR%\Scripts\activate.bat" (
    set ERROR_MSG=Virtual environment activation script not found
    call :ERROR_EXIT
)

REM === Activate venv ===
call "%VENV_DIR%\Scripts\activate.bat"
if errorlevel 1 (
    set ERROR_MSG=Failed to activate virtual environment
    call :ERROR_EXIT
)

REM === Upgrade pip and install requirements ===
python -m pip install --upgrade pip
if errorlevel 1 (
    set ERROR_MSG=Failed to upgrade pip
    call :ERROR_EXIT
)

python -m pip install numpy<2

if exist requirements.txt (
    python -m pip install -r requirements.txt
    if errorlevel 1 (
        set ERROR_MSG=Failed to install requirements
        call :ERROR_EXIT
    )
)

REM === Install Spyder IDE (optional) ===
python -m pip install spyder
if errorlevel 1 (
    set ERROR_MSG=Failed to install Spyder
    call :ERROR_EXIT
)

REM === Set PYTHONPATH ===
set PYTHONPATH=%REPOROOT%\fds\Utilities\Python;%PYTHONPATH%

REM === Run test script ===
cd /d "%REPOROOT%\fds\Utilities\Python"
if not exist hello_world.py (
    set ERROR_MSG=hello_world.py not found
    call :ERROR_EXIT
)

python hello_world.py
if errorlevel 1 (
    set ERROR_MSG=hello_world.py failed
    call :ERROR_EXIT
)

REM === Return to original directory ===
cd /d "%CURDIR%"

echo.
echo Python environment setup complete.

REM === Function-like error handler using variable ===
:ERROR_EXIT
if "%ERROR_MSG%"=="" (
    echo 
) else (
    echo ***Error: %ERROR_MSG%
)
exit /b 1


endlocal
