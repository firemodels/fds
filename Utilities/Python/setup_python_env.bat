@echo off
REM Configure Python Virtual Environment for FDS (WINDOWS)
REM Usage (Windows Command Prompt):
REM    1. Install Python 3.7 or higher and add its "python.exe" folder to PATH
REM    2. Open a normal Windows Command Prompt
REM    3. call setup_python_env.bat [--batchmode]
REM
REM For PowerShell, use setup_python_env.ps1 instead.

REM === Parse command-line arguments ===
set "BATCHMODE=false"
for %%a in (%*) do (
    if /i "%%~a"=="--batchmode" set "BATCHMODE=true"
)

REM === Check for python ===
where python >nul 2>nul
if errorlevel 1 (
    set "ERROR_MSG=python is not installed or not in PATH"
    goto :ERROR_EXIT
)

REM === Ensure Python version >= 3.7 ===
set "PY_MAJOR="
set "PY_MINOR="
for /f "tokens=1,2" %%a in ('python -c "import sys; print(sys.version_info.major, sys.version_info.minor)"') do (
    set "PY_MAJOR=%%a"
    set "PY_MINOR=%%b"
)

if not defined PY_MAJOR (
    set "ERROR_MSG=Failed to determine Python version"
    goto :ERROR_EXIT
)

if %PY_MAJOR% LSS 3 (
    set "ERROR_MSG=Python 3.7 or higher is required. Found %PY_MAJOR%.%PY_MINOR%."
    goto :ERROR_EXIT
)
if %PY_MAJOR%==3 if %PY_MINOR% LSS 7 (
    set "ERROR_MSG=Python 3.7 or higher is required. Found %PY_MAJOR%.%PY_MINOR%."
    goto :ERROR_EXIT
)

echo Python version is OK: %PY_MAJOR%.%PY_MINOR%

REM === Save current directory and navigate to .github folder ===
set "CURDIR=%CD%"
cd /d "%~dp0..\..\.."
if errorlevel 1 (
    set "ERROR_MSG=Failed to navigate to repo root"
    goto :ERROR_EXIT
)
set "REPOROOT=%CD%"
set "GITHUB_DIR=%REPOROOT%\fds\.github"

cd /d "%GITHUB_DIR%"
if errorlevel 1 (
    set "ERROR_MSG=Directory not found: %GITHUB_DIR%"
    goto :ERROR_EXIT
)

REM === Setup virtual environment ===
set "VENV_DIR=%GITHUB_DIR%\fds_python_env"
set "INSTALL_REQUIREMENTS=false"

if not exist "%VENV_DIR%\" goto :CREATE_VENV

if /i "%BATCHMODE%"=="true" (
    REM Preserve existing Windows batch-mode behavior: do not delete an existing environment.
    echo Batch mode: activating existing virtual environment without prompts or deletion.
    goto :ACTIVATE_ENV
)

echo Virtual environment "%VENV_DIR%" already exists.
set "choice="
set /p "choice=Do you want to reinstall everything? (y/N): "
if /i "%choice%"=="y" goto :REINSTALL_VENV
if /i "%choice%"=="yes" goto :REINSTALL_VENV

echo Activating existing environment...
goto :ACTIVATE_ENV

:REINSTALL_VENV
REM Deactivate any currently active venv before removing/recreating this one.
if defined VIRTUAL_ENV (
    if exist "%VIRTUAL_ENV%\Scripts\deactivate.bat" call "%VIRTUAL_ENV%\Scripts\deactivate.bat"
)

echo Removing old environment...
rmdir /s /q "%VENV_DIR%"
if exist "%VENV_DIR%\" (
    set "ERROR_MSG=Failed to remove existing virtual environment"
    goto :ERROR_EXIT
)

:CREATE_VENV
echo Creating new virtual environment...
python -m venv "%VENV_DIR%"
if errorlevel 1 (
    set "ERROR_MSG=Failed to create virtual environment"
    goto :ERROR_EXIT
)
set "INSTALL_REQUIREMENTS=true"

:ACTIVATE_ENV
REM === Activate environment in the current cmd.exe session ===
if not exist "%VENV_DIR%\Scripts\activate.bat" (
    set "ERROR_MSG=Virtual environment activation script not found: %VENV_DIR%\Scripts\activate.bat"
    goto :ERROR_EXIT
)
call "%VENV_DIR%\Scripts\activate.bat"
if errorlevel 1 (
    set "ERROR_MSG=Failed to activate virtual environment"
    goto :ERROR_EXIT
)

REM === Upgrade pip and install requirements if flagged ===
if /i "%INSTALL_REQUIREMENTS%"=="true" (
    echo Installing/updating required Python packages...
    python -m pip install --upgrade pip
    if errorlevel 1 (
        set "ERROR_MSG=Failed to upgrade pip"
        goto :ERROR_EXIT
    )

    if exist "%GITHUB_DIR%\requirements.txt" (
        python -m pip install -r "%GITHUB_DIR%\requirements.txt"
        if errorlevel 1 (
            set "ERROR_MSG=Failed to install requirements"
            goto :ERROR_EXIT
        )
    )
)

REM === Set PYTHONPATH ===
set "FDS_PYTHON_PATH=%REPOROOT%\fds\Utilities\Python"
if defined PYTHONPATH (
    set "PYTHONPATH=%FDS_PYTHON_PATH%;%PYTHONPATH%"
) else (
    set "PYTHONPATH=%FDS_PYTHON_PATH%"
)

REM === Run test script ===
cd /d "%FDS_PYTHON_PATH%"
if errorlevel 1 (
    set "ERROR_MSG=Failed to find script directory: %FDS_PYTHON_PATH%"
    goto :ERROR_EXIT
)
if not exist "hello_world.py" (
    set "ERROR_MSG=hello_world.py not found"
    goto :ERROR_EXIT
)
python hello_world.py
if errorlevel 1 (
    set "ERROR_MSG=hello_world.py failed"
    goto :ERROR_EXIT
)

REM === Return to original directory ===
cd /d "%CURDIR%"

echo.
echo Python environment setup complete.

REM Clean up setup-only variables.  Keep PATH, VIRTUAL_ENV and PYTHONPATH.
set "BATCHMODE="
set "INSTALL_REQUIREMENTS="
set "PY_MAJOR="
set "PY_MINOR="
set "choice="
set "FDS_PYTHON_PATH="
set "GITHUB_DIR="
set "VENV_DIR="
set "REPOROOT="
set "CURDIR="
set "ERROR_MSG="
exit /b 0

:ERROR_EXIT
if defined CURDIR cd /d "%CURDIR%" >nul 2>nul
echo.
echo *** Error: %ERROR_MSG%
exit /b 1
