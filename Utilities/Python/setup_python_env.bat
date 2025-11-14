@echo off
REM Configure Spyder IDE for FDS (WINDOWS)
REM Usage (WINDOWS): 
REM    1. Install Python 3.11 or below and add its "python.exe" folder to the PATH
REM    2. Open normal Windows Command Prompt
REM    3. call ./setup_python_env.bat [--batchmode]

setlocal enabledelayedexpansion

REM === Parse command-line arguments ===
set BATCHMODE=false
for %%a in (%*) do (
    if "%%a"=="--batchmode" set BATCHMODE=true
)

REM === Check for python ===
where python >nul 2>nul
if errorlevel 1 (
    set ERROR_MSG=python is not installed or not in PATH
    call :ERROR_EXIT
)

REM === Ensure Python version >= 3.7 ===
for /f "delims=" %%v in ('python -c "import sys; print(sys.version_info[0])"') do set PY_MAJOR=%%v
for /f "delims=" %%v in ('python -c "import sys; print(sys.version_info[1])"') do set PY_MINOR=%%v

if %PY_MAJOR% LSS 3 (
    set ERROR_MSG=Python 3.7 or higher is required. Found %PY_MAJOR%.%PY_MINOR%.
    call :ERROR_EXIT
)
if %PY_MAJOR%==3 if %PY_MINOR% LSS 7 (
    set ERROR_MSG=Python 3.7 or higher is required. Found %PY_MAJOR%.%PY_MINOR%.
    call :ERROR_EXIT
)

echo Python version is OK: %PY_MAJOR%.%PY_MINOR%

REM === Save current directory and navigate to .github folder ===
set CURDIR=%cd%
cd /d "%~dp0..\..\.."
if errorlevel 1 (
    set ERROR_MSG=Failed to navigate to repo root
    call :ERROR_EXIT
)
set REPOROOT=%cd%
cd /d "%REPOROOT%\fds\.github"
if errorlevel 1 (
    set ERROR_MSG=Directory not found: %REPOROOT%\fds\.github
    call :ERROR_EXIT
)

REM === Setup virtual environment ===
set VENV_DIR=fds_python_env
set INSTALL_REQUIREMENTS=false

if exist "%VENV_DIR%" (
    if "%BATCHMODE%"=="true" (
        echo Batch mode: activating existing virtual environment without prompts or deletion.
    ) else (
        echo ⚠ Virtual environment "%VENV_DIR%" already exists.
        set /p choice="Do you want to reinstall everything? (y/N): "
        if /i "!choice!"=="y" (
            REM Deactivate if currently active (only works if environment was manually activated)
            if defined VIRTUAL_ENV (
                call deactivate
            )
            echo Removing old environment...
            rmdir /s /q "%VENV_DIR%" || (
                set ERROR_MSG=Failed to remove existing virtual environment
                call :ERROR_EXIT
            )
            echo Creating new virtual environment...
            python -m venv "%VENV_DIR%" || (
                set ERROR_MSG=Failed to create virtual environment
                call :ERROR_EXIT
            )
            set INSTALL_REQUIREMENTS=true
        ) else (
            echo Activating existing environment...
        )
    )
) else (
    echo Creating new virtual environment...
    python -m venv "%VENV_DIR%" || (
        set ERROR_MSG=Failed to create virtual environment
        call :ERROR_EXIT
    )
    set INSTALL_REQUIREMENTS=true
)

REM === Activate environment ===
if not exist "%VENV_DIR%\Scripts\activate.bat" (
    set ERROR_MSG=Virtual environment activation script not found
    call :ERROR_EXIT
)
call "%VENV_DIR%\Scripts\activate.bat" || (
    set ERROR_MSG=Failed to activate virtual environment
    call :ERROR_EXIT
)

REM === Upgrade pip and install requirements if flagged ===
if "%INSTALL_REQUIREMENTS%"=="true" (
    echo Installing/updating required Python packages...
    python -m pip install --upgrade pip || (
        set ERROR_MSG=Failed to upgrade pip
        call :ERROR_EXIT
    )
    if exist requirements.txt (
        python -m pip install -r requirements.txt || (
            set ERROR_MSG=Failed to install requirements
            call :ERROR_EXIT
        )
    )
)

REM === Set PYTHONPATH ===
set PYTHONPATH=%REPOROOT%\fds\Utilities\Python;%PYTHONPATH%

REM === Run test script ===
cd /d "%REPOROOT%\fds\Utilities\Python"
if not exist hello_world.py (
    set ERROR_MSG=hello_world.py not found
    call :ERROR_EXIT
)
python hello_world.py || (
    set ERROR_MSG=hello_world.py failed
    call :ERROR_EXIT
)

REM === Return to original directory ===
cd /d "%CURDIR%"

echo.
echo Python environment setup complete.

REM === Error handler ===
:ERROR_EXIT
if "%ERROR_MSG%"=="" (
    echo.
) else (
    echo *** Error: %ERROR_MSG%
)
exit /b 1

endlocal
