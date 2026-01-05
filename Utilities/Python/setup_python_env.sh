#!/bin/bash
# Configure Python Virtual Environment for FDS
# Usage (Linux and macOS):
#    source ./setup_python_env.sh

BATCHMODE=false

# Parse command-line arguments
for arg in "$@"; do
    case $arg in
        --batchmode)
        BATCHMODE=true
        shift
        ;;
    esac
done

function error_exit {
    echo "*** Error: $1"
    return 1
}


# Check Python version > 3.7
major=$(python3 -c 'import sys; print(sys.version_info.major)')
minor=$(python3 -c 'import sys; print(sys.version_info.minor)')

if (( major < 3 || (major == 3 && minor < 7) )); then
    echo "Python 3.7 or higher is required. Found $major.$minor."
    echo "Stopping python environment setup."
    return
else
    echo "Python version is OK: $major.$minor"
fi


# Save current directory to return later
setup_python_env_pwd=$(pwd)

# Go to .github folder
cd "$(dirname "${BASH_SOURCE[0]}")/../../.." || error_exit "Failed to locate repo root"
reporoot=$(pwd)
cd "$reporoot/fds/.github" || error_exit "Directory not found: $reporoot/fds/.github"

VENV_DIR="fds_python_env"
VENV_BIN="$VENV_DIR/bin"
INSTALL_REQUIREMENTS=false

REINSTALL=false

# Check if virtual environment exists
if [ -d "$VENV_DIR" ]; then
    if [ "$BATCHMODE" = true ]; then
        echo "Batch mode: will reinstall virtual environment."
        REINSTALL=true
    else
        echo "Virtual environment '$VENV_DIR' already exists."
        read -p "Do you want to reinstall everything? (y/N): " choice
        case "$choice" in
            [yY]|[yY][eE][sS])
                REINSTALL=true
                ;;
            *)
                echo "Activating existing environment..."
                ;;
        esac
    fi
else
    REINSTALL=true
fi

# Reinstall if requested or required
if [ "$REINSTALL" = true ]; then
    if [ -d "$VENV_DIR" ]; then
        echo "Removing old environment..."
        rm -rf "$VENV_DIR" || error_exit "Failed to remove existing virtual environment"
    fi

    echo "Creating new virtual environment..."
    python3 -m venv "$VENV_DIR" || error_exit "Failed to create virtual environment"
    INSTALL_REQUIREMENTS=true
fi

source "$VENV_BIN/activate" || error_exit "Failed to activate virtual environment"

# Install requirements only if flagged
if [ "$INSTALL_REQUIREMENTS" = true ]; then
    echo "Installing/updating required Python packages..."
    pip install --upgrade pip
    if [ -f "requirements.txt" ]; then
        pip install -r requirements.txt || error_exit "Failed to install requirements"
    fi
fi

export PYTHONPATH="$reporoot/fds/Utilities/Python:$PYTHONPATH"

# Run test script to verify environment
cd "$reporoot/fds/Utilities/Python" || error_exit "Failed to find script directory"
python hello_world.py || error_exit "hello_world.py failed"

# Return to original directory
cd "$setup_python_env_pwd"

echo "Python environment setup complete."

