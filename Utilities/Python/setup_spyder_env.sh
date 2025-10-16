#!/bin/bash
# Configure spyder IDE for FDS
# Usage (Linux and MacOS): 
#    1. source ./setup_spyder_env.sh
#    2. spyder   

function error_exit {
  echo "***Error: $1"
  return 1
}

# Check for python3
if ! command -v python3 >/dev/null 2>&1; then
  error_exit "python3 is not installed"
fi

curdir=$(pwd)
cd "$(dirname "${BASH_SOURCE[0]}")/../../.." || error_exit "Failed to locate repo root"
reporoot=$(pwd)

# Create and activate virtual environment
cd "$reporoot/fds/.github" || error_exit "Directory not found: $reporoot/fds/.github"
VENV_DIR="fds_python_env"
VENV_BIN="$VENV_DIR/bin"
if [ ! -d "$VENV_DIR" ]; then
  python3 -m venv "$VENV_DIR" || error_exit "Failed to create virtual environment"
fi
source "$VENV_BIN/activate"

# Install requirements
pip install --upgrade pip
pip install -r requirements.txt || error_exit "Failed to install requirements"
pip install spyder

# Set PYTHONPATH
export PYTHONPATH="$reporoot/fds/Utilities/Python:$PYTHONPATH"

# Run test script
cd "$reporoot/fds/Utilities/Python" || error_exit "Failed to find script directory"
python hello_world.py || error_exit "hello_world.py failed"

# Return to original directory
cd "$curdir"

echo "Python environment setup complete."

