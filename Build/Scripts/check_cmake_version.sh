#!/bin/bash

# Required minimum version
REQUIRED_VERSION="3.21.0"

# Function to compare versions
version_ge() {
  [ "$(printf '%s\n' "$2" "$1" | sort -V | head -n1)" = "$2" ]
}

# Check if cmake is installed
if ! command -v cmake &> /dev/null; then
  echo "Error: CMake is not installed. Please install CMake version $REQUIRED_VERSION or newer."
  exit 1
fi

# Get installed cmake version
INSTALLED_VERSION=$(cmake --version | head -n1 | awk '{print $3}')

# Compare versions
if version_ge "$INSTALLED_VERSION" "$REQUIRED_VERSION"; then
  echo "CMake version $INSTALLED_VERSION is sufficient."
else
  echo "Error: Installed CMake version is $INSTALLED_VERSION. Version $REQUIRED_VERSION or newer is required."
  exit 1
fi
