#!/bin/bash

echo "Checking for hypre repository..."
dir=`pwd`

if [ -d "$FIREMODELS/hypre" ]; then
  echo "Hypre repository exists. Building hypre library."
  cd $FIREMODELS/hypre/src
  export HYPRE_VERSION=$(git describe)
  cp $FIREMODELS/fds/Build/Scripts/HYPRE/confmake_impi_intel_linux.sh .
  ./confmake_impi_intel_linux.sh
  cd $dir
  export HYPRE_HOME=$FIREMODELS/libs/hypre/$HYPRE_VERSION
else
  echo "Hypre repository does not exist."
fi
