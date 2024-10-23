#!/bin/bash
CONFMAKE=$1

echo "Checking for hypre repository..."
dir=`pwd`

if [ -d "$FIREMODELS/hypre" ]; then
  echo "Hypre repository exists. Building hypre library."
  cd $FIREMODELS/hypre/src
  export HYPRE_VERSION=$(git describe)
  cp $FIREMODELS/fds/Build/Scripts/HYPRE/$CONFMAKE .
  ./$CONFMAKE
  cd $dir
  export HYPRE_HOME=$FIREMODELS/libs/hypre/$HYPRE_VERSION
else
  echo "Hypre repository does not exist."
fi
