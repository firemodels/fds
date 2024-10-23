#!/bin/bash
CONFMAKE=$1
CLEAN_HYPRE=$2

dir=`pwd`

echo "CLEAN_HYPRE = $CLEAN_HYPRE"
if [ "$CLEAN_HYPRE" = true ]; then
  echo "Removing hypre library ..."
  rm -r $FIREMODELS/libs/hypre
fi

echo "Checking for hypre library..."

if [ -d "$FIREMODELS/libs/hypre" ]; then
  echo "Hypre library exists.  Skipping hypre build."
  # List all directories under $FIREMODELS/libs/hypre
  hypre_lib_dir=$(ls -d $FIREMODELS/libs/hypre/*/)
  # Extract the version part (removes the leading path)
  HYPRE_VERSION=$(basename $hypre_lib_dir)
  export HYPRE_HOME=$FIREMODELS/libs/hypre/$HYPRE_VERSION
  echo "Hypre library:" $FIREMODELS/libs/hypre/$HYPRE_VERSION
  return 0
else
  echo "Hypre library does not exist."
fi

echo "Checking for hypre repository..."

if [ -d "$FIREMODELS/hypre" ]; then
  echo "Hypre repository exists. Building hypre library."
  cd $FIREMODELS/hypre/src
  export HYPRE_VERSION=$(git describe)
  echo "Cleaning hypre repository..."
  make distclean
  cp $FIREMODELS/fds/Build/Scripts/HYPRE/$CONFMAKE .
  ./$CONFMAKE
  cd $dir
  export HYPRE_HOME=$FIREMODELS/libs/hypre/$HYPRE_VERSION
  echo "Hypre library:" $FIREMODELS/libs/hypre/$HYPRE_VERSION
  return 0
else
  echo "Hypre repository does not exist."
fi
