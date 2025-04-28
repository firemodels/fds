#!/bin/bash
HYPRE_LIB_TAG=v2.33.0

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
  cd $FIREMODELS/hypre
  # Handle possible corrupted state of repository
  if git branch | grep -q "* master"; then
    echo "On master branch"
  else
    git checkout master
  fi
  git checkout .
  git clean -dxf
  if [[ "$(git tag -l $HYPRE_LIB_TAG)" == $HYPRE_LIB_TAG ]]; then
    echo "Checking out $HYPRE_LIB_TAG"
    git checkout $HYPRE_LIB_TAG
  else
    echo "Your HYPRE repository is not up to date with the required tag: $HYPRE_LIB_TAG."
    echo "The FDS build requires HYPRE version $HYPRE_LIB_TAG. Please update your HYPRE repository."
    exit 1
  fi
  cd $FIREMODELS/hypre/build
  export HYPRE_VERSION=$(git describe)
  echo "Cleaning hypre repository..."
  rm -r $FIREMODELS/hypre/build/*
  cp $FIREMODELS/fds/Build/Scripts/HYPRE/$CONFMAKE .
  ./$CONFMAKE
  # get back from detached HEAD state
  cd $FIREMODELS/hypre
  git checkout master
  cd $dir
  export HYPRE_HOME=$FIREMODELS/libs/hypre/$HYPRE_VERSION
  echo "Hypre library:" $FIREMODELS/libs/hypre/$HYPRE_VERSION
  return 0
else
  echo "Hypre repository does not exist."
fi
