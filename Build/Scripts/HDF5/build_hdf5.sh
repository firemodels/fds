#!/bin/bash
HDF5_LIB_TAG=hdf5_1.14.5

CONFMAKE=$1
CLEAN_HDF5=$2

dir=`pwd`

echo "CLEAN_HDF5 = $CLEAN_HDF5"
if [ "$CLEAN_HDF5" = true ]; then
  echo "Removing hdf5 library ..."
  rm -r $FIREMODELS/libs/hdf5
  rm -r $FIREMODELS/hdf5/BUILDDIR
fi

echo "Checking for hdf5 library..."

if [ -d "$FIREMODELS/libs/hdf5" ]; then
  echo "HDF5 library exists.  Skipping HDF5 build."
  # List all directories under $FIREMODELS/libs/hdf5
  hdf5_lib_dir=$(ls -d $FIREMODELS/libs/hdf5/*/)
  # Extract the version part (removes the leading path)
  HDF5_VERSION=$(basename $hdf5_lib_dir)
  export HDF5_HOME=$FIREMODELS/libs/hdf5/$HDF5_VERSION
  echo "HDF5 library:" $FIREMODELS/libs/hdf5/$HDF5_VERSION
  return 0
else
  echo "HDF5 library does not exist."
fi

echo "Checking for hdf5 repository..."

if [ -d "$FIREMODELS/hdf5" ]; then
  echo "hdf5 repository exists. Building hdf5 library."
  cd $FIREMODELS/hdf5
  # Handle possible corrupted state of repository
  if git branch | grep -q "* develop"; then
    echo "On develop branch"
  else
    git checkout develop
  fi
  git checkout .
  git clean -dxf
  if [[ "$(git tag -l $HDF5_LIB_TAG)" == $HDF5_LIB_TAG ]]; then
    echo "Checking out $HDF5_LIB_TAG"
    git checkout $HDF5_LIB_TAG
  fi 

  mkdir $FIREMODELS/hdf5/BUILDDIR
  cd $FIREMODELS/hdf5/BUILDDIR
  echo "Creating library directry..."
  export HDF5_VERSION=$(git describe)
  mkdir $FIREMODELS/libs/hdf5/$HDF5_VERSION
  echo "Cleaning hdf5 repository..."
  rm -r $FIREMODELS/hdf5/BUILDDIR/*
  cp $FIREMODELS/fds/Build/Scripts/HDF5/$CONFMAKE .
  ./$CONFMAKE
  # get back from detached HEAD state
  cd $FIREMODELS/hdf5
  git checkout develop
  cd $dir
  export HDF5_HOME=$FIREMODELS/libs/hdf5/$HDF5_VERSION
  echo "HDF5 library:" $FIREMODELS/libs/hdf5/$HDF5_VERSION
  return 0
else
  echo "HDF5 repository does not exist."
fi
