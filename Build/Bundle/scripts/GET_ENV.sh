#!/bin/bash

BUNDLE_HOME=$HOME/.bundle

mkdir -p $BUNDLE_HOME/pubs
mkdir -p $BUNDLE_HOME/BUNDLE
mkdir -p $BUNDLE_HOME/OPENMPMI

# this script is run from a script in Build/Bundle/linux or
# Build/Bundle/osx

if [ -e $BUNDLE_HOME/FDS_SMV_ENVpc.sh ]; then
  source $BUNDLE_HOME/FDS_SMV_ENVpc.sh
else
  if [ -e $BUNDLE_HOME/FDS_SMV_ENV.sh ]; then
    source $BUNDLE_HOME/FDS_SMV_ENV.sh
  else
    source ../scripts/FDS_SMV_ENV.sh
  fi
fi

if [ "$FDS_VERSION" != "" ]; then
  export fds_version=$FDS_VERSION
fi
if [ "$SMV_VERSION" != "" ]; then
  export smv_version=$SMV_VERSION
fi
../scripts/CHECK_VARS.sh

