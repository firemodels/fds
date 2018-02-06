#!/bin/bash

# this script is run from a script in Build/Bundle/linux or
# Build/Bundle/osx

if [ -e ~/FDS_SMV_ENVpc.sh ]; then
  source ~/FDS_SMV_ENVpc.sh
else
  if [ -e ~/FDS_SMV_ENV.sh ]; then
    source ~/FDS_SMV_ENV.sh
  else
    source ../scripts/FDS_SMV_ENV.sh
  fi
fi

