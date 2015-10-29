#/bin/bash

# MPI distribution location

export MPIDIST=$1

if [[ "$MPIDIST" == "" ]]; then
  echo "*** Warning: the MPI distribution location is not defined."
  echo "Make sure the MPIDIST_ETH and/or MPIDIST_IB environment"
  echo "variables are defined in your startup file"
  exit
fi

if [[ ! -d $MPIDIST ]]; then
  echo "*** Warning: the MPI distribution, $MPIDIST, does not exist"
  echo "*** Compilation aborted"
  MPIDIST=
  exit
fi

# Update LD_LIBRARY_PATH and PATH
LD_LIBRARY_PATH=$MPIDIST/lib:$LD_LIBRARY_PATH
PATH=$MPIDIST/bin:$PATH
export LD_LIBRARY_PATH PATH
