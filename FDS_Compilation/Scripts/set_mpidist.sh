#/bin/bash

# MPI distribution location

export MPIDIST=$1

if [[ "$MPIDIST" != "" && ! -d $MPIDIST ]]; then
  echo "*** Warning: the MPI distribution, $MPIDIST, does not exist"
  echo "*** Compilation aborted"
  MPIDIST=
  exit
fi

# Update LD_LIBRARY_PATH and PATH
echo Building FDS with the MPI distribution: $MPIDIST
LD_LIBRARY_PATH=$MPIDIST/lib:$LD_LIBRARY_PATH
PATH=$MPIDIST/bin:$PATH
export LD_LIBRARY_PATH PATH
