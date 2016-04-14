#/bin/bash

# MPI distribution location

export MPITYPE=$1
export MPIDIST=$2

if [[ "$MPIDIST" == "" ]]; then
  if [[ "$MPITYPE" == "ib" ]]; then
    if [[ -d /shared/openmpi_64ib ]]; then
      MPIDIST=/shared/openmpi_64ib
    fi
  else
    if [[ -d /shared/openmpi_64 ]]; then
      MPIDIST=/shared/openmpi_64
    fi
  fi
fi

if [[ "$MPIDIST" == "" ]]; then
  echo "*** Warning: the MPI distribution location is not defined."
  echo "Make sure that the MPIDIST_ETH and/or MPIDIST_IB environment"
  echo "variables are defined in your startup file"
  exit
fi

if [[ ! -d $MPIDIST ]]; then
  echo "*** Warning: the MPI distribution, $MPIDIST, does not exist."
  echo "*** Compilation aborted"
  MPIDIST=
  exit
fi

# Update LD_LIBRARY_PATH and PATH
LD_LIBRARY_PATH=$MPIDIST/lib:$LD_LIBRARY_PATH
PATH=$MPIDIST/bin:$PATH
export LD_LIBRARY_PATH PATH
