#/bin/bash

# MPI distribution location

export MPITYPE=$1
export MPIDIST=$2

if [ "$MPIFORT" == "" ]; then
  MPIFORT=mpifort
fi
export MPIFORT

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
  if [[ "$MPITYPE" == "ib" ]]; then
    echo "*** Warning: the infiniband MPI library location is not defined."
    echo "Make sure that MPIDIST_IB is defined in your startup file."
  else
    echo "*** Warning: the ethernet MPI library location is not defined."
    echo "Make sure that MPIDIST_ETH is defined in your startup file."
  fi
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
