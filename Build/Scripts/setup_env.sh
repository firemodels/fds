#!/bin/bash
network=$1

if [ "$OPENMPI_MODULE" == "" ]; then
  MPINETWORK=$MPIDIST_ETH
  if [ "network" == "ib" ]; then
    MPINETWORK=$MPIDIST_IB
  fi
  if [ "$IFORT_COMPILER" != "" ]; then
    source $IFORT_COMPILER/bin/compilervars.sh $platform
  fi
  source ../Scripts/set_mpidist.sh $nework $MPINETWORK
  if [ "$MPIDIST" == "" ]; then
# if MPIDIST was not defined above, abort
    echo 1
  fi
  echo 0
fi
if [[ "$OPENMPI_MODULE" == "eth" ]] && [[ "$network" == "ib" ]]; then
  echo 1
fi
if [[ "$OPENMPI_MODULE" == "ib" ]] && [[ "$network" == "eth" ]]; then
  echo 1
fi
echo 0
