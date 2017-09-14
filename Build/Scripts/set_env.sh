#!/bin/bash
network=$1

# setup environment on systems without modules

if [ "$OPENMPI_MODULE" == "" ]; then
  MPINETWORK=$MPIDIST_ETH
  if [ "network" == "ib" ]; then
    MPINETWORK=$MPIDIST_IB
  fi
  if [ "$IFORT_COMPILER" != "" ]; then
    source $IFORT_COMPILER/bin/compilervars.sh intel64
  fi
  SCRIPTPATH=$( cd $(dirname $0) ; pwd -P )
  echo SCRIPTPATH=$SCRIPTPATH
  source $SCRIPTPATH/../Scripts/set_mpidist.sh $network $MPINETWORK
  if [ "$MPIDIST" == "" ]; then
# if MPIDIST was not defined above, abort
  echo "*** Warning: unable to define an openmpi environment"
    return 1
  fi
  return 0
fi

# make sure the right type of module is loaded (infiniband openmpi 
# for infiniband fds builds, ethernet openmpi for ethernet fds builds)

if [[ "$OPENMPI_MODULE" != "ib" ]] && [[ "$network" == "ib" ]]; then
  echo "*** Warning: An infiniband openmpi module was not detected."
  echo "             Set OPENMPI_MODULE to ib in your startup file"
  echo "             or module definition to remove this warning."
  return 1
fi
if [[ "$OPENMPI_MODULE" != "eth" ]] && [[ "$network" == "eth" ]]; then
  echo "*** Warning: An ethernet openmpi module was not detected."
  echo "             Set OPENMPI_MODULE to eth in your startup file"
  echo "             or module definition to remove this warning."
  return 1
fi
return 0
