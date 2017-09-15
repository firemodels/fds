#!/bin/bash
network=$1

if [ "$OPENMPI_MODULE" == "" ]; then
# not using modules
  if [ "$network" == "eth" ]; then
# using ethernet openmpi make sure MPIDIST_ETH is defined
    if [ "$MPIDIST_ETH" != "" ]; then
      MPINETWORK=$MPIDIST_ETH
    else
      echo "***warning: trying to build an ethernet fds but"
      echo "   MPIDIST_ETH is not defined."
      echo "   Define MPIDIST_ETH to point to an ethnernet"
      echo "   version of openmpi"
      return 1
    fi
  fi
  if [ "$network" == "ib" ]; then
# using infiniband openmpi make sure MPIDIST_IB is defined
    if [ "$MPIDIST_IB" != "" ]; then
      MPINETWORK=$MPIDIST_IB
    else
      echo "***warning: trying to build an infiniband fds but"
      echo "   MPIDIST_IB is not defined."
      echo "   Define MPIDIST_IB to point to an infiniband"
      echo "   version of openmpi"
      return 1
    fi
  fi
# openmpi is defined now define compiler environment
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

# using modules

if [[ "$OPENMPI_MODULE" != "ib" ]] && [[ "$network" == "ib" ]]; then
  # using infiniband, make sure an infiniband openmpi module is loaded
  if [ "$OPENMPI_MODULE" == "eth" ]; then
    echo "*** Error: An ethernet openmpi module is loaded."
    echo "    But you are trying to build and infiniband fds."
  else
    echo "*** Warning: An infiniband openmpi module was not detected."
    echo "    If an infiniband openmpi module is loaded, load an infiniband"
    echo "    openmpi module and set OPENMPI_MODULE to ib in your startup"
    echo "    file or module definition to remove this warning."
  fi
  return 1
fi
if [[ "$OPENMPI_MODULE" != "eth" ]] && [[ "$network" == "eth" ]]; then
  # using ethernet, make sure an ethernet openmpi module is loaded
  if [ "$OPENMPI_MODULE" == "ib" ]; then
    echo "*** Error: An infiniband openmpi module is loaded but"
    echo "    you are trying to build and ethernet fds."
  else
    echo "*** Warning: An ethernet openmpi module was not detected."
    echo "    If an ethernet openmpi module is loaded, load an ethernet"
    echo "    openmpi module and set OPENMPI_MODULE to eth in your startup"
    echo "    file or module definition to remove this warning."
  fi
  return 1
fi
# passed all checks
return 0
