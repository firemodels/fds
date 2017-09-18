#!/bin/bash
network=$1

export setup_fortran=0
if [ "$MPITYPE" == "" ]; then
# CASE 1 - NOT USING MODULES
#          using variables defined in startup script - not using modules
  if [ "$network" == "eth" ]; then
    if [ "$MPIDIST_ETH" == "" ];  then
       MPIDIST_ETH=/shared/openmpi_64
    fi
    if [ ! -d $MPIDIST_ETH ]; then
      echo "***error: The ethernet openmpi library $MPIDIST_ETH does not exist."
      echo "   Define MPIDIST_ETH in a startup file to point to "
      echo "   an existing ethnernet version of openmpi."
      echo "   Example: export MPIDIST_ETH=/shared/openmpi_64"
      export mpi_error=1
      return
    fi
    MPINETWORK=$MPIDIST_ETH
  fi
  if [ "$network" == "ib" ]; then
    if [ "$MPIDIST_IB" == "" ];  then
       MPIDIST_IB=/shared/openmpi_64ib
    fi
    if [ ! -d $MPIDIST_IB ]; then
      echo "***error: The infiniband openmpi library $MPIDIST_IB does not exist."
      echo "   Define MPIDIST_IB in a startup file to point to "
      echo "   an existing infiniband version of openmpi."
      echo "   Example: export MPIDIST_IB/shared/openmpi_64"
      export mpi_error=1
      return
    fi
    MPINETWORK=$MPIDIST_IB
  fi
 
  export not_installed=`ifort -v |& head -1 | grep "not found" | wc -l`

  if [ $not_installed -eq 1 ]; then
    if [ "$IFORT_COMPILER" != "" ]; then
      if [ ! -d $IFORT_COMPILER/bin/compilervars.sh ]; then
         echo "***error: An Intel compiler could not be found at"
         echo "   $IFORT_COMPILER"
         echo "   Define IFORT_COMPILER in a startup file to point to your"
         echo "   compiler location"
         export mpi_error=1
         return
      fi
    fi
    if [ "$IFORT_COMPILER" == "" ]; then
      echo "***error: An Intel compiler could not be found."
      echo "   Define IFORT_COMPILER in a startup file to point to your compiler"
      echo "   location or load a module defining the compiler environmnt."
      export mpi_error=1
      return
    fi
    export setup_fortran=1
  fi

  export LD_LIBRARY_PATH=$MPINETWORK/lib:$LD_LIBRARY_PATH
  export PATH=$MPINETWORK/bin:$PATH
  export MPIDIST=$MPINETWORK
  export MPIFORT=mpifort
  export mpi_error=0
  return
fi

# CASE 2 - USING MODULES
  
export not_installed=`ifort -v |& head -1 | grep "not found" | wc -l`
if [ $not_installed -eq 1 ]; then
  echo "***error: An Intel compiler could not be found."
  echo "   Load a module defining the compiler enviornment or"
  echo "   define the environment in your startup file"
  export mpi_error=1
  return
fi

if [[ "$MPITYPE" != "ib" ]] && [[ "$network" == "ib" ]]; then
  # using infiniband, make sure an infiniband openmpi module is loaded
  if [ "$MPITYPE" == "eth" ]; then
    echo "*** Error: An ethernet openmpi module is loaded but"
    echo "    you are trying to build an infiniband fds."
    echo "    Unload the ethernet openmpi module then load an "
    echo "    infiniband one or build an ethernet version of fds."
  else
    echo "*** Error: An infiniband openmpi module was not detected."
    echo "    If an infiniband openmpi module is loaded, load an infiniband"
    echo "    openmpi module and set MPITYPE to ib in your startup"
    echo "    file or module definition to remove this warning."
  fi
  export mpi_error=1
  return
fi
if [[ "$MPITYPE" != "eth" ]] && [[ "$network" == "eth" ]]; then
  # using ethernet, make sure an ethernet openmpi module is loaded
  if [ "$MPITYPE" == "ib" ]; then
    echo "*** Error: An infiniband openmpi module is loaded but"
    echo "    you are trying to build an ethernet fds."
    echo "    Unload the infiniband openmpi module then load an "
    echo "    ethernet one or build an infiniband version of fds."
  else
    echo "*** Warning: An ethernet openmpi module was not detected."
    echo "    If an ethernet openmpi module is loaded, load an ethernet"
    echo "    openmpi module and set MPITYPE to eth in your startup"
    echo "    file or module definition to remove this warning."
  fi
  export mpi_error=1
  return
fi

# passed all checks
export mpi_error=0
return
