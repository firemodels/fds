#!/bin/bash

# PARSE OPTIONS FOR CLEAN LIBRARY BUILDS ####################################

# Parse the long options first using getopt
OPTIONS=$(getopt -o "" --long clean-hypre,clean-sundials -- "$@")

# Check if getopt parsed successfully
if [ $? -ne 0 ]; then
    echo "Error parsing options."
    exit 1
fi

# Evaluate the parsed options
eval set -- "$OPTIONS"

# Initialize variables for options
clean_hypre=false
clean_sundials=false
ARG=""

# Loop through the options
while true; do
    case "$1" in
        --clean-hypre)
            clean_hypre=true  # Set the flag to true when --clean-hypre is used
            shift
            ;;
        --clean-sundials)
            clean_sundials=true
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Invalid option."
            exit 1
            ;;
    esac
done

# After all options are processed, check for any remaining positional argument (ARG)
if [ -n "$1" ]; then
    ARG=$1
fi

# Use ARG and the options
#echo "ARG is: $ARG"

if [ "$clean_hypre" = true ]; then
    echo "Option --clean-hypre is set."
fi

if [ "$clean_sundials" = true ]; then
    echo "Option --clean-sundials is set."
fi

# FINISHED WITH CLEANING OPTIONS ###########################################

source ../Scripts/set_intel_compiler.sh $ARG

dir=`pwd`
target=${dir##*/}

# build hypre
source ../Scripts/HYPRE/build_hypre.sh confmake_impi_intel_linux.sh $clean_hypre

## build sundials
source ../Scripts/SUNDIALS/build_sundials.sh confmake_impi_intel_linux.sh $clean_sundials

# build fds
echo Building $target with Intel MPI and $INTEL_IFORT
make VPATH="../../Source" -f ../makefile $target
