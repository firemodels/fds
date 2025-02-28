#!/bin/bash

# Decide compilers
source ../Scripts/set_compilers.sh


# PARSE OPTIONS FOR CLEAN LIBRARY BUILDS ####################################

# Set FIREMODELS environment variable if it is not already exists.
if [ -z "${FIREMODELS}" ]; then
    export FIREMODELS="$(readlink -f "$(pwd)/../../../")"
fi

echo "FIREMODELS=$FIREMODELS"

clean_fds=false
clean_hypre=false
clean_sundials=false
clean_hdf5=false
with_gpu=false
gpu_arch="-"
no_libs=false
ARG=""

# Define an array of allowed GPU architectures
allowed_gpus=("hip" "cuda" "sycl")

# Function to check if a value exists in the allowed list
is_valid_gpu() {
    local value="$1"
    for gpu in "${allowed_gpus[@]}"; do
        if [[ "$gpu" == "$value" ]]; then
            return 0  # Found a match, return success
        fi
    done
    return 1  # No match found, return failure
}


# Loop through the options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --clean-all)
            clean_fds=true
            clean_hypre=true
            clean_sundials=true
            clean_hdf5=true
            shift
            ;;
        --clean-fds)
            clean_fds=true
            shift
            ;;
        --clean-hypre)
            clean_fds=true
            clean_hypre=true  # Set the flag to true when --clean-hypre is used
            shift
            ;;
        --clean-sundials)
            clean_fds=true
            clean_sundials=true
            shift
            ;;
        --clean-hdf5)
            clean_fds=true
            clean_hdf5=true
            shift
            ;;
        --with-gpu=*)
            with_gpu=true
            gpu_arch="${1#*=}"
	    # Check if the value is in the allowed list
            if ! is_valid_gpu "$gpu_arch"; then
                echo "Error: Invalid option for --with-gpu. Allowed values: ${allowed_gpus[*]}." >&2
                exit 1
            fi
	    export BUILD_WITH_GPU=ON
	    export GPU_ARCH="$gpu_arch"
            shift
            ;;
	    --with-gpu)  # Error if --with-gpu is provided without a value
            echo "Error: --with-gpu requires a value (e.g., --with-gpu=hip, --with-gpu=cuda, or --with-gpu=sycl)" >&2
            exit 1
            ;;
	    --no-libs)
            no_libs=true
            clean_fds=true
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
        ARG="${ARG} $1"  # Append unrecognized arguments to ARG
            shift
            ;;
    esac
done

# Trim leading spaces from ARG, if necessary
ARG="${ARG#"${ARG%%[![:space:]]*}"}"

if [ "$clean_fds" = true ]; then
    echo "Option --clean-fds is set."
    rm *.o *.mod >& /dev/null
fi

if [ "$clean_hypre" = true ]; then
    echo "Option --clean-hypre is set."
fi

if [ "$clean_sundials" = true ]; then
    echo "Option --clean-sundials is set."
fi

if [ "$clean_hdf5" = true ]; then
    echo "Option --clean-hdf5 is set."
fi

# FINISHED WITH CLEANING OPTIONS ###########################################


if [ "$no_libs" == false ]; then
   # build hypre
   source ../Scripts/HYPRE/build_hypre.sh confmake.sh $clean_hypre

   # build sundials
   source ../Scripts/SUNDIALS/build_sundials.sh confmake.sh $clean_sundials

   # build hdf5
   source ../Scripts/HDF5/build_hdf5.sh confmake.sh $clean_hdf5
else
   unset SUNDIALS_HOME
   unset HYPRE_HOME
   unset HDF5_HOME
   echo "Building FDS without third-party libraries."
fi


