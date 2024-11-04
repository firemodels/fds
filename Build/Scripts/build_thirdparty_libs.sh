# PARSE OPTIONS FOR CLEAN LIBRARY BUILDS ####################################

# Set FIREMODELS environment variable if it is not already exists.
if [ -z "${FIREMODELS}" ]; then
    export FIREMODELS="$(readlink -f "$(pwd)/../../../")"
fi 

echo "FIREMODELS=$FIREMODELS"

clean_hypre=false
clean_sundials=false
clean_hdf5=false
ARG=""

# Loop through the options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --clean-all)
            clean_hypre=true
            clean_sundials=true
            shift
            ;;
        --clean-hypre)
            clean_hypre=true  # Set the flag to true when --clean-hypre is used
            shift
            ;;
        --clean-sundials)
            clean_sundials=true
            shift
            ;;
        --clean-hdf5)
            clean_hdf5=true
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

# Decide compilers
source ../Scripts/set_thirdparty_compilers.sh

# build hypre
source ../Scripts/HYPRE/build_hypre.sh confmake.sh $clean_hypre

## build sundials
source ../Scripts/SUNDIALS/build_sundials.sh confmake.sh $clean_sundials

## build hdf5
source ../Scripts/HDF5/build_hdf5.sh confmake.sh $clean_hdf5

