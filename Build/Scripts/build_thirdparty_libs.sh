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

if [ "$clean_hypre" = true ]; then
    echo "Option --clean-hypre is set."
fi

if [ "$clean_sundials" = true ]; then
    echo "Option --clean-sundials is set."
fi

# FINISHED WITH CLEANING OPTIONS ###########################################

# Decide compilers
echo "Before setting tpt comp"
source ../Scripts/set_thirdparty_compilers.sh

# build hypre
source ../Scripts/HYPRE/build_hypre.sh confmake.sh $clean_hypre

## build sundials
source ../Scripts/SUNDIALS/build_sundials.sh confmake.sh $clean_sundials


# Use ARG and the options
#echo "ARG is: $ARG"
if [ "$SOURCE_INTEL_IFORT" -eq 1 ]; then
   source ../Scripts/set_intel_compiler.sh $ARG
fi   
