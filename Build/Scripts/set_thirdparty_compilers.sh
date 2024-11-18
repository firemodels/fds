echo "FDS build target = $FDS_BUILD_TARGET"

# Initialize variables and check environment variables
set_compiler_var() {
    local var_name=$1    # Variable name to set (e.g., CC, CXX, FC)
    local env_var=$2     # Environment variable to check (e.g., FIREMODELS_LIBS_CC)
    local set_flag_var=$3 # Flag variable name (e.g., set_CC, set_CXX, set_FC)

    if [ -n "${!env_var}" ]; then
        eval "$var_name=${!env_var}"
        if command -v "${!var_name}" &> /dev/null; then
            eval "$set_flag_var=1"
        else
            echo "Warning: ${!env_var} specified by $env_var is not available. Searching for an alternative."
        fi
    fi
}


# Set compilers based on the build target
select_compiler() {
    local var_name=$1       # Variable to set (CC, CXX, FC)
    shift
    local compilers=("$@")  # List of compilers to check in order
    local set_flag_var="set_$var_name"

    # Only proceed if the compiler flag is not set
    if [ "$(eval echo \$$set_flag_var)" -eq 0 ]; then
        for compiler in "${compilers[@]}"; do
            if command -v "$compiler" &> /dev/null; then
                # Set the compiler variable
                eval "$var_name=$compiler"
                # Set the flag variable to 1 (indicating the compiler was found)
                eval "$set_flag_var=1"
                return
            fi
        done
        echo "Error: None of the specified compilers (${compilers[*]}) are available for $var_name on this system."
        exit 1
    fi
}


# Following variables indicate if compilers are set using environment variables FIREMODELS_LIBS_XXX.
set_CC=0
set_CXX=0
set_FC=0

# Check environment variables for compilers
set_compiler_var CC FIREMODELS_LIBS_CC set_CC
set_compiler_var CXX FIREMODELS_LIBS_CXX set_CXX
set_compiler_var FC FIREMODELS_LIBS_FC set_FC

# Determine compiler list based on build target
if [[ "$FDS_BUILD_TARGET" == *"osx"* ]]; then
    select_compiler CC mpicc clang gcc
    select_compiler CXX mpicxx clang++ g++
    select_compiler FC mpifort gfortran
elif [[ "$FDS_BUILD_TARGET" == *"intel"* ]]; then
    select_compiler CC mpiicx icx mpiicc icc
    select_compiler CXX mpiicpx icpx mpiicpc icpc
    select_compiler FC ifort mpiifx ifx
else  # Default to GNU compilers
    select_compiler CC mpicc gcc
    select_compiler CXX mpicxx g++
    select_compiler FC mpifort gfortran
fi


echo "Third-party libs C Compiler=$CC"
echo "Third-party libs C++ compiler=$CXX"
echo "Third-party libs Fortran compiler=$FC"

export CC=$CC
export CXX=$CXX
export FC=$FC



