export HDF5_INSTALL_PREFIX=$FIREMODELS/libs/hdf5/$HDF5_VERSION

cmake_args=(
  -DCMAKE_INSTALL_PREFIX="$HDF5_INSTALL_PREFIX"
  -DCMAKE_C_COMPILER="$COMP_CC"
  -DCMAKE_CXX_COMPILER="$COMP_CXX"
  -DCMAKE_Fortran_COMPILER="$COMP_FC"
  -DCMAKE_BUILD_TYPE="Release"
  -DBUILD_SHARED_LIBS="OFF"
  -DBUILD_TESTING="OFF"
  -DHDF5_BUILD_TOOLS="OFF"
  -DHDF5_BUILD_FORTRAN="ON"
  -DHDF5_ENABLE_PARALLEL="ON"
)


# Add OSX deployment target if building for macOS
if [[ "$FDS_BUILD_TARGET" == *"osx"* ]]; then
   if [ "$(uname -m)" == "x86_64" ]; then
      cmake_args+=(-DCMAKE_OSX_DEPLOYMENT_TARGET="10.6")
   else
      cmake_args+=(-DCMAKE_OSX_DEPLOYMENT_TARGET="13.0")
   fi
fi

# Run cmake with the arguments
cmake ../ "${cmake_args[@]}"

make install

