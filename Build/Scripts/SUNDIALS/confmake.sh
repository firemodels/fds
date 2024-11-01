export SUNDIALS_INSTALL_PREFIX=$FIREMODELS/libs/sundials/$SUNDIALS_VERSION

cmake_args=(
  -DCMAKE_INSTALL_PREFIX="$SUNDIALS_INSTALL_PREFIX"
  -DEXAMPLES_INSTALL_PATH="$SUNDIALS_INSTALL_PREFIX/examples"
  -DCMAKE_C_COMPILER="$CC"
  -DCMAKE_CXX_COMPILER="$CXX"
  -DCMAKE_Fortran_COMPILER="$FC"
  -DBUILD_FORTRAN_MODULE_INTERFACE=ON
  -DEXAMPLES_ENABLE_CXX=OFF
  -DEXAMPLES_ENABLE_CUDA=OFF
  -DEXAMPLES_ENABLE_F2003=OFF
  -DENABLE_OPENMP=OFF
  -DCMAKE_INSTALL_LIBDIR="lib"
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
