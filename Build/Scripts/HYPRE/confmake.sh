export HYPRE_INSTALL_PREFIX=$FIREMODELS/libs/hypre/$HYPRE_VERSION

if [[ "$FDS_BUILD_TARGET" == *"osx"* ]]; then
   C_FLAGS="-O3"
elif [[ "$FDS_BUILD_TARGET" == *"intel"* ]]; then
   C_FLAGS="-O3 -fno-unsafe-math-optimizations -fp-model=precise"
else
   C_FLAGS="-O3"
fi

cmake_args=(
  -DCMAKE_INSTALL_PREFIX="$HYPRE_INSTALL_PREFIX"
  -DCMAKE_C_COMPILER="$CC"
  -DCMAKE_C_FLAGS="$C_FLAGS"
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

# ./configure CC=$CC FC=mpiifort CFLAGS="-O3 -fno-unsafe-math-optimizations -fp-model=precise" FFLAGS="-O3 -fno-unsafe-math-optimizations -fp-model=precise" \
#       --prefix=$FIREMODELS/libs/hypre/$HYPRE_VERSION

make install
