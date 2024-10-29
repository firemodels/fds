export HYPRE_INSTALL_PREFIX=$FIREMODELS/libs/hypre/$HYPRE_VERSION

if [[ "$FDS_BUILD_TARGET" == *"intel"* ]]; then
   C_FLAGS="-O3 -fno-unsafe-math-optimizations -fp-model=precise"
else
   C_FLAGS="-O3"
fi

cmake ../ \
-DCMAKE_INSTALL_PREFIX=$HYPRE_INSTALL_PREFIX \
-DCMAKE_C_COMPILER=$CC \
-DCMAKE_C_FLAGS=$C_FLAGS \
-DCMAKE_INSTALL_LIBDIR="lib" 
#-DCMAKE_OSX_DEPLOYMENT_TARGET="14.0"

# ./configure CC=$CC FC=mpiifort CFLAGS="-O3 -fno-unsafe-math-optimizations -fp-model=precise" FFLAGS="-O3 -fno-unsafe-math-optimizations -fp-model=precise" \
#       --prefix=$FIREMODELS/libs/hypre/$HYPRE_VERSION

make install
