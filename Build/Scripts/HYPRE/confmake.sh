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

make install
