./configure CC=mpiicc FC=mpiifort CFLAGS="-O3 -fno-unsafe-math-optimizations -fp-model=precise" FFLAGS="-O3 -fno-unsafe-math-optimizations -fp-model=precise" \
            --prefix=$FIREMODELS/libs/hypre/$HYPRE_VERSION
make install
