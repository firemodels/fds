# Check for mpiicx or mpiicc
if command -v mpiicx &> /dev/null; then
    CC=mpiicx
elif command -v mpiicc &> /dev/null; then
    CC=mpiicc
else
    echo "Error: Neither mpiicx nor mpiicc is available on this system."
    exit 1
fi

./configure CC=$CC FC=mpiifort CFLAGS="-O3 -fno-unsafe-math-optimizations -fp-model=precise" FFLAGS="-O3 -fno-unsafe-math-optimizations -fp-model=precise" \
            --prefix=$FIREMODELS/libs/hypre/$HYPRE_VERSION
make install
