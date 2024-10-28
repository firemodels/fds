export SUNDIALS_INSTALL_PREFIX=$FIREMODELS/libs/sundials/$SUNDIALS_VERSION

# Check for Intel C compiler (icx or icc)
if command -v icx &> /dev/null; then
    CC=icx
elif command -v icc &> /dev/null; then
    CC=icc
else
    echo "Error: Neither icx nor icc is available on this system."
    exit 1
fi

# Check for Intel C++ compiler (icpx or icpc)
if command -v icpx &> /dev/null; then
    CXX=icpx
elif command -v icpc &> /dev/null; then
    CXX=icpc
else
    echo "Error: Neither icpx nor icpc is available on this system."
    exit 1
fi

cmake ../ \
-DCMAKE_INSTALL_PREFIX=$SUNDIALS_INSTALL_PREFIX \
-DEXAMPLES_INSTALL_PATH=$SUNDIALS_INSTALL_PREFIX/examples \
-DCMAKE_C_COMPILER=$CC \
-DCMAKE_CXX_COMPILER=$CXX \
-DCMAKE_Fortran_COMPILER=ifort \
-DBUILD_FORTRAN_MODULE_INTERFACE=ON \
-DEXAMPLES_ENABLE_CXX=OFF \
-DEXAMPLES_ENABLE_CUDA=OFF \
-DEXAMPLES_ENABLE_F2003=OFF \
-DENABLE_OPENMP=ON

make install
