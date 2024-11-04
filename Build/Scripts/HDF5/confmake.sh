export HDF5_INSTALL_PREFIX=$FIREMODELS/libs/hdf5/$HDF5_VERSION

../configure  --enable-fortran --enable-parallel CXX="$(which mpiicpc) -cc=$(which icpx)" CC="$(which mpiicc) -cc=$(which icx)" FC="$(which mpiifort) -fc=$(which ifort)" --prefix=$HDF5_INSTALL_PREFIX

make install
