echo "FDS build target = $FDS_BUILD_TARGET"

# For following variables 1- indicate availble through 
# environment variable FIREMODELS_LIBS_CC, FIREMODELS_LIBS_CXX, and FIREMODELS_LIBS_FC
set_CC=0 
set_CXX=0
set_FC=0

if [ -n "$FIREMODELS_LIBS_CC" ]; then
   CC=$FIREMODELS_LIBS_CC
   if command -v $CC &> /dev/null; then
      set_CC=1
   else
      echo "The compiler specified by the FIREMODELS_LIBS_CC environment variable ($CC) is not available on this system. Searching for an alternative compiler."
   fi
fi

if [ -n "$FIREMODELS_LIBS_CXX" ]; then
   CXX=$FIREMODELS_LIBS_CXX
   if command -v $CXX &> /dev/null; then
      set_CXX=1
   else
      echo "The compiler specified by the FIREMODELS_LIBS_CXX environment variable ($CXX) is not available on this system. Searching for an alternative compiler."
   fi
fi

if [ -n "$FIREMODELS_LIBS_FC" ]; then
   FC=$FIREMODELS_LIBS_FC
   if command -v $FC &> /dev/null; then
      set_FC=1
   else
      echo "The compiler specified by the FIREMODELS_LIBS_FC environment variable ($FC) is not available on this system. Searching for an alternative compiler."
   fi
fi


if [[ "$FDS_BUILD_TARGET" == *"osx"* ]]; then
   # Check for C compiler (mpicc, gcc, clang)
   if [ $set_CC -eq 0 ]; then
      if command -v mpicc &> /dev/null; then
         CC=mpicc
      elif command -v clang &> /dev/null; then
         CC=clang
      elif command -v gcc &> /dev/null; then
         CC=gcc
      else
         echo "Error: Any of mpicc, clang, or gcc compiler is not available on this system."
         exit 1
      fi
   fi

   # Check for clang C++ compiler (mpicxx, g++, clang ++)
   if [ $set_CXX -eq 0 ]; then
      if command -v mpicxx &> /dev/null; then
         CXX=mpicxx
      elif command -v clang++ &> /dev/null; then
         CXX=clang++
      elif command -v g++ &> /dev/null; then
         CXX=g++
      else
         echo "Error: Any of mpicxx, clang++, or g++ compiler is not available on this system."
         exit 1
      fi
   fi

   # Check for Fortran compiler (gfortran)
   if [ $set_FC -eq 0 ]; then
      if command -v mpifort &> /dev/null; then
         FC=mpifort
      elif command -v gfortran &> /dev/null; then
         FC=gfortran
      else
         echo "Error: Any of mpifort or gfortran compiler is not available on this system."
         exit 1
      fi
   fi

elif [[ "$FDS_BUILD_TARGET" == *"intel"* ]]; then
   # Check for Intel C compiler
   if [ $set_CC -eq 0 ]; then
      if command -v mpiicx &> /dev/null; then
         CC=mpiicx
      elif command -v icx &> /dev/null; then
         CC=icx
      elif command -v mpiicc &> /dev/null; then
         CC=mpiicc
      elif command -v icc &> /dev/null; then
         CC=icc
      else
         echo "Error: Any of mpiicx, icx, mpiicc, or icc compiler is not available on this system."
         exit 1
      fi
   fi

   # Check for Intel C++ compiler
   if [ $set_CXX -eq 0 ]; then
      if command -v mpiicpx &> /dev/null; then
         CXX=mpiicpx
      elif command -v icpx &> /dev/null; then
         CXX=icpx
      elif command -v mpiicpc &> /dev/null; then
         CXX=mpiicpc
      elif command -v icpc &> /dev/null; then
         CXX=icpc
      else
         echo "Error: Any of mpiicpx, icpx, mpiicpc, or icpc compiler is not available on this system."
         exit 1
      fi
   fi

   # Check for Intel Fortran compiler (ifort). ifx will be added in future.
   if [ $set_FC -eq 0 ]; then
      if command -v mpiifort &> /dev/null; then
         FC=mpiifort
      elif command -v ifort &> /dev/null; then
         FC=ifort
      else
         echo "Error: Any of mpiifort, or ifort compiler is not available on this system."
         exit 1
      fi
   fi
else #gnu
   # Check for gnu C compiler (gcc)
   if [ $set_CC -eq 0 ]; then
      if command -v mpicc &> /dev/null; then
         CC=mpicc
      elif command -v gcc &> /dev/null; then
         CC=gcc
      else
         echo "Error: Any of mpicc, gcc compiler is not available on this system."
         exit 1
      fi
   fi

   # Check for Intel C++ compiler (g++)
   if [ $set_CXX -eq 0 ]; then
      if command -v mpicxx &> /dev/null; then
         CXX=mpicxx
      elif command -v g++  &> /dev/null; then
         CXX=g++
      else
         echo "Error: Any of mpicxx, g++ compiler is not available on this system."
         exit 1
      fi
   fi

   # Check for Fortran compiler (gfortran)
   if [ $set_FC -eq 0 ]; then
      if command -v mpifort &> /dev/null; then
         FC=mpifort
      elif command -v gfortran &> /dev/null; then
         FC=gfortran
      else
         echo "Error: Any of mpifort, gfortran compiler is not available on this system."
         exit 1
      fi
   fi
fi

echo "Third-party libs C Compiler=$CC"
echo "Third-party libs C++ compiler=$CXX"
echo "Third-party libs Fortran compiler=$FC"

export CC=$CC
export CXX=$CXX
export FC=$FC
