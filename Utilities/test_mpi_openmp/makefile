# Makefile for test_mpi_openmp

# A make utility for windows may be downloaded from http://gnuwin32.sourceforge.net/packages/make.htm

# The variable VPATH defines where the source code is relative to the current directory

VPATH = cd

# Definition of the non-MPI object variables

obj_serial = test_mpi_openmp.o

objwin_serial = $(obj_serial:.o=.obj)

# Definition of the MPI object variables

obj_mpi = test_mpi_openmp.o

objwin_mpi = $(obj_mpi:.o=.obj)

# General Purpose Rules

no_target:
	@echo \******** You did not specify a make target \********
	@echo Please read the comments at the top of the makefile

setup:
%.o : %.mod

setup_win:
%.obj : %.mod

.SUFFIXES: .f90 .o .obj

.f90.o:
	$(FCOMPL) -c $(FFLAGS) $<
.f90.obj:
	$(FCOMPL) -c $(FFLAGS) $<

#Compiler Rules

openmp_mpi_intel_win: MPILIB = "C:\Program Files (x86)\Intel\MPI\5.0.0.028\intel64\lib\impi.lib"
openmp_mpi_intel_win: MPIINCLUDE = "C:\Program Files (x86)\Intel\MPI\5.0.0.028\intel64\include"
openmp_mpi_intel_win: FFLAGS = /Qipo /O2 /Qvec_report0 /openmp /I$(MPIINCLUDE) 
openmp_mpi_intel_win: FCOMPL = mpiifort
openmp_mpi_intel_win: obj = test_mpi_openmp
openmp_mpi_intel_win: setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F1000000000 $(objwin_mpi) $(MPILIB)

openmp_mpi_intel_linux : FFLAGS = -m64 -O2 -vec_report0 
openmp_mpi_intel_linux : LFLAGS = -static-intel
openmp_mpi_intel_linux : FCOMPL = mpif90
openmp_mpi_intel_linux : FOPENMPFLAGS = -openmp -openmp-link static -liomp5 
openmp_mpi_intel_linux : obj = test_mpi_openmp
openmp_mpi_intel_linux : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi)

#*** Clean Target to remove Object and Module files ***

.PHONY : clean
clean:
	-rm -f *.o *.mod *.obj 

.PHONY : winclean
winclean:
	-erase *.mod *.obj 
