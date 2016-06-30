# Makefile for test_mpi_openmp

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

openmp_impi_intel_win: MPILIB = "C:\Program Files (x86)\Intel\MPI\5.0.0.028\intel64\lib\impi.lib"
openmp_impi_intel_win: MPIINCLUDE = "C:\Program Files (x86)\Intel\MPI\5.0.0.028\intel64\include"
openmp_impi_intel_win: FFLAGS = /Qipo /O2 /openmp /I$(MPIINCLUDE) -mt_mpi
openmp_impi_intel_win: FCOMPL = mpiifort
openmp_impi_intel_win: obj = test_mpi_openmp
openmp_impi_intel_win: setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F1000000000 $(objwin_mpi) $(MPILIB)

openmp_impi_intel_linux : FFLAGS = -m64 -O2 -openmp
openmp_impi_intel_linux : LFLAGS = -static-intel
openmp_impi_intel_linux : FCOMPL = mpiifort
openmp_impi_intel_linux : FOPENMPFLAGS = -openmp-link static -liomp5 
openmp_impi_intel_linux : obj = test_mpi_openmp
openmp_impi_intel_linux : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi)

openmp_mpi_intel_linux : FFLAGS = -m64 -O2 -openmp
openmp_mpi_intel_linux : LFLAGS = -static-intel
openmp_mpi_intel_linux : FCOMPL = mpifort
openmp_mpi_intel_linux : FOPENMPFLAGS = -openmp-link static -liomp5
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
