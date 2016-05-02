# Makefile for test_mpi_write

# A make utility for windows may be downloaded from http://gnuwin32.sourceforge.net/packages/make.htm

# The variable VPATH defines where the source code is relative to the current directory

VPATH = cd

# Definition of the non-MPI object variables

obj_serial = test_mpi_write.o

objwin_serial = $(obj_serial:.o=.obj)

# Definition of the MPI object variables

obj_mpi = test_mpi_write.o

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

impi_intel_win: MPILIB = "$(I_MPI_ROOT)\intel64\lib\impi.lib"
impi_intel_win: MPIINCLUDE = "$(I_MPI_ROOT)\intel64\include"
impi_intel_win: FFLAGS = /Qipo /O2 /traceback /I$(MPIINCLUDE)
impi_intel_win: FCOMPL = ifort
impi_intel_win: obj = test_mpi_write
impi_intel_win: setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F1000000000 $(objwin_mpi) $(MPILIB)

	
	
#*** Clean Target to remove Object and Module files ***

.PHONY : clean
clean:
	-rm -f *.o *.mod *.obj 

.PHONY : winclean
winclean:
	-erase *.mod *.obj 
