# Makefile for test_mpi

# A make utility for windows may be downloaded from http://gnuwin32.sourceforge.net/packages/make.htm

# The variable VPATH defines where the source code is relative to the current directory

VPATH = ..

# Definition of the non-MPI object variables

obj_serial = test_mpi.o

objwin_serial = $(obj_serial:.o=.obj)

# Definition of the MPI object variables

obj_mpi = test_mpi.o

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

# windows

impi_intel_win: SHELL=C:/Windows/System32/cmd.exe
impi_intel_win: MPILIB = "$(I_MPI_RELEASE_ROOT)\impi.lib"
impi_intel_win: MPIINCLUDE = "$(IMPI_INCLUDE)"
impi_intel_win: FFLAGS = /Qipo /O2 /I$(MPIINCLUDE)
impi_intel_win: FCOMPL = mpiifx
impi_intel_win: obj = test_mpi
impi_intel_win: setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F1000000000 $(objwin_mpi) $(MPILIB)

# Linux

impi_intel_linux: FFLAGS = -m64 -O2 -traceback -no-wrap-margin
impi_intel_linux: FCOMPL = mpiifx
impi_intel_linux: obj = test_mpi
impi_intel_linux: $(obj_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(obj_mpi)

# OSX

ompi_intel_osx: FFLAGS = -m64 -O2 -traceback -no-wrap-margin
ompi_intel_osx: LFLAGS = -static-intel -ld_classic
ompi_intel_osx: FCOMPL = mpifort
ompi_intel_osx: obj = test_mpi
ompi_intel_osx: $(obj_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(LFLAGS) $(obj_mpi)

ompi_gnu_osx: FFLAGS = -m64 -O2
ompi_gnu_osx: LFLAGS =
ompi_gnu_osx: FCOMPL = mpifort
ompi_gnu_osx: obj = test_mpi
ompi_gnu_osx: $(obj_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(LFLAGS) $(obj_mpi)

#*** Clean Target to remove Object and Module files ***

.PHONY : clean
clean:
	-rm -f *.o *.mod *.obj 

.PHONY : winclean
winclean:
	-erase *.mod *.obj 
