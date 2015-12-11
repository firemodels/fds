SOURCE_DIR = ../..
FILTERC =
FILTERF =

VPATH = ../../Data_Processing
bin = .

# Definition of the object variables

obj = fds2ascii.o
objwin = fds2ascii.obj

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .F90 .c .obj .o .f90 .cpp

.F90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.obj:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<

# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : FFLAGS    = -O2 /nologo /fpscomp:general
intel_win_64 : FC        = ifort
intel_win_64 : exe       = fds2ascii_win_64.exe

intel_win_64 : $(objwin)
	$(FC) -o $(bin)/$(exe) $(objwin)

# ------------- gnu linux 64 ----------------
#  64 bit gcc linux

gnu_linux_64 : FFLAGS    = -O0 -m64 -ffree-form -frecord-marker=4
gnu_linux_64 : FC        = gfortran
gnu_linux_64 : exe       = fds2ascii_linux_64

gnu_linux_64 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- gnu osx 64 ----------------
#  64 bit gcc linux

gnu_osx_64 : FFLAGS    = -O0 -m64 -ffree-form -frecord-marker=4
gnu_osx_64 : FC        = gfortran
gnu_osx_64 : exe       = fds2ascii_osx_64

gnu_osx_64 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)



# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : FFLAGS    = -O2 -m64
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = fds2ascii_linux_64

intel_linux_64 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : FFLAGS    = -O2 -m64
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = fds2ascii_osx_64

intel_osx_64 : $(obj)
	$(FC) -o $(bin)/$(exe)  -mmacosx-version-min=10.4 $(obj)

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm *.o
