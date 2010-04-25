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


# ------------- Win 32 ----------------
#  32 bit windows

intel_win_32 : FFLAGS    = -O2 /Qvec_report0 /nologo 
intel_win_32 : FC        = ifort
intel_win_32 : exe       = fds2ascii_win_32.exe

intel_win_32 : $(objwin)
	$(FC) -o $(bin)/$(exe) $(objwin)


# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : FFLAGS    = -O2 /Qvec_report0 /nologo 
intel_win_64 : FC        = ifort
intel_win_64 : exe       = fds2ascii_win_64.exe

intel_win_64 : $(objwin)
	$(FC) -o $(bin)/$(exe) $(objwin)

# ------------- Linux 32 ----------------
#  32 bit linux

intel_linux_32 : FFLAGS    = -O2 -m32 -vec_report0 
intel_linux_32 : FC        = ifort
intel_linux_32 : exe       = fds2ascii_linux_32

intel_linux_32 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : FFLAGS    = -O2 -m64 -vec_report0 
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = fds2ascii_linux_64

intel_linux_64 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- OSX 32 ----------------
#  32 bit OSX

intel_osx_32 : FFLAGS    = -O2 -m32 -vec_report0 
intel_osx_32 : FC        = ifort
intel_osx_32 : exe       = fds2ascii_osx_32

intel_osx_32 : $(obj)
	$(FC) -o $(bin)/$(exe)  -mmacosx-version-min=10.4 $(obj)

# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : FFLAGS    = -O2 -m64 -vec_report0 
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = fds2ascii_osx_64

intel_osx_64 : $(obj)
	$(FC) -o $(bin)/$(exe)  -mmacosx-version-min=10.4 $(obj)

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm *.o
