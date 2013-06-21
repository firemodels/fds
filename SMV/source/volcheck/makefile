FILTERF =

VPATH = .:../smokeview
bin = .

# Definition of the object variables

obj = getdata.o volcheck.o
objwin = getdata.obj volcheck.obj

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

intel_win_64 : FFLAGS    = -O2 /Qvec_report0 /nologo /fpp
intel_win_64 : FC        = ifort
intel_win_64 : exe       = volcheck.exe

intel_win_64 : $(objwin)
	$(FC) -o $(bin)/$(exe) $(objwin)

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = volcheck_linux_64

intel_linux_64 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = volcheck_osx_64

intel_osx_64 : $(obj)
	$(FC) -o $(bin)/$(exe)  -mmacosx-version-min=10.4 $(obj)

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm *.o
