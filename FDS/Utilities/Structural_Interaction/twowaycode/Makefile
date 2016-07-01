FILTERC =
FILTERF =

VPATH = ../source
bin = .

# Definition of the object variables

obj = 2waycode.o interface.o be2ftmi.o write_gc.o
objwin = 2waycode.obj interface.obj be2ftmi.obj write_gc.obj

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

setup_win:
%.obj : %.mod		

.SUFFIXES: .f90 .o .obj

.f90.obj:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<


# ------------- Win 32 ----------------
#  32 bit windows

intel_win_32 : FFLAGS    = -O3 /fp:source /check:bounds /fpp /DSYS=2
intel_win_32 : FC        = ifort
intel_win_32 : exe       = twowaycode_win_32.exe

intel_win_32 : setup_win $(objwin)
	$(FC) -o $(exe) $(FFLAGS) $(objwin)


# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : FFLAGS    = -O3 /fp:source /CB /fpp /DSYS=2
intel_win_64 : FC        = ifort
intel_win_64 : exe       = twowaycode_win_64.exe

intel_win_64 : setup_win $(objwin)
	$(FC) -o $(exe) $(FFLAGS) $(objwin)

# ------------- Linux 32 ----------------
#  32 bit linux

intel_linux_32 : FFLAGS    = -O3 -check bounds -fp-model source -fpp -DSYS=1
intel_linux_32 : FC        = ifort
intel_linux_32 : exe       = twowaycode_linux_32

intel_linux_32 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : FFLAGS    = -O3 -check bounds -fp-model source -fpp -DSYS=1
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = twowaycode_linux_64

intel_linux_64 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- debug Linux 64 ---------
#  64 bit linux

intel_linux_64_db : FFLAGS    = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fltconsistency -stand:f08 -fpp -DSYS=1
intel_linux_64_db : FC        = ifort
intel_linux_64_db : exe       = twowaycode_linux_64_db

intel_linux_64_db : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- OSX 32 ----------------
#  32 bit OSX

intel_osx_32 : FFLAGS    = -O3 -check bounds -fp-model source -fpp -DSYS=1
intel_osx_32 : FC        = ifort
intel_osx_32 : exe       = twowaycode_osx_32

intel_osx_32 : $(obj)
	$(FC) -o $(bin)/$(exe)  -mmacosx-version-min=10.4 $(obj)

# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : FFLAGS    = -O3 -check bounds -fp-model source -fpp -DSYS=1 
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = twowaycode_osx_64

intel_osx_64 : $(obj)
	$(FC) -o $(bin)/$(exe)  -mmacosx-version-min=10.4 $(obj)

#*** End Compiler Specific Rules ***

# *** Object Dependencies ***

2waycode.o : interface.o be2ftmi.o write_gc.o
interface.o : be2ftmi.o
write_gc.o : be2ftmi.o

2waycode.obj : interface.obj be2ftmi.obj write_gc.obj
interface.obj : be2ftmi.obj
write_gc.obj : be2ftmi.obj

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	-rm *.o

.PHONY : winclean
winclean:
	-erase *.mod *.obj *o