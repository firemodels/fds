FILTERC =
FILTERF =

VPATH = ../source
bin = .

# Definition of the object variables

obj = interface.o fds2ast.o fds2lsdyna.o
objgnu = interface.o fds2ast.o fds2lsdyna.o
objwin = interface.obj fds2ast.obj fds2lsdyna.obj

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


# ------------- Win 64 ----------------
#  64 bit windows

intel_win : FFLAGS    = -O3 /fp:source /CB
intel_win : FC        = ifort
intel_win : obj       = fds2ftmi_win

intel_win : setup_win $(objwin)
	$(FC) -o $(obj) $(FFLAGS) $(objwin)

# ----------- debug Win 64 --------------
#  64 bit windows

intel_win_db : FFLAGS    = /CB -auto -WB -traceback -fltconsistency -stand:f08  
intel_win_db : FC        = ifort
intel_win_db : obj       = fds2ftmi_win_db

intel_win_db : setup_win $(objwin)
	$(FC) -o $(obj) $(FFLAGS) $(objwin)

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux : FFLAGS    = -O3 -check bounds -fp-model source 
intel_linux : FC        = ifort
intel_linux : exe       = fds2ftmi_linux

intel_linux : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- debug Linux 64 ---------
#  64 bit linux

intel_linux_db : FFLAGS    = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fltconsistency -stand:f08
intel_linux_db : FC        = ifort
intel_linux_db : exe       = fds2ftmi_linux_db

intel_linux_db : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx : FFLAGS    = -O3 -check bounds -fp-model source  
intel_osx : FC        = ifort
intel_osx : exe       = fds2ftmi_osx

intel_osx : $(obj)
	$(FC) -o $(bin)/$(exe)  -mmacosx-version-min=10.4 $(obj)

# -------------GFORTRAN---------------
# ------------- Win 64 ----------------
#  64 bit GNU windows

gnu_win : FFLAGS    = -O3 -fbounds-check -fno-range-check -ffree-form -frecord-marker=4
gnu_win : FC        = gfortran
gnu_win : obj       = fds2ftmi_gnu_win
gnu_win : setup_win $(objgnu)
	$(FC) -o $(obj) $(FFLAGS) $(objgnu)	

#*** End Compiler Specific Rules ***

# *** Object Dependencies ***

interface.o : fds2ast.o fds2lsdyna.o	
fds2lsdyna.o : fds2ast.o	

interface.obj : fds2ast.obj fds2lsdyna.obj		
fds2lsdyna.obj : fds2ast.obj

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	-rm *.o

.PHONY : winclean
winclean:
	-erase *.mod *.obj *o
