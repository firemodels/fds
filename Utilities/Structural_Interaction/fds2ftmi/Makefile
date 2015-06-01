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


# ------------- Win 32 ----------------
#  32 bit windows

intel_win_32 : FFLAGS    = -O3 /fp:source /CB
intel_win_32 : FC        = ifort
intel_win_32 : obj       = fds2ftmi_win_32
intel_win_32 : setup_win $(objwin)
	$(FC) -o $(obj) $(FFLAGS) $(objwin)


# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : FFLAGS    = -O3 /fp:source /CB
intel_win_64 : FC        = ifort
intel_win_64 : obj       = fds2ftmi_win_64

intel_win_64 : setup_win $(objwin)
	$(FC) -o $(obj) $(FFLAGS) $(objwin)

# ----------- debug Win 64 --------------
#  64 bit windows

intel_win_64_db : FFLAGS    = /CB -auto -WB -traceback -fltconsistency -stand:f08
intel_win_64_db : FC        = ifort
intel_win_64_db : obj       = fds2ftmi_win_64_db

intel_win_64_db : setup_win $(objwin)
	$(FC) -o $(obj) $(FFLAGS) $(objwin)

# ------------- Linux 32 ----------------
#  32 bit linux

intel_linux_32 : FFLAGS    = -O3 -check bounds -fp-model source
intel_linux_32 : FC        = ifort
intel_linux_32 : exe       = fds2ftmi_linux_32

intel_linux_32 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : FFLAGS    = -O3 -check bounds -fp-model source 
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = fds2ftmi_linux_64

intel_linux_64 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- debug Linux 64 ---------
#  64 bit linux

intel_linux_64_db : FFLAGS    = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fltconsistency -stand:f08
intel_linux_64_db : FC        = ifort
intel_linux_64_db : exe       = fds2ftmi_linux_64_db

intel_linux_64_db : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- OSX 32 ----------------
#  32 bit OSX

intel_osx_32 : FFLAGS    = -O3 -check bounds -fp-model source 
intel_osx_32 : FC        = ifort
intel_osx_32 : exe       = fds2ftmi_osx_32

intel_osx_32 : $(obj)
	$(FC) -o $(bin)/$(exe)  -mmacosx-version-min=10.4 $(obj)

# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : FFLAGS    = -O3 -check bounds -fp-model source  
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = fds2ftmi_osx_64

intel_osx_64 : $(obj)
	$(FC) -o $(bin)/$(exe)  -mmacosx-version-min=10.4 $(obj)

# -------------GFORTRAN---------------
# ------------- Win 64 ----------------
#  64 bit GNU windows

gnu_win_64 : FFLAGS    = -O3 -fbounds-check -fno-range-check -ffree-form -frecord-marker=4
gnu_win_64 : FC        = gfortran
gnu_win_64 : obj       = fds2ftmi_gnu_win_64
gnu_win_64 : setup_win $(objgnu)
	$(FC) -o $(obj) $(FFLAGS) $(objgnu)	

# ------------- Win 32 ----------------
#  32 bit GNU windows

gnu_win_32 : FFLAGS    = -O3 -fbounds-check -fno-range-check
gnu_win_32 : FC        = gfortran
gnu_win_32 : obj       = fds2ftmi_gnu_win_32
gnu_win_32 : setup_win $(objgnu)
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
