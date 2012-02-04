SOURCE_DIR = ../../source/
INC_DIR = -I $(SOURCE_DIR)/mergepo -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/include2
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/mergepo:$(SOURCE_DIR)/shared
bin = .

LINUX64FORTLIBDIR = $(IFORT_COMPILER)/lib/intel64

# Definition of the object variables

obj = main.o dmalloc.o assert.o translate.o string_util.o
objwin = $(obj:.o=.obj)

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .F90 .c .obj .o .f90 .cpp

.c.obj:
	$(CC) -c $(CFLAGS) $(FILTERC) $(INC_DIR) $<
.c.o:
	$(CC) -c $(CFLAGS) $(FILTERC) $(INC_DIR) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.F90.obj:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.obj:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<



# ********  rules for each platform supported *************

# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : CFLAGS    = -O2 /Qvec_report0 /nologo -D X64 -D WIN32
intel_win_64 : FFLAGS    = -O2 /Qvec_report0 /nologo -fpp -D pp_cvf -D X64 /iface:stdref
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : FC        = ifort
intel_win_64 : exe       = mergepo_win_64.exe ../../Build/LIBS/lib_win_intel_64/pthreadVC2_x64.lib

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin)

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_noappend -D pp_LINUX64
intel_linux_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_LINUX64
intel_linux_64 : FORTLIBDIR = $(LINUX64FORTLIBDIR)
intel_linux_64 : CC        = icc
intel_linux_64 : CPP       = icc
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = mergepo_linux_64

intel_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) -lifcore -L $(FORTLIBDIR)

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm *.o
