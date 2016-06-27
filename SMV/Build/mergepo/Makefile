SOURCE_DIR = ../../Source/
INC_DIR = -I $(SOURCE_DIR)/mergepo -I $(SOURCE_DIR)/shared 
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

intel_win_64 : CFLAGS    = -O2 -I ../../source/pthreads /Qvec_report0 /nologo -D X64 -D WIN32 -D PTW32_STATIC_LIB
intel_win_64 : LFLAGS    = ../../Build/LIBS/lib_win_intel_64/pthreads.lib
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : exe       = mergepo_win_64.exe

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) -l $(LFLAGS) $(objwin)

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_LINUX64
intel_linux_64 : CC        = icc
intel_linux_64 : CPP       = icc
intel_linux_64 : exe       = mergepo_linux_64

intel_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj)

# ------------- Linux 64 debug ----------------
#  64 bit linux debug

intel_linux_64_db : CFLAGS    = -O0 -g -traceback -check-uninit -ftrapuv -m64 -vec_report0 -D pp_LINUX64
intel_linux_64_db : LFLAGS    = -traceback
intel_linux_64_db : CC        = icc
intel_linux_64_db : CPP       = icc
intel_linux_64_db : exe       = mergepo_linux_64_db

intel_linux_64_db : $(obj)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(obj)
