SOURCE_DIR = ../../source/
INC_DIR = -I $(SOURCE_DIR)/makepo  -I $(SOURCE_DIR)/shared 
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/makepo:$(SOURCE_DIR)/shared
bin = .

LINUX64FORTLIBDIR = $(IFORT_COMPILER)/lib/intel64

# Definition of the object variables

obj = main.o string_util.o dmalloc.o assert.o
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
intel_win_64 : exe       = makepo_win_64.exe 

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) $(LFLAGS)

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_LINUX64
intel_linux_64 : LFLAGS    = -lpthread
intel_linux_64 : CC        = icc
intel_linux_64 : CPP       = icc
intel_linux_64 : exe       = makepo_linux_64

intel_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj)

# ------------- gcc Linux 64 ----------------
#  64 bit linux

gcc_linux_64 : CFLAGS    = -O2 -m64 -D pp_LINUX64
gcc_linux_64 : CC        = gcc
gcc_linux_64 : CPP       = g++
gcc_linux_64 : FC        = ifort
gcc_linux_64 : exe       = makepo_linux_64

gcc_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) 


#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm *.o
