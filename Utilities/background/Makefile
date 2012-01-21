SOURCE_DIR = ../../../SMV/source/
INC_DIR = -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/background
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/background $(SOURCE_DIR)/shared
bin = .

# Definition of the object variables

objwin = main.obj string_util.obj dmalloc.obj
obj = main.o string_util.o dmalloc.o

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

# ------------- Win 32 ----------------
#  32 bit windows

intel_win_32 : CFLAGS    = -O1 /Qvec_report0 /nologo -D WIN32 -I $(SOURCE_DIR)/shared
intel_win_32 : FFLAGS    = 
intel_win_32 : CC        = icl
intel_win_32 : CPP       = icl
intel_win_32 : FC        = ifort
intel_win_32 : exe       = background.exe

intel_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) 

# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : CFLAGS    = -O2 /Qvec_report0 /nologo -D WIN32
intel_win_64 : FFLAGS    = 
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : FC        = ifort
intel_win_64 : exe       = background.exe

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) 

# ------------- Linux 32 ----------------
#  32 bit linux

intel_linux_32 : FFLAGS    = 
intel_linux_32 : CFLAGS    = -O2 -m32 -vec_report0 -D pp_LINUX
intel_linux_32 : CC        = icc
intel_linux_32 : CPP       = icc
intel_linux_32 : FC        = ifort
intel_linux_32 : exe       = background

intel_linux_32 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) 

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm *.o
