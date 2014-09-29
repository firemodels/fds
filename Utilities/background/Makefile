SOURCE_DIR = ../../../SMV/source
INC_DIR = -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/background
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/background $(SOURCE_DIR)/shared
bin = .

# Definition of the object variables

obj = main.o string_util.o dmalloc.o
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

# ********  rules for each platform supported *************

# ------------- Win 32 ----------------
#  32 bit windows

intel_win_32 : CFLAGS    = -O1 /nologo -D WIN32 -I $(SOURCE_DIR)/shared
intel_win_32 : CC        = icl
intel_win_32 : exe       = background.exe

intel_win_32 : $(objwin)
	$(CC) -o $(bin)/$(exe) $(objwin) 

# ------------- Linux 32 ----------------
#  32 bit linux

intel_linux_32 : CFLAGS    = -O2 -m32 -D pp_LINUX
intel_linux_32 : CC        = icc
intel_linux_32 : exe       = background

intel_linux_32 : $(obj)
	$(CC) -o $(bin)/$(exe) $(obj) 

# ------------- gcc Linux 32 ----------------

gcc_linux_32 : CFLAGS    = -O2 -m32 -D pp_LINUX
gcc_linux_32 : CC        = gcc
gcc_linux_32 : exe       = background

gcc_linux_32 : $(obj)
	$(CC) -m32 -o $(bin)/$(exe) $(obj) -lm


# ------------- OSX 32 ----------------
#  32 bit osx

intel_osx_32 : CFLAGS    = -O2 -m32 -D pp_OSX
intel_osx_32 : CC        = icc
intel_osx_32 : exe       = background

intel_osx_32 : $(obj)
	$(CC) -o $(bin)/$(exe) -m32 $(obj) 

# ------------- gcc OSX 32 ----------------

gcc_osx_32 : CFLAGS    = -O2 -m32 -D pp_OSX
gcc_osx_32 : CC        = gcc
gcc_osx_32 : exe       = background

gcc_osx_32 : $(obj)
	$(CC) -m32 -o $(bin)/$(exe) $(obj) -lm
