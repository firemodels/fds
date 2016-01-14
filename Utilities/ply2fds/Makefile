SOURCE_DIR = ../../../SMV/source/
INC2_DIR = -I $(SOURCE_DIR)/include2  -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/ply2fds
INC_DIR =  -I $(SOURCE_DIR)/include -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/ply2fds

FILTERC =

VPATH = $(SOURCE_DIR)/ply2fds $(SOURCE_DIR)/shared
bin = .
ifeq ($(shell echo "check_quotes"),"check_quotes")
  GIT_HASH := $(shell ..\..\..\Utilities\Scripts\githash)
  GIT_DATE := $(shell ..\..\..\Utilities\Scripts\gitlog)
else
  GIT_HASH := $(shell git describe --long --dirty)
  GIT_DATE := $(shell git log -1 --format=%cd)
endif
GITINFO=-Dpp_GITHASH=\"$(GIT_HASH)\" -Dpp_GITDATE=\""$(GIT_DATE)\""


# Definition of the object variables

obj = main.o ply.o string_util.o dmalloc.o
objwin = $(obj:.o=.obj)

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .c .obj .o .cpp

.c.obj:
	$(CC) -c $(CFLAGS) $(FILTERC) $(INC_DIR) $<
.c.o:
	$(CC) -c $(CFLAGS) $(FILTERC) $(INC2_DIR) $<

# ********  rules for each platform supported *************

intel_win_64 : CFLAGS    = -Od /nologo -D WIN32 -D PTW32_STATIC_LIB -D X64 -D pp_WIN_INTEL $(GITINFO)
intel_win_64 : FFLAGS    = -Od /nologo -fpp -D pp_cvf -D X64 /iface:stdref /fpscomp:general
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : FC        = ifort
intel_win_64 : exe       = smokezip_win_64.exe

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) $(LIB_DIR)/lib_win_intel_64/pthreads.lib $(LIB_DIR)/lib_win_intel_64/zlib.lib

# ------------- gnu Linux ----------------
#  gnu linux

gnu_linux_64 : FFLAGS    = -O0 -m64 -x f95-cpp-input -D pp_GCC -D pp_noappend -ffree-form -frecord-marker=4
gnu_linux_64 : CFLAGS    = -O0 -m64 -D pp_LINUX $(GITINFO)
gnu_linux_64 : LFLAGS    = -m64
gnu_linux_64 : CC        = gcc
gnu_linux_64 : CPP       = g++
gnu_linux_64 : FC        = gfortran
gnu_linux_64 : exe       = smokezip_linux_64

gnu_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -lgfortran -L $(LIB_DIR)/lib_linux_gnu_64 -lz -lpthread

# ------------- gnu OSX ----------------
#  gnu linux

gnu_osx_64 : FFLAGS    = -O0 -m64 -x f95-cpp-input -D pp_GCC -D pp_noappend -ffree-form -frecord-marker=4
gnu_osx_64 : CFLAGS    = -O0 -m64 -D pp_LINUX -D pp_OSX $(GITINFO)
gnu_osx_64 : LFLAGS    = -m64
gnu_osx_64 : CC        = gcc
gnu_osx_64 : CPP       = g++
gnu_osx_64 : FC        = gfortran
gnu_osx_64 : exe       = smokezip_osx_64

gnu_osx_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L/usr/local/lib/gcc/5 -lgfortran -L $(LIB_DIR)/lib_linux_gnu_64 -lz -lpthread


# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : CFLAGS    = -O2 -m64 -D pp_LINUX $(GITINFO)
intel_linux_64 : CC        = icc
intel_linux_64 : CPP       = icc
intel_linux_64 : exe       = ply2fds_linux_64

intel_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) 

	

# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : CFLAGS    = -O2 -m64 -D pp_OSX $(GITINFO)
intel_osx_64 : CC        = icc
intel_osx_64 : CPP       = icc
intel_osx_64 : exe       = smokezip_osx_64

intel_osx_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) -m64 $(obj)  

#-------------- compile sring_util.c so revision strings are updated -----------

string_util.o:  .FORCE
string_util.obj:        .FORCE

.FORCE:

