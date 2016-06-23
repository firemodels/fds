SOURCE_DIR = ../../../Source/
INC_DIR =  -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/ply2fds

FILTERC =

VPATH = $(SOURCE_DIR)/ply2fds $(SOURCE_DIR)/shared
bin = .
ifeq ($(shell echo "check_quotes"),"check_quotes")
  GIT_HASH := $(shell ..\..\..\..\Utilities\Scripts\githash)
  GIT_DATE := $(shell ..\..\..\..\Utilities\Scripts\gitlog)
else
  GIT_HASH := $(shell git describe --long --dirty)
  GIT_DATE := $(shell git log -1 --format=%cd)
endif
GITINFO=-Dpp_GITHASH=\"$(GIT_HASH)\" -Dpp_GITDATE=\""$(GIT_DATE)\""


# Definition of the object variables

obj = main.o ply.o string_util.o dmalloc.o assert.o file_util.o
objwin = $(obj:.o=.obj)

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .c .obj .o .cpp

.c.obj:
	$(CC) -c $(CFLAGS) $(FILTERC) $(INC_DIR) $<
.c.o:
	$(CC) -c $(CFLAGS) $(FILTERC) $(INC_DIR) $<

# ********  rules for each platform supported *************

intel_win_64 : CFLAGS    = -Od /nologo -D WIN32 -D PTW32_STATIC_LIB -D X64 -D pp_WIN_INTEL $(GITINFO)
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : exe       = ply2fds_win_64.exe

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin)

# ------------- gnu Linux ----------------
#  gnu linux

gnu_linux_64 : CFLAGS    = -O0 -m64 -D pp_LINUX $(GITINFO)
gnu_linux_64 : LFLAGS    = -m64
gnu_linux_64 : CC        = gcc
gnu_linux_64 : CPP       = g++
gnu_linux_64 : exe       = ply2fds_linux_64

gnu_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS)

# ------------- gnu OSX ----------------
#  gnu linux

gnu_osx_64 : CFLAGS    = -O0 -m64 -D pp_LINUX -D pp_OSX $(GITINFO)
gnu_osx_64 : LFLAGS    = -m64
gnu_osx_64 : CC        = gcc
gnu_osx_64 : CPP       = g++
gnu_osx_64 : exe       = ply2fds_osx_64

gnu_osx_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS)


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
intel_osx_64 : exe       = ply2fds_osx_64

intel_osx_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) -m64 $(obj)  

#-------------- compile sring_util.c so revision strings are updated -----------

string_util.o:  .FORCE
string_util.obj:        .FORCE

.FORCE:

