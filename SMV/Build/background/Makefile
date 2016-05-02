SOURCE_DIR = ../../../source
INC_DIR = -I $(SOURCE_DIR)/zlib128 -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/background
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/background $(SOURCE_DIR)/shared
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

obj = main.o string_util.o dmalloc.o file_util.o
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

# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : CFLAGS    = -O1 /nologo -D WIN32 -D pp_INTEL -I $(SOURCE_DIR)/shared $(GITINFO)
intel_win_64 : CC        = icl
intel_win_64 : exe       = background.exe

intel_win_64 : $(objwin)
	$(CC) -o $(bin)/$(exe) $(objwin) 

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : CFLAGS    = -O2 -m64 -D pp_INTEL -D pp_LINUX $(GITINFO)
intel_linux_64 : CC        = icc
intel_linux_64 : exe       = background

intel_linux_64 : $(obj)
	$(CC) -o $(bin)/$(exe) $(obj) 

# ------------- gnu Linux 64 ----------------

gnu_linux_64 : CFLAGS    = -O2 -m64 -D pp_LINUX $(GITINFO)
gnu_linux_64 : CC        = gcc
gnu_linux_64 : exe       = background

gnu_linux_64: $(obj)
	$(CC) -m64 -o $(bin)/$(exe) $(obj) -lm

# ------------- OSX 64 ----------------
#  64 bit osx

intel_osx_64 : CFLAGS    = -O2 -m64 -D pp_INTEL -D pp_OSX $(GITINFO)
intel_osx_64 : CC        = icc
intel_osx_64 : exe       = background

intel_osx_64 : $(obj)
	$(CC) -o $(bin)/$(exe) -m64 $(obj) 

# ------------- gcc OSX 64 ----------------
gnu_osx_64 : CFLAGS    = -O2 -m64 -D pp_OSX $(GITINFO)
gnu_osx_64 : CC        = gcc
gnu_osx_64 : exe       = background

gnu_osx_64 : $(obj)
	$(CC) -m64 -o $(bin)/$(exe) $(obj)  -lm

#-------------- compile sring_util.c so revision strings are updated -----------

string_util.o:  .FORCE
string_util.obj:        .FORCE

.FORCE:

