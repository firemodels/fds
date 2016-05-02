SOURCE_DIR = ../../../source/
INC_DIR = -I $(SOURCE_DIR)/zlib128 -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/sh2bat
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/sh2bat:$(SOURCE_DIR)/shared
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

obj = sh2bat.o string_util.o dmalloc.o assert.o file_util.o
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

intel_win_64 : CFLAGS    = -D WIN32 -O2 -D pp_INTEL /nologo $(GITINFO)
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : exe       = sh2bat_win_64.exe

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin)

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	erase *.obj

#-------------- compile sring_util.c so revision strings are updated -----------

string_util.o:  .FORCE
string_util.obj:        .FORCE

.FORCE:

