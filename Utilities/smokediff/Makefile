SOURCE_DIR = ../../SMV_5/source/
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/setfdspath
bin = .

# Definition of the object variables

obj = main.obj dmalloc.o assert.o utilities.o readsmv.o IOdslice.o IOdplot.o
objwin = $(obj:.o=.obj)

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .F90 .c .obj .o .f90 .cpp

.c.obj:
	$(CC) -c $(CFLAGS) $(FILTERC) $<
.c.o:
	$(CC) -c $(CFLAGS) $(FILTERC) $<

# ********  rules for each platform supported *************

# ------------- Win 32 ----------------
#  32 bit windows

CFLAGS    = -O2 /Qvec_report0 /nologo -D WIN32
CC        = icl
CPP       = icl
FC        = ifort
exe       = smokediff.exe

intel_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) smokeview_fort.lib

install	:	$(exe)
	copy $(exe) ..\Scripts\bundle_setup\$(exe)

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	erase *.obj
