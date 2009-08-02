SOURCE_DIR = ../../SMV_5/source/
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/setfdspath
bin = .

# Definition of the object variables

obj = main.obj dmalloc.o assert.o
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

intel_win_32 : CFLAGS    = -O2 /Qvec_report0
intel_win_32 : CC        = icl
intel_win_32 : CPP       = icl
intel_win_32 : FC        = ifort
intel_win_32 : exe       = set_fds_path.exe

intel_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) 

install	:	set_fds_path.exe
	copy set_fds_path.exe ..\Scripts\bundle_setup\set_fds_path.exe

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	erase *.obj
