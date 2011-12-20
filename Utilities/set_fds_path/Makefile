SOURCE_DIR = ../../../SMV/source/
INC_DIR = -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/set_path
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/set_path $(SOURCE_DIR)/shared
bin = .

# Definition of the object variables

obj = main.o string_util.o
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

intel_win_32 : CFLAGS    = -D WIN32 -O2 /Qvec_report0 /nologo 
intel_win_32 : CC        = icl
intel_win_32 : CPP       = icl
intel_win_32 : FC        = ifort
intel_win_32 : exe       = set_path32.exe

intel_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) Advapi32.lib 

# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : CFLAGS    = -D WIN32 -O2 /Qvec_report0 /nologo 
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : FC        = ifort
intel_win_64 : exe       = set_path64.exe

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) Advapi32.lib 

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	erase *.obj
