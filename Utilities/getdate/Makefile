SOURCE_DIR = ../../../SMV/source/
BUILD_DIR = ../../../SMV/Build
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/getdate $(SOURCE_DIR)/shared  
bin = .


# Definition of the object variables

obj =    main.o
objwin = $(obj:.o=.obj)

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .c .obj .o 

.c.obj:
	$(CC) -c $(CFLAGS) $(FILTERC) $<
.c.o:
	$(CC) -c $(CFLAGS) $(FILTERC) $<



# ********  rules for each platform supported *************

intel_win_64 : CFLAGS    = -Od /nologo
intel_win_64 : CC        = icl
intel_win_64 : exe       = getdate_win_64.exe

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin)

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : CFLAGS    = -O2 -m64
intel_linux_64 : CC        = icc
intel_linux_64 : exe       = getdate_linux_64

intel_linux_64 : $(obj)
	$(CC) -o $(bin)/$(exe) $(obj) 

	
# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : CFLAGS    = -O2 -m64
intel_osx_64 : CC        = icc
intel_osx_64 : exe       = getdate_osx_64

intel_osx_64 : $(obj)
	$(CC) -o $(bin)/$(exe) -m64 $(obj)
