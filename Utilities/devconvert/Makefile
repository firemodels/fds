SOURCE_DIR = ../../../SMV/source/
BUILD_DIR = $(SOURCE_DIR)/../Build
LIB_DIR = $(BUILD_DIR)/LIBS
INC2_DIR = -I $(SOURCE_DIR)/include2  -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/smokezip 
INC_DIR =  -I $(SOURCE_DIR)/include -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/smokezip 
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/devconvert $(SOURCE_DIR)/shared 
bin = .

#LINUX32FORTLIBDIR = /opt/intel/fc/10.1.021/lib
LINUX32FORTLIBDIR = /opt/intel/11/lib/ia32
LINUX64FORTLIBDIR = /opt/intel/11/lib/intel64
OSX32FORTLIBDIR = /opt/intel/Compiler/11.1/080/lib
OSX64FORTLIBDIR = /opt/intel/Compiler/11.1/080/lib

# Definition of the object variables

obj =    main.o assert.o dmalloc.o string_util.o
objwin = $(obj:.o=.obj)

# main has to be compiled without optimization (on a Windows PC) otherwise errors occcur
main.obj: CFLAGS = -Od /Qvec_report0 /nologo -D WIN32 -D pp_WIN_INTEL

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .F90 .c .obj .o .f90 .cpp

.c.obj:
	$(CC) -c $(CFLAGS) $(FILTERC) $(INC_DIR) $<
.c.o:
	$(CC) -c $(CFLAGS) $(FILTERC) $(INC2_DIR) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.obj:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<



# ********  rules for each platform supported *************


# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : CFLAGS    = -Od /Qvec_report0 /nologo -D WIN32 -D X64 -D pp_WIN_INTEL
intel_win_64 : FFLAGS    = -Od /Qvec_report0 /nologo -fpp -D pp_cvf -D X64 /iface:stdref
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : FC        = ifort
intel_win_64 : exe       = devconvert_win_64.exe

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) $(LIB_DIR)/lib_win_intel_64/pthreadVC2_x64.lib

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_noappend -D pp_LINUX64
intel_linux_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_LINUX64
intel_linux_64 : FORTLIBDIR = $(LINUX64FORTLIBDIR)
intel_linux_64 : CC        = icc
intel_linux_64 : CPP       = icc
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = devconvert_linux_64

intel_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj)


# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_OSX64
intel_osx_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_OSX64
intel_osx_64 : FORTLIBDIR = $(OSX32FORTLIBDIR)
intel_osx_64 : CC        = icc
intel_osx_64 : CPP       = icc
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = smokezip_osx_64

intel_osx_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) -m64 $(obj)

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm *.o
