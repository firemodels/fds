SOURCE_DIR = ..\..\..\Source\get_time
INC_DIR =  -I $(SOURCE_DIR)
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)
bin = .

# Definition of the object variables

obj =    get_time.o 
objwin = $(obj:.o=.obj)

# main has to be compiled without optimization (on a Windows PC) otherwise errors occcur
main.obj: CFLAGS = -Od /nologo -D WIN32

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .F90 .c .obj .o .f90 .cpp

.c.obj:
	$(CC) -c $(CFLAGS) $(FILTERC) $(INC_DIR) $<
.c.o:
	$(CC) -c $(CFLAGS) $(FILTERC) $(INC_DIR) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.obj:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<



# ********  rules for each platform supported *************

# ------------- Win 32 ----------------
#  32 bit windows

intel_win_32 : CFLAGS    = -O2 /nologo -D WIN32 -D PTW32_STATIC_LIB
intel_win_32 : FFLAGS    = -O2 /nologo -fpp -D pp_cvf /iface:stdref
intel_win_32 : CC        = icl
intel_win_32 : CPP       = icl
intel_win_32 : FC        = ifort
intel_win_32 : exe       = smokezip_win_32.exe

intel_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) $(LIB_DIR)/lib_win_intel_32/pthreads.lib $(LIB_DIR)/lib_win_intel_32/zlib.lib


intel_win_64 : CFLAGS    = -Od /nologo -D WIN32 -D X64
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : FC        = ifort
intel_win_64 : exe       = get_time_64.exe

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) 

# ------------- Linux 32 profile ----------------
#  32 bit linux

intel_linux_32p : FFLAGS    = -O2 -m32 -fpp -D pp_noappend
intel_linux_32p : CFLAGS    = -O2 -m32 -D pp_LINUX -p
intel_linux_32p : CC        = icc
intel_linux_32p : CPP       = icc
intel_linux_32p : FC        = ifort
intel_linux_32p : exe       = smokezip_linux_32p

intel_linux_32p : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) -p -L $(FORTLIBDIR)/ia32 -lifcore -L $(LIB_DIR)/lib_linux_intel_32 -lz
	
LINUXFORTLIBS_32=$(FORTLIBDIR)/ia32/libifcore.a
LINUXFORTLIBS_64=$(FORTLIBDIR)/intel64/libifcore.a

# ------------- Linux 32 ----------------
#  32 bit linux

intel_linux_32 : FFLAGS    = -O2 -m32 -fpp -D pp_noappend
intel_linux_32 : CFLAGS    = -O2 -m32 -D pp_LINUX 
intel_linux_32 : CC        = icc
intel_linux_32 : CPP       = icc
intel_linux_32 : FC        = ifort
intel_linux_32 : exe       = smokezip_linux_32

intel_linux_32 : $(obj)
	$(CPP) -o $(bin)/$(exe) -static-intel $(obj) $(LINUXFORTLIBS_32) -L $(LIB_DIR)/lib_linux_intel_32 -lz -lpthread

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : FFLAGS    = -O2 -m64 -fpp -D pp_noappend -D pp_LINUX
intel_linux_64 : CFLAGS    = -O2 -m64 -D pp_LINUX
intel_linux_64 : CC        = icc
intel_linux_64 : CPP       = icc
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = smokezip_linux_64

intel_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) -static-intel $(obj) $(LINUXFORTLIBS_64) -L $(LIB_DIR)/lib_linux_intel_64 -lz -lpthread
	
OSXFORTLIBS=$(FORTLIBDIR)/libifcoremt.a $(FORTLIBDIR)/libimf.a $(FORTLIBDIR)/libirc.a $(FORTLIBDIR)/libsvml.a

# ------------- OSX 32 ----------------
#  32 bit OSX

intel_osx_32 : FFLAGS    = -O2 -m32 -fpp -D pp_OSX
intel_osx_32 : CFLAGS    = -O2 -m32 -D pp_OSX
intel_osx_32 : CC        = icc
intel_osx_32 : CPP       = icc
intel_osx_32 : FC        = ifort
intel_osx_32 : exe       = smokezip_osx_32

intel_osx_32 : $(obj)
	$(CPP) -o $(bin)/$(exe) -m32 -static-intel $(obj) -mmacosx-version-min=10.4 $(OSXFORTLIBS) -L $(LIB_DIR)/lib_osx_intel_32 -lz -lpthread

# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : FFLAGS    = -O2 -m64 -fpp -D pp_OSX
intel_osx_64 : CFLAGS    = -O2 -m64 -D pp_OSX
intel_osx_64 : CC        = icc
intel_osx_64 : CPP       = icc
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = smokezip_osx_64

intel_osx_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) -m64 -static-intel $(obj)  -mmacosx-version-min=10.4 $(OSXFORTLIBS) -L $(LIB_DIR)/lib_osx_intel_64 -lz -lpthread
