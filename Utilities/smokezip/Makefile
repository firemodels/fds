SOURCE_DIR = ../../../SMV_5/source/
BUILD_DIR = $(SOURCE_DIR)/../Build
LIB_DIR = $(BUILD_DIR)/LIBS
INC_DIR = -I $(SOURCE_DIR)/include -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/smokezip 
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/smokezip $(SOURCE_DIR)/shared $(SOURCE_DIR)/smokeview 
bin = .

#LINUX32FORTLIBDIR = /opt/intel/fc/10.1.021/lib
LINUX32FORTLIBDIR = /opt/intel/11/lib/ia32
LINUX64FORTLIBDIR = /opt/intel/11/lib/intel64
OSX32FORTLIBDIR = /opt/intel/Compiler/11.1/080/lib
OSX64FORTLIBDIR = /opt/intel/Compiler/11.1/080/lib

# Definition of the object variables

obj =    main.o CNV3dsmoke.o CNVboundary.o CNViso.o CNVpart.o CNVplot3d.o CNVslice.o assert.o csphere.o dmalloc.o egz_stdio.o endian.o histogram.o readfiles.o utilities.o isobox.o threader.o getdata.o getsizes.o
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
	$(CC) -c $(CFLAGS) $(FILTERC) $(INC_DIR) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.obj:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<



# ********  rules for each platform supported *************

# ------------- MS Win 32 ----------------
#  32 bit windows (ms visual studio C compiler)

ms_win_32 : CFLAGS    = -O2 -D WIN32
ms_win_32 : FFLAGS    = -O2 -fpp -D pp_cvf -D pp_PART5 /iface:stdref
ms_win_32 : CC        = cl
ms_win_32 : CPP       = cl
ms_win_32 : FC        = ifort
ms_win_32 : exe       = smokezip.exe
ms_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) $(LIB_DIR)/lib_win_intel_32/zlib.lib  $(LIB_DIR)/lib_win_intel_32/pthreadVC2.lib

# ------------- Win 32 ----------------
#  32 bit windows

intel_win_32 : CFLAGS    = -O2 /Qvec_report0 /nologo -D WIN32 -D pp_WIN_INTEL
intel_win_32 : FFLAGS    = -O2 /Qvec_report0 /nologo -fpp -D pp_cvf /iface:stdref
intel_win_32 : CC        = icl
intel_win_32 : CPP       = icl
intel_win_32 : FC        = ifort
intel_win_32 : exe       = smokezip.exe

intel_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) $(LIB_DIR)/lib_win_intel_32/zlib.lib $(LIB_DIR)/lib_win_intel_32/pthreadVC2.lib

# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : CFLAGS    = -Od /Qvec_report0 /nologo -D WIN32 -D X64 -D pp_WIN_INTEL
intel_win_64 : FFLAGS    = -Od /Qvec_report0 /nologo -fpp -D pp_cvf -D X64 /iface:stdref
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : FC        = ifort
intel_win_64 : exe       = smokezip_64.exe

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) $(LIB_DIR)/lib_win_intel_64/zlib.lib $(LIB_DIR)/lib_win_intel_64/pthreadVC2_x64.lib

# ------------- Linux 32 profile ----------------
#  32 bit linux

intel_linux_32p : FFLAGS    = -O2 -m32 -vec_report0 -fpp -D pp_noappend
intel_linux_32p : CFLAGS    = -O2 -m32 -vec_report0 -D pp_LINUX -p
intel_linux_32p : FORTLIBDIR = $(LINUX32FORTLIBDIR)
intel_linux_32p : CC        = icc
intel_linux_32p : CPP       = icc
intel_linux_32p : FC        = ifort
intel_linux_32p : exe       = smokezip_linux_32p

intel_linux_32p : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) -p -lifcore -L $(FORTLIBDIR) ../linux_lib32/libz.a

# ------------- Linux 32 ----------------
#  32 bit linux

intel_linux_32 : FFLAGS    = -O2 -m32 -vec_report0 -fpp -D pp_noappend
intel_linux_32 : CFLAGS    = -O2 -m32 -vec_report0 -D pp_LINUX 
intel_linux_32 : FORTLIBDIR = $(LINUX32FORTLIBDIR)
intel_linux_32 : CC        = icc
intel_linux_32 : CPP       = icc
intel_linux_32 : FC        = ifort
intel_linux_32 : exe       = smokezip_linux_32

intel_linux_32 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) -p -lifcore -L $(FORTLIBDIR) ../linux_lib32/libz.a

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_noappend
intel_linux_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_LINUX
intel_linux_64 : FORTLIBDIR = $(LINUX64FORTLIBDIR)
intel_linux_64 : CC        = icc
intel_linux_64 : CPP       = icc
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = smokezip_linux_64

intel_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) -lifcore -L $(FORTLIBDIR) ../linux_lib64/libz.a

# ------------- OSX 32 ----------------
#  32 bit OSX

intel_osx_32 : FFLAGS    = -O2 -m32 -vec_report0 -fpp -D pp_OSX
intel_osx_32 : CFLAGS    = -O2 -m32 -vec_report0 -D pp_OSX
intel_osx_32 : FORTLIBDIR = $(OSX32FORTLIBDIR)
intel_osx_32 : CC        = icc
intel_osx_32 : CPP       = icc
intel_osx_32 : FC        = ifort
intel_osx_32 : exe       = smokezip_osx_32

intel_osx_32 : $(obj)
	$(CPP) -o $(bin)/$(exe) -m32 $(obj)  -mmacosx-version-min=10.4 $(FORTLIBDIR)/libifcoremt.a $(FORTLIBDIR)/libimf.a $(FORTLIBDIR)/libirc.a $(FORTLIBDIR)/libsvml.a ../osx_lib32/libz.a

# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_OSX
intel_osx_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_OSX
intel_osx_64 : FORTLIBDIR = $(OSX32FORTLIBDIR)
intel_osx_64 : CC        = icc
intel_osx_64 : CPP       = icc
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = smokezip_osx_64

intel_osx_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) -m64 $(obj)  -mmacosx-version-min=10.4 $(FORTLIBDIR)/libifcoremt.a $(FORTLIBDIR)/libimf.a $(FORTLIBDIR)/libirc.a $(FORTLIBDIR)/libsvml.a ../osx_lib64/libz.a

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm *.o
