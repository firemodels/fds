SOURCE_DIR = ../../../SMV_5/source/
INC_DIR = -I $(SOURCE_DIR)/shared
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

obj = assert.o CNV3dsmoke.o CNVboundary.o CNViso.o CNVpart.o CNVplot3d.o CNVslice.o csphere.o dmalloc.o egz_stdio.o endian.o histogram.o main.o readfiles.o utilities.o getdata.o getsizes.o
objwin = main.obj dmalloc.obj assert.obj utilities.obj readsmv.obj IOdboundary.obj IOdslice.obj IOdplot.obj getdata.obj getsizes.obj histogram.obj

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
ms_win_32 : FFLAGS    = -O2 -fpp -D pp_cvf -D pp_PART5
ms_win_32 : CC        = cl
ms_win_32 : CPP       = cl
ms_win_32 : FC        = ifort
ms_win_32 : exe       = smokezip.exe
ms_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin)

# ------------- Win 32 ----------------
#  32 bit windows

intel_win_32 : CFLAGS    = -O2 /Qvec_report0 /nologo -D WIN32
intel_win_32 : FFLAGS    = -O2 /Qvec_report0 /nologo -fpp -D pp_cvf -D pp_PART5
intel_win_32 : CC        = icl
intel_win_32 : CPP       = icl
intel_win_32 : FC        = ifort
intel_win_32 : exe       = smokezip.exe

intel_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin)

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
	$(CPP) -o $(bin)/$(exe) $(obj) -lifcore -L $(FORTLIBDIR) ../linux_lib32/libz.a

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
	$(CPP) -o $(bin)/$(exe) -m32 $(obj) $(FORTLIBDIR)/libifcoremt.a $(FORTLIBDIR)/libimf.a $(FORTLIBDIR)/libirc.a $(FORTLIBDIR)/libsvml.a ../osx_lib32/libz.a

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
	$(CPP) -o $(bin)/$(exe) -m64 $(obj) $(FORTLIBDIR)/libifcoremt.a $(FORTLIBDIR)/libimf.a $(FORTLIBDIR)/libirc.a $(FORTLIBDIR)/libsvml.a ../osx_lib64/libz.a

#-------------- Clean Target to remove Object and Module files -----------
#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	erase *.obj
