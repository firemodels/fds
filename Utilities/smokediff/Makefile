SOURCE_DIR = ../../../SMV_5/source/
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/smokediff $(SOURCE_DIR)/smokeview
bin = .

LINUX32FORTLIBDIR = /opt/intel/fc/10.1.021/lib
OSX32FORTLIBDIR = /opt/intel/fc/10.1.006/lib

# Definition of the object variables

obj = main.o dmalloc.o assert.o utilities.o readsmv.o IOdboundary.o IOdslice.o IOdplot.o getsizes.o getdata.o
objwin = main.obj dmalloc.obj assert.obj utilities.obj readsmv.obj IOdboundary.obj IOdslice.obj IOdplot.obj getdata.obj getsizes.obj

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .F90 .c .obj .o .f90 .cpp

.c.obj:
	$(CC) -c $(CFLAGS) $(FILTERC) $<
.c.o:
	$(CC) -c $(CFLAGS) $(FILTERC) $<
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
ms_win_32 : exe       = smokediff.exe
ms_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin)

# ------------- Win 32 ----------------
#  32 bit windows

intel_win_32 : CFLAGS    = -O2 /Qvec_report0 /nologo -D WIN32
intel_win_32 : FFLAGS    = -O2 /Qvec_report0 /nologo -fpp -D pp_cvf -D pp_PART5
intel_win_32 : CC        = icl
intel_win_32 : CPP       = icl
intel_win_32 : FC        = ifort
intel_win_32 : exe       = smokediff.exe

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
intel_linux_32 : exe       = smokediff_linux_32

intel_linux_32 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) -lifcore -L $(FORTLIBDIR)

# ------------- OSX 32 ----------------
#  32 bit OSX

intel_osx_32 : FFLAGS    = -O2 -m32 -vec_report0 -fpp -D pp_OSX
intel_osx_32 : CFLAGS    = -O2 -m32 -vec_report0 -D pp_OSX
intel_osx_32 : FORTLIBDIR = $(OSX32FORTLIBDIR)
intel_osx_32 : CC        = icc
intel_osx_32 : CPP       = icc
intel_osx_32 : FC        = ifort
intel_osx_32 : exe       = smokediff_osx_32

intel_osx_32 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(FORTLIBDIR)/libifcoremt.a $(FORTLIBDIR)/libimf.a $(FORTLIBDIR)/libirc.a $(FORTLIBDIR)/libsvml.a

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	erase *.obj
