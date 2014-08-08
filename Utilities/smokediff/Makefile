SOURCE_DIR = ../../../SMV/source
INC_DIR = -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/smokediff
FILTERC =
FILTERF =

VPATH = $(SOURCE_DIR)/smokediff ../../../FDS_Source $(SOURCE_DIR)/shared $(SOURCE_DIR)/smokeview
bin = .

FORTLIBDIR = /opt/intel/composerxe/lib

# Definition of the object variables

obj = assert.o dmalloc.o histogram.o IOdboundary.o IOdplot.o IOdslice.o main.o readsmv.o smv_endian.o file_util.o string_util.o utilities.o fdsmodules.o gsmv.o getdata.o getsizes.o
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
.F90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.F90.obj:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.obj:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<

# ------------- Win 64 debug ----------------
#  64 bit debug windows

intel_win_64_db : CFLAGS    =  -D X64 -D WIN32 /Qvec_report0 /nologo /debug:full /Od /Qtrapuv /fp:strict /fp:except 
intel_win_64_db : FFLAGS    =  -fpp -D pp_cvf -D X64 /iface:stdref /Qvec_report0 /nologo /debug:full /Od /extend_source:132 /warn:unused /warn:nointerfaces /Qtrapuv /fp:strict /fp:except
intel_win_64_db : CC        = icl
intel_win_64_db : CPP       = icl
intel_win_64_db : FC        = ifort
intel_win_64_db : exe       = smokediff_win_64_db.exe

intel_win_64_db : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin)

# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : CFLAGS    = -O2 /Qvec_report0 /nologo -D X64 -D WIN32
intel_win_64 : FFLAGS    = -O2 /Qvec_report0 /nologo -fpp -D pp_cvf -D X64 /iface:stdref
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : FC        = ifort
intel_win_64 : exe       = smokediff_win_64.exe

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin)
	
LINUXFORTLIBS_32=$(FORTLIBDIR)/ia32/libifcore.a
LINUXFORTLIBS_64=$(FORTLIBDIR)/intel64/libifcore.a

# ------------- Linux gcc ----------------
#  gcc linux

gcc_linux : FFLAGS    = -O0 -m64 -x f95-cpp-input -D pp_GCC -D pp_noappend -ffree-form -frecord-marker=4
gcc_linux : CFLAGS    = -O0 -m64 -D pp_LINUX64
gcc_linux : LFLAGS    = -m64
gcc_linux : CC        = gcc
gcc_linux : CPP       = g++
gcc_linux : FC        = gfortran
gcc_linux : exe       = smokediff_linux_64

gcc_linux : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -lgfortran

# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_noappend -D pp_LINUX64
intel_linux_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_LINUX64
intel_linux_64 : CC        = icc
intel_linux_64 : CPP       = icc
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = smokediff_linux_64

intel_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) -static-intel $(obj) $(LINUXFORTLIBS_64)
	

MACFORTLIBS=$(FORTLIBDIR)/libifcoremt.a $(FORTLIBDIR)/libimf.a $(FORTLIBDIR)/libirc.a $(FORTLIBDIR)/libsvml.a

# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_OSX64
intel_osx_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_OSX64
intel_osx_64 : CC        = icc
intel_osx_64 : CPP       = icc
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = smokediff_osx_64

intel_osx_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) -m64 -static-intel $(obj) -mmacosx-version-min=10.4 $(MACFORTLIBS)
