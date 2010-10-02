#*** Smokeview, Version 5, Makefile ***

#To use this Makefile, you will need to provide a make target at the command line.
#For example, to make the 32 bit version in the current directory with the Intel
#compiler, you would need to type 'make intel_linux_32' at the prompt.
#To make it in a sub-directory type 'make -f ../Makefile intel_linux_32'

# use the following line if building Smokeview in the directory containing this Makefile
#SOURCE_DIR = ../source/

# use the following line if building Smokeview in a directory 'below' Makefile

SOURCE_DIR = ../../../SMV_5/source/

# use FILTERC and FILTERF to customize smokeview builds (ie -D pp_BETA)

FILTERC =
FILTERF =

#----------------------------------------------------------------------------
# Should not need to edit lines below except to add or 'debug' target entries

# The variable VPATH defines the source code directory relative to the current directory

VPATH = $(SOURCE_DIR)/thread_test
bin = .

# BUILD_DIR defined directory containing this Makefile

BUILD_DIR = $(SOURCE_DIR)/../Build
LIB_DIR = $(BUILD_DIR)/LIBS
INC_DIR = $(SOURCE_DIR)/include2

# SYSTEM libraries

SYSTEMLIBS = -lpthread 

# FORTRAN Library directories for various platforms

#OSX32FORTLIBDIR  = /opt/intel/fc/10.1.006/lib
OSX32FORTLIBDIR  = /opt/intel/Compiler/11.1/080/lib
OSX64FORTLIBDIR  = /opt/intel/Compiler/11.1/080/lib
#LINUX32FORTLIBDIR  = /opt/intel/fc/10.1.021/lib
LINUX32FORTLIBDIR  = /opt/intel/11/lib/ia32
LINUX64FORTLIBDIR  = /opt/intel/11/lib/intel64

# Definition of the object variables

obj = main.o 
objwin = $(obj:.o=.obj)

# include file directory locations (should be the same for all builds)

INC = -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/smokeview

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .F90 .c .obj .o .f90 .cpp

.c.obj:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<
.c.o:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<

.cpp.obj:
	$(CPP) -c $(CFLAGS) $(INC) $(FILTERC) $<
.cpp.o:
	$(CPP) -c $(CFLAGS) $(INC) $(FILTERC) $<

.F90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.obj:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<

# ********  rules for each platform supported *************

# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : INC = -I ../../../SMV_5/source/include 
intel_win_64 : IMAGELIBDIR = $(LIB_DIR)/lib_win_intel_64
intel_win_64 : CFLAGS    = -O2 -D WIN32 -D _CONSOLE -D pp_WIN_INTEL -D X64 /Qvec_report0
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : exe       = thread_test

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) $(IMAGELIBDIR)\pthreadVC2_x64.lib


# ------------- Win 32 ----------------
#  32 bit windows

intel_win_32 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview
intel_win_32 : IMAGELIBDIR = $(LIB_DIR)/lib_win_intel_32
intel_win_32 : WINLIBS = kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib
intel_win_32 : FFLAGS    = -O2 -fpp -D pp_cvf /Qvec_report0 /iface:stdref
intel_win_32 : CFLAGS    = -O2 -D WIN32 -D _CONSOLE -D pp_WIN_INTEL /Qvec_report0
intel_win_32 : CC        = icl
intel_win_32 : CPP       = icl
intel_win_32 : FC        = ifort
intel_win_32 : exe       = smv5_win_32

intel_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) opengl32.lib $(WINLIBS) $(IMAGELIBDIR)\pthreadVC.lib $(IMAGELIBDIR)\glui.lib $(IMAGELIBDIR)\glut32.lib $(IMAGELIBDIR)\glew32.lib $(IMAGELIBDIR)\gd.lib $(IMAGELIBDIR)\jpeg.lib $(IMAGELIBDIR)\png.lib $(IMAGELIBDIR)\zlib.lib   /link /machine:x64


#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm *.o
