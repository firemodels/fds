# $Date$ 
# $Revision$
# $Author$

#*** Smokeview, VMakefile ***

#To use this Makefile cd to a sub-directory and type make_smv.sh or make_smv.bat

SOURCE_DIR = ../../source/
SMV_TESTFLAG = 
SMV_TESTSTRING =
SMV_PROFILEFLAG = 
SMV_PROFILESTRING =
ifeq ($(shell echo "check_quotes"),"check_quotes")
  GIT_HASH := $(shell ..\..\..\Utilities\Scripts\githash)
  GIT_DATE := $(shell ..\..\..\Utilities\Scripts\gitlog)
else
  GIT_HASH := $(shell git describe --long --dirty)
  GIT_DATE := $(shell git log -1 --format=%cd)
endif
GITINFO=-Dpp_GITHASH=\"$(GIT_HASH)\" -Dpp_GITDATE=\""$(GIT_DATE)\""

#----------------------------------------------------------------------------
# Should not need to edit lines below except to add or 'debug' target entries

# The variable VPATH defines the source code directory relative to the current directory

VPATH = $(SOURCE_DIR)/smokeview:../../../FDS_Source:$(SOURCE_DIR)/shared:$(SOURCE_DIR)/glew
bin = .

# BUILD_DIR defined directory containing this Makefile

BUILD_DIR = $(SOURCE_DIR)/../Build
LIB_DIR = $(BUILD_DIR)/LIBS
INC_DIR = $(SOURCE_DIR)/include2

# Definition of the object variables

obj = main.o menus.o IOscript.o IOshooter.o glui_devices.o glui_shooter.o glui_wui.o\
      csphere.o colortimebar.o smv_endian.o camera.o IOgeometry.o IOhrr.o IOwui.o\
      IOobject.o IOtour.o glui_display.o glui_tour.o getdatacolors.o smokeview.o\
      output.o renderfile.o isobox.o getdatabounds.o readsmv.o scontour2d.o\
      glui_3dsmoke.o glui_clip.o glui_stereo.o glui_geometry.o glui_motion.o\
      glui_bounds.o dmalloc.o interp.o assert.o compress3dc.o\
      compress.o IOvolsmoke.o IOsmoke.o IOplot3d.o IOslice.o IOboundary.o\
      IOpart.o IOzone.o IOiso.o IOtarget.o callbacks.o drawGeometry.o\
      glui_colorbar.o skybox.o file_util.o string_util.o startup.o glui_trainer.o\
      shaders.o unit.o threader.o histogram.o translate.o update.o viewports.o\
      smv_geometry.o geometry.o showscene.o fdsmodules.o\
      gsmv.o getdata.o getsizes.o glew.o
objwin = $(obj:.o=.obj)

# include file directory locations (should be the same for all builds)

INC = -I $(INC_DIR)/glut -I $(INC_DIR)/glui -I $(INC_DIR)/gd -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/smokeview -I $(SOURCE_DIR)/glew -I $(INC_DIR)/zlib

# libraries used by smokeview

SMV_LIBS = -lglui -lgd -ljpeg -lpng -lz
FULL_SMV_LIBS = $(SMV_LIBS) -lglut

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .F90 .c .obj .o .f90 .cpp

.c.obj:
	$(CC) -c $(CFLAGS) $(INC) $<
.c.o:
	$(CC) -c $(CFLAGS) $(INC) $<

.cpp.obj:
	$(CPP) -c $(CFLAGS) $(INC) $<
.cpp.o:
	$(CPP) -c $(CFLAGS) $(INC) $<

.F90.o:
	$(FC) -c $(FFLAGS) $<
.f90.obj:
	$(FC) -c $(FFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $<

# ********  rules for each platform supported *************

# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV  Windows VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

WIN_LIBS = user32.lib gdi32.lib comdlg32.lib shell32.lib

WIN_LIBDIR_64 = $(LIB_DIR)/lib_win_intel_64
SMV_LIBS_64 = $(WIN_LIBDIR_64)\glui.lib $(WIN_LIBDIR_64)\glut32.lib $(WIN_LIBDIR_64)\gd.lib $(WIN_LIBDIR_64)\jpeg.lib $(WIN_LIBDIR_64)\png.lib $(WIN_LIBDIR_64)\zlib.lib $(WIN_LIBDIR_64)\pthreads.lib

# ------------- ms_win_64 ----------------

ms_win_64 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
ms_win_64 : FFLAGS    = -O2 /iface:stdref /fpp
ms_win_64 : CFLAGS    = -O1 -D WIN32 -D _CONSOLE -D X64 -D GLEW_STATIC -D PTW32_STATIC_LIB $(SMV_TESTFLAG)
ms_win_64 : LFLAGS    = /F32768000 
ms_win_64 : CC        = cl
ms_win_64 : CPP       = cl
ms_win_64 : FC        = ifort
ms_win_64 : exe       = smokeview_win_$(SMV_TESTSTRING)64

ms_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(objwin) $(WIN_LIBS) $(SMV_LIBS_64)

# note: use /fpscomp:general option to allow smokeview to read files being currently written to by FDS
#       (instead of the non-standard SHARED keyword in a FORTRAN OPEN statement)

# ------------- intel_win_64_db ----------------

intel_win_64_db : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
intel_win_64_db : FFLAGS    = /Od /iface:stdref /fpp -D WIN32 /Z7 /nologo /fpscomp:general /debug:full /extend_source:132 /warn:unused /warn:nointerfaces /Qtrapuv /fp:strict /fp:except /traceback /check:all /stand:f08
intel_win_64_db : CFLAGS    = /Od /W5 /debug:full /Z7 -D WIN32 -D pp_INTEL -D _CONSOLE -D X64 -D GLEW_STATIC -D PTW32_STATIC_LIB $(SMV_TESTFLAG) $(GITINFO)
intel_win_64_db : LFLAGS    = /F32768000 /DEBUG
intel_win_64_db : CC        = icl
intel_win_64_db : CPP       = icl
intel_win_64_db : FC        = ifort
intel_win_64_db : exe       = smokeview_win_$(SMV_TESTSTRING)64_db

intel_win_64_db : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(objwin) $(WIN_LIBS) $(SMV_LIBS_64)

# ------------- intel_win_64 ----------------

intel_win_64 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
intel_win_64 : FFLAGS    = -O2 /iface:stdref /fpp -D WIN32 /fpscomp:general 
intel_win_64 : CFLAGS    = -O1 -D WIN32 -D pp_INTEL -D _CONSOLE -D X64 -D GLEW_STATIC -D PTW32_STATIC_LIB $(SMV_TESTFLAG) $(GITINFO)
intel_win_64 : LFLAGS    = /F32768000 
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : FC        = ifort
intel_win_64 : exe       = smokeview_win_$(SMV_TESTSTRING)64

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(objwin) $(WIN_LIBS) $(SMV_LIBS_64)

# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV  Linux VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

LINUX_SYSTEMLIBS = -lpthread -lX11 -lXmu -lGLU -lGL -lm

LINUX_INTEL_LIBS_64 = $(IFORT_COMPILER_LIB)/intel64/libifcore.a  $(IFORT_COMPILER_LIB)/intel64/libirc.a

# ------------- intel_linux_64 ----------------

intel_linux_64 : FFLAGS    = -O0 -traceback -m64 -static-intel -fpp -D pp_noappend
intel_linux_64 : CFLAGS    = -O0 -traceback -m64 -static-intel -D pp_LINUX -D pp_INTEL $(SMV_TESTFLAG) $(SMV_PROFILEFLAG) $(GITINFO)
intel_linux_64 : LFLAGS    = -m64 -static-intel
intel_linux_64 : CC        = icc 
intel_linux_64 : CPP       = icpc
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = smokeview_linux_$(SMV_TESTSTRING)64$(SMV_PROFILESTRING)

intel_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LIB_DIR)/lib_linux_intel_64 $(FULL_SMV_LIBS) $(LINUX_INTEL_LIBS_64) $(LINUX_SYSTEMLIBS)

# ------------- intel_linux_64_db ----------------

intel_linux_64_db : FFLAGS    = -O0 -m64 -static-intel -traceback -g -fpe0 -fltconsistency -WB -fpp -D pp_noappend -stand:f08
intel_linux_64_db : CFLAGS    = -O0 -g -m64 -static-intel -D pp_LINUX -D pp_INTEL $(SMV_PROFILEFLAG) -Wall -Wextra -check=stack,uninit -fp-stack-check -fp-trap-all=divzero,invalid,overflow -ftrapuv -Wuninitialized -Wunused-function -Wunused-variable $(GITINFO)
intel_linux_64_db : LFLAGS    = -m64 -static-intel
intel_linux_64_db : CC        = icc
intel_linux_64_db : CPP       = icpc
intel_linux_64_db : FC        = ifort
intel_linux_64_db : exe       = smokeview_linux_64_db$(SMV_PROFILESTRING)

intel_linux_64_db : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LIB_DIR)/lib_linux_intel_64 $(FULL_SMV_LIBS) $(LINUX_INTEL_LIBS_64) $(LINUX_SYSTEMLIBS)


# ------------- intel_linux_64_profile ----------------

intel_linux_64_profile : FFLAGS    = -O0 -m64 -static-intel -fpp -D pp_noappend -p
intel_linux_64_profile : CFLAGS    = -O0 -m64 -static-intel -D pp_LINUX -D pp_INTEL -p $(GITINFO)
intel_linux_64_profile : LFLAGS    = -m64 -static-intel
intel_linux_64_profile : CC        = icc 
intel_linux_64_profile : CPP       = icpc
intel_linux_64_profile : FC        = ifort
intel_linux_64_profile : exe       = smokeview_linux_64p

intel_linux_64_profile : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LIB_DIR)/lib_linux_intel_64 $(FULL_SMV_LIBS) $(LINUX_INTEL_LIBS_64) $(LINUX_SYSTEMLIBS) -p

# ------------- gnu_linux_64_db ----------------

gnu_linux_64_db : FFLAGS    = -O0 -m64 -x f95-cpp-input -D pp_GCC -D pp_noappend -ffree-form -frecord-marker=4
gnu_linux_64_db : CFLAGS    = -O0 -m64 -D pp_LINUX
gnu_linux_64_db : LFLAGS    = -m64
gnu_linux_64_db : CC        = gcc
gnu_linux_64_db : CPP       = g++
gnu_linux_64_db : FC        = gfortran
gnu_linux_64_db : exe       = smokeview_linux_$(SMV_TESTSTRING)64

gnu_linux_64_db : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L $(LIB_DIR)/lib_linux_gnu_64 $(FULL_SMV_LIBS) -lgfortran $(LINUX_SYSTEMLIBS)

# ------------- gnu_linux_64 ----------------

gnu_linux_64 : FFLAGS    = -O0 -m64 -x f95-cpp-input -D pp_GCC -D pp_noappend -ffree-form -frecord-marker=4
gnu_linux_64 : CFLAGS    = -O0 -m64 -D pp_LINUX
gnu_linux_64 : LFLAGS    = -m64
gnu_linux_64 : CC        = gcc
gnu_linux_64 : CPP       = g++
gnu_linux_64 : FC        = gfortran
gnu_linux_64 : exe       = smokeview_linux_$(SMV_TESTSTRING)64

gnu_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L $(LIB_DIR)/lib_linux_gnu_64 $(FULL_SMV_LIBS) -lgfortran $(LINUX_SYSTEMLIBS)

# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV   OSX   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

OSX_INTEL_LIBS = $(IFORT_COMPILER_LIB)/libifcoremt.a $(IFORT_COMPILER_LIB)/libimf.a $(IFORT_COMPILER_LIB)/libirc.a $(IFORT_COMPILER_LIB)/libsvml.a

# ------------- intel_osx_64 ----------------

intel_osx_64 : FFLAGS    = -O0 -m64 -static-intel -fpp -D pp_OSX -mmacosx-version-min=10.4
intel_osx_64 : CFLAGS    = -O0 -m64 -static-intel -D pp_OSX  -D pp_INTEL -mmacosx-version-min=10.4 $(SMV_TESTFLAG) $(SMV_PROFILEFLAG) $(GITINFO)
intel_osx_64 : LFLAGS    = -m64 -static-intel -framework OpenGL -framework GLUT -mmacosx-version-min=10.4
intel_osx_64 : CC        = icc 
intel_osx_64 : CPP       = icpc
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = smokeview_osx_$(SMV_TESTSTRING)64 $(SMV_PROFILESTRING)

intel_osx_64 : $(obj)
	icpc -o $(bin)/$(exe) $(LFLAGS) $(obj)  -L $(LIB_DIR)/lib_osx_intel_64 $(SMV_LIBS) $(OSX_INTEL_LIBS)

# ------------- intel_osx_64_fg ----------------

intel_osx_64_fg : FFLAGS    = -O0 -m64 -static-intel -fpp -D pp_OSX -mmacosx-version-min=10.4
intel_osx_64_fg : CFLAGS    = -O0 -m64 -I /shared/freeglut/include -static-intel -D pp_OSX  -D pp_INTEL -mmacosx-version-min=10.4 $(SMV_TESTFLAG) $(SMV_PROFILEFLAG) $(GITINFO)
intel_osx_64_fg : LFLAGS    = -m64 -static-intel -framework OpenGL -L/shared/freeglut/lib -lglut -mmacosx-version-min=10.4
intel_osx_64_fg : CC        = icc 
intel_osx_64_fg : CPP       = icpc
intel_osx_64_fg : FC        = ifort
intel_osx_64_fg : exe       = smokeview_osx_$(SMV_TESTSTRING)64_fg $(SMV_PROFILESTRING)

intel_osx_64_fg : $(obj)
	icpc -o $(bin)/$(exe) $(LFLAGS) $(obj)  -L $(LIB_DIR)/lib_osx_intel_64 $(SMV_LIBS) $(OSX_INTEL_LIBS)

# ------------- intel_osx_64_db ----------------

intel_osx_64_db : FFLAGS    = -O0 -m64 -static-intel -fpp -D pp_OSX -mmacosx-version-min=10.4 -stand:f08
intel_osx_64_db : CFLAGS    = -O0 -m64 -static-intel -D pp_OSX  -D pp_INTEL -mmacosx-version-min=10.4 $(SMV_TESTFLAG) $(SMV_PROFILEFLAG) $(GITINFO)
intel_osx_64_db : LFLAGS    = -m64 -static-intel -framework OpenGL -framework GLUT -mmacosx-version-min=10.4
intel_osx_64_db : CC        = icc 
intel_osx_64_db : CPP       = icpc
intel_osx_64_db : FC        = ifort
intel_osx_64_db : exe       = smokeview_osx_$(SMV_TESTSTRING)64_db $(SMV_PROFILESTRING)

intel_osx_64_db : $(obj)
	icpc -o $(bin)/$(exe) $(LFLAGS) $(obj)  -L $(LIB_DIR)/lib_osx_intel_64 $(SMV_LIBS) $(OSX_INTEL_LIBS)

# ------------- gnu_osx_64_db ----------------

gnu_osx_64_db : FFLAGS    = -O0 -m64 -x f95-cpp-input -D pp_GCC -D pp_noappend -ffree-form
gnu_osx_64_db : CFLAGS    = -O0 -m64 -D pp_OSX -Wno-write-strings
gnu_osx_64_db : LFLAGS    = -m64 -framework OpenGL -framework GLUT
gnu_osx_64_db : CC        = gcc
gnu_osx_64_db : CPP       = g++
gnu_osx_64_db : FC        = gfortran
gnu_osx_64_db : exe       = smokeview_osx_64_db

gnu_osx_64_db : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LIB_DIR)/lib_osx_gnu_64 $(SMV_LIBS) -L $(GLIBDIR) -lgfortran

# ------------- gnu_osx_64 ----------------

gnu_osx_64 : FFLAGS    = -O0 -m64 -x f95-cpp-input -D pp_GCC -D pp_noappend -ffree-form
gnu_osx_64 : CFLAGS    = -O0 -m64 -D pp_OSX -Wno-write-strings
gnu_osx_64 : LFLAGS    = -m64 -framework OpenGL -framework GLUT
gnu_osx_64 : CC        = gcc
gnu_osx_64 : CPP       = g++
gnu_osx_64 : FC        = gfortran
gnu_osx_64 : exe       = smokeview_osx_64

gnu_osx_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LIB_DIR)/lib_osx_gnu_64 $(SMV_LIBS) -L $(GLIBDIR) -lgfortran

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm -f *.o *.mod

#-------------- force compilation of string_util.c -----------

string_util.o:	.FORCE
string_util.obj:	.FORCE

.FORCE:
