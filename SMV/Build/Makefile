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

#----------------------------------------------------------------------------
# Should not need to edit lines below except to add or 'debug' target entries

# The variable VPATH defines the source code directory relative to the current directory

VPATH = $(SOURCE_DIR)/smokeview:$(SOURCE_DIR)/shared:$(SOURCE_DIR)/glew
bin = .

# BUILD_DIR defined directory containing this Makefile

BUILD_DIR = $(SOURCE_DIR)/../Build
LIB_DIR = $(BUILD_DIR)/LIBS
INC_DIR = $(SOURCE_DIR)/include2

# Definition of the object variables

obj = main.o menus.o IOscript.o IOshooter.o glui_devices.o glui_shooter.o glui_wui.o csphere.o colortimebar.o smv_endian.o camera.o IOembed.o IOhrr.o IOwui.o IOobject.o IOtour.o glui_display.o glui_tour.o getdatacolors.o smokeview.o output.o renderfile.o isobox.o getdatabounds.o readsmv.o scontour2d.o glui_3dsmoke.o glui_clip.o glui_stereo.o glui_blockedit.o glui_motion.o glui_bounds.o egz_stdio.o dmalloc.o interp.o assert.o compress3dc.o compress.o IOvolsmoke.o IOsmoke.o IOplot3d.o IOslice.o IOboundary.o IOpart.o IOzone.o IOiso.o IOtarget.o callbacks.o drawGeometry.o glui_colorbar.o skybox.o file_util.o string_util.o startup.o glui_trainer.o shaders.o unit.o threader.o histogram.o translate.o update.o viewports.o smv_geometry.o geometry.o showscene.o getdata.o getsizes.o glew.o
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

WIN_LIBDIR_32 = $(LIB_DIR)/lib_win_intel_32
WIN_LIBDIR_64 = $(LIB_DIR)/lib_win_intel_64
SMV_LIBS_32 = $(WIN_LIBDIR_32)\glui.lib $(WIN_LIBDIR_32)\glut32.lib $(WIN_LIBDIR_32)\gd.lib $(WIN_LIBDIR_32)\jpeg.lib $(WIN_LIBDIR_32)\png.lib $(WIN_LIBDIR_32)\zlib.lib $(WIN_LIBDIR_32)\pthreads.lib 
SMV_LIBS_64 = $(WIN_LIBDIR_64)\glui.lib $(WIN_LIBDIR_64)\glut32.lib $(WIN_LIBDIR_64)\gd.lib $(WIN_LIBDIR_64)\jpeg.lib $(WIN_LIBDIR_64)\png.lib $(WIN_LIBDIR_64)\zlib.lib $(WIN_LIBDIR_64)\pthreads.lib

# ------------- intel_win_32 ----------------

intel_win_32 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
intel_win_32 : FFLAGS    = -O2 -fpp -D pp_cvf /iface:stdref /arch:IA32
intel_win_32 : CFLAGS    = -O2 -D WIN32 -D pp_INTEL -D _CONSOLE  -D GLEW_STATIC -D PTW32_STATIC_LIB $(SMV_TESTFLAG) /arch:IA32
intel_win_32 : LFLAGS    = /F32768000 
intel_win_32 : CC        = icl
intel_win_32 : CPP       = icl
intel_win_32 : FC        = ifort
intel_win_32 : exe       = smokeview_win_$(SMV_TESTSTRING)32

intel_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(objwin) $(WIN_LIBS) $(SMV_LIBS_32)

# ------------- intel_win_64_db ----------------

intel_win_64_db : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
intel_win_64_db : FFLAGS    = -O2 -fpp -D pp_cvf -D X64 /iface:stdref
intel_win_64_db : CFLAGS    = -O1 -D WIN32 -D pp_INTEL -D _CONSOLE -D X64 -D GLEW_STATIC -D PTW32_STATIC_LIB $(SMV_TESTFLAG)
intel_win_64_db : LFLAGS    = /F32768000 
intel_win_64_db : CC        = icl
intel_win_64_db : CPP       = icl
intel_win_64_db : FC        = ifort
intel_win_64_db : exe       = smokeview_win_64_db

intel_win_64_db : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(objwin) $(WIN_LIBS) $(SMV_LIBS_64)

# ------------- intel_win_64 ----------------

# intel_win_64 : INC = -I ../../source/include3 -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
intel_win_64 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
intel_win_64 : FFLAGS    = -O2 -fpp -D pp_cvf -D X64 /iface:stdref
intel_win_64 : CFLAGS    = -O1 -D WIN32 -D pp_INTEL -D _CONSOLE -D X64 -D GLEW_STATIC -D PTW32_STATIC_LIB $(SMV_TESTFLAG)
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

LINUX_INTEL_LIBS_32 = $(IFORT_COMPILER_LIB)/ia32/libifcore.a  $(IFORT_COMPILER_LIB)/ia32/libirc.a
LINUX_INTEL_LIBS_64 = $(IFORT_COMPILER_LIB)/intel64/libifcore.a  $(IFORT_COMPILER_LIB)/intel64/libirc.a

# ------------- intel_linux_32 ----------------

intel_linux_32 : FFLAGS    = -O2 -m32 -static-intel -vec_report0 -fpp -D pp_noappend
intel_linux_32 : CFLAGS    = -O2 -m32 -static-intel -vec_report0 -D pp_LINUX -D pp_INTEL $(SMV_TESTFLAG) $(SMV_PROFILEFLAG)
intel_linux_32 : LFLAGS    = -m32 -static-intel
intel_linux_32 : CC        = icc 
intel_linux_32 : CPP       = icpc
intel_linux_32 : FC        = ifort
intel_linux_32 : exe       = smokeview_linux_$(SMV_TESTSTRING)32$(SMV_PROFILESTRING)

intel_linux_32 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LIB_DIR)/lib_linux_intel_32 $(FULL_SMV_LIBS) $(LINUX_INTEL_LIBS_32) -L/usr/lib $(LINUX_SYSTEMLIBS)

# ------------- intel_linux_64 ----------------

intel_linux_64 : FFLAGS    = -O2 -m64 -static-intel -vec_report0 -fpp -D pp_noappend
intel_linux_64 : CFLAGS    = -O2 -m64 -static-intel -vec_report0 -D pp_LINUX64 -D pp_INTEL $(SMV_TESTFLAG) $(SMV_PROFILEFLAG)
intel_linux_64 : LFLAGS    = -m64 -static-intel
intel_linux_64 : CC        = icc 
intel_linux_64 : CPP       = icpc
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = smokeview_linux_$(SMV_TESTSTRING)64$(SMV_PROFILESTRING)

intel_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LIB_DIR)/lib_linux_intel_64 $(FULL_SMV_LIBS) $(LINUX_INTEL_LIBS_64) $(LINUX_SYSTEMLIBS)

# ------------- intel_linux_64_db ----------------

intel_linux_64_db : FFLAGS    = -O0 -m64 -static-intel -traceback -g -fpe0 -fltconsistency -WB -vec_report0 -fpp -D pp_noappend
intel_linux_64_db : CFLAGS    = -O0 -g -m64 -static-intel -vec_report0 -D pp_LINUX64 -D pp_INTEL $(SMV_PROFILEFLAG)
intel_linux_64_db : LFLAGS    = -m64 -static-intel
intel_linux_64_db : CC        = icc
intel_linux_64_db : CPP       = icpc
intel_linux_64_db : FC        = ifort
intel_linux_64_db : exe       = smokeview_linux_64_db$(SMV_PROFILESTRING)

intel_linux_64_db : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LIB_DIR)/lib_linux_intel_64 $(FULL_SMV_LIBS) $(LINUX_INTEL_LIBS_64) $(LINUX_SYSTEMLIBS)


# ------------- intel_linux_64_profile ----------------

intel_linux_64_profile : FFLAGS    = -O2 -m64 -static-intel -vec_report0 -fpp -D pp_noappend -p
intel_linux_64_profile : CFLAGS    = -O2 -m64 -static-intel -vec_report0 -D pp_LINUX64 -D pp_INTEL -p
intel_linux_64_profile : LFLAGS    = -m64 -static-intel
intel_linux_64_profile : CC        = icc 
intel_linux_64_profile : CPP       = icpc
intel_linux_64_profile : FC        = ifort
intel_linux_64_profile : exe       = smokeview_linux_64p

intel_linux_64_profile : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LIB_DIR)/lib_linux_intel_64 $(FULL_SMV_LIBS) $(LINUX_INTEL_LIBS_64) $(LINUX_SYSTEMLIBS) -p

# ------------- gcc_linux ----------------

gcc_linux : FFLAGS    = -O2 -m64 -x f95-cpp-input -D pp_noappend -ffree-form -frecord-marker=4
gcc_linux : CFLAGS    = -O2 -m64 -D pp_LINUX64
gcc_linux : LFLAGS    = -m64
gcc_linux : CC        = gcc
gcc_linux : CPP       = g++
gcc_linux : FC        = gfortran
gcc_linux : exe       = smokeview_gcc_$(SMV_TESTSTRING)64

gcc_linux : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L $(LIB_DIR)/lib_linux_gcc_64 $(FULL_SMV_LIBS) -lgfortran $(LINUX_SYSTEMLIBS)

# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV   OSX   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

OSX_INTEL_LIBS = $(IFORT_COMPILER)/lib/libifcoremt.a $(IFORT_COMPILER)/lib/libimf.a $(IFORT_COMPILER)/lib/libirc.a $(IFORT_COMPILER)/lib/libsvml.a

# ------------- intel_osx_32 ----------------

intel_osx_32 : FFLAGS    = -O2 -m32 -static-intel -vec_report0 -fpp -D pp_OSX
intel_osx_32 : CFLAGS    = -O2 -m32 -static-intel -vec_report0 -D pp_OSX -D pp_INTEL $(SMV_TESTFLAG) $(SMV_PROFILEFLAG)
intel_osx_32 : LFLAGS    = -m32 -static-intel -mmacosx-version-min=10.4 -framework OpenGL -framework GLUT
intel_osx_32 : CC        = icc 
intel_osx_32 : CPP       = icpc
intel_osx_32 : FC        = ifort
intel_osx_32 : exe       = smokeview_osx_$(SMV_TESTSTRING)32$(SMV_PROFILESTRING)

intel_osx_32 : $(obj)
	icpc -o $(bin)/$(exe) $(LFLAGS) $(obj) -L $(LIB_DIR)/lib_osx_intel_32 $(SMV_LIBS) $(OSX_INTEL_LIBS)

# ------------- intel_osx_64 ----------------

intel_osx_64 : FFLAGS    = -O2 -m64 -static-intel -vec_report0 -fpp -D pp_OSX64 -mmacosx-version-min=10.4
intel_osx_64 : CFLAGS    = -O2 -m64 -static-intel -vec_report0 -D pp_OSX64  -D pp_INTEL -mmacosx-version-min=10.4 $(SMV_TESTFLAG) $(SMV_PROFILEFLAG)
intel_osx_64 : LFLAGS    = -m64 -static-intel -framework OpenGL -framework GLUT -mmacosx-version-min=10.4
intel_osx_64 : CC        = icc 
intel_osx_64 : CPP       = icpc
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = smokeview_osx_$(SMV_TESTSTRING)64 $(SMV_PROFILESTRING)

intel_osx_64 : $(obj)
	icpc -o $(bin)/$(exe) $(LFLAGS) $(obj)  -L $(LIB_DIR)/lib_osx_intel_64 $(SMV_LIBS) $(OSX_INTEL_LIBS)

# ------------- gcc_osx ----------------

gcc_osx : FFLAGS    = -O2 -m64 -x f95-cpp-input -D pp_noappend -ffree-form
gcc_osx : CFLAGS    = -O2 -m64 -D pp_OSX -Wno-write-strings
gcc_osx : LFLAGS    = -m64 -framework OpenGL -framework GLUT
gcc_osx : CC        = gcc
gcc_osx : CPP       = g++
gcc_osx : FC        = gfortran
gcc_osx : exe       = smokeview_gcc_64

gcc_osx : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LIB_DIR)/lib_osx_gcc_64 $(SMV_LIBS) -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin10/4.6.2 -lgfortran

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm -f *.o *.mod
