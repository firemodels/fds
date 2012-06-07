# $Date$ 
# $Revision$
# $Author$

#*** Smokeview, Version 5, Makefile ***

#To use this Makefile, you will need to provide a make target at the command line.
#For example, to make the 32 bit version in the current directory with the Intel
#compiler, you would need to type 'make intel_linux_32' at the prompt.
#To make it in a sub-directory type 'make -f ../Makefile intel_linux_32'

# use the following line if building Smokeview in a directory 'below' Makefile

SOURCE_DIR = ../../source/

# use FILTERC and FILTERF to customize smokeview builds (ie -D pp_BETA)

FILTERC =
FILTERF =

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

obj = main.o menus.o IOscript.o IOshooter.o glui_devices.o glui_message.o glui_shooter.o glui_wui.o csphere.o colortimebar.o smv_endian.o camera.o IOembed.o IOhrr.o IOwui.o IOobject.o IOtour.o glui_labels.o glui_tour.o getdatacolors.o smokeview.o output.o renderfile.o isobox.o getdatabounds.o readsmv.o scontour2d.o glui_3dsmoke.o glui_clip.o glui_stereo.o glui_blockedit.o glui_motion.o glui_bounds.o egz_stdio.o dmalloc.o interp.o assert.o compress3dc.o compress.o IOvolsmoke.o IOsmoke.o IOplot3d.o IOslice.o IOboundary.o IOpart.o IOzone.o IOiso.o IOtarget.o callbacks.o drawGeometry.o glui_colorbar.o skybox.o file_util.o string_util.o startup.o glui_trainer.o shaders.o unit.o threader.o histogram.o translate.o update.o viewports.o smv_geometry.o geometry.o showscene.o getdata.o getsizes.o glew.o
objwin = $(obj:.o=.obj)

# include file directory locations (should be the same for all builds)

INC = -I $(INC_DIR)/glut -I $(INC_DIR)/glui -I $(INC_DIR)/gd -I $(SOURCE_DIR)/shared -I $(SOURCE_DIR)/smokeview -I $(SOURCE_DIR)/glew

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


# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV  Windows VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

WINLIBS = opengl32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib

LIBWIN32 = $(LIB_DIR)/lib_win_intel_32
SMVLIBS32 = $(LIBWIN32)\pthreadVC.lib $(LIBWIN32)\glui.lib $(LIBWIN32)\glut32.lib $(LIBWIN32)\gd.lib $(LIBWIN32)\jpeg.lib $(LIBWIN32)\png.lib $(LIBWIN32)\zlib.lib

LIBWIN64 = $(LIB_DIR)/lib_win_intel_64
SMVLIBS64 = $(LIBWIN64)\glui.lib $(LIBWIN64)\glut32.lib $(LIBWIN64)\gd.lib $(LIBWIN64)\jpeg.lib $(LIBWIN64)\png.lib $(LIBWIN64)\zlib.lib $(LIBWIN64)\pthreadVC2_x64.lib


# ------------- ms_win_32 ----------------

ms_win_32 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
ms_win_32 : FFLAGS    = -O2 -fpp -D pp_cvf /iface:stdref
ms_win_32 : CFLAGS    = -O2 -D WIN32 -D _CONSOLE  -D GLEW_STATIC
ms_win_32 : LFLAGS    = /F32768000 
ms_win_32 : CC        = cl
ms_win_32 : CPP       = cl
ms_win_32 : FC        = ifort
ms_win_32 : exe       = smokeview_win_32

ms_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(objwin) $(WINLIBS) $(SMVLIBS32)

# ------------- gcc_win_32 ----------------

gcc_win_32 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
gcc_win_32 : FFLAGS    = -O2 -m32 -D pp_noappend -ffree-form
gcc_win_32 : CFLAGS    = -O2 -m32 -D WIN32 -D GLEW_STATIC
gcc_win_32 : CC        = gcc
gcc_win_32 : CPP       = g++
gcc_win_32 : FC        = gfortran
gcc_win_32 : exe       = smokeview_win_32

gcc_win_32 : $(obj)
	$(CPP) -m32 -o $(bin)/$(exe) $(obj) $(SMVLIBS32)

# ------------- intel_win_32 ----------------

intel_win_32 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
intel_win_32 : FFLAGS    = -O2 -fpp -D pp_cvf /iface:stdref /arch:IA32
intel_win_32 : CFLAGS    = -O2 -D WIN32 -D _CONSOLE  -D GLEW_STATIC /arch:IA32
intel_win_32 : LFLAGS    = /F32768000 
intel_win_32 : CC        = icl
intel_win_32 : CPP       = icl
intel_win_32 : FC        = ifort
intel_win_32 : exe       = smokeview_win_32

intel_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(objwin) $(WINLIBS) $(SMVLIBS32)

# ------------- intel_win_test_32 ----------------

intel_win_test_32 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
intel_win_test_32 : FFLAGS    = -O2 -fpp -D pp_cvf -D pp_BETA   /iface:stdref /arch:IA32
intel_win_test_32 : CFLAGS    = -O2 -D WIN32 -D _CONSOLE -D pp_BETA -D GLEW_STATIC /arch:IA32
intel_win_test_32 : LFLAGS    = /F32768000 
intel_win_test_32 : CC        = icl
intel_win_test_32 : CPP       = icl
intel_win_test_32 : FC        = ifort
intel_win_test_32 : exe       = smokeview_win_test_32

intel_win_test_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(objwin) $(WINLIBS) $(SMVLIBS32)

# ------------- intel_win_64 ----------------

intel_win_64 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
intel_win_64 : FFLAGS    = -O2 -fpp -D pp_cvf -D X64 /iface:stdref
intel_win_64 : CFLAGS    = -O2 -D WIN32 -D _CONSOLE -D X64 -D GLEW_STATIC
intel_win_64 : LFLAGS    = /F32768000 
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : FC        = ifort
intel_win_64 : exe       = smokeview_win_64

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(objwin) $(WINLIBS) $(SMVLIBS64)

# ------------- intel_win_test_64 ----------------

intel_win_test_64 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
intel_win_test_64 : FFLAGS    = -O2 -traceback -fpp -D pp_cvf -D pp_BETA -D X64 /iface:stdref
intel_win_test_64 : CFLAGS    = -O2 -traceback -D WIN32 -D _CONSOLE -D X64 -D pp_BETA -D GLEW_STATIC
intel_win_test_64 : LFLAGS    = /F32768000
intel_win_test_64 : CC        = icl
intel_win_test_64 : CPP       = icl
intel_win_test_64 : FC        = ifort
intel_win_test_64 : exe       = smokeview_win_test_64

intel_win_test_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(objwin) $(WINLIBS) $(SMVLIBS64)

# ------------- ms_win_64 ----------------

ms_win_64 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
ms_win_64 : FFLAGS    = -O2 -fpp -D pp_cvf -D X64 /iface:stdref
ms_win_64 : CFLAGS    = -O2 -D WIN32 -D _CONSOLE -D X64 -D GLEW_STATIC
ms_win_64 : LFLAGS    = /F32768000 
ms_win_64 : CC        = cl
ms_win_64 : CPP       = cl
ms_win_64 : FC        = ifort
ms_win_64 : exe       = smokeview_win_64

ms_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(objwin) $(WINLIBS) $(SMVLIBS64)

# ------------- ms_win_test_64 ----------------

ms_win_test_64 : INC = -I ../../source/include -I ../../source/shared -I ../../source/smokeview -I ../../source/glew
ms_win_test_64 : FFLAGS    = -O2 -fpp -D pp_cvf -D X64 -D pp_BETA /iface:stdref
ms_win_test_64 : CFLAGS    = -O2 -D WIN32 -D _CONSOLE -D X64 -D pp_BETA -D GLEW_STATIC
ms_win_test_64 : LFLAGS    = /F32768000 
ms_win_test_64 : CC        = cl
ms_win_test_64 : CPP       = cl
ms_win_test_64 : FC        = ifort
ms_win_test_64 : exe       = smokeview_win_64

ms_win_test_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(LFLAGS) $(objwin) $(WINLIBS) $(SMVLIBS64)

# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV  Linux VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

LINUXSYSTEMLIBS = -lpthread -lXmu -lGLU -lGL -lm

LINUXFORTLIBDIR32  = $IFORT_COMPILER/lib/ia32
LINUXLIBDIR32 = $(LIB_DIR)/lib_linux_intel_32

LINUXFORTLIBDIR64  = $IFORT_COMPILER/lib/intel64
LINUXLIBDIR64 = $(LIB_DIR)/lib_linux_intel_64

LINUXSMVLIBS = -lglui -lglut -lgd -ljpeg -lpng -lz
LINUXINTELLIBS = -lifcore  -lintlc


# ------------- intel_linux_32 ----------------

intel_linux_32 : FORTLIBDIR  = $(LINUXFORTLIBDIR32)
intel_linux_32 : FFLAGS    = -O2 -m32 -static-intel -vec_report0 -fpp -D pp_noappend
intel_linux_32 : CFLAGS    = -O2 -m32 -static-intel -vec_report0 -D pp_LINUX
intel_linux_32 : LFLAGS    = -m32 -static-intel
intel_linux_32 : CC        = icc 
intel_linux_32 : CPP       = icpc
intel_linux_32 : FC        = ifort
intel_linux_32 : exe       = smokeview_linux_32

intel_linux_32 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LINUXLIBDIR32) -L $(FORTLIBDIR) $(LINUXSMVLIBS) $(LINUXINTELLIBS) -L/usr/lib $(LINUXSYSTEMLIBS)

# ------------- intel_linux_test_32 ----------------

intel_linux_test_32 : FFLAGS    = -O2 -m32 -static-intel -vec_report0 -fpp -D pp_noappend
intel_linux_test_32 : CFLAGS    = -O2 -m32 -static-intel -vec_report0 -D pp_BETA -D pp_LINUX
intel_linux_test_32 : LFLAGS    = -m32 -static-intel
intel_linux_test_32 : CC        = icc 
intel_linux_test_32 : CPP       = icpc
intel_linux_test_32 : FC        = ifort
intel_linux_test_32 : exe       = smokeview_linux_test_32

intel_linux_test_32 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LINUXLIBDIR32) -L $(LINUXFORTLIBDIR32) $(LINUXSMVLIBS) $(LINUXINTELLIBS) -L/usr/lib $(LINUXSYSTEMLIBS)

# ------------- intel_linux_32_dbg ----------------

intel_linux_32_dbg : FFLAGS    = -O0 -m32 -static-intel -traceback -g -fpe0 -fltconsistency -WB -vec_report0 -fpp -D pp_noappend
intel_linux_32_dbg : CFLAGS    = -O0 -g -m32 -static-intel -vec_report0 -D pp_LINUX
intel_linux_32_dbg : LFLAGS    = -m32 -static-intel
intel_linux_32_dbg : CC        = icc 
intel_linux_32_dbg : CPP       = icpc
intel_linux_32_dbg : FC        = ifort
intel_linux_32_dbg : exe       = smokeview_linux_32_dbg

intel_linux_32_dbg : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LINUXLIBDIR32) -L $(LINUXFORTLIBDIR32) $(LINUXSMVLIBS) $(LINUXINTELLIBS) $(LINUXSYSTEMLIBS)

# ------------- intel_linux_64 ----------------

intel_linux_64 : FFLAGS    = -O2 -m64 -static-intel -vec_report0 -fpp -D pp_noappend
intel_linux_64 : CFLAGS    = -O2 -m64 -static-intel -vec_report0 -D pp_LINUX64
intel_linux_64 : LFLAGS    = -m64 -static-intel
intel_linux_64 : CC        = icc 
intel_linux_64 : CPP       = icpc
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = smokeview_linux_64

intel_linux_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LINUXLIBDIR64) -L $(LINUXFORTLIBDIR64) $(LINUXSMVLIBS) $(LINUXINTELLIBS) $(LINUXSYSTEMLIBS)

# ------------- intel_linux_test_64 ----------------

intel_linux_test_64 : FFLAGS    = -O2 -m64 -static-intel -traceback -vec_report0 -fpp -D pp_noappend
intel_linux_test_64 : CFLAGS    = -O2 -m64 -static-intel -traceback -vec_report0 -D pp_LINUX64 -D pp_BETA
intel_linux_test_64 : LFLAGS    = -m64 -static-intel -traceback
intel_linux_test_64 : CC        = icc 
intel_linux_test_64 : CPP       = icpc
intel_linux_test_64 : FC        = ifort
intel_linux_test_64 : exe       = smokeview_linux_test_64

intel_linux_test_64 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LINUXLIBDIR64) -L $(LINUXFORTLIBDIR64) $(LINUXSMVLIBS) $(LINUXINTELLIBS) $(LINUXSYSTEMLIBS)

# ------------- intel_linux_test_64_dbg ----------------

intel_linux_test_64_dbg : FFLAGS    = -O0 -m64 -static-intel -traceback -g -fpe0 -fltconsistency -WB -vec_report0 -fpp -D pp_noappend
intel_linux_test_64_dbg : CFLAGS    = -O0 -g -m64 -static-intel -vec_report0 -D pp_LINUX64 -D pp_BETA
intel_linux_test_64_dbg : LFLAGS    = -m64 -static-intel
intel_linux_test_64_dbg : CC        = icc 
intel_linux_test_64_dbg : CPP       = icpc
intel_linux_test_64_dbg : FC        = ifort
intel_linux_test_64_dbg : exe       = smokeview_linux_test_64_dbg

intel_linux_test_64_dbg : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LINUXLIBDIR64) -L $(LINUXFORTLIBDIR64) $(LINUXSMVLIBS) $(LINUXINTELLIBS) $(LINUXSYSTEMLIBS)


# ------------- intel_linux_64_dbg ----------------

intel_linux_64_dbg : FFLAGS    = -O0 -m64 -static-intel -traceback -g -fpe0 -fltconsistency -WB -vec_report0 -fpp -D pp_noappend
intel_linux_64_dbg : CFLAGS    = -O0 -g -m64 -static-intel -vec_report0 -D pp_LINUX64
intel_linux_64_dbg : LFLAGS    = -m64 -static-intel
intel_linux_64_dbg : CC        = icc 
intel_linux_64_dbg : CPP       = icpc
intel_linux_64_dbg : FC        = ifort
intel_linux_64_dbg : exe       = smokeview_linux_64_dbg

intel_linux_64_dbg : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LINUXLIBDIR64) -L $(LINUXFORTLIBDIR64) $(LINUXSMVLIBS) $(LINUXINTELLIBS) $(LINUXSYSTEMLIBS)

# ------------- intel_linux_64_profile ----------------

intel_linux_64_profile : FFLAGS    = -O2 -m64 -static-intel -vec_report0 -fpp -D pp_noappend -p
intel_linux_64_profile : CFLAGS    = -O2 -m64 -static-intel -vec_report0 -D pp_LINUX64 -p
intel_linux_64_profile : LFLAGS    = -m64 -static-intel
intel_linux_64_profile : CC        = icc 
intel_linux_64_profile : CPP       = icpc
intel_linux_64_profile : FC        = ifort
intel_linux_64_profile : exe       = smokeview_linux_64p

intel_linux_64_profile : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(LINUXLIBDIR64) -L $(LINUXFORTLIBDIR64) $(LINUXSMVLIBS) $(LINUXINTELLIBS) $(LINUXSYSTEMLIBS) -p

# ------------- gcc_linux_32 ----------------

gcc_linux_32 : FORTLIBDIR  = $(LINUXFORTLIBDIR32)
gcc_linux_32 : FFLAGS    = -O2 -m32 -D pp_noappend -ffree-form
gcc_linux_32 : CFLAGS    = -O2 -m32 -D pp_LINUX
gcc_linux_32 : LFLAGS    = -m32
gcc_linux_32 : CC        = gcc
gcc_linux_32 : CPP       = g++
gcc_linux_32 : FC        = gfortran
gcc_linux_32 : exe       = smokeview_linux_32

gcc_linux_32 : $(obj)
	$(CPP) -o $(bin)/$(exe) $(obj) $(LFLAGS) -L$(IMAGELIBDIR) -L $(LIB_DIR)/lib_linux_gcc_32 $(LINUXSMVLIBS) -lgfortran $(LINUXSYSTEMLIBS)

# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV   OSX   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

OSXINTELLIB32 = $(LIB_DIR)/lib_osx_intel_32
OSXSMVLIBS32 = -lz $(OSXINTELLIB32)/libpng.a $(OSXINTELLIB32)/libgd.a $(OSXINTELLIB32)/libglui.a

OSXINTELLIB64 = $(LIB_DIR)/lib_osx_intel_64
OSXSMVLIBS64 = -lz $(OSXINTELLIB64)/libpng.a $(OSXINTELLIB64)/libgd.a $(OSXINTELLIB64)/libglui.a

OSXFORTLIBDIR32  = $(IFORT_COMPILER2)/lib
OSXINTELLIBS32 = $(OSXFORTLIBDIR32)/libifcoremt.a $(OSXFORTLIBDIR32)/libimf.a $(OSXFORTLIBDIR32)/libirc.a $(OSXFORTLIBDIR32)/libsvml.a

OSXFORTLIBDIR64  = $(IFORT_COMPILER2)/lib
OSXINTELLIBS32 = $(OSXFORTLIBDIR64)/libifcoremt.a $(OSXFORTLIBDIR64)/libimf.a $(OSXFORTLIBDIR64)/libirc.a $(OSXFORTLIBDIR64)/libsvml.a

# ------------- intel_osx_32 ----------------

intel_osx_32 : FFLAGS    = -O2 -m32 -vec_report0 -fpp -D pp_OSX
intel_osx_32 : CFLAGS    = -O2 -m32 -vec_report0 -D pp_OSX
intel_osx_32 : LFLAGS    = -m32 -static-intel -mmacosx-version-min=10.4 -framework OpenGL -framework GLUT
intel_osx_32 : CC        = icc 
intel_osx_32 : CPP       = icpc
intel_osx_32 : FC        = ifort
intel_osx_32 : exe       = smokeview_osx_32

intel_osx_32 : $(obj)
	icpc -o $(bin)/$(exe) $(LFLAGS) $(obj)  $(OSXSMVLIBS32) $(OSXINTELLIBS32)

# ------------- intel_osx_test_32 ----------------

intel_osx_test_32 : FFLAGS    = -O2 -m32 -vec_report0 -fpp -D pp_OSX -D pp_BETA
intel_osx_test_32 : CFLAGS    = -O2 -m32 -vec_report0 -D pp_OSX -D pp_BETA
intel_osx_test_32 : LFLAGS    = -m32 -static-intel -mmacosx-version-min=10.4 -framework OpenGL -framework GLUT
intel_osx_test_32 : CC        = icc 
intel_osx_test_32 : CPP       = icpc
intel_osx_test_32 : FC        = ifort
intel_osx_test_32 : exe       = smokeview_osx_test_32

intel_osx_test_32 : $(obj)
	icpc -o $(bin)/$(exe) $(LFLAGS) $(obj) $(OSXSMVLIBS32) $(OSXINTELLIBS32)

# ------------- intel_osx_64 ----------------

intel_osx_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_OSX64 -mmaccox-version-min=10.4
intel_osx_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_OSX64 -mmacosx-version-min=10.4
intel_osx_64 : LFLAGS    = -m64 -static-intel -framework OpenGL -framework GLUT -mmacosx-version-min=10.4
intel_osx_64 : CC        = icc 
intel_osx_64 : CPP       = icpc
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = smokeview_osx_64

intel_osx_64 : $(obj)
	icpc -o $(bin)/$(exe) $(LFLAGS) $(obj)  $(OSXSMVLIBS64) $(OSXINTELLIBS64)

# ------------- intel_osx_test_64 ----------------

intel_osx_test_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_OSX64 -D pp_BETA
intel_osx_test_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_OSX64 -D pp_BETA
intel_osx_test_64 : LFLAGS    = -m64 -static-intel -framework OpenGL -framework GLUT -mmacosx-version-min=10.4
intel_osx_test_64 : CC        = icc 
intel_osx_test_64 : CPP       = icpc
intel_osx_test_64 : FC        = ifort
intel_osx_test_64 : exe       = smokeview_osx_test_64

intel_osx_test_64 : $(obj)
	icpc -o $(bin)/$(exe) $(LFLAGS) $(obj) $(OSXSMVLIBS64) $(OSXINTELLIBS64)

# ------------- intel_osx_test_64_db ----------------

intel_osx_test_64_dbg : FFLAGS    = -O0 -m64 -traceback -g -fpe0 -fltconsistency -WB -vec_report0 -fpp -D pp_BETA -D pp_OSX64 -D pp_noappend -D _DEBUG
intel_osx_test_64_dbg : CFLAGS    = -O0 -g -m64 -vec_report0 -D pp_OSX64 -D pp_BETA -D _DEBUG
intel_osx_test_64_dbg : LFLAGS    = -m64 -static-intel -framework OpenGL -framework GLUT -mmacosx-version-min=10.4
intel_osx_test_64_dbg : CC        = icc 
intel_osx_test_64_dbg : CPP       = icpc
intel_osx_test_64_dbg : FC        = ifort
intel_osx_test_64_dbg : exe       = smokeview_osx_test_64_dbg

intel_osx_test_64_dbg : $(obj)
	icpc -o $(bin)/$(exe) $(LFAGS) $(obj) $(OSXSMVLIBS64) $(OSXINTELLIBS64)


# ------------- make depend ----------------
# ------------- make depend ----------------

depend : FFLAGS    = -MM 
depend : CFLAGS    = -MM
depend : CC        = gcc
depend : CPP       = g++
depend : FC        = gfortran
depend : exe       = smokeview_linux_32

depend : $(obj)

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm -f *.o
