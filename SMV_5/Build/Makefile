#*** Smokeview, Version 5, Makefile ***

#To use this Makefile, you will need to provide a make target at the command line.
#For example, to make the 32 bit version in the current directory with the Intel
#compiler, you would need to type 'make intel_linux_32' at the prompt.
#To make it in a sub-directory type 'make -f ../Makefile intel_linux_32'

# use the following line if building Smokeview in the directory containing this Makefile
#SOURCE_DIR = ../source/

# use the following line if building Smokeview in a directory 'below' Makefile

SOURCE_DIR = ../../source/

# use FILTERC and FILTERF to customize smokeview builds (ie -D pp_BETA)

FILTERC =
FILTERF =

#----------------------------------------------------------------------------
# Should not need to edit lines below except to add or 'debug' target entries

# The variable VPATH defines the source code directory relative to the current directory

VPATH = $(SOURCE_DIR)/smokeview:$(SOURCE_DIR)/shared
bin = $(SOURCE_DIR)/../bin

# BUILD_DIR defined directory containing this Makefile

BUILD_DIR = $(SOURCE_DIR)/../Build
LIB_DIR = $(BUILD_DIR)/LIBS
INC_DIR = $(SOURCE_DIR)/include2

# SYSTEM libraries

SYSTEMLIBS = -lpthread -lXmu -lGLU -lGL -lm

# FORTRAN Library directories for various platforms

#OSX32FORTLIBDIR  = /opt/intel/fc/10.1.006/lib
OSX32FORTLIBDIR  = /opt/intel/Compiler/11.1/080/lib
OSX64FORTLIBDIR  = /opt/intel/Compiler/11.1/080/lib
#LINUX32FORTLIBDIR  = /opt/intel/fc/10.1.021/lib
LINUX32FORTLIBDIR  = /opt/intel/11/lib/ia32
LINUX64FORTLIBDIR  = /opt/intel/11/lib/intel64

# Definition of the object variables

obj = main.o IOscript.o IOshooter.o glui_message.o glui_shooter.o glui_wui.o csphere.o colorbar.o endian.o camera.o IOhrr.o IOwui.o IOobject.o IOtour.o glui_labels.o glui_tour.o smoothlabel.o getdatacolors.o smokeview.o output.o renderfile.o isobox.o getdatabounds.o readsmv.o scontour2d.o sv_api.o glui_3dsmoke.o glui_clip.o glui_stereo.o glui_blockedit.o glui_motion.o glui_bounds.o egz_stdio.o parseobst.o dmalloc.o assert.o compress3dc.o IOsmoke.o IOplot3d.o IOslice.o IOboundary.o IOpart.o IOzone.o IOiso.o IOtarget.o callbacks.o drawGeometry.o glui_colorbar.o skybox.o blockage_test.o segtri.o string_util.o startup.o glui_trainer.o shaders.o unit.o getdata.o getsizes.o histogram.o
objwin = $(obj:.o=.obj)

# include file directory locations (should be the same for all builds)

INC = -I $(INC_DIR)/glut -I $(INC_DIR)/glui -I $(INC_DIR)/gd -I $(SOURCE_DIR)/shared

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

intel_win_64 : INC = -I ../../source/include -I ../../source/shared
intel_win_64 : IMAGELIBDIR = $(LIB_DIR)/lib_win_intel_64
intel_win_64 : WINLIBS = kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib
intel_win_64 : FFLAGS    = -O2 -fpp -D pp_PART5 -D pp_cvf /Qvec_report0
intel_win_64 : CFLAGS    = -O2 -D WIN32 -D _CONSOLE -D pp_WIN_INTEL -D X64 /Qvec_report0
intel_win_64 : CC        = icl
intel_win_64 : CPP       = icl
intel_win_64 : FC        = ifort
intel_win_64 : exe       = smv5_win_64

intel_win_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) opengl32.lib $(WINLIBS) $(IMAGELIBDIR)\glui.lib $(IMAGELIBDIR)\glut32.lib $(IMAGELIBDIR)\gd.lib $(IMAGELIBDIR)\jpeg.lib $(IMAGELIBDIR)\png.lib $(IMAGELIBDIR)\zlib.lib $(IMAGELIBDIR)/smokeview_fort_64.lib


# ------------- Win 64 ----------------
#  64 bit windows

intel_win_test_64 : INC = -I ../../source/include -I ../../source/shared
intel_win_test_64 : IMAGELIBDIR = $(LIB_DIR)/lib_win_intel_64
intel_win_test_64 : WINLIBS = kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib
intel_win_test_64 : FFLAGS    = -O2 -fpp -D pp_PART5 -D pp_cvf -D pp_BETA /Qvec_report0
intel_win_test_64 : CFLAGS    = -O2 -D WIN32 -D _CONSOLE -D pp_WIN_INTEL -D X64 -D pp_BETA /Qvec_report0
intel_win_test_64 : CC        = icl
intel_win_test_64 : CPP       = icl
intel_win_test_64 : FC        = ifort
intel_win_test_64 : exe       = smv5_win_test_64

intel_win_test_64 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) opengl32.lib $(WINLIBS) $(IMAGELIBDIR)\glui.lib $(IMAGELIBDIR)\glut32.lib $(IMAGELIBDIR)\gd.lib $(IMAGELIBDIR)\jpeg.lib $(IMAGELIBDIR)\png.lib $(IMAGELIBDIR)\zlib.lib $(IMAGELIBDIR)/smokeview_fort_64.lib

# ------------- Win 32 ----------------
#  32 bit windows

intel_win_32 : INC = -I ../../source/include -I ../../source/shared
intel_win_32 : IMAGELIBDIR = $(LIB_DIR)/lib_win_intel_32
intel_win_32 : WINLIBS = kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib
intel_win_32 : FFLAGS    = -O2 -fpp -D pp_cvf -D pp_PART5 /Qvec_report0
intel_win_32 : CFLAGS    = -O2 -D WIN32 -D _CONSOLE -D pp_WIN_INTEL /Qvec_report0
intel_win_32 : CC        = icl
intel_win_32 : CPP       = icl
intel_win_32 : FC        = ifort
intel_win_32 : exe       = smv5_win_32

intel_win_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) opengl32.lib $(WINLIBS) $(IMAGELIBDIR)\pthreadVC.lib $(IMAGELIBDIR)\glui.lib $(IMAGELIBDIR)\glut32.lib $(IMAGELIBDIR)\glew32.lib $(IMAGELIBDIR)\gd.lib $(IMAGELIBDIR)\jpeg.lib $(IMAGELIBDIR)\png.lib $(IMAGELIBDIR)\zlib.lib


# ------------- Win test 32 ----------------
#  32 bit windows

intel_win_test_32 : INC = -I ../../source/include -I ../../source/shared
intel_win_test_32 : IMAGELIBDIR = $(LIB_DIR)/lib_win_intel_32
intel_win_test_32 : WINLIBS = kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib
intel_win_test_32 : FFLAGS    = -O2 -fpp -D pp_cvf -D pp_PART5  -D pp_BETA /Qvec_report0
intel_win_test_32 : CFLAGS    = -O2 -D WIN32 -D _CONSOLE -D pp_WIN_INTEL  -D pp_BETA /Qvec_report0
intel_win_test_32 : CC        = icl
intel_win_test_32 : CPP       = icl
intel_win_test_32 : FC        = ifort
intel_win_test_32 : exe       = smv5_win_test_32

intel_win_test_32 : $(objwin)
	$(CPP) -o $(bin)/$(exe) $(objwin) opengl32.lib $(WINLIBS) $(IMAGELIBDIR)\pthreadVC.lib $(IMAGELIBDIR)\glui.lib $(IMAGELIBDIR)\glut32.lib $(IMAGELIBDIR)\glew32.lib $(IMAGELIBDIR)\gd.lib $(IMAGELIBDIR)\jpeg.lib $(IMAGELIBDIR)\png.lib $(IMAGELIBDIR)\zlib.lib

# ------------- OSX test 32 ----------------
#  OSX running snow leopard (bluesky)

intel_osx_test_32 : IMAGELIBDIR = $(LIB_DIR)/lib_osx_intel_32
intel_osx_test_32 : FORTLIBDIR  = $(OSX32FORTLIBDIR)
intel_osx_test_32 : FFLAGS    = -O2 -m32 -vec_report0 -fpp -D pp_OSX -D pp_BETA
intel_osx_test_32 : CFLAGS    = -O2 -m32 -vec_report0 -D pp_OSX -D pp_BETA
intel_osx_test_32 : CC        = icc 
intel_osx_test_32 : CPP       = icpc
intel_osx_test_32 : FC        = ifort
intel_osx_test_32 : exe       = smv5_osx_test_32

intel_osx_test_32 : $(obj)
	icpc -m32 -o $(bin)/$(exe) $(obj) -mmacosx-version-min=10.4 -framework OpenGL -framework GLUT -lz $(IMAGELIBDIR)/libpng.a $(IMAGELIBDIR)/libgd.a $(IMAGELIBDIR)/libglui.a $(FORTLIBDIR)/libifcoremt.a $(FORTLIBDIR)/libimf.a $(FORTLIBDIR)/libirc.a $(FORTLIBDIR)/libsvml.a

# ------------- OSX 32 ----------------
#  OSX running snow leopard (bluesky)

intel_osx_32 : IMAGELIBDIR = $(LIB_DIR)/lib_osx_intel_32
intel_osx_32 : FORTLIBDIR  = $(OSX32FORTLIBDIR)
intel_osx_32 : FFLAGS    = -O2 -m32 -vec_report0 -fpp -D pp_OSX
intel_osx_32 : CFLAGS    = -O2 -m32 -vec_report0 -D pp_OSX
intel_osx_32 : CC        = icc 
intel_osx_32 : CPP       = icpc
intel_osx_32 : FC        = ifort
intel_osx_32 : exe       = smv5_osx_32

intel_osx_32 : $(obj)
	icpc -m32 -o $(bin)/$(exe) $(obj)  -mmacosx-version-min=10.4 -framework OpenGL -framework GLUT -lz $(IMAGELIBDIR)/libpng.a $(IMAGELIBDIR)/libgd.a $(IMAGELIBDIR)/libglui.a $(FORTLIBDIR)/libifcoremt.a $(FORTLIBDIR)/libimf.a $(FORTLIBDIR)/libirc.a $(FORTLIBDIR)/libsvml.a


# ------------- OSX 64 ----------------
#  OSX running test leopard (bluesky)

intel_osx_test_64 : IMAGELIBDIR = $(LIB_DIR)/lib_osx_intel_64
intel_osx_test_64 : FORTLIBDIR  = $(OSX64FORTLIBDIR)
intel_osx_test_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_OSX64 -D pp_BETA
intel_osx_test_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_OSX64 -D pp_BETA
intel_osx_test_64 : CC        = icc 
intel_osx_test_64 : CPP       = icpc
intel_osx_test_64 : FC        = ifort
intel_osx_test_64 : exe       = smv5_osx_test_64

intel_osx_test_64 : $(obj)
	icpc -m64 -o $(bin)/$(exe) $(obj)  -mmacosx-version-min=10.4 -framework OpenGL -framework GLUT -lz $(IMAGELIBDIR)/libpng.a $(IMAGELIBDIR)/libgd.a $(IMAGELIBDIR)/libglui.a $(FORTLIBDIR)/libifcoremt.a $(FORTLIBDIR)/libimf.a $(FORTLIBDIR)/libirc.a $(FORTLIBDIR)/libsvml.a

# ------------- OSX 64 ----------------
#  OSX running leopard (bluesky)

intel_osx_64 : IMAGELIBDIR = $(LIB_DIR)/lib_osx_intel_64
intel_osx_64 : FORTLIBDIR  = $(OSX64FORTLIBDIR)
intel_osx_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_OSX64
intel_osx_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_OSX64
intel_osx_64 : CC        = icc 
intel_osx_64 : CPP       = icpc
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = smv5_osx_64

intel_osx_64 : $(obj)
	icpc -m64 -o $(bin)/$(exe) $(obj)  -mmacosx-version-min=10.4 -framework OpenGL -framework GLUT -lz $(IMAGELIBDIR)/libpng.a $(IMAGELIBDIR)/libgd.a $(IMAGELIBDIR)/libglui.a $(FORTLIBDIR)/libifcoremt.a $(FORTLIBDIR)/libimf.a $(FORTLIBDIR)/libirc.a $(FORTLIBDIR)/libsvml.a

# ------------- Linux 32 ----------------
#  32 bit linux (acrux) using gnu compilers 

gcc_linux_32 : IMAGELIBDIR = $(LIB_DIR)/lib_linux_gcc_32
gcc_linux_32 : FORTLIBDIR  = $(LINUX32FORTLIBDIR)
gcc_linux_32 : FFLAGS    = -O2 -m32 -D pp_noappend -ffree-form
gcc_linux_32 : CFLAGS    = -O2 -m32 -D pp_LINUX
gcc_linux_32 : CC        = gcc
gcc_linux_32 : CPP       = g++
gcc_linux_32 : FC        = gfortran
gcc_linux_32 : exe       = smv5_linux_32

gcc_linux_32 : $(obj)
	$(CPP) -m32 -o $(bin)/$(exe) $(obj) -L$(IMAGELIBDIR) -L $(FORTLIBDIR) -lglui -lglut -lgd -ljpeg -lpng -lz -lgfortran $(SYSTEMLIBS)

# ------------- Linux 32 DEBUG ----------------
#  32 bit debug linux (acrux)

intel_linux_32_dbg : IMAGELIBDIR = $(LIB_DIR)/lib_linux_intel_32
intel_linux_32_dbg : FORTLIBDIR  = $(LINUX32FORTLIBDIR)
intel_linux_32_dbg : FFLAGS    = -O0 -m32 -traceback -g -fpe0 -fltconsistency -WB -vec_report0 -fpp -D pp_noappend
intel_linux_32_dbg : CFLAGS    = -O0 -g -m32 -vec_report0 -D pp_LINUX
intel_linux_32_dbg : CC        = icc 
intel_linux_32_dbg : CPP       = icpc
intel_linux_32_dbg : FC        = ifort
intel_linux_32_dbg : exe       = smv5_linux_32_dbg

intel_linux_32_dbg : $(obj)
	$(CPP) -m32 -o $(bin)/$(exe) $(obj) -L$(IMAGELIBDIR) -L $(FORTLIBDIR) -lglui -lglut -lgd -ljpeg -lpng -lz -lifcore $(SYSTEMLIBS)

# ------------- Linux TEST 32 ----------------
#  32 bit test linux (acrux)

intel_linux_test_32 : IMAGELIBDIR = $(LIB_DIR)/lib_linux_intel_32
intel_linux_test_32 : FORTLIBDIR  = $(LINUX32FORTLIBDIR)
intel_linux_test_32 : FFLAGS    = -O2 -m32 -vec_report0 -fpp -D pp_noappend
intel_linux_test_32 : CFLAGS    = -O2 -m32 -vec_report0 -D pp_BETA -D pp_LINUX
intel_linux_test_32 : CC        = icc 
intel_linux_test_32 : CPP       = icpc
intel_linux_test_32 : FC        = ifort
intel_linux_test_32 : exe       = smv5_linux_test_32

intel_linux_test_32 : $(obj)
	$(CPP) -m32 -o $(bin)/$(exe) $(obj) -L$(IMAGELIBDIR) -L $(FORTLIBDIR) -lglui -lglut -lgd -ljpeg -lpng -lz -lifcore -L/usr/lib $(SYSTEMLIBS)


# ------------- Linux 32 ----------------
#  32 bit linux (fire78)

intel_linux_32 : IMAGELIBDIR = $(LIB_DIR)/lib_linux_intel_32
intel_linux_32 : FORTLIBDIR  = $(LINUX32FORTLIBDIR)
intel_linux_32 : FFLAGS    = -O2 -m32 -vec_report0 -fpp -D pp_noappend
intel_linux_32 : CFLAGS    = -O2 -m32 -vec_report0 -D pp_LINUX
intel_linux_32 : CC        = icc 
intel_linux_32 : CPP       = icpc
intel_linux_32 : FC        = ifort
intel_linux_32 : exe       = smv5_linux_32

intel_linux_32 : $(obj)
	$(CPP) -m32 -o $(bin)/$(exe) $(obj) -L$(IMAGELIBDIR) -L $(FORTLIBDIR) -lglui -lglut -lgd -ljpeg -lpng -lz -lifcore -L/usr/lib $(SYSTEMLIBS)

# ------------- Linux 64 DEBUG ----------------
# 64 bit debug linux (fire79)

intel_linux_64_dbg : IMAGELIBDIR = $(LIB_DIR)/lib_linux_intel_64
intel_linux_64_dbg : FORTLIBDIR  = $(LINUX64FORTLIBDIR)
intel_linux_64_dbg : FFLAGS    = -O0 -m64 -traceback -g -fpe0 -fltconsistency -WB -vec_report0 -fpp -D pp_noappend
intel_linux_64_dbg : CFLAGS    = -O0 -g -m64 -vec_report0 -D pp_LINUX64
intel_linux_64_dbg : CC        = icc 
intel_linux_64_dbg : CPP       = icpc
intel_linux_64_dbg : FC        = ifort
intel_linux_64_dbg : exe       = smv5_linux_64_dbg

intel_linux_64_dbg : $(obj)
	$(CPP) -m64 -o $(bin)/$(exe) $(obj) -L$(IMAGELIBDIR) -L $(FORTLIBDIR) -lglui -lglut -lgd -ljpeg -lpng -lz -lifcore  -lintlc $(SYSTEMLIBS)

# ------------- Linux 64 Test ----------------
# 64 bit test linux (fire79)

intel_linux_test_64 : IMAGELIBDIR = $(LIB_DIR)/lib_linux_intel_64
intel_linux_test_64 : FORTLIBDIR  = $(LINUX64FORTLIBDIR)
intel_linux_test_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_noappend
intel_linux_test_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_LINUX64 -D pp_BETA
intel_linux_test_64 : CC        = icc 
intel_linux_test_64 : CPP       = icpc
intel_linux_test_64 : FC        = ifort
intel_linux_test_64 : exe       = smv5_linux_test_64

intel_linux_test_64 : $(obj)
	$(CPP) -m64 -o $(bin)/$(exe) $(obj) -L$(IMAGELIBDIR) -L $(FORTLIBDIR) -lglui -lglut -lgd -ljpeg -lpng -lz -lifcore  -lintlc $(SYSTEMLIBS)

# ------------- Linux 64 ----------------
# 64 bit linux (fire79)

intel_linux_64 : IMAGELIBDIR = $(LIB_DIR)/lib_linux_intel_64
intel_linux_64 : FORTLIBDIR  = $(LINUX64FORTLIBDIR)
intel_linux_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_noappend
intel_linux_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_LINUX64
intel_linux_64 : CC        = icc 
intel_linux_64 : CPP       = icpc
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = smv5_linux_64

intel_linux_64 : $(obj)
	$(CPP) -m64 -o $(bin)/$(exe) $(obj) -L$(IMAGELIBDIR) -L $(FORTLIBDIR) -lglui -lglut -lgd -ljpeg -lpng -lz -lifcore  -lintlc $(SYSTEMLIBS)

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm *.o
