#*** Smokeview, Version 5, Makefile ***

#To use this makefile, you will need to provide a make target at the
#command line.  For example, to make the 32 bit version with the Intel
#compiler, you would need to type 'make intel_linux' at the prompt.

VPATH = ../source/smokeview
#VPATH = ../../source/smokeview

obj = main.o IOscript.o glui_wui.o csphere.o colorbar.o endian.o camera.o IOhrr.o IOwui.o IOobject.o IOtour.o glui_labels.o glui_tour.o smoothlabel.o getdatacolors.o smokeview.o output.o renderfile.o isobox.o getdatabounds.o readsmv.o scontour2d.o sv_api.o getsizes.o getdata.o  glui_3dsmoke.o glui_clip.o glui_stereo.o glui_blockedit.o glui_motion.o glui_bounds.o egz_stdio.o parseobst.o dmalloc.o objectedit.o assert.o compress3dc.o IOsmoke.o IOplot3d.o IOslice.o IOboundary.o IOpart.o IOzone.o IOiso.o IOtarget.o callbacks.o drawGeometry.o glui_colorbar.o skybox.o blockage_test.o segtri.o string_util.o startup.o glui_trainer.o shaders.o

#*** General Purpose Rules ***
no_target:
	@echo \******** You did not specify a make target \********
	@echo Please read the comments at the top of the makefile

FILTERC = 
FILTERC =
INC = -I ../source/include2/glut -I ../source/include2/glui
bin = ../bin

# FFLAGS contains the compiler options to use on a compile

.SUFFIXES: .c .o .f90 .cpp

.c.o:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<
.cpp.o:
	icpc -c $(CFLAGS) $(INC) $(FILTERC) $<

.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<

# Linux 32

intel_linux_32 : IMAGELIBS = -Llib_linux_intel_32
intel_linux_32 : FORTLIBS  = -L/opt/intel/fc/10.1.021/lib
intel_linux_32 : FFLAGS    = -O2 -m32 -vec_report0 -fpp -D pp_noappend
intel_linux_32 : CFLAGS    = -O2 -m32 -vec_report0 -D pp_LINUX
intel_linux_32 : CC        = icc 
intel_linux_32 : FC        = ifort
intel_linux_32 : exe       = smv5_linux_32

intel_linux_32 : $(obj)
	icpc -m32 -o $(bin)/$(exe) $(obj) $(IMAGELIBS) $(FORTLIBS) -lglui -lglut -lgd -ljpeg -lpng -lz -lifcore -lXmu -lGLU -lGL -lm

# Linux 64

intel_linux_64 : IMAGELIBS = -Llib_linux_intel_64
intel_linux_64 : FORTLIBS  = -L/opt/intel/Compiler/11.0/074/lib/intel64
intel_linux_64 : FFLAGS    = -O2 -m64 -vec_report0 -fpp -D pp_noappend
intel_linux_64 : CFLAGS    = -O2 -m64 -vec_report0 -D pp_LINUX
intel_linux_64 : CC        = icc 
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = smv5_linux_64

intel_linux_64 : $(obj)
	icpc -m64 -o $(bin)/$(exe) $(obj) $(IMAGELIBS) $(FORTLIBS) -lglui -lglut -lgd -ljpeg -lpng -lz -lifcore -lXmu -lGLU -lGL -lm


#*** Clean Target to remove Object and Module files ***

.PHONY : clean
clean:
	rm *.o
