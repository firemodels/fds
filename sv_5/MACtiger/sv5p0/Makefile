# FILTER contains the keys to use when filtering source files.  

FILTERF = -D pp_OSX -vec_report0
FILTERC = -D pp_OSX -vec_report0 -D pp_SV5P0

# FFLAGS contains the compiler options to use on a compile

CC        = icc
CPP        = icc
CFLAGS = 
FC        = ifort
FFLAGS =  -fpp
INCLUDE = -I ../glut/include -I ../glui -I ../gd
bin = ../../bin
exe = sv5p0_osx_intel

VPATH = ../../source/smokeview
.SUFFIXES:
.SUFFIXES: .c .o .f90 .cpp


.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $(FILTERC) $<
.cpp.o:
	$(CPP) -c $(CFLAGS) $(INCLUDE) $(FILTERC) $<

.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<

obj = main.o colorbar.o endian.o camera.o IOhrr.o IOwui.o IOobject.o IOtour.o glui_labels.o glui_tour.o smoothlabel.o getdatacolors.o smokeview.o output.o renderfile.o isobox.o getdatabounds.o readsmv.o scontour2d.o sv_api.o getsizes.o getdata.o  glui_3dsmoke.o glui_clip.o glui_stereo.o glui_blockedit.o glui_motion.o glui_bounds.o egz_stdio.o parseobst.o dmalloc.o objectedit.o assert.o compress3dc.o IOsmoke.o IOplot3d.o IOslice.o IOboundary.o IOpart.o IOzone.o IOiso.o IOtarget.o callbacks.o drawGeometry.o glui_colorbar.o skybox.o blockage_test.o segtri.o string_util.o startup.o glui_trainer.o

$(bin)/$(exe) : $(obj)
	g++ -o $(bin)/$(exe) $(obj) -framework OpenGL -framework GLUT ../lib//libz.a ../lib/libjpeg.a ../lib/libpng.a ../lib/libgd.a ../lib/libglui.a ../fortlib/libifcoremt.a ../fortlib/libimf.a ../fortlib/libirc.a

clean	:
	rm *.o
