# FILTER contains the keys to use when filtering source files.  

FILTERF = -Cpp -D pp_noappend -D pp_LAHEY -D pp_PART5
FILTERC = -D pp_SV5P0
#FILTERC = -D pp_aloft -D pp_DRAWISO -D EGZ -D USE_ZLIB -D pp_glui -fno-gnu-linker -D pp_NISTREVIEW

# FFLAGS contains the compiler options to use on a compile

CC        = gcc
CFLAGS = 
FC        = lf95
FFLAGS = -O --tpp
INCLUDE = -I ../glut/include -I ../glui
bin = ../../bin
exe = sv5p0_linux

VPATH = ../../source/smokeview
.SUFFIXES:
.SUFFIXES: .c .o .f90 .cpp


.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $(FILTERC) $<
.cpp.o:
	g++ -c $(CFLAGS) $(INCLUDE) $(FILTERC) $<

.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<

obj = main.o colorbar.o endian.o camera.o IOhrr.o IOwui.o IOobject.o IOtour.o glui_labels.o glui_tour.o smoothlabel.o getdatacolors.o smokeview.o output.o renderfile.o isobox.o getdatabounds.o readsmv.o scontour2d.o sv_api.o getsizes.o getdata.o  glui_3dsmoke.o glui_clip.o glui_stereo.o glui_blockedit.o glui_motion.o glui_bounds.o egz_stdio.o parseobst.o dmalloc.o objectedit.o assert.o compress3dc.o IOsmoke.o IOplot3d.o IOslice.o IOboundary.o IOpart.o IOzone.o IOiso.o IOtarget.o callbacks.o drawGeometry.o glui_colorbar.o skybox.o blockage_test.o segtri.o string_util.o startup.o glui_trainer.o

$(bin)/$(exe) : $(obj)
	gcc -o $(bin)/$(exe) $(obj) ../lib/libglui.a ../lib/libglut.a /usr/X11R6/lib/libXmu.a -lGLU -lGL -lm  ../lib/libgd.a ../lib/libjpeg.a ../lib/libpng.a ../lib/libz.a ./libfst.a ./libfj9fpp.a ./libfj9i6.a ./libfj9ipp.a ./libfj9e6.a ./libfj9f6.a 

clean	:
	rm *.o
