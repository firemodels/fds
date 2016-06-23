#If the ar command fails on your system, consult the ar manpage
#for your system. 
AR=ar

CFLAGS=-g  -DHAVE_LIBPNG -DHAVE_LIBZ -DHAVE_LIBJPEG -DHAVE_LIBX11

INCLUDEDIRS=-I. -I../png-1.6.21 -I../jpeg-9b -I../zlib128

#
#
# Changes should not be required below here.
#
#

# Update these with each release!

MAJOR_VERSION=2
VERSION=2.0.4

COMPILER=gcc -fPIC
#COMPILER=icc
#SIZE=-m32
SIZE=-m64

CC=$(COMPILER) $(SIZE) $(INCLUDEDIRS)

.SUFFIXES:      .c .o .obj

.c.obj:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<
.c.o:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<

LIBOBJS=gd.o gd_gd.o gd_gd2.o gd_io.o gd_io_dp.o \
		gd_io_file.o gd_ss.o gd_io_ss.o gd_png.o gd_jpeg.o gdxpm.o \
		gdfontt.o gdfonts.o gdfontmb.o gdfontl.o gdfontg.o \
		gdtables.o gdft.o gdcache.o gdkanji.o wbmp.o \
		gd_wbmp.o gdhelpers.o gd_topal.o 
LIBOBJSWIN = $(LIBOBJS:.o=.obj)

#Static library, if you really need one for some reason.
libgd.a: ${LIBOBJS}
	$(RM) libgd.a
	$(AR) rc libgd.a ${LIBOBJS}
	-ranlib libgd.a

libgd.lib: ${LIBOBJSWIN}
	lib ${LIBOBJSWIN} -OUT:libgd.lib 
