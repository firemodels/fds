# Generated automatically from makefile.cfg by configure.
# Makefile for Independent JPEG Group's software

# makefile.cfg is edited by configure to produce a custom Makefile.

# Read installation instructions before saying "make" !!

# For compiling with source and object files in different directories.
srcdir = .

# The name of your C compiler:
COMPILER=icc
SIZE=-m32
CC=$(COMPILER) $(SIZE)

# You may need to adjust these cc options:
CFLAGS= -O2  -I$(srcdir)

.c.obj:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<
.c.o:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<

# Link-time cc options:
LDFLAGS= 

# Library version ID; libtool uses this for the shared library version number.
# Note: we suggest this match the macro of the same name in jpeglib.h.
JPEG_LIB_VERSION = 62

# Put here the object file name for the correct system-dependent memory
# manager file.  For Unix this is usually jmemnobs.o, but you may want
# to use jmemansi.o or jmemname.o if you have limited swap space.
SYSDEPMEM= jmemnobs.o

# miscellaneous OS-dependent stuff
SHELL= /bin/sh
# linker
LN= $(CC)
# file deletion command
RM= rm -f
# library (.a) file creation command
AR= ar rc
# second step in .a creation (use "touch" if not needed)
AR2= ranlib

# End of configurable options.


# library object files common to compression and decompression
COMOBJECTS= jcomapi.o jutils.o jerror.o jmemmgr.o $(SYSDEPMEM)
# compression library object files
CLIBOBJECTS= jcapimin.o jcapistd.o jctrans.o jcparam.o \
        jdatadst.o jcinit.o jcmaster.o jcmarker.o jcmainct.o \
        jcprepct.o jccoefct.o jccolor.o jcsample.o jchuff.o \
        jcphuff.o jcdctmgr.o jfdctfst.o jfdctflt.o \
        jfdctint.o
# decompression library object files
DLIBOBJECTS= jdapimin.o jdapistd.o jdtrans.o jdatasrc.o \
        jdmaster.o jdinput.o jdmarker.o jdhuff.o jdphuff.o \
        jdmainct.o jdcoefct.o jdpostct.o jddctmgr.o \
        jidctfst.o jidctflt.o jidctint.o jidctred.o \
        jdsample.o jdcolor.o jquant1.o jquant2.o jdmerge.o
# These objectfiles are included in libjpeg.a
LIBOBJECTS= $(CLIBOBJECTS) $(DLIBOBJECTS) $(COMOBJECTS)

LIBOBJECTSWIN = $(LIBOBJECTS:.o=.obj)


.SUFFIXES: .o .obj

# the library:

libjpeg.a:  $(LIBOBJECTS)
	$(RM) libjpeg.a
	$(AR) libjpeg.a  $(LIBOBJECTS)
	$(AR2) libjpeg.a

libjpeg.lib:  $(LIBOBJECTSWIN)
	lib   $(LIBOBJECTSWIN) -OUT:libjpeg.lib
