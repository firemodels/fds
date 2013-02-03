# makefile for libpng.a and libpng12.so on Linux ELF with gcc
# Copyright (C) 1998, 1999, 2002 Greg Roelofs and Glenn Randers-Pehrson
# Copyright (C) 1996, 1997 Andreas Dilger
# For conditions of distribution and use, see copyright notice in png.h

LIBNAME = libpng12
PNGMAJ = 0
PNGMIN = 1.2.5
PNGVER = $(PNGMAJ).$(PNGMIN)

COMPILER=icc
SIZE=-m32
CC=$(COMPILER) $(SIZE)

# where "make install" puts libpng12.a, libpng12.so*,
# libpng12/png.h and libpng12/pngconf.h
# Prefix must be a full pathname.
prefix=/usr/local

# Where the zlib library and include files are located.
ZLIBLIB=../zlib114
ZLIBINC=../zlib114

ALIGN=
# for i386:
#ALIGN=-malign-loops=2 -malign-functions=2

WARNMORE=-Wwrite-strings -Wpointer-arith -Wshadow \
	-Wmissing-declarations -Wtraditional -Wcast-align \
	-Wstrict-prototypes -Wmissing-prototypes #-Wconversion

# for pgcc version 2.95.1, -O3 is buggy; don't use it.

CFLAGS=-I$(ZLIBINC) -Wall -O3 -funroll-loops \
	$(ALIGN) # $(WARNMORE) -g -DPNG_DEBUG=5

LDFLAGS=-L. -Wl,-rpath,. -L$(ZLIBLIB) -Wl,-rpath,$(ZLIBLIB) -lpng12 -lz -lm
LDFLAGS_A=-L$(ZLIBLIB) -Wl,-rpath,$(ZLIBLIB) libpng.a -lz -lm

.c.obj:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<
.c.o:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<

RANLIB=ranlib
RM=/bin/rm -f
OBJS = png.o pngset.o pngget.o pngrutil.o pngtrans.o pngwutil.o \
	pngread.o pngrio.o pngwio.o pngwrite.o pngrtran.o \
	pngwtran.o pngmem.o pngerror.o pngpread.o

.SUFFIXES:      .c .o .pic.o

.c.pic.o:
	$(CC) -c $(CFLAGS) -fPIC -o $@ $*.c

all: libpng.a $(LIBNAME).so pngtest pngtest-static libpng.pc libpng-config

libpng.a: $(OBJS)
	ar rc $@ $(OBJS)
	$(RANLIB) $@

clean:
	$(RM) *.o libpng.a
