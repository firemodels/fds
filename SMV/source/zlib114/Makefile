# Makefile for zlib
# Copyright (C) 1995-2002 Jean-loup Gailly.
# For conditions of distribution and use, see copyright notice in zlib.h 

# To compile and test, type:
#   ./configure; make test
# The call of configure is optional if you don't have special requirements
# If you wish to build zlib as a shared library, use: ./configure -s

# To install /usr/local/lib/libz.* and /usr/local/include/zlib.h, type:
#    make install
# To install in $HOME instead of /usr/local, use:
#    make install prefix=$HOME

COMPILER=icc
SIZE=-m64
CC=$(COMPILER)

CFLAGS=-O2 $(SIZE)

VER=1.1.4

AR=ar rc
RANLIB=ranlib

.SUFFIXES:      .c .o .obj

.c.obj:
	$(CC) -c $(CFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<

OBJS = adler32.o compress.o crc32.o gzio.o uncompr.o deflate.o trees.o \
       zutil.o inflate.o infblock.o inftrees.o infcodes.o infutil.o inffast.o
OBJSWIN = $(OBJS:.o=.obj)

libz.a: $(OBJS)
	ar rc $@ $(OBJS)
	$(RANLIB) $@

libz.lib: $(OBJSWIN)
	lib $(OBJSWIN) -OUT:$@ 
