# makefile for libpng.a and libpng12.so on Linux ELF with gcc
# Copyright (C) 1998, 1999, 2002 Greg Roelofs and Glenn Randers-Pehrson
# Copyright (C) 1996, 1997 Andreas Dilger
# For conditions of distribution and use, see copyright notice in png.h

LIBNAME = libglut
STD_DEFINES = -Dlinux -D__i386__ -D_POSIX_SOURCE -D_BSD_SOURCE -D_GNU_SOURCE -DX_LOCALE



#CC=icc
CC=gcc

CFLAGS= -O2 $(STD_DEFINES) -I ../include

RANLIB=ranlib

DESTDIR=

OBJS = glut_8x13.o glut_9x15.o glut_bitmap.o glut_bwidth.o glut_cindex.o glut_cmap.o glut_cursor.o glut_dials.o glut_dstr.o glut_event.o\
glut_ext.o

.SUFFIXES:      .c .o .pic.o

libglut.a: $(OBJS)
	ar rc $@ $(OBJS)
	$(RANLIB) $@

