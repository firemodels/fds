# makefile for libpng.a and libpng12.so on Linux ELF with gcc
# Copyright (C) 1998, 1999, 2002 Greg Roelofs and Glenn Randers-Pehrson
# Copyright (C) 1996, 1997 Andreas Dilger
# For conditions of distribution and use, see copyright notice in png.h

SIZE=-m32
COMPILER=icc
CC=$(COMPILER) $(SIZE)

FILTERC=
# FILTERC= -D WIN32 -D _LIB

.c.obj:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<
.c.o:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<


CFLAGS= -O2 $(PLATFORM) -I ../GLINC

RANLIB=ranlib

# use for Linux and OSX
OBJS =  glut_8x13.o glut_9x15.o glut_bitmap.o glut_bwidth.o glut_cindex.o glut_cmap.o glut_cursor.o glut_dials.o glut_dstr.o glut_event.o \
glut_ext.o glut_fullscrn.o glut_gamemode.o glut_get.o glut_glxext.o glut_hel10.o glut_hel12.o glut_hel18.o glut_init.o  glut_input.o \
glut_joy.o glut_key.o glut_keyctrl.o glut_keyup.o glut_menu.o glut_menu2.o glut_modifier.o glut_mroman.o glut_overlay.o glut_roman.o \
glut_shapes.o glut_space.o glut_stroke.o   glut_swap.o glut_swidth.o glut_tablet.o glut_teapot.o glut_tr10.o glut_tr24.o glut_util.o \
glut_win.o glut_winmisc.o layerutil.o glut_mesa.o glut_warp.o glut_vidresize.o

# use (for now) for the PC
OBJWINS =  glut_8x13.o glut_9x15.o glut_bitmap.o glut_bwidth.o glut_cindex.o glut_cmap.o glut_cursor.o glut_dials.o glut_dstr.o glut_event.o \
glut_ext.o glut_fullscrn.o glut_gamemode.o glut_get.o glut_glxext.o glut_hel10.o glut_hel12.o glut_hel18.o glut_init.o  glut_input.o \
glut_joy.o glut_key.o glut_keyctrl.o glut_keyup.o glut_modifier.o glut_mroman.o glut_overlay.o glut_roman.o \
glut_shapes.o glut_space.o glut_stroke.o   glut_swap.o glut_swidth.o glut_tablet.o glut_teapot.o glut_tr10.o glut_tr24.o glut_util.o \
glut_win.o glut_winmisc.o glut_mesa.o glut_warp.o glut_vidresize.o win32_glx.o win32_menu.o win32_util.o win32_x11.o win32_winproc.o

OBJWINSWIN = $(OBJWINS:.o=.obj)

.SUFFIXES:      .c .o .obj .pic.o

libglut.a: $(OBJS)
	ar rc $@ $(OBJS)
	$(RANLIB) $@

libglutwin.a: $(OBJWINS)
	ar rc $@ $(OBJWINS)
	$(RANLIB) $@

libglutwin.lib: $(OBJWINSWIN)
	lib  $(OBJWINSWIN) -OUT:$@
