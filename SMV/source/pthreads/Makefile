COMPILER=icc
SIZE=-m32
CC=$(COMPILER) $(SIZE)

INC=-I.
CFLAGS=-o2
FILTERC=-D PTW32_BUILD -D PTW32_STATIC_LIB

.SUFFIXES:      .c .o .obj

.c.obj:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<
.c.o:
	$(CC) -c $(CFLAGS) $(INC) $(FILTERC) $<

RANLIB=ranlib
RM=/bin/rm -f
OBJS = pthread.o
OBJSWIN = $(OBJS:.o=.obj)

libpthreads.a: $(OBJS)
	ar rc $@ $(OBJS)
	$(RANLIB) $@

libpthreads.lib: $(OBJSWIN)
	lib $(OBJSWIN) -OUT:$@
