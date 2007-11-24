# FILTER contains the keys to use when filtering source files.  

FILTERC = 

# FFLAGS contains the compiler options to use on a compile

CC        = gcc
CFLAGS = 
FC        = lf95
FFLAGS =
INCLUDE = 
bin = ../../bin
exe = blockaid_linux

VPATH = ../../source/blockaid
.SUFFIXES:
.SUFFIXES: .c .o 


.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $(FILTERC) $<
.cpp.o:
	g++ -c $(CFLAGS) $(INCLUDE) $(FILTERC) $<

.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<

obj = main.o assert.o dmalloc.o readfds.o utilities.o

$(bin)/$(exe) : $(obj)
	gcc -o $(bin)/$(exe) -lm $(obj)
\c
