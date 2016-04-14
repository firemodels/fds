LIBNAME = lpeg
LUADIR = ../../../../

DLLEXT = .dll

COPT = -O2
# COPT = -DLPEG_DEBUG -g

CWARNS = -Wall -Wextra -pedantic \
	-Waggregate-return \
	-Wcast-align \
	-Wcast-qual \
	-Wdisabled-optimization \
	-Wpointer-arith \
	-Wshadow \
	-Wsign-compare \
	-Wundef \
	-Wwrite-strings \
	-Wbad-function-cast \
	-Wdeclaration-after-statement \
	-Wmissing-prototypes \
	-Wnested-externs \
	-Wstrict-prototypes \
# -Wunreachable-code \

RANLIB=ranlib

CFLAGS = $(CWARNS) $(COPT) -std=c99 -I../lua-5.3.1/src -fPIC
CC = gcc

FILES = lpvm.o lpcap.o lptree.o lpcode.o lpprint.o

DLLFLAGS = -shared -fPIC
# DLLFLAGS =

linux:
	make lpeg.so "DLLFLAGS = -shared -fPIC"
mingw:
	make lpeg$(DLLEXT)  "DLLFLAGS = -shared -fPIC"

# For Mac OS
macosx:
	make lpeg.so "DLLFLAGS = -bundle -undefined dynamic_lookup"

lpeg.dll: $(FILES)
	$(CC) $(DLLFLAGS) $(FILES) -o $@ -L ../lua-5.3.1/src -llua53

lpeg.lib: $(FILES)
	lib $(FILES) -OUT:$@

lpeg.so: $(FILES)
	$(CC) $(DLLFLAGS) $(FILES) -o $@ -L ../lua-5.3.1/src

$(FILES): makefile

test: test.lua re.lua lpeg$(DLLEXT)
	./test.lua

clean:
	rm -f $(FILES) lpeg$(DLLEXT)


lpcap.o: lpcap.c lpcap.h lptypes.h
lpcode.o: lpcode.c lptypes.h lpcode.h lptree.h lpvm.h lpcap.h
lpprint.o: lpprint.c lptypes.h lpprint.h lptree.h lpvm.h lpcap.h
lptree.o: lptree.c lptypes.h lpcap.h lpcode.h lptree.h lpvm.h lpprint.h
lpvm.o: lpvm.c lpcap.h lptypes.h lpvm.h lpprint.h lptree.h
