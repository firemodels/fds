# Makefile for Fire Dynamics Simulator (FDS)

VPATH = .
.SUFFIXES: .f90 .o

.f90.o:
	$(FCOMPL) -c $(FFLAGS) $<
obj_serial = fisosmoke.o main.o

intel_linux_64 : FFLAGS = -m64 -O3 -vec_report0
intel_linux_64 : CFLAGS = -m64 -O -Dpp_noappend -vec_report0
intel_linux_64 : LFLAGS = -static-intel
intel_linux_64 : FCOMPL = ifort
intel_linux_64 : CCOMPL = icc
intel_linux_64 : obj = geomtest
intel_linux_64 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)
