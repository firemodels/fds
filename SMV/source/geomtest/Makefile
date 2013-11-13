# Makefile for Fire Dynamics Simulator (FDS)

VPATH = .
.SUFFIXES: .f90 .o .obj

.f90.o:
	$(FCOMPL) -c $(FFLAGS) $<
.f90.obj:
	$(FCOMPL) -c $(FFLAGS) $<

obj_serial = main.o
objwin_serial = $(obj_serial:.o=.obj)

intel_win_64 : FFLAGS = /O2 /Qvec_report0 /Qipo
intel_win_64 : LFLAGS = -static-intel
intel_win_64 : FCOMPL = ifort
intel_win_64 : CCOMPL = icc
intel_win_64 : obj = geomtest
intel_win_64 : $(objwin_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(objwin_serial)
