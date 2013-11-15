# Makefile for Fire Dynamics Simulator (FDS)

VPATH = .
.SUFFIXES: .f90 .o .obj

.f90.o:
	$(FCOMPL) -c $(FFLAGS) $<
.f90.obj:
	$(FCOMPL) -c $(FFLAGS) $<

obj_serial = geom.o main.o
objwin_serial = $(obj_serial:.o=.obj)

intel_win_64 : FFLAGS = /O2 /Qvec_report0 /Qipo
intel_win_64 : LFLAGS = -static-intel
intel_win_64 : FCOMPL = ifort
intel_win_64 : obj = geomtest
intel_win_64 : $(objwin_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(objwin_serial)

intel_linux_64 : FFLAGS = -m64 -O2 -vec_report0 
intel_linux_64 : LFLAGS = -static-intel
intel_linux_64 : FCOMPL = ifort
intel_linux_64 : obj = geomtest
intel_linux_64 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

intel_osx_64 : FFLAGS = -m64 -O2 -vec_report0 
intel_osx_64 : LFLAGS = -static-intel
intel_osx_64 : FCOMPL = ifort
intel_osx_64 : obj = geomtest
intel_osx_64 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)
