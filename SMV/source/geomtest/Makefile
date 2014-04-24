# Makefile for Fire Dynamics Simulator (FDS)

VPATH = .
.SUFFIXES: .f90 .o .obj

.f90.o:
	$(FCOMPL) -c $(FFLAGS) $<
.f90.obj:
	$(FCOMPL) -c $(FFLAGS) $<

obj_serial = temp.o geom.o main.o
objwin_serial = $(obj_serial:.o=.obj)

intel_win_64 : FFLAGS = /O2 /Qvec_report0 /Qipo
intel_win_64 : LFLAGS = -static-intel
intel_win_64 : FCOMPL = ifort
intel_win_64 : obj = geomtest_win_64
intel_win_64 : $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(LFLAGS) $(objwin_serial)

intel_linux_64 : FFLAGS = -m64 -O2 -vec_report0 
intel_linux_64 : LFLAGS = -static-intel
intel_linux_64 : FCOMPL = ifort
intel_linux_64 : obj = geomtest_linux_64
intel_linux_64 : $(obj_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(LFLAGS) $(obj_serial)

intel_linux_64p : FFLAGS = -m64 -O2 -vec_report0 -p
intel_linux_64p : LFLAGS = -static-intel
intel_linux_64p : FCOMPL = ifort
intel_linux_64p : obj = geomtestp
intel_linux_64p : $(obj_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(LFLAGS) $(obj_serial)

intel_osx_64 : FFLAGS = -m64 -O2 -vec_report0 
intel_osx_64 : LFLAGS = -static-intel
intel_osx_64 : FCOMPL = ifort
intel_osx_64 : obj = geomtest_osx_64
intel_osx_64 : $(obj_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(LFLAGS) $(obj_serial)
