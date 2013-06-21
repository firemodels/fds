# Makefile for volcheck

# The variable VPATH defines where the source code is relative to the current directory

VPATH = .;../smokeview

obj_serial = getdata.o volcheck.o

objwin_serial = $(obj_serial:.o=.obj)

# General Purpose Rules

no_target:
	@echo \******** You did not specify a make target \********
	@echo Please read the comments at the top of the makefile

setup:
%.o : %.mod

setup_win:
%.obj : %.mod

.SUFFIXES: .f90 .o .obj

.f90.o:
	$(FCOMPL) -c $(FFLAGS) $<
.f90.obj:
	$(FCOMPL) -c $(FFLAGS) $<


intel_win_64 : FFLAGS = /O3 /Qvec_report0 /Qipo -fpp
intel_win_64 : FCOMPL = ifort
intel_win_64 : obj = fds_win_64
intel_win_64 : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F1000000000 $(objwin_serial)

intel_linux_64_profile : FFLAGS = -m64 -O2 -vec_report0 -p
intel_linux_64_profile : LFLAGS = -static-intel
intel_linux_64_profile : FCOMPL = ifort
intel_linux_64_profile : obj = fds_intel_linux_64_profile
intel_linux_64_profile : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

intel_linux_64 : FFLAGS = -m64 -O3 -vec_report0 -ipo
intel_linux_64 : LFLAGS = -static-intel
intel_linux_64 : FCOMPL = ifort
intel_linux_64 : obj = fds_intel_linux_64
intel_linux_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)


intel_osx_64 : FFLAGS = -O3 -m64 -heap-arrays -ipo -mmacosx-version-min=10.4
intel_osx_64 : LFLAGS = -static-intel
intel_osx_64 : FCOMPL = ifort
intel_osx_64 : obj = fds_intel_osx_64
intel_osx_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

#*** Clean Target to remove Object and Module files ***

.PHONY : clean
clean:
	-rm -f *.o *.mod *.obj 

.PHONY : winclean
winclean:
	-erase *.mod *.obj 
