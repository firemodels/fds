# Makefile for Fire Dynamics Simulator (FDS)

# To use this makefile, select the appropriate compiler and OS from the list below. For example, to
# compile the 64 bit parallel version (MPI) of FDS under Linux using the Intel compilers type:
#
# make mpi_intel_linux_64

# A make utility for windows may be downloaded from http://gnuwin32.sourceforge.net/packages/make.htm

# The variable VPATH defines where the source code is relative to the current directory

VPATH = ../FDS_Source

# Definition of the non-MPI object variables

obj_serial = prec.o mpis.o cons.o devc.o data.o type.o mesh.o func.o smvv.o irad.o \
             ieva.o pois.o scrc.o radi.o evac.o geom.o part.o vege.o ctrl.o turb.o dump.o hvac.o mass.o \
             read.o wall.o fire.o divg.o init.o velo.o pres.o main.o

objwin_serial = $(obj_serial:.o=.obj)

# Definition of the MPI object variables

obj_mpi = prec.o mpip.o cons.o devc.o data.o type.o mesh.o func.o smvv.o irad.o \
          ieva.o pois.o scrc.o radi.o evac.o geom.o part.o vege.o ctrl.o turb.o dump.o hvac.o mass.o \
          read.o wall.o fire.o divg.o init.o velo.o pres.o main.o

objwin_mpi = $(obj_mpi:.o=.obj)

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

# Files which have OpenMP compiler directives are compiled using the
# FOPENMPFLAG and can be enabled/disabled on a per-file basis below.
# If the target name includes "openmp", then the OpenMP flags are appended.
ifeq ($(findstring openmp, $(MAKECMDGOALS)), openmp)
# divg.o: FFLAGS += $(FOPENMPFLAGS)
# fire.o: FFLAGS += $(FOPENMPFLAGS)
func.o: FFLAGS += $(FOPENMPFLAGS)
# mass.o: FFLAGS += $(FOPENMPFLAGS)
# part.o: FFLAGS += $(FOPENMPFLAGS)
# pois.o: FFLAGS += $(FOPENMPFLAGS)
pres.o: FFLAGS += $(FOPENMPFLAGS)
# turb.o: FFLAGS += $(FOPENMPFLAGS)
# velo.o: FFLAGS += $(FOPENMPFLAGS)

# divg.obj: FFLAGS += $(FOPENMPFLAGS)
# fire.obj: FFLAGS += $(FOPENMPFLAGS)
func.obj: FFLAGS += $(FOPENMPFLAGS)
# mass.obj: FFLAGS += $(FOPENMPFLAGS)
# part.obj: FFLAGS += $(FOPENMPFLAGS)
# pois.obj: FFLAGS += $(FOPENMPFLAGS)
pres.obj: FFLAGS += $(FOPENMPFLAGS)
# turb.obj: FFLAGS += $(FOPENMPFLAGS)
# velo.obj: FFLAGS += $(FOPENMPFLAGS)
endif

# Windows Intel Compiler, version 11

#*** note: the following script (located in ..\Scripts directory) must be run before running make with this target
#    ..\Scripts\ifortvars ia32

intel_win_32_db : FFLAGS = /Qvec_report0 /nologo /debug:full /Od /extend_source:132 /warn:all /warn:nointerfaces /Qtrapuv /fp:strict /fp:except /traceback /check:all /stand:f03
intel_win_32_db : FCOMPL = ifort
intel_win_32_db : obj = fds_win_32_db
intel_win_32_db : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F100000000 $(objwin_serial)

openmp_intel_win_32_db : FFLAGS = /Qvec_report0 /nologo /debug:full /Od /extend_source:132 /warn:all /warn:nointerfaces /Qtrapuv /fp:strict /fp:except /traceback /check:all /stand:f03
openmp_intel_win_32_db : FCOMPL = ifort
openmp_intel_win_32_db : FOPENMPFLAGS = /Qopenmp
openmp_intel_win_32_db : obj = fds_openmp_win_32_db
openmp_intel_win_32_db : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(FOPENMPFLAGS) /F100000000 $(objwin_serial)

intel_win_32 : FFLAGS = /O2 /Qvec_report0 /Qipo
intel_win_32 : FCOMPL = ifort
intel_win_32 : obj = fds_win_32
intel_win_32 : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F100000000 $(objwin_serial)

intel_win_32_ps : FFLAGS = /O3 /Qipo  /Qvec_report0 
intel_win_32_ps : FCOMPL = ifort
intel_win_32_ps : obj = fds_win_32_ps
intel_win_32_ps : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F100000000 $(objwin_serial) /link /SUBSYSTEM:CONSOLE,5.01

openmp_intel_win_32 : FFLAGS = /O2 /Qvec_report0 /Qipo
openmp_intel_win_32 : FCOMPL = ifort
openmp_intel_win_32 : FOPENMPFLAGS = /Qopenmp
openmp_intel_win_32 : obj = fds_openmp_win_32
openmp_intel_win_32 : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(FOPENMPFLAGS) /F100000000 $(objwin_serial)

mpi_intel_win_32 : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
mpi_intel_win_32 : MPIINCLUDE = "c:\program files\mpich2\include"
mpi_intel_win_32 : FFLAGS = /O2 /Qipo /I $(MPIINCLUDE)
mpi_intel_win_32 : FCOMPL = ifort
mpi_intel_win_32 : obj = fds_mpi_win_32
mpi_intel_win_32 : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F100000000 $(objwin_mpi) $(MPILIB)

mpi_intel_win_32_ps : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
mpi_intel_win_32_ps : MPIINCLUDE = "c:\program files\mpich2\include"
mpi_intel_win_32_ps : FFLAGS = /O2 /Qipo /Qvec_report0  /QxHost  /I $(MPIINCLUDE)
mpi_intel_win_32_ps : FCOMPL = ifort
mpi_intel_win_32_ps : obj = fds_mpi_win_32_ps
mpi_intel_win_32_ps : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F100000000 $(objwin_mpi) $(MPILIB) /link /SUBSYSTEM:CONSOLE,5.01

openmp_mpi_intel_win_32 : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
openmp_mpi_intel_win_32 : MPIINCLUDE = "c:\program files\mpich2\include"
openmp_mpi_intel_win_32 : FFLAGS = /O2 /Qipo /I $(MPIINCLUDE)
openmp_mpi_intel_win_32 : FCOMPL = ifort
openmp_mpi_intel_win_32 : FOPENMPFLAGS = /Qopenmp
openmp_mpi_intel_win_32 : obj = fds_openmp_mpi_win_32
openmp_mpi_intel_win_32 : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(FOPENMPFLAGS) /F100000000 $(objwin_mpi) $(MPILIB)

mpi_intel_win_32_db : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
mpi_intel_win_32_db : MPIINCLUDE = "c:\program files\mpich2\include"
mpi_intel_win_32_db : FFLAGS = /O2 /I $(MPIINCLUDE) /Qvec_report0 /nologo /debug:full /Od /extend_source:132 /warn:all /warn:nointerfaces /Qtrapuv /fp:strict /fp:except /traceback /check:all /stand:f03
mpi_intel_win_32_db : FCOMPL = ifort
mpi_intel_win_32_db : obj = fds_mpi_win_32_db
mpi_intel_win_32_db : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F100000000 $(objwin_mpi) $(MPILIB)

openmp_mpi_intel_win_32_db : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
openmp_mpi_intel_win_32_db : MPIINCLUDE = "c:\program files\mpich2\include"
openmp_mpi_intel_win_32_db : FFLAGS = /O2 /I $(MPIINCLUDE) /Qvec_report0 /nologo /debug:full /Od /extend_source:132 /warn:all /warn:nointerfaces /Qtrapuv /fp:strict /fp:except /traceback /check:all /stand:f03
openmp_mpi_intel_win_32_db : FCOMPL = ifort
openmp_mpi_intel_win_32_db : FOPENMPFLAGS = /Qopenmp
openmp_mpi_intel_win_32_db : obj = fds_openmp_mpi_win_32_db
openmp_mpi_intel_win_32_db : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(FOPENMPFLAGS) /F100000000 $(objwin_mpi) $(MPILIB)

 #*** note: the following scripts (located in ..\Scripts directory) must be run before running make with this target
 #    ..\Scripts\ifortvars intel64

intel_win_64 : FFLAGS = /O2 /Qvec_report0 /Qipo
intel_win_64 : FCOMPL = ifort
intel_win_64 : obj = fds_win_64
intel_win_64 : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F1000000000 $(objwin_serial)

intel_win_64_ps : FFLAGS = -O2 /Qipo /Qvec_report0  
intel_win_64_ps : FCOMPL = ifort
intel_win_64_ps : obj = fds_win_64_ps
intel_win_64_ps : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F1000000000 $(objwin_serial) /link /SUBSYSTEM:CONSOLE,5.02

openmp_intel_win_64 : FFLAGS = /O2 /Qvec_report0 /Qipo
openmp_intel_win_64 : FCOMPL = ifort
openmp_intel_win_64 : FOPENMPFLAGS = /Qopenmp
openmp_intel_win_64 : obj = fds_openmp_win_64
openmp_intel_win_64 : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(FOPENMPFLAGS) /F1000000000 $(objwin_serial)

intel_win_64_dv : FFLAGS = /O2 /Qvec_report0
intel_win_64_dv : FCOMPL = ifort
intel_win_64_dv : obj = fds_win_64_dv
intel_win_64_dv : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F1000000000 $(objwin_serial)

openmp_intel_win_64_dv : FFLAGS = /O2 /Qvec_report0
openmp_intel_win_64_dv : FCOMPL = ifort
openmp_intel_win_64_dv : FOPENMPFLAGS = /Qopenmp
openmp_intel_win_64_dv : obj = fds_openmp_win_64_dv
openmp_intel_win_64_dv : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(FOPENMPFLAGS) /F1000000000 $(objwin_serial)

intel_win_64_db : FFLAGS = /Qvec_report0 /nologo /debug:full /Od /extend_source:132 /warn:all /warn:nointerfaces /Qtrapuv /fp:strict /fp:except /traceback /check:all /stand:f03
intel_win_64_db : FCOMPL = ifort
intel_win_64_db : obj = fds_win_64_db
intel_win_64_db : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F1000000000 $(objwin_serial)

openmp_intel_win_64_db : FFLAGS = $(FOPTIMIZATIONFLAGS) /Qvec_report0 /nologo /debug:full /extend_source:132 /traceback /stand:f03 /Zi /MDd
openmp_intel_win_64_db : FCOMPL = ifort
openmp_intel_win_64_db : FOPENMPFLAGS = /Qopenmp
openmp_intel_win_64_db : FOPTIMIZATIONFLAGS = /Od
openmp_intel_win_64_db : obj = fds_openmp_win_64_db
openmp_intel_win_64_db : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(FOPENMPFLAGS) /F1000000000 $(objwin_serial)

mpi_intel_win_64 : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
mpi_intel_win_64 : MPIINCLUDE = "c:\program files\mpich2\include"
mpi_intel_win_64 : FFLAGS = /O2 /Qipo /Qvec_report0 /I $(MPIINCLUDE)
mpi_intel_win_64 : FCOMPL = ifort
mpi_intel_win_64 : obj = fds_mpi_win_64
mpi_intel_win_64 : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F1000000000 $(objwin_mpi) $(MPILIB)

mpi_intel_win_64_ps : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
mpi_intel_win_64_ps : MPIINCLUDE = "c:\program files\mpich2\include"
mpi_intel_win_64_ps : FFLAGS = -O2 /Qipo /Qvec_report0  /QxHost  /I $(MPIINCLUDE)
mpi_intel_win_64_ps : FCOMPL = ifort
mpi_intel_win_64_ps : obj = fds_mpi_win_64_ps
mpi_intel_win_64_ps : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F1000000000 $(objwin_mpi) $(MPILIB) /link /SUBSYSTEM:CONSOLE,5.02

openmp_mpi_intel_win_64 : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
openmp_mpi_intel_win_64 : MPIINCLUDE = "c:\program files\mpich2\include"
openmp_mpi_intel_win_64 : FFLAGS = /O2 /Qipo /Qvec_report0 /I $(MPIINCLUDE)
openmp_mpi_intel_win_64 : FCOMPL = ifort
openmp_mpi_intel_win_64 : FOPENMPFLAGS = /Qopenmp
openmp_mpi_intel_win_64 : obj = fds_openmp_mpi_win_64
openmp_mpi_intel_win_64 : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(FOPENMPFLAGS) /F1000000000 $(objwin_mpi) $(MPILIB)

mpi_intel_win_64_db : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
mpi_intel_win_64_db : MPIINCLUDE = "c:\program files\mpich2\include"
mpi_intel_win_64_db : FFLAGS = /O2 /I $(MPIINCLUDE) /Qvec_report0 /nologo /debug:full /Od /extend_source:132 /warn:all /warn:nointerfaces /Qtrapuv /fp:strict /fp:except /traceback /check:all /stand:f03
mpi_intel_win_64_db : FCOMPL = ifort
mpi_intel_win_64_db : obj = fds_mpi_win_64_db
mpi_intel_win_64_db : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F1000000000 $(objwin_mpi) $(MPILIB)

openmp_mpi_intel_win_64_db : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
openmp_mpi_intel_win_64_db : MPIINCLUDE = "c:\program files\mpich2\include"
openmp_mpi_intel_win_64_db : FFLAGS = /O2 /I $(MPIINCLUDE) /Qvec_report0 /nologo /debug:full /Od /extend_source:132 /warn:all /warn:nointerfaces /Qtrapuv /fp:strict /fp:except /traceback /check:all /stand:f03
openmp_mpi_intel_win_64_db : FCOMPL = ifort
openmp_mpi_intel_win_64_db : FOPENMPFLAGS = /Qopenmp
openmp_mpi_intel_win_64_db : obj = fds_openmp_mpi_win_64_db
openmp_mpi_intel_win_64_db : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(FOPENMPFLAGS) /F1000000000 $(objwin_mpi) $(MPILIB)

# Intel Compiler for Linux, version 11

intel_linux_32 : FFLAGS = -m32 -O2 -vec_report0 -ipo
intel_linux_32 : LFLAGS = -static-intel
intel_linux_32 : FCOMPL = ifort
intel_linux_32 : obj = fds_intel_linux_32
intel_linux_32 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

openmp_intel_linux_32 : FFLAGS = -m32 -O2 -vec_report0 -ipo
openmp_intel_linux_32 : LFLAGS = -static-intel
openmp_intel_linux_32 : FCOMPL = ifort
openmp_intel_linux_32 : FOPENMPFLAGS = -openmp -openmp-link static
openmp_intel_linux_32 : obj = fds_openmp_intel_linux_32
openmp_intel_linux_32 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_serial)

intel_linux_32_dv : FFLAGS = -m32 -O2 -vec_report0 -traceback
intel_linux_32_dv : LFLAGS = -static-intel
intel_linux_32_dv : FCOMPL = ifort
intel_linux_32_dv : obj = fds_intel_linux_32_dv
intel_linux_32_dv : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

openmp_intel_linux_32_dv : FFLAGS = -m32 -O2 -vec_report0 -traceback
openmp_intel_linux_32_dv : LFLAGS = -static-intel
openmp_intel_linux_32_dv : FCOMPL = ifort
openmp_intel_linux_32_dv : FOPENMPFLAGS = -openmp -openmp-link static
openmp_intel_linux_32_dv : obj = fds_openmp_intel_linux_32_dv
openmp_intel_linux_32_dv : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_serial)

#
# to profile FDS:
# 1.  build FDS using intel_linux_32_profile target
# 2.  run case as usual
# 3.  type: gprof fds_linux_32_profile > results.out

intel_linux_32_profile : FFLAGS = -m32 -O2 -vec_report0 -p
intel_linux_32_profile : LFLAGS = -static-intel
intel_linux_32_profile : FCOMPL = ifort
intel_linux_32_profile : obj = fds_intel_linux_32_profile
intel_linux_32_profile : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

openmp_intel_linux_32_profile : FFLAGS = -m32 -O2 -vec_report0 -p
openmp_intel_linux_32_profile : LFLAGS = -static-intel
openmp_intel_linux_32_profile : FCOMPL = ifort
openmp_intel_linux_32_profile : FOPENMPFLAGS = -openmp -openmp-link static
openmp_intel_linux_32_profile : obj = fds_openmp_intel_linux_32_profile
openmp_intel_linux_32_profile : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_serial)

intel_linux_64_profile : FFLAGS = -m64 -O2 -vec_report0 -p
intel_linux_64_profile : LFLAGS = -static-intel
intel_linux_64_profile : FCOMPL = ifort
intel_linux_64_profile : obj = fds_intel_linux_64_profile
intel_linux_64_profile : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

openmp_intel_linux_64_profile : FFLAGS = -m64 -O2 -vec_report0 -p
openmp_intel_linux_64_profile : LFLAGS = -static-intel
openmp_intel_linux_64_profile : FCOMPL = ifort
openmp_intel_linux_64_profile : FOPENMPFLAGS = -openmp -openmp-link static
openmp_intel_linux_64_profile : obj = fds_openmp_intel_linux_64_profile
openmp_intel_linux_64_profile : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_serial)

intel_linux_64 : FFLAGS = -m64 -O2 -vec_report0 -ipo
intel_linux_64 : LFLAGS = -static-intel
intel_linux_64 : FCOMPL = ifort
intel_linux_64 : obj = fds_intel_linux_64
intel_linux_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

openmp_intel_linux_64 : FFLAGS = -m64 -O2 -vec_report0 -ipo
openmp_intel_linux_64 : LFLAGS = -static-intel
openmp_intel_linux_64 : FCOMPL = ifort
openmp_intel_linux_64 : FOPENMPFLAGS = -openmp -openmp-link static
openmp_intel_linux_64 : obj = fds_openmp_intel_linux_64
openmp_intel_linux_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_serial)

intel_linux_64_dv : FFLAGS = -m64 -O2 -vec_report0 -traceback
intel_linux_64_dv : LFLAGS = -static-intel
intel_linux_64_dv : FCOMPL = ifort
intel_linux_64_dv : obj = fds_intel_linux_64_dv
intel_linux_64_dv : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

openmp_intel_linux_64_dv : FFLAGS = -m64 -O2 -vec_report0 -traceback
openmp_intel_linux_64_dv : LFLAGS = -static-intel
openmp_intel_linux_64_dv : FCOMPL = ifort
openmp_intel_linux_64_dv : FOPENMPFLAGS = -openmp -openmp-link static
openmp_intel_linux_64_dv : obj = fds_openmp_intel_linux_64_dv
openmp_intel_linux_64_dv : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_serial)

mpi_intel_linux_32 : FFLAGS = -m32 -O2 -vec_report0 -ipo
mpi_intel_linux_32 : LFLAGS = -static-intel
mpi_intel_linux_32 : FCOMPL = mpif77
mpi_intel_linux_32 : obj = fds_mpi_intel_linux_32
mpi_intel_linux_32 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_mpi)

openmp_mpi_intel_linux_32 : FFLAGS = -m32 -O2 -vec_report0 -ipo
openmp_mpi_intel_linux_32 : LFLAGS = -static-intel
openmp_mpi_intel_linux_32 : FCOMPL = mpif77
openmp_mpi_intel_linux_32 : FOPENMPFLAGS = -openmp -openmp-link static
openmp_mpi_intel_linux_32 : obj = fds_openmp_mpi_intel_linux_32
openmp_mpi_intel_linux_32 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi)

mpi_intel_linux_64 : FFLAGS = -m64 -O2 -vec_report0 -ipo
mpi_intel_linux_64 : LFLAGS = -static-intel
mpi_intel_linux_64 : FCOMPL = mpif77
mpi_intel_linux_64 : obj = fds_mpi_intel_linux_64
mpi_intel_linux_64 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_mpi)

mpi_intel_linux_64ib : FFLAGS = -m64 -O2 -vec_report0 -ipo
mpi_intel_linux_64ib : LFLAGS = -static-intel
mpi_intel_linux_64ib : FCOMPL = mpif77
mpi_intel_linux_64ib : obj = fds_mpi_intel_linux_64ib
mpi_intel_linux_64ib : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_mpi)

openmp_mpi_intel_linux_64 : FFLAGS = -m64 -O2 -vec_report0 -ipo
openmp_mpi_intel_linux_64 : LFLAGS = -static-intel
openmp_mpi_intel_linux_64 : FCOMPL = mpif77
openmp_mpi_intel_linux_64 : FOPENMPFLAGS = -openmp -openmp-link static
openmp_mpi_intel_linux_64 : obj = fds_openmp_mpi_intel_linux_64
openmp_mpi_intel_linux_64 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi)

intel_linux_32_db : FFLAGS = -m32 -check all -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -vec_report0 -ftrapuv -debug full  -stand:f03
intel_linux_32_db : LFLAGS = -static-intel
intel_linux_32_db : FCOMPL = ifort
intel_linux_32_db : obj = fds_intel_linux_32_db
intel_linux_32_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

openmp_intel_linux_32_db : FFLAGS = -m32 -check all -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -vec_report0 -ftrapuv -debug full  -stand:f03
openmp_intel_linux_32_db : LFLAGS = -static-intel
openmp_intel_linux_32_db : FCOMPL = ifort
openmp_intel_linux_32_db : FOPENMPFLAGS = -openmp -openmp-link static
openmp_intel_linux_32_db : obj = fds_openmp_intel_linux_32_db
openmp_intel_linux_32_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS $(FOPENMPFLAGS) -o $(obj) $(obj_serial)

intel_linux_64_db : FFLAGS = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -vec_report0 -stand:f03
intel_linux_64_db : LFLAGS = -static-intel
intel_linux_64_db : FCOMPL = ifort
intel_linux_64_db : obj = fds_intel_linux_64_db
intel_linux_64_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

openmp_intel_linux_64_db : FFLAGS = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -vec_report0 -stand:f03
openmp_intel_linux_64_db : LFLAGS = -static-intel
openmp_intel_linux_64_db : FCOMPL = ifort
openmp_intel_linux_64_db : FOPENMPFLAGS = -openmp -openmp-link static
openmp_intel_linux_64_db : obj = fds_openmp_intel_linux_64_db
openmp_intel_linux_64_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_serial)

mpi_intel_linux_32_db : FFLAGS = -m32 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency  -vec_report0
mpi_intel_linux_32_db : LFLAGS = -static-intel
mpi_intel_linux_32_db : FCOMPL = mpif77
mpi_intel_linux_32_db : obj = fds_mpi_intel_linux_32_db
mpi_intel_linux_32_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_mpi)

openmp_mpi_intel_linux_32_db : FFLAGS = -m32 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency  -vec_report0
openmp_mpi_intel_linux_32_db : LFLAGS = -static-intel
openmp_mpi_intel_linux_32_db : FCOMPL = mpif77
openmp_mpi_intel_linux_32_db : FOPENMPFLAGS = -openmp -openmp-link static
openmp_mpi_intel_linux_32_db : obj = fds_openmp_mpi_intel_linux_32_db
openmp_mpi_intel_linux_32_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi)

mpi_intel_linux_64_db : FFLAGS = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -vec_report0
mpi_intel_linux_64_db : LFLAGS = -static-intel
mpi_intel_linux_64_db : FCOMPL = mpif77
mpi_intel_linux_64_db : obj = fds_mpi_intel_linux_64_db
mpi_intel_linux_64_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_mpi)

mpi_intel_linux_64ib_db : FFLAGS = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -vec_report0
mpi_intel_linux_64ib_db : LFLAGS = -static-intel
mpi_intel_linux_64ib_db : FCOMPL = mpif77
mpi_intel_linux_64ib_db : obj = fds_mpi_intel_linux_64ib_db
mpi_intel_linux_64ib_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_mpi)

openmp_mpi_intel_linux_64_db : FFLAGS = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -vec_report0
openmp_mpi_intel_linux_64_db : LFLAGS = -static-intel
openmp_mpi_intel_linux_64_db : FCOMPL = mpif77
openmp_mpi_intel_linux_64_db : FOPENMPFLAGS = -openmp -openmp-link static
openmp_mpi_intel_linux_64_db : obj = fds_openmp_mpi_intel_linux_64_db
openmp_mpi_intel_linux_64_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi)

# 32 bit

intel_osx_32 : FFLAGS = -O2  -m32 -heap-arrays -vec_report0 -ipo -mmacosx-version-min=10.4
intel_osx_32 : LFLAGS = -static-intel
intel_osx_32 : FCOMPL = ifort
intel_osx_32 : obj = fds_intel_osx_32
intel_osx_32 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

mpi_intel_osx_32 : FFLAGS  = -O2 -m32 -heap-arrays -vec_report0 -ipo -mmacosx-version-min=10.4
mpi_intel_osx_32 : LFLAGS  = -L$(MPIDIST)/lib -static-intel
mpi_intel_osx_32 : FCOMPL  = mpif77
mpi_intel_osx_32 : obj = fds_mpi_intel_osx_32
mpi_intel_osx_32 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_mpi)

intel_osx_32_db : FFLAGS = -m32 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -vec_report0 -stand:f03 -mmacosx-version-min=10.4
intel_osx_32_db : FCOMPL = ifort
intel_osx_32_db : obj = fds_intel_osx_32_db
intel_osx_32_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

mpi_intel_osx_32_db : FFLAGS  = -check -warn unused -O0 -m32 -g -traceback -fpe:0 -vec_report0 -mmacosx-version-min=10.4
mpi_intel_osx_32_db : LFLAGS  = -L$(MPIDIST)/lib -static-intel
mpi_intel_osx_32_db : FCOMPL  = mpif77
mpi_intel_osx_32_db : obj = fds_mpi_intel_osx_32_db
mpi_intel_osx_32_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_mpi)

# 64 bit

intel_osx_64 : FFLAGS = -O2 -m64 -heap-arrays -ipo -mmacosx-version-min=10.4
intel_osx_64 : LFLAGS = -static-intel
intel_osx_64 : FCOMPL = ifort
intel_osx_64 : obj = fds_intel_osx_64
intel_osx_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_serial)

mpi_intel_osx_64 : FFLAGS  = -O2 -m64 -heap-arrays -vec_report0 -ipo -mmacosx-version-min=10.4
mpi_intel_osx_64 : LFLAGS  = -L$(MPIDIST)/lib -static-intel
mpi_intel_osx_64 : FCOMPL  = mpif77
mpi_intel_osx_64 : obj = fds_mpi_intel_osx_64
mpi_intel_osx_64 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj_mpi)

openmp_intel_osx_64 : FFLAGS = -m64 -O2 -vec_report0 -ipo
openmp_intel_osx_64 : LFLAGS = -static-intel
openmp_intel_osx_64 : FCOMPL = ifort
openmp_intel_osx_64 : FOPENMPFLAGS = -openmp -openmp-link static
openmp_intel_osx_64 : obj = fds_openmp_intel_osx_64
openmp_intel_osx_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_serial)

intel_osx_64_db : FFLAGS = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -vec_report0 -stand:f03 -mmacosx-version-min=10.4
intel_osx_64_db : LFLAGS = -static-intel
intel_osx_64_db : FCOMPL = ifort
intel_osx_64_db : obj = fds_intel_osx_64_db
intel_osx_64_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

mpi_intel_osx_64_db : FFLAGS = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -vec_report0 -stand:f03 -mmacosx-version-min=10.4
mpi_intel_osx_64_db : LFLAGS = -L$(MPIDIST)/lib -static-intel
mpi_intel_osx_64_db : FCOMPL = mpif77
mpi_intel_osx_64_db : obj = fds_mpi_intel_osx_64_db
mpi_intel_osx_64_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

openmp_intel_osx_64_db : FFLAGS = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -vec_report0 -stand:f03
openmp_intel_osx_64_db : LFLAGS = -static-intel
openmp_intel_osx_64_db : FCOMPL = ifort
openmp_intel_osx_64_db : FOPENMPFLAGS = -openmp -openmp-link static
openmp_intel_osx_64_db : obj = fds_openmp_intel_osx_64_db
openmp_intel_osx_64_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_serial)

intel_osx_64_dv : FFLAGS = -m64 -O2 -vec_report0 -traceback -stand:f03 -mmacosx-version-min=10.4
intel_osx_64_dv : LFLAGS = -static-intel
intel_osx_64_dv : FCOMPL = ifort
intel_osx_64_dv : obj = fds_intel_osx_64_dv
intel_osx_64_dv : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

#*** IBM XLF Compiler ***

ibm_osx_ppc : FFLAGS = -O4 -qstrict -qextname -qarch=auto -qtune=auto
ibm_osx_ppc : FCOMPL = /opt/ibmcmp/xlf/8.1/bin/xlf95
ibm_osx_ppc : obj = fds_osx_ppc_ibm
ibm_osx_ppc : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

#*** Lahey Compiler ***

lahey_linux : FFLAGS = -O --tpp --staticlink --nap --nchk --npca --nsav --ntrace
lahey_linux : FCOMPL = lf95
lahey_linux : obj = fds_lahey
lahey_linux : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

lahey_linux_db : FFLAGS = --swm 3909 --wo --tpp --chk --f95 -g --warn --staticlink
lahey_linux_db : FCOMPL = lf95
lahey_linux_db : obj = fds_lahey_db
lahey_linux_db : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

lahey_linux_mpi : FFLAGS = -O --tpp --staticlink --nap --nchk --npca --nsav --ntrace
lahey_linux_mpi : FCOMPL = hf77
lahey_linux_mpi : INCLUDEDIR = /usr/local/include
lahey_linux_mpi : obj = fds_mpi_lahey
lahey_linux_mpi : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

#*** Sun Solaris 10 f95 Compiler ***

f95_solaris10 : FFLAGS   = -O2 -m64 -Bstatic
f95_solaris10 : FCOMPL   = f95
f95_solaris10 : obj = fds_solaris
f95_solaris10 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

f95_solaris10_mpi : INC_FLAGS = -I$(OPT)/SUNWhpc/HPC7.0/include/amd64
f95_solaris10_mpi : LIB_FLAGS = -Bstatic -lfui -Bdynamic -L$(OPT)/SUNWspro/prod/lib/amd64 -L$(OPT)/SUNWhpc/HPC7.0/lib/amd64
f95_solaris10_mpi : FFLAGS = -m64 $(INC_FLAGS)
f95_solaris10_mpi : FCOMPL = mpif90
f95_solaris10_mpi : obj = fds_mpi_solaris
f95_solaris10_mpi : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

#*** Sun Solaris (SunOS 5.11) f95 Compiler - AMD64 Opteron(tm) Processor *** 

f95_solaris11 : FFLAGS   = -xarch=sse2a -m64 -xchip=native -Bstatic 
f95_solaris11 : FCOMPL   = f95 
f95_solaris11 : obj = fds_solaris11 
f95_solaris11 : $(obj_serial) 
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial) 

f95_solaris11_mpi : INC_FLAGS=-I/opt/SUNWhpc/HPC8.2/sun/include/amd64 
f95_solaris11_mpi : LIB_FLAGS=-Bstatic -L/usr/lib/amd64/ld.so.1 -L/opt/sunstudio12.1/lib/amd64 -L/opt/SUNWhpc/HPC8.2/sun/lib/amd64 
f95_solaris11_mpi : FFLAGS = -xarch=sse2a -m64 -xchip=native $(INC_FLAGS) 
f95_solaris11_mpi : FCOMPL = mpif90 
f95_solaris11_mpi : obj = fds_mpi_solaris11 
f95_solaris11_mpi : $(obj_mpi) 
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi) 

sun_linux_64 : FFLAGS   = -O2 -m64 -Bstatic
sun_linux_64 : FCOMPL   = /opt/oracle/solstudio12.2/bin/sunf95
sun_linux_64 : obj = fds_sun_linux_64
sun_linux_64 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)


#
# to profile FDS using the sun compiler:
# 1.  build FDS using sun_linux_64_profile target
# 2.  run case as usual
# 3.  type: 
# /opt/oracle/solstudio12.2/bin/xprof_btoa fds_sun_linux_64_profile.profile profile_results.out

sun_linux_64_profile : FFLAGS   = -O2 -m64 -Bstatic -xprofile=collect
sun_linux_64_profile : FCOMPL   = /opt/oracle/solstudio12.2/bin/sunf95
sun_linux_64_profile : obj = fds_sun_linux_64_profile
sun_linux_64_profile : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

#*** Portland Group Compiler***

# AMD Opteron 64Bit PC
pgi_linux_64 : FFLAGS = -fast -tp k8-64 -Mstandard
pgi_linux_64 : FCOMPL = pgf90
pgi_linux_64 : LFLAGS  = -fast -tp k8-64 -byteswapio -Bstatic -lc -lgcc_eh
pgi_linux_64 : obj = fds_pgi_Opteron64
pgi_linux_64 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

# AMD Opteron 64Bit PC
mpi_pgi_linux_64: FFLAGS = -fast -tp k8-64 -Mstandard
mpi_pgi_linux_64: FCOMPL = mpif77
mpi_pgi_linux_64: LFLAGS  = -fast -tp k8-64 -byteswapio -Bstatic -lc -lgcc_eh
mpi_pgi_linux_64: obj = fds_mpi_pgi_Opteron64
mpi_pgi_linux_64: $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

#XLF-Compiler, IBM Power6 Core, OpenMP, full optimized settings
xlf_openmp : FFLAGS  = -O2 -qhot -q64 -qtune=pwr6 -qarch=pwr6 -qmaxmem=-1 -bdatapsize:64K -bstackpsize:64K -btextpsize:64K -qsmp=omp
xlf_openmp : FCOMPL  = xlf90_r
xlf_openmp : obj     = fds_xlf_openmp
xlf_openmp : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial) 

#XLF-Compiler, IBM Power6 Core, OpenMP, non-optimized setting
xlf_openmp_nooptimization : FFLAGS  = -O0 -q64 -qsmp=omp
xlf_openmp_nooptimization : FCOMPL  = xlf90_r
xlf_openmp_nooptimization : obj     = fds_xlf_openmp_nooptimization
xlf_openmp_nooptimization : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial) 

#AIX, JUMP, MPI-Version
aix_mpi : FFLAGS  = -O2 -q64
aix_mpi : FCOMPL  = mpxlf90
aix_mpi : obj     = fds_jump_mpi_64
aix_mpi : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

#*** GNU Compilers ***

#GCC and gfortran Version

gnu_win_32 : FFLAGS  = -O1 -static -ffree-line-length-none -ffree-form
gnu_win_32 : FCOMPL  = gfortran
gnu_win_32 : obj     = fds_gnu_win_32
gnu_win_32 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

gnu_linux_32 : FFLAGS  = -O1 -static -ffree-line-length-none -ffree-form
gnu_linux_32 : FCOMPL  = gfortran
gnu_linux_32 : obj     = fds_gnu_linux_32
gnu_linux_32 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

gnu_osx_32 : FFLAGS  = -O1 -ffree-line-length-none -ffree-form
gnu_osx_32 : FCOMPL  = gfortran
gnu_osx_32 : obj     = fds_gnu_osx_32
gnu_osx_32 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

gnu_mpi_32 : FFLAGS  = -O1
gnu_mpi_32 : FCOMPL  = mpif90
gnu_mpi_32 : obj = fds_gnu_mpi_32
gnu_mpi_32 : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

gnu_linux_64 : FFLAGS = -O2 -unroll -static -frecord-marker=4
gnu_linux_64 : FCOMPL = gfortran
gnu_linux_64 : obj = fds_gnu_linux_64
gnu_linux_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

xxxmpi_gnu_linux_64 : FFLAGS = -O2 -unroll -static -frecord-marker=4 
xxxmpi_gnu_linux_64 : FCOMPL = mpif77
xxxmpi_gnu_linux_64 : obj = fds_mpi_gnu_linux_64 
xxxmpi_gnu_linux_64 : setup $(obj_mpi) 

mpi_gnu_linux_64 : FFLAGS = -m64 -O2 -unroll
mpi_gnu_linux_64 : FCOMPL = mpif77
mpi_gnu_linux_64 : obj = fds_mpi_gnu_linux_64 
mpi_gnu_linux_64 : setup $(obj_mpi) 
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi) 


mpi_gnu_osx_32_db : FFLAGS = -O0 -ffree-line-length-256 -g -Wall -Wsurprising -W 
mpi_gnu_osx_32_db : FCOMPL = mpif90
mpi_gnu_osx_32_db : obj = fds_mpi_gnu_osx_32_db
mpi_gnu_osx_32_db : setup $(obj_mpi) 
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi) 

mpi_gnu_osx_32 : FFLAGS = -O2 -ffree-line-length-256
mpi_gnu_osx_32 : FCOMPL = mpif90
mpi_gnu_osx_32 : obj = fds_mpi_gnu_osx_32 
mpi_gnu_osx_32 : setup $(obj_mpi) 
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi) 


#*** End Compiler Specific Rules ***

# *** Object Dependencies ***

cons.o : prec.o
mpis.o : prec.o
type.o : prec.o cons.o
devc.o : prec.o
pois.o : prec.o cons.o
mesh.o : prec.o type.o
data.o : prec.o func.o cons.o
func.o : prec.o cons.o type.o mesh.o devc.o
smvv.o : prec.o func.o
ctrl.o : prec.o cons.o type.o mesh.o func.o
turb.o : func.o prec.o cons.o mesh.o
geom.o : func.o prec.o cons.o mesh.o
scrc.o : prec.o cons.o mesh.o func.o pois.o
irad.o : func.o prec.o cons.o type.o mesh.o
ieva.o : func.o prec.o cons.o type.o mesh.o
fire.o : func.o prec.o cons.o type.o mesh.o
hvac.o : func.o prec.o cons.o type.o mesh.o ctrl.o
wall.o : func.o prec.o cons.o type.o mesh.o hvac.o mass.o
velo.o : func.o prec.o cons.o type.o mesh.o turb.o geom.o
evac.o : func.o prec.o cons.o type.o mesh.o ieva.o data.o
divg.o : func.o prec.o cons.o type.o mesh.o scrc.o evac.o
mass.o : func.o prec.o cons.o type.o mesh.o turb.o
radi.o : func.o prec.o cons.o type.o mesh.o irad.o
pres.o : func.o prec.o cons.o type.o mesh.o pois.o velo.o scrc.o
part.o : func.o prec.o cons.o type.o mesh.o devc.o data.o geom.o
vege.o : func.o prec.o cons.o type.o mesh.o devc.o part.o
read.o : func.o prec.o cons.o type.o mesh.o devc.o ctrl.o irad.o evac.o hvac.o scrc.o data.o geom.o mass.o
init.o : func.o prec.o cons.o type.o mesh.o devc.o irad.o pois.o data.o evac.o
dump.o : func.o prec.o cons.o type.o mesh.o devc.o evac.o smvv.o turb.o scrc.o data.o radi.o
main.o : func.o prec.o cons.o type.o mesh.o devc.o smvv.o mass.o divg.o velo.o wall.o fire.o irad.o ieva.o radi.o evac.o part.o vege.o dump.o read.o init.o pres.o pois.o ctrl.o turb.o hvac.o scrc.o data.o geom.o

cons.obj : prec.obj
mpis.obj : prec.obj
type.obj : prec.obj cons.obj
devc.obj : prec.obj
pois.obj : prec.obj cons.obj
mesh.obj : prec.obj type.obj
data.obj : prec.obj func.obj cons.obj
func.obj : prec.obj cons.obj type.obj mesh.obj devc.obj
smvv.obj : prec.obj func.obj
ctrl.obj : prec.obj cons.obj type.obj mesh.obj func.obj
turb.obj : func.obj prec.obj cons.obj mesh.obj
geom.obj : func.obj prec.obj cons.obj mesh.obj
scrc.obj : prec.obj cons.obj mesh.obj func.obj pois.obj
irad.obj : func.obj prec.obj cons.obj type.obj mesh.obj
ieva.obj : func.obj prec.obj cons.obj type.obj mesh.obj
fire.obj : func.obj prec.obj cons.obj type.obj mesh.obj
hvac.obj : func.obj prec.obj cons.obj type.obj mesh.obj ctrl.obj
wall.obj : func.obj prec.obj cons.obj type.obj mesh.obj hvac.obj mass.obj
velo.obj : func.obj prec.obj cons.obj type.obj mesh.obj turb.obj geom.obj
evac.obj : func.obj prec.obj cons.obj type.obj mesh.obj ieva.obj data.obj 
divg.obj : func.obj prec.obj cons.obj type.obj mesh.obj scrc.obj evac.obj
mass.obj : func.obj prec.obj cons.obj type.obj mesh.obj turb.obj
radi.obj : func.obj prec.obj cons.obj type.obj mesh.obj irad.obj
pres.obj : func.obj prec.obj cons.obj type.obj mesh.obj pois.obj velo.obj scrc.obj
part.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj data.obj geom.obj
vege.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj part.obj
read.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj ctrl.obj irad.obj evac.obj hvac.obj scrc.obj data.obj geom.obj mass.obj
init.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj irad.obj pois.obj data.obj evac.obj
dump.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj evac.obj smvv.obj turb.obj scrc.obj data.obj radi.obj
main.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj smvv.obj mass.obj divg.obj velo.obj wall.obj fire.obj irad.obj ieva.obj radi.obj evac.obj part.obj vege.obj dump.obj read.obj init.obj pres.obj pois.obj ctrl.obj turb.obj hvac.obj scrc.obj data.obj geom.obj


#*** Clean Target to remove Object and Module files ***

.PHONY : clean
clean:
	-rm -f *.o *.mod *.obj 

.PHONY : winclean
winclean:
	-erase *.mod *.obj 
