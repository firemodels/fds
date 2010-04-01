# Makefile for Fire Dynamics Simulator (FDS), Version 5

# To use this makefile, select the appropriate compiler and OS from the list below. For example, to
# compile the 64 bit parallel version (MPI) of FDS under Linux using the Intel compilers type:
#
# make intel_linux_mpi_64

# A make utility for windows may be downloaded from http://gnuwin32.sourceforge.net/packages/make.htm

# The variable VPATH defines where the source code is relative to the current directory

VPATH = ../FDS_Source

# Definition of the non-MPI or MPI object variables

# note: objf_serial and objf_mpi contain only Fortran routines ie not isob.c

objf_serial = prec.o smvv.o cons.o devc.o type.o mesh.o func.o irad.o \
             ieva.o pois.o radi.o evac.o part.o vege.o ctrl.o turb.o dump.o hvac.o read.o \
             mass.o wall.o fire.o divg.o init.o velo.o pres.o mpis.o main.o
obj_serial = isob.o $(objf_serial)

objwinf_serial = $(objf_serial:.o=.obj)
objwin_serial = $(obj_serial:.o=.obj)

objf_mpi = prec.o smvv.o cons.o devc.o type.o mesh.o func.o irad.o \
          ieva.o pois.o radi.o evac.o part.o vege.o ctrl.o turb.o dump.o hvac.o read.o \
          mass.o wall.o fire.o divg.o init.o velo.o pres.o mpip.o main.o
obj_mpi = isob.o $(objf_mpi)

objwinf_mpi = $(objf_mpi:.o=.obj)
objwin_mpi = $(obj_mpi:.o=.obj)

# General Purpose Rules

no_target:
	@echo \******** You did not specify a make target \********
	@echo Please read the comments at the top of the makefile

setup:
%.o : %.mod

setup_win:
%.obj : %.mod

.SUFFIXES: .c .f90 .o .obj

.f90.o:
	$(FCOMPL) -c $(FFLAGS) $<
.f90.obj:
	$(FCOMPL) -c $(FFLAGS) $<
.c.o:
	$(CCOMPL) -c $(CFLAGS) $<
.c.obj:
	$(CCOMPL) -c $(CFLAGS) $<

# Windows Intel Compiler, version 11

 #*** note: the following scripts (located in ..\Scripts directory) must be run before running make with this target
 #    ..\Scripts\iclvars ia32
 #    ..\Scripts\ifortvars ia32

intel_win_32_db : FFLAGS = /Qvec_report0 /nologo /debug:full /Od /extend_source:132 /warn:all /warn:nointerfaces /Qtrapuv /fp:strict /fp:except /traceback /check:all
intel_win_32_db : CFLAGS = -Dpp_noappend
intel_win_32_db : FCOMPL = ifort
intel_win_32_db : CCOMPL = cl
intel_win_32_db : obj = fds5_win_32_db
intel_win_32_db : setup_win $(objwinf_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F100000000 $(objwinf_serial) ../win32lib/isob.obj

intel_win_32 : FFLAGS = -O2 /Qvec_report0
intel_win_32 : CFLAGS = -O2 -Dpp_noappend
intel_win_32 : FCOMPL = ifort
intel_win_32 : CCOMPL = cl
intel_win_32 : obj = fds5_win_32
intel_win_32 : setup_win $(objwinf_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F100000000 $(objwinf_serial) ../win32lib/isob.obj

intel_win_32all : FFLAGS = -O2 /Qvec_report0
intel_win_32all : CFLAGS = -O2 -Dpp_noappend
intel_win_32all : FCOMPL = ifort
intel_win_32all : CCOMPL = cl
intel_win_32all : obj = fds5_win_32
intel_win_32all : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F100000000 $(objwin_serial)
	
mpi_intel_win_32 : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
mpi_intel_win_32 : MPIINCLUDE = "c:\program files\mpich2\include"
mpi_intel_win_32 : FFLAGS = -O2 /I $(MPIINCLUDE)
mpi_intel_win_32 : CFLAGS = -O2 -Dpp_noappend
mpi_intel_win_32 : FCOMPL = ifort
mpi_intel_win_32 : CCOMPL = cl
mpi_intel_win_32 : obj = fds5_mpi_win_32
mpi_intel_win_32 : setup_win $(objwinf_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F100000000 $(objwinf_mpi) ../win32lib/isob.obj $(MPILIB)

 #*** note: the following scripts (located in ..\Scripts directory) must be run before running make with this target
 #    ..\Scripts\iclvars intel64
 #    ..\Scripts\ifortvars intel64

intel_win_64 : FFLAGS = -O2 /Qvec_report0
intel_win_64 : CFLAGS = -O2 -Dpp_noappend
intel_win_64 : FCOMPL = ifort
intel_win_64 : CCOMPL = cl
intel_win_64 : obj = fds5_win_64
intel_win_64 : setup_win $(objwinf_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F1000000000 $(objwinf_serial) ../win64lib/isob.obj

mpi_intel_win_64 : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
mpi_intel_win_64 : MPIINCLUDE = "c:\program files\mpich2\include"
mpi_intel_win_64 : FFLAGS = -O2 /Qvec_report0 /I $(MPIINCLUDE)
mpi_intel_win_64 : CFLAGS = -O2 -Dpp_noappend 
mpi_intel_win_64 : FCOMPL = ifort
mpi_intel_win_64 : CCOMPL = cl
mpi_intel_win_64 : obj = fds5_mpi_win_64
mpi_intel_win_64 : setup_win $(objwinf_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F1000000000 $(objwinf_mpi) ../win64lib/isob.obj $(MPILIB)

# Intel Compiler for Linux, version 11

intel_linux_32 : FFLAGS = -m32 -O3 -static -vec_report0
intel_linux_32 : CFLAGS = -m32 -O -Dpp_noappend -vec_report0
intel_linux_32 : FCOMPL = ifort
intel_linux_32 : CCOMPL = icc
intel_linux_32 : obj = fds5_intel_linux_32
intel_linux_32 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_linux_64 : FFLAGS = -m64 -O3 -static -vec_report0
intel_linux_64 : CFLAGS = -m64 -O -Dpp_noappend -vec_report0
intel_linux_64 : FCOMPL = ifort
intel_linux_64 : CCOMPL = icc
intel_linux_64 : obj = fds5_intel_linux_64
intel_linux_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

mpi_intel_linux_32 : FFLAGS = -m32 -O3 -vec_report0
mpi_intel_linux_32 : CFLAGS = -m32 -O -Dpp_noappend -vec_report0
mpi_intel_linux_32 : FCOMPL = mpif77
mpi_intel_linux_32 : CCOMPL = icc
mpi_intel_linux_32 : obj = fds5_mpi_intel_linux_32
mpi_intel_linux_32 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

mpi_intel_linux_64 : FFLAGS = -m64 -O3 -vec_report0
mpi_intel_linux_64 : CFLAGS = -m64 -O -Dpp_noappend -vec_report0
mpi_intel_linux_64 : FCOMPL = mpif77
mpi_intel_linux_64 : CCOMPL = icc
mpi_intel_linux_64 : obj = fds5_mpi_intel_linux_64
mpi_intel_linux_64 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

intel_linux_32_db : FFLAGS = -m32 -check all -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -static -vec_report0 -ftrapuv -debug full
intel_linux_32_db : CFLAGS = -m32 -Dpp_noappend -vec_report0
intel_linux_32_db : FCOMPL = ifort
intel_linux_32_db : CCOMPL = icc
intel_linux_32_db : obj = fds5_intel_linux_32_db
intel_linux_32_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_linux_64_db : FFLAGS = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -static -vec_report0
intel_linux_64_db : CFLAGS = -m64 -Dpp_noappend -vec_report0
intel_linux_64_db : FCOMPL = ifort
intel_linux_64_db : CCOMPL = icc
intel_linux_64_db : obj = fds5_intel_linux_64_db
intel_linux_64_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

mpi_intel_linux_32_db : FFLAGS = -m32 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency  -vec_report0
mpi_intel_linux_32_db : CFLAGS = -m32 -Dpp_noappend -vec_report0
mpi_intel_linux_32_db : FCOMPL = mpif77
mpi_intel_linux_32_db : CCOMPL = icc
mpi_intel_linux_32_db : obj = fds5_mpi_intel_linux_32_db
mpi_intel_linux_32_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

mpi_intel_linux_64_db : FFLAGS = -m64 -check -warn unused -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -vec_report0
mpi_intel_linux_64_db : CFLAGS = -m64 -O -Dpp_noappend -vec_report0
mpi_intel_linux_64_db : FCOMPL = mpif77
mpi_intel_linux_64_db : CCOMPL = icc
mpi_intel_linux_64_db : obj = fds5_mpi_intel_linux_64_db
mpi_intel_linux_64_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

# 32 bit
intel_osx_32 : FFLAGS = -O3  -m32 -axSSSE3 -heap-arrays -static-intel -vec_report0
intel_osx_32 : CFLAGS = -O3 -m32 -Dpp_noappend -Dpp_OSX -vec_report0
intel_osx_32 : FCOMPL = ifort
intel_osx_32 : CCOMPL = icc
intel_osx_32 : obj = fds5_intel_osx_32
intel_osx_32 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

mpi_intel_osx_32 : FFLAGS  = -O3 -m32 -axSSSE3 -heap-arrays -static-intel -L/opt/openmpi32/lib -vec_report0
mpi_intel_osx_32 : CFLAGS  = -O3 -m32 -Dpp_noappend -Dpp_OSX -vec_report0
mpi_intel_osx_32 : FCOMPL  = /opt/openmpi32/bin/mpif77
mpi_intel_osx_32 : CCOMPL  = icc
mpi_intel_osx_32 : obj = fds5_mpi_intel_osx_32
mpi_intel_osx_32 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

intel_osx_32_db : FFLAGS = -check -warn unused -O0 -m32 -g -traceback -fpe:0 -fltconsistency -FR -WB -vec_report0
intel_osx_32_db : CFLAGS = -O0 -m32 -g -traceback -Dpp_noappend -Dpp_OSX -vec_report0
intel_osx_32_db : FCOMPL = ifort
intel_osx_32_db : CCOMPL = icc
intel_osx_32_db : obj = fds5_intel_osx_32_db
intel_osx_32_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_osx_mpi_32_db : FFLAGS  = -check -warn unused -O0 -m32 -g -traceback -fpe:0 -vec_report0
intel_osx_mpi_32_db : CFLAGS  = -O0 -m32 -g -traceback -Dpp_noappend -Dpp_OSX
intel_osx_mpi_32_db : FCOMPL  = /opt/openmpi32/bin/mpif77
intel_osx_mpi_32_db : CCOMPL  = icc
intel_osx_mpi_32_db : obj = fds5_mpi_intel_osx_32_db
intel_osx_mpi_32_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

# 64 bit
intel_osx_64 : FFLAGS = -O3 -m64 -heap-arrays -axSSSE3 -static-intel
intel_osx_64 : CFLAGS = -O3 -m64 -Dpp_noappend -Dpp_OSX -vec_report0
intel_osx_64 : FCOMPL = ifort
intel_osx_64 : CCOMPL = icc
intel_osx_64 : obj = fds5_intel_osx_64
intel_osx_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

mpi_intel_osx_64 : FFLAGS  = -O3 -m64 -heap-arrays -axSSSE3 -static-intel -L/opt/openmpi64/lib -vec_report0
mpi_intel_osx_64 : CFLAGS  = -O3 -m64 -Dpp_noappend -Dpp_OSX -vec_report0
mpi_intel_osx_64 : FCOMPL  = /opt/openmpi64/bin/mpif77
mpi_intel_osx_64 : CCOMPL  = icc
mpi_intel_osx_64 : obj = fds5_mpi_intel_osx_64
mpi_intel_osx_64 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

intel_osx_64_db : FFLAGS = -check -warn unused -O0 -m64 -g -traceback -save-temps -fpe:0 -fltconsistency -fp-stack-check -FR -WB -vec_report0
intel_osx_64_db : CFLAGS = -O0 -g -traceback -m64 -Dpp_noappend -Dpp_OSX -vec_report0
intel_osx_64_db : FCOMPL = ifort
intel_osx_64_db : CCOMPL = icc
intel_osx_64_db : obj = fds5_intel_osx_64_db
intel_osx_64_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_osx_mpi_64_db : FFLAGS  = -check -warn unused -O0 -g -traceback -m64 -fpe:0 -fp-stack-check -save-temps -vec_report0
intel_osx_mpi_64_db : CFLAGS  = -O0 -g -traceback -m64 -Dpp_noappend -Dpp_OSX -vec_report0
intel_osx_mpi_64_db : FCOMPL  = /opt/openmpi64/bin/mpif77
intel_osx_mpi_64_db : CCOMPL  = icc
intel_osx_mpi_64_db : obj = fds5_mpi_intel_osx_64_db
intel_osx_mpi_64_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

#*** IBM XLF Compiler ***

ibm_osx_ppc : FFLAGS = -O4 -qstrict -qextname -qarch=auto -qtune=auto
ibm_osx_ppc : CFLAGS = -O3 -D pp_OSX
ibm_osx_ppc : FCOMPL = /opt/ibmcmp/xlf/8.1/bin/xlf95
ibm_osx_ppc : CCOMPL = gcc
ibm_osx_ppc : obj = fds5_osx_ppc_ibm
ibm_osx_ppc : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

#*** Lahey Compiler ***

lahey_linux : FFLAGS = -O --tpp --staticlink --nap --nchk --npca --nsav --ntrace
lahey_linux : CFLAGS = -O
lahey_linux : FCOMPL = lf95
lahey_linux : CCOMPL = gcc
lahey_linux : obj = fds5_lahey
lahey_linux : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

lahey_linux_db : FFLAGS = --swm 3909 --wo --tpp --chk --f95 -g --warn --staticlink
lahey_linux_db : CFLAGS = -O
lahey_linux_db : FCOMPL = lf95
lahey_linux_db : CCOMPL = gcc
lahey_linux_db : obj = fds5_lahey_db
lahey_linux_db : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

lahey_linux_mpi : FFLAGS = -O --tpp --staticlink --nap --nchk --npca --nsav --ntrace
lahey_linux_mpi : CFLAGS = -O
lahey_linux_mpi : FCOMPL = hf77
lahey_linux_mpi : CCOMPL = gcc
lahey_linux_mpi : INCLUDEDIR = /usr/local/include
lahey_linux_mpi : obj = fds5_mpi_lahey
lahey_linux_mpi : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

#*** Sun Solaris 10 f95 Compiler ***

f95_solaris10 : FFLAGS   = -O2 -m64 -Bstatic
f95_solaris10 : CFLAGS   = -O2 -m64
f95_solaris10 : FCOMPL   = f95
f95_solaris10 : CCOMPL   = cc
f95_solaris10 : obj = fds5_solaris
f95_solaris10 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

f95_solaris10_mpi : INC_FLAGS = -I$(OPT)/SUNWhpc/HPC7.0/include/amd64
f95_solaris10_mpi : LIB_FLAGS = -Bstatic -lfui -Bdynamic -L$(OPT)/SUNWspro/prod/lib/amd64 -L$(OPT)/SUNWhpc/HPC7.0/lib/amd64
f95_solaris10_mpi : FFLAGS = -m64 $(INC_FLAGS)
f95_solaris10_mpi : CFLAGS = -m64
f95_solaris10_mpi : FCOMPL = mpif90
f95_solaris10_mpi : CCOMPL = cc
f95_solaris10_mpi : obj = fds5_mpi_solaris
f95_solaris10_mpi : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

#*** Sun Solaris (SunOS 5.11) f95 Compiler - AMD64 Opteron(tm) Processor *** 

f95_solaris11 : FFLAGS   = -xarch=sse2a -m64 -xchip=native -Bstatic 
f95_solaris11 : CFLAGS   = -xarch=sse2a -m64 -xchip=native 
f95_solaris11 : FCOMPL   = f95 
f95_solaris11 : CCOMPL   = cc 
f95_solaris11 : obj = fds5_solaris11 
f95_solaris11 : $(obj_serial) 
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial) 

f95_solaris11_mpi : INC_FLAGS=-I/opt/SUNWhpc/HPC8.2/sun/include/amd64 
f95_solaris11_mpi : LIB_FLAGS=-Bstatic -L/usr/lib/amd64/ld.so.1 -L/opt/sunstudio12.1/lib/amd64 -L/opt/SUNWhpc/HPC8.2/sun/lib/amd64 
f95_solaris11_mpi : FFLAGS = -xarch=sse2a -m64 -xchip=native $(INC_FLAGS) 
f95_solaris11_mpi : CFLAGS = -xarch=sse2a -m64 -xchip=native 
f95_solaris11_mpi : FCOMPL = mpif90 
f95_solaris11_mpi : CCOMPL = cc 
f95_solaris11_mpi : obj = fds5_mpi_solaris11 
f95_solaris11_mpi : $(obj_mpi) 
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi) 

#*** Portland Group Compiler***

# AMD Opteron 64Bit PC
pgi_linux_64 : FFLAGS = -fast -tp k8-64 -Mstandard
pgi_linux_64 : CFLAGS = -fast -tp k8-64
pgi_linux_64 : FCOMPL = pgf90
pgi_linux_64 : CCOMPL = pgcc
pgi_linux_64 : LFLAGS  = -fast -tp k8-64 -byteswapio -Bstatic -lc -lgcc_eh
pgi_linux_64 : obj = fds5_pgi_Opteron64
pgi_linux_64 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

# AMD Opteron 64Bit PC
mpi_pgi_linux_64: FFLAGS = -fast -tp k8-64 -Mstandard
mpi_pgi_linux_64: CFLAGS = -fast -tp k8-64
mpi_pgi_linux_64: FCOMPL = mpif77
mpi_pgi_linux_64: CCOMPL = mpicc
mpi_pgi_linux_64: LFLAGS  = -fast -tp k8-64 -byteswapio -Bstatic -lc -lgcc_eh
mpi_pgi_linux_64: obj = fds5_mpi_pgi_Opteron64
mpi_pgi_linux_64: $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

#XLF-Compiler, IBM Power6 Core, OpenMP, full optimized settings
xlf_openmp : FFLAGS  = -O3 -qhot -q64 -qtune=pwr6 -qarch=pwr6 -qmaxmem=-1 -bdatapsize:64K -bstackpsize:64K -btextpsize:64K -qsmp=omp
xlf_openmp : CFLAGS  = -O3 -qhot -Dpp_noappend -q64 -qtune=pwr6 -qarch=pwr6 -qmaxmem=-1 -bdatapsize:64K -bstackpsize:64K -btextpsize:64K -qsmp=omp
xlf_openmp : FCOMPL  = xlf90_r
xlf_openmp : CCOMPL  = xlc_r
xlf_openmp : obj     = fds5_xlf_openmp
xlf_openmp : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial) 

#XLF-Compiler, IBM Power6 Core, OpenMP, non-optimized setting
xlf_openmp_nooptimization : FFLAGS  = -O0 -q64 -qsmp=omp
xlf_openmp_nooptimization : CFLAGS  = -O0 -Dpp_noappend -q64 -qsmp=omp
xlf_openmp_nooptimization : FCOMPL  = xlf90_r
xlf_openmp_nooptimization : CCOMPL  = xlc_r
xlf_openmp_nooptimization : obj     = fds5_xlf_openmp_nooptimization
xlf_openmp_nooptimization : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial) 

#AIX, JUMP, MPI-Version
aix_mpi : FFLAGS  = -O3 -q64
aix_mpi : CFLAGS  = -O3 -Dpp_noappend -q64
aix_mpi : FCOMPL  = mpxlf90
aix_mpi : CCOMPL  = mpcc
aix_mpi : obj     = fds5_jump_mpi_64
aix_mpi : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

#*** GNU Compilers ***

#GCC and gfortran Version

gnu_win_32 : FFLAGS  = -O1 -ffree-line-length-none -ffree-form
gnu_win_32 : CFLAGS  = -O1
gnu_win_32 : FCOMPL  = gfortran
gnu_win_32 : CCOMPL  = gcc
gnu_win_32 : obj     = fds5_gnu_win_32
gnu_win_32 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

gnu_linux_32 : FFLAGS  = -O1 -ffree-line-length-none -ffree-form
gnu_linux_32 : CFLAGS  = -O1
gnu_linux_32 : FCOMPL  = gfortran
gnu_linux_32 : CCOMPL  = gcc
gnu_linux_32 : obj     = fds5_gnu_linux_32
gnu_linux_32 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

gnu_mpi_32 : FFLAGS  = -O1
gnu_mpi_32 : CFLAGS  = -O1 
gnu_mpi_32 : FCOMPL  = mpif90
gnu_mpi_32 : CCOMPL  = mpicc
gnu_mpi_32 : obj = fds5_gnu_mpi_32
gnu_mpi_32 : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

gnu_linux_64 : FFLAGS = -O3 -unroll -static -frecord-marker=4
gnu_linux_64 : CFLAGS = -O3 -unroll -static
gnu_linux_64 : FCOMPL = gfortran
gnu_linux_64 : CCOMPL = gcc
gnu_linux_64 : obj = fds5_gnu_linux_64
gnu_linux_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

xxxmpi_gnu_linux_64 : FFLAGS = -O3 -unroll -static -frecord-marker=4 
xxxmpi_gnu_linux_64 : CFLAGS = -O3 -unroll -static 
xxxmpi_gnu_linux_64 : FCOMPL = mpif77
xxxmpi_gnu_linux_64 : CCOMPL = mpicc 
xxxmpi_gnu_linux_64 : obj = fds5_mpi_gnu_linux_64 
xxxmpi_gnu_linux_64 : setup $(obj_mpi) 

mpi_gnu_linux_64 : FFLAGS = -m64 -O2 -unroll
mpi_gnu_linux_64 : CFLAGS = -m64 -O2 -unroll
mpi_gnu_linux_64 : FCOMPL = mpif77
mpi_gnu_linux_64 : CCOMPL = mpicc 
mpi_gnu_linux_64 : obj = fds5_mpi_gnu_linux_64 
mpi_gnu_linux_64 : setup $(obj_mpi) 
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi) 


#*** End Compiler Specific Rules ***

# *** Object Dependencies ***

smvv.o : isob.o
cons.o : prec.o
type.o : prec.o cons.o
devc.o : prec.o
pois.o : prec.o
mesh.o : prec.o type.o
func.o : prec.o cons.o type.o mesh.o
turb.o : func.o prec.o cons.o mesh.o
ctrl.o : prec.o cons.o type.o mesh.o func.o
irad.o : func.o prec.o cons.o type.o mesh.o
ieva.o : func.o prec.o cons.o type.o mesh.o
fire.o : func.o prec.o cons.o type.o mesh.o
wall.o : func.o prec.o cons.o type.o mesh.o
velo.o : func.o prec.o cons.o type.o mesh.o turb.o
divg.o : func.o prec.o cons.o type.o mesh.o
mass.o : func.o prec.o cons.o type.o mesh.o turb.o
radi.o : func.o prec.o cons.o type.o mesh.o irad.o
evac.o : func.o prec.o cons.o type.o mesh.o ieva.o
pres.o : func.o prec.o cons.o type.o mesh.o pois.o velo.o
part.o : func.o prec.o cons.o type.o mesh.o devc.o
vege.o : func.o prec.o cons.o type.o mesh.o devc.o part.o
hvac.o : func.o prec.o cons.o type.o mesh.o ctrl.o
read.o : func.o prec.o cons.o type.o mesh.o devc.o irad.o evac.o hvac.o
init.o : func.o prec.o cons.o type.o mesh.o devc.o irad.o pois.o
dump.o : func.o prec.o cons.o type.o mesh.o devc.o evac.o smvv.o isob.c turb.o
main.o : func.o prec.o cons.o type.o mesh.o devc.o smvv.o isob.c mass.o divg.o velo.o wall.o fire.o irad.o ieva.o radi.o evac.o part.o vege.o dump.o read.o init.o pres.o pois.o ctrl.o turb.o hvac.o 

cons.obj : prec.obj
type.obj : prec.obj cons.obj
devc.obj : prec.obj
pois.obj : prec.obj
mesh.obj : prec.obj type.obj
func.obj : prec.obj cons.obj type.obj mesh.obj
turb.obj : func.obj prec.obj cons.obj mesh.obj
ctrl.obj : prec.obj cons.obj type.obj mesh.obj func.obj
irad.obj : func.obj prec.obj cons.obj type.obj mesh.obj
ieva.obj : func.obj prec.obj cons.obj type.obj mesh.obj
fire.obj : func.obj prec.obj cons.obj type.obj mesh.obj
wall.obj : func.obj prec.obj cons.obj type.obj mesh.obj
velo.obj : func.obj prec.obj cons.obj type.obj mesh.obj turb.obj
divg.obj : func.obj prec.obj cons.obj type.obj mesh.obj
mass.obj : func.obj prec.obj cons.obj type.obj mesh.obj turb.obj
radi.obj : func.obj prec.obj cons.obj type.obj mesh.obj irad.obj
evac.obj : func.obj prec.obj cons.obj type.obj mesh.obj ieva.obj
pres.obj : func.obj prec.obj cons.obj type.obj mesh.obj pois.obj velo.obj
part.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj
vege.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj part.obj
hvac.obj : func.obj prec.obj cons.obj type.obj mesh.obj ctrl.obj
read.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj ctrl.obj irad.obj evac.obj hvac.obj
init.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj irad.obj pois.obj
dump.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj evac.obj smvv.obj isob.c turb.obj
main.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj smvv.obj isob.c mass.obj divg.obj velo.obj wall.obj fire.obj irad.obj ieva.obj radi.obj evac.obj part.obj vege.obj dump.obj read.obj init.obj pres.obj pois.obj ctrl.obj turb.obj hvac.obj


#*** Clean Target to remove Object and Module files ***

.PHONY : clean
clean:
	-rm -f *.o *.mod *.obj 

.PHONY : winclean
winclean:
	-erase *.mod *.obj 
