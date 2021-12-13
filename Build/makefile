# Makefile for Fire Dynamics Simulator (FDS)

# To use this makefile, select the appropriate compiler and OS from the list below. For example, to
# compile the 64 bit parallel version (MPI) of FDS under Linux using the Intel compilers type:
#
# make mpi_intel_linux_64

# A make utility for windows may be downloaded from http://gnuwin32.sourceforge.net/packages/make.htm

# The variable VPATH defines where the source code is relative to the current directory

VPATH = ../Source

ifeq ($(shell echo "check_quotes"),"check_quotes")
# windows
  GIT_HASH := $(shell ..\Scripts\githash)
  GIT_STAT := $(shell ..\Scripts\gitstat)
  GIT_DATE := $(shell ..\Scripts\gitlog)
  GIT_BRANCH := $(shell ..\Scripts\gitbranch)
  BUILD_DATE := $(shell date /t) $(shell time /t)
else
# linux/osx
  SHELL:=/bin/bash
  GIT_HASH := $(shell git describe --long)
  ifeq ("$(shell git diff --shortstat $(VPATH)/*.f90)","")
    GIT_STAT :=
  else
    GIT_STAT := -dirty
  endif
  GIT_DATE := $(shell git log -1 --format=%cd)
  GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD)
  BUILD_DATE := $(shell date "+%b %d, %Y  %T")
# The preprocessor used by xlf uses comma as separation between preprocessing commands.
# Commas needs to be escaped if they are to be used as part of a value in a preprocessing variable.
  BUILD_DATE_XLF := $(shell date "+%b %d\, %Y  %T")
endif

GITINFO_BASE=-DGITHASH_PP=\"$(GIT_HASH)$(GIT_STAT)-$(GIT_BRANCH)\" -DGITDATE_PP=\""$(GIT_DATE)\""
GITINFO=-fpp $(GITINFO_BASE) -DBUILDDATE_PP=\""$(BUILD_DATE)\""
GITINFOGNU=-cpp $(GITINFO_BASE) -DBUILDDATE_PP=\""$(BUILD_DATE)\""
GITINFOXLF=-qsuffix=cpp=f90 '-WF,-DGITHASH_PP=\"$(GIT_HASH)\",-DGITDATE_PP=\"$(GIT_DATE)\",-DBUILDDATE_PP=\"$(BUILD_DATE_XLF)\"'
GITINFOPGI=-Mpreprocess $(GITINFO_BASE) -DBUILDDATE_PP=\""$(BUILD_DATE)\""
GFORTRAN_OPTIONS=

ifeq ($(shell echo "check_quotes"),"check_quotes")
# windows
  INTELMPI_COMPVERSION := $(shell ..\Scripts\intel_compversion)
  GNU_COMPVERSION := $(shell ..\Scripts\gnu_compversion)
else
# linux/osx
  INTELMPI_COMPVERSION := $(shell ../Scripts/intelmpi_compversion.sh)
  INTELMAJOR_COMPVERSION := $(shell ../Scripts/intelmajor_compversion.sh)
  OPENMPI_COMPVERSION := $(shell ../Scripts/openmpi_compversion.sh)
  GNU_COMPVERSION := $(shell ../Scripts/gnu_compversion.sh)
  ifeq ($(shell ../Scripts/gnu_compmajorversion.sh),10)
    GFORTRAN_OPTIONS=-fallow-argument-mismatch
  endif
endif

INTELMPI_COMPINFO=-DCOMPVER_PP=\"$(INTELMPI_COMPVERSION)\"
OPENMPI_COMPINFO=-DCOMPVER_PP=\"$(OPENMPI_COMPVERSION)\"
GNU_COMPINFO=-DCOMPVER_PP=\"$(GNU_COMPVERSION)\"

FFLAGSMKL_INTEL =
LFLAGSMKL_INTEL =
FFLAGSMKL_OPENMPI =
LFLAGSMKL_OPENMPI =
FFLAGSMKL_CUSTOM =
LFLAGSMKL_CUSTOM =
ifdef MKLROOT # This assumes the MKL library is available.
ifeq ($(shell echo "check_quotes"),"check_quotes")
# windows
ifeq ($(shell ..\Scripts\havemkl),1)
# Use cluster_sparse_solver with intelmpi blacs
  FFLAGSMKL_INTEL = -DWITH_MKL /I"$(MKLROOT)"\include
  LFLAGSMKL_INTEL = "$(MKLROOT)"\lib\intel64\mkl_intel_lp64.lib "$(MKLROOT)"\lib\intel64\mkl_intel_thread.lib "$(MKLROOT)"\lib\intel64\mkl_core.lib "$(MKLROOT)"\lib\intel64\mkl_blacs_intelmpi_lp64.lib
endif
else
# linux/osx
ifneq ("$(wildcard ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a)","")
# Use cluster_sparse_solver with intelmpi blacs
  FFLAGSMKL_INTEL = -DWITH_MKL -I$(MKLROOT)/include
  LFLAGSMKL_INTEL = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl
endif
ifneq ("$(wildcard ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a)","")
# Use cluster_sparse_solver with openmpi blacs
  FFLAGSMKL_OPENMPI = -DWITH_MKL -I$(MKLROOT)/include
  LFLAGSMKL_OPENMPI = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl
  FFLAGSMKL_GNU_OPENMPI = -DWITH_MKL -I$(MKLROOT)/include
  LFLAGSMKL_GNU_OPENMPI =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lgomp -lpthread -lm -ldl
endif
ifneq ("$(wildcard ${MKLROOT}/lib/libmkl_blacs_custom_lp64.a)","")
# Use cluster_sparse_solver with custom blacs
  FFLAGSMKL_CUSTOM = -DWITH_MKL -I$(MKLROOT)/include
# MKL_LIBS needs to occur 3 times in order to avoid unresolved external reference warnings
  MKL_LIBS = ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_core.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_blacs_custom_lp64.a
  LFLAGSMKL_CUSTOM = ${MKL_LIBS} ${MKL_LIBS} ${MKL_LIBS} -lpthread -lm -ldl
  FFLAGSMKL_GNU_CUSTOM = -DWITH_MKL -I$(MKLROOT)/include
  LFLAGSMKL_GNU_CUSTOM =  ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a ${MKLROOT}/lib/libmkl_blacs_custom_lp64.a -lpthread -lm -ldl
endif
endif # OS test
endif # MKLROOT test

obj_mpi = prec.o cons.o devc.o type.o data.o mesh.o func.o gsmv.o smvv.o rcal.o turb.o soot.o \
          pois.o geom.o ccib.o radi.o part.o vege.o ctrl.o hvac.o mass.o \
          wall.o fire.o velo.o scrc.o pres.o init.o dump.o read.o divg.o main.o

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

# Files which have OpenMP compiler directives are compiled using the FOPENMPFLAG and can be enabled/disabled on a per-file basis.

divg.o: FFLAGS += $(FOPENMPFLAGS)
func.o: FFLAGS += $(FOPENMPFLAGS)
mass.o: FFLAGS += $(FOPENMPFLAGS)
pois.o: FFLAGS += $(FOPENMPFLAGS)
pres.o: FFLAGS += $(FOPENMPFLAGS)
turb.o: FFLAGS += $(FOPENMPFLAGS)
velo.o: FFLAGS += $(FOPENMPFLAGS)
rcal.o: FFLAGS += $(FOPENMPFLAGS)
radi.o: FFLAGS += $(FOPENMPFLAGS)
scrc.o: FFLAGS += $(FOPENMPFLAGS)
main.o: FFLAGS += $(FOPENMPFLAGS)

divg.obj: FFLAGS += $(FOPENMPFLAGS)
func.obj: FFLAGS += $(FOPENMPFLAGS)
mass.obj: FFLAGS += $(FOPENMPFLAGS)
pois.obj: FFLAGS += $(FOPENMPFLAGS)
pres.obj: FFLAGS += $(FOPENMPFLAGS)
turb.obj: FFLAGS += $(FOPENMPFLAGS)
velo.obj: FFLAGS += $(FOPENMPFLAGS)
rcal.obj: FFLAGS += $(FOPENMPFLAGS)
radi.obj: FFLAGS += $(FOPENMPFLAGS)
scrc.obj: FFLAGS += $(FOPENMPFLAGS)
main.obj: FFLAGS += $(FOPENMPFLAGS)

# Windows Intel Fortran Compiler and Intel MPI

#*** Note: the script, Build\Scripts\setup_intel_compilers.bat, must be run before running make with this target

impi_intel_win_64 : MPILIB = "$(IMPI_RELEASE_ROOT)\impi.lib"
impi_intel_win_64 : MPIINCLUDE = "$(IMPI_INCLUDE)"
impi_intel_win_64 : FFLAGS = /Qipo /O2 /I$(MPIINCLUDE) /wrap-margin- $(GITINFO) $(INTELMPI_COMPINFO) $(FFLAGSMKL_INTEL)
impi_intel_win_64 : FOPENMPFLAGS = /Qopenmp
impi_intel_win_64 : FCOMPL = ifort
impi_intel_win_64 : obj = fds_impi_win_64
impi_intel_win_64 : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(FOPENMPFLAGS) /F1000000000 $(objwin_mpi) $(MPILIB) $(LFLAGSMKL_INTEL)

impi_intel_win_64_dv : MPILIB = "$(IMPI_RELEASE_ROOT)\impi.lib"
impi_intel_win_64_dv : MPIINCLUDE = "$(IMPI_INCLUDE)"
impi_intel_win_64_dv : FFLAGS = /O1 /I$(MPIINCLUDE) /traceback /debug:full /wrap-margin- $(GITINFO) $(INTELMPI_COMPINFO) $(FFLAGSMKL_INTEL)
impi_intel_win_64_dv : FOPENMPFLAGS = /Qopenmp
impi_intel_win_64_dv : FCOMPL = ifort
impi_intel_win_64_dv : obj = fds_impi_win_64_dv
impi_intel_win_64_dv : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(FOPENMPFLAGS) /F1000000000 $(objwin_mpi) $(MPILIB) $(LFLAGSMKL_INTEL)

impi_intel_win_64_db : MPILIB = "$(IMPI_DEBUG_ROOT)\impi.lib"
impi_intel_win_64_db : MPIINCLUDE = "$(IMPI_INCLUDE)"
impi_intel_win_64_db : FFLAGS = /Od /nologo /debug:all /I$(MPIINCLUDE) /Z7 /extend-source:132 /warn:unused /warn:nointerfaces /Qtrapuv /fp:strict /fp:except /traceback /check:all /stand:f18 /wrap-margin- $(GITINFO) $(INTELMPI_COMPINFO) $(FFLAGSMKL_INTEL)
impi_intel_win_64_db : FOPENMPFLAGS = /Qopenmp
impi_intel_win_64_db : FCOMPL = ifort
impi_intel_win_64_db : obj = fds_impi_win_64_db
impi_intel_win_64_db : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS) $(FOPENMPFLAGS) /F1000000000 $(objwin_mpi) $(MPILIB) $(LFLAGSMKL_INTEL)

# Linux Intel Fortran Compiler and Intel MPI
 
impi_intel_linux_64 : FFLAGS = -m64 -O2 -ipo -no-wrap-margin $(GITINFO) $(INTELMPI_COMPINFO) $(FFLAGSMKL_INTEL) -DUSE_IFPORT
impi_intel_linux_64 : LFLAGSMKL = $(LFLAGSMKL_INTEL)
impi_intel_linux_64 : LFLAGS =
impi_intel_linux_64 : FCOMPL = mpiifort
impi_intel_linux_64 : FOPENMPFLAGS = -qopenmp
impi_intel_linux_64 : obj = fds_impi_intel_linux_64
impi_intel_linux_64 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

impi_intel_linux_64_db : FFLAGS = -m64 -check all -warn all -O0 -g -traceback -fpe0 -nofltconsistency -stand:f18 -no-wrap-margin $(GITINFO) $(INTELMPI_COMPINFO) $(FFLAGSMKL_INTEL) -DUSE_IFPORT
impi_intel_linux_64_db : LFLAGSMKL = $(LFLAGSMKL_INTEL)
impi_intel_linux_64_db : LFLAGS =
impi_intel_linux_64_db : FCOMPL = mpiifort
impi_intel_linux_64_db : FOPENMPFLAGS = -qopenmp
impi_intel_linux_64_db : obj = fds_impi_intel_linux_64_db
impi_intel_linux_64_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

impi_intel_linux_64_dv : FFLAGS = -m64 -warn unused -O1 -g -traceback -no-wrap-margin $(GITINFO) $(INTELMPI_COMPINFO) $(FFLAGSMKL_INTEL) -DUSE_IFPORT
impi_intel_linux_64_dv : LFLAGSMKL = $(LFLAGSMKL_INTEL)
impi_intel_linux_64_dv : LFLAGS =
impi_intel_linux_64_dv : FCOMPL = mpiifort
impi_intel_linux_64_dv : FOPENMPFLAGS = -qopenmp
impi_intel_linux_64_dv : obj = fds_impi_intel_linux_64_dv
impi_intel_linux_64_dv : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

# Linux Intel Fortran Compiler and Open MPI
 
mpi_intel_linux_64 : FFLAGS = -m64 -O2 -ipo -no-wrap-margin $(GITINFO) $(OPENMPI_COMPINFO) $(FFLAGSMKL_OPENMPI)
mpi_intel_linux_64 : LFLAGSMKL = $(LFLAGSMKL_OPENMPI)
mpi_intel_linux_64 : LFLAGS =
mpi_intel_linux_64 : FCOMPL = mpifort
mpi_intel_linux_64 : FOPENMPFLAGS = -qopenmp
mpi_intel_linux_64 : obj = fds_mpi_intel_linux_64
mpi_intel_linux_64 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

mpi_intel_linux_64_dv : FFLAGS = -m64 -warn unused -O1 -g -traceback -no-wrap-margin $(GITINFO) $(OPENMPI_COMPINFO) $(FFLAGSMKL_OPENMPI)
mpi_intel_linux_64_dv : LFLAGSMKL = $(LFLAGSMKL_OPENMPI)
mpi_intel_linux_64_dv : LFLAGS =
mpi_intel_linux_64_dv : FCOMPL = mpifort
mpi_intel_linux_64_dv : FOPENMPFLAGS = -qopenmp
mpi_intel_linux_64_dv : obj = fds_mpi_intel_linux_64_dv
mpi_intel_linux_64_dv : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

mpi_intel_linux_64_db : FFLAGS = -m64 -check -warn all -O0 -auto -WB -g -traceback -fpe0 -nofltconsistency -no-wrap-margin $(GITINFO) $(OPENMPI_COMPINFO) $(FFLAGSMKL_OPENMPI)
mpi_intel_linux_64_db : LFLAGSMKL = $(LFLAGSMKL_OPENMPI)
mpi_intel_linux_64_db : LFLAGS =
mpi_intel_linux_64_db : FCOMPL = mpifort
mpi_intel_linux_64_db : FOPENMPFLAGS = -qopenmp
mpi_intel_linux_64_db : obj = fds_mpi_intel_linux_64_db
mpi_intel_linux_64_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

# OSX Intel Fortran and Open MPI
 
mpi_intel_osx_64 : FFLAGS = -m64 -O2 -ipo -stand:f18 -no-wrap-margin $(GITINFO) $(OPENMPI_COMPINFO) $(FFLAGSMKL_CUSTOM)
mpi_intel_osx_64 : LFLAGSMKL = $(LFLAGSMKL_CUSTOM)
mpi_intel_osx_64 : LFLAGS = -static-intel
mpi_intel_osx_64 : FOPENMPFLAGS = -qopenmp -qopenmp-link static
mpi_intel_osx_64 : FCOMPL = mpifort
mpi_intel_osx_64 : obj = fds_mpi_intel_osx_64
mpi_intel_osx_64 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

mpi_intel_osx_64_db : FFLAGS = -m64 -check all -fp-stack-check -warn unused -O0 -auto -WB -g -traceback -fpe0 -nofltconsistency -stand:f18 -no-wrap-margin $(GITINFO) $(OPENMPI_COMPINFO) $(FFLAGSMKL_CUSTOM)
mpi_intel_osx_64_db : LFLAGSMKL = $(LFLAGSMKL_CUSTOM)
mpi_intel_osx_64_db : LFLAGS = -static-intel
mpi_intel_osx_64_db : FOPENMPFLAGS = -qopenmp -qopenmp-link static
mpi_intel_osx_64_db : FCOMPL = mpifort
mpi_intel_osx_64_db : obj = fds_mpi_intel_osx_64_db
mpi_intel_osx_64_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

mpi_intel_osx_64_dv : FFLAGS = -m64 -warn unused -O1 -g -traceback -stand:f18 -no-wrap-margin  $(GITINFO) $(OPENMPI_COMPINFO) $(FFLAGSMKL_CUSTOM)
mpi_intel_osx_64_dv : LFLAGSMKL = $(LFLAGSMKL_CUSTOM)
mpi_intel_osx_64_dv : LFLAGS = -static-intel
mpi_intel_osx_64_dv : FOPENMPFLAGS = -qopenmp -qopenmp-link static
mpi_intel_osx_64_dv : FCOMPL = mpifort
mpi_intel_osx_64_dv : obj = fds_mpi_intel_osx_64_dv
mpi_intel_osx_64_dv : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

lahey_linux_mpi : FFLAGS = -O --tpp --staticlink --nap --nchk --npca --nsav --ntrace
lahey_linux_mpi : FCOMPL = hf77
lahey_linux_mpi : INCLUDEDIR = /usr/local/include
lahey_linux_mpi : obj = fds_mpi_lahey
lahey_linux_mpi : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

f95_solaris10_mpi : INC_FLAGS = -I$(OPT)/SUNWhpc/HPC7.0/include/amd64
f95_solaris10_mpi : LIB_FLAGS = -Bstatic -lfui -Bdynamic -L$(OPT)/SUNWspro/prod/lib/amd64 -L$(OPT)/SUNWhpc/HPC7.0/lib/amd64
f95_solaris10_mpi : FFLAGS = -m64 $(INC_FLAGS)
f95_solaris10_mpi : FCOMPL = mpif90
f95_solaris10_mpi : obj = fds_mpi_solaris
f95_solaris10_mpi : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

f95_solaris11_mpi : INC_FLAGS=-I/opt/SUNWhpc/HPC8.2/sun/include/amd64
f95_solaris11_mpi : LIB_FLAGS=-Bstatic -L/usr/lib/amd64/ld.so.1 -L/opt/sunstudio12.1/lib/amd64 -L/opt/SUNWhpc/HPC8.2/sun/lib/amd64
f95_solaris11_mpi : FFLAGS = -xarch=sse2a -m64 -xchip=native $(INC_FLAGS)
f95_solaris11_mpi : FCOMPL = mpif90
f95_solaris11_mpi : obj = fds_mpi_solaris11
f95_solaris11_mpi : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

#*** Portland Group Compiler***

# AMD Opteron 64Bit PC
mpi_pgi_linux_64 : FFLAGS = -fast -tp k8-64 -Mstandard $(GITINFOPGI)
mpi_pgi_linux_64 : FCOMPL = mpifort
mpi_pgi_linux_64 : LFLAGS  = -fast -tp k8-64 -byteswapio -Bstatic -lc -lgcc_eh
mpi_pgi_linux_64 : obj = fds_mpi_pgi_Opteron64
mpi_pgi_linux_64 : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

#AIX, JUMP, MPI-Version
mpi_xlf_aix : FFLAGS  = -O2 -q64 -qtune=auto -qarch=auto $(GITINFOXLF)
mpi_xlf_aix : FCOMPL  = mpxlf90
mpi_xlf_aix : obj     = fds_jump_mpi_64
mpi_xlf_aix : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)


#*** GNU Compilers ***

mpi_gnu_linux_64 : FFLAGS = -m64 -O2 -std=f2018 -frecursive -ffpe-summary=none -fall-intrinsics $(GITINFOGNU) $(GNU_COMPINFO) $(FFLAGSMKL_GNU_OPENMPI) $(GFORTRAN_OPTIONS)
mpi_gnu_linux_64 : LFLAGSMKL = $(LFLAGSMKL_GNU_OPENMPI)
mpi_gnu_linux_64 : FCOMPL = mpifort
mpi_gnu_linux_64 : FOPENMPFLAGS = -fopenmp
mpi_gnu_linux_64 : obj = fds_mpi_gnu_linux_64
mpi_gnu_linux_64 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

mpi_gnu_linux_64_db : FFLAGS = -m64 -O0 -std=f2018 -ggdb -Wall -Wcharacter-truncation -Wno-target-lifetime -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -frecursive -ffpe-summary=none -fall-intrinsics $(GITINFOGNU) $(GNU_COMPINFO) $(FFLAGSMKL_GNU_OPENMPI) $(GFORTRAN_OPTIONS)
mpi_gnu_linux_64_db : LFLAGSMKL = $(LFLAGSMKL_GNU_OPENMPI)
mpi_gnu_linux_64_db : FCOMPL = mpifort
mpi_gnu_linux_64_db : FOPENMPFLAGS = -fopenmp
mpi_gnu_linux_64_db : obj = fds_mpi_gnu_linux_64_db
mpi_gnu_linux_64_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

mpi_gnu_osx_64 : FFLAGS = -m64 -O2 -std=f2018 -frecursive -ffpe-summary=none -fall-intrinsics $(GITINFOGNU) $(GNU_COMPINFO) $(FFLAGSMKL_GNU_CUSTOM) $(GFORTRAN_OPTIONS)
mpi_gnu_osx_64 : LFLAGSMKL = $(LFLAGSMKL_GNU_CUSTOM)
mpi_gnu_osx_64 : FCOMPL = mpifort
mpi_gnu_osx_64 : FOPENMPFLAGS =
mpi_gnu_osx_64 : obj = fds_mpi_gnu_osx_64
mpi_gnu_osx_64 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

mpi_gnu_osx_64_db : FFLAGS = -m64 -O0 -std=f2018 -ggdb -Wall -Wcharacter-truncation -Wno-target-lifetime -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -frecursive -ffpe-summary=none -fall-intrinsics $(GITINFOGNU) $(GNU_COMPINFO) $(FFLAGSMKL_GNU_CUSTOM) $(GFORTRAN_OPTIONS)
mpi_gnu_osx_64_db : LFLAGSMKL = $(LFLAGSMKL_GNU_CUSTOM)
mpi_gnu_osx_64_db : FCOMPL = mpifort
mpi_gnu_osx_64_db : FOPENMPFLAGS =
mpi_gnu_osx_64_db : obj = fds_mpi_gnu_osx_64_db
mpi_gnu_osx_64_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

mpi_gnu_osx_64_dv : FFLAGS = -m64 -O1 -fbacktrace -std=f2018 -frecursive -ffpe-summary=none -fall-intrinsics $(GITINFOGNU) $(GNU_COMPINFO) $(FFLAGSMKL_GNU_CUSTOM) $(GFORTRAN_OPTIONS)
mpi_gnu_osx_64_dv : LFLAGSMKL = $(LFLAGSMKL_GNU_CUSTOM)
mpi_gnu_osx_64_dv : FCOMPL = mpifort
mpi_gnu_osx_64_dv : FOPENMPFLAGS =
mpi_gnu_osx_64_dv : obj = fds_mpi_gnu_osx_64_dv
mpi_gnu_osx_64_dv : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) $(FOPENMPFLAGS) -o $(obj) $(obj_mpi) $(LFLAGSMKL)

# *** Object Dependencies ***

prec.o :
imkl.o :
cons.o : prec.o
devc.o : prec.o
type.o : prec.o cons.o imkl.o
data.o : prec.o func.o cons.o
rcal.o : func.o prec.o cons.o
mesh.o : prec.o type.o
func.o : prec.o cons.o type.o mesh.o devc.o
gsmv.o : prec.o func.o cons.o mesh.o
smvv.o : geom.o prec.o func.o
turb.o : prec.o func.o cons.o mesh.o
soot.o : prec.o mesh.o cons.o func.o turb.o
pois.o : prec.o cons.o
geom.o : gsmv.o func.o prec.o cons.o type.o mesh.o turb.o imkl.o
ccib.o : gsmv.o func.o prec.o cons.o type.o mesh.o turb.o geom.o
radi.o : func.o prec.o cons.o type.o mesh.o rcal.o geom.o
part.o : func.o prec.o cons.o type.o mesh.o devc.o data.o geom.o soot.o
vege.o : func.o prec.o cons.o type.o mesh.o geom.o
ctrl.o : prec.o cons.o type.o mesh.o func.o
hvac.o : func.o prec.o cons.o type.o mesh.o ctrl.o
mass.o : func.o prec.o cons.o type.o mesh.o turb.o soot.o geom.o ccib.o
wall.o : func.o prec.o cons.o type.o mesh.o hvac.o mass.o soot.o geom.o ccib.o
fire.o : func.o prec.o cons.o type.o mesh.o mass.o soot.o radi.o
velo.o : func.o prec.o cons.o type.o mesh.o turb.o geom.o ccib.o
scrc.o : prec.o cons.o mesh.o func.o pois.o geom.o velo.o imkl.o
pres.o : func.o prec.o cons.o type.o mesh.o pois.o velo.o scrc.o geom.o ccib.o imkl.o
init.o : func.o prec.o cons.o type.o mesh.o devc.o radi.o pois.o data.o pres.o geom.o
dump.o : func.o prec.o cons.o type.o mesh.o devc.o smvv.o turb.o scrc.o data.o radi.o soot.o geom.o ccib.o
read.o : func.o prec.o cons.o type.o mesh.o devc.o ctrl.o radi.o hvac.o scrc.o data.o geom.o mass.o soot.o
divg.o : func.o prec.o cons.o type.o mesh.o scrc.o mass.o ccib.o
main.o : func.o prec.o cons.o type.o mesh.o devc.o smvv.o mass.o divg.o velo.o wall.o fire.o radi.o part.o vege.o dump.o read.o init.o pres.o pois.o ctrl.o turb.o hvac.o scrc.o data.o geom.o ccib.o soot.o

prec.obj :
imkl.obj :
cons.obj : prec.obj
devc.obj : prec.obj
type.obj : prec.obj cons.obj imkl.obj
data.obj : prec.obj func.obj cons.obj
rcal.obj : func.obj prec.obj cons.obj
mesh.obj : prec.obj type.obj
func.obj : prec.obj cons.obj type.obj mesh.obj devc.obj
gsmv.obj : prec.obj func.obj cons.obj mesh.obj
smvv.obj : geom.obj prec.obj func.obj
turb.obj : prec.obj func.obj cons.obj mesh.obj
soot.obj : prec.obj mesh.obj cons.obj func.obj turb.obj
pois.obj : prec.obj cons.obj
geom.obj : gsmv.obj func.obj prec.obj cons.obj type.obj mesh.obj turb.obj imkl.obj
ccib.obj : gsmv.obj func.obj prec.obj cons.obj type.obj mesh.obj turb.obj geom.obj
radi.obj : func.obj prec.obj cons.obj type.obj mesh.obj rcal.obj geom.obj
part.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj data.obj geom.obj soot.obj
vege.obj : func.obj prec.obj cons.obj type.obj mesh.obj geom.obj
ctrl.obj : prec.obj cons.obj type.obj mesh.obj func.obj
hvac.obj : func.obj prec.obj cons.obj type.obj mesh.obj ctrl.obj
mass.obj : func.obj prec.obj cons.obj type.obj mesh.obj turb.obj soot.obj geom.obj ccib.obj
wall.obj : func.obj prec.obj cons.obj type.obj mesh.obj hvac.obj mass.obj soot.obj geom.obj ccib.obj
fire.obj : func.obj prec.obj cons.obj type.obj mesh.obj mass.obj soot.obj radi.obj
velo.obj : func.obj prec.obj cons.obj type.obj mesh.obj turb.obj geom.obj ccib.obj 
scrc.obj : prec.obj cons.obj mesh.obj func.obj pois.obj geom.obj velo.obj imkl.obj
pres.obj : func.obj prec.obj cons.obj type.obj mesh.obj pois.obj velo.obj scrc.obj geom.obj ccib.obj imkl.obj
init.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj pois.obj data.obj pres.obj geom.obj
dump.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj smvv.obj turb.obj scrc.obj data.obj radi.obj soot.obj geom.obj ccib.obj
read.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj ctrl.obj radi.obj hvac.obj scrc.obj data.obj geom.obj mass.obj soot.obj
divg.obj : func.obj prec.obj cons.obj type.obj mesh.obj scrc.obj mass.obj ccib.obj
main.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj smvv.obj mass.obj divg.obj velo.obj wall.obj fire.obj radi.obj part.obj vege.obj dump.obj read.obj init.obj pres.obj pois.obj ctrl.obj turb.obj hvac.obj scrc.obj data.obj geom.obj ccib.obj soot.obj

#*** Clean Target to remove Object and Module files ***

.PHONY : clean
clean:
	-rm -f *.o *.mod *.obj

.PHONY : winclean
winclean:
	-erase *.mod *.obj
