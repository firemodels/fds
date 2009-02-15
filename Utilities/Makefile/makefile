# Makefile for Fire Dynamics Simulator (FDS), Version 5

# To use this makefile, select the appropriate compiler and OS from the list below. For example, to
# compile the 64 bit parallel version (MPI) of FDS under Linux using the Intel compilers type:
#
# make intel_linux_mpi_64

# The variable VPATH defines where the source code is relative to the current directory

VPATH = ../../FDS_Source

# Definition of the Serial or MPI object variables (same, except for main vs main_mpi)

obj_serial = isob.o prec.o smvv.o cons.o devc.o type.o mesh.o func.o irad.o \
             ieva.o pois.o radi.o evac.o part.o vege.o ctrl.o dump.o read.o turb.o mass.o \
             wall.o fire.o pres.o divg.o init.o velo.o main.o
objwin_serial = $(obj_serial:.o=.obj)

obj_mpi = isob.o prec.o smvv.o cons.o devc.o type.o mesh.o func.o irad.o \
          ieva.o pois.o radi.o evac.o part.o vege.o ctrl.o dump.o read.o turb.o mass.o \
          wall.o fire.o pres.o divg.o init.o velo.o main_mpi.o
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

intel_win_32 : FFLAGS = -O2 /Qvec_report0
intel_win_32 : CFLAGS = -O2 -Dpp_noappend /Qvec_report0
intel_win_32 : FCOMPL = ifort
intel_win_32 : CCOMPL = icl
intel_win_32 : obj = fds5_win_32
intel_win_32 : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F100000000 $(objwin_serial)

mpi_intel_win_32 : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
mpi_intel_win_32 : MPIINCLUDE = "c:\program files\mpich2\include"
mpi_intel_win_32 : FFLAGS = -O2 /Qvec_report0 /I $(MPIINCLUDE)
mpi_intel_win_32 : CFLAGS = -O2 -Dpp_noappend /Qvec_report0
mpi_intel_win_32 : FCOMPL = ifort
mpi_intel_win_32 : CCOMPL = icl
mpi_intel_win_32 : obj = fds5_win_mpi_32
mpi_intel_win_32 : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F100000000 $(objwin_mpi) $(MPILIB)

 #*** note: the following scripts (located in ..\Scripts directory) must be run before running make with this target
 #    ..\Scripts\iclvars intel64
 #    ..\Scripts\ifortvars intel64

intel_win_64 : FFLAGS = -O2 /Qvec_report0
intel_win_64 : CFLAGS = -O2 -Dpp_noappend /Qvec_report0
intel_win_64 : FCOMPL = ifort
intel_win_64 : CCOMPL = icl
intel_win_64 : obj = fds5_win_64
intel_win_64 : setup_win $(objwin_serial)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F100000000 $(objwin_serial)

mpi_intel_win_64 : MPILIB = "c:\program files\mpich2\lib\fmpich2.lib"
mpi_intel_win_64 : MPIINCLUDE = "c:\program files\mpich2\include"
mpi_intel_win_64 : FFLAGS = -O2 /Qvec_report0 /I $(MPIINCLUDE)
mpi_intel_win_64 : CFLAGS = -O2 -Dpp_noappend /Qvec_report0
mpi_intel_win_64 : FCOMPL = ifort
mpi_intel_win_64 : CCOMPL = icl
mpi_intel_win_64 : obj = fds5_win_mpi_64
mpi_intel_win_64 : setup_win $(objwin_mpi)
	$(FCOMPL) -o $(obj) $(FFLAGS)  /F100000000 $(objwin_mpi) $(MPILIB)

# Intel Compiler, version 11

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

mpi_intel_linux_32 : FFLAGS = -m32 -O3
mpi_intel_linux_32 : CFLAGS = -m32 -O -Dpp_noappend
mpi_intel_linux_32 : FCOMPL = mpif77
mpi_intel_linux_32 : CCOMPL = icc
mpi_intel_linux_32 : obj = fds5_mpi_intel_linux_32
mpi_intel_linux_32 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

mpi_intel_linux_64 : FFLAGS = -m64 -O3
mpi_intel_linux_64 : CFLAGS = -m64 -O -Dpp_noappend
mpi_intel_linux_64 : FCOMPL = mpif77
mpi_intel_linux_64 : CCOMPL = icc
mpi_intel_linux_64 : obj = fds5_mpi_intel_linux_64
mpi_intel_linux_64 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)
	
intel_linux_32_db : FFLAGS = -m32 -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -static
intel_linux_32_db : CFLAGS = -m32 -Dpp_noappend
intel_linux_32_db : FCOMPL = ifort
intel_linux_32_db : CCOMPL = icc
intel_linux_32_db : obj = fds5_intel_linux_32_db
intel_linux_32_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)
	
intel_linux_64_db : FFLAGS = -m64 -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -static
intel_linux_64_db : CFLAGS = -m64 -Dpp_noappend
intel_linux_64_db : FCOMPL = ifort
intel_linux_64_db : CCOMPL = icc
intel_linux_64_db : obj = fds5_intel_linux_64_db
intel_linux_64_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

mpi_intel_linux_32_db : FFLAGS = -m32 -O0 -auto -WB -traceback -g -fpe0 -fltconsistency 
mpi_intel_linux_32_db : CFLAGS = -m32 -Dpp_noappend
mpi_intel_linux_32_db : FCOMPL = mpif77
mpi_intel_linux_32_db : CCOMPL = icc
mpi_intel_linux_32_db : obj = fds5_mpi_intel_linux_32_db
mpi_intel_linux_32_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

mpi_intel_linux_64_db : FFLAGS = -m64 -O0 -auto -WB -traceback -g -fpe0 -fltconsistency
mpi_intel_linux_64_db : CFLAGS = -m64 -O -Dpp_noappend
mpi_intel_linux_64_db : FCOMPL = mpif77
mpi_intel_linux_64_db : CCOMPL = icc
mpi_intel_linux_64_db : obj = fds5_mpi_intel_linux_64_db
mpi_intel_linux_64_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)


# OS X
osx32_libdir = /opt/intel/Compiler/11.0/056/lib
osx64_libdir = /opt/intel/Compiler/11.0/056/lib

# 32 bit
intel_osx_32 : FFLAGS = -O3  -m32 -axSSSE3 -heap-arrays -static-intel -L$(osx32_libdir) -vec_report0
intel_osx_32 : CFLAGS = -O3 -m32 -Dpp_noappend -Dpp_OSX -vec_report0
intel_osx_32 : FCOMPL = ifort
intel_osx_32 : CCOMPL = icc
intel_osx_32 : obj = fds5_osx_32
intel_osx_32 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_osx_mpi_32 : FFLAGS  = -O3 -m32 -axSSSE3 -heap-arrays -static-intel -L$(osx32_libdir)
intel_osx_mpi_32 : CFLAGS  = -O3 -m32 -Dpp_noappend -Dpp_OSX
intel_osx_mpi_32 : FCOMPL  = /usr/local/bin/mpif77
intel_osx_mpi_32 : CCOMPL  = icc
intel_osx_mpi_32 : obj = fds5_mpi_osx_32
intel_osx_mpi_32 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

intel_osx_32_db : FFLAGS = -O0 -m32 -g -traceback -fpe:0 -fltconsistency -FR -WB
intel_osx_32_db : CFLAGS = -O0 -m32 -g -traceback -Dpp_noappend -Dpp_OSX
intel_osx_32_db : FCOMPL = ifort
intel_osx_32_db : CCOMPL = icc
intel_osx_32_db : obj = fds5_osx_32_db
intel_osx_32_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_osx_mpi_32_db : FFLAGS  = -O0 -m32 -g -traceback -fpe:0
intel_osx_mpi_32_db : CFLAGS  = -O0 -m32 -g -traceback -Dpp_noappend -Dpp_OSX
intel_osx_mpi_32_db : FCOMPL  = /usr/local/bin/mpif77
intel_osx_mpi_32_db : CCOMPL  = icc
intel_osx_mpi_32_db : obj = fds5_mpi_osx_32_db
intel_osx_mpi_32_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

# 64 bit
intel_osx_64 : FFLAGS = -O3 -m64 -heap-arrays -axSSSE3 -static-intel -L$(osx64_libdir)
intel_osx_64 : CFLAGS = -O3 -m64 -Dpp_noappend -Dpp_OSX
intel_osx_64 : FCOMPL = ifort
intel_osx_64 : CCOMPL = icc
intel_osx_64 : obj = fds5_osx_64
intel_osx_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_osx_mpi_64 : FFLAGS  = -O3 -m64 -heap-arrays -axSSSE3 -static-intel -L$(osx64_libdir)
intel_osx_mpi_64 : CFLAGS  = -O3 -m64 -Dpp_noappend -Dpp_OSX
intel_osx_mpi_64 : FCOMPL  = /usr/local/bin/mpif77
intel_osx_mpi_64 : CCOMPL  = icc
intel_osx_mpi_64 : obj = fds5_mpi_osx_64
intel_osx_mpi_64 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

intel_osx_64_db : FFLAGS = -O0 -m64 -g -traceback -save-temps -fpe:0 -fltconsistency -fp-stack-check -FR -WB
intel_osx_64_db : CFLAGS = -O0 -g -traceback -m64 -Dpp_noappend -Dpp_OSX
intel_osx_64_db : FCOMPL = ifort
intel_osx_64_db : CCOMPL = icc
intel_osx_64_db : obj = fds5_osx_64_db
intel_osx_64_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_osx_mpi_64_db : FFLAGS  = -O0 -g -traceback -m64 -fpe:0 -fp-stack-check -check all -save-temps
intel_osx_mpi_64_db : CFLAGS  = -O0 -g -traceback -m64 -Dpp_noappend -Dpp_OSX
intel_osx_mpi_64_db : FCOMPL  = /usr/local/bin/mpif77
intel_osx_mpi_64_db : CCOMPL  = icc
intel_osx_mpi_64_db : obj = fds5_mpi_osx_64_db
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

#*** Portland Group Compiler***

# AMD Opteron 64Bit PC
pgi_linux : FFLAGS = -fast -tp k8-64 -Mstandard
pgi_linux : CFLAGS = -fast -tp k8-64
pgi_linux : FCOMPL = pgf90
pgi_linux : CCOMPL = pgcc
pgi_linux : LFLAGS  = -fast -tp k8-64 -byteswapio -Bstatic -lc -lgcc_eh
pgi_linux : obj = fds5_pgi_Opteron64
pgi_linux : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

# AMD Opteron 64Bit PC
pgi_linux_mpi : FFLAGS = -fast -tp k8-64 -Mstandard
pgi_linux_mpi : CFLAGS = -fast -tp k8-64
pgi_linux_mpi : FCOMPL = mpif77
pgi_linux_mpi : CCOMPL = mpicc
pgi_linux_mpi : LFLAGS  = -fast -tp k8-64 -byteswapio -Bstatic -lc -lgcc_eh
pgi_linux_mpi : obj = fds5_mpi_pgi_Opteron64
pgi_linux_mpi : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

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
gnu_32 : FFLAGS  = -O1
gnu_32 : CFLAGS  = -O1 -D pp_noappend
gnu_32 : FCOMPL  = gfortran
gnu_32 : CCOMPL  = gcc
gnu_32 : obj     = fds5_gnu_32
gnu_32 : $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

gnu_mpi_32 : FFLAGS  = -O1
gnu_mpi_32 : CFLAGS  = -O1 -D pp_noappend
gnu_mpi_32 : FCOMPL  = mpif90
gnu_mpi_32 : CCOMPL  = mpicc
gnu_mpi_32 : obj = fds5_gnu_mpi_32
gnu_mpi_32 : $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

gnu_linux_mpi_64 : FFLAGS = -O3 -unroll -static -frecord-marker=4 
gnu_linux_mpi_64 : CFLAGS = -O3 -unroll -static 
gnu_linux_mpi_64 : FCOMPL = mpif90 
gnu_linux_mpi_64 : CCOMPL = mpicc 
gnu_linux_mpi_64 : obj = fds5_mpi_linux_64 
gnu_linux_mpi_64 : setup $(obj_mpi) 
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi) 


#*** End Compiler Specific Rules ***

# *** Object Dependencies ***

isob.o : isob.c
smvv.o : smvv.f90 isob.c
prec.o : prec.f90
cons.o : cons.f90 prec.f90
type.o : type.f90 prec.f90 cons.f90
devc.o : devc.f90 prec.f90
pois.o : pois.f90 prec.f90
mesh.o : mesh.f90 prec.f90 type.f90
func.o : func.f90 prec.f90 cons.f90 type.f90 mesh.f90
turb.o : turb.f90 func.f90 prec.f90 cons.f90 mesh.f90
ctrl.o : ctrl.f90 prec.f90 cons.f90 type.f90 mesh.f90 func.f90
irad.o : irad.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
ieva.o : ieva.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
fire.o : fire.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
wall.o : wall.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
velo.o : velo.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 turb.f90
divg.o : divg.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
mass.o : mass.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 turb.f90
radi.o : radi.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 irad.f90
evac.o : evac.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 ieva.f90
pres.o : pres.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 pois.f90
part.o : part.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90
vege.o : vege.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90
read.o : read.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 irad.f90
init.o : init.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 irad.f90 pois.f90
dump.o : dump.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 smvv.f90 isob.c 
main.o : main.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 smvv.f90 isob.c mass.f90 divg.f90 velo.f90 wall.f90 fire.f90 irad.f90 ieva.f90 radi.f90 evac.f90 part.f90 vege.f90 dump.f90 read.f90 init.f90 pres.f90 pois.f90 ctrl.f90 turb.f90
main_mpi.o : main_mpi.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 smvv.f90 isob.c mass.f90 divg.f90 velo.f90 wall.f90 fire.f90 irad.f90 ieva.f90 radi.f90 evac.f90 part.f90 vege.f90 dump.f90 read.f90 init.f90 pres.f90 pois.f90 ctrl.f90 turb.f90

isob.obj : isob.c
smvv.obj : smvv.f90 isob.c
prec.obj : prec.f90
cons.obj : cons.f90 prec.f90
type.obj : type.f90 prec.f90 cons.f90
devc.obj : devc.f90 prec.f90
pois.obj : pois.f90 prec.f90
mesh.obj : mesh.f90 prec.f90 type.f90
func.obj : func.f90 prec.f90 cons.f90 type.f90 mesh.f90
turb.obj : turb.f90 func.f90 prec.f90 cons.f90 mesh.f90
ctrl.obj : ctrl.f90 prec.f90 cons.f90 type.f90 mesh.f90 func.f90
irad.obj : irad.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
ieva.obj : ieva.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
fire.obj : fire.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
wall.obj : wall.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
velo.obj : velo.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 turb.f90
divg.obj : divg.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
mass.obj : mass.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 turb.f90
radi.obj : radi.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 irad.f90
evac.obj : evac.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 ieva.f90
pres.obj : pres.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 pois.f90
part.obj : part.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90
vege.obj : vege.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90
read.obj : read.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 irad.f90
init.obj : init.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 irad.f90 pois.f90
dump.obj : dump.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 smvv.f90 isob.c 
main.obj : main.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 smvv.f90 isob.c mass.f90 divg.f90 velo.f90 wall.f90 fire.f90 irad.f90 ieva.f90 radi.f90 evac.f90 part.f90 vege.f90 dump.f90 read.f90 init.f90 pres.f90 pois.f90 ctrl.f90 turb.f90
main_mpi.obj : main_mpi.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 smvv.f90 isob.c mass.f90 divg.f90 velo.f90 wall.f90 fire.f90 irad.f90 ieva.f90 radi.f90 evac.f90 part.f90 vege.f90 dump.f90 read.f90 init.f90 pres.f90 pois.f90 ctrl.f90 turb.f90


#*** Clean Target to remove Object and Module files ***

.PHONY : clean
clean:
	-rm -f *.o *.mod *.obj
