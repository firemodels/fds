#*** Fire Dynamics Simulator (FDS), Version 5, Makefile ***

#To use this makefile, you will need to provide a make target
#at the command line.  For example, to make the MPI version with
#the intel compiler and the debugging flags set, you would need
#to type 'make intel_linux_mpi_db' at the prompt.
#Below this comment section there are make file Variables that 
#can be modified by the user. Please make the changes that 
#are necessary for your system configuration.

#*** Begin User Defined Variables Section ***
VPATH = ../../FDS_5/SOURCE_FDS5

#*** End User Defined Variables Section***

#*** Serial or MPI object variables ***
obj_serial = isob.o prec.o smvv.o cons.o devc.o type.o mesh.o func.o irad.o \
             ieva.o pois.o radi.o evac.o part.o vege.o ctrl.o dump.o read.o mass.o \
             wall.o fire.o pres.o divg.o init.o velo.o main.o

obj_mpi = isob.o prec.o smvv.o cons.o devc.o type.o mesh.o func.o irad.o \
          ieva.o pois.o radi.o evac.o part.o vege.o ctrl.o dump.o read.o mass.o \
          wall.o fire.o pres.o divg.o init.o velo.o main_mpi.o

#*** General Purpose Rules ***
no_target:
	@echo \******** You did not specify a make target \********
	@echo Please read the comments at the top of the makefile

setup:
%.o : %.mod

.SUFFIXES: .c .f90 .o

.f90.o:
	$(FCOMPL) -c $(FFLAGS) $<
.c.o:
	$(CCOMPL) -c $(CFLAGS) $<

#*** Begin Compiler Specific Rules ***

#*** Intel Compiler ***

# Linux
intel_linux : FFLAGS = -O3 -axP -ip -static -vec_report0
intel_linux : CFLAGS = -O -Dpp_noappend
intel_linux : FCOMPL = ifort
intel_linux : CCOMPL = icc
intel_linux : obj = fds5_intel
intel_linux : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)
	
intel_linux_db : $(obj_serial)
intel_linux_db : FFLAGS = -FR  -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -static
intel_linux_db : CFLAGS = -O -Dpp_noappend
intel_linux_db : FCOMPL = ifort
intel_linux_db : CCOMPL = icc
intel_linux_db : obj = fds5_intel_db
intel_linux_db : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)
	
intel_linux_mpi : FFLAGS = -O3 -axP -ip -vec_report0
intel_linux_mpi : CFLAGS = -O -Dpp_noappend
intel_linux_mpi : FCOMPL = mpif77
intel_linux_mpi : CCOMPL = icc
intel_linux_mpi : obj = fds5_mpi_intel
intel_linux_mpi : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

intel_linux_mpi_db : FFLAGS = -FR -check all -O0 -auto -WB -traceback -g -fpe0 -fltconsistency
intel_linux_mpi_db : CFLAGS = -O -Dpp_noappend
intel_linux_mpi_db : FCOMPL = mpif77
intel_linux_mpi_db : CCOMPL = icc
intel_linux_mpi_db : obj = fds5_mpi_intel_db
intel_linux_mpi_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

# OS X
intel_osx_32 : FFLAGS = -O3 -axTP -heap-arrays -ip -vec_report0 -m32
intel_osx_32 : CFLAGS = -O -Dpp_noappend -Dpp_OSX -m32
intel_osx_32 : FCOMPL = ifort
intel_osx_32 : CCOMPL = icc
intel_osx_32 : obj = fds5_osx_intel_32
intel_osx_32 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_osx_mpi_32 : FFLAGS  = -O3 -axTP -ip -vec_report0 -m32
intel_osx_mpi_32 : CFLAGS  = -Dpp_noappend -Dpp_OSX
intel_osx_mpi_32 : FCOMPL  = /usr/local/bin/mpif77
intel_osx_mpi_32 : CCOMPL  = icc
intel_osx_mpi_32 : obj = fds5_mpi_osx_intel_32
intel_osx_mpi_32 : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

intel_osx_db_32 : FFLAGS = -FR -O0 -m32 -heap-arrays -WB -traceback -g -fpe3 -fltconsistency
intel_osx_db_32 : CFLAGS = -m32 -Dpp_noappend -Dpp_OSX
intel_osx_db_32 : FCOMPL = ifort
intel_osx_db_32 : CCOMPL = icc
intel_osx_db_32 : obj = fds5_osx_intel_db_32
intel_osx_db_32 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_osx_mpi_db : FFLAGS  = -O -g
intel_osx_mpi_db : CFLAGS  = -O -g -Dpp_noappend -Dpp_OSX
intel_osx_mpi_db : FCOMPL  = /usr/local/bin/mpif77
intel_osx_mpi_db : CCOMPL  = icc
intel_osx_mpi_db : obj = fds5_osx_intel_mpi_db
intel_osx_mpi_db : setup $(obj_mpi)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_mpi)

intel_osx_64 : FFLAGS = -O3 -axTP -heap-arrays -ip -vec_report0 -m64
intel_osx_64 : CFLAGS = -O2 -Dpp_noappend -Dpp_OSX
intel_osx_64 : FCOMPL = ifort
intel_osx_64 : CCOMPL = icc
intel_osx_64 : obj = fds5_osx_intel_64
intel_osx_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_osx_64_mc : FFLAGS = -O3 -axTP -heap-arrays -ip -parallel -m64
intel_osx_64_mc : CFLAGS = -O2 -Dpp_noappend -Dpp_OSX
intel_osx_64_mc : FCOMPL = ifort
intel_osx_64_mc : CCOMPL = icc
intel_osx_64_mc : obj = fds5_osx_intel_64_mc
intel_osx_64_mc : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_osx_db_64 : FFLAGS = -FR  -O0 -m64 -heap-arrays -WB -traceback -g -fpe:0 -fltconsistency -fpstkchk
intel_osx_db_64 : CFLAGS = -m64 -Dpp_noappend -Dpp_OSX
intel_osx_db_64 : FCOMPL = ifort
intel_osx_db_64 : CCOMPL = icc
intel_osx_db_64 : obj = fds5_osx_intel_db_64
intel_osx_db_64 : setup $(obj_serial)
	$(FCOMPL) $(FFLAGS) -o $(obj) $(obj_serial)

intel_osx_mpi_64 : FFLAGS  = -O3 -axP -ip -parallel -m64
intel_osx_mpi_64 : CFLAGS  = -Dpp_noappend -Dpp_OSX
intel_osx_mpi_64 : FCOMPL  = /usr/local/bin/mpif77
intel_osx_mpi_64 : CCOMPL  = icc
intel_osx_mpi_64 : obj = fds5_mpi_osx_intel_64
intel_osx_mpi_64 : setup $(obj_mpi)
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
ctrl.o : ctrl.f90 prec.f90 cons.f90 type.f90 mesh.f90 func.f90
irad.o : irad.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
ieva.o : ieva.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
fire.o : fire.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
wall.o : wall.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
velo.o : velo.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
divg.o : divg.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
mass.o : mass.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90
radi.o : radi.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 irad.f90
evac.o : evac.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 ieva.f90
pres.o : pres.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 pois.f90
part.o : part.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90
vege.o : vege.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90
read.o : read.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 irad.f90
init.o : init.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 irad.f90 pois.f90
dump.o : dump.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 smvv.f90 isob.c
main.o : main.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 smvv.f90 isob.c mass.f90 divg.f90 velo.f90 wall.f90 fire.f90 irad.f90 ieva.f90 radi.f90 evac.f90 part.f90 vege.f90 dump.f90 read.f90 init.f90 pres.f90 pois.f90 ctrl.f90
main_mpi.o : main_mpi.f90 func.f90 prec.f90 cons.f90 type.f90 mesh.f90 devc.f90 smvv.f90 isob.c mass.f90 divg.f90 velo.f90 wall.f90 fire.f90 irad.f90 ieva.f90 radi.f90 evac.f90 part.f90 vege.f90 dump.f90 read.f90 init.f90 pres.f90 pois.f90 ctrl.f90

#*** Clean Target to remove Object and Module files ***

.PHONY : clean
clean:
	-rm -f *.o *.mod *.obj
