# Makefile for geom_test

VPATH = ../Source;../../Source

ifeq ($(shell echo "check_quotes"),"check_quotes")
# windows
  BUILD_DATE := $(shell date /t) $(shell time /t)
  GIT_HASH := $(shell ..\..\Utilities\Scripts\githash)
  GIT_DATE := $(shell ..\..\Utilities\Scripts\gitlog)
else
# linux/osx
  GIT_HASH := $(shell git describe --long --dirty)
  GIT_DATE := $(shell git log -1 --format=%cd)
  BUILD_DATE := $(shell date "+%b %d, %Y  %T")
# The preprocessor used by xlf uses comma as separation between preprocessing commands.
# Commas needs to be escaped if they are to be used as part of a value in a preprocessing variable.
  BUILD_DATE_XLF := $(shell date "+%b %d\, %Y  %T")
endif

GITINFO_BASE=-DGITHASH_PP=\"$(GIT_HASH)\" -DGITDATE_PP=\""$(GIT_DATE)\""
GITINFO=-fpp $(GITINFO_BASE) -DBUILDDATE_PP=\""$(BUILD_DATE)\""
GITINFOGNU=-cpp $(GITINFO_BASE) -DBUILDDATE_PP=\""$(BUILD_DATE)\""
GITINFOXLF=-qsuffix=cpp=f90 '-WF,-DGITHASH_PP=\"$(GIT_HASH)\",-DGITDATE_PP=\"$(GIT_DATE)\",-DBUILDDATE_PP=\"$(BUILD_DATE_XLF)\"'
GITINFOPGI=-Mpreprocess $(GITINFO_BASE) -DBUILDDATE_PP=\""$(BUILD_DATE)\""

# Definition of the non-MPI object variables

obj = prec.o cons.o devc.o type.o mesh.o func.o geom.o gsmv.o geom_main.o

objwin = $(obj:.o=.obj)

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

# Windows Intel Compiler

#*** note: the following scripts (located in ..\Scripts directory) must be run before running make with this target
#    ..\Scripts\ifortvars intel64

intel_win_64 : FFLAGS = /O /traceback $(GITINFO)
intel_win_64 : FCOMPL = ifort
intel_win_64 : obj = fds_win_64
intel_win_64 : setup_win $(objwin)
	$(FCOMPL) -o $(obj) $(FFLAGS) /F1000000000 $(objwin)

# linux

intel_linux_64 : FFLAGS = -m64 -O2 -traceback -no-wrap-margin -assume realloc_lhs $(GITINFO)
intel_linux_64 : LFLAGS = -static-intel
intel_linux_64 : FCOMPL = ifort
intel_linux_64 : obj = fds_intel_linux_64
intel_linux_64 : setup $(obj)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj)

# OSX

intel_osx_64 : FFLAGS = -m64 -no-wrap-margin $(GITINFO)
intel_osx_64 : LFLAGS = -static-intel
intel_osx_64 : FCOMPL = ifort
intel_osx_64 : obj = fds_intel_osx_64
intel_osx_64 : setup $(obj)
	$(FCOMPL) $(FFLAGS) $(LFLAGS) -o $(obj) $(obj)


# *** Object Dependencies ***

cons.o : prec.o
type.o : prec.o cons.o
devc.o : prec.o
mesh.o : prec.o type.o
func.o : prec.o cons.o type.o mesh.o devc.o
gsmv.o : func.o prec.o cons.o mesh.o
geom.o : gsmv.o func.o prec.o cons.o mesh.o
geom_main.o : func.o prec.o cons.o type.o mesh.o devc.o geom.o

cons.obj : prec.obj
type.obj : prec.obj cons.obj
devc.obj : prec.obj
mesh.obj : prec.obj type.obj
func.obj : prec.obj cons.obj type.obj mesh.obj devc.obj
gsmv.obj : func.obj prec.obj cons.obj mesh.obj
geom.obj : func.obj prec.obj cons.obj mesh.obj
geom_main.obj : func.obj prec.obj cons.obj type.obj mesh.obj devc.obj geom.obj

#*** Clean Target to remove Object and Module files ***

.PHONY : clean
clean:
	-rm -f *.o *.mod *.obj

.PHONY : winclean
winclean:
	-erase *.mod *.obj
