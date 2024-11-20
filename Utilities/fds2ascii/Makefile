SOURCE_DIR = ../..
FILTERC =
FILTERF =

VPATH = ../
bin = .

ifeq ($(shell echo "check_quotes"),"check_quotes")
# windows
  BUILD_DATE := $(shell date /t) $(shell time /t)
  GIT_HASH := $(shell ..\..\Scripts\githash)
  GIT_DATE := $(shell ..\..\Scripts\gitlog)
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



# Definition of the object variables

obj = fds2ascii.o
objwin = fds2ascii.obj

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .F90 .c .obj .o .f90 .cpp

.F90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.obj:
	$(FC) -c $(FFLAGS) $(FILTERF) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FILTERF) $<

intel_win : FFLAGS    = -O2 /nologo /fpscomp:general $(GITINFO)
intel_win : FC        = ifx
intel_win : exe       = fds2ascii_intel_win.exe
intel_win : $(objwin)
	$(FC) -o $(bin)/$(exe) $(objwin)

gnu_linux : FFLAGS    = -O0 -m64 -ffree-form -frecord-marker=4  $(GITINFOGNU)
gnu_linux : FC        = gfortran
gnu_linux : exe       = fds2ascii_gnu_linux
gnu_linux : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

gnu_osx : FFLAGS    = -O0 -m64 -ffree-form -frecord-marker=4  $(GITINFOGNU)
gnu_osx : FC        = gfortran
gnu_osx : exe       = fds2ascii_gnu_osx
gnu_osx : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

intel_linux : FFLAGS    = -O2 -m64 $(GITINFO)
intel_linux : FC        = ifx
intel_linux : exe       = fds2ascii_intel_linux
intel_linux : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

intel_osx : FFLAGS    = -O2 -m64 $(GITINFO)
intel_osx : FC        = ifort
intel_osx : exe       = fds2ascii_intel_osx
intel_osx : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm *.o
