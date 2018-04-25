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

# ------------- Win 64 ----------------
#  64 bit windows

intel_win_64 : FFLAGS    = -O2 /nologo /fpscomp:general $(GITINFO)
intel_win_64 : FC        = ifort
intel_win_64 : exe       = fds2ascii_win_64.exe

intel_win_64 : $(objwin)
	$(FC) -o $(bin)/$(exe) $(objwin)

# ------------- gnu linux 64 ----------------
#  64 bit gcc linux

gnu_linux_64 : FFLAGS    = -O0 -m64 -ffree-form -frecord-marker=4  $(GITINFOGNU)
gnu_linux_64 : FC        = gfortran
gnu_linux_64 : exe       = fds2ascii_linux_64

gnu_linux_64 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- gnu osx 64 ----------------
#  64 bit gcc linux

gnu_osx_64 : FFLAGS    = -O0 -m64 -ffree-form -frecord-marker=4  $(GITINFOGNU)
gnu_osx_64 : FC        = gfortran
gnu_osx_64 : exe       = fds2ascii_osx_64

gnu_osx_64 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)



# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : FFLAGS    = -O2 -m64 $(GITINFO)
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = fds2ascii_linux_64

intel_linux_64 : $(obj)
	$(FC) -o $(bin)/$(exe) $(obj)

# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : FFLAGS    = -O2 -m64 $(GITINFO)
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = fds2ascii_osx_64

intel_osx_64 : $(obj)
	$(FC) -o $(bin)/$(exe)  -mmacosx-version-min=10.7 $(obj)

#-------------- Clean Target to remove Object and Module files -----------

.PHONY : clean
clean:
	rm *.o
