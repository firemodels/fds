#!/bin/csh -f
# make tools on the fire linux cluster
icc -o ~/bin/catcol catcol.c
icc -o ~/bin/headskip headskip.c
ifort -o ~/bin/fds2ascii_linux fds2ascii.f90

# make tools on hercules
#pgcc -o ~/bin/catcol catcol.c
#pgcc -o ~/bin/headskip headskip.c
#pgf90 -o ~/bin/fds2ascii_linux fds2ascii.f

#
# 
