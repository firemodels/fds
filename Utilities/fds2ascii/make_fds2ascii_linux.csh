#!/bin/csh -f

#call %intelbin%\ifortvars ia32

ifort -o fds2ascii_linux ../Data_Processing/fds2ascii.f90
