@echo off

Title Building FDS for 64 bit Windows

make SHELL="%ComSpec%" VPATH="../../FDS_Source" -f ..\makefile gnu_win_64
pause

