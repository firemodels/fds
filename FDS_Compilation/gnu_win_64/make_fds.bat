@echo off

Title Building FDS for 64 bit Windows

make VPATH="../../FDS_Source" -f ..\makefile gnu_win_64
pause

