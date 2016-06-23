@echo off

Title Building FDS for 64 bit Windows

make SHELL="%ComSpec%" VPATH="../../Source" -f ..\makefile gnu_win_64_db
pause

