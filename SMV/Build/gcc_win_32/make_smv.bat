@echo off

Rem erase *.o
make -j4 -f ..\Makefile gcc_win_32
pause