@echo off

call ..\setopts %OPTS%
erase *.o
make COMPILER=%COMPILER% SIZE=%SIZE% libgd.a -f makefile_win