@echo off

call ..\setopts %OPTS%
erase *.o
make COMPILER=%COMPILER% SIZE=%SIZE% RM=erase libpng.a