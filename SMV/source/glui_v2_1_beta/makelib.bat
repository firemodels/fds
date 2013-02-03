@echo off

call ..\setopts %OPTS%
erase *.o
make COMPILER=%COMPILER% COMPILER2=%COMPILER2% SIZE=%SIZE% RM=erase libglui.a