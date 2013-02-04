@echo off

call ..\setopts %OPTS%
erase *.o
make COMPILER=%COMPILER% SIZE=%SIZE% FILTERC="-D WIN32" RM=erase libglutwin.a