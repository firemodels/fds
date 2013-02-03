@echo off

erase *.o
make COMPILER=gcc RM=erase libjpeg.a