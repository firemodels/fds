@echo off
set prog=%1
set in=%2
set out=%3

mpiexec -localonly -n 1 %prog% %in% > %out% 2>&1
exit