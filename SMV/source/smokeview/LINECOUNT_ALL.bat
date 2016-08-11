@echo off
cat ..\..\..\FDS_Source\gsmv.f90 *.c *.cpp *.h *.f90 ..\shared\* ..\background\*.c ..\background\*.h ..\smokediff\*.h ..\smokediff\*.c ..\smokezip\*.c ..\smokezip\*.h ..\smokezip\*.f90 ..\wind2fds\*.c ..\wind2fds\*.h| wc
