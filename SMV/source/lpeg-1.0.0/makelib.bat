@echo off
setlocal
call "%VS_COMPILER%\vcvarsall" x86_amd64
erase *.o *.lib
make lpeg.lib
endlocal
