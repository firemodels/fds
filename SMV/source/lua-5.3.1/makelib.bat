@echo off
setlocal
cd src
erase *.o *.lib *.def *.dll
call "%VS_COMPILER%\vcvarsall" x86_amd64
cd ..
make windows
cd src
lib /machine:x64 /def:liblua.def
endlocal
