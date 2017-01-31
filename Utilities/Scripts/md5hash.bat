@echo off
certutil -hashfile %1 MD5|find /v "hash of file"|find /v "CertUtil" > hash.out
set /p hash=<hash.out
set hash=%hash: =%
echo %hash%   %1
erase hash.out
