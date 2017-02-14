@echo off
set hashfile=%temp%\hash.out
certutil -hashfile %1 MD5|find /v "hash of file"|find /v "CertUtil" > %hashfile%
set /p hash=<%hashfile%
set hash=%hash: =%
echo %hash%   %1
erase %hashfile%
