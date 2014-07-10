@echo off

IF not EXIST placeholder.txt goto dircheck
echo ***error: This script is running in the wrong directory.
pause
exit
:dircheck

echo.
echo Wrapping up 64 bit Smokeview update
echo.

echo.
echo Associating the .smv file extension with smokeview.exe

ftype smvDoc="%CD%\smokeview.exe" "%%1" >Nul
assoc .smv=smvDoc>Nul

:: adding path

echo.
echo Adding %CD% to the system path
call "%CD%\set_path.exe" -s -m -a "%CD%"

:: add firewall rules

:: netsh advfirewall firewall add rule dir=in action=allow name="IMPI-MPIEXEC.SMPD" program="%CD%\mpiexec.smpd.exe" profile=domain protocol=tcp
:: netsh advfirewall firewall add rule dir=in action=allow name="IMPI-MPIEXEC.SMPD" program="%CD%\mpiexec.smpd.exe" profile=domain protocol=udp

:: netsh advfirewall firewall add rule dir=in action=allow name="IMPI-SMPD" program="%CD%\smpd.exe" profile=domain protocol=tcp
:: netsh advfirewall firewall add rule dir=in action=allow name="IMPI-SMPD" program="%CD%\smpd.exe" profile=domain protocol=udp

erase "%CD%\set_path.exe"
echo Press any key to complete update
pause>NUL
