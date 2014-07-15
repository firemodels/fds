@echo off
echo in setup_fds_fireall
pause

:: netsh advfirewall firewall add rule dir=in action=allow name="IMPI-MPIEXEC.SMPD" program="%CD%\mpiexec.smpd.exe" profile=domain protocol=tcp
:: netsh advfirewall firewall add rule dir=in action=allow name="IMPI-MPIEXEC.SMPD" program="%CD%\mpiexec.smpd.exe" profile=domain protocol=udp

:: netsh advfirewall firewall add rule dir=in action=allow name="IMPI-SMPD" program="%CD%\smpd.exe" profile=domain protocol=tcp
:: netsh advfirewall firewall add rule dir=in action=allow name="IMPI-SMPD" program="%CD%\smpd.exe" profile=domain protocol=udp

echo FDS directory=%CD%\bin

echo press any key to continue
pause>Nul
