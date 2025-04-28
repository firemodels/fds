@echo off
set bindir=%1
:: This script deletes previous Intel MPI firewall exceptions and creates two
:: new rules opening ports specified by the environment variable PORT_VAL.

set PORT_VAL=8670-8690
netsh advfirewall firewall delete rule name="Intel MPI Port for FDS" > Nul
netsh advfirewall firewall add rule dir=in action=allow name="Intel MPI Port for FDS" profile=domain,private protocol=tcp localport=%PORT_VAL% > Nul
netsh advfirewall firewall add rule dir=in action=allow name="Intel MPI Port for FDS" profile=domain,private protocol=udp localport=%PORT_VAL% > Nul

:: setx -m MPIEXEC_PORT_RANGE 8670:8690 > Nul
:: setx -m MPICH_PORT_RANGE   8670:8690 > Nul

%bindir%\hydra_service -install > Nul

