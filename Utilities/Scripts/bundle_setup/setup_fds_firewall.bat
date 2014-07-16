:: This script deletes previous Intel MPI firewall exceptions
:: and creates 2 new rules opening ports equivalent to PORT_VAL.
:: It then sets env variable MPICH_PORT_RANGE to match PORT_VAL.


@echo off
echo in setup_fds_fireall
pause

set PORT_VAL=8670-8690
netsh advfirewall firewall delete rule name="Intel MPI Port for FDS"
netsh advfirewall firewall add rule dir=in action=allow name="Intel MPI Port for FDS" profile=domain protocol=tcp localport=%PORT_VAL%
netsh advfirewall firewall add rule dir=in action=allow name="Intel MPI Port for FDS" profile=domain protocol=udp localport=%PORT_VAL%

setx MPIEXEC_PORT_RANGE 8670:8690
setx MPICH_PORT_RANGE 8670:8690

smpd -install

echo FDS directory=%CD%\bin

echo press any key to continue
pause>Nul