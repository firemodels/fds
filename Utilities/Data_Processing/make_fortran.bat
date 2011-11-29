set intelbin="%IFORT_COMPILER11%\bin"

call %intelbin%\ifortvars ia32

ifort -o fds2ascii fds2ascii.f90
ifort -o Hamins_CH4_Average Hamins_CH4_Average.f90
ifort -o layer_height layer_height.f
ifort -o NIST_RSE_1994 NIST_RSE_1994.f90
pause
