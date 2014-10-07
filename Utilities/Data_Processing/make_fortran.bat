set intelbin="%IFORT_COMPILER14%\bin"

call %intelbin%\ifortvars intel64

ifort -o fds2ascii fds2ascii.f90
ifort -o layer_height layer_height.f
ifort -o coupling_emulator -fpp -D pp_WINDOWS coupling_emulator.f90
pause
