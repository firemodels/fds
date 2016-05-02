set intelbin="%IFORT_COMPILER16%\bin"

call %intelbin%\ifortvars intel64

echo ***warning: use script in Utilities\fds2ascii\intel_win_64 to build fds2ascii
ifort -o layer_height layer_height.f
ifort -o coupling_emulator -fpp -D pp_WINDOWS coupling_emulator.f90
pause
