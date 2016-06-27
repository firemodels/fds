@echo off

:: run program that emulates coupling between FDS and code that models geometry changes

set coupler=..\..\Utilities\coupling_emulator\intel_win_64\coupling_emulator_win_64.exe
set fds=..\..\FDS_Compilation\intel_win_64\fds_win_64.exe
%coupler% %fds% W
