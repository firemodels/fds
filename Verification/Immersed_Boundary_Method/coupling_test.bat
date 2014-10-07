@echo off
set coupler=..\..\Utilities\Data_Processing\coupling_emulator.exe
set fds=..\..\FDS_Compilation\intel_win_64\fds_win_64.exe
%coupler% %fds%
