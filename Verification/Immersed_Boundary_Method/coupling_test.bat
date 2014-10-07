@echo off
set coupler=..\..\utilities\data_processing\coupling_emulator.exe
set fds=..\..\fds_compilation\intel_win_64\fds_win_64.exe
%coupler% %fds% tetra_1.fds
