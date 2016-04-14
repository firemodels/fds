@echo off
set fds=fds
::set fds=..\..\..\..\fds_compilation\mpi_intel_win_64\fds_mpi_win_64

set case=blodget
%fds% %case%.fds

set case=nist
%fds% %case%.fds

set case=sugarloaf
%fds% %case%.fds

set case=trails
%fds% %case%.fds
