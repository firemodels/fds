@echo off
set fds=fds
::set fds=..\..\..\..\FDS_Compilation\intel_win_64_db\fds_win_64_db
set fds=..\..\..\..\FDS_Compilation\mpi_intel_win_64\fds_mpi_win_64

set case=blodget
%fds% %case%.fds

set case=nist
%fds% %case%.fds

set case=sugarloaf
%fds% %case%.fds

set case=trails
%fds% %case%.fds

set case=test1x1
%fds% %case%.fds

set case=test1x2
%fds% %case%.fds

set case=test3x3
%fds% %case%.fds
