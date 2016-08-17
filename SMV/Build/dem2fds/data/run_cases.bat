@echo off
set fds=fds
::set fds=..\..\..\..\FDS_Compilation\intel_win_64_db\fds_win_64_db
::set fds=..\..\..\..\FDS_Compilation\mpi_intel_win_64\fds_mpi_win_64

%fds% blodget.fds 

%fds% demtest1.fds 
%fds% demtest2.fds 

%fds% nist.fds 

%fds% sugarloaf.fds 

%fds% tower.fds 

%fds% trails.fds 
%fds% trails2.fds 
