@echo off
set convert=..\..\utilities\wind2fds\intel_win_64\wind2fds_win_64

%convert% -prefix sd11 -offset "100.0 100.0 0.0" wind_test1a.csv
Rem %convert% -prefix sd12 -offset " 50.0 150.0 0.0" wind_test1b.csv
Rem %convert% -prefix sd21 -offset "150.0  50.0 0.0" wind_test1c.csv
Rem %convert% -prefix sd22 -offset "150.0 150.0 0.0" wind_test1d.csv
