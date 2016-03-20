set convert=..\..\..\..\utilities\wind2fds\intel_win_64\wind2fds_win_64

%convert% -prefix sd11 -offset " 50.0  50.0 0.0" a111024.csv
%convert% -prefix sd12 -offset " 50.0 150.0 0.0" a111025.csv
%convert% -prefix sd21 -offset "150.0  50.0 0.0" a111026.csv
%convert% -prefix sd22 -offset "150.0 150.0 0.0" a111027.csv
%convert% -mindate "12/25/2011 15:37:32" -maxdate "1/1/2012 10:01:30" -wv -prefix wv11 -offset "100.0 100.0 50.0" site17_1214_0102c_ORIG wv11.csv