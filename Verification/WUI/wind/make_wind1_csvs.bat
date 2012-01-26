set convert=..\..\..\utilities\wind2fds\intel_win_64\wind2fds_win_64

%convert% -prefix sd11 -offset " 50.0  50.0 0.0" a111024.csv
%convert% -prefix sd12 -offset " 50.0 150.0 0.0" a111025.csv
%convert% -prefix sd21 -offset "150.0  50.0 0.0" a111026.csv
%convert% -prefix sd22 -offset "150.0 150.0 0.0" a111027.csv