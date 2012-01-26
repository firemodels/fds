set convert=..\..\..\utilities\windsodar2fds\intel_win_64\wind2fds_win_64

%convert% -prefix sd11 -offset "   1.0    1.0 184.0" a111024.csv t111024_exp.csv
%convert% -prefix sd21 -offset "7999.0    1.0 254.0" a111026.csv t111026_exp.csv
%convert% -prefix sd22 -offset "7999.0 7999.0 124.0" a111027.csv t111027_exp.csv
%convert% -prefix sd12 -offset "   1.0 7999.0 244.0" a111025.csv t111025_exp.csv
