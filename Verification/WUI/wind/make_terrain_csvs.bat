set convert=..\..\..\utilities\sodar2fds\intel_win_64\sodar2fds_win_64

%convert% -prefix sd11 -offset "   0.0    0.0 184.0" a111024.csv t111024_exp.csv
%convert% -prefix sd21 -offset "8000.0    0.0 254.0" a111026.csv t111026_exp.csv
%convert% -prefix sd22 -offset "8000.0 8000.0 124.0" a111027.csv t111027_exp.csv
%convert% -prefix sd12 -offset "   0.0 8000.0 244.0" a111025.csv t111025_exp.csv
