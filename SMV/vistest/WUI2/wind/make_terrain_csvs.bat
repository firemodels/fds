set convert=..\..\..\utilities\wind2fds\intel_win_64\wind2fds_win_64

%convert% -prefix sd11 -offset " 1000  1000.0 184.0" a111024.csv t111024_exp.csv
%convert% -prefix sd21 -offset "7000.0 1000.0 254.0" a111026.csv t111026_exp.csv
%convert% -prefix sd22 -offset "7000.0 7000.0 124.0" a111027.csv t111027_exp.csv
%convert% -prefix sd12 -offset "1000.0 7000.0 244.0" a111025.csv t111025_exp.csv
%convert% -date 12/25/2011 -wv -prefix wv_1 -offset "4130.542 3924.751 250" site17_1214_0102c_ORIG wv_01.csv
%convert% -date 12/25/2011 -wv -prefix wv_2 -offset "3955.698 3697.627 250" site17_1214_0102c_ORIG wv_02.csv
%convert% -date 12/25/2011 -wv -prefix wv_3 -offset "4542.663 3815.957 250" site17_1214_0102c_ORIG wv_03.csv
%convert% -date 12/25/2011 -wv -prefix wv_4 -offset "4486.409  3974.691 250" site17_1214_0102c_ORIG wv_04.csv
%convert% -date 12/25/2011 -wv -prefix wv_5 -offset "3974.607 3804.412 250" site17_1214_0102c_ORIG wv_05.csv