% Charles Luo
% 4-20-2012
% cutcell_area.m
%

close all
clear all

dir = '../../Verification/Immersed_Boundary_Method/Demo/';
infile{1} = 'tri_cube_cut_cell_test_01_cc.csv';
infile{2} = 'tri_cube_cut_cell_test_02_cc.csv';
infile{3} = 'tri_cube_cut_cell_test_03_cc.csv';
infile{4} = 'tri_cube_cut_cell_test_04_cc.csv';
infile{5} = 'tri_cube_cut_cell_test_05_cc.csv';
infile{6} = 'tri_cube_cut_cell_test_06_cc.csv';
infile{7} = 'tri_cube_cut_cell_test_07_cc.csv';
infile{8} = 'tri_cube_cut_cell_test_08_cc.csv';
infile{9} = 'tri_cube_cut_cell_test_09_cc.csv';
infile{10} = 'tri_cube_cut_cell_test_10_cc.csv';
infile{11} = 'tri_cube_cut_cell_test_11_cc.csv';
infile{12} = 'tri_cube_cut_cell_test_12_cc.csv';
for kcase=1:12
    M = csvread([dir,infile{kcase}],1,0);
    [nrow,ncol] = size(M);
    for i=1:nrow
        icell = M(i,2);
        if (icell == 14)
            cut_area(kcase) = M(i,3);
            break
        end
    end
end

tri_area(1) = 0.433;
tri_area(2) = 0.2165;
tri_area(3) = 0.375;
tri_area(4) = 0.2179;
tri_area(5) = 0.866;
tri_area(6) = 0.559;
tri_area(7) = 0.7071;
tri_area(8) = 0.275;
tri_area(9) = 0.7071;
tri_area(10) = 0.3766;
tri_area(11) = 1.215;
tri_area(12) = 1.4614;

filename = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/cutcell_area.tex';
fid = fopen(filename,'wt');
fprintf(fid,'%s\n','\begin{table}[ht]');
fprintf(fid,'%s\n','\begin{center}');
fprintf(fid,'%s\n','\caption{Summary of testing examples for the cut cell algorithm.}');
fprintf(fid,'%s\n','\begin{tabular}{|c|c|c|c|}\hline');
fprintf(fid,'%s\n','Case number  &  Triangle Area  &  NXP  &  Cutting Area  \\ \hline');
fprintf(fid,'%s %6.4f %s %6.4f %s\n', '1    &',  tri_area(1), '&  3     &', cut_area(1), '\\');
fprintf(fid,'%s %6.4f %s %6.4f %s\n', '2    &',  tri_area(2), '&  3     &', cut_area(2), '\\');
fprintf(fid,'%s %6.4f %s %6.4f %s\n', '3    &',  tri_area(3), '&  3     &', cut_area(3), '\\');
fprintf(fid,'%s %6.4f %s %6.4f %s\n', '4    &',  tri_area(4), '&  3     &', cut_area(4), '\\');
fprintf(fid,'%s %6.4f %s %6.4f %s\n', '5    &',  tri_area(5), '&  1     &', cut_area(5), '\\');
fprintf(fid,'%s %6.4f %s %6.4f %s\n', '6    &',  tri_area(6), '&  3     &', cut_area(6), '\\');
fprintf(fid,'%s %6.4f %s %6.4f %s\n', '7    &',  tri_area(7), '&  3     &', cut_area(7), '\\');
fprintf(fid,'%s %6.4f %s %6.4f %s\n', '8    &',  tri_area(8), '&  3     &', cut_area(8), '\\');
fprintf(fid,'%s %6.4f %s %6.4f %s\n', '9    &',  tri_area(9), '&  2     &', cut_area(9), '\\');
fprintf(fid,'%s %6.4f %s %6.4f %s\n', '10   &',  tri_area(10), '&  5     &', cut_area(10), '\\');
fprintf(fid,'%s %6.4f %s %6.4f %s\n', '11   &',  tri_area(11), '&  7     &', cut_area(11), '\\');
fprintf(fid,'%s %6.4f %s %6.4f %s\n', '12   &',  tri_area(11), '&  9     &', cut_area(12), '\\');
fprintf(fid,'%s\n','\hline');
fprintf(fid,'%s\n','\end{tabular}');
fprintf(fid,'%s\n','\end{center}');
fprintf(fid,'%s\n','\label{tab_cutcell_areas}');
fprintf(fid,'%s\n','\end{table}');
fclose(fid);


