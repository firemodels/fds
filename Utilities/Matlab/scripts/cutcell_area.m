% Charles Luo
% 4-20-2012
% cutcell_area.m
%

close all
clear all

dir = '../../Verification/Immersed_Boundary_Method/';
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
    total_area(kcase) = 0.;
    switch kcase
        case 10
            target_cell = 4;
        case 11
            target_cell = 4;
        case 12
            target_cell = 6;
        otherwise
            target_cell = 1;
    end
    if ~exist([dir,infile{kcase}])
        display(['Error: File ',[dir,infile{kcase}],' does not exist. Skipping case.'])
        return
    end
    M = csvread([dir,infile{kcase}],1,0);
    [nrow,ncol] = size(M);
    for i=1:nrow
        icell = M(i,2);
        if (icell == target_cell)
            cut_area(kcase) = M(i,3);
        end
        total_area(kcase) = total_area(kcase)+M(i,3);
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
fprintf(fid,'%s\n','\caption[Summary of examples for the cut cell algorithm]{\label{tab_cutcell_areas} Summary of examples for the cut cell algorithm.}');
fprintf(fid,'%s\n','\begin{tabular}{|c|c|c|c|c|c|}\hline');
fprintf(fid,'%s\n','Case number  &  Illustration  &  Triangle Area  &  NXP  &  Cutcell Area  & Total Cutcell Area \\ \hline');
fprintf(fid,'%s %s %6.4f %s %6.4f %s %6.4f %s\n', '1    &',  '\includegraphics[totalheight=2.in]{SCRIPT_FIGURES/tri_cube_cut_cell_test_01_image}   &',   tri_area(1), '&  3     &', cut_area(1), '&    ', total_area(1), '\\ \noalign{\smallskip}\hline\noalign{\smallskip}');
fprintf(fid,'%s %s %6.4f %s %6.4f %s %6.4f %s\n', '2    &',  '\includegraphics[totalheight=2.in]{SCRIPT_FIGURES/tri_cube_cut_cell_test_02_image}   &',   tri_area(2), '&  3     &', cut_area(2), '&    ', total_area(2), '\\ \noalign{\smallskip}\hline\noalign{\smallskip}');
fprintf(fid,'%s %s %6.4f %s %6.4f %s %6.4f %s\n', '3    &',  '\includegraphics[totalheight=2.in]{SCRIPT_FIGURES/tri_cube_cut_cell_test_03_image}   &',   tri_area(3), '&  3     &', cut_area(3), '&    ', total_area(3), '\\ \noalign{\smallskip}\hline\noalign{\smallskip}');
fprintf(fid,'%s %s %6.4f %s %6.4f %s %6.4f %s\n', '4    &',  '\includegraphics[totalheight=2.in]{SCRIPT_FIGURES/tri_cube_cut_cell_test_04_image}   &',   tri_area(4), '&  3     &', cut_area(4), '&    ', total_area(4), '\\ \noalign{\smallskip}\hline\noalign{\smallskip}');
fprintf(fid,'%s\n','\end{tabular}');
fprintf(fid,'%s\n','\end{center}');
fprintf(fid,'%s\n','\end{table}');
fprintf(fid,'\n');

fprintf(fid,'%s\n','\begin{table}[ht]');
fprintf(fid,'%s\n','\begin{center}');
fprintf(fid,'%s\n','\begin{tabular}{|c|c|c|c|c|c|}\hline');
fprintf(fid,'%s\n','Case number  &  Illustration  &  Triangle Area  &  NXP  &  Cutcell Area  & Total Cutcell Area \\ \hline');
fprintf(fid,'%s %s %6.4f %s %6.4f %s %6.4f %s\n', '5    &',  '\includegraphics[totalheight=2.in]{SCRIPT_FIGURES/tri_cube_cut_cell_test_05_image}   &',   tri_area(5), '&  1     &', cut_area(5), '&    ', total_area(5), '\\ \noalign{\smallskip}\hline\noalign{\smallskip}');
fprintf(fid,'%s %s %6.4f %s %6.4f %s %6.4f %s\n', '6    &',  '\includegraphics[totalheight=2.in]{SCRIPT_FIGURES/tri_cube_cut_cell_test_06_image}   &',   tri_area(6), '&  3     &', cut_area(6), '&    ', total_area(6), '\\ \noalign{\smallskip}\hline\noalign{\smallskip}');
fprintf(fid,'%s %s %6.4f %s %6.4f %s %6.4f %s\n', '7    &',  '\includegraphics[totalheight=2.in]{SCRIPT_FIGURES/tri_cube_cut_cell_test_07_image}   &',   tri_area(7), '&  3     &', cut_area(7), '&    ', total_area(7), '\\ \noalign{\smallskip}\hline\noalign{\smallskip}');
fprintf(fid,'%s %s %6.4f %s %6.4f %s %6.4f %s\n', '8    &',  '\includegraphics[totalheight=2.in]{SCRIPT_FIGURES/tri_cube_cut_cell_test_08_image}   &',   tri_area(8), '&  3     &', cut_area(8), '&    ', total_area(8), '\\ \noalign{\smallskip}\hline\noalign{\smallskip}');
fprintf(fid,'%s\n','\end{tabular}');
fprintf(fid,'%s\n','\end{center}');
fprintf(fid,'%s\n','\end{table}');
fprintf(fid,'\n');

fprintf(fid,'%s\n','\begin{table}[ht]');
fprintf(fid,'%s\n','\begin{center}');
fprintf(fid,'%s\n','\begin{tabular}{|c|c|c|c|c|c|}\hline');
fprintf(fid,'%s\n','Case number  &  Illustration  &  Triangle Area  &  NXP  &  Cutcell Area  & Total Cutcell Area \\ \hline');
fprintf(fid,'%s %s %6.4f %s %6.4f %s %6.4f %s\n', '9    &',  '\includegraphics[totalheight=2.in]{SCRIPT_FIGURES/tri_cube_cut_cell_test_09_image}   &',   tri_area(9),  '&  2     &', cut_area(9), '&    ', total_area(9), '\\ \noalign{\smallskip}\hline\noalign{\smallskip}');
fprintf(fid,'%s %s %6.4f %s %6.4f %s %6.4f %s\n', '10   &',  '\includegraphics[totalheight=2.in]{SCRIPT_FIGURES/tri_cube_cut_cell_test_10_image}   &',   tri_area(10), '&  5     &', cut_area(10), '&    ', total_area(10), '\\ \noalign{\smallskip}\hline\noalign{\smallskip}');
fprintf(fid,'%s %s %6.4f %s %6.4f %s %6.4f %s\n', '11   &',  '\includegraphics[totalheight=2.in]{SCRIPT_FIGURES/tri_cube_cut_cell_test_11_image}   &',   tri_area(11), '&  7     &', cut_area(11), '&    ', total_area(11), '\\ \noalign{\smallskip}\hline\noalign{\smallskip}');
fprintf(fid,'%s %s %6.4f %s %6.4f %s %6.4f %s\n', '12   &',  '\includegraphics[totalheight=2.in]{SCRIPT_FIGURES/tri_cube_cut_cell_test_12_image}   &',   tri_area(12), '&  9     &', cut_area(12), '&    ', total_area(12), '\\ \noalign{\smallskip}\hline\noalign{\smallskip}');
fprintf(fid,'%s\n','\end{tabular}');
fprintf(fid,'%s\n','\end{center}');
fprintf(fid,'%s\n','\end{table}');

fclose(fid);


