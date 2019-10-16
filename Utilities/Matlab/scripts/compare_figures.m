% McGrattan
% 10-16-19
% compare_figures.m

close all
clear all

ref_dir = '../../../../fig/fds/Reference_Figures/';
use_dir = '../../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/';
ver_dir = '../../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

files = dir([ref_dir '*.png']);

fid = fopen('compare_figures.txt','w');

for i=1:length(files)
    
    if isfile([use_dir files(i).name])
        A = imread([use_dir files(i).name]);
    else
        A = imread([ver_dir files(i).name]);
    end
    
    B = imread([ref_dir files(i).name]);

    C = immse(A,B);
    fprintf(fid,'%-50s %0.4f\n',files(i).name,C);

clear A B C

end


