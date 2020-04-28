% McDermott
% 5-21-14
% openmp_timing_benchmarks.m

close all
clear all

plot_style

dir = '../../Verification/Timing_Benchmarks/';

ncores = [1,2,3,4,5,6,7,8];
a = {'a','b','c','d','e','f','g','h'};

for i=1:length(a)
    filename = [dir,'openmp_test64',a{i},'_devc.csv'];
    if ~exist(filename)
        display(['Error: File ' filename ' does not exist. Skipping case.'])
        return
    end
    M = importdata(filename,',',2);
    j = find(strcmp(M.colheaders,'"clock time"'));
    time64(i) = M.data(end,j);
end
time64 = time64/time64(1) * 100;

for i=1:length(a)
    filename = [dir,'openmp_test128',a{i},'_devc.csv'];
    if ~exist(filename)
        display(['Error: File ' filename ' does not exist. Skipping case.'])
        return
    end
    M = importdata(filename,',',2);
    j = find(strcmp(M.colheaders,'"clock time"'));
    time128(i) = M.data(end,j);
end
time128 = time128/time128(1) * 100;

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H(1)=plot(ncores,time64,'b^--'); hold on
H(2)=plot(ncores,time128,'rsq--');

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([1 8 0 100])

xlabel('Number of OpenMP threads')
ylabel('Relative clock time (%)')

h=legend(H,'64^3','128^3');
set(h,'Interpreter',Font_Interpreter)
set(h,'FontSize',Key_Font_Size)

% add version string if file is available

git_file = [dir,'openmp_test64a_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/openmp_timing_benchmarks')

% check errors
if time64(4) > 55.
   display(['Matlab Warning: Timing for openmp_test64 out of tolerance.',num2str(time64(4))])
end
if time128(4) > 55.
   display(['Matlab Warning: Timing for openmp_test128 out of tolerance.',num2str(time128(4))])
end

