% McDermott
% 5-21-14
% openmp_timing_benchmarks.m

close all
clear all

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

H(1)=plot(ncores,time64,'b^--'); hold on
H(2)=plot(ncores,time128,'rsq--');

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([1 8 0 100])

xlabel('Number of OpenMP threads')
ylabel('Relative clock time (%)')

h=legend(H,'64^3','128^3');
set(h,'Interpreter',Font_Interpreter)

% add SVN if file is available

svn_file = [dir,'openmp_test64a_git.txt'];
addverstr(gca,svn_file,'linear')
% if exist(svn_file,'file')
%     SVN = importdata(svn_file);
%     x_lim = get(gca,'XLim');
%     y_lim = get(gca,'YLim');
%     X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
%     Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
%     text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
%         'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
% end

% print to pdf
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/openmp_timing_benchmarks')

% check errors
if time64(4) > 56.
   display(['Matlab Warning: Timing for openmp_test64 out of tolerance.',num2str(time64(4))])
end
if time128(4) > 56.
   display(['Matlab Warning: Timing for openmp_test128 out of tolerance.',num2str(time128(4))])
end

