% McDermott
% 4-22-2001 (modified 12-9-2010)
% plotspec.m
%
% Read and plot energy spectra for Comte-Bellot and Corrsin data
%
% [] = plotspec(chid,N)
%
% chid = CHID from FDS input file
% N = number of cells in 1D
%
% Example: >> plotspec('csmag_32',32)

function [] = plotspec(chid,N)

close all
addpath('../../Verification/Turbulence')

% set FDS standard plot format
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gcf,'DefaultLineLineWidth',Line_Width)

L = 9*2*pi/100; % box length (m)
k0 = 2*pi/L;
kc = 1/2*N*k0;

uvw_file1 = [chid,'_uvw_001.csv'];
uvw_file2 = [chid,'_uvw_002.csv'];
uvw_file3 = [chid,'_uvw_003.csv'];

skip_case = 0;
if ~exist(uvw_file1)
    display(['Error: File ' uvw_file1 ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist(uvw_file2)
    display(['Error: File ' uvw_file2 ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist(uvw_file3)
    display(['Error: File ' uvw_file3 ' does not exist. Skipping case.'])
    skip_case = 1;
end
if skip_case
    return
end

% Plot the FDS data
H(1) = plotspec_uvw(uvw_file1,'k.-'); hold on
H(2) = plotspec_uvw(uvw_file2,'r.-');
H(3) = plotspec_uvw(uvw_file3,'b.-');

% Gather the Comte-Bellot/Corrsin data
CBC = load('cbcdata.txt');
k  = CBC(:,1)*1e2;
E1 = CBC(:,2)/1e6;
E2 = CBC(:,3)/1e6;
E3 = CBC(:,4)/1e6;

% Plot the CBC data
loglog(k,E1,'k-');
loglog(k,E2,'k-');
loglog(k,E3,'k-');

% format axes
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Key_Font_Size)
axis([10 1000 1/1e6 1000/1e6])
ylabel('{\it E(k)} (m^3/s^2)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
xlabel('{\it k} (1/m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

% Add lines to indicate cutoff wavenumbers
loglog( [k0*N/2,k0*N/2],[1e-10,1e-2],'k--'); % LES Nyquist limit

anno_obj_handle=annotation('textarrow',[.56,.44],[.85,.55],'String','Time');
set(anno_obj_handle,'FontName',Font_Name)

% % Apply transfer function (optional)
% delta = pi/kc;
% for j = 1:length(k)
%     if k(j) <= kc
%         %G(j) = sin(1/2*k(j)*delta)/(1/2*k(j)*delta); % box filter
%         %G(j) = exp(-k(j)^2*delta^2/24);              % Gaussian filter
%         G(j) = 1;                                     % spectral cutoff
%     elseif k(j) > kc & k(j) <= sqrt(2)*kc
%         G(j) = sqrt(3*kc/k(j) - 2);                   % RJM implicit filter
%     elseif k(j) > sqrt(2)*kc
%         G(j) = 0;
%     end
%     E1(j) = G(j)*G(j)*E1(j);
%     E2(j) = G(j)*G(j)*E2(j);
%     E3(j) = G(j)*G(j)*E3(j);
% end
% 
% % Plot the filtered CBC data
% loglog(k,E1,'ko','Linewidth',2)
% loglog(k,E2,'ro','Linewidth',2)
% loglog(k,E3,'bo','Linewidth',2)

% add SVN if file is available

SVN_Filename = [chid,'_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/',chid,'_spectra'])
    
