% McDermott
% 8-2-2018
% tunnel_pressure_drop.m
% Test of tunnel friction factor
%
% Case matrix
% -----------
%
% Case A: Vel = 2 m/s, Roughness = 1.E-4 m
% Case B: Vel = 2 m/s, Roughness = 1.E-1 m
% Case C: Vel = 10 m/s, Roughness = 1.E-4 m
% Case D: Vel = 10 m/s, Roughness = 1.E-1 m
% Case E: Vel = 4 m/s, Roughness = 1.E-2 m

close all
clear all

datadir = '../../../out/Moody_Chart/';
pltdir = '../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/';

plot_style

chid='tunnel_pressure_drop';
cases={'a','b','c','d','e'};
CASES={'A','B','C','D','E'};
res={'10','20'};
markers={'ko','k+'};
lines={'k--','k:'};

VEL = [2 2 10 10 4];
s = [1.0E-4 1.0E-1 1.0E-4 1.0E-1 1.0E-2];  % sand grain roughness (m) from input file
H = [10 10 10 10 7.2];      % tunnel height (m) from input file
L = [1600 1600 1600 1600 1600];    % tunnel length (m)

f_save = zeros(1,length(cases));
f_fds_save = zeros(length(res),length(cases));

for i=1:length(cases)

    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
    nh= 1; % legend handle index

    x = [10:10:90,100:100:L(i)]';

    for j=1:length(res)

        % check steady state

        M = importdata([datadir,chid,'_',cases{i},'_',res{j},'_devc.csv'],',',2);
        t = M.data(:,find(strcmp(M.colheaders,'Time')));
        U = M.data(:,find(strcmp(M.colheaders,'"UBAR"')));
        % plot(t,U)
        % xlabel('time (s)')
        % ylabel('u velocity (m/s)')

        % compute friction factor (f) from Colebrook equation
        mu = M.data(end,find(strcmp(M.colheaders,'"MU"'))); % dynamic viscosity of AIR
        rho = M.data(end,find(strcmp(M.colheaders,'"RHO"'))); % density of AIR

        Re = rho*H(i)*VEL(i)/mu;
        [f,error,iter] = colebrook(Re,s(i)/H(i),.001,1e-9);
        dpdx_exact = -f/H(i) * 0.5 * rho*VEL(i)^2;
        f_save(i) = f;

        % pressure drop

        P = M.data(end,find(strcmp(M.colheaders,'"P10"')):end)';

        hh(nh)=plot(x,P,markers{j}); hold on;
        nh=nh+1;

        set(gca,'FontName',Font_Name)
        set(gca,'FontSize',Label_Font_Size)

        % least squares to get slope (pressure drop)

        rsub = find(x>=10);
        xsub = x(rsub);
        A = [ones(length(xsub),1), xsub];
        y = inv(A'*A)*(A'*P(rsub));
        hh(nh)=plot(x,[y(1)+x*y(2)],lines{j}); nh=nh+1;

        % compute friction factor from DEVC pressure drop

        dpdx = y(2); % pressure drop (Pa/m)
        f_fds = 2*(-dpdx)*H(i)/(rho*U(end)^2);  % f from FDS
        f_fds_save(j,i) = f_fds;

    end

    hh(nh)=plot(x,[(x-L(i))*dpdx_exact],'k-');

    lh=legend(hh,'DEVC Pressure 10','Least squares fit 10','DEVC Pressure 20','Least squares fit 20','Exact pressure drop');
    set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

    % add plot title

    Min_Ind = 0;
    Max_Ind = L(i);
    YLimits = get(gca,'YLim');
    Min_Dep = YLimits(1);
    Max_Dep = 1.2*YLimits(2);
    Title_Position = [0.05 0.90];
    X_Title_Position = Min_Ind+Title_Position(1)*(Max_Ind-Min_Ind);
    Y_Title_Position = Min_Dep+Title_Position(2)*(Max_Dep-Min_Dep);
    Plot_Title = ['Case ',CASES{i}];
    text(X_Title_Position,Y_Title_Position,Plot_Title,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter);
    Plot_Title = ['{\itU} = ',num2str(VEL(i)),' m/s, {\its} = ',num2str(s(i)),' m'];
    text(X_Title_Position,0.9*Y_Title_Position,Plot_Title,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter);

    % axis labels

    axis([Min_Ind Max_Ind Min_Dep Max_Dep])
    xlabel('Distance (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('Gauge Pressure (Pa)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    set(gca,'XTick',[0 400 800 1200 1600])

    % add version string if file is available

    Git_Filename = [datadir,chid,'_',cases{i},'_10_git.txt'];
    addverstr(gca,Git_Filename,'linear')

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[pltdir,chid,'_',cases{i}])

end

% compute errors

for i=1:length(cases)
    err10 = abs((f_save(i)-f_fds_save(1,i))/f_save(i));
    err20 = abs((f_save(i)-f_fds_save(2,i))/f_save(i));
    max_error(i) = max(err10,err20)*100;
end

% write friction factors to latex

fid = fopen([pltdir,'tunnel_pressure_drop.tex'],'wt');
fprintf(fid, '%s\n', '\scriptsize');
fprintf(fid, '%s\n', '\caption[Friction factors in tunnels]{Friction factors for {\ct tunnel\_pressure\_drop} cases.}');
fprintf(fid, '%s\n', '\label{tab:tunnel_pressure_drop}');
fprintf(fid, '%s\n', '\centering');
fprintf(fid, '%s\n', '\begin{tabular}{lccccccc}');
fprintf(fid, '%s\n', '\hline');
fprintf(fid, '%s\n', 'Case    & Velocity (m/s) & Roughness (m) & Hydraulic Dia. (m) & $f$ Colebrook & $f$ FDS 10 & $f$ FDS 20 & Max Rel. Error (\%)\\');
fprintf(fid, '%s\n', '\hline');
fprintf(fid, '%s\n', ['A      & ',num2str(VEL(1)),' & ',num2str(s(1)),' & ',num2str(H(1)),' & ',num2str(f_save(1),3),' & ',num2str(f_fds_save(1,1),3),' & ',num2str(f_fds_save(2,1),3),' & ',num2str(max_error(1),2),' \\']);
fprintf(fid, '%s\n', ['B      & ',num2str(VEL(2)),' & ',num2str(s(2)),' & ',num2str(H(2)),' & ',num2str(f_save(2),3),' & ',num2str(f_fds_save(1,2),3),' & ',num2str(f_fds_save(2,2),3),' & ',num2str(max_error(2),2),' \\']);
fprintf(fid, '%s\n', ['C      & ',num2str(VEL(3)),' & ',num2str(s(3)),' & ',num2str(H(3)),' & ',num2str(f_save(3),3),' & ',num2str(f_fds_save(1,3),3),' & ',num2str(f_fds_save(2,3),3),' & ',num2str(max_error(3),2),' \\']);
fprintf(fid, '%s\n', ['D      & ',num2str(VEL(4)),' & ',num2str(s(4)),' & ',num2str(H(4)),' & ',num2str(f_save(4),3),' & ',num2str(f_fds_save(1,4),3),' & ',num2str(f_fds_save(2,4),3),' & ',num2str(max_error(4),2),' \\']);
fprintf(fid, '%s\n', ['E      & ',num2str(VEL(5)),' & ',num2str(s(5)),' & ',num2str(H(5)),' & ',num2str(f_save(5),3),' & ',num2str(f_fds_save(1,5),3),' & ',num2str(f_fds_save(2,5),3),' & ',num2str(max_error(5),2),' \\']);
fprintf(fid, '%s\n', '\hline');
fprintf(fid, '%s\n', '\end{tabular}');
fprintf(fid, '%s\n', '\normalsize');
fclose(fid);














