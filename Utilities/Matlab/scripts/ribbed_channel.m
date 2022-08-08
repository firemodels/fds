% McDermott
% 9-18-12
% ribbed_channel.m

close all
clear all

outdir = '../../../out/Casara_Arts_Ribbed_Channel/';
expdir = '../../../exp/Casara_Arts_Ribbed_Channel/';
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

plot_style

nx = [20 40 80 160]; % not correct at the moment, just referencing old file names
lnx = length(nx);
Ub=6.2; % exp bulk velocity
L= 0.3; % channel length
D = 0.1; % channel height
h = 0.03;
dx = L./[10 20 40 80];
fds_marker = {'r+-' 'c^-' 'g>-' 'k-'};
fds_key = {'FDS {\it h/\deltax}=3' 'FDS {\it h/\deltax}=6' 'FDS {\it h/\deltax}=12' 'FDS {\it h/\deltax}=24'};
geom = {'_' '_geom_'};

if ~exist([expdir,'ribbed_channel_data.csv'])
    display(['Error: File ' [expdir,'ribbed_channel_data.csv'] ' does not exist. Skipping case.'])
    return
end

DATA = importdata([expdir,'ribbed_channel_data.csv'],',',1);

for ii=1:length(geom)

    % check existence of device files

    for i=1:lnx
        if ~exist([outdir,'ribbed_channel',geom{ii},num2str(nx(i)),'_devc.csv'])
            display(['Error: File ' [outdir,'ribbed_channel',geom{ii},num2str(nx(i)),'_devc.csv'] ' does not exist. Skipping case.'])
            return
        end
    end

    % check bulk velocity is 6.2 m/s

    for i=1:lnx
        M{i} = importdata([outdir,'ribbed_channel',geom{ii},num2str(nx(i)),'_devc.csv'],',',2);
    end

    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    H(1)=plot([0 100],Ub*ones(1,2),'k-','linewidth',2); hold on

    for i=1:lnx
        t_fds = M{i}.data(:,1);
        Ub_fds = M{i}.data(:,2);
        H(1+i)=plot(t_fds,Ub_fds,fds_marker{i});
        t_range = find(t_fds>2);
        if abs(mean(Ub_fds(t_range))-Ub)/Ub > 0.01
            disp(['Matlab Warning: Ub mean nx ',geom{ii},num2str(nx(i)),' = ',num2str(mean(Ub_fds(t_range)))])
        end
    end

    xlabel('Time (s)','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
    ylabel('Bulk Velocity (m/s)','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)

    axis([0 10 0 1.5*Ub])

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Title_Font_Size)

    lh = legend(H,'Exp Bulk Velocity',fds_key{1:lnx},'Location','SouthEast');
    set(lh,'Interpreter',Font_Interpreter)

    % add Git revision if file is available

    Git_Filename = [outdir,'ribbed_channel',geom{ii},num2str(nx(1)),'_git.txt'];
    addverstr(gca,Git_Filename,'linear')

    % print to pdf

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[plotdir,'ribbed_channel',geom{ii},'Ubulk'])

    % check existence of line device files

    for i=1:lnx
        if ~exist([outdir,'ribbed_channel',geom{ii},num2str(nx(i)),'_line.csv'])
            display(['Error: File ' [outdir,'ribbed_channel',geom{ii},num2str(nx(i)),'_line.csv'] ' does not exist. Skipping case.'])
            return
        end
    end

    % read experimental and FDS data files

    for i=1:lnx
        M{i} = importdata([outdir,'ribbed_channel',geom{ii},num2str(nx(i)),'_line.csv'],',',2);
    end

    % organize and plot streamwise U along bottom of channel

    j = find(strcmp(DATA.colheaders,'x/h U strm'));
    xoh = DATA.data(:,j);
    j = find(strcmp(DATA.colheaders,'U strm'));
    u_data = DATA.data(:,j);

    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    H(1)=plot(xoh,u_data,'bo'); hold on

    for i=1:lnx
        j = find(strcmp(M{i}.colheaders,'u_strm_bot-x'));
        x = M{i}.data(:,j);
        I = find(x<0);
        x(I) = x(I) + L;
        [x I] = sort(x);

        j = find(strcmp(M{i}.colheaders,'u_strm_bot'));
        u_fds = M{i}.data(:,j);
        u_fds = u_fds(I);

        H(1+i)=plot(x/h,u_fds/Ub,fds_marker{i});
    end

    xlabel('\it x/h','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
    ylabel('\it U/U_b','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)

    axis([0 10 -0.6 0.6])
    set(gca,'YTick',-0.6:0.3:0.6)
    set(gca,'YTickLabel',-0.6:0.3:0.6)

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Title_Font_Size)

    lh = legend(H,'PIV data',fds_key{1:lnx},'Location','Northwest');
    set(lh,'Interpreter',Font_Interpreter)

    % add Git revision if file is available

    Git_Filename = [outdir,'ribbed_channel',geom{ii},num2str(nx(1)),'_git.txt'];
    addverstr(gca,Git_Filename,'linear')

    % print to pdf

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[plotdir,'ribbed_channel',geom{ii},'u_strm'])

    % streamwise urms along bottom of channel

    j = find(strcmp(DATA.colheaders,'x/h urms strm'));
    xoh = DATA.data(:,j);
    j = find(strcmp(DATA.colheaders,'urms strm'));
    urms_data = DATA.data(:,j);

    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    H(1)=plot(xoh,urms_data,'bo'); hold on

    for i=1:lnx
        j = find(strcmp(M{i}.colheaders,'urms_strm_bot-x'));
        x = M{i}.data(:,j);
        I = find(x<0);
        x(I) = x(I) + L;
        [x I] = sort(x);

        j = find(strcmp(M{i}.colheaders,'urms_strm_bot'));
        urms_fds = M{i}.data(:,j);
        urms_fds = urms_fds(I);

        H(1+i)=plot(x/h,urms_fds/Ub,fds_marker{i});
    end

    xlabel('\it x/h','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
    yh = ylabel('\it U_{rms}/U_b','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
    set(yh,'Position',[-0.8194-.2    0.2963-.025    1.0001])

    axis([0 10 0 0.6])

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Title_Font_Size)

    lh = legend(H,'PIV data',fds_key{1:lnx},'Location','Northwest');
    set(lh,'Interpreter',Font_Interpreter)

    % add Git revision if file is available

    Git_Filename = [outdir,'ribbed_channel',geom{ii},num2str(nx(1)),'_git.txt'];
    addverstr(gca,Git_Filename,'linear')

    % print to pdf

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[plotdir,'ribbed_channel',geom{ii},'urms_strm'])

    % streamwise U profile at x/h=0 (center of rib)

    j = find(strcmp(DATA.colheaders,'y/h U prof'));
    yoh = DATA.data(:,j);
    j = find(strcmp(DATA.colheaders,'U prof'));
    u_data = DATA.data(:,j);

    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    H(1)=plot(u_data,yoh,'bo'); hold on

    for i=1:lnx
        j = find(strcmp(M{i}.colheaders,'u_prof_rib-z'));
        I = find( M{i}.data(:,j)>h & M{i}.data(:,j)<D );
        y = M{i}.data(I,j);

        j = find(strcmp(M{i}.colheaders,'u_prof_rib'));
        u_fds = M{i}.data(I,j);

        H(1+i)=plot(u_fds/Ub,y/h,fds_marker{i});
    end

    ylabel('\it z/h','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
    xlabel('\it U/U_b','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)

    axis([0 2 1 D/h])

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Title_Font_Size)

    lh = legend(H,'PIV data',fds_key{1:lnx},'Location','Northwest');
    set(lh,'Interpreter',Font_Interpreter)

    % add Git revision if file is available

    Git_Filename = [outdir,'ribbed_channel',geom{ii},num2str(nx(1)),'_git.txt'];
    addverstr(gca,Git_Filename,'linear')

    % print to pdf

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

    print(gcf,'-dpdf',[plotdir,'ribbed_channel',geom{ii},'u_prof'])

    % streamwise urms profile at x/h=0 (center of rib)

    j = find(strcmp(DATA.colheaders,'y/h urms prof'));
    yoh = DATA.data(:,j);
    j = find(strcmp(DATA.colheaders,'urms prof'));
    urms_data = DATA.data(:,j);

    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    H(1)=plot(urms_data,yoh,'bo'); hold on

    for i=1:lnx
        j = find(strcmp(M{i}.colheaders,'urms_prof_rib-z'));
        I = find( M{i}.data(:,j)>h & M{i}.data(:,j)<D );
        y = M{i}.data(I,j);

        j = find(strcmp(M{i}.colheaders,'urms_prof_rib'));
        urms_fds = M{i}.data(I,j);

        H(1+i)=plot(urms_fds/Ub,y/h,fds_marker{i});
    end

    ylabel('\it z/h','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
    xlabel('\it U_{rms}/U_b','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)

    axis([0 1 1 D/h])

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Title_Font_Size)

    lh = legend(H,'PIV data',fds_key{1:lnx},'Location','Northeast');
    set(lh,'Interpreter',Font_Interpreter)

    % add Git revision if file is available

    Git_Filename = [outdir,'ribbed_channel',geom{ii},num2str(nx(1)),'_git.txt'];
    addverstr(gca,Git_Filename,'linear')

    % print to pdf

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[plotdir,'ribbed_channel',geom{ii},'urms_prof'])

end % ii




