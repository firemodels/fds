%Topi Sikanen
% 3.8.2016


expdir = '../../Validation/Submodules/macfp-db/Liquid_Pool_Fires/Waterloo_Methanol/Experimental_Data/';
fdsdir = '../../Validation/Waterloo_Methanol/FDS_Output_Files/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Waterloo_Methanol/';

for predc={'Prescribed','Predicted'},
pred = predc{1};
F3 = importdata([fdsdir,'Waterloo_Methanol_',pred,'_0p25cm_line.csv'],',',2);
F2 = importdata([fdsdir,'Waterloo_Methanol_',pred,'_0p5cm_line.csv'],',',2);
F1 = importdata([fdsdir,'Waterloo_Methanol_',pred,'_1cm_line.csv'],',',2);
git_file = [fdsdir,'Waterloo_Methanol_',pred,'_0p25cm_git.txt'];


%% Mean and RMS velocities at different z-offsets

E1 = importdata([expdir,'Mean_and_rms_values.csv'],',',1);

for z=[2:2:20,30],
    %Mean radial velocity
    figure;
    plot_style;
    % load data
    var = ['Up',num2str(z)];
    r = F1.data(:,1)*100;
    ubar1 = F1.data(:,find(strcmp(F1.colheaders,var)));
    ubar2 = F2.data(:,find(strcmp(F2.colheaders,var)));
    ubar3 = F3.data(:,find(strcmp(F3.colheaders,var)));

    ecol = find(strncmp(E1.colheaders,'v-bar',5));
    eind = E1.data(:,find(strcmp(E1.colheaders,'z')))==z;
    re   = E1.data(eind,find(strcmp(E1.colheaders,'r')));
    ye   = E1.data(eind,ecol);
    % plot
    H(1)=plot(re,ye,'ksq','MarkerSize',Marker_Size); hold on
    H(2) = plot(r,ubar1,'r-.','LineWidth',Line_Width); % dx = 1 cm
    H(3) = plot(r,ubar2,'m--','LineWidth',Line_Width); % dx = 0.5 cm
    H(4) = plot(r,ubar3,'b-','LineWidth',Line_Width);  % dx = 0.25 cm
    
    xlabel('Radial position (cm)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    ylabel('Radial velocity ( m/s )','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    lh = legend(H,'Exp','FDS 1 cm','FDS 0.5 cm','FDS 0.25 cm');
    set(lh,'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    % add Git revision if file is available
    addverstr(gca,git_file,'linear')

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
    %text(.03,.24,'Waterloo Methanol','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Units','normalized')
    text(.03,.08,['{\it z} = ',num2str(z),' cm'],'FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Units','normalized')

    % print to pdf
    outfile = ['Waterloo_Methanol_',pred,'_',var];
    print(gcf,'-dpdf',[pltdir,outfile])
    close
    
    %Mean horizontal velocity
    figure;
    plot_style;
    % load data
    var = ['Wp',num2str(z)];
    r = F1.data(:,1)*100;
    ubar1 = F1.data(:,find(strcmp(F1.colheaders,var)));
    ubar2 = F2.data(:,find(strcmp(F2.colheaders,var)));
    ubar3 = F3.data(:,find(strcmp(F3.colheaders,var)));

    ecol = find(strncmp(E1.colheaders,'u-bar',5));
    eind = E1.data(:,find(strcmp(E1.colheaders,'z')))==z;
    re   = E1.data(eind,find(strcmp(E1.colheaders,'r')));
    ye   = E1.data(eind,ecol);
    % plot
    H(1)=plot(re,ye,'ksq','MarkerSize',Marker_Size); hold on
    H(2) = plot(r,ubar1,'r-.','LineWidth',Line_Width); % dx = 1 cm
    H(3) = plot(r,ubar2,'m--','LineWidth',Line_Width); % dx = 0.5 cm
    H(4) = plot(r,ubar3,'b-','LineWidth',Line_Width);  % dx = 0.25 cm
    
    xlabel('Radial position (cm)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    ylabel('Horizontal velocity ( m/s )','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    lh = legend(H,'Exp','FDS 1 cm','FDS 0.5 cm','FDS 0.25 cm');
    set(lh,'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    % add Git revision if file is available
    addverstr(gca,git_file,'linear')

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
    %text(.03,.24,'Waterloo Methanol','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Units','normalized')
    text(.03,.08,['{\it z} = ',num2str(z),' cm'],'FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Units','normalized')

    % print to pdf
    outfile = ['Waterloo_Methanol_',pred,'_',var];
    print(gcf,'-dpdf',[pltdir,outfile])
    close
    
    %Mean temperature
    figure;
    plot_style;
    % load data
    var = ['Tp',num2str(z)];
    r = F1.data(:,1)*100;
    ubar1 = F1.data(:,find(strcmp(F1.colheaders,var)));
    ubar2 = F2.data(:,find(strcmp(F2.colheaders,var)));
    ubar3 = F3.data(:,find(strcmp(F3.colheaders,var)));

    ecol = find(strncmp(E1.colheaders,'T-bar',5));
    eind = E1.data(:,find(strcmp(E1.colheaders,'z')))==z;
    re   = E1.data(eind,find(strcmp(E1.colheaders,'r')));
    ye   = E1.data(eind,ecol)-273.15;
    % plot
    H(1)=plot(re,ye,'ksq','MarkerSize',Marker_Size); hold on
    H(2) = plot(r,ubar1,'r-.','LineWidth',Line_Width); % dx = 1 cm
    H(3) = plot(r,ubar2,'m--','LineWidth',Line_Width); % dx = 0.5 cm
    H(4) = plot(r,ubar3,'b-','LineWidth',Line_Width);  % dx = 0.25 cm
    
    xlabel('Radial position (cm)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    ylabel('Temperature ( \circC)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    lh = legend(H,'Exp','FDS 1 cm','FDS 0.5 cm','FDS 0.25 cm');
    set(lh,'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    % add Git revision if file is available
    addverstr(gca,git_file,'linear')

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
    %text(.03,.24,'Waterloo Methanol','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Units','normalized')
    text(.03,.08,['{\it z} = ',num2str(z),' cm'],'FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Units','normalized')
    % print to pdf
    outfile = ['Waterloo_Methanol_',pred,'_',var];
    print(gcf,'-dpdf',[pltdir,outfile])
    close
    
     %%%%%%%%%%
     %%%% RMS values
     %%%%%%%%%%
     
     %RMS radial velocity
    figure;
    plot_style;
    % load data
    var = ['Urms',num2str(z)];
    r = F1.data(:,1)*100;
    ubar1 = F1.data(:,find(strcmp(F1.colheaders,var)));
    ubar2 = F2.data(:,find(strcmp(F2.colheaders,var)));
    ubar3 = F3.data(:,find(strcmp(F3.colheaders,var)));

    ecol = find(strncmp(E1.colheaders,'v-rms',5));
    eind = E1.data(:,find(strcmp(E1.colheaders,'z')))==z;
    re   = E1.data(eind,find(strcmp(E1.colheaders,'r')));
    ye   = E1.data(eind,ecol);
    % plot
    H(1)=plot(re,ye,'ksq','MarkerSize',Marker_Size); hold on
    H(2) = plot(r,ubar1,'r-.','LineWidth',Line_Width); % dx = 1 cm
    H(3) = plot(r,ubar2,'m--','LineWidth',Line_Width); % dx = 0.5 cm
    H(4) = plot(r,ubar3,'b-','LineWidth',Line_Width);  % dx = 0.25 cm
    
    xlabel('Radial position (cm)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    ylabel('RMS radial velocity ( m/s )','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    lh = legend(H,'Exp','FDS 1 cm','FDS 0.5 cm','FDS 0.25 cm');
    set(lh,'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    % add Git revision if file is available
    addverstr(gca,git_file,'linear')

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
    %text(.03,.24,'Waterloo Methanol','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Units','normalized')
    text(.03,.08,['{\it z} = ',num2str(z),' cm'],'FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Units','normalized')

    % print to pdf
    outfile = ['Waterloo_Methanol_',pred,'_',var];
    print(gcf,'-dpdf',[pltdir,outfile])
    close
    
    %RMS horizontal velocity
    figure;
    plot_style;
    % load data
    var = ['Wrms',num2str(z)];
    r = F1.data(:,1)*100;
    ubar1 = F1.data(:,find(strcmp(F1.colheaders,var)));
    ubar2 = F2.data(:,find(strcmp(F2.colheaders,var)));
    ubar3 = F3.data(:,find(strcmp(F3.colheaders,var)));

    ecol = find(strncmp(E1.colheaders,'u-rms',5));
    eind = E1.data(:,find(strcmp(E1.colheaders,'z')))==z;
    re   = E1.data(eind,find(strcmp(E1.colheaders,'r')));
    ye   = E1.data(eind,ecol);
    % plot
    H(1)=plot(re,ye,'ksq','MarkerSize',Marker_Size); hold on
    H(2) = plot(r,ubar1,'r-.','LineWidth',Line_Width); % dx = 1 cm
    H(3) = plot(r,ubar2,'m--','LineWidth',Line_Width); % dx = 0.5 cm
    H(4) = plot(r,ubar3,'b-','LineWidth',Line_Width);  % dx = 0.25 cm
    
    xlabel('Radial position (cm)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    ylabel('RMS horizontal velocity ( m/s )','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    lh = legend(H,'Exp','FDS 1 cm','FDS 0.5 cm','FDS 0.25 cm');
    set(lh,'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    % add Git revision if file is available
    addverstr(gca,git_file,'linear')

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
    %text(.03,.24,'Waterloo Methanol','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Units','normalized')
    text(.03,.08,['{\it z} = ',num2str(z),' cm'],'FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Units','normalized')

    % print to pdf
    outfile = ['Waterloo_Methanol_',pred,'_',var];
    print(gcf,'-dpdf',[pltdir,outfile])
    close
    
    %RMS temperature
    figure;
    plot_style;
    % load data
    var = ['Trms',num2str(z)];
    r = F1.data(:,1)*100;
    ubar1 = F1.data(:,find(strcmp(F1.colheaders,var)));
    ubar2 = F2.data(:,find(strcmp(F2.colheaders,var)));
    ubar3 = F3.data(:,find(strcmp(F3.colheaders,var)));

    ecol = find(strncmp(E1.colheaders,'T-rms',5));
    eind = E1.data(:,find(strcmp(E1.colheaders,'z')))==z;
    re   = E1.data(eind,find(strcmp(E1.colheaders,'r')));
    ye   = E1.data(eind,ecol);
    % plot
    H(1)=plot(re,ye,'ksq','MarkerSize',Marker_Size); hold on
    H(2) = plot(r,ubar1,'r-.','LineWidth',Line_Width); % dx = 1 cm
    H(3) = plot(r,ubar2,'m--','LineWidth',Line_Width); % dx = 0.5 cm
    H(4) = plot(r,ubar3,'b-','LineWidth',Line_Width);  % dx = 0.25 cm
    
    xlabel('Radial position (cm)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    ylabel('RMS temperature ( \circC)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    lh = legend(H,'Exp','FDS 1 cm','FDS 0.5 cm','FDS 0.25 cm');
    set(lh,'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
    % add Git revision if file is available
    addverstr(gca,git_file,'linear')

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
    %text(.03,.24,'Waterloo Methanol','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Units','normalized')
    text(.03,.08,['{\it z} = ',num2str(z),' cm'],'FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Units','normalized')
    % print to pdf
    outfile = ['Waterloo_Methanol_',pred,'_',var];
    print(gcf,'-dpdf',[pltdir,outfile])
    close
end
end


