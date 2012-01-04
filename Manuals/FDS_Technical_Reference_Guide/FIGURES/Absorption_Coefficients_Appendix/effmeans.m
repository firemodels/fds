function effmeans

addpath Two_flux
addpath Data

% BEGIN PLOT_STYLE
% font properties
Font_Name       = 'Times';
Key_Font_Size   = 12;
Title_Font_Size = 16;
Label_Font_Size = 16;

% line properties
Line_Width      = 1.0;

% plot properties
Plot_Units      = 'inches';
Plot_Width      = 5.0;
Plot_Height     = 3.4;
Plot_X          = 1.2;
Plot_Y          = 0.8;

% paper properties
Paper_Units     = 'inches';
Paper_Width     = 6.5;
Paper_Height    = 4.5;

% print properties
Figure_Visibility = 'on';

% svn text position
SVN_Scale_X = 0.80;
SVN_Scale_Y = 1.05;

set(gcf,'DefaultLineLineWidth',Line_Width)
WPos = get(gcf,'Position');
set(gcf,'Position',[WPos(1) WPos(2) 640,420]);
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

%END PLOT_STYLE

x=0.0001:0.0001:0.005;
Kappa=zeros(length(x),5);

N=20;

sigma=5.67040040*10^-8; 

qin=sigma*1450^4;

temp=csvread('abscoeffs.csv',2);

v=temp(:,1);
for i=2:size(temp,2)
    a=temp(:,i);
    %ind=a>0;
    fprintf('%d\n',i)
    Kappa(:,i-1)=effmean3(a,v,x,1450);
end

plot(x*1e3,Kappa(:,[3 5]),'LineWidth',2,'LineSmoothing','on');
%semilogy(x*1e3,Kappa,'LineWidth',2,'LineSmoothing','on');
xlabel('Path length (mm)','FontSize',14,'FontName',Font_Name);
ylabel('Absorption coefficient (1/m)','FontSize',14,'FontName',Font_Name);
legend({'Water','Diesel'},'FontSize',14,'FontName',Font_Name);
%set(gca,'FontSize',14);
%axis([0 5 500 2500]);
% print to pdf
        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[Paper_Width Paper_Height]);
        set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
        display(['Printing plot ',num2str(j),'...'])
%        print(gcf,'-dpdf','../KAPPA_PATH_LENGTH')

exportfig(gcf,'../KAPPA_VS_PATH_LENGTH_WD.pdf','Renderer','painters', 'width',6.5,'height',5 ,'fontsize',1.2,...
            'Color','CMYK','Format','pdf','Resolution',900);
figure;
plot(x*1e3,Kappa(:,[1:2 4]),'LineWidth',2,'LineSmoothing','on');
%semilogy(x*1e3,Kappa,'LineWidth',2,'LineSmoothing','on');
xlabel('Path length (mm)','FontSize',14,'FontName',Font_Name);
ylabel('Absorption coefficient (1/m)','FontSize',14,'FontName',Font_Name);
legend({'Toluene','Methanol','Benzene'},'FontSize',14,'FontName',Font_Name);
%set(gca,'FontSize',14);
axis([0 5 0 600]);
% print to pdf
        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[Paper_Width Paper_Height]);
        set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
        display(['Printing plot ',num2str(j),'...'])
%        print(gcf,'-dpdf','../KAPPA_PATH_LENGTH')

exportfig(gcf,'../KAPPA_VS_PATH_LENGTH_TMB.pdf','Renderer','painters', 'width',6.5,'height',5 ,'fontsize',1.2,...
            'Color','CMYK','Format','pdf','Resolution',900);
        
clear i
Temp=linspace(800,1500,length(x));
x=0.003;
for i=1:length(Temp),
 
    Kappa(i,1)=effmean3(temp(:,2),v,x,Temp(i)); 
    Kappa(i,2)=effmean3(temp(:,3),v,x,Temp(i));
    Kappa(i,3)=effmean3(temp(:,4),v,x,Temp(i));
    Kappa(i,4)=effmean3(temp(:,5),v,x,Temp(i));
    Kappa(i,5)=effmean3(temp(:,6),v,x,Temp(i));
end
figure;
plot(Temp,Kappa(:,[3 5]),'LineWidth',2,'LineSmoothing','on');
%semilogy(Temp,Kappa,'LineWidth',2,'LineSmoothing','on');
xlabel('Blackbody temperature (K)','FontSize',14,'FontName',Font_Name);
ylabel('Absorption coefficient (1/m)','FontSize',14,'FontName',Font_Name);
legend({'Water','Diesel'},'FontSize',14,'FontName',Font_Name);
%set(gca,'FontSize',14);
% print to pdf
        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[Paper_Width Paper_Height]);
        set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
        display(['Printing plot ',num2str(i),'...'])
        %print(gcf,'-dpdf','../KAPPA_TEMP')
exportfig(gcf,'../KAPPA_VS_TEMP_WD.pdf','Renderer','painters', 'width',6.5,'height',4.8 ,'fontsize',1.2,...
            'Color','CMYK','Format','pdf','Resolution',900);

figure;
plot(Temp,Kappa(:,[1:2 4]),'LineWidth',2,'LineSmoothing','on');
%semilogy(Temp,Kappa,'LineWidth',2,'LineSmoothing','on');
xlabel('Blackbody temperature (K)','FontSize',14,'FontName',Font_Name);
ylabel('Absorption coefficient (1/m)','FontSize',14,'FontName',Font_Name);
legend({'Toluene','Methanol','Benzene'},'FontSize',14,'FontName',Font_Name);
%set(gca,'FontSize',14);
% print to pdf
        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[Paper_Width Paper_Height]);
        set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
        display(['Printing plot ',num2str(i),'...'])
        %print(gcf,'-dpdf','../KAPPA_TEMP')
exportfig(gcf,'../KAPPA_VS_TEMP_TMB.pdf','Renderer','painters', 'width',6.5,'height',4.8 ,'fontsize',1.2,...
            'Color','CMYK','Format','pdf','Resolution',900);
end
