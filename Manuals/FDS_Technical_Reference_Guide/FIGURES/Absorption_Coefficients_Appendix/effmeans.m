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
Tau=zeros(length(x),5);

N=20;

sigma=5.67040040*10^-8; 

qin=sigma*1450^4;

% Lasketaan myös Keefe et al. datoista
temp=importdata('Data/toluene/Keefe_toluene_2/tla105em.y'); % in L /molcm
vstart=temp(1)*10^2;
vend=temp(2)*10^2;
npt=temp(3);
v=linspace(vstart,vend,npt);
a=temp(4:end)*0.1*9339; % L/mol/cm=0.1 m^3/mol/m, Tolueenin tiheys 9339 mol/m^3
Tau(:,1)=effmean3(a,v,x,1450);


% testataan myös metanoli
[a v]=readyspec('Data/Methanol/Bertie/CH3OH/MTHEM93.y');
Tau(:,2)=effmean3(a*2.4719e+004,v,x,1450);
% molar mass 32.04 g mol?1 density 792 kg/m^3
% => density = 2.4719e+004 mol/m^3


% Vesi
[a v]=readyspec('Data/Water/Bertie/WTEREM95.y');
Tau(:,3)=effmean3(a*5.5556e+004,v,x,1450);
% molar mass 18.0153 g/mol density 1000 kg/m^3
% => density =5.5556e+004 mol/m^3

% Bentseeni Data/Benzene/Keefe_benzene/bzh6emf2.y
[a v]=readyspec('Data/Benzene/Keefe_benzene/bzh6emf2.y');
Tau(:,4)=effmean3(a*100*11161,v,x,1450);
% Density 11161 mol/m^3 at T=300K ja P = 1 atm
%kappa0=-log(Tau).'/(x*1e-3).';

% Diesel
A=csvread('Data/Diesel/Sazhin2004.csv',1,0);
[C,IA,IB]=UNION(A(:,1),A(:,1));
A=A(IB,:);
l=A(:,1)*1e-6;
a=A(:,2)*2*pi/l;
a=A(:,2)*2*pi./l;

v=1./l;
Tau(:,5)=effmean3(a,v,x,1450);
plot(x*1e3,Tau,'LineWidth',2,'LineSmoothing','on');
%semilogy(x*1e3,Tau,'LineWidth',2,'LineSmoothing','on');
xlabel('Path length (mm)','FontSize',14,'FontName',Font_Name);
ylabel('Absorption coefficient (1/m)','FontSize',14,'FontName',Font_Name);
legend({'Toluene','Methanol','Water','Benzene','Diesel'},'FontSize',14,'FontName',Font_Name);
%set(gca,'FontSize',14);
axis([0 5 0 2000]);
% print to pdf
        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[Paper_Width Paper_Height]);
        set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
        display(['Printing plot ',num2str(j),'...'])
        print(gcf,'-dpdf','../KAPPA_PATH_LENGTH')

%exportfig(gcf,'KAPPA_VS_PATH_LENGTH.png','Renderer','painters', 'width',14,'height',12 ,'fontsize',1.2,...
%            'Color','bw','Format','png','Resolution',900);
        
clear i
Temp=linspace(800,1500,length(x));
x=0.003;
for i=1:length(Temp),
    % Lasketaan myös Keefe et al. datoista
temp=importdata('Data/Toluene/Keefe_toluene_2/tla105em.y'); % in L /molcm
vstart=temp(1)*10^2;
vend=temp(2)*10^2;
npt=temp(3);
v=linspace(vstart,vend,npt);
a=temp(4:end)*0.1*9339; % L/mol/cm=0.1 m^3/mol/m, Tolueenin tiheys 9339 mol/m^3
Tau(:,1)=effmean3(a,v,x,1450);


% testataan myös metanoli
[a v]=readyspec('Data/Methanol/Bertie/CH3OH/MTHEM93.y');
Tau(i,2)=effmean3(a*2.4719e+004,v,x,Temp(i));
% molar mass 32.04 g mol?1 density 792 kg/m^3
% => density = 2.4719e+004 mol/m^3


% Vesi
[a v]=readyspec('Data/Water/Bertie/WTEREM95.y');
Tau(i,3)=effmean3(a*5.5556e+004,v,x,Temp(i));
% molar mass 18.0153 g/mol density 1000 kg/m^3
% => density =5.5556e+004 mol/m^3

% Bentseeni Data/Benzene/Keefe_benzene/bzh6emf2.y
[a v]=readyspec('Data/Benzene/Keefe_benzene/bzh6emf2.y');
Tau(i,4)=effmean3(a*100*11161,v,x,Temp(i));
% Density 11161 mol/m^3 at T=300K ja P = 1 atm
%kappa0=-log(Tau).'/(x*1e-3).';

% Diesel
A=csvread('Data/Diesel/Sazhin2004.csv',1,0);
[C,IA,IB]=UNION(A(:,1),A(:,1));
A=A(IB,:);
l=A(:,1)*1e-6;
a=A(:,2)*2*pi/l;
a=A(:,2)*2*pi./l;

v=1./l;
Tau(i,5)=effmean3(a,v,x,Temp(i));
end
figure;
plot(Temp,Tau,'LineWidth',2,'LineSmoothing','on');
%semilogy(Temp,Tau,'LineWidth',2,'LineSmoothing','on');
xlabel('Blackbody temperature (K)','FontSize',14,'FontName',Font_Name);
ylabel('Absorption coefficient (1/m)','FontSize',14,'FontName',Font_Name);
legend({'Toluene','Methanol','Water','Benzene','Diesel'},'FontSize',14,'FontName',Font_Name);
%set(gca,'FontSize',14);
% print to pdf
        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[Paper_Width Paper_Height]);
        set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
        display(['Printing plot ',num2str(i),'...'])
        print(gcf,'-dpdf','../KAPPA_TEMP')
%exportfig(gcf,'KAPPA_VS_TEMP.png','Renderer','painters', 'width',14,'height',12 ,'fontsize',1.2,...
%            'Color','bw','Format','png','Resolution',900);
        
end
