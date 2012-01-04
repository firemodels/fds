% make a semilogy containing the spectra used in determining the absorption
% coefficients

function plotpectra

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


subplot(3,2,1);
% Lasketaan myös Keefe et al. datoista
temp=importdata('Data/toluene/Keefe_toluene_2/tla105em.y'); % in L /molcm
vstart=temp(1)*10^2;
vend=temp(2)*10^2;
npt=temp(3);
v=linspace(vstart,vend,npt);
a=temp(4:end)*0.1*9339; % L/mol/cm=0.1 m^3/mol/m, Tolueenin tiheys 9339 mol/m^3
semilogy(v,a);
axis([0 4e5 0 1e5]);
title('Toluene','FontName',Font_Name)
xlabel('\nu (m^{-1})','FontName',Font_Name);
ylabel('a_{\lambda} (m^{-1})','FontName',Font_Name);


subplot(3,2,2);
% testataan myös metanoli
[a v]=readyspec('Data/Methanol/Bertie/CH3OH/MTHEM93.y');
semilogy(v,a*2.4719e+004);
title('Methanol','FontName',Font_Name)
axis([0 7e5 0 4e3]);
xlabel('\nu (m^{-1})','FontName',Font_Name);
ylabel('a_{\lambda} (m^{-1})','FontName',Font_Name);
% molar mass 32.04 g mol?1 density 792 kg/m^3
% => density = 2.4719e+004 mol/m^3

subplot(3,2,3);
% Vesi
[a v]=readyspec('Data/Water/Bertie/WTEREM95.y');
semilogy(v,a*5.5556e+004);
title('Water','FontName',Font_Name);
axis([0 1.5e6 0 5e5]);
xlabel('\nu (m^{-1})','FontName',Font_Name);
ylabel('a_{\lambda} (m^{-1})','FontName',Font_Name);
% molar mass 18.0153 g/mol density 1000 kg/m^3
% => density =5.5556e+004 mol/m^3

subplot(3,2,4);
% Bentseeni Data/Benzene/Keefe_benzene/bzh6emf2.y
[a v]=readyspec('Data/Benzene/Keefe_benzene/bzh6emf2.y');
semilogy(v,a*11161);
title('Benzene','FontName',Font_Name);
axis([0 6e5 0 4e3]);
xlabel('\nu (m^{-1})','FontName',Font_Name);
ylabel('a_{\lambda} (m^{-1})','FontName',Font_Name);
% Density 11161 mol/m^3 at T=300K ja P = 1 atm
%kappa0=-log(Tau).'/(x*1e-3).';

subplot(3,2,5);
% Diesel
A=csvread('Data/Diesel/Sazhin2004.csv',1,0);
[C,IA,IB]=UNION(A(:,1),A(:,1));
A=A(IB,:);
l=A(:,1)*1e-6;
a=A(:,2)*2*pi./l;
v=1./l;
semilogy(v,a);
title('Diesel','FontName',Font_Name);
xlabel('\nu (m^{-1})','FontName',Font_Name);
ylabel('a_{\lambda} (m^{-1})','FontName',Font_Name);

% print to pdf
        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[Paper_Width Paper_Height]);
        set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
        display(['Printing plot ',num2str(i),'...'])
        print(gcf,'-dpdf',['Plots','Absorption_spectra'])
        
end
