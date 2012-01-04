
function kappa=effmean_Suo_Anttila

addpath Two_flux

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


foo=csvread('Suo-anttila.csv',1);
x=foo(:,1)*1e-3;
Tau=zeros(length(x),9);
Tau(:,1:4)=foo(:,2:5);
N=20;

sigma=5.67040040*10^-8; 

qin=sigma*1450^4;

temp=csvread('abscoeffs.csv',2);

v=temp(:,1);
for i=2:size(temp,2)
    a=temp(:,i);
    %ind=a>0;
    Tau(:,i+3)=transmission(a,v,x,1450);
end

Tau(:,9)=transmission(a,v,x,1450);
for i=1:size(Tau,2),
    kappa0(i)=fminbnd(@(K) max((Tau(:,i)-exp(-K*x)).^2),0,1e4);
end
kappa0=kappa0.';

T=20;
kappa=zeros(size(kappa0));
I=zeros(length(kappa0),N);
xx=zeros(length(kappa0),N);
for i=1:length(kappa0),
    kappa(i)=fminbnd(@(X) targetfun(X,x,Tau(:,i),qin),0,2*kappa0(i));
    [q X] = irad_twoflux(max(x),qin,kappa(i),T,N);
    I(i,:)=q/qin;
    xx(i,:)=X;
end

I0=qin;



%[q,X] = irad_twoflux(max(x)*1e-3,qin,kappa,T,N);
figure;
plot(x*1e3,Tau(:,1),'*');
hold on;
plot(x*1e3,Tau(:,2),'x');
plot(x*1e3,Tau(:,3),'o');
plot(x*1e3,Tau(:,4),'d');
plot(x*1e3,Tau(:,5),'+');
plot(xx(1:5,:).'*1e3,I(1:5,:).','-k','LineWidth',2,'LineSmoothing','on');
legend({'Jp-8','EthTol','Ethanol','Toluene(Suo-Anttila et al.)','Toluene (Keefe et al.)'},'FontName',Font_Name);
% print to pdf
        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[Paper_Width Paper_Height]);
        set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
        display(['Printing plot ',num2str(i),'...'])
        %print(gcf,'-dpdf','../Suo_anttila_effective_FDS1')

exportfig(gcf,'../Suo_anttila_effective_FDS1.pdf','Renderer','painters', 'width',6.5,'height',5 ,'fontsize',1.2,...
            'Color','CMYK','Format','pdf','Resolution',900);
        figure;
plot(x*1e3,Tau(:,6),'s');
hold on;
plot(x*1e3,Tau(:,7),'^');
plot(x*1e3,Tau(:,8),'p');
plot(x*1e3,Tau(:,9),'v');
%plot(x*1e3,exp(-kappa*x.'))
%plot(X,q/qin,'r');
plot(xx(6:9,:).'*1e3,I(6:9,:).','-k','LineWidth',2,'LineSmoothing','on');
axis
legend('SouthEast',{'Methanol','Water','Benzene','Diesel'},'FontName',Font_Name);
% print to pdf
        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[Paper_Width Paper_Height]);
        set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
        display(['Printing plot ',num2str(i),'...'])
        %print(gcf,'-dpdf','../Suo_anttila_effective_FDS2','Resolution',900)

exportfig(gcf,'../Suo_anttila_effective_FDS2.pdf','Renderer','painters', 'width',6.5,'height',5 ,'fontsize',1.2,...
            'Color','CMYK','Format','pdf','Resolution',900);

function mse=targetfun(kappa,expx,expy,qin)
        [q,X] = irad_twoflux(max(expx),qin,kappa,20,20);
        q=interp1(X,q,expx,'linear',qin)/qin;
        mse=norm(expy-q);
end

end