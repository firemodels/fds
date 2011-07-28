
function [planck Temps]=planckmean(data,Temps)
if(isstruct(data))
    if(isfield(data,'type'))
        if(strcmp(data.type,'jcamp_struct'))
            D=data;
        else
            error('Wrong type struct');
        end
    else
        error('Wrong type struct');
    end
else
    if(ischar(data))
        D=readjcamp(data);
    else
        error('First argument needs to be jcamp_struct or filename');
    end
end


% 1/cm=1/(m*10^-2)=1/m*10^2
v=D.x*10^2;  % wave number 1/cm->1/m


sigma=5.67040040*10^-8; 
S=D.path_length*1e-2; % path length (cm=>m)
foo=regexp(D.state,'\w*','match');
state=foo{1};
switch(upper(state))
    case 'SOLUTION'
        c=str2double(foo(2))*10; % wt/vol% (g/100ml => 10 kg/m^3)
        aunits='$m^3 \cdot kg^{-1} \cdot m^{-1}$';
    case 'GAS'
        if(length(foo)>1)
            c=str2double(foo(2))*133.322; %mmHG -> Pa
            aunits='$(Pa\cdot m)^{-1}$';
        else
            c=1;
            aunits='$ m^{-1}$';
        end
    case 'LIQUID'
        c=1;
        aunits='$ m^{-1}$';
    otherwise
        error(strcat('Unknown state: ',foo(1)));
end

% Radiation constants
C1=1.1910428220e-16;
C2=1.438775225e-2; 


switch(upper(D.yunits))
    case 'TRANSMITTANCE'
%Absorbance
        T=D.y; % Transmittance
        %T(T<=0)=eps;
        T(T<=0)=realmin;
        
        A=-log(T);
    case 'ABSORBANCE'
        A=D.y;
end
%Spectral absorption coeff.
a=A./(S*c);



% figure;
% plot(10^-2*v,a);
% xlabel('Wave number $1/cm$','Interpreter','Latex');
% ylabel(['Absorption coefficient ' aunits ],'Interpreter','Latex');
% hold on
% ax1=gca;
% set(ax1,'XDir','reverse');
% ax2 = axes('Position',get(ax1,'Position'),...
%            'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XtickLabel',{},...
%            'XDir','reverse');
% Iv=(C1*v.^3)./(exp((C2*v)/300)-1);
% line(v,Iv,'Parent',ax2,'Color','r');
% ylabel(ax2,'Blackbody intensity $1/cm$','Interpreter','Latex');
% Iv=(C1*v.^3)./(exp((C2*v)/600)-1);
% line(v,Iv,'Parent',ax2,'Color','k');
% Iv=(C1*v.^3)./(exp((C2*v)/800)-1);
% line(v,Iv,'Parent',ax2,'Color','g');
% legend('300','600','800');
% 



k=1;
planck=zeros(size(Temps));
eta=linspace(min(v),max(v),length(v)*2);
ai=interp1(v,a,eta);
for Temp=Temps,
    Iv=(C1*eta.^3)./(exp((C2*eta)/Temp)-1);
    Ib=pi/(sigma*Temp^4);
    tmp=ai.*Iv;
    planck(k)=Ib*trapz(eta,tmp);
    k=k+1;
end
% figure;
% plot(300:1500,planck);
% xlabel('Temperature K','Interpreter','Latex');
% ylabel(['Plank mean absorption coefficient ',aunits],'Interpreter','Latex');
% title(D.title);
% Temps=300:1500;
% figure;
% plot(300:1500,planck*800);
% xlabel('Temperature K','Interpreter','Latex');
% ylabel(['Plank mean absorption coefficient ','$1/m$'],'Interpreter','Latex');
% title(D.title);
% Temps=300:1500;
end