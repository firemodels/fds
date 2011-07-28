function [amean amean0]=effmean3(a,v,path_length,varargin)
% Radiation constants
C1=1.1910428220e-16;
C2=1.438775225e-2; 
sigma=5.67040040*10^-8; 

optargin = size(varargin,2);
switch(optargin)
    case 0
        Temp=1450;
    case 1
        Temp=varargin{1};
    case 2
        Vin=varargin{1};
        Iin=varargin{2};   
    otherwise
       error('Wrong number of arguments'); 
end
        



% incident radiation
%dv=abs(min(diff(v)))/4;
%eta=1:dv:1.2e6;
eta=v;
%eta=linspace(min(v),max(v),length(v)*10);
ai=(interp1(v,a,eta,'linear',a(end))).';
ai(ai<0)=0;
if(optargin<=1)
    I0=((C1*eta.^3)./(exp((C2*eta)/Temp)-1)).';
else
    I0=interp1(Vin,Iin,eta,'linear',0).';
end
for j=1:length(path_length),
S=linspace(0,path_length(j),100);
Itot=zeros(size(S));
for i=1:length(S)
    I=I0.*exp(-ai*S(i));
    Itot(i)=trapz(eta,I)*pi;
end
%figure;
%plot(S,Itot/Itot(1));
%hold
I0=Itot(1);
T=Itot/I0;

[ncol nrow]=size(S);
if(ncol>1) S=S.'; end
[ncol nrow]=size(T);
if(ncol>1) T=T.'; end

amean0(j)=-log(T)/S; % Least squares fit.

% fit using the two-flux radiation model
[amean(j),fval,exitflag,output]=fminbnd(@(X) targetfun(X,S,Itot,I0),0,2*amean0(j))

%[q,X] = irad_twoflux(max(S),I0,amean,20,100);
%plot(S,exp(-amean0*S),'r');
%plot(X,q/I0,'g');
end

 function mse=targetfun(kappa,expx,expy,qin)
        [q,X] = irad_twoflux(max(expx),qin,kappa,20,1000);
        q=interp1(X,q,expx,'linear');
        mse=sum((expy-q).^2);
  end
end
