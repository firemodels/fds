function Tau=transmission(a,v,S,T)
C1=1.1910428220e-16;
C2=1.438775225e-2; 
sigma=5.67040040*10^-8; 
if(nargin<3) T=1450; end

dv=abs(min(diff(v)))/4;
eta=1:dv:1.2e6;
%eta=linspace(min(v),max(v),length(v)*10);
ai=(interp1(v,a,eta,'linear',a(end))).';
ai(ai<0)=0;

I0=((C1*eta.^3)./(exp((C2*eta)/T)-1)).';

Qin=trapz(eta,I0)*pi;
Itot=zeros(size(S));
for i=1:length(S)
    I=I0.*exp(-ai*S(i));
    Itot(i)=trapz(eta,I)*pi;
end

Tau=Itot./Qin;
end