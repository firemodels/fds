% webbook datoista laskettu effective mean

function kappa=effmean_webbook
x=(0:0.001:0.004).';
Tau=zeros(length(x),4);
N=20;

sigma=5.67040040*10^-8; 
qin=sigma*1450^4;

[a v]=abscoeffs('Data\N-heptane\142-82-5-IR.jdx');
Tau(:,1)=transmission(a,v,x,1450);
[a v]=abscoeffs('Data\N-hexane\110-54-3-IR.jdx');
Tau(:,2)=transmission(a,v,x,1450);
[a v]=abscoeffs('Data\Ethanol\64-17-5-IR(SOLUTION).jdx');
Tau(:,3)=transmission(a,v,x,1450);
[a v]=abscoeffs('Data\Methanol\67-56-1-IR.jdx');
Tau(:,4)=transmission(a,v,x,1450);
[a v]=abscoeffs('Data\dioxane\123-91-1-IR-2.jdx');
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

plot(x*1e3,Tau(:,1),'*');
hold on;
plot(x*1e3,Tau(:,2),'x');
plot(x*1e3,Tau(:,3),'o');
plot(x*1e3,Tau(:,4),'d');

%plot(x*1e3,exp(-kappa*x.'))
%plot(X,q/qin,'r');
plot(xx.'*1e3,I.','-k','LineWidth',2,'LineSmoothing','on');
legend();
exportfig(gcf,'Suo_webbook_effective_FDS.png','Renderer','painters', 'width',12,'height',12 ,'fontsize',1.2,...
            'Color','bw','Format','png','Resolution',600);
function mse=targetfun(kappa,expx,expy,qin)
        [q,X] = irad_twoflux(max(expx),qin,kappa,20,20);
        q=interp1(X,q,expx,'linear',qin)/qin;
        mse=mean((expy-q).^2);
end

end