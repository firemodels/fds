% McDermott
% 9-22-08
% energy_decay.m

close all
clear all

N = 32;
L = 9*(2*pi)/100;

M = dlmread(['cbc_',num2str(N),'_devc.csv'],',');
t = M(:,1);
ke = M(:,2)/(L^3);

figure(1)
H(3)=plot(t,ke,'k-'); hold on
xlabel('time (s)')
ylabel('kinetic energy (m^2/s^2)')

M = dlmread('cbc_32_nu0_devc.csv',',');
t = M(:,1);
ke = M(:,2)/(L^3);
H(1)=plot(t,ke,'b-'); hold on

M = dlmread('cbc_32_molnu_devc.csv',',');
t = M(:,1);
ke = M(:,2)/(L^3);
H(2)=plot(t,ke,'r-'); hold on


% Gather the Comte-Bellot/Corrsin data for comparison
% Open and read the grid information
infile = fopen('cbcdata.txt','r');

% Initialize the field
kCBC = zeros(20,1);
ECBC42 = zeros(20,1);
ECBC98 = zeros(20,1);
ECBC171 = zeros(20,1);

% Read the data
for j = 1:20   
   k(j) = fscanf(infile,'%f',1)*1e2;
   E_42(j) = fscanf(infile,'%f',1)/1e6;
   E_98(j) = fscanf(infile,'%f',1)/1e6;
   E_171(j) = fscanf(infile,'%f',1)/1e6;
end

% Close the CBC dat file
fclose(infile);

% Apply transfer function
k0 = 2*pi/L;
kmax = 0.5*N*k0;
delta = pi/kmax;
%G = sin(1/2*k.*delta)./(1/2*k.*delta);  % box filter
G = ones(1,length(k));
for j = 1:length(k)
  if k(j)>kmax
     tmp = 3*kmax/k(j) - 2;
     if tmp > 0
        G(j) = sqrt(tmp); % RJM filter
        %G(j) = 0; % spectral cutoff
     else
        G(j) = 0;
     end
  end
end

E_42   = G.*G.*E_42;
E_98   = G.*G.*E_98;
E_171  = G.*G.*E_171;

% % Use this to check that spectra are correct
% figure(2)
% loglog(k,E_42);hold on
% loglog(k,E_98)
% loglog(k,E_171)
% xlabel('\kappa (1/m)')
% ylabel('E(\kappa) (m^3/s^2)')

% Now integrate
Ebar_42 = 0;
Ebar_98 = 0;
Ebar_171 = 0;
for j = 1:length(k)-1
   dk = k(j+1)-k(j);
   %Ebar_42 = Ebar_42 + sqrt(E_42(j)*E_42(j+1))*dk;
   %Ebar_98 = Ebar_98 + sqrt(E_98(j)*E_98(j+1))*dk;
   %Ebar_171 = Ebar_171 + sqrt(E_171(j)*E_171(j+1))*dk;
   Ebar_42 = Ebar_42 + 0.5*(E_42(j)+E_42(j+1))*dk;
   Ebar_98 = Ebar_98 + 0.5*(E_98(j)+E_98(j+1))*dk;
   Ebar_171 = Ebar_171 + 0.5*(E_171(j)+E_171(j+1))*dk;
end
figure(1)
H(4)=plot(0.0,Ebar_42,'ro');
plot(0.28,Ebar_98,'ro')
plot(0.66,Ebar_171,'ro')
legend(H,'FDS no visc','FDS mol visc','FDS Smag','filtered CBC data','Location','East')

