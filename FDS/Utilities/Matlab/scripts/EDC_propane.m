%-----------------------
% C Weinschenk
% Verification of EDM (Species,Temp,Pres)
% Fuel: Propane
% 1/2012
%-----------------------

clear all
close all

%-----------------------
% Initialize species
%-----------------------

%-------------
% vector key
% (1) = nitrogen
% (2) = propane
% (3) = oxygen
% (4) = carbon dioxide
% (5) = water vapor
%-------------

temp_0 = 298.15;        % intial temperature [K]
volume = 0.001;         % volume [m^3]
rho=1.2043424;          % density [kg/m^3]
R = 8.3145;             % gas constnat [J/mol K]
mass = rho* volume;     % mass [kg]

yf0 = [0.720828; 0.060321; 0.218851; 0.0; 0.0];  % initial mass fraction

y_MW = [28.0134; 44.09562; 31.9988; 44.0095; 18.01528]; % [g/mol]

y_hf = [0.0; -103850; 0.0; -393522; -241826]; % [J/mol]

N0 = (1000*volume*rho*yf0)./(y_MW);  % initial moles

pres_0 = temp_0*sum(N0)*R/volume;      % initial pressure

%-----------------------
% Setup time vector for integration
%-----------------------
tspan=0:0.1:2;

%-----------------------
% Setup vector of moles for integrator
%-----------------------
y0=[N0(1) N0(2) N0(3) N0(4) N0(5)];

%-----------------------
% Pass information to integrator
%-----------------------
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[T,Y]=ode113(@EDC_spec_propane,tspan,y0, options);

Nf = Y(end,:);

%-----------------------
% Convert back to mass fraction
%-----------------------
yff(:,1)=(y_MW(1)*Y(:,1))/(1000*rho*volume);
yff(:,2)=(y_MW(2)*Y(:,2))/(1000*rho*volume);
yff(:,3)=(y_MW(3)*Y(:,3))/(1000*rho*volume);
yff(:,4)=(y_MW(4)*Y(:,4))/(1000*rho*volume);
yff(:,5)=(y_MW(5)*Y(:,5))/(1000*rho*volume);

%-----------------------
% Determine final temperature and pressure
%-----------------------
tol0 = 1e5;
tol = 0.1;
Tf_guess = [1000;2000;3000];    % final temperature guess [K]
count = 0;
while abs(tol0)>tol
    %---------
    %NOTE: coefficients were found from NIST Webbook
    %---------    
    
    %initial temperature coeffs
    coeff0(:,1) = [28.98641;1.853978;-9.647459;16.63537;0.000117;-8.671914;0.0];
    coeff0(:,2) = [1;1;1;1;1;1;1]; %placeholder for propane
    coeff0(:,3) = [31.32234;-20.23531;57.86644;-36.50624;-0.007374;-8.903471;0.0];
    coeff0(:,4) = [24.99735;55.18696;-33.69137;7.948387;-0.136638;-403.6075;-393.5224];
    coeff0(:,5) = [30.09200;6.832514;6.793435;-2.534480;0.082139;-250.8810;-241.8264];
    
    %nitrogen cp coeffs [J/mol K]
    if Tf_guess <=500
        coeff(:,1) = [28.98641;1.853978;-9.647459;16.63537;0.000117;-8.671914;0.0]; 
    elseif 500 < Tf_guess <= 2000
        coeff(:,1) = [19.50583;19.88705;-8.598535;1.369784;0.527601;-4.935202;0.0]; 
    else 
        coeff(:,1) = [35.51872;1.128728;-0.196103;0.014662;-4.553760;-18.97091;0.0]; 
    end
    %propane cp coeffs [J/mol K]
    if Tf_guess <= 1300
        coeff(:,2) = [1;1;1;1;1;1;1]; %placeholder for propane
    else
        coeff(:,2) = [1;1;1;1;1;1;1]; %placeholder for propane
    end
    %oxygen cp coeffs [J/mol K]
    if Tf_guess <=700
        coeff(:,3) = [31.32234;-20.23531;57.86644;-36.50624;-0.007374;-8.903471;0.0]; 
    elseif 700 < Tf_guess <= 2000
        coeff(:,3) = [30.03235;8.772972;-3.988133;0.788313;-0.741599;-11.32468;0.0]; 
    else 
        coeff(:,3) = [20.91111;10.72071;-2.020498;0.146449;9.245722;5.337651;0.0]; 
    end
    %carbon dooxide cp coeffs [J/mol K]
    if Tf_guess <= 1200
        coeff(:,4) = [24.99735;55.18696;-33.69137;7.948387;-0.136638;-403.6075;-393.5224];    
    else    
        coeff(:,4) = [58.16639;2.720074;-0.492289;0.038844;-6.447293;-425.9186;-393.5224];
    end    
    %water vapor cp coeffs [J/mol K] 
    if Tf_guess <= 1700
        coeff(:,5) = [30.09200;6.832514;6.793435;-2.534480;0.082139;-250.8810;-241.8264];
    else
        coeff(:,5) = [41.96246;8.622053;-1.499780;0.098199;-11.15764;-272.1797;-241.8264];
    end

    t=Tf_guess./1000;
    t0=temp_0/1000;
    for j=1:3
        for i=1:5
            del_h(i,j) = ((coeff(1,i)*t(j) + (1/2)*coeff(2,i)*t(j)^2 + (1/3)*coeff(3,i)*t(j)^3 + (1/4)*coeff(4,i)*t(j)^4 - coeff(5,i)/t(j) + coeff(6,i) - coeff(7,i))...
                -(coeff0(1,i)*t0 + (1/2)*coeff0(2,i)*t0^2 + (1/3)*coeff0(3,i)*t0^3 + (1/4)*coeff0(4,i)*t0^4 - coeff0(5,i)/t0 + coeff0(6,i) - coeff0(7,i)))...
                * 1000;
        end
    end
    
    %overwrite del_h(2,:) for propane
    coeff2(:,1) = [1.67536e-9;-5.46675e-6;4.38029e-3;2.7949;6.11728e2];
    for j=1:3
        del_h(2,j)=((1/5)*coeff2(1,1)*Tf_guess(j)^5+(1/4)*coeff2(2,1)*Tf_guess(j)^4+(1/3)*coeff2(3,1)*Tf_guess(j)^3+(1/2)*coeff2(4,1)*Tf_guess(j)^2+coeff2(5,1)*Tf_guess(j))...
            -((1/5)*coeff2(1,1)*temp_0^5+(1/4)*coeff2(2,1)*temp_0^4+(1/3)*coeff2(3,1)*temp_0^3+(1/2)*coeff2(4,1)*temp_0^2+coeff2(5,1)*temp_0)...
            .* (Tf_guess(j) - temp_0);
    end
   
    h_fc = Nf(4)*y_hf(4)+Nf(5)*y_hf(5)-N0(2)*y_hf(2);
    del_hN = Nf(1).*del_h(1,:)+Nf(4).*del_h(4,:)+Nf(5).*del_h(5,:);
    RT0 = sum(N0)*R*temp_0;
    RTf = -sum(Nf)*R*Tf_guess;

    Tf_solve = h_fc + del_hN + RT0 + RTf';
    tol0 = Tf_guess(1) - Tf_guess(3);
    
    if Tf_solve(2) < 0
        if Tf_solve(1) > 0
            mp = 0.5*(Tf_guess(1)+Tf_guess(2));
            Tf_guess = [Tf_guess(1);mp;Tf_guess(2)];
        else
            mp = 0.5*(Tf_guess(2)+Tf_guess(3));
            Tf_guess = [Tf_guess(2);mp;Tf_guess(3)];
        end
    else
        if Tf_solve(1) < 0
            mp = 0.5*(Tf_guess(1)+Tf_guess(2));
            Tf_guess = [Tf_guess(1);mp;Tf_guess(2)];
        else
            mp = 0.5*(Tf_guess(2)+Tf_guess(3));
            Tf_guess = [Tf_guess(2);mp;Tf_guess(3)];
        end
    end
 
    count = count+1;
end
Tf=Tf_guess(2);
Pf=(Tf_guess(2)*sum(Nf)*R)/volume;
dP = Pf-pres_0;
Tf = Tf - 273.15;

%-----------------------
% Create SS Vector of T and P
%-----------------------

tss=[1;1.25;1.5;1.75;2];
Tf=[Tf;Tf;Tf;Tf;Tf];
dP=[dP;dP;dP;dP;dP];

%------------------------------------------------
% Write Expected Data CSV File
%------------------------------------------------
yff(:,1) = tspan;

header1 = {'Time','O2','C3H8','CO2','H2O'};
filename1 = '../../Verification/Species/reactionrate_EDC_flim_1step_C3H8_spec.csv';
fid = fopen(filename1,'wt');
fprintf(fid,'%s, %s, %s, %s, %s\n',header1{:});
for j=1:length(tspan)
    fprintf(fid,'%f, %f, %f, %f, %f\n',yff(j,:));
end
fclose(fid);

header1 = {'Time','TEMP','PRES'};
filename1 = '../../Verification/Species/reactionrate_EDC_flim_1step_C3H8_temppres.csv';
fid = fopen(filename1,'wt');
fprintf(fid,'%s, %s, %s\n',header1{:});
for j=1:length(tss)
    fprintf(fid,'%f, %f, %f\n',tss(j),Tf(j),dP(j));
end
fclose(fid);


