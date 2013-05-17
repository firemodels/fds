%-----------------------
% C Weinschenk
% Verification of Extinction Model
% Fuel: Methane
% 6/2012
%-----------------------

clear all
close all

addpath('../../Verification/Extinction/')

%-----------------------
% Initialize species
%-----------------------

%-------------
% vector key
% (1) = nitrogen
% (2) = methane
% (3) = oxygen
% (4) = carbon dioxide
% (5) = water vapor
%-------------
T_crit = 1600;          % critical flame temperature [K]
mass = 0.1;             % mass [kg]
volume = 0.001;         % volume [m^3]
R = 8.3145;             % gas constnat [J/mol K]

seed = 0.05:0.1:.95;
for ii = 1:10    
    for jj = 1:10
        fuel_seed(jj+10*(ii-1)) = seed(ii);
    end
end

seed2 = 300:175:2035;
for ii = 1:10    
    for jj = 1:10
        temp_seed(jj+10*(ii-1)) = seed2(jj);
    end
end

num_samp = 100;         % number of random input sets
ignite=zeros(num_samp,2);
extinct_o2=zeros(num_samp,2);
extinct_fuel=zeros(num_samp,2);

for i=1:num_samp
    phi_fuel(i) = fuel_seed(i); % initial mass fraction of fuel
    phi(i,:) = [0.77*(1-fuel_seed(i)); fuel_seed(i); 0.23*(1-fuel_seed(i)); 0.0; 0.0];  % initial mass fraction
    T0(i) = temp_seed(i); % initial temperature [K]
end

nu = [0; -1; -2; 1; 2]; % stoichiometric coefficients
y_MW = [28.0134; 16.042460; 31.9988; 44.0095; 18.01528]; % [g/mol]
y_hf = [0.0; -74873; 0.0; -393522; -241826]; % [J/mol]
pres_0 = 1.013253421185575e+05; % initial presure [Pa]
    
%-----------------------
% Determine max change in species
%-----------------------
for jj = 1:num_samp
    nu_mw_sum = 0;
    phi_sum = 0;

    for i=1:length(nu)
        if nu(i) < 0
            nu_mw_sum = nu_mw_sum + abs(nu(i))*y_MW(i);
            phi_sum = phi_sum + phi(jj,i);
        end
    end

    for i=1:length(nu)
        if nu(i) < 0
            phi_st(i) = abs(nu(i))*y_MW(i)/nu_mw_sum;
            phi_r(i) = phi(jj,i)/phi_sum;
        end
    end

    RR = phi_r./phi_st;
    ILR = 2;
    for i=3:length(nu)
        if nu(i) < 0
            if RR(i) < RR(ILR)
                ILR = i;
            end
        end
    end

    for i=1:length(nu)
        d_phi(jj,i) = nu(i)*y_MW(i)/(abs(nu(ILR))*y_MW(ILR))*phi(jj,ILR);
    end
    phi_new(jj,:) = phi(jj,:) + d_phi(jj,:);
    
    %-----------------------
    % Determine Extinction
    %-----------------------
        %---------
        %NOTE: coefficients were found from NIST Webbook
        %---------
        %reference temperature coeffs
        coeff_ref(:,1) = [28.98641;1.853978;-9.647459;16.63537;0.000117;-8.671914;0.0];
        coeff_ref(:,2) = [-0.703029;108.4773;-42.52157;5.862788;0.678565;-76.84376;-74.87310];
        coeff_ref(:,3) = [31.32234;-20.23531;57.86644;-36.50624;-0.007374;-8.903471;0.0];
        coeff_ref(:,4) = [24.99735;55.18696;-33.69137;7.948387;-0.136638;-403.6075;-393.5224];
        coeff_ref(:,5) = [30.09200;6.832514;6.793435;-2.534480;0.082139;-250.8810;-241.8264];   
        
        %initial temperature coeffs
        %nitrogen cp coeffs [J/mol K]
        if T0(jj) <=500
            coeff0(:,1) = [28.98641;1.853978;-9.647459;16.63537;0.000117;-8.671914;0.0]; 
        elseif 500 < T0(jj) <= 2000
            coeff0(:,1) = [19.50583;19.88705;-8.598535;1.369784;0.527601;-4.935202;0.0]; 
        else 
            coeff0(:,1) = [35.51872;1.128728;-0.196103;0.014662;-4.553760;-18.97091;0.0]; 
        end
        %methane cp coeffs [J/mol K]
        if T0(jj) <= 1300
            coeff0(:,2) = [-0.703029;108.4773;-42.52157;5.862788;0.678565;-76.84376;-74.87310];
        else
            coeff0(:,2) = [85.81217;11.26467;-2.114146;0.138190;-26.42221;-153.5327;-74.87310];
        end
        %oxygen cp coeffs [J/mol K]
        if T0(jj) <=700
            coeff0(:,3) = [31.32234;-20.23531;57.86644;-36.50624;-0.007374;-8.903471;0.0]; 
        elseif 700 < T0(jj) <= 2000
            coeff0(:,3) = [30.03235;8.772972;-3.988133;0.788313;-0.741599;-11.32468;0.0]; 
        else 
            coeff0(:,3) = [20.91111;10.72071;-2.020498;0.146449;9.245722;5.337651;0.0]; 
        end
        %carbon dooxide cp coeffs [J/mol K]
        if T0(jj) <= 1200
            coeff0(:,4) = [24.99735;55.18696;-33.69137;7.948387;-0.136638;-403.6075;-393.5224];    
        else    
            coeff0(:,4) = [58.16639;2.720074;-0.492289;0.038844;-6.447293;-425.9186;-393.5224];
        end    
        %water vapor cp coeffs [J/mol K] 
        if T0(jj) <= 1700
            coeff0(:,5) = [30.09200;6.832514;6.793435;-2.534480;0.082139;-250.8810;-241.8264];
        else
            coeff0(:,5) = [41.96246;8.622053;-1.499780;0.098199;-11.15764;-272.1797;-241.8264];
        end
        
        %critical temperature coeffs
        %nitrogen cp coeffs [J/mol K]
        if T_crit <=500
            coeff(:,1) = [28.98641;1.853978;-9.647459;16.63537;0.000117;-8.671914;0.0]; 
        elseif 500 < T_crit <= 2000
            coeff(:,1) = [19.50583;19.88705;-8.598535;1.369784;0.527601;-4.935202;0.0]; 
        else 
            coeff(:,1) = [35.51872;1.128728;-0.196103;0.014662;-4.553760;-18.97091;0.0]; 
        end
        %methane cp coeffs [J/mol K]
        if T_crit <= 1300
            coeff(:,2) = [-0.703029;108.4773;-42.52157;5.862788;0.678565;-76.84376;-74.87310];
        else
            coeff(:,2) = [85.81217;11.26467;-2.114146;0.138190;-26.42221;-153.5327;-74.87310];
        end
        %oxygen cp coeffs [J/mol K]
        if T_crit <=700
            coeff(:,3) = [31.32234;-20.23531;57.86644;-36.50624;-0.007374;-8.903471;0.0]; 
        elseif 700 < T_crit <= 2000
            coeff(:,3) = [30.03235;8.772972;-3.988133;0.788313;-0.741599;-11.32468;0.0]; 
        else 
            coeff(:,3) = [20.91111;10.72071;-2.020498;0.146449;9.245722;5.337651;0.0]; 
        end
        %carbon dooxide cp coeffs [J/mol K]
        if T_crit <= 1200
            coeff(:,4) = [24.99735;55.18696;-33.69137;7.948387;-0.136638;-403.6075;-393.5224];    
        else    
            coeff(:,4) = [58.16639;2.720074;-0.492289;0.038844;-6.447293;-425.9186;-393.5224];
        end    
        %water vapor cp coeffs [J/mol K] 
        if T_crit <= 1700
            coeff(:,5) = [30.09200;6.832514;6.793435;-2.534480;0.082139;-250.8810;-241.8264];
        else
            coeff(:,5) = [41.96246;8.622053;-1.499780;0.098199;-11.15764;-272.1797;-241.8264];
        end
        t=T_crit/1000;
        t0=T0(jj)/1000;
        for i=1:5
            del_h_init(i) = 1000*(coeff0(1,i)*t0 + (1/2)*coeff0(2,i)*t0^2 + (1/3)*coeff0(3,i)*t0^3 + (1/4)*coeff0(4,i)*t0^4 - coeff0(5,i)/t0 + coeff0(6,i) - coeff0(7,i));
            del_h(i) = 1000*(coeff(1,i)*t + (1/2)*coeff(2,i)*t^2 + (1/3)*coeff(3,i)*t^3 + (1/4)*coeff(4,i)*t^4 - coeff(5,i)/t + coeff(6,i) - coeff(7,i));
        end

        Q(jj) = 0;
        h = 0;
        h0 = 0;
        for i=1:length(y_hf)
            Q(jj) = Q(jj) - y_hf(i)*d_phi(jj,i)/y_MW(i);
            h0 = h0 +phi(jj,i).*del_h_init(i)/y_MW(i);
            h = h + phi(jj,i).*del_h(i)/y_MW(i);
        end           
    
    if h0 + Q(jj)  > h
        ignite(jj,:) = [T0(jj);phi(jj,3)];
    else
        extinct_o2(jj,:) = [T0(jj);phi(jj,3)];  
    end
end

%-----------------------
% Import FDS results
%-----------------------
epsilon = 1e-10;

if ~exist('extinction_devc.csv')
    display('Error: File extinction_devc.csv does not exist. Skipping case.')
    return
end

extinction_2=importdata('extinction_devc.csv');
extinct_2(:,:)=extinction_2.data; % data

for i=1:100
   hrr_ext(:,i) = extinct_2(:,4*i-2);
   temp_ext(:,i) = extinct_2(:,4*i-1)+273.15;
   o2_ext(:,i) = extinct_2(:,4*i);
   fu_ext(:,i) = extinct_2(:,4*i+1);
end

for i=1:100
    if sum(hrr_ext(:,i)) > epsilon
        fds_ignite(i,:) = [temp_ext(1,i);o2_ext(1,i)];
    else
        fds_ext_o2(i,:) = [temp_ext(1,i);o2_ext(1,i)];
    end
end

%-----------------------
% Simple Extinction Model
%-----------------------

simple_o2 = [0.15 0];
simple_temp = [273.15 1700];


%-----------------------
% Plotting
%-----------------------
figure(1)
h=plot(simple_temp,simple_o2,'k',ignite(:,1),ignite(:,2),'rs',fds_ignite(:,1),fds_ignite(:,2),'r+',extinct_o2(:,1),extinct_o2(:,2),'bo',fds_ext_o2(:,1),fds_ext_o2(:,2),'b*','LineWidth',0.5,'MarkerSize',4);
set(h([1]),'LineWidth',1)
set(h([2 4]),'MarkerSize',7)
axis([273.15 1900 0 0.23])
plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Scat_Plot_X,Scat_Plot_Y,Scat_Plot_Width,Scat_Plot_Height])
xlabel('Temperature (K)','Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size,'FontName',Font_Name)
ylabel('Mass Fraction Oxygen','Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size,'FontName',Font_Name)
legend('Simple Model','Expected Burning','FDS Burning','Expected Extinction','FDS Extinction','Location','NorthEast')

% add SVN if file is available

svn_file = 'extinction_svn.txt';

if exist(svn_file,'file')
    SVN = importdata(svn_file);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/extinction');
