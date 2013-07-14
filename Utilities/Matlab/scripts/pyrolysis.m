% McGrattan
% 8-10-2009
% pyrolysis.m
%
% This script solves the ODE
%
% dY/dt = -A*Y*exp(-E/(R*T)) where Y(0)=1
%
% Both Y and T are functions of time, t, but dT/dt is a specified constant.
% Typical TGA heating rates are 5, 10 and 20 degrees C per minute.


close all
clear all

% Set global reaction rate parameters

global dTdt
global R0
global E
global A
global residue

addpath('../../Verification/Pyrolysis')

show_fds = 1;

for i_plot=1:2
    
    close all
    
    plot_style
    set(gcf,'DefaultLineLineWidth',Line_Width)
    set(gca,'FontName',Font_Name)
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X-.15,Plot_Y,Plot_Width,Plot_Height])

    dTdt = 5./60.;
    R0 = 8314.3;

    if i_plot==1
        n_components = 3;
        T_p(1) = 100.+273.;
        T_p(2) = 300.+273.;
        Y_0(1)=0.00;
        Y_0(2)=1.00;
        Y_0(3)=0.00;
        delta_T(1) = 10.;
        delta_T(2) = -80.;
        r_p(2)     = 0.002;
        residue(1) = 0.0;
        residue(2) = 0.0; 
    else
        n_components = 3;
        T_p(1) = 100.+273.;
        T_p(2) = 300.+273.;
        Y_0(1)=0.1;
        Y_0(2)=0.9;
        Y_0(3)=0.0;
        delta_T(1) = 10.;
        delta_T(2) = 80.;
        residue(1) = 0.0;
        residue(2) = 0.2;
    end

    % Calculate A and E
    
    for i=1:n_components-1
        if delta_T(i)>0 
            r_p(i)=2*dTdt*(1.-residue(i))/delta_T(i);
        end
        E(i)=exp(1.)*r_p(i)*R0*T_p(i)*T_p(i)/dTdt;
        A(i)=exp(1.)*r_p(i)*exp(E(i)./(R0.*T_p(i)));
    end

    % Solve the ODE dY/dT=f(T,Y)

    options = odeset('RelTol',1e-10,'AbsTol',1e-8);
    [T,Y] = ode15s(@reaction_rate,[273 673],Y_0,options);

    % Compute dm/dT

    m = sum(Y,2);
    for i=2:length(Y(:,1))-1
        dmdT(i) = (m(i+1,1)-m(i-1,1))/(T(i+1,1)-T(i-1,1));
    end
    dmdT(1) = 0.;
    dmdT(length(Y(:,1))) = 0.;

    % Plot the solution

    TC = T-273;
    [AX,H1,H2] = plotyy(TC,m,TC,-dmdT*dTdt*1000,'plot'); hold on

    % Add the FDS solution

    if show_fds==1
        if i_plot==1
            if ~exist('pyrolysis_1_devc.csv')
                display('Error: File pyrolysis_1_devc.csv does not exist. Skipping case.')
                return
            end
            FDS = csvread('pyrolysis_1_devc.csv',2,0);
        else
            if ~exist('pyrolysis_2_devc.csv')
                display('Error: File pyrolysis_2_devc.csv does not exist. Skipping case.')
                return
            end
            FDS = csvread('pyrolysis_2_devc.csv',2,0);
        end
        h=plot(FDS(:,4),FDS(:,2),'b--');
        h=plot(FDS(:,4),FDS(:,3)*500,'g--'); % Multiply by 500 for right axis
    end

    % Plot attributes
    
    set(gca,'FontName',Font_Name)
    set(AX(1),'FontName',Font_Name)
    set(AX(2),'FontName',Font_Name)
    xlabel('Temperature (\circC)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    set(get(AX(1),'Ylabel'),'String','Mass Fraction','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    set(get(AX(2),'Ylabel'),'String','Reaction Rate (s^{-1}) \times 10^3','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    set(AX(1),'YLim',[0 1.1])
    set(AX(2),'YLim',[0 2.2])
    set(AX(1),'YTickMode','manual')
    set(AX(2),'YTickMode','manual')
    set(AX(1),'YTick',[0 0.2 0.4 0.6 0.8 1.0])
    set(AX(2),'YTick',[0 0.4 0.8 1.2 1.6 2.0])
    set(H1,'LineStyle','-')
    set(H2,'LineStyle','-')

    % Add vertical lines to indicate the temperature peaks
    
    for i=1:n_components-1
        axes(AX(2));
        line([T_p(i)-273 T_p(i)-273],[0.00 Y_0(i)*r_p(i)*(1-residue(i))*1000],'LineStyle','-','Color','black','LineWidth',1)
    end
    
    % add SVN if file is available
    
    if i_plot==1
        chid = 'pyrolysis_1';
    else
        chid = 'pyrolysis_2';
    end
    
    SVN_Filename = [chid,'_svn.txt'];
    if exist(SVN_Filename,'file')
        SVN = importdata(SVN_Filename);
        x_lim = get(gca,'XLim');
        y_lim = get(gca,'YLim');
        X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
        Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
        text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
            'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
    end
    
    % Create the PDF files
    
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width*1.1 Paper_Height]);
    set(gcf,'PaperPosition',[0 0 Paper_Width*1.1 Paper_Height]);
    if i_plot==1
        print(gcf,'-dpdf','../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/pyrolysis_1')
    else
        print(gcf,'-dpdf','../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/pyrolysis_2')
    end
    
end




