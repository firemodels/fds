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

plot_style
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X-.15,Plot_Y,Plot_Width,Plot_Height])

% Set global reaction rate parameters

global dTdt
global R0
global E
global A

dTdt = 5./60.;
R0 = 8314.3;
T_0 = 300.+273.;
delta_T = 80.;
r_0=2*dTdt/delta_T;
E=r_0*R0*T_0*T_0/(dTdt*0.4);
A=(E./(R0.*T_0.*T_0))*dTdt*exp(E./(R0.*T_0));

% Solve the ODE dY/dT=f(T,Y)

options = odeset('RelTol',1e-10,'AbsTol',1e-8);
[T,Y] = ode15s(@reaction_rate,[373 673],1.,options);

% Compute dY/dT

for i=2:length(Y)-1
    dYdT(i) = (Y(i+1)-Y(i-1))/(T(i+1)-T(i-1));
end
dYdT(1) = 0.;
dYdT(length(Y)) = 0.;

% Plot the solution

dTdt = 5./60.;
TC = T-273;
[AX,H1,H2] = plotyy(TC,Y,TC,-dYdT*dTdt*1000,'plot'); hold on

set(gca,'FontName',Font_Name)
set(AX(1),'FontName',Font_Name)
set(AX(2),'FontName',Font_Name)
xlabel('Temperature ($^\circ$C)','Interpreter','LaTeX','FontSize',Label_Font_Size)
set(get(AX(1),'Ylabel'),'String','Mass Fraction','Interpreter','LaTeX','FontSize',Label_Font_Size) 
set(get(AX(2),'Ylabel'),'String','Reaction Rate (s$^{-1}$) $\times$ 10$^3$','Interpreter','LaTeX','FontSize',Label_Font_Size) 
set(AX(1),'YLim',[0 1.1])
set(AX(2),'YLim',[0 2.2])
set(AX(1),'YTickMode','manual')
set(AX(2),'YTickMode','manual')
set(AX(1),'YTick',[0 0.2 0.4 0.6 0.8 1.0])
set(AX(2),'YTick',[0 0.4 0.8 1.2 1.6 2.0])
set(H1,'LineStyle','-')
set(H2,'LineStyle','-')

line([300. 300.],[0.00 0.96],'LineStyle','--','Color','black','LineWidth',1)
%line([300. 400.],[0.96 0.96],'LineStyle','--','Color','black','LineWidth',1)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../../Manuals/FDS_5_User_Guide/FIGURES/pyrolysis')




