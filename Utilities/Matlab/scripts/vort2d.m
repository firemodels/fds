%-------------------------------------------------------------------------%
% FILE: vort2d.m
% AUTHOR: Max Gould
% LAST DATE EDITED: 26 June 2012
%-------------------------------------------------------------------------%
close all
clear all

% Open directory with input files
input_dir = '../../Verification/NS_Analytical_Solution/';

%-------------------------------------------------------------------------%
% Define Variables
%-------------------------------------------------------------------------%

% Define basic flow parameters
L = 0.3112;
U0 = 35.0;
Rc = 0.01556;
BigGamma = 0.0359157;

% Define mesh sizes to be analyzed
meshsize = [40 80 160 320];

% Set time step of DEVC output
DT = 1.1114*10.^(-4);

% Set the flow-through time
    % The time period required for the vortex to pass through the
    % entire computational domain and return to where it was initialized
FLOWTIME = int16(0.3112/(U0*DT));
    
% Set the number of time steps to be analyzed
    % Do not set below 3 (see 'Export Error for Plotting')
TIMESTEPS = 4;

% Allocate space for arrays
MatIn = cell(4);
DEVCNUM = zeros(1,4);
STEPSIZE = zeros(1,4);
AxisPlotLabels = cell(1,TIMESTEPS+1);
COUNT = zeros(4,TIMESTEPS);
RMS = zeros(4,TIMESTEPS);
    
% Define step sizes for each mesh
for m = 1:4
    DEVCNUM(m) = meshsize(m)/5;
    STEPSIZE(m) = 0.311200/meshsize(m);
end
    
%-------------------------------------------------------------------------%
% Import Data
%-------------------------------------------------------------------------%

for m = 1:4
    meshname = num2str(40*2.^(m-1));
    endrow = FLOWTIME*TIMESTEPS+3;
    endcol = DEVCNUM(m);
    MatIn{m}  = csvread([input_dir,'vort2d_',meshname,'_devc.csv'],2,0,...
                       [2,0,endrow,endcol]);
end
    
%-------------------------------------------------------------------------%
% Plot Data
%-------------------------------------------------------------------------%

% Open directory for plot output
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

% Define standard plotting parameters
plot_style

% Run for each mesh size
for m = 1:4
        
    % Create mesh-size string for variable file-name output
    meshname = num2str(40*2.^(m-1));

    %-----------------------------------------------------------------%
    % Plot U-Velocity at X=0 Data for a Range of Times
    % Calculate RMS Error Values
    %-----------------------------------------------------------------%
    
    % Define plot range
    Z0 = STEPSIZE(m)/2;
    Z  = -2*Rc+Z0:STEPSIZE(m):2*Rc-Z0;
    
    % Plot analytical solution
    AxisPlotExact = plot(Z,U0-BigGamma.*Z.*...
                         exp(-Z.^2./(2*Rc.^2))./Rc.^2,'--');
    
    % Define plot color as black
    set(AxisPlotExact,'Color',[0 0 0],'LineWidth',1.4);
    
    % Hold analytical solution plot to be included in future plot
    hold on
    
    % Plot FDS simulation values with
    % different colors for each time step
    for t = 1:TIMESTEPS+1
        
        % Set time steps to be integer multiples
        % of the flow-through time
        TIME = FLOWTIME*(t-1)+1;
        
        % Initialize error arrays
        COUNT(m,t) = 0.0;
        RMS(m,t) = 0.0;
        
        % Allocate space for array
        UArray = zeros(1,DEVCNUM(m));
        
        for k = 2:DEVCNUM(m)+1
            % Create array of u-velocity values
            UArray(k-1) = MatIn{m}(TIME,k);
            
            % Define actual z-coordinate (see UZEXACT)
            ZVAL = -2*Rc+Z0+(k-2)*STEPSIZE(m);
            
            % Calculate analytical u-velocity values
            AxisExact = U0-BigGamma.*ZVAL.*...
                            exp(-ZVAL.^2./(2*Rc.^2))./Rc.^2;
            
            % Calculate initial error values
            RMS(m,t) = RMS(m,t)+((UArray(k-1)-AxisExact)).^2;
            COUNT(m,t) = COUNT(m,t)+1.0;
        end
        
        % Calculate averaged error values
        RMS(m,t) = sqrt(RMS(m,t)/COUNT(m,t));
        
        % Create plot of u-velocity along a line at x=0
        AxisPlotFDS = plot(Z,UArray);
        
        % Define color functions
        CFRed = real(0.5*(tanh(6*t/TIMESTEPS-3)+1));
        CFGreen = real(-0.5*(tanh(7*t/TIMESTEPS-6)-1));
        
        % Define plot colors
        set(AxisPlotFDS,'Color',[CFRed CFGreen 0]);
    end
    
    % Plot all plots together (from hold on)
    hold off
    
    % Specify axes for above plots
    axis([-2*Rc+Z0 2*Rc-Z0 33.4 36.6]);

    % Set standard plotting parameters
    set(gca,'Units',Plot_Units)
    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Title_Font_Size)
    set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
    set(gcf,'DefaultLineLineWidth',Line_Width)
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
    
    % Label axes
    xlabel('Z-Coordinate (m)','FontSize',Title_Font_Size,...
           'Interpreter',Font_Interpreter);
    ylabel('U-Velocity (m/s)','FontSize',Title_Font_Size,...
           'Interpreter',Font_Interpreter);
    
    % Define legend-label strings
    AxisPlotLabels{1} = 'Analytical';
    for t = 2:TIMESTEPS+2
        TIME_INDEX = num2str(double((t-2)*FLOWTIME)*DT,'%6.3f');
        AxisPlotLabels{t} = ['FDS (t=',TIME_INDEX,' s)'];
    end
    
    % Create legend
    AxisPlotLeg = legend(AxisPlotLabels,'Location','NorthEast');
    set(AxisPlotLeg,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter);

    % Save plot to file
    print(gcf,'-dpdf',[plot_dir,'vort2d_',meshname,'_uzgraph']);

    %-----------------------------------------------------------------%
    % Plot U-Velocity at a Point as a Function of Time
    %-----------------------------------------------------------------%

    % Define plot range
    T = 0.0:DT:DT*double(FLOWTIME*TIMESTEPS+1);
        
    % Define actual coordinates of point for analytical solution
        % Must correspond to device coordinates
        % in numerical simulation (see UPGRAPH)
    XP = 0.0;
    ZP = -2*Rc+Z0;
        
    % Plot analytical solution (periodic to 6 loops)
    PointPlotExact = plot(T,U0-BigGamma.*ZP.*...
                            (exp(-(XP.^2+ZP.^2+(U0.^2).*(T.^2)-...
                             2.*U0.*XP.*T)./(2.*Rc.^2))+...
                             exp(-((XP+1*L).^2+ZP.^2+(U0.^2).*(T.^2)-...
                             2.*U0.*(XP+1*L).*T)./(2*Rc.^2))+...
                             exp(-((XP+2*L).^2+ZP.^2+(U0.^2).*(T.^2)-...
                             2.*U0.*(XP+2*L).*T)./(2*Rc.^2))+...
                             exp(-((XP+3*L).^2+ZP.^2+(U0.^2).*(T.^2)-...
                             2.*U0.*(XP+3*L).*T)./(2*Rc.^2))+...
                             exp(-((XP+4*L).^2+ZP.^2+(U0.^2).*(T.^2)-...
                             2.*U0.*(XP+4*L).*T)./(2*Rc.^2))+...
                             exp(-((XP+5*L).^2+ZP.^2+(U0.^2).*(T.^2)-...
                             2.*U0.*(XP+5*L).*T)./(2*Rc.^2))+...
                             exp(-((XP+6*L).^2+ZP.^2+(U0.^2).*(T.^2)-...
                             2.*U0.*(XP+6*L).*T)./(2*Rc.^2)))./Rc.^2);
        
    % Hold analytical solution plot to be included in future plot
    hold on
                
    % Create plot of FDS simulation values
        % Device specification must correspond to
        % point location for analytical solution
        % ---------------------------v
    PointPlotFDS = plot(T,MatIn{m}(:,2),'--');
        
    % Plot all plots together (from hold on)
    hold off
            
    % Define plot colors
    set(PointPlotExact,'Color',[0.3 0.3 0.3]);
    set(PointPlotFDS,'Color',[1 0 0]);

    % Specify axes for above plots
    axis([0.0 0.037 34.6 36.0]);
    
    % Set standard plotting parameters
    set(gca,'Units',Plot_Units)
    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Title_Font_Size)
    set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
    set(gcf,'DefaultLineLineWidth',Line_Width)
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
    
    % Label axes
    xlabel('Time (s)','FontSize',Title_Font_Size,...
           'Interpreter',Font_Interpreter);
    ylabel('U-Velocity (m/s)','FontSize',Title_Font_Size,...
           'Interpreter',Font_Interpreter);
    
    % Create legend
    PointPlotLeg = legend('Analytical','FDS','Location','SouthEast');
    set(PointPlotLeg,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter);

    % Save plot to file
    print(gcf,'-dpdf',[plot_dir,'vort2d_',meshname,'_upgraph']);
end

%-------------------------------------------------------------------------%
% Export Error Data for Plotting
%-------------------------------------------------------------------------%

% Define dx and dx^2 values
dx =[0.0077800 0.0038900 0.0019450 0.0009725];
dx2=[0.0000605 0.0000151 0.0000038 0.0000009];

% Open data file for reference from FDS_verification_script.m
RMSFile = fopen([input_dir,'vort2d_error.csv'],'wt');
    
% Write column headers
fprintf(RMSFile,'%s\n','dx, dxmod, dx^2, rms1, rms2, rms3');
    
% Write error data for all meshes
for m = 1:4
    fprintf(RMSFile,'%12.9f, %12.9f, %12.9f, %12.9f, %12.9f, %12.9f\n',...
            dx(m),120*dx(m),4500*dx2(m),RMS(m,2),RMS(m,3),RMS(m,4));
end


    
    
    
    
    
    
    