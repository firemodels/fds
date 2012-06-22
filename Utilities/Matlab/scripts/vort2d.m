%-------------------------------------------------------------------------%
% FILE: vort2d.m
% AUTHOR: Max Gould
% LAST DATE EDITED: 21 June 2012
%-------------------------------------------------------------------------%
close all
clear all

% Open directory with input files
input_dir = '../../Verification/NS_Analytical_Solution/';

%-------------------------------------------------------------------------%
% Define Variables
%-------------------------------------------------------------------------%

% Define basic flow parameters
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
TIMESTEPS = 4;
    
% Allocate space for STEP, STEPSIZE, and DEVCNUM
STEP = zeros(1,4);
STEPSIZE = zeros(1,4);
DEVCNUM = zeros(1,4);
    
% Define step sizes and total device count for each mesh
for m = 1:4
    STEP(m) = meshsize(m)/5+1;
    STEPSIZE(m) = 0.311200/meshsize(m);
    DEVCNUM(m) = STEP(m)*STEP(m);
end
    
% Allocate space for arrays
MatU = cell(4);
COUNT = zeros(4,TIMESTEPS);
ERROR = zeros(4,TIMESTEPS);
    
%-------------------------------------------------------------------------%
% Import Data
%-------------------------------------------------------------------------%
disp('Importing Data...');

for m = 1:4
    meshname = num2str(40*2.^(m-1));
    endrow = FLOWTIME*TIMESTEPS+3;
    endcol = STEP(m)*STEP(m);
    MatU{m}  = csvread([input_dir,'vort2d_',meshname,'_devc.csv'],2,0,...
                       [2,0,endrow,endcol]);
end
    
%-------------------------------------------------------------------------%
% Plot Data
%-------------------------------------------------------------------------%
disp('Plotting...');

% Open directory for plot output
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

% Define standard plotting parameters
plot_style

% Begin figure
figure

% Set standard plotting parameters
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);

% Run for each mesh size
for m = 1:4
        
    % Create mesh-size string for variable file-name output
    meshname = num2str(40*2.^(m-1));
        
    % Define plot range for analytical solutions
    Z0 = -STEPSIZE(m)/2;
    X  = -2*Rc:STEPSIZE(m):2*Rc;
    Z  = -2*Rc+Z0:STEPSIZE(m):2*Rc+Z0;

    %-----------------------------------------------------------------%
    % Construct U-Velocity Matrices for Each Time Step
    %-----------------------------------------------------------------%
        
    for t = 1:TIMESTEPS+1
            
        % Set time steps to be integer multiples
        % of the flow-through time
        TIME = FLOWTIME*(t-1)+1;
            
        for i = 1:STEP(m)
                
            % Specify devices corresponding to current
            % i-coordinate from .csv file
            JBEGIN = (i-1)*STEP(m)+2;
            JEND   =  i   *STEP(m)+1;
                
            for j = JBEGIN:JEND
                    
                % Create array of u-velocity values as
                % measured by devices from simulation
                UARRAY(STEP(m)-i+1,j-JBEGIN+1,t) = MatU{m}(TIME,j);
            end
        end
    end
        
    %-----------------------------------------------------------------%
    % Produce U-Velocity Contour Plots
    %-----------------------------------------------------------------%
        
    for t = 1:TIMESTEPS+1
            
        % Create time-step string for variable file-name output
        NUM = num2str(t);
            
        % Define file name for output
        filename = ['vort2d_',meshname,'_ugraph',NUM];
            
        % Create time-specific array of u-velocity values
        UARRAYPART = UARRAY(:,:,t);
            
        % Create contour plot of simulated u-velocities
        contour(UARRAYPART,50);
            
        % Save contour plot to file
        print(gcf,'-dpdf',[plot_dir,filename]);
    end
        
    %-----------------------------------------------------------------%
    % Plot U-Velocity at Z=0 Data for a Range of Times
    %-----------------------------------------------------------------%
            
    % Plot analytical solution
    UXEXACT = plot(X,U0,'--');
        
    % Define plot color as black
    set(UXEXACT,'Color',[0 0 0],'LineWidth',1.4);
        
    % Hold analytical solution plot to be included in future plot
    hold on
        
    % Plot FDS simulation values with
    % different colors for each time step
    for t = 1:TIMESTEPS+1
            
        % Create arrays of u-velocity values
        % above and below the Z=0 line
            % As a consequence of using multiple mesh sizes,
            % an average must be done to acquire u-velocity
            % values from the same point in space from each mesh
        UXLOW = UARRAY((STEP(m)+1)/2-1,:,t);
        UXHIGH = UARRAY((STEP(m)+1)/2,:,t);
            
        UXARRAY = (UXLOW+UXHIGH)./2;
            
        % Create plot of u-velocity along a line at z=0
        UXGRAPH = plot(X,UXARRAY);
            
        % Define color functions
        CFRed = real(0.5*(tanh(6*t/TIMESTEPS-3)+1));
        CFGreen = real(-0.5*(tanh(7*t/TIMESTEPS-6)-1));
            
        % Define plot colors
        set(UXGRAPH,'Color',[CFRed CFGreen 0]);
    end
        
    % Plot all plots together (from hold on)
    hold off

    % Specify axes for above plots
    axis([-0.031 0.031 34.5 35.5]);

    % Save plot to file
    print(gcf,'-dpdf',[plot_dir,'vort2d_',meshname,'_uxgraph']);
        
    %-----------------------------------------------------------------%
    % Plot U-Velocity at X=0 Data for a Range of Times
    %-----------------------------------------------------------------%
        
    % Plot analytical solution
    UZEXACT = plot(Z,U0-BigGamma.*Z.*...
                   exp(-Z.^2./(2*Rc.^2))./Rc.^2,'--');
         
    % Define plot color as black
    set(UZEXACT,'Color',[0 0 0],'LineWidth',1.4);
        
    % Hold analytical solution plot to be included in future plot
    hold on
        
    % Plot FDS simulation values with
    % different colors for each time step
    for t = 1:TIMESTEPS+1
            
        % Create time-specific array of u-velocity values
        UZARRAY = flipdim(UARRAY(:,(STEP(m)+1)/2,t),1);

        % Create plot of u-velocity along a line at x=0
        UZGRAPH = plot(Z,UZARRAY);
            
        % Define color functions
        CFRed = real(0.5*(tanh(6*t/TIMESTEPS-3)+1));
        CFGreen = real(-0.5*(tanh(7*t/TIMESTEPS-6)-1));
            
        % Define plot colors
        set(UZGRAPH,'Color',[CFRed CFGreen 0]);
    end

    % Plot all plots together (from hold on)
    hold off

    % Specify axes for above plots
    axis([-0.031 0.031 33.4 36.6]);

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
    XP = -2*Rc;
    ZP = -2*Rc+Z0;
        
    % Define periodicity of vortex
    X0 = 0.3112;
        
    % Plot analytical solution (periodic to 6 loops)
    UPEXACT = plot(T,U0-BigGamma.*ZP.*...
                     (exp(-(XP.^2+ZP.^2+(U0.^2).*(T.^2)-...
                      2.*U0.*XP.*T)./(2.*Rc.^2))+...
                      exp(-((XP+1*X0).^2+ZP.^2+(U0.^2).*(T.^2)-...
                      2.*U0.*(XP+1*X0).*T)./(2*Rc.^2))+...
                      exp(-((XP+2*X0).^2+ZP.^2+(U0.^2).*(T.^2)-...
                      2.*U0.*(XP+2*X0).*T)./(2*Rc.^2))+...
                      exp(-((XP+3*X0).^2+ZP.^2+(U0.^2).*(T.^2)-...
                      2.*U0.*(XP+3*X0).*T)./(2*Rc.^2))+...
                      exp(-((XP+4*X0).^2+ZP.^2+(U0.^2).*(T.^2)-...
                      2.*U0.*(XP+4*X0).*T)./(2*Rc.^2))+...
                      exp(-((XP+5*X0).^2+ZP.^2+(U0.^2).*(T.^2)-...
                      2.*U0.*(XP+5*X0).*T)./(2*Rc.^2))+...
                      exp(-((XP+6*X0).^2+ZP.^2+(U0.^2).*(T.^2)-...
                      2.*U0.*(XP+6*X0).*T)./(2*Rc.^2)))./...
                     Rc.^2);
        
    % Hold analytical solution plot to be included in future plot
    hold on
                
    % Create plot of FDS simulation values
        % Device specification must correspond to
        % oint location for analytical solution
        % ---------------------v
    UPGRAPH = plot(T,MatU{m}(:,2),'--');
        
    % Plot all plots together (from hold on)
    hold off
            
    % Define plot colors
    set(UPEXACT,'Color',[0.3 0.3 0.3]);
    set(UPGRAPH,'Color',[1 0 0]);

    % Specify axes for above plots
    axis([0.0 0.045 34.6 35.8]);

    % Save plot to file
    print(gcf,'-dpdf',[plot_dir,'vort2d_',meshname,'_upgraph']);
    
    %-----------------------------------------------------------------%
    % Calculate Error Values for X=0
    %-----------------------------------------------------------------%

    for t = 1:TIMESTEPS
            
        % Create time-specific array of u-velocity values
        UZARRAY = flipdim(UARRAY(:,(STEP(m)+1)/2,t+1),1);
            
        % Initialize error arrays
        COUNT(m,t) = 0.0;
        ERROR(m,t) = 0.0;
            
        for j = 1:STEP(m)
                
            % Define actual z-coordinate (see UZEXACT)
            ZVAL = -2*Rc+Z0+(j-1)*STEPSIZE(m);
                
            % Calculate analytical u-velocity values
            UZEXACT = U0-BigGamma.*ZVAL.*...
                          exp(-ZVAL.^2./(2*Rc.^2))./Rc.^2;
                
            % Calculate initial error values
            ERROR(m,t) = ERROR(m,t)+abs((UZARRAY(j)-UZEXACT)/UZEXACT);
            COUNT(m,t) = COUNT(m,t)+1.0;
        end
            
        % Calculate averaged error values
        ERROR(m,t) = ERROR(m,t)/COUNT(m,t);
    end
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
fprintf(RMSFile,'%s\n','dx, dxmod, dx^2, error1, error2, error3');
    
% Write error data for all meshes
for m = 1:4
    fprintf(RMSFile,'%12.8f, %12.8f, %12.8f, %12.8f, %12.8f, %12.8f\n',...
            dx(m),8*dx(m),120*dx2(m),ERROR(m,1),ERROR(m,2),ERROR(m,3));
end

disp('Done');
    
    
    
    
    
    
    