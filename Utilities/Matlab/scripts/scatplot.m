% R. McDermott
% 7-06-2009
% scatplot.m
%
% Generate scatter plots.  Must first run dataplot.m to generate
% saved_data and drange.
%
% [] = scatplot(saved_data,drange,[qrange])
%
% Arguments:
%    saved_data - cell array of packed data, may be obtained by running
%    the dataplot.m function
%
%    drange - obtained from (or input to) function dataplot.m
%
%    [optional] qrange - specify the range of lines to read from the
%    scatterplot configuration file
%
% Dependencies:
%    ../scripts/define_qrow_variables.m
%
% Example:
%    From the Matlab/functions/ directory, type
%    >> scatplot(saved_data,drange,3:4)

function [] = scatplot(saved_data,drange,varargin)

% unpack data
Save_Quantity         = saved_data{:,1};
Save_Group_Style      = saved_data{:,2};
Save_Fill_Color       = saved_data{:,3};
Save_Group_Key_Label  = saved_data{:,4};
Save_Measured_Metric  = saved_data{:,5};
Save_Predicted_Metric = saved_data{:,6};

qfil = [pwd,'/scatterplot_config_matlab.csv'];

if nargin==3
    qrange = varargin{1};
else
    qrange = [2:100];
end

paper_width  = 6.0; % inches
paper_height = 6.0; % inches

Q = importdata(qfil);
H = textscan(Q{1},'%q','delimiter',',');
headers = H{:}'; clear H

for j=qrange
    if j>length(Q); break; end
    
    define_qrow_variables
    
    clear Measured_Metric
    clear Predicted_Metric
    
    k = 0;
    for i=drange
        if i>length(Save_Quantity); break; end
        if strcmp(Save_Quantity(i),Scatter_Plot_Title)
            k = k+1;
            Measured_Metric(k)  = Save_Measured_Metric(i);
            Predicted_Metric(k) = Save_Predicted_Metric(i);
            Group_Key_Label(k)  = Save_Group_Key_Label(i);
            K(k) = plot(Measured_Metric(k),Predicted_Metric(k),...
                char(Save_Group_Style(i)),'MarkerFaceColor',char(Save_Fill_Color(i))); hold on
        end
    end
    
    if strcmp(Scatter_Plot_Title,'Verification')
        k = 1000;
        %Measured_Metric = normrnd(1:1000,(Sigma_2_E/200)*(1:1000),[1 1000]);
        %Predicted_Metric = normrnd(0.5*(1:1000),0.02*(1:1000),[1 1000]);
        Measured_Metric = [1:k].*(1 + (Sigma_2_E/200)*randn(1,k));
        Predicted_Metric = 0.5*[1:k].*(1 + 0.02*randn(1,k));
        K = plot(Measured_Metric,Predicted_Metric,'ko'); hold on    
    end
    
    if k>0
        
        % statistics
        
        E_bar = mean(log(Measured_Metric));
        M_bar = mean(log(Predicted_Metric));
        M_hat = M_bar - E_bar + log(Measured_Metric);
        u  = sqrt( sum((log(Predicted_Metric)-M_hat).*(log(Predicted_Metric)-M_hat))/(k-1) );
        Sigma_E = Sigma_2_E/200;
        omega_M = sqrt( max(0,u*u - Sigma_E.^2) );
        delta = exp(M_bar-E_bar+0.5*omega_M.^2);
        Sigma_M = delta*omega_M;
        
        plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max],'k-')                    
        plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max*(1+2*Sigma_E)],'k--') 
        plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max*(1-2*Sigma_E)],'k--') 
        
        if strcmp(Model_Error,'yes')
            plot([Plot_Min,Plot_Max],[Plot_Min,delta*Plot_Max],'r-')
            plot([Plot_Min,Plot_Max],[Plot_Min,delta*Plot_Max*(1+2*Sigma_M)],'r--')
            plot([Plot_Min,Plot_Max],[Plot_Min,delta*Plot_Max*(1-2*Sigma_M)],'r--')
        end
        
        % format the legend and axis labels
        xlabel(Ind_Title,'Interpreter','LaTeX','FontSize',14)
        ylabel(Dep_Title,'Interpreter','LaTeX','FontSize',14)
        axis([Plot_Min Plot_Max Plot_Min Plot_Max])
        
        set(gca,'Units','inches')
        set(gca,'FontName','Times')
        set(gca,'FontSize',12)
        set(gca,'YTick',get(gca,'XTick'))
        set(gca,'Position',[1,1,4.5,4.5])
        
        text(Plot_Min+Title_Position(1)*(Plot_Max-Plot_Min),Plot_Min+Title_Position(2)*(Plot_Max-Plot_Min),...
            Scatter_Plot_Title,'FontSize',14,'FontName','Times','Interpreter','LaTeX')
  
        text(Plot_Min+(Title_Position(1)+0.05)*(Plot_Max-Plot_Min),Plot_Min+(Title_Position(2)-0.05)*(Plot_Max-Plot_Min),...
            ['$2 \, \sigma_E$=',num2str(2*Sigma_E,'%4.2f')],'FontSize',12,'FontName','Times','Interpreter','LaTeX')
        
        if strcmp(Model_Error,'yes')
            text(Plot_Min+(Title_Position(1)+0.05)*(Plot_Max-Plot_Min),Plot_Min+(Title_Position(2)-0.10)*(Plot_Max-Plot_Min),...
                ['$2 \, \sigma_M$=',num2str(2*Sigma_M,'%4.2f')],'FontSize',12,'FontName','Times','Interpreter','LaTeX')
            
            text(Plot_Min+(Title_Position(1)+0.05)*(Plot_Max-Plot_Min),Plot_Min+(Title_Position(2)-0.15)*(Plot_Max-Plot_Min),...
                ['Bias =',num2str(delta,'%4.2f')],'FontSize',12,'FontName','Times','Interpreter','LaTeX')
        end
        
        if strcmp(Scatter_Plot_Title,'Verification')==0
            C = stripcell(Group_Key_Label);
            [B I] = unique(C);
            legend(K(I),C(I),'Location',Key_Position,'FontSize',12)
            legend boxon
        end
        
        hold off
        
        % print to pdf
        set(gcf,'Visible','on');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperSize',[paper_width paper_height]);
        set(gcf,'PaperPosition',[0 0 paper_width paper_height]);
        display(['Printing scatter plot ',num2str(j),'...'])
        print(gcf,'-dpdf',[pwd,'/../../Manuals/',Plot_Filename])
        
    else
        display(['No data for scatter plot ',Scatter_Plot_Title])
    end
    
    clear Measured_Metric Predicted_Metric Group_Key_Label K
end




