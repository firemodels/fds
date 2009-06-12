% R. McDermott
% 6-11-2009
% scatterplot.m
%
% Generate scatter plots.  Must first run read_dline.m to generate
% measured_metric and predicted_metric.
%
% Dependencies:
%    define_qrow_variables.m
%    ../scatterplot_config_matlab.csv

close all

qfil = ['../scatterplot_config_matlab.csv'];
qrange = 2:2;

addpath('../functions')
paper_width  = 4.5; % inches
paper_height = 4.5; % inches

Q = importdata(qfil);
H = textscan(Q{1},'%q','delimiter',',');
headers = H{:}'; clear H

for j=qrange
    if j>length(Q); break; end
    
    define_qrow_variables
    
    for i=drange
        if strcmp(Save_Quantity(i),Scatter_Plot_Title)
            K(i) = plot(Save_Measured_Metric(i),Save_Predicted_Metric(i),...
                char(Save_Group_Style(i)),'MarkerFaceColor',char(Save_Fill_Color(i))); hold on
        end
    end
    hold off
    
    % format the legend and axis labels
    xlabel(Ind_Title,'Interpreter','LaTeX','FontSize',16)
    ylabel(Dep_Title,'Interpreter','LaTeX','FontSize',16)
    axis([Plot_Min Plot_Max Plot_Min Plot_Max])
    
    set(gca,'FontName','Times')
    set(gca,'FontSize',14)
    set(gca,'YTick',get(gca,'XTick'))
    
    text(Title_Position(1)*(Plot_Max-Plot_Min),Title_Position(2)*(Plot_Max-Plot_Min),...
        Scatter_Plot_Title,'FontSize',16,'FontName','Times','Interpreter','LaTeX')
    
    C = stripcell(Save_Group_Key_Label);
    [B I] = unique(C);
    legend(K(I),Save_Group_Key_Label(I),'Location',Key_Position,'FontSize',12)
    legend boxoff
end