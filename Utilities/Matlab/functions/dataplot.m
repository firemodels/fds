% R. McDermott and C. Cruz
% 7-06-2009
% dataplot.m
%
% [saved_data,drange] = dataplot(data_type,[drange])
%
% Output:
%    saved_data - cell array containing data needed in scatplot.m
%
%    drange - data range needed for scatplot.m, must be commensurate with
%    saved_data
%
% Input:
%    data_type - a character string with one of two options:
%    'verification' or 'validation'; this argument is just used to simplify
%    the path for the output files and configuration files.
%
%    [optional] drange - a vector for the 'd' lines you want to read from the
%    config file.  For example, [2:5,7:8,10,12].
%
% Dependencies:
%    ../scripts/define_drow_variables.m
%    dvcread.m
%    parse.m
%
% Example: From the command line within the Matlab/functions/ directory,
%    type
%    >> [saved_data,drange] = dataplot('validation',[2:4,6:8]);

function [saved_data,drange] = dataplot(varargin)

if nargin<1|nargin>2; 
    display('Error in argument list')
end
if nargin>=1
    if strcmp(varargin{1},'verification')|strcmp(varargin{1},'Verification')
        cfil = ['../verification_data_config_matlab.csv'];
        vdir = ['../../../Verification/'];
    elseif strcmp(varargin{1},'validation')|strcmp(varargin{1},'Validation')
        cfil = ['../validation_data_config_matlab.csv'];
        vdir = ['../../../Validation/'];
    end
end
if nargin==2
    drange = varargin{2};
else
    drange = 2:1055;
end

addpath('../scripts')
paper_width  = 6.0; % inches
paper_height = 4.5; % inches

A = importdata(cfil);
H = textscan(A{1},'%q','delimiter',',');
headers = H{:}'; clear H

for i=drange
    if i>length(A); break; end
    
    P = textscan(A{i},'%q','delimiter',',');
    parameters = P{:}';
    
    if strcmp(parameters(find(strcmp(headers,'switch_id'))),'d')
        
        define_drow_variables
        
        % save for scatter plots
        
        Save_Quantity(i)        = Quantity;
        Save_Group_Style(i)     = Group_Style;
        Save_Fill_Color(i)      = Fill_Color;
        Save_Group_Key_Label(i) = Group_Key_Label;
        
        % plot the experimental data or analytical solution (d1)
        
        [H M] = dvcread(d1_Filename,d1_Col_Name_Row);
        d1_Ind_Col = find(strcmp(H,d1_Ind_Col_Name));
        S1 = parse(d1_Dep_Col_Name);
        style = parse(d1_Style);
        for j=1:length(S1)
            d1_Dep_Col = find(strcmp(H,S1(j)));
            clear indices
            indices = find(d1_Comp_Start<=M(:,d1_Ind_Col) & M(:,d1_Ind_Col)<=d1_Comp_End);
            if Metric=='max'
                Save_Measured_Metric(i) = max(M(indices,d1_Dep_Col))-d1_Initial_Value;
            end
            if Metric=='min'
                Save_Measured_Metric(i) = d1_Initial_Value-min(M(indices,d1_Dep_Col));
            end
            clear indices
            indices = find(d1_Start<=M(:,d1_Ind_Col) & M(:,d1_Ind_Col)<=d1_End);
            if strcmp(Flip_Axis,'no')
                X = M(indices,d1_Ind_Col)/Scale_Ind;
                Y = M(indices,d1_Dep_Col)/Scale_Dep;
            else
                X = M(indices,d1_Dep_Col)/Scale_Dep;
                Y = M(indices,d1_Ind_Col)/Scale_Ind;
            end
            if Plot_Type=='linear'
                K(j) = plot(X,Y,char(style(j))); hold on
            elseif Plot_Type=='loglog'
                K(j) = loglog(X,Y,char(style(j))); hold on
            end
        end

        % plot the FDS or model data (d2)
        
        [H M] = dvcread(d2_Filename,d2_Col_Name_Row);
        d2_Ind_Col = find(strcmp(H,d2_Ind_Col_Name));
        S2 = parse(d2_Dep_Col_Name);
        style = parse(d2_Style);
        for j=1:length(S2)
            d2_Dep_Col = find(strcmp(H,S2(j)));
            clear indices
            indices = find(d2_Comp_Start<=M(:,d2_Ind_Col) & M(:,d2_Ind_Col)<=d2_Comp_End);
            if Metric=='max'
                Save_Predicted_Metric(i) = max(M(indices,d2_Dep_Col))-d2_Initial_Value;
            end
            if Metric=='min'
                Save_Predicted_Metric(i) = d2_Initial_Value-min(M(indices,d2_Dep_Col));
            end
            clear indices
            indices = find(d2_Start<=M(:,d2_Ind_Col) & M(:,d2_Ind_Col)<=d2_End);
            if strcmp(Flip_Axis,'no')
                X = M(indices,d2_Ind_Col)/Scale_Ind;
                Y = M(indices,d2_Dep_Col)/Scale_Dep;
            else
                X = M(indices,d2_Dep_Col)/Scale_Dep;
                Y = M(indices,d2_Ind_Col)/Scale_Ind;
            end
            if Plot_Type=='linear'
                K(length(S1)+j) = plot(X,Y,char(style(j)));
            elseif Plot_Type=='loglog'
                K(length(S1)+j) = loglog(X,Y,char(style(j)));
            end
        end
        hold off
        
        % format the legend and axis labels
        
        set(gca,'Units','inches')
        set(gca,'FontName','Times')
        set(gca,'FontSize',14)
        %set(gca,'PlotBoxAspectRatio',[1,0.7,1])
        set(gca,'Position',[1,0.75,4.5,3.15])
        
        if strcmp(Flip_Axis,'no')
            xlabel(Ind_Title,'Interpreter','LaTeX','FontSize',14)
            ylabel(Dep_Title,'Interpreter','LaTeX','FontSize',14)
            axis([Min_Ind Max_Ind Min_Dep Max_Dep])
            text(Min_Ind+Title_Position(1)*(Max_Ind-Min_Ind),Min_Dep+Title_Position(2)*(Max_Dep-Min_Dep),...
                Plot_Title,'FontSize',14,'FontName','Times','Interpreter','LaTeX')
        else
            xlabel(Dep_Title,'Interpreter','LaTeX','FontSize',14)
            ylabel(Ind_Title,'Interpreter','LaTeX','FontSize',14)
            axis([Min_Dep Max_Dep Min_Ind Max_Ind])
            text(Min_Dep+Title_Position(1)*(Max_Dep-Min_Dep),Min_Ind+Title_Position(2)*(Max_Ind-Min_Ind),...
                Plot_Title,'FontSize',14,'FontName','Times','Interpreter','LaTeX')
        end
        if size(Key_Position)>0
            legend(K,[parse(d1_Key),parse(d2_Key)],'Location',Key_Position,'Interpreter','LaTeX','FontSize',12)
            legend boxon
        end
        
        % print to pdf
        set(gcf,'Visible','on');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperSize',[paper_width paper_height]);
        set(gcf,'PaperPosition',[0 0 paper_width paper_height]); 
        display(['Printing plot ',num2str(i),'...'])
        print(gcf,'-dpdf',['../../../Manuals/',Plot_Filename])
        
    end
    clear S1 S2 K style H M X Y P parameters
end
clear A

% pack data for use in scatter plot
saved_data = [{Save_Quantity'},...
              {Save_Group_Style'},...
              {Save_Fill_Color'},...
              {Save_Group_Key_Label'},...
              {Save_Measured_Metric'},...
              {Save_Predicted_Metric'}];
              
display('Done!')
display('Why? Because...')
why


