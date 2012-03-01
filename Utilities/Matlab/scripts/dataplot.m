% R. McDermott and C. Cruz
% 7-06-2009
% dataplot.m
%
% [saved_data,drange] = dataplot(cfil,vdir,plotdir,[drange])
%
% Output:
%    saved_data - cell array containing data needed in scatplot.m
%
%    drange - data range needed for scatplot.m, must be commensurate with
%    saved_data
%
% Input:
%
%    cfil - base configuration file (set in master script)
%
%    vfil - base input file directory (set in master)
%
%    plotdir - base plot directory (set in master)
%
%    [optional] drange - a vector for the 'd' lines you want to read from the
%    config file.  For example, [2:5,7:8,10,12].
%    drange may also be a text string matching the 'Dataname' column in the
%    configuration file.  For example, [currently] the 'WTC' text string
%    identifies drange = [13:18], but is much easier to remember if you
%    happen to work with this set of validation cases frequently.
%
% Dependencies:
%    ../scripts/define_drow_variables.m
%    dvcread.m
%    parse.m
%
% Example: From the command line within the Matlab/functions/ directory,
%    type
%
%    >> [saved_data,drange] = dataplot(cfil,vdir,plotdir,[2:4,6:8]);
%
%    >> [saved_data,drange] = dataplot(cfil,vdir,plotdir,'WTC');

function [saved_data,drange] = dataplot(varargin)

if nargin<3||nargin>4; 
    display('Error in argument list')
end
if nargin>=3
    cfil = varargin{1};
    vdir = varargin{2};
    plotdir = varargin{3};
end
if nargin==4
    drange = varargin{4};
else
    drange = 2:2000;
end

% set the plot style parameters

plot_style
set(gcf,'DefaultLineLineWidth',Line_Width)
WPos = get(gcf,'Position');
set(gcf,'Position',[WPos(1) WPos(2) 640,420]);
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

% read the configuration file

A = importdata(cfil);
H = textscan(A{1},'%q','delimiter',',');
headers = H{:}'; clear H

if ~isnumeric(drange)
    dataname_col = find(strcmp(headers,'Dataname'));
    dstring = drange;
    drange_index = 0;
    clear drange
else
    dstring = 'null';
end

% process the "d" lines one by one

for i=2:2000
    
    if i>length(A); break; end
    
    P = textscan(A{i},'%q','delimiter',',');
    parameters = P{:}';
    
    % check for shortname specification instead of numeric drange
    
    if strcmp(dstring,'null')
        itest = ismember(i,drange);
    else
        itest = strcmp(parameters(dataname_col),dstring);
        if itest
            drange_index = drange_index + 1;
            drange(drange_index) = i;
        end
    end
    
    % check to see if d line has been activated in configuration file
    
    dtest = strcmp(parameters(find(strcmp(headers,'switch_id'))),'d');
    
    if itest & dtest
        
        define_drow_variables
        
        % save for scatter plots
        
        Save_Quantity(i)        = Quantity;
        Save_Group_Style(i)     = Group_Style;
        Save_Fill_Color(i)      = Fill_Color;
        Save_Group_Key_Label(i) = Group_Key_Label;
                
        % plot the experimental data or analytical solution (d1)
        
        if ~exist(d1_Filename,'file')
           disp(['Warning: File ' d1_Filename ' does not exist.'])
		   hold off
           continue
        end
        [H M] = dvcread(d1_Filename,d1_Col_Name_Row);
        R1 = parse(d1_Ind_Col_Name);
        S1 = parse(d1_Dep_Col_Name);
        style = parse(d1_Style);
        for j=1:length(S1)
            d1_Ind_Col = find(strcmp(H,R1(min(j,length(R1)))));
            d1_Dep_Col = find(strcmp(H,S1(j)));
            clear indices
            indices = find(d1_Comp_Start<=M(:,d1_Ind_Col) & M(:,d1_Ind_Col)<=d1_Comp_End);
            if strcmp(Metric,'max')
                Save_Measured_Metric(i,j) = max(M(indices,d1_Dep_Col))-d1_Initial_Value;
            elseif strcmp(Metric,'min')
                Save_Measured_Metric(i,j) = d1_Initial_Value-min(M(indices,d1_Dep_Col));
            elseif strcmp(Metric,'mean')
                Save_Measured_Metric(i,j) = mean(M(indices,d1_Dep_Col));
            else
                Save_Measured_Metric(i,j) = 0;
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
       
        if ~exist(d2_Filename,'file')
           disp(['Warning: File ' d2_Filename ' does not exist.'])
		   hold off
           continue
        end
        [H M] = dvcread(d2_Filename,d2_Col_Name_Row);
        R2 = parse(d2_Ind_Col_Name);
        S2 = parse(d2_Dep_Col_Name);
        style = parse(d2_Style);
        for j=1:length(S2)
            d2_Ind_Col = find(strcmp(H,R2(min(j,length(R2)))));
            d2_Dep_Col = find(strcmp(H,S2(j)));
            clear indices
            indices = find(d2_Comp_Start<=M(:,d2_Ind_Col) & M(:,d2_Ind_Col)<=d2_Comp_End);
            if strcmp(Metric,'max')
                Save_Predicted_Metric(i,j) = max(M(indices,d2_Dep_Col))-d2_Initial_Value;
            elseif strcmp(Metric,'min')
                Save_Predicted_Metric(i,j) = d2_Initial_Value-min(M(indices,d2_Dep_Col));
            elseif strcmp(Metric,'mean')
                Save_Predicted_Metric(i,j) = mean(M(indices,d2_Dep_Col));
            else
                Save_Predicted_Metric(i,j) = 0;
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
        
        if Plot_Type=='linear' & strcmp(Flip_Axis,'no')
            X_Title_Position = Min_Ind+Title_Position(1)*(Max_Ind-Min_Ind);
            Y_Title_Position = Min_Dep+Title_Position(2)*(Max_Dep-Min_Dep);
        elseif Plot_Type=='linear' & strcmp(Flip_Axis,'yes')
            X_Title_Position = Min_Dep+Title_Position(1)*(Max_Dep-Min_Dep);
            Y_Title_Position = Min_Ind+Title_Position(2)*(Max_Ind-Min_Ind);
        elseif Plot_Type=='loglog' & strcmp(Flip_Axis,'no')
            X_Title_Position = 10^(log10(Min_Ind)+Title_Position(1)*(log10(Max_Ind)-log10(Min_Ind)));
            Y_Title_Position = 10^(log10(Min_Dep)+Title_Position(2)*(log10(Max_Dep)-log10(Min_Dep)));
        elseif Plot_Type=='loglog' & strcmp(Flip_Axis,'yes')
            X_Title_Position = 10^(log10(Min_Dep)+Title_Position(1)*(log10(Max_Dep)-log10(Min_Dep)));
            Y_Title_Position = 10^(log10(Min_Ind)+Title_Position(2)*(log10(Max_Ind)-log10(Min_Ind)));
        end
        
        set(gca,'FontName',Font_Name)
        set(gca,'FontSize',Label_Font_Size)
     
        if strcmp(Flip_Axis,'no')
            xlabel(Ind_Title,'Interpreter','LaTeX','FontSize',Label_Font_Size)
            ylabel(Dep_Title,'Interpreter','LaTeX','FontSize',Label_Font_Size)
            axis([Min_Ind Max_Ind Min_Dep Max_Dep])
            text(X_Title_Position,Y_Title_Position,...
                Plot_Title,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter','LaTeX')
        else
            xlabel(Dep_Title,'Interpreter','LaTeX','FontSize',Label_Font_Size)
            ylabel(Ind_Title,'Interpreter','LaTeX','FontSize',Label_Font_Size)
            axis([Min_Dep Max_Dep Min_Ind Max_Ind])
            text(X_Title_Position,Y_Title_Position,...
                Plot_Title,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter','LaTeX')
        end
        if size(Key_Position)>0
            legend_handle = legend(K,[parse(d1_Key),parse(d2_Key)],'Location',Key_Position);
            set(legend_handle,'Interpreter','LaTeX');
            set(legend_handle,'Fontsize',Key_Font_Size);
            set(legend_handle,'Box','on');
            if size(d1_Tick)>0
               set(gca,'XTick',d1_Tick)
            end
            if size(d2_Tick)>0
               set(gca,'YTick',d2_Tick)
            end
            if size(Legend_XYWidthHeight)>0
               legend_position=get(legend_handle,'Position');
               if Legend_XYWidthHeight(1)>0; legend_position(1)=Legend_XYWidthHeight(1); end % X
               if Legend_XYWidthHeight(2)>0; legend_position(2)=Legend_XYWidthHeight(2); end % Y
               if Legend_XYWidthHeight(3)>0; legend_position(3)=Legend_XYWidthHeight(3); end % Width
               if Legend_XYWidthHeight(4)>0; legend_position(4)=Legend_XYWidthHeight(4); end % Height
               set(legend_handle,'Position',legend_position)
            end
        end
        
        % add SVN if file is available
        
        if exist(SVN_Filename,'file')
            SVN = importdata(SVN_Filename);
            x_lim = get(gca,'XLim');
            y_lim = get(gca,'YLim');
            if strcmp(Plot_Type,'loglog')
                X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
                Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
            else
                X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
                Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
            end
            text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
                'FontSize',10,'FontName',Font_Name,'Interpreter','LaTeX')
        end
        
        % print to pdf
		
		PDF_Paper_Width = Paper_Width_Factor*Paper_Width;
        
        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
        set(gcf,'PaperPosition',[0 0 PDF_Paper_Width Paper_Height]); 
        display(['Printing plot ',num2str(i),'...'])
        print(gcf,'-dpdf',[plotdir,Plot_Filename])
        
    end
    clear S1 S2 K style H M X Y P parameters
end

clear A

% pack data for use in scatter plot

saved_data = [{Save_Quantity'},...
              {Save_Group_Style'},...
              {Save_Fill_Color'},...
              {Save_Group_Key_Label'},...
              {Save_Measured_Metric},...
              {Save_Predicted_Metric}];

              
display('dataplot completed successfully!')


