%define surfaces to calculate distances
ab.BL = 0;
ab.SR = 0.4572;
ab.SL = -0.4572;
ab.TF = 2.2860;
ab.F = 2.0574;
bdStraight.AB = 0.24;
bdStraight.BL = -0.24;
bdStraight.BK = 0.32;
bdStraight.FR = -0.32;
bdTee.AB = 0.72;
bdTee.BL = -0.24;
bdTee.FL = -0.32;
bdTee.FU = -0.32;
bdElbow.BL = -0.24;
bdElbow.FR = -0.32;
bdElbow.BK = 0.32;
bdElbow.RT = 4.24;
plotShapes = {'-s','-o','-*','-^','-x','-p','-.'};

%[filename,path] = uigetfile('*.csv');
path = uigetdir;
cd(path);
files = dir(fullfile(path, '*_line.csv'));

for f=1:length(files)
    figure('visible','off');
    filename = files(f).name;
    
    %use filename to choose geometry
    config = regexp(filename,'c(.*)_', 'tokens');
    if size(config)==0
        config{1,1}{1,1}='';
    end
    geometry = eval(strcat(filename(1:2),config{1,1}{1,1}));
    
    table = readtable(strcat(path,'/',filename));
    data = table2array(table);
    [rows, columns] = size(data);
    
    %get list of DEVC locations
    locations = fieldnames(geometry);
    
    %iterate through locations
    for i = 1:length(locations)
        X=zeros(1,1);
        Y=zeros(1,1);
        %get a logical array of columns relevant to location
        idx = strncmp(data(1,:),locations{i},length(locations{1}));
        %ignore the columns with distances
        idx(1:2:end)=0;
        %get the indices of the logical array
        indices = find(idx(1,:)==1);
        
        
        %iterate through the planes of DEVCs at each distance and take max
        %incident energy
        for j = 2:rows
            if not(strcmp(data{j,(indices(1)-1)},'') | strcmp(data{j,(indices(1)-1)},'NaN'))
                X(j-1) = abs(cellfun(@str2num,data(j,(indices(1)-1)))-geometry.(locations{i}));
                Y(j-1) = max(cellfun(@str2num,data(j,idx)))/1000;
            end
        end
        plot(X,Y,plotShapes{i})
        hold on
        syms x;
        curve = fit(X',Y','Exp2');
        %plot(curve);
        TScurve = @(x)30-curve(x);
        TPcurve = @(x)15-curve(x);
        tsZOI(i) = fzero(TScurve,0);
        tpZOI(i) = fzero(TPcurve,0);
        
        
    end
    plotTitle = regexprep(filename,'_line\.csv','');
    plotTitle = regexprep(plotTitle,'_',' ');
    title(plotTitle)
    xlabel('Distance from Surface (m)')
    ylabel('Incident Energy (MJ/m^{2})')
    plot(xlim, 15*[1 1])
    plot(xlim, 30*[1 1])
    legend(locations)
    set(gca,'FontSize',15)
    yl = ylim;
    for i=1:length(locations)
        text(1.5,(0.8*yl(2))-(0.08*yl(2)*i),strcat(locations{i},': TS: ',num2str(tsZOI(i)),'m',' TP: ',num2str(tpZOI(i)),'m'));
    end
    print(regexprep(filename,'.csv','.pdf'),'-dpdf')
    
end