% McDermott
% 9-23-13
% layer_height.m
%
% replaces old layer_height.f (by Kevin McGrattan) used to post-process HGL cases

close all
clear all

% this is horrendous, but the ls command was giving file permission problems with Samba
list_dir = dir('../../Validation/*');
k = 0;
for i=1:length(list_dir)
    Output_Directory = (['../../../out/',list_dir(i).name,'/FDS_Output_Files']);
    if exist(Output_Directory)
        list_files = dir([Output_Directory,'/*HGL.input']);
        if size(list_files)>0
            for j=1:length(list_files)
                if ~strcmp(list_files(j).name(1),'.') % ignore hidden files
                    k=k+1;
                    output_dir{k} = ['../../../out/',list_dir(i).name,'/FDS_Output_Files/'];
                    input_file{k} = list_files(j).name;
                end
            end
        end
    end
end

%return % uncomment to just list the files for testing purposes

ntd=20000;
ncd=500;

for i=1:length(input_file) % input_file loop

    clear M
    tmp=zeros(ntd,ncd);

    fid = fopen([output_dir{i},input_file{i}],'r');

    % read number of trees
    ntrees = str2num(fgetl(fid));

    % read number of TCs in the tree
    ntc = str2num(fgetl(fid));
    for n=1:ntc
        tline = fgetl(fid);
        cell_array = textscan(tline,'%s','delimiter',' ','MultipleDelimsAsOne',1);
        S = cell_array{:};
        ztc(n)=str2num(S{1});
        for nn=1:ntrees
            icol(n,nn)=str2num(S{nn+1});
        end
    end

    % read weight of each tree
    tline = fgetl(fid);
    cell_array = textscan(tline,'%s','delimiter',' ','MultipleDelimsAsOne',1);
    S = cell_array{:};
    for nn=1:ntrees
        wgt(nn)=str2num(S{nn});
    end

    % read data file name
    infile = fgetl(fid);

    % read number of columns in data file
    nc = str2num(fgetl(fid));

    % read row number where data starts
    nr = str2num(fgetl(fid));

    % read ceiling height
    z_h = str2num(fgetl(fid));

    % read starting time
    t_start = str2num(fgetl(fid));

    % read name of output file
    outfile = fgetl(fid);
    fclose(fid);

    % read data from file
    M = importdata([output_dir{i},infile],',',nr-1);
    t = M.data(:,1);
    d = M.data(:,2:nc);

    z_0 = 0;
    for n=1:ntc-1
        z(n) = (ztc(n)+ztc(n+1))/2;
    end
    z(ntc) = z_h;

    fout = fopen([output_dir{i},outfile],'w');

    fprintf(fout,'%s, %s, %s, %s\n','Time','Height','T_lower','T_upper');

    for i=1:length(t) % time loop
        if t(i)<t_start
            continue
        end
        tmp(i,n) = 0;
        for nn=1:ntrees
            for n=1:ntc
                tmp(i,n) = tmp(i,n) + (273+d(i,icol(n,nn)-1))*wgt(nn);
            end
        end
        i1 = 0;
        i2 = 0;
        for n=1:ntc
            if n==1
                i1 = i1 + tmp(i,n)*(z(n)-z_0);
                i2 = i2 + (1./tmp(i,n))*(z(n)-z_0);
            else
                i1 = i1 + tmp(i,n)*(z(n)-z(n-1));
                i2 = i2 + (1./tmp(i,n))*(z(n)-z(n-1));
            end
        end
        zint(i)=tmp(i,1)*(i1*i2-z_h^2)/(i1+i2*tmp(i,1)^2-2*tmp(i,1)*z_h);
        tmpl(i)=tmp(i,1);
        i1 = 0;
        for n=1:ntc
            if n==1
                if z(n)>zint(i)
                    if z_0>=zint(i)
                        i1 = i1 + tmp(i,n)*(z(n)-z_0);
                    end
                    if z_0<zint(i)
                        i1 = i1 + tmp(i,n)*(z(n)-zint(i));
                    end
                end
            else
                if z(n)>zint(i)
                    if z(n-1)>=zint(i)
                        i1 = i1 + tmp(i,n)*(z(n)-z(n-1));
                    end
                    if z(n-1)<zint(i)
                        i1 = i1 + tmp(i,n)*(z(n)-zint(i));
                    end
                end
            end
        end
        tmph(i) = max(tmpl(i),(1./(z_h-zint(i)))*i1);

        fprintf(fout,'%f, %f, %f, %f\n',t(i),zint(i),tmpl(i)-273,tmph(i)-273);

    end % time_loop

    fclose(fout);

end % input_file loop




