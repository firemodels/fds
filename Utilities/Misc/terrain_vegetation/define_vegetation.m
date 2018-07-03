% Main script for vegetation definition in terrain:
%
% Instructions:
%
% 1. Go to smv repo and update it.
%
% 2. Navigate to the build directory of dem2fds: smv/Build/dem2fds and 
% select your compilation target.
% In burn to compile with Intel ifort select intel_linux_64 and execute
% ./make_dem2fds.sh
% This will provide the executable dem2fds_linux_64.
%
% 3. Execute the shell script that will invoke dem2fds. In our case we are
% interested in the Gatlinburg 1000mx1000m terrain case, defined in the
% dem2fds input script Gatlinburg_1000m.in. The script is
% Make_Gatlinburg_fds_input.sh
% This script runs dem2fds and produces both the OBST and GEOM version of
% the terrain of interest and also image data in .png format.
% For the Gatlinburg test case the name of the fds file using ZVALS (GEOM)
% is Gatlinburg_1000mg.fds
%
% Note: Several examples on the execution of dem2fds can be found in:
% smv/Build/dem2fds/data
%
% 4. Use the fds input file name (without the .fds extension) for geoms 
% on this script, and execute it.
% See below tree terrain and size distribution parameters. 
%
% Note: If you rotate the final matlab plot you can see the tree
% distribution in elevation.
%
% -------------------------------------------------------------------------
close all
clear all
clc

% Basedir and file name:
basedir  = './';
casename = 'Gatlinburg_1000m_g';

fds_file =[casename '.fds'];
png_file =[casename '.png'];

disp(['Loading Zvals and image data for case : ' casename '...'])

% Tree distribution over terrain:
% Number of trees to distribute:
N_TREES        =10000;

% Tree cone distance from terrain provided by Zvals:
H_FROM_TERRAIN = 3.00; % meters.

% Tree type : eoter block or cone:
TREE_TYPE = 'cone';

% Trees Properties: We produce a uniform distribution using rand().
% Mean Tree cone base radius:
RADIUS  = 2.00; % meters

% Maximum variation of cone base radius:
DRADIUS = 1.00; % meters

% Parameter to define a buffer layer from domain boundaries with no trees.
BUFF_DIST = 4*(RADIUS+DRADIUS);

% Mean cone height:
HEIGHT  = 10.0; % meters

% Maximum height variation
DHEIGHT =  2.5; % meters

% Mass per volume for the trees:
MASS_PER_VOLUME = 1.5;  % kg/m^3 

% Particles per fds fluid cell:
N_PARTICLES_PER_CELL = 1;

% Load Zvals from fds file:
t=cputime;
disp(' ')
fprintf('Loading information from .fds file ..');
[fid]=fopen([basedir fds_file],'r');
while(1)
   line=fgetl(fid);
   if(length(line) < 5)
       continue
   end
   % Check for the end:
   if(strcmp(line(1:5),'&TAIL'))
       break
   end
   % Check if we arrived at &GEOM:
   if(strcmp(line(1:5),'&GEOM'))
        line=fgetl(fid);
        if(strcmp(line(1:3),'IJK'))
            % Find XB:
            for ix=5:length(line)
                if(strcmp(line(ix),'X'))
                    break
                end
            end
            % Now do IJK, XB:
            size_zvals.IJK=str2num(line(5:ix-1));
            size_zvals.XB =str2num(line(ix+3:end));
        end
        line=fgetl(fid);
        if(strcmp(line(1:5),'ZVALS'))
           zvalsv=[];
           while(1)
              line=fgetl(fid);
              if(strcmp(line(end),'/'))
                  zvalsv = [zvalsv str2num(line(1:end-1))];
                  break
              end
              zvalsv = [zvalsv str2num(line)];
           end
           % Now build zvals 2D array:
           ij=0;
           zvals=zeros(size_zvals.IJK(1),size_zvals.IJK(2));
           for j=1:size_zvals.IJK(2)
               for i=1:size_zvals.IJK(1)
                  ij=ij+1;
                  zvals(i,j)=zvalsv(ij);
               end
           end
        end
   end   
end
fclose(fid);
dx = (size_zvals.XB(2)-size_zvals.XB(1))/size_zvals.IJK(1);
dy = (size_zvals.XB(4)-size_zvals.XB(3))/size_zvals.IJK(2);
x=linspace(size_zvals.XB(1),size_zvals.XB(2),size_zvals.IJK(1));
y=linspace(size_zvals.XB(3),size_zvals.XB(4),size_zvals.IJK(2));
fprintf(' done. Time taken : %6.4f sec.\n',cputime-t)

% Load png file, figure out dx, dy per pixel, rgb values:
t=cputime;
fprintf('Loading information from .png file ..');
IMG=imread([basedir png_file]);
nxp=length(IMG(1,:,1));
nyp=length(IMG(:,1,1));
dxp=(x(end)-x(1))/nxp;
dyp=abs((y(end)-y(1)))/nyp;
fprintf(' done. Time taken : %6.4f sec.\n',cputime-t)

% Now do Dart throwing for pine trees:
% Threshold for G value respect to R and B in pixel:
G_thresh       = 1.20;
TREE_XYZ = zeros(N_TREES,3);
t=cputime;
itry=0;
itree=0;
disp(' ')
disp('Producing tree distribution ..')
while (1)
    itry=itry+1;
    if mod(itry,500) == 0
        disp(['Number of trials : ' num2str(itry) ...
              ', itree :' num2str(itree)])
        rng('shuffle');
    end
    xi = (size_zvals.XB(1)+BUFF_DIST) + ...
         (size_zvals.XB(2)-size_zvals.XB(1)-2*BUFF_DIST)*rand;
    yj = (size_zvals.XB(3)+BUFF_DIST) + ...
         (size_zvals.XB(4)-size_zvals.XB(3)-2*BUFF_DIST)*rand;
    
    % Find RGB -> G value:
    jp=max(1,floor((xi-size_zvals.XB(1))/dxp));
    ip=nyp-max(1,floor((yj-size_zvals.XB(3))/dyp));
    
    R = IMG(ip,jp,1);
    G = IMG(ip,jp,2);
    B = IMG(ip,jp,3);
    flg = (G > G_thresh*R) & (G > G_thresh*B);
    if(~flg) 
        continue
    end  
    % Here this point gets a zval.
    ig=max(1,floor((xi-size_zvals.XB(1))/dx));
    jg=max(1,floor((yj-size_zvals.XB(3))/dy));
    
    ilo=max(1,ig-1);
    ihi=min(size_zvals.IJK(1),ig+1);
    jlo=max(1,jg-1);
    jhi=min(size_zvals.IJK(2),jg+1);     
    zk = max(max(zvals(ilo:ihi,jlo:jhi)));
    
    itree=itree+1;
    TREE_XYZ(itree,1:3) = [xi yj zk+H_FROM_TERRAIN];
    
    if(itree == N_TREES)
        disp(['Defined location for ' num2str(N_TREES) ...
              ' trees. Time taken :' num2str(cputime-t) ' sec.'])
        break
    end
end
disp('Producing tree distribution done.')
disp(' ')

% Write &INIT PART_ID file for all trees:
dist_file=[TREE_TYPE '_tree_init_NT_' num2str(N_TREES) ...
           '_NPC_' num2str(N_PARTICLES_PER_CELL) '.txt'];
t=cputime;
fprintf(['Writing tree distribution file :' dist_file])
[fid]=fopen([basedir dist_file],'w');
if strcmp(TREE_TYPE,'cone')
    for itree=1:N_TREES
        rnum=rand;
        RAD = (RADIUS-DRADIUS) + 2*DRADIUS*rnum;
        HGT = (HEIGHT-DHEIGHT) + 2*DHEIGHT*rnum;
        % Case Cones:
        fprintf(fid,'&INIT PART_ID=''foliage'', XYZ=%6.2f,%6.2f,%6.2f, RADIUS=%6.2f, HEIGHT=%6.2f, SHAPE=''CONE'', N_PARTICLES_PER_CELL=%d, MASS_PER_VOLUME=%6.2f /\n',...
            TREE_XYZ(itree,1),TREE_XYZ(itree,2),TREE_XYZ(itree,3),RAD,HGT,N_PARTICLES_PER_CELL,MASS_PER_VOLUME);
    end
elseif strcmp(TREE_TYPE,'block')
    for itree=1:N_TREES
        rnum=rand;
        RAD = (RADIUS-DRADIUS) + 2*DRADIUS*rnum;
        HGT = (HEIGHT-DHEIGHT) + 2*DHEIGHT*rnum;
        % Case Blocks:
        fprintf(fid,'&INIT PART_ID=''foliage'', XB=%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f, N_PARTICLES_PER_CELL=%d, MASS_PER_VOLUME=%6.2f /\n',...
            TREE_XYZ(itree,1)-RAD,TREE_XYZ(itree,1)+RAD,...
            TREE_XYZ(itree,2)-RAD,TREE_XYZ(itree,2)+RAD,...
            TREE_XYZ(itree,3)    ,TREE_XYZ(itree,3)+HGT,...
            N_PARTICLES_PER_CELL,MASS_PER_VOLUME);
    end
end    
fclose(fid);
fprintf(' done. Time taken : %6.4f sec.\n',cputime-t)

% Make Figure:
ttp=['^g'; 'sg'];
if(strcmp(TREE_TYPE,'cone'))
    tti=ttp(1,:);
else
    tti=ttp(2,:);
end
figure
hold on
surf(x,y,zvals','FaceAlpha',0.4)
image(size_zvals.XB([1 2]),size_zvals.XB([4 3]),IMG)
for itree=1:N_TREES
    plot3(TREE_XYZ(itree,1),TREE_XYZ(itree,2),TREE_XYZ(itree,3),tti)
end
axis equal