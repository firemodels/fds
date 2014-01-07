
% This program reads heat transfer FEA model data e.g., nodal coordinates
% and element connectivity to determine which nodes are on the boundary of
% the structural components.
%
% It then reads six sets of FDS boundary file data (in +x, -x, +y, -y, +z, -z )
% directions containing temperatures
% at gridpoints located on the exposed surface of structural objects.
% It then maps FDS boundary tempeature data to the heat transfer FEA model
% surface nodes.

% The algorithm for finding nodes on the surface of the FEA model 
% works as follows. Using element connectivity data in the FEA model, 
% it obtains (1) a column vector containing the number of elements that use each node, 
%(2) a column vector pointing to the beginning of the element list for each node, and 
%(3) finally, a column vector containing a list of elements surrounding each node.
% Next, nodes on faces of each element are compared with nodes of other elements 
% that are connected to the same nodes to determine if all nodes along the 
% faces of each pair of elements match. If all nodes along common faces 
% of a pair of elements match, then these nodes are not on free surface. 
% When this match is not found, then nodes of the element face are on free surface 
% and these nodes are written to a column vector. A unique list of such nodes 
% can be obtained using a simple MATLAB functionality that removes any duplication 
% of surface nodes in the column vector containing surface nodes.

% The program creates two output files in csv format:
%
% 1. "testout.csv" - It contains FEA model surface nodes and their
% coordinates (node #, x, y, z)
%
% 2. "testout.csv" - It contains the FEA model surface nodes, their
% coordinates, and their mapped temperatures (from the FDS model) in 
% ANSYS table read format as shown below:
%
% First line: Header info (e.g. 0, 1, 2, 3, 4, 5)
% All other lines: Record #, node #, x, y, z, Temp)
% 
% At the end of this program, there are commands that can be used to
% generate MATLAB slice plots, MATLAB scatter plots. MATLAB slice plots can
% be really useful. However, these plots may have to be done in several
% steps. For example, for an  I-beam slice plots, 3 separete plots should
% be done: one for top flange, one for the web section, and the other one
% for the bottom flange.



function fds_to_ansys_data_mapping

clear all
close all

nuf=[];
nup=[];
nu=[];
kount=[];


free_node = [];
ff_node = [];
tansys = [];
tt = [];

% Read FEA heat transfer model nodal coordinates in csv format
%P = importdata('nodal_coord_ibeam_new.csv', ',');

fname = input('Enter FEA heat transfer model nodal coordinates file name without extension( [.csv extension assumed]): ', 's');

fname1 = strcat(fname, '.csv');
P = importdata(fname1, ',');


nodes = P(:,1);
xb = P(:,2);
yb = P(:,3);
zb = P(:,4);


nnod = length(nodes);

ns = 0;

nuf = zeros(nnod,1);
nuf1 = zeros(nnod,1);
nup = zeros(nnod,1);

n_stat = zeros(nnod,1);
kount = zeros(nnod,1);

% Read FEA heat transfer model element connecvity data in csv format
% ncon = importdata('elem_conn_ibeam_new.csv', ',');

fname2 = input('Enter FEA heat transfer model element connecvity file name without extension( [.csv extension assumed]): ', 's');
fname3 = strcat(fname2, '.csv');

ncon = importdata(fname3, ',');

nel = length(ncon);
el_stat = zeros(nel,1);


% The following section of the code develops a list of elements that surrounds each node  
% The list is in vector "nu".  The number of elements that uses a node "i" 
% is "nuf[i]".  "nup[i]" points to the beginning of the element list for  
% for each node i. 


face = zeros(6,5);
face(1,:)= [1 1 4 3 2];
face(2,:)= [2 1 2 6 5];
face(3,:)= [3 2 3 7 6];
face(4,:)= [4 3 4 8 7];
face(5,:)= [5 4 1 5 8];
face(6,:)= [6 5 6 7 8];

fpe = 6;
npf = 4;
npe = 8;

iindex = 1;

for i = 1:nel
    for j = 2:npe+1
       node = ncon(i,j);
       nuf(node) = nuf(node)+1;
    end
end

nup(1) = 1;

for i = ns+2:ns+nnod
    nup(i) = nup(i-1)+ nuf(i-1);
end

dim = nup(nnod) + nuf(nnod);
nu = zeros(dim,1);

for i = 1:nnod
    nuf(i) = 0;
end


for i = 1:nel
   for j = 2:npe+1 
      node = ncon(i,j);
      nu( ns+ nup( node ) + nuf( node ) ) = i;
%      kount(node) = nup( node ) + nuf( node );
      nuf(node) = nuf(node) + 1;
   end
end

for i = 1:nel
    for j = 1:fpe
        for k =2:npf+1
            kk = face(j,k);
            ll = ncon(i,kk+1);
            n_stat(ll) = 1;
        end
        
        kk = face(j,2);
        n = ncon(i,kk+1);
        ifind = 0;
        
        for l = 1:nuf(n)
            e = nu( nup(n) + l - 1);
            aa = nup(n) + l;
            
            if ( e == i )
                continue;
            end
                      
            mm = 0;
            
            for m = 2:npe+1
                s1 = ncon(e,m);
                if( n_stat(ncon(e,m)) == 1)
                    mm = mm + 1;
                end
            end
            
           if ( mm == npf )
                ifind = 1;
                break;
            end
            
        end
        
        if ( ifind == 0 )
         
            for m = 2:npe+1
                if( n_stat(ncon(i,m)) == 1)
                    free_node(iindex) = ncon(i,m);
                    iindex = iindex + 1;
                end
                
            end
        end
        
        for k =2:npf+1
            kk = face(j,k);
            ll = ncon(i,kk+1);
            n_stat(ll) = 0;
        end
        
    end
    
end    

ff = unique(free_node);

kk1 = length(ff);

% Write the FEA surface nodes and the nodal coordinates
fw = fopen(['test_out.csv'],'w' );

for i = 1:kk1
   jjj = ff(i);
   bound_node(i) = jjj;
   jindex = find( P(:,1) == jjj);
   xansys(i) = xb(jindex);
   yansys(i) = yb(jindex);
   zansys(i) = zb(jindex);
   fprintf(fw, '%s\n', [ num2str(jindex) ',' num2str(xb(jindex)) ',' num2str(yb(jindex)), ',' num2str(zb(jindex)) ]);
end

fclose(fw);


% READ FDS DATA
% Read FDS boundary data for all six directions
% FDS problem name "case_5_v1_test" in this example
% Input file fields should be separated by comma

% +x direction

fname4 = input('Enter FDS boundary file [+x direction] name without extension( [.csv extension assumed]): ', 's');
fname5 = strcat(fname4, '.csv');

%XPOS = importdata('case_5_v1_test_xp.prn', ' ');
XPOS = importdata(fname5, ',');

xp = XPOS(:,1);
yp = XPOS(:,2);
zp = XPOS(:,3);
tp = XPOS(:,4);

xp_length = length(XPOS);
index = 0;

for i=1:xp_length
    if( tp(i) > 20.0 )
        index = index + 1;
        xp1(index) = xp(i);
        yp1(index) = yp(i);
        zp1(index) = zp(i); 
        r_in(index) = tp(i);
     end
end

% -x direction
%
fname6 = input('Enter FDS boundary file [-x direction] name without extension( [.csv extension assumed]): ', 's');
fname7 = strcat(fname6, '.csv');

%XPOS = importdata('case_5_v1_test_xm.prn', ' ');
XPOS = importdata(fname7, ',');

xp = XPOS(:,1);
yp = XPOS(:,2);
zp = XPOS(:,3);
tp = XPOS(:,4);

xp_length = length(XPOS);

for i=1:xp_length
    if( tp(i) > 20.0 )
        index = index + 1;
        xp1(index) = xp(i);
        yp1(index) = yp(i);
        zp1(index) = zp(i); 
        r_in(index) = tp(i);
     end
end

% +y direction
%
fname8 = input('Enter FDS boundary file [+y direction] name without extension( [.csv extension assumed]): ', 's');
fname9 = strcat(fname8, '.csv');

%XPOS = importdata('case_5_v1_test_yp.prn', ' ');
XPOS = importdata(fname9, ',');

xp = XPOS(:,1);
yp = XPOS(:,2);
zp = XPOS(:,3);
tp = XPOS(:,4);

xp_length = length(XPOS);

for i=1:xp_length
    if( tp(i) > 20.0 )
        index = index + 1;
        xp1(index) = xp(i);
        yp1(index) = yp(i);
        zp1(index) = zp(i); 
        r_in(index) = tp(i);
     end
end

% -y direction
%
fname10 = input('Enter FDS boundary file [-y direction] name without extension( [.csv extension assumed]): ', 's');
fname11 = strcat(fname10, '.csv');

%XPOS = importdata('case_5_v1_test_ym.prn', ' ');
XPOS = importdata(fname11, ',');

xp = XPOS(:,1);
yp = XPOS(:,2);
zp = XPOS(:,3);
tp = XPOS(:,4);

xp_length = length(XPOS);

for i=1:xp_length
    if( tp(i) > 20.0 )
        index = index + 1;
        xp1(index) = xp(i);
        yp1(index) = yp(i);
        zp1(index) = zp(i); 
        r_in(index) = tp(i);
     end
end

% +z direction
%
fname12 = input('Enter FDS boundary file [+z direction] name without extension( [.csv extension assumed]): ', 's');
fname13 = strcat(fname12, '.csv');

XPOS = importdata(fname13, ',');
%XPOS = importdata('case_5_v1_test_zp.prn', ' ');

xp = XPOS(:,1);
yp = XPOS(:,2);
zp = XPOS(:,3);
tp = XPOS(:,4);

xp_length = length(XPOS);

for i=1:xp_length
    if( tp(i) > 20.0 )
        index = index + 1;
        xp1(index) = xp(i);
        yp1(index) = yp(i);
        zp1(index) = zp(i); 
        r_in(index) = tp(i);
     end
end

% -z direction
%
fname14 = input('Enter FDS boundary file [-z direction] name without extension( [.csv extension assumed]): ', 's');
fname15 = strcat(fname14, '.csv');

XPOS = importdata(fname15, ',');
%XPOS = importdata('case_5_v1_test_zm.prn', ' ');

xp = XPOS(:,1);
yp = XPOS(:,2);
zp = XPOS(:,3);
tp = XPOS(:,4);

xp_length = length(XPOS);

for i=1:xp_length
    if( tp(i) > 20.0 )
        index = index + 1;
        xp1(index) = xp(i);
        yp1(index) = yp(i);
        zp1(index) = zp(i); 
        r_in(index) = tp(i);
     end
end

xp2 = transpose(xp1);
yp2 = transpose(yp1);
zp2 = transpose(zp1);
r_in2 = transpose(r_in);

xansys1 = transpose(xansys);
yansys1 = transpose(yansys);
zansys1 = transpose(zansys);

xmin = min(xansys1);
xmax = max(xansys1);

ymin = min(yansys1);
ymax = max(yansys1);

zmin = min(zansys1);
zmax = max(zansys1);

% nodiv is number of divisions used to generate MATLAB griddata
% A default value of 30 is chosen.

nodiv = 30;

xx = linspace(xmin,xmax,nodiv);
yy = linspace(ymin,ymax,nodiv);
zz = linspace(zmin,zmax,nodiv);


[xinterp,yinterp,zinterp] = meshgrid(xx,yy,zz);

tansys=griddata(xp2,yp2,zp2,r_in2,xinterp,yinterp,zinterp,'nearest');


for i = 1:kk1
    xcoord = xansys1(i);
    ycoord = yansys1(i);
    zcoord = zansys1(i);
    
    [d1, p1] = min( abs(xx-xcoord) );
    [d2, p2] = min( abs(yy-ycoord) );
    [d3, p3] = min( abs(zz-zcoord) );
    
    tt(i) = tansys(p2,p1,p3);
    
end

% tt1 contains the mapped data to the FEA boundary nodes

tt1 = transpose(tt);
fw1 = fopen('test_out1.csv','w' );

% data written such that they can be read by ANSYS in TABLE format using
% TREAD command
%
% header info for ANSYS 
% Additional first column for TREAD command to work
%
fprintf(fw1, '%s\n', [ num2str(0) ',' num2str(1) ',' num2str(2) ',' num2str(3), ',' num2str(4), ',' num2str(5)]);
%
% Write the data
%
for i = 1:kk1
   fprintf(fw1, '%s\n', [ num2str(i) ',' num2str(bound_node(i)) ',' num2str(xansys1(i)) ',' num2str(yansys1(i)), ',' num2str(zansys1(i)), ',' num2str(tt1(i))]);
end
fclose(fw1);


% Generate plot in MATLAB (scatter plot)
% figure
% plot3(xansys1, yansys1, tt1 ,'.')

%view(3)
%axis on
%grid on
%light
%lighting phong
%camlight('left')
%shading interp


% Generate slice plots in MATLAB
% Sample example is shown below for a beam. 3 slice plots are created: 1
% for top flange, 1 for web, and the last one for the bottom flange.

% Example 1
%figure
% Top flange 
%xslice = [3.0, 3.4];
%yslice = [0.5, 4.5];
%zslice = [2.925,3.1];
%slice(xinterp,yinterp,zinterp,tansys,xslice,yslice,zslice);

%hold on;

% web
%xslice = [3.125, 3.275];
%yslice = [0.5, 4.5];
%zslice = [];
%slice(xinterp,yinterp,zinterp,tansys,xslice,yslice,zslice);

%hold on;

% bottom flange
%xslice = [3.0, 3.4];
%yslice = [0.5, 4.5];
%zslice = [2.0,2.175];
%slice(xinterp,yinterp,zinterp,tansys,xslice,yslice,zslice)
%%%%
%
% Example 2
% Slice example

%xmin1 = 3;
%xmax1 = 3.4;
%ymin1 = 0.5;
%ymax1 = 4.5;
%zmin1 = 2.925;
%zmax1 = 3.1;

%xmin2 = 3.125;
%xmax2 = 3.275;
%ymin2 = 0.5;
%ymax2 = 4.5;
%zmin2 = 2.175;
%zmax2 = 2.925;

%xmin3 = 3;
%xmax3 = 3.4;
%ymin3 = 0.5;
%ymax3 = 4.5;
%zmin3 = 2.0;
%zmax3 = 2.175;

%nodiv = 30;

% top flange
%xx1 = linspace(xmin1,xmax1,nodiv);
%yy1 = linspace(ymin1,ymax1,nodiv);
%zz1 = linspace(zmin1,zmax1,nodiv);

%[xinterp1,yinterp1,zinterp1] = meshgrid(xx1,yy1,zz1);
%tansys1=griddata(xp2,yp2,zp2,r_in2,xinterp1,yinterp1,zinterp1,'nearest');
%figure

%xslice = [3.0,3.4];
%yslice = [0.5,4.5];
%zslice = [2.925,3.1];


%slice(xinterp1,yinterp1,zinterp1,tansys1,xslice,yslice,zslice);

%hold on;

% Web

%xx2 = linspace(xmin2,xmax2,nodiv);
%yy2 = linspace(ymin2,ymax2,nodiv);
%zz2 = linspace(zmin2,zmax2,nodiv);

%[xinterp2,yinterp2,zinterp2] = meshgrid(xx2,yy2,zz2);
%tansys2=griddata(xp2,yp2,zp2,r_in2,xinterp2,yinterp2,zinterp2,'nearest');


%xslice = [3.125,3.275];
%yslice = [0.5,4.5];
%zslice = [];

%slice(xinterp2,yinterp2,zinterp2,tansys2,xslice,yslice,zslice);

%hold on;

% Bottom flange
%xx3 = linspace(xmin3,xmax3,nodiv);
%yy3 = linspace(ymin3,ymax3,nodiv);
%zz3 = linspace(zmin3,zmax3,nodiv);

%[xinterp3,yinterp3,zinterp3] = meshgrid(xx3,yy3,zz3);

%tansys3=griddata(xp2,yp2,zp2,r_in2,xinterp3,yinterp3,zinterp3,'nearest');


%xslice = [3.0,3.4];
%yslice = [0.5,4.5];
%zslice = [2.0,2.175];

%slice(xinterp3,yinterp3,zinterp3,tansys3,xslice,yslice,zslice);
%shading interp;
%colorbar;









