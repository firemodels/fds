
% This routine maps temperatures from a solid thermal FEA model to a 
% structural beam element model. This is is the second of the two routines for mapping
% results for beam elements. In this approach, transfer points coordinates are 
% chosen along the beam cross section and  the transfer point coordinates are entered in a file by the user.
% Output files for each temperature data sets are created. Output
% temperatures are written out for each transfer point belonging to a node
% Node #, Temp(TP1), Temp(TP2),......
%

% Input files (examples):
%
% 1. Structural Model (beam elements):
% a. beam_nodal_coordinates.csv
% b. transfer_coordinates.csv
% c. beam_element_connectivity.csv
%
% 2. Heat Transfer Model (solid elements):
% a. solid_element_connectivity.csv
% b. solid_nodal_coordinates.csv
%

function beam_script_neutral

clear all
close all

b_nodes = [];
x_coord = [];
y_coord = [];
z_coord = [];
tp1 = [];

x_cen = [];
y_cen = [];
z_cen = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!!!!!!!!!!!! Read Structural Model Data !!!!!!!!!!!!!!!!!!!!!!!!!!!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% read axial coordinates of nodes in the beam model in ascending order
% format: node#, x_coordinates
%
% y and z coordinates for all structural (beam model) nodes should be the same. Only
% x-coordinates for structural nodes will vary.

fname = input('Enter structural model nodal coordinates file name without extension( [.csv extension assumed]): ', 's');
fname1 = strcat(fname, '.csv');
P = importdata(fname1, ',');

%P = importdata('beam_nodal_coordinates.csv', ',');

b_nodes = P(:,1);
x_coord = P(:,2);

nnod = length(b_nodes);
tp_mean = [];


% beam elements are assumed to be linking two cosecutive beam nodes
% format: element#, node#1, node#2
% element numbering should follow nodal list sequences, e.g.:
% e1, n1, n2
% e2, n2, n3
% e3, n3, n4
% ......

fname2 = input('Enter structural model element connecvity file name without extension( [.csv extension assumed]): ', 's');
fname3 = strcat(fname2, '.csv');
R = importdata(fname3, ',');

%R = importdata('beam_element_connectivity.csv', ',');
elems = R(:,1);
bn1 = R(:,2);
bn2 = R(:,3);

nels = length(elems);

% Enter Transfer Point information for each structural node
%
% It is assumed that all structural nodes have same number of transfer
% points with same y and z-coordinates. The only difference is that the
% axial dimensions vary from node to node. Axial dimension of each transfer
% point is the x-coordinate value of the corresponding node, which is read
% from a separate file as shown earlier.
%
% Format: Transfer Point #, y_coordinate, z_coordinate


fname4 = input('Enter structural node transfer point coordinate file name without extension( [.csv extension assumed]): ', 's');
fname5 = strcat(fname4, '.csv');

S = importdata(fname5, ',');
tpnum = S(:,1);
tpy = S(:,2);
tpz = S(:,3);

% example file name: transfer_coordinates.csv


% np1 is the number of transfer points for each structural node. 
% The number of transfer points is same for all nodes.

np1 = length( S(:,1) );


% x_coord, y_coord, and z_coord arrays have the coordinates of the transfer points 


for i = 1:np1
   y_coord(i) = tpy(i);
   z_coord(i) = tpz(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% !!!!!!!!!!!! Read Heat Transfer Thermal Model Data !!!!!!!!!!!!!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enter heat transfer FEA model info
% Enter element connectivities: elem#, node1, node2,....,node8 (brick element)

fname6 = input('Enter solid HT FEA element connectivity file name without extension( [.csv extension assumed]): ', 's');
fname7 = strcat(fname6, '.csv');
M = importdata(fname7, ',');

%M = importdata('solid_element_connectivity.csv', ',');

el_conn_length = length(M);
nel = M(el_conn_length,1);


%Enter FEA nodal coordiates for the HT model: node #, x, y, z

fname8 = input('Enter solid HT FEA nodal coordinates file name without extension( [.csv extension assumed]): ', 's');
fname9 = strcat(fname8, '.csv');
N = importdata(fname9, ',');

%N = importdata('solid_nodal_coordinates.csv', ',');

node = N(:, 1);
xx = N(:, 2);
yy = N(:, 3);
zz = N(:, 4);

%
%Enter nodal temperatures from HT model for all time steps:
%node #, temp for step 1, temp for step 2,........,temp for step nt
%Temp1 = importdata('solid_thermal_temperatures.csv', ',');

fname10 = input('Enter solid HT FEA nodal temperatures file name without extension( [.csv extension assumed]): ', 's');
fname11 = strcat(fname10, '.csv');
Temp1 = importdata(fname11, ',');

nt = size(Temp1,2) - 1;

%
% t_int is the computed temperatures for transfer points for each of the
% "nnod" nodes in the structural beam model for each of the "nt"
% temperature data sets

t_int = zeros(np1,nnod,nt);


%determine centroids of each FEA heat transfer model elements

[x_cen, y_cen, z_cen] = centroid(M, N);


for jjj = 1:np1
   tp1(jjj) = 0.;
end

for jnode = 2:nt+1
    Temp = [Temp1(:,1),Temp1(:,jnode)];
    
    for inode = 1:nnod
        [tp1]=beamnodes(M,N,node,Temp,x_coord,y_coord,z_coord,inode,x_cen, y_cen, z_cen);
        t_int(:,inode,jnode) = tp1;
    end
    
end



for jnode = 2:nt+1
    % output written to a separate file for each set of temperatures entered
    % Format: Node#, Temp_pt_1, Temp_pt_2,........,Temp_pt_n
    
    out_filename = num2str(jnode-1);
    fw = fopen([out_filename 'out.csv'],'w' );
    
        
    for inode = 1:nnod
        
        tint1 = t_int(:,inode,jnode);
    
        fprintf(fw, '%s', num2str(inode), ',' );
        
        for i = 1:np1
            fprintf(fw, '%s', num2str(tint1(i)), ',' );
        end
        
        fprintf(fw, '\n');
        
    end
    
    fclose(fw);
end



function [tp1] = beamnodes(M,N,node,Temp,xp1,yp1,zp1,inode,x_cen,y_cen,z_cen)

% Finds the hexahedral elements which contain the transfer points.
% Shape functions of hexahedral elements are used in natural coordinates.
% Uses finite element concepts to locate the element that contains the transfer point.
% Average temperature at each transfer point is computed with shape
% functions of the HT FEA model element and nodal temperatures for that
% element.

wt = [];

g = [];
dgdr = [];
dgds = [];
dgdt = [];

found = 0;
j = 1;
jj = 1;
mean_temp = 0.;

el_conn_length = length(M);
nel = M(el_conn_length,1);

connectedData = zeros(8,3);

index = 0;

np = length(yp1);

for jjj = 1:np
   tp1(jjj) = 0.;
end

for kk = 1:np
    
    found = 0;
    x_pt = xp1(inode);
    y_pt = yp1(kk);
    z_pt = zp1(kk);
    
    
% select a few elements in the vicinity

   xp = repmat(x_pt,1,nel);
   yp = repmat(y_pt,1,nel);
   zp = repmat(z_pt,1,nel);
    
   r_FE = sqrt( (x_cen - xp).^2 + (y_cen - yp).^2 + (z_cen - zp).^2);
   
   % search for those elements that are within the "r_FE" distance
   % from the transfer point.
   % A default value of 30 is used. This means that those elements whose
   % centroids are located within a distance of 30 from the transfer point
   % will be searched.
   % 
   
   ind_FE_keep = find(r_FE<= 25 );
   
   if isempty(ind_FE_keep)
      continue;
   end
   
   
   ne = length(ind_FE_keep);
   
    
    while ( found == 0 ) && ( j <= ne )
        
        r = 0.;
        s = 0.;
        t = 0.;
        sum = 1.;
        g_iter = 0;
        
        connectedData = zeros(8,3);
        
        index1 = ind_FE_keep(j);
        
        if isempty(index1)
            j = j + 1;
            continue;
        end
        
        index = find(M(:,1) == index1);
        
        if isempty(index)
            j = j + 1;
            continue;
        end
        
        for k = 1:8
            x(k)=0.;
            y(k)=0.;
            z(k)=0.;
        end
        
        for k = 1:8
            connectedNode = M(index,k+1);
            mask = node == connectedNode;
            connectedData(k,1:3) = N(mask,2:4);
            x(k) = connectedData(k,1);
            y(k) = connectedData(k,2);
            z(k) = connectedData(k,3);
            el_conn(k) = connectedNode;
        end
        
        while( sum > 0.001) && ( g_iter < 10 )
            
            g(1) = ( 1. - r ) * ( 1. - s ) * ( 1. - t ) / 8.;
            g(2) = ( 1. + r ) * ( 1. - s ) * ( 1. - t ) / 8.;
            g(3) = ( 1. + r ) * ( 1. + s ) * ( 1. - t ) / 8.;
            g(4) = ( 1. - r ) * ( 1. + s ) * ( 1. - t ) / 8.;
            g(5) = ( 1. - r ) * ( 1. - s ) * ( 1. + t ) / 8.;
            g(6) = ( 1. + r ) * ( 1. - s ) * ( 1. + t ) / 8.;
            g(7) = ( 1. + r ) * ( 1. + s ) * ( 1. + t ) / 8.;
            g(8) = ( 1. - r ) * ( 1. + s ) * ( 1. + t ) / 8.;
            
            dgdr(1) = - ( 1. - s ) * ( 1. - t ) / 8.;
            dgdr(2) =   ( 1. - s ) * ( 1. - t ) / 8.;
            dgdr(3) =   ( 1. + s ) * ( 1. - t ) / 8.;
            dgdr(4) = - ( 1. + s ) * ( 1. - t ) / 8.;
            dgdr(5) = - ( 1. - s ) * ( 1. + t ) / 8.;
            dgdr(6) =   ( 1. - s ) * ( 1. + t ) / 8.;
            dgdr(7) =   ( 1. + s ) * ( 1. + t ) / 8.;
            dgdr(8) = - ( 1. + s ) * ( 1. + t ) / 8.;
            
            dgds(1) = - ( 1. - r ) * ( 1. - t ) / 8.;
            dgds(2) = - ( 1. + r ) * ( 1. - t ) / 8.;
            dgds(3) =   ( 1. + r ) * ( 1. - t ) / 8.;
            dgds(4) =   ( 1. - r ) * ( 1. - t ) / 8.;
            dgds(5) = - ( 1. - r ) * ( 1. + t ) / 8.;
            dgds(6) = - ( 1. + r ) * ( 1. + t ) / 8.;
            dgds(7) =   ( 1. + r ) * ( 1. + t ) / 8.;
            dgds(8) =   ( 1. - r ) * ( 1. + t ) / 8.;
            
            dgdt(1) = - ( 1. - r ) * ( 1. - s ) / 8.;
            dgdt(2) = - ( 1. + r ) * ( 1. - s ) / 8.;
            dgdt(3) = - ( 1. + r ) * ( 1. + s ) / 8.;
            dgdt(4) = - ( 1. - r ) * ( 1. + s ) / 8.;
            dgdt(5) =   ( 1. - r ) * ( 1. - s ) / 8.;
            dgdt(6) =   ( 1. + r ) * ( 1. - s ) / 8.;
            dgdt(7) =   ( 1. + r ) * ( 1. + s ) / 8.;
            dgdt(8) =   ( 1. - r ) * ( 1. + s ) / 8.;
            
            rsx = x_pt - ( g(1) * x(1) + g(2) * x(2) + g(3) * x(3) + g(4) * x(4) + g(5) * x(5) + g(6) * x(6) + g(7) * x(7) + g(8) * x(8) );
            rsy = y_pt - ( g(1) * y(1) + g(2) * y(2) + g(3) * y(3) + g(4) * y(4) + g(5) * y(5) + g(6) * y(6) + g(7) * y(7) + g(8) * y(8) );
            rsz = z_pt - ( g(1) * z(1) + g(2) * z(2) + g(3) * z(3) + g(4) * z(4) + g(5) * z(5) + g(6) * z(6) + g(7) * z(7) + g(8) * z(8) );
            
            d11 = dgdr(1) * x(1) + dgdr(2) * x(2) + dgdr(3) * x(3) + dgdr(4) * x(4) + dgdr(5) * x(5) + dgdr(6) * x(6) + dgdr(7) * x(7) + dgdr(8) * x(8);
            d12 = dgds(1) * x(1) + dgds(2) * x(2) + dgds(3) * x(3) + dgds(4) * x(4) + dgds(5) * x(5) + dgds(6) * x(6) + dgds(7) * x(7) + dgds(8) * x(8);
            d13 = dgdt(1) * x(1) + dgdt(2) * x(2) + dgdt(3) * x(3) + dgdt(4) * x(4) + dgdt(5) * x(5) + dgdt(6) * x(6) + dgdt(7) * x(7) + dgdt(8) * x(8);
            d21 = dgdr(1) * y(1) + dgdr(2) * y(2) + dgdr(3) * y(3) + dgdr(4) * y(4) + dgdr(5) * y(5) + dgdr(6) * y(6) + dgdr(7) * y(7) + dgdr(8) * y(8);
            d22 = dgds(1) * y(1) + dgds(2) * y(2) + dgds(3) * y(3) + dgds(4) * y(4) + dgds(5) * y(5) + dgds(6) * y(6) + dgds(7) * y(7) + dgds(8) * y(8);
            d23 = dgdt(1) * y(1) + dgdt(2) * y(2) + dgdt(3) * y(3) + dgdt(4) * y(4) + dgdt(5) * y(5) + dgdt(6) * y(6) + dgdt(7) * y(7) + dgdt(8) * y(8);
            d31 = dgdr(1) * z(1) + dgdr(2) * z(2) + dgdr(3) * z(3) + dgdr(4) * z(4) + dgdr(5) * z(5) + dgdr(6) * z(6) + dgdr(7) * z(7) + dgdr(8) * z(8);
            d32 = dgds(1) * z(1) + dgds(2) * z(2) + dgds(3) * z(3) + dgds(4) * z(4) + dgds(5) * z(5) + dgds(6) * z(6) + dgds(7) * z(7) + dgds(8) * z(8);
            d33 = dgdt(1) * z(1) + dgdt(2) * z(2) + dgdt(3) * z(3) + dgdt(4) * z(4) + dgdt(5) * z(5) + dgdt(6) * z(6) + dgdt(7) * z(7) + dgdt(8) * z(8);
            
            det = d11 * ( d22 * d33 - d32 * d23 ) - d12 * ( d21 * d33 - d31 * d23 ) + d13 * ( d21 * d32 - d31 * d22 );

            % Using Cramer rule for solving matrix equations     
            
            if ( abs( det ) > 1.0E-10 )
                deltar = ( rsx * ( d22 * d33 - d32 * d23 ) - d12 * ( rsy * d33 - rsz * d23 ) + d13 * ( rsy * d32 - rsz * d22 ) ) / det;
                deltas = ( d11 * ( rsy * d33 - rsz * d23 ) - rsx * ( d21 * d33 - d31 * d23 ) + d13 * ( d21 * rsz - d31 * rsy ) ) / det;
                deltat = ( d11 * ( d22 * rsz - d32 * rsy ) - d12 * ( d21 * rsz - d31 * rsy ) + rsx * ( d21 * d32 - d31 * d22 ) ) / det;
            end
            
            r = r + deltar;
            s = s + deltas;
            t = t + deltat;
            
            sum = abs( deltar ) + abs( deltas ) + abs( deltat );
            
            if ( sum > 5. )
                sum = 0.;
            end
            
            g_iter = g_iter + 1;
            
        end
        
        if ( r > -1.01 && r < 1.01 && s > -1.01 && s < 1.01 && t > -1.01 && t < 1.01 )

            found = 1;
            elem = j;
            
            if ( r < -1. )
                r = -1.;
            end
            
            if ( r > 1. )
                r = 1.;
            end
            
            if ( s < -1. )
                s = -1.;
            end
            
            if ( s > 1. )
                s = 1.;
            end
            
            if ( t < -1. )
                t = -1.;
            end
            
            if ( t > 1. )
                t = 1.;
            end
            
            g(1) = ( 1. - r ) * ( 1. - s ) * ( 1. - t ) / 8.;
            g(2) = ( 1. + r ) * ( 1. - s ) * ( 1. - t ) / 8.;
            g(3) = ( 1. + r ) * ( 1. + s ) * ( 1. - t ) / 8.;
            g(4) = ( 1. - r ) * ( 1. + s ) * ( 1. - t ) / 8.;
            g(5) = ( 1. - r ) * ( 1. - s ) * ( 1. + t ) / 8.;
            g(6) = ( 1. + r ) * ( 1. - s ) * ( 1. + t ) / 8.;
            g(7) = ( 1. + r ) * ( 1. + s ) * ( 1. + t ) / 8.;
            g(8) = ( 1. - r ) * ( 1. + s ) * ( 1. + t ) / 8.;
            
            mean_temp = 0.;
            
            for k = 1:8
                index1 = find(Temp(:,1) == el_conn(k));
                mean_temp = mean_temp + g(k)*Temp(index1, 2);
            end
            
            tp1(kk) = mean_temp;
            j = 0;
            
        end
        
        j = j + 1;
        
    end
    
    if ( found == 0 )
        % When no element is found that contains the transfer point 
        
       XX = ['For beam node ', num2str(inode),', the transfer point ' num2str(kk), ' was not found within any solid element.'];
       disp(XX);     
       
    end
    
    
end


function [x_cen, y_cen, z_cen] = centroid(M,N)
 
    el_conn_length = length(M);
    nel = M(el_conn_length,1);
    node = N(:, 1);
    
    j = 1;
    r = 0.;
    s = 0.;
    t = 0.;
    
    connectedData = zeros(8,3);
    
    while ( j <= nel )
    
        index = find(M(:,1) == j);
        
        if isempty(index)
            j = j + 1;
            continue;
        end
        
        for k = 1:8
            x(k)=0.;
            y(k)=0.;
            z(k)=0.;
        end
        
        connectedData = zeros(8,3);
        
        for k = 1:8
            connectedNode = M(index,k+1);
            nindex = find(N(:,1) == connectedNode);
            connectedData(k,1:3) = N(nindex,2:4);
            x(k) = connectedData(k,1);
            y(k) = connectedData(k,2);
            z(k) = connectedData(k,3);

        end
        
% shape functions of hexahedral elements are used in natural coordinates 

        g(1) = ( 1. - r ) * ( 1. - s ) * ( 1. - t ) / 8.;
        g(2) = ( 1. + r ) * ( 1. - s ) * ( 1. - t ) / 8.;
        g(3) = ( 1. + r ) * ( 1. + s ) * ( 1. - t ) / 8.;
        g(4) = ( 1. - r ) * ( 1. + s ) * ( 1. - t ) / 8.;
        g(5) = ( 1. - r ) * ( 1. - s ) * ( 1. + t ) / 8.;
        g(6) = ( 1. + r ) * ( 1. - s ) * ( 1. + t ) / 8.;
        g(7) = ( 1. + r ) * ( 1. + s ) * ( 1. + t ) / 8.;
        g(8) = ( 1. - r ) * ( 1. + s ) * ( 1. + t ) / 8.;
        
        x_cen(j) = 0.;
        y_cen(j) = 0.;
        z_cen(j) = 0.;
        
        for k = 1:8
            x_cen(j) = x_cen(j) + g(k)*x(k);
            y_cen(j) = y_cen(j) + g(k)*y(k);
            z_cen(j) = z_cen(j) + g(k)*z(k);
        end
        
        j = j + 1;
        
    end
