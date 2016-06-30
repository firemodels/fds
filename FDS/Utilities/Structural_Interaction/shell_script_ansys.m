
% This routine maps temperatures from a solid thermal FEA model to a 
% structural shell element model. At, present output data are written out in 
% ANSYS format. A sample of the output file written by this program is shown below for an element
% number 25383, which is a 3-layer shell element with 4 nodes:
%
% BFE,25383,temp,1,179.16,178.9,178.4,178.91
% BFE,25383,temp,5,83.81,83.72,83.54,83.7
% BFE,25383,temp,9,51.31,51.27,51.18,51.26
%
% In this approach, there is one transfer point per structural node, whose
% coordinate is the same as that of the corresponding structural node.
% However, for each structural node, there is one transfer point at the mid
% plane of each shell layer. In other nodes for a structural node A (coordinate x1,
% y1, z1) belonging to a 4-node, 3-layer shell element with a slab thickness of h, 
% there are 3 transfer points whose coordinates are (x1,y1,z1+h/6), 
% (x1,y1,z1+h/2), (x1,y1,z1+5h/6). These are the coordinates of the mid plane
% of each of the three layers.
%
% 
% Input files:
%
% 1. Structural Model:
% a. structural_model_slab_nodal_coord.csv
% b. shell_dimensions.txt
% c. structural_model_slab_elem_conn.csv
%
% 2. Heat Transfer Model:
% a. solid_thermal_slab_elem_conn.csv
% b. solid_thermal_slab_nodal_coord.csv
%

function shell_script_ansys

clear all
close all

s_nodes = [];
z_coord = [];
x_coord = [];
y_coord = [];
tp1 = [];

x_cen = [];
y_cen = [];
z_cen = [];

xb = [];
yb = [];
zb = [];

elb = [];
n1b = [];
n2b = [];
n3b = [];
n4b = [];

eltn = [];
MM = [];

% !!!!!!!!!!!! Read Structural Model Data !!!!!!!!!!!!!
%
% read nodal coordinates of nodes in the shell model in ascending order
% format should be: node#, x_coord, y_coord, z_coord
% shell elements are assumed to be quadrilaterals with 4 nodes
% The height of shell elements should correspond to the height of the
% bottom surface of the slab in the equivalent thermal model

fname = input('Enter structural shell model nodal coordinates file name without extension( [.csv extension assumed]): ', 's');
fname1 = strcat(fname, '.csv');
P = importdata(fname1, ',');

%P = importdata('structural_model_slab_nodal_coord.csv', ',');

s_nodes = P(:,1);
xb = P(:,2);
yb = P(:,3);
zb = P(:,4);

nnod = length(s_nodes);
grad = zeros(3,nnod);


% Enter shell thickness, number of layers, number of transfer points per
% shell node, and number of temperature steps
% nlayer = number of layers in shell element model (typical number is 3 or 4)
%
% npoints = is the number of transfer points per structural shell element node 
% nt = number of time steps for which temp data are available

% Sample values of these parameters are shown
% sl_thick = 4.0   e.g. 4 inch thick slab;
% nt=6             e.g. 6 sets of temperature data;
% npoints = 1      e.g. only 1 set of transfer point per node per layer;
% nlayer = 3       e.g. 3-layer structural shell element;

fname2 = input('Enter structural shell model section dimensions file name without extension( [.csv extension assumed]): ', 's');
fname3 = strcat(fname2, '.csv');
T = importdata(fname3, ',');

%[sl_thick, nlayer, npoints, nt] = textread('shell_dimensions.txt','%f %f %f %d',1);
sl_thick = T(1,1);
nlayer = T(1,2);
npoints  = T(1,3);
nt = T(1,4);

% Read structural (shell) model elem connectivities, elem#, node1, node 2, node 3, node 4
fname4 = input('Enter structural shell model element connecvity file name without extension( [.csv extension assumed]): ', 's');
fname5 = strcat(fname4, '.csv');
Q = importdata(fname5, ',');

%Q = importdata('structural_model_slab_elem_conn.csv', ',');

elb = Q(:,1);
n1b = Q(:,2);
n2b = Q(:,3);
n3b = Q(:,4);
n4b = Q(:,5);

nels = length(elb);

% t_int is the computed temperatures for each transfer point for each of the
% "nnod" nodes in the structural shell model for each of the "nt"
% temperature data sets and "nlayers" number of layers.

t_int = zeros(npoints,nnod,nt,nlayer);


% !!!!!!!!!!!! Read Heat Transfer Thermal Model Data !!!!!!!!!!!!!
%
%
% Enter solid element connectivities for the heat transfer finite element model: 
% elem#, node1, node2,....,node8 (brick element)

fname6 = input('Enter solid HT FEA element connectivity file name without extension( [.csv extension assumed]): ', 's');
fname7 = strcat(fname6, '.csv');
M = importdata(fname7, ',');

%M = importdata('solid_thermal_slab_elem_conn.csv', ',');

el_conn_length = length(M);
nel = M(el_conn_length,1);


%Enter nodal coordinates for the heat transfer FEA model: 
% node#, x, y, z
fname8 = input('Enter solid HT FEA nodal coordinates file name without extension( [.csv extension assumed]): ', 's');
fname9 = strcat(fname8, '.csv');
N = importdata(fname9, ',');

%N = importdata('solid_thermal_slab_nodal_coord.csv', ',');
node = N(:, 1);
xx = N(:, 2);
yy = N(:, 3);
zz = N(:, 4);


%Enter nodal temperatures from HT model for all the time steps:
%node #, temp for step 1, temp for step 2,........,temp for step nt

fname10 = input('Enter solid HT FEA nodal temperature file name without extension( [.csv extension assumed]): ', 's');
fname11 = strcat(fname10, '.csv');
Temp1 = importdata(fname11, ',');

%Temp1 = importdata('solid_thermal_slab_temp.csv', ',');


for iii = 1:el_conn_length
    if (M(iii,1) > 0)
       eltn(iii,:) = iii;
    end
end


eltn = [eltn, M];


%determine centroids of each Heat Transfer FEA element

[x_cen, y_cen, z_cen] = centroid(M, N);


for jjj = 1:npoints
   tp1(jjj) = 0.;
end

% "sli" array is used to compute the height of the mid-plane of each layer
% in the shell model; naturally mid-plane offset is assumed.
% e.g., for a 3-layer element with mid-plane offset, sli values are:
%sli(1) = 1;
%sli(2) = 3;
%sli(3) = 5;

for i = 1:nlayer
    sli(i) = 2*i - 1;
end

zb_base = zb(1);

layer_thick = sl_thick/nlayer;

for klayer = 1:nlayer
%    zb1 = zb_base + sli(klayer)*(sl_thick/6.0);
    zb1 = zb_base + klayer*layer_thick - layer_thick/2.0;
    
    for jnode = 2:nt+1
        
        Temp = [Temp1(:,1),Temp1(:,jnode)];
        
        for inode = 1:nnod
            xcoord = xb(inode);
            ycoord = yb(inode);
            zcoord = zb1;
            [tp1]=shellnodes(M,N,node,Temp,xcoord,ycoord,zcoord,x_cen,y_cen,z_cen,eltn,npoints,inode);
            t_int(:,inode,jnode,klayer) = tp1;
        end
        
    end
    
end

for jnode = 2:nt+1
    % output written to a separate file for each set of temperatures entered
    
    out_filename = num2str(jnode-1);
    fw = fopen([out_filename 'slout.csv'],'w' );
    
    for jj = 1:nels
        k = 1;
        n1 = n1b(jj);
        n2 = n2b(jj);
        n3 = n3b(jj);
        n4 = n4b(jj);
        
        index11 = find( s_nodes == n1 );
        index12 = find( s_nodes == n2 );
        index13 = find( s_nodes == n3 );
        index14 = find( s_nodes == n4 );
       
        for ilayer = 1:nlayer
            tint1 = t_int(:,index11,jnode,ilayer);
            tint2 = t_int(:,index12,jnode,ilayer);
            tint3 = t_int(:,index13,jnode,ilayer);
            tint4 = t_int(:,index14,jnode,ilayer);
            
% Write mapped temperature data in ANSYS body force format          

            fprintf(fw, '%s\n', ['BFE,' num2str(elb(jj)) ',temp,' num2str(k) ','  num2str(tint1), ',' num2str(tint2), ',' num2str(tint3), ',' num2str(tint4) ]);
            k = k + 4;
        end
    end
    
    fclose(fw);
end





function [tp1] = shellnodes(M,N,node,Temp,xp1,yp1,zp1,x_cen,y_cen,z_cen,eltn,npoints,inode)

% Finds the hexahedral elements which contain the transfer points.
% Shape functions of hexahedral elements are used in natural coordinates
% Uses finite element concepts to locate the element that contains the transfer point.
%
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

j = 1;
index = 0;

np = npoints;

for jjj = 1:np
   tp1(jjj) = 0.;
end

for kk = 1:np
    found = 0;
    j = 1;
    
    x_pt = xp1;
    y_pt = yp1;
    z_pt = zp1;
 
    
% select a few elements in the vicinity
   xp = repmat(x_pt,1,el_conn_length);
   yp = repmat(y_pt,1,el_conn_length);
   zp = repmat(z_pt,1,el_conn_length);
    
   r_FE = sqrt( (x_cen - xp).^2 + (y_cen - yp).^2 + (z_cen - zp).^2);
   
   % search for those elements that are within the "r_FE" distance
   % from the transfer point.
   % A default value of 25 is used. This means that those elements whose
   % centroids are located within a distance of 25 from the transfer point
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
        index3 = eltn(index1,2);
        
        if isempty(index1)
            j = j + 1;
            continue;
        end
        
        if isempty(index3)
            j = j + 1;
            continue;
        end
        
        for k = 1:8
            x(k)=0.;
            y(k)=0.;
            z(k)=0.;
        end
        
        for k = 1:8
            connectedNode = eltn(index1,k+2);
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
                index4 = find(Temp(:,1) == el_conn(k));
                mean_temp = mean_temp + g(k)*Temp(index4, 2);
            end
            
            tp1(kk) = mean_temp;
            j = 0;
            
        end
        
        j = j + 1;
        
    end
    
    if ( found == 0 )
        
% When no element is found that contains the transfer point 

       XX = ['For shell node ', num2str(inode),', the transfer point ' num2str(kk), ' was not found within any solid element.'];
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
    
    while ( j <= el_conn_length )
        
        index = M(:,1);
        
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
            connectedNode = M(j,k+1);
            nindex = find(N(:,1) == connectedNode);
            connectedData(k,1:3) = N(nindex,2:4);
            
            x(k) = connectedData(k,1);
            y(k) = connectedData(k,2);
            z(k) = connectedData(k,3);
            
            el_conn(k) = connectedNode;
            
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
    