
% This routine maps temperatures from a solid thermal FEA model to a 
% structural beam element model. At, present the data are written out in 
% ANSYS format. This is is the fist of the two routines for mapping
% results for beam elements. In this approach, 15 transfer points are 
% chosen along the beam cross section. Their coordinates are calculated 
% from the cross-sectional dimensions entered. In the second approach,
% the transfer point coordinates are entered in a file by the user.
% Output files for each temperature data sets are created.
%
% At, present output data are written out in 
% ANSYS format. A sample of the output file is shown below for element
% number 21580, which is a 2-node beam element. The temperatures at each
% beam element end takes 3 temperatures, e.g. mean temperature, temperature
% at (x,1,0), and temperature at (x,0,1):
%
% BFE,21580,temp,1,363.55,363.18,349.54,355.16
% BFE,21580,temp,5,354.78,340.17

% Input files:
%
% 1. Structural Model:
% a. beam_nodal_coordinates.csv
% b. beam_dimensions.txt
% c. beam_element_connectivity.csv
%
% 2. Heat Transfer Model:
% a. solid_element_connectivity.csv
% b. solid_nodal_coordinates.csv
%

function beam_script_ansys

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

% !!!!!!!!!!!! Read Structural Model Data !!!!!!!!!!!!!
%
% read axial coordinates of nodes in the beam model in ascending order
% format: node#, x_coordinates
%
% y and z coordinates for all structural nodes should be the same. Only
% x-coordinates for structural nodes will vary.

fname = input('Enter structural model nodal coordinates file name without extension( [.csv extension assumed]): ', 's');
fname1 = strcat(fname, '.csv');
P = importdata(fname1, ',');
%P = importdata('beam_nodal_coordinates.csv', ',');

b_nodes = P(:,1);
x_coord = P(:,2);

nnod = length(b_nodes);
grad = zeros(3,nnod);


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


% Enter beam dimensions and number of steps for temperature data
% nt = number of time steps for which temp data are available
% W2 = flange thickness; W3 = beam height; t2 = flange thickness
% t3 = web thickness
% format: W2, W3, t2, t3, nt

fname4 = input('Enter structural model section dimensions file name without extension( [.csv extension assumed]): ', 's');
fname5 = strcat(fname4, '.csv');

S = importdata(fname5, ',');

%[W2, W3, t2, t3 nt] = textread('beam_dimensions.txt','%f %f %f %f %d',1);
W2 = S(1,1);
W3 = S(1,2);
t2 = S(1,3);
t3 = S(1,4);
nt = S(1,5);

W4 = W3 - 2*t2;
W1 = W2;
t1 = t2;

% t_int is the computed temperatures for 15 transfer points for each of the
% "nnod" nodes in the structural beam model for each of the "nt"
% temperature data sets

% np1 is the number of transfer points for each structural node. In this
% example the number of transfer points is chosen to be 15.

np1 = 15;

t_int = zeros(np1,nnod,nt);

y_base = P(1,3);
z_base = P(1,4);

for i = 1:np1
   y_coord(i) = 0.;
   z_coord(i) = 0.;
end

% define x and y-coordinates for 15 beam section transfer points using 
% data obtained from beam section dimension as entered in 
% "beam_dimensions.txt" file

y_coord(1) = y_base + W2/4+W2/8;
z_coord(1) = z_base + W4/2+t2/2+t2/4;

y_coord(2) = y_base + 0.0;
z_coord(2) = z_base + W4/2+t2/2+t2/4;

y_coord(3) = y_base -(W2/4+W2/8);
z_coord(3) = z_base + W4/2+t2/2+t2/4;

y_coord(4) = y_base + W2/4+W2/8;
z_coord(4) = z_base + W4/2+t2/4;

y_coord(5) = y_base + 0.0;
z_coord(5) = z_base + W4/2+t2/4;

y_coord(6) = y_base -(W2/4+W2/8);
z_coord(6) = z_base + W4/2+t2/4;

y_coord(7) = y_base + 0.0;
z_coord(7) = z_base + 1/2*W4/3+1/2*W4/3;

y_coord(8) = y_base + 0.0;
z_coord(8) = z_base + 0.0;

y_coord(9) = y_base + 0.0;
z_coord(9) = z_base -(1/2*W4/3+1/2*W4/3);

y_coord(10) = y_base + W2/4+W2/8;
z_coord(10) = z_base -(W4/2+t2/4);

y_coord(11) = y_base + 0.0;
z_coord(11) = z_base -(W4/2+t2/4);

y_coord(12) = y_base -(W2/4+W2/8);
z_coord(12) = z_base -(W4/2+t2/4);

y_coord(13) = y_base + W2/4+W2/8;
z_coord(13) = z_base -(W4/2+t2/2+t2/4);

y_coord(14) = y_base + 0.0;
z_coord(14) = z_base -(W4/2+t2/2+t2/4);

y_coord(15) = y_base -(W2/4+W2/8);
z_coord(15) = z_base -(W4/2+t2/2+t2/4);


% Write out transfer point coordinates (y, z)
fw2 = fopen('transfer_coordinates.csv','w' );

for i = 1:np1
        
   fprintf(fw2, '%s\n', [num2str(i) ',' num2str(y_coord(i)) ',' num2str(z_coord(i)) ] );

end
    
fclose(fw2);

% !!!!!!!!!!!! Read Heat Transfer Thermal Model Data !!!!!!!!!!!!!
%
%

% Enter heat transfer FEA model info
% Enter element connectivities: elem#, node1, node2,....,node8 (brick element)

fname6 = input('Enter solid HT FEA element connectivity file name without extension( [.csv extension assumed]): ', 's');
fname7 = strcat(fname6, '.csv');
M = importdata(fname7, ',');

%M = importdata('solid_element_connectivity.csv', ',');

el_conn_length = length(M);
nel = M(el_conn_length,1);


%Enter nodal coordiates for the HT model: node #, x, y, z
fname8 = input('Enter solid HT FEA nodal coordinates file name without extension( [.csv extension assumed]): ', 's');
fname9 = strcat(fname8, '.csv');
N = importdata(fname9, ',');

%N = importdata('solid_nodal_coordinates.csv', ',');

node = N(:, 1);
xx = N(:, 2);
yy = N(:, 3);
zz = N(:, 4);

%Enter nodal temperatures from HT model for all the time steps:
%node #, temp for step 1, temp for step 2,........,temp for step nt

fname10 = input('Enter solid HT FEA nodal temperatures file name without extension( [.csv extension assumed]): ', 's');
fname11 = strcat(fname10, '.csv');
Temp1 = importdata(fname11, ',');

%Temp1 = importdata('solid_thermal_temperatures.csv', ',');


%determine centroids of each element

[x_cen, y_cen, z_cen] = centroid(M, N);


for jjj = 1:15
   tp1(jjj) = 0.;
end

for jnode = 2:nt+1
    Temp = [Temp1(:,1),Temp1(:,jnode)];
    
    for inode = 1:nnod
        [tp1]=beamnodes(M,N,node,Temp,x_coord,y_coord,z_coord,inode,x_cen, y_cen, z_cen);
        t_int(:,inode,jnode) = tp1;
    end
    
end


Area = 	W2*t2+W4*t3+W2*t1;
wt(1) = W2*t2/8;
wt(3) = wt(1);
wt(4) = wt(1);
wt(6) = wt(1);
wt(2) = 2*wt(1);
wt(5) = 2*wt(1);
wt(7) = W4*t3/3;
wt(8) = wt(7);
wt(9) = wt(7);

wt(10) = W1*t1/8;
wt(12) = wt(10);
wt(13) = wt(10);
wt(15) = wt(10);
wt(11) = 2*wt(10);
wt(14) = wt(11);

for jnode = 2:nt+1
    % output written to a separate file for each set of temperatures entered
    
    out_filename = num2str(jnode-1);
    fw = fopen([out_filename 'out.csv'],'w' );
    
    for kk1 = 1:3
        for kk2 = 1:nnod
            grad(kk1,kk2) = 0.0;
        end
    end
        
    for inode = 1:nnod
        
        tint1 = t_int(:,inode,jnode);
    
        TY_term1 = (tint1(1)+tint1(3)+tint1(4)+tint1(6)+2*(tint1(2)+tint1(5)))/(4*(W3+W4));
        TY_term2 = (tint1(10)+tint1(12)+tint1(13)+tint1(15)+2*(tint1(11)+tint1(14)))/(4*(W3+W4));
        TY = TY_term1 - TY_term2;

        TX_term1 = 2*t2*(tint1(1)-tint1(3)+tint1(4)-tint1(6))/(3*(W1*t1+W2*t2));
        TX_term2 = 2*t1*(tint1(10)-tint1(12)+tint1(13)-tint1(15))/(3*(W1*t1+W2*t2));
%        TX_term1 = 2*t2*(-tint1(1)+tint1(3)-tint1(4)+tint1(6))/(3*(W1*t1+W2*t2));
%        TX_term2 = 2*t1*(-tint1(10)+tint1(12)-tint1(13)+tint1(15))/(3*(W1*t1+W2*t2));
        TX = TX_term1 + TX_term2;
    
        sum = 0.;

        for i = 1:15
            sum = sum + wt(i)*tint1(i);
        end

        T_mean = sum/Area;
        T_mean_y = T_mean+TY;
        T_mean_x = T_mean+TX;
    
        grad(1,inode) = T_mean;
        grad(2,inode) = T_mean_x;
        grad(3,inode) = T_mean_y;
    end
    
    kk = nnod-1;

    % Write mapped temperature data in ANSYS body force format
    
    for i = 1:kk
        fprintf(fw, '%s\n', ['BFE,' num2str(elems(i)) ',temp,1,'  num2str(grad(1,i)) ',' num2str(grad(2,i)) ',' num2str(grad(3,i)) ',' num2str(grad(1,i+1))]);
        fprintf(fw, '%s\n', ['BFE,' num2str(elems(i)) ',temp,5,' num2str(grad(2,i+1)) ',' num2str(grad(3,i+1))]);
    end
    
    fclose(fw);
    
end



function [tp1] = beamnodes(M,N,node,Temp,xp1,yp1,zp1,inode,x_cen,y_cen,z_cen)

% Finds the hexahedral elements which contain the transfer points.
% Shape functions of hexahedral elements are used in natural coordinates.
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
   
   ind_FE_keep = find(r_FE<= 30 );
   
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
