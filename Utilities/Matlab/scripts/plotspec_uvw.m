% McDermott
% 9-10-10
% plotspec_uvw.m

function [handle]=plotspec_uvw(filename,marker_format)

if ~exist(filename)
    display(['Error: File ' filename ' does not exist. Skipping case.'])
    return
end
M = importdata(filename,',',1);
s = str2num(M.textdata{:});
n = s(2);

% convert to 3D array
p=0;
for k=1:n
    for j=1:n
        for i=1:n
            p=p+1;
            u(i,j,k) = M.data(p,1);
            v(i,j,k) = M.data(p,2);
            w(i,j,k) = M.data(p,3);
            tke(i,j,k) = 0.5*(u(i,j,k)^2+v(i,j,k)^2+w(i,j,k)^2);
        end
    end
end

u_hat = fftn(u)/n^3;
v_hat = fftn(v)/n^3;
w_hat = fftn(w)/n^3;
for k=1:n
    for j=1:n
        for i=1:n
            tke_hat(i,j,k) = 0.5*( u_hat(i,j,k)*conj(u_hat(i,j,k)) + ...
                                   v_hat(i,j,k)*conj(v_hat(i,j,k)) + ...
                                   w_hat(i,j,k)*conj(w_hat(i,j,k)) );
        end
    end
end

% spectrum

L = 9*2*pi/100;
k0 = 2*pi/L;
kmax = n/2;
wn = k0*[0:n]; % wavenumber array
vt = zeros(size(wn));

for kx=1:n
    rkx = kx-1;
    if (kx>kmax+1); rkx=rkx-n; end
    
    for ky=1:n
        rky = ky-1;
        if (ky>kmax+1); rky=rky-n; end
        
        for kz=1:n
            rkz = kz-1;
            if (kz>kmax+1); rkz=rkz-n; end
            
            rk = sqrt(rkx^2+rky^2+rkz^2);
            k = round(rk);
        
            vt(k+1) = vt(k+1) + tke_hat(kx,ky,kz)/k0;
            
        end
    end
end
        
% plot the energy spectrum

handle=loglog(wn(2:n),vt(2:n),marker_format,'MarkerSize',15);
















