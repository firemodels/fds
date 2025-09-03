function [fid]=plot_TG_field(t,Nx,Ny,Xc,Yc,u,p)

Ux = reshape(u(:,1), [Nx,Ny])';      % Ny x Nx
Uy = reshape(u(:,2), [Nx,Ny])';
Pmat = reshape(p, [Nx,Ny])';

fid=0; clf

nexttile; contourf(Xc, Yc, Ux, 20, 'LineStyle','none'); axis equal tight;
colorbar; title(['u_x t=' num2str(t)]); xlabel('x'); ylabel('y');

nexttile; contourf(Xc, Yc, Uy, 20, 'LineStyle','none'); axis equal tight;
colorbar; title('u_y'); xlabel('x'); ylabel('y');

nexttile; contourf(Xc, Yc, Pmat, 20, 'LineStyle','none'); axis equal tight;
colorbar; title('pressure p'); xlabel('x'); ylabel('y');

end