%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier-Fourier-Chebyshev calculation of vorticity
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Wx, Wy, Wz] = vorticityFFC(Vx, Vy, Vz, Lx, Lz, D, mode_x_3D, mode_z_3D)

    kx_3D = 2*pi*mode_x_3D/Lx;
    kz_3D = 2*pi*mode_z_3D/Lz;

    dVxdz = 1j*kz_3D.*Vx;
    dVxdy = ddyCheb(Vx, D);
    dVzdx = 1j*Vz.*kx_3D;
    dVzdy = ddyCheb(Vz, D);
    dVydx = 1j*Vy.*kx_3D;
    dVydz = 1j*kz_3D.*Vy;
    
    Wx = dVzdy - dVydz;
    Wy = -dVzdx + dVxdz;
    Wz = dVydx - dVxdy;
    
end