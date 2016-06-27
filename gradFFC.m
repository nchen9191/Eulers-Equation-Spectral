%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier-Fourier-Chebyshev calculation of gradient
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gradx, grady, gradz] = gradFFC(V, Lx, Lz, mode_x_3D, mode_z_3D, D)

    kx = 2*pi*mode_x_3D/Lx;
    kz = 2*pi*mode_z_3D/Lz;
    
    gradx = 1j*kx.*V;
    gradz = 1j*kz.*V;
    grady = ddyCheb(V, D);
end