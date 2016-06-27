%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier-Fourier-Chebyshev transform to Physical-Physical-Physical
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lap = laplacianFFC(V, Lx, Ly, mode_x, mode_y, D)

    [gradx, grady, gradz] = gradFFC(V, Lx, Ly, mode_x, mode_y, D);
    lap = FFCDiv(gradx, grady, gradz, Lx, Ly, mode_x, mode_y, D);
    
end
