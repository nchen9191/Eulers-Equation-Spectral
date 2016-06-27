%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier-Fourier-Chebyshev transform to Physical-Physical-Physical
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FFC = iFFCT(A)
    
    Nx = size(A,1);
    Ny = size(A,2);
    Nz = size(A,3);
    N = Ny - 1;
    
    FFC = iChebTrans(A);
    FFC = ifft(ifft(FFC,[],3),[],1)*Nx*Nz;
end