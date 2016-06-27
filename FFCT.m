%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical-Physical-Physical to Fourier-Fourier-Chebyshev transform
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FFC = FFCT(A)

    Nx = size(A,1);
    Nz = size(A,3);

    FFC = fft(fft(A,[],1),[],3)/Nx/Nz;
    
    FFC = ChebTrans(FFC);
    
    if (mod(Nx,2) == 0)
        FFC(Nx/2+1,:,:) = 0;
    end
    
    if (mod(Nz,2) == 0)
        FFC(:,:,Nz/2+1) = 0;
    end
   
end