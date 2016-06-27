%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier-Fourier-Chebyshev transform to Fourier-Fourier-Physical
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FFP = iChebTrans(A)
    
    Nx = size(A,1);
    Ny = size(A,2);
    Nz = size(A,3);
    N = Ny - 1;
    
    if (mod(Nx,2) == 0)
        A(Nx/2+1,:,:) = 0;
    end
    
    if (mod(Nz,2) == 0)
        A(:,:,Nz/2+1) = 0;
    end
    
    A(:,1,:) = A(:,1,:)*2; 
    A(:,2:N,:) = A(:,2:N,:); 
    A(:,end,:) = A(:,end,:)*2;
    
    M = zeros(Nx,2*Ny-2,Nz);
    M(:,1:Ny,:) = A;
    M(:,Ny+1:end,:) = A(:,end-1:-1:2,:);
    v = ifft(M,[],2);
    FFP = v(:,N+1:-1:1,:)*N;
end