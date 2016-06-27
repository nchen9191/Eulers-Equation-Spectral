%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Chebyshev (cosine) transform
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FFP = ChebTrans(A)

    Nx = size(A,1);
    Ny = size(A,2);
    Nz = size(A,3);

    FFP = A;
    
    N = Ny-1;
    v = zeros(Nx,2*Ny - 2,Nz);
    v(:,1:Ny,:) = FFP(:,end:-1:1,:);
    v(:,Ny+1:end,:)= FFP(:,2:N,:);

    a = fft(v,[],2)/N;
    FFP(:,1,:) = a(:,1,:)/2; 
    FFP(:,2:N,:) = a(:,2:N,:); 
    FFP(:,N+1,:) = a(:,N+1,:)/2;
    
    if (mod(Nx,2) == 0)
        FFP(Nx/2+1,:,:) = 0;
    end
    
    if (mod(Nz,2) == 0)
        FFP(:,:,Nz/2+1) = 0;
    end
   
end