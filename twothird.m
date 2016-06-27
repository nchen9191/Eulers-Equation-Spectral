%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2/3 rule to de-alias nonlinear project
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Anew = twothird(A)

    Nx = size(A,1);
    Nz = size(A,3);
    numx = ceil(Nx/3);
    numz = ceil(Nz/3);

    Anew = A;
    
    if(mod(numx,2) == 1)
        numx = numx + 1;
    end
    
    if(mod(numz,2) == 1)
        numz = numz + 1;
    end
    
    Anew(Nx/2+1-numx/2:Nx/2+1+numx/2,:,:) = 0;
    Anew(:,:,Nz/2+1-numz/2:Nz/2+1+numz/2) = 0;
    
end