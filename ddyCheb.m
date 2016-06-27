%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Chebyshev derivative, need to supply D matrix
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dvdy = ddyCheb(A, D)

    Nz = size(A,3);
    dvdy = zeros(size(A));

    for i = 1:Nz
        dvdy(:,:,i) = A(:,:,i)*D.';
    end
end