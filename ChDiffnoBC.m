%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building the Chebyshev derivative Matrix scaled by L
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = ChDiffnoBC(N,L)

    D = zeros(N+1,N+1);
    for i = 1:2:N
        D = D + diag(i:N,i);
    end
    D(1,:) = D(1,:)/2;
    D = 2*D/L;
end
        
        
        