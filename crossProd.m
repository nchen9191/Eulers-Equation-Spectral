%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute cross-product of two vectors in 3D arrays
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CPPx, CPPy, CPPz] = crossProd(AxP,AyP,AzP,BxP,ByP,BzP)

    CPPx = AyP.*BzP - AzP.*ByP;
    CPPy = AzP.*BxP - AxP.*BzP;
    CPPz = AxP.*ByP - AyP.*BxP;
    
end