%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix to compute tau's with pressure boundary conditions built in
% analytically
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = tau_mat3(lap_G1,lap_G2)

    global o e 

    B = zeros(2,2);
    N = length(lap_G1)-1;
    
    %div(V) Nz mode
    B(2,1) = -lap_G1(end);
    B(2,2) = -lap_G2(end);
    
    %div(V) Nz-1 mode
    B(1,1) = -lap_G1(end-1);
    B(1,2) = 2*(N)-lap_G2(end-1);
    
end