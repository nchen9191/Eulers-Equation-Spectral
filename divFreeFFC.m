%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove divergence from n+1/2 steps by calculating pressure and taus
% Update velocities mode by mode
% Analytically computed pressure boundary conditions in term of taus
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vx_new, vy_new, vz_new] = divFreeFFC(D,kx,kz,nx,nz,...
            vx,vy,vz,div,green_cell,TH_inv)

    global o e g3 g4 I D2 tau_3D128 tau_matrix128
    
    %modify last two elements of divergence
    div_vhalf_mod = [div(1:end-2);o*vy;e*vy];    
        
    %calculate the homogenous and green functions of pi
    pi_h = TH_inv*div_vhalf_mod;

    G1 = green_cell{1,1}(nx,:,nz).'; 
    G2 = green_cell{1,2}(nx,:,nz).'; 
    
    gpi_h_y = D*pi_h;
    
    lap_pih = (D2 - (kx^2+kz^2)*I)*pi_h;    
    
    lap_G1 = green_cell{5,1}(nx,:,nz).';
    lap_G2 = green_cell{5,2}(nx,:,nz).';

    %RHS of matrix equation
    div_v_BC1 = lap_pih(end-1)-div(end-1); %Ny-1 mode
    div_v_BC2 = lap_pih(end)-div(end); %Ny mode
    g = [div_v_BC1; div_v_BC2];
    
    %LHS matrix
    B = tau_mat3(lap_G1,lap_G2);
    tau = B\g; %solve for taus
    
    tau_matrix128{nx,nz} = B;
    tau_3D128{nx,nz} = tau;
    
    gradx_pi = 1j*kx*(pi_h + tau(1)*G1 + tau(2)*G2);
    gradz_pi = 1j*kz*(pi_h + tau(1)*G1 + tau(2)*G2);
    grady_pi = gpi_h_y + D*(tau(1)*G1 + tau(2)*G2);
    
    vx_new = vx - gradx_pi;
    vy_new = vy - grady_pi + tau(1)*g3 + tau(2)*g4;
    vz_new = vz - gradz_pi;
end