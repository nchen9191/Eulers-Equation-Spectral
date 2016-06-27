%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precompute green functions for tau calculation
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [green_cell, TH_inv_cell] = Gen_green_euler(Nx, Nz, Lx, Lz, mode_x, mode_z, D)

    global I r1 r2 zero g1 g2 g3 g4 D2

    green_cell = cell(5,4);
    
    %Green's functions
    G1_3D = zero;
    G2_3D = zero;
    G3_3D = zero;
    G4_3D = zero;

    %gradient of green's functions
    dG1_x_3D = zero; dG2_x_3D = zero; dG3_x_3D = zero; dG4_x_3D = zero;
    dG1_y_3D = zero; dG2_y_3D = zero; dG3_y_3D = zero; dG4_y_3D = zero;
    dG1_z_3D = zero; dG2_z_3D = zero; dG3_z_3D = zero; dG4_z_3D = zero;
    
    %Laplacian of green's functions
    lap_G1_3D = zero;
    lap_G2_3D = zero;
    lap_G3_3D = zero;
    lap_G4_3D = zero;

    %Pre comupte inverses of laplacian matrices
    TH_inv_cell = cell(Nx,Nz);

    for nx = 2:Nx
        kx = 2*pi*mode_x(nx)/Lx;
        for nz = 1
            kz = 2*pi*mode_z(nz)/Lz;
            TH = D2 - (kx^2+kz^2)*I;
            TH(end-1,:) = r1;       %Neumann Boundary conditions
            TH(end,:) = r2;

            TH_inv = TH^-1;
            TH_inv_cell(nx,nz) = {TH_inv};  

            G1 = TH_inv*g1;
            G2 = TH_inv*g2;
            G3 = TH_inv*g3;
            G4 = TH_inv*g4;

            dG1_x = 1j*kx*G1; dG2_x = 1j*kx*G2; dG3_x = 1j*kx*G3; dG4_x = 1j*kx*G4;
            dG1_y = D*G1;     dG2_y = D*G2;     dG3_y = D*G3;     dG4_y = D*G4;
            dG1_z = 1j*kz*G1; dG2_z = 1j*kz*G2; dG3_z = 1j*kz*G3; dG4_z = 1j*kz*G4;

            lap_G1 = (D2 - (kx^2+kz^2)*I)*G1;
            lap_G2 = (D2 - (kx^2+kz^2)*I)*G2;
            lap_G3 = (D2 - (kx^2+kz^2)*I)*G3;
            lap_G4 = (D2 - (kx^2+kz^2)*I)*G4;      

            G1_3D(nx,:,nz) = G1.';
            G2_3D(nx,:,nz) = G2.';
            G3_3D(nx,:,nz) = G3.';
            G4_3D(nx,:,nz) = G4.';
            
            dG1_x_3D(nx,:,nz) = dG1_x.'; dG1_y_3D(nx,:,nz) = dG1_y.'; dG1_z_3D(nx,:,nz) = dG1_z.';
            dG2_x_3D(nx,:,nz) = dG2_x.'; dG2_y_3D(nx,:,nz) = dG2_y.'; dG2_z_3D(nx,:,nz) = dG2_z.';
            dG3_x_3D(nx,:,nz) = dG3_x.'; dG3_y_3D(nx,:,nz) = dG3_y.'; dG3_z_3D(nx,:,nz) = dG3_z.';
            dG4_x_3D(nx,:,nz) = dG4_x.'; dG4_y_3D(nx,:,nz) = dG4_y.'; dG4_z_3D(nx,:,nz) = dG4_z.';

            lap_G1_3D(nx,:,nz) = lap_G1.';
            lap_G2_3D(nx,:,nz) = lap_G2.';
            lap_G3_3D(nx,:,nz) = lap_G3.';
            lap_G4_3D(nx,:,nz) = lap_G4.';

        end
    end

    for nx = 1:Nx
        kx = 2*pi*mode_x(nx)/Lx;
        for nz = 2:Nz
            kz = 2*pi*mode_z(nz)/Lz;
            TH = D2 - (kx^2+kz^2)*I;
            TH(end-1,:) = r1;
            TH(end,:) = r2;

            TH_inv = TH^-1;
            TH_inv_cell(nx,nz) = {TH_inv};  

            G1 = TH_inv*g1;
            G2 = TH_inv*g2;
            G3 = TH_inv*g3;
            G4 = TH_inv*g4;

            dG1_x = 1j*kx*G1; dG2_x = 1j*kx*G2; dG3_x = 1j*kx*G3; dG4_x = 1j*kx*G4;
            dG1_y = D*G1;     dG2_y = D*G2;     dG3_y = D*G3;     dG4_y = D*G4;
            dG1_z = 1j*kz*G1; dG2_z = 1j*kz*G2; dG3_z = 1j*kz*G3; dG4_z = 1j*kz*G4;

            lap_G1 = (D2 - (kx^2+kz^2)*I)*G1;
            lap_G2 = (D2 - (kx^2+kz^2)*I)*G2;
            lap_G3 = (D2 - (kx^2+kz^2)*I)*G3;
            lap_G4 = (D2 - (kx^2+kz^2)*I)*G4;      

            G1_3D(nx,:,nz) = G1.';
            G2_3D(nx,:,nz) = G2.';
            G3_3D(nx,:,nz) = G3.';
            G4_3D(nx,:,nz) = G4.';
            
            dG1_x_3D(nx,:,nz) = dG1_x.'; dG1_y_3D(nx,:,nz) = dG1_y.'; dG1_z_3D(nx,:,nz) = dG1_z.';
            dG2_x_3D(nx,:,nz) = dG2_x.'; dG2_y_3D(nx,:,nz) = dG2_y.'; dG2_z_3D(nx,:,nz) = dG2_z.';
            dG3_x_3D(nx,:,nz) = dG3_x.'; dG3_y_3D(nx,:,nz) = dG3_y.'; dG3_z_3D(nx,:,nz) = dG3_z.';
            dG4_x_3D(nx,:,nz) = dG4_x.'; dG4_y_3D(nx,:,nz) = dG4_y.'; dG4_z_3D(nx,:,nz) = dG4_z.';

            lap_G1_3D(nx,:,nz) = lap_G1.';
            lap_G2_3D(nx,:,nz) = lap_G2.';
            lap_G3_3D(nx,:,nz) = lap_G3.';
            lap_G4_3D(nx,:,nz) = lap_G4.';

        end
    end
    
    green_cell(1,:) = {G1_3D, G2_3D, G3_3D, G4_3D};
    green_cell(2,:) = {dG1_x_3D, dG2_x_3D, dG3_x_3D, dG4_x_3D};
    green_cell(3,:) = {dG1_y_3D, dG2_y_3D, dG3_y_3D, dG4_y_3D};   
    green_cell(4,:) = {dG1_z_3D, dG2_z_3D, dG3_z_3D, dG4_z_3D};
    green_cell(5,:) = {lap_G1_3D, lap_G2_3D, lap_G3_3D, lap_G4_3D};
    
end