function [vx_new, vy_new, vz_new] = divFreeFFC_test(D,kx,kz,nx,nz,...
            vx,vy,vz,div,green_cell,TH_inv)

    global o e r1 r2 g1 g2 g3 g4 I nu dt D2 tau_3D tau_matrix
    
    %modify last two elements of divergence
    div_vhalf_mod = [div(1:end-2);0;0];
    
    %homogenous part of pi
    pi_h = TH_inv*div_vhalf_mod;
%     pi_h = TH_inv*div_vhalf_mod; G1 = TH_inv*dTnz0; G2 = TH_inv*dTnz1; G3 = TH_inv*a; G4 = TH_inv*b;
    
    %green functions
    G1 = green_cell{1,1}(nx,:,nz).'; 
    G2 = green_cell{1,2}(nx,:,nz).'; 
    G3 = green_cell{1,3}(nx,:,nz).'; 
    G4 = green_cell{1,4}(nx,:,nz).'; 

    %dpi_h/dy
    gpi_h_y = D*pi_h;
    
    %laplacian of pi_h
    lap_pih = (D2 - (kx^2+kz^2)*I)*pi_h;
    
%     G1_x = green_cell{2,1}(nx,:,nz)'; 
%     G2_x = green_cell{2,2}(nx,:,nz)'; 
%     G3_x = green_cell{2,3}(nx,:,nz)'; 
%     G4_x = green_cell{2,4}(nx,:,nz)';
%     
%     G1_y = green_cell{3,1}(nx,:,nz)';     
%     G2_y = green_cell{3,2}(nx,:,nz)';     
%     G3_y = green_cell{3,3}(nx,:,nz)';     
%     G4_y = green_cell{3,4}(nx,:,nz)';
%     
%     G1_z = green_cell{4,1}(nx,:,nz)'; 
%     G2_z = green_cell{4,2}(nx,:,nz)'; 
%     G3_z = green_cell{4,3}(nx,:,nz)'; 
%     G4_z = green_cell{4,4}(nx,:,nz)';
    
    %laplacian of green's function
    lap_G1 = green_cell{5,1}(nx,:,nz).';
    lap_G2 = green_cell{5,2}(nx,:,nz).';
    lap_G3 = green_cell{5,3}(nx,:,nz).';
    lap_G4 = green_cell{5,4}(nx,:,nz).';

    %RHS of matrix equation
    vy_h_BC1 = o*(gpi_h_y - vy); %V(1)
    vy_h_BC2 = e*(gpi_h_y - vy); %V(-1)
    div_v_BC1 = lap_pih(end-1)-div(end-1); %Ny-1 mode
    div_v_BC2 = lap_pih(end)-div(end); %Ny mode
    g = [vy_h_BC1; vy_h_BC2; div_v_BC1; div_v_BC2];
    
    %LHS matrix
    B = tau_mat2(G1,G2,G3,G4,lap_G1,lap_G2,lap_G3,lap_G4,kx,kz,D);
    tau = B\g; %solve for taus
    
    tau_3D{nx,nz} = tau;
    tau_matrix{nx,nz} = B;
    
%     %gradient of Green functions
%     gradx_G1 = 1j*kx*G1; grady_G1 = 1j*ky*G1; gradz_G1 = D*G1;
%     gradx_G2 = 1j*kx*G2; grady_G2 = 1j*ky*G2; gradz_G2 = D*G2;
%     gradx_G3 = 1j*kx*G3; grady_G3 = 1j*ky*G3; gradz_G3 = D*G3;
%     gradx_G4 = 1j*kx*G4; grady_G4 = 1j*ky*G4; gradz_G4 = D*G4;
    
    %gradient of pi
%     gradx_pi = dpih_dx + tau(1)*gradx_G1 + tau(2)*gradx_G2 + tau(3)*gradx_G3 + tau(4)*gradx_G4;
%     grady_pi = dpih_dy + tau(1)*grady_G1 + tau(2)*grady_G2 + tau(3)*grady_G3 + tau(4)*grady_G4;
%     gradz_pi = dpih_dz + tau(1)*gradz_G1 + tau(2)*gradz_G2 + tau(3)*gradz_G3 + tau(4)*gradz_G4;
    
    
    gradx_pi = 1j*kx*(pi_h + tau(1)*G1 + tau(2)*G2 + tau(3)*G3 + tau(4)*G4);
    gradz_pi = 1j*kz*(pi_h + tau(1)*G1 + tau(2)*G2 + tau(3)*G3 + tau(4)*G4);
    grady_pi = gpi_h_y + D*(tau(1)*G1 + tau(2)*G2 + tau(3)*G3 + tau(4)*G4);
    
%     max([max(abs(gradx_pi)), max(abs(grady_pi)), max(abs(gradz_pi))])
%     grady_pi
%     tau;
    
    vx_new = vx - gradx_pi;
    vy_new = vy - grady_pi + tau(1)*g3 + tau(2)*g4;
    vz_new = vz - gradz_pi;
end