function [mode_x, mode_y, mode_z] = fModes(Nx, Ny, Nz)

    mode_x = zeros(Nx,Ny,Nz);
    mode_y = mode_x;
    mode_z = mode_x;

    Nx_upper = ceil((Nx-1)/2); %highest modes (goes from -N_lower:N_upper with 2N+1 modes)
    Ny_upper = ceil((Ny-1)/2);
    Nz_upper = ceil((Nz-1)/2);
    
    Nx_lower = floor((Nx-1)/2); 
    Ny_lower = floor((Ny-1)/2);
    Nz_lower = floor((Nz-1)/2);
    
    mode_xvec = [0:Nx_upper, -Nx_lower:-1];
    mode_yvec = [0:Ny_upper, -Ny_lower:-1];
    mode_zvec = [0:Nz_upper, -Nz_lower:-1];
    
    mode_xmat = (ones(Ny,1)*mode_xvec).';
    mode_ymat = ones(Nx,1)*mode_yvec;
    mode_zmat = ones(Ny,1)*mode_zvec;
    
    for i = 1:Nz
        mode_x(:,:,i) = mode_xmat;
        mode_y(:,:,i) = mode_ymat;
    end
    
    for j = 1:Nx
        mode_z(j,:,:) = mode_zmat;
    end
end