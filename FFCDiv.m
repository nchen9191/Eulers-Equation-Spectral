function div = FFCDiv(Vx, Vy, Vz, Lx, Lz, mode_x_3D, mode_z_3D, D)

    kx = 2*pi*mode_x_3D/Lx;
    kz = 2*pi*mode_z_3D/Lz;
    
    dVx = 1j*kx.*Vx;
    dVy = ddyCheb(Vy, D);
    dVz = 1j*kz.*Vz;
    
    div = dVx + dVy + dVz;
end