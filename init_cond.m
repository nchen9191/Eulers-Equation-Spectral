function [VxP0, VyP0, VzP0] = init_cond(L,X,Y,Z,mode_x_3D,mode_z_3D,int_case)

    Nx = size(mode_x_3D,1);
    Ny = size(mode_x_3D,2);
    Nz = size(mode_x_3D,3);
    
    Lx = L(1);
    Ly = L(2);
    Lz = L(3);
    
    kx_3D = 2*pi*mode_x_3D/Lx;
    kz_3D = 2*pi*mode_z_3D/Lz;
    
    if(int_case == 1)
    
        %column vortex initial flow
        a = Lx/4;
        b = Lx/9;
        c = Lx/5;
        w0 = 1;
        vort = w0*exp(-(X-Lx/2).^2/a^2 - (Y-Ly/2).^2/b^2 - (Z - Lz/2).^2/c^2);

        psi = -1*ColumnVortexCheb(vort, Lx,Lz, mode_x_3D, mode_z_3D);
        psi_FFC = FFCT(psi);
        VxFFC0 = 1j*kz_3D.*psi_FFC;
        VzFFC0 = -1j*kx_3D.*psi_FFC;
        VyFFC0 = zeros(Nx,Ny,Nz);
        VxP0 = iFFCT(VxFFC0);
        VyP0 = iFFCT(VyFFC0);
        VzP0 = iFFCT(VzFFC0);
        
    elseif(int_case == 2)

        %Kelvin-Helmholtz instability
        %Warning: will become first order accurate because does not satisfy
        %divergence perfectly
        VxP0 = tanh(5*Y);
        VyP0 = 0.1*sin(2*X).*(-exp(-8)+exp(-8*(Y.^2)));
        VzP0 = zeros(Nx,Ny,Nz);
        
    else

        %Random Analytical function that satisfies vy = 0 at boundaries
        %And divergence free condition
        VxP0 = (2*Y).*cos(X);
        VyP0 = sin(X).*(Y.^2-1);
        VzP0 = zero;
    end
    
end