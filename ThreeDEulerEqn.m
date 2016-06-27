%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Eulers-Equation solver Periodic in x,z and chebyshev in y
% Spectral Methods - Fourier Series in x,z and Chebyshev Polynomials in y
% Special Tau methods with Green's Functions to enforce divergence free
%
% Nelson Chen 
% University of California, Berkeley
% Computational Fluid Dynamics Lab
% nchen9191@berkeley.edu
% Last revision: 6/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Global Variables initialization
    global n o e r1 r2 g1 g2 g3 g4 I D D2 zero tau_3D tau_matrix

%% initialize Parameters

    %Save parameters
    save_count = 2;
    save_num = 10;

    % Number of modes
    Nx = 128;
    Ny = 129;
    Nz = 4;

    %Domain size
    Lx = 2*pi;
    Lz = 2*pi;
    Ly = 2;
    
    %Time
    T = 8;                    %Final time
    dt = 0.001;                %time step-size
    Nt = floor(T/dt)+1;        %Number of time points
    
    %% Initialize Variables

    %increments
    inc_x = Lx/Nx;
    inc_z = Lz/Nz;
    inc_theta = pi/(Ny-1);
    
    %domain grids
    xpts = 0:inc_x:2*pi-inc_x;
    zpts = 0:inc_z:2*pi-inc_z;
    
    %Chebyshev domain grids, equally spaced in theta
    theta = pi:-inc_theta:0;
    ypts = (cos(theta))*Ly/2;

    %3D points
    [X, Y, Z] = ndgrid(xpts, ypts, zpts);

    %fourier modes in x and z
    mode_x = [0:Nx/2, -Nx/2+1:-1];
    mode_z = [0:Nz/2, -Nz/2+1:-1];
    
    %3D array of modes correspond to x,y,z points
    [mode_x_3D, mode_y_3D, mode_z_3D] = fModes(Nx, Ny, Nz);
    
    %global tau-arrays for diagnostics
    tau_3D = cell(Nx,Nz);
    tau_matrix = cell(Nx,Nz);
    
    %k = 2*pi*n/L, 3D arrays
    kx_3D = 2*pi*mode_x_3D/Lx;
    kz_3D = 2*pi*mode_z_3D/Lz;
    
    VxFFC_2 = zeros(Nx,Ny,Nz);
    VyFFC_2 = zeros(Nx,Ny,Nz);
    VzFFC_2 = zeros(Nx,Ny,Nz);

    Data = cell(floor((Nt-1)/save_num), 4);
    DataFFC = cell(floor((Nt-1)/save_num),4);

    %Chebychev derivative operators
    D = ChDiffnoBC(Ny-1,Ly/2);        %First derivative
    D2 = D^2;                         %Second derivative
    
    %dTN/dy
    g1 = [D(1:end-2,end-1);1;(-1)^(Ny-2)];
    g2 = [D(1:end-2,end);1;(-1)^(Ny-1)];
    
    %identity matrix and zero array for convenience
    I = eye(Ny,Ny);
    zero = zeros(Nx,Ny,Nz);

    %useful vectors for convenience
    g3 = [zeros(Ny-2,1); 1;0];
    g4 = [zeros(Ny-1,1); 1];

    %Boundary vectors Chebyshev 
    n = 0:Ny-1;
    o = ones(1,Ny);
    e = (-o).^n;
    r1 = o*D;
    r2 = e*D;

    %% Green's Functions pre-processing

    [green_cell, TH_inv_cell] = Gen_green_euler(Nx, Nz, Lx, Lz, mode_x, mode_z, D);

    'green'

    %% Initial flows
    
    % case 1: Column Vortex
    % case 2: Shear flow for Kelvin-Helmholtz Instability
    % case 3: Random analytical functions 
    
    int_case = 2;
    
    [VxP0, VyP0, VzP0] = init_cond([Lx,Ly,Lz],X,Y,Z,mode_x_3D,mode_z_3D,int_case);
    
    %Transform velocity field to spectral space
    VxFFC0 = FFCT(VxP0);
    VyFFC0 = FFCT(VyP0);
    VzFFC0 = FFCT(VzP0);
    
    %Save initial condition in data arrays
    DataFFC(1,:) = {VxFFC0, VyFFC0, VzFFC0, 0};
    Data(1,:) = {VxP0, VyP0, VzP0, 0};

    %% First Step Forward Euler
    
    %initialize velocity
    VxFFC = DataFFC{1,1};
    VyFFC = DataFFC{1,2};
    VzFFC = DataFFC{1,3};

    %physical space
    VxP = Data{1,1};
    VyP = Data{1,2};
    VzP = Data{1,3};

    %Calculate vorticity
    [WxFFC, WyFFC, WzFFC] = vorticityFFC(VxFFC, VyFFC, VzFFC, Lx, Lz, D, mode_x_3D, mode_z_3D);
    
    %Remove aliasing from nonlinear projection
    BxFFC = twothird(WxFFC); ByFFC = twothird(WyFFC); BzFFC = twothird(WzFFC);
    AxFFC = twothird(VxFFC); AyFFC = twothird(VyFFC); AzFFC = twothird(VzFFC);
    
    AxP = iFFCT(AxFFC); AyP = iFFCT(AyFFC); AzP = iFFCT(AzFFC);
    BxP = iFFCT(BxFFC); ByP = iFFCT(ByFFC); BzP = iFFCT(BzFFC);
    
    %calculate velocity x vorticity in physical space
    [VxWxP, VxWyP, VxWzP] = crossProd(AxP,AyP,AzP,BxP,ByP,BzP);
    VxWxFFC_1 = FFCT(VxWxP);
    VxWyFFC_1 = FFCT(VxWyP);
    VxWzFFC_1 = FFCT(VxWzP);
    
    %RHS of diff eq at current time-step
    Rx = VxWxFFC_1;
    Ry = VxWyFFC_1;
    Rz = VxWzFFC_1;

    %n+1/2 step
    Vx_1half = VxFFC + dt*Rx;
    Vy_1half = VyFFC + dt*Ry;
    Vz_1half = VzFFC + dt*Rz;
    
    %divergence of half step
    DIV = FFCDiv(Vx_1half, Vy_1half, Vz_1half, Lx, Lz, mode_x_3D, mode_z_3D, D);

    %special treatment for 0-0 modes
    nx = 1; nz = 1;
    VxFFC_2(nx,:,nz) = Vx_1half(nx,:,nz);
    VyFFC_2(nx,:,nz) = 0;
    VzFFC_2(nx,:,nz) = Vz_1half(nx,:,nz);
    
    %calculating pressure and taus to satisfy boundary conditions and
    %divergence-free velocity fields mode by mode
    nz = 1;
    kz = 2*pi*mode_z(nz)/Lz;
    for nx = 2:Nx
        kx = 2*pi*mode_x(nx)/Lx;
        
        %Pressure Laplacian Matrix
        TH_inv = TH_inv_cell{nx,nz};
        
        vx = Vx_1half(nx,:,nz).';
        vy = Vy_1half(nx,:,nz).';
        vz = Vz_1half(nx,:,nz).';
        
        div = DIV(nx,:,nz).';

        [vx_new, vy_new, vz_new] = divFreeFFC(D,kx,kz,nx,nz,...
                vx,vy,vz,div,green_cell,TH_inv);    

        VxFFC_2(nx,:,nz) = vx_new.';
        VyFFC_2(nx,:,nz) = vy_new.';
        VzFFC_2(nx,:,nz) = vz_new.';

    end

    for nx = 1:Nx
        kx = 2*pi*mode_x(nx)/Lx;
        for nz = 2:Nz
            kz = 2*pi*mode_z(nz)/Lz;

            TH_inv = TH_inv_cell{nx,nz};

            vx = Vx_1half(nx,:,nz).';
            vy = Vy_1half(nx,:,nz).';
            vz = Vz_1half(nx,:,nz).';

            div = DIV(nx,:,nz).';

            [vx_new, vy_new, vz_new] = divFreeFFC(D,kx,kz,nx,nz,...
                vx,vy,vz,div,green_cell,TH_inv);

            VxFFC_2(nx,:,nz) = vx_new.';
            VyFFC_2(nx,:,nz) = vy_new.';
            VzFFC_2(nx,:,nz) = vz_new.';

        end
    end
    
    %Convert to real velocities
    VxP_2 = iFFCT(VxFFC_2);
    VyP_2 = iFFCT(VyFFC_2);
    VzP_2 = iFFCT(VzFFC_2);
    
    if (save_num == 1)
        %Save first step
        DataFFC(2,:) = {VxFFC_2, VyFFC_2, VzFFC_2, dt};
        Data(2,:) = {VxP_2, VyP_2, VzP_2, dt};
        save_count = save_count + 1;
    end
    %% Adam-BashForth2 for remaining steps

    for t = 3:Nt

        %initialize velocity
        VxFFC = VxFFC_2;
        VyFFC = VyFFC_2;
        VzFFC = VzFFC_2;

        %physical space
        VxP = VxP_2;
        VyP = VyP_2;
        VzP = VzP_2;

         %Calculate vorticity
        [WxFFC_2, WyFFC_2, WzFFC_2] = vorticityFFC(VxFFC, VyFFC, VzFFC, Lx, Lz, D, mode_x_3D, mode_z_3D);

        %2/3 rule to de-alias nonlinear part
        BxFFC = twothird(WxFFC_2); ByFFC = twothird(WyFFC_2); BzFFC = twothird(WzFFC_2);
        AxFFC = twothird(VxFFC); AyFFC = twothird(VyFFC); AzFFC = twothird(VzFFC);
        AxP = iFFCT(AxFFC); AyP = iFFCT(AyFFC); AzP = iFFCT(AzFFC);
        BxP = iFFCT(BxFFC); ByP = iFFCT(ByFFC); BzP = iFFCT(BzFFC);
        
        %calculate velocity x vorticity in physical space
        [VxWxP_2, VxWyP_2, VxWzP_2] = crossProd(AxP,AyP,AzP,BxP,ByP,BzP);
        VxWxFFC_2 = FFCT(VxWxP_2);
        VxWyFFC_2 = FFCT(VxWyP_2);
        VxWzFFC_2 = FFCT(VxWzP_2);
        
        Rx = 3*VxWxFFC_2 - VxWxFFC_1;
        Ry = 3*VxWyFFC_2 - VxWyFFC_1;
        Rz = 3*VxWzFFC_2 - VxWzFFC_1;

        %n+1/3 step
        Vx_1half = VxFFC + 0.5*dt*Rx;
        Vy_1half = VyFFC + 0.5*dt*Ry;
        Vz_1half = VzFFC + 0.5*dt*Rz;

        DIV = FFCDiv(Vx_1half, Vy_1half, Vz_1half, Lx, Lz, mode_x_3D, mode_z_3D, D);

        VxFFC_2(1,:,1) = Vx_1half(1,:,1);
        VyFFC_2(1,:,1) = 0;
        VzFFC_2(1,:,1) = Vz_1half(1,:,1);

        for nx = 2:Nx
            nz = 1;
            kx = 2*pi*mode_x(nx)/Lx;
            kz = 2*pi*mode_z(nz)/Lz;

            TH_inv = TH_inv_cell{nx,nz};

            vx = Vx_1half(nx,:,nz).';
            vy = Vy_1half(nx,:,nz).';
            vz = Vz_1half(nx,:,nz).';

            div = DIV(nx,:,nz).';
        
            [vx_new, vy_new, vz_new] = divFreeFFC(D,kx,kz,nx,nz,...
                vx,vy,vz,div,green_cell,TH_inv);

            VxFFC_2(nx,:,nz) = vx_new.';
            VyFFC_2(nx,:,nz) = vy_new.';
            VzFFC_2(nx,:,nz) = vz_new.';

        end

        for nx = 1:Nx
            kx = 2*pi*mode_x(nx)/Lx;
            for nz = 2:Nz
                kz = 2*pi*mode_z(nz)/Lz;

                TH_inv = TH_inv_cell{nx,nz};

                vx = Vx_1half(nx,:,nz).';
                vy = Vy_1half(nx,:,nz).';
                vz = Vz_1half(nx,:,nz).';

                div = DIV(nx,:,nz).';

                [vx_new, vy_new, vz_new] = divFreeFFC(D,kx,kz,nx,nz,...
                vx,vy,vz,div,green_cell,TH_inv);
            
                VxFFC_2(nx,:,nz) = vx_new.';
                VyFFC_2(nx,:,nz) = vy_new.';
                VzFFC_2(nx,:,nz) = vz_new.';

            end
        end

        VxP_2 = iFFCT(VxFFC_2);
        VyP_2 = iFFCT(VyFFC_2);
        VzP_2 = iFFCT(VzFFC_2);
        
        %Saving old advection step for Adam-Bashforth
        VxWxFFC_1 = VxWxFFC_2;
        VxWyFFC_1 = VxWyFFC_2;
        VxWzFFC_1 = VxWzFFC_2;       
        
        %Print out time step
        if (mod(t-1,save_num) == 0)
            t
            %saving velocities into cell array
            DataFFC(save_count,:) = {VxFFC_2, VyFFC_2, VzFFC_2, (t-1)*dt};
            Data(save_count,:) = {VxP_2, VyP_2, VzP_2, (t-1)*dt};
            save_count = save_count + 1;
            DIV_diagFFC = FFCDiv(VxFFC_2, VyFFC_2, VzFFC_2, Lx, Lz, mode_x_3D, mode_z_3D, D);
            DIV_diagPPP = iFFCT(DIV_diagFFC);
            max(max(max(abs(DIV_diagPPP))))
        end

    end
