%{
Use waveguidemeshfull to construct the matrix of indices
%}

% Refractive indices:
n1 = 3.34;          % Lower cladding
n2 = 3.44;          % Core
n3 = 1.00;          % Upper cladding (air)

% Layer heights:
h1 = 2.0;           % Lower cladding
h2 = 1.3;           % Core thickness
h3 = 0.5;           % Upper cladding
rh = 1.1;           % Ridge height

% Horizontal dimensions:
rw = 1.0;           % Ridge half-width
d = 2.5;            % center-to-center separation
side = 2.5;         % Space on side

% Grid size:
dx = 0.0125;        % grid size (horizontal)
dy = 0.0125;        % grid size (vertical)

lambda = 1.55;      % wavelength
nmodes = 1;         % number of modes to compute

k_0 = (2*pi)/lambda; % Incident Wave Number
omega = k_0 * 3 * 10^8;
mu_0 = 4*pi*10^(-7);
epsilon_0 = 8.85 * 10E-12;
mu = mu_0;

%Construct the mesh of epsilon_xx coefficients. We'll feed it into
%waveguidemesh() just like for the refractive index matrix. For our
%purposes, all parameters will be the same except for [n1, n2, n3]. Note
%that waveguidemesh() squares the given inputs.
epsxx = 1;
epsyy = 1;
epszz = 1;

[x,y,xc,yc,nx,ny,epsxx] = waveguidemeshfull([epsxx, epsxx, epsxx],[h1, h2, h3] ,rh,rw,side,dx,dy);
[x,y,xc,yc,nx,ny,epsyy] = waveguidemeshfull([epsyy, epsyy, epsyy],[h1, h2, h3] ,rh,rw,side,dx,dy);
[x,y,xc,yc,nx,ny,epszz] = waveguidemeshfull([epszz, epszz, epszz],[h1, h2, h3] ,rh,rw,side,dx,dy);
[x,y,xc,yc,nx,ny,eps] = waveguidemeshfull([n1, n2, n3],[h1, h2, h3] ,rh,rw,side,dx,dy);

%Create new variables to be consistent with the naming convention of M rows x N columns
    M = ny; %rows    
    N = nx; %columns

%% Indexing and Coordinates
Hx_coord = zeros(M,N); Hx_coord(:) = (1:M*N);
Hy_coord = zeros(M,N); Hy_coord(:) = (M*N+1:2*M*N);

Hxp = ones(1, M*N); Hxp(:) = Hx_coord(1:M,1:N);
Hyp = ones(1, M*N); Hyp(:) = Hy_coord(1:M,1:N);

Hxw = ones(1, M*(N-1)); Hxw(:) = Hx_coord(1:M, 1:N-1);
Hyw = ones(1, M*(N-1)); Hyw(:) = Hy_coord(1:M, 1:N-1);

Hxe = ones(1, M*(N-1)); Hxe(:) = Hx_coord(1:M, 2:N);
Hye = ones(1, M*(N-1)); Hye(:) = Hy_coord(1:M, 2:N);

Hxn = ones(1, (M-1)*N); Hxn(:) = Hx_coord(1:M-1, 1:N);
Hyn = ones(1, (M-1)*N); Hyn(:) = Hy_coord(1:M-1, 1:N);

Hxs = ones(1, (M-1)*N); Hxs(:) = Hx_coord(2:M, 1:N);
Hys = ones(1, (M-1)*N); Hys(:) = Hy_coord(2:M, 1:N);

diag_coeff = (k_0)^2 ./ eps; diag_coeff = transpose(diag_coeff(:)); diag_coeff = diag_coeff - (4/(dx^2));
od_ew = (1/(dx^2)).*ones(1, M*(N-1)); % off diag east-west
od_ns = (1/(dx^2)).*ones(1, (M-1)*N); % off diag north-south

shift = (n2*k_0)^2;

arg1 = [Hxp, Hxe, Hxw, Hxn, Hxs, Hyp, Hye, Hyw, Hyn, Hys];
arg2 = [Hxp, Hxw, Hxe, Hxs, Hxn, Hyp, Hyw, Hye, Hys, Hyn];
arg3 = [diag_coeff, od_ew, od_ew, od_ns, od_ns, diag_coeff, od_ew, od_ew, od_ns, od_ns];

A = sparse(arg1, arg2, arg3);


options.tol = 1e-7;
options.disp = 0;						% suppress output
options.isreal = isreal(A);

[v,d] = eigs(A, speye(size(A)),10);
%[v,d] = eigs(A,speye(size(A)),20,shift,options);
