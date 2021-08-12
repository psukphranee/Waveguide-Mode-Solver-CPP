%{
This script constructs the cross section of a an optical waveguide, each
layer is specified by it's refractive index. The cross section will be
represented by a MxN (row x columns) matrix of refractive indices. For
simplification, one of the material principle axes will coincide with the
propogation direction z.
%}
clc
clear
%Define the slab thicknesses in the verticle direction
    d_top = 0.5E-7;
    d_center_vertical = 1.3E-7;
    d_bottom = 2E-7;
    d_vertical_sum = d_top + d_center_vertical + d_bottom;

%Define the width of the middle slab. d_left, d_right are the left and right wall thicknesses, respectively.
    d_center_horizontal = 8E-7;
    d_left = 3E-7;
    d_right = 3E-7;
    d_horizontal_sum = d_left + d_right + d_center_horizontal;

%Resolution value, res, is a scaling factor of how many entries will be in the matrix of our waveguide
    res = 4;

%Light wavelength in free space
    lambda_0 = 4E-7;
    k_0 = (2*pi)/lambda_0; % Incident Wave Number
    omega = k_0 * 3 * 10^8;
    mu_0 = 4*pi*10^(-7);

%refractive index for each layer
    n_top = 3.3;
    n_mid = 2.7 ;
    n_bottom = 1.2;
    n_left = 2.7;
    n_right = 3.3;

%Need to normalize thickness of dielectrics to get ratio for calculating number of points
    d_min = min([d_top, d_center_vertical, d_bottom, d_left, d_right, d_center_horizontal]);
    order_d_min = 10^(abs(order(d_min)));

%Different method of scaling layer thickness to matrix points; 4/24/18: I
%dont know what this comment means. 
    layer_top = d_top * order_d_min * res;
    layer_center_vertical = d_center_vertical * order_d_min * res;
    layer_bottom = d_bottom * order_d_min * res;
    layer_sum_vertical = layer_top + layer_center_vertical + layer_bottom;

    layer_left = d_left * order_d_min * res;
    layer_center_horizontal = d_center_horizontal * order_d_min * res;
    layer_right = d_right * order_d_min * res;
    layer_sum_horizontal = layer_left + layer_center_horizontal + layer_right;

%Number of points in waveguide, not including the boundary. 4/28/18: Why
%does it not include the boundary points? Looks like it does.
    L = layer_top + layer_center_vertical + layer_bottom
    W = layer_left + layer_center_horizontal + layer_right

%Create new variables to be consistent with the naming convention of M rows x N columns
    M = L %Vertical gridpoints;
    N = W %Horizontal gridpoints;

%Refractive index array, contains refractive index for each of the L x W points (L rows, W columns)
    n = zeros(L, W);

%Step size
    dx = d_horizontal_sum/res; %Take the smallest thickness and divide by resolution

%Initialize array of refractive indices n
    progress_bar = waitbar(0, sprintf('Refractive Index Array : Vertical'));
    disp('Initializing Refractive Index array, n')
    for k = 1:L
        if k <= layer_top
            for i = 1:W
                n(k,i) = n_top;
            end
        elseif k <= layer_top + layer_center_vertical
            for i = 1:W
                n(k,i) = n_mid;
            end
         else
            for i = 1:W
                n(k,i) = n_bottom;
            end
        end

        waitbar(k/L, progress_bar);    
    end

    pause(.5);
    waitbar(0, progress_bar, sprintf('Refractive Index Array : Horizontal'));
    pause(.5);

    for k = (layer_top+1):(layer_top+layer_center_vertical) %iterate through rows
        for i = 1:W
            if i <= layer_left
                n(k,i) = n_left;
            elseif i <= layer_left + layer_center_horizontal 
                0;
            else
                n(k,i) = n_right;
            end
        end
        waitbar(i/W, progress_bar);
    end
    close(progress_bar)

%Initialize the electric permitivitty tensor. Refer to Fallahkair/Murphy
%equation 1 for entry labels.
    exx = 1;
    exy = 0;
    eyx = 0;
    eyy = 1;
    ezz = 1;
    e_tensor = [exx, exy, 0; 
                eyx, eyy, 0;
                0  , 0  , ezz];
    e_tensor = e_tensor * (8.85E-12);

%Refer to notebook as to why we take the inverse and consider the diagonal
%elements. 
    e_tilde = inv(e_tensor);
    a1 = e_tilde(1,1);
    a2 = e_tilde(2,2);
    a3 = e_tilde(3,3);
%Consolidate constants, refer to notes
    k1 = a3/a2;
    k2 = (1 - a3/a2);
    k3 = omega^2*mu_0/a2;
    k = k3 - 2*k2 -2;
    
    l1 = a3/a1;
    l2 = (1 - a3/a1);
    l3 = omega^2*mu_0/a1;
    l = l3 - 2*l1 - 2;

%Create eigenmatrix for H, the transverse magnetic field
disp('Creating EigenMatrix');
eigenMatrix_dimension = M*N;
eigenMatrix = eye(2*M*N); %Size 2MN squared, Hx and Hy are each MN long

%Start by initializing the relationship for Hx. We start at (i,j) = (1,2)
%and stop when p=M(j-1) + i reaches the limit. See notes. 
Hx_iteration_limit = M*N - M - 1;
progress_bar = waitbar(0, sprintf('eigenMatrix entries for Hx'));
for i=1:M
    for j=2:N
        p = M*(j-1) + i;
        
        if p > Hx_iteration_limit
            break;
        end    
        eigenMatrix(p,p-M) = k1;
        eigenMatrix(p,p-1) = 1;
        eigenMatrix(p,p) = k;
        eigenMatrix(p,p+1) = 1;
        eigenMatrix(p,p+M) = k1;
        eigenMatrix(p,p+M*N+M+1) = k2/4;
        eigenMatrix(p,p+M*N-M-1) = k2/4;
        eigenMatrix(p,p+M*N-M+1) = k2/4;
        eigenMatrix(p,p+M*N+M-1) = k2/4;
        
        waitbar((M*i + j)/(M*N), progress_bar);
    end
end
waitbar(100, progress_bar);
close(progress_bar);

Hy_iteration_limit = M*N - M;
progress_bar = waitbar(0, sprintf('eigenMatrix entries for Hy'));
for i=2:M
    for j=2:N
        p = M*(j-1) + i;
        
        if p > Hy_iteration_limit
            break;
        end
        eigenMatrix(p+M*N,p+M*N-M) = 1;
        eigenMatrix(p+M*N,p+M*N-1) = l1;
        eigenMatrix(p+M*N,p+M*N) = l;
        eigenMatrix(p+M*N,p+M*N+1) = l1;
        eigenMatrix(p+M*N,p+M*N+M) = 1;
        eigenMatrix(p+M*N,p+M+1) = l2/4;
        eigenMatrix(p+M*N,p-M-1) = l2/4;
        eigenMatrix(p+M*N,p-M+1) = l2/4;
        eigenMatrix(p+M*N,p+M-1) = l2/4;
        
        waitbar((M*i + j)/(M*N), progress_bar);
    end
end
waitbar(100, progress_bar);
close(progress_bar);

disp('Finding eigenvalues/eigenvectors')
[eigenvectors, eigenvalues] = eig(eigenMatrix);
eigenvalues = transpose(diag(eigenvalues));
e_vects_back_up = eigenvectors;
e_vals_back_up = eigenvalues;
pos_e_vals = find(eigenvalues > 0);
fprintf('Amount of positive e_vals: %d\n', length(pos_e_vals));
disp('Done finding eigenvectors/eigenvalues')
