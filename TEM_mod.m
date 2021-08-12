%%%%%%%%%%%%%%%%%%%%% Modify these values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
%{ 
Modify TEM.m to keep complex eigenvalues. Eqn 9.1.21 suggests that beta
should be complex for transverse components to be non-zero. 
%}

%First we define the slab thicknesses in the verticle direction
d_top = 2E-7;
d_center_vertical = 6E-7;
d_bottom = 2E-7;

d_vertical_sum = d_top + d_center_vertical + d_bottom;

%We then define the width of the middle slab. d_left, d_right are the left and right wall thicknesses, respectively.
d_center_horizontal = 7E-7;
d_left = 1E-7;
d_right = 2E-7;

d_horizontal_sum = d_left + d_right + d_center_horizontal;

%Resolution value, res, is a scaling factor of how many entries will be in the matrix of our waveguide
res = 7;

%Light wavelength in free space
lambda_0 = 700E-9;
k_0 = (2*pi)/lambda_0; % Incident Wave Number
omega = k_0 * 3 * 10^8;
mu_0 = 4*pi*10^(-7);

%refractive index for each layer
n_top = 1.2;
n_mid = 1.8 ;
n_bottom = .5;
n_left = 1.3;
n_right = .8;
% need to normalize thickness of dielectrics to get ratio for calculating number of points
d_min = min([d_top, d_center_vertical, d_bottom, d_left, d_right, d_center_horizontal]);
order_d_min = 10^(abs(order(d_min)));

%Different method of scaling layer thickness to matrix points
layer_top = d_top * order_d_min * res;
layer_center_vertical = d_center_vertical * order_d_min * res;
layer_bottom = d_bottom * order_d_min * res;
layer_sum_vertical = layer_top + layer_center_vertical + layer_bottom;

layer_left = d_left * order_d_min * res;
layer_center_horizontal = d_center_horizontal * order_d_min * res;
layer_right = d_right * order_d_min * res;
layer_sum_horizontal = layer_left + layer_center_horizontal + layer_right;

%Number of points in waveguid, not including the boundary. 
L = layer_top + layer_center_vertical + layer_bottom
W = layer_left + layer_center_horizontal + layer_right

%Create new variables to be consistent with the convention of M rows x N columns
M = L %Vertical gridpoints;
N = W %Horizontal gridpoints;

%Refractive index array, contains refractive index for each of the L x W points (L rows, W columns)
n = zeros(L, W);

%Step size
dx_1 = d_vertical_sum/res; % test to see if it matches dx
dx = d_horizontal_sum/res; %Take the smallest thickness and divide by resolution

%Initialize array n
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
end

%printf(%s, "Refractive index array, n, initialized vertically");

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
end

disp('Creating eigenMatrix')
eigenMatrix_dimension = M*(N-2);
eigenMatrix = zeros(M*(N-2), M*(N-2), 'single'); %the M(N-2) x M(N-2) matrix we need to solve

%first create matrix "alpha" from matrix n, refer to documentation 
disp('Initializing alpha matrix')
alpha = zeros(L,W);
for i = 1:L
    for j = 1:W
        alpha(i,j) = (n(i,j)*dx*k_0)^2 - 4; 
    end
end
        
%Initialize eigenMatrix diagonals %
%%dummy_index_k = 1; %for indexing diagonals from (1,1) to (m(n-2), m(n-2))
disp('Initializing eigenMatrix diagonals')
pause(3)

for j = 2:(N-1) %i and j are swapped ... ends up being wrong why?
    for i = 1:M %these range limits are determined by equations 4 and 5 in the documentation
        diag_index = M*(j-1) + i - M
        eigenMatrix(diag_index, diag_index) = alpha(i,j);
    end
end

%Initialize eigenMatrix off diagonals (by 1)
disp('Initializing eigenMatrix off diagonals, 1')
pause(3)
for i = 1:eigenMatrix_dimension-1
    eigenMatrix(i, i + 1) = 1;
    eigenMatrix(i + 1, i) = 1;
    i;
end
%Initialize eigenMatrix off diagonals (by m)
disp('Initializing eigenMatrix off diagonals, 2')
pause(3)
for i = 1:eigenMatrix_dimension-M
    eigenMatrix(i, i + M) = 1;
    eigenMatrix(i + M, i) = 1;
    i;
end

disp('Finding eigenvalues/eigenvectors')
[e_vects, e_vals] = eig(eigenMatrix);
e_vals = transpose(diag(e_vals));
e_vects_back_up = e_vects;
e_vals_back_up = e_vals;



%{
pos_e_vals = find(e_vals > 0);
fprintf('Amount of positive e_vals: %d\n', length(pos_e_vals));
disp('Done finding eigenvectors/eigenvalues')

%**************************************************************************
%Extract Edge Values of e_vects matrix. 
%Do this by selecting the rows which correspond to the edges
row_length = length(pos_e_vals); %number of columns
edge_values = zeros(1, row_length);

%Extract the upper edge values for all possible matrices U by iterating
%through all the eigenvector rows that correspont to the upper edges.
disp('Getting edge values')
for i = 1:N-2  %There are N-2 columns of U embedded in the columns of the e_vects
    edge_values = [edge_values; e_vects(M*(i - 1) + 1,:)];
end
%Same thing happens as above except with lower edges of all posible matrices U
for i = 1:N-2
    edge_values = [edge_values; e_vects(M*(i - 1) + M,:)];
end

%delete the first row of edge_values because those are just 0's from the
%initial declaration
edge_values(1,:) = [];
edge_values_backup = edge_values;

%This script creates all the variables needed. Proceed to workspace_run.m
%to process the data
disp('TEM.m done executing')
%}