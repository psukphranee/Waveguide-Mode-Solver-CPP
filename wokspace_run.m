%This script is copied from the end of TEM.m . It is meant to process
%eigenvalues/eigenvectors saved from previous calculations. 
%{
%Reset matrices/arrays with their backup values
e_vects = e_vects_back_up;
e_vals = e_vals_back_up;
n_eff = n_eff_backup;
edge_values = edge_values_backup;
%}
%%
pos_e_vals = find(e_vals > 0);
fprintf('Amount of positive e_vals: %d\n', length(pos_e_vals));
 
%Delete negative eigenvalues and their correstponding eigenvectors
disp('Deleting Negative Eigenvalues/Eigenvectors');
e_vects(:, find(e_vals < 0)) = [];
e_vals(e_vals < 0) = [];

%Sort eigenvalues/eigenvectors
[e_vals , sort_order] = sort(e_vals);
e_vects = e_vects(:, sort_order);
disp('Eigenvectors sorted according to eigenvalues');

%Construct the U matrix using one of the eigenvectors
%Let V be some eigenvector corresponding to a positive eigenvalue
disp('Calculating n_eff');
n_eff = sqrt(e_vals) ./ (k_0 * dx);
%%
%Extract Edge Values of e_vects matrix. 
%Do this by selecting the rows which correspond to the edges
row_length = length(pos_e_vals); %number of columns
edge_values = zeros(1, row_length);

%Extract the upper edge values for all possible matrices U by iterating
%through all the eigenvector rows that correspont to the upper edges.
disp('Extracting Edge Values');
for i = 1:N-2 %There are N-2 columns
    edge_values = [edge_values; e_vects(M*(i - 1) + 1,:)];
end
%Same thing happens as above except with lower edges of all posible matrices U
for i = 1:N-2
    edge_values = [edge_values; e_vects(M*(i - 1) + M,:)];
end

%delete the first row of edge_values because those are just 0's from the
%initial declaration
edge_values(1,:) = [];

%*************************************************************************%
%We want to filter out eigevectors who's edgevalues are too large. We do so
%by summing up the edge values columnwise and delete if they are greater
%than some arbitrary "upper_bound"
%**************************************************************************

max_edge_values = max(abs(edge_values));
upper_bound = 1E-7;

e_vects(:, find(max_edge_values > upper_bound)) = [];
e_vals(:, find(max_edge_values > upper_bound)) = [];
n_eff(:, find(max_edge_values > upper_bound)) = [];

L = length(e_vals);

fprintf('There are %d eigenvectors after filtering.\n', L);

for i=L:-1:(L-100) %Use e_vals length to iterate through the eigenvectors
    V = e_vects(:,i);
    U = create_U(V, M, N); %V is vector to "unfold", M and N dimensions of matrix unfolding to
    %figure(i);
    mesh(U);
    info = sprintf('Center: %.2E W x %.2E L | Eigenvalue: %.4f | n effective: %.4f', d_center_horizontal, d_center_vertical, e_vals(i), n_eff(i));
    title(info);
    filename = sprintf('%d', i)
    %saveas(gcf, filename, 'png');
    pause(.2);
    if i<(2*L/3)
        break;
    end
end