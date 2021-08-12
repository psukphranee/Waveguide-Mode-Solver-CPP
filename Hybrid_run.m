%We can use TEM to calculate E_z and H_z, from there we can find the
%transverse components using 9.1.21. 
%Recall that Beta can be extracted from the eigenvalue by dividing by
%delta-x

%E_z and H_z has the same form. We can use the resulting modes found from
%TEM.m to calculate transverse components. The H_z and E_z can be taken from e_vects.
%Unwind columns in e_vects using create_U. dont for get Beta with e_vals.

omega = k_0 * 3 * 10^8;
mu_0 = 4*pi*10^(-7);

%Find betas
beta_matrix = sqrt(e_vals ./ (dx^2));
beta_matrix_backup = beta_matrix;

t = 1;
z = 1;

%Choose a column for H_z and the same for E_z. They must be the same
%because columns correspond to propogation constant beta
i = length(e_vals)-2; %column number
for(t=0:.02:100)
    t
H_z = create_U(e_vects(:,i), M, N);
%E_z = create_U(e_vects(:,i), M, N);
E_z = H_z;
beta = e_vals(i);
n_effective = n_eff(i); %n_eff is a matrix

%k_c should be a matrix
k_c_2 = ((k_0)^2 * (n .^ 2)) - beta^2;

%Use MATLAB's gradient() function to get partials 
[dx_Hz, dy_Hz] = gradient(E_z);
dx_Ez = dx_Hz;
dy_Ez = dy_Hz;

%impedance TE, TM. should be a matrix eqn 9.1.21
eta_TE = (omega * mu_0 / beta) * ones(M, N);
eta_TM = (n.^2);
eta_TM = 1 ./ eta_TM;
eta_TM = (beta / omega) * eta_TM;
% But we use 1/eta_TM
eta_TM_inv = (beta / omega) .* (n .^ 2);

%need to multiple the matrices by the exponential factor 
exp_factor = exp((1i) * (omega * t - beta * z));

%Calculate 9.1.21
leading_factor = -(1i)*beta ./ k_c_2;
Ex = leading_factor .* (dx_Ez + eta_TE .* dy_Hz) * exp_factor;
Ey = leading_factor .* (dy_Ez - eta_TE .* dx_Hz) * exp_factor;
Hx = leading_factor .* (dx_Hz - eta_TM_inv .* dy_Ez) * exp_factor;
Hy = leading_factor .* (dy_Hz + eta_TM_inv .* dx_Ez) * exp_factor;





quiver(Ex, Ey, 6);
hold on
contour(H_z);
info = sprintf('Center: %.2E W x %.2E L | Eigenvalue: %.4f | n effective: %.4f | exp_factor = %.3f', d_center_horizontal, d_center_vertical, e_vals(i), n_eff(i), exp_factor);
title(info);
pause(.2);
clf;
end






