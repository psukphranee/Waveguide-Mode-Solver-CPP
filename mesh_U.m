function null = mesh_U(i)
i

V = e_vects(:,i);
U = create_U(V, M, N); %V is vector to "unfold", M and N dimensions of matrix unfolding to
mesh(U);
title(['Column ', num2str(i),' | eigenvalue: ', num2str(e_vals(i)), ' | n effective: ' , num2str(n_eff(i))]);