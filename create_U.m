%{
Turns a column vector into a grid matrix U
%}
function U = create_U(column_vect, M, N)

U = zeros(M,N);
k = 1;
for j = 2:N-1
    for i = 1:M
        U(i, j) = column_vect(k);
        k = k+1;
    end
end
