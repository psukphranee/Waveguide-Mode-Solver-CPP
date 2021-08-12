%Construct Hx and Hy from some eigenvectors column
[rows, columns] = size(eigenvectors);
[X,Y] = meshgrid(1:N, 1:M);

for column_choice = 1:columns

    Hx = zeros(M,N);
    Hy = zeros(M,N);

    for i=1:M
        for j=1:N
            p_index = M*(j-1) + i;
            Hx(i,j) = eigenvectors(p_index, column_choice);
            Hy(i,j) = eigenvectors(p_index + M*N, column_choice);
        end
    end
    %fig = quiver(X, Y, Hx, Hy)
    %fig = contour(Hx)
    %pause(.2);
    
    figure
    subplot(2,1,1)       % add first plot in 2 x 1 grid
    title('Hx')
    subplot(2,1,2)       % add second plot in 2 x 1 grid
    plot(x,y2,'+')       % plot using + markers
    title('Hy')
end
