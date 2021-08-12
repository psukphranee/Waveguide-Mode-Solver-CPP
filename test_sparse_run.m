%graph test
make_x = zeros(M,N); 
make_y = make_x;

for i=1:10
    
    pick = v(:,i);
    make_x(:) = pick(1:M*N);
    make_y(:) = pick((M*N)+1:M*2*N);
    
    g = figure(1);
    %subplot(2,1,1);
    mesh(make_x);
    %subplot(2,1,2);
    h = figure(2);
    mesh(transpose(eps));
    
    pause(1);
    
    clf(g);
end
