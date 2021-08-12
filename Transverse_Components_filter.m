%We need to filter out the results to H having continous parallel
%components at the boundaries. The idea is to sweep the refractive index
%matrix (ie computational window matrix) horizontally (and then vertically)
%to see which indices have different adjacent indices. 

%for the horizontal sweep, we'll just subract the right columns from left
%columns of n and store it in the corresponding left column of
%horizontal_sweep. 
horizontal_sweep = zeros(M,N);
for i=1:N-1
    horizontal_sweep(:,i) = n(:,i+1) - n(:,i);
end

vertical_sweep = zeros(M,N);
for i=1:M-1
    vertical_sweep(i,:) = n(i+1,:) - n(i,:);
end

%We now go through the sweep matrices and see which entries are non-zero.
%from there we find out which entries in the eigenvector columns need to be
%"equal". Equal meaning the difference being within some epsilon. 
%We do this by finding non-zero entries and converting those entries into
%the corresponding entry 'p' of the [Hx Hy] vector. So a nonzero entry at
%(i,j) will now have an entry of M*(j-1) + i.
for i_column=1:N
    for i_row=1:M
        if horizontal_sweep(i_row,i_column) ~= 0
            horizontal_sweep(i_row,i_column) = M*(i_column-1) + i_row;
        end
        if vertical_sweep(i_row,i_column) ~= 0
            vertical_sweep(i_row,i_column) = M*(i_column-1) + i_row;
        end
    end
end
%We then reshape the above into a vector and delete all zero entries.
horizontal_sweep = horizontal_sweep(:);
horizontal_sweep = nonzeros(horizontal_sweep);
vertical_sweep = vertical_sweep(:);
vertical_sweep = nonzeros(vertical_sweep);

epsilon = .01;

%I want to find out which eigenvectors dont have
%continuous parallel boundary components. We need to use the entries given 
%in the "sweep" matrices to determine which entries in H need to be within
%epsilon of its neighbor. We'll save the original eigenvector matrix and
%make a copy to modify.

eigenvectors_copy = eigenvectors;

%We use the horizontal sweep entries to filter continuity of Hy at
%boundaries. Note that boundary points between two horizontal entries are M
%elements apart.
progress_bar = waitbar(0, sprintf('Filtering non continuous boundaries for Hy'));
[rows, columns] = size(eigenvectors_copy);
%for current_H=1:columns %iterate through columns, which are the H=[Hx Hy] arrays
current_column = 1;
while(current_column < columns)
    for i=1:length(horizontal_sweep) %now we iterate through the iterate through the "sweep" array entries to see which H(p) to compare
        %For the horizontal sweep, horizontally adjacent entries are M
        %entries apart in the H array
        p = horizontal_sweep(i);
        sprintf('Current Column: %d, H(%d) = %.3f' , current_column, p, eigenvectors_copy(p,current_column))

        if(abs(eigenvectors_copy(p,current_column) - eigenvectors_copy(p + M, current_column)) < epsilon)
            eigenvectors_copy(:,current_column) = [];
            current_column = current_column - 1;
            break;
        end
    end
    current_column = current_column + 1;
    [rows, columns] = size(eigenvectors_copy);
    waitbar(current_column/columns, progress_bar);
end
close(progress_bar);


%Perform horizontal component boundary continuity condition using the vertical sweep. Vertically adjacent Hx
%entries need to be continuous. They are 1 entry apart in the H array
progress_bar = waitbar(0, sprintf('Filtering non continuous boundaries for Hx'));
current_column = 1;
[rows, columns] = size(eigenvectors_copy);
while(current_column < columns)
%for current_H=1:columns %iterate through columns, which are the H=[Hx Hy] arrays
    for i=1:length(vertical_sweep) %now we iterate through the iterate through the "sweep" array entries to see which H(p) to compare
        %For the vertical sweep, vertically adjacent entries are 1
        %entries apart in the H array
        p = vertical_sweep(i);
        if(abs(eigenvectors_copy(p,current_column) - eigenvectors_copy(p + 1,current_column)) < epsilon)
            eigenvectors_copy(:,current_column) = [];
            current_column = current_column - 1;
            break;
        end
    end
    current_column = current_column + 1;
    waitbar(current_column/columns, progress_bar);
end
close(progress_bar);

