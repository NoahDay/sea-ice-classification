function perimeter = find_perimeter(mat,R,C)
    % Returns the sum of perimeter of shapes formed with 1s
    perimeter = 0;
    for i = 1:R
        for j = 1:C
            if mat(i,j)
                perimeter = perimeter + (4 - num_of_neighbour(mat, i, j));
            end

        end
    end
end

function count = num_of_neighbour(mat, i, j)
% Find the number of covered side for mat(i,j)
    count = 0;
    [R, C] = size(mat);
    % Up
    if i > 1 && mat(i-1,j) == 1
        count = count + 1;
    end
    
    % Left
    if j > 1 && mat(i,j-1)
        count = count + 1;
    end

    % Down
    if i < R-1 && mat(i+1,j)
        count = count + 1;
    end

    % Right
    if j < C-1 && mat(i, j+1)
        count = count + 1;
    end

end
