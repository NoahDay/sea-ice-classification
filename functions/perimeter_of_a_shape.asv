% # Python3 program to find perimeter of area
% # covered by 1 in 2D matrix consists of 0's and 1's.
% https://www.geeksforgeeks.org/find-perimeter-shapes-formed-1s-binary-matrix/
% 
% R = 3
% C = 5
% 
% # Find the number of covered side for mat[i][j].
% def numofneighbour(mat, i, j):
% 
% 	count = 0;
% 
% 	# UP
% 	if (i > 0 and mat[i - 1][j]):
% 		count+= 1;
% 
% 	# LEFT
% 	if (j > 0 and mat[i][j - 1]):
% 		count+= 1;
% 
% 	# DOWN
% 	if (i < R-1 and mat[i + 1][j]):
% 		count+= 1
% 
% 	# RIGHT
% 	if (j < C-1 and mat[i][j + 1]):
% 		count+= 1;
% 
% 	return count;
% 
% # Returns sum of perimeter of shapes formed with 1s
% def findperimeter(mat):
% 
% 	perimeter = 0;
% 
% 	# Traversing the matrix and finding ones to
% 	# calculate their contribution.
% 	for i in range(0, R):
% 		for j in range(0, C):
% 			if (mat[i][j]):
% 				perimeter += (4 - numofneighbour(mat, i, j));
% 
% 	return perimeter;
% 
% # Driver Code
% mat = [ [0, 1, 0, 0, 0],
% 		[1, 1, 1, 0, 0],
% 		[1, 0, 0, 0, 0] ]
% 
% print(findperimeter(mat), end="\n");
% 
% # This code is contributed by Akanksha Rai

R = 3;
C = 5;

mat = [0, 1, 0, 0, 0;
       1, 1, 1, 0, 0;
       1, 0, 0, 0, 0];
find_perimeter(mat,R,C)

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





