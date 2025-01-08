function X = tridiagonalSolverMatrix(A, C)

    % check that A is a square matrix
    [rows, cols] = size(A);
    if rows ~= cols
        fprintf('error: A must be a square matrix\n');
        return;
    else
        n = rows;
    end

    % for every diagonal element A(i, i) starting from the top left
    for i = 1 : n-1
        
        % convert the element directly below A(i, i) to 0 using row operations
            
            % A(i+1, i) - factor*A(i, i) = 0
            factor = A(i+1, i)/A(i, i);
        
            % row j = row j - factor*row i  
            A(i+1, :) = A(i+1, :) - factor*A(i, :);
            C(i+1) = C(i+1) - factor*C(i);
 
    end

    % solve for X starting from the bottom
    X = zeros(n, 1);
    X(n) = C(n)/A(n, n);
    for i = n-1 : -1 : 1
        X(i) = (C(i) - A(i, i+1)*X(i+1))/A(i, i);
    end

end