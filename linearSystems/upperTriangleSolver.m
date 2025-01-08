function X = upperTriangleSolver(A, C)
    
    % check that A is a square matrix
    [rows, cols] = size(A);
    if rows ~= cols
        fprintf('error: A must be a square matrix\n');
        X = NaN;
        return;
    else
        n = rows;
    end

    % initialise X
    X = zeros(n, 1);

    % calculate xn
    X(n) = C(n)/A(n, n);

    % calculate xn-1 to x1
    for i = n-1 : -1 : 1
        
        sum = 0;
        for j = i + 1 : n
            sum = sum + A(i, j)*X(j);
        end
    
        % ax + sum = c  =>  x = (c - sum)/a 
        X(i) = (C(i) - sum)/A(i, i);  % error if A(i, i) = 0

    end

end