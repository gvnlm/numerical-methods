function X = lowerTriangleSolver(A, C)
    
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

    % calculate x1
    X(1) = C(1)/A(1, 1);

    % calculate x2 to xn
    for i = 2 : n

        sum = 0;
        for j = 1 : i - 1
            sum = sum + A(i, j)*X(j);
        end
        
        % ax + sum = c  =>  x = (c - sum) / a 
        X(i) = (C(i) - sum)/A(i, i);  % error if A(i, i) = 0

    end

end