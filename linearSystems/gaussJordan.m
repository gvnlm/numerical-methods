function [APrime, CPrime] = gaussJordan(A, C)

    % check that A is a square matrix
    [rows, cols] = size(A);
    if rows ~= cols
        fprintf('error: A must be a square matrix\n');
        APrime = NaN;
        CPrime = NaN;
        return;
    else
        n = rows;
    end
    
    % convert A to an upper triangle using gauss elimination
    [A, C] = gaussElimination(A, C);

    % for every diagonal element A(i, i) starting from the bottom right
    for i = n : -1 : 2
        
        % for every element above A(i, i)
        for j = i-1 : -1 : 1
            
            % A(j, i) - factor*A(i, i) = 0
            factor = A(j, i)/A(i, i);

            % row j = row j - factor*row i  
            A(j, :) = A(j, :) - factor*A(i, :);
            C(j) = C(j) - factor*C(i);

        end

    end

    APrime = A;
    CPrime = C;

end

function [APrime, CPrime] = gaussElimination(A, C)
    
    % check that A is a square matrix
    [rows, cols] = size(A);
    if rows ~= cols
        fprintf('error: A must be a square matrix\n');
        APrime = NaN;
        CPrime = NaN;
        return;
    else
        n = rows;
    end

    % for every diagonal element A(i, i) starting from the top left
    for i = 1 : n-1
        
        % for every element below A(i, i)
        for j = i+1 : n

            % A(j, i) - factor*A(i, i) = 0
            factor = A(j, i)/A(i, i);

            % row j = row j - factor*row i  
            A(j, :) = A(j, :) - factor*A(i, :);
            C(j) = C(j) - factor*C(i);

        end
            
    end
    
    APrime = A;
    CPrime = C;

end
 