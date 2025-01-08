function X = pointJacobiMethod(A, C, initialGuess, relaxationParameter, ...
    tolerance, maxIterations)

    % check that A is a square matrix
    [rows, cols] = size(A);
    if rows ~= cols
        fprintf('error: A must be a square matrix\n');
        X = NaN;
        return;
    else
        n = rows;
    end
    
    % A = M - N
    % M = diagonal(A)
    M = diag(diag(A));
    N = M - A;

    % predict convergence
    P = M\N;
    if max(abs(eig(P))) < 1
        fprintf('convergence is expected\n');
    else
        fprintf('convergence is not expected\n');
    end

    % initial guess
    Xi = initialGuess.*ones(n, 1);
    residual = A*Xi - C;

    for k = 0 : maxIterations-1

        % display
        disp(table(Xi, residual));
    
        % check if [A]{Xi} ≈ {C}
        if max(abs(residual)) <= tolerance
            fprintf('a solution was found in %d iterations\n', k);
            X = Xi;
            return;
        end
        
        % calculate new guess
        % [M]{Xi} = [N][Xim1] + {C}
        Xim1 = Xi;
        Xi = diagonalSolver(M, N*Xi + C);   % error if A(i, i) = 0

        % no relaxation: w = 1
        % under-relaxation: 0 < w < 1
        % over-relaxation: 1 < w < 2
        Xi = relaxationParameter*Xi + (1 - relaxationParameter)*Xim1;

        % update residual
        residual = A*Xi - C;

    end

    fprintf('no solution was found');
    X = NaN;

end