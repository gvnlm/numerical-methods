function X = diagonalSolver(A, C)
    
    X = C./diag(A);   % error if A(i, i) = 0

end