function x = tridiagonalSolverVectors(a, B, y, Q)

    n = length(Q);
    
    % calculate and store Bi* and Qi*
    BPrime = zeros(1, n);
    QPrime = zeros(1, n);
    BPrime(1) = B(1);
    QPrime(1) = Q(1);
    for i = 2 : n
        BPrime(i) = B(i) - a(i-1)*y(i-1)/BPrime(i-1);
        QPrime(i) = Q(i) - a(i-1)*QPrime(i-1)/BPrime(i-1);
    end

    % calculate and store xi bottom-up (starting from xn)
    x = zeros(n, 1);
    x(n) = QPrime(n)/BPrime(n);
    for i = n-1 : -1 : 1
        x(i) = (QPrime(i) - y(i)*x(i+1))/BPrime(i);
    end

end