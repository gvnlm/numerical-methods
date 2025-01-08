function yPoints = linearDirichletBVP(p, q, r, xPoints, y0, yn)

    % number of points
    n = length(xPoints);

    f = @(xim1, xi, xip1) -p(xi)/(xip1-xim1) + 2/((xip1-xim1)*(xi-xim1));
    g = @(xim1, xi, xip1) -2/((xip1-xi)*(xi-xim1)) + q(xi);
    h = @(xim1, xi, xip1) p(xi)/(xip1-xim1) + 2/((xip1-xim1)*(xip1-xi));

    a = zeros(1, n-1);
    for i = 1 : length(a)-1
        a(i) = f(xPoints(i), xPoints(i+1), xPoints(i+2));
    end

    B = ones(1, n);
    for i = 2 : length(B)-1
        B(i) = g(xPoints(i-1), xPoints(i), xPoints(i+1));
    end

    y = zeros(1, n-1);
    for i = 2 : length(y)
        y(i) = h(xPoints(i-1), xPoints(i), xPoints(i+1));
    end

    Q = zeros(1, n);
    Q(1) = y0;
    Q(n) = yn;
    for i = 2 : n-1
        Q(i) = r(xPoints(i));
    end

    a
    B
    y
    Q

    yPoints = tridiagonalSolverVectors(a, B, y, Q);

end

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