% uses thomas' algorithm to solve c (more efficient than V1)

function [S, coefficients] = splineInterpolationCubicNaturalV2(xData, yData, domain)
    
    % check that xData and yData are the same length
    if length(xData) ~= length(yData)
        fprintf('error: xData and yData must be the same length\n');
        return;
    end

    % number of splines
    n = length(xData) - 1;

    % for Si(x), calculate then store ai
    a = zeros(1, n+1);    % need to store an+1 to calculate bn and cn
    for i = 1 : n+1
        a(i) = yData(i);
    end

    % for Si(x), calculate then store hi
    h = zeros(1, n);
    for i = 1 : n
        h(i) = xData(i+1) - xData(i);
    end

    % for Si(x), calculate then store ci
        alpha = [h(1 : end-1), 0];
        gamma = [0, h(1 : end)];
        beta = ones(1, n+1);
        for i = 2 : n
            beta(i) = 2 * (h(i-1) + h(i));
        end

        Q = zeros(1, n+1);
        for i = 2 : n
            Q(i) = 3/h(i) * (a(i+1) - a(i)) + 3/h(i-1) * (a(i-1) - a(i));
        end
    
        % [A]{c} = {C}
        c = (tridiagonalSolverVectors(alpha, beta, gamma, Q))';

    % for Si(x), calculate then store bi
    b = zeros(1, n);
    for i = 1 : n
        b(i) = 1/h(i) * (a(i+1) - a(i)) - h(i)/3 * (2*c(i) + c(i+1));
    end

    % for Si(x), calculate then store di
    d = zeros(1, n);
    for i = 1 : n
        d(i) = (c(i+1) - c(i))/(3*h(i));
    end

    % find the leftmost spline over the given domain
    % S(x) will be calculated starting from this spline
    currentSpline = 1;
    for i = 2 : n
        if domain(1) < xData(i)
            currentSpline = i - 1;
            break;
        end
    end

    % calculate S(x)
    S = zeros(1, length(domain));
    for i = 1 : length(domain)

        % move onto the next spline if we encounter the next x data point
        % if we are on the last spline, we do not move since there are no more splines
        if domain(i) >= xData(currentSpline + 1) && currentSpline < n 
            currentSpline = currentSpline + 1;
        end

        x = domain(i) - xData(currentSpline);
        ai = a(currentSpline);
        bi = b(currentSpline);
        ci = c(currentSpline);
        di = d(currentSpline);
        S(i) = ai + bi*x + ci*x^2 + di*x^3;

    end

    coefficients = [a(1 : end-1)', b', c(1 : end-1)', d']
    matlabCoefficients = flip(csape(xData,yData,'variational').coefs, 2)

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