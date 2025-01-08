% uses gauss elimination to solve c

function [S, coefficients] = splineInterpolationCubicNaturalV1(xData, yData, domain)
    
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
        A = zeros(n+1); 
        A(1, 1) = 1;
        A(end, end) = 1;
        for row = 2 : n
            A(row, row-1) = h(row-1);
            A(row, row) = 2 * (h(row-1) + h(row));
            A(row, row+1) = h(row);
        end
    
        C = zeros(n+1, 1);
        for row = 2 : n
            C(row, 1) = 3/h(row) * (a(row+1) - a(row)) + 3/h(row-1) * (a(row-1) - a(row));
        end

        % [A]{c} = {C}
        [APrime, CPrime] = gaussElimination(A, C);
        c = (upperTriangleSolver(APrime, CPrime))';

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