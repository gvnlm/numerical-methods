function [S, coefficients] = splineInterpolationQuadraticCiEq0(xData, yData, domain, i)
    
    % check that xData and yData are the same length
    if length(xData) ~= length(yData)
        fprintf('error: xData and yData must be the same length\n');
        return;
    end

    % number of splines
    n = length(xData) - 1;

    if i < 1 || i > n
        fprintf('error: i must be an int such that 1 <= i <= n\n');
        return;
    end

    % for Si(x), calculate then store ai
    a = zeros(1, n + 1);    % need to store an+1 to calculate bn and cn
    for j = 1 : n + 1
        a(j) = yData(j);
    end

    % for Si(x), calculate then store hi
    h = zeros(1, n);
    for j = 1 : n
        h(j) = xData(j + 1) - xData(j);
    end

    % for Si(x), calculate then store ci
    c = zeros(1, n);
    
        % right of Ci = 0
        for j = i+1 : n
            c(j) = (a(j + 1)/h(j) - a(j)*(1/h(j - 1) + 1/h(j)) + a(j - 1)/h(j - 1) - c(j - 1)*h(j - 1))/h(j);
        end
    
        % left of Ci = 0
        for j = i : -1 : 2
            c(j - 1) = 1/h(j - 1) * (a(j + 1)/h(j) - a(j)*(1/h(j - 1) + 1/h(j)) + a(j - 1)/h(j - 1) - c(j)*h(j));
        end

    % for Si(x), calculate then store bi
    b = zeros(1, n);
    for j = 1 : n
        b(j) = (a(j + 1) - a(j))/h(j) - c(j)*h(j);
    end

    % find the leftmost spline over the given domain
    % S(x) will be calculated starting from this spline
    currentSpline = 1;
    for j = 2 : n
        if domain(1) < xData(j)
            currentSpline = j - 1;
            break;
        end
    end

    % calculate S(x)
    S = zeros(1, length(domain));
    for j = 1 : length(domain)

        % move onto the next spline if we encounter the next x data point
        % if we are on the last spline, we do not move since there are no more splines
        if domain(j) >= xData(currentSpline + 1) && currentSpline < n 
            currentSpline = currentSpline + 1;
        end

        x = domain(j) - xData(currentSpline);
        ai = a(currentSpline);
        bi = b(currentSpline);
        ci = c(currentSpline);
        S(j) = ai + bi*x + ci*x^2;

    end

    coefficients = [a(1 : end-1)', b', c']

end