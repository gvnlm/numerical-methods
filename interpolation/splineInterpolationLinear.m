function [S, coefficients] = splineInterpolationLinear(xData, yData, domain)

    % check that xData and yData are the same length
    if length(xData) ~= length(yData)
        fprintf('error: xData and yData must be the same length\n');
        return;
    end

    % number of splines
    n = length(xData) - 1;

    % for Si(x), calculate then store ai
    a = zeros(1, n);
    for i = 1 : n
        a(i) = yData(i);
    end

    % for Si(x), calculate then store bi
    b = zeros(1, n);
    for i = 1 : n
        b(i) = (yData(i + 1) - yData(i))/(xData(i + 1) - xData(i));
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

        S(i) = a(currentSpline) + b(currentSpline)*(domain(i) - xData(currentSpline));

    end

    coefficients = [a', b'];

end