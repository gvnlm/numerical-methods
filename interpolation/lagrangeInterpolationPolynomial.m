function f = lagrangeInterpolationPolynomial(xData, yData, domain)
    
    % check that xData and yData are the same length
    if length(xData) ~= length(yData)
        fprintf('error: xData and yData must be the same length\n');
        return;
    else
        % order of the lagrange interpolation polynomial
        n = length(xData) - 1;
    end

    % calculate fn(x)
    f = zeros(size(domain));
    for i = 0 : n

        % calculate Li(x)
        Li = ones(size(domain));
        for j = 0 : n
            if j ~= i
                Li = Li .* (domain - xData(j + 1))./(xData(i + 1) - xData(j + 1));
            end
        end

        f = f + Li.*yData(i + 1);

    end
    
    % plot
    figure;
    hold on;

    plot(xData, yData, 'ro', 'markerFaceColor', 'r');
    plot(domain, f, 'b');
    
    title('lagrange interpolating polynomial');
    xlabel('x');
    ylabel('y');
    legend('data points', 'fn(x)');
    
end