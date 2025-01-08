function B = newtonInterpolationCoefficients(xData, yData)
    
    % check that xData and yData are the same length
    if length(xData) ~= length(yData)
        fprintf('error: xData and yData must be the same length\n');
        return;
    else
        % order of newton interpolation polynomial
        n = length(xData) - 1;
    end
    
    % store the divided differences  
    dividedDifferences = zeros(n + 1);

    % calculate col 1
    for row = 1 : n+1
        dividedDifferences(row, 1) = yData(row);
    end

    % calculate col 2, ... col n+1
    for col = 2 : n+1
        for row = 1 : n+2 - col
            dividedDifferences(row, col) = (dividedDifferences(row + 1, col - 1) - dividedDifferences(row, col - 1))/(xData(row + col - 1) - xData(row));
        end 
    end

    B = dividedDifferences(1, :);

end