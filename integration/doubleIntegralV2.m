function I = doubleIntegralV2(f, xDelta, yDelta, yStart, yEnd)

    yData = yStart : yDelta : yEnd;

    % calculate the x integral by letting y be constant
    Ix = zeros(1, length(yData));
    for i = 1 : length(yData)
        Ix(i) = simpsonsRule(f(xData, yData(i)), xDelta);
    end

    I = simpsonsRule(Ix, yDelta);
    
end

function I = simpsonsRule(yData, delta)

    I = delta/3 * (yData(1) + 4*sum(yData(2 : 2 : end-1)) + 2*sum(yData(3 : 2 : end - 2)) + yData(end));
    
end