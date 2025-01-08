function I = rombergIntegration(f, xStart, xEnd, delta1, delta2)

    [I1, ~] = trapezoidalRuleDelta(f, xStart, xEnd, delta1);
    [I2, ~] = trapezoidalRuleDelta(f, xStart, xEnd, delta2);
    I = I2*(1 + 1/((delta1/delta2)^2 - 1)) - I1*(1/((delta1/delta2)^2 - 1));

end

function [I, intervals] = trapezoidalRuleDelta(f, xStart, xEnd, delta)

    % data points seperating each interval
    xData = xStart : delta : xEnd;
    yData = f(xData);

    % number of intervals
    intervals = length(xData) - 1;

    % calculate integral
    I = delta/2 * (2*sum(yData) - yData(1) - yData(end));
    matlabI = delta * trapz(yData);

end