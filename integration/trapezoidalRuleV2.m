function [I, intervals] = trapezoidalRuleV2(f, xStart, xEnd, delta)

    % equally spaced data points seperating each interval
    xData = xStart : delta : xEnd;
    yData = f(xData);

    % number of intervals
    intervals = length(xData) - 1;

    % calculate integral
    I = delta/2 * (2*sum(yData) - yData(1) - yData(end));

end