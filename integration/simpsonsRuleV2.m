function [I, intervals] = simpsonsRuleV2(f, xStart, xEnd, delta)

    % equally spaced data points seperating each interval
    xData = xStart : delta : xEnd;
    yData = f(xData);

    % number of intervals
    intervals = length(xData) - 1;

    % check that the number of intervals is even
    if mod(intervals, 2) ~= 0
        fprintf('error: number of intervals must be even\n');
        return;
    end

    % calculate integral
    I = delta/3 * (yData(1) + 4*sum(yData(2 : 2 : end-1)) + 2*sum(yData(3 : 2 : end - 2)) + yData(end));

end