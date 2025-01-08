function I = trapezoidalRule(yData, delta)

    I = delta/2 * (2*sum(yData) - yData(1) - yData(end));
    
end