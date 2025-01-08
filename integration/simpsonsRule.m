function I = simpsonsRule(yData, delta)

    I = delta/3 * (yData(1) + 4*sum(yData(2 : 2 : end-1)) + 2*sum(yData(3 : 2 : end - 2)) + yData(end));
    
end