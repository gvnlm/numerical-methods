function xPoints = explicitEuler(f, tPoints, x0) 
% assumes tPoints are equally spaced

    deltaT = tPoints(2) - tPoints(1);
    n = length(tPoints);
    
    xPoints = zeros(1, n);
    xPoints(1) = x0; 
    for i = 2:n

        tn = tPoints(i-1);
        xn = xPoints(i-1);
        xPoints(i) = xn + deltaT.*f(tn, xn);
        
    end

end