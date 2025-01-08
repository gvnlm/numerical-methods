function xPoints = rungeKuttaOrder2(f, tPoints, x0)
% assumes tPoints are equally spaced
% assumes a2 = 1/2

    deltaT = tPoints(2) - tPoints(1);
    n = length(tPoints);
    
    xPoints = zeros(1, n);
    xPoints(1) = x0; 
    for i = 2:n

        tn = tPoints(i-1);
        xn = xPoints(i-1);
        k1 = f(tn, xn);
        k2 = f(tn + deltaT, xn + deltaT.*k1);
        xPoints(i) = xn + deltaT.*(k1./2 + k2./2);
        
    end
end