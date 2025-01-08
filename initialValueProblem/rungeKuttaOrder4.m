function xPoints = rungeKuttaOrder4(f, tPoints, x0)
% assumes tPoints are equally spaced

    deltaT = tPoints(2) - tPoints(1);
    n = length(tPoints);
    
    xPoints = zeros(1, n);
    xPoints(1) = x0; 
    for i = 2:n

        tn = tPoints(i-1);
        xn = xPoints(i-1);
        k1 = f(tn, xn);
        k2 = f(tn + deltaT./2, xn + (deltaT.*k1)./2);
        k3 = f(tn + deltaT./2, xn + (deltaT.*k2)./2);
        k4 = f(tn + deltaT, xn + deltaT.*k3);
        xPoints(i) = xn + deltaT.*(k1./6 + k2./3 + k3./3 + k4./6);
        
    end
end