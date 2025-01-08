function xPoints = crankNicolsonLinear(xnp1Formula, tPoints, x0) 
% assumes tPoints are equally spaced

% xnp1Formula(tn, xn, deltaT) - crank-nicolson formula rearranged
% such that xnp1 is a function of tn, xn and delta t

% function is the same as implicitEulerMethodLinear 

    deltaT = tPoints(2) - tPoints(1);
    n = length(tPoints);
    
    xPoints = zeros(1, n);
    xPoints(1) = x0; 
    for i = 2:n

        tn = tPoints(i-1);
        xn = xPoints(i-1);
        xPoints(i) = xnp1Formula(tn, xn, deltaT);
        
    end

end