function xPoints = rungeKuttaOrder2System2(f0, f1, tPoints, initialX0, initialX1) 
% assumes tPoints are equally spaced

% dx0dt = f0(t, x0, x1) 
% dx1dt = f1(t, x0, x1) 

    deltaT = tPoints(2) - tPoints(1);
    n = length(tPoints);
    
    xPoints = zeros(n, 2);
    xPoints(1, 1) = initialX0;
    xPoints(1, 2) = initialX1;
    for i = 2:n

        tn = tPoints(i-1);
        x0n = xPoints(i-1, 1);
        x1n = xPoints(i-1, 2);

        % k1 associated with x0 and x1 respectively
        k10 = f0(tn, x0n, x1n);
        k11 = f1(tn, x0n, x1n);
        

        % k2 associated with x0 and x1 respectively
        k20 = f0(tn + deltaT, x0n + deltaT.*k10, x1n + deltaT.*k11);
        k21 = f1(tn + deltaT, x0n + deltaT.*k10, x1n + deltaT.*k11);

        xPoints(i, 1) = x0n + deltaT.*(k10./2 + k20./2);
        xPoints(i, 2) = x1n + deltaT.*(k11./2 + k21./2);

    end

end