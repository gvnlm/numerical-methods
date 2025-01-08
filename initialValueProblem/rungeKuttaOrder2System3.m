function xPoints = rungeKuttaOrder2System3(f0, f1, f2, tPoints, initialX0, initialX1, initialX2) 
% assumes tPoints are equally spaced

% dx0dt = f0(t, x0, x1, x2) 
% dx1dt = f1(t, x0, x1, x2) 
% dx2dt = f2(t, x0, x1, x2) 

    deltaT = tPoints(2) - tPoints(1);
    n = length(tPoints);
    
    xPoints = zeros(n, 3);
    xPoints(1, 1) = initialX0;
    xPoints(1, 2) = initialX1;
    xPoints(1, 3) = initialX2;
    for i = 2:n

        tn = tPoints(i-1);
        x0n = xPoints(i-1, 1);
        x1n = xPoints(i-1, 2);
        x2n = xPoints(i-1, 3);

        % k1 associated with x0, x1 and x2 respectively
        k10 = f0(tn, x0n, x1n, x2n);
        k11 = f1(tn, x0n, x1n, x2n);
        k12 = f2(tn, x0n, x1n, x2n);
        

        % k2 associated with x0, x1 and x2 respectively
        k20 = f0(tn + deltaT, x0n + deltaT.*k10, x1n + deltaT.*k11, x2n + deltaT.*k12);
        k21 = f1(tn + deltaT, x0n + deltaT.*k10, x1n + deltaT.*k11, x2n + deltaT.*k12);
        k22 = f2(tn + deltaT, x0n + deltaT.*k10, x1n + deltaT.*k11, x2n + deltaT.*k12);

        xPoints(i, 1) = x0n + deltaT.*(k10./2 + k20./2);
        xPoints(i, 2) = x1n + deltaT.*(k11./2 + k21./2);
        xPoints(i, 3) = x2n + deltaT.*(k12./2 + k22./2);

    end

end