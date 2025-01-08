function xPoints = implicitEulerLinearSystem2(K, tPoints, initialX0, initialX1) 
% assumes tPoints are equally spaced

    deltaT = tPoints(2) - tPoints(1);
    n = length(tPoints);
    
    xPoints = zeros(n, 2);
    xPoints(1, 1) = initialX0;
    xPoints(1, 2) = initialX1;
    for i = 2:n

        xPoints(i, :) = inv(eye(2) - K.*deltaT) * xPoints(i-1, :)';
        
    end

end