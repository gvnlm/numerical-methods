function xPoints = taylorOrder2(f, dfdt, dfdx, tPoints, x0) 
% assumes tPoints are equally spaced
% dfdt = partial derivative of f with respect to t
% dfdxConstant = the constant partial derivative of f with respect to x

% local error: O(deltaT^3)
% global error: O(deltaT^2)

    deltaT = tPoints(2) - tPoints(1);
    n = length(tPoints);
    
    xPoints = zeros(1, n);
    xPoints(1) = x0; 
    for i = 2:n

        tn = tPoints(i-1);
        xn = xPoints(i-1);
        xPoints(i) = xn + deltaT.*f(tn, xn) + deltaT.^2./2 .* ...
            (dfdt(tn, xn) + dfdx(tn, xn).*f(tn, xn));
        
    end

end