function xPoints = implicitEulerNonlinear(f, dfdx, tPoints, x0, ...
    tolerance, maxIterations) 
% assumes tPoints are equally spaced

% dfdx(t, x) - partial derivative of f with respect to x

    deltaT = tPoints(2) - tPoints(1);
    n = length(tPoints);
    
    % starting from x0, compute xnp1
    xPoints = NaN .* ones(1, n);
    xPoints(1) = x0; 
    for i = 1:n-1
        
        % constants
        tnp1 = tPoints(i+1);
        xn = xPoints(i);

        % rearrange implicit euler formula such that it is a root finding
        % problem
        % => xnp1 - xn - deltaT*f(tnp1, xnp1) = 0
        g = @(xnp1) xnp1 - xn- deltaT.*f(tnp1, xnp1); % = 0
        dg = @(xnp1) 1 - deltaT.*dfdx(tnp1, xnp1);

        % use newton raphson root finding method to solve xnp1 with initial
        % guess xn
        xPoints(i+1) = newtonRaphsonMethodManualDf(g, dg, xn, ...
            tolerance, maxIterations);

        if isnan(xPoints(i))
            fprintf(['error: could not find a solution for x%d. ' ...
                'try increasing maxIterations\n'], i-1);
            return;
        end
        
    end

end

function root = newtonRaphsonMethodManualDf(f, df, initialGuess, tolerance, maxIterations)
    
    % initialise iteration using initial guess
    xi = initialGuess;
    fxi = f(xi);
    dfxi = df(xi);
    
    for i = 1 : maxIterations

        if isnan(xi)
            root = NaN;
            return;
        end

        % check if f(xi) ≈ 0
        if abs(fxi) <= tolerance
            root = xi;
            return;
        end

        % update
        xi = xi - fxi/dfxi;
        fxi = f(xi);
        dfxi = df(xi);

    end

    root = NaN;
    
end