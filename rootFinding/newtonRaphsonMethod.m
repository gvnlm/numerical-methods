function root = newtonRaphsonMethod(f, initialGuess, tolerance, maxIterations)
    
    % get f'(x) as an anonymous function
    df = matlabFunction(diff(sym(f)));
    
    % initialise iteration using initial guess
    xi = initialGuess;
    fxi = f(xi);
    dfxi = df(xi);
    
    for i = 1 : maxIterations
        
        % display
        fprintf('x%d = %.6f, f(x%d) = %.6f\n', i, xi, i, fxi);

        % check if xi = NaN
        if isnan(xi)
            break;
        end

        % check if f(xi) ≈ 0
        if abs(fxi) <= tolerance
            fprintf('a root was found at %.6f after %d iterations\n', ...
                xi, i);
            root = xi;
            return;
        end

        % update
        xi = xi - fxi/dfxi;
        fxi = f(xi);
        dfxi = df(xi);

    end

    fprintf('a root could not be found within %d iterations\n', maxIterations);
    root = NaN;
    
end