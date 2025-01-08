function root = methodOfFalsePosition(f, xLower, xUpper, tolerance, maxIterations)
    
    % initialise iteration using given bounds
    xRootGuess = xUpper - f(xUpper)*(xLower - xUpper)/(f(xLower) - f(xUpper));
    fxLower = f(xLower);
    fxRootGuess = f(xRootGuess);
    fxUpper = f(xUpper);

    for i = 1 : maxIterations

        % display
        fprintf('xLower = %.6f, xRootGuess = %.6f, xUpper = %.6f\n', ...
            xLower, xRootGuess, xUpper);
        fprintf('f(xLower) = %.6f, f(RootGuess) = %.6f, f(xUpper) = %.6f\n\n', ...
            fxLower, fxRootGuess, fxUpper);
        
        % check if f(xRootGuess) ≈ 0
        if abs(fxRootGuess) <= tolerance
            fprintf('a root was found at %.6f after %d iterations\n', ...
                xRootGuess, i);
            root = xRootGuess;
            return;
        end
    
        % update lower/upper bound
        if fxLower * fxRootGuess < 0
            xUpper = xRootGuess;
            fxUpper = f(xUpper);
        elseif fxUpper * fxRootGuess < 0 
            xLower = xRootGuess;
            fxLower = f(xLower);
        else
            fprintf('error: xLower, xRootGuess and xUpper share the same sign\n');
            root = NaN;
            return;
        end

        % calculate new xRootGuess
        xRootGuess = xUpper - f(xUpper)*(xLower - xUpper)/(f(xLower) - f(xUpper));
        fxRootGuess = f(xRootGuess);

    end

    fprintf('a root could not be found within %d iterations\n', maxIterations);
    root = NaN;
    
end