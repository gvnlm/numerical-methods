function root = bisectionMethod(f, xLower, xUpper, tolerance, maxIterations)

    % initialise iteration using given bounds
    xMiddle = (xUpper + xLower)/2;
    fxLower = f(xLower);
    fxMiddle = f(xMiddle);
    fxUpper = f(xUpper);

    for i = 1 : maxIterations

        % display 
        fprintf('xLower = %.6f, xMiddle = %.6f, xUpper = %.6f\n', ...
            xLower, xMiddle, xUpper);
        fprintf('f(xLower) = %.6f, f(xMiddle) = %.6f, f(xUpper) = %.6f\n\n', ...
            fxLower, fxMiddle, fxUpper);
        
        % check if f(xMiddle) ≈ 0
        if abs(fxMiddle) <= tolerance
            fprintf('a root was found at %.6f after %d iterations\n', ...
                xMiddle, i);
            root = xMiddle;
            return;
        end
    
        % update lower/upper bound
        if fxLower * fxMiddle < 0
            xUpper = xMiddle;
            fxUpper = f(xUpper);
        elseif fxUpper * fxMiddle < 0 
            xLower = xMiddle;
            fxLower = f(xLower);
        else
            fprintf('error: xLower, xMiddle and xUpper share the same sign\n');
            root = NaN;
            return;
        end
        
        % calculate new xMiddle
        xMiddle = (xUpper + xLower)/2;
        fxMiddle = f(xMiddle);

    end

    fprintf('a root could not be found within %d iterations\n', maxIterations);
    root = NaN;

end