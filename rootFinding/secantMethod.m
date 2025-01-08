function root = secantMethod(f, initialGuess0, initialGuess1, tolerance, maxIterations)
    
    % initialise iteration using initial guesses
    xCurr = initialGuess1;
    xPrev = initialGuess0;
    fxCurr = f(xCurr);
    fxPrev = f(xPrev);
    
    for i = 1 : maxIterations

        % display
        fprintf('x%d = %.6f, f(x%d) = %.6f\n', i, xCurr, i, fxCurr);

        % check if f(xCurr) ≈ 0
        if abs(fxCurr) <= tolerance
            fprintf('a root was found at %.6f after %d iterations\n', ...
                xCurr, i);
            root = xCurr;
            return;
        end

        % update
        temp = xCurr;
        xCurr = xCurr - (xCurr - xPrev)*fxCurr/(fxCurr - fxPrev);
        xPrev = temp;
        fxPrev = fxCurr;
        fxCurr = f(xCurr);
        
    end

    fprintf('a root could not be found within %d iterations\n', maxIterations);
    root = NaN;
    
end