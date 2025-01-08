function root = fixedPointIteration(g, initialGuess, tolerance, maxIterations)

    % initialise iteration using initial guess
    xi = initialGuess;
    gxi = g(xi);

    for i = 1 : maxIterations

        % display
        fprintf('x%d = %.6f, g(x%d) = %.6f, difference = %.6f\n', ...
            i, xi, i, gxi, abs(xi - gxi));

        % check if xi ≈ gxi
        if abs(xi - gxi) <= tolerance
            fprintf('a root was found at %.6f after %d iterations\n', ...
                xi, i);
            root = xi;
            return;
        end

        % update
        xi = gxi;
        gxi = g(xi);

    end

    fprintf('a root could not be found within %d iterations\n', maxIterations);
    root = NaN;

end