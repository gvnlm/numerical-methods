function I = doubleIntegral(f, xDelta, yDelta, integralSolver)

    [~, cols] = size(f);
   
    % % calculate the x intergral when Y IS CONSTANT
    % xIntegral = zeros(1, rows);
    % for i = 1 : rows
    %     xIntegral(i) = integralSolverFunction(f(i, :), xDelta);
    % end

    % calculate the y intergral when X IS CONSTANT
    yIntegral = zeros(1, cols);
    for i = 1 : cols
        yIntegral(i) = integralSolver(f(:, i), yDelta);
    end
    
    I = integralSolver(yIntegral, xDelta);

end