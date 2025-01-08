function df = BDA1stOrder3Point(f, xPoints)
% uses 1st order (error) 3-point backward difference scheme to 
% differentiate f over the generally spaced grid xPoints

    df = zeros(1, length(xPoints));
    
    % since we have no information on xi-1 at x1, we use 1st order 3-point 
    % FDA to calculate df(x1)
    xi = xPoints(1);
    xip1 = xPoints(2);
    df(1) = (f(xip1) - f(xi))./(xip1 - xi);

    % BDA
    xi = xPoints(2:end);
    xim1 = xPoints(1:end-1);
    df(2:end) = (f(xi) - f(xim1))./(xi - xim1);

end