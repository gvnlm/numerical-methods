function df = FDA1stOrder3Point(f, xPoints)
% uses 1st order (error) 3-point forward difference scheme to differentiate 
% f over the generally spaced  grid xPoints

    df = zeros(1, length(xPoints));
    
    % since we have no information on xi+1 at xn, we use 1st order 3-point 
    % BDA to calculate df(xn)
    xi = xPoints(end);
    xim1 = xPoints(end-1);
    df(end) = (f(xi) - f(xim1))./(xi - xim1);

    % FDA
    xi = xPoints(1:end-1);
    xip1 = xPoints(2:end);
    df(1:end-1 ) = (f(xip1) - f(xi))./(xip1 - xi);

end