function df = CDA2ndOrder3Point(f, xPoints)
% uses 2nd order (error) 3-point centred difference scheme to differentiate 
% f over the generally spaced grid xPoints

    df = zeros(1, length(xPoints));

    % since we have no information on xi-1 at x1, we use 1st order 3-point 
    % FDA to calculate df(x1)
    xi = xPoints(1);
    xip1 = xPoints(2);
    df(1) = (f(xip1) - f(xi))./(xip1 - xi);
    
    % CDA
    xim1 = xPoints(1:end-2);
    xip1 = xPoints(3:end);
    df(2:end-1)= (f(xip1)-f(xim1))./(xip1-xim1);

    % since we have no information on xi+1 at xn, we use 1st oder 3-point 
    % BDA to calculate df(xn)
    xi = xPoints(end);
    xim1 = xPoints(end-1);
    df(end) = (f(xi) - f(xim1))./(xi - xim1);

end