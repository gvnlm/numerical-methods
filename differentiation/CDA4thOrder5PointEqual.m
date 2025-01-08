function df = CDA4thOrder5PointEqual(f, xPoints)
% uses 4th order (error) 5-point centred difference scheme to differentiate 
% f over the equally spaced grid xPoints

    delta = xPoints(2) - xPoints(1);

    df = zeros(1, length(xPoints));

    % since we have no information on xi-1 and xi-2 at x1, we use 2nd order 
    % 5-point FDA to calculate df(x1)
    xi = xPoints(1);
    xip1 = xPoints(2);
    xip2 = xPoints(3);
    df(1) = (-3*f(xi) + 4*f(xip1) - f(xip2))/(2*delta);

    % since we have no information on xi-2 at x2, we use 3rd order 
    % 5-point partial FDA to calculate df(x2)
    xim1 = xPoints(1);
    xi = xPoints(2);
    xip1 = xPoints(3);
    xip2 = xPoints(4);
    df(2) = (-2*f(xim1) - 3*f(xi) + 6*f(xip1) - f(xip2))/(6*delta);
    
    % 5-point CDA
    xim2 = xPoints(1:end-4);
    xim1 = xPoints(2:end-3);
    xip1 = xPoints(4:end-1);
    xip2 = xPoints(5:end);
    df(3:end-2) = (f(xim2) - 8*f(xim1) + 8*f(xip1) - f(xip2))/(12*delta); 

    % since we have no information on xi+1 and xi+2 at xn, we use 2nd order 
    % 5-point BDA to calculate df(xn)
    xi = xPoints(end);
    xim1 = xPoints(end-1);
    xim2 = xPoints(end-2);
    df(end) = (3*f(xi) - 4*f(xim1) + f(xim2))/(2*delta);

    % since we have no information on xi+2 at xn-1, we use 3rd order
    % 5-point partial BDA to calculate df(xn-1)
    xip1 = xPoints(end);
    xi = xPoints(end-1);
    xim1 = xPoints(end-2);
    xim2 = xPoints(end-3);
    df(end-1) = (f(xim2) - 6*f(xim1) + 3*f(xi) + 2*f(xip1))/(6*delta);

end