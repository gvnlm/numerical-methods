function I = twoPointGaussLegendre(g, xStart, xEnd)

    % optimal xi
    x0 = -1/sqrt(3);
    x1 = 1/sqrt(3);

    t = @(x) x*(xEnd-xStart)/2 + (xEnd+xStart)/2;
    f = @(x) g(t(x)) * ((xEnd-xStart)/2);
    
    I = f(x0) + f(x1);  % note that the optimal weights wi are all equal to 1

end