function I = threePointGaussLegendre(g, xStart, xEnd)
    
    % optimal xi
    x0 = -(sqrt(3/5));
    x1 = 0;
    x2 = sqrt(3/5);

    % optimal wi
    w0 = 5/9;
    w1 = 8/9;
    w2 = 5/9;

    t = @(x) x*(xEnd-xStart)/2 + (xEnd+xStart)/2;
    f = @(x) g(t(x)) * ((xEnd-xStart)/2);
    
    I = w0*f(x0) + w1*f(x1) + w2*f(x2);  
    
end