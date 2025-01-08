function I = kPointGaussLegendre(g, xStart, xEnd, k)
    
    [w, x] = gaussLegendrePoints(k);

    t = @(x) x*(xEnd-xStart)/2 + (xEnd+xStart)/2;
    f = @(x) g(t(x)) * ((xEnd-xStart)/2);
    
    I = 0;
    for i = 1 : length(w)
        I = I + w(i)*f(x(i));
    end

end

function [w,x] = gaussLegendrePoints(k)
    
    syms t
    x = double(vpasolve(legendreP(k,t) == 0));
    w = 2*(1-x.^2)./((k+1)^2.*(legendreP(k+1,x).^2));
    
end