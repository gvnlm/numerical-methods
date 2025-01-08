function f = newtonInterpolationEvaluation(xData, B, domain)

    % check that xData and B are the same length
    if length(xData) ~= length(B)
        fprintf('error: xData and B must be the same length\n');
        return;
    end
    
    f = zeros(size(domain));
    product = 1;
    for i = 1 : length(B)
        f = f + B(i).*product;
        product = product.*(domain - xData(i));
    end

end