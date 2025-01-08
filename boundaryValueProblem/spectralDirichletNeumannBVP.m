function yPoints = spectralDirichletNeumannBVP(f, xPoints, y0, ynDf)
% y'' = f(xPoints)
% dirichlet boundary condition at y0
% neumann boundary condition at yn

    % number of points
    nPoints = length(xPoints);

    % [A]{y} = {c}
        
        % [A]
        D = spectralDfMatrix(xPoints);
        D2 = D^2;
        A = D2(2:end, 2:end);
        A(end, :) = D(end, 2:end);

        % {y}
        yPoints = zeros(nPoints, 1);
        yPoints(1) = y0;
    
        % {c}
        c = zeros(nPoints-1, 1);
        c(1:end-1) = f(xPoints(2:end-1)) - y0.*D2(2:end-1, 1);
        c(end) = ynDf - D(end, 1);
       
        yPoints(2:end) = linsolve(A, c);

end

function D = spectralDfMatrix(tPoints)
% creates the 1st order derivative matrix for f(tPoints), evaluated from
% tPoints

% {f'} = [D]{f}
    
    a = tPoints(1);
    b = tPoints(end);

    % fit tPoints over the domain [-1, 1]
    xPoints = (2.*tPoints)./(b-a) - (b+a)/(b-a);
    Dx = spectralDfMatrixNeg1To1(xPoints);

    D = 2/(b-a) * Dx;
    
end

function D = spectralDfMatrixNeg1To1(xPoints)
% creates the 1st order derivative matrix for f(xPoints), evaluated from
% xPoints where xPoints = [-1, 1]

% xPoints must be from -1 to 1

% {f'} = [D]{f}
    
    % number of points
    nPoints = length(xPoints);

    D = zeros(nPoints);
    
    numerator = ones(1, nPoints);
    for i = 1:nPoints
        for k = 1:nPoints
            if k ~= i
                numerator(i) = numerator(i)*(xPoints(i)-xPoints(k));
            end
        end
    end

    denominator = ones(1, nPoints);
    for j = 1:nPoints
        for k = 1:nPoints
            if k~= j
                denominator(j) = denominator(j)*(xPoints(j)-xPoints(k));
            end
        end
    end

    % calculate all D(i, j), where i ~= j (all non-diagonal elements)
    for i = 1:nPoints
        for j = 1:nPoints
            if i ~= j
                D(i, j) = numerator(i)/((xPoints(i)-xPoints(j))*denominator(j));
            end
        end
    end

    % calculate all D(i, i) (all diagonal elements)
    for i = 1:nPoints
        D(i, i) = 0;
        for k = 1:nPoints
            if k ~= i
                D(i, i) = D(i, i) + 1/(xPoints(i) - xPoints(k));
            end
        end
    end

end