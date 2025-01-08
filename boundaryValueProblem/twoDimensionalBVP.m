function phi = twoDimensionalBVP(phi, Q, delta, tolerance)
% Q can be either a matrix or function of x, y
% remember to divide Q by k if k ~= 1

% x and y spacing must be equal
        
    [rows, cols] = size(phi);

    % whilst any approximation of ϕ(j, i) has an error > tolerance
    maxError = Inf;
    while maxError > tolerance

        % approximate ϕ(j, i)
        for j = 2 : rows-1
            for i = 2 : cols-1
                phi(j, i) = (phi(j, i-1) + phi(j, i+1) + phi(j-1, i) + ...
                    phi(j+1, i) - delta^2*Q(j, i))/4;
            end
        end

        % find the max error of any approximation
        maxError = 0;
        for j = 2 : rows-1
            for i = 2 : cols-1
                
                % error = current approx - next approx
                error = abs(phi(j, i) - (phi(j, i-1) + phi(j, i+1) + ...
                    phi(j-1, i) + phi(j+1, i) - delta^2*Q(j, i))/4);
                if error > maxError
                    maxError = error;
                end

            end
        end

    end

end