function X = leastSquaresPolynomial(xData, yData, n)

    % check that xData and yData are the same length
    if length(xData) ~= length(yData)
        fprintf('error: length(xData) ~= length(yData)\n');
        X = NaN;
        return;
    end

    % store sum(xData^i) to prevent repeating the same calculations
    xSums = zeros(1, 2*n + 1);
    for i = 0 : 2*n
        xSums(i + 1) = sum(xData.^i);
    end

    % [A]
    A = zeros(n + 1);
    for row = 1 : n+1
        for col = 1:n + 1
            A(row, col) = xSums(row + col - 1);
        end
    end

    % {C} 
    C = zeros(n + 1, 1);
    for row = 1 : n+1
        C(row, 1) = sum(xData.^(row - 1).*yData);
    end

    % [A]{X} = {C}
    X = LUDecompositionSolver(A, C);

end

function X = LUDecompositionSolver(A, C)

    % [L][U]{X} = {C}
    [L, U] = croutsAlgorithm(A);

    % [L]{R} = {C}
    R = lowerTriangleSolver(L, C);

    % [U]{X} = {R}
    X = upperTriangleSolver(U, R);
    
end

function [L,U]=croutsAlgorithm(A)

[n,n]=size(A);
L=zeros(n,n);
U=zeros(n,n);

for i=1:n
    L(i,1)=A(i,1);
end
for j=2:n
    U(1,j)=A(1,j)/L(1,1);
end

for j=2:n
    for i=j:n
        sum=0;
        for k=1:j-1
            sum=sum+L(i,k)*U(k,j);
        end
        L(i,j)=A(i,j)-sum;
    end
        
    for k=j+1:n
        sum=0;
        for i=1:j-1
            sum=sum+L(j,i)*U(i,k);
        end
        U(j,k)=(A(j,k)-sum)/L(j,j);
    end
end

for i=1:n
   U(i,i)=1.0; 
end
    
end    

function X = lowerTriangleSolver(A, C)
    
    % check that A is a square matrix
    [rows, cols] = size(A);
    if rows ~= cols
        fprintf('error: A must be a square matrix\n');
        X = NaN;
        return;
    else
        n = rows;
    end

    % initialise X
    X = zeros(n, 1);

    % calculate x1
    X(1) = C(1)/A(1, 1);

    % calculate x2 to xn
    for i = 2 : n

        sum = 0;
        for j = 1 : i - 1
            sum = sum + A(i, j)*X(j);
        end
        
        % ax + sum = c  =>  x = (c - sum) / a 
        X(i) = (C(i) - sum)/A(i, i);  % error if A(i, i) = 0

    end

end

function X = upperTriangleSolver(A, C)
    
    % check that A is a square matrix
    [rows, cols] = size(A);
    if rows ~= cols
        fprintf('error: A must be a square matrix\n');
        X = NaN;
        return;
    else
        n = rows;
    end

    % initialise X
    X = zeros(n, 1);

    % calculate xn
    X(n) = C(n)/A(n, n);

    % calculate xn-1 to x1
    for i = n-1 : -1 : 1
        
        sum = 0;
        for j = i + 1 : n
            sum = sum + A(i, j)*X(j);
        end
    
        % ax + sum = c  =>  x = (c - sum)/a 
        X(i) = (C(i) - sum)/A(i, i);  % error if A(i, i) = 0

    end

end