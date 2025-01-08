close all;
clear;
clc;

initialGuess = [-10; 10];
tolerance = 1e-6;

% SOLVING USING FSOLVE() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fsolve() without jacobian option
xMatlabSol1 = fsolve(@nonLinearSystem, initialGuess)

% fsolve() with jacobian option
options = optimoptions(@fsolve, 'specifyObjectiveGradient', true);
xMatlabSol2 = fsolve(@nonLinearSystemWithJacobian, initialGuess, options)

% SOLVING USING NEWTON-RAPHSON METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [J]{Xi - Xim1} = {-f}
% Xi = new guess
% Xim1 = old guess

% initialise iteration using inital guess
Xim1 = initialGuess;
[f, J] = nonLinearSystemWithJacobian(Xim1);

while max(abs(f)) > tolerance

    % calculate the new guess using a linear system solver
    Xi = LUDecompositionSolver(J, -f) + Xim1;
    Xim1 = Xi;
    
    % update f and J using the new guesses
    [f, J] = nonLinearSystemWithJacobian(Xi);

end

xNewtonRaphsonSol = Xi

% INPUT EQUATIONS AND PARTIAL DERIVATIVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = nonLinearSystem(x)
    f(1) = x(1)*exp(x(2)) - x(1)^5 + x(2) - 3;
    f(2) = x(1) + x(2) + tan(x(1)) - sin(x(2));
end

function [f, J] = nonLinearSystemWithJacobian(x)
    f(1) = x(1)*exp(x(2)) - x(1)^5 + x(2) - 3;
    f(2) = x(1) + x(2) + tan(x(1)) - sin(x(2));

    % partial derivatives of f1
    J(1, 1) = exp(x(2)) - 5*x(1)^4;
    J(1, 2) = x(1)*exp(x(2)) + 1;

    % partial derivatives of f2
    J(2, 1) = 1 + sec(x(1))*sec(x(1));
    J(2, 2) = 1  - cos(x(2));
end

% LINEAR SYSTEMS SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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