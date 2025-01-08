function L = choleskyDecomposition(A)

[n,n]=size(A);
L=zeros(n,n);

% column by column
for j = 1:n
    
    % start with diagonal value L(j, j)
    sum = 0;
    for k = 1 : j - 1
        sum = sum + L(j, k)^2;
    end
    L(j, j) = sqrt(A(j, j) - sum);

    % then do values directly below L(i, j);
    for i = j + 1 : n
        sum = 0;
        for k = 1 : j - 1
            sum = sum + L(i, k)*L(j, k);
        end
        L(i, j) = (A(i, j) - sum)/L(j, j);
    end
    
    % move on to next column
end

end     