function D = CDA2ndOrder3PointDfMatrix(n, delta)
% 2nd order (error) CDA derivative matrix for n equally points

    D = zeros(n);
    D(1, 1) = -2;
    D(1, 2) = 2;
    D(n, n-1) = -2;
    D(n, n) = 2;
    for row = 2 : n-1
        D(row, row-1) = -1;
        D(row, row+1) = 1;
    end
    D = D./(2*delta);

end