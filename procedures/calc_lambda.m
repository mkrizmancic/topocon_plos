function [lambda, fiedler_ok, lambda_n] = calc_lambda(A, failed)
%CALC_LAMBDA Calculate value of lambda2 from adjacency matrix.
%
%   Input:
%       * A - Adjacency matrix [n x n]
%   Output:
%       * lambda - Value of lambda2.
%       * fiedler - Eigenvector corresponding to lambda2
%       * lambda_n - Value of lambdaN

    D = diag(sum(A, 2));
    L = D - A;

    [F, V] = eig(L);

    [v, ind] = sort(diag(V));
    f = F(:, ind);

    lambda = v(2);
    fiedler = f(:, 2);
    lambda_n = v(end);

    fiedler_ok = failed';
    fiedler_ok(failed==1) = NaN;
    fiedler_ok(failed==0) = fiedler;
end

