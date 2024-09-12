% This function calculates the maximum value of discretization step size
% according to (Kia, 2019, S17). The input is the system matrix A (not the
% adjacency matrix A).
%
% For a lot of systems, this procedure does not return good values (10e15)

function [sigma] = kia_s17_stable(A)
    [~, V] = eig(A);
    v = diag(V);

    sigma = min(-2 * real(v) ./ abs(v).^2, [], 'all');
end