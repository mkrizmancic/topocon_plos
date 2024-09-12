function A = extract_init_topo(A0)
%EXTRACT_INIT_TOPO Extract local views on the initial topology for each agent.
    n = size(A0, 1);
    A = zeros(n, n, n);

    for i=1:n
        A(i, :, i) = A0(i, :);
        A(:, i, i) = A0(:, i);
    end
end