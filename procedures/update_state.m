function A_new = update_state(l, A, tau, neighbors, sigma, params)
%UPDATE_STEP Update the current estimate of the state, i.e. current topology.

    n = params.n;
    A_new = zeros(n);

    for i=1:n
        for j=1:n
            if i ~= j
                A_new(i, j) = A(i, j, l) + sigma * calc_update_diff(A, tau, neighbors, l, i, j, n);
                A_new(i, j) = max(A_new(i, j), 0);
            end
        end
    end

end


function delta = calc_update_diff(A, tau, neighbors, l, i, j, n)
    delta = 0;
    % This is done only once per i-j pair like this, not for every p.
    if l == i
        d = tau(j) - A(i, j, l);
    elseif l == j
        d = tau(i) - A(i, j, l);
    else
        d = 0;
    end

    for p=1:n
        if neighbors(p)
            delta = delta + A(l, p, l) * (A(i, j, p) - A(i, j, l));
        end
    end

    delta = delta + d;
end