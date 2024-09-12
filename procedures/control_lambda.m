function [neighbors, no_disconnect, no_reconnect, changed] = control_lambda( ...
    k, l, l2, l2_ref, K_l2, f, A, neighbors, no_disconnect, no_reconnect, params)
%CONTROL_LAMBDA Add or remove links to achieve desired value of lambda2.

    f = round(f, 3);
    num_active = params.n - sum(isnan(f));
    changed = [];

    % TODO: it would be better to have a parameter for multiple of tau_window

    % Adding links
    if l2 < l2_ref - K_l2
        changed = add_link(f, A, l);
        if ~isempty(changed)
            neighbors(changed) = 1;
            no_disconnect(l, changed) = k + 3*params.tau_window;
        end

    % Removing links
    elseif l2 > l2_ref + K_l2
        changed = remove_link(f, A, l);
        if ~isempty(changed)
            d_min_desired = (num_active - 1) * l2_ref / num_active;
            d_neighbor = sum(A(changed,:));
            d_me = sum(A(l,:));
            if d_neighbor - A(changed, l) > d_min_desired && d_me - A(l, changed) > d_min_desired && sum(neighbors) > 1
                neighbors(changed) = 0;
                no_reconnect(l, changed) = k + 3*params.tau_window;
            end
        end

    % Removing redundant links
    elseif params.remove_zero_links
        changed = remove_zero_link(f, A, l);

        if ~isempty(changed)
            d_min_desired = (num_active - 1) * l2_ref / num_active;
            d_neighbor = sum(A(changed,:));
            d_me = sum(A(l,:));
            if d_neighbor - A(changed, l) > d_min_desired && d_me - A(l, changed) > d_min_desired && sum(neighbors) > 1
                neighbors(changed) = 0;
                no_reconnect(l, changed) = k + 3*params.tau_window;
            end
        end
    end

end


function added = add_link(f, A, l)
    added = [];

    n = size(A, 1);

    found = false;
    combs = nchoosek(1:n, 2);
    search = zeros(size(combs, 1), 1);
    for k=1:size(combs, 1)
        i = combs(k, 1);
        j = combs(k, 2);
        if A(i, j) < 0.1
            found = true;
            search(k) = (f(i) - f(j))^2;
        end
    end

    if found
        [M, I] = max(search);
        if combs(I, 1) == l
            added = combs(I, 2);
        elseif combs(I, 2) == l
            added = combs(I, 1);
        end
    end
end

function removed = remove_link(f, A, l)
    removed = [];

    n = size(A, 1);

    found = false;
    combs = nchoosek(1:n, 2);
    search = 1000 * ones(size(combs, 1), 1);
    for k=1:size(combs, 1)
        i = combs(k, 1);
        j = combs(k, 2);
        if abs(f(i) - f(j)) > 10e-3 && A(i, j) > 0.9
            found = true;
            search(k) = (f(i) - f(j))^2;
        end
    end

    if found
        [M, I] = min(search);
        if combs(I, 1) == l
            removed = combs(I, 2);
        elseif combs(I, 2) == l
            removed = combs(I, 1);
        end
    end
end

function removed = remove_zero_link(f, A, l)
    removed = [];

    n = size(A, 1);

    found = false;
    combs = nchoosek(1:n, 2);
    search = 1000 * ones(size(combs, 1), 1);
    for k=1:size(combs, 1)
        i = combs(k, 1);
        j = combs(k, 2);
        if abs(f(i) - f(j)) < 10e-7 && A(i, j) > 0.5
            found = true;
            search(k) = (f(i) - f(j))^2;
        end
    end

    if found
        [M, I] = min(search);
        if combs(I, 1) == l
            removed = combs(I, 2);
        elseif combs(I, 2) == l
            removed = combs(I, 1);
        end
    end
end