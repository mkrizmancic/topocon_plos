function [K_l2] = calc_Kl2(l, A, f, K_l2_old, ready, params)
% CALC_KL2 Calculate the error threshold for lambda tracking.

    persistent buffer;

    if isempty(buffer)
        buffer = zeros(params.n, 5);
    end

    if params.adapt_Kl2 > 0
        if ready
            K_l2 = buffer(l,1);
            buffer(l,1:end-1) = buffer(l,2:end);
            buffer(l,end) = params.adapt_Kl2 * calculate(A, f);
        else
            K_l2 = K_l2_old;
        end
    else
        K_l2 = params.K_l2;
    end
end


function max_diff = calculate(A, f)
    [existing_link_1, existing_link_2] = find(A > 0.2);
    max_diff = 0;
    for i=1:length(existing_link_1)
        max_diff = max(max_diff, (f(existing_link_1(i)) - f(existing_link_2(i)))^2);
        % d = f(existing_link_1(i)) - f(existing_link_2(i))^2;
        % if d > max_diff
        %     max_diff = d;
        %     most_sensitive = [existing_link_1(i), existing_link_2(i)];
        % end
    end
end
