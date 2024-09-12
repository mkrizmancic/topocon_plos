%   This script will compare the time it takes to find the most sensitive
% link from the Fiedler vector using two approaches.
%
%   As a reminder, the most sensitive link is the one for which the
% absolute difference of the corresponding vector elements is the biggest.
%   The first method considered uses 'find' to get row and column indices
% of the matrix where a link exists. Each row and col pair is a link to
% consider. Then it iterates over all such links and finds the maximum.
% The downside of this approach is that each link is considered twice
% since the matrix is symmetrical.
%   The second approach uses 'nchoosek' to get all possible link
% combinations without repetition and then iterates over them to find the
% maximum.
%
% Result: 'find' is faster by almost double.
% Note: I also tested different versions of max, and this one is the
% quickest.

tf = @() test_find(A_ok, fiedler, params, steps);
tc = @() test_comb(A_ok, fiedler, params, steps);

time_find = timeit(tf)
time_comb = timeit(tc)


function test_find(A_ok, fiedler, params, steps)
    l=1;
    for k=2:steps
        test = calc_find(A_ok{l}, fiedler(:,l,k), params);
    end
end

function [K_l2] = calc_find(A, f, params)
% Calculate the error threshold for lambda tracking.
    
    if params.adapt_Kl2 > 0
        [existing_link_1, existing_link_2] = find(A > 0.2);
        max_diff = 0;
        for i=1:length(existing_link_1)
            max_diff = max(max_diff, (f(existing_link_1(i)) - f(existing_link_2(i)))^2);
        end
        K_l2 = params.adapt_Kl2 * max_diff;
    else
        K_l2 = params.K_l2;
    end
end

function test_comb(A_ok, fiedler, params, steps)
    l=1;
    for k=2:steps
        test = calc_comb(A_ok{l}, fiedler(:,l,k), params);
    end
end

function [K_l2] = calc_comb(A, f, params)
% Calculate the error threshold for lambda tracking.
    
    if params.adapt_Kl2 > 0
        combs = nchoosek(1:params.n, 2);
        max_diff = 0;
        for k=1:size(combs, 1)
            i = combs(k, 1);
            j = combs(k, 2);
            if A(i, j) > 0.2
                max_diff = max(max_diff, (f(i) - f(j))^2);
            end
        end
        K_l2 = params.adapt_Kl2 * max_diff;
    else
        K_l2 = params.K_l2;
    end
end