function [failed,neighbors,A_ok] = detect_failed(k,l,A,failed,neighbors,params)
% DETECT_FAILED Detect failed agents and remove them from the network.

    for j=1:params.n
        if l ~= j && sum(A(j,:)) < 0.5 && ~failed(j)
            failed(j) = 1;
            neighbors(j) = 0;
            fprintf("Agent %d: Failed agent %d in step %d\n",l, j, k)
        elseif l ~= j && sum(A(j,:)) > 0.5 && failed(j)
            failed(j) = 0;
            fprintf("Agent %d: Agent %d rejoined in step %d\n",l, j, k)
        end
    end

    A_ok = A(~failed, ~failed);

end

