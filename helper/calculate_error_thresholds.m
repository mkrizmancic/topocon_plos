function [k,thresholds] = calculate_error_thresholds(l2_refs, K_l2)
%CALCULATE_ERROR_THRESHOLDS Calculate allowable error thresholds for each reference point based on K_l2.

j = 2;
k = l2_refs(1,2):l2_refs(1,end);
thresholds = zeros(2, length(k));
for i=1:length(k)
    if l2_refs(1, j) == k(i)
        l2_ref = l2_refs(2, j);
        j = j + 1;
    end
    thresholds(1,i) = l2_ref + K_l2(k(i));
    thresholds(2,i) = l2_ref - K_l2(k(i));
end

end

