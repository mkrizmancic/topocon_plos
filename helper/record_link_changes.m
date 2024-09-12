function [new] = record_link_changes(k, neighbors, A0, option, already_recorded)
%RECORD_LINK_CHANGES Detect which links were changed and print them.

changes = neighbors(:,:,k) - neighbors(:,:,k-1);
[r,c] = ind2sub(size(A0), find(changes));
%     new = unique(sort([r,c], 2), 'rows');
new = [r,c];
if exist('already_recorded', 'var') && ~isempty(already_recorded)
    new = new(~ismember(new,already_recorded,"rows"), :);  % Remove pairs that were already recorded earlier.
end
for row=1:size(new, 1)
    if changes(new(row, 1), new(row, 2)) > 0
        fprintf("%s  Added link %d - %d in step %d\n", option, new(row, 1), new(row, 2), k)
    else
        fprintf("%s Removed link %d - %d in step %d\n", option, new(row, 1), new(row, 2), k)
    end
end

end