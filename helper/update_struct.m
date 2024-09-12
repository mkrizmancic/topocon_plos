function [old_struct] = update_struct(old_struct,new_struct)
%UPDATE_STRUCT Update old_struct with values from new_struct for common fields.

old_fields = fieldnames(old_struct);

for i = 1:numel(old_fields)
    % Check if struct2 also contains the field
    if isfield(new_struct, old_fields{i})
        % Update the value in struct1 with the value from struct2
        old_struct.(old_fields{i}) = new_struct.(old_fields{i});
    end
end
end