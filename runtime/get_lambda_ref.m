function r = get_lambda_ref(k, refs)
%GET_LAMBDA_REF Get the next connectivity reference.

    persistent index
    persistent ref
    if isempty(index)
        index = 1;
    end

    if ~isempty(refs)
        if k == refs(1, index)
            ref = refs(2, index);
            index = min(index + 1, size(refs, 2));
        end
    else
        ref = -1;
    end

    r = ref;

end

