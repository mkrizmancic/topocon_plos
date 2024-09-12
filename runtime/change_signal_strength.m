function [signal_strength] = change_signal_strength(k, signal_strength, signal_refs)
%CHANGE_SIGNAL_STRENGTH Simulate a change in signal strength of selected links.
    persistent index
    if isempty(index)
        index = 1;
    end

    if ~isempty(signal_refs)
        if k == signal_refs{index, 1}
            for c=1:size(signal_refs{index, 2}, 1)
                instr = signal_refs{index, 2}(c,:);
                signal_strength(instr{1}, instr{2}) = instr{3};
            end
            index = min(index + 1, size(signal_refs, 1));
        end
    end

end

