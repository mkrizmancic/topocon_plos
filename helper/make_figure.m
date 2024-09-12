function [enabled, figures] = make_figure(flags, figures, name)
%MAKE_FIGURE Create a figure if the corresponding flag is set.

enabled = false;
if isfield(flags, name) && flags.(name)
    enabled = true;
    if isfield(figures, name) && isgraphics(figures.(name))
        figure(figures.(name));
    else
        figures.(name) = figure();
    end
end

end