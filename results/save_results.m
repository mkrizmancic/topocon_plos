z_baseFilename = get_filename();

if ~isempty(z_baseFilename)
    z_filename = strcat(z_baseFilename, '.mat');
    save(z_filename{1}, 'A0', 'l2_refs', 'signal_refs', 'params', 'estimation_period', 'steps', 'rand_stream')

    z_tempFields = fieldnames(z_figures);
    for i = 1:numel(z_tempFields)
        z_tempFigHandle = z_figures.(z_tempFields{i});
        z_filename = strcat(z_baseFilename, '_', z_tempFields{i});

        % Save figure as SVG
        saveas(z_tempFigHandle, z_filename{1}, 'svg');
    end
end

% Put stuff in function to avoid messing up the workspace.
function [base_filename] = get_filename()
    name = inputdlg({'Enter file name:'},'Input',[1 45],{'matlab'});
    script_path = mfilename('fullpath');
    [script_folder, ~, ~] = fileparts(script_path);
    base_filename = fullfile(script_folder, name);
end