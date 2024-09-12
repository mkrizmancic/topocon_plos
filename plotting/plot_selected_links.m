function [fig, tiles] = plot_selected_links(fig, A, is, js, options)
%PLOT_SELECTED_LINKS Plot quality values of selected links.

n = size(A, 1);
subplots = length(is);

% fig = figure('Position', [50, 50, 1400, 900]);

tiles = tiledlayout(subplots,1,'TileSpacing','Compact','Padding','Compact');

for l=1:subplots
    nexttile
    plot_one_link(A, is(l), js(l), n)
    ylabel(sprintf("$a_{%d%d}$", is(l), js(l)), options.labels{:})
    xlabel('Step (k)', options.labels{:})
    legend('Location', 'eastoutside')
    grid on
end

end


function [] = plot_one_link(A,i,j,n)

% color = zeros(5,3);
% color(1,:) = [1,0,0];
% color(2,:) = [0,1,0];
% color(3,:) = [0,0,1];
% color(4,:) = [240,150,60] ./ 255;
% color(5,:) = [1,0,1];

for l=1:n
    plot(squeeze(A(i, j, l, :)), ...
        'DisplayName', sprintf("a^{%d}_{%d%d}", l, i, j))
    hold on
end
hold off

end