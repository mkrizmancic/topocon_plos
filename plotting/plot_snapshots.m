function plot_snapshots(lambda, A, select, moments, graph_layout)
% PLOT_SNAPSHOTS Plot snapshots of the graph at given moments.

n = size(A, 1);

if length(moments) == 1
    num_moments = moments;
    moments = round(linspace(1, size(A, 4), moments));
else
    num_moments = length(moments);
end

if isa(graph_layout, 'matlab.graphics.chart.primitive.GraphPlot')
    vertices = [graph_layout.XData', graph_layout.YData'];
elseif strcmp(graph_layout, 'auto')
    vertices = nsidedpoly(n).Vertices;
else
    error('Unknown graph layout');
end

[p, ~] = numSubplots(num_moments);  % https://www.mathworks.com/matlabcentral/fileexchange/26310-numsubplots-neatly-arrange-subplots

tiles = tiledlayout(min(p), max(p), 'TileSpacing','Compact','Padding','Compact');

for m=moments
    if select == 0
        A_ = squeeze(mean(A(:,:,:,m), 3));
        l = mean(lambda(:,m));
    else
        A_ = A(:,:,select,m);
        l = lambda(select,m);
    end

    A_ = round((A_ + A_')/2, 1);

    nexttile;
    G = graph(A_);
    minWidth = 2; maxWidth = 4;
    LWidths = minWidth + (maxWidth-minWidth)*G.Edges.Weight/max(G.Edges.Weight);
    p = plot(G,'XData',vertices(:,1),'YData',vertices(:,2),'MarkerSize',18,'LineWidth',LWidths,'EdgeColor','k');
    text(p.XData, p.YData, p.NodeLabel, ...
        'VerticalAlignment','middle',...
        'HorizontalAlignment', 'center',...
        'Color','k','FontWeight','bold',...
        'FontSize', 10)
    p.NodeLabel = {};
    set(gca,'xtick',[])
    set(gca,'ytick',[])

    title(sprintf("\nk = %d,  \\lambda_2 = %.4f", m, l), 'FontSize', 14);
    % annotations = [annotations; m, lambda];
end

title(tiles, '.');

end