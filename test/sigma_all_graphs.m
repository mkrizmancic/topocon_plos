% Calculate stable and optimal values of sigma and convergence rate for each
% graph example. Graphs are loaded from the standard list in `get_topology.m`.
% Save the results in a table for easier reading.

graphs = {
    'dogtail', 0;
    'karlo5', 0;
    'goodit_rpi', 0;
    'goodit_15_ok', 0;
    'goodit_15_bad', 0;
    'old_6', 0;
    'star', 5;
    'star', 10;
    'house', 0;
    'line', 5;
    'line', 10;
    'circle', 5;
    'circle', 10;
    'bridge_2', 0;
};

T = table();

for i=1:size(graphs, 1)
    [A0, ~, ~] = get_topology(graphs{i, 1}, graphs{i, 2});    % Initial topology
    n = size(A0, 1);

    D = diag(sum(A0, 2));
    L = D - A0;

    [~, V] = eig(L);
    [v, ~] = sort(diag(V));

    standard.sigma_o = 2 / (v(end) + v(2)) * 1.0;
    standard.sigma_s = 2 / (v(end));
    standard.ro = 1 - standard.sigma_o * v(2);

    extended.sigma_o = 2 / (v(end) + 2);
    extended.sigma_s = 2 / (v(end) + 1);
    extended.ro_o = 1 - extended.sigma_o;

    T = [T; table(struct2table(standard), struct2table(extended))];
end

T = round(T, 4);

for i=1:size(graphs, 1)
    if graphs{i, 2} ~= 0
        T.Properties.RowNames{i} = [graphs{i, 1}, '_', num2str(graphs{i, 2})];
    else
        T.Properties.RowNames{i} = graphs{i, 1};
    end
end

T.Properties.VariableNames = {'standard', 'extended'};