% Make a plot for comparing connectivity with and without sigma adaptation.

% 1. Run the main script with sigma adaptation set to 'false'.
% 2. Rename the resulting variable 'lambda' to 'lambda_old'.
% 3. Run the main script with sigma adaptation set to 'true'.
% 4. Run this script to plot the results.

fig = figure(99);
do_plot.subSigma = true;

% Set up axis.
if isfield(do_plot, 'subSigma') && do_plot.subSigma
    fig.Position(4) = fig.Position(4) * 1.5;
    tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
    ax = nexttile;
else
    ax = gca;
end

% Plot old values of lambda.
h = plot(lambda_old', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Reduce opacity of old lines.
for i=1:numel(h)
    c = get(h(i), 'Color');
    set(h(i), 'Color', [c 0.3]);
end

% Reset colors.
hold on
set(gca,'ColorOrderIndex',1)

% Plot new values of lambda.
l2_refs = [l2_refs, [steps; l2_refs(2,end)]];
h = plot(lambda', 'LineWidth', 1.5);
set(h, {'DisplayName'}, z_dispNames')
stairs(l2_refs(1,2:end), l2_refs(2,2:end), 'k--', 'DisplayName','ref')
hold off
grid on
xlim([0 steps])
% ylim([0 3.1])
xlabel('Step (k)', z_labelOptions{:})
ylabel('Algebraic connectivity $\lambda_2$', z_labelOptions{:})
% title('\lambda_2 through steps')
legend('Location','best','NumColumns',floor(n/5),'FontSize',10)

% Plot values of sigma.
if isfield(do_plot, 'subSigma') && do_plot.subSigma
    nexttile
    h = plot(sigma', 'LineWidth', 1.5);
    set(h, {'DisplayName'}, z_dispNames')
    grid on
    xlabel('Step (k)', z_labelOptions{:})
    ylabel('Update step size $\sigma$', z_labelOptions{:})
    % title('\\sigma through steps')
    legend('Location','best','NumColumns',floor(n/5),'FontSize',10)
    % ylim([0.05, 0.45])
end