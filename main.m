% Here goes the intro.

% Note to self: in the code, it makes more sense to use (k, k-1) than
% (k+1, k) because of the logical explanation of the algorithm. We
% calculate the current (k) estimate of the adjacency matrix based on
% former values and calculate the current (k) value of lambda. We are not
% interested in future (k+1) values.

clearvars -except z_figures lambda_old;
clear get_lambda_ref change_signal_strength simulate_tau calc_Kl2; % Clear persistant variables from functions that use them.

%% Setup
%%% Choose this option if you want to recreate figures from the paper.
[A0, l2_refs, signal_refs, override_params, z_doPlot] = get_setup_for_fig(11);
%%% Choose this option if you want to explore and play around.
% [A0, l2_refs, signal_refs, override_params] = get_topology_new('follow', 6);

n = size(A0, 1);                    % Number of agents

%%% User settings
params.adapt_sigma  = true;         % Use adaptive sigma.
params.sigma        = 0.15;         % Fixed value of sigma if adaption is disabled.
params.adapt_Kl2    = 0.20;         % K_l2 spread. If zero, adaptation is disabled.
params.K_l2         = 0.50;         % Fixed value of K_l2 if adaptation is disabled.
params.kappa        = 2 * n;        % Connectivity control period. 'n' is enough for all agents to detect convergence.
params.norms_limit  = 0.01;         % TODO: this should be a function of the topology.
handshake.neighbors = 'active';     % Immediately update neighbors only for the agent initiating the link change ('active') or for 'both' of them.
handshake.taus      = 'none';       % Immediately update tau only for the agent initiating the link change ('active'), 'both' of them, or 'none'.
handshake.As        = 'none';       % Immediately update A matrix only for the agent initiating the link change ('active'), 'both' of them, or 'none'.
params.handshake = handshake;
params.tau_method = 'simple';       % Method of calculating tau average - 'simple', 'ema', or 'force'.
params.tau_window = 20;             % Window for averaging number of received messages.
params.n = n;
params.remove_zero_links = false;   % Enable energy conserving feature.

%%% Select random number generator.
% rand_stream = RandStream.getGlobalStream;
rand_stream = qrandstream('halton',1,'Skip',11,'Leap',1e2);

%%% Simulation params
steps = 1000;                   % Total number of steps
estimation_period = 100;        % Allow a period for initial estimation
signal_strength = ones(n);      % Set the signal strength to mimic the link quality between agents.

%%% Override params
params = update_struct(params, override_params);
if isfield(override_params, 'steps') && ~isempty(override_params.steps)
    steps = override_params.steps;
end

%%% Declare variables. Matlab is faster with preallocation.
% Essential variables for consensus
A = zeros(n, n, n, steps);
lambda = zeros(n, steps);
lambda_n = zeros(n, steps);
fiedler = zeros(n, n, steps);
tau = zeros(n, n, steps);
neighbors = zeros(n, n, steps);
failed = zeros(n, n);
A_ok = cell(1, n);
sigma = zeros(n, steps);
K_l2 = zeros(n, steps);
% Helper variables for connectivity control
no_disconnect = zeros(n, n);
no_reconnect = zeros(n, n);
% Helper variables for convergence detection
local_norms = zeros(n, steps);
ready_flags = zeros(n, n, steps);
ready_all = zeros(n, steps);
ready_delay = zeros(n, 1);
ready_mask = ones(n) - eye(n);
estimation_done = zeros(n, 1);

%%% Initialize variables.
A(:, :, :, 1) = extract_init_topo(A0);
neighbors(:, :, 1) = A0;
tau(:, :, 1) = A0;
sigma(:,1) = params.sigma;


%% Experiment loop
tic;

for k=2:steps

    % 0) Preparation actions not part of the algorithm
    l2_ref = get_lambda_ref(k, l2_refs);  % Take the new reference value.
    signal_strength = change_signal_strength(k, signal_strength, signal_refs);
    temp_norms = zeros(n, 1);

    % 1) Send the current estimate to neighbors.
    % Not needed.

    % 2) Wait for neighbors' responses.
    % Not needed.

    % 3) Measure link qualities.
    randoms = ones(n, n) * rand(rand_stream, 1);
    for l=1:n
        tau(l,:,k) = simulate_tau(l, neighbors(:,:,k-1), tau(l,:,k-1), params, signal_strength, randoms);
    end

    % 3.5) Respond to link changes.
    for l=1:n
        neighbors(l,:,k) = respond_to_changes(k, l, neighbors(:,:,k-1), tau(l,:,k), no_disconnect, no_reconnect, params);
    end
    % --> If we started receiving messages and tau is big enough to add a new neighbor,
    % --> that means we have the state from this new neighbor and we should start including
    % --> it in our calculations. If tau reduces so we disconnect from the neighbor, that
    % --> means we have his (very) old state and we should not include it in calculations
    % --> any more.
    % --> Conclusion: use neighbors(k) in this update and for matrix norm calculation.
    r_changes = record_link_changes(k, neighbors, A0, 'R');

    % 4) Based on received information, update the current estimate.
    for l=1:n
        A(:,:,l,k) = update_state(l, A(:,:,:,k-1), tau(l,:,k), neighbors(l,:,k), sigma(l,k-1), params);
    end

    % 4.1) Calculate local matrix norms.
    for l=1:n
        for j=1:n
            if neighbors(l,j,k)
                local_norms(l, k) = max(local_norms(l, k), norm(A(:,:,l,k-1) - A(:,:,j,k-1), "fro"));
            end
        end
    end
    % --> We only know our own A(k). For all other neighbors, we have A(k-1). This is why
    % --> we must compare matrices from the previous step when calculating the norm.

    % 4.2) Monitor matrix norms to detect convergence and run second-level consensus.
    for l=1:n
        if k > estimation_period && local_norms(l,k) < params.norms_limit
            ready_flags(l,l,k) = 1;

            for j=1:n
                if neighbors(l,j,k-1)
                    ready_flags(l,:,k) = ready_flags(l,:,k) | (ready_flags(j,:,k-1) & ready_mask(l,:)) | failed(l,:);
                end                                            % ^ Don't allow agent j to set the the ready flag of agent l.
            end
            % --> We based our decision about readiness on data from previous step. So,
            % --> we also need to use ready flags and neighbors from the previous step.
        end
        % --> Judging from the last step, we are ready for new control period.

        if all(ready_flags(l,:,k)) || ~any(neighbors(l,:,k))
            ready_all(l,k) = 1;
            estimation_done(l) = 1;
        end
    end

    % TESTME: This was moved from after step 6. It shouldn't affect the
    %   execution in any way, but if there is a problem, this would be the
    %   first place to look.
    % 5) Detect any agents that may have failed.
    for l=1:n
        if k > estimation_period && estimation_done(l) == 1
            [failed(l,:), neighbors(l,:,k), A_ok{l}] = detect_failed(k,l,A(:,:,l,k),failed(l,:),neighbors(l,:,k),params);
        else
            A_ok{l} = A(:,:,l,k);
        end
    end

    % 6) Calculate the current value of lambda_2, sigma, and K_l2
    for l=1:n
        [lambda(l,k), fiedler(:,l,k), lambda_n(l,k)] = calc_lambda(A_ok{l}, failed(l,:));
        sigma(l,k) = calc_sigma(lambda_n(l,k), params);
        K_l2(l,k) = calc_Kl2(l, A_ok{l}, fiedler(:,l,k), K_l2(l,k-1), ready_all(l,k), params);
    end

    % 7) Add or remove links if necesary.
    for l=1:n
        % Trigger connectivity control when agent is ready, but only in specified discrete moments.
        % This helps with the problem when one agent makes a change in step k, and the other one
        % remains ready until step k + x when it makes its change. This can result in removing two
        % links that are important for the connectivity. When using fixed moments, a minimum of n
        % steps must pass between two actions which should be enough time to change the norms and
        % react to already made changes. n is equal to the maximum diameter of the graph.
        if k > estimation_period && k > ready_delay(l) && l2_ref > 0 && ready_all(l,k) && mod((k - estimation_period), params.kappa) == 0
            [neighbors(l,:,k), no_disconnect, no_reconnect, changed] = control_lambda( ...
                k, l, lambda(l,k), l2_ref, K_l2(l,k), fiedler(:,l,k), A(:,:,l,k), neighbors(l,:,k-1), no_disconnect, no_reconnect, params);
            % If link was added or removed, reset your ready flag after 10
            % steps (to allow other agents to become ready as well).
            % I guess this is only in case when local norm doesn't leave
            % the specified threshold.
            if ~isempty(changed)
                ready_delay(l) = k + 10;
            end
        end

        if k == ready_delay(l)
            ready_flags(l,:,k) = zeros(1,n);
        end
    end

    % Record the link changes.
    record_link_changes(k, neighbors, A0, 'D', r_changes);

    % Simulate different handshake options.
    [neighbors(:,:,k), tau(:,:,k), A(:,:,:,k)] = simulate_handshake(k, neighbors, tau(:,:,k), A(:,:,:,k), params);

end

elapsed = toc

%% Plotting and anaysis

% Override plotting settings
% --------------------------
% z_doPlot.Graph = true;       % Visual representation of the graph.

% z_doPlot.Lambda = true;      % Lambda2 through time.
% z_doPlot.K_l2 = true;        % K_l2 thresholds on the lambda2 plot.
% z_doPlot.subSigma = true;    % Value of sigma below the lambda2 plot.
% z_doPlot.avgMatrix = true;   % Average link qualities below the lambda2 plot.
% z_doPlot.subLink = true;     % Selected link quality below the lambda2 plot.

% z_doPlot.Tau = true;         % Tau through time.
% z_doPlot.Sigma = true;       % Value of sigma on a separate plot.
% z_doPlot.LambdaN = true;     % LambdaN through time.

% z_doPlot.Links = true;       % Individual link qualities through time - selectable.
% z_doPlot.Matrix = true;      % All links of a given matrix through time.
% z_doPlot.Snapshots = true;   % Snapshots of the graph through time.

% z_doPlot.LocalNorm = true;   % Local norms through time.
% z_doPlot.ReadyFlags = true;  % Ready flags through time.

z_shouldPlot = @(var) ~isempty(var) && isfield(z_doPlot, var) && z_doPlot.(var);

z_dispNames = sprintfc('Agent {%d}',1:n);
if ~exist('z_figures', 'var')
    z_figures = struct();
end

z_labelOptions = {'interpreter','latex','FontSize',14};
z_figureOptions.labels = z_labelOptions;

% Plot initial topology graph.
[z_show, z_figures] = make_figure(z_doPlot, z_figures, 'Graph');
if z_show
    z_G0 = graph(A0);
    z_plotHandle = plot(z_G0,'Layout','auto','MarkerSize',18,'LineWidth',1,'EdgeColor','k','EdgeAlpha',1);
    text(z_plotHandle.XData, z_plotHandle.YData, z_plotHandle.NodeLabel, ...
        'VerticalAlignment','middle',...
        'HorizontalAlignment', 'center',...
        'Color','k','FontWeight','bold',...
        'FontSize', 10)
    % title('Initial topology')
    z_plotHandle.NodeLabel = {};
    z_axisHandle = gca;
    set(z_axisHandle, 'XTick', [], 'XTickLabel', []);
    set(z_axisHandle, 'YTick', [], 'YTickLabel', []);
    set(get(z_axisHandle, 'XAxis'), 'Visible', 'off');
    set(get(z_axisHandle, 'YAxis'), 'Visible', 'off');
end

% Plot values of lambda2 through time.
[z_show, z_figures] = make_figure(z_doPlot, z_figures, 'Lambda');
if z_show
    if z_shouldPlot('subSigma') || z_shouldPlot('avgMatrix') || z_shouldPlot('subLink')
        z_defaultFigPosition = get(groot, 'defaultfigureposition');
        if z_figures.Lambda.Position(4) == z_defaultFigPosition(4)
            z_figures.Lambda.Position(4) = z_figures.Lambda.Position(4) * 1.5;
        end
        tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
        nexttile;
    end
    l2_refs = [l2_refs, [steps; l2_refs(2,end)]];
    z_LineObject = plot(lambda', 'LineWidth', 1.5);
    set(z_LineObject, {'DisplayName'}, z_dispNames')
    hold on
    stairs(l2_refs(1,2:end), l2_refs(2,2:end), 'k--', 'DisplayName','ref')
    % Additionally plot K_l2 thresholds on the same axis.
    z_legendEntries = n + 1;
    if z_shouldPlot('K_l2')
        [z_etx, z_error_threshold] = calculate_error_thresholds(l2_refs, K_l2(1,:));
        plot(z_etx, z_error_threshold(1,:), 'r:', 'DisplayName','+K_{l2}')
        plot(z_etx, z_error_threshold(2,:), 'b:', 'DisplayName','-K_{l2}')
        z_legendEntries = z_legendEntries + 2;
    end
    hold off
    grid on
    xlim([0 steps])
    ylim([0 1.05 * max(lambda, [], 'all')])
    xlabel('Step (k)', z_labelOptions{:})
    ylabel('Algebraic connectivity $\lambda_2$', z_labelOptions{:})
    % title('\lambda_2 through steps')
    legend('Location','best','NumColumns',ceil(z_legendEntries/6),'FontSize',10)

    % Additionally plot sigma through time on the second subplot.
    if z_shouldPlot('subSigma')
        nexttile
        z_LineObject = plot(sigma', 'LineWidth', 1.5);
        set(z_LineObject, {'DisplayName'}, z_dispNames')
        grid on
        xlabel('Step (k)', z_labelOptions{:})
        ylabel('Update step size $\sigma$', z_labelOptions{:})
        % title('\\sigma through steps')
        legend('Location','southeast','NumColumns',floor(n/5),'FontSize',10)
        % ylim([0.05, 0.45])
    end

    % Additionally plot average link qualities through time on the second subplot.
    if z_shouldPlot('avgMatrix')
        nexttile
        for i=1:n
            for j=(i+1):n
                plot(squeeze(mean(A(i,j,:,:), 3)), 'DisplayName', sprintf("$\\bar{a}_{%d%d}$", i, j))
                hold on
            end
        end
        hold off
        grid on
        ylim([-0.1, 1.1])
        xlabel("Step (k)", z_labelOptions{:})
        ylabel("Average link weight $\bar{a}_{ij}$", z_labelOptions{:})
        legend('Location','best','Interpreter','latex','NumColumns',floor((n*n-n)/2/5),'FontSize',12)
    end

    % Additionally plot individual link quality through time on the second subplot.
    if z_shouldPlot('subLink')
        nexttile
        i = 1; j = 4;  % Select links.
        for l=1:n
            plot(squeeze(A(i, j, l, :)), 'DisplayName', sprintf("$a^%d_{%d%d}$", l, i, j))
            hold on
        end
        hold off
        grid on
        xlabel("Step (k)", z_labelOptions{:})
        ylabel("Link weight $a_{ij}$", z_labelOptions{:})
        legend('Location','best','Interpreter','latex','NumColumns',floor(n/5),'FontSize',12)
    end
end

% Plot tau through time
[z_show, z_figures] = make_figure(z_doPlot, z_figures, 'Tau');
if z_show
    t = tiledlayout(n,1,'TileSpacing','Compact','Padding','Compact');
    for i=1:n
        nexttile
        plot(squeeze(tau(i,:,:)).')
        ylim([-0.1, 1.1])
        legend('Location','eastoutside','NumColumns',ceil(n/5))
        ylabel(strcat('$\tau_', num2str(i), '$'), z_labelOptions{:})
    end
    title(t, 'Link measurements through time')
    xlabel(t, 'Step (k)', z_labelOptions{:})
end

% Plot sigma through time
[z_show, z_figures] = make_figure(z_doPlot, z_figures, 'Sigma');
if z_show
    plot(sigma')
    grid on
    xlabel('Step (k)', z_labelOptions{:})
    ylabel('Update step size $\sigma$', z_labelOptions{:})
    title('\sigma through steps')
    legend('Location','best','NumColumns',floor(n/5))
    % ylim([0.05, 0.45])
end

% Plot values of lambdaN through time.
[z_show, z_figures] = make_figure(z_doPlot, z_figures, 'LambdaN');
if z_show
    plot(lambda_n')
    grid on
    xlim([0 steps])
    xlabel('Step (k)', z_labelOptions{:})
    ylabel('$\lambda_n$', z_labelOptions{:})
    title('\lambda_n through steps')
    legend('Location','best','NumColumns',floor(n/5))
end

% Plot individual link qualities through time.
[z_show, z_figures] = make_figure(z_doPlot, z_figures, 'Links');
if z_show
    z_selectedLink_i = [1, 1];
    z_selectedLink_j = [4, 5];
    plot_selected_links(z_figures.Links, A, z_selectedLink_i, z_selectedLink_j, z_figureOptions);
end

% Plot all links of a given matrix through time.
[z_show, z_figures] = make_figure(z_doPlot, z_figures, 'Matrix');
if z_show
    l = 1;
    for i=1:n
        for j=(i+1):n
            plot(squeeze(A(i,j,l,:)), 'DisplayName', sprintf("a^{%d}_{%d,%d}", l, i, j))
            hold on
        end
    end
    hold off
    grid on
    ylim([-0.1, 1.1])
    xlabel('Step(k)', z_labelOptions{:})
    ylabel('Link weight $a_{ij}$', z_labelOptions{:})
    legend('Location','southoutside','NumColumns',ceil((n*n-n)/2/5))
end

% Plot snapshots of the graph through time.
[z_show, z_figures] = make_figure(z_doPlot, z_figures, 'Snapshots');
if z_show
    z_defaultFigPosition = get(groot, 'defaultfigureposition');
    if z_figures.Snapshots.Position(3:4) == z_defaultFigPosition(3:4)
        z_figures.Snapshots.Position(3) = z_figures.Snapshots.Position(3) * 1.5;
        z_figures.Snapshots.Position(4) = z_figures.Snapshots.Position(4) * 1.5;
    end
    z_GraphLayout = matlab.graphics.chart.primitive.GraphPlot('BasicGraph',MLGraph(graph(A0)), 'Layout', 'force');
    % z_GraphLayout = 'auto';
    plot_snapshots(lambda, A, 0, [150, 250, 350, 450, 550, 700], z_GraphLayout);
end

% Plot local norms through time
[z_show, z_figures] = make_figure(z_doPlot, z_figures, 'LocalNorm');
if z_show
    z_LineObject = plot(local_norms');
    hold on
    plot([1, steps], [params.norms_limit, params.norms_limit], 'k--', 'DisplayName', 'threshold')
    hold off
    grid on
    xlabel('Step (k)', z_labelOptions{:})
    ylabel('Norm', z_labelOptions{:})
    title('Local maximum of matrix norms')
    legend('Location','best','NumColumns',floor(n/5))
    set(z_LineObject, {'DisplayName'}, z_dispNames')

%     within_threshold = zeros(n,1);
%     for i=1:n
%         within_threshold(i) = find(local_norms(i,600:end) < 0.0001, 1);
%     end
%     fprintf("Max difference in step of reaching consensus %d\n", max(within_threshold) - min(within_threshold))
end

% Plot ready flags through time
[z_show, z_figures] = make_figure(z_doPlot, z_figures, 'ReadyFlags');
if z_show
    z_LineObject = plot(ready_all');
    xlabel('Step (k)', z_labelOptions{:})
    ylabel('Flag value', z_labelOptions{:})
    title('Ready flags')
    legend('Location','best','NumColumns',floor(n/5))
    set(z_LineObject, {'DisplayName'}, z_dispNames')
end

%% Save the run parameters
save_results;
