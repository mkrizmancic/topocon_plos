function tau = simulate_tau(me, neighbors, tau, params, signal_strength, randoms)
%SIMULATE_TAU Simulate link quality changes.
%
%   Input:
%       * me - id of the current agent
%       * neighbors - complete neighbors matrix [n x n]
%       * tau - tau vector of the current agent [1 x n]
%       * params - various algorithm parameters
%   Output:
%       * tau - updated tau vector

    if strcmp(params.tau_method, 'simple')
        tau = simulate_tau_simple(me, neighbors, tau, params, signal_strength, randoms);
    elseif strcmp(params.tau_method, 'ema')
        tau = simulate_tau_ema(me, neighbors, tau, params, signal_strength, randoms);
    elseif strcmp(params.tau_method, 'force')
        tau = simulate_tau_force(me, neighbors, tau, params, signal_strength);
    else
        error('Unknown tau method');
    end
end

function tau = simulate_tau_simple(me, neighbors, tau, params, signal_strength, randoms)
    persistent history;
    persistent old_tau;

    if isempty(history)
        history = -1 * ones(params.n, params.n, params.tau_window);
        old_tau = -1 * ones(params.n, params.n);
    end
    if history(me, 1, 1) == -1
        history(me, :, :) = repmat(tau, 1, 1, params.tau_window);
    end

    for j=1:params.n
        % If the output tau from the last call changed, simulate_handshake function overwrote it.
        % We have to update the history to match that change.
        if old_tau(me, j) ~= -1 && old_tau(me, j) ~= tau(j)
            history(me, j, :) = tau(j);
        end

        s = neighbors(j, me) == 1 && randoms(j,me) < signal_strength(j,me);
        history(me, j, :) = cat(3, history(me, j, 2:end), s);
        tau(j) = sum(history(me, j, :)) / params.tau_window;

        old_tau(me, j) = tau(j);
    end

end

function tau = simulate_tau_ema(me, neighbors, tau, params, signal_strength, randoms)
    alpha = 1 - exp(-5/params.tau_window);

    for j=1:params.n
        s = neighbors(j, me) == 1 && randoms(j,me) < signal_strength(j,me);
        tau(j) = min(max(alpha*s + (1-alpha)*tau(j), 0), 1);
    end
end

function tau = simulate_tau_force(me, neighbors, tau, params, signal_strength)
    for j=1:params.n
        tau(j) = neighbors(j, me) * signal_strength(j, me);
    end
end