function [new_neighbors, new_tau, new_A] = simulate_handshake(k, neighbors, tau, A, params)
%SIMULATE_HANDSHAKE Simulate a handshake between two agents.
%   Handshake defines how the two agents that are modifying their connection
%   share the information about the change between them. Options are:
%       * 'both' - Both agents update their neighbors, tau, or A.
%       * 'active' - Only the agent who made the change updates its neighbors, tau, or A.
%       * 'none' - No information is shared between the agents.

new_neighbors = neighbors(:,:,k);
new_tau = tau;
new_A = A;

options = params.handshake;

for l=1:params.n
    changes = neighbors(l, :, k) - neighbors(l, :, k-1);
    new = find(changes);
    if length(new) > 1
        warning("Found two changes in neighbors vector in step %d", k)
    end


    if ~isempty(new)
        val = changes(new) > 0;

        %% Updating neighbors.
        % Update neighbors lists of of the agent who made a change
        % and the other agent.
        if strcmp(options.neighbors, 'both')
            new_neighbors(new, l) = val;
        end

        % Updating only the agent who made the change is default behavior.

        %% Updating tau.
        % Immediately update tau towards the newly (dis)connected agent.
        if (strcmp(options.taus, 'both') || strcmp(options.taus, 'active'))
            new_tau(l, new) = val;
        end

        % Immediately update tau of other agent that was newly (dis)connected.
        if strcmp(options.taus, 'both')
            new_tau(new, l) = val;
        end

        %% Updating adjacency matrix.
        % Immediately update the link towards the newly (dis)connected agent.
        if (strcmp(options.As, 'both') || strcmp(options.As, 'active'))
            new_A(l, new, l) = val;
            new_A(new, l, l) = val;
        end

        % Immediately update the link of other agent that was newly (dis)connected.
        if strcmp(options.As, 'both')
            new_A(l, new, new) = val;
            new_A(new, l, new) = val;
        end
    end

end


end