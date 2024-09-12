function new_neighbors = respond_to_changes(k, me, neighbors, tau, no_disconnect, no_reconnect, params)
%RESPOND_TO_CHANGES Simulate link quality changes and respond to them accordingly.
%
%   Input:
%       * k - current step
%       * me - id of the current agent
%       * neighbors - complete neighbors matrix [n x n]
%       * tau - tau vector of the current agent [1 x n]
%       * no_disconnect - matrix of steps before which agents should not disconnect (because they recently connected)
%       * no_reconnect - matrix of steps before which agents should not reconnect (because they recently disconnected)
%       * params - simulation parameters
%   Output:
%       * new_neighbors - updated neighbors matrix [1 x n]

    n = params.n;

    % Respond to changes in link qualities.
    for j=1:n
        % This agent was a neighbor of 'j', but link got very bad so it should disconnect.
        if neighbors(me, j) == 1 && tau(j) < 0.2 && k > no_disconnect(me, j)
            neighbors(me, j) = 0;
        % This agent was not a neighbor of 'j', but suddenly started receiving messages so it should connect.
        elseif neighbors(me, j) == 0 && tau(j) > 0.2 && k > no_reconnect(me, j)
            neighbors(me, j) = 1;
        end
    end

    new_neighbors = neighbors(me, :);
end