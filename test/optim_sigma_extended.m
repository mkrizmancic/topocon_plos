% This script tests our claim that the optimal value of sigma for the 
% system given with G = I - sigma * (I+L), corresponding to additional
% measurements is equal to sigma = 2 / (lambda_n + 2).
%
% By changing the value of sigma to be larger or smaller than the optimal
% value, this script shows that consensus with given G can work with 
% different values, but converges in least number of steps when sigma is
% set to the optimal value.

[A0, ~, ~] = get_topology('dogtail', 10);    % Initial topology

n = size(A0, 1);
steps = 1000;

x = zeros(n, steps);
rng(123)
x(:, 1) = rand(n, 1);
x_avg = mean(x(:, 1));
u_avg = x_avg * 2;

D = diag(sum(A0, 2));
L = D - A0;

[~, V] = eig(L);
[v, ~] = sort(diag(V));
sigma_o(1) = 2 / (v(end) + 2);
sigma_s(1) = 2 / (v(end) + 1);
ro(1) = v(end) / (v(end) + 2);

sigma_stable = sigma_s(1);
sigma_optimal = sigma_o(1);
sigma = sigma_stable * 1.1;

G = eye(n) - sigma * (eye(n) + L);
done = 0;

for k=1:steps-1
    % Only the first agent can measure something.
    x(:, k+1) = G * x(:,k) + sigma * [u_avg; x(2:end, k)];
    % All agents measure the same thing.
    % x(:, k+1) = G * x(:,k) + sigma * u_avg * ones(n, 1);
    if all(abs(x(:,k+1) - u_avg) <= 0.01 * u_avg) && done == 0
        fprintf("Consensus reached in step %d\n", k)
        done = 1;
    end
end

figure(1)
plot(x')
grid on
xlabel('step (k)')
ylabel('x')
title('values through steps')
legend('Location','best','NumColumns',floor(n/5))

% T = table(sigma_s', sigma_o', ro', RowNames={'eigen', 'hurwitz', 'fmincon'}, VariableNames={'stable', 'optimal', 'rate'})
