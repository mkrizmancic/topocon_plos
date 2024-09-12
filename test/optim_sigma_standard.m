% This script tests the claim from (Bertrand, 2013) that the optimal value
% of sigma for the system given with G = I - sigma * L is equal to
% sigma = 2 / (lambda_n + lambda_2).
%
% By changing the value of sigma to be larger or smaller than the optimal
% value, this script shows that consensus with given G can work with 
% different values, but converges in least number of steps when sigma is
% set to the optimal value.
%
% However, this value of sigma is not optimal when additional measurements
% are included in consensus calculation.

[A0, ~, ~] = get_topology('goodit_15_ok', 10);    % Initial topology

n = size(A0, 1);
steps = 3000;

x = zeros(n, steps);
rng(123)
x(:, 1) = rand(n, 1);
x_avg = mean(x(:, 1));

D = diag(sum(A0, 2));
L = D - A0;

[~, V] = eig(L);
[v, ~] = sort(diag(V));

% Values of sigma derived from spectral radius of L.
sigma_o(1) = 2 / (v(end) + v(2)) * 1.0;
sigma_s(1) = 2 / (v(end));
ro(1) = 1 - sigma_o(1) * v(2);

% According to (Kia, 2019, eq. S17), discrete system is in form of
%   x(k+1) = (I + sigma*A) * x(k).
%   * A is a system matrix, not an adjacency matrix.
% According to (Bertrand, 2013), the consensus system is in form of
%   x(k+1) = (I - sigma*L) * x(k).
% That is:
%   x(k+1) = (I - sigma*(D-A)) * x(k)
%   x(k+1) = (I + sigma*(A-D)) * x(k)
sigma_s(2) = kia_s17_stable(A0-D);
sigma_o(2) = -1;
ro(2) = -1;

% According to (Kia, 2019, Theorem 1), discrete static consensus converges
% if ||In - sigma*L - 1n*1n^T|| < 1.
[sigma_s(3), sigma_o(3), ro(3)] = kia_T1_stable(L);

% According to Kristian Hengster-Movric, a true discrete form of consensus
% could be used instead of calculating optimal sigma. It looks like this:
%   x(k+1) = (I+D)^-1 * (I+A) * x(k)
%   x(k+1) = (I - (I+D)^-1 * L) * x(k)
sigma_o_kristian = (eye(n)+D)^-1;

% Select which result to use in the experiment.
sigma_stable = sigma_s(1) * 0.999;
sigma_optimal = sigma_o(3);
G = eye(n) - sigma_optimal * L;
% G = (eye(n) + D)^-1 * (eye(n) + A0);
done = 0;

for k=1:steps-1
    x(:, k+1) = G * x(:,k);
    % if all(abs(x(:,k+1) - x_avg) <= 0.01 * x_avg) && done == 0
    %     fprintf("Consensus reached in step %d\n", k)
    %     done = 1;
    % end
    if ((max(x(:,k+1)) - min(x(:,k+1))) / max(x(:,k+1)) <= 0.01) && done == 0
        gap = abs(mean(x(:,k+1)) - x_avg) * 100;
        fprintf("Consensus reached in step %d. Value: %f (%f, %f %%)\n", k, mean(x(:,k+1)), x_avg, gap)
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

T = table(sigma_s', sigma_o', ro', RowNames={'eigen', 'hurwitz', 'fmincon'}, VariableNames={'stable', 'optimal', 'rate'})
