% This script tests our claim that the optimal value of sigma for the 
% system given with G = I - sigma * (I+L), corresponding to additional
% measurements is equal to sigma = 2 / (lambda_n + 2). Additionally, it tests
% this on a system with state modeled as vector.
%
% By changing the value of sigma to be larger or smaller than the optimal
% value, this script shows that consensus with given G can work with 
% different values, but converges in least number of steps when sigma is
% set to the optimal value.

[A0, ~, ~] = get_topology('goodit_15_ok', 10);    % Initial topology

n = size(A0, 1);
m = 2;
steps = 100;


x = zeros(n, m, steps);
rng(123)
x(:, :, 1) = rand(n, m);
x_avg = mean(x(:, :, 1));
a = 0.9; b = 1.1;
u_avg = (a + (b-a)*rand(1,m)) .* randsample([-1, 1], 2) .* x_avg;


D = diag(sum(A0, 2));
L = D - A0;

[~, V] = eig(L);
[v, ~] = sort(diag(V));
sigma_o(1) = 2 / (v(end) + 2);
sigma_s(1) = 2 / (v(end) + 1);
ro(1) = v(end) / (v(end) + 2);


sigma_stable = sigma_s(1);
sigma_optimal = sigma_o(1);
sigma = sigma_optimal * 1;

G = eye(n) - sigma * (eye(n) + L);

done = 0;

for k=1:steps-1
    x(:,:,k+1) = G * x(:,:,k) + sigma * u_avg .* ones(n, 2);
    if all(abs(x(:,:,k+1) - u_avg) <= 0.01 * abs(u_avg), 'all') && done == 0
        fprintf("Consensus reached in step %d\n", k)
        done = 1;
    end
end

figure(1)
subplot(211)
plot(squeeze(x(:,1,:))')
grid on
xlabel('step (k)')
ylabel('x')

subplot(212)
plot(squeeze(x(:,2,:))')
grid on
xlabel('step (k)')
ylabel('x')

title('values through steps')
legend('Location','best','NumColumns',floor(n/5))

% T = table(sigma_s', sigma_o', ro', RowNames={'eigen', 'hurwitz', 'fmincon'}, VariableNames={'stable', 'optimal', 'rate'})
