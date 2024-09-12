% This script analyzes and plots the root locus of systems given with 
% a) G = I - sigma * L 
% b) G = I - sigma * (I+L). 
% The systems are stable if all of their poles given by the eigenvalues of L are
% in the unit circle. The eigenvalues are modified with update parameter sigma. 
% Depending on the system, different expressions define optimal and stable
% values of sigma.
%
% All eigenvalues of system b) are larger by 1 than the eigenvalues of a).
% The smallest eigenvalue of a) is not plotted because it is always zero.

all_fig = findall(0, 'type', 'figure');
close(all_fig)

[A0, ~, ~] = get_topology('dogtail', 15);    % Initial topology
n = size(A0, 1);
D = diag(sum(A0, 2));
L = D - A0;

[~, V] = eig(L);
[v, ~] = sort(diag(V));

% Regular system
% sigma_o = 2 / (v(end) + v(2)) * 1.0
% sigma_s = 2 / (v(end))
% ro_o = 1 - sigma_o * v(2)
% lambdas = v(2:end);

% Modified system
sigma_o = 2 / (v(end) + 2)
sigma_s = 2 / (v(end) + 1)
ro_o = 1 - sigma_o
lambdas = v + 1;

fig = uifigure('Name', 'Root Locus Slider');
fig.Position(3:4) = [1000 100];
g1 = uigridlayout(fig, [1 1]);
g1.RowHeight = {'1x'};
g1.ColumnWidth = {'1x'};

sl = uislider(g1);

sl.Layout.Row = 1;
sl.Layout.Column = 1;

sl.Limits = [0, 1.2];
sl.Value = sigma_o/sigma_s;
sl.MajorTicks = [0:0.2:1.2 sigma_o/sigma_s];
sl.MinorTicks = 0:0.02:1.2;

ko = find(sigma_o < sl.MajorTicks * sigma_s, 1) - 1;
sl.MajorTickLabels{ko} = '*';

sl.ValueChangedFcn = @(src, event) updateSigma(src, event, lambdas, sigma_s);

dummyE.Value = sl.Value;
updateSigma([], dummyE, lambdas, sigma_s)


function updateSigma(src, event, eigvals, sigma)

    sigma_value = event.Value * sigma
    den = [1 -1];
    sys = tf(sigma_value, den, -1);

    r = rlocus(sys, eigvals);
    ro2 = max(abs(r));

    figure(1)
    rlocus(sys, eigvals)
    hold on
    plot(real(r), imag(r), 'ro')
    for k=1:size(r, 2)
        text(real(r(k))+0.03,imag(r(k))-0.03,num2str(k),'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle')
    end
    rectangle('Position', [0-ro2 0-ro2 2*ro2 2*ro2],'Curvature', [1 1], 'EdgeColor','g')
    hold off
    legend('sys', 'sys poles @ \lambda_i')
    % zgrid
    axis equal
end

