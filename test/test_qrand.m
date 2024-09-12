n=1; % Looks good for 3, but 4 uncovers the error.
steps=1000;

%% Wrong approach
wrong_result = zeros(n, n, steps);
wrong_stream = qrandstream('halton',1,'Skip',1e4,'Leap',1e3);

for k=1:steps
    for i=1:n
        for j=1:n
            r = rand();
            wrong_result(i,j,k) = r;
        end
    end
end

figure(1)
plot_results(wrong_result, n, steps)
figure(2)
plot_histogram(wrong_result, n)

%% Correct approach
correct_result = zeros(n, n, steps);
correct_stream = qrandstream('halton',n,'Skip',1e6,'Leap',1e3);

for k=1:steps
    randoms = qrand(correct_stream, n);
    % randoms = reshape(randoms(randperm(n*n)), n, n);
    for i=1:n
        for j=1:n
            correct_result(i,j,k) = randoms(i,j);
        end
    end
    % or better
    % correct_result(:,:,k) = randoms;
end

figure(3)
plot_results(correct_result, n, steps)
figure(4)
plot_histogram(correct_result, n)

%%
figure(99)
plot_results(test_save, n, steps)


function plot_results(result, n, steps)
reshaped = squeeze(reshape(result, [], 1, steps));
plot(reshaped')

legend_entries = cell(1, n^2);
for i = 1:n
    for j = 1:n
        legend_entries{(i-1)*n+j} = [num2str(i), num2str(j)];
    end
end
legend(legend_entries, 'Location', 'best');
end

function plot_histogram(result, n)
tiledlayout(n, n, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n
    for j = 1:n
        nexttile
        histogram(squeeze(result(i,j,:)))
        title([num2str(i), num2str(j)])
    end
end

end

