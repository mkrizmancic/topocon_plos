figure(10)

q = qrandstream('halton', 1,'Skip',1e3,'Leap',1e2);
r = rand(q,1000,1);

x = zeros(1, length(r));
x(r < 0.9) = 1;
y1 = zeros(1, length(x));
y2 = zeros(1, length(x));
w = 20;
alpha = 1 - exp(-5/w);
buffer = zeros(1, w);


for k=2:length(x)
    y1(k) = alpha * x(k) + (1-alpha) * y1(k-1);
    buffer = [buffer(2:end) x(k)];
    y2(k) = sum(buffer) / w;
end

scatter(1:length(x), x, '.')
hold on
plot(y1)
plot(y2)
hold off
legend('', 'wema', 'ma')

% Percentage of the original signal achieved after i time constants.
% i=1 63.21 
% i=2 86.47
% i=3 95.02
% i=4 98.17
% i=5 99.33
