function [stable_sigma, optim_sigma, ro] = kia_T1_stable(L)

n = size(L, 1);

% valid_sigma = 0;
% sigma = 0.01;
% ro = 0;
% 
% while ro < 1
%     valid_sigma = sigma;
%     sigma = sigma + 0.01;
%     ro = norm(eye(n) - sigma * L - ones(n)/n);
% end

fun = @(sigma)(-sigma);
x0 = 0.001;
options = optimoptions('fmincon','Display','off');
stable_sigma = fmincon(fun,x0,[],[],[],[],[],[],@stability,options);

    function [c,ceq] = stability(x)
    c = norm(eye(n) - x * L - ones(n)/n) - 1;
    ceq = [];
    end


fun = @(x)(norm(eye(n) - x * L - ones(n)/n));
x0 = 0.001;
options = optimoptions('fmincon','Display','off');
[optim_sigma, ro] = fmincon(fun,x0,[],[],[],[],[],[],[],options);


end


