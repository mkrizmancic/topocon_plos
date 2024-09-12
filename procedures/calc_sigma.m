function [sigma] = calc_sigma(lambda_n, params)
%CALC_SIGMA Calculate the optimal value of consensus update step (sigma).
%
%   Input:
%       * lambda_n - The largest eigenvalue of the Laplacian matrix
%       * params - Struct containing parameters
%   Output:
%       * sigma - Optimal value of consensus update step

    if ~params.adapt_sigma
        sigma = params.sigma;
        return
    end

    sigma = 2 / (lambda_n + 2);
end