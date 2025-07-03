function [mu] = karcher_mean_sphere(X, varargin)
%KARCHER_MEAN_SPHERE Computes the Karcher (intrinsic) mean on the 2-sphere
%
%   mu = karcher_mean_sphere(X)
%   mu = karcher_mean_sphere(X, W)
%   mu = karcher_mean_sphere(X, W, max_iter)
%   mu = karcher_mean_sphere(X, W, max_iter, lambda)
%
%   Inputs:
%       X         - 3xN matrix, each column is a point on the unit sphere (||X(:,i)|| = 1)
%       W         - 1xN weight vector (optional, default: uniform weights)
%       max_iter  - maximum number of iterations (optional, default: 10000)
%       lambda    - learning rate (optional, default: 0.7)
%
%   Output:
%       mu - estimated Karcher mean on the unit sphere
%
%   This function uses a gradient descent scheme based on the Riemannian
%   exponential and logarithmic maps to compute the intrinsic mean.

% Dimensions
[n, N] = size(X);

% Default parameters
lambda = 0.7;
max_iter = 10000;
W = ones(1, N) / N;

% Optional input arguments
if nargin >= 4
    lambda = varargin{3};
end
if nargin >= 3
    max_iter = varargin{2};
end
if nargin >= 2
    W = varargin{1};
    W = W / sum(W);  % normalize weights
end

% Initialization
iter = 0;
error = realmax;
mu = X(:, 1);  % initial mean

% Gradient descent
while true
    mu_old = mu;
    error_old = error;
    gradient = zeros(n, 1);
    
    % Compute weighted gradient
    for i = 1:N
        gradient = gradient + W(i) * logmap_sphere(mu, X(:, i));
    end
    gradient = 2 * gradient / N;
    
    % Update mean using exponential map
    mu = expmap_sphere(mu, lambda * gradient);
    
    % Compute error (intrinsic variance)
    error = 0;
    for i = 1:N
        error = error + W(i) * acos(mu' * X(:, i))^2;
    end
    error = error / N;
    
    % Accept or reject the update
    if error < error_old
        iter = iter + 1;
    else
        mu = mu_old;
        error = error_old;
        lambda = 0.95 * lambda;  % reduce learning rate
    end
    
    % Stopping conditions
    if (iter >= max_iter) || (norm(gradient) < 1e-4) || (lambda < 0.9)
        break;
    end
end

end
