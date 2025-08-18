function [u] = generateEmbedding(H, i, weights, factor, sparse_deg, options)
%GENERATEEMBEDDING Given the dual certificate and arguments creates a real
% valued embedding for the nodes.
%
%   Detailed explanation goes here

%% Argument validation
arguments
    H (:, :) {mustBeGraph}
    i (1, 1) double {mustBePositive}
    weights (1, :) int64 {mustBeNonnegative}
    factor (:, :) {mustBeGraph}
    sparse_deg (:, :) {mustBeGraph}
    options.init (1, 1) double {mustBeNonnegative} = 1
    options.rate (1, :) char {mustBeMember(options.rate, {'d', 'n', 'infty', 'KL'})} = 'n'
    options.eta  (1, 1) double {mustBePositive} = 0.5
    options.pwr_k (1, 1) int64 {mustBePositive} = 1
    options.embedding_dim (1, 1) int64 {mustBePositive} = 1
end

dweights = double(weights);
init = options.init;
rate = options.rate;
eta = options.eta;
pwr_k = options.pwr_k;
embedding_dim = options.embedding_dim;


%% Computing current eta
n = size(H, 1);

% LEARNING RATE INITIALIZATION
if(strcmp(rate,'d'))
    current_eta = eta*sqrt(8*log(n)/i);
elseif (strcmp(rate, 'KL'))
    current_eta = eta*sqrt(8*entr/i);
else    
    current_eta = eta;
end

% SPECTRAL PARTITIONING
half_mu = diag(diag(factor).^(-1));
%% SECOND EIGENVALUE
if(strcmp(rate,'infty'))
    opts.tol = 1e-6;
    ddweights = diag(sparse(double(weights)));
    % [u, ~] = eigs(@(x) (factor * ((init + i - 1) .* sparse_deg - H) * factor * x + sum(half_mu * x) * half_mu * ones(size(x))), n, sparse_deg, pwr_k, 'SA', opts);
    [u, ~] = eigs(@(x) (((init + i - 1) .* sparse_deg - H) * x + (ddweights * x)' * ddweights * ones(size(x))), n, ddweights, pwr_k, 'SA', opts);
    u(:) = factor * u;
end
%% Parallel vector cut/matching
if(~strcmp(rate, 'infty'))
    % RANDOM BISECTION INITIALIZATION;
    s = round(rand(n, pwr_k * embedding_dim));
    s = s - mean(s);
    
    % s = factor * s;

    M = factor * ((init + i - 1) .* sparse_deg - H) * factor;
    for step=1:pwr_k * embedding_dim
        %%%  RANDOM WALK STEP 
        u(:, step) = factor * expv((-1)*current_eta, M, s(:, step), 1e-4);
        % u(:, step) = factor * expmv(M, s(:, step), (-1)*current_eta);
    end
    u = reshape(u, n, pwr_k, embedding_dim);

    %%% CENTER
    u = u - mean(u, 1, Weights=dweights);

    %%% RESCALE V for better tolerance
    
    for i=1:pwr_k
        u(:, i, :) = u(:, i, :) / norm(half_mu * u(:, i, :), 'fro');
    end
end

