function [A, B] = cutPlayer(H, i, weight, factor, sparse_deg, options)
%CUTPLAYER Cut player - Produces an S-T cut to route a matching over
%   

%% Argument validation
arguments
    H (:, :) {mustBeGraph}
    i (1, 1) int64 {mustBePositive}
    weight (1, :) int64 {mustBeNonnegative}
    factor (:, :) {mustBeGraph}
    sparse_deg (:, :) {mustBeGraph}
    options.init (1, 1) double {mustBeNonnegative} = 1
    options.rate (1, :) char {mustBeMember(options.rate, {'d', 'n', 'infty', 'KL'})} = 'n'
    options.eta  (1, 1) double {mustBePositive} = 0.5
    options.pwr_k (1, 1) int64 {mustBePositive} = 1
    options.embedding_dim (1, 1) int64 {mustBePositive} = 1
    options.directed (1, 1) logical = false
    options.t (1, 1) double = 3
    options.b (1, 1) double = 1/10;
end

init = options.init;
rate = options.rate;
eta = options.eta;
pwr_k = options.pwr_k;
embedding_dim = options.embedding_dim;
directed = options.directed;
t = options.t;
b = options.b;

%% Generate embedding, round the cut and if certificate is directed execute directed round cut

v = generateEmbedding(H, i, weight, factor, sparse_deg, init=init, rate=rate, eta=eta, pwr_k=pwr_k, embedding_dim=embedding_dim);
[S, T] = roundCut(v, weight, t, b);
if directed
    [A, B] = directedRoundCut(S, T, v, weight);
else
    A = S;
    B = T;
end
end

