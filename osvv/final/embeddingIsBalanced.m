function [ret] = embeddingIsBalanced(v, weights, t, b)
%EMBEDDINGISBALANCED Check if the embedding is (t, b)-balanced
%   In order for `v` to be (t, b)-balanced the weighted variance of
%   all the node embeddings whose ||v_i||_^2 <= t must account for at least
%   a `b` fraction of the total variance of the embedding.

%% Argument processing

arguments
    v (:, :) double
    weights (1, :) double
    t (1, 1) double
    b (1, 1) double
end

%% Checking condition

r = sum(v .* v, 2);
ind = r <= t;
ret = sum(weights(ind) * r(ind) >= b * weights * r);
end
