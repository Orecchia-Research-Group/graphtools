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
% R_t
r = sum(v .* v, 2);
ind = r <= t;
if sum(ind) == 0
    ret = false;
    return;
end

% Weighted variance inside inside R_t
v_bar_rt = mean(v(ind), 1, Weights=weights(ind));
rt_var = weights(ind) * sum((v - v_bar_rt).^2, 2);

% mu(R_t^2) >= mu(V^2)
ret = rt_var >= b * weights * r;
end
