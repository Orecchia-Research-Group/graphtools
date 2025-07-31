function [S, T] = roundCut(v, weights, t, b)
%ROUNDCUT Given an (n, d) embedding of the nodes 
%   Detailed explanation goes here

%% Argument processing
arguments
    v (:, :) double
    weights (1, :) double
    t (1, 1) double
    b (1, 1) double
end

[n, d] = size(v);
mu_V = sum(weights);

if embeddingIsBalanced(v, weights, t, b)
    %% Balanced embedding
    % Create random projection
    g = rand(d, 1);
    r = v * g;
    [~, idx] = sort(r, "descend");
    mu_S_up = 0;
    mu_S_down = mu_V;
    S_up_mask = zeros(n, 1, "logical");
    S_down_mask = ones(n, 1, "logical");
    best_S_index = -1;
    for i=1:n
        u = idx(i);
        mu_u = weights(u);
        mu_S_up = mu_S_up + mu_u;
        S_up_mask(u) = true;
        if (t * mu_S_up >= mu_V) && (best_S_index < 0)
            best_S_index = i;
        end
        if t * (mu_S_down - mu_u) < mu_V
            T = find(S_down_mask);
            break;
        end
        mu_S_down = mu_S_down - mu_u;
        S_down_mask(u) = false;
    end
    S = sort(idx(1:best_S_index));
else
    %% Unbalanced embedding
    % Separate based on measure
    r = sum(v .^ 2, 2);
    [~, idx] = sort(r, "descend");
    T = find(r <= t / 2);
    mu_S = 0;
    S_mask = false(n, 1);
    thr = 24 * mu_V / (100 * log(mu_V));    % Compute threshold
    for i=1:n
        u = idx(i);
        mu_S = mu_S + weights(u);
        S_mask(u) = true;
        if mu_S * r(u) >= thr
            S = find(S_mask);
            break;
        end
        if r(u) < t^2
            printf(2, "This wasn't supposed to happen. Refer to Algorithm 4 of https://arxiv.org/pdf/2301.08920.")
            break;
        end
    end
    if exist('S', 'var') == 0
        S = find(r > t / 2);
    end
end
end

