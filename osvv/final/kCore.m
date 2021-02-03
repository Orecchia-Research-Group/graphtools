function [new_G, new_n, new_m] =kCore(G, k)
%KCORE Iteratively eliminate all nodes with degree less than or equal to k
%
% Keeps only the largest component.

new_G = G(:, :);
weight = full(sum(G));
id = find(weight > k);
new_n = length(G);
new_m = nnz(G);
while length(id) ~= new_n
    new_G = new_G(id, id);
    new_n = length(id);
    new_m = nnz(G);
    weight = full(sum(new_G));
    id = find(weight > k);
end

[~, C] = graphconncomp(new_G);
id = (C == 1);
new_G = new_G(id, id);
