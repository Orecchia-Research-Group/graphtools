function [ edgesCut, L, R] = kernighanLin(filename, eta)
%KERNIGHANLIN Kernighan-Lin greedy algorithm
%   Given a graph. produce a bisection and iteratively move the node with
%   the largest crossing degree to the other side.

    [G, n, ~] = loadeg2graph(filename);
overlappingNodes = [];
D = diag(sum(G));
v = round(rand(n,1));
u = expv((-1) * eta, (D - G), v);
[~, index] = sort(u);
Ltemp = zeros(n, 1, 'logical');
Ltemp(index(1:floor(n / 2))) = true;
Rtemp = ~Ltemp;
edgeCut = sum(sum(G(Ltemp, Rtemp)));
crossing_edges = zeros(n, 1);
crossing_edges(Ltemp, :) = sum(G(Rtemp, Ltemp));
crossing_edges(Rtemp, :) = sum(G(Ltemp, Rtemp));
non_crossing_edges = zeros(n, 1);
non_crossing_edges(Ltemp) = sum(G(Ltemp, Ltemp));
non_crossing_edges(Rtemp) = sum(G(Rtemp, Rtemp));
difference_edges = crossing_edges - non_crossing_edges;
i = 0;
while true
    [m, idx] = max(difference_edges);
    if m <= 0
        break;
    end
    i = i + 1;
    edgeCut = edgeCut - difference_edges(idx);
    same = Ltemp == Ltemp(idx);
    diff = ~same;
    difference_edges(same, :) = difference_edges(same, :) + 2 * G(same, idx);
    difference_edges(diff, :) = difference_edges(diff, :) - 2 * G(diff, idx);
    Ltemp(idx) = ~Ltemp(idx);
    Rtemp(idx) = ~Rtemp(idx);
    difference_edges(idx) = -difference_edges(idx);
    % fprintf('%6d: %d %6d %3d %6d %6d\n', i, Ltemp(idx), idx, -difference_edges(idx), full(cutedges), full(sum(sum(G(Ltemp, Rtemp)))));
end
Lmask = Ltemp;
Rmask = Rtemp;

L = sort(int64(find(Lmask)));
R = sort(int64(find(Rmask)));
edgesCut = sum(sum(G(L, R)));

end
