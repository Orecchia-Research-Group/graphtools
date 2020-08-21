function [ edgesCut, L, R] = sweepCut( filename, eta)
%SWEEPCUT Sweep Cut greedy algorithm
%   Given a node ordering produce an locally optimal sweep cut.
%   Iteratively choose the node with the most crosses to include in the overlap


[G, n, ~] = loadeg2graph(filename);
start = 1;
overlappingNodes = [];
degree = sum(G);
D = diag(degree);
% opts.tol = 0.001;
% [temp, ~ ] = eigs(D-G, 2, 'SA', opts);
v = round(rand(n,1));
u = expv((-1)*eta, (D - G), v);
[~, index] = sort(u);
Lsize = start;
Rsize = n - Lsize;
Ltemp = zeros(n, 1, 'logical');
Ltemp(index(1:Lsize), 1) = true;
Rtemp = ~Ltemp;
cutedges = sum(sum(G(Ltemp, Rtemp)));
best_expansion = cutedges;
Lmask = Ltemp;
Rmask = Rtemp;
for i= start + 1:n - start
    Lsize = Lsize + degree(index(i));
    Rsize = Rsize - degree(index(i));
    cutedges = cutedges - sum(sum(G(index(i), Ltemp)));
    Ltemp(index(i), 1) = true;
    Rtemp(index(i), 1) = false;
    cutedges = cutedges + sum(sum(G(index(i), Rtemp)));
    expansion = cutedges / min(Lsize, Rsize);
%     if (i < 100) || (i > n - 100) || ((i > n /2 -50) && (i < n/2 + 50))
%         fprintf('%d %f\n', i, full(expansion));
%     end
    if expansion < best_expansion
        best_expansion = expansion;
        Lmask = Ltemp;
        Rmask = Rtemp;
    end
end

L = sort(int64(find(Lmask)));
R = sort(int64(find(Rmask)));
edgesCut = sum(sum(G(L, R)));

end

