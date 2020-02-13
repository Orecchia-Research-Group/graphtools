function [ edgeCut, overlappingNodes, L, R] = KernighanLin( filename, overlapNodes)
%KERNIGHANLIN Kernighan Lin greedy algorithm
%   Iteratively choose the node with the most crosses to include in the
%   overlap


[G, n, ~] = loadeg2graph(filename);
overlappingNodes = [];
D = diag(sum(G));
opts.tol = 0.001;
[temp, ~ ] = eigs(D-G, 2, 'SA', opts);
u=temp(1:n, 2);
[~, index] = sort(u);
Lsize = 0;
Rsize = n;
Ltemp = sparse(index(1), 1, true, n, 1);
Rtemp = ~Ltemp;
cutedges = sum(sum(G(Ltemp, Rtemp)));
best_expansion = cutedges;
Lmask = Ltemp;
Rmask = Rtemp;
for i=2:n-1
    Lsize = Lsize + 1;
    Rsize = Rsize - 1;
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

edgeCut = [edgesCut];

degrees = zeros(n, 1);
degrees(L) = sum(G(R, L));
degrees(R) = sum(G(L, R));

Cmask = Lmask & Rmask;

Lsize = full(sum(Lmask));
Rsize = full(sum(Rmask));

for i=1:overlapNodes
    [~, id] = max(degrees ./ (Lmask * Lsize + Rmask * Rsize));
    maxd = degrees(id);
    if maxd == 0
        break
    end

    degrees(id) = 0;
    edgesCut = edgesCut - maxd;
    edgeCut(end + 1) = edgesCut;
    Cmask(id) = 1;
    if Lmask(id) == 1
        mask = Rmask & ~Cmask;
        degrees(mask) = degrees(mask) - full(G(mask, id));
        Rmask(id) = 1;
        Lsize = Lsize - 1;
    else
        mask = Lmask & ~Cmask;
        degrees(mask) = degrees(mask) - full(G(mask, id));
        Lmask(id) = 1;
        Rsize = Rsize - 1;
    end
    overlappingNodes(end + 1) = double(id);
end

end

