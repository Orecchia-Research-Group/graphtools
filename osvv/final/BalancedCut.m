%% MATLAB FUNCTION: BalancedCut
%
% PURPOSE:  Run cutfind iteratively, breaking up the largest component every
%           time until cluster_count is reached.
%
% INPUTS: 
%   The rest are as described in cutfind
% 
% OUTPUTS:
%       edgesCut - number of edges in best cuts found
%       L - list of vertices composing best cut found
%       R- total time taken

function [edgesCut, L, R] = BalancedCut(c, FileToRead, outputfile, suffix, t, stop,  eta, init, seed, p, rate, lwbd, certificatespec, varargin)

% Mixed cut or edge cut?
if (size(varargin, 2) > 0)
    lamda = varargin{1};
else
    lamda = -1.0;
end

if(ischar(FileToRead))
    [G, n, ~] = loadeg2graph(FileToRead);
else
    G = FileToRead;
    n = size(G, 1);
end

shared = false(n, 1);
clusters = {[1:n]'};
clusterSizes = [n];

clusterExpansions = [0];
edgesCut = 0;

clusterCountIndex = 1;
[largestClusterSize, largestClusterIndex] = max(clusterSizes);
while largestClusterSize > c * n
    largestClusterNodes = clusters{largestClusterIndex};
    
    largestMask = sparse(largestClusterNodes, 1, true, n, 1);
    sharedBefore = largestMask & shared;
    largestNotShared = largestMask & (~shared);
    % fprintf('shared Before: %d. largestNotShared: %d.\n', sum(sharedBefore'), sum(largestNotShared'));
    flowgraph = G(largestNotShared, largestNotShared);
    nodes = find(largestNotShared);
    
    % Check to see if largestSubgraph is disconnected.
    grph = graph(flowgraph);
    bins = conncomp(grph);
    if (length(unique(bins)) > 1)
        a = hist(bins,unique(bins));
        fprintf(2, 'Size: %d\n', a);
    end
    
    if lamda > 0
        [expansionFound, edgeCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, lower] = cutfind(flowgraph,  outputfile, suffix, t, stop,  eta, init, seed, p, rate, lwbd, certificatespec, lamda);
    else
        [expansionFound, edgeCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, lower] = cutfind(flowgraph,  outputfile, suffix, t, stop,  eta, init, seed, p, rate, lwbd, certificatespec);
    end
    
    Lnodes = nodes(L, 1);
    Rnodes = nodes(R, 1);
    Lmask = sparse(Lnodes, 1, true, n, 1);
    Rmask = sparse(Rnodes, 1, true, n, 1);
    Cmask = Lmask & Rmask;
    
    shared = shared | Cmask;

    clusters{end+1, 1} = find(Lmask | sharedBefore);
    clusters{largestClusterIndex, 1} = find(Rmask | sharedBefore);
    
    clusterSizes = [clusterSizes; size(clusters{end, 1}, 1)];
    clusterSizes(largestClusterIndex, 1) = size(clusters{largestClusterIndex, 1}, 1);
    % fprintf('%d %d %d %d\n', full(sum(Lmask & ~shared)), full(sum(~(Lmask | shared))), full(sum(Lmask & ~shared & ~(Lmask | shared))), nnz(G(Lmask & ~shared, ~(Lmask | shared))));
    % fprintf('%d %d %d %d\n', full(sum(Rmask & ~shared)), full(sum(~(Rmask | shared))), full(sum(Rmask & ~shared & ~(Rmask | shared))), nnz(G(Rmask & ~shared, ~(Rmask | shared))));

    clusterExpansions(end+1, 1) = nnz(G(Lmask & ~shared, ~(Lmask | shared)));
    clusterExpansions(largestClusterIndex, 1) = nnz(G(Rmask & ~shared, ~(Rmask | shared)));
    if lamda > 0
        clusterExpansions(end, 1) = clusterExpansions(end, 1) + lamda * sum(Cmask | sharedBefore);
        clusterExpansions(largestClusterIndex, 1) = clusterExpansions(largestClusterIndex, 1) + lamda * sum(Cmask | sharedBefore);
    end
    
    edgesCut = edgesCut + nnz(flowgraph(L, R));

    [largestClusterSize, largestClusterIndex] = max(clusterSizes);
end

fprintf('%d\n', clusterSizes);
[~, sizeIndex] = sort(clusterSizes, 'descend');
L = clusters{sizeIndex(1)};
R = clusters{sizeIndex(2)};
Lmask = sparse(L, 1, true, n, 1);
Rmask = sparse(R, 1, true, n, 1);
if size(sizeIndex, 1) > 3
    for i = sizeIndex(3:end)
        newCluster = clusters{i};
        newMask = sparse(newCluster, 1, true, n, 1);
        Lcount = nnz(G(Lmask & (~shared), newMask & (~shared)));
        Rcount = nnz(G(Rmask & (~shared), newMask & (~shared)));
        if Rcount > Lcount
            Rmask = Rmask | newMask;
        else
            Lmask = Lmask | newMask;
        end
    end
end

end

