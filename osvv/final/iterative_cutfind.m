%% MATLAB FUNCTION: iterative_cutfind
%
% PURPOSE:  Run cutfind iteratively, breaking up the largest component every
%           time until cluster_count is reached.
%
% INPUTS: 
%       (int) cluster_count - Desired number of clusters
%   The rest are as described in cutfind
% 
% OUTPUTS:
%       edgesCut - number of edges in best cuts found
%       cutFound - list of vertices composing best cut found
%       endtime - total time taken
%       cuttime - time taken by cutfind calls
%       inittime - time taken by initializations
%       spectime - time taken by spectral computations
%       flowtime - time taken by flow computations

function [edgesCut, cutsFound, sizes, expansions, endtime, cuttimes, inittimes, spectimes, flowtimes] = iterative_cutfind(clusterCount, FileToRead, outputfile, suffix, t, stop,  eta, init, seed, p, rate, lwbd, certificatespec, varargin)

clusterCount = sort(clusterCount);
if(any(floor(clusterCount) ~= clusterCount))
    error('cluster_count must be integer');
end

% Mixed cut or edge cut?
if (size(varargin, 2) > 0)
    lamda = varargin{1};
else
    lamda = -1.0;
end

edgesCut = 0;
cuttimes = 0;
inittimes = 0;
spectimes = 0;
flowtimes = 0;

[G, n, ~] = loadeg2graph(FileToRead);

shared = false(n, 1);
clusters = {[1:n]'};
clusterSizes = [n];
clusterExpansions = [0];

clusterCountIndex = 1;
for k=2:max(clusterCount)
    [~, largestClusterIndex] = max(clusterSizes);
    largestClusterNodes = clusters{largestClusterIndex};
    
    largestMask = full(sparse(largestClusterNodes, 1, true, n, 1));
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
    % Update time counters
    cuttimes = cuttimes + endtime;
    innitimes = inittimes + inittime;
    spectimes = spectimes + spectime;
    flowtimes = flowtimes + flowtime;

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

    if k == clusterCount(clusterCountIndex)
        cutsFound{clusterCountIndex} = clusters;
        expansions{clusterCountIndex} = clusterExpansions ./ clusterSizes;
        sizes{clusterCountIndex} = clusterSizes;
        clusterCountIndex = clusterCountIndex + 1;
    end
end

end

    