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

function [edgesCut, cutsFound, endtime, cuttimes, inittimes, spectimes, flowtimes] = iterative_cutfind(cluster_count, FileToRead, outputfile, suffix, t, stop,  eta, init, seed, p, rate, lwbd, certificatespec, varargin)

if(floor(cluster_count) ~= cluster_count)
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

subgraphs = {G};
clusters = {[1:n]};
clusterSizes = [n];

for k=1:cluster_count-1
    [largestClusterSize, largestClusterIndex] = max(clusterSizes);
    largestSubgraph = subgraphs{largestClusterIndex};
    largestClusterNodes = clusters{largestClusterIndex};
    if (~isreal(largestSubgraph))
        fprintf(2, 'Graph is not real!');
    end
    
    
    grph = graph(largestSubgraph);
    bins = conncomp(grph);
    
    if (length(unique(bins)) > 1)
        a = hist(bins,unique(bins));
        fprintf(2, 'Size: %d\n', a);
        % fprintf(2, 'Graph disconnected: %d', length(unique(bins)));
    end
    
    if lamda > 0
        [expansionFound, edgeCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, lower] = cutfind(largestSubgraph,  outputfile, suffix, t, stop,  eta, init, seed, p, rate, lwbd, certificatespec, lamda);
    else
        [expansionFound, edgeCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, lower] = cutfind(largestSubgraph,  outputfile, suffix, t, stop,  eta, init, seed, p, rate, lwbd, certificatespec);
    end
    
    % Update time counters
    cuttimes = cuttimes + endtime;
    innitimes = inittimes + inittime;
    spectimes = spectimes + spectime;
    flowtimes = flowtimes + flowtime;
    
    RGraph = largestSubgraph(R, R);
    LGraph = largestSubgraph(L, L);
    
    fprintf(2, 'size of largest subgraph: %d\n', size(largestSubgraph, 1));
    subgraphs{size(subgraphs, 1)+1, 1} = largestSubgraph(R, R);
    subgraphs{largestClusterIndex, 1} = largestSubgraph(L, L);
    
    clusters{size(clusters, 1)+1, 1} = largestClusterNodes(1, R);
    clusters{largestClusterIndex, 1} = largestClusterNodes(1, L);
    
    clusterSizes = [clusterSizes; size(R, 1)];
    clusterSizes(largestClusterIndex, 1) = size(L, 1);
    
    edgesCut = edgesCut + nnz(largestSubgraph(L, R));
    
    for i=1:size(subgraphs, 1)
        fprintf(2, 'Size of subgraph %d: %d\n', i, size(subgraphs{i}, 1));
    end
end
cutsFound = clusters;


end

