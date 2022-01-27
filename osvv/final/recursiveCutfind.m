%% MATLAB FUNCTION: recursiveCutfind
%
% PURPOSE:  Run cutfind recursively, breaking up the largest component every
%           time until cluster_count is reached.
%
% INPUTS: 
%   The rest are as described in cutfind
% 
% OUTPUTS:
%       edgesCut - number of edges in best cuts found
%       L - list of vertices composing best cut found
%       R- total time taken



function [score, clusters] = recursiveCutfind(clusterCount, FileToRead, outputfile, suffix, t, stop,  eta, init, seed, p, pwr_k, rate, lwbd, certificatespec, ufactor, varargin)

% Mixed cut or edge cut?
if (size(varargin, 2) > 0)
    if length(varargin) < 2
        error(error_string, 'lambda');
    end
    lamda_num = int64(varargin{1});
    lamda_den = int64(varargin{2});
    if lamda_num > lamda_den
        error('Lambda needs to be less than or equal to 1');
    end
else
    lamda_num = int64(-1);
    lamda_den = int64(1);
end

if(ischar(FileToRead))
    [G, weight] = loadMetisGraph(FileToRead);
    n = size(G, 1);
    weight = int64(weight);
else
    G = FileToRead;
    n = size(G, 1);
    degree = int64(full(sum(G)));
    weight = ones(1, n, 'int64');
    weight(:) = degree;
end
totalWeight = sum(weight);

clusters = {[1:n]'};
clusterSizes = [n];

for c=2:clusterCount
    [~, largestClusterIndex] = max(clusterSizes);
    largestClusterNodes = clusters{largestClusterIndex};
    flowgraph = G(largestClusterNodes, largestClusterNodes);
    currentWeight = sum(weight(largestClusterNodes));
    currentUfactor = ufactor * totalWeight / currentWeight;
    
    if lamda_num > 0
        [expansionFound, edgeCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, lower] = cutfind(flowgraph,  outputfile, suffix, t, stop,  eta, init, seed, p, pwr_k, rate, lwbd, certificatespec, currentUfactor, lamda_num, lamda_den);
    else
        [expansionFound, edgeCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, lower] = cutfind(flowgraph,  outputfile, suffix, t, stop,  eta, init, seed, p, pwr_k, rate, lwbd, certificatespec, currentUfactor);
    end
    grph = graph(G(largestClusterNodes(R), largestClusterNodes(R)));
    bins = conncomp(grph);
    comp = length(unique(bins));
    if (comp > 1)
        a = hist(bins, unique(bins));
        fprintf(2, 'Size: %d\n', a);
        [~, connIndex] = max(a);
        L = [L; R(bins ~= connIndex)];
        R = R(bins == connIndex);
    end
    grph = graph(G(largestClusterNodes(L), largestClusterNodes(L)));
    bins = conncomp(grph);
    comp = length(unique(bins));
    if (comp > 1)
        a = hist(bins, unique(bins));
        fprintf(2, 'Size: %d\n', a);
        [~, connIndex] = max(a);
        R = [R; L(bins ~= connIndex)];
        L = L(bins == connIndex);
    end
    L = sort(L);
    R = sort(R);
    clusters{end+1, 1} = largestClusterNodes(L);
    clusters{largestClusterIndex, 1} = largestClusterNodes(R);
    clusterSizes = [clusterSizes; size(clusters{end, 1}, 1)];
    clusterSizes(largestClusterIndex, 1) = size(clusters{largestClusterIndex, 1}, 1);
    fprintf(2, 'clusterSize[%d] = %d\n', [1:length(clusterSizes); clusterSizes']);
end

for c=1:clusterCount
    L = clusters{c};
    Rmask = false(n, 1);
    for cr=1:clusterCount
        if cr == c
            continue
        end
        newMask = sparse(clusters{cr}, 1, true, n, 1);
        Rmask = Rmask | newMask;
    end
    R = find(Rmask);
    [~, ~, clusterExpansion(c)] = cutexp(G, int64(lamda_num), int64(lamda_den), int64(weight), int64(L), int64(R));
end
score = min(clusterExpansion);

end

