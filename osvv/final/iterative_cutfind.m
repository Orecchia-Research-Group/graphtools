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

function [edgesCut, cutsFound, sizes, expansions, endtime, cuttimes, inittimes, spectimes, flowtimes] = ...
    iterative_cutfind(clusterCount, fileToRead, options)

arguments
    clusterCount (1, :) int64 {mustBePositive}
    fileToRead (1, :) char {mustBeFileOrGraph}
    options.outputfile (1, :) {mustBeFileOrID(options.outputfile, 0, 1)} = 1
    options.suffix (1, :) char = ''
    options.t (1, 1) int16 {mustBePositive} = 100
    options.stop (1, :) int64 {mustBeNumeric, mustBeGreaterThanOrEqual(options.stop, 1)} = 10
    options.eta (1, 1) double {mustBePositive} = 0.5
    options.init (1, 1) double {mustBeNonnegative} = 1
    options.seed (1, 1) {mustBeNumeric} = 0
    options.p (1, 1) int64 {mustBeGreaterThanOrEqual(options.p, 1)} = 1000
    options.pwr_k (1, 1) int64 {mustBePositive} = 1
    options.rate (1, :) char {mustBeMember(options.rate, {'d', 'n', 'infty', 'KL'})} = 'n'
    options.lwbd (1, :) char {mustBeMember(options.lwbd, {'y', 'n', 'ylast'})} = 'n'
    options.matchingAlgorithm (1, :) char {mustBeMember(options.matchingAlgorithm, {'dinic', 'dynamic'})} = 'dinic'
    options.certificateSpec (1, 1) {mustBeNumericOrLogical, mustBeInRange(options.certificateSpec, 0, 1)} = 0
    options.lambda_num  (1, 1) int64 = 1
    options.lambda_den  (1, 1) int64 {mustBePositive} = 1
end

outputfile = options.outputfile;
suffix = options.suffix;
t = int16(options.t);
stop = options.stop;
eta = options.eta;
init = options.init;
seed = options.seed;
p = int64(options.p);
pwr_k = options.pwr_k;
rate = options.rate;
lwbd = options.lwbd;
matchingAlgorithm = options.matchingAlgorithm;
certificateSpec = options.certificateSpec;
lambda_num = options.lambda_num;
lambda_den = options.lambda_den;

clusterCount = sort(clusterCount);

edgesCut = 0;
cuttimes = 0;
inittimes = 0;
spectimes = 0;
flowtimes = 0;

[G, weight] = loadMetisGraph(fileToRead);
n = size(G, 1);
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
    
    [expansionFound, edgesCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, lower] = ...
                    cutfind(flowgraph, outputFile=1, suffix=suffix, t=t, stop=stop, eta=eta, init=init, seed=seed, ...
                    p=p, pwr_k=pwr_k, rate=rate, lwbd=lwbd, matchingAlgorithm=matchingAlgorithm, certificateSpec=certificateSpec, ...
                    lambda_num=lambda_num, lambda_den=lambda_den);
    
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
    if lambda_num > 0
        clusterExpansions(end, 1) = clusterExpansions(end, 1) + lambda_num * sum(Cmask | sharedBefore) / lambda_den;
        clusterExpansions(largestClusterIndex, 1) = clusterExpansions(largestClusterIndex, 1) + lambda_num * sum(Cmask | sharedBefore) / lambda_den;
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

    