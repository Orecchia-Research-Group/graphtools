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



function [score, clusters] = recursiveCutfind(clusterCount, fileToRead, options)

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
    options.ufactor (1, 1) double {mustBeLessThanOrEqual(options.ufactor, 0.5)} = 0
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
ufactor = options.ufactor;
lambda_num = options.lambda_num;
lambda_den = options.lambda_den;

if(ischar(fileToRead))
    [G, weight] = loadMetisGraph(fileToRead);
    n = size(G, 1);
    weight = int64(weight);
else
    G = fileToRead;
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
    
    [expansionFound, edgesCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, lower] = ...
                    cutfind(flowgraph, outputFile=1, suffix=suffix, t=t, stop=stop, eta=eta, init=init, seed=seed, ...
                    p=p, pwr_k=pwr_k, rate=rate, lwbd=lwbd, matchingAlgorithm=matchingAlgorithm, certificateSpec=certificateSpec, ...
                    ufactor=currentUfactor, lambda_num=lambda_num, lambda_den=lambda_den);
    
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
        Rmask(clusters{cr}) = true;
    end
    R = find(Rmask);
    [~, ~, clusterExpansion(c)] = cutexp(G, int64(lambda_num), int64(lambda_den), int64(weight), int64(L), int64(R));
end
score = min(clusterExpansion);

end

