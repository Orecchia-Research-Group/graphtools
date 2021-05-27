function [L, R] = RunOnce(graph, partition, outfile, lamda_num, lamda_den)
%RUNONCE Run cutfind once on a graph and its partition
%
% INPUT:
%   (char) graph - a eg2 undirected graph file to read - must be a valid graph
%   (char) partition - a starting partition to find an overlapping partition
%   (char) outfile - .ptn file to write the produced partition
%   (int) lamda_num, lamda_den - fraction of the internal edge to use
%
% OUTPUT:
%   (double) L - One side of the overlapping partitions
%   (double) R - Other side of the overlapping partitions

error_string = 'Error in parameter %s. See README file for usage.\n';

if (~ischar(graph) && (~isnumeric(graph) || (size(graph, 1) ~= size(graph, 2))))
    error(error_string, 'eg2 input file');
end

if (~ischar(partition) && (~isnumeric(partition{1}) || ~isnumeric(partition{2})))
    error(error_string, 'partition input file');
end

if (ischar(graph))
    [G, n, ~] = loadeg2graph(graph);
else
    G = graph;
    n = size(G, 1);
end

lamda_num = int64(lamda_num);
lamda_den = int64(lamda_den);

if (ischar(partition))
    partitions = readPtn(partition);
    cut = partitions{1};
    reciprocalCut = partitions{2};
else
    cut = partition{1};
    reciprocalCut = partition{2};
end

minweirdrat_num = double(10);
minweirdrat_den = int64(1);
minweirdrat = double(10);
p = int64(5000);
nomatching = 0;

m = nnz(G)/2;
degree = int64(full(sum(G)));
weight = ones(1, n, 'int64');
weight(:) = degree;

[weirdrat_num, weirdrat_den, weirdrat] =  cutweird(G, cut, reciprocalCut, cut, weight, lamda_num, lamda_den); % COMPUTE NEW WEIRDRAT
[oldex_num, oldex_den, oldex] = cutexp(G, lamda_num, lamda_den, weight, cut, reciprocalCut);
[minweirdrat_num, minweirdrat_den, minweirdrat, ex_num, ex_den, ex, cut, reciprocalCut, matching, matchrat, iterflownumber] =  RunFlow(G, cut, weight, minweirdrat_num, minweirdrat_den, minweirdrat, p, nomatching, 0, lamda_num, lamda_den);
fprintf('Name: %s. lamda: %d / %d = %.3f. Old score: %d / %d = %f. New score: %d / %d = %f\n', outfile, lamda_num, lamda_den, double(lamda_num) / double(lamda_den), oldex_num, oldex_den, oldex, ex_num, ex_den, ex);
partitions{1} = cut';
partitions{2} = reciprocalCut';
ptnFile = fopen(outfile, 'w');
toPtn(ptnFile, partitions);
fclose(ptnFile);
L = cut;
R = reciprocalCut;
end

