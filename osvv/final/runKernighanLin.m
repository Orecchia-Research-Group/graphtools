function [ ] = runKernighanLin( inputDirectory, outputDirectory)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

directorySearch = fullfile(inputDirectory, '*.eg2');
files = dir(directorySearch);
for f=1:length(files)
    fprintf('%s\n', files(f).name);
    [~, dataset, ~] = fileparts(files(f).name);
    if any(strcmp(dataset, ["bcsstk29" "bcsstk31" "fe_body" "fe_pwt"]))
        continue;                                                     
    end
    inputFilename = fullfile(inputDirectory, files(f).name);
    [~, n, ~] = loadeg2graph(inputFilename);
    [edgeCut, overlappingNodes, L, R] = KernighanLin(inputFilename, n);
    ptnFilename = fullfile(outputDirectory, sprintf('%s_KernighanLin.ptn', dataset));
    partitions{1} = L';
    partitions{2} = R';
    toPtn(ptnFilename, partitions);
    overlapNodesFilename = fullfile(outputDirectory, sprintf('%s_KernighanLin_nodes.txt', dataset));
    f_overlap = fopen(overlapNodesFilename, 'w');
    for node=overlappingNodes
        fprintf(f_overlap, '%d\n', node);
    end
    
end

