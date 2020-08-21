function [ ] = runAlgorithm( inputDirectory, outputDirectory, func, name)
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
    [edgeCut, L, R] = func(inputFilename, 5);
    ptnFilename = fullfile(outputDirectory, sprintf('%s_%s.ptn', dataset, name));
    partitions{1} = L';
    partitions{2} = R';
    toPtn(ptnFilename, partitions);
    
end

