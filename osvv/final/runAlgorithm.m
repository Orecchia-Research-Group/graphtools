function [ ] = runAlgorithm( inputDirectory, outputDirectory, func, name, eta, ufactor)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

directorySearch = fullfile(inputDirectory, '*.eg2');
files = dir(directorySearch);
for f=1:length(files)
    fprintf('%s\n', files(f).name);
    [~, dataset, ~] = fileparts(files(f).name);
    if ~any(strcmp(dataset, {'amazon', 'dcExtractedDblp', 'youtube'}))
        continue;                                                     
    end
    inputFilename = fullfile(inputDirectory, files(f).name);
    [~, n, ~] = loadeg2graph(inputFilename);
    [edgeCut, L, R] = func(inputFilename, eta, ufactor);
    ptnFilename = fullfile(outputDirectory, sprintf('%s_%s.ptn', dataset, name));
    partitions{1} = L';
    partitions{2} = R';
    toPtn(ptnFilename, partitions);
    
end

