function [  ] = RunDirectoryOnce(graphDirectory, resultDirectory, lamdas)
%RUNDIRECTORYGRAPHS Runs all graphs in the input directory and saves
%   results in the output directory
%
%   Finds all .eg2 files in the input directory and saves all results in
%   the outputDirectory.
%
% INPUTS:
%   (char) graphDirectory - Directory to search for .eg2 graph files
%   (char) resultDirectory - Directory to read algorithm results and write the new ones
%   (int) lamdas - list of lamdas for which to run and save results
%
% OUTPUTS:
%
% DESCRIPTION:
%   Orchestration function designed to take in non-overlapping partitions
%   and produce overlapping ones using different lambdas by running FI just once
%   Preferrably call from `screen` on a terminal or use some kind of cloud (eg SLURM).

algs = {'KernighanLin', 'SweepCut', 'metis'};
% algs = {'metis'};
directorySearch = fullfile(graphDirectory, '*.eg2');
files = dir(directorySearch);
for f=1:length(files)
    fprintf('%s\n', files(f).name);
    [~, dataset, ~] = fileparts(files(f).name);
    inputFilename = fullfile(graphDirectory, files(f).name);
    [G, ~, ~] = loadeg2graph(inputFilename);
    for alg=algs
        inputFilename = fullfile(resultDirectory, sprintf('%s_%s.ptn', dataset, alg{1}));
        fprintf('%s\n', inputFilename);
        if ~isfile(inputFilename)
            continue;
        end
        partitions = readPtn(inputFilename);
        for lamda=lamdas
                outputFilename = fullfile(resultDirectory, sprintf('%s_%s_runOnce_%.0f.ptn', dataset, alg{1}, 100 / lamda));
            if isfile(outputFilename)
                continue;
            end
            try
                [~, ~] = RunOnce(G, partitions, outputFilename, lamda);
            catch
                fprintf(2, 'Failed %s for lambda = %f', dataset, lamda)
                continue;
            end
        end
    end
end

