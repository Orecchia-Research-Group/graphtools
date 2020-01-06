function [  ] = runDirectoryGraphs(inputDirectory, outputDirectory, clusterCounts, lamdas)
%RUNDIRECTORYGRAPHS Runs all graphs in the input directory and saves
%   results in the output directory
%
%   Finds all .eg2 files in the input directory and saves all results in
%   the outputDirectory.
%
% INPUTS:
%   (char) inputDirectory - Directory to search for .eg2 graph files
%   (char) outputDirectory - Directory to write results
%   (int) clusterCounts - number of cluster for which to save results
%   (int) lamdas - list of lamdas for which to run and save results
%
% OUTPUTS:
%
% DESCRIPTION:
%   Orchestration function designed to be called preferrably from `screen`
%   on a terminal for many small experiments.

directorySearch = fullfile(inputDirectory, '*.eg2');
files = dir(directorySearch);
for f=1:length(files)
    fprintf('%s\n', files(f).name);
    [~, dataset, ~] = fileparts(files(f).name);
    inputFilename = fullfile(inputDirectory, files(f).name);
    for lamda=lamdas
        [edgesCut, cutsFound, endtime, cuttimes, inittimes, spectimes, flowtimes] = iterative_cutfind(clusterCounts, inputFilename, 1, '', 1000, 10, 5, 1, 1, 10, 'infty', 'n', 1, lamda);
        fprintf('% 15s end: %9.2f. cut: %9.2f. init: %9.2f. spec: %9.2f. flow: %9.2f\n', dataset, endtime, cuttimes, inittimes, spectimes, flowtimes);
        for i=1:length(clusterCounts)
            if lamda > 0
                outputFilename = fullfile(outputDirectory, sprintf('%s.%d.%d.ptn', dataset, lamda, clusterCounts(i)));
            else
                outputFilename = fullfile(outputDirectory, sprintf('%s.%d.ptn', dataset, clusterCounts(i)));
            end
            toPtn(outputFilename, cutsFound{i});
        end
    end
end

end

