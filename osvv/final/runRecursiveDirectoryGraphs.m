function [ ] = runRecursiveDirectoryGraphs(inputDirectory, outputDirectory, lamdas_num, lamdas_den, balances, clusterCount)
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

directorySearch = fullfile(inputDirectory, '*.metis');
files = dir(directorySearch);
for f=1:length(files)
    fprintf('%s\n', files(f).name);
    [~, dataset, ~] = fileparts(files(f).name);
    if ~any(strcmp(dataset, ['dcExtractedDblp']))
        continue;
    end
    inputFilename = fullfile(inputDirectory, files(f).name);
    [G, weights] = loadMetisGraph(inputFilename);
    if ~exist(fullfile(outputDirectory, sprintf('%s.mat', dataset)), 'file') 
        [vec, ~] = eigs(diag(sum(G)) - G, 3, 'SA');
        save(fullfile(outputDirectory, sprintf('%s.mat', dataset)), 'vec');
    end
    for balance=balances
        for l=1:length(lamdas_num)
            lamda_num = lamdas_num(l);
            lamda_den = lamdas_den(l);
            ptnFilename = fullfile(outputDirectory, sprintf('%s_recursive_%d_%d_%d_%d.ptn', dataset, lamda_num, lamda_den, balance, clusterCount));
            if exist(ptnFilename)
                continue;
            end
            %try
                [score, clusters] = recursiveCutfind(clusterCount, G, 1, '', 1000, 5, 10, 1, 42, 1000, 4, 'KL', 'n', 1, balance/1000, lamda_num, lamda_den);
            %catch
            %    fprintf(2, 'Failed lambda=%d / %d\n', lamda_num, lambda_den);
            %end
            for i=1:length(clusters)
                clusters{i} = clusters{i}';
            end
            ptnFile = fopen(ptnFilename, 'w');
            toPtn(ptnFile, clusters);
            fclose(ptnFile);
        end
    end
end

end



