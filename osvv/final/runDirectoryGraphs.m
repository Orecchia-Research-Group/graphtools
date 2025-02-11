function [ ] = runDirectoryGraphs(inputDirectory, outputDirectory, lambdas_num, lambdas_den, balances, pwr_k, eta)
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
    if ~any(strcmp(dataset, ['dcExtractedDblp']))       % TODO: Remove working only for youtube.
       continue;
    end
    inputFilename = fullfile(inputDirectory, files(f).name);
    [G, weights] = loadMetisGraph(inputFilename);
    n = size(G, 1);
    %[G, n, m] = loadeg2graph(inputFilename);
    %if ~exist(fullfile(outputDirectory, sprintf('%s.mat', dataset)), 'file') 
    %    [vec, ~] = eigs(diag(sum(G)) - G, 3, 'SA');
    %    save(fullfile(outputDirectory, sprintf('%s.mat', dataset)), 'vec');
    %end
    for balance=balances
        for l=1:length(lambdas_num)
            lambda_num = lambdas_num(l);
            lambda_den = lambdas_den(l);
            ptnFilename = fullfile(outputDirectory, sprintf('%s_parallel_%d_%d_%d_%d_%d.ptn', dataset, lambda_num, lambda_den, balance, pwr_k, eta));
            %if exist(ptnFilename)
            %    continue;
            %end
            %try
                [expansionFound, edgesCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, lower] = ...
                    cutfind(inputFilename, outputFile=1, suffix='', t=100, stop=1:5, eta=eta, init=1, seed=0, ...
                    p=100000, pwr_k=pwr_k, rate='KL', lwbd='y', matchingAlgorithm='dinic', certificateSpec=1, ...
                    ufactor=balance/1000, lambda_num=lambda_num, lambda_den=lambda_den);
            %catch
            %    fprintf(2, 'Failed lambda= %d / %d = %.2f\n', lambda_num , lambda_den, lambda_num / lambda_den);
            %end
            fprintf('% 15s end: %9.2f. init: %9.2f. spec: %9.2f. flow: %9.2f\n', dataset, endtime, inittime, spectime, flowtime);
            Lmask = sparse(double(L), 1, true, n, 1);
            Rmask = sparse(double(R), 1, true, n, 1);
            Cmask = Lmask & Rmask;
            ptnFile = fopen(ptnFilename, 'w');
            partitions{1} = L';
            partitions{2} = R';
            toPtn(ptnFile, partitions);
            fclose(ptnFile);
            if sum(sum(G(Lmask & ~Cmask, Rmask & ~Cmask))) == 0
                break;
            end
            lower_filename = fullfile(outputDirectory, sprintf('%s_parallel_%d_%d_%d_%d_%d.lower', dataset, lambda_num, lambda_den, balance, pwr_k, eta));
            lower_file = fopen(lower_filename, 'w');
            for i=1:size(lower, 1)
                fprintf(lower_file, '%d %d %f %f %f\n', lower(i, 1), lower(i, 2), lower(i, 3), lower(i, 4), lower(i, 5));
            end
        end
    end
end

end

