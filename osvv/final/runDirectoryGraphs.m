function [ ] = runDirectoryGraphs(inputDirectory, outputDirectory, lamdas_num, lamdas_den, balances)
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
%     if ~any(strcmp(dataset, ["auto" "wave"]))
%         continue;
%     end
    inputFilename = fullfile(inputDirectory, files(f).name);
    [G, n, m] = loadeg2graph(inputFilename);
    if ~exist(fullfile(outputDirectory, sprintf('%s.mat', dataset)), 'file') 
        [vec, ~] = eigs(diag(sum(G)) - G, 3, 'SA');
        save(fullfile(outputDirectory, sprintf('%s.mat', dataset)), 'vec');
    end
    for balance=balances
        for l=1:length(lamdas_num)
            lamda_num = lamdas_num(l);
            lamda_den = lamdas_den(l);
            ptnFilename = fullfile(outputDirectory, sprintf('%s_balanced_%d_%d_%d.ptn', dataset, lamda_num, lamda_den, balance));
            if exist(ptnFilename)
                continue;
            end
            try
                [expansionFound, edgesCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, lower] = cutfind(G, 1, '', 1000, 4, 5, 1, 42, 10, 'infty', 'n', 1, balance/1000, lamda_num, lamda_den);
            catch
                fprintf(2, 'Failed lambda=%.2f\n', lamda);
            end
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
        end
    end
end

end

