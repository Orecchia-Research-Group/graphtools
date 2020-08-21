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
    if any(strcmp(dataset, ["bcsstk29" "bcsstk31" "fe_body" "fe_pwt"]))
        continue;                                                     
    end
%     if ~(strcmp(dataset, "minedDBLP"))
%         continue;
%     end
%     if ~any(strcmp(dataset, ["auto" "wave"]))
%         continue;
%     end
    %if exist(fullfile(outputDirectory, sprintf('%s.gif', dataset))) == 2
    %    continue;
    %end
    inputFilename = fullfile(inputDirectory, files(f).name);
    [G, n, m] = loadeg2graph(inputFilename);
    if ~exist(fullfile(outputDirectory, sprintf('%s.mat', dataset)), 'file') 
        [vec, ~] = eigs(diag(sum(G)) - G, 3, 'SA');
        save(fullfile(outputDirectory, sprintf('%s.mat', dataset)), 'vec');
    end
    for lamda=lamdas
%         plotFile = fullfile(outputDirectory, sprintf('%s.%02d.png', dataset, lamda));
%         if exist(plotFile, 'file') == 2
%             continue;
%         end
        try
            [expansionFound, edgesCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, lower] = cutfind(G, 1, '', 1000, 4, 5, 1, 42, 10, 'infty', 'n', 1, lamda);
            fprintf('% 15s end: %9.2f. init: %9.2f. spec: %9.2f. flow: %9.2f\n', dataset, endtime, inittime, spectime, flowtime);
            Lmask = sparse(double(L), 1, true, n, 1);
            Rmask = sparse(double(R), 1, true, n, 1);
            Cmask = Lmask & Rmask;
            HFilename = fullfile(outputDirectory, sprintf('%s_final_H_%.0f.mat', dataset, 100 / lamda));
            save(HFilename, 'H');
            ptnFilename = fullfile(outputDirectory, sprintf('%s_final_%.0f.ptn', dataset, 100 / lamda));
            ptnFile = fopen(ptnFilename, 'w');
            partitions{1} = L';
            partitions{2} = R';
            toPtn(ptnFile, partitions);
            fclose(ptnFile);
            if sum(sum(G(Lmask & ~Cmask, Rmask & ~Cmask))) == 0
                break;
            end
        catch
            fprintf(2, 'Failed lambda=%.2f\n', lamda);
        end
    end
end

end

