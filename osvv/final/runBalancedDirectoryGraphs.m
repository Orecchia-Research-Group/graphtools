function [  ] = runBalancedDirectoryGraphs(inputDirectory, outputDirectory, clusterCounts, lamdas, bal)
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
    if any(strcmp(dataset, ["144" "bcsstk29" "bcsstk31" "fe_body" "fe_pwt"]))
        continue;                                                     
    end
    %if exist(fullfile(outputDirectory, sprintf('%s.gif', dataset))) == 2
    %    continue;
    %end
    inputFilename = fullfile(inputDirectory, files(f).name);
    [G, n, m] = loadeg2graph(inputFilename);
    %[vec, ~] = eigs(diag(sum(G)) - G, 3, 'SA');
    %save(fullfile(outputDirectory, sprintf('%s.mat', dataset)), 'vec');
    for lamda=lamdas
        [edgesCut, L, R] = BalancedCut(bal, G, 1, '', 1000, 4, 5, 1, 42, 10, 'infty', 'n', 1, lamda);
        Lmask = sparse(double(L), 1, true, n, 1);
        Rmask = sparse(double(R), 1, true, n, 1);
        Cmask = Lmask & Rmask;
        
        ptnFilename = fullfile(outputDirectory, sprintf('bal_%s_%02d.ptn', dataset, lamda));
        ptnFile = fopen(ptnFilename, 'w');
        partitions{1} = L';
        partitions{2} = R';
        toPtn(ptnFile, partitions);
        fclose(ptnFile);
        if sum(Cmask) == 0
            break;
        end
        
    end
end

end

