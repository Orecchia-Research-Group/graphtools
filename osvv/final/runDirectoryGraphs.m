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
%     if ~any(strcmp(dataset, ["auto" "wave"]))
%         continue;
%     end
    %if exist(fullfile(outputDirectory, sprintf('%s.gif', dataset))) == 2
    %    continue;
    %end
    inputFilename = fullfile(inputDirectory, files(f).name);
    [G, n, m] = loadeg2graph(inputFilename);
    [vec, ~] = eigs(diag(sum(G)) - G, 3, 'SA');
    save(fullfile(outputDirectory, sprintf('%s.mat', dataset)), 'vec');
    for lamda=lamdas
        plotFile = fullfile(outputDirectory, sprintf('%s.%02d.png', dataset, lamda));
        if exist(plotFile) == 2
            continue;
        end
        [expansionFound, edgesCut, L, R, H, endtime, inittime, spectime, flowtime, iterations, lower] = cutfind(G, 1, '', 1000, 4, 5, 1, 42, 10, 'infty', 'n', 1, lamda);
        fprintf('% 15s end: %9.2f. init: %9.2f. spec: %9.2f. flow: %9.2f\n', dataset, endtime, inittime, spectime, flowtime);
        Lmask = sparse(double(L), 1, true, n, 1);
        Rmask = sparse(double(R), 1, true, n, 1);
        Cmask = Lmask & Rmask;
%         colorIndex = ones(n, 1);
%         colorIndex(Lmask & ~Cmask) = 2;
%         colorIndex(Rmask & ~Cmask) = 3;
%         colors = [0 1 0; 1 0 0; 0 0 1];
%         c = colors(colorIndex, 1:end);
%         figure;
%         hold on;
%         scatter(vec(1:end, 2), vec(1:end, 3), 10, c, 'filled');
%         axis off;
%         title(sprintf('Graph: %s. Lambda: %2d. Overlap: %3d nodes.', dataset, lamda, full(sum(Cmask))));
        
%         cutnodes = full(sum(Cmask));
%         cutnodes_perc = full(cutnodes / n);
%         cutedges = edgesCut - lamda * cutnodes;
%         cutedges_perc = full(cutedges / m);
        
        ptnFilename = fullfile(outputDirectory, sprintf('%s_%02d.ptn', dataset, lamda));
        ptnFile = fopen(ptnFilename, 'w');
%         fprintf(ptnFile, 'cutedges %d,cutedges_perc %f,expansion %f,overlap %d,overlap_perc %f\n', cutedges, cutedges_perc, expansionFound, cutnodes, cutnodes_perc);
        partitions{1} = L';
        partitions{2} = R';
        toPtn(ptnFile, partitions);
        fclose(ptnFile);
%         
%         str = {
%             sprintf('Cutedges: %d', cutedges),
%             sprintf('Overlap: %d', cutnodes),
%             sprintf('Expansion: %.3f', expansionFound)
%         };
%         text(-0.1, 0.85, str, 'Units', 'normalized');
%         saveas(gcf, plotFile);
        if sum(Cmask) == 0
            break;
        end
        
    end
end

end

