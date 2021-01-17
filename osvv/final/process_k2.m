function [dataset, algorithm, lambda, score, overlap] = process_k2(graphDirectory, inputDirectory, version, embedding)
%PROCESS_K2 Process results for k=2 communities
%   
% INPUTS:
%   (char) graphDirectory - Directory that holds the graphs
%   (char) inputDirectory - Directory where the graph partitions are saved
%   (char) version - Which results you want processed
%   (char) embedding - Which embedding to use for plotting
%
% OUTPUTS:
%   (char) dataset - Names for each dataset in results
%   (char) algorithm - Which one was used to produce the 
%   (double) lambda - lambda used to produce each result
%   (double) score - for each dataset saves the sparsest cut
%   (double) overlap - Weight in the overlap for each partition


%% Parameter checking
error_string = 'Error in parameter %s. See README file for usage.\n';

if ~ischar(graphDirectory)
    error(error_string, 'graph directory');
end

if ~ischar(inputDirectory)
    error(error_string, 'input directory');
end

if ~ischar(embedding)
    error(error_string, 'embedding');
end


%% Read in data and compute results

algos = {'Sweep Cut', 'Kernighan-Lin', 'METIS'};
partitionExtensions = {'SweepCut', 'KernighanLin', 'metis'};

i = 1;
graphFiles = dir(fullfile(graphDirectory, '*.eg2'));
for f=1:length(graphFiles)
    % Load graph;
    graphFilename = graphFiles(f).name;
    [~, datasetName, ~] = fileparts(graphFilename);
    fprintf('Processing %s.', datasetName);
    [G, n, m] = loadeg2graph(fullfile(graphDirectory, graphFilename));
    weight = full(sum(G));
    volume = sum(weight);
    
    % Load embedding
    if embedding == 'spectral'
        vectorFilename = fullfile(inputDirectory, sprintf('%s.mat', datasetName));
        try
            vectors = load(vectorFilename).vec;
        catch
            fprintf(2, 'Unable to open vector file %s. Skipping dataset %s.\n', vectorFilename, datasetName);
            continue;
        end 
    end
    
    resultFiles = dir(fullfile(inputDirectory, sprintf('%s_%s_*.ptn', datasetName, version)));
    for f_res=1:length(resultFiles)
        resultFilename = fullfile(inputDirectory, resultFiles(f_res).name);
        [~, filename, ~] = fileparts(resultFilename);
        filenamePieces = split(filename, '_');
        partitions = readPtn(resultFilename);

        dataset{i} = datasetName;
        algorithm{i} = 'ORC-SDP';
        lambda{i} = 100 / str2double(filenamePieces{3});
        [~, ~, score{i}] = cutexp(G, 0, int64(weight), partitions{1}, partitions{2});
        overlappingNodes = intersect(partitions{1}, partitions{2});
        overlap{i} = 100 * sum(weight(overlappingNodes)) / volume;
        
        i = i + 1;
    end
    
    for j=1:length(algos)
        partitionFilename = fullfile(inputDirectory, ...
            sprintf('%s_%s_runOnce_*.ptn', datasetName, partitionExtensions{j}));
        competingResultFiles = dir(partitionFilename);
        for f_cmp=1:length(competingResultFiles)
            resultFilename = fullfile(resultFiles(f_cmp).folder, resultFiles(f_cmp).name);
            [~, filename, ~] = fileparts(resultFilename);
            filenamePieces = split(filename, '_');        
            partitions = readPtn(resultFilename);
            
            dataset{i} = datasetName;
            algorithm{i} = filenamePieces{2};
            lambda{i} = 100 / str2double(filenamePieces{4});
            [~, ~, score{i}] = cutexp(G, 0, int64(weight), partitions{1}, partitions{2});
            overlappingNodes = intersect(partitions{1}, partitions{2});
            overlap{i} = 100 * sum(weight(overlappingNodes)) / volume;

            i = i + 1;
        end
    end
       
    fprintf(' Total results read %d.\n', i-1);
end


%% Plot results
algs = {'ORC-SDP', algos{:}};
for ds=unique(dataset)
    figure;
    hold on;
    
    algMask = false(1, length(algs));
    k = 1;
    for alg=algs
        mask = strcmp(dataset, ds{1}) & strcmp(algorithm, alg);
        if sum(mask) == 0
            k = k + 1;
            continue;
        end
        algMask(k) = 1;
        sc = [score{mask}];
        ol = [overlap{mask}];
        la = [lambda{mask}];
        [~, id] = sort(la, 'descend');
        plot(ol(id), sc(id));
    end
    
    yl = ylim();
    xlabel('Overlap (%)');
    ylabel('Conductance');
    t = title(sprintf('Conductance vs Overlap (%%) for %s', ds{1}));
    %t.Position(2) = t.Position(2) + 0.05 * (yl(2) - yl(1));
    legend({algs{algMask}});
    exportgraphics(gcf, fullfile(inputDirectory, sprintf('%s_lambda.png', ds{1})));
    hold off;

end
end