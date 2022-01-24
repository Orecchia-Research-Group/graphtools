function [] = processSocial(graphDirectory, inputDirectory)
%ProcessGreedySocial
%
% INPUTS:
%   (char) graphDirectory - Directory that holds the graphs
%   (char) inputDirectory - Directory where the graph partitions are saved
%
% OUTPUTS:

%% Parameter checking
error_string = 'Error in parameter %s. See README file for usage.\n';

if ~ischar(graphDirectory)
    error(error_string, 'graph directory');
end

if ~ischar(inputDirectory)
    error(error_string, 'input directory');
end

algorithms = {'${\tt cm + improve}$', '${\tt METIS} + {\tt OverlapImprove}$', '${\tt METIS} + {\tt GreedyImprove}$', '${\tt Spectral} + {\tt GreedyImprove}$'};

dataset = {};
algorithm = {};
C = {};
p = {};
Lvol = [];
Rvol = [];
Cvol = [];
edgeConductance = [];
mixedConductance = [];
lambda_num = [];
lambda_den = [];
lambda = [];
balances = [];

i = 1;
graphFiles = dir(fullfile(graphDirectory, '*.metis'));
for f=1:length(graphFiles)
    % Load graph;
    graphFilename = graphFiles(f).name;
    [~, datasetName, ~] = fileparts(graphFilename);
    datasetPieces = split(graphFilename, '_');
        
    % ORC
    fprintf('Processing %s.\n', datasetName);
    if any(strcmp(datasetName, {'orkut', 'livejournal'}))
       continue
    end
    [G, weight] = loadMetisGraph(fullfile(graphDirectory, graphFilename));
    volume = sum(weight);
    
    try
        temp = importdata(fullfile(graphDirectory, sprintf('%s_authorAreas.txt', datasetName)));
        category_names = split(temp.textdata{1});
        category_names = upper(category_names([2 4 6 3 5 7]));
        categories = temp.data(:, [1 3 5 2 4 6]);
    catch
    end
    resultFiles = dir(fullfile(inputDirectory, sprintf('%s_parallel_*_*_*.ptn', datasetName)));
    for f_res=1:length(resultFiles)
        resultFilename = fullfile(inputDirectory, resultFiles(f_res).name);
        [~, filename, ~] = fileparts(resultFilename);
        filenamePieces = split(filename, '_');
        partitions = readPtn(resultFilename);
        if partitions{1}(1) == 1
            L{i} = partitions{1};
            R{i} = partitions{2};
        else
            L{i} = partitions{2};
            R{i} = partitions{1};
        end
        inter = intersect(partitions{1}, partitions{2});
        Cvol(i) = sum(weight(inter));
        Lvol(i) = sum(weight(partitions{1}));
        Rvol(i) = sum(weight(partitions{2}));
        if strcmp(datasetName, 'dcExtractedDblp')
            for p=1:2
                category_counts(i, p, :) = sum(categories(partitions{p}, :)) - sum(categories(inter, :));
            end
            category_counts(i, 3, :) = sum(categories(inter, :));
        end
        dataset{i} = datasetName;
        algorithm{i} = '${\tt cm + improve}$';
        lambda_num(i) = str2double(filenamePieces{3});
        lambda_den(i) = str2double(filenamePieces{4});
        lambda(i) = lambda_num(i) / lambda_den(i);
        balances(i) = str2double(filenamePieces{5});
        [~, ~, edgeConductance(i)] = cutexp(G, int64(-1), int64(1), int64(weight), partitions{1}, partitions{2});
        [~, ~, mixedConductance(i)] = cutexp(G, int64(lambda_num(i)), int64(lambda_den(i)), int64(weight), partitions{1}, partitions{2});
        i = i + 1;
    end
    cm_count = i - 1;
    fprintf('Read %d results for cm + improve\n', cm_count);
    
    % METIS + OverlapImprove
    resultFiles = dir(fullfile(inputDirectory, sprintf('%s_metis_runOnce_*_*_*.ptn', datasetName)));
    for f_res=1:length(resultFiles)
        resultFilename = fullfile(inputDirectory, resultFiles(f_res).name);
        [~, filename, ~] = fileparts(resultFilename);
        filenamePieces = split(filename, '_');
        partitions = readPtn(resultFilename);
        if partitions{1}(1) == 1
            L{i} = partitions{1};
            R{i} = partitions{2};
        else
            L{i} = partitions{2};
            R{i} = partitions{1};
        end
        inter = intersect(partitions{1}, partitions{2});
        Cvol(i) = sum(weight(inter));
        Lvol(i) = sum(weight(partitions{1}));
        Rvol(i) = sum(weight(partitions{2}));
        if strcmp(datasetName, 'dcExtractedDblp')
            for p=1:2
                category_counts(i, p, :) = sum(categories(partitions{p}, :)) - sum(categories(inter, :));
            end
            category_counts(i, 3, :) = sum(categories(inter, :));
        end
        dataset{i} = datasetName;
        algorithm{i} = '${\tt METIS} + {\tt OverlapImprove}$';
        lambda_num(i) = str2double(filenamePieces{4});
        lambda_den(i) = str2double(filenamePieces{5});
        lambda(i) = lambda_num(i) / lambda_den(i);
        balances(i) = str2double(filenamePieces{6});
        [~, ~, edgeConductance(i)] = cutexp(G, int64(-1), int64(1), int64(weight), partitions{1}, partitions{2});
        [~, ~, mixedConductance(i)] = cutexp(G, int64(lambda_num(i)), int64(lambda_den(i)), int64(weight), partitions{1}, partitions{2});
        i = i + 1;
    end
    runOnce_count = i - 1 - cm_count;
    fprintf('Read %d results for metis + runOnce\n', runOnce_count);
    
    % METIS + greedyImprove
    resultFiles = dir(fullfile(inputDirectory, sprintf('%s_metis_*.nodes', datasetName)));
    LvolIndex = 2;
    RvolIndex = 3;
    CvolIndex = 4;
    precissionIndex = 5;
    recallIndex = 6;
    f1scoreIndex = 7;
    edgeConductanceIndex = 8;
    mixedConductanceIndex = 9;
    lambdaIndex = 10;
    for f_res=1:length(resultFiles)
        resultFilename = fullfile(inputDirectory, resultFiles(f_res).name);
        [~, filename, ~] = fileparts(resultFilename);
        filenamePieces = split(filename, '_');
        data = importdata(resultFilename);
        currentLength = length(data);
        
        for j=i:i+currentLength-1
            dataset{j} = datasetName;
            algorithm{j} = '${\tt METIS} + {\tt GreedyImprove}$';
            balances(j) = 500 - str2double(filenamePieces{3}) / 2;
        end
        currentIndex = i:i+currentLength-1;
        edgeConductance(currentIndex) = data(:, edgeConductanceIndex);
        mixedConductance(currentIndex) = data(:, mixedConductanceIndex);
        lambda(currentIndex) = data(:, lambdaIndex);
        i = i + currentLength;
    end
    greedy_count = i - 1 - runOnce_count;
    fprintf('Read %d results for metis + greedy\n', greedy_count);
    
    % SweepCut + greedyImprove
    resultFiles = dir(fullfile(inputDirectory, sprintf('%s_SweepCut_450.nodes', datasetName)));
    LvolIndex = 2;
    RvolIndex = 3;
    CvolIndex = 4;
    precissionIndex = 5;
    recallIndex = 6;
    f1scoreIndex = 7;
    edgeConductanceIndex = 8;
    mixedConductanceIndex = 9;
    lambdaIndex = 10;
    for f_res=1:length(resultFiles)
        resultFilename = fullfile(inputDirectory, resultFiles(f_res).name);
        [~, filename, ~] = fileparts(resultFilename);
        filenamePieces = split(filename, '_');
        data = importdata(resultFilename);
        currentLength = length(data);
        
        for j=i:i+currentLength-1
            dataset{j} = datasetName;
            algorithm{j} = '${\tt Spectral} + {\tt GreedyImprove}$';
            balances(j) = str2double(filenamePieces{3});
        end
        currentIndex = i:i+currentLength-1;
        edgeConductance(currentIndex) = data(:, edgeConductanceIndex);
        mixedConductance(currentIndex) = data(:, mixedConductanceIndex);
        lambda(currentIndex) = data(:, lambdaIndex);
        i = i + currentLength;
    end
    sweepCut_count = i - 1 - greedy_count;
    fprintf('Read %d results for metis + greedy\n', sweepCut_count);
end


%% Plot
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

colors = {'b', 'r', 'g', 'k'};
markers = {'o', 'o', 'none', 'none'};
resolution = 500;

unique_balances = sort(unique(balances));

datasets = unique(dataset);
for d=1:length(datasets)
    d_index = strcmp(dataset, datasets{d});
    for b=1:length(unique_balances)
        b_index = d_index & (unique_balances(b) == balances);
        if unique_balances(b) ~= 450
            continue;
        end
        figure;
        hold on;
        for a=1:length(algorithms)-1
            a_index = b_index & strcmp(algorithm, algorithms{a});
            if sum(a_index) == 0
                continue
            end
            currLambda = lambda(a_index);
            currEdgeConductance = edgeConductance(a_index);
            currMixedConductance = mixedConductance(a_index);
            [~, index] = sort(currLambda);
            
            plot(1 ./ currLambda(index), currMixedConductance(index), 'Color', colors{a}, 'Marker', markers{a}, 'MarkerFaceColor', colors{a}, 'MarkerSize', 3);
        end
        set(gca, 'XScale', 'log');
        % title({'Conductance relative to $1/\lambda$ in different ways to get an overlapping cut from', sprintf('SweepCut for synthetic dataset with %s communities and %s overlap', balance, o_title)});
        axis([0 inf 0 inf]);
        xlabel('$1/\lambda$');
        ylabel('$q_{G, \lambda}\left(\left[S, T\right]\right)$');
        legend(algorithms);
        exportgraphics(gcf, fullfile(inputDirectory, sprintf('%s_lambda_%d.png', datasets{d}, unique_balances(b))), 'Resolution', resolution);
    end
end


%% Radarplots
lambdas = unique(lambda);
for d=1:length(datasets)
    if ~strcmp(datasets{d}, 'dcExtractedDblp')
        continue
    end
    d_index = strcmp(dataset, datasets{d});
    for b=1:length(unique_balances)
        b_index = d_index & (unique_balances(b) == balances);
        if unique_balances(b) ~= 450
            continue;
        end
        for a=1:2
            a_index = b_index & strcmp(algorithm, algorithms{a});
            if sum(a_index) == 0
                continue;
            end
            for l=1:length(lambdas)
                l_index = a_index & (lambda == lambdas(l));
                if sum(l_index) == 0
                    continue;
                end
                figure;
                hold on;
                to_plot = squeeze(category_counts(l_index, :, :)) ./ sum(categories);
                sp = sum(to_plot(1, :));
                pr = to_plot(1, :) / sp;
                fprintf('Balance = %.3f Lambda = %.3f: %.3f %.3f %.3f %.3f %.3f %.3f\n', balances(l_index), lambda(l_index), pr(1), pr(2), pr(3), pr(4), pr(5), pr(6));
                hnd = radarplot(to_plot, category_names);
                legend(hnd, {'L', 'R', 'C'});
                % title({'Percentage of papers in each category that belong', sprintf('to each community when target balance = %d', tb(i))})
                saveas(gcf, fullfile(inputDirectory, sprintf('%s_spider_%d_%d_%d_%d.png', datasets{d}, a, lambda_num(l_index), lambda_den(l_index), balances(l_index))));
                close;
            end
        end
    end
end


end