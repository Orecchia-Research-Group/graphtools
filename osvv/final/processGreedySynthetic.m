function [] = processGreedySynthetic(graphDirectory, inputDirectory)
%ProcessGreedySynthetic
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

%% Read in data and compute results

balances = {'balanced', 'unbalanced'};  % {'balanced', 'unbalanced'}
overs = {'log', 'invSqrt'};                        % {'log', 'invSqrt'}
ps = {'normal', };                        % {'normal', 'sqrt'}
qs = {'sparse', 'normal', 'dense', 'constant'};  % {'sparse', 'normal', 'dense', 'constant'}
algs = {'SweepCut', 'KernighanLin'};

dataset = {};
algorithm = {};
overlapAlgorithm = {};
balanced = {};
C = {};
p = {};
q = {};
runIndex = [];
Lvol = [];
Rvol = [];
Cvol = [];
precission = [];
recall = [];
f1score = [];
edgeConductance = [];
mixedConductance = [];
lambda = [];

i = 0;
graphFiles = dir(fullfile(graphDirectory, '*.eg2'));
for f=1:length(graphFiles)
    % Load graph;
    graphFilename = graphFiles(f).name;
    [~, datasetName, ~] = fileparts(graphFilename);
    datasetPieces = split(graphFilename, '_');
    if ~any(strcmp(datasetPieces{4}, overs)) || ~any(strcmp(datasetPieces{5}, ps))
        continue;
    end
    
    fprintf('Processing %s.\n', datasetName);
    [G, n, m] = loadeg2graph(fullfile(graphDirectory, graphFilename));
    weight = full(sum(G));
    volume = sum(weight);
    
    vectorFilename = fullfile(inputDirectory, sprintf('%s.mat', datasetName));
    try
        vectors = load(vectorFilename).vec;
    catch
        fprintf(2, 'Unable to open vector file %s. Skipping dataset %s.\n', vectorFilename, datasetName);
        continue;
    end
    resultFiles = dir(fullfile(inputDirectory, sprintf('%s_*.nodes', datasetName)));
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
        
        for j=i+1:i+currentLength
            dataset{j} = datasetName;
            algorithm{j} = filenamePieces{8};
            overlapAlgorithm = 'greedy';
            balanced{j} = filenamePieces{3};
            C{j} = filenamePieces{4};
            p{j} = filenamePieces{5};
            q{j} = filenamePieces{6};
            runIndex(j) = str2double(filenamePieces{7});
        end
        currentIndex = i+1:i+currentLength;
        Lvol(currentIndex) = data(:, LvolIndex);
        Rvol(currentIndex) = data(:, RvolIndex);
        Cvol(currentIndex) = data(:, CvolIndex);
        precission(currentIndex) = data(:, precissionIndex);
        recall(currentIndex) = data(:, recallIndex);
        f1score(currentIndex) = data(:, f1scoreIndex);
        edgeConductance(currentIndex) = data(:, edgeConductanceIndex);
        mixedConductance(currentIndex) = data(:, mixedConductanceIndex);
        lambda(currentIndex) = data(:, lambdaIndex);
        i = i + currentLength;
    end
end

%% Aggregate results and plot
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

resolution = 500;

for p_ind = 1:1
    p_prob = ps{p_ind};
    p_index = strcmp(p, p_prob);
    for b=1:2
        balance = balances{b};
        b_index = p_index & strcmp(balanced, balance);
        for o=1:2
            if strcmp(overs{o}, 'log')
                o_title = '$\log(n)$';
            else
                o_title = '$\sqrt{n}$';
            end
            o_index = b_index & strcmp(C, overs{o});
            for a=2:2
                a_title = algs{a};
                a_index = o_index & strcmp(algorithm, algs{a});
                for j=1:4
                    legd{j} = sprintf('C edges are %s', qs{j});
                    j_index = a_index & strcmp(q, qs{j});
                    l_values = sort(unique(lambda(j_index)));
                    points = length(l_values);
                    currPrecission = zeros(5, points);
                    currRecall = zeros(5, points);
                    currF1score = zeros(5, points);
                    currEdgeConductance = zeros(5, points);
                    currMixedConductance = zeros(5, points);
                    for r=1:5
                        r_index = j_index & (runIndex == r);
                        noise = rand(1, sum(r_index))*1e-9;
                        pre = precission(r_index);
                        currPrecission(r, :) = interp1(lambda(r_index) + noise, pre, l_values, 'linear', pre(end));
                        rec = recall(r_index);
                        currRecall(r, :) = interp1(lambda(r_index) + noise, rec, l_values, 'linear', rec(end));
                        f1 = f1score(r_index);
                        currF1score(r, :) = interp1(lambda(r_index) + noise, f1, l_values, 'linear', f1(end));
                        ec = edgeConductance(r_index);
                        currEdgeConductance(r, :) = interp1(lambda(r_index) + noise, ec, l_values, 'linear', ec(end));
                        mc = mixedConductance(r_index);
                        currMixedConductance(r, :) = interp1(lambda(r_index) + noise, mc, l_values, 'linear', mc(end));
                    end
                    currPrecission(isnan(currPrecission)) = 0;
                    currRecall(isnan(currRecall)) = 0;
                    currF1score(isnan(currF1score)) = 0;
                    currEdgeConductance(isnan(currEdgeConductance)) = 0;
                    currMixedConductance(isnan(currMixedConductance)) = 0;
                    meanPrecission = mean(currPrecission);
                    meanRecall = mean(currRecall);
                    meanF1score = mean(currF1score);
                    meanEdgeConductance = mean(currEdgeConductance);
                    meanMixedConductance = mean(currMixedConductance);
                    figure;
                    hold on;

                    X = [1./l_values fliplr(1./l_values)];
                    pre_plot = plot(1./l_values, meanPrecission, 'g-');
                    Y = [min(currPrecission) fliplr(max(currPrecission))];
                    patch(X, Y, 'g', 'FaceAlpha', 0.2);
                    rec_plot = plot(1./l_values, meanRecall, 'r-');
                    Y = [min(currRecall) fliplr(max(currRecall))];
                    patch(X, Y, 'r', 'FaceAlpha', '0.2');
                    f1_plot = plot(1./l_values, meanF1score, 'y-');
                    Y = [min(currF1score) fliplr(max(currF1score))];
                    patch(X, Y, 'y', 'FaceAlpha', '0.2');
                    ylabel('Precision/Recall/F1 Score');
                    xlabel('$1/\lambda$');
                    axis([1 inf 0 inf]);
                    yyaxis right;
                    cond_plot = plot(1./l_values, meanEdgeConductance, 'b-', 'LineWidth', 1);
                    Y = [min(currEdgeConductance) fliplr(max(currEdgeConductance))];
                    patch(X, Y, 'blue', 'FaceAlpha', 0.2);
                    set(gca, 'XScale', 'log');
                    axis([0 inf 0 inf]);
                    ylabel('$\frac{E(L, R)}{\min(\pi(L)), \pi(R)) + \pi(C)}$', 'Color', [0.15, 0.15, 0.15], 'FontSize', 15);
                    legend([pre_plot rec_plot f1_plot cond_plot], {'Precision', 'Recall', 'F1-Score', 'ORC Objective'}, 'Location', 'best');
                    title({'Conductance and reconstruction metrics for synthetic dataset', sprintf('with %s communities with %s overlap starting from the %s cut', balance, o_title, a_title)});
                    exportgraphics(gcf, fullfile(inputDirectory, sprintf('greedy_overlap_fscore_%s_%s_%s_%s_%s.png', balance, p_prob, overs{o}, qs{j}, a_title)), 'Resolution', resolution)
                end
            end
        end
    end
end
    
end

