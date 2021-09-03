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
algs = {'SweepCut', 'KernighanLin', 'ORC'};
oas = {'greedy', 'runOnce', 'ORC'};

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
lambda_num = [];
lambda_den = [];
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
    % Greedy overlap
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
            overlapAlgorithm{j} = 'greedy';
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
    
    % runOnce Overlap
    i = i + 1;
    current_true_partitions = readPtn(fullfile(graphDirectory, sprintf('%s.partition', datasetName)));
    current_y_true = zeros(n, 1);
    current_y_true(current_true_partitions{2}) = 1;
    current_y_true(current_true_partitions{3}) = 2;
    true_Cvol = sum(weight(current_true_partitions{3}));
    resultFiles = dir(fullfile(inputDirectory, sprintf('%s_SweepCut_runOnce_*.ptn', datasetName)));
    for f_res=1:length(resultFiles)
        resultFilename = fullfile(inputDirectory, resultFiles(f_res).name);
        [~, filename, ~] = fileparts(resultFilename);
        filenamePieces = split(filename, '_');
        lambda_num(i) = str2double(filenamePieces(10));
        lambda_den(i) = str2double(filenamePieces(11));
        lambda(i) = lambda_num(i) / lambda_den(i);
        if lambda(i) < 1/8
            continue;
        end
        partitions = readPtn(resultFilename);
        if partitions{1}(1) == 1
            L{i} = partitions{1};
            R{i} = partitions{2};
        else
            L{i} = partitions{2};
            R{i} = partitions{1};
        end
        overlap = intersect(partitions{1}, partitions{2});
        precission(i) = sum(weight(intersect(overlap, current_true_partitions{3}))) / true_Cvol;
        Cvol(i) = sum(weight(overlap));
        Lvol(i) = sum(weight(partitions{1}));
        Rvol(i) = sum(weight(partitions{2}));dataset{i} = datasetName;
        algorithm{i} = filenamePieces{8};
        overlapAlgorithm{i} = filenamePieces{9};
        balanced{i} = filenamePieces{3};
        C{i} = filenamePieces{4};
        p{i} = filenamePieces{5};
        q{i} = filenamePieces{6};
        runIndex(i) = str2double(filenamePieces{7});
        [~, ~, edgeConductance(i)] = cutexp(G, int64(-1), int64(1), int64(weight), partitions{1}, partitions{2});
        [~, ~, mixedConductance(i)] = cutexp(G, int64(lambda_num(i)), int64(lambda_den(i)), int64(weight), partitions{1}, partitions{2});
        i = i + 1;
    end
   
    i = i + 1;
    resultFiles = dir(fullfile(inputDirectory, sprintf('%s_infty_*.ptn', datasetName)));
    for f_res=1:length(resultFiles)
        resultFilename = fullfile(inputDirectory, resultFiles(f_res).name);
        [~, filename, ~] = fileparts(resultFilename);
        filenamePieces = split(filename, '_');
        if str2double(filenamePieces{9}) > 8
            continue;
        end
        partitions = readPtn(resultFilename);
        if partitions{1}(1) == 1
            L{i} = partitions{1};
            R{i} = partitions{2};
        else
            L{i} = partitions{2};
            R{i} = partitions{1};
        end
        Cvol(i) = sum(weight(intersect(partitions{1}, partitions{2})));
        Lvol(i) = sum(weight(partitions{1}));
        Rvol(i) = sum(weight(partitions{2}));
        dataset{i} = datasetName;
        algorithm{i} = 'ORC';
        overlapAlgorithm{i} = 'ORC';
        balanced{i} = '';
        balanced{i} = filenamePieces{3};
        C{i} = filenamePieces{4};
        p{i} = filenamePieces{5};
        q{i} = filenamePieces{6};
        runIndex(i) = str2double(filenamePieces{7});
        lambda_num(i) = 1;
        lambda_den(i) = str2double(filenamePieces{9});
        lambda(i) = lambda_num(i) / lambda_den(i);
        [~, ~, edgeConductance(i)] = cutexp(G, int64(-1), int64(1), int64(weight), partitions{1}, partitions{2});
        [~, ~, mixedConductance(i)] = cutexp(G, int64(lambda_num(i)), int64(lambda_den(i)), int64(weight), partitions{1}, partitions{2});
        i = i + 1;
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
                    for oa=1:1
                        oa_index = j_index & strcmp(overlapAlgorithm, oas{oa});
                        l_values = sort(unique(lambda(oa_index)));
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
                            cl = Lvol(r_index);
                            cr = Rvol(r_index);
                            cc = Cvol(r_index);
                            currC = interp1(lambda(r_index) + noise, cc, l_values, 'linear', cc(end));
                            currL = interp1(lambda(r_index) + noise, cl, l_values, 'linear', cl(end));
                            currR = interp1(lambda(r_index) + noise, cr, l_values, 'linear', cr(end));
                            currMixedConductance(r, :) = currC .* l_values ./min(currL, currR) + currEdgeConductance(r, :);
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
                        cond_plot = plot(1./l_values, meanMixedConductance, 'b-', 'LineWidth', 1);
                        Y = [min(currMixedConductance) fliplr(max(currMixedConductance))];
                        patch(X, Y, 'blue', 'FaceAlpha', 0.2);
                        set(gca, 'XScale', 'log');
                        axis([0 inf 0 inf]);
                        % title({'Conductance and reconstruction metrics for synthetic dataset', sprintf('with %s communities with %s overlap starting from the %s cut', balance, o_title, a_title)});
                        ylabel('$q_{G, \lambda}\left(\left[S, T\right]\right)$', 'Color', [0.15, 0.15, 0.15], 'FontSize', 15);
                        legend([pre_plot rec_plot f1_plot cond_plot], {'Precision', 'Recall', 'F1-Score', 'ORC Objective'}, 'Location', 'best');
                        exportgraphics(gcf, fullfile(inputDirectory, sprintf('greedy_overlap_fscore_%s_%s_%s_%s_%s.png', balance, p_prob, overs{o}, qs{j}, a_title)), 'Resolution', resolution)
                    end
                end
            end
        end
    end
end

%% SweepCut greedy vs runOnce
colors = {'r', 'g', 'b'};

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

            for j=1:4
                legd{j} = sprintf('C edges are %s', qs{j});
                j_index = o_index & strcmp(q, qs{j});
                figure;
                hold on;
                plt = [];
                lgd = {};
                for a=[1 3]
                    a_title = algs{a};
                    a_index = j_index & strcmp(algorithm, algs{a});
                    for oa=1:3
                        oa_index = a_index & strcmp(overlapAlgorithm, oas{oa});
                        l_values = sort(unique(lambda(oa_index)));
                        points = length(l_values);
                        currPrecission = zeros(5, points);
                        currEdgeConductance = zeros(5, points);
                        currMixedConductance = zeros(5, points);
                        
                        for r=1:5
                            r_index = oa_index & (runIndex == r);
                            if sum(r_index) == 0
                                continue;
                            end
                            noise = rand(1, sum(r_index))*1e-6;
                            if a == 1
                                pre = precission(r_index);
                                currPrecission(r, :) = interp1(lambda(r_index) + noise, pre, l_values, 'linear', pre(end));
                            end
                            ec = edgeConductance(r_index);
                            currEdgeConductance(r, :) = interp1(lambda(r_index) + noise, ec, l_values, 'linear', ec(end));
                            cl = Lvol(r_index);
                            cr = Rvol(r_index);
                            cc = Cvol(r_index);
                            currC = interp1(lambda(r_index) + noise, cc, l_values, 'linear', cc(end));
                            currL = interp1(lambda(r_index) + noise, cl, l_values, 'linear', cl(end));
                            currR = interp1(lambda(r_index) + noise, cr, l_values, 'linear', cr(end));
                            currMixedConductance(r, :) = currC .* l_values ./min(currL, currR) + currEdgeConductance(r, :);
                        end
                        if a == 1
                            currPrecission(isnan(currPrecission)) = 0;
                            meanPrecission = mean(currPrecission);
                        end
                        currEdgeConductance(isnan(currEdgeConductance)) = 0;
                        currMixedConductance(isnan(currMixedConductance)) = 0;
                        meanEdgeConductance = mean(currEdgeConductance);
                        meanMixedConductance = mean(currMixedConductance);
                        if length(l_values) == 0
                            continue;
                        end
                        if (a == 1) && (max(max(currPrecission)) > 0)
                            fprintf('Prec = %f\n', max(max(currPrecission)));
                        end
                        X = [1./l_values fliplr(1./l_values)];
                        plt(oa) = plot(1./l_values, meanMixedConductance, colors{oa}, 'LineWidth', 1);
                        lgd{oa} = oas{oa};
                        Y = [min(currMixedConductance) fliplr(max(currMixedConductance))];
                        patch(X, Y, colors{oa}, 'FaceAlpha', 0.2);                        
                    end                    
                end
                set(gca, 'XScale', 'log');
                % title({'Conductance relative to $1/\lambda$ in different ways to get an overlapping cut from', sprintf('SweepCut for synthetic dataset with %s communities and %s overlap', balance, o_title)});
                axis([0 inf 0 inf]);
                xlabel('$1/\lambda$');
                ylabel('$q_{G, \lambda}\left(\left[S, T\right]\right)$', 'FontSize', 15);
                legend(plt, {'${\tt SweepCut} + {\tt GreedyImprove}$', '${\tt SweepCut} + {\tt OverlapImprove}$', '${\tt cm + improve}$'}, 'Location', 'best');
                exportgraphics(gcf, fullfile(inputDirectory, sprintf('sweepCut_%s_%s_%s_%s.png', balance, p_prob, overs{o}, qs{j})), 'Resolution', resolution)
            end
        end
    end
end
    
end

