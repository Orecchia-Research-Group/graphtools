function [dataset, algorithm, balanced, C, p, q, index, lambda, score, volume_overlap] = processSynthetic(graphDirectory, inputDirectory, version)
%PROCESSSYNTHETIC Process synthetic results for k=2 communities
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
color = [0 0 1; 1 0 0; 0 1 0];
balances = {'balanced', 'unbalanced'};  % {'balanced', 'unbalanced'}
overs = {'log', 'invSqrt'};                        % {'log', 'invSqrt'}
ps = {'normal', };                        % {'normal', 'sqrt'}
qs = {'sparse', 'normal', 'dense', 'constant'};  % {'sparse', 'normal', 'dense', 'constant'}

dataset = {};
algorithm = {};
balanced = {};
C = {};
p = {};
q = {};
index = [];
lambda = [];
score = [];
volume_overlap = [];
node_overlap = [];
true_volume_overlap = [];
true_node_overlap = [];
acc = [];
precision = [];
recall = [];
f_score = [];

i = 1;
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
    current_true_partitions = readPtn(fullfile(graphDirectory, sprintf('%s.partition', datasetName)));
    current_y_true = zeros(n, 1);
    current_y_true(current_true_partitions{2}) = 1;
    current_y_true(current_true_partitions{3}) = 2;
    try
        vectors = load(vectorFilename).vec;
    catch
        fprintf(2, 'Unable to open vector file %s. Skipping dataset %s.\n', vectorFilename, datasetName);
        continue;
    end 

    resultFiles = dir(fullfile(inputDirectory, sprintf('%s_%s_*.ptn', datasetName, version)));
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
        dataset{i} = datasetName;
        algorithm{i} = 'ORC-SDP';
        balanced{i} = filenamePieces{3};
        C{i} = filenamePieces{4};
        p{i} = filenamePieces{5};
        q{i} = filenamePieces{6};
        index(i) = str2double(filenamePieces{7});
        lambda(i) = 1 / str2double(filenamePieces{9});
        
        [~, ~, edgescore(i)] = cutexp(G, int64(-1), int64(1), int64(weight), partitions{1}, partitions{2});
        [~, ~, mixedscore(i)] = cutexp(G, int64(1), int64(str2double(filenamePieces{9})), int64(weight), partitions{1}, partitions{2});
        overlappingNodes = intersect(partitions{1}, partitions{2});
        volume_overlap(i) = 100 * sum(weight(overlappingNodes)) / volume;
        node_overlap(i) = 100 * length(overlappingNodes) / n;
        
        y_true{i} = current_y_true;
        true_overlapping_nodes = current_true_partitions{3};
        true_volume_overlap(i) = 100 * sum(weight(true_overlapping_nodes)) / volume;
        true_node_overlap(i) = 100 * length(true_overlapping_nodes) / n;
        
        y_exp = zeros(n, 1);
        y_exp(R{i}) = 1;
        y_exp(overlappingNodes) = 2;
        y{i} = y_exp;
        correctness = y_true{i} == y{i};
        acc(i) = sum(correctness) / n;
        precision(i) = sum(correctness(overlappingNodes)) / length(overlappingNodes);
        recall(i) = sum(correctness(true_overlapping_nodes)) / length(true_overlapping_nodes);
        f_score(i) = 2 * precision(i) * recall(i) / (precision(i) + recall(i));
        fprintf('accuracy = %.2f precision = %.2f recall = %.2f F1 = %.2f\n', acc(i), precision(i), recall(i), f_score(i));
        i = i + 1;
    end
end
precision(isnan(precision)) = 0;
fprintf('Read in total %d results\n', length(dataset));

%% Aggregate results
xtl = {'\begin{tabular}{c} $Pr[L] = 50\%$ \\ $Pr[C] = \Theta(\frac{1}{\sqrt{n}})$ \end{tabular}', ...
       '\begin{tabular}{c} $Pr[L] = 25\%$ \\ $Pr[C] = \Theta(\frac{1}{\sqrt{n}})$ \end{tabular}', ...
       '\begin{tabular}{c} $Pr[L] = 50\%$ \\ $Pr[C] = \Theta(\frac{1}{\log(n)})$ \end{tabular}', ...
       '\begin{tabular}{c} $Pr[L] = 25\%$ \\ $Pr[C] = \Theta(\frac{1}{\log(n)})$ \end{tabular}'};

% lgd = {'$q = \Theta(1/|C|)$', '$q = \Theta(\log(|C|)/ |C|)$', '$q = \Theta(1/\sqrt{|C|})$', '$q = \Theta(1)$'};
lgd = {'$q = 3/2|C|$', '$q = 3\log(|C|) / 4|C|)$', '$q = 3 / 5\sqrt{|C|}$', '$q = 0.03$'};

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560]};
markers = {'o', 'd', 'h', 's'};
resolution = 500;
logWidth = 0.97;
lams = sort(1./unique(lambda));
lams = lams(1:21);
%{
for p_ind=1:1
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
                for l_inv=lams
                    l_index = j_index & (lambda == 1/l_inv);
                    node_acc(l_inv, j) = mean(acc(l_index));
                    node_acc_err(l_inv, j) = std(acc(l_index));
                    volm_over(l_inv, j) = mean(volume_overlap(l_index) ./ true_volume_overlap(l_index));
                    volm_over_err(l_inv, j) = std(volume_overlap(l_index) ./ true_volume_overlap(l_index));
                    node_over(l_inv, j) = mean(node_overlap(l_index) ./ true_node_overlap(l_index));
                    node_over_err(l_inv, j) = std(node_overlap(l_index) ./ true_node_overlap(l_index));
                    f_score_mean(l_inv, j) = mean(f_score(l_index));
                    f_score_err(l_inv, j) = std(f_score(l_index));
                end
            end
            f_score_mean(isnan(f_score_mean)) = 0;
            f_score_err(isnan(f_score_err)) = 0;
            figure;
            hold on;
            % plot(data, 'Marker', 'o')
            for j=1:4
                X = lams * (logWidth^(1 - 2 * (j - 1) / 3));
                errorbar(X, node_acc(lams, j), node_acc_err(lams, j), 'Color', colors{j}, 'Marker', markers{j}, 'MarkerFaceColor', 'auto', 'MarkerSize', 3);
            end
            title(sprintf('Accuracy of predicted blocks with %s communities and %s overlap', balance, o_title));
            xlabel('$1/\lambda$');
            xticks(lams);
            set(gca, 'XScale', 'log');
            axis([lams(1) * logWidth, lams(end) / logWidth, 0.95, 1.005]); 
            ylabel('Accuracy');
            legend(lgd, 'Location', 'southoutside', 'Orientation', 'horizontal');
            exportgraphics(gcf, fullfile(inputDirectory, sprintf('accuracy_%s_%s_%s.png', balance, p_prob, overs{o})), 'Resolution', resolution);

            figure;
            hold on;
            for j=1:4
                X = lams * (logWidth^(1 - 2 * (j - 1) / 3));
                errorbar(X, node_over(lams, j), node_over_err(lams, j), 'Color', colors{j}, 'Marker', markers{j}, 'MarkerFaceColor', 'auto', 'MarkerSize', 3);
            end
            set(gca, 'XScale', 'log');
            yline(1, 'r--');
            title({'Volume overlap relative to groundtruth'  sprintf('with %s communities and %s overlap', balance, o_title)});
            xlabel('$1/\lambda$');
            xticks(lams);
            ylabel('Volume overlap relative to groundtruth');
            legend(lgd, 'Location', 'southoutside', 'Orientation', 'horizontal');
            exportgraphics(gcf, fullfile(inputDirectory, sprintf('volume_overlap_%s_%s_%s.png', balance, p_prob, overs{o})), 'Resolution', resolution);

            figure;
            hold on;
            for j=1:4
                X = lams * (logWidth^(1 - 2 * (j - 1) / 3));
                errorbar(X, volm_over(lams, j), volm_over_err(lams, j), 'Color', colors{j}, 'Marker', markers{j}, 'MarkerFaceColor', 'auto', 'MarkerSize', 3);
            end
            set(gca, 'XScale', 'log');
            yline(1, 'r--');
            title({'Overlap size relative to groundtruth'  sprintf('with %s communities and %s overlap', balance, o_title)});
            xlabel('$1/\lambda$');
            xticks(lams);
            axis([lams(1) * logWidth, lams(end) / logWidth, 0, 1.6]); 
            ylabel('Overlap size relative to groundtruth');
            legend(lgd, 'Location', 'southoutside', 'Orientation', 'horizontal');
            exportgraphics(gcf, fullfile(inputDirectory, sprintf('node_overlap_%s_%s_%s.png', balance, p_prob, overs{o})), 'Resolution', resolution);
            
            figure;
            hold on;
            for j=1:4
                X = lams * (logWidth^(1 - 2 * (j - 1) / 3));
                errorbar(X, f_score_mean(lams, j), f_score_err(lams, j), 'Color', colors{j}, 'Marker', markers{j}, 'MarkerFaceColor', 'auto', 'MarkerSize', 3);
            end
            set(gca, 'XScale', 'log');
            yline(1, 'r--', 'Alpha', 0.5);
            title({'F1-score for overlap reconstruction'  sprintf('in %s communities with %s overlap', balance, o_title)});
            xlabel('$1/\lambda$');
            xticks(lams);
            axis([lams(1) * logWidth, lams(end) / logWidth, 0, 1.1]); 
            ylabel('F1-score');
            legend(lgd, 'Location', 'southoutside', 'Orientation', 'horizontal');
            exportgraphics(gcf, fullfile(inputDirectory, sprintf('overlap_fscore_%s_%s_%s.png', balance, p_prob, overs{o})), 'Resolution', resolution);

        end
    end
end
%}
%% Single F1-score with precision and recall

p_ind = 1;
p_prob = ps{p_ind};
p_index = strcmp(p, p_prob);

for b = 1:2
balance = balances{b};
b_index = p_index & strcmp(balanced, balance);

o = 2;
if strcmp(overs{o}, 'log')
    o_title = '$\log(n)$';
else
    o_title = '$\sqrt{n}$';
end
o_index = b_index & strcmp(C, overs{o});

j = 2;
j_index = b_index & strcmp(q, qs{j});
if strcmp(overs{o}, 'log')
    o_title = '$\log(n)$';
else
    o_title = '$\sqrt{n}$';
end
o_index = j_index & strcmp(C, overs{o});
for l_inv=lams
    l_index = o_index & (lambda == 1/l_inv);
    node_pre(l_inv, j) = mean(precision(l_index));
    node_pre_err(l_inv, j) = std(precision(l_index));
    node_rec(l_inv, j) = mean(recall(l_index));
    node_rec_err(l_inv, j) = std(recall(l_index));
    volm_over(l_inv, j) = mean(volume_overlap(l_index));
    volm_over_err(l_inv, j) = std(volume_overlap(l_index));
    node_over(l_inv, j) = mean(node_overlap(l_index) ./ true_node_overlap(l_index));
    node_over_err(l_inv, j) = std(node_overlap(l_index) ./ true_node_overlap(l_index));
    f_score_mean(l_inv, j) = mean(f_score(l_index));
    f_score_err(l_inv, j) = std(f_score(l_index));
    over_sc(l_inv, j) = mean(edgescore(l_index));
    over_sc_err(l_inv, j) = std(edgescore(l_index));
    over_msc(l_inv, j) = mean(mixedscore(l_index));
    over_msc_err(l_inv, j) = std(mixedscore(l_index));
    
end

f_score_mean(isnan(f_score_mean)) = 0;
f_score_err(isnan(f_score_err)) = 0;
figure;
hold on;
X = lams;
ind = lams;

set(gca, 'XScale', 'log');
plot_pre = errorbar(X, node_pre(ind, j), node_pre_err(ind, j), 'Color', colors{2}, 'Marker', markers{2}, 'MarkerFaceColor', 'auto', 'MarkerSize', 3);
plot_rec = errorbar(X, node_rec(ind, j), node_rec_err(ind, j), 'Color', colors{3}, 'Marker', markers{3}, 'MarkerFaceColor', 'auto', 'MarkerSize', 3);
plot_f1 = errorbar(X, f_score_mean(ind, j), f_score_err(ind, j), 'Color', [0, 0.5, 0], 'Marker', markers{1}, 'MarkerFaceColor', 'auto', 'MarkerSize', 3);
yline(1, 'r--', 'Alpha', 0.5);
% title({'F1-score for overlap reconstruction'  sprintf('in %s communities with %s overlap', balance, o_title)});
xlabel('$1/\lambda$');
axis([X(1) * logWidth, X(end) / logWidth, 0, 1.1]); 
ylabel('Precision/Recall/F1 Score');

% Add textarrow annotation
[y_end, l_end_ind] = max(f_score_mean(ind, j));
lim = get(gca,'Position');
y_end = lim(2) + y_end / 1.1 * lim(4);
x_end = lim(1) + log(X(l_end_ind)) / log(max(X)) * lim(3);
ta_x = [x_end + 0.075, x_end];
ta_y = [y_end - 0.2, y_end];

ta = annotation('textarrow', ta_x, ta_y, 'String', {'Close to perfect', 'reconstruction'}, 'FontWeight', 'bold');

% ORC score plot
yyaxis right;
plot_orc = errorbar(X, over_sc(ind, j), over_sc_err(ind, j), '-', 'Color', colors{1}, 'MarkerFaceColor', 'auto', 'MarkerSize', 6, 'LineWidth', 3);
%plot_msc = errorbar(X, over_msc(ind, j), over_msc_err(ind, j), '-', 'Color', 'black', 'MarkerFaceColor', 'auto', 'MarkerSize', 6, 'LineWidth', 3);
axis([1, inf, 0, inf])
% Elbow pont textarrow
lim = get(gca, 'Position');
y_end = over_sc(l_end_ind, j);
y_end = lim(2) + y_end / (max(over_sc(ind, j)) * 1.1) * lim(4);
ta_y = [y_end + 0.06, y_end];
ta = annotation('textarrow', ta_x, ta_y, 'String', {'Elbow', 'point'}, 'FontWeight', 'bold');
ylabel('$q_{G, 0}\left(\left[S, T\right]\right)$', 'Color', [0.15, 0.15, 0.15], 'FontSize', 15);

% Legend and save
legend([plot_pre, plot_rec, plot_f1, plot_orc], {'Precision', 'Recall', 'F1-Score', 'ORC Objective'}, 'Location', 'southeast');
exportgraphics(gcf, fullfile(inputDirectory, sprintf('lambda_prf_%s_%s_%s.png', balance, p_prob, overs{o})), 'Resolution', resolution);

%% Same plot against overlap
figure;
hold on;
X = volm_over(lams, j);
ind = lams;

plot_pre = errorbar(X, node_pre(ind, j), node_pre_err(ind, j), 'Color', colors{2}, 'Marker', markers{2}, 'MarkerFaceColor', 'auto', 'MarkerSize', 3);
plot_rec = errorbar(X, node_rec(ind, j), node_rec_err(ind, j), 'Color', colors{3}, 'Marker', markers{3}, 'MarkerFaceColor', 'auto', 'MarkerSize', 3);
plot_f1 = errorbar(X, f_score_mean(ind, j), f_score_err(ind, j), 'Color', [0, 0.5, 0], 'Marker', markers{1}, 'MarkerFaceColor', 'auto', 'MarkerSize', 3);
yline(1, 'r--', 'Alpha', 0.5);
title({'F1-score for overlap reconstruction'  sprintf('in %s communities with %s overlap', balance, o_title)});
xlabel('Overlap (%)');
axis([0, inf, 0, inf]); 
ylabel('Precision/Recall/F1 Score');

% Add textarrow annotation
[y_end, l_end_ind] = max(f_score_mean(ind, j));
y_end = lim(2) + y_end / 1.1 * lim(4);
lim = get(gca,'Position');
x_end = lim(1) + X(l_end_ind) / (max(X) / logWidth) * lim(3);
ta_x = [x_end - 0.1, x_end];
ta_y = [y_end - 0.3, y_end];

ta = annotation('textarrow', ta_x, ta_y, 'String', {'Close to perfect', 'reconstruction'}, 'FontWeight', 'bold');

% ORC score plot
yyaxis right;
plot_orc = errorbar(X, over_sc(ind, j), over_sc_err(ind, j), '-', 'Color', colors{1}, 'MarkerFaceColor', 'auto', 'MarkerSize', 6, 'LineWidth', 3);
axis([-inf, inf, 0, 1.2*max(over_sc(ind, j))])
% Elbow pont textarrow
lim = get(gca, 'Position');
y_end = over_sc(l_end_ind, j);
y_end = lim(2) + y_end / (max(over_sc(ind, j)) * 1.1) * lim(4);
ta_x = [x_end- 0.01, x_end];
ta_y = [y_end + 0.1, y_end];
ta = annotation('textarrow', ta_x, ta_y, 'String', {'Elbow', 'point'}, 'FontWeight', 'bold');
ylabel('$q_{G, \lambda}\left(\left[S, T\right]\right)$', 'Color', [0.15, 0.15, 0.15], 'FontSize', 15);

% Legend and save
legend([plot_pre, plot_rec, plot_f1, plot_orc], {'Precision', 'Recall', 'F1-Score', 'ORC Objective'}, 'Location', 'south');
exportgraphics(gcf, fullfile(inputDirectory, sprintf('overlap_prf_%s_%s_%s.png', balance, p_prob, overs{o})), 'Resolution', resolution);
end