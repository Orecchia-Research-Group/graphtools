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
        part_mask = zeros(n, 1);
        part_mask(partitions{1}) = 1;
        if sum(part_mask(1:1000)) > 700
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
        
        [~, ~, score(i)] = cutexp(G, 0, int64(weight), partitions{1}, partitions{2});
        overlappingNodes = intersect(partitions{1}, partitions{2});
        volume_overlap(i) = 100 * sum(weight(overlappingNodes)) / volume;
        node_overlap(i) = 100 * length(overlappingNodes) / n;
        
        L_sum = full(sum(G(1:4000, :)));
        if strcmp(balanced{i}, 'balanced')
            R_sum = full(sum(G(5000:9700, :)));
        else
            R_sum = full(sum(G(7500:9700, :)));
        end
        L_avg = mean(L_sum);
        R_avg = mean(R_sum);
        L_mask = L_sum > (L_avg / 2);
        R_mask = R_sum > (R_avg / 2);
        if strcmp(balanced{i}, 'balanced')
            fprintf('high L = %d low R = %d ', sum(L_mask(5000:9700)), sum(R_mask(1:4500)));
            L_mask(5100:9700) = false;
        else
            fprintf('high L = %d low R = %d ', sum(L_mask(7500:9700)), sum(R_mask(1:4500)));
            L_mask(7600:9700) = false;
        end
        R_mask(1:4500) = false;
        curr_y = zeros(n, 1);
        curr_y(R_mask) = 1;
        curr_y(L_mask & R_mask) = 2;
        y_true{i} = curr_y;
        true_volume_overlap(i) = 100 * sum(weight(L_mask & R_mask)) / volume;
        true_node_overlap(i) = 100 * sum(L_mask & R_mask) / n;
        
        y_exp = zeros(n, 1);
        y_exp(R{i}) = 1;
        y_exp(overlappingNodes) = 2;
        y{i} = y_exp;
        acc(i) = sum(y_true{i} == y{i}) / n;
        fprintf('Lmask = %d Rmask = %d L = %d R = %d i = %d acc = %f\n', sum(L_mask), sum(R_mask), length(L{i}), length(R{i}), index(i), acc(i));
        
        i = i + 1;
    end
end
fprintf('Read in total %d results\n', length(dataset));

% %% Aggregate results
xtl = {'\begin{tabular}{c} $Pr[L] = 50\%$ \\ $Pr[C] = \Theta(\frac{1}{\sqrt{n}})$ \end{tabular}', ...
       '\begin{tabular}{c} $Pr[L] = 25\%$ \\ $Pr[C] = \Theta(\frac{1}{\sqrt{n}})$ \end{tabular}', ...
       '\begin{tabular}{c} $Pr[L] = 50\%$ \\ $Pr[C] = \Theta(\frac{1}{\log(n)})$ \end{tabular}', ...
       '\begin{tabular}{c} $Pr[L] = 25\%$ \\ $Pr[C] = \Theta(\frac{1}{\log(n)})$ \end{tabular}'};

lgd = {'q=$\Theta(1/|C|)$', '$q = \Theta(\log(|C|)/ |C|)$', '$q = \Theta(1/\sqrt{|C|})$', '$q = \Theta(1)$'};

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
 
% for l=[1/2 1/3 1/4]
%     l_index = lambda == l;
%     for p_prob=ps
%         p_title = 'dense';
%         if strcmp(p_prob{1}, 'normal')
%             p_title = 'normal';
%         end
%         p_index = l_index & strcmp(p, p_prob{1});
%         %fprintf('p index = %d\n', sum(p_index));
%         toBar = zeros(4, 4);
%         eb = zeros(4, 4);
%         for j=1:4
%             j_index = p_index & strcmp(q, qs{j});
%             %fprintf('j index = %d\n', sum(j_index));
%             for o=1:2
%                 o_index = j_index & strcmp(C, overs{o});
%                 %fprintf('o index = %d\n', sum(o_index));
%                 for b=1:2
%                     index =  o_index & strcmp(balanced, balances{b});
%                     % fprintf('Files fitting lambda=%f p=%s q=%s overlap=%s balance=%s: %d\n', ...
%                     %     l, p_prob{1}, qs{j}, overs{o}, balances{b}, sum(index));
%                     toBar(j, (o - 1) * 2 + b) = mean(volume_overlap(index) ./ true_volume_overlap(index))  ;
%                     eb(j, (o - 1) * 2 + b) = std(volume_overlap(index) ./ true_volume_overlap(index));
%                     % pause;
%                 end
%             end   
%         end
% 
%         figure;
%         hold on;
%         pl = bar(toBar);
%         [ngroups, nbars] = size(toBar);
%         groupwidth = min(0.8, nbars/(nbars + 1.5));
%         for i = 1:nbars
%             % Calculate center of each bar
%             x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%             errorbar(x, toBar(:,i), eb(:,i), 'k', 'linestyle', 'none');
%         end
% 
%         set(gca, 'TickLabelInterpreter', 'latex');
%         xticks(1:4);
%         xticklabels(xtl);
%         title(sprintf('Overlap (%%) in %s communities for lambda = %.3f', p_title, l));
%         %xtickangle(15);
%         legend(pl, lgd, 'Interpreter','latex', 'Location', 'southoutside', 'Orientation', 'horizontal');
%         exportgraphics(gcf, fullfile(inputDirectory, sprintf('%s_%.0f.png', p_prob{1}, 1/l)));
%         hold off;
%     end
% end

%% Plot Accuracy and Overlap
p_prob = 'normal';
logWidth = 0.97;
lams = sort(1./unique(lambda));
lams = lams(1:21);
for b=1:2
    balance = balances{b};
    for o=1:2
        if strcmp(overs{o}, 'log')
            o_title = '$\log(n)$';
        else
            o_title = '$\sqrt{n}$';
        end
        o_index = strcmp(balanced, balance) & strcmp(p, p_prob) & strcmp(C, overs{o});
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
            end
        end
        figure;
        hold on;
        % plot(data, 'Marker', 'o')
        for j=1:4
            X = lams * (logWidth^(1 - 2 * (j - 1) / 3));
            errorbar(X, node_acc(lams, j), node_acc_err(lams, j), 'Marker', 'o', 'MarkerFaceColor', 'auto');
        end
        xticks(lams);
        title(sprintf('Accuracy of predicted blocks with %s communities and %s overlap', balance, o_title));
        xlabel('$1/\lambda$');
        xticks(lams);
        set(gca, 'XScale', 'log');
        axis([lams(1)-0.3 lams(end)+0.3 0.95, 1]); 
        ylabel('Accuracy');
        legend(lgd, 'Location', 'southoutside', 'Orientation', 'horizontal');
        exportgraphics(gcf, fullfile(inputDirectory, sprintf('accuracy_%s_%s_%s.png', balance, p_prob, overs{o})));
        
        figure;
        hold on;
        for j=1:4
            X = lams * (logWidth^(1 - 2 * (j - 1) / 3));
            errorbar(X, volm_over(lams, j), volm_over_err(lams, j), 'Marker', 'o', 'MarkerFaceColor', 'auto');
        end
        xticks(lams);
        set(gca, 'XScale', 'log');
        yline(1, 'r--');
        title(sprintf('Volume of predicted blocks relative to groundtruth with %s communities and %s overlap', balance, o_title));
        xlabel('$1/\lambda$');
        xticks(lams);
        axis([lams(1)-0.3 lams(end)+0.3 0 1.6]); 
        ylabel('Volume of overlap relative to groundtruth');
        legend(lgd, 'Location', 'southoutside', 'Orientation', 'horizontal');
        exportgraphics(gcf, fullfile(inputDirectory, sprintf('volume_overlap_%s_%s_%s.png', balance, p_prob, overs{o})));
                
    end
end
end

