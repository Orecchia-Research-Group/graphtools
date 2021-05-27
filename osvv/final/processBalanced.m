function [] = processBalanced(graphDirectory, inputDirectory)
%PROCESSBALANCED Plot balanced results

%% Parameter checking
error_string = 'Error in parameter %s. See README file for usage.\n';

if ~ischar(graphDirectory)
    error(error_string, 'graph directory');
end

if ~ischar(inputDirectory)
    error(error_string, 'input directory');
end

%% Read in data and compute results
graphFiles = dir(fullfile(graphDirectory, '*.eg2'));
i = 1;
for f=1:length(graphFiles)
    % Load graph;
    graphFilename = graphFiles(f).name;
    [~, datasetName, ~] = fileparts(graphFilename);
    datasetPieces = split(graphFilename, '_');
    
    fprintf('Processing %s.\n', datasetName);
    
    vectorFilename = fullfile(inputDirectory, sprintf('%s.mat', datasetName));
    try
        vectors = load(vectorFilename).vec;
        temp = importdata(fullfile(graphDirectory, sprintf('%s_authorAreas.txt', datasetName)));
        category_names = split(temp.textdata{1});
        category_names = upper(category_names(2:end));
        categories = temp.data;
        
    catch
        fprintf(2, 'Unable to open vector or areas file. Skipping dataset %s.\n', vectorFilename, datasetName);
        continue;
    end 
    
    resultFiles = dir(fullfile(inputDirectory, sprintf('%s_balanced_01_*.ptn', datasetName)));
    if ~isempty(resultFiles)
        [G, n, m] = loadeg2graph(fullfile(graphDirectory, graphFilename));
        weight = full(sum(G));
        volume = sum(weight);
    end
    for f_res=1:length(resultFiles)
        resultFilename = fullfile(inputDirectory, resultFiles(f_res).name);
        [~, filename, ~] = fileparts(resultFilename);
        filenamePieces = split(filename, '_');
        partitions = readPtn(resultFilename);
        dataset{i} = datasetName;
        algorithm{i} = 'ORC-SDP';
        lambda(i) = 1 / str2double(filenamePieces{3});
        target_balance(i) = str2double(filenamePieces{4});
        balance(i) = min(sum(weight(partitions{1})), sum(weight(partitions{2}))) / volume;
        [~, ~, score(i)] = cutexp(G, 1, int64(weight), partitions{1}, partitions{2});
        overlappingNodes = intersect(partitions{1}, partitions{2});
        volume_overlap(i) = 100 * sum(weight(overlappingNodes)) / volume;
        node_overlap(i) = 100 * length(overlappingNodes) / n;
        for a=1:2
            category_counts(i, a, :) = sum(categories(partitions{a}, :));
        end
        i = i + 1;
    end
end
fprintf('Read in total %d results\n', length(dataset));

%% Plot Results
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

for f = unique(dataset)
    figure;
    ind = strcmp(f{1}, dataset);
    if sum(ind) == 0
        continue;
    end
    [~, index] = sort(target_balance(ind));
    b = balance(ind);
    s = score(ind);
    plot(b(index), s(index), 'Color', 'b', 'Marker', 'o', 'MarkerFaceColor', 'b');
    title(sprintf('Score vs Balance for the %s dataset', f{1}));
    xlabel('Balance = $\frac{\min\big(\pi(L), \pi(R)\big) + \pi(C)}{\pi(V)}$');
    ylabel('Score = $\frac{E(L, R) + \pi(C)}{\min\big(\pi(L), \pi(R)\big) + \pi(C)}$');
    %axes('Position',[.18 .68 .2 .2], 'YAxisLocation','right')
    %box on
    %plot(b(index(1:2)), s(index(1:2)), 'Color', 'b', 'Marker', 'o', 'MarkerFaceColor', 'b');
    saveas(gcf, fullfile(inputDirectory, sprintf('%s_balance.png', f{1})));
    
    % Spider plot for category
    cat_count = category_counts(ind, :, :);
    tb = target_balance(ind);
    entropy = zeros(length(index), 1);
    for i=index
        figure;
        to_plot = squeeze(cat_count(i, :, :)) ./ sum(categories);
        sp = sum(to_plot(1, :));
        pr = to_plot(1, :) / sp;
        for p=1:length(to_plot)
            if to_plot(1, p) > 0
                entropy(i) = entropy(i) - to_plot(1, p) / sp  * log2(to_plot(1, p) / sp);
            end
            if to_plot(1, p) < 1
               entropy(i) = entropy(i) - (1 - to_plot(1, p) / sp) * log2(1 - to_plot(1, p) / sp);
            end
        end
        entropy(i) = entropy(i) / length(to_plot);
        fprintf('Balance = %.3f: %.3f %.3f %.3f %.3f %.3f %.3f\n', b(i), pr(1), pr(2), pr(3), pr(4), pr(5), pr(6));
        hnd = radarplot(to_plot, category_names);
        legend(hnd, {'L', 'R'});
        title({'Percentage of papers in each category that belong', sprintf('to each community when target balance = %d', tb(i))})
        saveas(gcf, fullfile(inputDirectory, sprintf('%s_spider_%d.png', f{1}, tb(i))));
    end
    %s
    figure();
    plot(b(index), entropy(index), 'Color', 'b', 'Marker', 'o', 'MarkerFaceColor', 'b');
    title('Shannon Entropy relative to balance');
    xlabel('Balance = $\frac{\min\big(\pi(L), \pi(R)\big) + \pi(C)}{\pi(V)}$');
    ylabel('$H(p)$');
    saveas(gcf, fullfile(inputDirectory, sprintf('%s_communityEntropy.png', f{1})));
    
end
end
