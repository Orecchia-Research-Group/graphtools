function [] = processHypergraph(graphDirectory, inputDirectory)

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

i = 1;
graphFiles = dir(fullfile(graphDirectory, '*.heg2'));
for f=1:length(graphFiles)
    graphFilename = graphFiles(f).name;
    [~, datasetName, ~] = fileparts(graphFilename);
    [G, weights] = loadMetisGraph(fullfile(graphDirectory, graphFilename));
    n = length(G);
    row_sum = sum(G, 1)';
    col_sum = sum(G, 2);
    volume = sum(weights);
    categories = importdata('../../graphs/social/hyperExtractedDblp_authorAreas.txt') + 1;
    category_names = {'THEORY', 'DMW', 'IVC', 'ML', 'IDM', 'NET'};
    total_cat_count = sum(categories == 1:6)' / 2;
    resultFiles = dir(fullfile(inputDirectory, sprintf('%s_balanced_*_*_*.ptn', datasetName)));
    for f_res=1:length(resultFiles)
        resultFilename = fullfile(inputDirectory, resultFiles(f_res).name);
        [~, filename, ~] = fileparts(resultFilename);
        filenamePieces = split(filename, '_');
        lambda_num(i) = str2double(filenamePieces{3});
        lambda_den(i) = str2double(filenamePieces{4});
        lambda(i) = lambda_num(i) / lambda_den(i);
        balances(i) = str2double(filenamePieces{5});
        partitions = readPtn(resultFilename);
        L{i} = partitions{1};
        R{i} = partitions{2};
        inter = intersect(L{i}, R{i});
        Lmask = zeros(n, 1, 'logical');
        Rmask = zeros(n, 1, 'logical');
        Cmask = zeros(n, 1, 'logical');
        Lmask(L{i}) = true;
        Rmask(R{i}) = true;
        Cmask(inter) = true;
        authorCount = sum(weights > 0);
        ind = 1:(n-authorCount)/2;
        node = authorCount + 2 * ind - 1;
        cat_id{1} = Lmask(node) .* Lmask(node + 1) .* (sum(G(node + 1, L{i}), 2) == col_sum(node + 1));
        cat_id{2} = Rmask(node) .* Rmask(node + 1) .* (sum(G(node + 1, R{i}), 2) == col_sum(node + 1));
        for p=1:2
            category_counts(i, p, :) = sum(categories(2 * find(cat_id{p})) == 1:6);
        end
        category_counts(i, 3, :) = total_cat_count - squeeze(sum(category_counts(i, :, :)));
        to_plot = squeeze(category_counts(i, :, :)) ./ total_cat_count';
        figure;
        hold on;
        hnd = radarplot(to_plot, category_names);
        legend(hnd, {'L', 'R', 'Cut'});
        title({'Category distribution of paper hyperedges for', sprintf('target balance = %.3f and lambda = %.0f / %.0f = %.6f', balances(i) / 1000, lambda_num(i), lambda_den(i), lambda_num(i) / lambda_den(i))});
        saveas(gcf, fullfile(inputDirectory, sprintf('%s_spider_%d_%d_%d.png', datasetName, lambda_num(i), lambda_den(i), balances(i))));
        close;
        i = i + 1;
end
end
end

