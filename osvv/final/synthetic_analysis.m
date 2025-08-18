function [result_cell] = synthetic_analysis(folder, dataset, lambda_num, lambda_den, options)
arguments
    folder char
    dataset char
    lambda_num (1, :) int64 = 1
    lambda_den (1, :) int64 = 1
    options.which {mustBeMember(options.which, {'a', 'C', 'non overlapping'})} = 'non overlapping'
    options.q (1, :) double = -1
end
    if nargin < 2
        dataset = 'balanced';
    end

    if nargin < 3
        lambda_num = 1;
        lambda_den = 1;
    end

    % List all .mtx files in the folder
    mtxFiles = dir(fullfile(folder, sprintf('*%s*.mtx', dataset)));

    if isempty(mtxFiles)
        error('No .mtx files found in %s', folder);
    end
    result_column_names = {'Name', 'lambda', 'lambda_num', 'lambda_den', 'bestExpNum', 'bestExpDen', 'bestExp', 'p', 'q', 'repeat', 'Ground-truth', 'Time (s)', 'Best found', 'Lower bound', 'Accuracy', 'Precision', 'Recall', 'F1-score', 'LL', 'LC', 'LR', 'CL', 'CC', 'CR', 'RL', 'RC', 'RR'};
    result_cell = result_column_names;

    % Start a parpool for the parfor
    p = gcp('nocreate'); % Attempt to get the current pool without creating one
    if isempty(p)
        % No parallel pool exists, so create one
        parpool(4); % This will create a parallel pool with 4 workers
    end

    % Loop through each graph
    for k = 1:length(mtxFiles)
        mtxName = mtxFiles(k).name;
        mtxPath = fullfile(folder, mtxFiles(k).name);

        % Parse corresponding .ptn filename from .mtx filename
        % e.g., synthetic_balanced_50_70_3.mtx => synthetic_balanced_50_70_3.ptn
        [~, baseName, ~] = fileparts(mtxFiles(k).name);
        ptnPath = fullfile(folder, [baseName '.ptn']);

        if ~isfile(ptnPath)
            warning('Missing partition file: %s', ptnPath);
            continue;
        end

        % Extract metadata from filename: synthetic_{dataset}_{p}_{q}_{r}.mtx
        tokens = regexp(mtxName, 'synthetic_([^_]+)_([0-9]+)_([0-9]+)_([0-9]+)\.mtx', 'tokens');
        if isempty(tokens)
            warning('Could not parse filename: %s', mtxName);
            continue;
        end
        tokens = tokens{1};
        datasetName = tokens{1};
        p = str2double(tokens{2});
        q = str2double(tokens{3});
        r = str2double(tokens{4});

        % Check if processing this file
        if (options.q ~= -1 && ~ismember(q, options.q))
            fprintf('Skipping: %s\n', mtxFiles(k).name)
            continue;
        end

        fprintf('Processing: %s\n', mtxFiles(k).name);

        % Run cutfind
        [G, weight] = loadGraph(mtxPath);
        number_of_edges = nnz(G);

        % Read ground truth partitions
        partitions = readPtn(ptnPath);
        [A, B] = partitions{:};
        n = max(max(A), max(B));

        for i=1:length(lambda_num)
            lam_num = int64(lambda_num(i));
            lam_den = int64(lambda_den(i));
            lam = double(lam_num) / double(lam_den);
            [expansionFound, edgesCut, S, T, ~, endtime, ~, ~, ~, ~, ~, lower] = cutfind(mtxPath, stop=40, pwr_k=4, eta=1, lambda_num=lam_num, lambda_den=lam_den);

            % Check if prediction is flipped
            if length(intersect(A, S)) + length(intersect(B, T)) < length(intersect(A, T)) + length(intersect(B, S))
                [S, T] = deal(T, S);
            end

            % Compute metrics (micro)
            accuracy = [];
            precision = [];
            recall = [];
            f1 = [];

            if ~strcmp(options.which, 'C')
                [accuracy(end + 1), precision(end + 1), recall(end + 1), f1(end + 1)] = evalMetrics(setdiff(A, B), setdiff(S, T), n);
                [accuracy(end + 1), precision(end + 1), recall(end + 1), f1(end + 1)] = evalMetrics(setdiff(B, A), setdiff(T, S), n);
            end
            if ~strcmp(options.which, 'non overlapping')
                [accuracy(end + 1), precision(end + 1), recall(end + 1), f1(end + 1)] = evalMetrics(intersect(A, B), intersect(S, T), n);
            end

            accuracy = mean(accuracy);
            precision = mean(precision);
            recall = mean(recall);
            f1 = mean(f1);

            [~, ~, realExp] = cutexp(G, int64(lam_num), int64(lam_den), int64(weight), A, B);
            [bestExpNum, bestExpDen, bestExp] = cutexp(G, int64(lam_num), int64(lam_den), int64(weight), S, T);


            LL = sum(weight(intersect(setdiff(A, B), setdiff(S, T))));
            LC = sum(weight(intersect(setdiff(A, B), intersect(S, T))));
            LR = sum(weight(intersect(setdiff(A, B), setdiff(T, S))));
            CL = sum(weight(intersect(intersect(A, B), setdiff(S, T))));
            CC = sum(weight(intersect(intersect(A, B), intersect(S, T))));
            CR = sum(weight(intersect(intersect(A, B), setdiff(T, S))));
            RL = sum(weight(intersect(setdiff(B, A), setdiff(S, T))));
            RC = sum(weight(intersect(setdiff(B, A), intersect(S, T))));
            RR = sum(weight(intersect(setdiff(B, A), setdiff(T, S))));
            % Append to results: [p, q, r, accuracy, precision, recall, f1]
            result_cell(end + 1, :) = {datasetName, lam, lam_num, lam_den, bestExpNum, bestExpDen, bestExp, p, q, r, realExp, endtime, expansionFound, lower, accuracy, precision, recall, f1, LL, LC, LR, CL, CC, CR, RL, RC, RR};
        end
    end
    % result_column_names = {'p', 'q', 'repeat', 'accuracy', 'precision', 'recall', 'f1'};
    % result_cell = [result_column_names; num2cell(results)];
    writecell(result_cell, fullfile(folder, sprintf('../%s_%03d_results.csv', dataset, options.q)));
end