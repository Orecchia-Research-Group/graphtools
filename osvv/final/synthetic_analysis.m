function [result_cell] = synthetic_analysis(folder, dataset)
    if nargin < 2
        dataset = 'balanced';
    end

    lambda_num = 1;
    lambda_den = 1;
    % List all .mtx files in the folder
    mtxFiles = dir(fullfile(folder, sprintf('*%s*.mtx', dataset)));

    if isempty(mtxFiles)
        error('No .mtx files found in %s', folder);
    end
    result_column_names = {'Name', 'Edges', 'p', 'q', 'repeat', 'Ground-truth', 'Time (s)', 'Best found', 'Lower bound', 'Accuracy', 'Precision', 'Recall', 'F1-score'};
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
        fprintf('Processing: %s\n', mtxFiles(k).name);

        % Parse corresponding .ptn filename from .mtx filename
        % e.g., synthetic_dsname_50_70_3.mtx => synthetic_dsname_50_70_3.ptn
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

        % Run cutfind
        [G, weight] = loadGraph(mtxPath);
        number_of_edges = nnz(G);
        try
            [expansionFound, ~, L, R, ~, endtime, ~, ~, ~, ~, ~, lower] = cutfind(mtxPath, stop=40, pwr_k=4, eta=1);
        

            % Read ground truth partitions
            partitions = readPtn(ptnPath);
            [A, B] = partitions{:};
    
            % Convert to logical sets
            n = max(max(A), max(B));
            A_label = false(n, 1); % ground truth for A
            A_label(A) = true;
            B_label = false(n, 1); % ground truth for B
            B_label(B) = true;
    
            L_label = false(n, 1); % predicted for L
            L_label(L) = true;
            R_label = false(n, 1); % predicted for R
            R_label(R) = true;
    
            % Check if prediction is flipped
            if sum(A_label == L_label) + sum(B_label == R_label) > sum(A_label == R_label) + sum(B_label == L_label)
                A_pred_label = L_label;
                B_pred_label = R_label;
            else
                A_pred_label = R_label;
                B_pred_label = L_label;
            end
    
            % Compute metrics (micro)
            TP = sum(A_pred_label & A_label) + sum(B_pred_label & B_label);
            FP = sum(A_pred_label & ~A_label) + sum(B_pred_label & ~B_label);
            FN = sum(~A_pred_label & A_label) + sum(~B_pred_label & B_label);
            TN = sum(~A_pred_label & ~A_label) + sum(~B_pred_label & ~B_label);
    
            accuracy = (TP + TN) / double(2 * n);
            recall = TP / (TP + FN);
            precision = TP / (TP + FP);
            f1 = 2 * precision * recall / (precision + recall + eps);
    
            [~, ~, realExp] = cutexp(G, int64(lambda_num), int64(lambda_den), int64(weight), A, B);

            % Append to results: [p, q, r, accuracy, precision, recall, f1]
            result_cell(end + 1, :) = {datasetName, number_of_edges, p, q, r, realExp, endtime, expansionFound, lower, accuracy, precision, recall, f1};
        catch ME
            warning('%s\n', getReport(ME, 'extended'));
            continue;
        end
    end
    % result_column_names = {'p', 'q', 'repeat', 'accuracy', 'precision', 'recall', 'f1'};
    % result_cell = [result_column_names; num2cell(results)];
    writecell(result_cell, fullfile(folder, sprintf('../%s_results.csv', dataset)));
end