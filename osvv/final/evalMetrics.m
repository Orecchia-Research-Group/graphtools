function [accuracy, precision, recall, f1] = evalMetrics(trueSet, predSet, n)
%EVALMETRICS Computes accuracy, precision, recall, and F1-score for two sets.
%
%   Inputs:
%       trueSet   - Array or set of ground-truth elements
%       predSet   - Array or set of predicted elements
%       n         - (Optional) Number of elements
%
%   Outputs:
%       accuracy  - (TP + TN) / n
%       precision - TP / (TP + FP)
%       recall    - TP / (TP + FN)
%       f1        - 2 * precision * recall / (precision + recall)

    if nargin < 3
        n = length(union(trueSet, predSet));
    end

    n = double(n);

    % Convert to row vectors
    trueSet = unique(trueSet(:)');
    predSet = unique(predSet(:)');

    % True positives, false positives, false negatives, true negatives
    TP = double(numel(intersect(trueSet, predSet)));
    FP = double(numel(setdiff(predSet, trueSet)));
    FN = double(numel(setdiff(trueSet, predSet)));
    TN = n - double(numel(union(trueSet, predSet)));

    % Compute metrics
    if TP + FP > 0
        precision = TP / (TP + FP);
    else
        precision = 0;
    end

    if TP + FN > 0
        recall = TP / (TP + FN);
    else
        recall = 0;
    end

    if precision + recall > 0
        f1 = 2 * precision * recall / (precision + recall);
    else
        f1 = 0;
    end

    accuracy = (TP + TN) / n;
end   
