function [] = mustBeSquare(A)
%MUSTBESQUARE Checks that a matrix is square
    msgType = 'Input must be a two dimensional, square matrix.';
    if length(size(A)) ~= 2
        eidType = 'mustBeSquare:notTwoDimensional';
        error(eidType, msgType);
    elseif size(A, 1) ~= size(A, 2)
        eidType = 'mustBeSquare:notSquare';
        error(eidType, msgType);
    end
end

