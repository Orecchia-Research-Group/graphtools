function [] = mustBeGraph(G)
%MUSTBEGRAPH Checks that it is a sparse square matrix

mustBeNumeric(G);
mustBeSparse(G);
mustBeSquare(G);

end

