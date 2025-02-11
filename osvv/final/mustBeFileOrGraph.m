function [] = mustBeFileOrGraph(filename)
%MUSTBEFILEORGRAPH Checks that either mustBeFile or mustBeGraph checks out.

if isa(filename, 'char') || (size(filename, 1) == 1)
    mustBeFileOrID(filename);
else
    mustBeGraph(filename)
end

