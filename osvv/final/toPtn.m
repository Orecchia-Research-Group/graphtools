function [ ] = toPtn(outputFilename, partitions)
% toPtn  Takes as input some (possibly overlapping) partitions and writes 
%   them to the filename provided
%
% INPUTS:
%   (char) outputFilename - filename for the partitions to be written
%   (int cell) partitions - cell i contains the i-th sparse partition of
%   the graph
%
% OUTPUTS:
%
% DESCRIPTION:
%   Every line contains space separated the communities that the
%   corresponding node belongs to. For non overlapping communitiees each
%   line will contain exactly one number between 0 and number of clusters
%   minus 1.
%   The filename is encuraged to be of the form dataset.clusterNumber.ptn

n = max([partitions{:}]);
p = length(partitions);
if ischar(outputFilename)
    outputFile = fopen(outputFilename, 'w');
else
    outputFile = outputFilename;
end
if outputFile == -1
    fprintf(2, 'Failed to open %s\n', outputFilename);
end
partitionIndex = ones(p, 1);
for node=1:n
    clusterList = [];
    for i=1:p
        if (partitionIndex(i) <= length(partitions{i})) && (partitions{i}(partitionIndex(i)) < node)
            partitionIndex(i) = partitionIndex(i) + 1;
        end
        
        if (partitionIndex(i) <= length(partitions{i})) && (partitions{i}(partitionIndex(i)) == node)
            clusterList = [clusterList i-1];
        end
    end
    if isempty(clusterList)
        fprintf(2, 'Node %d does not belong to any cluster\n', node);
    end
    fprintf(outputFile, '%d', clusterList(1));
    for c=clusterList(2:end)
        fprintf(outputFile, ' %d', c);
    end
    fprintf(outputFile, '\n');
end
if ischar(outputFilename)
    fclose(outputFile);
end

end

