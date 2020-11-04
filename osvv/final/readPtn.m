function [partitions] = readPtn(Filename)
%READPTN Read a ptn file and return the partitions 
%        as an array of lists of nodes
%
% INPUT:
%   (char) Filename - name of the .ptn file containing the partition
%
% OUTPUT:
%   partitions: Array of lists of nodes of the 2 partitions

if (~ischar(Filename))
    error('readPtn: Error in parameter filename. Not a char.');
end

file = fopen(Filename, 'r');

partitions{1} = [];
partitions{2} = [];
i = int64(1);
line = fgetl(file);
while(line ~= -1)
    c = textscan(line, repmat('%f ',1,2));
    m = cellfun(@numel, c);
    k = arrayfun(@(x,y)[x{1};nan(max(m)-y,1)],c,m,'un',0);
    out = cat(2,k{:});
    for p=out
        if isnan(p)
            continue;
        end
        partitions{p + 1} = [partitions{p + 1}; i];
    end
    i = i + 1;
    line = fgetl(file);
end

fclose(file);

end

