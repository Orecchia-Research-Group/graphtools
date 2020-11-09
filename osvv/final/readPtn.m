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

data = importdata(Filename);
partitions{1} = find(data);
partitions{2} = find(~data);

end
