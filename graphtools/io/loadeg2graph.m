% MATLAB FUNCTION: loadeg2graph
%
% PURPOSE: loads graph from file in eg2 format to a matlab sparse array
%
% INPUTS: 
%     (char) FileToRead - a eg2 graph file to read - must be a valid graph
%
% OUTPUTS:
%     sparse (double) G - matrix holding the loaded graph
%     (double) n - number of vertices 
%     (double) m - number of arcs
%
% VARIABLES;
%    (double) data - holds temporary data from importdata
%
% NOTE: the base for the numbering of vertices in eg2 files is 0. In matlab that is
% converted to 1.


function [G, n, m] = loadeg2graph(FileToRead)

% IMPORT DATA
data = importdata(FileToRead);

% READ DATA
n = data(1,1);
m = data(1,3);

% ERROR CHECK
if(~isnumeric(n) || ~isnumeric(m) || n < 1 || m < 0)
  error('Invalid eg2 file.\n');
end

% ADD +1 - SEE NOTES
data(:,1) = data(:,1) + 1;
data(:,2) = data(:,2) + 1;


% CHECK NUMER OF EDGES IS CORRECT
if(size(data, 1) ~= m + 1) 
  error('Invalid eg2 file.\n');
end

% BUILD SPARSE - ERROR CHECKING ON THE INDICES HAPPENS HERE
G = sparse(data(2:m+1,1), data(2:m+1,2), data(2:m+1,3),n,n);

% ERROR CHECKING ON THE WEIGHTS
if(min(min(G)) < 0)
  error('Invalid eg2 file.\n');
end

%SHOULD CHECK IF SYMMETRIC

m = nnz(G)/2;


end
