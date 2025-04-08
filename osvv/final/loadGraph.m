function [G, weight] = loadGraph(fileToRead)
%LOADGRAPH Factory function to get G and node weights


loaders = dictionary([".eg2", ".metis", ".mtx"], {@loadeg2graph, @loadMetisGraph, @loadMtxGraph});
if(ischar(fileToRead))
    [~, ~, ext] = fileparts(fileToRead);
    ext = lower(ext);
    if(~isKey(loaders, ext))
        validExts = strjoin(loaders.keys, ', ');
        error('loadGraph:InvalidFileExtension', 'Invalid graph file type %s. Supported graph file extensions are %s', fileToRead, validExts);
    end
    loader = loaders(ext);
    loader = loader{1};
    [G, weight] = loader(fileToRead);
    weight = int64(weight);
else
    G = fileToRead;
    degree = int64(full(sum(G)));
    weight = ones(1, n, 'int64');
    weight(:) = degree;
end

% Check that the graph is strongly connected
[S, ~] = graphconncomp(G);
if (S > 1)
    error('Graph is not strongly connected. S = %d', S);
end

end

