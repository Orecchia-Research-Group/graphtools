function [partitions] = readGexf(inputFilename)
%READGEXF Reads the communities from the godawful .gexf xml file
%         that Jaewon Yang of Stanford thought was "a good idea at a time"*
%   *Not an actual quote.



partitions = {};
nodeCount = 0;
graph = xmlread(inputFilename);
nodes = graph.item(0).item(1).item(3);
nodeLength = nodes.getLength;
for i=1:2:nodeLength - 2
    node = nodes.item(i);
    nodeId = str2double(node.getAttributes.item(1).getValue) + 1;
    nodeCount = nodeCount + 1;
    if node.getLength < 6
        continue;
    end
    communities = node.item(5);
    for j=1:2:communities.getLength-2
        community = str2double(communities.item(j).getAttributes.item(0).getValue) + 1;
        if length(partitions) < community
            partitions{community} = [];
        end
        partitions{community} = [partitions{community}; nodeId];
    end
end
end
