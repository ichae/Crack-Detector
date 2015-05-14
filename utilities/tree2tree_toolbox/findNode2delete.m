function [node] = findNode2delete(Graph)
%FINDNODE2DELETE delete the node with the least neuronness

nodeList = Graph.nodeList ;

numNodes = length(nodeList);

for i = 1 : numNodes
    nodeVal(i) = nodeList{i}.NodeValue;
end

[val pos] = min(nodeVal);

node = pos;

end

