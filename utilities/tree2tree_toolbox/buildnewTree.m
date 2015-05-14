function [filtMST,newConnCompNodes,new_leaf_connectivity,new_pathBetweenLeaves] = buildnewTree(node2keep,MST,ConnCompNodes,leaf_connectivity,pathBetweenLeaves)
%BUILDNEWTREE build the tree keeping only the specified nodes (ie the heaviest connected graph)

nodes_in_tree = length(node2keep(:));

filtMST = zeros(nodes_in_tree);
newConnCompNodes = cell(1,nodes_in_tree);
new_leaf_connectivity = cell(nodes_in_tree);
new_pathBetweenLeaves = cell(nodes_in_tree);

% Build the new MST
for i = 1 : nodes_in_tree
    node1 = node2keep(i);
    newConnCompNodes{i} = ConnCompNodes{node1};
    for j = 1 : nodes_in_tree
       node2 = node2keep(j);
       filtMST(i,j) = MST(node1,node2);
       new_leaf_connectivity{i,j} = leaf_connectivity{node1,node2};
       new_pathBetweenLeaves{i,j} = pathBetweenLeaves{node1,node2};
    end
end

% filtMST = sparse(filtMST);

end

