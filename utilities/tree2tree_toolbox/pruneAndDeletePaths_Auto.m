function [filtMST,newConnCompNodes,new_leaf_connectivity,new_pathBetweenLeaves] = pruneAndDeletePaths_Auto(MST, connCompNodes , leafConnectivity , alpha, numCC,pathcost, pathBetweenLeaves)



if alpha == 0  % no change is tolerated
    filtMST = MST+MST';
    newConnCompNodes = connCompNodes;
    new_leaf_connectivity = leafConnectivity;
    new_pathBetweenLeaves = pathBetweenLeaves;
    filtMST = sparse(filtMST);
    fprintf('\n Not deleting any edges');
else
    
    % we keep deleting till the change in neurone-ness is more than alpha

    initNurenoness = findTotalNeuroneness(connCompNodes);
    
    for i = 1 : numCC-1
        topK = i;
        echo = 0;
        [~,pruned_newConnCompNodes,~,pruned_new_pathBetweenLeaves] = pruneAndDeletePaths(MST, connCompNodes , leafConnectivity , topK, numCC,pathcost, pathBetweenLeaves,echo);
        postNurenoness = findTotalNeuroneness(pruned_newConnCompNodes);
        change = (initNurenoness-postNurenoness)/initNurenoness;
        fprintf('\n change from deleting %d to %d is %f ',i-1,i,change);
        if change >= alpha
            topK = topK-1;
            break;
        end
        initNurenoness = postNurenoness;
        
    end
    
    fprintf('\n Automatically selected %d nodes to delete',topK);
    echo=1;
    [filtMST,newConnCompNodes,new_leaf_connectivity,new_pathBetweenLeaves] = pruneAndDeletePaths(MST, connCompNodes , leafConnectivity , topK, numCC,pathcost, pathBetweenLeaves,echo);
    
end


end


function val = findTotalNeuroneness(CC_nodes)

n_nodes = length(CC_nodes(:));

val=0;
for j = 1:n_nodes
    val=val+ CC_nodes{j}.NodeValue;
end


end