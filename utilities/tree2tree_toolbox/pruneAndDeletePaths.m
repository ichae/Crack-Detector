function [filtMST,newConnCompNodes,new_leaf_connectivity,new_pathBetweenLeaves] = pruneAndDeletePaths(MST, connCompNodes , leafConnectivity , topK, numCC,pathcost, pathBetweenLeaves,echo)
%DELETEPATHS Prune the graph and delete paths
%   filtMST --  The final post deletion MST
%   newConnCompNodes -- the final set of CC structure
%   new_leaf_connectivity --    the final connectivity between the leaves
%   new_pathBetweenLeaves --    the paths between the remaining leaves

if topK >= numCC - 1
    topK = numCC -1;
    if echo
        fprintf('\n Too many nodes to delete, will preserve only one component');
    end
end

mstf = full(MST);
mstf = mstf+mstf';

n = length(connCompNodes); % number of connected components
nzmask = (tril(mstf) > 0);

% ------------------THE DELETION STARTS HERE---------------------------
% -------------------------------------------------------------------------
% Sort according to the maximum path cost for deletion
if topK > 0 % code for deletion
    
    [val , IX] = sort(pathcost(:),'descend');
    [leftCC rightCC] = ind2sub(size(pathcost),IX);  % get the index of the connected component
    
    nodeMarker = zeros(n,1);
    num_deleted = 0;
    
    for i = 1 : topK
        CC1 = leftCC(i);
        CC2 = rightCC(i);
        
        % now remove the smaller of the two connected components
        
        CC1_wt = connCompNodes{CC1}.NodeValue;
        CC2_wt = connCompNodes{CC2}.NodeValue;
        
        % Mark which nodes to delete
        if CC1_wt <= CC2_wt       % delete CC1
            if nodeMarker(CC1) == 0
                mstf(CC1,:) = 0 ; % remove all other connections with CC1
                mstf(:,CC1) = 0 ; % remove all other connected to CC1
                nodeMarker(CC1) = 1 ; % indicate that this node is removed
                if echo
                    fprintf('\n Node %d to be deleted',CC1);
                end
                num_deleted = num_deleted+1;
            end
        else
            if nodeMarker(CC2) == 0
                mstf(CC2,:) = 0 ; % remove all other connections with CC1
                mstf(:,CC2) = 0 ; % remove all other connected to CC1
                nodeMarker(CC2) = 1 ; % indicate that this node is removed
                if echo
                    fprintf('\n Node %d to be deleted',CC2);
                end
                num_deleted = num_deleted+1;
            end
        end
        
    end
    
    % We now need to remove the nodes from 1. leaf_connectivity and 2.
    % connCompNodes and 3. MST
    
    newLength = n-num_deleted;
    [filtMST,newConnCompNodes,new_leaf_connectivity,new_pathBetweenLeaves] = deleteNode(newLength,nodeMarker,mstf,connCompNodes,leafConnectivity,pathBetweenLeaves);
    
    % WE NEED TO CHECK IF THE RESULTING GRAPH IS DISCONNECTED. IF
    % DISCONNECTED, WE PRESERVE THE LARGEST CONNECTED COMPONENT
    
    [numGraphComp nodeLabel] = graphconncomp(sparse(filtMST),'Weak','true');
    
    if numGraphComp > 1
        if echo
            fprintf('\n The tree has been disconnected, will preserve the largest subtree');
        end
        Component = struct('Nodes',{},'Weight',{});
        for i = 1 : numGraphComp
            Component(i).Nodes = find(nodeLabel == i);
            num_nodes = length(Component(i).Nodes(:)) ;
            wt = 0;
            for j = 1 : num_nodes
                node = Component(i).Nodes(j);
                wt = wt + newConnCompNodes{node}.NodeValue;
                % if we want to do by tree2tree paper.........
                % len = length(newConnCompNodes{node}.NodePosition.x);
                % wt = wt + (newConnCompNodes{node}.NodeValue)^.33 + len;
            end
            Component(i).Weight = wt;
        end
        
        CC_weight = zeros(numGraphComp,1);
        for i = 1 : numGraphComp
            CC_weight(i) =  Component(i).Weight;
        end
        
        [maxWt wtInd] = max(CC_weight(:));
        if echo
            fprintf('\n Preserving the component with weight %d',maxWt);
        end
        node2keep = Component(wtInd).Nodes;
        
        [filtMST,newConnCompNodes,new_leaf_connectivity,new_pathBetweenLeaves] = buildnewTree(node2keep,filtMST,newConnCompNodes,...
            new_leaf_connectivity,new_pathBetweenLeaves);
    end
end

if topK == 0
    filtMST = MST+MST';
    newConnCompNodes = connCompNodes;
    new_leaf_connectivity = leafConnectivity;
    new_pathBetweenLeaves = pathBetweenLeaves;
    if echo
        fprintf('\n Not deleting any edges');
    end
end

filtMST = sparse(filtMST);





end

