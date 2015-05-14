function [filtMST,newConnCompNodes,new_leaf_connectivity,new_pathBetweenLeaves] = deleteNode(newLength,Label,mstf,connCompNodes, leafConnectivity,pathBetweenLeaves )
%DELETENODE Summary of this function goes here

    n = length(connCompNodes);
    newConnCompNodes = cell(1,newLength);
    new_leaf_connectivity = cell(newLength);
    new_pathBetweenLeaves = cell(newLength);
    filtMST = zeros(newLength);
    
    inner_index = 1;
    for i = 1 : n
        if Label(i) ~= 1
            newConnCompNodes{inner_index} = connCompNodes{i};
            inner_index = inner_index+1;
        end
    end
    
    
    inner_index = 1 ;
    for i = 1 : n
        for j = 1 : n
            if Label(i) ~= 1 && Label(j) ~= 1
%                 [X Y] = ind2sub(size(filtMST),inner_index)
                new_leaf_connectivity{inner_index} = leafConnectivity{j,i}; % column major order
                new_pathBetweenLeaves{inner_index} = pathBetweenLeaves{j,i} ;
                filtMST(inner_index) = mstf(j,i);
                inner_index = inner_index+1;
            end
        end
    end




end

