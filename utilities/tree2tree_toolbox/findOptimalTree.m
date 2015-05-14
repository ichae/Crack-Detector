function [optTree,score] = findOptimalTree(treeSet , t1, t2)
%FINDOPTIMALTREE find the optimal subtree of the subtrees
%   min sum{A_i*v_i} + sum{B_j*e_j}

numTrees = length(treeSet);
score = zeros(numTrees,1);

for i = 1 : numTrees
    
    thistreeNodes = treeSet{i}.nodeList ;
    thistreePathcost = treeSet{i}.pathCost;
    score(i) = getScore(thistreeNodes,thistreePathcost,t1,t2);
    
end

[~,pos] = min(score);
optTree = treeSet(pos);

end

function val = getScore(V,E,t1,t2)

    numNodes = length(V);
    nodeVal = zeros(numNodes,1);
    
    for i = 1 : numNodes
        nodeVal(i) = V{i}.NodeValue;
    end
    
    % vectorize the edges and nodes
    % currently, neuroness is given by the node volume only...will consider
    % other later

    vectE = [];
    vectV = [];
%     nodeVal = nodeVal/(max(nodeVal(:)));
    
    for i = 1 : length(E(:))
       if E(i) ~= 0
          vectE = [vectE;E(i)]; 
       end
    end

    for i = 1 : numNodes
       if nodeVal(i) ~= 0
          vectV = [vectV;nodeVal(i)]; 
       end
    end
    
    vectE = vectE/max(vectE(:));
    vectV = vectV/max(vectV(:));
    
    A = -log((1 - exp(-(vectV)/t1))./(exp(-(vectV)/t1)));
    B = log((1-exp(-(vectE)/t2))./(exp(-(vectE)/t2)));
    
    v1 = sum(A);
    v2 = sum(B);
    val = (v1+v2);
end