function [modifiedGraph,modifiedNodeList,modifiedleafConnectivity,modifiedPathBetweenLeaves] = modifyGraph(node , G)
%MODIFYGRAPH delete 'node' from G

Mat = G.graphMat ;
nList = G.nodeList ;
leafConn   = G.leafConnectivity;
path_bet_leaves = G.pathBetweenLeaves;


num_nodes = length(nList);

% Remove the node from the adjoint matrix
Mat(node,:) = [];
Mat(:,node) = []; 
modifiedGraph = Mat;

% Remove from the node list
nList(node) = [];
modifiedNodeList = nList ;

% Remove from leaf connectivity
leafConn(node,:) = [];
leafConn(:,node) = [];
modifiedleafConnectivity = leafConn;

% Remove from pathBetweenLeaves
path_bet_leaves(node,:) = [];
path_bet_leaves(:,node) = [];
modifiedPathBetweenLeaves = path_bet_leaves;

end


