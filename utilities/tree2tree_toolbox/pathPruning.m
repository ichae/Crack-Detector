function [filtMST,newConnCompNodes,new_leaf_connectivity,new_pathBetweenLeaves] = pathPruning(Img, MST, connCompNodes , leafConnectivity , topK, numCC, eImg,offset,zoffset)
%PATHPRUNING Detect the presence of a path
%       Img:              Original Image
%       MST:              The generated Minimum Spanning Tree
%       connCompNodes:    A cell of structures for each connected components
%       leafConnectivity: The connectivity of the leaves, by considering the
%                         neuron as a collection of connected components
%       topK :            The number of paths to remove
%       new_pathBetweenLeaves: Gives the shortest path between two CC


if topK >= numCC - 1
    topK = numCC -1;
    fprintf('\n Too many nodes to delete, will preserve only one component');
end

mstf = full(MST);
mstf = mstf+mstf';

n = length(connCompNodes); % number of connected components
init_num_edges = n-1;      % initial number of edges
nzmask = (tril(mstf) > 0);
linIndex=find(nzmask);     % gives the linear indices of the position of edges

pathBetweenLeaves = cell(n) ; % store the shortest path between the leaves

%-------------------------------------------------------------------------
% Find the Largest connected component
prev_wt = 0;
nodeWt = zeros(n,1);
for i = 1 : n
    nodeWt(i,1) = connCompNodes{i}.NodeValue ; % weight of each node
end
[biggestValue ind] = max(nodeWt);
myBiggestComponent = ind ; % this is the biggest component

%-------------------------------------------------------------------------
%--------CALCULATE THE SHORTEST PATH COST FOR EACH SUCH PATH ----------

pathCost = zeros(n);
%  Get the Path Cost for each edge
for i = 1 : n
    for j = i : n
        Connectivity = leafConnectivity{i,j};
        if ~isempty(Connectivity) && mstf(i,j)       % if an edge is present between the leaves
            Leaf_of_i = Connectivity(1);
            Leaf_of_j = Connectivity(2);
            startNode = connCompNodes{i};
            endNode = connCompNodes{j};
            
            start_x = startNode.LeafData.LeafPosition.y(Leaf_of_i);
            start_y = startNode.LeafData.LeafPosition.x(Leaf_of_i);
            start_z = startNode.LeafData.LeafPosition.z(Leaf_of_i);
            
            end_x = endNode.LeafData.LeafPosition.y(Leaf_of_j);
            end_y = endNode.LeafData.LeafPosition.x(Leaf_of_j);
            end_z = endNode.LeafData.LeafPosition.z(Leaf_of_j);
            
            startVoxel = round([start_x;start_y;start_z]);
            endVoxel = round([end_x;end_y;end_z]);
            
            X = [start_x;end_x]; Y = [start_y;end_y]; Z = [start_z;end_z];
%             offset = 20;
            
            % Testing with the vesselness img
            [subImg ,start , finish , xL,yL,zL] = createSubImg3D(Img,offset,zoffset,round(min(start_x,end_x)),round(min(start_y,end_y)),round(min(start_z,end_z)),...
                round(max(start_x,end_x)),round(max(start_y,end_y)),round(max(start_z,end_z)),X,Y,Z); % create a subimage to compute the graph
            
            start = round(start); finish = round(finish);
%             disp([start' finish'])
            
            ImageGraph = createGraph3D_new(subImg);
            
            % -------------GET THE SHORTEST PATH----------------------
            S = sub2ind(size(subImg),start(1),start(2),start(3));
            T = sub2ind(size(subImg),finish(1),finish(2),finish(3));
            S = round(S);
            T = round(T);
            [pathcost(i,j), shortestpath, predec]=graphshortestpath(ImageGraph,round(S),round(T),'Directed','false');
            
            % Get the path coords in reference to the original image
            [x y z] = ind2sub(size(subImg),shortestpath);  % (x,y,z) corresponds to subimg
            
            if isempty(x) == 0
                x_offset = zeros(1,length(x(:))); x_offset(1,:) = xL ;
                y_offset = zeros(1,length(y(:))); y_offset(1,:) = yL ;
                z_offset = zeros(1,length(z(:))); z_offset(1,:) = zL ;
                
                x_act = x+x_offset; y_act = y+y_offset ; z_act = z + z_offset; % get the actual coordinates as row vector
                actualPath = [x_act; y_act; z_act] ;  % sace the actual path as row vectors of x,y,z
                pathBetweenLeaves{i,j} = actualPath;  % the image path between the leaves
            end
            fprintf('\n Examined Path (%d,%d)',i,j);
            
        end
        
    end
end

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
                fprintf('\n Node %d to be deleted',CC1);
                num_deleted = num_deleted+1;
            end
        else
            if nodeMarker(CC2) == 0
                mstf(CC2,:) = 0 ; % remove all other connections with CC1
                mstf(:,CC2) = 0 ; % remove all other connected to CC1
                nodeMarker(CC2) = 1 ; % indicate that this node is removed
                fprintf('\n Node %d to be deleted',CC2);
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
        
        fprintf('\n The tree has been disconnected, will preserve the largest subtree');
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
        fprintf('\n Preserving the component with weight %d',maxWt);
        node2keep = Component(wtInd).Nodes;
        
        [filtMST,newConnCompNodes,new_leaf_connectivity,new_pathBetweenLeaves] = buildnewTree(node2keep,filtMST,newConnCompNodes,...
            new_leaf_connectivity,new_pathBetweenLeaves);
    end
end
%--------------------------------------------------------------------------
% OLD CODE ......... no longer valid ...........
%  isOrphan = zeros(newLength,1);
%  for i = 1 : newLength
%      rowVect = filtMST(i,:);
%      if any(rowVect) == 0 && newConnCompNodes{i}.NodeValue < biggestValue % if a node is orphan and it is not the biggest CC
%         isOrphan(i) = 1;
%         fprintf('\n Node %d is orphaned,will be deleted',i);
%      end
%  end
%  t = find(isOrphan == 1);
%  latestLength = newLength - length(t(:));
%  [filtMST,newConnCompNodes,new_leaf_connectivity] = deleteNode(latestLength,isOrphan,filtMST,newConnCompNodes,new_leaf_connectivity );
%--------------------------------------------------------------------------
if topK == 0
    filtMST = MST+MST';
    newConnCompNodes = connCompNodes;
    new_leaf_connectivity = leafConnectivity;
    new_pathBetweenLeaves = pathBetweenLeaves;
end

filtMST = sparse(filtMST);
end

