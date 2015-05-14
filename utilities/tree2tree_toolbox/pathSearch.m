function [pathcost, pathBetweenLeaves] = pathSearch(Img, MST, connCompNodes , leafConnectivity , eImg,offset,zoffset)
%GETPATHS Find the shortest paths between the connected components
%   1. Create the subimage and subgraph
%   2. Find the paths between each component
%   3. return the data structure: pathBetweenLeaves and pathcost.
%   pathcost -- gets the cost of the paths between each pair of CC
%   pathBeweenLeaves -- gets the actual path between each pair of CC


% ---------- Data Structures ------------------------------------------

mstf = full(MST);
mstf = mstf+mstf';

n = length(connCompNodes); % number of connected components
init_num_edges = n-1;      % initial number of edges
nzmask = (tril(mstf) > 0);
linIndex=find(nzmask);     % gives the linear indices of the position of edges
pathBetweenLeaves = cell(n) ; % store the shortest path between the leaves


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




end

