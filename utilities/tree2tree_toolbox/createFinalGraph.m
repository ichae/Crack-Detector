function [finalGraph,NodePositions,startNode,endNode] = createFinalGraph(connCompNodes,connectivity,pathBetweenCC)
%CREATEFINALGRAPH


num_cc = length(connCompNodes); % no of conn comps (nodes) in MST

% Intra component connection
leafStack = [];
startNode = [];
endNode = [];
nodePos_X = [];
nodePos_Y = [];
nodePos_Z = [];

node_ind = 1;
CC_done = zeros(1,num_cc);

nodeStack = [];
predecStack = [];

for i = 1 : num_cc
    
    currentCC = connCompNodes{i};
    CC_graph = full(currentCC.MedialGraph);  % get the connected comp graph
    num_pts = length(currentCC.NodePositions.x); % this connected comp. has num_pts coordinates
    
%     child_connected = connectivity(i,:);
    
    for j = 1 : num_pts % fill up each node value and its parent-son relation
        startNode = [startNode ; node_ind];
        endNode = [endNode;node_ind+1];
        
        startNode = [startNode ; node_ind+1];
        endNode = [endNode;node_ind];
        node_ind = node_ind+1 ;
        xval = currentCC.NodePositions.x(j);
        yval = currentCC.NodePositions.y(j);
        zval = currentCC.NodePositions.z(j);
        nodePos_X = [nodePos_X ; xval];
        nodePos_Y = [nodePos_Y ; yval];
        nodePos_Z = [nodePos_Z ; zval];
    end
    CC_done(i) = 1 ;    % completed this current CC
    
    % Now that we have recorded the nodes for a single connected component,
    % look for its connections to other CC
    
    % push in stack the children that have not been connected
    for  j = i : num_cc
        child = connectivity(i,j);
        if ~isempty(child) && CC_done(j) == 0  % if such a connection exists and it has not been examined
            nodeStack = [nodeStack;j]; % insert that node in the queue
            predecStack = [predecStack;i];  % stack the id of its predec
        end
    end
    
    %     thisChild = nodeStack(1);   % get the first
    
    while isempty(nodeStack) == 0 % do while all the children are visited
        % Bridge the gap between the child and parent
        thisChild = nodeStack(1);
        predec = predecStack(1);
        nodeStack(1) = [];  % pop
        predecStack(1) = [];
        
        parentID = predec;
        child =  thisChild;
        path = pathBetweenCC{i,j};
        path_x = path(1,:);
        path_y = path(2,:);
        path_z = path(3,:);
        
        childID = max(max(endNode(:)),max(startNode(:)));
        node_ind = childID;
        
        % Connect the parent and child
        startNode = [startNode;node_ind+1];
        endNode = [endNode;parentID];
        startNode = [startNode;parentID];
        endNode = [endNode;node_ind+1];
        % now connect the points inside the path
        for  j = 1 : length(path_x(:))
            startNode = [startNode ; node_ind];
            endNode = [endNode;node_ind+1];
            node_ind = node_ind+1 ;
            nodePos_X = [nodePos_X ; path_x(j)];
            nodePos_Y = [nodePos_Y ; path_y(j)];
            nodePos_Z = [nodePos_Z ; path_z(j)];
        end
        
        % now take care of this child (i.e inter CC connection)
        childCC = connCompNodes{child};  % get the child cc
        child_num_pts = length(childCC.NodePositions.x); % this connected comp. has child_num_pts coordinates
        
        for j = 1 : child_num_pts % fill up each node value and its parent-son relation
            startNode = [startNode ; node_ind];
            endNode = [endNode;node_ind+1];
            node_ind = node_ind+1 ;
            xval = childCC.NodePositions.x(j);
            yval = childCC.NodePositions.y(j);
            zval = childCC.NodePositions.z(j);
            nodePos_X = [nodePos_X ; xval];
            nodePos_Y = [nodePos_Y ; yval];
            nodePos_Z = [nodePos_Z ; zval];
        end
        
        CC_done(thisChild) = 1 ;    % completed this current child
        
        
        % push in stack the children that have not been connected
        for  j = thisChild : num_cc
%               this_child_connected = connectivity{thisChild,:};
%             child = child_connected{j};
            if ~isempty(connectivity(thisChild,j)) && CC_done(j) == 0  % if such a connection exists and it has not been examined
                nodeStack = [nodeStack;j]; % insert that node in the stack
                predecStack = [predecStack;childID]; % insert the parent
            end
        end
        
    end
    
    
    
end

% disp(length(startNode));
% disp(length(endNode));
m = max(max(startNode(:)),max(endNode(:)));
array = zeros(m);
for i = 1 : length(startNode)
   x = startNode(i);
   y = endNode(i);
   array(x,y)=1;
end
% array=array+array';

finalGraph = sparse(array);

% edgeweight = ones(1,length(startNode));
% finalGraph = sparse(startNode,endNode,edgeweight);

% 
NodePositions.x = [];
NodePositions.y = [];
NodePositions.z = [];

% NodePositions = struct('x',{},'y',{},'z',{});

NodePositions.x = nodePos_X;
NodePositions.y = nodePos_Y;
NodePositions.z = nodePos_Z;


end

