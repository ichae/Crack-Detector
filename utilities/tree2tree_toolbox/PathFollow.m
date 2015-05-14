

function [splinedMedialGraph,spNodePositions,startnode,endnode]=PathFollow(MST,connCompNodes,connectivity,pathBetweenLeaves,resolution)

% HERE THE SHORTEST PATH IS STICHED BETWEEN TWO CONNECTED COMPONENTS


n=length(connCompNodes); % no of conn comps (nodes) in MST

nodeOffset=zeros(1,n+1);% n+1th position will never be used
newNodePositions.x=[];
newNodePositions.y=[];
newNodePositions.z=[];
startnode=[];
endnode=[];

path = [];
% stack the connections inside each individual conn comp
for I=1:n,
    
    nodeFull=full(connCompNodes{I}.MedialGraph);
    
    m=length(connCompNodes{I}.NodePositions.x); % no of nodes in connCompNodes(I)
    
    newNodePositions.x=[newNodePositions.x connCompNodes{I}.NodePositions.x];
    newNodePositions.y=[newNodePositions.y connCompNodes{I}.NodePositions.y];
    newNodePositions.z=[newNodePositions.z connCompNodes{I}.NodePositions.z];
    
    for J=1:m,
        for K=1:J-1,
            if nodeFull(J,K)==1,
                startnode(end+1)=nodeOffset(I)+J;
                endnode(end+1)=nodeOffset(I)+K;
                startnode(end+1)=nodeOffset(I)+K;
                endnode(end+1)=nodeOffset(I)+J;
            end
        end
    end
    
    nodeOffset(I+1)=nodeOffset(I)+m; % new nodeOffset in the next iteration
    fprintf('\n Connected points inside a CC(%d)',I);
    
end

% Now take care of the connectivity between the conn comps

nodeID = max(max(startnode(:)),max(endnode(:)))+1; % new node number to start numbering

if n > 1
    MST_full=full(MST)+full(MST)';

    for I=1:n,
        
        for J=I+1:n,
                 
            if MST_full(I,J)>0, % edge exists between I and J
                 
                leaf_I_index=connectivity{I,J}(1); % leaf numbers for the connection
                leaf_J_index=connectivity{I,J}(2);
                
                leaf_I=connCompNodes{I}.LeafData.LeafNodes(leaf_I_index);
                leaf_J=connCompNodes{J}.LeafData.LeafNodes(leaf_J_index);
                
                startLeaf = leaf_I + nodeOffset(I);
                endLeaf = leaf_J + nodeOffset(J);
                % We need to stich a path between startLeaf and endLeaf
                
                prunedpath = pathBetweenLeaves{I,J};
                % deformed_path = pathDeform(prunedpath);  %to be
                % implemented
                
                if isempty(prunedpath)    % connect these two leaves directly in case of no path
                    startnode(end+1)=nodeOffset(I)+leaf_I;
                    endnode(end+1)=nodeOffset(J)+leaf_J;
                    startnode(end+1)=nodeOffset(J)+leaf_J;
                    endnode(end+1)=nodeOffset(I)+leaf_I;
                    fprintf('\n Force joining the nodes...no prunedpath');
                else
                    xval = prunedpath(1,:);
                    yval = prunedpath(2,:);
                    zval = prunedpath(3,:);
                    %                     if isempty(xval)==0 && isempty(yval)==0 && isempty(zval)==0
                    fprintf('\n inserting a path');
                    % Use saurav's convention
                    skip_pts = ceil((length(xval)-1)/10); % CHANGE TO 1 IN CASE OF BUGS
                    
                    % BUGFIX........
                    % ----- start from yval(2:skip)...changed from (1:skip)
                    newNodePositions.x=[newNodePositions.x yval(2:skip_pts:end)];
                    newNodePositions.y=[newNodePositions.y xval(2:skip_pts:end)];
                    newNodePositions.z=[newNodePositions.z zval(2:skip_pts:end)];
                                        
                    startnode(end+1) = startLeaf;
                    endnode(end+1) = nodeID;   % connect the leaf and the start of path
                    startnode(end+1) = nodeID;
                    endnode(end+1) = startLeaf;   % connect the leaf and the start of path
                    fprintf('\n Start Stich between (%d,%d)',startLeaf,nodeID)
                    
                    % ................BUGFIX.........................
                    % ----- start from yval(2:skip)...changed from (1:skip)
                    for k = 1 : length(xval(2:skip_pts:end))-1
                        startnode(end+1)=nodeID;
                        endnode(end+1)=nodeID+1;
                        startnode(end+1)=nodeID+1;
                        endnode(end+1)=nodeID;
                        fprintf('\n Stich between (%d,%d)',nodeID,nodeID+1);
                        nodeID = nodeID+1;
                    end
                    % stich the end leaf
                    startnode(end+1)=nodeID;
                    endnode(end+1)=endLeaf;
                    startnode(end+1)=endLeaf;
                    endnode(end+1)=nodeID;
                    fprintf('\n End Stich between (%d,%d)',nodeID,endLeaf);
                    nodeID = nodeID+1;
                end
            end
        end
    end
end

edgeweight=ones(1,length(startnode));

newMedialGraph=sparse(startnode,endnode,edgeweight);

splinedMedialGraph = newMedialGraph;
spNodePositions = newNodePositions;

% [splinedMedialGraph,spNodePositions]=graphSpline3D(newMedialGraph,newNodePositions,resolution);




