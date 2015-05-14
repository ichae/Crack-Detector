
function [splinedMedialGraph,spNodePositions]=fitSplineToMST3D(MST,connCompNodes,connectivity,resolution)

n=length(connCompNodes); % no of conn comps (nodes) in MST

nodeOffset=zeros(1,n+1);% n+1th position will never be used
newNodePositions.x=[];
newNodePositions.y=[];
newNodePositions.z=[];
startnode=[];
endnode=[];


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
        
    end
    
    % Now take care of the connectivity between the conn comps
    
    nodeID = max(max(startnode(:)),max(endnode(:)));
    
if n > 1    
    MST_full=full(MST);
    
    for I=1:n,
        
        for J=1:I-1, 
            
            if MST_full(I,J)>0, % edge exists between I and J
                
                leaf_I_index=connectivity{I,J}(1); % leaf numbers for the connection
                leaf_J_index=connectivity{I,J}(2);
                
                leaf_I=connCompNodes{I}.LeafData.LeafNodes(leaf_I_index);
                leaf_J=connCompNodes{J}.LeafData.LeafNodes(leaf_J_index);
                
%                 startLeaf = leaf_I + nodeOffset(I);
%                 endLeaf = leaf_J + nodeOffset(J);
                
                % now place the path points between them
                
                
                
                startnode(end+1)=nodeOffset(I)+leaf_I;
                endnode(end+1)=nodeOffset(J)+leaf_J;
                startnode(end+1)=nodeOffset(J)+leaf_J;
                endnode(end+1)=nodeOffset(I)+leaf_I;
                
            end
            
        end
        
    end
%else    % n == 1, one connected component
    
    
    
end

edgeweight=ones(1,length(startnode));

newMedialGraph=sparse(startnode,endnode,edgeweight);

[splinedMedialGraph,spNodePositions]=graphSpline3D(newMedialGraph,newNodePositions,resolution);




