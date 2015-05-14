
function [newMedialGraph,newNodePositions,indexMap]=...
           mergeAndSmoothNodes(oldMedialGraph,oldNodePositions,resolution)

OMG_full=full(oldMedialGraph);

no_of_nodes_old=length(oldNodePositions.x);

node_flag=zeros(1,no_of_nodes_old);
node_flag(1)=1; % start numbering from the 1st node

new_node_count=1; % first node already accounted for

% oldNodePositions.x
% length(oldNodePositions.x)
% 'no_of_old_nodes=',no_of_nodes_old

for I=2:no_of_nodes_old, 
    
    for J=1:I-1,
    
        if OMG_full(I,J)==1, % I and J are connected
            
            dist_I_to_J=sqrt((oldNodePositions.x(I)-oldNodePositions.x(J)).^2 +...
                (oldNodePositions.y(I)-oldNodePositions.y(J)).^2);
                    
            if dist_I_to_J<resolution, % have to be merged
                
                if node_flag(I)==0 && node_flag(J)==0, % none of the nodes have been merged
                    
                    new_node_count=new_node_count+1;
                    node_flag(I)=new_node_count;
                    node_flag(J)=node_flag(I);
                    
                elseif node_flag(I)~=0 && node_flag(J)==0 % startnode has been merged
                    
                    node_flag(J)=node_flag(I);
                    
                elseif node_flag(I)==0 && node_flag(J)~=0 % endnode has been merged
                    
                    node_flag(I)=node_flag(J);
                    
                else % both nodes have been merged 
                    
                    where_is_J=find(node_flag==node_flag(J));
                    node_flag(where_is_J)=node_flag(I); % set a common flag
                    
                end
             
            else % not to be merged
                
                if node_flag(I)==0, % assign a new flag if it has none
                    
                    new_node_count=new_node_count+1;
                    node_flag(I)=new_node_count;
                    
                end
                
                if node_flag(J)==0, % assign a new flag if it has none
                    
                    new_node_count=new_node_count+1;
                    node_flag(J)=new_node_count;
                    
                end 
                    
            end
                        
        end
        
    end
    
end

% get a unique node numbering for the repeated values in
% node_flag

[unique_nodes,ii,jj]=unique(node_flag);

n=length(unique_nodes);
% n

% get the new node positions of the merged nodes

newNodePositions=struct('x',zeros(1,n),'y',zeros(1,n));

% newNodePositions
% oldNodePositions.x
% oldNodePositions.y

for I=1:n,
    
    % collect nodes with node_flag from unique_nodes
    
    where_is_I=find(node_flag==unique_nodes(I)); 
    
%     I
%     
%     where_is_I
%     oldNodePositions.x(where_is_I)
    
    % set the new node position to mean of the merged nodes
    
    newNodePositions.x(I)=sum(oldNodePositions.x(where_is_I))/length(where_is_I);
    newNodePositions.y(I)=sum(oldNodePositions.y(where_is_I))/length(where_is_I);
    
end
    
indexMap=zeros(1,no_of_nodes_old);

% create new adjacency matrix

NMG_full=zeros(n);

for I=2:no_of_nodes_old, 
    
    for J=1:I-1,
        
        % check if I and J are connected and have different labels/flags
        
        if OMG_full(I,J)==1 && (node_flag(I)~=node_flag(J)),
            
            I_mapped_to_this=find(unique_nodes==node_flag(I));
            J_mapped_to_this=find(unique_nodes==node_flag(J));
            
            NMG_full(I_mapped_to_this(1),J_mapped_to_this(1))=1;
            NMG_full(J_mapped_to_this(1),I_mapped_to_this(1))=1; % symmetric
            indexMap(I)=I_mapped_to_this(1);
            indexMap(J)=J_mapped_to_this(1);
            
        end
        
    end
    
end

newMedialGraph=sparse(NMG_full);


    
    

