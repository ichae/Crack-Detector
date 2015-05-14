
function [filtMST,newConnCompNodes,new_leaf_connectivity]=alphaBetaFilterGraph(...
                                                          MST,connCompNodes,...
                                                          leaf_connectivity,alpha,beta)

mstf=full(MST);
mstf=mstf+mstf'; % make it symmetric
initialGraphWeight=0;
n=length(connCompNodes);
init_num_edges=n-1;
edgeWeights=mstf(tril(mstf)>0);
[ew_sorted,sortindex]=sort(edgeWeights,'ascend');
ew_sorted_norm=ew_sorted/sum(ew_sorted);
mstf_current=mstf/sum(ew_sorted);
check_matrix=(mstf_current>0);

% 0-not visited, 1-deleted, 2-visited and not deleted
edge_delete_flag=zeros(1,init_num_edges);

for I=1:n,

    initialGraphWeight=initialGraphWeight+connCompNodes{I}.NodeValue;
    
end
    
stop='no';
current_edge_index=init_num_edges;
sum_of_deleted_edges=0;
currentLargestGraphWeight=initialGraphWeight;

sauravcheck=0

while strcmp(stop,'no') && (current_edge_index > 0),
    
    sauravcheck=sauravcheck+1
    current_edge_weight=ew_sorted_norm(current_edge_index);
    
%     current_edge_weight
        
    if (current_edge_weight+sum_of_deleted_edges) >= alpha, % within noise removal limit
        
        % find the edges with the given weight that has not been visited
        
        checkindex=(tril(mstf_current)==current_edge_weight...
                                   & tril(check_matrix)~=2);
        [nonemptyI,nonemptyJ]=find(checkindex==1);
        nonemptyI=nonemptyI(1);
        nonemptyJ=nonemptyJ(1);
        
%         nonemptyI
%         nonemptyJ
        
        % do all checks with a copy of the mstf
        mstf_copy=mstf_current;
        
        % disconnect this edge for checking
        mstf_current(nonemptyI,nonemptyJ)=0;
        mstf_current(nonemptyJ,nonemptyI)=0;
        
        % get the isolated trees
        [S,C] = graphconncomp(sparse(mstf_current));
        
%         S
%         C
        
        biggest_weight=0;
        biggest_component=0;
        
        % get the tree with the biggest node weight
        for gr_count=1:S, 
            
            gr_index=find(C==gr_count);
            graphWeight=0;
            
            for node_count=1:length(gr_index),
                
                graphWeight=graphWeight+connCompNodes{gr_index(node_count)}.NodeValue;
                
            end
            
            if graphWeight>biggest_weight,
                
                biggest_weight=graphWeight;
                biggest_component=gr_count;
                
            end
            
        end
        
        if (currentLargestGraphWeight-biggest_weight) ...
                >= beta*initialGraphWeight, % allowable loss of conn comps
            
            currentLargestGraphWeight=currentLargestGraphWeight-biggest_weight;
            edge_delete_flag(current_edge_index)=1; % deleted
            current_edge_index=current_edge_index-1; % move to next edge
            sum_of_deleted_edges=sum_of_deleted_edges+...
                                 current_edge_weight; 
            
            % indicate that this edge has been visited 
            check_matrix(sub2ind(size(check_matrix),nonemptyI,nonemptyJ))=2;
            check_matrix(sub2ind(size(check_matrix),nonemptyJ,nonemptyI))=2;
            
            % disconnect all the remaining nodes
            other_nodes=find(C~=biggest_component);
                        
            for other_node_count=1:length(other_nodes),
                
                mstf_copy(other_nodes(other_node_count),:)=0;
                
            end
            
            % copy back, since mstf has changed
            mstf_current=mstf_copy;
            
        elseif ew_sorted_norm(current_edge_index)...
                ==ew_sorted_norm(current_edge_index-1), % 2 edges of same length
            
            edge_delete_flag(current_edge_index)=2; % visited but not deleted
            current_edge_index=current_edge_index-1;
            
            % indicate that this edge has been visited
            check_matrix(sub2ind(size(check_matrix),nonemptyI,nonemptyJ))=2;
            check_matrix(sub2ind(size(check_matrix),nonemptyJ,nonemptyI))=2;
            
        else
            
            stop='yes'; % we can stop the deletion now
            
        end
        
    else
        
        % stop deletion because edge weight has fallen below specified
        % value
        stop='yes'; 
        
    end
    
end

new_node_count=0;
new_node_mapping=zeros(1,n);

% finally, get a new node mapping for the result tree
for I=1:n,
    
    if sum(mstf_current(I,:))>0,
        
        new_node_count=new_node_count+1;
        new_node_mapping(I)=new_node_count;
        
    end
    
end

new_mstf=zeros(new_node_count);
new_leaf_connectivity=cell(new_node_count);

% re-number the conn comp connectivity matrix
for I=1:n,
    
    for J=1:I-1,
        
        if mstf_current(I,J)>0,
%             
%             I
%             J
%             new_node_mapping(I)
%             new_node_mapping(J)
%             pause
            
            new_mstf(new_node_mapping(I),new_node_mapping(J))=mstf_current(I,J);
            new_mstf(new_node_mapping(J),new_node_mapping(I))=mstf_current(I,J);
            new_leaf_connectivity(new_node_mapping(I),new_node_mapping(J))=...
                                                         leaf_connectivity(I,J);
            new_leaf_connectivity(new_node_mapping(J),new_node_mapping(I))=...
                                                         leaf_connectivity(J,I);
        end
        
    end
    
end

newConnCompNodes=[];

% renumber the remaining connected components
for I=1:new_node_count,
    
    newConnCompNodes{end+1}=connCompNodes{new_node_mapping==I};
    
end

filtMST=sparse(new_mstf);
                
    

