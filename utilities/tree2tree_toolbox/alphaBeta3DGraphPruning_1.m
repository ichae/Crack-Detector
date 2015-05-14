
function [filtMST,newConnCompNodes,...
          new_leaf_connectivity]=...
          alphaBeta3DGraphPruning_1(MST,connCompNodes,...
                                  leaf_connectivity,alpha,beta)

mstf=full(MST);
mstf=mstf+mstf'; % make it symmetric

n=length(connCompNodes); % number of connected components
init_num_edges= n-1;
nzmask = (tril(mstf)>0);
linIndex=find(nzmask);

% if length(linIndex)==init_num_edges,
%     disp('hoorah!');
%     pause;
% end

% find the biggest component: SUVADIP
prev_wt = 0;
for i = 1 : n
     nodeWt(i,1) = connCompNodes{i}.NodeValue ;
   
end
[val ind] = max(nodeWt);
myBiggestComponent = ind ; % this is the biggest component

% -------------------------------------------------------------
    
edgeWeights=mstf(nzmask);
[ew_sorted,sortindex]=sort(edgeWeights,'descend');
ew_sorted_norm=ew_sorted/sum(ew_sorted);

mstf_current=mstf/sum(ew_sorted);
init_mstf=mstf_current;

% 0-not visited, 1-deleted, 2-visited and not deleted
edge_delete_flag=zeros(1,init_num_edges);

initialGraphWeight=0;
for I=1:n,

    initialGraphWeight=initialGraphWeight+connCompNodes{I}.NodeValue;
    
end

initial_edge_wt=sum(ew_sorted_norm);
residual_edge_wt=initial_edge_wt;
% current_edge_index=init_num_edges;

step=1;
% sum(sum(init_mstf-init_mstf'))
% spy(sparse(init_mstf))
% pause

while residual_edge_wt > beta*initial_edge_wt ...
                      && step <= init_num_edges,
                  
    step
    
    if edge_delete_flag(step)==0, % not visited this edge

        [nonemptyI,nonemptyJ]=ind2sub(size(init_mstf),linIndex(sortindex(step)));

        % keep copy of the mstf in case we have to revert back
        mstf_copy=mstf_current;
        
        % disconnect this edge for checking
        mstf_current(nonemptyI,nonemptyJ)=0;
        mstf_current(nonemptyJ,nonemptyI)=0;
       
        % get the isolated trees: S- # CC, C- the nodes
        [S,C] = graphconncomp(sparse(mstf_current));
        
        bigger_weight=-Inf;
        bigger_component=0;
        
        % get the tree with the biggest node weight
        for gr_count=1:S, 
            
            gr_index=find(C==gr_count);
            graphWeight=0;
            
            for node_count=1:length(gr_index),
                
                graphWeight=graphWeight+connCompNodes{gr_index(node_count)}.NodeValue;
                
            end
            
            if graphWeight > bigger_weight,
                
                bigger_weight=graphWeight;
                bigger_component=gr_count;
                
            end
            
        end
        
        % disconnect the smaller trees
         
        mstf_current=init_mstf;
        mstf_current((C~=bigger_component),:)=0;
        mstf_current(:,(C~=bigger_component))=0;
        nzmask1=nzmask;
        nzmask1((C==bigger_component),:)=0;
        linIndex1=find(nzmask1);
%         linIndex1
%         linIndex
        
        edge_del_index=zeros(size(linIndex1));
        
        for I=1:length(linIndex1),
            
            edge_del_index(I)=find(linIndex==linIndex1(I));
            
        end
        
        residual_edge_wt = sum(mstf_current(:))/2;
        
        % check whether edge weight has not fallen below allowable limits
        if (bigger_weight >= alpha*initialGraphWeight) 
            
            disp('success!');

            % delete the edges of the smaller trees
            for I=1:length(edge_del_index),
                
                edge_delete_flag(sortindex==edge_del_index(I))=1;
                
            end
            edge_delete_flag(step)=1;
            
        else
            
            disp('cannot delete !');
            edge_delete_flag(step)=2; % cannot delete,
            mstf_current=mstf_copy;
            residual_edge_wt = sum(mstf_current(:))/2;
            
        end
        
    end
    
%     step
%     edge_del_index
%     edge_delete_flag'
%     pause
    
    step=step+1;
    
end

edge_delete_flag'

new_node_count=0;
new_node_mapping=zeros(1,n);

% [S1,C1]=graphconncomp(sparse(mstf_current))
% sum(sum(mstf_current-mstf_current'))
% spy(sparse(mstf_current))
% sparse(mstf_current)
% pause


% finally, get a new node mapping for the result tree
for I=1:n,
    
    if sum(mstf_current(I,:))>0,
        
        new_node_count=new_node_count+1;
        new_node_mapping(I)=new_node_count;
        
    end
    
end
fprintf('\n Number of remaining edges = %d',new_node_count);

if new_node_count > 0  % added by SUVADIP
    
    new_mstf=zeros(new_node_count);
    new_leaf_connectivity=cell(new_node_count);
    
    % re-number the conn comp connectivity matrix
    for I=1:n,
        
        for J=1:I-1,
            
            if mstf_current(I,J)>0,
                
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
    
else % when only one node remains: SUVADIP
    
    filtMST = [];
    newConnCompNodes{1} = connCompNodes{myBiggestComponent} ;
    new_leaf_connectivity = [];
end
    

