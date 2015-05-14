
function [prunedMedialGraph,prunedNodePositions]=...
                    removeUnlikelyBranches(splinedMedialGraph,...
                    spNodePositions,cl_dist,grayImage,bwImage,resolution)
        
% get all leaf and bifurcation nodes

G=full(splinedMedialGraph);
node_degree=sum(G,1);
critical_nodelist=find(node_degree~=2);
S=critical_nodelist(1); % start from a leaf node for DFS search

% create an image that is just a tubular neighborhood 
% of the splined graph

nbhd_image=zeros(size(grayImage));
graph_x=max(1,floor(spNodePositions.x));
graph_y=max(1,floor(spNodePositions.y));
ind=sub2ind(size(grayImage),graph_y,graph_x);
nbhd_image(ind)=1;
R=3; % radius of the tubular neighborhood
se=strel('disk',R,0);
nbhd_image=imdilate(nbhd_image,se);

% create a adjacency matrix for inter-critical_node 
% conenctivity (or the branches)

n_crit_nodes=length(critical_nodelist);
branchCheck=zeros(n_crit_nodes);

% DFS search of the old graph
[disc, pred, closed]= graphtraverse(splinedMedialGraph,...
                      S, 'Directed', 'false');

n=length(disc); % no of nodes in the graph

% Initialize buffers and sentinels
last_crit_node=S;
last_crit_node_index=1;
current_crit_node=[];
current_crit_node_index=[];
branch_buffer=S;
branch_strength=0;

%% Get a validity score for each branch in the splined graph

for I=2:n,
    
%     I
    
    if disc(I-1)==pred(disc(I)), % traversing branch in same path
        
        
        index=find(critical_nodelist==disc(I)); % check to see whether critical node
        
        if  ~isempty(index), % found a critical node
            
            branch_buffer(end+1)=disc(I);
            current_crit_node=disc(I);
            current_crit_node_index=index;
            
            % get the node positions for this branch
            
            branch_x=graph_x(branch_buffer);
            branch_y=graph_y(branch_buffer);
            branchCheck(current_crit_node_index,last_crit_node_index)=...
                getBranchScore(branch_buffer,branch_x,branch_y,...
                grayImage,bwImage,nbhd_image,resolution);
%             
%             current_crit_node_index
%             last_crit_node_index
%             branchCheck(current_crit_node_index,last_crit_node_index)
% %             'boom'
%             pause
            
            % reset all sentinels and buffers
            
            last_crit_node=current_crit_node; % reset the last critical node
            last_crit_node_index=current_crit_node_index; % reset the index as well
            branch_buffer=last_crit_node; % start a new branch buffer
            
        else % keep on adding to branch buffer
            
            branch_buffer(end+1)=disc(I);
            
        end
        
    else % changed over to branch in another path
        
        % first reset the last critical index and branch buffer
        
        last_crit_node=pred(disc(I));
        last_crit_node_index=find(critical_nodelist==last_crit_node);
        branch_buffer=last_crit_node; % start a new branch buffer
        
        index=find(critical_nodelist==disc(I)); % check to see whether critical node
        
        
        if  ~isempty(index), % found a critical node
            
            branch_buffer(end+1)=disc(I);
            current_crit_node=disc(I);
            current_crit_node_index=index;
            
            % get the node positions for this branch
            
            branch_x=graph_x(branch_buffer);
            branch_y=graph_y(branch_buffer);
            branchCheck(current_crit_node_index,last_crit_node_index)=...
                getBranchScore(branch_buffer,branch_x,branch_y,...
                grayImage,bwImage,nbhd_image,resolution);
            
%             'boom'
            
            % reset all sentinels and buffers
            
            last_crit_node=current_crit_node; % reset the last critical node
            last_crit_node_index=current_crit_node_index; % reset the index as well
            branch_buffer=last_crit_node; % start a new branch buffer
            
        else % keep on adding to branch buffer
            
            branch_buffer(end+1)=disc(I);
            
        end
        
    end
    
end

branchCheck1=branchCheck+branchCheck';
branchCheck=branchCheck1;

branchCheck

%% cluster critical nodes based on given clustering distance

node_flag=zeros(1,n_crit_nodes);
node_flag(1)=1; % start numbering from the 1st node
new_node_count=1; % first node already accounted for

for I=2:n_crit_nodes, 
    
    for J=1:I-1,
    
        if branchCheck(I,J)>0, % crit nodes I and J are connected
            
            if branchCheck(I,J) < cl_dist, % have to be merged
                
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
 
node_flag

%% delete branches which are not within cl_dist-path

biggest_cluster=1;
biggest_cluster_size=0;

clusters=unique(node_flag);

% find the biggest cluster

for cluster_no=1:length(clusters); 
    
    cluster_size=sum(node_flag==clusters(cluster_no));
    
    if cluster_size>biggest_cluster_size,
        
        biggest_cluster_size=cluster_size;
        biggest_cluster=clusters(cluster_no);
        
    end
    
end

biggest_cluster

% now delete edges which came from a branch with unacceptable validity
% score

% Initialize buffers and sentinels

last_crit_node=S;
last_crit_node_index=1;
current_crit_node=[];
current_crit_node_index=[];
branch_buffer=S;
branch_strength=0;
vertex1=[];
vertex2=[];

for I=2:n,
    
        if disc(I-1)==pred(disc(I)), % traversing branch in same path
           
                    
            index=find(critical_nodelist==disc(I)); % check to see whether critical node         
            
            if  ~isempty(index), % found a critical node
                
                branch_buffer(end+1)=disc(I);
                current_crit_node=disc(I);
                current_crit_node_index=index;
                
                % both endpoints belong to the biggest cluster
                
                if node_flag(current_crit_node_index)==biggest_cluster && ...
                   node_flag(last_crit_node_index)==biggest_cluster,
               
                    % put this branch into the edge connectivity list
                    
                    vertex1=[vertex1 branch_buffer(1:end-1)];
                    vertex2=[vertex2 branch_buffer(2:end)];
                    
                end
                
                % reset all sentinels and buffers
                
                last_crit_node=current_crit_node; % reset the last critical node
                last_crit_node_index=current_crit_node_index; % reset the index as well
                branch_buffer=last_crit_node; % start a new branch buffer
                               
            else % keep on adding to branch buffer
                
                branch_buffer(end+1)=disc(I);
                
            end
               
       else % changed over to branch in another path
           
           % first reset the last critical index and branch buffer
           
           last_crit_node=pred(disc(I));
           last_crit_node_index=find(critical_nodelist==last_crit_node);
           branch_buffer=last_crit_node; % start a new branch buffer
   
           index=find(critical_nodelist==disc(I)); % check to see whether critical node
         
           if  ~isempty(index), % found a critical node
    
                branch_buffer(end+1)=disc(I);
                current_crit_node=disc(I);
                current_crit_node_index=index;
                     
                % both endpoints belong to the biggest cluster
                
                if node_flag(current_crit_node_index)==biggest_cluster && ...
                   node_flag(last_crit_node_index)==biggest_cluster,
               
                    % put this branch into the edge connectivity list
                    
                    vertex1=[vertex1 branch_buffer(1:end-1)];
                    vertex2=[vertex2 branch_buffer(2:end)];
                    
                end
            
                % reset all sentinels and buffers
                
                last_crit_node=current_crit_node; % reset the last critical node
                last_crit_node_index=current_crit_node_index; % reset the index as well
                branch_buffer=last_crit_node; % start a new branch buffer
                               
            else % keep on adding to branch buffer
                
                branch_buffer(end+1)=disc(I);
                
            end
                       
        end
           
end

vertex1_list=[vertex1 vertex2]; % undirected graph

[vertex1_uniq,ii1,jj1]=unique(vertex1_list); % get the unique nodes and correspondnig numbering in jj1

shiftamount=length(vertex1);

jj2=circshift(jj1,[0 shiftamount]); % create the ending vertex list

edge_weights=ones(1,length(vertex1_list));

if isempty(vertex1), % the medial graph had only one node
    
    prunedMedialGraph=0; % G=[0] in this case 
    prunedNodePositions.x=spNodePositions.x;
    prunedNodePositions.y=spNodePositions.y;
    
else % otherwise treat as a general case
    
    prunedMedialGraph=sparse(jj1,jj2,edge_weights);
    prunedNodePositions=struct('x',[],'y',[]);
    prunedNodePositions.x=spNodePositions.x(vertex1_uniq);
    prunedNodePositions.y=spNodePositions.y(vertex1_uniq);

end






 