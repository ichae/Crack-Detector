
function [splinedMedialGraph,spNodePositions]=...
    graphSpline3D(oldMedialGraph,oldNodePositions,resolution)


resolution=abs(resolution); % check against accidental negative resolution

% [oldNodePositions] = checkUniqueness(oldNodePositions);

% get all leaf and bifurcation nodes
G=full(oldMedialGraph);
node_degree=sum(G,1);
critical_nodelist=find(node_degree~=2);

% This if statement is added by suvadip
if numel(critical_nodelist) == 0
    splinedMedialGraph = oldMedialGraph ;
    spNodePositions = oldNodePositions;
else
    
    S=critical_nodelist(1); % start from a leaf node for DFS 
    % critical_nodelist
    % pause
    
    % start new node numbering after critical nodes
    last_new_node=length(critical_nodelist);
    
    % fill up the list initially with critical nodes
    spNodePositions.x=oldNodePositions.x(critical_nodelist);
    spNodePositions.y=oldNodePositions.y(critical_nodelist);
    spNodePositions.z=oldNodePositions.z(critical_nodelist);
    
    % DFS search of the old graph
    [disc, pred, closed]= graphtraverse(oldMedialGraph, S, 'Directed', 'false');
    
    n=length(disc); % no of nodes in the graph
    
    % Initialize buffers and sentinels
    last_crit_node=S;
    last_crit_node_index=1;
    current_crit_node=[];
    current_crit_node_index=[];
    old_branch_buffer=S;
    branch_length=0;
    
    vertex1_list=[];
    vertex2_list=[];
    
    for I=2:n,
        
        if disc(I-1)== pred(disc(I)), % traversing branch in same path
            
            branch_length(end+1)=branch_length(end) +...
                sqrt((oldNodePositions.x(disc(I))-oldNodePositions.x(disc(I-1))).^2 + ...
                (oldNodePositions.y(disc(I))-oldNodePositions.y(disc(I-1))).^2 + ...
                (oldNodePositions.z(disc(I))-oldNodePositions.z(disc(I-1))).^2);
            
            index=find(critical_nodelist==disc(I)); % check to see whether critical node
            
            if  ~isempty(index), % found a critical node
                
                old_branch_buffer(end+1)=disc(I);
                current_crit_node=disc(I);
                current_crit_node_index=index;
                
                % special care for just a one node branch
                
                if resolution > branch_length(end),
                    
                    new_index=0;
                    
                else
                    
                    new_index=(0:resolution:branch_length(end));
                    
                end
                
                % always make sure the start and end points are critical
                % nodes
                
                if new_index(end)~=branch_length(end),
                    
                    new_index(end+1)=branch_length(end);
                    
                end
                
                % get the splined nodes for this branch with uniform
                % resolution
                
                branch_x=spline(branch_length,oldNodePositions.x(old_branch_buffer),new_index);
                branch_y=spline(branch_length,oldNodePositions.y(old_branch_buffer),new_index);
                branch_z=spline(branch_length,oldNodePositions.z(old_branch_buffer),new_index);
                
                new_branch_length=length(branch_x);
                
                branch_nodelist=[last_crit_node_index ...
                    ((last_new_node+1):1:(last_new_node+new_branch_length-2))...
                    current_crit_node_index];
                
                %                 branch_nodelist
                %                 scatter(branch_x,branch_y,'y','filled');pause;
                
                for J=1:new_branch_length-1, % connect the nodes in the branch
                    
                    vertex1_list(end+1)=branch_nodelist(J);
                    vertex2_list(end+1)=branch_nodelist(J+1);
                    vertex1_list(end+1)=branch_nodelist(J+1);
                    vertex2_list(end+1)=branch_nodelist(J);
                    
                end
                
                % add the new node coordinates to spNodePositions
                spNodePositions.x=[spNodePositions.x branch_x(2:end-1)];
                spNodePositions.y=[spNodePositions.y branch_y(2:end-1)];
                spNodePositions.z=[spNodePositions.z branch_z(2:end-1)];
                
                % reset all sentinels and buffers
                
                last_crit_node=current_crit_node; % reset the last critical node
                last_crit_node_index=current_crit_node_index; % reset the index as well
                old_branch_buffer=last_crit_node; % start a new branch buffer
                branch_length=0; % a new branch begins again
                last_new_node=last_new_node+new_branch_length-2; % prepare for next branch
                
            else % keep on adding to branch buffer
                
                old_branch_buffer(end+1)=disc(I);
                
            end
            
        else % changed over to branch in another path
            
            % first reset the last critical index and branch buffer
            
            last_crit_node=pred(disc(I));
            last_crit_node_index=find(critical_nodelist==last_crit_node);
            old_branch_buffer=last_crit_node; % start a new branch buffer
            branch_length=0; % a new branch begins again
            
            % now start measuring branch length again
            
            branch_length(end+1)=branch_length(end) +...
                sqrt((oldNodePositions.x(disc(I))-oldNodePositions.x(last_crit_node)).^2 + ...
                (oldNodePositions.y(disc(I))-oldNodePositions.y(last_crit_node)).^2 +...
                (oldNodePositions.z(disc(I))-oldNodePositions.z(last_crit_node)).^2);
            
            
            index=find(critical_nodelist==disc(I)); % check to see whether critical node
            
            
            if  ~isempty(index), % found a critical node
                
                old_branch_buffer(end+1)=disc(I);
                current_crit_node=disc(I);
                current_crit_node_index=index;
                
                % special care for just a one node branch
                
                if resolution>branch_length(end),
                    
                    new_index=0;
                    
                else
                    
                    new_index=(0:resolution:branch_length(end));
                    
                end
                
                % always make sure the start and end points are critical
                % nodes
                
                if new_index(end)~=branch_length(end),
                    
                    new_index(end+1)=branch_length(end);
                    
                end
                
                % get the splined nodes for this branch with uniform
                % resolution
                
                branch_x=spline(branch_length,oldNodePositions.x(old_branch_buffer),new_index);
                branch_y=spline(branch_length,oldNodePositions.y(old_branch_buffer),new_index);
                branch_z=spline(branch_length,oldNodePositions.z(old_branch_buffer),new_index);
                %                 scatter(branch_x,branch_y,'y','filled');pause;
                
                new_branch_length=length(branch_x);
                
                branch_nodelist=[last_crit_node_index ...
                    ((last_new_node+1):1:(last_new_node+new_branch_length-2))...
                    current_crit_node_index];
                
                %                 branch_nodelist
                %                 scatter(branch_x,branch_y,'y','filled');pause;
                
                for J=1:new_branch_length-1, % connect the nodes in the branch
                    
                    vertex1_list(end+1)=branch_nodelist(J);
                    vertex2_list(end+1)=branch_nodelist(J+1);
                    vertex1_list(end+1)=branch_nodelist(J+1);
                    vertex2_list(end+1)=branch_nodelist(J);
                    
                end
                
                % add the new node coordinates to spNodePositions
                spNodePositions.x=[spNodePositions.x branch_x(2:end-1)];
                spNodePositions.y=[spNodePositions.y branch_y(2:end-1)];
                spNodePositions.z=[spNodePositions.z branch_z(2:end-1)];
                
                % reset all sentinels and buffers
                
                last_crit_node=current_crit_node; % reset the last critical node
                last_crit_node_index=current_crit_node_index; % reset the index as well
                old_branch_buffer=last_crit_node; % start a new branch buffer
                branch_length=0; % a new branch begins again
                last_new_node=last_new_node+new_branch_length-2; % prepare for next branch
                
            else % keep on adding to branch buffer
                
                old_branch_buffer(end+1)=disc(I);
                
            end
            
        end
        
    end
    
    
    edges=ones(1,length(vertex1_list));
    
    splinedMedialGraph=sparse(vertex1_list,vertex2_list,edges);
    
    if n==1,
        
        splinedMedialGraph=1;
        spNodePositions.x=oldNodePositions.x;
        spNodePositions.y=oldNodePositions.y;
        spNodePositions.y=oldNodePositions.z;
        
    end
end

