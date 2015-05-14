
function [uMedialGraph,uNodePositions]=uniformlyResolveGraph(oldMedialGraph,oldNodePositions,resolution)

%% %%%%%%%%%%%%%%%%%% critical node list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

node_degree=sum(full(oldMedialGraph),1);
critical_nodelist=find(node_degree~=2);

%% %%%%% create uniformly sampled graph %%%%%%%

S=critical_nodelist(1); % start node for DFS search

[disc,pred,closed]= graphtraverse(oldMedialGraph,S,'Directed','false');

startnode=S;
% last_visited_crit_node=S;
% endnode=0;
dist_from_last_crit_node=0;

vertex1_list=[];
vertex2_list=[];

for I=2:length(disc),
    
        if disc(I-1)==pred(disc(I)), % traversing same branch
           
            index=find(critical_nodelist==disc(I)); % check to see whether critical node
            
            if  ~isempty(index), % found a critical node
                
                dist_from_last_crit_node=0;
                last_visited_crit_node=critical_nodelist(index);
                endnode=last_visited_crit_node;
 
                vertex1_list(end+1)=startnode;
                vertex2_list(end+1)=endnode;
                
                startnode=endnode; % reset the starting node
                
            elseif mod(dist_from_last_crit_node+1,resolution)==0 % node at sampling resolution
                
                dist_from_last_crit_node=dist_from_last_crit_node+1;
                endnode=disc(I);
                
                vertex1_list(end+1)=startnode;
                vertex2_list(end+1)=endnode;

                startnode=endnode; % reset the starting node
                
            else % unimportant node, just skip over
                
                dist_from_last_crit_node=dist_from_last_crit_node+1;
                
            end
               
       else % changed over to another branch
           
           gobackindex=find(disc==pred(disc(I)));
           last_visited_crit_node=disc(gobackindex);
           dist_from_last_crit_node=0;
           startnode=last_visited_crit_node;
           
%            gobackindex
           
           index=find(critical_nodelist==disc(I)); % check to see whether critical node

            
           if  ~isempty(index), % found a critical node
    
               dist_from_last_crit_node=0;
               last_visited_crit_node=critical_nodelist(index);
               endnode=last_visited_crit_node;

               vertex1_list(end+1)=startnode;
               vertex2_list(end+1)=endnode;
                
               startnode=endnode; % reset the starting node
    
           elseif mod(dist_from_last_crit_node+1,resolution)==0 % node at sampling resolution
               
               dist_from_last_crit_node=dist_from_last_crit_node+1;
               endnode=disc(I);
               
               vertex1_list(end+1)=startnode;
               vertex2_list(end+1)=endnode;

               startnode=endnode; % reset the starting node
    
           else % unimportant node, just skip over
               
               dist_from_last_crit_node=dist_from_last_crit_node+1;
               
           end
           
       end
end

% we will create a vertex list with consecutive nodes

vertex1=[vertex1_list vertex2_list]; % undirected graph

[vertex1_uniq,ii1,jj1]=unique(vertex1); % get the unique nodes and correspondnig numbering in jj1

shiftamount=length(vertex1_list);

jj2=circshift(jj1,[0 shiftamount]); % create the ending vertex list

edge_weights=ones(1,length(vertex1));

if isempty(vertex1_list), % the medial graph had only one node
    
    uMedialGraph=oldMedialGraph; % oldMedialGraph=[0] in this case 
    uNodePositions.x=oldNodePositions.x;
    uNodePositions.y=oldNodePositions.y;
    
else % otherwise treat as a general case
    
    uMedialGraph=sparse(jj1,jj2,edge_weights);
    uNodePositions.x=oldNodePositions.x(vertex1_uniq);
    uNodePositions.y=oldNodePositions.y(vertex1_uniq);
    
end

