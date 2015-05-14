
function connCompNode=give3DMedialGraph(connComp,resolution)

[r,c,s]=size(connComp.binImg);
offset=connComp.Offset;

MG_image = find3DSkeleton(connComp.binImg); % raw medial axis image

%%%%%%%%%%%% Create an adjacency matrix from medial axis %%%%%%%%%%

labelImg=zeros(r,c,s);
I=find(MG_image==1);

sz=length(I); % sz = number of nodes
% AdjMat=diag(ones(1,sz));
AdjMat=zeros(sz);


% labelImg(I)=I; % give unique node number to each white node
% size(labelImg)

count=0;
for y=1:r,  % number the nodes of the graph in this loop
    for x=1:c,
        for z=1:s,
            if MG_image(y,x,z)==1,
                count=count+1;
                labelImg(y,x,z)=count;  % ordering the bright pixels
            end
        end
    end
end

Pos.x=zeros(1,sz);
Pos.y=zeros(1,sz);
Pos.z=zeros(1,sz);

for y=2:r-1,
    for x=2:c-1,
        for z=2:s-1,
            
            if MG_image(y,x,z)==1,
                
                % This is the current position
                Pos.x(labelImg(y,x,z))=x;
                Pos.y(labelImg(y,x,z))=y;
                Pos.z(labelImg(y,x,z))=z;
                
                % find the children of the node 'pos'
                if MG_image(y,x,z+1)==1, % U
                    AdjMat(labelImg(y,x,z),labelImg(y,x,z+1))=1; 
                    AdjMat(labelImg(y,x,z+1),labelImg(y,x,z))=1;
                end
                if MG_image(y,x,z-1)==1, % D
                    AdjMat(labelImg(y,x,z),labelImg(y,x,z-1))=1; 
                    AdjMat(labelImg(y,x,z-1),labelImg(y,x,z))=1;
                end
                if MG_image(y+1,x,z)==1, % N
                    AdjMat(labelImg(y,x,z),labelImg(y+1,x,z))=1; 
                    AdjMat(labelImg(y+1,x,z),labelImg(y,x,z))=1;
                end
                if MG_image(y-1,x,z)==1, % S
                    AdjMat(labelImg(y,x,z),labelImg(y-1,x,z))=1; 
                    AdjMat(labelImg(y-1,x,z),labelImg(y,x,z))=1;
                end             
                if MG_image(y,x+1,z)==1, % E
%                     labelImg(y,x,z)
%                     labelImg(y,x+1,z)
                    AdjMat(labelImg(y,x,z),labelImg(y,x+1,z))=1; 
                    AdjMat(labelImg(y,x+1,z),labelImg(y,x,z))=1;
                end                 
                if MG_image(y,x-1,z)==1, % W
                    AdjMat(labelImg(y,x,z),labelImg(y,x-1,z))=1; 
                    AdjMat(labelImg(y,x-1,z),labelImg(y,x,z))=1;
                end 
                if MG_image(y+1,x+1,z)==1, % NE
                    AdjMat(labelImg(y,x,z),labelImg(y+1,x+1,z))=1.4142; % sqrt(2) 
                    AdjMat(labelImg(y+1,x+1,z),labelImg(y,x,z))=1.4142; % sqrt(2)
                end                 
                if MG_image(y-1,x+1,z)==1, % SE
                    AdjMat(labelImg(y,x,z),labelImg(y-1,x+1,z))=1.4142; % sqrt(2) 
                    AdjMat(labelImg(y-1,x+1,z),labelImg(y,x,z))=1.4142; % sqrt(2)
                end                
                if MG_image(y+1,x-1,z)==1, % NW
                    AdjMat(labelImg(y,x,z),labelImg(y+1,x-1,z))=1.4142; % sqrt(2) 
                    AdjMat(labelImg(y+1,x-1,z),labelImg(y,x,z))=1.4142; % sqrt(2)
                end                 
                if MG_image(y-1,x-1,z)==1, % SW
                    AdjMat(labelImg(y,x,z),labelImg(y-1,x-1,z))=1.4142; % sqrt(2) 
                    AdjMat(labelImg(y-1,x-1,z),labelImg(y,x,z))=1.4142; % sqrt(2)
                end               
                if MG_image(y+1,x,z+1)==1, % UN
                    AdjMat(labelImg(y,x,z),labelImg(y+1,x,z+1))=1.4142; % sqrt(2) 
                    AdjMat(labelImg(y+1,x,z+1),labelImg(y,x,z))=1.4142; % sqrt(2)
                end
                if MG_image(y-1,x,z+1)==1, % US
                    AdjMat(labelImg(y,x,z),labelImg(y-1,x,z+1))=1.4142; % sqrt(2) 
                    AdjMat(labelImg(y-1,x,z+1),labelImg(y,x,z))=1.4142; % sqrt(2)
                end              
                if MG_image(y+1,x,z-1)==1, % DN
                    AdjMat(labelImg(y,x,z),labelImg(y+1,x,z-1))=1.4142; % sqrt(2) 
                    AdjMat(labelImg(y+1,x,z-1),labelImg(y,x,z))=1.4142; % sqrt(2)
                end               
                if MG_image(y-1,x,z-1)==1, % DS
                    AdjMat(labelImg(y,x,z),labelImg(y-1,x,z-1))=1.4142; % sqrt(2) 
                    AdjMat(labelImg(y-1,x,z-1),labelImg(y,x,z))=1.4142; % sqrt(2)
                end                
                if MG_image(y,x+1,z+1)==1, % UE
                    AdjMat(labelImg(y,x,z),labelImg(y,x+1,z+1))=1.4142; % sqrt(2) 
                    AdjMat(labelImg(y,x+1,z+1),labelImg(y,x,z))=1.4142; % sqrt(2)
                end
                if MG_image(y,x+1,z-1)==1, % DE
                    AdjMat(labelImg(y,x,z),labelImg(y,x+1,z-1))=1.4142; % sqrt(2) 
                    AdjMat(labelImg(y,x+1,z-1),labelImg(y,x,z))=1.4142; % sqrt(2)
                end               
                if MG_image(y,x-1,z+1)==1, % UW
                    AdjMat(labelImg(y,x,z),labelImg(y,x-1,z+1))=1.4142; % sqrt(2) 
                    AdjMat(labelImg(y,x-1,z+1),labelImg(y,x,z))=1.4142; % sqrt(2)
                end               
                if MG_image(y,x-1,z-1)==1, % DW
                    AdjMat(labelImg(y,x,z),labelImg(y,x-1,z-1))=1.4142; % sqrt(2) 
                    AdjMat(labelImg(y,x-1,z-1),labelImg(y,x,z))=1.4142; % sqrt(2)
                end                  
                if MG_image(y+1,x+1,z+1)==1, % UNE
                    AdjMat(labelImg(y,x,z),labelImg(y+1,x+1,z+1))=1.7321; % sqrt(3) 
                    AdjMat(labelImg(y+1,x+1,z+1),labelImg(y,x,z))=1.7321; % sqrt(3)
                end               
                if MG_image(y+1,x+1,z-1)==1, % DNE
                    AdjMat(labelImg(y,x,z),labelImg(y+1,x+1,z-1))=1.7321; % sqrt(3) 
                    AdjMat(labelImg(y+1,x+1,z-1),labelImg(y,x,z))=1.7321; % sqrt(3)
                end                 
                if MG_image(y+1,x-1,z+1)==1, % UNW
                    AdjMat(labelImg(y,x,z),labelImg(y+1,x-1,z+1))=1.7321; % sqrt(3) 
                    AdjMat(labelImg(y+1,x-1,z+1),labelImg(y,x,z))=1.7321; % sqrt(3)
                end 
                if MG_image(y+1,x-1,z-1)==1, % DNW
                    AdjMat(labelImg(y,x,z),labelImg(y+1,x-1,z-1))=1.7321; % sqrt(3) 
                    AdjMat(labelImg(y+1,x-1,z-1),labelImg(y,x,z))=1.7321; % sqrt(3)
                end                 
                if MG_image(y-1,x-1,z+1)==1, % USW
                    AdjMat(labelImg(y,x,z),labelImg(y-1,x-1,z+1))=1.7321; % sqrt(3) 
                    AdjMat(labelImg(y-1,x-1,z+1),labelImg(y,x,z))=1.7321; % sqrt(3)
                end
                if MG_image(y-1,x+1,z+1)==1, % USE
                    AdjMat(labelImg(y,x,z),labelImg(y-1,x+1,z+1))=1.7321; % sqrt(3)
                    AdjMat(labelImg(y-1,x+1,z+1),labelImg(y,x,z))=1.7321; % sqrt(3)
                end
                if MG_image(y-1,x-1,z-1)==1, % DSW
                    AdjMat(labelImg(y,x,z),labelImg(y-1,x-1,z-1))=1.7321; % sqrt(3)
                    AdjMat(labelImg(y-1,x-1,z-1),labelImg(y,x,z))=1.7321; % sqrt(3)
                end
                if MG_image(y-1,x+1,z-1)==1, % DSE
                    AdjMat(labelImg(y,x,z),labelImg(y-1,x+1,z-1))=1.7321; % sqrt(3)
                    AdjMat(labelImg(y-1,x+1,z-1),labelImg(y,x,z))=1.7321; % sqrt(3)
                end
                 
            end
        end
    end
end

% spy(sparse(AdjMat))
% sparse(AdjMat)
% pause
% sum(sum(AdjMat-AdjMat'))
% [L,Num]=bwlabeln(MG_image);
% Num
% [how_many_conn_comps,C] = graphconncomp(sparse(AdjMat),'DIRECTED',false);
% DIST = GRAPHALLSHORTESTPATHS(sparse(AdjMat));
% DIST(1,:)
% spy(sparse(DIST));
% sum(sum(isinf(DIST)))
% how_many_conn_comps
% pause

%% %%%%%%%%%%%%%%%%%% critical node list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the number of critical nodes (leaf or bifurcation)
% the above algo produces the graph from the points
[G1,g_pred]=graphminspantree(sparse(AdjMat)); % To make sure there are no cycles

% G1 is assymetric and directed, make it symmetric

G=sparse(triu(full(G1))+triu(full(G1))'+tril(full(G1))+tril(full(G1))');

% nz_index=find(G~=0);
AdjMat1=zeros(sz);
AdjMat1(G~=0)=1; % just the connection information
node_degree=sum(AdjMat1,1);
critical_nodelist=find(node_degree~=2);


%% %%%%% create subsampled graph %%%%%%%

S=critical_nodelist(1); % start node for DFS search

% S
% sparse(AdjMat1)
% [how_many_conn_comps,C] = graphconncomp(sparse(AdjMat1));
% how_many_conn_comps
% pause

[disc, pred, closed]= graphtraverse(sparse(AdjMat1), S, 'Directed', 'false');

% disc
% pause

startnode=S;
last_visited_crit_node=S;
endnode=0;
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

%%%%%%%%%%%%%%% sub sampled graph created %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 %% %%% data structure for output data %%%%%%
 
connCompNode = struct('MedialGraph',[],...
                    'NodePositions',[],...
                    'LeafData',[],...
                    'NodeValue',[]);
                
LeafData = struct('no_of_leaves',[],...
                    'LeafNodes',[],...
                    'LeafPosition',[],...
                    'LeafTangents',[]);

% we will create a vertex list with consecutive nodes

vertex1 = [vertex1_list vertex2_list]; % undirected graph

[vertex1_uniq,ii1,jj1]=unique(vertex1); % get the unique nodes and correspondnig numbering in jj1

shiftamount=length(vertex1_list);

jj2=circshift(jj1,[0 shiftamount]); % create the ending vertex list

edge_weights=ones(1,length(vertex1));

% size(jj1)
% size(jj2)
% size(edge_weights)
% vertex1_list
% pause

if isempty(vertex1_list), % the medial graph had only one node
    
    connCompNode.MedialGraph=G; % G=[0] in this case 
%     connCompNode.NodeValue=0;
    connCompNode.NodePositions.x=Pos.x+offset(1);
    connCompNode.NodePositions.y=Pos.y+offset(2);
    connCompNode.NodePositions.z=Pos.z+offset(3);
    
    
else % otherwise treat as a general case
    
    connCompNode.MedialGraph=sparse(jj1,jj2,edge_weights);
%     connCompNode.NodeValue=sum(edge_weights)/2; % higher for more elongated conn comps
%     jj1
%     jj2
%     pause
    connCompNode.NodePositions=struct('x',[],'y',[],'z',[]);
    connCompNode.NodePositions.x=Pos.x(vertex1_uniq)+offset(1);
    connCompNode.NodePositions.y=Pos.y(vertex1_uniq)+offset(2);
    connCompNode.NodePositions.z=Pos.z(vertex1_uniq)+offset(3);
 
end

% Merge and smooth the nodes which are nearer than resolution

[connCompNode.MedialGraph,connCompNode.NodePositions,indexMap]=...
        mergeAndSmoothNodes3D(connCompNode.MedialGraph,...
        connCompNode.NodePositions,resolution);
    

% clear AdjMat G PosMat;

dummyMG=full(connCompNode.MedialGraph);

new_nodedegree=sum(dummyMG,1);
% num_of_branches=sum(new_nodedegree(new_nodedegree~=2))/2;
connCompArea=sum(connComp.binImg(:));
connCompLength=length(connCompNode.NodePositions.x)-1;
% connCompNode.NodeValue=(connCompLength/(num_of_branches+0.01))...
%                          + sqrt(connCompArea);

% connCompNode.NodeValue=connCompLength + sqrt(connCompArea);
connCompNode.NodeValue=connCompLength + ((connCompArea)^(1/3));

leaf_nodelist=find(new_nodedegree==1); % position of leaves in vertex list

LeafData.no_of_leaves=length(leaf_nodelist);
LeafData.LeafNodes=leaf_nodelist;

LeafData.LeafPosition=struct('x',[],'y',[],'z',[]);
LeafData.LeafTangents=struct('tx',[],'ty',[],'tz',[]);

LeafData.LeafPosition.x=connCompNode.NodePositions.x(leaf_nodelist);
LeafData.LeafPosition.y=connCompNode.NodePositions.y(leaf_nodelist);
LeafData.LeafPosition.z=connCompNode.NodePositions.z(leaf_nodelist);

for count=1:LeafData.no_of_leaves,
    
    connectivity=dummyMG(LeafData.LeafNodes(count),:);
    prev_node_index=find(connectivity==1);
    prev_node=prev_node_index(1); % redundant, since there is only one
                                  % connecting node to a leaf
    prev_node_x=connCompNode.NodePositions.x(prev_node);
    prev_node_y=connCompNode.NodePositions.y(prev_node);
    prev_node_z=connCompNode.NodePositions.z(prev_node);
    
    LeafData.LeafTangents.tx(end+1)=...
               LeafData.LeafPosition.x(count)-prev_node_x;
           
    LeafData.LeafTangents.ty(end+1)=...
               LeafData.LeafPosition.y(count)-prev_node_y;  
           
    LeafData.LeafTangents.tz(end+1)=...
               LeafData.LeafPosition.z(count)-prev_node_z;
           
end

% set leaf node to graph node if graph consists of 1 node only

if isempty(leaf_nodelist),
    
    LeafData.no_of_leaves=1;
    LeafData.LeafNodes=1;
    LeafData.LeafPosition.x=connCompNode.NodePositions.x(1);
    LeafData.LeafPosition.y=connCompNode.NodePositions.y(1);
    LeafData.LeafPosition.z=connCompNode.NodePositions.z(1);
    
end

connCompNode.LeafData=LeafData;
