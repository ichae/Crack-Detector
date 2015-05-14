
function [connCompNode]=giveMedialGraph(connCompImage,topLeftOffset,gaussKer,resolution)

[r,c]=size(connCompImage);

expand=20;

% expand the component in both direction to get rid of edge effects

connCompImage1=zeros(r+2*expand,c+2*expand);
connCompImage1(expand+1:expand+r,expand+1:expand+c)=connCompImage;

% diffuse the conn comp to smooth out jagged edges

connCompBlurred = imfilter(connCompImage1,gaussKer,'replicate','same');

connCompNew=zeros(size(connCompBlurred));

connCompNew(find(connCompBlurred>=0.2))=1; % a smoother conn comp

MG_image= bwmorph(connCompNew,'skel',Inf); % raw medial axis image

% shrink back to original size

% MG_image=MG_image1(expand+1:expand+r,expand+1:expand+c);

% imshow(MG_image,[]); hold on; pause;


%% %%%%%%%%%% Create an adjacency matrix from medial axis %%%%%%%%%%

r=r+2*expand;
c=c+2*expand;

clean1=zeros(r,c);
I=find(MG_image==1);
sz=length(I); % sz = number of nodes
AdjMat=diag(ones(1,sz));
count=0;

for row=1:r,  % number the nodes of the graph in this loop
    for col=1:c,
        if MG_image(row,col)==1,
            count=count+1;
            clean1(row,col)=count;
        end
    end
end

PosMat=zeros(sz,2);

for row=2:r-1,
    for col=2:c-1,
 
        if MG_image(row,col)==1,
            
            PosMat(clean1(row,col),1)=col-expand;
            PosMat(clean1(row,col),2)=row-expand;
      
            if MG_image(row-1,col-1)==1,
                AdjMat(clean1(row,col),clean1(row-1,col-1))=1.4142; % sqrt(2)
                AdjMat(clean1(row-1,col-1),clean1(row,col))=1.4142;
            end
            if MG_image(row-1,col)==1,
                AdjMat(clean1(row,col),clean1(row-1,col))=1;
                AdjMat(clean1(row-1,col),clean1(row,col))=1;
            end
            if MG_image(row,col-1)==1,
                AdjMat(clean1(row,col),clean1(row,col-1))=1;
                AdjMat(clean1(row,col-1),clean1(row,col))=1;
            end
            if MG_image(row+1,col+1)==1,
                AdjMat(clean1(row,col),clean1(row+1,col+1))=1.4142;
                AdjMat(clean1(row+1,col+1),clean1(row,col))=1.4142;
            end
            if MG_image(row+1,col)==1,
                AdjMat(clean1(row,col),clean1(row+1,col))=1;
                AdjMat(clean1(row+1,col),clean1(row,col))=1;
            end
            if MG_image(row,col+1)==1,
                AdjMat(clean1(row,col),clean1(row,col+1))=1;
                AdjMat(clean1(row,col+1),clean1(row,col))=1;
            end
            if MG_image(row-1,col+1)==1,
                AdjMat(clean1(row,col),clean1(row-1,col+1))=1.4142;
                AdjMat(clean1(row-1,col+1),clean1(row,col))=1.4142;
            end
            if MG_image(row+1,col-1)==1,
                AdjMat(clean1(row,col),clean1(row+1,col-1))=1.4142;
                AdjMat(clean1(row+1,col-1),clean1(row,col))=1.4142;
            end
        end
    end
end

%% %%%%%%%%%%%%%%%%%% critical node list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the number of critical nodes (leaf or bifurcation)

[G1,g_pred]=graphminspantree(sparse(AdjMat)); % To make sure there are no cycles

% G1 is assymetric and directed, make it symmetric

G=sparse(triu(full(G1))+triu(full(G1))'+tril(full(G1))+tril(full(G1))');

nz_index=find(full(G)~=0);
AdjMat1=zeros(sz);
AdjMat1(nz_index)=1; % just the connection information
node_degree=sum(AdjMat1,1);
critical_nodelist=find(node_degree~=2);

%% %%%%% create subsampled graph %%%%%%%

S=critical_nodelist(1); % start node for DFS search

[disc, pred, closed]= graphtraverse(G, S, 'Directed', 'false');

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
 
connCompNode=struct('MedialGraph',[],...
                    'NodePositions',[],...
                    'LeafData',[],...
                    'NodeValue',[]);
                
LeafData=struct('no_of_leaves',[],...
                    'LeafNodes',[],...
                    'LeafPosition',[],...
                    'LeafTangents',[]);

% we will create a vertex list with consecutive nodes

vertex1=[vertex1_list vertex2_list]; % undirected graph

[vertex1_uniq,ii1,jj1]=unique(vertex1); % get the unique nodes and correspondnig numbering in jj1

shiftamount=length(vertex1_list);

jj2=circshift(jj1,[0 shiftamount]); % create the ending vertex list

edge_weights=ones(1,length(vertex1));

% size(jj1)
% size(jj2)
% size(edge_weights)

if isempty(vertex1_list), % the medial graph had only one node
    
    connCompNode.MedialGraph=G; % G=[0] in this case 
%     connCompNode.NodeValue=0;
    connCompNode.NodePositions.x=PosMat(1,1)+topLeftOffset.x;
    connCompNode.NodePositions.y=PosMat(1,2)+topLeftOffset.y;
    
else % otherwise treat as a general case
    
    connCompNode.MedialGraph=sparse(jj1,jj2,edge_weights);
%     connCompNode.NodeValue=sum(edge_weights)/2; % higher for more elongated conn comps
    
    connCompNode.NodePositions=struct('x',[],'y',[]);
    
    connCompNode.NodePositions.x=PosMat(vertex1_uniq,1)+topLeftOffset.x;
    connCompNode.NodePositions.y=PosMat(vertex1_uniq,2)+topLeftOffset.y;

end

% G
% vertex1
% vertex1_uniq
% connCompNode.NodePositions.x
% r
% c

% Merge and smooth the nodes which are nearer than resolution

[connCompNode.MedialGraph,connCompNode.NodePositions,indexMap]=...
        mergeAndSmoothNodes(connCompNode.MedialGraph,...
        connCompNode.NodePositions,resolution);
    

% clear AdjMat G PosMat;

dummyMG=full(connCompNode.MedialGraph);

new_nodedegree=sum(dummyMG,1);
% num_of_branches=sum(new_nodedegree(new_nodedegree~=2))/2;
connCompArea=sum(connCompImage(:));
connCompLength=length(connCompNode.NodePositions.x)-1;
% connCompNode.NodeValue=(connCompLength/(num_of_branches+0.01))...
%                          + sqrt(connCompArea);

connCompNode.NodeValue=connCompLength + sqrt(connCompArea);

leaf_nodelist=find(new_nodedegree==1); % position of leaves in vertex list

LeafData.no_of_leaves=length(leaf_nodelist);
LeafData.LeafNodes=leaf_nodelist;

LeafData.LeafPosition=struct('x',[],'y',[]);
LeafData.LeafTangents=struct('tx',[],'ty',[]);

LeafData.LeafPosition.x=connCompNode.NodePositions.x(leaf_nodelist);
LeafData.LeafPosition.y=connCompNode.NodePositions.y(leaf_nodelist);

for count=1:LeafData.no_of_leaves,
    
    connectivity=dummyMG(LeafData.LeafNodes(count),:);
    prev_node_index=find(connectivity==1);
    prev_node=prev_node_index(1); % redundant, since there is only one
                                  % connecting node to a leaf
    prev_node_x=connCompNode.NodePositions.x(prev_node);
    prev_node_y=connCompNode.NodePositions.y(prev_node);
    
    LeafData.LeafTangents.tx(end+1)=...
               LeafData.LeafPosition.x(count)-prev_node_x;
           
    LeafData.LeafTangents.ty(end+1)=...
               LeafData.LeafPosition.y(count)-prev_node_y;       
           
end

% set leaf node to graph node if graph consists of 1 node only

if isempty(leaf_nodelist),
    
    LeafData.no_of_leaves=1;
    LeafData.LeafNodes=1;
    LeafData.LeafPosition.x=connCompNode.NodePositions.x(1);
    LeafData.LeafPosition.y=connCompNode.NodePositions.y(1);
    
end

connCompNode.LeafData=LeafData;

% scatter(connCompNode.NodePositions.x-topLeftOffset.x,connCompNode.NodePositions.y-topLeftOffset.y,'r');


