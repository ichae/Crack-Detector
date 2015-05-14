
function [neuronPrimitiveGraph,connCompNodes]=make_kNN_graph(thresholded_image,k,resolution)

% neuronPrimitiveGraph=struct('adjmat',[],'clusterMat',[],'leaf_connectivity',[]);
neuronPrimitiveGraph=struct('adjmat',[],'leaf_connectivity',[]);
dirtSize=10;

SetCells = getCellParts(thresholded_image,dirtSize);

n = size(SetCells,1); % the number of connected components
gaussKer=fspecial('gaussian',[7 7]);
% resolution=5;
topLeftOffset = struct('x',[],'y',[]);
connCompNodes = cell(1,n);

for I=1:n, % Create a property list for the n conn comps 
    
    topLeftOffset.x=SetCells{I,2}(1);
    topLeftOffset.y=SetCells{I,2}(2);
    connCompNodes{I}=giveMedialGraph(SetCells{I,1},...
                     topLeftOffset,gaussKer,resolution);
    
end

distance_matrix=zeros(n); % matrix giving distance between each conn comp
cluster_matrix=zeros(n);
connectivity_info=cell(n);

for I=1:n, % create the distance matrix based on
           % (1) conn comp proximity and (2) tangent alignment
    
    num_I_leafnodes=connCompNodes{I}.LeafData.no_of_leaves;
    nodeweight_I=connCompNodes{I}.NodeValue;
    
    for J=1:I-1,
    
        num_J_leafnodes=connCompNodes{J}.LeafData.no_of_leaves;
        nodeweight_J=connCompNodes{J}.NodeValue;
        
%         num_I_leafnodes
%         num_J_leafnodes
        
        dist_options=zeros(num_I_leafnodes,num_J_leafnodes);
        
        for I_leaf=1:num_I_leafnodes,
            
            for J_leaf=1:num_J_leafnodes,
                
                euclid_dist=sqrt((connCompNodes{J}.LeafData.LeafPosition.x(J_leaf)...
                            -connCompNodes{I}.LeafData.LeafPosition.x(I_leaf)).^2 ...
                            +(connCompNodes{J}.LeafData.LeafPosition.y(J_leaf)...
                            -connCompNodes{I}.LeafData.LeafPosition.y(I_leaf)).^2);
                
                if ~isempty(connCompNodes{J}.LeafData.LeafTangents.tx),
                    
                    angleJleaf=givePrincipalArgumentAngle(...
                        connCompNodes{J}.LeafData.LeafTangents.tx(J_leaf),...
                        connCompNodes{J}.LeafData.LeafTangents.ty(J_leaf));
                else
                   
                    angleJleaf=[];
                    
                end
                
                if ~isempty(connCompNodes{I}.LeafData.LeafTangents.tx),
                    
                    angleIleaf=givePrincipalArgumentAngle(...
                        connCompNodes{I}.LeafData.LeafTangents.tx(I_leaf),...
                        connCompNodes{I}.LeafData.LeafTangents.ty(I_leaf));
                         
                else
                    
                    angleIleaf=[];
                    
                end
                    
                if ~isempty(angleIleaf) && ~isempty(angleJleaf),
                    
                    % acute angle between the tangents
                    angle=mod(abs(angleIleaf-angleJleaf),2*pi); %SB
%                       angle =  pi- abs(angleIleaf-angleJleaf); % Suvadip
                
                else
                    
                    angle=0; % bad match in angle if one medial graph is only a point
                    
                end
                
                % distance metric that is larger for small angles (min at
                % angle=pi), smaller for smaller euclidean distances,
                % and smaller for connecting to highly probable nodes
                
%                 dist_options(I_leaf,J_leaf)=euclid_dist/(10*angle*nodeweight_I*nodeweight_J+0.01);  
                
%                 dist_options(I_leaf,J_leaf)=(euclid_dist).^2 + (pi-angle);

%                 dist_options(I_leaf,J_leaf)=(euclid_dist).^2 + 180*(1-(angle/pi));% Saurav
%                     dist_options(I_leaf,J_leaf)=abs(pi-angle) ; % Suvadip
                  dist_options(I_leaf,J_leaf)=(abs(pi-angle))*(euclid_dist).^5;% Suvadip2
            end
            
        end
        
%         size(dist_options)
%         dist_options
        
        [choice_I,choice_J]=find(dist_options==min(dist_options(:)));
        choice_I=choice_I(1);
        choice_J=choice_J(1);
        
        distance_matrix(I,J)=dist_options(choice_I,choice_J); % closest distance possible
        connectivity_info{I,J}=[choice_I choice_J];
        distance_matrix(J,I)=distance_matrix(I,J); % symmetric
        connectivity_info{J,I}=[choice_J choice_I]; % symmetric
        
%         cluster_matrix(I,J)=dist_options(choice_I,choice_J)...
%             /(connCompNodes{I}.NodeValue+connCompNodes{J}.NodeValue);
%         cluster_matrix(J,I)=cluster_matrix(I,J); % symmetric
        
    end
    
end
    
if k>n, k=n;end % check for valid k


for I=1:n,

    nodelist=distance_matrix(I,:); % consider the distances from node I
    [sorted_list,sort_index]=sort(nodelist); % sort the distances
%     sort_index_k=sort_index(1:k); % take the nearest k nodes
    
    if k<n,
        
        sort_index_from_k=sort_index(k+1:n); % the farther neighbors
        
        for J=1:length(sort_index_from_k),
             
            distance_matrix(I,sort_index_from_k(J))=0; % disconnect the edges farther than kNN
%             cluster_matrix(I,sort_index_from_k(J))=0;
            
        end
        
    end
    
end

neuronPrimitiveGraph.adjmat=sparse(distance_matrix);
% neuronPrimitiveGraph.clusterMat=sparse(cluster_matrix);
neuronPrimitiveGraph.leaf_connectivity=connectivity_info;

    
