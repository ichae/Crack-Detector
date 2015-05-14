
function [neuronPrimitiveGraph,connCompNodes] = make_3D_kNN_graph(binImg3D,k,resolution)

neuronPrimitiveGraph = struct('adjmat',[],'leaf_connectivity',[]);
connComps = get3DConnComps(binImg3D);

n=size(connComps,2); % the number of connected components
fprintf('\n Number of conncted components = %d',n);
connCompNodes=cell(1,n);

for I=1:n, % Create a property list for the n conn comps 
    
    connCompNodes{I} = give3DMedialGraph(connComps{I},resolution);                
    
end

if n > 1 % more than one CC.. changed by SUVADIP
    
distance_matrix=zeros(n); % matrix giving distance between each conn comp
connectivity_info=cell(n);

for I=1:n, % create the distance matrix based on
           % (1) conn comp proximity (2) tangent twist and 
           % (3) tangent bend
    
    num_I_leafnodes=connCompNodes{I}.LeafData.no_of_leaves;
    
    for J=1:I-1,
    
        num_J_leafnodes=connCompNodes{J}.LeafData.no_of_leaves;
        
%         num_I_leafnodes
%         num_J_leafnodes
        
        dist_options=zeros(num_I_leafnodes,num_J_leafnodes);
        
        for I_leaf=1:num_I_leafnodes,
            
            for J_leaf=1:num_J_leafnodes,
                
                euclid_dist_sq=(connCompNodes{J}.LeafData.LeafPosition.x(J_leaf)...
                            -connCompNodes{I}.LeafData.LeafPosition.x(I_leaf)).^2 ...
                            +(connCompNodes{J}.LeafData.LeafPosition.y(J_leaf)...
                            -connCompNodes{I}.LeafData.LeafPosition.y(I_leaf)).^2 ...
                           +(connCompNodes{J}.LeafData.LeafPosition.z(J_leaf)...
                            -connCompNodes{I}.LeafData.LeafPosition.z(I_leaf)).^2;
                        

                
                if ~isempty(connCompNodes{J}.LeafData.LeafTangents.tx) && ...
                   ~isempty(connCompNodes{I}.LeafData.LeafTangents.tx),
                   
               
                  tanJ_vector=[connCompNodes{J}.LeafData.LeafTangents.tx(J_leaf);
                      connCompNodes{J}.LeafData.LeafTangents.ty(J_leaf);
                      connCompNodes{J}.LeafData.LeafTangents.tz(J_leaf)];

               
                  tanI_vector=[connCompNodes{I}.LeafData.LeafTangents.tx(I_leaf);
                      connCompNodes{I}.LeafData.LeafTangents.ty(I_leaf);
                      connCompNodes{I}.LeafData.LeafTangents.tz(I_leaf)];

               
                  connect_vector=[(connCompNodes{J}.LeafData.LeafPosition.x(J_leaf)...
                      -connCompNodes{I}.LeafData.LeafPosition.x(I_leaf));...
                      (connCompNodes{J}.LeafData.LeafPosition.y(J_leaf)...
                      -connCompNodes{I}.LeafData.LeafPosition.y(I_leaf));...
                      (connCompNodes{J}.LeafData.LeafPosition.z(J_leaf)...
                      -connCompNodes{I}.LeafData.LeafPosition.z(I_leaf))];
                  connect_vector=connect_vector/norm(connect_vector);
                  
                  if isnan(connect_vector),
                      
                      error('connect_vector is nan');
                      
                  end
               
                  rotated_tanI_vector_axial=dot(tanI_vector,connect_vector)*connect_vector;
                  twist_vec_1=tanI_vector-rotated_tanI_vector_axial;
                  normal_vec=cross(tanJ_vector,connect_vector);
                  dummy=normal_vec;
                  normal_vec=normal_vec/(norm(normal_vec)+eps);
                  
                    if isnan(normal_vec),
                      
                        norm(dummy)
                        tanJ_vector
                        connect_vector
                      error('normal_vec is nan');
                      
                    end
                  
                  twist_vec_2=twist_vec_1-dot(normal_vec,twist_vec_1)*normal_vec;
                  twist_vec_2=(norm(twist_vec_1)/(norm(twist_vec_2)+eps))*twist_vec_2;
                  twist_amount=norm((twist_vec_2/(norm(twist_vec_2)+eps))-...
                      (twist_vec_1)/(norm(twist_vec_1)+eps));
                  
                  if isnan(twist_vec_2),
                      
                      error('twist_vec_2 is nan');
                      
                  end
                  
                  if isnan(twist_vec_1),
                      
                      error('twist_vec_1 is nan');
                      
                  end
                  
                  rotated_tanI_vector=rotated_tanI_vector_axial+twist_vec_2;
                  tanJ_vector_norm=tanJ_vector/norm(tanJ_vector);
                  
                  if isnan(tanJ_vector) ,
                      
                      error('tanJ_vector is nan');
                      
                  end
                  
                  rotated_tanI_vector_norm=rotated_tanI_vector/norm(rotated_tanI_vector);
                  
                  if isnan(rotated_tanI_vector),
                      
                      error('rotated_tanI_vector is nan');
                      
                  end
                  
                  bend_amount=norm(-tanJ_vector_norm-rotated_tanI_vector_norm); % first minus is necessary !!
                  
                  if isnan(twist_amount) ,
                      
                      error('twist_amount is nan');
                      
                  end
                  if isnan(bend_amount) ,
                      
%                       error('bend_amount is nan');
                        bend_amount = 0;    % changed by suvadip
                      
                  end
     
                  angle=twist_amount+bend_amount;
                  
                else
                    
                  angle = 8; % 2^2+2^2, bad match
                    
                end
                
                
                % distance metric that is larger for smaller for small bends and twists,
                % smaller for smaller euclidean distances,

                dist_options(I_leaf,J_leaf)= euclid_dist_sq + angle.^2;
                
                if isnan(euclid_dist_sq),
                    error('dist is nan');
                end
                if isnan(angle),
                    error('angle is nan');
                end
                
            end
            
        end
        
%         size(dist_options)
        dist_options
        
        [choice_I,choice_J]=find(dist_options==min(dist_options(:)));
        choice_I=choice_I(1);
        choice_J=choice_J(1);
        
        distance_matrix(I,J)=dist_options(choice_I,choice_J); % closest distance possible
        connectivity_info{I,J}=[choice_I choice_J];
        distance_matrix(J,I)=distance_matrix(I,J); % symmetric
        connectivity_info{J,I}=[choice_J choice_I]; % symmetric
       
    end
    
end
    
if k>n, k=n;end % check for valid k


for I=1:n,

    nodelist=distance_matrix(I,:); % consider the distances from node I
    [sorted_list,sort_index]=sort(nodelist); % sort the distances
    
    if k<n,
        
        sort_index_from_k=sort_index(k+1:n); % the farther neighbors
        
        for J=1:length(sort_index_from_k),
             
            distance_matrix(I,sort_index_from_k(J))=0; % disconnect the edges farther than kNN
%             cluster_matrix(I,sort_index_from_k(J))=0;
            
        end
        
    end
    
end

neuronPrimitiveGraph.adjmat=sparse(distance_matrix);
neuronPrimitiveGraph.leaf_connectivity=connectivity_info;

else
    
neuronPrimitiveGraph.adjmat=[];
neuronPrimitiveGraph.leaf_connectivity=[];

    
end
