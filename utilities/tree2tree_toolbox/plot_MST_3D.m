
function plot_MST_3D(MST,connCompNodes,connectivity,bounds)

% figure;imshow(image,[],'init', 'fit');hold on;
figure;axis([0 bounds.x 0 bounds.y 0 bounds.z]); hold on;
no_of_nodes=length(connCompNodes);

%first plot the connected components

for I=1:no_of_nodes,
    
    mstf=full(connCompNodes{I}.MedialGraph);
    n=size(mstf,1);
    for J=1:n,
        
        for K=1:J-1,
            
            if mstf(J,K)>0,
                
                plot3([connCompNodes{I}.NodePositions.x(J) connCompNodes{I}.NodePositions.x(K)],...
                     [connCompNodes{I}.NodePositions.y(J) connCompNodes{I}.NodePositions.y(K)],...
                     [connCompNodes{I}.NodePositions.z(J) connCompNodes{I}.NodePositions.z(K)],...
                     'b','Linewidth',3);drawnow;
                 
            end
            
        end
        
    end
   
    scatter3(connCompNodes{I}.NodePositions.x,connCompNodes{I}.NodePositions.y,...
            connCompNodes{I}.NodePositions.z,5,'r','filled');
    drawnow;
    
%     gplot(connCompNodes{I}.MedialGraph,...
%           [connCompNodes{I}.NodePositions.x' connCompNodes{I}.NodePositions.y']);drawnow;
      
end

MST_full=full(MST);

for I=1:no_of_nodes,
    
    for J=1:no_of_nodes,
        
        if MST_full(I,J)>0, % edge exists between I and J
            
            leaf_I=connectivity{I,J}(1); % leaf numbers for the connection
            leaf_J=connectivity{I,J}(2);
            
            leaf_I_x=connCompNodes{I}.LeafData.LeafPosition.x(leaf_I);
            leaf_I_y=connCompNodes{I}.LeafData.LeafPosition.y(leaf_I);
            leaf_I_z=connCompNodes{I}.LeafData.LeafPosition.z(leaf_I);
            leaf_J_x=connCompNodes{J}.LeafData.LeafPosition.x(leaf_J);
            leaf_J_y=connCompNodes{J}.LeafData.LeafPosition.y(leaf_J);
            leaf_J_z=connCompNodes{J}.LeafData.LeafPosition.z(leaf_J);
            
            % plot the edges of the MST between the connected component
            % leaves
            
            plot3([leaf_I_x,leaf_J_x],[leaf_I_y,leaf_J_y],...
                 [leaf_I_z,leaf_J_z],'g', 'Linewidth',3);drawnow;
            
        end
        
    end
    
end

hold off;