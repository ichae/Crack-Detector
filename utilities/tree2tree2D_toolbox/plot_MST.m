
function plot_MST(MST,connCompNodes,connectivity,image)

figure;imshow(image,[],'init', 'fit');hold on;
no_of_nodes=length(connCompNodes);

%first plot the connected components

for I=1:no_of_nodes,
    
    mstf=full(connCompNodes{I}.MedialGraph);
    n=size(mstf,1);
    for J=1:n,
        
        for K=1:J-1,
            
            if mstf(J,K)>0,
                
                plot([connCompNodes{I}.NodePositions.x(J) connCompNodes{I}.NodePositions.x(K)],...
                     [connCompNodes{I}.NodePositions.y(J) connCompNodes{I}.NodePositions.y(K)],...
                     'b','Linewidth',3);drawnow;
                 
            end
            
        end
        
    end
   
%     scatter(connCompNodes{I}.NodePositions.x,connCompNodes{I}.NodePositions.y,5,'r','filled');
%     drawnow;
    
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
            leaf_J_x=connCompNodes{J}.LeafData.LeafPosition.x(leaf_J);
            leaf_J_y=connCompNodes{J}.LeafData.LeafPosition.y(leaf_J);
            
            % plot the edges of the MST between the connected component
            % leaves
            
            plot([leaf_I_x,leaf_J_x],[leaf_I_y,leaf_J_y],'g', 'Linewidth',3);drawnow;
            
        end
        
    end
    
end

hold off;