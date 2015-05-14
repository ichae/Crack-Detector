
function plot3DGraph(adjMat,NodePositions,bounds)

figure;axis([0 bounds.x 0 bounds.y 0 bounds.z]); hold on;
no_of_nodes=size(adjMat,1);

for I=1:no_of_nodes,
    
    for J=1:I-1,
        
        if adjMat(I,J)>0, % edge exists between I and J
            
            % plot the edges of the MST between the connected component
            % leaves
            
            plot3([NodePositions.x(I),NodePositions.x(J)],...
                  [NodePositions.y(I),NodePositions.y(J)],...
                  [NodePositions.z(I),NodePositions.z(J)],...
                  'm', 'Linewidth',4);drawnow;
            
        end
        
    end
    
end