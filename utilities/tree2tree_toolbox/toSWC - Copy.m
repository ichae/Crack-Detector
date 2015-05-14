
function toSWC(nGraph,nodePos,rootNode,row,col,depth)

nNodes = size(nGraph,1);

% make sure the nodes are in dfs order in the swc file
% otherwise there will be a parsing error

[d,p,c]=graphtraverse(nGraph,rootNode,'Directed', 'false');
G=full(nGraph);
newG=zeros(nNodes);

for I=1:nNodes,
    for J=1:nNodes,
        
        if G(d(I),d(J))==1,
            
            newG(I,J)=1;
            
        end
        
    end
end
nGraph=sparse(newG);    
nodePos.x=nodePos.x(d);
nodePos.y=nodePos.y(d);
nodePos.z=nodePos.z(d);

% nodePos.x=nodePos.x(d);
% nodePos.y=nodePos.y(d);
% nodePos.z=nodePos.z(d);
% 

[fname,pname]=uiputfile('*.swc','Save SWC File as');
p_write = strcat(pname,fname);
fid=fopen(p_write,'wt');

header='##n,type,x,y,z,radius,parent';
fprintf(fid,'%s\n\n',header);

fprintf(fid,'1 3 %6.3f %6.3f %6.3f 5 -1\n',...
        (nodePos.x(1)),(nodePos.y(1)),(nodePos.z(1)));

value=zeros(1,7);

[disc, pred, closed]= graphtraverse(nGraph, 1, 'Directed', 'false');

for I=2:nNodes,
    
    parent=pred(I);
    value(1)=I;
    value(2)=3; % default
    value(3)=nodePos.x(I);
    value(4)=nodePos.y(I);
    value(5)=nodePos.z(I);
    value(6)=2; % default
    value(7)=parent;
    newline=num2str(value);
    fprintf(fid,'%s\n',newline);

end

fclose(fid);


    


