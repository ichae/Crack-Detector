
function toSWC(Img,nGraph,nodePos,rootNode,Rmin,Rmax,flip)




nNodes = size(nGraph,1);

% make sure the nodes are in dfs order in the swc file
% otherwise there will be a parsing error
[row col depth] = size(Img);

[d,p,c]=graphtraverse(nGraph,rootNode,'Directed', 'false');
G=full(nGraph);
newG=zeros(nNodes);

m = numel(d);
% fprintf('\n number of nodes = %d',length(nodePos.x));
% fprintf('\n Size that is shown = %d',m);
%m = nNodes; % changed by suvadip

for I=1:m,
    for J=1:m,
        %         disp([d(I) d(J)])
        if G(d(I),d(J))==1,
            
            newG(I,J)=1;
            
        end
        
    end
end
nGraph=sparse(newG);
nodePos.x=nodePos.x(d);
nodePos.y=nodePos.y(d);
nodePos.z=nodePos.z(d);


testImg = zeros(size(Img));
% testImg(round(nodePos.x),round(nodePos.y),round(nodePos.z)) = 1;


% nodePos.x=nodePos.x(d);
% nodePos.y=nodePos.y(d);
% nodePos.z=nodePos.z(d);
%

[fname,pname]=uiputfile('*.swc','Save SWC File as');
p_write = strcat(pname,fname);
fid=fopen(p_write,'wt');

header='##n,type,x,y,z,radius,parent';
fprintf(fid,'%s\n\n',header);

if flip
    fprintf(fid,'1 3 %6.3f %6.3f %6.3f 2 -1\n',...
        (nodePos.x(1)),(mod(-nodePos.y(1),row+1)),(nodePos.z(1)));
else
    fprintf(fid,'1 3 %6.3f %6.3f %6.3f 2 -1\n',...
        (nodePos.x(1)),nodePos.y(1),(nodePos.z(1)));
end
% value=zeros(1,7);

[disc, pred, closed]= graphtraverse(nGraph, 1, 'Directed', 'false');

for I=2:m,
    %    [ nodePos.x(I)    nodePos.y(I)     nodePos.z(I)]
    
    parent=pred(I);
    value(I,1)=I;
    value(I,2)=3; % default
    value(I,3)=nodePos.x(I);
    if flip
        value(I,4)=mod(-round(nodePos.y(I)),row+1);
    else
        value(I,4)=nodePos.y(I);
    end
    value(I,5)=nodePos.z(I);
    %     value(6)=2; % default
    value(I,6) = findRadius(Img,round(nodePos.y(I)),round(nodePos.x(I)),round(nodePos.z(I)),Rmin,Rmax);
    % % need to implement this function
    value(I,7)=parent;
    %     newline=num2str(value);
    %     fprintf(fid,'%s\n',newline);
    
    %     testImg(round(value(4)),round(value(3)),round(value(5))) = 1;
end
% write3Dstack(testImg,1,1,'check');

%% Do a Gaussian smoothing of the radii
Radii = value(:,6);
smoothRad = imfilter(Radii,1/7*[1 1 1 1 1 1 1],'symmetric','same');
value(:,6) = smoothRad;
%%
for I=2:m,
    newline=num2str(value(I,:));
    fprintf(fid,'%s\n',newline);
end


fclose(fid);
end


% function yval = flipY(y,I)
%
%  [row col depth] = size(I);
%
%
% end




