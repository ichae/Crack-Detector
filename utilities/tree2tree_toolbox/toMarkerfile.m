function [] = toMarkerfile(skel3D)
%TOMARKERFILE Create the markerfile for V3D
%   skel3D: The 3D skeleton of the binarized image. The points serve as the
%   markers


[row, col, depth] = size(skel3D);

numPts = length(find(skel3D == 1));

nodePos = getMarkers(skel3D,row,col,depth,numPts);

% Create the markerfile
[fname,pname]=uiputfile('*.marker','Save MarkerFile as');
p_write = strcat(pname,fname);
fid=fopen(p_write,'w');

header='##x,y,z,radius,shape,name,comment, color_r,color_g,color_b';
fprintf(fid,'%s\n\n',header);

% value = zeros(1,10);
for i = 1:numPts
   value1 = nodePos(i,2); % X
   value2 = nodePos(i,1); % Y
   value3 = nodePos(i,3); % Z
   value4 = 0 ; % radius
   value5 = 1; % shape
   value6 = strcat('landmark ',num2str(i)) ; %name
   value7 =' ' ; %comment
   value8 = round(1+254*rand(1)); %R
   value9 = round(1+254*rand(1)); %G
   value10 = round(1+254*rand(1)); %B
   toPrint = [num2str(value1) ',' num2str(value2) ',' num2str(value3) ','...
               num2str(value4) ',' num2str(value5) ',' value6 ','...
               value7 ',' num2str(value8) ',' num2str(value9) ','...
               num2str(value10) ','];
   fprintf(fid,'%s\n',toPrint);
end




end


function [nodes] = getMarkers(Img,r,c,d,pts)

nodes = zeros(pts,3);

index = 1;

for i = 1 : r
    for j = 1 : c
        for k = 1 : d
            if Img(i,j,k) == 1
               nodes(index,1) = i;
               nodes(index,2) = j;
               nodes(index,3) = k;
               index = index+1;
            end
            
        end
    end
end

end