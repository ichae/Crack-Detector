function [] = viewMarkers( coord )
%VIEWMARKERS 

xval = coord(1,:);
yval = coord(2,:);
zval = coord(3,:);

numPts = length(xval(:));

[fname,pname]=uiputfile('*.marker','Save MarkerFile as');
p_write = strcat(pname,fname);
fid=fopen(p_write,'w');

header='##x,y,z,radius,shape,name,comment, color_r,color_g,color_b';
fprintf(fid,'%s\n\n',header);

% value = zeros(1,10);
for i = 1:numPts
   value1 = yval(i); % X
   value2 = xval(i); % Y
   value3 = zval(i); % Z
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

