function [A] = readSWC(fileName)
%READSWC returns the coordinate vector of all points from the SWC
%        the output file has the coordinate vector in each row
fid = fopen(fileName);
if (fid == -1)
    fprintf('File = %s Not found\n',fileName);
    coordSwc = -1 ;
    return;
end

start = fgetl(fid);
count = 0;
while(start(1)=='#')
    start = fgetl(fid);
%     next = fgetl(fid);
%     disp(start+1)
count = count+1;
    if numel(start)== 0 
        break
    elseif( start(1)~= '#' )
        break;
    end
    
%     ftell(fid)
%     disp(start)
end
% fprintf(fid,'%s',start);
% fseek(fid,-2,'cof')
frewind(fid);
for i = 1:count
    start = fgetl(fid);
end

% line1 = textscan(start,'%d %f %f %f %f %f %d');
% x1 = line1{1,3}(1);
% y1 = line1{1,4}(1);
% z1 = line1{1,5}(1);
% p1 = line1{1,7}(1);
% coord = [x1 y1 z1 ];
A = fscanf(fid,'%d %f %f %f %f %f %d',[7 inf]);

fclose(fid);

% Each column of the coordSwc matrix gives the coordinates of the points
% coordSwc = zeros(3,size(A,2));
% parent = zeros(1,size(A,2));
% coordSwc(1:3,:) = A(3:5,:);
% coordSwc = coordSwc' ;
% coordSwc = [coord;coordSwc];
% 
% parent = A(7,:);
% parent = parent';
% parent = [p1;parent];


%         

% node_num = A(1,:);
% node_num = node_num';
% coordSwc has been transposed....each row is the position vector


end

