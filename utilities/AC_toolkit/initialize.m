function [ pts ] = initialize(I)
%INITIALIZE the level curve

[nrow,ncol] = size(I);
[y_c,x_c] = ginput(1);

plot(y_c,x_c,'*r');
% 
% [y_r,x_r] = ginput(1);
% 
% plot(y_r,x_r,'*g');
% 
% rad = sqrt((x_c-x_r).^2 + (y_c-y_r).^2);
% 
% theta = 0:res:360;
% 
% pts.x = x_c + rad*cosd(theta);
% pts.y = y_c + rad*sind(theta);
% 
% pts.x = [pts.x pts.x(1)];
% pts.y = [pts.y pts.y(1)];
% 
% t1 = find(pts.x < 1);
% pts.x(t1)=1;
% t1 = find(pts.y < 1);
% pts.y(t1)=1;
% 
% t2 = find(pts.x >= nrow);
% pts.x(t2)=nrow;
% t2 = find(pts.y >= ncol);
% pts.y(t2)=ncol;


pts.x = x_c;
pts.y = y_c;
end

