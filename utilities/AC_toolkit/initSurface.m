function [surface] = initSurface(I,pts , rad)
%INITSURFACE Initialize the signed distance surface based on the contour.
%   surface < 0 --> INSIDE
%   surface > 0 --> OUTSIDE

surface=zeros(size(I));
if length(size(I)) == 3
    surface(pts.x- rad:pts.x+rad , pts.y- rad:pts.y+rad,pts.z - rad:pts.z+rad) = 1;
else
    pts.x = round(pts.x);
    pts.y = round(pts.y);
    surface(pts.x- rad:pts.x+rad , pts.y- rad:pts.y+rad) = 1;
end
surface = logical(surface);
% surface = bwdist(surface) - bwdist(~surface);

% [nrow,ncol,nslice]            = size(I);
% 
% [x,y,z]                       = meshgrid(1:ncol,1:nrow,1:nslice);
% 
% d = sqrt((x - pts.x).^2 + (y - pts.y).^2 + (z - pts.z).^2);
% surface                        = d - rad;             % Signed distance function, <0 for inside
    

end

