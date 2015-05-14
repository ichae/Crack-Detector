function [phi] = findZeroCrossing(Img,inside)
%FINDZEROCROSSING : Find the object contour and find the re_initialized
%surface by hole filling

[nrow ncol] = size(Img);

if inside == 0
    bwI = Img <= 0 ;
else
    bwI = Img >= 0;
end

se_close = strel('disk',3);
se_open = strel('disk',1);
area_frac = 0.8;
size_cc = sum(bwI(:));
area_openSz = round(area_frac*size_cc);

% % phi = zeros(size(bwI));
bwI = imclose(bwI,se_close);
% bwI = imopen(bwI,se_open);
% bwI = bwareaopen(bwI,area_openSz,8);
% CC = bwconncomp(bwI);
% largest_cc = max();
% bwI = bwareaopen(bwI,20,8);

% We want to add to the previous positive surface, not reiniitalize completely. But remove the negitive part. Let that be reinitialized

masked_prev_phi = Img.*bwI ;    
phi = double(bwdist(bwI) - bwdist(~bwI))+im2double(bwI)-0.5;     % The signed distance transform

if inside
%     phi = -phi + masked_prev_phi;
    phi = -phi;
end


end


