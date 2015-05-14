function [phi] = reInitPhi(Img,Sz, areaOpenSz, is_inside)
%REINITPHI Re-initialize the contour. Try to join nearby locations using
%the open-close filter. Img is the signed distance function image.

if is_inside == 0
    bwI = Img <= 0 ;
else
    bwI = Img >= 0;
end

phi             = zeros(size(bwI));
strel           = ones(Sz,Sz,Sz);


% bwI             = imopen(bwI,strel);
bwI             = imclose(bwI,strel);
% 
bwI             = bwareaopen(bwI,areaOpenSz,26);

phi             = double(bwdist(bwI) - bwdist(~bwI));     % The signed distance transform

masked_prev_phi = Img.*bwI ;    

if is_inside
   phi = -phi  ; 
end

end

