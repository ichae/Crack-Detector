function [g] = getImageTerm(Img,sigma)
%GETIMAGETERM Just 1/1+gradI


hsize = 8*sigma;
h = fspecial('gaussian',hsize,sigma);
filtImg = conv2(Img,h,'same');
% Img = edge(Img,'canny');
[gx gy] = gradient(filtImg);
grad = sqrt(gx.^2+gy.^2);

% grad = double(edge(Img,'canny'));

g = 1./(.1+grad.^2);
g = (g - min(g(:)))/(max(g(:))-min(g(:)));


end

